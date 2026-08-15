"""Frozen three-zone baseline configuration and WNTR/EPANET helpers.

This module fixes every choice of the main calibration experiment (originally in
03_pipeline_net3_02.ipynb) so that all follow-up experiments (threshold sensitivity,
displaced prior, structural error, noise sweep, Fisher information, sensor bias) run
against exactly the same synthetic truth, monitoring array, seeds and priors.

sigma = 0.1 mg/L is interpreted as ONE standard deviation of the Gaussian
observation error.

The network file is read from a FROZEN copy under models/net3_frozen/ rather than from the
installed WNTR package, so that upgrading WNTR cannot silently change the model. The SHA-256 is
checked on import; a mismatch is a hard error, because every cached result is tied to this file.
"""
import hashlib
import os
import numpy as np
import wntr

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NET3_INP = os.path.join(REPO_ROOT, "models", "net3_frozen", "Net3.inp")
NET3_INP_SHA256 = "ea3e825c4fef0b5cba47fb06301bc85253f18b6364dc96c44d9fb492c40faa52"


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def _check_net3_inp():
    if not os.path.exists(NET3_INP):
        raise RuntimeError(
            f"frozen network file missing: {NET3_INP}\n"
            "It is the provenance anchor for every cached result; restore it from git rather "
            "than falling back to the installed WNTR copy.")
    got = sha256_of(NET3_INP)
    if got != NET3_INP_SHA256:
        raise RuntimeError(
            f"frozen network file changed:\n  expected {NET3_INP_SHA256}\n  got      {got}\n"
            "Every cached result in baseline_cache/ was produced against the expected file.")
    return got


_check_net3_inp()

# ---- monitoring array: two nodes per zone (new 107/113 | old 15/145 | average 209/231) ----
MONITOR_NODES = ["107", "113", "15", "145", "209", "231"]

# ---- boundary / initial conditions ----
INLET_CHLORINE_MGL = 1.0
TANK_INIT_MGL = 0.5
SECONDS_PER_DAY = 24 * 3600

# ---- simulation timing ----
# Set by the Step 0 convergence test, not by convention. The chlorine field first becomes
# cyclostationary to within 0.005 mg/L at the monitors after a 120 h warm-up; 120 + 48 = 168 h is
# also exactly the model horizon, because pump 10 runs on absolute-time controls enumerated only to
# 159 h. The assessment window stays 48 h, so N_RESID is unchanged at 294 and every quantity derived
# from N (the behavioural threshold, the objective sampling sd, the profile cut-offs) carries over.
DURATION_H = 168
WARMUP_H = 120
HYDRAULIC_TIMESTEP_S = 3600
REPORT_TIMESTEP_S = 3600
QUALITY_TIMESTEP_S = 300

# ---- fixed / true parameters (m/day for wall, 1/day for bulk) ----
KB_FIXED = -0.5
KW_OLD_TRUE, KW_AVG_TRUE, KW_NEW_TRUE = -1.0, -0.1, -0.05

# Finite-difference steps for the Jacobian, one per coefficient rather than one shared absolute
# step. The three truths span a factor of 20, so a single 0.02 m/day is 2% of the old coefficient
# but 40% of the new one — for `new` that is a secant over a large arc, not a derivative. Each step
# is ~2-5% of its own truth and each is a member of that coefficient's convergence sweep in Step 7.
FD_STEP = {"old": 0.02, "average": 0.005, "new": 0.0025}
FD_STEP_KB = 0.05

# ---- GLUE configuration ----
PRIOR = {"old": (-1.5, -0.2), "avg": (-0.2, -0.04), "new": (-0.10, -0.005)}
# Sobol rather than pseudo-random: a scrambled Sobol set fills the 3-D prior box more evenly, and
# because every leading 2^k subset is itself balanced, the sampling-convergence diagnostic can
# compare 1024/2048/4096/8192 draws exactly rather than by ad-hoc thinning. 8192 = 2^13 is needed
# because the formal likelihood is far sharper than the informal score: at 2000 prior draws its
# effective sample size is only ~37.
N_MC_LOG2 = 13
N_MC = 2 ** N_MC_LOG2    # 8192
SIGMA_OBS = 0.1          # one standard deviation of the Gaussian observation error (mg/L)
# Behavioural thresholds. The draft used 0.12; the revised analysis uses a principled
# ~95% acceptance band tied to the sampling distribution of the RMSE objective (see below).
# Both thresholds belong to the INFORMAL GLUE comparator only. The primary analysis
# (PRIMARY_WEIGHTING below) is a formal likelihood and carries no threshold at all, because a hard
# acceptance cut-off is a feature of behavioural weighting, not of a likelihood.
RMSE_THR_DRAFT = 0.12    # the draft's loose threshold, kept so its configuration is reproducible
RMSE_THR = 0.107         # the defensible one: one-sided 95% band of the objective at sigma = 0.1
NOISE_SEED = 42          # seed for the baseline noisy observation set
SAMPLE_SEED = 0          # seed for the Sobol scramble of the prior draws

# number of residuals per candidate evaluation: monitors x post-warm-up hours
N_RESID = len(MONITOR_NODES) * ((DURATION_H - WARMUP_H) + 1)   # 6 x 49 = 294


def sobol_draws(box, n=None, seed=None):
    """Scrambled-Sobol draws from an arbitrary 3-D box {"old": (a,b), "avg": ..., "new": ...}.

    Displaced-prior designs must be sampled the same way as the baseline, otherwise a difference
    between them mixes a change of prior with a change of sampler.
    """
    from scipy.stats import qmc
    n = N_MC if n is None else n
    seed = SAMPLE_SEED if seed is None else seed
    m = int(np.log2(n))
    if 2 ** m != n:
        raise ValueError(f"Sobol designs must be a power of two, got {n}")
    u = qmc.Sobol(d=3, scramble=True, seed=seed).random_base2(m)     # (n, 3) in [0, 1)
    lo = np.array([box["old"][0], box["avg"][0], box["new"][0]])
    hi = np.array([box["old"][1], box["avg"][1], box["new"][1]])
    x = qmc.scale(u, lo, hi)
    return {"old": x[:, 0], "avg": x[:, 1], "new": x[:, 2]}


def prior_draws(n=None, seed=None):
    """Scrambled-Sobol draws from the baseline uniform prior box, as a dict of arrays.

    Leading 2^k subsets of a Sobol sequence are themselves balanced, so `prior_draws(n)[:m]` is a
    valid smaller design for m a power of two — that is what the sampling-convergence check uses.
    """
    return sobol_draws(PRIOR, n=n, seed=seed)


def threshold_for_sigma(sigma, z=1.645):
    """Principled behavioural threshold = one-sided acceptance band of the RMSE objective.

    The RMSE at the truth has mean ~ sigma and sampling sd ~ sigma / sqrt(2 N_RESID); the
    threshold accepts parameter sets whose RMSE is within z sampling-sd of the noise floor.
    Default z = 1.645 gives a ~95% one-sided band (0.107 mg/L at sigma = 0.1).

    This is a statement about the sampling distribution of the objective AT THE TRUTH (~95% of noise
    realisations there would be accepted), NOT a 95% credible interval for the parameters.
    """
    return float(sigma * (1.0 + z / (2.0 * N_RESID) ** 0.5))


def objective_sampling_sd(sigma=None, n_resid=None):
    """sd of the RMSE objective at the truth: SSE/sigma^2 ~ chi2(N) => sd(RMSE) ~ sigma/sqrt(2N)."""
    sigma = SIGMA_OBS if sigma is None else sigma
    n_resid = N_RESID if n_resid is None else n_resid
    return float(sigma / np.sqrt(2.0 * n_resid))


# ==================== likelihoods and weighting schemes ====================
# Three weightings are kept side by side deliberately. The formal ones are the primary analysis;
# the informal GLUE score is retained as a COMPARATOR because the draft used it and because the
# contrast is itself a result (Stedinger 2008; Mantovan & Todini 2006).
#
#   formal censored  PRIMARY. Gaussian density on uncensored points, Phi(-mu/sigma) on points the
#                    sensor floor clipped to zero. Nests the iid case when nothing is clipped.
#   formal iid       Gaussian on every point, treating a clipped zero as an exact measurement.
#   informal GLUE    exp(-0.5 (RMSE/sigma)^2) x 1[RMSE < threshold]. NOT a Gaussian likelihood: it
#                    drops the factor N, i.e. it is a Gaussian with sigma_eff = sigma*sqrt(N)
#                    = 1.7 mg/L at the baseline, 17x the sensor noise and larger than the inlet
#                    concentration, which is why it is nearly flat inside the behavioural set.

WEIGHTINGS = ("formal_censored", "formal_iid", "informal_glue")
PRIMARY_WEIGHTING = "formal_censored"
COMPARATOR_WEIGHTINGS = ("formal_iid", "informal_glue")


def weighting_provenance(comparators=COMPARATOR_WEIGHTINGS, **extra):
    """The stanza every result artifact must carry so a table can never be read under the wrong rule.

    Without it a JSON file records numbers whose weighting is only implied by the script that wrote
    them, which is how an informal-GLUE table ends up quoted under a formal heading.
    """
    return {"primary_weighting": PRIMARY_WEIGHTING, "comparators": list(comparators), **extra}


def _sum_over_data_axes(a):
    """Sum every axis except the leading candidate axis."""
    a = np.asarray(a)
    return a.reshape(a.shape[0], -1).sum(axis=1)


def log_gaussian(pred, obs, sigma=None):
    """Formal iid Gaussian log-likelihood per candidate, up to an additive constant.

    pred: (N, ...) candidate predictions; obs: (...) observations broadcast against them.
    """
    sigma = SIGMA_OBS if sigma is None else sigma
    resid = np.asarray(pred, dtype=np.float64) - np.asarray(obs, dtype=np.float64)[None]
    return -0.5 * _sum_over_data_axes((resid / sigma) ** 2)


def log_censored(pred, obs, sigma=None, floor=0.0):
    """Formal censored Gaussian log-likelihood per candidate, up to an additive constant.

    Observations equal to `floor` are treated as left-censored ("at most the floor") and contribute
    log Phi((floor - mu)/sigma) instead of a squared residual. Treating them as exact zeros instead
    is the naive alternative and is what log_gaussian does; Step 9 compares the two.
    """
    from scipy.special import log_ndtr
    sigma = SIGMA_OBS if sigma is None else sigma
    pred = np.asarray(pred, dtype=np.float64)
    obs = np.asarray(obs, dtype=np.float64)
    censored = (obs <= floor)[None]
    resid = pred - obs[None]
    dens = -0.5 * (resid / sigma) ** 2
    tail = log_ndtr((floor - pred) / sigma)
    return _sum_over_data_axes(np.where(censored, tail, dens))


def glue_score(rmse, sigma=None):
    """The informal GLUE score of the draft, as a LOG score so it composes with the formal ones.

    Returned as a log so every scheme can go through weights_from_loglik; the value is
    -0.5 (RMSE/sigma)^2, which is the formal Gaussian log-likelihood DIVIDED BY N.
    """
    sigma = SIGMA_OBS if sigma is None else sigma
    return -0.5 * (np.asarray(rmse, dtype=np.float64) / sigma) ** 2


def weights_from_loglik(loglik, mask=None):
    """Normalised weights plus the diagnostics needed to judge whether they are usable.

    mask: optional boolean acceptance mask (the behavioural filter of the informal scheme). The
    formal schemes take no threshold; a hard cut-off is part of the GLUE comparator, not of a
    likelihood.
    """
    loglik = np.asarray(loglik, dtype=np.float64)
    w = np.exp(loglik - loglik.max())
    if mask is not None:
        w = w * np.asarray(mask, dtype=np.float64)
    total = w.sum()
    if total <= 0:
        raise ValueError("all weights are zero: no candidate is acceptable")
    w = w / total
    nz = w > 0
    ess = float(1.0 / np.sum(w ** 2))
    entropy = float(-np.sum(w[nz] * np.log2(w[nz])))
    return w, {
        "n_candidates": int(w.size),
        "n_nonzero": int(nz.sum()),
        "ess": ess,
        "ess_frac": ess / w.size,
        "max_weight": float(w.max()),
        "entropy_bits": entropy,
        "entropy_bits_if_uniform": float(np.log2(w.size)),
    }


def rmse_of(pred, obs):
    """Root-mean-square residual per candidate, the objective the informal comparator scores."""
    resid = np.asarray(pred, dtype=np.float64) - np.asarray(obs, dtype=np.float64)[None]
    return np.sqrt((resid.reshape(resid.shape[0], -1) ** 2).mean(axis=1))


def all_weightings(pred, obs, sigma=None, threshold=None, schemes=WEIGHTINGS):
    """{scheme: (weights, diagnostics)} for one set of predictions and observations.

    Every step that compares the primary rule with the comparator goes through this function, so the
    two can never drift apart in the details (which sigma scales the score, whether the behavioural
    mask is applied, how the effective sample size is defined). A scheme whose weights are all zero
    — possible only for the thresholded comparator — comes back as (None, None) rather than raising,
    because a realisation with an empty behavioural set is a result about the comparator.
    """
    sigma = SIGMA_OBS if sigma is None else sigma
    threshold = threshold_for_sigma(sigma) if threshold is None else threshold
    out = {}
    rmse = rmse_of(pred, obs) if "informal_glue" in schemes else None
    for s in schemes:
        if s == "formal_censored":
            ll, mask = log_censored(pred, obs, sigma), None
        elif s == "formal_iid":
            ll, mask = log_gaussian(pred, obs, sigma), None
        elif s == "informal_glue":
            ll, mask = glue_score(rmse, sigma), rmse < threshold
        else:
            raise ValueError(f"unknown weighting scheme {s!r}")
        try:
            out[s] = weights_from_loglik(ll, mask)
        except ValueError:
            out[s] = (None, None)
    return out


def weighted_mean_sd(w, x):
    m = float(np.sum(w * x))
    return m, float(np.sqrt(np.sum(w * (np.asarray(x) - m) ** 2)))


def window_hours(C):
    """Length of the assessment window in hours: T reporting points bound T-1 one-hour intervals."""
    return np.asarray(C).shape[-2] - 1


def below_threshold_duration(C, c_min):
    """Hours below c_min, trapezoidal over the one-hour intervals. Time is the second-last axis.

    Section 2.8.1 defines the duration D as an integral over the T-1 one-hour intervals, not as a
    mean over the T reporting points: the two endpoints carry half weight and the divisor is the
    window length rather than the point count. Step 10 and paper_figs computed it this way while
    the perturbation scripts averaged the indicator over the 49 points, so the baseline field and
    the arms it was compared against were not the same quantity. Every risk field now goes through
    here, and tests/test_verification.py pins the two forms apart so the split cannot reopen.
    """
    return np.trapezoid((np.asarray(C) < c_min).astype(np.float64), dx=1.0, axis=-2)


def pbar_field(C, c_min, w=None):
    """Time-averaged below-threshold probability, P_bar = E[D] / T_window.

    With w, C is a candidate cube (N_MC, T, nodes) and the ensemble is weighted first. Without w,
    C is a single (T, nodes) trajectory, so the result is that trajectory's own time-average.
    """
    D = below_threshold_duration(C, c_min)
    if w is not None:
        D = np.tensordot(w, D, axes=(0, 0))
    return D / window_hours(C)


def weighted_quantile(x, w, q):
    """Weighted quantile with Hazen plotting positions (matches the convention used in Step 10)."""
    x = np.asarray(x, dtype=np.float64)
    o = np.argsort(x)
    xs, ws = x[o], np.asarray(w, dtype=np.float64)[o]
    cw = (np.cumsum(ws) - 0.5 * ws) / ws.sum()
    return np.interp(np.atleast_1d(q), cw, xs)

# ---- three contiguous zones by node coordinates ----
ZONE_Y_LOW = 10.0
ZONE_X_MID = 26.0
ZONE_RANK = {"new": 0, "average": 1, "old": 2}


def per_day_to_per_second(v):
    return v / SECONDS_PER_DAY


def zone_of(x, y):
    """Bottom band (y <= 10) = average; otherwise left (x <= 26) = new, right = old."""
    if y <= ZONE_Y_LOW:
        return "average"
    return "new" if x <= ZONE_X_MID else "old"


def assign_materials_zones(inp_file=NET3_INP):
    """Assign each pipe to a zone; a cross-zone pipe is given the weaker (newer) side."""
    wn = wntr.network.WaterNetworkModel(inp_file)
    mat = {}
    for p in wn.pipe_name_list:
        lk = wn.get_link(p)
        zs = zone_of(*wn.get_node(lk.start_node_name).coordinates)
        ze = zone_of(*wn.get_node(lk.end_node_name).coordinates)
        mat[p] = zs if zs == ze else (zs if ZONE_RANK[zs] <= ZONE_RANK[ze] else ze)
    return mat


MATERIAL_ZONES = assign_materials_zones()


def make_kw_hook(kw_old, kw_avg, kw_new, material=None):
    """Return a pre_run hook that sets per-pipe wall_coeff by zone."""
    mat = MATERIAL_ZONES if material is None else material
    kwz = {"old": kw_old, "average": kw_avg, "new": kw_new}

    def _hook(wn):
        for p in wn.pipe_name_list:
            wn.get_link(p).wall_coeff = per_day_to_per_second(kwz[mat[p]])
    return _hook


# ==================== concentration units ====================
# WNTR's Python API stores CONCENTRATION internally in kg/m^3, not in the units named by
# `options.quality.inpfile_units` — that field only governs how the .inp is read and written. So a
# value assigned to `initial_quality` is a kg/m^3 value, and `run_sim()` returns kg/m^3, whatever the
# .inp header says. Assigning 1.0 intending "1 mg/L" therefore simulates 1000 mg/L, and reading the
# result as mg/L repeats the error on the way out.
#
# This is the same factor the project already knew about in one place: the supervisor's notebook
# measured that a ZERO-ORDER global bulk coefficient needs x1000, because a zero-order coefficient
# carries a mass/volume dimension. It is the same kg/m^3 fact, and it applies to concentrations
# themselves as well. Under FIRST-order kinetics it is invisible — the field scales exactly with the
# source, and sigma and the risk threshold were scaled with it — which is why nothing broke and no
# known-answer test failed. The first-order coefficients are unaffected: k_b is 1/time and k_w is
# length/time, so neither carries a mass unit.
#
# Convention from here: everything OUTSIDE this module is in mg/L. The conversion happens at the
# model boundary, in build_model on the way in and in the two run helpers on the way out.
MG_L_PER_KG_M3 = 1000.0


def mgl_to_internal(v):
    """mg/L -> kg/m^3, the unit WNTR actually stores."""
    return v / MG_L_PER_KG_M3


def internal_to_mgl(v):
    """kg/m^3 -> mg/L."""
    return v * MG_L_PER_KG_M3


# A DELIBERATELY STRICT SOLVER TOLERANCE, in mg/L. Read the unit twice: unlike concentrations, this
# number is NOT converted. WNTR writes options.quality.tolerance into the .inp verbatim, and EPANET
# reads it in the file's own quality unit, so the .inp receives `TOLERANCE 1e-05` and EPANET treats
# it as 1e-5 MG/L -- a thousand times stricter than its 0.01 mg/L default. Verified by writing the
# .inp and reading the [OPTIONS] block back; the same probe shows a 0.001 kg/m^3 source coming out
# as `River 1.0`, i.e. concentrations converted, tolerance not.
#
# WHY THIS VALUE. The tolerance is an ABSOLUTE concentration below which EPANET treats two water
# parcels as identical, so what it means depends entirely on the concentration scale beside it. The
# superseded 1000x-too-high configuration ran the EPANET default 0.01 against a 1000 mg/L source, a
# ratio of 1e-5. Correcting the concentrations while leaving 0.01 in place would have left 0.01
# against 1 mg/L -- a ratio of 1e-2, i.e. a solver a thousand times coarser than every result
# published before it. 1e-5 mg/L holds the ratio fixed, which is what lets step15 compare the
# corrected and superseded runs as the SAME numerical experiment and attribute their agreement
# (2.3e-07) to the units alone. That is a choice made for comparability, not a unit conversion, and
# it is not the EPANET default -- state it as a chosen setting wherever the numerics are described.
QUALITY_TOLERANCE = 1e-5      # mg/L, as EPANET receives it; 1e-5 of a 1 mg/L source


def build_model(kb_per_day, kw_per_day, inp_file=NET3_INP, inlet_mgl=INLET_CHLORINE_MGL,
                tank_mgl=TANK_INIT_MGL, duration_hours=DURATION_H, bulk_order=1, wall_order=1,
                quality="CHEMICAL", pre_run=None):
    """The one place the baseline model is configured; both callers below use it.

    `inlet_mgl` and `tank_mgl` are in mg/L and are converted to WNTR's internal kg/m^3 here.
    """
    wn = wntr.network.WaterNetworkModel(inp_file)
    wn.options.time.duration = duration_hours * 3600
    wn.options.time.hydraulic_timestep = HYDRAULIC_TIMESTEP_S
    wn.options.time.report_timestep = REPORT_TIMESTEP_S
    wn.options.time.quality_timestep = QUALITY_TIMESTEP_S
    wn.options.quality.parameter = quality
    wn.options.quality.chemical_name = "Chlorine"
    wn.options.quality.inpfile_units = "mg/L"
    wn.options.quality.tolerance = QUALITY_TOLERANCE
    wn.options.reaction.bulk_order = bulk_order
    wn.options.reaction.wall_order = wall_order
    wn.options.reaction.bulk_coeff = per_day_to_per_second(kb_per_day)
    wn.options.reaction.wall_coeff = per_day_to_per_second(kw_per_day)
    # AGE carries no concentration: its "initial quality" is an age in seconds, so the conversion
    # must not be applied there.
    conv = mgl_to_internal if quality == "CHEMICAL" else (lambda v: v)
    for r in wn.reservoir_name_list:
        wn.get_node(r).initial_quality = conv(inlet_mgl)
    for t in wn.tank_name_list:
        wn.get_node(t).initial_quality = conv(tank_mgl)
    if pre_run is not None:
        pre_run(wn)
    return wn


def run_model(wn):
    """Run and return (results, hours index). Use when link or head results are needed too.

    For a CHEMICAL run the node quality is converted to mg/L; for AGE it is left alone (seconds),
    because the caller converts it to hours itself.
    """
    res = wntr.sim.EpanetSimulator(wn).run_sim()
    if wn.options.quality.parameter == "CHEMICAL":
        res.node["quality"] = internal_to_mgl(res.node["quality"])
    return res, np.asarray(res.node["quality"].index, dtype=float) / 3600.0


def simulate_chlorine(kb_per_day, kw_per_day, monitor_nodes=None, inp_file=NET3_INP,
                      inlet_mgl=INLET_CHLORINE_MGL, tank_mgl=TANK_INIT_MGL,
                      duration_hours=DURATION_H, bulk_order=1, wall_order=1, pre_run=None):
    """Single chlorine simulation. Returns a DataFrame (index = hours, cols = nodes)."""
    if monitor_nodes is None:
        monitor_nodes = MONITOR_NODES
    wn = build_model(kb_per_day, kw_per_day, inp_file=inp_file, inlet_mgl=inlet_mgl,
                     tank_mgl=tank_mgl, duration_hours=duration_hours, bulk_order=bulk_order,
                     wall_order=wall_order, quality="CHEMICAL", pre_run=pre_run)
    res = wntr.sim.EpanetSimulator(wn).run_sim()
    q = internal_to_mgl(res.node["quality"][monitor_nodes])     # kg/m^3 -> mg/L at the boundary
    q.index = q.index / 3600.0
    q.index.name = "hours"
    return q


def all_junctions(inp_file=NET3_INP):
    return wntr.network.WaterNetworkModel(inp_file).junction_name_list
