"""Step 12: operational temperature / ageing scenario projection of the calibrated ensemble.

Transplants the supervisor's WSP risk-scenario framework onto THIS project's three-zone calibration
(not the homogeneous Gaussian posterior of the enclosed notebook). The ensemble is weighted by the
PRIMARY formal censored likelihood; the informal GLUE score plays no part in any number here, so the
word GLUE does not describe this step's weights.

Pipeline
--------
  the calibrated kinetics ensemble, weighted by the formal censored likelihood (Step 1)
    -> Arrhenius temperature scaling of k_b and the three zonal k_w
    -> optional illustrative ageing-reactivity multipliers on the three zones
    -> network-wide chlorine prediction
    -> window-breach probability, time-averaged probability, duration, deficit
    -> likelihood x consequence risk bands, maps, risk register
    -> control-measure evaluation: heatwave source dosing 1.00 / 1.15 / 1.30 mg/L

Three sources of scenario uncertainty are propagated jointly, with COMMON RANDOM NUMBERS
(one draw per behavioural member, reused across every scenario and dose) so that scenario
differences are physical rather than Monte-Carlo noise:
  1. kinetic coefficients          -- the formal-likelihood-weighted ensemble itself
  2. activation energies           -- Ea_bulk ~ N(45, 8^2), Ea_wall ~ N(35, 10^2) kJ/mol
  3. water temperature actually reached -- dT ~ N(0, 1^2) degC added to the scenario mean

Two probability definitions are reported and kept DISTINCT (they are not interchangeable):
  P_min  = sum_i w_i * 1[ min_t C_i(t) < C_crit ]   window-breach probability
  P_bar  = E[D] / T_window                          time-averaged below-threshold
                                                     probability (the Step-10 quantity)
The assessment window is the post-warm-up record (48 one-hour intervals), so P_min is a 48-hour
window minimum -- systematically higher than a 24-hour "daily minimum" and therefore not
directly comparable with the supervisor's notebook figures.

Numerical bound: zonal k_w are clipped to CLIP_LO purely as a solver guard. CLIP_LO is set
far outside the physical range reached by any scenario, and the number of clipped members
is reported and asserted to be zero, so no reported mean is distorted by the bound.

Caveats for the thesis: T_ref = 12 degC and the ageing alphas are illustrative planning
assumptions, not Net3 asset records; calibrated coefficients are effective parameters at an
assumed reference regime; the product is a calibration-conditioned scenario projection, not
a sensor nowcast, and scenario maps cannot be verified against data.
"""
from __future__ import annotations

import os
import json
import time

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import truncnorm
import wntr
from wntr.metrics.hydraulic import average_expected_demand

import wq_common as B

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures")
CACHEDIR = os.path.join(HERE, "baseline_cache")
os.makedirs(FIGDIR, exist_ok=True)

# ---- operational / scenario constants (illustrative) ----
C_CRIT = 0.2
T_REF_C = 12.0
R_GAS = 8.314
EA_BULK_MEAN, EA_BULK_SD = 45000.0, 8000.0     # J/mol
EA_WALL_MEAN, EA_WALL_SD = 35000.0, 10000.0    # J/mol
T_SD_C = 1.0                                   # water-temperature uncertainty (degC)
EA_MIN = 5000.0                                # J/mol, physical lower bound for Ea draws
T_VALID_C = (8.0, 24.0)                        # stated validity range of the assessment
CLIP_LO = -8.0                                 # inert solver guard (see module docstring)
DOSES = [1.00, 1.15, 1.30]
DRAW_SEED = 12

# Illustrative ageing-reactivity stress on the three material/age zones. The baseline model
# ALREADY distinguishes new/average/old k_w, so ageing is applied as an escalation-only
# stress (all alpha >= 1); alpha_new < 1 would weaken already-weak new pipes and
# double-count the zone structure. These are NOT asset measurements for Net3.
ALPHA_SETS = {
    "mild":    {"new": 1.00, "average": 1.15, "old": 1.40},
    "central": {"new": 1.00, "average": 1.35, "old": 1.85},
    "severe":  {"new": 1.00, "average": 1.50, "old": 2.20},
}
ALPHA_ZONE = ALPHA_SETS["central"]
ALPHA_NONE = {"new": 1.0, "average": 1.0, "old": 1.0}

SCENARIO_DEF = {
    "A_baseline": {"label": "A. Baseline 12 °C", "T": 12.0, "alpha": ALPHA_NONE},
    "B_warm":     {"label": "B. Warm season 16 °C", "T": 16.0, "alpha": ALPHA_NONE},
    "C_heatwave": {"label": "C. Heatwave 20 °C", "T": 20.0, "alpha": ALPHA_NONE},
    "D_heat_age": {"label": "D. Heatwave 20 °C + ageing stress", "T": 20.0, "alpha": ALPHA_ZONE},
}

LIK_LABEL = {1: "rare", 2: "unlikely", 3: "possible", 4: "likely", 5: "almost certain"}
CONS_LABEL = {0: "non-consumer", 1: "minor", 2: "moderate", 3: "major"}
# Severity axis: how LONG a node is expected to spend below the threshold, not whether it goes below
# at all. The bands are absolute hours in the 48 h assessment window, PRE-DECLARED rather than taken
# as quantiles of this network's own results, because a severity scale that moves with the data
# cannot be compared across scenarios (its own bands would shift with the heatwave). E[D] is the
# axis rather than E[A] because hours-below has an operational meaning that a mg/L*h integral does
# not; E[A] is carried alongside as the DEPTH diagnostic (E[A]/E[D] = mean depth while below), since
# a long shallow excursion and a short deep one give the same E[D].
SEV_LABEL = {1: "negligible", 2: "brief", 3: "sustained", 4: "prolonged", 5: "persistent"}
SEV_EDGES_H = (1.0, 6.0, 12.0, 24.0)     # < 1 | 1-6 | 6-12 | 12-24 | >= 24 (half the window)
CONTROL = {
    "very high": "Review dosing/booster strategy; consider flushing or turnover improvement",
    "high": "Confirmatory sampling; review local operation and demand assumptions",
    "medium": "Scheduled sampling within the monitoring programme; watch trend",
    "low": "No action beyond routine monitoring",
    "not applicable": "Not applicable",
}
# GOVERNANCE RULE, pre-declared. The register carries two risk products -- breach probability x
# consequence and severity x consequence -- and says that neither is the correct one. A register
# still has to drive a single action, so "neither is correct" has to be a RULE and not an evasion:
# the action is taken on the HIGHER of the two bands, i.e. neither axis alone is allowed to reduce
# an action. Both per-axis control measures are emitted alongside it so the reason stays visible.
#
# The rule is not the same thing as choosing the breach axis, even where the two coincide. In this
# baseline field the severity band never exceeds the breach band, so the governing band equals the
# breach band at every junction -- but that is measured (and reported as such), not assumed: the
# bound E[D] <= T_window * P_min permits an inversion of one band, so the rule can bind elsewhere.
BAND_ORDER = {"not applicable": 0, "low": 1, "medium": 2, "high": 3, "very high": 4}
GOVERNANCE_RULE = ("act on the higher of the breach band and the severity band; neither axis alone "
                   "may reduce an action")

trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz


def arrhenius_factor(t_c, ea, t_ref_c=T_REF_C):
    t = np.asarray(t_c, dtype=float) + 273.15
    t_ref = float(t_ref_c) + 273.15
    return np.exp(-(np.asarray(ea, dtype=float) / R_GAS) * (1.0 / t - 1.0 / t_ref))


LIK_EDGES = (0.05, 0.20, 0.50, 0.80)     # < 0.05 | 0.05-0.20 | 0.20-0.50 | 0.50-0.80 | >= 0.80


def likelihood_band(p):
    """Likelihood score 1-5 from P_min, the probability of ANY breach in the window."""
    for i, edge in enumerate(LIK_EDGES):
        if float(p) < edge:
            return i + 1
    return 5


def severity_reachability(t_window):
    """What E[D] <= T_window * P_min does and does NOT guarantee about the two banded scores.

    The inequality is real: D_i <= T_window * 1[member i breaches], so taking weighted expectations
    bounds E[D] by T_window * P_min. It is tempting -- and wrong -- to read that as "the severity
    score can never exceed the likelihood score". The two axes are banded on DIFFERENT and
    unaligned scales, so the bound constrains the inversion without forbidding it: with the edges
    used here, P_min = 0.03 with E[D] = 1.2 h satisfies the bound (1.2 <= 1.44) yet scores
    likelihood 1 and severity 2.

    What the bound does buy is a limit on how large an inversion can be, computed here from the
    edges themselves rather than asserted, so it stays true if either scale is ever changed.
    """
    sup_p = list(LIK_EDGES) + [1.0]
    rows = []
    for L in range(1, len(sup_p) + 1):
        d_sup = t_window * sup_p[L - 1]
        # the supremum is not attained for L < 5 (strict <), so probe just below it
        s_max = severity_band(d_sup if L == len(sup_p) else d_sup - 1e-12)
        rows.append({"likelihood_band": L, "sup_P_min": sup_p[L - 1],
                     "sup_E_duration_h": round(float(d_sup), 4),
                     "max_reachable_severity_band": int(s_max),
                     "max_inversion": int(max(0, s_max - L))})
    return rows


def severity_band(d_hours):
    """Severity score 1-5 from the expected below-threshold duration E[D], in hours."""
    d = float(d_hours)
    for i, edge in enumerate(SEV_EDGES_H):
        if d < edge:
            return i + 1
    return 5


def risk_band(score):
    if score == 0:
        return "not applicable"
    if score <= 3:
        return "low"
    if score <= 6:
        return "medium"
    if score <= 9:
        return "high"
    return "very high"


def demand_L_s(wn, nodes):
    """Pattern-aware average expected demand in L/s.

    NOT the base value. Four Net3 junctions (15, 35, 123, 203) encode their demand as 1 GPM times a
    large pattern, so a base-only reading makes them 0.063 L/s each instead of 16.7 / 108.4 / 75.3 /
    284.6, and the network total 192.6 instead of 690.7 L/s. Since demand is the CONSEQUENCE axis of
    the risk register, that error propagates into the terciles, the bands and the control measures -
    node 15 alone moves from minor/medium to major/very high. average_expected_demand applies the
    pattern and the global demand multiplier.

    WNTR parses the .inp and stores every demand internally in SI base units (m^3/s) regardless of
    the file's own unit declaration, so x1000 converts to L/s. The frozen Net3.inp declares
    `Units GPM` -- the conversion here is from WNTR's internal m^3/s, NOT from the file's units.
    """
    return average_expected_demand(wn)[list(nodes)] * 1000.0


def consequence_from_demand(demand):
    served = demand[demand > 0]
    q1, q2 = served.quantile([1 / 3, 2 / 3])

    def band(d):
        if d <= 0:
            return 0
        return 1 if d <= q1 else (2 if d <= q2 else 3)

    return demand.apply(band), float(q1), float(q2)


def metrics_from_C(C, w, t_window):
    """C: (n_members, n_t, n_nodes) -> the four node-level risk metrics."""
    D = trapz((C < C_CRIT).astype(np.float64), dx=1.0, axis=1)          # hours below
    A = trapz(np.maximum(0.0, C_CRIT - C), dx=1.0, axis=1)              # mg/L*h deficit
    M = C.min(axis=1)                                                   # window minimum
    P_min = w @ (M < C_CRIT).astype(np.float64)
    Dbar = w @ D
    return {"P_min": P_min, "P_bar": Dbar / t_window, "Dbar": Dbar, "Abar": w @ A, "M": M}


def run_ensemble(idx, S_old, S_avg, S_new, ea_b, ea_w, dT, T_mean, alpha, inlet,
                 nodes, tank_mgl=None, duration_h=B.DURATION_H, warmup_h=B.WARMUP_H,
                 members=None, quiet=False):
    """Forward every behavioural member under one scenario / dose."""
    sel = np.arange(len(idx)) if members is None else members
    n = len(sel)
    n_t = (duration_h - warmup_h) + 1
    C = np.empty((n, n_t, len(nodes)), dtype=np.float32)
    kb_used = np.empty(n)
    kw_used = {"old": np.empty(n), "average": np.empty(n), "new": np.empty(n)}
    n_clipped = 0
    tank = B.TANK_INIT_MGL if tank_mgl is None else tank_mgl
    t0 = time.time()
    for k, m in enumerate(sel):
        i = idx[m]
        t_water = T_mean + dT[m]
        fb = float(arrhenius_factor(t_water, ea_b[m]))
        fw = float(arrhenius_factor(t_water, ea_w[m]))
        raw = {"old": S_old[i] * fw * alpha["old"],
               "average": S_avg[i] * fw * alpha["average"],
               "new": S_new[i] * fw * alpha["new"]}
        n_clipped += sum(1 for v in raw.values() if v < CLIP_LO)
        kwz = {z: float(np.clip(v, CLIP_LO, 0.0)) for z, v in raw.items()}
        kb = B.KB_FIXED * fb
        kb_used[k] = kb
        for z in kw_used:
            kw_used[z][k] = kwz[z]
        C[k] = B.simulate_chlorine(
            kb, 0.0,
            pre_run=B.make_kw_hook(kwz["old"], kwz["average"], kwz["new"]),
            monitor_nodes=nodes, inlet_mgl=inlet, tank_mgl=tank,
            duration_hours=duration_h,
        ).values[warmup_h:].astype(np.float32)
        if not quiet and ((k + 1) % 300 == 0 or (k + 1) == n):
            print(f"    {k + 1}/{n}  ({time.time() - t0:.1f}s)")
    return C, kb_used, kw_used, n_clipped


def draw_network_background(ax, wn):
    for p in wn.pipe_name_list:
        lk = wn.get_link(p)
        x0, y0 = wn.get_node(lk.start_node_name).coordinates
        x1, y1 = wn.get_node(lk.end_node_name).coordinates
        ax.plot([x0, x1], [y0, y1], color="0.85", lw=0.6, zorder=1)


def scatter_nodes(ax, wn, nodes, values, vmin, vmax, cmap="viridis"):
    xs = [wn.get_node(n).coordinates[0] for n in nodes]
    ys = [wn.get_node(n).coordinates[1] for n in nodes]
    return ax.scatter(xs, ys, c=values, s=30, cmap=cmap, vmin=vmin, vmax=vmax,
                      edgecolors="0.3", linewidths=0.3, zorder=2)


def mark_monitors(ax, wn):
    xs = [wn.get_node(n).coordinates[0] for n in B.MONITOR_NODES]
    ys = [wn.get_node(n).coordinates[1] for n in B.MONITOR_NODES]
    ax.scatter(xs, ys, s=90, facecolors="none", edgecolors="navy",
               linewidths=1.6, zorder=5, label="monitors")


# ============================== setup ==============================
cache = np.load(os.path.join(CACHEDIR, "baseline.npz"), allow_pickle=True)
S_old = cache["S_old"].astype(np.float64)
S_avg = cache["S_avg"].astype(np.float64)
S_new = cache["S_new"].astype(np.float64)
RMSE = cache["RMSE"].astype(np.float64)
C_all = cache["C_all"].astype(np.float64)
ALL_NODES = [str(n) for n in cache["all_nodes"]]

# Weighting: the formal censored likelihood is primary. It carries no threshold, so the "ensemble"
# is every draw with non-negligible weight rather than a behavioural set; members whose weight is
# below WEIGHT_FLOOR of the maximum are dropped only to keep the scenario runs affordable, and the
# discarded mass is reported so the truncation is auditable.
WEIGHT_FLOOR = 1e-6
w_full, w_diag = B.weights_from_loglik(cache["loglik_censored"])
keep = w_full >= WEIGHT_FLOOR * w_full.max()
idx = np.where(keep)[0]
w = w_full[idx] / w_full[idx].sum()
n_beh = len(idx)
w_mass_dropped = float(1.0 - w_full[idx].sum())
T_WINDOW = C_all.shape[1] - 1

print(f"kinetics ensemble: {n_beh}/{len(RMSE)} draws retained above a {WEIGHT_FLOOR:g} relative "
      f"weight floor; ESS {w_diag['ess']:.1f}; discarded weight mass {w_mass_dropped:.2e}")
print(f"assessment window: {C_all.shape[1]} points = {T_WINDOW} h "
      f"(t = {B.WARMUP_H}..{B.DURATION_H})")
print(f"T_ref = {T_REF_C} °C, water-T uncertainty SD = {T_SD_C} °C, clip guard {CLIP_LO} m/day")
print(f"ageing stress (central) = {ALPHA_ZONE}")

# Activation energies are drawn from normals TRUNCATED at EA_MIN rather than clipped, so no
# draw can be non-physical and no probability mass piles up on the bound.
rng = np.random.default_rng(DRAW_SEED)


def draw_truncated_ea(mean, sd, n):
    a = (EA_MIN - mean) / sd
    return truncnorm.rvs(a, np.inf, loc=mean, scale=sd, size=n, random_state=rng)


T_MEAN_LO = min(s["T"] for s in SCENARIO_DEF.values())
T_MEAN_HI = max(s["T"] for s in SCENARIO_DEF.values())
# dT bounds that keep EVERY scenario mean inside the stated validity range (here ±4 °C = ±4 SD)
DT_LO, DT_HI = T_VALID_C[0] - T_MEAN_LO, T_VALID_C[1] - T_MEAN_HI


def draw_truncated_dT(sd, n, lo, hi):
    """Water-temperature offset, truncated so no draw can leave the stated validity range.

    Out-of-range draws are RESAMPLED rather than the distribution being replaced, so a run in
    which nothing falls outside (the case at the current seed and ensemble size) is bit-identical
    to an untruncated normal. Truncation matters once the ensemble grows: at ±4 SD the expected
    number of rejected draws is ~0.07 for 1000 members but ~1.3 for 20000.
    """
    d = rng.normal(0.0, sd, n)
    bad = (d < lo) | (d > hi)
    n_resampled = int(bad.sum())
    while bad.any():
        d[bad] = rng.normal(0.0, sd, int(bad.sum()))
        bad = (d < lo) | (d > hi)
    return d, n_resampled


ea_b = draw_truncated_ea(EA_BULK_MEAN, EA_BULK_SD, n_beh)
ea_w = draw_truncated_ea(EA_WALL_MEAN, EA_WALL_SD, n_beh)
dT, n_dT_resampled = draw_truncated_dT(T_SD_C, n_beh, DT_LO, DT_HI)   # common random numbers

assert ea_b.min() >= EA_MIN and ea_w.min() >= EA_MIN, "non-physical activation energy drawn"
T_lo = T_MEAN_LO + float(dT.min())
T_hi = T_MEAN_HI + float(dT.max())
assert T_VALID_C[0] <= T_lo and T_hi <= T_VALID_C[1], (
    f"sampled water temperature {T_lo:.2f}–{T_hi:.2f} °C leaves the stated "
    f"{T_VALID_C[0]}–{T_VALID_C[1]} °C validity range")
print(f"Ea draws (kJ/mol): bulk {ea_b.min()/1000:.1f}–{ea_b.max()/1000:.1f}, "
      f"wall {ea_w.min()/1000:.1f}–{ea_w.max()/1000:.1f} (truncated at {EA_MIN/1000:.0f})")
print(f"sampled water temperature spans {T_lo:.2f}–{T_hi:.2f} °C, "
      f"inside the {T_VALID_C[0]}–{T_VALID_C[1]} °C validity range "
      f"(dT truncated to [{DT_LO:+.1f}, {DT_HI:+.1f}] °C, {n_dT_resampled} draw(s) resampled)")

wn0 = wntr.network.WaterNetworkModel(B.NET3_INP)
DEMAND = demand_L_s(wn0, ALL_NODES)
CONSEQUENCE, q1, q2 = consequence_from_demand(DEMAND)
DEM = DEMAND.values
DEM_TOT = float(DEM.sum())
print(f"consequence terciles (L/s): {q1:.2f}, {q2:.2f}; total demand {DEM_TOT:.1f} L/s "
      f"over {int((DEM > 0).sum())} consumer nodes\n")


def risk_bands_for(P):
    scores = np.array([likelihood_band(float(p)) * int(CONSEQUENCE[n])
                       for n, p in zip(ALL_NODES, P)])
    return np.array([risk_band(int(s)) for s in scores]), scores


def severity_bands_for(Dbar):
    """The severity product on the same consequence axis. Costs nothing per scenario: E[D] is
    already computed for every scenario, so banding it is arithmetic on an array in memory."""
    scores = np.array([severity_band(float(d)) * int(CONSEQUENCE[n])
                       for n, d in zip(ALL_NODES, Dbar)])
    return np.array([risk_band(int(s)) for s in scores]), scores


def governing_bands_for(breach_bands, sev_bands):
    """The pre-declared governance rule: neither axis alone may reduce an action."""
    return np.array([b if BAND_ORDER[b] >= BAND_ORDER[s] else s
                     for b, s in zip(breach_bands, sev_bands)])


# ============================== scenarios ==============================
t_all = time.time()
results = {}
clip_log = {}
for key, spec in SCENARIO_DEF.items():
    print(f"=== {spec['label']} ===")
    C, kb_used, kw_used, n_clip = run_ensemble(
        idx, S_old, S_avg, S_new, ea_b, ea_w, dT,
        T_mean=spec["T"], alpha=spec["alpha"], inlet=B.INLET_CHLORINE_MGL, nodes=ALL_NODES)
    met = metrics_from_C(C, w, T_WINDOW)
    results[key] = {**met, "label": spec["label"], "T": spec["T"], "C": C,
                    "kb_mean": float(w @ kb_used),
                    "kw_old_mean": float(w @ kw_used["old"]),
                    "kw_avg_mean": float(w @ kw_used["average"]),
                    "kw_new_mean": float(w @ kw_used["new"])}
    clip_log[key] = n_clip
    at = met["P_min"] > 0.5
    print(f"  weighted mean kb {results[key]['kb_mean']:.3f} /day | "
          f"kw_old {results[key]['kw_old_mean']:.3f} m/day | clipped {n_clip}")
    print(f"  P_min>0.5 nodes {int(at.sum())} | demand at risk {DEM[at].sum():.1f} L/s\n")

# ratio check: scenario D vs C wall coefficient must equal alpha_old once the guard is inert
ratio_DC = results["D_heat_age"]["kw_old_mean"] / results["C_heatwave"]["kw_old_mean"]
print(f"verification  mean kw_old(D)/mean kw_old(C) = {ratio_DC:.4f}  "
      f"(alpha_old = {ALPHA_ZONE['old']:.2f})")
assert sum(clip_log.values()) == 0, f"clip guard active: {clip_log}"
assert abs(ratio_DC - ALPHA_ZONE["old"]) < 1e-6, "ageing multiplier not applied cleanly"
print("clip guard inert in every scenario (0 clipped draws)\n")

# reference row at exactly T_ref (cached forward runs) for continuity with Steps 1-11
ref = metrics_from_C(C_all[idx], w, T_WINDOW)

# ---- ageing-stress sensitivity: is scenario D driven by alpha_old = 1.85? ----
print("=== ageing-stress sensitivity (mild / central / severe) ===")
alpha_rows = []
for name, alpha in ALPHA_SETS.items():
    if name == "central":
        met = results["D_heat_age"]
    else:
        C_a, _, kwa, n_clip = run_ensemble(idx, S_old, S_avg, S_new, ea_b, ea_w, dT,
                                           T_mean=20.0, alpha=alpha,
                                           inlet=B.INLET_CHLORINE_MGL, nodes=ALL_NODES,
                                           quiet=True)
        assert n_clip == 0, f"clip guard active in ageing sensitivity ({name})"
        met = metrics_from_C(C_a, w, T_WINDOW)
        del C_a
    P = met["P_min"]
    at = P > 0.5
    bands, _ = risk_bands_for(P)
    alpha_rows.append({
        "ageing_set": name, "alpha_avg": alpha["average"], "alpha_old": alpha["old"],
        "P_min_gt_0.5_nodes": int(at.sum()),
        "demand_at_risk_L_s": round(float(DEM[at].sum()), 1),
        "high_or_very_high": int(np.isin(bands, ["high", "very high"]).sum()),
        "indeterminate": int(((P > 0.05) & (P < 0.95)).sum()),
        "net_mean_E_duration_h": round(float(met["Dbar"].mean()), 3),
        "net_mean_E_deficit": round(float(met["Abar"].mean()), 3),
    })
    print(f"  {name:>7}: P_min>0.5 {int(at.sum()):3d} | demand {DEM[at].sum():5.1f} L/s "
          f"| high/v-high {int(np.isin(bands, ['high', 'very high']).sum()):2d} "
          f"| indeterminate {int(((P > 0.05) & (P < 0.95)).sum()):2d} "
          f"| mean E[A] {met['Abar'].mean():.3f}")
print()

# ============================== corrective dosing ==============================
# Inlet dose scales BOTH reservoir source quality and tank initial quality, so the whole
# source regime is raised consistently. Leaving tanks at their un-dosed initial value would
# confound the result with a fixed boundary condition.
print("=== corrective dosing under heatwave (control-measure evaluation) ===")
SUB_STRIDE = 8
sub = np.arange(0, n_beh, SUB_STRIDE)
w_sub = w[sub] / w[sub].sum()
dose_results = {}
dose_C_sub = {}          # short-horizon trajectories of the SAME subset, for the paired test
for dose in DOSES:
    print(f"  inlet {dose:.2f} mg/L (tank initial {B.TANK_INIT_MGL * dose:.2f} mg/L)")
    if abs(dose - 1.0) < 1e-12:
        met = {k: results["C_heatwave"][k] for k in ("P_min", "P_bar", "Dbar", "Abar", "M")}
        dose_C_sub[dose] = results["C_heatwave"]["C"][sub].copy()
    else:
        C_d, _, _, n_clip = run_ensemble(idx, S_old, S_avg, S_new, ea_b, ea_w, dT,
                                         T_mean=20.0, alpha=ALPHA_NONE, inlet=dose,
                                         nodes=ALL_NODES, tank_mgl=B.TANK_INIT_MGL * dose,
                                         quiet=True)
        assert n_clip == 0
        met = metrics_from_C(C_d, w, T_WINDOW)
        dose_C_sub[dose] = C_d[sub].copy()
        del C_d
    dose_results[dose] = met

# linearity known-answer check: with all sources scaled, first-order kinetics give C ∝ dose
m0 = int(np.argmax(w))
C1, _, _, _ = run_ensemble(idx, S_old, S_avg, S_new, ea_b, ea_w, dT, T_mean=20.0,
                           alpha=ALPHA_NONE, inlet=1.0, nodes=ALL_NODES,
                           tank_mgl=B.TANK_INIT_MGL, members=[m0], quiet=True)
C13, _, _, _ = run_ensemble(idx, S_old, S_avg, S_new, ea_b, ea_w, dT, T_mean=20.0,
                            alpha=ALPHA_NONE, inlet=1.3, nodes=ALL_NODES,
                            tank_mgl=B.TANK_INIT_MGL * 1.3, members=[m0], quiet=True)
lin_err = float(np.max(np.abs(C13 - 1.3 * C1)))
print(f"  linearity check max|C(1.3) - 1.3*C(1.0)| = {lin_err:.2e} mg/L "
      "(first-order kinetics scale with the source regime)")

# Paired warm-up test: IDENTICAL members, weights and scenario draws; only the simulation horizon
# changes, so any difference is attributable to warm-up length alone. The roles are now reversed
# relative to the earlier version of this script: the CURRENT configuration is the long one that
# Step 0 justified, and the comparator is the SHORT 72 h / 24 h setting the draft used. The test is
# therefore a measurement of what the draft's warm-up cost, not a check on the current one.
SHORT_DUR, SHORT_WARM = 72, 24
paired_rows = []
print(f"  paired warm-up test on {len(sub)} identical members "
      f"(draft {SHORT_DUR} h / {SHORT_WARM} h warm-up  vs  current "
      f"{B.DURATION_H} h / {B.WARMUP_H} h)")
for dose in DOSES:
    met_l = metrics_from_C(dose_C_sub[dose], w_sub, T_WINDOW)          # current configuration
    C_s, _, _, n_clip = run_ensemble(idx, S_old, S_avg, S_new, ea_b, ea_w, dT, T_mean=20.0,
                                     alpha=ALPHA_NONE, inlet=dose, nodes=ALL_NODES,
                                     tank_mgl=B.TANK_INIT_MGL * dose, duration_h=SHORT_DUR,
                                     warmup_h=SHORT_WARM, members=sub, quiet=True)
    assert n_clip == 0
    met_s = metrics_from_C(C_s, w_sub, SHORT_DUR - SHORT_WARM)
    del C_s
    at_s, at_l = met_s["P_min"] > 0.5, met_l["P_min"] > 0.5
    paired_rows.append({
        "inlet_dose_mgl": dose,
        "short_P_min_gt_0.5_nodes": int(at_s.sum()),
        "long_P_min_gt_0.5_nodes": int(at_l.sum()),
        "short_demand_at_risk_L_s": round(float(DEM[at_s].sum()), 1),
        "long_demand_at_risk_L_s": round(float(DEM[at_l].sum()), 1),
        "short_net_mean_E_duration_h": round(float(met_s["Dbar"].mean()), 3),
        "long_net_mean_E_duration_h": round(float(met_l["Dbar"].mean()), 3),
        "short_net_mean_E_deficit": round(float(met_s["Abar"].mean()), 4),
        "long_net_mean_E_deficit": round(float(met_l["Abar"].mean()), 4),
        "long_minus_short_rel_E_deficit": round(
            float(met_l["Abar"].mean() / met_s["Abar"].mean() - 1.0), 4),
    })
del dose_C_sub
print(f"\ntotal runtime {time.time() - t_all:.1f}s\n")

# ============================== tables ==============================
# Severity and governing bands are computed for EVERY scenario, not only the baseline. The absolute
# pre-declared duration edges were justified by being portable across scenarios; that justification
# is only worth anything if the portability is exercised, and E[D] is already in memory for all four,
# so it costs no simulation. It also removes an inconsistency: the register acts on the governing
# band, so escalation must be judged on the governing band too.
SCEN_BANDS = {}
for key, r in results.items():
    b_breach, _ = risk_bands_for(r["P_min"])
    b_sev, _ = severity_bands_for(r["Dbar"])
    SCEN_BANDS[key] = {"breach": b_breach, "severity": b_sev,
                       "governing": governing_bands_for(b_breach, b_sev),
                       "severity_score": np.array([severity_band(float(d)) for d in r["Dbar"]])}

summary_rows = []
for key, r in results.items():
    P = r["P_min"]
    bands = SCEN_BANDS[key]["breach"]
    sev_bands = SCEN_BANDS[key]["severity"]
    gov_bands = SCEN_BANDS[key]["governing"]
    sev_sc = SCEN_BANDS[key]["severity_score"]
    at = P > 0.5
    summary_rows.append({
        "scenario": r["label"],
        "mean_kb": round(r["kb_mean"], 3),
        "mean_kw_old": round(r["kw_old_mean"], 3),
        "P_min_gt_0.5_nodes": int(at.sum()),
        "demand_at_risk_L_s": round(float(DEM[at].sum()), 1),
        "pct_demand_at_risk": round(100 * float(DEM[at].sum()) / DEM_TOT, 1),
        "high_or_very_high": int(np.isin(bands, ["high", "very high"]).sum()),
        "high_or_very_high_severity": int(np.isin(sev_bands, ["high", "very high"]).sum()),
        "high_or_very_high_governing": int(np.isin(gov_bands, ["high", "very high"]).sum()),
        "severity_band_counts": {SEV_LABEL[k]: int((sev_sc[CONSEQUENCE.values > 0] == k).sum())
                                 for k in range(1, 6)},
        "n_persistent_consumers": int((sev_sc[CONSEQUENCE.values > 0] == 5).sum()),
        "demand_at_risk_severity_L_s": round(
            float(DEM[sev_sc >= 3].sum()), 1),      # sustained or worse
        "indeterminate": int(((P > 0.05) & (P < 0.95)).sum()),
        "net_mean_P_bar": round(float(r["P_bar"].mean()), 4),
        "net_mean_E_duration_h": round(float(r["Dbar"].mean()), 3),
        "net_mean_E_deficit": round(float(r["Abar"].mean()), 4),
    })
summary_df = pd.DataFrame(summary_rows)
print("Scenario risk summary (formal censored-likelihood ensemble; network means are unweighted "
      "arithmetic means over all 92 junctions):\n")
print(summary_df.to_string(index=False))
print(f"\nreference (T = T_ref exactly, cached): P_min>0.5 nodes "
      f"{int((ref['P_min'] > 0.5).sum())}, demand at risk "
      f"{DEM[ref['P_min'] > 0.5].sum():.1f} L/s, net-mean E[A] {ref['Abar'].mean():.4f}")

dose_rows = []
base_at = results["A_baseline"]["P_min"] > 0.5
base_at_risk = float(DEM[base_at].sum())
for dose, met in dose_results.items():
    at = met["P_min"] > 0.5
    dose_rows.append({
        "inlet_dose_mgl": dose,
        "P_min_gt_0.5_nodes": int(at.sum()),
        "demand_at_risk_L_s": round(float(DEM[at].sum()), 1),
        "pct_demand_at_risk": round(100 * float(DEM[at].sum()) / DEM_TOT, 1),
        "net_mean_P_min": round(float(met["P_min"].mean()), 4),
        "net_mean_E_duration_h": round(float(met["Dbar"].mean()), 3),
        "net_mean_E_deficit": round(float(met["Abar"].mean()), 4),
        "median_over_nodes_of_mean_window_min_mgl": round(float(np.median(w @ met["M"])), 3),
    })
dose_df = pd.DataFrame(dose_rows)
print(f"\nBaseline (12 °C, 1.0 mg/L) demand at risk: {base_at_risk:.1f} L/s")
print("Heatwave dosing evaluation (reservoir AND tank source quality scaled):\n")
print(dose_df.to_string(index=False))
restored = dose_df["demand_at_risk_L_s"].min() <= base_at_risk + 1e-9
print(f"\nsource dosing alone {'restores' if restored else 'does NOT restore'} the "
      "pre-heatwave demand-at-risk position (control-measure evaluation, not a "
      "dosing recommendation).")
print(f"\nPaired warm-up test — identical {len(sub)} members, weights and draws; "
      "only the horizon differs:\n")
print(pd.DataFrame(paired_rows).to_string(index=False))

# ---- escalation and ageing increment ----
P_A = results["A_baseline"]["P_min"]
P_C = results["C_heatwave"]["P_min"]
P_D = results["D_heat_age"]["P_min"]
dP_age = P_D - P_C
bands_A = SCEN_BANDS["A_baseline"]["breach"]
bands_D = SCEN_BANDS["D_heat_age"]["breach"]
gov_A = SCEN_BANDS["A_baseline"]["governing"]
gov_D = SCEN_BANDS["D_heat_age"]["governing"]
sev_A = SCEN_BANDS["A_baseline"]["severity"]
sev_D = SCEN_BANDS["D_heat_age"]["severity"]

# A breach-only escalation test is STRUCTURALLY BLIND to any node already at P_min ~ 1 in both
# scenarios: its likelihood band cannot change, so it can never be flagged however much its
# below-threshold duration grows. Those nodes are exactly the ones the severity axis exists for.
esc = np.where((gov_A != gov_D) & (DEM > 0))[0]
esc_breach_only = np.where((bands_A != bands_D) & (DEM > 0))[0]
blind = np.where((P_A > 0.999) & (P_D > 0.999) & (DEM > 0))[0]
esc = esc[np.argsort(-DEM[esc])]
esc_rows = [{
    "node": ALL_NODES[j],
    "P_min_baseline": round(float(P_A[j]), 3),
    "P_min_heatwave": round(float(P_C[j]), 3),
    "P_min_heat_ageing": round(float(P_D[j]), 3),
    "dP_ageing": round(float(dP_age[j]), 3),
    "demand_L_s": round(float(DEM[j]), 2),
    "band_current": gov_A[j],
    "band_heat_ageing": gov_D[j],
    "band_current_breach": bands_A[j],
    "band_heat_ageing_breach": bands_D[j],
    "band_current_severity": sev_A[j],
    "band_heat_ageing_severity": sev_D[j],
    "E_duration_h_A": round(float(results["A_baseline"]["Dbar"][j]), 2),
    "E_duration_h_D": round(float(results["D_heat_age"]["Dbar"][j]), 2),
    "escalates_on_breach_only": bool(bands_A[j] != bands_D[j]),
    "monitored": ALL_NODES[j] in B.MONITOR_NODES,
} for j in esc]
esc_df = pd.DataFrame(esc_rows)
esc_compare = {
    "test": "risk band changes between scenario A and scenario D",
    "on_governing_band": {
        "n_consumer_junctions": int(len(esc)),
        "demand_L_s": round(float(DEM[esc].sum()), 2),
        "n_unmonitored": int(sum(1 for j in esc if ALL_NODES[j] not in B.MONITOR_NODES)),
    },
    "on_breach_band_only": {
        "n_consumer_junctions": int(len(esc_breach_only)),
        "demand_L_s": round(float(DEM[esc_breach_only].sum()), 2),
        "n_unmonitored": int(sum(1 for j in esc_breach_only
                                 if ALL_NODES[j] not in B.MONITOR_NODES)),
    },
    "found_only_by_the_governing_test": [ALL_NODES[j] for j in esc
                                         if j not in set(esc_breach_only.tolist())],
    "the_governing_rule_does_NOT_remove_the_blind_spot": {
        "why": "max(breach, severity) SATURATES. Where the breach band is already 'very high' -- "
               "which is what P_min ~ 1 at a consuming junction gives -- no rise in the severity "
               "band can move the governing band, so a band-change escalation test stays blind "
               "however far E[D] grows. The rule fixes which axis drives an ACTION; it does not "
               "make a band comparison sensitive to change at the ceiling",
        "nodes": [{
            "node": ALL_NODES[j],
            "demand_L_s": round(float(DEM[j]), 2),
            "E_duration_h_A": round(float(results["A_baseline"]["Dbar"][j]), 2),
            "E_duration_h_D": round(float(results["D_heat_age"]["Dbar"][j]), 2),
            "delta_E_duration_h": round(float(results["D_heat_age"]["Dbar"][j]
                                              - results["A_baseline"]["Dbar"][j]), 2),
            "severity_band_A": sev_A[j],
            "severity_band_D": sev_D[j],
            "severity_band_changes": bool(sev_A[j] != sev_D[j]),
            "governing_band_A": gov_A[j],
            "governing_band_D": gov_D[j],
            "flagged_by_either_escalation_test": bool(gov_A[j] != gov_D[j]),
        } for j in sorted(blind, key=lambda k: -(results["D_heat_age"]["Dbar"][k]
                                                 - results["A_baseline"]["Dbar"][k]))],
        "n_changing_severity_band_but_flagged_by_neither_test": int(sum(
            1 for j in blind if sev_A[j] != sev_D[j] and gov_A[j] == gov_D[j])),
        "demand_L_s_of_those": round(float(sum(
            DEM[j] for j in blind if sev_A[j] != sev_D[j] and gov_A[j] == gov_D[j])), 2),
        "recommendation": "for junctions whose band cannot move, escalation has to be read on the "
                          "CONTINUOUS metric (delta E[D], delta E[A]), not on a band change. The "
                          "band-change columns are necessary but not sufficient",
    },
    "two_distinct_reasons_a_band_cannot_move": {
        "1_saturation_at_the_top": "a junction already at 'very high' cannot rise. This is what "
                                   "P_min ~ 1 at a moderate or major consumer produces, and it is "
                                   "why max(breach, severity) does not help there",
        "2_ceiling_of_the_multiplicative_matrix": "score = axis score (1-5) x consequence score "
                                                  "(0-3), so a MINOR consumer (consequence 1) tops "
                                                  "out at score 5 -> 'medium' however severe it "
                                                  "gets. The register cannot express 'this small "
                                                  "consumer has no chlorine for the whole window' "
                                                  "as anything above medium",
        "max_attainable_band_by_consequence_score": {
            str(c): risk_band(5 * c) for c in (0, 1, 2, 3)},
        "n_junctions_capped_below_high": int((CONSEQUENCE.values == 1).sum()),
        "demand_L_s_capped_below_high": round(
            float(DEM[CONSEQUENCE.values == 1].sum()), 2),
        "worked_example": {
            "node": "243",
            "E_duration_h_A": round(float(results["A_baseline"]["Dbar"][ALL_NODES.index("243")]), 2),
            "E_duration_h_D": round(float(results["D_heat_age"]["Dbar"][ALL_NODES.index("243")]), 2),
            "demand_L_s": round(float(DEM[ALL_NODES.index("243")]), 2),
            "consequence_score": int(CONSEQUENCE["243"]),
            "note": "goes to the FULL window below the threshold under scenario D and stays "
                    "'medium' in both, because its consequence score of 1 caps the product",
        },
    },
    "structural_blind_spot_of_the_breach_test": {
        "n_consumer_junctions_with_P_min_1_in_both_A_and_D": int(len(blind)),
        "demand_L_s": round(float(DEM[blind].sum()), 2),
        "why": "their likelihood band cannot change between the two scenarios, so a breach-only "
               "escalation test can never flag them however far their below-threshold duration "
               "grows; this is the conflation the severity axis exists to prevent, present in the "
               "escalation column",
        "max_E_duration_growth_h": round(float(max(
            (results["D_heat_age"]["Dbar"][j] - results["A_baseline"]["Dbar"][j] for j in blind),
            default=0.0)), 2),
    },
}
esc_demand = float(esc_df["demand_L_s"].sum()) if len(esc_df) else 0.0
print(f"\n{len(esc_df)} consumer junctions change risk band A -> D; demand {esc_demand:.1f} L/s "
      f"({100 * esc_demand / DEM_TOT:.0f}% of network)\n")
if len(esc_df):
    print(esc_df.head(12).to_string(index=False))

# ---- risk register ----
# P_A is a weighted indicator average, so floating-point summation can put it a few ULP outside
# [0, 1]; clip before P(1-P) or a certain node reports a negative-zero sampling priority.
unc = np.clip(P_A, 0.0, 1.0) * (1.0 - np.clip(P_A, 0.0, 1.0))
prio_raw = CONSEQUENCE.values * unc
prio = np.where(DEM > 0, prio_raw / (prio_raw.max() + 1e-12), 0.0)

# ---- severity axis (parallel to the breach-probability axis, NOT a replacement) ----
# The breach product answers "does this node go below the threshold at all, and who does it serve".
# The severity product answers "for how long, and who does it serve". They are different questions
# and they do not have to agree; where they disagree is itself a reportable result, so both scores
# are carried and the disagreement is counted rather than resolved by picking one.
DBAR_A = results["A_baseline"]["Dbar"]
ABAR_A = results["A_baseline"]["Abar"]
SEV_SCORE = np.array([severity_band(d) for d in DBAR_A])
SEV_RISK = SEV_SCORE * CONSEQUENCE.values
SEV_BANDS = np.array([risk_band(int(s)) for s in SEV_RISK])
# mean depth below the threshold while below it; undefined (0) where the node never goes below
DEPTH_A = np.where(DBAR_A > 1e-9, ABAR_A / np.maximum(DBAR_A, 1e-9), 0.0)

shift = np.array([BAND_ORDER[SEV_BANDS[j]] - BAND_ORDER[bands_A[j]] for j in range(len(ALL_NODES))])
GOV_BANDS = np.array([bands_A[j] if BAND_ORDER[bands_A[j]] >= BAND_ORDER[SEV_BANDS[j]]
                      else SEV_BANDS[j] for j in range(len(ALL_NODES))])

register = pd.DataFrame([{
    "node": n,
    "P_min_current": round(float(P_A[j]), 3),
    "P_min_heatwave": round(float(P_C[j]), 3),
    "P_min_heat_ageing": round(float(P_D[j]), 3),
    "P_bar_current": round(float(results["A_baseline"]["P_bar"][j]), 3),
    "E_duration_current_h": round(float(DBAR_A[j]), 2),
    "E_deficit_current_mgL_h": round(float(ABAR_A[j]), 3),
    "E_depth_while_below_mgL": round(float(DEPTH_A[j]), 4),
    "demand_L_s": round(float(DEM[j]), 2),
    "likelihood": LIK_LABEL[likelihood_band(float(P_A[j]))],
    "severity": SEV_LABEL[int(SEV_SCORE[j])],
    "consequence": CONS_LABEL[int(CONSEQUENCE[n])],
    "risk_score_breach": int(likelihood_band(float(P_A[j])) * int(CONSEQUENCE[n])),
    "risk_score_severity": int(SEV_RISK[j]),
    "risk_band_current": bands_A[j],
    "risk_band_severity": SEV_BANDS[j],
    "risk_band_governing": GOV_BANDS[j],
    "risk_score_governing": int(max(int(likelihood_band(float(P_A[j])) * int(CONSEQUENCE[n])),
                                    int(SEV_RISK[j]))),
    "band_shift_severity_minus_breach": int(shift[j]),
    "risk_band_heat_ageing": bands_D[j],
    "escalates_under_heat": bool(gov_A[j] != gov_D[j]),                 # the governing rule
    "escalates_under_heat_breach_only": bool(bands_A[j] != bands_D[j]),
    "monitored": n in B.MONITOR_NODES,
    "sampling_priority": round(float(prio[j]), 3),
    "control_measure": CONTROL.get(GOV_BANDS[j], "Not applicable"),          # the governing rule
    "control_measure_breach": CONTROL.get(bands_A[j], "Not applicable"),
    "control_measure_severity": CONTROL.get(SEV_BANDS[j], "Not applicable"),
} for j, n in enumerate(ALL_NODES)])
register = register.sort_values(["risk_score_governing", "risk_score_breach", "P_min_current"],
                                ascending=False)

# how much the two axes disagree, and in which direction; a node is only interesting here if it
# serves demand, since a zero-demand node scores 0 on both by construction
cons_mask = CONSEQUENCE.values > 0
sev_counts = {SEV_LABEL[k]: int((SEV_SCORE[cons_mask] == k).sum()) for k in range(1, 6)}
agree = int((shift[cons_mask] == 0).sum())
sev_up = int((shift[cons_mask] > 0).sum())
sev_dn = int((shift[cons_mask] < 0).sum())
mov = np.where(cons_mask & (shift != 0))[0]
mov = mov[np.argsort(-np.abs(shift[mov]))]
# Every column the log tabulates for these nodes is written here, including the two derived ones,
# so the table has a single source. The depth spread is the point: "brief" is a statement about
# duration and says nothing about how far below the threshold a node sits while it is under it.
severity_moves = [{
    "node": ALL_NODES[j],
    "P_min": round(float(P_A[j]), 3),
    "E_duration_h": round(float(DBAR_A[j]), 2),
    "E_deficit_mgL_h": round(float(ABAR_A[j]), 4),
    "E_depth_while_below_mgL": round(float(DEPTH_A[j]), 4),
    "mean_C_while_below_mgL": round(float(C_CRIT - DEPTH_A[j]), 4),
    "demand_L_s": round(float(DEM[j]), 2),
    "band_breach": bands_A[j],
    "band_severity": SEV_BANDS[j],
    "direction": "severity higher" if shift[j] > 0 else "severity lower",
} for j in mov[:12]]
_depths = [m["E_depth_while_below_mgL"] for m in severity_moves]
_defs = [m["E_deficit_mgL_h"] for m in severity_moves]
# How far each consumer junction is from scoring HIGHER on severity than on likelihood. This is the
# empirical counterpart of the reachability table: the ordering is not guaranteed, so the question
# is not "can it invert" (it can) but "how close does this field come".
LIK_SCORE = np.array([likelihood_band(float(p_)) for p_ in P_A])
sev_minus_lik = SEV_SCORE - LIK_SCORE
hours_to_inversion = []
for j in range(len(ALL_NODES)):
    if not cons_mask[j]:
        continue
    if sev_minus_lik[j] > 0:
        hours_to_inversion.append(0.0)
        continue
    need = next((e for e in SEV_EDGES_H
                 if e > DBAR_A[j] and severity_band(e) > LIK_SCORE[j]), None)
    hours_to_inversion.append(float(need - DBAR_A[j]) if need is not None else float("inf"))
finite_gaps = [g for g in hours_to_inversion if g != float("inf")]

severity_axis = {
    "axis": "expected below-threshold duration E[D] in the assessment window, banded on absolute "
            "pre-declared hours",
    "band_edges_h": list(SEV_EDGES_H),
    "scenario": "A_baseline (the same scenario the breach product is scored on)",
    "consumer_junctions": int(cons_mask.sum()),
    "severity_band_counts_consumers": sev_counts,
    "risk_band_agreement_consumers": {"same band": agree, "severity band higher": sev_up,
                                      "severity band lower": sev_dn},
    "nodes_that_move": severity_moves,
    "depth_spread_over_movers": {
        "min_depth_mgL": min(_depths) if _depths else None,
        "max_depth_mgL": max(_depths) if _depths else None,
        "ratio": round(max(_depths) / min(_depths), 2) if _depths and min(_depths) > 0 else None,
        "min_mean_C_while_below_mgL": round(C_CRIT - max(_depths), 4) if _depths else None,
        "max_mean_C_while_below_mgL": round(C_CRIT - min(_depths), 4) if _depths else None,
        "note": "the movers share a probability and a limited duration, NOT a depth. Summarising "
                "them as 'marginally below' flattens a factor-of-several spread and repeats, one "
                "level down, the conflation the severity axis exists to prevent",
    },
    "deficit_spread_within_one_severity_band": {
        "band": "prolonged",
        "members": [m["node"] for m in severity_moves if m["band_severity"] and
                    m["E_duration_h"] >= SEV_EDGES_H[2] and m["E_duration_h"] < SEV_EDGES_H[3]],
        "E_deficit_mgL_h": [m["E_deficit_mgL_h"] for m in severity_moves
                            if m["E_duration_h"] >= SEV_EDGES_H[2]
                            and m["E_duration_h"] < SEV_EDGES_H[3]],
        "note": "E[D] bands cannot separate depth; this is the blind spot the depth column exists "
                "for, measured inside a single band",
    },
    "n_nodes_that_move": int(len(mov)),
    "reading": "the two products are not interchangeable: where the severity band is LOWER the node "
               "breaches reliably but only briefly, and where it is HIGHER the breach is less "
               "certain but long when it happens. Neither ordering is the correct one — the pair is "
               "the result",
    "ordering_is_not_guaranteed": {
        "bound": "E[D] <= T_window * P_min, from D_i <= T_window * 1[member i breaches]",
        "bound_violations_over_all_junctions": int(
            (DBAR_A > T_WINDOW * np.clip(P_A, 0.0, 1.0) + 1e-9).sum()),
        "what_it_does_NOT_imply": "an ordering of the two banded SCORES. The scales are not "
                                  "aligned, so the bound constrains an inversion without "
                                  "forbidding it",
        "counterexample": {"P_min": 0.03, "E_duration_h": 1.2,
                           "bound_at_this_P_min_h": 48 * 0.03,
                           "satisfies_bound": True,
                           "likelihood_band": likelihood_band(0.03),
                           "severity_band": severity_band(1.2)},
        "what_it_DOES_imply": "an inversion of at most one band, and none once P_min >= the top "
                              "likelihood edge — see reachability below",
        "reachability": severity_reachability(T_WINDOW),
        "max_inversion_permitted_by_the_bound": max(
            r["max_inversion"] for r in severity_reachability(T_WINDOW)),
    },
    "governance_rule": {
        "rule": GOVERNANCE_RULE,
        "applies_to": "the control_measure column and the register's sort order",
        "why_not_pick_an_axis": "the register says neither product is the correct one; that has to "
                                "be a rule rather than a stance, or a single control_measure "
                                "silently reinstates one axis as the operational one",
        "columns": {"control_measure": "from risk_band_governing (the higher of the two)",
                    "control_measure_breach": "from risk_band_current, for visibility",
                    "control_measure_severity": "from risk_band_severity, for visibility"},
        "n_junctions_where_severity_governs": int((shift > 0).sum()),
        "n_consumer_junctions_where_severity_governs": int((shift[cons_mask] > 0).sum()),
        "currently_equivalent_to_the_breach_axis": bool((shift > 0).sum() == 0),
        "equivalence_is_measured_not_assumed": "in this baseline field the severity band never "
                                               "exceeds the breach band, so the governing band "
                                               "equals the breach band everywhere. The bound "
                                               "E[D] <= T_window * P_min permits an inversion of "
                                               "one band, so the rule is not vacuous in general",
    },
    "observed_ordering": {
        "note": "EMPIRICAL, not implied by the bound",
        "consumer_junctions": int(cons_mask.sum()),
        "n_severity_score_above_likelihood_score": int((sev_minus_lik[cons_mask] > 0).sum()),
        "n_equal": int((sev_minus_lik[cons_mask] == 0).sum()),
        "n_severity_score_below": int((sev_minus_lik[cons_mask] < 0).sum()),
        "min_hours_of_E_duration_to_first_inversion": (round(min(finite_gaps), 3)
                                                       if finite_gaps else None),
        "n_consumers_that_can_never_invert_at_their_P_min": int(
            sum(1 for g in hours_to_inversion if g == float("inf"))),
    },
}
print(f"\nseverity axis on E[D] (consumer junctions only, n={int(cons_mask.sum())}):")
print("  band counts   " + ", ".join(f"{k} {v}" for k, v in sev_counts.items()))
print(f"  vs breach product: {agree} same band, {sev_up} severity higher, {sev_dn} severity lower")
if severity_moves:
    print(pd.DataFrame(severity_moves).to_string(index=False))
reg_path = os.path.join(CACHEDIR, "step12_risk_register.csv")
register.to_csv(reg_path, index=False)
print(f"\nrisk register -> {reg_path} ({len(register)} rows)")

# ============================== figures ==============================
fig, axes = plt.subplots(2, 2, figsize=(12.5, 10))
for ax, key in zip(axes.ravel(), ["A_baseline", "B_warm", "C_heatwave", "D_heat_age"]):
    draw_network_background(ax, wn0)
    sc = scatter_nodes(ax, wn0, ALL_NODES, results[key]["P_min"], 0.0, 1.0, cmap="inferno")
    mark_monitors(ax, wn0)
    ax.set_title(results[key]["label"])
    ax.set_aspect("equal")
    ax.axis("off")
fig.colorbar(sc, ax=axes.ravel().tolist(), shrink=0.85,
             label=rf"$P_{{\min}}$: P(min C over {B.WARMUP_H}–{B.DURATION_H} h < {C_CRIT} mg/L)")
fig.suptitle("Step 12 — posterior-propagated window-breach probability under temperature and "
             "ageing-stress scenarios (formal censored likelihood)", y=0.98)
fig.savefig(os.path.join(FIGDIR, "step12_scenario_maps.png"), dpi=140, bbox_inches="tight")
plt.close(fig)

fig, ax = plt.subplots(figsize=(8.5, 7))
draw_network_background(ax, wn0)
lim = max(float(np.max(np.abs(dP_age))), 0.05)
sc = scatter_nodes(ax, wn0, ALL_NODES, dP_age, -lim, lim, cmap="coolwarm")
mark_monitors(ax, wn0)
for j in np.argsort(-dP_age)[:8]:
    x, y = wn0.get_node(ALL_NODES[j]).coordinates
    ax.annotate(ALL_NODES[j], (x, y), textcoords="offset points", xytext=(4, 4),
                fontsize=7, color="0.15")
ax.set_aspect("equal")
ax.axis("off")
fig.colorbar(sc, ax=ax, shrink=0.85, label=r"$\Delta P_{\min}$ = P(D) − P(C)")
ax.set_title("Ageing-reactivity increment at 20 °C\n"
             "(extra window-breach probability from the illustrative zone multipliers)")
ax.legend(loc="lower left", fontsize=8)
fig.savefig(os.path.join(FIGDIR, "step12_ageing_delta.png"), dpi=140, bbox_inches="tight")
plt.close(fig)

fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.6))
labs = [r["scenario"].split(". ", 1)[-1] for r in summary_rows]
x = np.arange(len(labs))
axes[0].bar(x - 0.18, [r["P_min_gt_0.5_nodes"] for r in summary_rows], 0.36,
            color="firebrick", label=r"nodes $P_{\min}>0.5$")
axes[0].set_ylabel(r"junctions with $P_{\min}>0.5$")
ax0b = axes[0].twinx()
ax0b.plot(x + 0.18, [r["pct_demand_at_risk"] for r in summary_rows], "o-",
          color="steelblue", lw=2)
ax0b.set_ylabel("% of network demand at risk")
axes[0].set_xticks(x)
axes[0].set_xticklabels(labs, rotation=15, ha="right")
axes[0].set_title("(a) Scenario escalation")
axes[0].grid(alpha=0.3, axis="y")

axes[1].bar([f"{d:.2f}" for d in dose_df["inlet_dose_mgl"]],
            dose_df["demand_at_risk_L_s"], color="darkorange")
axes[1].axhline(base_at_risk, color="steelblue", ls="--", lw=1.5,
                label=f"baseline ({base_at_risk:.1f} L/s)")
ax1b = axes[1].twinx()
ax1b.plot([f"{d:.2f}" for d in dose_df["inlet_dose_mgl"]],
          dose_df["net_mean_E_deficit"], "s-", color="0.25", lw=1.8,
          label="mean E[A]")
ax1b.set_ylabel("network-mean E[A] (mg/L·h)")
axes[1].set_xlabel("heatwave inlet dose (mg/L)")
axes[1].set_ylabel("demand at risk (L/s)")
axes[1].set_title("(b) Control-measure evaluation")
axes[1].legend(fontsize=8, loc="lower left")
axes[1].grid(alpha=0.3, axis="y")

a_lab = [r["ageing_set"] for r in alpha_rows]
axes[2].bar(a_lab, [r["P_min_gt_0.5_nodes"] for r in alpha_rows], color="seagreen")
axes[2].axhline(summary_rows[2]["P_min_gt_0.5_nodes"], color="crimson", ls="--", lw=1.5,
                label="heatwave, no ageing")
for k, r in enumerate(alpha_rows):
    axes[2].annotate(f"α_old={r['alpha_old']:.2f}", (k, r["P_min_gt_0.5_nodes"]),
                     textcoords="offset points", xytext=(0, 4), ha="center", fontsize=8)
axes[2].set_ylabel(r"junctions with $P_{\min}>0.5$")
axes[2].set_title("(c) Ageing-stress sensitivity")
axes[2].legend(fontsize=8)
axes[2].grid(alpha=0.3, axis="y")

fig.tight_layout()
fig.savefig(os.path.join(FIGDIR, "step12_summary.png"), dpi=140, bbox_inches="tight")
plt.close(fig)
print("figures -> figures/step12_{scenario_maps,ageing_delta,summary}.png")

# ============================== json ==============================
report = {
    **B.weighting_provenance(comparators=[]),
    "description": "Step 12 operational temperature/ageing scenario projection of the three-zone "
                   "posterior ensemble, weighted by the primary formal censored likelihood",
    "informal_threshold_reference_only": B.RMSE_THR,
    "n_retained_draws": int(n_beh),
    "n_total_draws": int(len(RMSE)),
    "weight_floor_relative": WEIGHT_FLOOR,
    # The point of a numerical truncation is that its cost can be checked, so the discarded mass
    # belongs in the artifact and not only on stdout: a reader has to be able to see that dropping
    # 5996 of 8192 draws cost a fraction of the total weight far below the resolution of any number
    # reported here, without re-running the step.
    "retained_weight_mass": float(1.0 - w_mass_dropped),
    "discarded_weight_mass": w_mass_dropped,
    "retention_rule": "draws whose formal censored relative weight exceeds WEIGHT_FLOOR; this is a "
                      "numerical truncation of a likelihood, NOT a behavioural acceptance threshold",
    "C_CRIT": C_CRIT,
    "assessment_window_h": [B.WARMUP_H, B.DURATION_H],
    "T_window_intervals": int(T_WINDOW),
    "definitions": {
        "P_min": f"sum_i w_i * 1[min_t C_i(t) < C_crit] over t = {B.WARMUP_H}..{B.DURATION_H} h "
                 f"({T_WINDOW}-hour window minimum, NOT a 24-hour daily minimum)",
        "P_bar": "E[D]/T_window, the Step-10 time-averaged below-threshold probability",
        "E_duration_h": f"weighted expected hours below C_crit, trapezoid over {T_WINDOW} intervals",
        "E_deficit": "weighted expected cumulative deficit, mg/L*h",
        "network_mean": "unweighted arithmetic mean over all 92 junctions",
        "indeterminate": "0.05 < P_min < 0.95",
        "demand_at_risk": "sum of average expected demand over junctions with P_min > 0.5 "
                          "(zero-demand nodes contribute 0)",
        "median_over_nodes_of_mean_window_min_mgl":
            "per junction take the likelihood-weighted mean of the per-member window minimum, "
            "then the median of those 92 node values (NOT a pooled member-node median)",
        "likelihood_bands_on_P_min": {
            "rare": "P_min < 0.05", "unlikely": "0.05 <= P_min < 0.20",
            "possible": "0.20 <= P_min < 0.50", "likely": "0.50 <= P_min < 0.80",
            "almost certain": "P_min >= 0.80"},
        "likelihood_scores": {"rare": 1, "unlikely": 2, "possible": 3, "likely": 4,
                              "almost certain": 5},
        "consequence_scores": {"non-consumer (d = 0)": 0,
                               "minor (0 < d <= tercile1)": 1,
                               "moderate (tercile1 < d <= tercile2)": 2,
                               "major (d > tercile2)": 3},
        "consequence_terciles": "1/3 and 2/3 quantiles of average expected demand over the 59 "
                                "junctions with non-zero demand. d is the PATTERN-AWARE average "
                                "expected demand (wntr.metrics.hydraulic.average_expected_demand), "
                                "not the base_demand field: four Net3 junctions encode a large "
                                "demand as 1 GPM times a large pattern and a base-only reading "
                                "puts them in the minor band. WNTR stores demand internally in "
                                "m^3/s, so x1000 -> L/s; the .inp itself declares GPM",
        "risk_score_breach": "likelihood score (on P_min) x consequence score, range 0-15. "
                             "Answers: does this node breach at all, and who does it serve",
        "risk_score_severity": "severity score (on E[D]) x consequence score, range 0-15. "
                               "Answers: for how long, and who does it serve. It is a PARALLEL "
                               "axis, not a replacement — the two products answer different "
                               "questions and are reported together with their disagreement",
        "severity_bands_on_E_duration_h": {
            "negligible": "E[D] < 1 h", "brief": "1 <= E[D] < 6 h",
            "sustained": "6 <= E[D] < 12 h", "prolonged": "12 <= E[D] < 24 h",
            "persistent": "E[D] >= 24 h (half the assessment window or more)"},
        "severity_scores": {"negligible": 1, "brief": 2, "sustained": 3, "prolonged": 4,
                            "persistent": 5},
        "severity_bands_are_absolute": "the edges are pre-declared hours, NOT quantiles of this "
                                       "network's own results, so the same scale applies across "
                                       "scenarios; a data-driven scale would move with the "
                                       "heatwave and could not show escalation",
        "E_depth_while_below_mgL": "E[A]/E[D], the mean depth below C_crit while below it; 0 "
                                   "where E[D] = 0. E[D] alone cannot separate a long shallow "
                                   "excursion from a short deep one, so the depth is carried "
                                   "with it",
        "risk_band_mapping": {"not applicable": "score = 0", "low": "1 <= score <= 3",
                              "medium": "4 <= score <= 6", "high": "7 <= score <= 9",
                              "very high": "score >= 10"},
        "sampling_priority": "consequence SCORE (0-3, not raw demand) x P_min(1 - P_min), "
                             "normalised by its maximum; forced to 0 at zero-demand nodes. This "
                             "ranks where MEASUREMENT would reduce uncertainty most, so a node "
                             "with P_min = 1 scores 0: it is confidently at risk, not low "
                             "priority. Intervention priority is the separate risk_score_breach "
                             "and risk_score_severity columns; the three must not be merged.",
    },
    "uncertainty_sources": {
        "kinetics": f"formal censored-likelihood weights over {n_beh} retained draws "
                    f"(ESS {w_diag['ess']:.1f}); informal GLUE is a comparator, not the primary",
        # the drawn ranges, not just the parameters: the log quotes them as evidence that no draw
        # is non-physical, which is only checkable if they are stored
        "Ea_bulk_J_per_mol": {"mean": EA_BULK_MEAN, "sd": EA_BULK_SD,
                              "drawn_min": float(ea_b.min()), "drawn_max": float(ea_b.max()),
                              "truncated_at": EA_MIN},
        "Ea_wall_J_per_mol": {"mean": EA_WALL_MEAN, "sd": EA_WALL_SD,
                              "drawn_min": float(ea_w.min()), "drawn_max": float(ea_w.max()),
                              "truncated_at": EA_MIN},
        "water_temperature_span_C": [T_lo, T_hi],
        "water_temperature_C": f"dT ~ N(0, {T_SD_C}^2) added to the scenario mean, truncated to "
                               f"[{DT_LO:+.1f}, {DT_HI:+.1f}] degC so every scenario stays inside "
                               f"the {T_VALID_C[0]}-{T_VALID_C[1]} degC validity range",
        "water_temperature_draws_resampled": n_dT_resampled,
        "common_random_numbers": True,
    },
    "T_ref_C": T_REF_C,
    "clip_guard_m_per_day": CLIP_LO,
    "clipped_draws_per_scenario": clip_log,
    "kw_old_ratio_D_over_C": round(float(ratio_DC), 6),
    "ageing_alphas_illustrative": ALPHA_SETS,
    "ageing_sensitivity": alpha_rows,
    "sigma_convention": "observation sigma = 0.1 mg/L is one standard deviation",
    "product_statement": (
        "calibration-conditioned scenario projection from the network model; not a "
        "sensor nowcast, not a spatial measurement, and not a statement that water is safe"
    ),
    "consequence_terciles_L_s": [q1, q2],
    "severity_axis": severity_axis,
    "scenario_summary": summary_rows,
    "reference_at_T_ref_exact": {
        "P_min_gt_0.5_nodes": int((ref["P_min"] > 0.5).sum()),
        "demand_at_risk_L_s": round(float(DEM[ref["P_min"] > 0.5].sum()), 1),
        "net_mean_E_deficit": round(float(ref["Abar"].mean()), 4),
    },
    "dosing_heatwave": dose_rows,
    "dosing_boundary": "reservoir source quality AND tank initial quality both scaled by dose",
    "dosing_linearity_check": {
        "max_abs_err_mgL": lin_err,
        "conditions": "holds under fixed hydraulics and demands, first-order bulk and wall "
                      "kinetics, and proportional scaling of ALL source and initial chlorine "
                      "concentrations; a dose factor r is then equivalent to evaluating the "
                      "unscaled field against C_crit / r",
    },
    "dosing_paired_warmup_test": {
        "design": "identical member subset, weights and scenario draws; only the simulation "
                  "horizon differs, so the contrast isolates warm-up length. 'long' is the CURRENT "
                  "configuration justified by Step 0; 'short' is the draft setting, so the rows "
                  "measure what the draft's warm-up cost rather than testing the current one",
        "short": {"duration_h": SHORT_DUR, "warmup_h": SHORT_WARM, "role": "draft comparator"},
        "long": {"duration_h": B.DURATION_H, "warmup_h": B.WARMUP_H, "role": "current baseline"},
        "n_members": int(len(sub)), "subset_stride": SUB_STRIDE, "rows": paired_rows,
    },
    "baseline_demand_at_risk_L_s": base_at_risk,
    "dosing_restores_baseline": bool(restored),
    "escalation_test": esc_compare,
    "n_escalating_consumer_nodes": int(len(esc_df)),
    "escalating_demand_L_s": esc_demand,
    "top_escalating": esc_rows[:15],
    "top_ageing_delta": [
        {"node": ALL_NODES[j], "dP_min": round(float(dP_age[j]), 3),
         "P_C": round(float(P_C[j]), 3), "P_D": round(float(P_D[j]), 3),
         "demand_L_s": round(float(DEM[j]), 2)}
        for j in np.argsort(-dP_age)[:10]
    ],
    "review_triggers": [
        "sensor QA failure / reference-check disagreement / sensor service",
        "hydraulic model or demand allocation revised",
        "source/treatment regime changes inlet chlorine",
        "water temperature moves outside the 8-24 degC scenario range",
        "forecast heat episode (switch to the scenario-D register; pre-authorise response)",
        "mains rehabilitation or burst altering the assumed age profile",
        "calibration record exceeds its approved age",
    ],
    "runtime_s": round(time.time() - t_all, 1),
}
with open(os.path.join(CACHEDIR, "step12_scenarios.json"), "w") as f:
    json.dump(report, f, indent=2)
print("json -> baseline_cache/step12_scenarios.json")
print("\nDONE Step 12")
