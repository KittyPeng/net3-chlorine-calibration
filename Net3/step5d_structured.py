"""Step 5d: STRUCTURED within-zone heterogeneity — does the fit move off the arithmetic mean?

Location-consistent design (addresses the critique of the length-based five-group truth):
the truth keeps the three LOCATION zones, but WITHIN each zone the wall coefficient is
correlated with pipe length: longer pipes decay more strongly. The per-zone ARITHMETIC mean is
held exactly at the zone mean, while the LENGTH-weighted mean is shifted stronger, so the two
candidate "targets" of a three-zone homogeneous fit are separated and can be told apart.

Length weighting is only ONE candidate effective weighting, and it is NOT the hydraulically
effective one. The reaction weight of a pipe also depends on its flow, direction, diameter,
residence time and on how strongly the monitors see it, so "length-weighted" must never be read as
"residence-weighted" or as "the" effective mean — they are different quantities that happen to be
separated from the arithmetic mean in the same direction here. What this step establishes is a
DIRECTION of travel, not a target that is recovered: all three fitted coefficients move off their
arithmetic means in the length-weighted direction by a fraction of the gap whose median over 30 noise
realisations is near 1 (0.74-1.26 raw, every zone and every CORR level). A SINGLE realisation cannot
order the zones against each other: at CORR = 0.50 one draw reads 60/147/45% raw while the 30-draw
medians read 86/112/92%, and the 5-95 intervals overlap heavily, so which zone appears to pass the
proxy is a property of the draw and not of the network. Identifying the genuinely effective weighting needs a sensitivity/Jacobian-weighted mean, which this
step does not compute, so the claim available from it is bounded accordingly.

PRIMARY: formal censored likelihood. COMPARATOR: informal GLUE at the primary threshold. This
matters for the size of the effect as well as its direction: on the homogeneous baseline the
informal score's weighted mean already sits away from the truth, so part of any shift it reports is
its own artefact, and the baseline offset of each scheme is recorded here for exactly that reason.

Contrast with Step 5a/5c (random, uncorrelated jitter): there length-weighted ~ arithmetic, so the
two targets coincide and no separation is observable. Here the correlation is what separates them.

Reuses the frozen baseline candidate predictions (truth-independent).
"""
import os
import json
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import wntr

import wq_common as B

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures")
ZKEYS = ["old", "average", "new"]
ZONE_MEAN = {"old": B.KW_OLD_TRUE, "average": B.KW_AVG_TRUE, "new": B.KW_NEW_TRUE}
C_MIN = 0.2
CORR = 0.5                       # strength of the length-correlation (factor 1 + CORR*s, s in [-1,1])

cache = np.load(os.path.join(HERE, "baseline_cache", "baseline.npz"), allow_pickle=True)
C_all = cache["C_all"]
ALL_NODES = list(cache["all_nodes"])
mon_pos = list(cache["mon_pos"])
S = {"old": cache["S_old"], "average": cache["S_avg"], "new": cache["S_new"]}
C_all_mon = C_all[:, :, mon_pos]

wn0 = wntr.network.WaterNetworkModel(B.NET3_INP)
PIPE_LEN = {p: wn0.get_link(p).length for p in wn0.pipe_name_list}
zone_pipes = {z: [p for p in wn0.pipe_name_list if B.MATERIAL_ZONES[p] == z] for z in ZKEYS}

# ---- build the structured (length-correlated) truth ----
KW_PIPE = {}
for z in ZKEYS:
    ps = zone_pipes[z]
    lens = np.array([PIPE_LEN[p] for p in ps])
    order = np.argsort(lens)                     # shortest -> longest
    rank = np.empty(len(ps))
    rank[order] = np.arange(len(ps))
    s = 2.0 * (rank / (len(ps) - 1) - 0.5) if len(ps) > 1 else np.zeros(len(ps))  # in [-1,1], mean 0
    for p, sp in zip(ps, s):
        KW_PIPE[p] = ZONE_MEAN[z] * (1.0 + CORR * sp)   # longer pipe -> stronger decay


def hook(wn):
    for p in wn.pipe_name_list:
        wn.get_link(p).wall_coeff = B.per_day_to_per_second(KW_PIPE[p])


true_ref = {}
for z in ZKEYS:
    ps = zone_pipes[z]
    kws = np.array([KW_PIPE[p] for p in ps])
    lens = np.array([PIPE_LEN[p] for p in ps])
    true_ref[z] = {"n": len(ps), "arith": float(kws.mean()),
                   "lenwt": float((kws * lens).sum() / lens.sum()),
                   "min": float(kws.min()), "max": float(kws.max())}

truth_all = B.simulate_chlorine(B.KB_FIXED, 0.0, pre_run=hook, monitor_nodes=ALL_NODES).values
truth_mon = truth_all[:, mon_pos]
truth_post_mon = truth_mon[B.WARMUP_H:]

# grid library for structural residual + grid fit
kw_old_grid = np.round(np.linspace(-0.2, -1.5, 7), 3)
kw_avg_grid = np.round(np.linspace(-0.04, -0.2, 7), 3)
kw_new_grid = np.round(np.linspace(-0.005, -0.10, 7), 3)
GRID = {}
t0 = time.time()
for kwo in kw_old_grid:
    for kwa in kw_avg_grid:
        for kwn in kw_new_grid:
            GRID[(float(kwo), float(kwa), float(kwn))] = B.simulate_chlorine(
                B.KB_FIXED, 0.0, pre_run=B.make_kw_hook(kwo, kwa, kwn)).values[B.WARMUP_H:]
struct = min(np.sqrt(((sim - truth_post_mon) ** 2).mean()) for sim in GRID.values())
print(f"grid {len(GRID)} sims ({time.time()-t0:.1f}s); structural residual {struct:.4f}")

# noisy obs, then both weightings from the cached candidates
SCHEMES = [B.PRIMARY_WEIGHTING, "informal_glue"]
rng_n = np.random.default_rng(B.NOISE_SEED)
obs = np.clip(truth_mon + rng_n.normal(0, B.SIGMA_OBS, truth_mon.shape), 0, None)[B.WARMUP_H:]
RMSE = B.rmse_of(C_all_mon, obs)
beh = RMSE < B.RMSE_THR
wts = B.all_weightings(C_all_mon, obs, threshold=B.RMSE_THR, schemes=SCHEMES)
w = wts[B.PRIMARY_WEIGHTING][0]
GRIDFIT = None
best_rmse = np.inf
for kw, sim in GRID.items():
    r = np.sqrt(((sim - obs) ** 2).mean())
    if r < best_rmse:
        best_rmse, GRIDFIT = r, kw
gf = {"old": GRIDFIT[0], "average": GRIDFIT[1], "new": GRIDFIT[2]}

# The same seed on the HOMOGENEOUS truth, so each scheme's structural shift can be separated from
# the offset it already carries with no structural error at all.
obs_hom = cache["obs_glue"]
wts_hom = B.all_weightings(C_all_mon, obs_hom, threshold=B.RMSE_THR, schemes=SCHEMES)
fit = {s: {z: dict(zip(("mean", "sd"), B.weighted_mean_sd(wts[s][0], S[z]))) for z in ZKEYS}
       for s in SCHEMES}
base_offset = {s: {z: float(B.weighted_mean_sd(wts_hom[s][0], S[z])[0] - ZONE_MEAN[z])
                   for z in ZKEYS} for s in SCHEMES}

# risk under the primary rule
P_fit = B.pbar_field(C_all, C_MIN, w)
P_true = B.pbar_field(truth_all[B.WARMUP_H:], C_MIN)
rank_fit = [str(ALL_NODES[i]) for i in np.argsort(P_fit)[::-1][:6]]
rank_true = [str(ALL_NODES[i]) for i in np.argsort(P_true)[::-1][:6]]

# ---- dose-response over the correlation strength, with 30 noise realisations each ----
# The truth here is deterministic (length rank order), so there is no heterogeneity FIELD to average
# over as in Step 5c; the single-realisation risk is the noise draw instead. And if the effect is
# real it must vanish at CORR = 0, which is the homogeneous control, and grow with CORR. Both are
# tested together: dose-response plus a control is much stronger evidence than one point.
CORR_SWEEP = [0.0, 0.25, 0.50, 0.75]
N_NOISE_5D = 30
print(f"\n=== correlation dose-response ({N_NOISE_5D} noise realisations per level) ===")
print(f"{'CORR':>5} | " + " | ".join(f"{z} shift-frac med [5,95]" for z in ZKEYS))
dose = []
for corr in CORR_SWEEP:
    kwp, ref = {}, {}
    for z in ZKEYS:
        ps = zone_pipes[z]
        lens = np.array([PIPE_LEN[p] for p in ps])
        rk = np.empty(len(ps))
        rk[np.argsort(lens)] = np.arange(len(ps))
        sv = 2.0 * (rk / (len(ps) - 1) - 0.5) if len(ps) > 1 else np.zeros(len(ps))
        for p, sp in zip(ps, sv):
            kwp[p] = ZONE_MEAN[z] * (1.0 + corr * sp)
        kk = np.array([kwp[p] for p in ps])
        ref[z] = {"arith": float(kk.mean()),
                  "lenwt": float((kk * lens).sum() / lens.sum())}

    def hook_c(wn, _kw=kwp):
        for p in wn.pipe_name_list:
            wn.get_link(p).wall_coeff = B.per_day_to_per_second(_kw[p])

    tm = B.simulate_chlorine(B.KB_FIXED, 0.0, pre_run=hook_c,
                             monitor_nodes=ALL_NODES).values[:, mon_pos]
    fr = {s: {z: [] for z in ZKEYS} for s in SCHEMES}
    for seed in range(42, 42 + N_NOISE_5D):
        rr = np.random.default_rng(seed)
        ob = np.clip(tm + rr.normal(0, B.SIGMA_OBS, tm.shape), 0, None)[B.WARMUP_H:]
        wc = B.all_weightings(C_all_mon, ob, threshold=B.RMSE_THR, schemes=SCHEMES)
        for s in SCHEMES:
            if wc[s][0] is None:
                continue
            for z in ZKEYS:
                m = float(np.sum(wc[s][0] * S[z]))
                gap = ref[z]["lenwt"] - ref[z]["arith"]
                # fraction of the arithmetic -> length-weighted gap travelled; undefined at CORR = 0,
                # where the two targets coincide, so that row is the control instead
                fr[s][z].append((m - ref[z]["arith"]) / gap if abs(gap) > 1e-9 else np.nan)
    row = {"corr": corr, "gap": {z: ref[z]["lenwt"] - ref[z]["arith"] for z in ZKEYS},
           "by_scheme": {}}
    line = f"{corr:>5.2f} |"
    for s in SCHEMES:
        row["by_scheme"][s] = {}
        for z in ZKEYS:
            a = np.array(fr[s][z], dtype=float)
            if np.all(np.isnan(a)):
                row["by_scheme"][s][z] = {"shift_frac_med": None,
                                          "note": "targets coincide at CORR = 0"}
                if s == B.PRIMARY_WEIGHTING:
                    line += "        control        |"
                continue
            row["by_scheme"][s][z] = {"shift_frac_med": float(np.median(a)),
                                      "shift_frac_5_95": [float(np.percentile(a, 5)),
                                                          float(np.percentile(a, 95))]}
            if s == B.PRIMARY_WEIGHTING:
                line += (f" {np.median(a):+6.2f} [{np.percentile(a, 5):+.2f},"
                         f"{np.percentile(a, 95):+.2f}] |")
    dose.append(row)
    print(line)
print("shift-frac = (fitted mean − arithmetic) / (length-weighted − arithmetic). 0 means the fit sits")
print("on the arithmetic mean, 1 on the length-weighted value. At CORR = 0 the two coincide, so that")
print("row is the homogeneous control and no fraction is defined.")

GRID_AXES = {"old": kw_old_grid, "average": kw_avg_grid, "new": kw_new_grid}
report = {**B.weighting_provenance(comparators=["informal_glue"]),
          "design": "within-zone length-correlated heterogeneity (CORR=%.2f)" % CORR,
          "target_semantics": "the length-weighted value is an illustrative DIRECTIONAL "
                              "comparator, not an estimand or a bound; it is NOT the "
                              "residence-weighted or hydraulically effective coefficient. The "
                              "single-realisation gap fractions below order the zones only for "
                              "THIS noise draw; the 30-realisation dose-response is what supports "
                              "a median fraction near 1, and its 5-95 intervals overlap across "
                              "zones, so only the DIRECTION is established here",
          "informal_threshold": B.RMSE_THR, "structural_residual": float(struct),
          "behavioural_informal": int(beh.sum()), "rmse_min": float(RMSE.min()),
          "homogeneous_baseline_offset": base_offset,
          # the grid resolution bounds how close any grid fit can get, so it is recorded rather
          # than left to be recomputed by hand when the fit is discussed
          "grid": {z: {"nodes": [float(v) for v in ax],
                       "step": float(abs(ax[1] - ax[0])),
                       "half_step": float(abs(ax[1] - ax[0]) / 2),
                       # distance from the truth to the nearest node: the quantity that decides
                       # whether "recovery to the nearest grid node" means anything
                       "nearest_node_to_truth": float(ax[np.argmin(np.abs(ax - ZONE_MEAN[z]))]),
                       "distance_to_truth": float(np.min(np.abs(ax - ZONE_MEAN[z])))}
                   for z, ax in GRID_AXES.items()},
          "grid_fit": GRIDFIT, "zones": {}, "rank_fit": rank_fit, "rank_true": rank_true,
          "correlation_dose_response": {"corr_levels": CORR_SWEEP, "n_noise": N_NOISE_5D,
                                        "rows": dose}}
print("\n=== Step 5d: structured (length-correlated) within-zone heterogeneity ===")
print(f"structural residual {struct:.4f}; min RMSE over the candidate library {RMSE.min():.4f}")
for s in SCHEMES:
    tag = "PRIMARY" if s == B.PRIMARY_WEIGHTING else "comparator"
    print(f"\n--- {s} ({tag}) ---")
    print(f"{'zone':>8} | {'arith':>7} {'lenwt':>7} | {'fit mean±sd':>15} | {'gridfit':>8} | "
          f"{'shift':>8} {'hom. offset':>11} | {'lenwt-arith':>11} | frac of gap")
    for z in ZKEYS:
        tr = true_ref[z]
        m, sd = fit[s][z]["mean"], fit[s][z]["sd"]
        bias = m - tr["arith"]
        gap = tr["lenwt"] - tr["arith"]
        rec = {"mean": m, "sd": sd, "bias": bias,
               "homogeneous_baseline_offset": base_offset[s][z],
               # the shift net of what the scheme already does with no structural error at all
               "bias_net_of_baseline": bias - base_offset[s][z],
               "shift_frac_of_lenwt_gap": bias / gap if abs(gap) > 1e-12 else None,
               "shift_frac_net_of_baseline": ((bias - base_offset[s][z]) / gap
                                             if abs(gap) > 1e-12 else None)}
        if s == B.PRIMARY_WEIGHTING:
            report["zones"][z] = {**tr, "grid_fit": gf[z], "lenwt_minus_arith": gap,
                                  "by_scheme": {s: rec}}
        else:
            report["zones"][z]["by_scheme"][s] = rec
        print(f"{z:>8} | {tr['arith']:7.3f} {tr['lenwt']:7.3f} | {m:7.3f}±{sd:.3f} | "
              f"{gf[z]:8.3f} | {bias:>+8.3f} {base_offset[s][z]:>+11.3f} | {gap:>+11.3f} | "
              f"{100 * bias / gap:5.0f}%")
print(f"\n'hom. offset' is the same scheme's departure from the truth on the HOMOGENEOUS baseline:")
print("subtract it before reading any shift as structural.")
print(f"risk (primary) {rank_fit}\nrisk TRUE      {rank_true}")


def _jsafe(o):
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.integer):
        return int(o)
    return str(o)


with open(os.path.join(HERE, "baseline_cache", "step5d_structured.json"), "w") as f:
    json.dump(report, f, indent=2, default=_jsafe)

# figure: kw vs pipe length per zone + arithmetic / length-weighted / fitted lines
fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
for ax, z in zip(axes, ZKEYS):
    ps = zone_pipes[z]
    lens = np.array([PIPE_LEN[p] for p in ps])
    kws = np.array([KW_PIPE[p] for p in ps])
    m_prim = fit[B.PRIMARY_WEIGHTING][z]["mean"]
    ax.scatter(lens, kws, s=18, color="0.5", label="true pipe k_w")
    ax.axhline(true_ref[z]["arith"], color="black", lw=2, label="true arith. mean")
    ax.axhline(true_ref[z]["lenwt"], color="green", lw=2, ls="-.",
               label="length-weighted proxy (not residence-weighted)")
    ax.axhline(m_prim, color="steelblue", lw=2, ls="--",
               label="posterior mean, formal censored")
    ax.axhline(fit["informal_glue"][z]["mean"], color="0.35", lw=1.2, ls=":",
               label="posterior mean, informal GLUE")
    frac = (m_prim - true_ref[z]["arith"]) / (true_ref[z]["lenwt"] - true_ref[z]["arith"])
    # the defensible number is net of the homogeneous control; the raw one is kept beside it
    # so the size of that correction stays visible
    frac_net = report["zones"][z]["by_scheme"][B.PRIMARY_WEIGHTING][
        "shift_frac_net_of_baseline"]
    ax.set_xlabel("pipe length")
    ax.set_ylabel("k_w (m/day)")
    ax.set_title(f"{z} zone\nnet of control {100 * frac_net:.0f}% of the arith. -> "
                 f"length-weighted gap (raw {100 * frac:.0f}%)", fontsize=10)
    ax.grid(alpha=0.3)
axes[0].legend(fontsize=6.5)
fig.suptitle("Step 5d — structured (length-correlated) heterogeneity: the fitted coefficients move "
             "away from the arithmetic mean in the length-weighted direction\n"
             "(ONE noise draw — the zone ordering below is not resolved; see the 30-realisation "
             "dose-response, whose median fraction is near 1 in every zone)",
             y=1.06)
plt.tight_layout()
figpath = os.path.join(FIGDIR, "step5d_structured.png")
plt.savefig(figpath, dpi=130, bbox_inches="tight")
print("figure saved to", figpath)
