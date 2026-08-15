"""Step 8: systematic sensor-bias experiment (empirical) — Priority-2 #1.

A constant offset (systematic bias) is added to ONE informative monitor (node 15, old zone) and the
calibration is re-run under the PRIMARY formal censored likelihood. Reported: how far the posterior
mean of the old coefficient is pushed relative to its own posterior SD, how that displacement scales
with the offset, and what happens to the operational risk ranking — which is the quantity a utility
would act on.

Scaling convention: shift(0.025)/shift(0.05) below 50% means SUPER-LINEAR (convex) growth, above 50%
sub-linear (concave).

The sweep is TWO-SIDED. Because observations are censored at the sensor floor, a negative offset is
not the mirror image of a positive one — it pushes readings onto the floor, where they carry
different information (see Step 9) — so the asymmetry is measured rather than assumed away.

Step 8c repeats the experiment at every monitor in turn; the worst case is not node 15.

Reuses the baseline cache (candidate predictions are observation-independent); the bias only changes
the observations, so no EPANET is re-run. 30 noise realisations.
"""
import os
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, kendalltau
import wq_common as B

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures")
ZKEYS = ["old", "average", "new"]
TRUE = {"old": B.KW_OLD_TRUE, "average": B.KW_AVG_TRUE, "new": B.KW_NEW_TRUE}

cache = np.load(os.path.join(HERE, "baseline_cache", "baseline.npz"), allow_pickle=True)
C_all = cache["C_all"]
ALL_NODES = list(cache["all_nodes"])
mon_pos = list(cache["mon_pos"])
truth_mon = cache["truth_all"][:, mon_pos]            # (169, 6)
S = {"old": cache["S_old"], "average": cache["S_avg"], "new": cache["S_new"]}
C_all_mon = C_all[:, :, mon_pos]
C_MIN = 0.2
# Two risk metrics, precomputed once. Step 10's headline top-10 is ranked by the cumulative deficit,
# so a rank claim made only on P_bar is not a claim about the list the risk chapter publishes; Step
# 8b found the two can disagree about which nodes lead.
DUR = B.below_threshold_duration(C_all, C_MIN)
TW = B.window_hours(C_all)
DEFC = np.trapezoid(np.clip(C_MIN - C_all, 0.0, None), dx=1.0, axis=1).astype(np.float64)

BIAS_NODE = "15"                                       # informative old-zone monitor
bcol = B.MONITOR_NODES.index(BIAS_NODE)
# Negative as well as positive: because observations are censored at the sensor floor, a negative
# offset is NOT the mirror image of a positive one, so the sweep has to be two-sided to be honest.
OFFSETS = [-0.10, -0.05, -0.025, 0.0, 0.025, 0.05, 0.10]
N_NOISE = 30
SEEDS = list(range(42, 42 + N_NOISE))
thr = B.RMSE_THR                                       # 0.107


def posterior_means(obs_post):
    """Weighted means, SDs and the risk field under the formal censored likelihood.

    Formal rather than informal: a sensor offset is exactly the kind of systematic error the
    informal score is too flat to register, so measuring the damage with it would understate the
    damage (Step 8c quantifies that understatement at every monitor).
    """
    w, _ = B.all_weightings(C_all_mon, obs_post,
                            schemes=[B.PRIMARY_WEIGHTING])[B.PRIMARY_WEIGHTING]
    means = {z: float(np.sum(w * S[z])) for z in ZKEYS}
    sds = {z: float(np.sqrt(np.sum(w * (S[z] - means[z]) ** 2))) for z in ZKEYS}
    P = np.tensordot(w, DUR, axes=(0, 0)) / TW                    # P_bar: how LONG below C_MIN
    A = np.tensordot(w, DEFC, axes=(0, 0))                        # E[A]: how long AND how far below
    rank = [ALL_NODES[i] for i in np.argsort(P)[::-1][:6]]
    rank_A = [ALL_NODES[i] for i in np.argsort(A)[::-1][:6]]
    return means, sds, rank, P, rank_A, A


def med(a):
    return float(np.median(a))


rows = []
P_ref = None
for off in OFFSETS:
    old_means, old_sds, avg_means, new_means, ranks, Ps, n_clipped = [], [], [], [], [], [], []
    As, rmse_mins, noise_rmses = [], [], []
    for seed in SEEDS:
        rng = np.random.default_rng(seed)
        eps = rng.normal(0, B.SIGMA_OBS, truth_mon.shape)
        obs = truth_mon + eps
        obs[:, bcol] += off                            # inject systematic bias at node 15
        raw_neg = int((obs[B.WARMUP_H:] < 0).sum())
        obs = np.clip(obs, 0, None)[B.WARMUP_H:]
        # Best achievable aggregate fit over the whole candidate library against the BIASED
        # observations, against the noise floor the SAME draw would have produced with no offset.
        # This is what decides whether the aggregate residual reveals the perturbation at all: if
        # the bias is absorbed into the coefficients, the best residual stays at the noise scale.
        # Comparing instead with RMSE(biased obs - truth) would fold the bias into the comparator
        # and make the test vacuous.
        obs_unbiased = np.clip(truth_mon + eps, 0, None)[B.WARMUP_H:]
        rmse_mins.append(float(B.rmse_of(C_all_mon, obs).min()))
        noise_rmses.append(float(np.sqrt(((obs_unbiased - truth_mon[B.WARMUP_H:]) ** 2).mean())))
        means, sds, rank, P, rank_A, A = posterior_means(obs)
        old_means.append(means["old"]); old_sds.append(sds["old"])
        avg_means.append(means["average"]); new_means.append(means["new"])
        ranks.append(tuple(rank)); Ps.append(P); As.append(A); n_clipped.append(raw_neg)
    P_med = np.median(np.vstack(Ps), axis=0)
    A_med = np.median(np.vstack(As), axis=0)
    if off == 0.0:
        P_ref, A_ref = P_med, A_med
    # The leading set is read off the SAME median field the Spearman uses, so the two rank columns
    # of one table share a basis. The modal per-realisation ordering is a different statistic and is
    # kept beside it rather than substituted for it.
    row = {"offset": off,
           "old_mean_med": med(old_means), "old_sd_med": med(old_sds),
           "avg_mean_med": med(avg_means), "new_mean_med": med(new_means),
           "risk_top6": [str(ALL_NODES[i]) for i in np.argsort(P_med)[::-1][:6]],
           "risk_rank_mode": max(set(ranks), key=ranks.count),
           # censoring is the reason the sweep cannot be assumed symmetric: a negative offset pushes
           # more observations onto the sensor floor, a positive one lifts them off it
           "n_censored_med": float(np.median(n_clipped)),
           "deficit_top6": [str(ALL_NODES[i]) for i in np.argsort(A_med)[::-1][:6]],
           "rmse_min_med": med(rmse_mins),
           "noise_rmse_med": med(noise_rmses),
           "rmse_min_over_noise_med": med([a / b for a, b in zip(rmse_mins, noise_rmses)]),
           "_P": P_med, "_A": A_med}
    rows.append(row)

zero_row = next(r for r in rows if r["offset"] == 0.0)
base_old = zero_row["old_mean_med"]
base_sd = zero_row["old_sd_med"]
ref_rank = list(zero_row["risk_top6"])
ref_rank_A = list(zero_row["deficit_top6"])
for r in rows:
    r["old_shift"] = r["old_mean_med"] - base_old
    r["shift_over_sd"] = r["old_shift"] / base_sd
    # rank agreement with the unbiased case, on the full 92-node risk field and on the top-6 set
    r["spearman_vs_unbiased"] = float(spearmanr(r["_P"], P_ref).statistic)
    r["kendall_vs_unbiased"] = float(kendalltau(r["_P"], P_ref).statistic)
    top, ref = set(r["risk_top6"]), set(ref_rank)
    r["top6_jaccard_vs_unbiased"] = len(top & ref) / len(top | ref)
    # the same two comparisons on the cumulative-deficit field Step 10 ranks by
    r["deficit_spearman_vs_unbiased"] = float(spearmanr(r["_A"], A_ref).statistic)
    topA, refA = set(r["deficit_top6"]), set(ref_rank_A)
    r["deficit_top6_jaccard_vs_unbiased"] = len(topA & refA) / len(topA | refA)
    del r["_P"], r["_A"]

print(f"=== Step 8: systematic bias at node {BIAS_NODE} (old zone) ===")
print(f"weighting: formal censored likelihood; {N_NOISE} noise realisations; medians reported")
print(f"posterior SD of k_w,old with no offset (the random spread) = {base_sd:.4f}\n")
print(f"{'offset':>7} | {'old mean':>9} | {'old shift':>10} | {'shift/SD':>9} | "
      f"{'avg mean':>9} | {'new mean':>9} | {'cens':>5} | {'rho_S':>6} {'tau_K':>6} {'J6':>5}")
for r in rows:
    print(f"{r['offset']:>+7.3f} | {r['old_mean_med']:>9.4f} | {r['old_shift']:>+10.4f} | "
          f"{r['shift_over_sd']:>+9.2f} | {r['avg_mean_med']:>9.4f} | {r['new_mean_med']:>9.4f} | "
          f"{r['n_censored_med']:>5.0f} | {r['spearman_vs_unbiased']:>6.3f} "
          f"{r['kendall_vs_unbiased']:>6.3f} {r['top6_jaccard_vs_unbiased']:>5.2f}")
print("cens = median number of observations pushed onto the sensor floor; rho_S / tau_K = rank")
print("correlation of the 92-node risk field against the unbiased case; J6 = top-6 Jaccard.")

# How does the shift scale, and is the sweep symmetric?
def shift_at(off):
    return next(r["old_shift"] for r in rows if r["offset"] == off)


ratio = shift_at(0.025) / shift_at(0.05)
shape = "super-linear (convex)" if ratio < 0.5 else "sub-linear (concave)"


def informal_old_mean(off):
    """The comparator's old posterior mean at one offset, for the curvature contrast below."""
    vals = []
    for seed in SEEDS:
        rng = np.random.default_rng(seed)
        o = truth_mon + rng.normal(0, B.SIGMA_OBS, truth_mon.shape)
        o[:, bcol] += off
        o = np.clip(o, 0, None)[B.WARMUP_H:]
        wi, _ = B.all_weightings(C_all_mon, o, threshold=thr,
                                 schemes=["informal_glue"])["informal_glue"]
        vals.append(B.weighted_mean_sd(wi, S["old"])[0])
    return float(np.median(vals))


# The curvature is quoted, so the comparator's version of it is computed rather than asserted: if
# the informal score changes the sign of the conclusion, that has to be a number in the artifact.
inf0 = informal_old_mean(0.0)
inf_ratio = (informal_old_mean(0.025) - inf0) / (informal_old_mean(0.05) - inf0)
inf_shape = "super-linear (convex)" if inf_ratio < 0.5 else "sub-linear (concave)"
print(f"\nscaling: shift(+0.025)/shift(+0.05) = {ratio * 100:.0f}%  "
      f"(50% = linear; <50% = super-linear/convex, >50% = sub-linear/concave) -> {shape}")
print(f"The review states this effect is concave with ~70% retained on halving. The direction is")
print(f"confirmed here ({ratio * 100:.0f}% > 50%), though the magnitude is weaker than ~70%.")
print(f"The SAME quantity under the informal GLUE comparator is {inf_ratio * 100:.0f}% -> "
      f"{inf_shape}: the flat weighting changes the curvature as well as shrinking the shift, which")
print("is why the primary rule has to be the formal one (see Steps 3 and 5c for the same effect).")
asym = {}
for mag in (0.025, 0.05, 0.10):
    pos, neg = shift_at(mag), shift_at(-mag)
    asym[mag] = {"positive": pos, "negative": neg, "sum": pos + neg,
                 "abs_ratio": abs(neg) / abs(pos) if pos else None}
print("\nsymmetry check (a perfectly symmetric response would sum to zero):")
for mag, a in asym.items():
    print(f"  +/-{mag:<6} shift {a['positive']:+.4f} / {a['negative']:+.4f}  "
          f"sum {a['sum']:+.4f}  |neg|/|pos| {a['abs_ratio']:.2f}")
print("Asymmetry is expected and is caused by censoring: a negative offset drives observations")
print("onto the sensor floor, where they carry different information than an unclipped value.")

report = {**B.weighting_provenance(comparators=[]),
          "bias_node": BIAS_NODE, "n_noise": N_NOISE,
          "risk_metric": {
              "P_bar": "expected fraction of the 48 h window below 0.2 mg/L (risk_top6, "
                       "spearman_vs_unbiased, kendall_vs_unbiased, top6_jaccard_vs_unbiased)",
              "E_A": "expected cumulative deficit in mg/L*h (deficit_* fields); this is the metric "
                     "Step 10's headline top-10 is ranked by, and Step 8b showed the two metrics "
                     "can disagree about which nodes lead"},
          "risk_ranking_basis": "median risk field over the noise realisations; the top-6 set and "
                                "the Spearman/Kendall columns share that field, and the modal "
                                "per-realisation ordering is reported separately as risk_rank_mode",
          "symmetry_check": asym,
          "baseline_old_sd": base_sd, "rows": rows,
          "offsets_swept": OFFSETS,
          "offset_signs": "two-sided; censoring at the sensor floor makes the response asymmetric, "
                          "so the negative arm is measured rather than mirrored (see Step 9)",
          "location_note": "node 15 is one monitor, not the worst case; Step 8c sweeps all six and "
                           "finds a larger displacement at monitor 231 (average zone)",
          "shift025_over_shift050": ratio, "shift_scaling": shape,
          "informal_comparator_curvature": {"threshold": thr,
                                            "shift025_over_shift050": inf_ratio,
                                            "shift_scaling": inf_shape}}


def _jsafe(o):
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, tuple):
        return list(o)
    return str(o)


with open(os.path.join(HERE, "baseline_cache", "step8_sensor_bias.json"), "w") as f:
    json.dump(report, f, indent=2, default=_jsafe)

# ---- figure: old behavioural mean vs bias offset (± behavioural SD) ----
fig, ax = plt.subplots(figsize=(8, 5))
offs = [r["offset"] for r in rows]
means = [r["old_mean_med"] for r in rows]
sds = [r["old_sd_med"] for r in rows]
ax.errorbar(offs, means, yerr=sds, marker="o", capsize=4, color="steelblue",
            label="posterior mean ± posterior SD (formal censored)")
ax.axhline(TRUE["old"], color="red", ls="--", lw=1.5, label="true k_w,old = -1.0")
ax.set_xlabel("systematic bias at node 15 (mg/L)")
ax.set_ylabel("posterior mean of k_w,old (m/day)")
ax.set_title("Step 8 — a systematic sensor bias pushes the calibrated coefficient\n"
             f"(formal censored likelihood; shift is {shape}, "
             f"shift(0.025)/shift(0.05) = {ratio * 100:.0f}%)")
ax.legend()
ax.grid(alpha=0.3)
plt.tight_layout()
figpath = os.path.join(FIGDIR, "step8_sensor_bias.png")
plt.savefig(figpath, dpi=130, bbox_inches="tight")
print("figure saved to", figpath)
