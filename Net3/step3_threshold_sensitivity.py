"""Step 3: behavioural-threshold sensitivity of the INFORMAL GLUE comparator.

This section is deliberately and entirely informal: the threshold is a feature of behavioural
weighting, the formal likelihood has none, so "how much does the answer depend on the threshold?" is
a question that can only be asked of the comparator. What it establishes is how sensitive the
draft's conclusions were to an analyst's choice — not what the data support, which is the primary
formal analysis in Steps 1 and 7.

Uses the frozen baseline cache only (no EPANET re-run). For each RMSE threshold it recomputes the
behavioural set, the normalised GLUE weights, the weighted mean/SD of the three coefficients, the
prior-to-behavioural width ratio, the 5-95% predictive band at a near-source and a downstream node,
and the top low-chlorine risk nodes.
"""
import os
import json
import numpy as np
import wq_common as B

HERE = os.path.dirname(os.path.abspath(__file__))
d = np.load(os.path.join(HERE, "baseline_cache", "baseline.npz"), allow_pickle=True)

RMSE = d["RMSE"]
S = {"old": d["S_old"], "avg": d["S_avg"], "new": d["S_new"]}
C_all = d["C_all"]                      # (8192, 49, 92)
all_nodes = list(d["all_nodes"])
TRUE = {"old": B.KW_OLD_TRUE, "avg": B.KW_AVG_TRUE, "new": B.KW_NEW_TRUE}
PRIOR = B.PRIOR
C_MIN = 0.2

N_RESID = len(B.MONITOR_NODES) * d["obs_glue"].shape[0]     # 6 x 49 = 294
sd_rmse = B.SIGMA_OBS / np.sqrt(2 * N_RESID)
print(f"N residuals = {N_RESID}, sd(RMSE) at truth ~ sigma/sqrt(2N) = {sd_rmse:.4f} mg/L")
print(f"min RMSE observed = {RMSE.min():.4f} mg/L\n")


def prior_sd(g):
    a, b = PRIOR[g]
    return (b - a) / np.sqrt(12)


def weighted_quantile(vals, w, q):
    o = np.argsort(vals)
    v = vals[o]
    cw = np.cumsum(w[o])
    cw /= cw[-1]
    return v[np.searchsorted(cw, min(q, 1.0))]


def band_width(bidx, wb, node):
    k = all_nodes.index(node)
    widths = []
    for t in range(C_all.shape[1]):
        vals = C_all[bidx, t, k]
        lo = weighted_quantile(vals, wb, 0.05)
        hi = weighted_quantile(vals, wb, 0.95)
        widths.append(hi - lo)
    return float(np.mean(widths))


THRESHOLDS = [0.107, 0.11, 0.12]
rows = []
for thr in THRESHOLDS:
    beh = RMSE < thr
    L = np.exp(-0.5 * (RMSE / B.SIGMA_OBS) ** 2)
    w = L * beh
    w = w / w.sum()

    def wmean(x):
        return float(np.sum(w * x))

    def wstd(x):
        m = wmean(x)
        return float(np.sqrt(np.sum(w * (x - m) ** 2)))

    stats = {}
    for g in ["old", "avg", "new"]:
        stats[g] = {
            "mean": wmean(S[g]),
            "sd": wstd(S[g]),
            "sd_retained": wstd(S[g]) / prior_sd(g),
        }

    bidx = np.where(beh)[0]
    wb = w[bidx] / w[bidx].sum()
    band15 = band_width(bidx, wb, "15")
    band107 = band_width(bidx, wb, "107")

    P_node = B.pbar_field(C_all, C_MIN, w)
    order = np.argsort(P_node)[::-1]
    top = [(all_nodes[i], round(float(P_node[i]), 3)) for i in order[:6]]
    n_nonzero = int((P_node > 0).sum())

    rows.append({
        "threshold": thr,
        "count": int(beh.sum()),
        "retention": round(float(beh.mean()), 4),
        "sd_above_floor": round((thr - B.SIGMA_OBS) / sd_rmse, 2),
        "old": stats["old"], "avg": stats["avg"], "new": stats["new"],
        "band_node15": round(band15, 4),
        "band_node107": round(band107, 4),
        "risk_nonzero_nodes": n_nonzero,
        "top_risk": top,
    })

print(f"{'thr':>6} {'count':>6} {'retain':>7} {'SDabove':>8} "
      f"{'old_ret':>8} {'avg_ret':>8} {'new_ret':>8} {'band15':>7} {'band107':>8} {'riskNodes':>9}")
for r in rows:
    print(f"{r['threshold']:>6} {r['count']:>6} {r['retention']:>7.3f} {r['sd_above_floor']:>8} "
          f"{r['old']['sd_retained']*100:>7.1f}% {r['avg']['sd_retained']*100:>7.1f}% "
          f"{r['new']['sd_retained']*100:>7.1f}% {r['band_node15']:>7.3f} {r['band_node107']:>8.3f} "
          f"{r['risk_nonzero_nodes']:>9}")

print("\nweighted means (old/avg/new) per threshold:")
for r in rows:
    print(f"  thr={r['threshold']}: "
          f"old {r['old']['mean']:+.3f}+/-{r['old']['sd']:.3f} | "
          f"avg {r['avg']['mean']:+.3f}+/-{r['avg']['sd']:.3f} | "
          f"new {r['new']['mean']:+.3f}+/-{r['new']['sd']:.3f}")

print("\ntop-6 risk nodes per threshold:")
for r in rows:
    print(f"  thr={r['threshold']}: {r['top_risk']}")

with open(os.path.join(HERE, "baseline_cache", "step3_threshold.json"), "w") as f:
    json.dump({"weighting": "informal_glue",
               "weighting_role": "comparator only — a threshold sweep is a property of behavioural "
                                 "weighting; the primary formal likelihood has no threshold, so no "
                                 "row here may be quoted as a headline result",
               "primary_weighting_elsewhere": B.PRIMARY_WEIGHTING,
               "thresholds": THRESHOLDS, "rows": rows}, f, indent=2)
print("\nsaved step3_threshold.json")
