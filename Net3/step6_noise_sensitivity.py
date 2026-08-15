"""Step 6: noise / sensor-accuracy sensitivity in the three-zone setup.

Re-runs the calibration against the SAME three-zone truth for a range of observation noise levels
σ = {0.02, 0.05, 0.10, 0.15} mg/L (the foundation's 0.02/0.05/0.10 plus the supervisor's requested
±0.05 / ±0.10 / ±0.15). σ is interpreted as ONE standard deviation.

PRIMARY: the formal censored Gaussian likelihood at the same σ that generated the observations.
COMPARATOR: the informal GLUE score with its behavioural threshold scaling as the ~95% acceptance
band of the objective, threshold(σ) = σ(1 + 1.645/√(2N)) (0.107 at σ = 0.1).

The distinction decides the answer. "How accurate must a sensor be?" is a question about what the
data can support, so it has to be asked of an efficient likelihood; asked of the informal score it
returns the score's own inefficiency (a 2.8-3.1x wider posterior at σ = 0.1, Step 7) as if it were a
property of the sensor, and therefore demands a better instrument than the measurement actually
needs. Both are reported so the size of that distortion is visible at every noise level.

Everything reuses the frozen baseline cache (candidate predictions are noise-independent), so no
EPANET is re-run. Each σ is repeated over 30 independent noise realisations; results are reported as
median [IQR]. This answers §3.5 (how the posterior width of each grouped coefficient scales with σ)
and the email question (required sensor accuracy for useful chlorine predictions).
"""
import os
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import wq_common as B

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures")

cache = np.load(os.path.join(HERE, "baseline_cache", "baseline.npz"), allow_pickle=True)
C_all = cache["C_all"]                              # (N_MC, 49, 92), noise-independent candidates
ALL_NODES = list(cache["all_nodes"])
mon_pos = list(cache["mon_pos"])
truth_mon = cache["truth_all"][:, mon_pos]          # three-zone truth at the monitors
S = {"old": cache["S_old"], "average": cache["S_avg"], "new": cache["S_new"]}
C_all_mon = C_all[:, :, mon_pos].astype(np.float64)

ZKEYS = ["old", "average", "new"]
PRIOR_OF = {"old": B.PRIOR["old"], "average": B.PRIOR["avg"], "new": B.PRIOR["new"]}
PRIOR_SD = {z: (PRIOR_OF[z][1] - PRIOR_OF[z][0]) / np.sqrt(12) for z in ZKEYS}
SIGMAS = [0.02, 0.05, 0.10, 0.15]
N_NOISE = 30
NOISE_SEEDS = list(range(42, 42 + N_NOISE))
NODE_BAND = "15"                                    # old-zone downstream node for the predictive band
kk = ALL_NODES.index(NODE_BAND)
SCHEMES = [B.PRIMARY_WEIGHTING, "informal_glue"]
# One effective-sample-size rule, used by the JSON, the printout and the figure alike. A 5-95%
# interval needs roughly 100 effective members before each tail rests on more than a handful.
ESS_MIN = 100


def med_iqr(a):
    a = np.asarray(a, float)
    if a.size == 0:
        return (float("nan"),) * 3
    return float(np.median(a)), float(np.percentile(a, 25)), float(np.percentile(a, 75))


def wquantile(vals, w, q):
    o = np.argsort(vals)
    cw = np.cumsum(w[o])
    cw /= cw[-1]
    return vals[o][np.searchsorted(cw, min(q, 1.0))]


rows = []
for sigma in SIGMAS:
    thr = B.threshold_for_sigma(sigma)
    acc = {s: {"ess": [], "band": [], "retention": [],
               **{z: {"sd": [], "sd_ret": [], "mean": []} for z in ZKEYS}} for s in SCHEMES}
    n_empty = {s: 0 for s in SCHEMES}
    for seed in NOISE_SEEDS:
        rng = np.random.default_rng(seed)
        obs = np.clip(truth_mon + rng.normal(0, sigma, truth_mon.shape), 0, None)[B.WARMUP_H:]
        wts = B.all_weightings(C_all_mon, obs, sigma=sigma, threshold=thr, schemes=SCHEMES)
        rmse = B.rmse_of(C_all_mon, obs)
        for s in SCHEMES:
            w, diag = wts[s]
            if w is None:
                # An empty behavioural set used to be skipped silently, which turns the reported
                # median into a median over the realisations that happened to sample a good
                # parameter set — a selection bias that grows as sigma shrinks. Now it is counted.
                n_empty[s] += 1
                continue
            acc[s]["ess"].append(diag["ess"])
            acc[s]["retention"].append(float((rmse < thr).mean()) if s == "informal_glue"
                                       else float(diag["ess_frac"]))
            for z in ZKEYS:
                m, sd = B.weighted_mean_sd(w, S[z])
                acc[s][z]["mean"].append(m)
                acc[s][z]["sd"].append(sd)
                acc[s][z]["sd_ret"].append(sd / PRIOR_SD[z])
            lo = np.array([wquantile(C_all[:, t, kk], w, 0.05) for t in range(C_all.shape[1])])
            hi = np.array([wquantile(C_all[:, t, kk], w, 0.95) for t in range(C_all.shape[1])])
            acc[s]["band"].append(float(np.mean(hi - lo)))

    row = {"sigma": sigma, "threshold_informal": thr, "n_realisations": len(NOISE_SEEDS),
           "primary_weighting": B.PRIMARY_WEIGHTING, "comparators": ["informal_glue"],
           "ess_criterion": f"median ESS < {ESS_MIN} or any realisation with empty weights",
           "by_scheme": {}}
    for s in SCHEMES:
        ess = np.asarray(acc[s]["ess"], float)
        d = {"n_valid": int(ess.size), "n_empty": n_empty[s],
             "retention_med": med_iqr(acc[s]["retention"])[0],
             "band_node15_med": med_iqr(acc[s]["band"])[0],
             # the full ESS distribution, not just its median: at sigma = 0.05 the median is in the
             # hundreds while the worst realisation is in the tens, and it is the worst realisation
             # that decides whether a 5-95% interval from this row can be quoted
             "ess_med": float(np.median(ess)) if ess.size else float("nan"),
             "ess_p5": float(np.percentile(ess, 5)) if ess.size else float("nan"),
             "ess_min": float(ess.min()) if ess.size else float("nan"),
             "ess_frac_below_min": float((ess < ESS_MIN).mean()) if ess.size else 1.0,
             "sampling_limited": bool(n_empty[s] > 0 or (ess.size and np.median(ess) < ESS_MIN))}
        for z in ZKEYS:
            m, lo, hi = med_iqr(acc[s][z]["sd_ret"])
            d[z] = {"sd_ret_med": m, "sd_ret_iqr": [lo, hi],
                    "sd_med": med_iqr(acc[s][z]["sd"])[0],
                    "mean_med": med_iqr(acc[s][z]["mean"])[0]}
        row["by_scheme"][s] = d
    prim, comp = row["by_scheme"][B.PRIMARY_WEIGHTING], row["by_scheme"]["informal_glue"]
    row["informal_over_formal_sd_ratio"] = {
        z: (comp[z]["sd_med"] / prim[z]["sd_med"] if prim[z]["sd_med"] else None) for z in ZKEYS}
    rows.append(row)

# ---- report ----
print("=== Step 6: sensor-accuracy sensitivity (three-zone) ===")
print(f"primary {B.PRIMARY_WEIGHTING}; comparator informal GLUE (threshold scales with σ)")
print(f"sampling-limited rule: median ESS < {ESS_MIN} or any empty-weight realisation\n")
for s in SCHEMES:
    print(f"--- {s} ---")
    print(f"{'σ':>5} {'ESSmed':>7} {'ESSp5':>7} {'ESSmin':>7} {'<100':>5} {'band15':>7} | "
          f"{'old SDret':>18} {'avg SDret':>18} {'new SDret':>18}")
    for r in rows:
        d = r["by_scheme"][s]

        def f(z, _d=d):
            q = _d[z]
            return (f"{q['sd_ret_med'] * 100:4.0f}% "
                    f"[{q['sd_ret_iqr'][0] * 100:.0f}-{q['sd_ret_iqr'][1] * 100:.0f}]")

        flag = " <- sampling-limited" if d["sampling_limited"] else ""
        print(f"{r['sigma']:>5} {d['ess_med']:>7.1f} {d['ess_p5']:>7.1f} {d['ess_min']:>7.1f} "
              f"{d['ess_frac_below_min'] * 100:>4.0f}% {d['band_node15_med']:>7.3f} | "
              f"{f('old'):>18} {f('average'):>18} {f('new'):>18}{flag}")
    print()
print("SD retained = posterior SD as a fraction of the prior SD; band15 = mean width of the 5-95%")
print("predictive band at node 15. ESS columns are the median, 5th percentile and minimum over the")
print(f"{N_NOISE} realisations, plus the share of realisations below ESS = {ESS_MIN}.\n")
print("formal / informal posterior SD ratio (how much of the required accuracy is the score's):")
for r in rows:
    print(f"  σ = {r['sigma']:<5} " + ", ".join(
        f"{z} x{r['informal_over_formal_sd_ratio'][z]:.2f}" for z in ZKEYS))
prim_50 = next(r for r in rows if r["sigma"] == 0.05)["by_scheme"][B.PRIMARY_WEIGHTING]
prim_100 = next(r for r in rows if r["sigma"] == 0.10)["by_scheme"][B.PRIMARY_WEIGHTING]
print(f"\nUnder the primary likelihood a σ = 0.10 sensor already retains only "
      f"{prim_100['old']['sd_ret_med'] * 100:.0f}/{prim_100['average']['sd_ret_med'] * 100:.0f}/"
      f"{prim_100['new']['sd_ret_med'] * 100:.0f}% of the prior width (old/avg/new) against "
      f"{prim_50['old']['sd_ret_med'] * 100:.0f}/{prim_50['average']['sd_ret_med'] * 100:.0f}/"
      f"{prim_50['new']['sd_ret_med'] * 100:.0f}% at σ = 0.05, so the accuracy the measurement "
      f"needs is milder than the informal comparator implies.")


def _jsafe(o):
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.integer):
        return int(o)
    return str(o)


with open(os.path.join(HERE, "baseline_cache", "step6_noise_sensitivity.json"), "w") as f:
    json.dump({**B.weighting_provenance(comparators=["informal_glue"]),
               "ess_min_criterion": ESS_MIN, "n_noise": N_NOISE, "rows": rows},
              f, indent=2, default=_jsafe)

# ---- figure: SD retained (%) vs σ, primary vs comparator, and the predictive band ----
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.6))
sig = [r["sigma"] for r in rows]
colors = {"old": "tab:red", "average": "tab:orange", "new": "tab:green"}
for z in ZKEYS:
    prim = [r["by_scheme"][B.PRIMARY_WEIGHTING][z] for r in rows]
    ax1.plot(sig, [p["sd_ret_med"] * 100 for p in prim], marker="o", color=colors[z],
             label=f"{z} (formal censored)")
    ax1.fill_between(sig, [p["sd_ret_iqr"][0] * 100 for p in prim],
                     [p["sd_ret_iqr"][1] * 100 for p in prim], color=colors[z], alpha=0.15)
    ax1.plot(sig, [r["by_scheme"]["informal_glue"][z]["sd_ret_med"] * 100 for r in rows],
             marker="s", ms=4, ls=":", lw=1.2, color=colors[z], alpha=0.8,
             label=f"{z} (informal GLUE)")
ax1.axhline(70, color="gray", ls="--", lw=1, label="70% of prior retained")
for r in rows:                        # flag on the SAME ESS rule the JSON uses, not on retention
    if r["by_scheme"][B.PRIMARY_WEIGHTING]["sampling_limited"]:
        ax1.axvspan(r["sigma"] - 0.006, r["sigma"] + 0.006, color="red", alpha=0.08)
        ax1.annotate(f"ESS < {ESS_MIN}", (r["sigma"], 8), fontsize=7, ha="center", color="crimson")
ax1.set_xlabel("observation σ (mg/L)")
ax1.set_ylabel("posterior SD retained (% of prior)")
ax1.set_title("Coefficient identifiability vs sensor noise\n"
              "solid = formal censored (primary), dotted = informal GLUE (comparator)", fontsize=10)
ax1.legend(fontsize=6.5, ncol=2)
ax1.grid(alpha=0.3)

for s, mk, ls in ((B.PRIMARY_WEIGHTING, "s", "-"), ("informal_glue", "^", ":")):
    ax2.plot(sig, [r["by_scheme"][s]["band_node15_med"] for r in rows], marker=mk, ls=ls,
             color="steelblue" if s == B.PRIMARY_WEIGHTING else "0.5",
             label="formal censored" if s == B.PRIMARY_WEIGHTING else "informal GLUE")
ax2.set_xlabel("observation σ (mg/L)")
ax2.set_ylabel("5–95% predictive band at node 15 (mg/L)")
ax2.set_title("Prediction uncertainty vs sensor noise")
ax2.legend(fontsize=8)
ax2.grid(alpha=0.3)
fig.suptitle("Step 6 — sensor-accuracy sensitivity: the required accuracy depends on the inference "
             "rule, not only on the sensor", y=1.02)
plt.tight_layout()
figpath = os.path.join(FIGDIR, "step6_noise_sensitivity.png")
plt.savefig(figpath, dpi=130, bbox_inches="tight")
print("figure saved to", figpath)
