"""Step 4b: single-parameter identifiability curves (why old is one-sided, avg/new flat).

For each coefficient, sweep it across (and a little beyond) its prior range while holding
the other two at their TRUE values, and compute the monitored RMSE against the baseline
noisy observations. This visualises the shape of the GLUE objective:
  - old  : an asymmetric valley (steep on the weak side, flat/saturated on the strong side)
  - avg  : a shallow valley, small on the mg/L scale of the threshold
  - new  : a shallow valley, small on the mg/L scale of the threshold

"Shallow" is reported against BOTH reference scales, because they disagree and the choice drives
the identifiability claim:
  - the behavioural threshold band (thr - noise floor ~ 0.009 mg/L), on which avg/new look flat;
  - the sampling SD of the objective itself, sigma/sqrt(2N) = 0.0041 mg/L, on which the avg/new
    swings are still ~1.7 and ~2.2 SD and therefore NOT beyond detection.
The second scale is the one a formal likelihood uses, which is why the formal profile of Step 7b
constrains avg/new while the informal GLUE score of Step 1 does not.

Cheap (~240 monitor-only EPANET runs). Outputs a figure + a JSON summary.
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
os.makedirs(FIGDIR, exist_ok=True)

cache = np.load(os.path.join(HERE, "baseline_cache", "baseline.npz"), allow_pickle=True)
obs_glue = cache["obs_glue"]                       # (49, 6) baseline noisy observations
TRUE = {"old": B.KW_OLD_TRUE, "avg": B.KW_AVG_TRUE, "new": B.KW_NEW_TRUE}
PRIOR = {"old": B.PRIOR["old"], "avg": B.PRIOR["avg"], "new": B.PRIOR["new"]}

# RMSE at the truth (noise floor reference)
truth_sim = B.simulate_chlorine(
    B.KB_FIXED, 0.0,
    pre_run=B.make_kw_hook(TRUE["old"], TRUE["avg"], TRUE["new"]),
).values[B.WARMUP_H:]
NOISE_FLOOR = float(np.sqrt(((truth_sim - obs_glue) ** 2).mean()))

# sampling SD of the RMSE objective at the truth: SSE/sigma^2 ~ chi2(N) => sd(RMSE) ~ sigma/sqrt(2N)
OBJ_SAMPLING_SD = B.SIGMA_OBS / np.sqrt(2.0 * B.N_RESID)


def rmse_for(kw_old, kw_avg, kw_new):
    sim = B.simulate_chlorine(
        B.KB_FIXED, 0.0,
        pre_run=B.make_kw_hook(kw_old, kw_avg, kw_new),
    ).values[B.WARMUP_H:]
    return float(np.sqrt(((sim - obs_glue) ** 2).mean()))


# display sweep ranges (wider than the prior so the valley shape is visible)
SWEEP = {
    "old": np.linspace(-2.2, -0.05, 80),
    "avg": np.linspace(-0.30, -0.005, 80),
    "new": np.linspace(-0.15, -0.002, 80),
}

results = {}
for g in ["old", "avg", "new"]:
    xs = SWEEP[g]
    rmse = np.empty(len(xs))
    for i, v in enumerate(xs):
        kw = dict(TRUE)
        kw[g] = float(v)
        rmse[i] = rmse_for(kw["old"], kw["avg"], kw["new"])
    a, b = PRIOR[g]
    in_prior = (xs >= min(a, b)) & (xs <= max(a, b))
    swing = float(rmse[in_prior].max() - rmse[in_prior].min())
    results[g] = {
        "x": xs.tolist(),
        "rmse": rmse.tolist(),
        "prior_range": [a, b],
        "truth": TRUE[g],
        "prior_rmse_min": float(rmse[in_prior].min()),
        "prior_rmse_max": float(rmse[in_prior].max()),
        "prior_rmse_swing": swing,
        "prior_rmse_swing_in_sampling_sd": swing / OBJ_SAMPLING_SD,
    }

# ---- figure ----
fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
labels = {"old": "k_w,old", "avg": "k_w,avg", "new": "k_w,new"}
for ax, g in zip(axes, ["old", "avg", "new"]):
    xs = np.array(results[g]["x"])
    rmse = np.array(results[g]["rmse"])
    a, b = results[g]["prior_range"]
    ax.plot(xs, rmse, color="steelblue", lw=2)
    ax.axvspan(min(a, b), max(a, b), color="0.85", label="prior range")
    ax.axvline(results[g]["truth"], color="red", lw=2, label="true value")
    ax.axhline(B.RMSE_THR_DRAFT, color="crimson", ls="--", lw=1,
               label=f"draft threshold {B.RMSE_THR_DRAFT}")
    ax.axhline(B.RMSE_THR, color="darkorange", ls=":", lw=1,
               label=f"primary threshold {B.RMSE_THR}")
    ax.axhline(NOISE_FLOOR, color="green", ls="-.", lw=1, label="noise floor")
    ax.set_xlabel(f"{labels[g]} (m/day)")
    ax.set_ylabel("monitored RMSE (mg/L)")
    ax.set_title(f"{labels[g]}: RMSE swing over prior\n"
                 f"{results[g]['prior_rmse_swing']:.3f} mg/L = "
                 f"{results[g]['prior_rmse_swing_in_sampling_sd']:.1f} x objective sampling SD",
                 fontsize=10)
    ax.grid(alpha=0.3)
axes[0].legend(fontsize=7, loc="upper center")
fig.suptitle("Single-parameter objective (other two held at truth): old is a deep asymmetric "
             "valley; avg/new are shallow\non the threshold scale but still a few objective "
             "sampling SD deep", y=1.04)
plt.tight_layout()
figpath = os.path.join(FIGDIR, "step4b_sensitivity_curves.png")
plt.savefig(figpath, dpi=130, bbox_inches="tight")

summary = {
    "weighting": "none",
    "weighting_role": "this step plots the RMSE objective itself, so no weighting is applied; the "
                      "two thresholds are drawn only to show where the informal comparator would "
                      "cut the curve",
    "primary_weighting_elsewhere": B.PRIMARY_WEIGHTING,
    "hold_others_at": "true values",
    "observations": "baseline noisy (seed 42)",
    "noise_floor_rmse": NOISE_FLOOR,
    "threshold_primary": B.RMSE_THR,
    "threshold_draft": B.RMSE_THR_DRAFT,
    "objective_sampling_sd": float(OBJ_SAMPLING_SD),
    "prior_rmse_swing": {g: results[g]["prior_rmse_swing"] for g in ["old", "avg", "new"]},
    "prior_rmse_swing_in_sampling_sd": {g: results[g]["prior_rmse_swing_in_sampling_sd"]
                                        for g in ["old", "avg", "new"]},
    "figure": os.path.relpath(figpath, HERE),
}
with open(os.path.join(HERE, "baseline_cache", "step4b_sensitivity.json"), "w") as f:
    json.dump({"summary": summary, "sweeps": results}, f, indent=2)

print("noise floor RMSE =", round(NOISE_FLOOR, 4))
print(f"objective sampling SD = sigma/sqrt(2N) = {OBJ_SAMPLING_SD:.4f} mg/L (N = {B.N_RESID})")
print("RMSE swing across each prior (others at truth):")
for g in ["old", "avg", "new"]:
    print(f"  {labels[g]:>9}: min {results[g]['prior_rmse_min']:.4f} -> "
          f"max {results[g]['prior_rmse_max']:.4f}  (swing {results[g]['prior_rmse_swing']:.4f} mg/L "
          f"= {results[g]['prior_rmse_swing_in_sampling_sd']:.1f} x sampling SD)")
print("figure saved to", figpath)
