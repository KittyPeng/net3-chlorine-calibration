"""Step 4: displaced-prior experiment, first look (informal GLUE, one noise realisation).

Same truth, same monitoring array, same baseline noisy observations. Only the prior ranges change:
each range keeps its width but is shifted so its midpoint sits one prior standard deviation from the
true value (all ranges stay non-positive and still contain the truth). Question: do the observations
pull the estimate back toward the truth?

SUPERSEDED BY STEP 4d for every reported conclusion. This version uses the informal comparator, one
noise realisation and one displacement direction, so it can neither separate a genuine pull-back from
a lucky noise draw nor separate the data's information from the score's flatness. It is kept because
it is the configuration the draft ran, and because the contrast with Step 4d is instructive: 4d
displaces in both directions, repeats over 30 realisations and reports the primary formal rule
alongside the comparator.
"""
import os
import json
import time
import numpy as np
import wq_common as B

HERE = os.path.dirname(os.path.abspath(__file__))
cache = np.load(os.path.join(HERE, "baseline_cache", "baseline.npz"), allow_pickle=True)
obs_glue = cache["obs_glue"]                       # identical baseline observations
ALL_NODES = list(cache["all_nodes"])
mon_pos = list(cache["mon_pos"])

TRUE = {"old": B.KW_OLD_TRUE, "avg": B.KW_AVG_TRUE, "new": B.KW_NEW_TRUE}


def prior_sd(rng_tuple):
    a, b = rng_tuple
    return (b - a) / np.sqrt(12)


def displaced_range(orig, truth):
    """Keep width; put midpoint at truth - 1*prior_SD (downward keeps values non-positive)."""
    a, b = orig
    width = b - a
    sd = width / np.sqrt(12)
    mid = truth - sd
    return (mid - width / 2, mid + width / 2)


ORIG = B.PRIOR                                    # {'old','avg','new'}
DISP = {g: displaced_range(ORIG[g], TRUE[{"old": "old", "avg": "avg", "new": "new"}[g]])
        for g in ["old", "avg", "new"]}
# guard: non-positive and truth inside
for g in ["old", "avg", "new"]:
    a, b = DISP[g]
    assert b <= 0, f"{g} upper bound positive"
    assert a <= TRUE[g] <= b, f"{g} truth not in displaced range"

rng = np.random.default_rng(B.SAMPLE_SEED)
S = {
    "old": rng.uniform(*DISP["old"], B.N_MC),
    "avg": rng.uniform(*DISP["avg"], B.N_MC),
    "new": rng.uniform(*DISP["new"], B.N_MC),
}

Tn = obs_glue.shape[0]
RMSE = np.empty(B.N_MC)
t0 = time.time()
for s in range(B.N_MC):
    cs = B.simulate_chlorine(
        B.KB_FIXED, 0.0,
        pre_run=B.make_kw_hook(S["old"][s], S["avg"][s], S["new"][s]),
        monitor_nodes=ALL_NODES,
    ).values[B.WARMUP_H:]
    d = cs[:, mon_pos] - obs_glue
    RMSE[s] = np.sqrt((d ** 2).mean())
    if (s + 1) % 500 == 0:
        print(f"  {s + 1}/{B.N_MC} ({time.time() - t0:.1f}s)")

L = np.exp(-0.5 * (RMSE / B.SIGMA_OBS) ** 2)
beh = RMSE < B.RMSE_THR
w = L * beh
w = w / w.sum()


def wmean(x):
    return float(np.sum(w * x))


def wstd(x):
    m = wmean(x)
    return float(np.sqrt(np.sum(w * (x - m) ** 2)))


rows = []
for g in ["old", "avg", "new"]:
    sd = prior_sd(DISP[g])                          # same as original width sd
    disp_mid = 0.5 * (DISP[g][0] + DISP[g][1])
    bm = wmean(S[g])
    bsd = wstd(S[g])
    move_toward_truth = (bm - disp_mid) / sd        # truth is +1 SD above disp_mid
    rows.append({
        "group": g,
        "orig_mid": round(0.5 * (ORIG[g][0] + ORIG[g][1]), 4),
        "disp_mid": round(disp_mid, 4),
        "disp_range": [round(DISP[g][0], 4), round(DISP[g][1], 4)],
        "true": TRUE[g],
        "beh_mean": round(bm, 4),
        "beh_sd": round(bsd, 4),
        "sd_retained": round(bsd / sd, 3),
        "move_toward_truth_priorSD": round(move_toward_truth, 3),
        "gap_closed_frac": round(move_toward_truth / 1.0, 3),  # gap = 1 prior SD
    })

summary = {
    "weighting": "informal_glue",
    "weighting_role": "comparator; SUPERSEDED by step4d, which reports the primary formal rule over "
                      "30 noise realisations and both displacement directions",
    "primary_weighting_elsewhere": B.PRIMARY_WEIGHTING,
    "n_noise_realisations": 1,
    "behavioural_count": int(beh.sum()),
    "retention": round(float(beh.mean()), 4),
    "rmse_min": round(float(RMSE.min()), 4),
    "rows": rows,
}
print("\n=== displaced-prior, informal GLUE, single realisation (superseded by Step 4d) ===")
print(f"behavioural {summary['behavioural_count']}/{B.N_MC} "
      f"({summary['retention']*100:.1f}%) | min RMSE {summary['rmse_min']}")
print(f"{'grp':>4} {'origMid':>8} {'dispMid':>8} {'true':>7} {'behMean':>8} "
      f"{'behSD':>7} {'SDret':>6} {'moveToTruth/SD':>15} {'gapClosed':>10}")
for r in rows:
    print(f"{r['group']:>4} {r['orig_mid']:>8} {r['disp_mid']:>8} {r['true']:>7} "
          f"{r['beh_mean']:>8} {r['beh_sd']:>7} {r['sd_retained']*100:>5.1f}% "
          f"{r['move_toward_truth_priorSD']:>15} {r['gap_closed_frac']*100:>9.0f}%")

with open(os.path.join(HERE, "baseline_cache", "step4_displaced_prior.json"), "w") as f:
    json.dump(summary, f, indent=2)
print("\nsaved step4_displaced_prior.json")
