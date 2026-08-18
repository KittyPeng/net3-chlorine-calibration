"""Step 16: is the metric-specific shortlist result an artefact of one assessment window?

Section 3.1.2 records that the water-age and cumulative-deficit criteria are never satisfied inside
the 168 h horizon, so the risk field is not cyclostationary at 120 h even though the concentration
field is. Section 3.6.2 then compares perturbed shortlists against reference shortlists, and every
one of those comparisons is computed on the single 120-168 h window. A window shared by both arms
cancels the horizon effect common to them, but it does not by itself establish that the RELATIVE
reordering survives a different window, and the duration shortlist is decided by near-ties.

This step tests that directly. Calibration is left exactly where the paper puts it, on the
120-168 h monitor series, because the concentration criteria that justify that choice do pass at
120 h; only the window over which D and A are integrated moves, to the adjacent 96-144 h cycle
pair. Everything else follows Steps 8b and 8c: 30 noise realisations, the risk field taken as the
median across them, top-6 Jaccard and Spearman against the reference field on the SAME window.
If the verdicts agree across the two windows, the Section 3.6.2 conclusion is a property of the
perturbation rather than of the window.

The forward model has to be re-run: baseline.npz and step8b_preds_kb*.npy store only the
post-warm-up slice, so no cached array reaches back to 96 h. About thirteen minutes.
"""
import os
import json
import time
import numpy as np
from scipy.stats import spearmanr
import wq_common as B

HERE = os.path.dirname(os.path.abspath(__file__))
CACHEDIR = os.path.join(HERE, "baseline_cache")
C_MIN = 0.2
TOP_K = 6
N_NOISE = 30
SEEDS = list(range(42, 42 + N_NOISE))
WINDOWS = {"120-168": (120, 169), "96-144": (96, 145)}     # paper window, adjacent cycle pair
PAPER_WIN = "120-168"
KB_ARMS = [-0.4, -0.6]
BIAS_ARMS = [("15", 0.10), ("231", -0.10)]                 # a shortlist mover, and the headline arm
SCHEME = B.PRIMARY_WEIGHTING

cache = np.load(os.path.join(CACHEDIR, "baseline.npz"), allow_pickle=True)
ALL_NODES = list(cache["all_nodes"])
mon_pos = list(cache["mon_pos"])
truth_mon = cache["truth_all"][:, mon_pos]                 # true field (k_b = -0.5) at the monitors
S = [cache["S_old"], cache["S_avg"], cache["S_new"]]


def library(kb):
    """Monitor predictions on the calibration window, and D and A per node on each risk window."""
    Cmon = np.empty((B.N_MC, B.DURATION_H - B.WARMUP_H + 1, len(mon_pos)), np.float32)
    D = {w: np.empty((B.N_MC, len(ALL_NODES)), np.float64) for w in WINDOWS}
    A = {w: np.empty((B.N_MC, len(ALL_NODES)), np.float64) for w in WINDOWS}
    t0 = time.time()
    for s in range(B.N_MC):
        full = B.simulate_chlorine(kb, 0.0, monitor_nodes=ALL_NODES,
                                   pre_run=B.make_kw_hook(S[0][s], S[1][s], S[2][s])).values
        Cmon[s] = full[B.WARMUP_H:][:, mon_pos]
        for w, (i, j) in WINDOWS.items():
            seg = full[i:j]
            D[w][s] = B.below_threshold_duration(seg, C_MIN)
            A[w][s] = np.trapezoid(np.clip(C_MIN - seg, 0.0, None), dx=1.0, axis=0)
        if (s + 1) % 2000 == 0:
            print(f"    {s + 1}/{B.N_MC} ({time.time() - t0:.0f}s)", flush=True)
    print(f"    done in {time.time() - t0:.0f}s")
    return Cmon, D, A


def median_fields(Cmon, D, A, bias_col=None, offset=0.0):
    """Median-over-noise risk fields per window, following the Step 8b/8c protocol."""
    acc = {w: {"P": [], "A": []} for w in WINDOWS}
    for seed in SEEDS:
        rng = np.random.default_rng(seed)
        obs = truth_mon + rng.normal(0, B.SIGMA_OBS, truth_mon.shape)
        if bias_col is not None:
            obs[:, bias_col] += offset
        obs = np.clip(obs, 0, None)[B.WARMUP_H:]
        w, _ = B.all_weightings(Cmon.astype(np.float64), obs, threshold=B.RMSE_THR,
                                schemes=[SCHEME])[SCHEME]
        for win in WINDOWS:
            acc[win]["P"].append(np.tensordot(w, D[win], axes=(0, 0)) / (WINDOWS[win][1]
                                                                        - WINDOWS[win][0] - 1))
            acc[win]["A"].append(np.tensordot(w, A[win], axes=(0, 0)))
    return {win: {k: np.median(np.vstack(v), axis=0) for k, v in acc[win].items()} for win in WINDOWS}


def compare(pert, ref):
    top_p = [ALL_NODES[i] for i in np.argsort(pert)[::-1][:TOP_K]]
    top_r = [ALL_NODES[i] for i in np.argsort(ref)[::-1][:TOP_K]]
    return {"spearman": float(spearmanr(pert, ref).statistic),
            f"top{TOP_K}_jaccard": len(set(top_p) & set(top_r)) / len(set(top_p) | set(top_r)),
            f"top{TOP_K}": top_p}


print("=== Step 16: assessment-window sensitivity of the Section 3.6.2 shortlist comparisons ===")
print(f"  calibration fixed on {B.WARMUP_H}-{B.DURATION_H} h; risk windows "
      f"{', '.join(f'{a}-{b - 1} h' for a, b in WINDOWS.values())}\n")

print("  reference library, k_b = -0.5")
Cb, Db, Ab = library(B.KB_FIXED)
ref = median_fields(Cb, Db, Ab)

rows = []
for kb in KB_ARMS:
    print(f"  library, k_b = {kb}")
    Ck, Dk, Ak = library(kb)
    fld = median_fields(Ck, Dk, Ak)
    for win in WINDOWS:
        rows.append({"arm": f"k_b={kb}", "window": win,
                     "duration": compare(fld[win]["P"], ref[win]["P"]),
                     "deficit": compare(fld[win]["A"], ref[win]["A"])})

for node, off in BIAS_ARMS:            # only the observations change, so the library is reused
    print(f"  sensor bias, node {node} {off:+.2f} mg/L")
    fld = median_fields(Cb, Db, Ab, bias_col=B.MONITOR_NODES.index(node), offset=off)
    for win in WINDOWS:
        rows.append({"arm": f"bias {node}{off:+.2f}", "window": win,
                     "duration": compare(fld[win]["P"], ref[win]["P"]),
                     "deficit": compare(fld[win]["A"], ref[win]["A"])})

hdr = f"\n  {'arm':>14} {'window':>9} | {'rho_D':>7} {'J6_D':>5} | {'rho_A':>7} {'J6_A':>5} | deficit steadier"
print(hdr)
for r in rows:
    steadier = r["deficit"][f"top{TOP_K}_jaccard"] >= r["duration"][f"top{TOP_K}_jaccard"]
    print(f"  {r['arm']:>14} {r['window']:>9} | {r['duration']['spearman']:>7.3f} "
          f"{r['duration'][f'top{TOP_K}_jaccard']:>5.2f} | {r['deficit']['spearman']:>7.3f} "
          f"{r['deficit'][f'top{TOP_K}_jaccard']:>5.2f} | {'yes' if steadier else 'NO'}")

# The claim Section 3.6.2 rests on is comparative and per arm: the deficit shortlist survives a
# perturbation that the duration shortlist does not. The window matters only if that verdict flips.
verdicts, flips = {}, []
for arm in dict.fromkeys(r["arm"] for r in rows):
    v = {r["window"]: r["deficit"][f"top{TOP_K}_jaccard"] >= r["duration"][f"top{TOP_K}_jaccard"]
         for r in rows if r["arm"] == arm}
    verdicts[arm] = v
    if len(set(v.values())) > 1:
        flips.append(arm)
deficit_intact = {arm: {w: next(r["deficit"][f"top{TOP_K}_jaccard"] for r in rows
                                if r["arm"] == arm and r["window"] == w) == 1.0 for w in WINDOWS}
                  for arm in verdicts}

report = {**B.weighting_provenance(comparators=[]),
          "purpose": "does the metric-specific shortlist result of Section 3.6.2 depend on the "
                     "assessment window, given that the risk field is not warm-up converged",
          "calibration_window_h": [B.WARMUP_H, B.DURATION_H],
          "risk_windows_h": {k: [v[0], v[1] - 1] for k, v in WINDOWS.items()},
          "paper_window": PAPER_WIN, "c_min": C_MIN, "top_k": TOP_K, "n_noise": N_NOISE,
          "rows": rows, "deficit_steadier_verdict": verdicts,
          "deficit_top6_fully_preserved": deficit_intact,
          "arms_whose_verdict_flips": flips,
          "verdict_stable_across_windows": len(flips) == 0}
with open(os.path.join(CACHEDIR, "step16_window_sensitivity.json"), "w") as f:
    json.dump(report, f, indent=2)
print(f"\n  arms whose metric-specific verdict flips between windows: {flips or 'none'}")
print("  saved step16_window_sensitivity.json")
