"""Step 8d: sensor DRIFT — a time-varying offset — Priority-2 #1, second half.

Step 8 and Step 8c model a *constant* sensor bias. The review's comment is headed "systematic sensor
error (bias, drift)", and a drifting sensor is the more common failure in practice: a probe fouls or
its calibration ages, so the offset grows over the record instead of sitting still. This step closes
that half.

DESIGN. Drift is a linear ramp across the 48 h assessment window:

    b(t) = D * (t - t0) / 48,    t in [120, 168] h

so the sensor reads correctly at the start of the window and is off by D at the end; its MEAN offset
over the window is D/2. Nothing is added before the window, which is irrelevant anyway because only
the post-warm-up slice enters the likelihood.

THE QUESTION THIS ANSWERS. A drift has a mean and a shape. If a drift did the same damage as a
constant bias of the same mean, then modelling constant bias alone would be sufficient and Step 8
would already cover this comment. So every drift arm is run against TWO constant-bias controls on
identical noise:

    const(D/2)   the MEAN-equivalent constant bias  — does the shape add anything?
    const(D)     the END-equivalent constant bias   — does the drift reach the same damage?

If drift(0->D) lands on const(D/2) the answer is "a drift is its own average" and constant-bias
results transfer; if it lands nearer const(D) the time structure matters on its own.

Two-sided, because observations are censored at the sensor floor and a negative offset is therefore
not the mirror image of a positive one (same reason as Step 8).

Two locations, both taken from what the earlier steps established: node 15 (old zone) is the Step 8
headline and the physically dominant coefficient; node 231 (average zone) is where Step 8c found the
largest NORMALISED corruption (-5.94 posterior SD at -0.10).

Weighting is the primary formal censored likelihood: Steps 3, 5c, 8 and 8b all show the informal
score is too flat to register a systematic error, so measuring drift with it would understate it.

Reuses the baseline cache — the candidate predictions do not depend on the observations, so no
EPANET is re-run. 30 noise realisations, medians reported.
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
ZONE_OF_NODE = {"107": "new", "113": "new", "15": "old", "145": "old",
                "209": "average", "231": "average"}

cache = np.load(os.path.join(HERE, "baseline_cache", "baseline.npz"), allow_pickle=True)
C_all = cache["C_all"]
ALL_NODES = list(cache["all_nodes"])
mon_pos = list(cache["mon_pos"])
truth_mon = cache["truth_all"][:, mon_pos]              # (169, 6) full record
S = {"old": cache["S_old"], "average": cache["S_avg"], "new": cache["S_new"]}
C_all_mon = C_all[:, :, mon_pos]
C_MIN = 0.2
# Two risk metrics: P_bar (how long below the threshold) and E[A] (how long AND how far below).
# Step 10's headline top-10 is ranked by E[A], and Step 8b showed the two can disagree about which
# nodes lead, so a rank claim made only on P_bar is not a claim about the published list.
DUR = B.below_threshold_duration(C_all, C_MIN)
TW = B.window_hours(C_all)
DEFC = np.trapezoid(np.clip(C_MIN - C_all, 0.0, None), dx=1.0, axis=1).astype(np.float64)

DRIFT_NODES = ["15", "231"]
MAGNITUDES = [-0.10, -0.05, 0.05, 0.10]
N_NOISE = 30
SEEDS = list(range(42, 42 + N_NOISE))
TN = C_all.shape[1]                                    # 49 reporting points in the window
RAMP = np.arange(TN, dtype=float) / (TN - 1)           # 0 -> 1 across the window; mean 0.5


def offset_profile(kind, D):
    """The additive error at each of the TN window points, in mg/L."""
    if kind == "drift":
        return D * RAMP                                # 0 -> D, mean D/2
    if kind == "const_mean":
        return np.full(TN, D / 2.0)                    # mean-equivalent control
    if kind == "const_end":
        return np.full(TN, D)                          # end-equivalent control
    raise ValueError(kind)


def evaluate(node, kind, D):
    """Median posterior summary and risk field over the 30 noise realisations."""
    col = B.MONITOR_NODES.index(node) if node is not None else None
    prof = None if node is None else offset_profile(kind, D)
    means_by_z = {z: [] for z in ZKEYS}
    sds_by_z = {z: [] for z in ZKEYS}
    ranks, Ps, As, n_clipped = [], [], [], []
    for seed in SEEDS:
        rng = np.random.default_rng(seed)
        obs = truth_mon + rng.normal(0, B.SIGMA_OBS, truth_mon.shape)
        obs = obs[B.WARMUP_H:].copy()                  # (TN, 6); only the window is used
        if prof is not None:
            obs[:, col] += prof
        n_clipped.append(int((obs < 0).sum()))
        obs = np.clip(obs, 0, None)
        w, _ = B.all_weightings(C_all_mon, obs, schemes=[B.PRIMARY_WEIGHTING])[B.PRIMARY_WEIGHTING]
        for z in ZKEYS:
            m, sd = B.weighted_mean_sd(w, S[z])
            means_by_z[z].append(m)
            sds_by_z[z].append(sd)
        P = np.tensordot(w, DUR, axes=(0, 0)) / TW
        Ps.append(P)
        As.append(np.tensordot(w, DEFC, axes=(0, 0)))
        ranks.append(tuple(ALL_NODES[i] for i in np.argsort(P)[::-1][:6]))
    P_med = np.median(np.vstack(Ps), axis=0)
    A_med = np.median(np.vstack(As), axis=0)
    # The leading set is read off the SAME median field the Spearman uses. The modal set over the
    # per-realisation orderings answers a different question ("which ordering occurs most often")
    # and is kept as its own field rather than mixed in: reporting one of them under a metadata line
    # that declares the other is how a single table ends up with two incompatible bases.
    return {"means": {z: float(np.median(means_by_z[z])) for z in ZKEYS},
            "sds": {z: float(np.median(sds_by_z[z])) for z in ZKEYS},
            "n_censored_med": float(np.median(n_clipped)),
            "top6": [str(ALL_NODES[i]) for i in np.argsort(P_med)[::-1][:6]],
            "modal_top6_across_realisations": [str(x) for x in max(set(ranks), key=ranks.count)],
            "deficit_top6": [str(ALL_NODES[i]) for i in np.argsort(A_med)[::-1][:6]],
            "_P": P_med, "_A": A_med}


print("=== Step 8d: sensor DRIFT (time-varying offset) vs its constant-bias controls ===")
print(f"weighting: {B.PRIMARY_WEIGHTING}; {N_NOISE} noise realisations; medians reported")
print(f"drift = linear ramp 0 -> D across the {TN - 1} h window (mean offset D/2)\n")

base = evaluate(None, None, 0.0)
P_ref, ref_top6 = base["_P"], set(base["top6"])
A_ref, ref_top6_A = base["_A"], set(base["deficit_top6"])
print("unbiased baseline: " + ", ".join(
    f"{z} {base['means'][z]:+.4f} (SD {base['sds'][z]:.4f})" for z in ZKEYS))
print()

rows = []
for node in DRIFT_NODES:
    own = ZONE_OF_NODE[node]
    print(f"--- monitor {node} ({own} zone); own coefficient = k_w,{own} ---")
    print(f"{'D':>7} | {'arm':>11} | {'own mean':>9} | {'shift':>8} | {'shift/SD':>8} | "
          f"{'cens':>5} | {'rho_S':>6} {'tau_K':>6} {'J6':>5}")
    for D in MAGNITUDES:
        for kind in ("drift", "const_mean", "const_end"):
            r = evaluate(node, kind, D)
            shift = r["means"][own] - base["means"][own]
            row = {"node": node, "zone": own, "D": D, "arm": kind,
                   "means": r["means"], "sds": r["sds"],
                   "own_coef": own,
                   "own_shift": shift,
                   "own_shift_over_sd": shift / base["sds"][own],
                   "shift_over_sd_all": {z: (r["means"][z] - base["means"][z]) / base["sds"][z]
                                         for z in ZKEYS},
                   "n_censored_med": r["n_censored_med"],
                   "risk_spearman_vs_unbiased": float(spearmanr(r["_P"], P_ref).statistic),
                   "risk_kendall_vs_unbiased": float(kendalltau(r["_P"], P_ref).statistic),
                   "risk_top6_jaccard_vs_unbiased":
                       len(set(r["top6"]) & ref_top6) / len(set(r["top6"]) | ref_top6),
                   "top6": r["top6"],
                   "modal_top6_across_realisations": r["modal_top6_across_realisations"],
                   "deficit_spearman_vs_unbiased": float(spearmanr(r["_A"], A_ref).statistic),
                   "deficit_top6_jaccard_vs_unbiased":
                       len(set(r["deficit_top6"]) & ref_top6_A)
                       / len(set(r["deficit_top6"]) | ref_top6_A),
                   "deficit_top6": r["deficit_top6"]}
            rows.append(row)
            print(f"{D:>+7.3f} | {kind:>11} | {r['means'][own]:>+9.4f} | {shift:>+8.4f} | "
                  f"{row['own_shift_over_sd']:>+8.2f} | {r['n_censored_med']:>5.0f} | "
                  f"{row['risk_spearman_vs_unbiased']:>6.3f} "
                  f"{row['risk_kendall_vs_unbiased']:>6.3f} "
                  f"{row['risk_top6_jaccard_vs_unbiased']:>5.2f}")
    print()


def get(node, D, arm):
    return next(r for r in rows if r["node"] == node and r["D"] == D and r["arm"] == arm)


# Is a drift equivalent to its mean-equivalent constant bias? The ratio is the whole point of the
# step, so it is computed rather than eyeballed off the table.
equivalence = {}
for node in DRIFT_NODES:
    own_sd = base["sds"][ZONE_OF_NODE[node]]
    per_D = {}
    for D in MAGNITUDES:
        d = get(node, D, "drift")["own_shift"]
        cm = get(node, D, "const_mean")["own_shift"]
        ce = get(node, D, "const_end")["own_shift"]
        per_D[f"{D:+.3f}"] = {
            "drift_shift": d, "const_mean_shift": cm, "const_end_shift": ce,
            "drift_over_const_mean": d / cm if cm else None,
            "drift_over_const_end": d / ce if ce else None,
            # the residual, in the units the rest of the log uses: how far the drift lands from its
            # mean-equivalent control measured against the unbiased posterior spread
            "drift_minus_const_mean_in_sd": (d - cm) / own_sd}
    equivalence[node] = per_D

print("Is a drift just its own mean? (drift shift / mean-equivalent constant shift)")
for node, per_D in equivalence.items():
    vals = [v["drift_over_const_mean"] for v in per_D.values() if v["drift_over_const_mean"]]
    print(f"  monitor {node}: " + ", ".join(
        f"D={k} {v['drift_over_const_mean']:.2f}" for k, v in per_D.items()) +
        f"   (mean {np.mean(vals):.2f})")
print("A ratio of 1.00 means the drift does exactly the damage of a constant bias at its own mean,")
print("i.e. only the mean of a drift matters and Step 8's constant-bias sweep already covers it.")

report = {**B.weighting_provenance(comparators=[]),
          "drift_nodes": DRIFT_NODES,
          "magnitudes": MAGNITUDES,
          "n_noise": N_NOISE,
          "n_window_points": int(TN),
          "drift_model": "linear ramp b(t) = D*(t-t0)/48 across the 48 h assessment window; the "
                         "sensor is correct at the window start and off by D at its end, so the "
                         "mean offset is D/2",
          "arms": {"drift": "linear ramp 0 -> D",
                   "const_mean": "constant bias at D/2, the mean-equivalent control",
                   "const_end": "constant bias at D, the end-equivalent control"},
          "risk_metric": {
              "P_bar": "expected fraction of the 48 h window below 0.2 mg/L (risk_* fields)",
              "E_A": "expected cumulative deficit in mg/L*h (deficit_* fields); the metric Step 10's "
                     "headline top-10 is ranked by"},
          "risk_ranking_basis": "median risk field over the 30 noise realisations; the top-6 "
                                "set and the Spearman are read off the SAME field, and the "
                                "modal per-realisation ordering is reported separately as "
                                "modal_top6_across_realisations",
          "baseline": {"means": base["means"], "sds": base["sds"],
                       "n_censored_med": base["n_censored_med"], "top6": base["top6"]},
          "equivalence": equivalence,
          "rows": [{k: v for k, v in r.items() if not k.startswith("_")} for r in rows]}

with open(os.path.join(HERE, "baseline_cache", "step8d_sensor_drift.json"), "w") as f:
    json.dump(report, f, indent=2)

# ---- figure: own-coefficient shift in posterior SD, drift against both controls ----
fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharex=True)
for ax, node in zip(axes, DRIFT_NODES):
    own = ZONE_OF_NODE[node]
    for kind, style in (("drift", dict(marker="o", color="steelblue", lw=2)),
                        ("const_mean", dict(marker="s", color="darkorange", ls="--")),
                        ("const_end", dict(marker="^", color="grey", ls=":"))):
        ys = [get(node, D, kind)["own_shift_over_sd"] for D in MAGNITUDES]
        ax.plot(MAGNITUDES, ys, label=kind.replace("_", " "), **style)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_title(f"monitor {node} ({own} zone)")
    ax.set_xlabel("drift magnitude D at the end of the window (mg/L)")
    ax.set_ylabel(f"shift of k_w,{own} (posterior SD)")
    ax.grid(alpha=0.3)
    ax.legend()
fig.suptitle("Step 8d — a drifting sensor against constant-bias controls of the same mean and end "
             "value\n(formal censored likelihood, 30 noise realisations)")
plt.tight_layout()
figpath = os.path.join(FIGDIR, "step8d_sensor_drift.png")
plt.savefig(figpath, dpi=130, bbox_inches="tight")
print("\nfigure saved to", figpath)
