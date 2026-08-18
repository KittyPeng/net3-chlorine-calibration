"""Step 8c: does the sensor-bias LOCATION matter? Inject the same offset at each of the six monitors
in turn and compare — this is what answers "why node 15, and would another node behave differently?"

Monitors are two per zone: new 107/113 | old 15/145 | average 209/231. A constant offset is added to
one monitor's column, the calibration is re-run from the cached predictions, and the displacement of
every coefficient is reported in units of its own unbiased posterior SD, together with what happens
to the risk field.

PRIMARY: formal censored likelihood. COMPARATOR: informal GLUE at the primary threshold.
Both signs are swept. Because observations are censored at the sensor floor, a negative offset is
not the mirror image of a positive one — it pushes readings onto the floor, where they carry
different information — so a one-sided sweep would report half of the behaviour as if it were all
of it. 30 noise realisations per cell; medians reported.
"""
import os
import json
import numpy as np
from scipy.stats import spearmanr
import wq_common as B

HERE = os.path.dirname(os.path.abspath(__file__))
ZKEYS = ["old", "average", "new"]
ESS_MIN = 100          # the same threshold Section 2.5.1 declares and step6 applies
ZONE_OF = {"107": "new", "113": "new", "15": "old", "145": "old", "209": "average", "231": "average"}

cache = np.load(os.path.join(HERE, "baseline_cache", "baseline.npz"), allow_pickle=True)
C_all = cache["C_all"].astype(np.float64)
ALL_NODES = list(cache["all_nodes"])
mon_pos = list(cache["mon_pos"])
truth_mon = cache["truth_all"][:, mon_pos]
S = {"old": cache["S_old"], "average": cache["S_avg"], "new": cache["S_new"]}
C_all_mon = C_all[:, :, mon_pos]
C_MIN = 0.2
DUR = B.below_threshold_duration(C_all, C_MIN)
TW = B.window_hours(C_all)
# see step8/step8b: P_bar and E[A] can rank the leading nodes differently, so both are carried
DEFC = np.trapezoid(np.clip(C_MIN - C_all, 0.0, None), dx=1.0, axis=1).astype(np.float64)
OFFSETS = [-0.10, -0.05, 0.05, 0.10]
SEEDS = list(range(42, 72))
SCHEMES = [B.PRIMARY_WEIGHTING, "informal_glue"]
TOP_K = 6


def run(bcol, off):
    """Median coefficient means/SDs, median risk field and censoring count over the noise seeds."""
    out = {s: {"means": {z: [] for z in ZKEYS}, "sds": {z: [] for z in ZKEYS}, "P": [], "A": [],
               "ess": [], "maxw": []}
           for s in SCHEMES}
    n_clipped = []
    rmse_mins, noise_rmses = [], []
    for seed in SEEDS:
        rng = np.random.default_rng(seed)
        eps = rng.normal(0, B.SIGMA_OBS, truth_mon.shape)
        obs = truth_mon + eps
        if bcol is not None:
            obs[:, bcol] += off
        n_clipped.append(int((obs[B.WARMUP_H:] < 0).sum()))
        obs = np.clip(obs, 0, None)[B.WARMUP_H:]
        # Aggregate fit against the noise floor, for every arm rather than for node 15 alone
        # (Step 8). The paper pairs the largest displacement, which comes from this sweep, with
        # a statement about the residual; that pairing is only honest if the residual is
        # measured over the same 24 arms. rmse_min is the best RMSE any candidate achieves
        # against the biased observations; noise_rmse is what the unbiased observations
        # themselves score against the truth.
        obs_unbiased = np.clip(truth_mon + eps, 0, None)[B.WARMUP_H:]
        rmse_mins.append(float(B.rmse_of(C_all_mon, obs).min()))
        noise_rmses.append(float(np.sqrt(((obs_unbiased - truth_mon[B.WARMUP_H:]) ** 2).mean())))
        wts = B.all_weightings(C_all_mon, obs, threshold=B.RMSE_THR, schemes=SCHEMES)
        for s in SCHEMES:
            # the diagnostics used to be discarded here. Section 2.5.1 declares median ESS < 100 a
            # sampling-limited row and Section 3.3.1 applies it to the noise sweep, so the same
            # criterion has to be checkable for these arms; without it a displacement quoted from a
            # concentrated arm cannot be told apart from one the library actually resolves.
            w, diag = wts[s]
            if w is None:
                continue
            out[s]["ess"].append(float(diag["ess"]))
            out[s]["maxw"].append(float(w.max()))
            for z in ZKEYS:
                m, sd = B.weighted_mean_sd(w, S[z])
                out[s]["means"][z].append(m)
                out[s]["sds"][z].append(sd)
            out[s]["P"].append(np.tensordot(w, DUR, axes=(0, 0)) / TW)
            out[s]["A"].append(np.tensordot(w, DEFC, axes=(0, 0)))
    res = {}
    for s in SCHEMES:
        P_med = np.median(np.vstack(out[s]["P"]), axis=0)
        A_med = np.median(np.vstack(out[s]["A"]), axis=0)
        ess = np.asarray(out[s]["ess"], float)
        res[s] = {"means": {z: float(np.median(out[s]["means"][z])) for z in ZKEYS},
                  "sds": {z: float(np.median(out[s]["sds"][z])) for z in ZKEYS},
                  "ess_med": float(np.median(ess)), "ess_p5": float(np.percentile(ess, 5)),
                  "ess_min": float(ess.min()), "maxw_med": float(np.median(out[s]["maxw"])),
                  "sampling_limited": bool(np.median(ess) < ESS_MIN),
                  "P": P_med, "A": A_med,
                  "top": [ALL_NODES[i] for i in np.argsort(P_med)[::-1][:TOP_K]],
                  "top_deficit": [ALL_NODES[i] for i in np.argsort(A_med)[::-1][:TOP_K]]}
    fit = {"rmse_min_med": float(np.median(rmse_mins)),
           "noise_rmse_med": float(np.median(noise_rmses)),
           "rmse_min_over_noise_med": float(np.median(np.array(rmse_mins)
                                                      / np.array(noise_rmses)))}
    return res, float(np.median(n_clipped)), fit


base, base_clip, base_fit = run(None, 0.0)

print("=== Step 8c: sensor-bias location sweep (both signs, 30 noise realisations) ===")
for s in SCHEMES:
    tag = "PRIMARY" if s == B.PRIMARY_WEIGHTING else "comparator"
    print(f"{s:>16} ({tag}) unbiased means: " +
          ", ".join(f"{z} {base[s]['means'][z]:+.4f}" for z in ZKEYS))
    print(f"{'':>16}       unbiased SDs  : " +
          ", ".join(f"{z} {base[s]['sds'][z]:.4f}" for z in ZKEYS))
print(f"median censored observations with no offset: {base_clip:.0f} of {6 * (C_all.shape[1])}\n")

report = {**B.weighting_provenance(comparators=["informal_glue"]),
          "informal_threshold": B.RMSE_THR, "offsets": OFFSETS, "n_noise": len(SEEDS),
          "offset_signs": "two-sided; censoring at the sensor floor makes the response asymmetric",
          "risk_ranking_basis": f"median risk field over the {len(SEEDS)} noise realisations",
          "risk_metric": {
              "P_bar": "expected fraction of the 48 h window below 0.2 mg/L (risk_* columns)",
              "E_A": "expected cumulative deficit in mg/L*h (deficit_* columns); this is the metric "
                     "Step 10's headline top-10 is ranked by, and Step 8b showed the two can "
                     "disagree about which nodes lead"},
          "baseline": {s: {"means": base[s]["means"], "sds": base[s]["sds"],
                           f"top{TOP_K}": base[s]["top"]} for s in SCHEMES},
          "baseline_censored_med": base_clip, "baseline_fit": base_fit, "rows": []}

for s in SCHEMES:
    tag = "PRIMARY" if s == B.PRIMARY_WEIGHTING else "comparator"
    print(f"--- {s} ({tag}) ---")
    print(f"{'node':>5} {'zone':>8} {'offset':>7} | {'Δold/SD':>8} {'Δavg/SD':>8} {'Δnew/SD':>8} | "
          f"{'own/SD':>7} | {'cens':>5} {'rho':>6} {'J6':>5}")
    for off in OFFSETS:
        for node in B.MONITOR_NODES:
            bcol = B.MONITOR_NODES.index(node)
            res, clip, fit = run(bcol, off)
            r = res[s]
            d = {z: r["means"][z] - base[s]["means"][z] for z in ZKEYS}
            dsd = {z: d[z] / base[s]["sds"][z] for z in ZKEYS}
            own = ZONE_OF[node]
            rho = float(spearmanr(r["P"], base[s]["P"]).statistic)
            jac = (len(set(r["top"]) & set(base[s]["top"]))
                   / len(set(r["top"]) | set(base[s]["top"])))
            rho_A = float(spearmanr(r["A"], base[s]["A"]).statistic)
            jac_A = (len(set(r["top_deficit"]) & set(base[s]["top_deficit"]))
                     / len(set(r["top_deficit"]) | set(base[s]["top_deficit"])))
            report["rows"].append({
                "scheme": s, "node": node, "zone": own, "offset": off,
                "d_old": d["old"], "d_avg": d["average"], "d_new": d["new"],
                "d_old_over_sd": dsd["old"], "d_avg_over_sd": dsd["average"],
                "d_new_over_sd": dsd["new"], "own_shift_over_sd": dsd[own],
                "ess_med": r["ess_med"], "ess_p5": r["ess_p5"], "ess_min": r["ess_min"],
                "maxw_med": r["maxw_med"], "sampling_limited": r["sampling_limited"],
                "sd_same_arm": r["sds"],
                "n_censored_med": clip, **fit, "risk_spearman_vs_unbiased": rho,
                f"risk_top{TOP_K}_jaccard_vs_unbiased": jac, f"top{TOP_K}": r["top"],
                "deficit_spearman_vs_unbiased": rho_A,
                f"deficit_top{TOP_K}_jaccard_vs_unbiased": jac_A,
                f"deficit_top{TOP_K}": r["top_deficit"]})
            print(f"{node:>5} {own:>8} {off:>+7.3f} | {dsd['old']:>+8.2f} {dsd['average']:>+8.2f} "
                  f"{dsd['new']:>+8.2f} | {dsd[own]:>+7.2f} | {clip:>5.0f} {rho:>6.3f} {jac:>5.2f}")
    print()

prim = [r for r in report["rows"] if r["scheme"] == B.PRIMARY_WEIGHTING]
worst = max(prim, key=lambda r: abs(r["own_shift_over_sd"]))
report["summary"] = {
    "max_own_shift_over_sd": {"node": worst["node"], "zone": worst["zone"],
                              "offset": worst["offset"],
                              "value": worst["own_shift_over_sd"]},
    # the paper quotes this as the bound on how far the aggregate residual rises under any bias
    # arm, so it is stored rather than left to be re-derived from the rows
    "max_rmse_min_over_noise": max(r["rmse_min_over_noise_med"] for r in prim),
    # the same bound as a percentage excess, which is the form the paper states it in
    "max_residual_excess_pct": 100 * (max(r["rmse_min_over_noise_med"] for r in prim) - 1),
    # how many arms fail the resolution criterion the study sets for itself, and which they are.
    # Reported rather than dropped, exactly as the sampling-limited noise rows are in Section 3.3.1.
    "n_arms_sampling_limited": sum(1 for r in prim if r["sampling_limited"]),
    "sampling_limited_arms": [f"{r['node']}@{r['offset']:+.2f}" for r in prim
                              if r["sampling_limited"]],
    "ess_med_min_over_arms": min(r["ess_med"] for r in prim),
    "ess_criterion": f"median ESS < {ESS_MIN}, the criterion of Section 2.5.1",
    "min_risk_spearman": min(r["risk_spearman_vs_unbiased"] for r in prim),
    f"min_risk_top{TOP_K}_jaccard": min(r[f"risk_top{TOP_K}_jaccard_vs_unbiased"] for r in prim),
    # The paper states how many of the 24 arms move the duration shortlist. That count was read
    # off the printed rows and then written into the prose, which left it with no path an
    # automated check could bind to; it is stored here so the claim registry can hold it.
    "n_arms_altering_risk_top6": sum(1 for r in prim
                                     if r[f"risk_top{TOP_K}_jaccard_vs_unbiased"] < 1.0),
    "n_arms_altering_deficit_top6": sum(1 for r in prim
                                        if r[f"deficit_top{TOP_K}_jaccard_vs_unbiased"] < 1.0),
    # the paper states the deficit side as the complement ("unchanged in 23 of 24"), so store the
    # complement rather than leave the sentence to be checked against a number it does not contain
    "n_arms_holding_deficit_top6": sum(1 for r in prim
                                       if r[f"deficit_top{TOP_K}_jaccard_vs_unbiased"] >= 1.0),
    "n_arms": len(prim),
    "min_deficit_spearman": min(r["deficit_spearman_vs_unbiased"] for r in prim),
    f"min_deficit_top{TOP_K}_jaccard": min(r[f"deficit_top{TOP_K}_jaccard_vs_unbiased"]
                                           for r in prim),
    "asymmetry_note": "compare the +0.05 and -0.05 rows of the same node: the censoring count "
                      "differs, so the two are not mirror images"}
print(f"largest own-coefficient displacement under the primary rule: node {worst['node']} "
      f"({worst['zone']}) at {worst['offset']:+.3f} -> {worst['own_shift_over_sd']:+.2f} posterior SD")
print(f"across every cell the 92-node risk ranking keeps Spearman >= "
      f"{report['summary']['min_risk_spearman']:.3f} and top-{TOP_K} Jaccard >= "
      f"{report['summary'][f'min_risk_top{TOP_K}_jaccard']:.2f}: the coefficients move far more "
      f"than the operational ranking does.")


def _jsafe(o):
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.integer):
        return int(o)
    return str(o)


with open(os.path.join(HERE, "baseline_cache", "step8c_bias_bynode.json"), "w") as f:
    json.dump(report, f, indent=2, default=_jsafe)
print("saved step8c_bias_bynode.json")
