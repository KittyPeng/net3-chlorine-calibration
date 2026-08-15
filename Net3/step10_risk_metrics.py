"""Step 10: operational risk metrics (duration + depth) and water-age corroboration.

A single probability P(C < C_MIN) conflates *how long* and *how far below* the threshold a node sits.
We therefore report, from the behavioural ensemble (weights recomputed at the primary threshold
0.107): the original probability re-expressed as a duration, two genuine severity metrics (minimum
concentration and cumulative deficit), and one reaction-independent hydraulic diagnostic (water age).

IMPORTANT — time axis: the post-warm-up record is t = WARMUP_H .. DURATION_H, which is one MORE
reporting point than it is hours (49 points, 48 one-hour intervals). Durations and cumulative
deficits are therefore trapezoidally integrated over the intervals, not summed over the points.
All window bounds are derived from the cached array shape, so they follow wq_common.

Per behavioural member i and node n:
    D_{i,n} = ∫ 1[C<C_MIN] dt        (below-threshold duration, trapezoid, h)
    A_{i,n} = ∫ max(0, C_MIN−C) dt   (cumulative deficit, trapezoid, mg/L·h)
    M_{i,n} = min_t C_{i,n}(t)        (minimum concentration; per-member min, then weighted — Method A)
Reported as weighted expectation AND weighted 5/50/95 % quantiles across the ensemble.
"""
import os
import json
import numpy as np
from scipy.stats import spearmanr, pearsonr
import wntr
from wntr.metrics.hydraulic import average_expected_demand
import wq_common as B

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures")
C_MIN = 0.2

cache = np.load(os.path.join(HERE, "baseline_cache", "baseline.npz"), allow_pickle=True)
C_all = cache["C_all"].astype(np.float64)                 # (N_MC, n_t, 92)
ALL_NODES = list(cache["all_nodes"])
RMSE = cache["RMSE"]
sigma = B.SIGMA_OBS
n_t = C_all.shape[1]                                       # reporting points in the window
T_WINDOW = n_t - 1                                         # one-hour intervals

# weighting: the formal censored likelihood is primary; the informal GLUE score is carried as a
# comparator so the question "is the operational pattern robust to the inference convention?" can be
# answered rather than assumed.
SCHEMES = {
    "formal_censored": (cache["loglik_censored"], None),
    "informal_glue": (B.glue_score(RMSE), RMSE < B.RMSE_THR),
}
weights, diags = {}, {}
for sname, (ll, mask) in SCHEMES.items():
    weights[sname], diags[sname] = B.weights_from_loglik(ll, mask)
w = weights[B.PRIMARY_WEIGHTING]
print(f"primary weighting: {B.PRIMARY_WEIGHTING} (ESS {diags[B.PRIMARY_WEIGHTING]['ess']:.1f} of "
      f"{len(RMSE)})")
print(f"comparator: informal GLUE at threshold {B.RMSE_THR} "
      f"({(RMSE < B.RMSE_THR).sum()}/{len(RMSE)} behavioural, ESS "
      f"{diags['informal_glue']['ess']:.1f})")
print(f"time window: {n_t} reporting points = {T_WINDOW} h "
      f"({B.WARMUP_H}->{B.DURATION_H})")

# ---- per-member metrics (trapezoidal over the one-hour intervals) ----
trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))
below = (C_all < C_MIN).astype(np.float64)
D = trapz(below, dx=1.0, axis=1)                           # (N_MC, 92) hours below, max T_WINDOW
deficit = np.maximum(0.0, C_MIN - C_all)
A = trapz(deficit, dx=1.0, axis=1)                         # (N_MC, 92) mg/L·h
M = C_all.min(axis=1)                                      # (N_MC, 92) per-member minimum

Dbar = w @ D                                               # expected duration
Abar = w @ A                                               # expected cumulative deficit
Mbar = w @ M
frac = Dbar / T_WINDOW                                     # = original time-averaged P(C<C_MIN)


def wquantile(vals, weights, qs):
    o = np.argsort(vals)
    v, cw = vals[o], np.cumsum(weights[o])
    cw = (cw - 0.5 * weights[o]) / weights.sum()           # Hazen plotting positions
    return [float(np.interp(q, cw, v)) for q in qs]


# ---- water age (reaction-independent) with time-window diagnostics ----
wn = wntr.network.WaterNetworkModel(B.NET3_INP)
wn.options.time.duration = B.DURATION_H * 3600
wn.options.time.hydraulic_timestep = B.HYDRAULIC_TIMESTEP_S
wn.options.time.report_timestep = B.REPORT_TIMESTEP_S
wn.options.time.quality_timestep = B.QUALITY_TIMESTEP_S
wn.options.quality.parameter = "AGE"
age_df = wntr.sim.EpanetSimulator(wn).run_sim().node["quality"][ALL_NODES] / 3600.0   # (169, 92) h
age_full = age_df.values                        # hourly reporting from 0 to 168 h
age_post = age_full[B.WARMUP_H:]                # post-warm-up window: 120..168 h, 49 reporting points
mean_age_post = age_post.mean(axis=0)           # primary definition: mean over the 120-168 h window
mean_age_last24 = age_full[B.DURATION_H - 24:].mean(axis=0)   # last 24 h interval: 144..168 h
# NOT a steady-state age: Step 0 shows the age field is still filling at 168 h and would need ~600 h,
# which the model's 168 h horizon cannot supply. The window is a choice and is recorded as one.
final_age = age_full[-1]                        # single final reporting hour, t = 168 h

T_A, T_B, T_C = B.WARMUP_H, (B.WARMUP_H + B.DURATION_H) // 2, B.DURATION_H
print("\nwater-age window diagnostics for the leading nodes:")
print(f"{'node':>5} | {'t=' + str(T_A):>6} {'t=' + str(T_B):>6} {'t=' + str(T_C):>6} | "
      f"{'mean' + str(T_A) + '-' + str(T_C):>11} "
      f"{'mean' + str(B.DURATION_H - 24) + '-' + str(T_C):>11}")
for node in ["243", "131", "166", "141"]:
    j = ALL_NODES.index(node)
    print(f"{node:>5} | {age_full[T_A, j]:6.1f} {age_full[T_B, j]:6.1f} {age_full[T_C, j]:6.1f} | "
          f"{mean_age_post[j]:11.1f} {mean_age_last24[j]:11.1f}")

# ---- association between water age and risk: DESCRIPTIVE ONLY ----
# No p-value and no ordinary bootstrap are reported. Both assume the 92 junctions are independent
# samples, and they are not: neighbouring junctions share pipes, tanks and the same demand patterns,
# so the effective number of independent units is a small fraction of 92. A p-value computed as if
# n = 92 was previously quoted as ~1e-16, which is a statement about the wrong null hypothesis rather
# than a small one. The rank correlation itself is still a legitimate DESCRIPTION of the pattern, and
# a spatial-BLOCK bootstrap is reported beside it as a deliberately conservative width: whole
# coordinate blocks of junctions are resampled together, so within-block dependence is preserved.
n_nodes = len(ALL_NODES)
rho_dur = float(spearmanr(mean_age_post, Dbar).statistic)
rho_def = float(spearmanr(mean_age_post, Abar).statistic)
r_dur = float(pearsonr(mean_age_post, Dbar).statistic)

N_BLOCKS = 10
coords = np.array([wn.get_node(n).coordinates for n in ALL_NODES], dtype=float)
cz = (coords - coords.mean(axis=0)) / coords.std(axis=0)
rng = np.random.default_rng(0)
centres = cz[rng.choice(n_nodes, N_BLOCKS, replace=False)]
for _ in range(50):                        # Lloyd's algorithm; spatial blocks, not a claim about clusters
    lab = np.argmin(((cz[:, None, :] - centres[None]) ** 2).sum(axis=2), axis=1)
    for b in range(N_BLOCKS):
        if np.any(lab == b):
            centres[b] = cz[lab == b].mean(axis=0)
blocks = [np.where(lab == b)[0] for b in range(N_BLOCKS) if np.any(lab == b)]
boot = []
for _ in range(2000):
    pick = np.concatenate([blocks[i] for i in rng.integers(0, len(blocks), len(blocks))])
    s = spearmanr(mean_age_post[pick], Dbar[pick]).statistic
    if np.isfinite(s):
        boot.append(s)
ci_block = (float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5)))
print(f"\nwater age vs risk over all n={n_nodes} junctions — DESCRIPTIVE association:")
print(f"  duration : Spearman {rho_dur:.2f}   Pearson {r_dur:.2f}")
print(f"  deficit  : Spearman {rho_def:.2f}")
print(f"  spatial-block bootstrap 95% (whole coordinate blocks resampled, {len(blocks)} blocks): "
      f"[{ci_block[0]:.2f}, {ci_block[1]:.2f}]")
print("No p-value is reported: the junctions are not independent samples, so a p-value computed at")
print("n = 92 tests a null hypothesis that does not apply. The block interval is conservative and is")
print("a width, not a significance test.")

# ---- top nodes by cumulative deficit, with 5-95% ensemble bands ----
order = np.argsort(Abar)[::-1]
print(f"\n=== top-10 risk nodes (ranked by expected cumulative deficit) ===")
print(f"{'node':>5} | {'dur E[5-95]':>18} | {'deficit E[5-95]':>20} | {'minC 50[5-95]':>18} | {'age':>5}")
top = []
for i in order[:10]:
    node = ALL_NODES[i]
    dq = wquantile(D[:, i], w, [0.05, 0.95])
    aq = wquantile(A[:, i], w, [0.05, 0.95])
    mq = wquantile(M[:, i], w, [0.05, 0.50, 0.95])
    top.append({"node": node, "dur_mean": float(Dbar[i]), "dur_5_95": dq,
                "deficit_mean": float(Abar[i]), "deficit_5_95": aq,
                "minC_5_50_95": mq, "frac": float(frac[i]), "mean_age_h": float(mean_age_post[i]),
                "mean_age_last24_h": float(mean_age_last24[i]),
                "age_final_h": float(final_age[i])})
    print(f"{node:>5} | {Dbar[i]:5.1f} [{dq[0]:4.1f},{dq[1]:4.1f}] | "
          f"{Abar[i]:5.2f} [{aq[0]:4.2f},{aq[1]:4.2f}] | {mq[1]:.3f}[{mq[0]:.2f},{mq[2]:.2f}] | {mean_age_post[i]:5.1f}")

# ---- three network averages, because they answer different questions ----
# An unweighted mean over 92 junctions counts a zero-demand node the same as the largest consumer.
# That is the right average for "how is the network doing hydraulically" and the wrong one for
# "how much service is affected", so all three are reported instead of one being chosen silently.
# base_demand is the BASE value only. Four Net3 junctions (15, 35, 123, 203) encode their real
# demand as 1 GPM times a large time-varying pattern (means 264 / 1718 / 1194 / 4512), so reading the
# base alone makes them 0.063 L/s each instead of 16.7 / 108.4 / 75.3 / 284.6 - together 70% of the
# network's actual demand, counted as four of its smallest users. average_expected_demand applies the
# pattern and the global multiplier, which is what a demand weight has to be.
dem = average_expected_demand(wn)[ALL_NODES].values * 1000.0                      # L/s
consumer = dem > 0
net_avgs = {}
for label, metric in (("E_duration_h", Dbar), ("E_deficit", Abar), ("min_C", Mbar)):
    net_avgs[label] = {
        "unweighted_all_junctions": float(metric.mean()),
        "consumer_only": float(metric[consumer].mean()),
        "demand_weighted": float(np.sum(metric * dem) / dem.sum()),
    }
print(f"\n=== network averages ({int(consumer.sum())} of {n_nodes} junctions have non-zero demand, "
      f"total {dem.sum():.1f} L/s) ===")
print(f"{'metric':>14} | {'unweighted':>11} {'consumer-only':>14} {'demand-weighted':>16}")
for label, v in net_avgs.items():
    print(f"{label:>14} | {v['unweighted_all_junctions']:11.4f} {v['consumer_only']:14.4f} "
          f"{v['demand_weighted']:16.4f}")
print("Demand-weighting is the service-relevant average; it is lower here than the unweighted mean,")
print("so quoting the unweighted figure is the conservative choice rather than a flattering one.")

# ---- is the operational pattern robust to the inference convention? ----
TOP_K = 10
xs = {}
print(f"\n=== robustness of the risk pattern to the weighting scheme ===")
print(f"{'scheme':>16} {'ESS':>8} | {'net-mean E[D] h':>15} {'net-mean E[A]':>13} | "
      f"top-{TOP_K} overlap with primary")
primary_top = None
for sname, ws in weights.items():
    dbar_s, abar_s = ws @ D, ws @ A
    top_s = [ALL_NODES[i] for i in np.argsort(abar_s)[::-1][:TOP_K]]
    if primary_top is None:
        primary_top = top_s
    jac = len(set(top_s) & set(primary_top)) / len(set(top_s) | set(primary_top))
    # Per-node 5-95% duration bands for the two schemes, not just the network means: the claim that
    # the formal rule tightens the risk bands is a NODE-level claim and needs a node-level number.
    j15 = ALL_NODES.index("15")
    band15 = B.weighted_quantile(D[:, j15], ws, [0.05, 0.95])
    xs[sname] = {"ess": diags[sname]["ess"],
                 "net_mean_E_duration_h": float(dbar_s.mean()),
                 "net_mean_E_deficit": float(abar_s.mean()),
                 "node15_E_duration_h": float(dbar_s[j15]),
                 "node15_duration_5_95": [float(band15[0]), float(band15[1])],
                 f"top{TOP_K}": top_s,
                 f"top{TOP_K}_jaccard_vs_primary": jac}
    print(f"{sname:>16} {diags[sname]['ess']:8.1f} | {dbar_s.mean():15.3f} {abar_s.mean():13.4f} | "
          f"{jac:.2f}  {'(primary)' if sname == B.PRIMARY_WEIGHTING else ''}")

# ---- does hourly output miss short excursions? ----
# The water-quality solver steps at QUALITY_TIMESTEP_S (300 s) but the risk metrics integrate the
# HOURLY report, so a dip that starts and ends inside one hour is invisible. Trapezoidal integration
# of a binary indicator also places any threshold crossing at the midpoint between two reports.
# Re-running the highest-weight members at finer reporting resolution measures the error directly.
REPORT_STEPS_S = [3600, 900, 300]
TS_MEMBERS = 12                                    # highest-weight members; enough for a bias check
ts_idx = np.argsort(w)[::-1][:TS_MEMBERS]
ts_w = w[ts_idx] / w[ts_idx].sum()
print(f"\n=== reporting-resolution sensitivity ({TS_MEMBERS} highest-weight members) ===")
print(f"{'report step':>12} {'points':>7} | {'net-mean E[D] (h)':>18} {'net-mean E[A]':>14} | "
      f"{'top-10 Jaccard':>14}")
ts_rows, ts_top_ref = [], None
for step_s in REPORT_STEPS_S:
    D_ts = np.empty((len(ts_idx), n_nodes))
    A_ts = np.empty((len(ts_idx), n_nodes))
    for k, i in enumerate(ts_idx):
        wn_ts = B.build_model(B.KB_FIXED, 0.0,
                              pre_run=B.make_kw_hook(float(cache["S_old"][i]),
                                                     float(cache["S_avg"][i]),
                                                     float(cache["S_new"][i])))
        wn_ts.options.time.report_timestep = step_s
        res_ts = wntr.sim.EpanetSimulator(wn_ts).run_sim()
        # this path bypasses simulate_chlorine, so the kg/m^3 -> mg/L conversion is applied here
        q_ts = B.internal_to_mgl(res_ts.node["quality"][ALL_NODES])
        hours = np.asarray(q_ts.index, dtype=float) / 3600.0
        keep = hours >= B.WARMUP_H
        c_ts = q_ts.values[keep]
        dt_h = step_s / 3600.0
        D_ts[k] = trapz((c_ts < C_MIN).astype(float), dx=dt_h, axis=0)
        A_ts[k] = trapz(np.maximum(0.0, C_MIN - c_ts), dx=dt_h, axis=0)
    dbar, abar = ts_w @ D_ts, ts_w @ A_ts
    ts_top = [ALL_NODES[i] for i in np.argsort(abar)[::-1][:10]]
    if ts_top_ref is None:
        ts_top_ref = ts_top
    jac = len(set(ts_top) & set(ts_top_ref)) / len(set(ts_top) | set(ts_top_ref))
    ts_rows.append({"report_timestep_s": step_s, "points_in_window": int(c_ts.shape[0]),
                    "net_mean_E_duration_h": float(dbar.mean()),
                    "net_mean_E_deficit": float(abar.mean()),
                    "top10_jaccard_vs_hourly": jac, "top10": ts_top})
    print(f"{step_s:>12} {c_ts.shape[0]:>7} | {dbar.mean():>18.3f} {abar.mean():>14.4f} | "
          f"{jac:>14.2f}")
h, f = ts_rows[0], ts_rows[-1]
rel_d = (f["net_mean_E_duration_h"] - h["net_mean_E_duration_h"]) / h["net_mean_E_duration_h"]
rel_a = (f["net_mean_E_deficit"] - h["net_mean_E_deficit"]) / h["net_mean_E_deficit"]
print(f"hourly -> {REPORT_STEPS_S[-1]} s: E[D] {rel_d * +100:+.1f}%, E[A] {rel_a * 100:+.1f}%, "
      f"top-10 Jaccard {f['top10_jaccard_vs_hourly']:.2f}")
print("Hourly reporting is adequate if these are small: the risk metrics are then integrals of a")
print("field that is smooth on the hour scale, not a sampling of spikes.")

report = {"C_MIN": C_MIN, "threshold": B.RMSE_THR, "T_window_h": T_WINDOW, "n_nodes": n_nodes,
          "reporting_resolution_sensitivity": {
              "n_members": TS_MEMBERS, "steps_s": REPORT_STEPS_S,
              "quality_timestep_s": B.QUALITY_TIMESTEP_S, "rows": ts_rows,
              "hourly_to_finest_rel_change": {"E_duration": float(rel_d), "E_deficit": float(rel_a)}},
          "primary_weighting": B.PRIMARY_WEIGHTING,
          "weighting_diagnostics": diags,
          "network_averages": {"n_consumer_junctions": int(consumer.sum()),
                               "total_demand_L_s": float(dem.sum()), "by_metric": net_avgs},
          "cross_scheme_robustness": xs,
          "age_risk_association": {
              "status": "descriptive only",
              "spearman_dur": rho_dur, "spearman_def": rho_def, "pearson_dur": r_dur,
              "n_junctions": n_nodes,
              "p_value": None,
              "p_value_omitted_because": "the 92 junctions share pipes, tanks and demand patterns, "
                                         "so they are not independent samples; a p-value computed "
                                         "at n = 92 (previously quoted as ~1e-16) tests a null "
                                         "hypothesis that does not apply here",
              "iid_bootstrap_omitted_because": "same reason: resampling individual junctions "
                                               "assumes the independence that fails",
              "block_bootstrap": {"n_blocks": len(blocks), "block_rule": "k-means on standardised "
                                  "node coordinates; whole blocks resampled with replacement",
                                  "n_resamples": len(boot), "ci95_spearman_dur": list(ci_block),
                                  "interpretation": "a conservative width for the rank "
                                                    "correlation, not a significance test"}},
          "age_windows": {
              "mean_age_h": f"primary: mean over the post-warm-up record ({B.WARMUP_H}-"
                            f"{B.DURATION_H} h)",
              "mean_age_last24_h": f"mean over the last diurnal cycle "
                                   f"({B.DURATION_H - 24}-{B.DURATION_H} h)",
              "age_final_h": f"single final reporting hour (t = {B.DURATION_H} h)",
              "why": "water age is still rising at the high-age nodes, so the window choice moves "
                     "the number; all three are stored to keep the log verifiable"},
          "top10": top}
with open(os.path.join(HERE, "baseline_cache", "step10_risk_metrics.json"), "w") as f:
    json.dump(report, f, indent=2)
print("\nsaved step10_risk_metrics.json")

# ---- figures ----
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.5, 5))

top5 = set(ALL_NODES[i] for i in order[:5])
colors = ["crimson" if n in top5 else "steelblue" for n in ALL_NODES]
axA.scatter(mean_age_post, Dbar, c=colors, s=28, alpha=0.75)
for i in order[:5]:
    axA.annotate(ALL_NODES[i], (mean_age_post[i], Dbar[i]), textcoords="offset points",
                 xytext=(5, 4), fontsize=9, color="crimson")
axA.set_xlabel("mean water age (h), post-warm-up — reaction-independent")
axA.set_ylabel(f"expected hours below {C_MIN} mg/L (of {T_WINDOW} h)")
axA.set_title(f"(a) Risk associated with water age\n"
              f"Spearman ρ = {rho_dur:.2f}, descriptive association (n={n_nodes} dependent junctions)")
axA.grid(alpha=0.3)

top10 = order[:10]
nodes10 = [ALL_NODES[i] for i in top10]
durs = [Dbar[i] for i in top10]
err = np.array([[Dbar[i] - wquantile(D[:, i], w, [0.05])[0] for i in top10],
                [wquantile(D[:, i], w, [0.95])[0] - Dbar[i] for i in top10]])
err = np.clip(err, 0, None)
axB.barh(range(len(top10)), durs, xerr=err, color="firebrick", capsize=3,
         error_kw={"ecolor": "0.3", "lw": 1})
axB.set_yticks(range(len(top10)))
axB.set_yticklabels(nodes10)
axB.invert_yaxis()
for k, i in enumerate(top10):
    axB.annotate(f"minC={Mbar[i]:.2f}, age={mean_age_post[i]:.0f}h",
                 (durs[k], k), textcoords="offset points", xytext=(6, 0),
                 va="center", fontsize=7.5, color="0.25")
axB.axvline(T_WINDOW, color="0.6", ls=":", lw=1)
axB.set_xlabel("expected hours below 0.2 mg/L (bars = 5-95% ensemble band)")
axB.set_title("(b) Top-10 nodes: expected duration (+5-95%), depth and age")
axB.grid(alpha=0.3, axis="x")
axB.set_xlim(0, T_WINDOW * 1.45)

plt.tight_layout()
figpath = os.path.join(FIGDIR, "step10_risk_metrics.png")
plt.savefig(figpath, dpi=130, bbox_inches="tight")
print("figure saved to", figpath)
