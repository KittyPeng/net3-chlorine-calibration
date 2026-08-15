"""Step 11: held-out validation at three levels of difficulty.

The draft's actual claim is that risk extrapolates to the 86 unmonitored junctions. Leave-one-monitor-
out alone cannot support that claim: every held-out monitor has a partner in the same zone, so the
zone stays observed. Three tests are therefore run in increasing order of severity, and they give
different answers.

  A  leave-one-MONITOR-out   drop 1 of 6; the held-out monitor's zone is still observed by its partner
  B  leave-one-ZONE-out      drop BOTH monitors of a zone; that zone becomes unobserved entirely
  C  unmonitored junctions   predict at junctions that never enter any calibration, against the truth

PRIMARY: formal censored likelihood. COMPARATOR: informal GLUE, whose threshold is rescaled to the
number of monitors actually used so that dropping a monitor does not also loosen the acceptance band.
Both are reported at every level, because the coverage of a predictive band is exactly the quantity
an inflated parameter spread distorts: the informal score buys coverage by being wide, which looks
like a well-calibrated model and is not one.

Reuses the cached candidate predictions (no EPANET); 30 noise realisations, median [IQR].

Two predictive quantities that must not be confused:

  predicting an OBSERVATION at a sensor needs parameter spread AND measurement noise
      band = mean ± 1.645·sqrt(Var_ensemble + σ²)
  predicting the TRUTH at an unmonitored junction needs parameter spread only
      band = mean ± 1.645·sqrt(Var_ensemble)

Test C uses the second form, because there is no sensor there to add noise. Both forms assume
normality; a weighted predictive quantile is reported alongside, since the ensemble is skewed and
censored at zero (the review asked for exactly this).
"""
import os
import json
import numpy as np
import wq_common as B

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures")
ZKEYS = ["old", "average", "new"]
TRUE = {"old": B.KW_OLD_TRUE, "average": B.KW_AVG_TRUE, "new": B.KW_NEW_TRUE}
ZONE_OF = {"107": "new", "113": "new", "15": "old", "145": "old", "209": "average", "231": "average"}
sigma = B.SIGMA_OBS

cache = np.load(os.path.join(HERE, "baseline_cache", "baseline.npz"), allow_pickle=True)
C_all = cache["C_all"].astype(np.float64)
mon_pos = list(cache["mon_pos"])
truth_full = cache["truth_all"][:, mon_pos]               # (169, 6)
S = {"old": cache["S_old"], "average": cache["S_avg"], "new": cache["S_new"]}
C_mon = C_all[:, :, mon_pos]                              # (N_MC, 49, 6)
NMON = len(B.MONITOR_NODES)
Z = 1.645
PRIMARY = B.PRIMARY_WEIGHTING
SCHEMES = [PRIMARY, "informal_glue"]


def thr_for(nmon):
    """The comparator's acceptance band for a subset of monitors.

    It has to be rescaled: the band is the sampling spread of the objective, which depends on the
    number of residuals, so keeping 0.107 while dropping a monitor would quietly widen the set as
    well as removing information and confound the two.
    """
    return sigma * (1.0 + Z / np.sqrt(2.0 * nmon * C_mon.shape[1]))


def calibrate(obs, cols, scheme=PRIMARY):
    """Weights from a monitor subset under one weighting scheme; None if the scheme accepts nothing."""
    w, _ = B.all_weightings(C_mon[:, :, cols], obs[:, cols], threshold=thr_for(len(cols)),
                            schemes=[scheme])[scheme]
    return w


def q(a):
    return float(np.median(a)), float(np.percentile(a, 25)), float(np.percentile(a, 75))


SEEDS = list(range(42, 72))
all_cols = list(range(NMON))

# ---- full-6-monitor reference (median over noise), per scheme ----
full = {s: {z: [] for z in ZKEYS} for s in SCHEMES}
for seed in SEEDS:
    rng = np.random.default_rng(seed)
    obs = np.clip(truth_full + rng.normal(0, sigma, truth_full.shape), 0, None)[B.WARMUP_H:]
    for s in SCHEMES:
        w = calibrate(obs, all_cols, s)
        if w is None:
            continue
        for z in ZKEYS:
            full[s][z].append(float(w @ S[z]))
print("=== Step 11: leave-one-monitor-out validation (30 noise, median) ===")
for s in SCHEMES:
    print(f"full-6 reference [{s}]: " +
          "  ".join(f"{z} {np.median(full[s][z]):.3f}" for z in ZKEYS))

report = {**B.weighting_provenance(comparators=["informal_glue"]),
          "sigma": sigma, "noise_floor": sigma,
          "full6": {s: {z: float(np.median(full[s][z])) for z in ZKEYS} for s in SCHEMES},
          "rows": []}
example = {}
for s in SCHEMES:
    tag = "PRIMARY" if s == PRIMARY else "comparator"
    print(f"\n--- A: leave-one-monitor-out, {s} ({tag}) ---")
    print(f"{'held-out':>12} | {'k_old':>7} {'k_avg':>7} {'k_new':>7} | {'pred RMSE@m':>11} | "
          f"{'90% cov':>7}")
    for m in range(NMON):
        node = B.MONITOR_NODES[m]
        cols = [c for c in all_cols if c != m]
        kw = {z: [] for z in ZKEYS}
        pred_rmse, cover = [], []
        for seed in SEEDS:
            rng = np.random.default_rng(seed)
            obs = np.clip(truth_full + rng.normal(0, sigma, truth_full.shape), 0, None)[B.WARMUP_H:]
            w = calibrate(obs, cols, s)
            if w is None:
                continue
            for z in ZKEYS:
                kw[z].append(float(w @ S[z]))
            pm = w @ C_mon[:, :, m]                        # predicted mean at the held-out monitor
            pv = w @ (C_mon[:, :, m] - pm[None]) ** 2      # ensemble variance
            psd = np.sqrt(pv + sigma ** 2)                 # + observation noise
            lo, hi = pm - Z * psd, pm + Z * psd
            o = obs[:, m]
            pred_rmse.append(float(np.sqrt(((pm - o) ** 2).mean())))
            cover.append(float(((o >= lo) & (o <= hi)).mean()))
            if seed == 42 and s == PRIMARY:
                example[node] = {"t": np.arange(B.WARMUP_H, B.DURATION_H + 1), "pm": pm,
                                 "lo": lo, "hi": hi, "obs": o}
        row = {"scheme": s, "node": node, "zone": ZONE_OF[node],
               "k_old": q(kw["old"]), "k_avg": q(kw["average"]), "k_new": q(kw["new"]),
               "pred_rmse": q(pred_rmse), "coverage90": q(cover)}
        report["rows"].append(row)
        print(f"{node + ' (' + ZONE_OF[node] + ')':>12} | {row['k_old'][0]:7.3f} "
              f"{row['k_avg'][0]:7.3f} {row['k_new'][0]:7.3f} | {row['pred_rmse'][0]:11.3f} | "
              f"{row['coverage90'][0]:6.2f}")

print(f"\nnoise floor σ = {sigma}; pred RMSE ≈ σ ⇒ the held-out monitor is predicted as well as noise "
      f"allows.\nBut its zone was still observed by the partner monitor, so this is the easy case.")
print("Compare the two coverage columns: the comparator reaches nominal coverage with a much wider")
print("band, so coverage alone cannot tell a calibrated model from an over-dispersed one.")

# ================= B: leave-one-ZONE-out =================
# Dropping both monitors of a zone removes that zone from the observations entirely. This is the
# test that speaks to the actual claim, and it is the one the draft was missing.
ZONE_COLS = {z: [i for i, n in enumerate(B.MONITOR_NODES) if ZONE_OF[n] == z] for z in ZKEYS}
print(f"\n=== B: leave-one-ZONE-out (both monitors of a zone dropped) ===")
print("Both weightings are run, so this differs from test D only in the truth being homogeneous.")
print(f"{'scheme':>16} {'zone dropped':>13} {'mons':>9} | {'k_old':>8} {'k_avg':>8} {'k_new':>8} | "
      f"{'own coef err':>12} {'own SD ret':>10} | {'pred RMSE':>9} {'90% cov':>7}")
zone_rows = []
for scheme in SCHEMES:
    for z in ZKEYS:
        drop = ZONE_COLS[z]
        cols = [c for c in all_cols if c not in drop]
        kw = {k: [] for k in ZKEYS}
        sd_ret, pred_rmse_z, cover_z = [], [], []
        prior_sd_z = (B.PRIOR[{"old": "old", "average": "avg", "new": "new"}[z]][1]
                      - B.PRIOR[{"old": "old", "average": "avg", "new": "new"}[z]][0]) / np.sqrt(12)
        for seed in SEEDS:
            rng = np.random.default_rng(seed)
            obs = np.clip(truth_full + rng.normal(0, sigma, truth_full.shape), 0, None)[B.WARMUP_H:]
            w = calibrate(obs, cols, scheme)
            if w is None:
                continue
            for k in ZKEYS:
                kw[k].append(float(w @ S[k]))
            m_own = float(w @ S[z])
            sd_ret.append(float(np.sqrt(w @ (S[z] - m_own) ** 2) / prior_sd_z))
            for mcol in drop:
                pm = w @ C_mon[:, :, mcol]
                pv = w @ (C_mon[:, :, mcol] - pm[None]) ** 2
                psd = np.sqrt(pv + sigma ** 2)
                o = obs[:, mcol]
                pred_rmse_z.append(float(np.sqrt(((pm - o) ** 2).mean())))
                cover_z.append(float(((o >= pm - Z * psd) & (o <= pm + Z * psd)).mean()))
        row = {"scheme": scheme, "zone_dropped": z,
               "monitors_dropped": [B.MONITOR_NODES[c] for c in drop],
               "k_old": q(kw["old"]), "k_avg": q(kw["average"]), "k_new": q(kw["new"]),
               "own_coef_error": float(np.median(kw[z]) - TRUE[z]),
               "own_sd_retained": q(sd_ret), "pred_rmse": q(pred_rmse_z),
               "coverage90": q(cover_z)}
        zone_rows.append(row)
        print(f"{scheme:>16} {z:>13} {','.join(row['monitors_dropped']):>9} | {row['k_old'][0]:8.3f} "
              f"{row['k_avg'][0]:8.3f} {row['k_new'][0]:8.3f} | {row['own_coef_error']:+12.3f} "
              f"{row['own_sd_retained'][0] * 100:9.0f}% | {row['pred_rmse'][0]:9.3f} "
              f"{row['coverage90'][0]:6.2f}")
report["leave_one_zone_out"] = zone_rows
print("`own coef err` is the dropped zone's coefficient error and `own SD ret` its posterior width")
print("as a fraction of the prior: together they say whether an unobserved zone is still constrained.")

# ================= C: unmonitored-junction validation =================
# Junctions that never enter any calibration, predicted against the noise-free truth. Predicting a
# truth value needs the parameter spread only, with no measurement-noise term added.
ALL_NODES = [str(n) for n in cache["all_nodes"]]
truth_all_post = cache["truth_all"][B.WARMUP_H:]                    # (Tn, 92)
unmon_idx = [i for i, n in enumerate(ALL_NODES) if n not in B.MONITOR_NODES]
rng_pick = np.random.default_rng(7)
VAL_N = 20
val_idx = sorted(rng_pick.choice(unmon_idx, size=VAL_N, replace=False).tolist())


def wquantile_cols(vals, w, qs):
    """Weighted quantiles down axis 0 for a (members, T) block."""
    out = np.empty((len(qs), vals.shape[1]))
    order = np.argsort(vals, axis=0)
    for t in range(vals.shape[1]):
        o = order[:, t]
        cw = np.cumsum(w[o])
        cw = (cw - 0.5 * w[o]) / w.sum()
        out[:, t] = np.interp(qs, cw, vals[o, t])
    return out


print(f"\n=== C: {VAL_N} unmonitored junctions, predicted against the noise-free truth ===")
print(f"{'weighting':>10} | {'90% cov (normal)':>17} | {'90% cov (quantile)':>19} | "
      f"{'mean |err|/|truth|':>18}")
unmon = {}
for scheme in SCHEMES:
    norm_cov, quant_cov, rel_err = [], [], []
    for seed in SEEDS:
        rng = np.random.default_rng(seed)
        obs = np.clip(truth_full + rng.normal(0, sigma, truth_full.shape), 0, None)[B.WARMUP_H:]
        w = calibrate(obs, all_cols, scheme)
        if w is None:
            continue
        for j in val_idx:
            block = C_all[:, :, j]                                  # (members, Tn)
            pm = w @ block
            pv = w @ (block - pm[None]) ** 2
            # no + sigma**2 here: the target is the TRUTH, not a noisy reading
            psd = np.sqrt(pv)
            tj = truth_all_post[:, j]
            norm_cov.append(float(((tj >= pm - Z * psd) & (tj <= pm + Z * psd)).mean()))
            qlo, qhi = wquantile_cols(block, w, [0.05, 0.95])
            quant_cov.append(float(((tj >= qlo) & (tj <= qhi)).mean()))
            denom = max(float(np.abs(tj).mean()), 1e-9)
            rel_err.append(float(np.abs(pm - tj).mean() / denom))
    unmon[scheme] = {"coverage90_normal_approx": q(norm_cov),
                     "coverage90_weighted_quantile": q(quant_cov),
                     "mean_abs_rel_error": q(rel_err)}
    print(f"{scheme:>10} | {np.median(norm_cov):17.2f} | {np.median(quant_cov):19.2f} | "
          f"{np.median(rel_err):18.3f}")
report["unmonitored_validation"] = {
    "n_junctions": VAL_N, "junctions": [ALL_NODES[j] for j in val_idx],
    "target": "noise-free truth (no measurement-noise term in the band)",
    "nominal_coverage": 0.90, "by_scheme": unmon,
}
print("Nominal coverage is 0.90. Over-coverage means the band is conservative rather than")
print("calibrated, which is what an inflated parameter spread produces; under-coverage would mean")
print("the band is too narrow to be trusted. Neither scheme is externally validated here — the")
print("truth is generated by the same model, so this measures internal consistency only.")

# ================= D: the hardest case — leave-one-zone-out on a HETEROGENEOUS truth =================
# Tests A-C all use a truth generated by the same three-zone model that is fitted, so a structural
# discrepancy cannot appear. Here the truth has per-pipe heterogeneity (the Step 5c design) AND a zone
# is unobserved, so extrapolation and structural error are tested together. This is the combination
# the risk map actually relies on and the one nothing else in this log covers.
import wntr

wn_h = wntr.network.WaterNetworkModel(B.NET3_INP)
ZONE_MEAN = {"old": B.KW_OLD_TRUE, "average": B.KW_AVG_TRUE, "new": B.KW_NEW_TRUE}
zone_pipes = {z: [p for p in wn_h.pipe_name_list if B.MATERIAL_ZONES[p] == z] for z in ZKEYS}
HET_JITTER = 0.20
HET_FIELDS = 8               # a few fields, so the answer is not one arrangement
print(f"\n=== D: leave-one-zone-out on a heterogeneous truth "
      f"(+/-{HET_JITTER:.0%} per-pipe, {HET_FIELDS} fields) ===")
print(f"{'zone dropped':>13} | {'own coef err':>12} {'own SD ret':>10} | {'pred RMSE':>9} "
      f"{'90% cov':>7} | {'unmon. rel err':>14}")
het_rows = []
for z in ZKEYS:
    drop = ZONE_COLS[z]
    cols = [c for c in all_cols if c not in drop]
    errs, sdr, prmse, cov, uerr = [], [], [], [], []
    prior_sd_z = (B.PRIOR[{"old": "old", "average": "avg", "new": "new"}[z]][1]
                  - B.PRIOR[{"old": "old", "average": "avg", "new": "new"}[z]][0]) / np.sqrt(12)
    for f in range(HET_FIELDS):
        rf = np.random.default_rng(90_000 + f)
        kwp = {p: ZONE_MEAN[B.MATERIAL_ZONES[p]] * (1.0 + rf.uniform(-HET_JITTER, HET_JITTER))
               for p in wn_h.pipe_name_list}

        def hk(wn, _kw=kwp):
            for p in wn.pipe_name_list:
                wn.get_link(p).wall_coeff = B.per_day_to_per_second(_kw[p])

        t_all = B.simulate_chlorine(B.KB_FIXED, 0.0, pre_run=hk,
                                    monitor_nodes=ALL_NODES).values
        t_mon, t_post = t_all[:, mon_pos], t_all[B.WARMUP_H:]
        arith_z = float(np.array([kwp[p] for p in zone_pipes[z]]).mean())
        rr = np.random.default_rng(B.NOISE_SEED + f)
        ob = np.clip(t_mon + rr.normal(0, sigma, t_mon.shape), 0, None)[B.WARMUP_H:]
        wz = calibrate(ob, cols, PRIMARY)
        m = float(wz @ S[z])
        errs.append(m - arith_z)
        sdr.append(float(np.sqrt(wz @ (S[z] - m) ** 2) / prior_sd_z))
        for mc in drop:
            pm = wz @ C_mon[:, :, mc]
            psd = np.sqrt(wz @ (C_mon[:, :, mc] - pm[None]) ** 2 + sigma ** 2)
            o = ob[:, mc]
            prmse.append(float(np.sqrt(((pm - o) ** 2).mean())))
            cov.append(float(((o >= pm - Z * psd) & (o <= pm + Z * psd)).mean()))
        for j in val_idx:
            pmj = wz @ C_all[:, :, j]
            tj = t_post[:, j]
            uerr.append(float(np.abs(pmj - tj).mean() / max(float(np.abs(tj).mean()), 1e-9)))
    row = {"zone_dropped": z, "jitter": HET_JITTER, "n_fields": HET_FIELDS,
           "own_coef_error_vs_field_arith": q(errs), "own_sd_retained": q(sdr),
           "pred_rmse_at_dropped_monitors": q(prmse), "coverage90": q(cov),
           "unmonitored_rel_error": q(uerr)}
    het_rows.append(row)
    print(f"{z:>13} | {row['own_coef_error_vs_field_arith'][0]:+12.3f} "
          f"{row['own_sd_retained'][0] * 100:9.0f}% | {row['pred_rmse_at_dropped_monitors'][0]:9.3f} "
          f"{row['coverage90'][0]:7.2f} | {row['unmonitored_rel_error'][0]:14.3f}")
report["heterogeneous_leave_one_zone_out"] = {
    "design": "per-pipe jitter truth (Step 5c design) with one zone unobserved; formal weighting; "
              "coefficient error measured against the field's own arithmetic zone mean",
    "rows": het_rows}
print("Coefficient error is measured against each field's own arithmetic zone mean, so it is not")
print("inflated by the heterogeneity itself. Compare with test B: if the numbers are similar, the")
print("structural discrepancy adds little on top of losing the sensors.")

with open(os.path.join(HERE, "baseline_cache", "step11_loo.json"), "w") as f:
    json.dump(report, f, indent=2, default=lambda o: o.tolist() if isinstance(o, np.ndarray) else float(o))
print("saved step11_loo.json")

# ---- figure ----
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.5, 4.8))

prim_rows = [r for r in report["rows"] if r["scheme"] == PRIMARY]
nodes = [r["node"] for r in prim_rows]
rmses = [r["pred_rmse"][0] for r in prim_rows]
errlo = [r["pred_rmse"][0] - r["pred_rmse"][1] for r in prim_rows]
errhi = [r["pred_rmse"][2] - r["pred_rmse"][0] for r in prim_rows]
cols = ["firebrick" if ZONE_OF[n] == "old" else ("goldenrod" if ZONE_OF[n] == "average" else "steelblue")
        for n in nodes]
axA.bar([f"{n}\n({ZONE_OF[n]})" for n in nodes], rmses, yerr=[errlo, errhi], color=cols, capsize=3,
        error_kw={"ecolor": "0.3"})
axA.axhline(sigma, color="k", ls="--", lw=1.3, label=f"noise floor σ = {sigma}")
axA.set_ylabel("held-out prediction RMSE (mg/L)")
axA.set_title("(a) LOO out-of-sample prediction error vs noise floor\n"
              "(formal censored weighting; median, bars = IQR over 30 noise)")
axA.legend(fontsize=8)
axA.grid(alpha=0.3, axis="y")

ex = example.get("15", example[nodes[0]])
node_ex = "15" if "15" in example else nodes[0]
axB.fill_between(ex["t"], np.clip(ex["lo"], 0, None), ex["hi"], color="steelblue", alpha=0.25,
                 label="90% predictive band")
axB.plot(ex["t"], ex["pm"], color="steelblue", lw=2, label="predicted mean (held out)")
axB.plot(ex["t"], ex["obs"], "o", color="crimson", ms=3, label=f"held-out obs (node {node_ex})")
axB.set_xlabel("time (h)")
axB.set_ylabel("chlorine (mg/L)")
axB.set_title(f"(b) Predicting held-out monitor {node_ex} from the other five")
axB.legend(fontsize=8)
axB.grid(alpha=0.3)

plt.tight_layout()
figpath = os.path.join(FIGDIR, "step11_loo.png")
plt.savefig(figpath, dpi=130, bbox_inches="tight")
print("figure saved to", figpath)
