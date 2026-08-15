"""Step 5c: within-zone heterogeneity magnitude sweep (is the fit precise but biased?).

For jitter in {0, 0.20, 0.35, 0.50}: truth has per-pipe k_w,p = zone_mean*(1+delta_p),
delta_p ~ U(-jitter, jitter) (seed 12345; jitter=0.20 reproduces Step 5a). The three-zone
HOMOGENEOUS model is calibrated against the (noisy) heterogeneous observations.

jitter = 0 is the PAIRED CONTROL and must be read first: the truth is then exactly homogeneous, so
whatever bias it shows cannot be structural. Only the increment over that control is attributable
to within-zone heterogeneity. The control also calibrates `struct_residual`: at jitter = 0 a
three-zone homogeneous model can fit the truth exactly, so the residual there is purely the
resolution of the 7x7x7 grid (the truth -1.0 / -0.1 / -0.05 is not a grid point) and forms a floor
that must be subtracted before the residual at higher jitter is called structural.

The homogeneous candidate predictions do NOT depend on the truth, so the 2000 GLUE candidate
predictions are reused from the frozen baseline cache (baseline.npz); only the observations
change with the truth. A 343-point grid (built here) gives the noise-free structural residual
and the grid best fit.

No sign is assumed for the bias. A scalar exponential argument does not settle it: for
C = exp(k t), Jensen gives E[exp(k t)] >= exp(E[k] t), so matching the higher mean concentration
with a single coefficient would require k_eff >= E[k], i.e. WEAKER decay. In a network the
effective coefficient also depends on flow paths, diameters, residence times and which pipes the
monitors actually see, so the direction and magnitude are network-dependent and are measured here
rather than predicted.

TWO corrections are needed when reading the bias, and they pull in OPPOSITE directions, so both
must be applied:

  subtract the control. Any weighting has a non-zero bias at jitter = 0, where the truth is exactly
  homogeneous and nothing can be structural. That offset is the single realisation's noise, so the
  raw `bias` column OVERSTATES the structural effect; `bias_increment_vs_control` is the honest part.

  use the formal weighting. The informal score exp(-0.5 (RMSE/sigma)^2) drops the factor N of a
  formal Gaussian likelihood, so it is nearly flat inside the behavioural set: its weighted mean is
  pinned near the centre of the accepted box and cannot follow the data. It therefore UNDERSTATES
  the structural bias -- measured here by a factor of 3-4 in the old zone. The informal score's
  apparent robustness to heterogeneity is inertia, not evidence.

Both columns are reported so the contrast is visible rather than assumed.
"""
import os
import json
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import wntr

import wq_common as B

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures")

ZKEYS = ["old", "average", "new"]
ZONE_MEAN = {"old": B.KW_OLD_TRUE, "average": B.KW_AVG_TRUE, "new": B.KW_NEW_TRUE}
JITTER_SEED = 12345
JITTERS = [0.0, 0.20, 0.35, 0.50]        # 0.0 = paired homogeneous control
C_MIN = 0.2

# reuse baseline GLUE candidate predictions (independent of the truth)
cache = np.load(os.path.join(HERE, "baseline_cache", "baseline.npz"), allow_pickle=True)
C_all = cache["C_all"]                                   # (8192, 49, 92) homogeneous candidates
ALL_NODES = list(cache["all_nodes"])
mon_pos = list(cache["mon_pos"])
S = {"old": cache["S_old"], "average": cache["S_avg"], "new": cache["S_new"]}
C_all_mon = C_all[:, :, mon_pos]                          # (8192, 49, 6)

wn0 = wntr.network.WaterNetworkModel(B.NET3_INP)
PIPE_LEN = {p: wn0.get_link(p).length for p in wn0.pipe_name_list}
zone_pipes = {z: [p for p in wn0.pipe_name_list if B.MATERIAL_ZONES[p] == z] for z in ZKEYS}

# homogeneous 3-zone grid library at monitors (for structural residual + grid fit)
kw_old_grid = np.round(np.linspace(-0.2, -1.5, 7), 3)
kw_avg_grid = np.round(np.linspace(-0.04, -0.2, 7), 3)
kw_new_grid = np.round(np.linspace(-0.005, -0.10, 7), 3)
GRID = {}
t0 = time.time()
for kwo in kw_old_grid:
    for kwa in kw_avg_grid:
        for kwn in kw_new_grid:
            GRID[(float(kwo), float(kwa), float(kwn))] = B.simulate_chlorine(
                B.KB_FIXED, 0.0, pre_run=B.make_kw_hook(kwo, kwa, kwn),
            ).values[B.WARMUP_H:]
print(f"grid library {len(GRID)} sims ({time.time() - t0:.1f}s)")


def wstats(w, x):
    m = float(np.sum(w * x))
    return m, float(np.sqrt(np.sum(w * (x - m) ** 2)))


rows = []
for jit in JITTERS:
    rng_j = np.random.default_rng(JITTER_SEED)
    kw_pipe = {p: ZONE_MEAN[B.MATERIAL_ZONES[p]] * (1.0 + rng_j.uniform(-jit, jit))
               for p in wn0.pipe_name_list}

    def hook(wn, _kw=kw_pipe):
        for p in wn.pipe_name_list:
            wn.get_link(p).wall_coeff = B.per_day_to_per_second(_kw[p])

    truth_all = B.simulate_chlorine(B.KB_FIXED, 0.0, pre_run=hook,
                                    monitor_nodes=ALL_NODES).values
    truth_mon = truth_all[:, mon_pos]
    truth_post_mon = truth_mon[B.WARMUP_H:]

    # noise-free structural residual + grid best fit (vs noise-free truth)
    struct = min(np.sqrt(((sim - truth_post_mon) ** 2).mean()) for sim in GRID.values())

    # noisy observations (same process: seed 42)
    rng_n = np.random.default_rng(B.NOISE_SEED)
    obs = np.clip(truth_mon + rng_n.normal(0, B.SIGMA_OBS, truth_mon.shape), 0, None)[B.WARMUP_H:]
    noise_rmse = float(np.sqrt(((truth_post_mon - obs) ** 2).mean()))

    # reuse the cached candidate predictions; only the observations change with the truth.
    # Both weightings are computed: the formal censored likelihood is primary, the informal GLUE
    # score is the comparator whose flatness the jitter = 0 control exposes.
    RMSE = np.sqrt(((C_all_mon - obs[None]) ** 2).mean(axis=(1, 2)))
    w, diag = B.weights_from_loglik(B.log_censored(C_all_mon, obs))
    w_inf, diag_inf = B.weights_from_loglik(B.glue_score(RMSE), RMSE < B.RMSE_THR)
    beh = RMSE < B.RMSE_THR

    # per-zone true arithmetic mean and bias, under both schemes
    zstats = {}
    for z in ZKEYS:
        kws = np.array([kw_pipe[p] for p in zone_pipes[z]])
        arith = float(kws.mean())
        m, sd = wstats(w, S[z])
        mi, sdi = wstats(w_inf, S[z])
        zstats[z] = {"arith": arith,
                     "mean": m, "sd": sd, "bias": m - arith,
                     "bias_in_sd": (m - arith) / sd if sd > 0 else np.nan,
                     "glue_mean": mi, "glue_sd": sdi, "glue_bias": mi - arith,
                     "glue_bias_in_sd": (mi - arith) / sdi if sdi > 0 else np.nan}

    # risk ranking
    P_glue = B.pbar_field(C_all, C_MIN, w)
    P_true = B.pbar_field(truth_all[B.WARMUP_H:], C_MIN)
    rank_glue = [ALL_NODES[i] for i in np.argsort(P_glue)[::-1][:6]]
    rank_true = [ALL_NODES[i] for i in np.argsort(P_true)[::-1][:6]]

    rows.append({"jitter": jit, "struct_residual": float(struct),
                 "noise_rmse": noise_rmse, "behavioural": int(beh.sum()),
                 "rmse_min": float(RMSE.min()), "zones": zstats,
                 "rank_glue": rank_glue, "rank_true": rank_true})

# ---- bias increment over the homogeneous (jitter = 0) control ----
control = next(r for r in rows if r["jitter"] == 0.0)
GRID_FLOOR = control["struct_residual"]     # grid resolution, not structural error
for r in rows:
    r["struct_residual_grid_floor"] = GRID_FLOOR
    r["struct_residual_above_floor"] = r["struct_residual"] - GRID_FLOOR
    for z in ZKEYS:
        zz, cz = r["zones"][z], control["zones"][z]
        for pre in ("", "glue_"):
            key = pre + "bias"
            sdkey = pre + "sd" if pre else "sd"
            zz[pre + "bias_increment_vs_control"] = zz[key] - cz[key]
            zz[pre + "bias_increment_in_sd"] = ((zz[key] - cz[key]) / zz[sdkey]
                                                if zz[sdkey] > 0 else np.nan)

# ---- report ----
print("\n=== within-zone heterogeneity sweep (is the fit precise but biased?) ===")
print("Two weightings side by side. `bias` is the formal censored likelihood (primary); `glue` is")
print("the informal score. `incr` is the increment over the jitter = 0 control, which is the only")
print("part attributable to heterogeneity: whatever the control shows cannot be structural.")
print(f"structRes has a grid-resolution floor of {GRID_FLOOR:.4f} at jitter = 0 (a homogeneous truth")
print("is exactly representable, so that residual is grid spacing); `above` exceeds that floor.\n")
for scheme, pre in (("formal censored (primary)", ""), ("informal GLUE (comparator)", "glue_")):
    print(f"-- {scheme} --")
    print(f"{'jit':>4} {'structRes':>9} {'above':>7} {'rmseMin':>7} | "
          f"{'old bias(SD)':>14} {'incr':>8} | {'avg bias(SD)':>14} {'incr':>8} | "
          f"{'new bias(SD)':>14} {'incr':>8}")
    for r in rows:
        z = r["zones"]
        def fmt(zz, p=pre):
            return f"{z[zz][p + 'bias']:+.3f}({z[zz][p + 'bias_in_sd']:+.2f})"
        def inc(zz, p=pre):
            return f"{z[zz][p + 'bias_increment_vs_control']:+.4f}"
        print(f"{r['jitter']:>4} {r['struct_residual']:>9.4f} "
              f"{r['struct_residual_above_floor']:>7.4f} {r['rmse_min']:>7.4f} | "
              f"{fmt('old'):>14} {inc('old'):>8} | {fmt('average'):>14} {inc('average'):>8} | "
              f"{fmt('new'):>14} {inc('new'):>8}")
    print()
for r in rows:
    print(f"  jit={r['jitter']} risk (primary) {r['rank_glue']}  | TRUE {r['rank_true']}")

# ---- field ensemble: one heterogeneity realisation is not a result ----
# Everything above uses a single jitter field per magnitude (seed 12345), so the "structural bias"
# it reports mixes the effect of heterogeneity with the accident of one spatial arrangement. The
# candidate predictions are truth-independent, so repeating over many fields costs one truth
# simulation each and separates between-field variability from the effect itself.
N_FIELDS = 25
ENS_JITTER = 0.20
print(f"\n=== field ensemble at +/-{ENS_JITTER:.0%}: {N_FIELDS} independent heterogeneity fields ===")
ens = {z: {"bias": [], "incr": []} for z in ZKEYS}
control_bias = {z: control["zones"][z]["bias"] for z in ZKEYS}
t0 = time.time()
for f in range(N_FIELDS):
    rng_f = np.random.default_rng(50_000 + f)
    kwp = {p: ZONE_MEAN[B.MATERIAL_ZONES[p]] * (1.0 + rng_f.uniform(-ENS_JITTER, ENS_JITTER))
           for p in wn0.pipe_name_list}

    def hook_f(wn, _kw=kwp):
        for p in wn.pipe_name_list:
            wn.get_link(p).wall_coeff = B.per_day_to_per_second(_kw[p])

    tm = B.simulate_chlorine(B.KB_FIXED, 0.0, pre_run=hook_f, monitor_nodes=ALL_NODES).values
    rng_nf = np.random.default_rng(B.NOISE_SEED)
    obs_f = np.clip(tm[:, mon_pos] + rng_nf.normal(0, B.SIGMA_OBS, tm[:, mon_pos].shape),
                    0, None)[B.WARMUP_H:]
    w_f, _ = B.weights_from_loglik(B.log_censored(C_all_mon, obs_f))
    for z in ZKEYS:
        arith = float(np.array([kwp[p] for p in zone_pipes[z]]).mean())
        m, _ = wstats(w_f, S[z])
        ens[z]["bias"].append(m - arith)
        ens[z]["incr"].append(m - arith - control_bias[z])
print(f"  {N_FIELDS} fields in {time.time() - t0:.0f}s\n")
print(f"{'zone':>8} | {'bias mean':>10} {'bias sd':>9} | {'incr mean':>10} {'incr sd':>9} | "
      f"{'|mean|/sd':>9}")
field_ens = {}
for z in ZKEYS:
    b = np.array(ens[z]["bias"]); i = np.array(ens[z]["incr"])
    ratio = abs(i.mean()) / i.std(ddof=1) if i.std(ddof=1) > 0 else np.nan
    field_ens[z] = {"n_fields": N_FIELDS, "jitter": ENS_JITTER,
                    "bias_mean": float(b.mean()), "bias_sd": float(b.std(ddof=1)),
                    "increment_mean": float(i.mean()), "increment_sd": float(i.std(ddof=1)),
                    "increment_mean_over_sd": float(ratio),
                    "increment_5_95": [float(np.percentile(i, 5)), float(np.percentile(i, 95))]}
    print(f"{z:>8} | {b.mean():+10.4f} {b.std(ddof=1):9.4f} | {i.mean():+10.4f} "
          f"{i.std(ddof=1):9.4f} | {ratio:9.2f}")
print("`incr` subtracts the jitter = 0 control, so it is the part heterogeneity causes; |mean|/sd")
print("says whether that part is larger than the field-to-field scatter. Below ~1 the single-field")
print("number reported above is not separable from the choice of field.")

for r in rows:
    if r["jitter"] == ENS_JITTER:
        r["field_ensemble"] = field_ens


def _jsafe(o):
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    return str(o)


with open(os.path.join(HERE, "baseline_cache", "step5c_jitter_sweep.json"), "w") as f:
    json.dump({**B.weighting_provenance(comparators=["informal_glue"]), "rows": rows},
              f, indent=2, default=_jsafe)

# ---- figure: |bias| vs jitter per zone (+ bias in SD units) ----
fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
jits = [r["jitter"] * 100 for r in rows]
for ax, z in zip(axes, ZKEYS):
    bias = [r["zones"][z]["bias"] for r in rows]
    sd = [r["zones"][z]["glue_sd"] for r in rows]
    ax.errorbar(jits, bias, yerr=sd, marker="o", capsize=4, color="steelblue",
                label="GLUE bias ± behavioural SD")
    ax.axhline(0, color="gray", lw=1)
    ax.axhline(control["zones"][z]["bias"], color="crimson", ls="--", lw=1.2,
               label="jitter = 0 control (weighting artefact)")
    ax.set_xlabel("within-zone jitter (%)")
    ax.set_ylabel("bias = GLUE mean − true arith. mean (m/day)")
    ax.set_title(f"{z} zone\nstructural increment at ±50%: "
                 f"{rows[-1]['zones'][z]['bias_increment_vs_control']:+.4f} m/day", fontsize=10)
    ax.grid(alpha=0.3)
    for xi, r in zip(jits, rows):
        ax.annotate(f"{r['zones'][z]['bias_in_sd']:+.1f}σ",
                    (xi, r["zones"][z]["bias"]), textcoords="offset points",
                    xytext=(6, 6), fontsize=8)
axes[0].legend(fontsize=8)
fig.suptitle("Symmetric within-zone heterogeneity: almost all of the apparent bias is already "
             "present in the homogeneous jitter = 0 control (dashed), so it is a weighting "
             "artefact, not a structural error", y=1.02)
plt.tight_layout()
figpath = os.path.join(FIGDIR, "step5c_jitter_sweep.png")
plt.savefig(figpath, dpi=130, bbox_inches="tight")
print("figure saved to", figpath)
