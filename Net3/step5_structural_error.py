"""Step 5: structural-error experiment (breaks the inverse crime).

Truth: pipe-level heterogeneous wall decay. Each pipe keeps its zone (old/avg/new) but its
true coefficient is  k_w,p = zone_mean * (1 + delta_p),  delta_p ~ U(-0.2, 0.2), fixed seed.
Calibration model: the SAME three-zone HOMOGENEOUS model (one coefficient per zone) — so it is
structurally simpler than the truth and cannot represent within-zone heterogeneity.

Reports:
  - noise-free structural misfit: best homogeneous fit vs the heterogeneous truth (how far the
    minimum RMSE sits above the noise floor -> irreducible structural residual);
  - grid-search best fit and GLUE behavioural means against the true field's arithmetic and
    length-weighted per-zone averages (is the fit precise but biased?);
  - risk-map ranking vs the true heterogeneous field.
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
os.makedirs(FIGDIR, exist_ok=True)

ZKEYS = ["old", "average", "new"]                       # match MATERIAL_ZONES labels
ZONE_MEAN = {"old": B.KW_OLD_TRUE, "average": B.KW_AVG_TRUE, "new": B.KW_NEW_TRUE}
PRIOR_OF = {"old": B.PRIOR["old"], "average": B.PRIOR["avg"], "new": B.PRIOR["new"]}
JITTER = 0.20
JITTER_SEED = 12345
C_MIN = 0.2

wn0 = wntr.network.WaterNetworkModel(B.NET3_INP)
ALL_NODES = wn0.junction_name_list
mon_pos = [ALL_NODES.index(m) for m in B.MONITOR_NODES]
PIPE_LEN = {p: wn0.get_link(p).length for p in wn0.pipe_name_list}

# ---- 1) build the pipe-level heterogeneous truth ----
rng_j = np.random.default_rng(JITTER_SEED)
KW_PIPE = {}
for p in wn0.pipe_name_list:
    z = B.MATERIAL_ZONES[p]
    KW_PIPE[p] = ZONE_MEAN[z] * (1.0 + rng_j.uniform(-JITTER, JITTER))


def hetero_hook(wn):
    for p in wn.pipe_name_list:
        wn.get_link(p).wall_coeff = B.per_day_to_per_second(KW_PIPE[p])


# per-zone reference "true effective" values
zone_pipes = {z: [p for p in wn0.pipe_name_list if B.MATERIAL_ZONES[p] == z]
              for z in ZKEYS}
true_ref = {}
for z, ps in zone_pipes.items():
    kws = np.array([KW_PIPE[p] for p in ps])
    lens = np.array([PIPE_LEN[p] for p in ps])
    lw = float((kws * lens).sum() / lens.sum()) if lens.sum() > 0 else float(kws.mean())
    true_ref[z] = {"n_pipes": len(ps), "arith_mean": float(kws.mean()),
                   "length_weighted": lw, "min": float(kws.min()), "max": float(kws.max()),
                   "zone_mean": ZONE_MEAN[z]}

# heterogeneous truth at all nodes + monitors
truth_all = B.simulate_chlorine(B.KB_FIXED, 0.0, pre_run=hetero_hook,
                                monitor_nodes=ALL_NODES).values
truth_mon = truth_all[:, mon_pos]

# ---- 2) homogeneous 3-zone grid library (monitors), built once ----
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

truth_post = truth_mon[B.WARMUP_H:]                       # (49, 6) noise-free
NOISE_FLOOR_STRUCT = None

# noise-free structural misfit (best homogeneous fit vs heterogeneous truth, no noise)
best_nf, best_nf_rmse = None, np.inf
for kw, sim in GRID.items():
    r = np.sqrt(((sim - truth_post) ** 2).mean())
    if r < best_nf_rmse:
        best_nf_rmse, best_nf = r, kw

# ---- 3) noisy observations (same process as baseline: seed 42, sigma, clip) ----
rng_n = np.random.default_rng(B.NOISE_SEED)
noisy = np.clip(truth_mon + rng_n.normal(0, B.SIGMA_OBS, truth_mon.shape), 0, None)
obs = noisy[B.WARMUP_H:]
# TWO noise floors, and only one of them is the right comparator here. Every RMSE below is computed
# on the 49x6 assessment window (`obs`), so the floor it must be read against is the realised noise
# RMSE on that same window. The full-record figure describes the observation set as a whole and is
# kept for that purpose only; quoting it beside a window RMSE compares two different samples, which
# is exactly the mistake this line used to make.
noise_rmse_full = float(np.sqrt(((noisy - truth_mon) ** 2).mean()))
noise_rmse = float(np.sqrt(((noisy - truth_mon)[B.WARMUP_H:] ** 2).mean()))

# grid-search best fit vs noisy obs
best_fit, best_rmse = None, np.inf
for kw, sim in GRID.items():
    r = np.sqrt(((sim - obs) ** 2).mean())
    if r < best_rmse:
        best_rmse, best_fit = r, kw

# ---- 4) GLUE (2000) homogeneous over all nodes, vs the heterogeneous-truth noisy obs ----
rng = np.random.default_rng(B.SAMPLE_SEED)
S_old = rng.uniform(*B.PRIOR["old"], B.N_MC)
S_avg = rng.uniform(*B.PRIOR["avg"], B.N_MC)
S_new = rng.uniform(*B.PRIOR["new"], B.N_MC)
Tn = obs.shape[0]
C_all = np.empty((B.N_MC, Tn, len(ALL_NODES)), dtype=np.float32)
RMSE = np.empty(B.N_MC)
t0 = time.time()
for s in range(B.N_MC):
    cs = B.simulate_chlorine(B.KB_FIXED, 0.0,
                             pre_run=B.make_kw_hook(S_old[s], S_avg[s], S_new[s]),
                             monitor_nodes=ALL_NODES).values[B.WARMUP_H:]
    C_all[s] = cs.astype(np.float32)
    RMSE[s] = np.sqrt(((cs[:, mon_pos] - obs) ** 2).mean())
    if (s + 1) % 500 == 0:
        print(f"  GLUE {s + 1}/{B.N_MC} ({time.time() - t0:.1f}s)")

L = np.exp(-0.5 * (RMSE / B.SIGMA_OBS) ** 2)
beh = RMSE < B.RMSE_THR
w = L * beh
w = w / w.sum()


def wmean(x):
    return float(np.sum(w * x))


def wstd(x):
    m = wmean(x)
    return float(np.sqrt(np.sum(w * (x - m) ** 2)))


glue = {"old": (wmean(S_old), wstd(S_old)), "average": (wmean(S_avg), wstd(S_avg)),
        "new": (wmean(S_new), wstd(S_new))}
GRIDFIT = {"old": best_fit[0], "average": best_fit[1], "new": best_fit[2]}

# ---- 5) risk maps: fitted GLUE ensemble vs true heterogeneous field ----
P_node_glue = B.pbar_field(C_all, C_MIN, w)
truth_post_all = truth_all[B.WARMUP_H:]
P_node_true = B.pbar_field(truth_post_all, C_MIN)
rank_glue = [ALL_NODES[i] for i in np.argsort(P_node_glue)[::-1][:6]]
rank_true = [ALL_NODES[i] for i in np.argsort(P_node_true)[::-1][:6]]

# ---- report ----
report = {
    "weighting": "informal_glue",
    "weighting_role": "comparator; SUPERSEDED by step5c (jitter ensemble, formal primary) and "
                      "step5d (structured heterogeneity, formal primary). A single jitter field and "
                      "the informal score cannot separate a structural bias from either the field "
                      "drawn or the score's own offset.",
    "primary_weighting_elsewhere": B.PRIMARY_WEIGHTING,
    "jitter": JITTER, "jitter_seed": JITTER_SEED,
    "noise_floor_rmse": noise_rmse,                # 49x6 window; compare with rmse_min
    "noise_floor_rmse_full_record": noise_rmse_full,   # 169x6, the observation set as a whole
    "structural_min_rmse_noisefree": best_nf_rmse,
    "structural_best_homogeneous": best_nf,
    "gridsearch_best_fit_noisy": best_fit, "gridsearch_best_rmse": best_rmse,
    "glue_behavioural": {z: {"mean": glue[z][0], "sd": glue[z][1]} for z in glue},
    "true_reference": true_ref,
    "behavioural_count": int(beh.sum()), "rmse_min": float(RMSE.min()),
    "risk_rank_glue": rank_glue, "risk_rank_true": rank_true,
    "risk_rank_baseline_reported": ["131", "243", "141", "139", "15", "143"],
}
def _jsafe(o):
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


with open(os.path.join(HERE, "baseline_cache", "step5_structural_error.json"), "w") as f:
    json.dump(report, f, indent=2, default=_jsafe)

print("\n=== structural error (+/-20% pipe-level jitter) ===")
print(f"noise floor RMSE            = {noise_rmse:.4f} on the calibration window "
      f"({noise_rmse_full:.4f} over the full record)")
print(f"structural min RMSE (noise-free best homogeneous) = {best_nf_rmse:.4f}"
      f"  -> irreducible structural residual = {best_nf_rmse:.4f} mg/L")
print(f"grid-search best fit (noisy) = {best_fit}, RMSE {best_rmse:.4f}")
print(f"GLUE behavioural {int(beh.sum())}/{B.N_MC}, min RMSE {RMSE.min():.4f}")
print(f"{'zone':>4} | {'arith':>7} {'lenwt':>7} {'range':>16} | {'GLUE mean±sd':>16} "
      f"| {'gridfit':>8} | bias(GLUE-arith)")
for z in ZKEYS:
    tr = true_ref[z]
    gf = GRIDFIT[z]
    bias = glue[z][0] - tr["arith_mean"]
    print(f"{z:>4} | {tr['arith_mean']:7.3f} {tr['length_weighted']:7.3f} "
          f"[{tr['min']:.2f},{tr['max']:.2f}] | {glue[z][0]:7.3f}±{glue[z][1]:.3f} "
          f"| {gf:8.3f} | {bias:+.3f}")
print(f"\nrisk top-6 GLUE  : {rank_glue}")
print(f"risk top-6 TRUE  : {rank_true}")
print(f"risk top-6 base  : {report['risk_rank_baseline_reported']}")

# ---- figure: true pipe-kw spread vs fitted effective values, per zone ----
fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
for ax, z in zip(axes, ZKEYS):
    kws = np.array([KW_PIPE[p] for p in zone_pipes[z]])
    ax.hist(kws, bins=12, color="0.8", edgecolor="0.5", label="true pipe k_w (spread)")
    ax.axvline(true_ref[z]["arith_mean"], color="black", lw=2, label="true arith. mean")
    ax.axvline(glue[z][0], color="steelblue", lw=2, ls="--", label="GLUE behavioural mean")
    ax.axvline(GRIDFIT[z], color="crimson", lw=1.5, ls=":", label="grid best fit")
    ax.set_title(f"{z}  (n={true_ref[z]['n_pipes']} pipes)")
    ax.set_xlabel("k_w (m/day)")
    ax.grid(alpha=0.3)
axes[0].set_ylabel("pipe count")
axes[0].legend(fontsize=7)
fig.suptitle(f"Structural error (±{int(JITTER*100)}% pipe-level jitter): fitted homogeneous "
             "coefficients vs the true heterogeneous field", y=1.02)
plt.tight_layout()
figpath = os.path.join(FIGDIR, "step5_structural_error.png")
plt.savefig(figpath, dpi=130, bbox_inches="tight")
print("figure saved to", figpath)
