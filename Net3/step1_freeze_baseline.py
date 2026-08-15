"""Step 1: freeze and cache the three-zone baseline.

Builds the synthetic truth, one noisy observation set, and the forward prediction of every prior
draw, then caches all of it so the later experiments never re-run EPANET.

Three things changed relative to the draft version of this script, each for a stated reason:

  warm-up      120 h, not 24 h. Set by the Step 0 convergence test; 120 + 48 = 168 h is also the
               model horizon. The 48 h assessment window is unchanged, so N_RESID stays 294.
  sampling     scrambled Sobol, 2^13 = 8192 draws, not 2000 pseudo-random ones. The formal
               likelihood is far sharper than the informal GLUE score, so 2000 prior draws gave it
               an effective sample size of only ~37. Every leading 2^k subset of a Sobol set is
               itself balanced, which is what makes the convergence check below exact.
  weighting    three schemes are computed and cached, not one. The formal censored Gaussian
               likelihood is the primary analysis; the formal iid one isolates the effect of
               treating a clipped zero as an exact measurement; the informal GLUE score is retained
               as a comparator because the draft used it and because the contrast is a result.

Outputs (in Net3/baseline_cache/):
  baseline.npz       - draws, RMSE, both formal log-likelihoods, C_all (N_MC x Tn x nodes),
                       truth_all, obs_glue, noisy, node lists, monitor positions
  baseline_meta.json - full config, per-scheme summaries, and the sampling-convergence table
"""
import os
import json
import time
import numpy as np

import wq_common as B

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "baseline_cache")
os.makedirs(OUT, exist_ok=True)

ALL_NODES = B.all_junctions()
mon_pos = [ALL_NODES.index(m) for m in B.MONITOR_NODES]

# 1) synthetic truth over the whole network (monitors are a subset)
truth_all_df = B.simulate_chlorine(
    B.KB_FIXED, 0.0,
    pre_run=B.make_kw_hook(B.KW_OLD_TRUE, B.KW_AVG_TRUE, B.KW_NEW_TRUE),
    monitor_nodes=ALL_NODES,
)
truth_all = truth_all_df.values                       # (DURATION_H + 1, 92)
truth_index = np.asarray(truth_all_df.index, dtype=float)
truth_mon = truth_all_df[B.MONITOR_NODES].values

# 2) baseline noisy observations: Gaussian(0, sigma) added to monitor truth, clipped at the sensor
#    floor of 0. The clipping is not cosmetic: those points are left-censored, which is what the
#    censored likelihood below accounts for and the iid one does not.
rng_noise = np.random.default_rng(B.NOISE_SEED)
raw = truth_mon + rng_noise.normal(0.0, B.SIGMA_OBS, size=truth_mon.shape)
n_clip = int((raw < 0).sum())
noisy = np.clip(raw, 0.0, None)
obs_glue = noisy[B.WARMUP_H:]                          # (Tn, 6)
n_clip_window = int((obs_glue == 0).sum())
# TWO noise floors, because they are not interchangeable. The candidate RMSE below is computed on
# the 49x6 assessment window, so that is the floor it must be compared with; the full-record figure
# is kept because it describes the observation set as a whole. Quoting the full-record floor beside a
# window RMSE compares two different samples.
noise_rmse_full = float(np.sqrt(((noisy - truth_mon) ** 2).mean()))
noise_rmse = float(np.sqrt(((noisy - truth_mon)[B.WARMUP_H:] ** 2).mean()))   # calibration window

# 3) prior draws: scrambled Sobol over the three-zone box
draws = B.prior_draws()
S_old, S_avg, S_new = draws["old"], draws["avg"], draws["new"]

print(f"=== Step 1: freezing the baseline ===")
print(f"  timing      {B.DURATION_H} h run, {B.WARMUP_H} h warm-up -> {obs_glue.shape[0]} reporting "
      f"points x {len(B.MONITOR_NODES)} monitors = {B.N_RESID} residuals")
print(f"  sampling    scrambled Sobol, 2^{B.N_MC_LOG2} = {B.N_MC} draws (seed {B.SAMPLE_SEED})")
print(f"  noise       sigma {B.SIGMA_OBS} mg/L, seed {B.NOISE_SEED}; {n_clip} of "
      f"{raw.size} raw points below zero, {n_clip_window} of {obs_glue.size} inside the window")
print(f"  noise floor RMSE {noise_rmse:.4f} mg/L on the calibration window "
      f"({noise_rmse_full:.4f} over the full record)\n")

# 4) forward each draw over the whole network; cache predictions, RMSE and both log-likelihoods
Tn = obs_glue.shape[0]
C_all = np.empty((B.N_MC, Tn, len(ALL_NODES)), dtype=np.float32)
RMSE = np.empty(B.N_MC)
t0 = time.time()
for s in range(B.N_MC):
    cs = B.simulate_chlorine(
        B.KB_FIXED, 0.0,
        pre_run=B.make_kw_hook(S_old[s], S_avg[s], S_new[s]),
        monitor_nodes=ALL_NODES,
    ).values[B.WARMUP_H:]
    C_all[s] = cs.astype(np.float32)
    d = cs[:, mon_pos] - obs_glue
    RMSE[s] = np.sqrt((d ** 2).mean())
    if (s + 1) % 1000 == 0:
        print(f"  {s + 1}/{B.N_MC} forward runs ({time.time() - t0:.0f}s)")
run_s = time.time() - t0
print(f"  forward runs done in {run_s:.0f}s\n")

C_mon = C_all[:, :, mon_pos].astype(np.float64)
loglik_iid = B.log_gaussian(C_mon, obs_glue)
loglik_cens = B.log_censored(C_mon, obs_glue)

# 5) the three weighting schemes, side by side
SCHEMES = {
    "formal_censored": (loglik_cens, None),
    "formal_iid": (loglik_iid, None),
    "informal_glue": (B.glue_score(RMSE), RMSE < B.RMSE_THR),
    "informal_glue_draft_thr": (B.glue_score(RMSE), RMSE < B.RMSE_THR_DRAFT),
}
TRUE = {"old": B.KW_OLD_TRUE, "avg": B.KW_AVG_TRUE, "new": B.KW_NEW_TRUE}
S = {"old": S_old, "avg": S_avg, "new": S_new}
PRIOR_SD = {g: (B.PRIOR[g][1] - B.PRIOR[g][0]) / np.sqrt(12) for g in S}

schemes_out = {}
print("weighting schemes (kw in m/day; SD retained = behavioural SD / prior SD):")
for name, (ll, mask) in SCHEMES.items():
    w, diag = B.weights_from_loglik(ll, mask)
    entry = {"diagnostics": diag, "coef": {}}
    print(f"  {name:>23} ESS {diag['ess']:8.1f} ({diag['ess_frac'] * 100:5.2f}% of {B.N_MC})  "
          f"max w {diag['max_weight']:.4f}  entropy {diag['entropy_bits']:5.2f} of "
          f"{diag['entropy_bits_if_uniform']:.2f} bits")
    for g in ("old", "avg", "new"):
        m, sd = B.weighted_mean_sd(w, S[g])
        q = B.weighted_quantile(S[g], w, [0.05, 0.5, 0.95])
        entry["coef"][g] = {"mean": m, "sd": sd, "q05_50_95": [float(v) for v in q],
                            "true": TRUE[g], "sd_retained": sd / PRIOR_SD[g],
                            "bias": m - TRUE[g]}
        print(f"      {g:>4} {m:+.4f} +/- {sd:.4f}  [{q[0]:+.4f}, {q[2]:+.4f}]  "
              f"true {TRUE[g]:+.3f}  SD retained {100 * sd / PRIOR_SD[g]:5.1f}%")
    schemes_out[name] = entry

# 6) sampling convergence: leading 2^m Sobol subsets are themselves balanced designs
print(f"\nsampling convergence under the primary scheme ({B.PRIMARY_WEIGHTING}):")
print(f"{'N':>6} {'ESS':>8} {'ESS%':>6} | " +
      " | ".join(f"{g} median [5,95]" for g in ("old", "avg", "new")))
conv = []
prev = None
for m in range(10, B.N_MC_LOG2 + 1):
    n = 2 ** m
    w, diag = B.weights_from_loglik(loglik_cens[:n])
    row = {"n": n, "ess": diag["ess"], "ess_frac": diag["ess_frac"], "coef": {}}
    line = f"{n:>6} {diag['ess']:8.1f} {diag['ess_frac'] * 100:5.2f}% |"
    for g in ("old", "avg", "new"):
        q = B.weighted_quantile(S[g][:n], w, [0.05, 0.5, 0.95])
        row["coef"][g] = {"q05": float(q[0]), "median": float(q[1]), "q95": float(q[2])}
        line += f" {q[1]:+.4f} [{q[0]:+.4f},{q[2]:+.4f}] |"
    if prev is not None:
        drifts = [abs(row["coef"][g]["median"] - prev["coef"][g]["median"]) / PRIOR_SD[g]
                  for g in ("old", "avg", "new")]
        row["median_drift_in_prior_sd"] = {g: d for g, d in zip(("old", "avg", "new"), drifts)}
        line += f"  max median drift {max(drifts):.3f} prior SD"
    print(line)
    conv.append(row)
    prev = row

summary = {
    "n_mc": B.N_MC, "n_resid": B.N_RESID, "Tn": int(Tn), "n_nodes": len(ALL_NODES),
    "rmse_min": float(RMSE.min()),
    "noise_rmse": noise_rmse,                       # 49x6 calibration window; compare with rmse_min
    "noise_rmse_full_record": noise_rmse_full,      # 169x6, the observation set as a whole
    "n_clip_raw": n_clip, "n_clip_window": n_clip_window,
    "behavioural_count_primary_thr": int((RMSE < B.RMSE_THR).sum()),
    "behavioural_count_draft_thr": int((RMSE < B.RMSE_THR_DRAFT).sum()),
    "true": [B.KW_OLD_TRUE, B.KW_AVG_TRUE, B.KW_NEW_TRUE],
    "prior_sd": PRIOR_SD,
    "schemes": schemes_out,
    "sampling_convergence": conv,
    "runtime_s": round(run_s, 1),
}

config = {
    "monitor_nodes": B.MONITOR_NODES,
    "inlet_mgl": B.INLET_CHLORINE_MGL, "tank_mgl": B.TANK_INIT_MGL,
    "duration_h": B.DURATION_H, "warmup_h": B.WARMUP_H,
    "warmup_basis": "step0_warmup_convergence.json",
    "hydraulic_timestep_s": B.HYDRAULIC_TIMESTEP_S,
    "report_timestep_s": B.REPORT_TIMESTEP_S,
    "quality_timestep_s": B.QUALITY_TIMESTEP_S,
    "kb_fixed": B.KB_FIXED,
    "kw_true": [B.KW_OLD_TRUE, B.KW_AVG_TRUE, B.KW_NEW_TRUE],
    "prior": B.PRIOR, "n_mc": B.N_MC, "sampling": f"scrambled Sobol 2^{B.N_MC_LOG2}",
    "sigma_obs": B.SIGMA_OBS,
    "rmse_thr_primary": B.RMSE_THR, "rmse_thr_draft": B.RMSE_THR_DRAFT,
    "noise_seed": B.NOISE_SEED, "sample_seed": B.SAMPLE_SEED,
    "primary_weighting": B.PRIMARY_WEIGHTING,
    "sigma_convention": "sigma = 0.1 mg/L is one standard deviation of the Gaussian observation "
                        "error",
    "weighting_convention": "the formal schemes carry no behavioural threshold; a hard cut-off is "
                            "part of the informal GLUE comparator, not of a likelihood",
}

np.savez_compressed(
    os.path.join(OUT, "baseline.npz"),
    S_old=S_old, S_avg=S_avg, S_new=S_new,
    RMSE=RMSE, loglik_iid=loglik_iid, loglik_censored=loglik_cens,
    C_all=C_all, truth_all=truth_all, truth_index=truth_index,
    obs_glue=obs_glue, noisy=noisy,
    all_nodes=np.asarray(ALL_NODES), monitor_nodes=np.asarray(B.MONITOR_NODES),
    mon_pos=np.asarray(mon_pos),
)
with open(os.path.join(OUT, "baseline_meta.json"), "w") as f:
    json.dump({**B.weighting_provenance(), "config": config, "summary": summary}, f, indent=2)

print(f"\nbehavioural at {B.RMSE_THR}: {summary['behavioural_count_primary_thr']}/{B.N_MC}; "
      f"at {B.RMSE_THR_DRAFT}: {summary['behavioural_count_draft_thr']}/{B.N_MC}; "
      f"min RMSE {RMSE.min():.4f}")
print("cache saved to", OUT)
