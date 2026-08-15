"""Step 7c: error autocorrelation (AR(1)) — Priority-2 #4, done properly.

Instead of mechanically multiplying every CRLB by the scalar ESS factor √[(1+ρ)/(1−ρ)], we build
the full AR(1) observation covariance and recompute the Fisher information:

  per-sensor 49×49 block   Σ_block[t,s] = σ² · ρ^|t−s|      (hourly AR(1), ρ = 0.4)
  full covariance          Σ = block_diag(6 sensor blocks)  (sensors independent)
  Fisher (AR(1))           F = Jᵀ Σ⁻¹ J     vs     independent  F = JᵀJ / σ²
  CRLB_j = sqrt[(F⁻¹)_jj]

Each coefficient's CRLB widens by its own factor (not a single 1.53), because J varies across
nodes and time. This is the idealised Case A (k_w only, k_b known, no bias) plus autocorrelation.
"""
import os
import json
import numpy as np
from scipy.linalg import block_diag
import wq_common as B

HERE = os.path.dirname(os.path.abspath(__file__))
ZKEYS = ["old", "average", "new"]
TRUTH = [B.KW_OLD_TRUE, B.KW_AVG_TRUE, B.KW_NEW_TRUE]
PRIOR_SD = {"old": (1.5 - 0.2) / np.sqrt(12), "average": (0.2 - 0.04) / np.sqrt(12),
            "new": (0.10 - 0.005) / np.sqrt(12)}
HS = [B.FD_STEP[z] for z in ZKEYS]          # scale-dependent, same as Step 7
SIGMA = 0.10
RHO = 0.40
NT = B.DURATION_H - B.WARMUP_H + 1          # 49 post-warm-up hours
NMON = len(B.MONITOR_NODES)


def sim_mon(kwo, kwa, kwn):
    return B.simulate_chlorine(B.KB_FIXED, 0.0,
                               pre_run=B.make_kw_hook(kwo, kwa, kwn)).values[B.WARMUP_H:]  # (49,6)


# ---- Jacobian in NODE-MAJOR order (so Σ is block-diagonal per sensor) ----
cols = []
for j in range(3):
    h = HS[j]
    kp, km = list(TRUTH), list(TRUTH)
    kp[j] += h
    km[j] -= h
    dC = (sim_mon(*kp) - sim_mon(*km)) / (2 * h)        # (49, 6)
    cols.append(dC.T.ravel())                            # node-major: node0 h0..48, node1 ...
J = np.column_stack(cols)                                # (294, 3)

# ---- AR(1) covariance ----
tt = np.arange(NT)
block = SIGMA ** 2 * RHO ** np.abs(tt[:, None] - tt[None, :])   # 49×49
Sigma = block_diag(*[block for _ in range(NMON)])               # 294×294
Sigma_inv = np.linalg.inv(Sigma)

# ---- Fisher + CRLB, independent vs AR(1) ----
F_indep = (J.T @ J) / SIGMA ** 2
F_ar1 = J.T @ Sigma_inv @ J
crlb_indep = np.sqrt(np.diag(np.linalg.inv(F_indep)))
crlb_ar1 = np.sqrt(np.diag(np.linalg.inv(F_ar1)))

ess_factor = np.sqrt((1 + RHO) / (1 - RHO))              # scalar reference 1.53
n_eff = NT * NMON * (1 - RHO) / (1 + RHO)

print(f"=== Step 7c: AR(1) autocorrelation (ρ={RHO}) — Case A (k_w only) ===")
print(f"scalar ESS reference: inflation √[(1+ρ)/(1−ρ)] = {ess_factor:.2f}; "
      f"N_eff ≈ {n_eff:.0f} of {NT*NMON}\n")
print(f"{'coef':>8} | {'CRLB indep':>10} | {'CRLB AR(1)':>10} | {'widening':>8} | "
      f"{'CRLB/prior indep':>16} | {'CRLB/prior AR(1)':>16}")
report = {"rho": RHO, "ess_factor": float(ess_factor), "n_eff": float(n_eff),
          "fd_default_steps": dict(B.FD_STEP), "coef": {}}
for j, z in enumerate(ZKEYS):
    wf = crlb_ar1[j] / crlb_indep[j]
    report["coef"][z] = {"crlb_indep": float(crlb_indep[j]), "crlb_ar1": float(crlb_ar1[j]),
                         "widening": float(wf),
                         "crlb_over_prior_indep": float(crlb_indep[j] / PRIOR_SD[z]),
                         "crlb_over_prior_ar1": float(crlb_ar1[j] / PRIOR_SD[z])}
    print(f"{z:>8} | {crlb_indep[j]:10.3f} | {crlb_ar1[j]:10.3f} | {wf:8.2f}x | "
          f"{crlb_indep[j]/PRIOR_SD[z]:16.2f} | {crlb_ar1[j]/PRIOR_SD[z]:16.2f}")

# ---- sweep over rho: rho is ASSUMED, not estimated, so its influence must be shown ----
# The baseline observations are generated iid, so no value of rho is supported by this data set.
# What can be reported honestly is how sensitive the intervals are to the assumption.
RHO_SWEEP = [0.0, 0.2, 0.4, 0.6, 0.8]
sweep = []
print(f"\n=== sensitivity to the ASSUMED autocorrelation (rho is not estimated here) ===")
print(f"{'rho':>5} {'N_eff':>7} | " + " | ".join(f"{z} CRLB (x indep)" for z in ZKEYS))
for rho in RHO_SWEEP:
    blk = SIGMA ** 2 * rho ** np.abs(tt[:, None] - tt[None, :])
    Sig = block_diag(*[blk for _ in range(NMON)])
    crlb = np.sqrt(np.diag(np.linalg.inv(J.T @ np.linalg.inv(Sig) @ J)))
    row = {"rho": rho, "n_eff": float(NT * NMON * (1 - rho) / (1 + rho)),
           "ess_factor": float(np.sqrt((1 + rho) / (1 - rho))) if rho < 1 else None,
           "coef": {z: {"crlb": float(crlb[j]),
                        "widening_vs_indep": float(crlb[j] / crlb_indep[j]),
                        "crlb_over_prior": float(crlb[j] / PRIOR_SD[z])}
                    for j, z in enumerate(ZKEYS)}}
    sweep.append(row)
    print(f"{rho:>5.1f} {row['n_eff']:>7.0f} | " +
          " | ".join(f"{crlb[j]:.4f} ({crlb[j] / crlb_indep[j]:.2f}x)" for j in range(3)))

worst = sweep[-1]["coef"]
print(f"\nAcross the tested range the CRLB inflates by up to "
      f"{max(v['widening_vs_indep'] for v in worst.values()):.1f}x, and at rho = 0.8 the ratio to "
      f"the prior SD reaches {max(v['crlb_over_prior'] for v in worst.values()):.2f}.")
print("So every interval in this project is a FLOOR whose height depends on an assumption that the")
print("synthetic data cannot test. Estimating rho needs real monitoring residuals (ACF/PACF or a")
print("likelihood in which rho is a free parameter); until then the thesis must write 'if the errors")
print("had AR(1) structure with rho = X' and never 'the errors are autocorrelated with rho = X'.")

report["rho_sweep"] = sweep
report["rho_is_assumed"] = ("the baseline observations are generated iid, so this is an "
                            "assumed-covariance sensitivity analysis, not an estimate of rho")
with open(os.path.join(HERE, "baseline_cache", "step7c_ar1.json"), "w") as f:
    json.dump(report, f, indent=2)
print("\nsaved step7c_ar1.json")
