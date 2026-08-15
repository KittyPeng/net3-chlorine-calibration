"""Step 7c (profile part): AR(1) profile likelihood — the REAL computation.

Step 7b's profile used the independent Gaussian NLL = (N/2σ²)·RMSE², which only needs the scalar
RMSE. Under AR(1) the objective is NLL = ½·eᵀ Σ⁻¹ e, which needs the full 294-vector residual e at
every grid point. The Step 7b grid cache stored only RMSE, so this script rebuilds the 21³ grid
storing the residual vectors (node-major, to match the block-diagonal AR(1) covariance of Step 7c)
and then reads off BOTH the independent profile (sanity check vs Step 7b) and the AR(1) profile.

  per-sensor 49×49 block   Σ_block[t,s] = σ²·ρ^|t−s|   (ρ = 0.4)
  Σ = block_diag(6 blocks) (294×294, sensors independent)
  NLL_indep(θ) = (1/2σ²)·eᵀe          profile = min over the other two coefficients
  NLL_ar1(θ)   = ½·eᵀ Σ⁻¹ e           ΔNLL ≤ 1.92 → 95% single-parameter interval
"""
import os
import json
import time
import numpy as np
from scipy.linalg import block_diag
import wq_common as B
import provenance

HERE = os.path.dirname(os.path.abspath(__file__))
ZKEYS = ["old", "average", "new"]
TRUE = [B.KW_OLD_TRUE, B.KW_AVG_TRUE, B.KW_NEW_TRUE]
NG = 21
SIGMA = B.SIGMA_OBS
RHO = 0.40
DNLL_95 = 1.92
NT, NMON = B.DURATION_H - B.WARMUP_H + 1, len(B.MONITOR_NODES)   # 49, 6

cache = np.load(os.path.join(HERE, "baseline_cache", "baseline.npz"), allow_pickle=True)
obs = cache["obs_glue"]                                          # (49, 6)
grids = {"old": np.linspace(*B.PRIOR["old"], NG),
         "average": np.linspace(*B.PRIOR["avg"], NG),
         "new": np.linspace(*B.PRIOR["new"], NG)}
go, ga, gn = grids["old"], grids["average"], grids["new"]

# ---- residual grid e[i,j,k,:] = (sim − obs), node-major (node0 h0..48, node1 ...) ----
RES_CACHE = os.path.join(HERE, "baseline_cache", "step7c_resid_grid.npy")
E = provenance.load_keyed_array(RES_CACHE, ng=NG)
if E is not None:
    print("loaded cached residual grid", E.shape, "(configuration matches)")
else:
    E = np.empty((NG, NG, NG, NT * NMON), dtype=np.float64)
    t0 = time.time()
    for i, kwo in enumerate(go):
        for j, kwa in enumerate(ga):
            for k, kwn in enumerate(gn):
                sim = B.simulate_chlorine(B.KB_FIXED, 0.0,
                                          pre_run=B.make_kw_hook(kwo, kwa, kwn)).values[B.WARMUP_H:]
                E[i, j, k] = (sim - obs).T.ravel()               # node-major
        if (i + 1) % 3 == 0:
            print(f"  grid slice {i+1}/{NG} ({time.time()-t0:.0f}s)", flush=True)
    provenance.save_keyed_array(RES_CACHE, E, ng=NG)
    print(f"built residual grid in {time.time()-t0:.0f}s")

# ---- AR(1) covariance and its inverse ----
tt = np.arange(NT)
block = SIGMA ** 2 * RHO ** np.abs(tt[:, None] - tt[None, :])
Sigma = block_diag(*[block for _ in range(NMON)])
Sinv = np.linalg.inv(Sigma)

# ---- objectives over the whole grid ----
ee = np.einsum("ijkr,ijkr->ijk", E, E)                           # eᵀe
NLL_indep = ee / (2 * SIGMA ** 2)                                # (N/2σ²)RMSE²
NLL_ar1 = 0.5 * np.einsum("ijkr,rs,ijks->ijk", E, Sinv, E)       # ½ eᵀΣ⁻¹e


def interval_95(x, nll_grid, axes):
    prof = nll_grid.min(axis=axes)
    dnll = prof - prof.min()
    below = np.where(dnll <= DNLL_95)[0]
    if len(below) == 0:
        return None, prof
    return {"lo": float(x[below.min()]), "hi": float(x[below.max()]),
            "open_low": bool(below.min() == 0), "open_high": bool(below.max() == len(x) - 1),
            "half_width": float((x[below.max()] - x[below.min()]) / 2)}, prof


AX = {"old": (1, 2), "average": (0, 2), "new": (0, 1)}
XG = {"old": go, "average": ga, "new": gn}
report = {"NG": NG, "sigma": SIGMA, "rho": RHO, "coef": {}}
print(f"\n=== Step 7c: AR(1) profile likelihood (ρ={RHO}, ΔNLL≤1.92=95%) ===")
print(f"{'coef':>8} | {'independent 95%':>24} | {'AR(1) 95%':>24} | {'hw indep':>9} | {'hw AR1':>7} | {'widen':>6}")
for z in ZKEYS:
    iv_i, _ = interval_95(XG[z], NLL_indep, AX[z])
    iv_a, _ = interval_95(XG[z], NLL_ar1, AX[z])
    wf = (iv_a["half_width"] / iv_i["half_width"]) if (iv_i and iv_a) else float("nan")
    report["coef"][z] = {"indep": iv_i, "ar1": iv_a, "widening": wf}
    si = f"[{iv_i['lo']:.3f},{iv_i['hi']:.3f}]" if iv_i else "flat"
    sa = f"[{iv_a['lo']:.3f},{iv_a['hi']:.3f}]" if iv_a else "flat"
    print(f"{z:>8} | {si:>24} | {sa:>24} | {iv_i['half_width']:9.3f} | "
          f"{iv_a['half_width']:7.3f} | {wf:5.2f}x")

with open(os.path.join(HERE, "baseline_cache", "step7c_profile_ar1.json"), "w") as f:
    json.dump(report, f, indent=2)
print("\nsaved step7c_profile_ar1.json")
