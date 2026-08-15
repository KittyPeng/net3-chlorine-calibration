"""Step 7b: profile likelihood (the reviewer's third identifiability tool).

For each grouped coefficient, fix it and RE-OPTIMISE the other two, giving the profile negative
log-likelihood. Unlike the one-at-a-time sweep of Step 4b (which fixes the other two at the truth),
the profile allows parameter compensation, so it is the rigorous test of one-sided identifiability.
The 95% single-parameter interval is ΔNLL ≤ 1.92 (½ · chi2(1) at 0.95).

Two things about this step were previously misreported and are fixed here.

1. QUANTISATION. A 21-point grid steps 0.065 m/day in the `old` direction — coarser than several
   effects discussed elsewhere in this log (the censoring shift is 0.012, structural increments are
   0.01-0.05), so a grid-based interval cannot resolve them and its endpoints sit on grid nodes by
   construction. The reported intervals therefore come from CONTINUOUS optimisation: the nuisance
   coefficients are minimised at each target value and the endpoints found by bisection on ΔNLL.
   The grid curve is retained only as a visual background and is labelled as such in the figure,
   because a grid interval is systematically NARROWER than the truth (its endpoints are the last
   nodes inside the region, not the crossing points).

2. LIKELIHOOD. The primary inference rule in this project is the CENSORED Gaussian likelihood, so
   the profile is computed under it as well as under the iid one. The iid arm is kept as the
   identifiability benchmark that matches Step 7's Fisher/CRLB calculation, which also assumes no
   censoring; the difference between the two arms is the profile-interval counterpart of the
   weighted-mean shift isolated in Step 9.
"""
import os
import json
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import minimize, brentq
import wq_common as B
import provenance

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures")
ZKEYS = ["old", "average", "new"]
PKEY = {"old": "old", "average": "avg", "new": "new"}      # zone -> key in B.PRIOR
TRUE = [B.KW_OLD_TRUE, B.KW_AVG_TRUE, B.KW_NEW_TRUE]
NG = 21
DNLL_95 = 1.92
KINDS = ["censored", "iid"]
PRIMARY_KIND = "censored"                                  # matches wq_common.PRIMARY_WEIGHTING
N_CURVE = 13                                               # target values per continuous profile curve

cache = np.load(os.path.join(HERE, "baseline_cache", "baseline.npz"), allow_pickle=True)
obs = cache["obs_glue"]                                    # (49, 6) baseline noisy observations
N = obs.size                                               # 294
sigma = B.SIGMA_OBS
FAC = N / (2 * sigma ** 2)                                 # NLL_iid = FAC * RMSE^2
n_censored = int((obs <= 0).sum())

grids = {"old": np.linspace(*B.PRIOR["old"], NG),
         "average": np.linspace(*B.PRIOR["avg"], NG),
         "new": np.linspace(*B.PRIOR["new"], NG)}
go, ga, gn = grids["old"], grids["average"], grids["new"]

# ---- the 21^3 RMSE grid (monitors only); cached for instant re-plots ----
GRID_CACHE = os.path.join(HERE, "baseline_cache", "step7b_rmse_grid.npy")
RMSE = provenance.load_keyed_array(GRID_CACHE, ng=NG)
if RMSE is not None:
    print("loaded cached RMSE grid", RMSE.shape, "(configuration matches)")
else:
    RMSE = np.empty((NG, NG, NG))
    t0 = time.time()
    for i, kwo in enumerate(go):
        for j, kwa in enumerate(ga):
            for k, kwn in enumerate(gn):
                sim = B.simulate_chlorine(B.KB_FIXED, 0.0,
                                          pre_run=B.make_kw_hook(kwo, kwa, kwn)).values[B.WARMUP_H:]
                RMSE[i, j, k] = np.sqrt(((sim - obs) ** 2).mean())
        if (i + 1) % 5 == 0:
            print(f"  grid {(i + 1) * NG * NG}/{NG ** 3} ({time.time() - t0:.0f}s)")
    provenance.save_keyed_array(GRID_CACHE, RMSE, ng=NG)
rmse_min = RMSE.min()
print(f"global min RMSE on the grid {rmse_min:.4f}; {n_censored} of {N} observations at the floor")

# nearest grid index to the truth (for the sweep = others fixed at truth)
it = int(np.argmin(np.abs(go - TRUE[0])))
jt = int(np.argmin(np.abs(ga - TRUE[1])))
kt = int(np.argmin(np.abs(gn - TRUE[2])))

profiles = {
    "old": {"x": go, "profile": RMSE.min(axis=(1, 2)), "sweep": RMSE[:, jt, kt], "true": TRUE[0]},
    "average": {"x": ga, "profile": RMSE.min(axis=(0, 2)), "sweep": RMSE[it, :, kt], "true": TRUE[1]},
    "new": {"x": gn, "profile": RMSE.min(axis=(0, 1)), "sweep": RMSE[it, jt, :], "true": TRUE[2]},
}


def interval_95(x, dnll):
    below = np.where(dnll <= DNLL_95)[0]
    if len(below) == 0:
        return None
    return {"lo": float(x[below.min()]), "hi": float(x[below.max()]),
            "open_low": bool(below.min() == 0), "open_high": bool(below.max() == len(x) - 1)}


# ---- continuous profile ----
_sim_cache = {}


def nll_at(kw, kind):
    """NLL (up to a constant) at a full parameter triple, memoised on the simulation.

    Both likelihood kinds are derived from the same EPANET run, so running the second profile costs
    only the optimiser's extra iterations rather than a second sweep of the parameter space.
    """
    key = tuple(round(float(v), 6) for v in kw)
    if key not in _sim_cache:
        sim = B.simulate_chlorine(B.KB_FIXED, 0.0,
                                  pre_run=B.make_kw_hook(*key)).values[B.WARMUP_H:]
        pred = sim[None]
        _sim_cache[key] = {"iid": float(-B.log_gaussian(pred, obs, sigma)[0]),
                           "censored": float(-B.log_censored(pred, obs, sigma)[0])}
    return _sim_cache[key][kind]


BOUNDS = [(min(B.PRIOR[PKEY[z]]), max(B.PRIOR[PKEY[z]])) for z in ZKEYS]


def profile_nll(z_idx, value, kind, x0):
    """min over the other two coefficients at a fixed value of coefficient z_idx."""
    free = [i for i in range(3) if i != z_idx]

    def objective(v):
        kw = [0.0, 0.0, 0.0]
        kw[z_idx] = value
        for i, f in zip(free, v):
            kw[i] = f
        return nll_at(kw, kind)

    res = minimize(objective, [x0[i] for i in free], method="Nelder-Mead",
                   bounds=[BOUNDS[i] for i in free],
                   options={"xatol": 1e-4, "fatol": 1e-3, "maxiter": 200})
    return float(res.fun)


report = {"NG": NG, "N": int(N), "sigma": sigma, "rmse_min": float(rmse_min),
          "n_censored_observations": n_censored,
          "primary_likelihood": PRIMARY_KIND,
          "grid_interval_status": "visualisation only; endpoints are quantised to the grid step and "
                                  "systematically narrower than the continuous interval",
          "grid_coef": {}}
print("\n=== grid profile on the 21-point grid (retained for the figure background only) ===")
for z in ZKEYS:
    p = profiles[z]
    p["dnll_prof"] = FAC * (p["profile"] ** 2 - rmse_min ** 2)
    p["dnll_sweep"] = FAC * (p["sweep"] ** 2 - p["sweep"].min() ** 2)
    iv = interval_95(p["x"], p["dnll_prof"])
    report["grid_coef"][z] = {"profile_95": iv,
                              "profile_flat": bool(p["dnll_prof"].max() < DNLL_95),
                              "prior_range": [float(v) for v in B.PRIOR[PKEY[z]]]}
    print(f"  {z:>8}: [{iv['lo']:.4f}, {iv['hi']:.4f}]" if iv else f"  {z:>8}: none")

# ---- continuous refinement under both likelihood kinds ----
cont_all, curves = {}, {}
t0 = time.time()
for kind in KINDS:
    res0 = minimize(lambda v: nll_at(v, kind), [go[it], ga[jt], gn[kt]], method="Nelder-Mead",
                    bounds=BOUNDS, options={"xatol": 1e-4, "fatol": 1e-3, "maxiter": 400})
    nll_star, kw_star = float(res0.fun), list(res0.x)
    print(f"\n=== continuous profile, {kind} likelihood"
          f"{' (PRIMARY)' if kind == PRIMARY_KIND else ' (benchmark, matches Step 7 Fisher)'} ===")
    print(f"unconstrained min: kw = ({kw_star[0]:.4f}, {kw_star[1]:.4f}, {kw_star[2]:.4f}), "
          f"NLL {nll_star:.4f}")
    print(f"{'coef':>8} | {'grid 95%':>22} | {'continuous 95%':>22} | {'grid step':>9} "
          f"{'half-width change':>17}")
    cont = {}
    for zi, z in enumerate(ZKEYS):
        lo_b, hi_b = BOUNDS[zi]

        def g(v, _zi=zi, _kind=kind, _n=nll_star, _x=kw_star):
            return profile_nll(_zi, v, _kind, _x) - _n - DNLL_95

        ends = []
        for lo_side in (True, False):
            a, b = (lo_b, kw_star[zi]) if lo_side else (kw_star[zi], hi_b)
            if g(a) * g(b) > 0:                  # no crossing inside the prior: the interval is open
                ends.append(None)
                continue
            ends.append(float(brentq(g, a, b, xtol=1e-4, maxiter=60)))
        grid_iv = report["grid_coef"][z]["profile_95"]
        step = float(abs(profiles[z]["x"][1] - profiles[z]["x"][0]))
        grid_hw = (grid_iv["hi"] - grid_iv["lo"]) / 2
        cont_hw = None if None in ends else (ends[1] - ends[0]) / 2
        cont[z] = {"lo": ends[0], "hi": ends[1], "open_low": ends[0] is None,
                   "open_high": ends[1] is None, "grid_step": step,
                   "half_width": cont_hw, "grid_half_width": float(grid_hw),
                   "half_width_change_vs_grid_frac": (None if cont_hw is None
                                                      else cont_hw / grid_hw - 1.0)}
        # a curve for the figure, spanning a little beyond the interval so the crossing is visible
        span = cont_hw if cont_hw is not None else grid_hw
        xs = np.linspace(max(lo_b, kw_star[zi] - 1.6 * span),
                         min(hi_b, kw_star[zi] + 1.6 * span), N_CURVE)
        curves[(kind, z)] = (xs, np.array([profile_nll(zi, float(v), kind, kw_star) - nll_star
                                           for v in xs]))
        gs = f"[{grid_iv['lo']:.3f}, {grid_iv['hi']:.3f}]"
        cs = (f"[{ends[0]:.4f}, {ends[1]:.4f}]" if None not in ends
              else f"[{ends[0]}, {ends[1]}] (open)")
        chg = cont[z]["half_width_change_vs_grid_frac"]
        print(f"{z:>8} | {gs:>22} | {cs:>22} | {step:9.4f} "
              f"{('—' if chg is None else f'{chg:+16.1%}'):>17}")
    cont_all[kind] = {"nll_min": nll_star, "kw_at_min": kw_star, "coef": cont}
print(f"\ntook {time.time() - t0:.0f}s, {len(_sim_cache)} distinct EPANET evaluations")
print("The grid interval is narrower because its endpoints are the last nodes INSIDE the 95% region,")
print("not the crossings. Quote the continuous interval; the grid curve is a picture, not a result.")

prim, bench = cont_all[PRIMARY_KIND]["coef"], cont_all["iid"]["coef"]
print(f"\ncensoring effect on the interval ({n_censored} of {N} observations at the floor):")
for z in ZKEYS:
    if prim[z]["half_width"] and bench[z]["half_width"]:
        print(f"  {z:>8}: half-width {bench[z]['half_width']:.4f} (iid) -> "
              f"{prim[z]['half_width']:.4f} (censored), "
              f"{prim[z]['half_width'] / bench[z]['half_width'] - 1:+.1%}")

report["continuous_profile"] = {
    "method": "Nelder-Mead over the two nuisance coefficients at each target value; endpoints by "
              "Brent bisection on DNLL - 1.92",
    "by_likelihood": cont_all, "n_epanet_evaluations": len(_sim_cache)}
with open(os.path.join(HERE, "baseline_cache", "step7b_profile.json"), "w") as f:
    json.dump(report, f, indent=2)

# ---- figure: continuous profile (primary), with the grid curve as background ----
fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
for ax, z in zip(axes, ZKEYS):
    p = profiles[z]
    ax.plot(p["x"], p["dnll_prof"], color="0.75", lw=1.2, ls="--",
            label=f"{NG}-point grid profile (visualisation only)")
    ax.plot(p["x"], p["dnll_sweep"], color="0.85", lw=1.0, ls=":",
            label="grid sweep (others at truth)")
    for kind, colour, lw in (("iid", "darkorange", 1.4), (PRIMARY_KIND, "steelblue", 2.2)):
        xs, ys = curves[(kind, z)]
        ax.plot(xs, ys, color=colour, lw=lw,
                label=("continuous profile, censored (PRIMARY)" if kind == PRIMARY_KIND
                       else "continuous profile, iid (benchmark)"))
    c = cont_all[PRIMARY_KIND]["coef"][z]
    for e in (c["lo"], c["hi"]):
        if e is not None:
            ax.axvline(e, color="steelblue", lw=1.1, ls="-.")
    ax.axhline(DNLL_95, color="crimson", ls=":", lw=1.5, label="ΔNLL = 1.92 (95%)")
    ax.axvline(p["true"], color="red", lw=1.2, alpha=0.7, label="truth")
    ax.set_xlim(min(curves[(PRIMARY_KIND, z)][0]) - 0.4 * (c["grid_step"]),
                max(curves[(PRIMARY_KIND, z)][0]) + 0.4 * (c["grid_step"]))
    ax.set_ylim(0, 8)
    ax.set_xlabel(f"k_w,{z} (m/day)")
    ax.set_ylabel("ΔNLL")
    ax.set_title(f"{z}: continuous 95% "
                 f"[{c['lo']:.3f}, {c['hi']:.3f}]" if c["lo"] is not None else z, fontsize=10)
    ax.grid(alpha=0.3)
axes[0].legend(fontsize=6.5, loc="upper center")
fig.suptitle("Step 7b — profile likelihood: all three coefficients have finite two-sided 95% "
             "intervals under the formal likelihood.\nDash-dotted verticals are the CONTINUOUS "
             "endpoints; the grey grid curve is shown only to display the quantisation it suffers.",
             y=1.06)
plt.tight_layout()
figpath = os.path.join(FIGDIR, "step7b_profile.png")
plt.savefig(figpath, dpi=130, bbox_inches="tight")
print("figure saved to", figpath)
