"""Step 7: a-priori identifiability via Fisher information / CRLB, with nuisance parameters.

Sensitivity (Jacobian) of monitored chlorine to each grouped wall coefficient at the truth, then
the Fisher information and the Cramér–Rao lower bound (CRLB = smallest achievable standard error).
The CRLB is computed for the three k_w under three nested models, adding the confounders the
reviewer flagged:

  Case A  k_w only              (k_b known, no sensor bias)          — idealised best case
  Case B  k_w + k_b             (k_b unknown; trades off with k_w)   — Priority-2 #5
  Case C  k_w + k_b + 6 offsets (per-monitor systematic bias)        — Priority-2 #1

Marginal CRLB for k_w uses the Schur complement  F_pp − F_pn F_nn⁻¹ F_np.
Independent Gaussian noise is assumed (AR(1) autocorrelation would inflate further, Priority-2 #4).

The idealised Case A shows the observations DO contain information about all three coefficients
(so GLUE's "avg/new uninformed" is partly the inefficiency of its informal likelihood — the
Stedinger 2008 / Mantovan & Todini 2006 critique). Adding the realistic confounders (Case C)
inflates the avg/new CRLB above their prior width, reconciling Fisher with GLUE.
"""
import os
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import wq_common as B

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures")
ZKEYS = ["old", "average", "new"]
TRUTH = [B.KW_OLD_TRUE, B.KW_AVG_TRUE, B.KW_NEW_TRUE]
KB = B.KB_FIXED
PRIOR_SD = {"old": (1.5 - 0.2) / np.sqrt(12), "average": (0.2 - 0.04) / np.sqrt(12),
            "new": (0.10 - 0.005) / np.sqrt(12)}
HS = [B.FD_STEP[z] for z in ZKEYS]        # one step per coefficient, ~2-5% of its own truth
HKB = B.FD_STEP_KB
SIGMA = 0.10
NMON = len(B.MONITOR_NODES)


def sim(kb, kwo, kwa, kwn):
    return B.simulate_chlorine(kb, 0.0, pre_run=B.make_kw_hook(kwo, kwa, kwn)).values[B.WARMUP_H:]


# ---- Jacobian columns: old, avg, new, kb (finite differences), then 6 sensor offsets ----
def col(idx, h):
    kwp, kwm, kbp, kbm = list(TRUTH), list(TRUTH), KB, KB
    if idx < 3:
        kwp[idx] += h; kwm[idx] -= h
    else:
        kbp += h; kbm -= h
    return ((sim(kbp, *kwp) - sim(kbm, *kwm)) / (2 * h)).ravel()   # (294,) row = hour*6 + node


phys = [col(0, HS[0]), col(1, HS[1]), col(2, HS[2]), col(3, HKB)]  # old,avg,new,kb
n_obs = phys[0].size
node_of_row = np.arange(n_obs) % NMON                              # row -> monitor index
offsets = [(node_of_row == m).astype(float) for m in range(NMON)] # per-monitor additive bias
J_full = np.column_stack(phys + offsets)                          # (294, 4+6)
JtJ = J_full.T @ J_full
LABELS = ["old", "average", "new", "kb"] + [f"off_{n}" for n in B.MONITOR_NODES]


def marginal(idx_all, sigma):
    """Marginal CRLB (σ units) and eigenvalues for the 3 k_w within the column subset idx_all."""
    sub = JtJ[np.ix_(idx_all, idx_all)]
    p = [0, 1, 2]                                                 # kw are always the first three
    n = [i for i in range(len(idx_all)) if i not in p]
    Fpp = sub[np.ix_(p, p)]
    if n:
        Fpn = sub[np.ix_(p, n)]
        Fnn = sub[np.ix_(n, n)]
        marg = Fpp - Fpn @ np.linalg.pinv(Fnn) @ Fpn.T
    else:
        marg = Fpp
    cov = np.linalg.pinv(marg)
    crlb = sigma * np.sqrt(np.clip(np.diag(cov), 0, None))
    eig = np.linalg.eigvalsh(marg) / sigma ** 2
    return crlb, eig


cases = {"A: kw only": [0, 1, 2],
         "B: kw+kb": [0, 1, 2, 3],
         "C: kw+kb+6 offsets": list(range(J_full.shape[1]))}

# ---- finite-difference step convergence ----
# The default step is SCALE-DEPENDENT (wq_common.FD_STEP): one shared absolute step would be 2% of
# the old truth but 40% of the new one, and at 40% a central difference is a secant over a large arc
# rather than a derivative. Each coefficient is re-differenced over a range of steps and the
# resulting Case-A CRLB compared; a converged derivative gives a flat column. Every sweep is
# asserted to contain its own default, so the reported spread really does bound the error the
# default carries.
STEP_SWEEP = {0: [0.005, 0.01, 0.02, 0.04],
              1: [0.002, 0.005, 0.01, 0.02],
              2: [0.001, 0.0025, 0.005, 0.01]}
for i, z in enumerate(ZKEYS):
    assert B.FD_STEP[z] in STEP_SWEEP[i], f"default step for {z} is outside its convergence sweep"
print("=== finite-difference step convergence (Case A CRLB, other columns at their default step) ===")
print("default steps: " + ", ".join(f"{z} {B.FD_STEP[z]:g}" for z in ZKEYS) +
      f" (each ~{100 * B.FD_STEP['old'] / abs(TRUTH[0]):.0f}-"
      f"{100 * B.FD_STEP['new'] / abs(TRUTH[2]):.0f}% of its own truth)")
print(f"{'coef':>8} {'step':>8} {'step/|truth|':>12} {'CRLB':>9} {'rel. change':>12}")
fd_conv = {}
for idx, steps in STEP_SWEEP.items():
    z = ZKEYS[idx]
    rows_fd, prev = [], None
    for h in steps:
        cols3 = [col(i, h if i == idx else HS[i]) for i in range(3)]
        J3 = np.column_stack(cols3)
        crlb3 = SIGMA * np.sqrt(np.diag(np.linalg.pinv(J3.T @ J3)))
        rel = None if prev is None else abs(crlb3[idx] - prev) / prev
        rows_fd.append({"step": h, "step_over_truth": h / abs(TRUTH[idx]),
                        "is_default": h == B.FD_STEP[z],
                        "crlb": float(crlb3[idx]),
                        "rel_change_vs_previous": None if rel is None else float(rel)})
        print(f"{z:>8} {h:>8.4f} {h / abs(TRUTH[idx]):>11.0%} {crlb3[idx]:>9.5f} "
              f"{'—' if rel is None else f'{rel:11.2%}'}"
              f"{'   <- default' if h == B.FD_STEP[z] else ''}")
        prev = crlb3[idx]
    spread = max(r["crlb"] for r in rows_fd) / min(r["crlb"] for r in rows_fd) - 1.0
    fd_conv[z] = {"default_step": B.FD_STEP[z],
                  "default_step_over_truth": B.FD_STEP[z] / abs(TRUTH[idx]),
                  "rows": rows_fd, "crlb_spread_over_sweep": float(spread)}
    print(f"{'':>8} {'spread over the sweep':>32}: {spread:.2%}\n")
print("A flat column means the derivative is converged and the step is not driving the CRLB. Each")
print("default is a member of its own sweep, so the spread above bounds the error it carries.")

# Empirical posterior SD from the frozen ensemble, under EVERY weighting scheme.
# The CRLB is a bound on the variance of an efficient estimator under the correct likelihood, so
# only the FORMAL schemes are comparable with it. Comparing the CRLB with the informal GLUE SD --
# as the earlier version of this script did -- puts an efficient bound next to an inefficient
# estimator and makes the data look less informative than it is.
cache = np.load(os.path.join(HERE, "baseline_cache", "baseline.npz"), allow_pickle=True)
RMSE = cache["RMSE"]
S = {"old": cache["S_old"], "average": cache["S_avg"], "new": cache["S_new"]}
SCHEMES = {
    "formal_censored": (cache["loglik_censored"], None),
    "formal_iid": (cache["loglik_iid"], None),
    "informal_glue": (B.glue_score(RMSE), RMSE < B.RMSE_THR),
}
ensemble_sd, ensemble_diag = {}, {}
for sname, (ll, mask) in SCHEMES.items():
    wS, diag = B.weights_from_loglik(ll, mask)
    ensemble_sd[sname] = {z: B.weighted_mean_sd(wS, S[z])[1] for z in ZKEYS}
    ensemble_diag[sname] = diag

# ---- prior-scaled Fisher: a condition number on the raw matrix is a units artefact ----
# The three coefficients differ in scale by a factor of 20, so eigenvalues of the unnormalised
# information matrix mix m/day with m/day and say as much about the units as about the geometry.
# Rescaling by the prior SD makes the parameters dimensionless and comparable, which is the form in
# which a condition number and a "sloppy direction" mean something.
S_scale = np.diag([PRIOR_SD[z] for z in ZKEYS])
F_raw = JtJ[np.ix_([0, 1, 2], [0, 1, 2])] / SIGMA ** 2
F_scaled = S_scale @ F_raw @ S_scale
ev_raw = np.linalg.eigvalsh(F_raw)[::-1]
ev_sc, evec_sc = np.linalg.eigh(F_scaled)
order = np.argsort(ev_sc)[::-1]
ev_sc, evec_sc = ev_sc[order], evec_sc[:, order]
sloppy = evec_sc[:, -1]
print(f"\n=== Fisher geometry: raw vs prior-scaled (Case A) ===")
print(f"raw eigenvalues        {np.array2string(ev_raw, precision=1)}  cond "
      f"{ev_raw.max() / ev_raw.min():.0f}")
print(f"prior-scaled           {np.array2string(ev_sc, precision=1)}  cond "
      f"{ev_sc.max() / ev_sc.min():.1f}")
print(f"sloppiest direction (prior-SD units): " +
      ", ".join(f"{z} {v:+.3f}" for z, v in zip(ZKEYS, sloppy)))
print("The raw condition number is inflated by the parameter scales; on the dimensionless matrix the")
print("spread is far smaller, so the three coefficients are much more comparably informed than the")
print("raw number suggests. The sloppiest direction names the combination the data constrain least.")

report = {**B.weighting_provenance(),
          "sigma": SIGMA, "prior_sd": PRIOR_SD,
          "fd_default_steps": {**B.FD_STEP, "kb": HKB},
          "fd_step_convergence": fd_conv,
          "fisher_geometry": {
              "note": "a condition number on the unnormalised matrix mixes units; the prior-scaled "
                      "form is the dimensionless one",
              "eigenvalues_raw": [float(v) for v in ev_raw],
              "condition_number_raw": float(ev_raw.max() / ev_raw.min()),
              "eigenvalues_prior_scaled": [float(v) for v in ev_sc],
              "condition_number_prior_scaled": float(ev_sc.max() / ev_sc.min()),
              "sloppiest_direction_prior_sd_units": {z: float(v) for z, v in zip(ZKEYS, sloppy)}},
          "ensemble_sd_by_scheme": ensemble_sd, "ensemble_diagnostics": ensemble_diag,
          "comparability_note": "the CRLB is a bound for an efficient estimator under the correct "
                                "likelihood, so only the formal schemes are like-for-like; the "
                                "informal GLUE SD is listed as a comparator, not as a benchmark",
          "crlb_agreement_claim": "the formal posterior spread is LOCALLY CONSISTENT with the "
                                  "Case-A CRLB. It is not a demonstration of frequentist "
                                  "efficiency: the CRLB bounds the variance of an unbiased "
                                  "estimator over repeated sampling, whereas the posterior SD here "
                                  "is the width from ONE noise realisation, and both use the same "
                                  "model, truth and error assumptions. Step 14 does the repeated-"
                                  "sampling test that efficiency and calibration actually require.",
          "cases": {}}
print("=== Step 7: Fisher / CRLB with nuisance parameters (σ = 0.10) ===\n")
print(f"{'case':>22} | {'coef':>7} | {'CRLB':>7} | {'priorSD':>7} | {'CRLB/prior':>10} | identif.? | cond#")
for cname, idx in cases.items():
    crlb, eig = marginal(idx, SIGMA)
    cond = float(eig.max() / eig.min())
    report["cases"][cname] = {"condition_number": cond, "eigenvalues": [float(e) for e in eig[::-1]],
                              "coef": {}}
    for j, z in enumerate(ZKEYS):
        ratio = crlb[j] / PRIOR_SD[z]
        report["cases"][cname]["coef"][z] = {"crlb": float(crlb[j]), "crlb_over_prior": float(ratio),
                                             "identifiable": bool(ratio < 1.0)}
        print(f"{cname:>22} | {z:>7} | {crlb[j]:7.3f} | {PRIOR_SD[z]:7.3f} | {ratio:10.2f} | "
              f"{str(ratio < 1.0):>9} | {cond:8.1f}")
    print()

print("empirical ensemble SD by weighting scheme (compare with CRLB above):")
for sname, sd in ensemble_sd.items():
    tag = "like-for-like with the CRLB" if sname.startswith("formal") else "comparator only"
    print(f"  {sname:>16} ESS {ensemble_diag[sname]['ess']:7.1f} | "
          + ", ".join(f"{z} {sd[z]:.3f}" for z in ZKEYS) + f"   ({tag})")
print("        prior SD                  | "
      + ", ".join(f"{z} {PRIOR_SD[z]:.3f}" for z in ZKEYS))
print("\nCase A is the like-for-like benchmark for the FORMAL schemes (both fix k_b and assume no")
print("sensor bias). The informal GLUE SD is 2-3x wider than the formal one on the same data, which")
print("is the inefficiency of the informal score, not a property of the observations.")
print("\nAgreement between the formal posterior SD and the Case-A CRLB is a LOCAL CONSISTENCY check,")
print("not proof of efficiency: the bound is over repeated sampling while the posterior width comes")
print("from one noise realisation, and the two share a model, a truth and an error assumption. The")
print("repeated-sampling version of the claim — bias, estimator variance, interval coverage — is")
print("Step 14.")

# CRLB(C) vs sigma (scales linearly)
crlbC_unit, _ = marginal(cases["C: kw+kb+6 offsets"], 1.0)
report["crlbC_vs_sigma"] = {f"{s:.2f}": {z: float(s * crlbC_unit[j]) for j, z in enumerate(ZKEYS)}
                            for s in [0.05, 0.10, 0.15]}

# ---- k_b confounding: two quantities that are easy to confuse and carry opposite signs ----
# The Jacobian-column cosine F_ij / sqrt(F_ii F_jj) says how collinear two SENSITIVITY DIRECTIONS
# are: a property of the experimental design. The estimator correlation, normalised from
# Cov(theta_hat) ~ sigma^2 (J'J)^-1, says how the two ESTIMATES co-vary. For a pair of parameters the
# two generally carry OPPOSITE SIGNS, and it is the negative estimator correlation that carries the
# bulk-wall compensation story Step 8b measures empirically. step7_verify.py prints the same two
# columns side by side; they are recorded here so the log can quote them against an artifact.
idxB = cases["B: kw+kb"]
FB = JtJ[np.ix_(idxB, idxB)]
dB = np.sqrt(np.diag(FB))
jac_cos = FB / np.outer(dB, dB)
covB = SIGMA ** 2 * np.linalg.pinv(FB)
sB = np.sqrt(np.diag(covB))
par_corr = covB / np.outer(sB, sB)
crlbA_kb, _ = marginal(cases["A: kw only"], SIGMA)
crlbB_kb, _ = marginal(idxB, SIGMA)
report["kb_confounding"] = {
    "note": "jacobian_cosine_with_kb is a design collinearity; estimator_correlation_with_kb is how "
            "the estimates co-vary. They differ in sign and must not be interchanged: only the "
            "second is a parameter correlation.",
    "by_coef": {z: {"jacobian_cosine_with_kb": float(jac_cos[j, 3]),
                    "estimator_correlation_with_kb": float(par_corr[j, 3]),
                    "crlb_inflation_when_kb_freed": float(crlbB_kb[j] / crlbA_kb[j])}
                for j, z in enumerate(ZKEYS)}}
print("\nk_b confounding (design collinearity vs estimator correlation — opposite signs):")
print(f"{'coef':>8} | {'cos(J_kw,J_kb)':>14} | {'corr(kw_hat,kb_hat)':>19} | CRLB inflation when k_b freed")
for j, z in enumerate(ZKEYS):
    print(f"{z:>8} | {jac_cos[j, 3]:+14.3f} | {par_corr[j, 3]:+19.3f} | "
          f"{crlbB_kb[j] / crlbA_kb[j]:.2f}x")


def _jsafe(o):
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.integer):
        return int(o)
    return str(o)


with open(os.path.join(HERE, "baseline_cache", "step7_fisher.json"), "w") as f:
    json.dump(report, f, indent=2, default=_jsafe)

# ---- figure: CRLB/prior for cases A/B/C per coefficient (log; >1 = unidentifiable) ----
fig, ax = plt.subplots(figsize=(9, 5))
x = np.arange(len(ZKEYS))
wid = 0.25
for k, cname in enumerate(cases):
    vals = [report["cases"][cname]["coef"][z]["crlb_over_prior"] for z in ZKEYS]
    ax.bar(x + (k - 1) * wid, vals, wid, label=cname)
ax.axhline(1.0, color="crimson", ls="--", lw=1.5, label="CRLB = prior (identifiable below)")
ax.set_yscale("log")
ax.set_xticks(x)
ax.set_xticklabels(ZKEYS)
ax.set_ylabel("CRLB / prior SD (log; <1 identifiable, >1 not)")
ax.set_title("Step 7 — Fisher/CRLB: adding k_b and sensor offsets pushes avg/new above the prior\n"
             "(idealised A: all identifiable → realistic C: only old survives)")
ax.legend(fontsize=8)
ax.grid(alpha=0.3, axis="y")
plt.tight_layout()
figpath = os.path.join(FIGDIR, "step7_fisher.png")
plt.savefig(figpath, dpi=130, bbox_inches="tight")
print("\nfigure saved to", figpath)
