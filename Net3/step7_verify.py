"""Verify the Step 7 Jacobian and reconcile Fisher vs GLUE.

(a) Directly print how much monitored C changes when each coefficient is perturbed, to confirm the
    Jacobian magnitudes are physical (not a bug).
(b) Compare CRLB with kb FIXED vs kb MARGINALISED (kb trades off with kw): marginalising kb should
    inflate the kw CRLBs toward the pessimistic (GLUE-like) view.
(c) Report the two quantities that are easy to confuse, separately:
      - Jacobian-column cosine, F_ij / sqrt(F_ii F_jj): how COLLINEAR two sensitivity directions
        are. This is a property of the design, and it is NOT a parameter correlation;
      - estimator correlation, normalised from Cov(theta_hat) ~ sigma^2 (J'J)^-1: how the two
        ESTIMATES co-vary. For two parameters the two differ in SIGN, so they must not be
        interchanged.
"""
import numpy as np
import wq_common as B

TRUTH = [B.KW_OLD_TRUE, B.KW_AVG_TRUE, B.KW_NEW_TRUE]
KB = B.KB_FIXED
ZKEYS = ["old", "average", "new"]
HS = [B.FD_STEP[z] for z in ZKEYS]          # scale-dependent, same as Step 7
HKB = B.FD_STEP_KB
PRIOR_SD = {"old": (1.5 - 0.2) / np.sqrt(12), "average": (0.2 - 0.04) / np.sqrt(12),
            "new": (0.10 - 0.005) / np.sqrt(12)}


def sim(kb, kwo, kwa, kwn):
    return B.simulate_chlorine(kb, 0.0, pre_run=B.make_kw_hook(kwo, kwa, kwn)).values[B.WARMUP_H:]


# ---- (a) direct C changes ----
base = sim(KB, *TRUTH)
print("=== direct ΔC at monitors (perturb one coeff, others at truth) ===")
print(f"{'perturb':>14} | " + " ".join(f"{n:>8}" for n in B.MONITOR_NODES))
for j, name in enumerate(ZKEYS):
    h = HS[j]
    kp = list(TRUTH); kp[j] += h
    km = list(TRUTH); km[j] -= h
    dC = sim(KB, *kp) - sim(KB, *km)                       # ΔC over 2h
    per_node_max = np.max(np.abs(dC), axis=0)
    print(f"{name + ' ±' + f'{h:g}':>14} | " + " ".join(f"{v:8.4f}" for v in per_node_max)
          + f"   (max|ΔC| over 2h={2 * h:g})")

# ---- (b) Fisher with kb fixed vs kb marginalised ----
def jac_col(perturb_idx, h):
    # perturb_idx: 0=old,1=avg,2=new,3=kb
    kwp, kwm = list(TRUTH), list(TRUTH)
    kbp, kbm = KB, KB
    if perturb_idx < 3:
        kwp[perturb_idx] += h; kwm[perturb_idx] -= h
    else:
        kbp += h; kbm -= h
    return ((sim(kbp, *kwp) - sim(kbm, *kwm)) / (2 * h)).ravel()


cols = [jac_col(0, HS[0]), jac_col(1, HS[1]), jac_col(2, HS[2]), jac_col(3, HKB)]
J4 = np.column_stack(cols)                                 # (294, 4): old,avg,new,kb
sigma = 0.10
names = ["old", "average", "new", "kb"]

# kb FIXED: use only the 3 kw columns
J3 = J4[:, :3]
crlb_fixed = sigma * np.sqrt(np.diag(np.linalg.inv(J3.T @ J3)))

# kb MARGINALISED: 4x4 inverse, take kw diagonals
crlb_marg = sigma * np.sqrt(np.diag(np.linalg.inv(J4.T @ J4)))[:3]

# collinearity of the sensitivity directions (design property) — NOT a parameter correlation
F4 = J4.T @ J4
d = np.sqrt(np.diag(F4))
jac_cos = F4 / np.outer(d, d)

# correlation of the ESTIMATORS, from the inverse Fisher covariance
cov4 = sigma ** 2 * np.linalg.inv(F4)
s4 = np.sqrt(np.diag(cov4))
par_corr = cov4 / np.outer(s4, s4)

print("\n=== CRLB (σ=0.10): kb fixed vs kb marginalised ===")
print(f"{'coef':>8} | {'priorSD':>8} | {'CRLB(kb fixed)':>14} | {'CRLB(kb marg)':>14} | "
      f"{'inflation':>9} | {'cos(J_kw,J_kb)':>14} | corr(kw_hat,kb_hat)")
for j, z in enumerate(ZKEYS):
    infl = crlb_marg[j] / crlb_fixed[j]
    print(f"{z:>8} | {PRIOR_SD[z]:8.3f} | {crlb_fixed[j]:14.3f} | {crlb_marg[j]:14.3f} | "
          f"{infl:9.1f}x | {jac_cos[j, 3]:+14.3f} | {par_corr[j, 3]:+.3f}")
print("cos(J_kw,J_kb) measures how collinear the two sensitivity directions are (design property);")
print("corr(kw_hat,kb_hat) is the correlation of the estimates. They are different quantities and")
print("generally carry opposite signs — do not report the first one as a parameter correlation.")
