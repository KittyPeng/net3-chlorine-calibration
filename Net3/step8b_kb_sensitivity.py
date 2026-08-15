"""Step 8b: calibration under a misspecified k_b (± 20%) — Priority-2 #5.

The observations are generated with the true k_b = -0.5. The calibration is then re-run while FIXING
k_b at -0.4, -0.5, -0.6 (±20%) to see how the three grouped k_w move to compensate: the bulk-wall
trade-off. This is the empirical partner to Fisher Case B, which predicts the same trade-off from
the Jacobian alone.

PRIMARY: formal censored likelihood. COMPARATOR: informal GLUE at the primary threshold. Running
both matters here because the whole quantity of interest is a systematic displacement, and a flat
score under-reports displacement — the same effect Step 8 finds for sensor bias.

k_b = -0.5 reuses the baseline cache; k_b = -0.4 / -0.6 are simulated here over all 92 nodes (the
risk field needs them) and cached, so a re-run costs nothing. 30 noise realisations; the reported
coefficient is the median across them, and the risk ranking is taken from the median risk field
across the same realisations rather than from one arbitrary seed.
"""
import os
import json
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import wq_common as B
import provenance

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures")
CACHEDIR = os.path.join(HERE, "baseline_cache")
ZKEYS = ["old", "average", "new"]
JKEY = {"old": "old", "average": "avg", "new": "new"}     # report key per zone
TRUE = {"old": B.KW_OLD_TRUE, "average": B.KW_AVG_TRUE, "new": B.KW_NEW_TRUE}
C_MIN = 0.2

cache = np.load(os.path.join(CACHEDIR, "baseline.npz"), allow_pickle=True)
ALL_NODES = list(cache["all_nodes"])
mon_pos = list(cache["mon_pos"])
truth_mon = cache["truth_all"][:, mon_pos]                # true field (k_b = -0.5) at the monitors
S = {"old": cache["S_old"], "average": cache["S_avg"], "new": cache["S_new"]}
N_MC = B.N_MC
Tn = cache["obs_glue"].shape[0]

KBS = [-0.4, -0.5, -0.6]
N_NOISE = 30
SEEDS = list(range(42, 42 + N_NOISE))
SCHEMES = [B.PRIMARY_WEIGHTING, "informal_glue"]
TOP_K = 6

# candidate predictions per k_b (all nodes); reuse the baseline cache for k_b = -0.5
preds = {-0.5: cache["C_all"]}
for kb in [-0.4, -0.6]:
    path = os.path.join(CACHEDIR, f"step8b_preds_kb{kb}.npy")
    Ck = provenance.load_keyed_array(path, kb=kb, n_mc=N_MC)
    if Ck is None:
        Ck = np.empty((N_MC, Tn, len(ALL_NODES)), dtype=np.float32)
        t0 = time.time()
        for s in range(N_MC):
            Ck[s] = B.simulate_chlorine(kb, 0.0,
                                        pre_run=B.make_kw_hook(S["old"][s], S["average"][s],
                                                               S["new"][s]),
                                        monitor_nodes=ALL_NODES).values[B.WARMUP_H:].astype(
                                            np.float32)
            if (s + 1) % 1000 == 0:
                print(f"  kb={kb}: {s + 1}/{N_MC} ({time.time() - t0:.0f}s)")
        provenance.save_keyed_array(path, Ck, kb=kb, n_mc=N_MC)
    else:
        print(f"  kb={kb}: reusing cached prediction library {Ck.shape}")
    preds[kb] = Ck


def med(a):
    return float(np.median(a))


rows = []
for kb in KBS:
    Call = preds[kb].astype(np.float64)
    Cmon = Call[:, :, mon_pos]
    dur = B.below_threshold_duration(Call, C_MIN)
    # TWO risk metrics, because they rank differently and the k_b conclusion is load-bearing.
    #   P_bar = expected fraction of the window below C_MIN  -> "how long"
    #   E[A]  = expected cumulative deficit, trapezoid in h  -> "how long AND how far below"
    # Step 10's headline top-10 is ranked by E[A], so a statement about "the hot-spot list" that
    # rested only on P_bar would not be about the list the risk chapter actually publishes.
    defc = np.trapezoid(np.clip(C_MIN - Call, 0.0, None), dx=1.0, axis=1)   # (N_MC, nodes), mg/L·h
    acc = {s: {z: [] for z in ZKEYS} for s in SCHEMES}
    sd_acc = {s: {z: [] for z in ZKEYS} for s in SCHEMES}
    ess_acc = {s: [] for s in SCHEMES}
    P_acc = {s: [] for s in SCHEMES}
    A_acc = {s: [] for s in SCHEMES}
    rmse_mins, noise_rmses = [], []
    for seed in SEEDS:
        rng = np.random.default_rng(seed)
        obs = np.clip(truth_mon + rng.normal(0, B.SIGMA_OBS, truth_mon.shape), 0, None)[B.WARMUP_H:]
        # Best achievable aggregate fit over the candidate library BUILT AT THIS k_b, against
        # observations generated at the true k_b, and the realised noise floor of the same draw.
        # This is what decides whether a wrong fixed bulk coefficient is visible in the aggregate
        # residual at all, or is absorbed into the wall coefficients without trace.
        rmse_mins.append(float(B.rmse_of(Cmon, obs).min()))
        noise_rmses.append(float(np.sqrt(((obs - truth_mon[B.WARMUP_H:]) ** 2).mean())))
        wts = B.all_weightings(Cmon, obs, threshold=B.RMSE_THR, schemes=SCHEMES)
        for s in SCHEMES:
            w, diag = wts[s]
            if w is None:
                continue
            ess_acc[s].append(diag["ess"])
            for z in ZKEYS:
                m, sd = B.weighted_mean_sd(w, S[z])
                acc[s][z].append(m)
                sd_acc[s][z].append(sd)
            # the risk field per realisation, so the ranking below is a median over the noise
            # rather than whatever one arbitrary seed happened to give
            P_acc[s].append(np.tensordot(w, dur, axes=(0, 0)) / B.window_hours(Call))
            A_acc[s].append(np.tensordot(w, defc, axes=(0, 0)))
    row = {"kb": kb,
           "rmse_min_med": med(rmse_mins),
           "noise_rmse_med": med(noise_rmses),
           "rmse_min_over_noise_med": med([a / b for a, b in zip(rmse_mins, noise_rmses)]),
           "by_scheme": {}}
    for s in SCHEMES:
        P_med = np.median(np.vstack(P_acc[s]), axis=0)
        A_med = np.median(np.vstack(A_acc[s]), axis=0)
        order = np.argsort(P_med)[::-1]
        order_A = np.argsort(A_med)[::-1]
        row["by_scheme"][s] = {
            **{JKEY[z]: med(acc[s][z]) for z in ZKEYS},
            **{JKEY[z] + "_sd": med(sd_acc[s][z]) for z in ZKEYS},
            "ess_med": med(ess_acc[s]),
            f"risk_top{TOP_K}": [ALL_NODES[i] for i in order[:TOP_K]],
            f"deficit_top{TOP_K}": [ALL_NODES[i] for i in order_A[:TOP_K]],
            "_P": P_med, "_A": A_med}
    rows.append(row)

# shift of each coefficient relative to the reference k_b, which is what the log tabulates
ref = next(r for r in rows if r["kb"] == B.KB_FIXED)
for r in rows:
    for s in SCHEMES:
        a, b = r["by_scheme"][s], ref["by_scheme"][s]
        a["shift_vs_kb_ref"] = {k: a[k] - b[k] for k in ("old", "avg", "new")}
        a["shift_over_own_sd"] = {k: (a[k] - b[k]) / b[k + "_sd"] if b[k + "_sd"] else None
                                  for k in ("old", "avg", "new")}
        a["risk_spearman_vs_kb_ref"] = float(spearmanr(a["_P"], b["_P"]).statistic)
        top, reftop = set(a[f"risk_top{TOP_K}"]), set(b[f"risk_top{TOP_K}"])
        a[f"risk_top{TOP_K}_jaccard_vs_kb_ref"] = len(top & reftop) / len(top | reftop)
        # the same two comparisons on the deficit field Step 10 actually ranks by
        a["deficit_spearman_vs_kb_ref"] = float(spearmanr(a["_A"], b["_A"]).statistic)
        topA, refA = set(a[f"deficit_top{TOP_K}"]), set(b[f"deficit_top{TOP_K}"])
        a[f"deficit_top{TOP_K}_jaccard_vs_kb_ref"] = len(topA & refA) / len(topA | refA)

# ---- how deep does the reordering go? -----------------------------------------------------------
# TOP_K = 6 is a chosen cut-off, so a Jaccard computed only there cannot distinguish "the leading
# nodes genuinely reshuffled" from "two nodes swapped across the 6th/7th boundary". Three diagnostics
# separate them: the same comparison at several k; the rank each node holds in both fields; and,
# decisively, the rank in the CORRECTLY SPECIFIED field of every node that enters the leading set. A
# node entering the top 6 from reference rank 7 is a cut-off effect; one entering from rank 20 is not.
TOP_KS = [3, 5, 6, 10, 15]


def ranks_of(P):
    """rank 1 = highest risk."""
    order = np.argsort(P)[::-1]
    r = np.empty(len(P), dtype=int)
    r[order] = np.arange(1, len(P) + 1)
    return r


def depth_diagnostics(Pa, Pb):
    """top-k Jaccard at several k, rank changes, and where each entering node came from."""
    ra, rb = ranks_of(Pa), ranks_of(Pb)
    oa, ob = np.argsort(Pa)[::-1], np.argsort(Pb)[::-1]
    d = np.abs(ra - rb)
    return {
        "topk_jaccard_vs_kb_ref": {
            f"top{k}": len(set(oa[:k]) & set(ob[:k])) / len(set(oa[:k]) | set(ob[:k]))
            for k in TOP_KS},
        "entered_top6": [{"node": str(ALL_NODES[i]), "rank_here": int(ra[i]),
                          "reference_rank": int(rb[i])}
                         for i in oa[:TOP_K] if i not in set(ob[:TOP_K])],
        "left_top6": [{"node": str(ALL_NODES[i]), "reference_rank": int(rb[i]),
                       "rank_here": int(ra[i])}
                      for i in ob[:TOP_K] if i not in set(oa[:TOP_K])],
        "rank_change": {"max_over_92": int(d.max()), "median_over_92": float(np.median(d)),
                        "max_within_reference_top15": int(d[ob[:15]].max())}}


for r in rows:
    for s in SCHEMES:
        a, b = r["by_scheme"][s], ref["by_scheme"][s]
        # both metrics get the same treatment; a depth claim made on P_bar is not a depth claim
        # about the list Step 10 publishes
        a["depth_deficit"] = depth_diagnostics(a["_A"], b["_A"])
        Pa, Pb = a["_P"], b["_P"]
        ra, rb = ranks_of(Pa), ranks_of(Pb)
        oa, ob = np.argsort(Pa)[::-1], np.argsort(Pb)[::-1]
        a["topk_jaccard_vs_kb_ref"] = {
            f"top{k}": len(set(oa[:k]) & set(ob[:k])) / len(set(oa[:k]) | set(ob[:k]))
            for k in TOP_KS}
        a["entered_top6"] = [{"node": str(ALL_NODES[i]), "rank_here": int(ra[i]),
                              "reference_rank": int(rb[i])}
                             for i in oa[:TOP_K] if i not in set(ob[:TOP_K])]
        a["left_top6"] = [{"node": str(ALL_NODES[i]), "reference_rank": int(rb[i]),
                           "rank_here": int(ra[i])}
                          for i in ob[:TOP_K] if i not in set(oa[:TOP_K])]
        d = np.abs(ra - rb)
        a["rank_change"] = {"max_over_92": int(d.max()), "median_over_92": float(np.median(d)),
                            "max_within_reference_top15": int(d[ob[:15]].max())}

# Is rank 6 a cliff or a plateau? If consecutive risk values near the cut-off are nearly equal, a
# small perturbation reorders them without changing the risk field in any material way.
refP = ref["by_scheme"][B.PRIMARY_WEIGHTING]["_P"]
o_ref = np.argsort(refP)[::-1][:12]
reference_profile = [{"rank": i + 1, "node": str(ALL_NODES[j]), "risk": float(refP[j])}
                     for i, j in enumerate(o_ref)]

print("\n=== Step 8b: calibration under a misspecified k_b (observations generated at k_b = -0.5) ===")
for s in SCHEMES:
    tag = "PRIMARY" if s == B.PRIMARY_WEIGHTING else "comparator"
    print(f"\n--- {s} ({tag}) ---")
    print(f"{'k_b':>5} {'ESS':>7} | {'old':>8} {'shift/SD':>9} | {'avg':>8} {'shift/SD':>9} | "
          f"{'new':>8} {'shift/SD':>9} | {'rho_risk':>8} {'J6':>5}")
    for r in rows:
        a = r["by_scheme"][s]
        sos = a["shift_over_own_sd"]
        print(f"{r['kb']:>5} {a['ess_med']:>7.0f} | {a['old']:>8.3f} {sos['old']:>+9.2f} | "
              f"{a['avg']:>8.3f} {sos['avg']:>+9.2f} | {a['new']:>8.3f} {sos['new']:>+9.2f} | "
              f"{a['risk_spearman_vs_kb_ref']:>8.3f} "
              f"{a[f'risk_top{TOP_K}_jaccard_vs_kb_ref']:>5.2f}")
    b = ref["by_scheme"][s]
    print(f"      posterior SD at the true k_b: " +
          ", ".join(f"{JKEY[z]} {b[JKEY[z] + '_sd']:.4f}" for z in ZKEYS))
print("\n--- the same comparison on the two risk metrics (PRIMARY scheme) ---")
print("P_bar = expected fraction of the window below 0.2; E[A] = expected cumulative deficit,")
print("which is what Step 10's headline top-10 is ranked by.")
print(f"{'k_b':>5} | {'rho(P_bar)':>10} {'J6(P_bar)':>9} | {'rho(E[A])':>10} {'J6(E[A])':>9}")
for r in rows:
    a = r["by_scheme"][B.PRIMARY_WEIGHTING]
    print(f"{r['kb']:>5} | {a['risk_spearman_vs_kb_ref']:>10.3f} "
          f"{a[f'risk_top{TOP_K}_jaccard_vs_kb_ref']:>9.2f} | "
          f"{a['deficit_spearman_vs_kb_ref']:>10.3f} "
          f"{a[f'deficit_top{TOP_K}_jaccard_vs_kb_ref']:>9.2f}")
for r in rows:
    a = r["by_scheme"][B.PRIMARY_WEIGHTING]
    print(f"  k_b={r['kb']}: P_bar top6 {a[f'risk_top{TOP_K}']}")
    print(f"  k_b={r['kb']}: E[A]  top6 {a[f'deficit_top{TOP_K}']}")

print("\nshift/SD = displacement of the coefficient caused by the wrong k_b, in units of that")
print(f"coefficient's own posterior SD at the true k_b; rho_risk / J{TOP_K} compare the 92-node risk")
print("field and its leading set against the correctly specified case.")

print("\n--- is the top-6 turnover a cut-off effect or a real reshuffle? (PRIMARY scheme) ---")
print("reference risk profile at the true k_b (rank, node, expected fraction of hours below 0.2):")
for e in reference_profile:
    mark = "   <- cut-off" if e["rank"] == TOP_K else ""
    print(f"  {e['rank']:>2}  node {e['node']:>4}  {e['risk']:.4f}{mark}")
print("\nJaccard against the correctly specified field at several k:")
hdr = "  ".join(f"top{k:<2}" for k in TOP_KS)
print(f"{'k_b':>5}  {hdr}   max|Δrank| (all 92 / ref top-15)")
for r in rows:
    a = r["by_scheme"][B.PRIMARY_WEIGHTING]
    js = "  ".join(f"{a['topk_jaccard_vs_kb_ref'][f'top{k}']:<5.2f}" for k in TOP_KS)
    print(f"{r['kb']:>5}  {js}   {a['rank_change']['max_over_92']:>3} / "
          f"{a['rank_change']['max_within_reference_top15']:>3}")
for r in rows:
    if r["kb"] == B.KB_FIXED:
        continue
    a = r["by_scheme"][B.PRIMARY_WEIGHTING]
    ent = ", ".join(f"{e['node']} (was rank {e['reference_rank']})" for e in a["entered_top6"])
    lef = ", ".join(f"{e['node']} (now rank {e['rank_here']})" for e in a["left_top6"])
    print(f"  k_b={r['kb']}: entered top-6: {ent or 'none'} | left: {lef or 'none'}")
print("A node entering from reference rank 7-8 is a cut-off effect; one entering from far down is a")
print("genuine reshuffle of the leading nodes. Read this together with the spacing table above.")

for r in rows:
    for s in SCHEMES:
        del r["by_scheme"][s]["_P"]
        del r["by_scheme"][s]["_A"]
report = {**B.weighting_provenance(comparators=["informal_glue"]),
          "kbs": KBS, "kb_true": B.KB_FIXED, "n_noise": N_NOISE,
          "informal_threshold": B.RMSE_THR,
          "risk_ranking_basis": f"median risk field over the {N_NOISE} noise realisations",
          "risk_metric": {
              "P_bar": "likelihood-weighted expected fraction of the 48 h window below 0.2 mg/L "
                       "(risk_top6, risk_spearman_*): how LONG a node sits below the threshold",
              "E_A": "likelihood-weighted expected cumulative deficit, trapezoid over the 48 h "
                     "window in mg/L*h (deficit_top6, deficit_spearman_*): how long AND how far "
                     "below. This is the metric Step 10's headline top-10 is ranked by, so a claim "
                     "about 'the hot-spot list' has to be made on it and not on P_bar alone."},
          "topk_note": "TOP_K = 6 is a chosen cut-off; the same comparison is reported at k = 3, 5, "
                       "6, 10, 15, together with per-node rank changes and the reference rank of "
                       "every node that enters the leading set, so a cut-off effect can be told "
                       "apart from a reshuffle of the leading nodes",
          "reference_risk_profile_top12": reference_profile,
          "rows": rows}


def _jsafe(o):
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.integer):
        return int(o)
    return str(o)


with open(os.path.join(CACHEDIR, "step8b_kb_sensitivity.json"), "w") as f:
    json.dump(report, f, indent=2, default=_jsafe)

# ---- figure: posterior mean of each k_w vs the fixed k_b, primary vs comparator ----
fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), sharey=False)
kbs = [r["kb"] for r in rows]
for ax, s in zip(axes, SCHEMES):
    for z, c in zip(ZKEYS, ["tab:red", "tab:orange", "tab:green"]):
        ax.plot(kbs, [r["by_scheme"][s][JKEY[z]] for r in rows], marker="o", color=c, label=z)
        ax.axhline(TRUE[z], color=c, ls=":", lw=1)
    ax.axvline(B.KB_FIXED, color="gray", ls="--", lw=1, label="true k_b")
    ax.set_xlabel("fixed k_b in calibration (1/day)")
    ax.set_ylabel("posterior mean of k_w (m/day)")
    ax.set_title(("formal censored likelihood (PRIMARY)" if s == B.PRIMARY_WEIGHTING
                  else f"informal GLUE, threshold {B.RMSE_THR} (comparator)"), fontsize=10)
    ax.grid(alpha=0.3)
axes[0].legend(fontsize=8)
fig.suptitle("Step 8b — bulk-wall compensation: misspecifying k_b by ±20% displaces the k_w "
             "estimates\n(dotted = true k_w; the risk ranking is checked separately)", y=1.03)
plt.tight_layout()
figpath = os.path.join(FIGDIR, "step8b_kb_sensitivity.png")
plt.savefig(figpath, dpi=130, bbox_inches="tight")
print("figure saved to", figpath)
