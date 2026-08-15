"""Step 4d: robust displaced-prior experiment (option d = a + b).

If a calibration only ever looks good because the prior happens to be centred on the truth, moving
the prior off the truth will expose it: an informative dataset pulls the estimate back, an
uninformative one leaves it at the displaced midpoint. Two displacement designs are used, each
keeping the prior width fixed and the truth inside the range:

  DOWN    : all three midpoints set to  truth - 1 prior SD  (into the strong-decay regime)
  OLDUP   : old midpoint set to  truth + 1 prior SD  (weaker/steep side; the upper bound is capped
            at -0.005 to stay non-positive); avg/new displaced DOWN as above.

Both designs displace all three coefficients, so gap_closed is well defined for every group; only
the DIRECTION of the old displacement differs between them.

PRIMARY: formal censored likelihood (no threshold). COMPARATOR: informal GLUE at both thresholds.
The comparison is the point of running both — the pull-back is a statement about how much
information the weighting extracts, so measuring it with the informal score answers a different
question from the one the section asks.

Each design's 8192 forward simulations are run ONCE (monitors only) and cached with the design box
as part of the cache key; then the gap-closing statistic is recomputed over N_NOISE independent
noisy observation sets. Robustness is reported as median [IQR] across the noise realisations.

gap_closed = (posterior_mean - displaced_mid) / (truth - displaced_mid)
  0  = the estimate stayed at the displaced prior midpoint (data uninformative)
  1  = the estimate was pulled all the way to the truth (data fully informative)
"""
import os
import json
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import wq_common as B
import provenance

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures")
CACHEDIR = os.path.join(HERE, "baseline_cache")
os.makedirs(FIGDIR, exist_ok=True)

cache = np.load(os.path.join(CACHEDIR, "baseline.npz"), allow_pickle=True)
mon_pos = list(cache["mon_pos"])
truth_mon = cache["truth_all"][:, mon_pos]          # noise-free monitor truth

TRUE = {"old": B.KW_OLD_TRUE, "avg": B.KW_AVG_TRUE, "new": B.KW_NEW_TRUE}
GROUPS = ["old", "avg", "new"]
THRESHOLDS = [B.RMSE_THR, B.RMSE_THR_DRAFT]
N_NOISE = 30
NOISE_SEEDS = list(range(42, 42 + N_NOISE))
UPPER_CAP = -0.005
PRIMARY = B.PRIMARY_WEIGHTING
COMPARATORS = [f"informal_glue@{t}" for t in THRESHOLDS]


def prior_sd(rng):
    a, b = rng
    return (b - a) / np.sqrt(12)


def down_range(orig, truth):
    a, b = orig
    w = b - a
    mid = truth - w / np.sqrt(12)
    return (mid - w / 2, mid + w / 2)


CAP_ACTIVATIONS = {}                   # group -> how far the cap had to shift the window, m/day


def up_range(orig, truth, cap=UPPER_CAP, group=None):
    """Displace upward by one prior SD, then shift down if that would go non-negative.

    The shift is a GUARD: it silently changes the design when it binds, so how much it moved and
    for which group is recorded in CAP_ACTIVATIONS and written to the artifact. A guard whose
    activation is not in the record cannot be audited from the artifact alone.
    """
    a, b = orig
    w = b - a
    mid = truth + w / np.sqrt(12)
    lo, hi = mid - w / 2, mid + w / 2
    if group is not None:
        CAP_ACTIVATIONS[group] = 0.0
    if hi > cap:                       # shift the window down to keep it non-positive
        shift = hi - cap
        if group is not None:
            CAP_ACTIVATIONS[group] = float(shift)
        lo, hi = lo - shift, hi - shift
    return (lo, hi)


def run_design(name, priors):
    """Sobol draws from the displaced box + cached monitor predictions + midpoints.

    Sobol rather than pseudo-random, and the same scramble seed as the baseline: otherwise a
    difference between this design and the baseline mixes a change of prior with a change of
    sampler. The prediction library is keyed on the box itself, so editing a design invalidates its
    own cache instead of silently reusing the previous one.
    """
    S = B.sobol_draws(priors)
    path = os.path.join(CACHEDIR, f"step4d_preds_{name}.npy")
    box = [float(v) for g in GROUPS for v in priors[g]]
    preds = provenance.load_keyed_array(path, design=name, box=box, n_mc=B.N_MC)
    if preds is None:
        preds = np.empty((B.N_MC, truth_mon.shape[0] - B.WARMUP_H, len(B.MONITOR_NODES)),
                         dtype=np.float32)
        t0 = time.time()
        for s in range(B.N_MC):
            preds[s] = B.simulate_chlorine(
                B.KB_FIXED, 0.0,
                pre_run=B.make_kw_hook(S["old"][s], S["avg"][s], S["new"][s]),
            ).values[B.WARMUP_H:].astype(np.float32)
            if (s + 1) % 1000 == 0:
                print(f"  [{name}] {s + 1}/{B.N_MC} ({time.time() - t0:.0f}s)")
        provenance.save_keyed_array(path, preds, design=name, box=box, n_mc=B.N_MC)
    else:
        print(f"  [{name}] reusing cached prediction library {preds.shape}")
    mids = {g: 0.5 * (priors[g][0] + priors[g][1]) for g in GROUPS}
    return S, preds.astype(np.float64), mids


def evaluate(S, preds, mids, priors):
    """Per-scheme, per-group arrays of gap_closed / sd_retained across the noisy observation sets."""
    schemes = [PRIMARY] + COMPARATORS
    out = {s: {g: {"gap": [], "sd_ret": []} for g in GROUPS} for s in schemes}
    diag = {s: {"ess": [], "retention": []} for s in schemes}
    psd = {g: prior_sd(priors[g]) for g in GROUPS}
    for seed in NOISE_SEEDS:
        rng = np.random.default_rng(seed)
        obs = np.clip(truth_mon + rng.normal(0, B.SIGMA_OBS, truth_mon.shape), 0, None)
        obs = obs[B.WARMUP_H:]
        rmse = B.rmse_of(preds, obs)
        w_primary, d_primary = B.all_weightings(preds, obs, schemes=[PRIMARY])[PRIMARY]
        per_scheme = {PRIMARY: (w_primary, d_primary, None)}
        for thr in THRESHOLDS:
            w, d = B.all_weightings(preds, obs, threshold=thr,
                                    schemes=["informal_glue"])["informal_glue"]
            per_scheme[f"informal_glue@{thr}"] = (w, d, float((rmse < thr).mean()))
        for s, (w, d, ret) in per_scheme.items():
            if w is None:
                continue
            diag[s]["ess"].append(d["ess"])
            diag[s]["retention"].append(d["ess_frac"] if ret is None else ret)
            for g in GROUPS:
                m, sd = B.weighted_mean_sd(w, S[g])
                denom = TRUE[g] - mids[g]
                out[s][g]["gap"].append((m - mids[g]) / denom if abs(denom) > 1e-9 else np.nan)
                out[s][g]["sd_ret"].append(sd / psd[g])
    return out, diag


def med_iqr(a):
    a = np.asarray(a, float)
    a = a[~np.isnan(a)]
    if len(a) == 0:
        return (np.nan, np.nan, np.nan)
    return (float(np.median(a)), float(np.percentile(a, 25)), float(np.percentile(a, 75)))


# ---- designs (both displace all three so every gap-closed is valid) ----
PRIORS_DOWN = {g: down_range(B.PRIOR[g], TRUE[g]) for g in GROUPS}
PRIORS_OLDUP = {"old": up_range(B.PRIOR["old"], TRUE["old"], group="old"),
                "avg": down_range(B.PRIOR["avg"], TRUE["avg"]),
                "new": down_range(B.PRIOR["new"], TRUE["new"])}
DESIGNS = [("DOWN", PRIORS_DOWN), ("OLDUP", PRIORS_OLDUP)]
for name, pr in DESIGNS:
    for g in GROUPS:
        a, b = pr[g]
        assert b <= 1e-9, f"{name}/{g} upper bound positive"
        assert a <= TRUE[g] <= b, f"{name}/{g} truth not in range"

report = {**B.weighting_provenance(comparators=COMPARATORS),
          # the upward-displacement guard, and whether it bound: a cap that silently moved the
          # design has to be visible in the artifact, not only in the resulting prior box
          "upper_cap_m_per_day": UPPER_CAP,
          "upper_cap_shift_applied": CAP_ACTIVATIONS,
          "upper_cap_note": "the OLDUP design displaces upward by one prior SD; if that would put "
                            "the window above the cap it is shifted down, which changes the design. "
                            "A non-zero shift here means the reported OLDUP prior is NOT simply "
                            "truth + 1 prior SD.",
          "n_noise": N_NOISE, "informal_thresholds": THRESHOLDS,
          "sampler": "scrambled Sobol from the displaced box, same seed as the baseline",
          "designs": {}}
box = {}   # for the figure: (scheme, design) -> per-group gap samples

for name, priors in DESIGNS:
    S, preds, mids = run_design(name, priors)
    ev, diag = evaluate(S, preds, mids, priors)
    d_out = {"priors": {g: list(priors[g]) for g in GROUPS}, "midpoints": mids, "by_scheme": {}}
    for s in [PRIMARY] + COMPARATORS:
        gstats = {}
        for g in GROUPS:
            gm, glo, ghi = med_iqr(ev[s][g]["gap"])
            sm, slo, shi = med_iqr(ev[s][g]["sd_ret"])
            gstats[g] = {"gap_med": gm, "gap_iqr": [glo, ghi],
                         "sd_ret_med": sm, "sd_ret_iqr": [slo, shi]}
            box[(s, name, g)] = [x for x in ev[s][g]["gap"] if not np.isnan(x)]
        rm, rlo, rhi = med_iqr(diag[s]["retention"])
        d_out["by_scheme"][s] = {"retention_med": rm, "retention_iqr": [rlo, rhi],
                                 "ess_med": med_iqr(diag[s]["ess"])[0], "groups": gstats}
    report["designs"][name] = d_out

with open(os.path.join(CACHEDIR, "step4d_displaced_robust.json"), "w") as f:
    json.dump(report, f, indent=2)

# ---- print summary ----
for name in ["DOWN", "OLDUP"]:
    d = report["designs"][name]
    print(f"\n=== {name} ===  priors: " +
          ", ".join(f"{g}[{d['priors'][g][0]:.3f},{d['priors'][g][1]:.3f}]" for g in GROUPS))
    for s in [PRIMARY] + COMPARATORS:
        b = d["by_scheme"][s]
        tag = "PRIMARY" if s == PRIMARY else "comparator"
        print(f"  {s} ({tag}) ESS {b['ess_med']:.0f}:")
        for g in GROUPS:
            q = b["groups"][g]
            print(f"    {g:>3}: gap closed {q['gap_med'] * 100:5.0f}% "
                  f"[{q['gap_iqr'][0] * 100:.0f}-{q['gap_iqr'][1] * 100:.0f}%]  "
                  f"SD retained {q['sd_ret_med'] * 100:4.0f}%")
print("\ngap closed = fraction of the displaced-midpoint -> truth distance recovered by the data.")
print("Read the primary row: the comparator rows show how much of any apparent prior dependence is")
print("the informal score's flatness rather than a limit of the observations.")

# ---- figure: gap_closed under the primary rule, with the comparator beside it ----
PANELS = [(PRIMARY, "formal censored likelihood (PRIMARY)"),
          (COMPARATORS[0], f"informal GLUE, threshold {THRESHOLDS[0]} (comparator)")]
fig, axes = plt.subplots(1, len(PANELS), figsize=(15, 4.8), sharey=True)
series = [(g, dn) for g in GROUPS for dn in ("DOWN", "OLDUP")]
for ax, (s, title) in zip(axes, PANELS):
    ax.boxplot([box[(s, dn, g)] for g, dn in series],
               tick_labels=[f"{g} ({dn})" for g, dn in series], showmeans=True)
    ax.axhline(0.0, color="gray", ls=":", lw=1, label="0% (stayed at prior mid)")
    ax.axhline(1.0, color="green", ls="--", lw=1, label="100% (reached truth)")
    ax.set_title(title, fontsize=10)
    ax.set_ylabel("gap to truth closed")
    ax.tick_params(axis="x", labelrotation=20)
    ax.grid(alpha=0.3, axis="y")
axes[0].legend(fontsize=8, loc="upper left")
fig.suptitle(f"Displaced-prior pull-back over {N_NOISE} noise realisations "
             "(all three coefficients displaced; primary rule vs the draft's comparator)", y=1.02)
plt.tight_layout()
figpath = os.path.join(FIGDIR, "step4d_displaced_robust.png")
plt.savefig(figpath, dpi=130, bbox_inches="tight")
print("\nfigure saved to", figpath)
