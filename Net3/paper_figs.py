"""Main-text figures and tables for the paper.

Every number is read from baseline_cache/ (JSON artifacts and baseline.npz). Nothing is
recomputed from EPANET and nothing is copied from the superseded Word draft, so a figure
can never drift from the artifact that produced it.

Style is fixed here, once, for all figures: one palette, one set of zone colours, one set
of inference-scheme colours, one grid weight.

    zones    #D55E00 old      #0072B2 average   #009E73 new
    schemes  #3F51B5 formal   #C25A00 informal@0.107   #ED9B4A informal@0.120

#ED9B4A is the lightest of the six and reads faintly against white, so every panel that
uses it also carries a legend and an in-panel direct label, and Table 3 gives the same
numbers as text.

Usage:
    python paper_figs.py all
    python paper_figs.py fig2 fig3
"""
import json
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "baseline_cache")
OUTDIR = os.path.join(HERE, "figures", "paper")

# ==================== locked style ====================

ZONE = {"old": "#D55E00", "average": "#0072B2", "new": "#009E73"}
ZONE_KEY = {"old": "old", "avg": "average", "average": "average", "new": "new"}
SCHEME = {
    "formal_censored": "#3F51B5",
    "informal_glue": "#C25A00",          # threshold 0.107
    "informal_glue_draft_thr": "#ED9B4A",  # threshold 0.120
}
SCHEME_LABEL = {
    "formal_censored": "Formal censored likelihood",
    "informal_glue": "Informal GLUE (RMSE < 0.107)",
    "informal_glue_draft_thr": "Informal GLUE (RMSE < 0.120)",
}
PRIOR_FILL, PRIOR_EDGE = "#E4E2DD", "#A8A59F"
INK, INK2, INK3 = "#1A1A19", "#55534E", "#8A8781"
GRID = "#E8E6E1"

# main text is 12 pt at 1.5 spacing on A4 with 25 mm margins, so the text column is
# 160 mm; a full-width figure is 160 mm = 6.30 in.
W_FULL, W_HALF = 6.30, 3.10

plt.rcParams.update({
    "figure.dpi": 110,
    "savefig.dpi": 400,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 7.5,
    "axes.titlesize": 8,
    "axes.labelsize": 7.5,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "axes.edgecolor": INK3,
    "axes.linewidth": 0.6,
    "axes.labelcolor": INK,
    "axes.titlecolor": INK,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "text.color": INK,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "legend.frameon": False,
    "lines.linewidth": 1.4,
    "lines.markersize": 4.5,
})

PARAMS = [("old", "$k_{w,\\mathrm{old}}$"), ("avg", "$k_{w,\\mathrm{average}}$"),
          ("new", "$k_{w,\\mathrm{new}}$")]
PRIOR_BOX = {"old": (-1.5, -0.2), "avg": (-0.2, -0.04), "new": (-0.10, -0.005)}
TRUTH = {"old": -1.0, "avg": -0.1, "new": -0.05}


def load(name):
    with open(os.path.join(CACHE, name)) as fh:
        return json.load(fh)


def npz():
    """The prediction library. Absent from the published repository, so say why rather than
    letting a fresh clone meet a bare FileNotFoundError from numpy."""
    path = os.path.join(CACHE, "baseline.npz")
    if not os.path.exists(path):
        raise SystemExit(
            "baseline.npz is not in this checkout.\n"
            "  It is 118 MB, past GitHub's per-file limit, so it is not published with the code.\n"
            "  Regenerate it, from the Net3 directory:\n"
            "      python step1_freeze_baseline.py\n"
            "  Then check it against the copy the published results came from:\n"
            "      shasum -a 256 baseline_cache/baseline.npz\n"
            "  and compare with baseline_npz_sha256 in baseline_cache/cache_manifest.json.")
    return np.load(path)


def note(fig, text, y=0.02, size=6.4, width=118):
    """Figure footnote, wrapped to the canvas width.

    savefig.bbox is "tight", so a single line that overruns the figure widens the saved image
    instead of being clipped. Figure 4's note stretched it to 8.7:1 against a designed 3.1:1, and
    the panels shrank to fit the page. Wrapping keeps the canvas the size it was laid out at.
    """
    import textwrap
    fig.text(0.5, y, "\n".join(textwrap.wrap(" ".join(text.split()), width)),
             ha="center", va="bottom", fontsize=size, color=INK3, linespacing=1.35)


def finish(fig, stem):
    os.makedirs(OUTDIR, exist_ok=True)
    path = os.path.join(OUTDIR, stem + ".png")
    fig.savefig(path)
    plt.close(fig)
    print("wrote", os.path.relpath(path, HERE))
    return path


def tidy(ax, grid_axis="y"):
    ax.grid(True, axis=grid_axis, zorder=0)
    ax.set_axisbelow(True)


def weighted_density(x, w, grid, bounds, bw_scale=1.0):
    """Gaussian-kernel weighted density with reflection at the prior bounds.

    Without the reflection an almost-uniform sample (which is what the informal score
    produces) is rendered as falling to zero at the edges of its own support, so the
    figure would show structure the posterior does not have.
    """
    w = np.asarray(w, dtype=np.float64)
    w = w / w.sum()
    mu = np.sum(w * x)
    sd = np.sqrt(max(np.sum(w * (x - mu) ** 2), 1e-30))
    ess = 1.0 / np.sum(w ** 2)
    h = bw_scale * 1.06 * sd * ess ** (-0.2)
    lo, hi = bounds
    dens = np.zeros_like(grid)
    for src in (x, 2 * lo - x, 2 * hi - x):
        z = (grid[:, None] - src[None, :]) / h
        dens += (np.exp(-0.5 * z ** 2) @ w) / (h * np.sqrt(2 * np.pi))
    return dens


def weights_from_loglik(ll, mask=None):
    ll = np.asarray(ll, dtype=np.float64).copy()
    if mask is not None:
        ll[~mask] = -np.inf
    ll -= ll.max()
    w = np.exp(ll)
    return w / w.sum()


# ==================== Figure 2 — inference-rule dependence ====================

def fig2():
    d = npz()
    meta = load("baseline_meta.json")["summary"]
    S = {"old": d["S_old"], "avg": d["S_avg"], "new": d["S_new"]}
    rmse = d["RMSE"]

    wts = {
        "formal_censored": weights_from_loglik(d["loglik_censored"]),
        "informal_glue": weights_from_loglik(-0.5 * (rmse / 0.1) ** 2, rmse < 0.107),
        "informal_glue_draft_thr": weights_from_loglik(-0.5 * (rmse / 0.1) ** 2, rmse < 0.12),
    }

    fig = plt.figure(figsize=(W_FULL, 4.35))
    gs = fig.add_gridspec(3, 2, width_ratios=[2.45, 1.0], hspace=0.55, wspace=0.30)

    for row, (key, tex) in enumerate(PARAMS):
        ax = fig.add_subplot(gs[row, 0])
        lo, hi = PRIOR_BOX[key]
        grid = np.linspace(lo, hi, 600)
        prior_h = 1.0 / (hi - lo)

        ax.fill_between(grid, 0, np.full_like(grid, prior_h), color=PRIOR_FILL,
                        edgecolor=PRIOR_EDGE, linewidth=0.6, zorder=1)
        peak = prior_h
        for scheme in ("informal_glue_draft_thr", "informal_glue", "formal_censored"):
            dens = weighted_density(S[key], wts[scheme], grid, (lo, hi))
            ax.plot(grid, dens, color=SCHEME[scheme], zorder=3,
                    linewidth=1.6 if scheme == "formal_censored" else 1.2)
            peak = max(peak, dens.max())

        ax.axvline(TRUTH[key], color=INK, linestyle=(0, (3, 2)), linewidth=1.0, zorder=4)
        ax.set_ylim(0, peak * 1.32)
        ax.set_xlim(lo, hi)
        ax.set_yticks([])
        ax.set_ylabel("density", color=INK2)
        ax.spines["left"].set_visible(False)
        ax.set_title(tex, loc="left", pad=3)
        if row == 2:
            ax.set_xlabel("wall-decay coefficient (m day$^{-1}$)")

        # direct labels: the lightest informal step carries a contrast WARN, so it is
        # never identified by colour alone
        ax.annotate("truth", xy=(TRUTH[key], peak * 1.32), xytext=(3, -1),
                    textcoords="offset points", ha="left", va="top",
                    fontsize=6.5, color=INK2)
        if row == 0:
            dens = weighted_density(S[key], wts["formal_censored"], grid, (lo, hi))
            imax = int(dens.argmax())
            # left shoulder, so the label cannot collide with the truth rule at the peak
            ishoulder = int(np.argmin(np.abs(dens[:imax] - 0.55 * dens.max())))
            ax.annotate("formal", xy=(grid[ishoulder], dens[ishoulder]), xytext=(-3, 1),
                        textcoords="offset points", ha="right", va="bottom", fontsize=6.5,
                        color=SCHEME["formal_censored"], fontweight="bold")
            ax.annotate("informal", xy=(lo + 0.04 * (hi - lo), prior_h * 1.45),
                        fontsize=6.5, color=SCHEME["informal_glue"], fontweight="bold")
            ax.annotate("prior", xy=(hi - 0.02 * (hi - lo), prior_h * 0.42),
                        ha="right", fontsize=6.5, color=INK3)

        # ---- SD-retained bars ----
        bx = fig.add_subplot(gs[row, 1])
        order = ["formal_censored", "informal_glue", "informal_glue_draft_thr"]
        vals = [meta["schemes"][s]["coef"][key]["sd_retained"] * 100 for s in order]
        ypos = np.arange(len(order))[::-1]
        bx.barh(ypos, vals, height=0.62, color=[SCHEME[s] for s in order], zorder=2)
        for y, v in zip(ypos, vals):
            bx.text(v + 2.5, y, f"{v:.0f}%", va="center", fontsize=6.5, color=INK2)
        bx.set_xlim(0, 118)
        bx.set_ylim(-0.7, len(order) - 0.3)
        bx.set_yticks(ypos)
        bx.set_yticklabels(["formal", "GLUE 0.107", "GLUE 0.120"], fontsize=6.2)
        bx.tick_params(axis="y", length=0, pad=1)
        bx.spines["left"].set_visible(False)
        bx.set_xticks([0, 50, 100])
        bx.grid(True, axis="x", zorder=0)
        bx.set_axisbelow(True)
        bx.set_title("prior SD retained", loc="left", pad=3, fontsize=7.5, color=INK2)
        if row == 2:
            bx.set_xlabel("% of prior SD")

    handles = [Patch(facecolor=PRIOR_FILL, edgecolor=PRIOR_EDGE, label="Prior (uniform)")]
    handles += [Line2D([], [], color=SCHEME[s], linewidth=1.6, label=SCHEME_LABEL[s])
                for s in ("formal_censored", "informal_glue", "informal_glue_draft_thr")]
    handles += [Line2D([], [], color=INK, linestyle=(0, (3, 2)), label="Truth")]
    fig.legend(handles=handles, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.055),
               columnspacing=1.4, handlelength=1.8)
    return finish(fig, "fig2_inference_rule_dependence")


# ==================== Figure 3 — triangulated identifiability ====================

def fig3():
    meta = load("baseline_meta.json")["summary"]
    fisher = load("step7_fisher.json")["cases"]["A: kw only"]["coef"]
    prof = load("step7b_profile.json")["continuous_profile"]["by_likelihood"]["censored"]["coef"]
    rep = load("step14_repeated_noise.json")["by_scheme"]["formal_censored"]["coef"]

    name = {"old": "old", "avg": "average", "new": "new"}
    rows = [
        ("Formal ensemble, 5–95%", "ens"),
        ("Profile likelihood, 95%", "prof"),
        ("Fisher bound, 95%", "fish"),
        ("Repeated noise, mean ± SD", "rep"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(W_FULL, 1.95))
    for ax, (key, tex) in zip(axes, PARAMS):
        nm = name[key]
        c = ZONE[nm]
        ens = meta["schemes"]["formal_censored"]["coef"][key]
        f = fisher[nm]
        p = prof[nm]
        r = rep[nm]
        spec = {
            "ens": (ens["q05_50_95"][0], ens["q05_50_95"][2], ens["mean"]),
            "prof": (p["lo"], p["hi"], None),
            "fish": (TRUTH[key] - 1.96 * f["crlb"], TRUTH[key] + 1.96 * f["crlb"], None),
            "rep": (r["mean_of_posterior_mean"] - r["empirical_sd_of_posterior_mean"],
                    r["mean_of_posterior_mean"] + r["empirical_sd_of_posterior_mean"],
                    r["mean_of_posterior_mean"]),
        }
        for i, (_, tag) in enumerate(rows):
            y = len(rows) - 1 - i
            lo, hi, mid = spec[tag]
            ax.plot([lo, hi], [y, y], color=c, linewidth=1.8, solid_capstyle="round", zorder=3)
            ax.plot([lo, lo, hi, hi], [y - .14, y + .14, y - .14, y + .14], ls="none",
                    marker="|", color=c, markersize=5, markeredgewidth=1.2, zorder=3)
            if mid is not None:
                ax.plot([mid], [y], marker="o", color=c, markersize=4,
                        markeredgecolor="white", markeredgewidth=0.7, zorder=4)
        ax.axvline(TRUTH[key], color=INK, linestyle=(0, (3, 2)), linewidth=1.0, zorder=2)
        ax.set_ylim(-0.6, len(rows) - 0.4)
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels([r[0] for r in rows][::-1] if ax is axes[0] else [])
        ax.set_title(tex, loc="left", pad=3)
        ax.grid(True, axis="x", zorder=0)
        ax.set_axisbelow(True)
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
        ax.set_xlabel("m day$^{-1}$")

    axes[0].annotate("truth", xy=(TRUTH["old"], len(rows) - 0.45), xytext=(3, 0),
                     textcoords="offset points", fontsize=6.5, color=INK2)
    fig.tight_layout(w_pad=1.2, rect=(0, 0.06, 1, 1))
    fig.text(0.5, 0.015, "Nominal levels differ by construction and are labelled "
             "individually; the four rows are not a like-for-like nesting.",
             ha="center", fontsize=6.5, color=INK3)
    return finish(fig, "fig3_triangulated_identifiability")


# ==================== Figure 6 — prediction without identification ====================

def fig6():
    loo = load("step11_loo.json")
    lozo = [r for r in loo["leave_one_zone_out"] if r["scheme"] == "formal_censored"]
    meta = load("baseline_meta.json")["summary"]
    key_of = {"old": "old", "average": "avg", "new": "new"}
    Z = 1.645  # 90% of a Gaussian, to match the ensemble 5-95% convention

    fig, axes = plt.subplots(1, 4, figsize=(W_FULL, 2.75),
                             gridspec_kw={"width_ratios": [1, 1, 1, 1.25], "wspace": 0.62})
    fig.subplots_adjust(left=0.095, right=0.985, top=0.845, bottom=0.275)

    for ax, row in zip(axes[:3], lozo):
        zone = row["zone_dropped"]
        key = key_of[zone]
        c = ZONE[zone]
        lo, hi = PRIOR_BOX[key]
        mid = (lo + hi) / 2
        est = row["k_" + key][0]
        # the artifact stores the RETENTION ratio, not the posterior width itself: the
        # width is recovered as retention x prior SD, which is the quantity the panel
        # is about. row["k_*"][1:] are the 5th/95th percentiles of the posterior MEAN
        # across noise realisations and are not a posterior interval.
        sd = row["own_sd_retained"][0] * meta["prior_sd"][key]

        ax.fill_between([0.35, 1.65], lo, hi, color=PRIOR_FILL, edgecolor=PRIOR_EDGE,
                        linewidth=0.6, zorder=1)
        ax.plot([0.35, 1.65], [mid] * 2, color=PRIOR_EDGE, linewidth=0.9,
                linestyle=(0, (2, 2)), zorder=2)
        ax.plot([1, 1], [est - Z * sd, est + Z * sd], color=c, linewidth=3.0,
                solid_capstyle="butt", zorder=3)
        ax.plot([1], [est], marker="o", color=c, markersize=5.5,
                markeredgecolor="white", markeredgewidth=0.9, zorder=4)
        ax.plot([0.35, 1.65], [TRUTH[key]] * 2, color=INK, linestyle=(0, (3, 2)),
                linewidth=1.0, zorder=5)

        span = hi - lo
        ax.set_title(f"{zone}", pad=3)
        ax.set_xlim(0.35, 1.65)
        ax.set_ylim(lo - 0.13 * span, hi + 0.04 * span)
        ax.set_xticks([])
        ax.set_ylabel("m day$^{-1}$" if ax is axes[0] else "")
        ax.grid(True, axis="y", zorder=0)
        ax.set_axisbelow(True)
        ax.spines["bottom"].set_visible(False)
        ax.text(0.5, 0.012, f"{row['own_sd_retained'][0] * 100:.0f}% of prior SD",
                transform=ax.transAxes, ha="center", va="bottom",
                fontsize=6.5, color=c, fontweight="bold")

    # held-out predictive skill, same three cases
    ax = axes[3]
    zones = [r["zone_dropped"] for r in lozo]
    rmse = [r["pred_rmse"][0] for r in lozo]
    lo_e = [r["pred_rmse"][0] - r["pred_rmse"][1] for r in lozo]
    hi_e = [r["pred_rmse"][2] - r["pred_rmse"][0] for r in lozo]
    x = np.arange(3)
    ax.bar(x, rmse, width=0.62, color=[ZONE[z] for z in zones], zorder=2)
    ax.errorbar(x, rmse, yerr=[lo_e, hi_e], fmt="none", ecolor=INK2, elinewidth=0.9,
                capsize=2.5, capthick=0.9, zorder=3)
    # NOT the truth dash of the other panels: a different style, because a dashed
    # black rule already means "truth" in this figure
    ax.axhline(loo["noise_floor"], color=INK2, linestyle=(0, (1, 1.8)), linewidth=1.1,
               zorder=4)
    ax.text(0.5, 0.985, "observation noise $\\sigma$ = 0.1", transform=ax.transAxes,
            ha="center", va="top", fontsize=6.5, color=INK2)
    ax.set_xticks(x)
    ax.set_xticklabels(zones)
    ax.set_ylim(0, 0.145)
    ax.set_title("held-out prediction", loc="left", pad=3)
    ax.set_ylabel("RMSE (mg L$^{-1}$)")
    tidy(ax)

    handles = [
        Patch(facecolor=PRIOR_FILL, edgecolor=PRIOR_EDGE, label="Prior support"),
        Line2D([], [], color=PRIOR_EDGE, linestyle=(0, (2, 2)), label="Prior midpoint"),
        Line2D([], [], color=INK, linestyle=(0, (3, 2)), label="Truth"),
        Line2D([], [], color=INK3, linewidth=3.0, label="Posterior mean, 90% width"),
        Line2D([], [], color=INK2, linestyle=(0, (1, 1.8)), label="Noise floor"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=5, bbox_to_anchor=(0.5, 0.085),
               columnspacing=1.6, handlelength=1.8)
    fig.text(0.5, 0.955, "zone left unobserved", ha="center", fontsize=8.5, color=INK)
    # Not "returns to its prior" without qualification: old and average retain ~100% of the
    # prior SD, but new retains 57%, so the third panel contradicts a universal statement.
    note(fig, "An unobserved zone's coefficient loses most or all of its constraint, fully for "
               "old and average and partly for new; prediction at the same dropped monitors "
               "stays at the noise floor throughout.")
    return finish(fig, "fig6_prediction_without_identification")


# ==================== Figure 5 — symmetric versus structured heterogeneity ====

# One standardising denominator for every "shift in posterior SD" reported in Figures
# 4 and 5: the baseline formal-censored posterior SD, median over the 30 noise
# realisations of Step 6. Mixing denominators across cells is what makes such a matrix
# unreadable, so it is fixed here once.
def base_sd():
    r = [x for x in load("step6_noise_sensitivity.json")["rows"] if x["sigma"] == 0.1][0]
    f = r["by_scheme"]["formal_censored"]
    return {"old": f["old"]["sd_med"], "avg": f["average"]["sd_med"],
            "new": f["new"]["sd_med"]}


def fig5():
    sd = base_sd()
    jit = load("step5c_jitter_sweep.json")
    st = load("step5d_structured.json")
    j20 = [r for r in jit["rows"] if r["jitter"] == 0.2][0]
    zmap = [("old", "old"), ("avg", "average"), ("new", "new")]

    fig, axes = plt.subplots(1, 4, figsize=(W_FULL, 2.70),
                             gridspec_kw={"width_ratios": [1.0, 1.0, 1.3, 1.15],
                                          "wspace": 0.46})
    fig.subplots_adjust(left=0.085, right=0.985, top=0.78, bottom=0.325)
    y = np.arange(3)[::-1]

    # ---- A: 25 symmetric fields, structural increment ----
    ax = axes[0]
    for i, (k, zn) in enumerate(zmap):
        fe = j20["field_ensemble"][zn]
        m, s = fe["increment_mean"] / sd[k], fe["increment_sd"] / sd[k]
        p5, p95 = [v / sd[k] for v in fe["increment_5_95"]]
        ax.plot([p5, p95], [y[i]] * 2, color=ZONE[zn], linewidth=1.2, alpha=0.55, zorder=2)
        ax.plot([m - s, m + s], [y[i]] * 2, color=ZONE[zn], linewidth=3.0,
                solid_capstyle="butt", zorder=3)
        ax.plot([m], [y[i]], marker="o", color=ZONE[zn], markersize=5,
                markeredgecolor="white", markeredgewidth=0.8, zorder=4)
    ax.axvline(0, color=INK, linewidth=0.9, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels([z for _, z in zmap])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(-0.6, 2.6)
    ax.set_xlim(-1.4, 1.4)
    ax.set_title("A  25 symmetric fields", loc="left", pad=4)
    ax.set_xlabel("structural increment (SD)")
    ax.grid(True, axis="x", zorder=0)
    ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)

    # ---- B: one field can point the wrong way ----
    ax = axes[1]
    for i, (k, zn) in enumerate(zmap):
        one = j20["zones"][zn]["bias_increment_vs_control"] / sd[k]
        many = j20["field_ensemble"][zn]["increment_mean"] / sd[k]
        ax.plot([one, many], [y[i]] * 2, color=ZONE[zn], linewidth=0.9, alpha=0.6, zorder=2)
        ax.plot([one], [y[i]], marker="X", color=ZONE[zn], markersize=6, zorder=3)
        ax.plot([many], [y[i]], marker="o", color=ZONE[zn], markersize=5,
                markeredgecolor="white", markeredgewidth=0.8, zorder=4)
    ax.axvline(0, color=INK, linewidth=0.9, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(-0.6, 2.6)
    ax.set_title("B  one field vs 25", loc="left", pad=4)
    ax.set_xlabel("structural increment (SD)")
    ax.grid(True, axis="x", zorder=0)
    ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)

    # ---- C: structured case, where the estimate lands ----
    ax = axes[2]
    # The reported quantity is the 30-realisation median of the RAW fraction f_j, exactly as
    # defined in 2.6.3 and quoted in the text and the Table 4 note (0.86/1.12/0.92). Do not
    # subtract the homogeneous-control offset here: that yields a different estimand, and
    # plotting it against text and table that quote the raw one made the figure contradict
    # the paper (it also flipped which zone sits past the proxy). The faint open marker is
    # the single paired reference realisation, which IS control-adjusted; it is kept only to
    # show that the spread across draws is real, and the caption says so.
    dose = [r for r in st["correlation_dose_response"]["rows"] if r["corr"] == 0.50][0]
    for i, (k, zn) in enumerate(zmap):
        s = st["zones"][zn]["by_scheme"]["formal_censored"]
        m = dose["by_scheme"]["formal_censored"][zn]
        med = m["shift_frac_med"]
        lo, hi = m["shift_frac_5_95"]
        ax.plot([lo, hi], [y[i]] * 2, color=ZONE[zn], linewidth=2.6,
                solid_capstyle="butt", alpha=0.85, zorder=3)
        ax.plot([med], [y[i]], marker="s", color=ZONE[zn], markersize=6.2,
                markeredgecolor="white", markeredgewidth=0.9, zorder=5)
        ax.plot([s["shift_frac_net_of_baseline"]], [y[i] + 0.26], marker="o",
                markerfacecolor="white", markeredgecolor=ZONE[zn], markeredgewidth=0.9,
                markersize=3.8, zorder=4)
        # Beside the bar's left end, not above the square: the single-draw circle sits at +0.26
        # and the raw median can land directly under it. Going sideways also keeps this panel's
        # row spacing identical to A and B, which carry the shared zone labels.
        ax.annotate(f"{med * 100:.0f}%", xy=(lo, y[i]), xytext=(-5, 0),
                    textcoords="offset points", ha="right", va="center",
                    fontsize=6.3, color=ZONE[zn], fontweight="bold")
    ax.axvline(0, color=INK, linewidth=0.9, zorder=1)
    ax.axvline(1, color=INK, linewidth=0.9, linestyle=(0, (3, 2)), zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.tick_params(axis="y", length=0)
    # left margin holds the percentage labels; the right end is set by the raw average-zone
    # 95th percentile, which reaches 2.43
    ax.set_xlim(-0.85, 2.60)
    ax.set_ylim(-1.15, 2.6)
    ax.set_title("C  length-correlated, 30 draws", loc="left", pad=4)
    ax.set_xlabel("position between the two targets")
    ax.set_xticks([0, 1.0, 2.0])
    ax.grid(True, axis="x", zorder=0)
    ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)
    # staggered, because at this x-range the two reference lines are close enough that
    # centred captions on one line would collide
    ax.text(0.0, -0.50, "arithmetic mean", ha="center", va="center",
            fontsize=6.2, color=INK2)
    ax.text(1.0, -0.95, "length-weighted\nproxy (illustrative)", ha="center", va="center",
            fontsize=6.2, color=INK2)

    # ---- D: residual buys no protection ----
    ax = axes[3]
    for r in jit["rows"]:
        worst = max(abs(r["zones"][zn]["bias_increment_vs_control"]) / sd[k]
                    for k, zn in zmap)
        ax.plot([r["struct_residual"]], [worst], marker="o", color=INK3,
                markersize=5, markeredgecolor="white", markeredgewidth=0.8, zorder=3)
        off = {0.0: (6, -1), 0.2: (-5, -3), 0.35: (0, 7), 0.5: (6, -1)}[r["jitter"]]
        ax.annotate(f"±{int(r['jitter'] * 100)}%" if r["jitter"] else "control",
                    xy=(r["struct_residual"], worst), xytext=off,
                    textcoords="offset points", fontsize=6.2, color=INK2,
                    ha="right" if r["jitter"] == 0.2 else
                       ("center" if r["jitter"] == 0.35 else "left"))
    worst_st = max(abs(st["zones"][zn]["by_scheme"]["formal_censored"]
                       ["bias_net_of_baseline"]) / sd[k] for k, zn in zmap)
    ax.plot([st["structural_residual"]], [worst_st], marker="D", color=ZONE["old"],
            markersize=6, markeredgecolor="white", markeredgewidth=0.8, zorder=4)
    ax.annotate("length-\ncorrelated", xy=(st["structural_residual"], worst_st),
                xytext=(4, -2), textcoords="offset points", va="top",
                fontsize=6.2, color=ZONE["old"], fontweight="bold")
    ax.set_xlim(0, 0.030)
    ax.set_ylim(-0.12, 2.35)
    ax.set_xticks([0, 0.01, 0.02, 0.03])
    ax.set_title("D  residual vs bias", loc="left", pad=4)
    ax.set_xlabel("minimum noise-free grid RMSE (mg L$^{-1}$)")
    ax.set_ylabel("largest |shift| (SD)")
    tidy(ax, "both")

    handles = [
        Line2D([], [], marker="X", color=INK3, ls="none", label="single field (B)"),
        Line2D([], [], marker="o", color=INK3, ls="none", label="mean of 25 fields (A, B)"),
        Line2D([], [], marker="o", markerfacecolor="white", markeredgecolor=INK3,
               ls="none", label="single noise draw (C)"),
        Line2D([], [], marker="s", color=INK3, ls="none",
               label="median of 30 draws, 5–95% (C, reported)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, bbox_to_anchor=(0.5, 0.10),
               columnspacing=1.5, handletextpad=0.3)
    note(fig, "Symmetric mean-zero heterogeneity leaves no bias resolvable against the "
               "field-to-field scatter. A length-correlated field displaces all three "
               "coefficients at the same residual, by a median fraction near 1; the intervals "
               "overlap, so the zones are not ordered.")
    return finish(fig, "fig5_symmetric_vs_structured")


# ==================== Figure 4 — standardised-effect matrix ====================

def fig4():
    sd = base_sd()
    s6 = {r["sigma"]: r["by_scheme"]["formal_censored"] for r in
          load("step6_noise_sensitivity.json")["rows"]}
    ar1 = load("step7c_ar1.json")["coef"]
    kb = {r["kb"]: r["by_scheme"]["formal_censored"] for r in
          load("step8b_kb_sensitivity.json")["rows"]}
    bias = [r for r in load("step8c_bias_bynode.json")["rows"]
            if r["scheme"] == "formal_censored"]
    zc = load("step9_zeroclip.json")["coef"]
    j20 = [r for r in load("step5c_jitter_sweep.json")["rows"] if r["jitter"] == 0.2][0]
    st = load("step5d_structured.json")["zones"]

    zmap = [("old", "old"), ("avg", "average"), ("new", "new")]
    bias_key = {"old": "d_old_over_sd", "avg": "d_avg_over_sd", "new": "d_new_over_sd"}

    wide_cols = ["$\\sigma$ 0.10 to 0.15", "AR(1) $\\rho$ = 0.4"]
    W = np.array([[s6[0.15][zn]["sd_ret_med"] / s6[0.1][zn]["sd_ret_med"],
                   ar1[zn]["widening"]] for _, zn in zmap])

    # The sensor-bias column takes the maximum SEPARATELY for each coefficient, so its three
    # cells generally come from three different arms (here node 15 +0.10, node 231 -0.10 and
    # node 113 -0.10). Labelling it "worst arm" implied one offset produced all three at once.
    shift_cols = ["$k_b$ = −0.4", "$k_b$ = −0.6", "sensor bias\n(max over arms,\nper coefficient)",
                  "zero\ncensoring", "symmetric\nheterogeneity", "structured\nheterogeneity"]
    S = np.array([[
        kb[-0.4]["shift_over_own_sd"][k],
        kb[-0.6]["shift_over_own_sd"][k],
        max((r[bias_key[k]] for r in bias), key=abs),
        zc[zn]["delta_median"] / sd[k],
        j20["field_ensemble"][zn]["increment_mean"] / sd[k],
        st[zn]["by_scheme"]["formal_censored"]["bias_net_of_baseline"] / sd[k],
    ] for k, zn in zmap])

    # tall enough for a two-line tick label under the matrix AND a wrapped note under that;
    # at 2.05 in the note ran into the labels
    fig, axes = plt.subplots(1, 2, figsize=(W_FULL, 2.55),
                             gridspec_kw={"width_ratios": [2.0, 6.0], "wspace": 0.09})
    fig.subplots_adjust(left=0.115, right=0.985, top=0.80, bottom=0.30)

    for ax, M, cols, cmap, vmax, title, vmin in (
            (axes[0], W, wide_cols, "Blues", 2.0, "A  interval widening (×)", 1.0),
            (axes[1], np.abs(S), shift_cols, "Oranges", 6.0,
             "B  displacement (baseline posterior SD)", 0.0)):
        ax.imshow(M, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels(cols, fontsize=6.6)
        ax.set_yticks(range(3))
        ax.set_yticklabels([z for _, z in zmap] if ax is axes[0] else [], fontsize=7)
        ax.tick_params(length=0)
        ax.set_title(title, loc="left", pad=5)
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_xticks(np.arange(-0.5, len(cols), 1), minor=True)
        ax.set_yticks(np.arange(-0.5, 3, 1), minor=True)
        ax.grid(which="minor", color="white", linewidth=1.6)
        ax.tick_params(which="minor", length=0)

    for i in range(3):
        for j in range(2):
            axes[0].text(j, i, f"{W[i, j]:.2f}", ha="center", va="center", fontsize=6.8,
                         color="white" if (W[i, j] - 1.0) > 0.55 else INK)
        for j in range(6):
            v = S[i, j]
            axes[1].text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=6.8,
                         color="white" if abs(v) / 6.0 > 0.55 else INK)

    note(fig, "Sign is printed, not encoded by colour. Column 3 maximises over the 24 bias arms "
               "separately for each coefficient, so its three cells need not share an arm. Drift "
               "is omitted: two monitors only, and it reproduces 0.89–0.99 of the mean-equivalent "
               "bias already in column 3. The structured column is a single noise draw (Fig. 5C).",
         y=0.015)
    return finish(fig, "fig4_standardised_effects")


# ==================== Figure 7 — risk robustness is metric-specific ==========

def per_node_risk():
    """Weighted E[D] and E[A] for all 92 junctions, from the cached prediction cube.

    Reproduces Step 10 exactly (checked against its stored top-10 and both network
    means), so the 92-node field used here and the headline risk artifact are the same
    quantity computed the same way.
    """
    d = npz()
    C = d["C_all"].astype(np.float64)
    w = weights_from_loglik(d["loglik_censored"])
    dur = np.trapezoid((C < 0.2).astype(np.float64), dx=1.0, axis=1)
    dfc = np.trapezoid(np.clip(0.2 - C, 0, None), dx=1.0, axis=1)
    return d["all_nodes"], w @ dur, w @ dfc


def fig7():
    nodes, Ed, Ea = per_node_risk()
    kb = {r["kb"]: r["by_scheme"]["formal_censored"] for r in
          load("step8b_kb_sensitivity.json")["rows"]}
    bias = [r for r in load("step8c_bias_bynode.json")["rows"]
            if r["scheme"] == "formal_censored"]

    fig, axes = plt.subplots(1, 3, figsize=(W_FULL, 2.75),
                             gridspec_kw={"width_ratios": [0.95, 1.25, 1.3],
                                          "wspace": 0.40})
    fig.subplots_adjust(left=0.055, right=0.985, top=0.80, bottom=0.30)

    # ---- A: the two metrics do not agree on the leading nodes ----
    ax = axes[0]
    rd = {n: i + 1 for i, n in enumerate(nodes[np.argsort(-Ed)])}
    ra = {n: i + 1 for i, n in enumerate(nodes[np.argsort(-Ea)])}
    shown = sorted({n for n in nodes if rd[n] <= 8 or ra[n] <= 8}, key=lambda n: rd[n])
    for n in shown:
        moved = rd[n] != ra[n]
        ax.plot([0, 1], [rd[n], ra[n]], color=ZONE["old"] if moved else INK3,
                linewidth=1.4 if moved else 0.8, alpha=0.95 if moved else 0.5, zorder=2)
        ax.plot([0, 1], [rd[n], ra[n]], ls="none", marker="o", markersize=3.6,
                color=ZONE["old"] if moved else INK3, zorder=3)
        ax.annotate(n, xy=(0, rd[n]), xytext=(-4, 0), textcoords="offset points",
                    ha="right", va="center", fontsize=6.2, color=INK2)
        ax.annotate(n, xy=(1, ra[n]), xytext=(4, 0), textcoords="offset points",
                    ha="left", va="center", fontsize=6.2, color=INK2)
    ax.set_xlim(-0.42, 1.42)
    ax.set_ylim(max(max(rd[n] for n in shown), max(ra[n] for n in shown)) + 0.6, 0.4)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["by $\\bar{P}$\n(duration)", "by E[A]\n(deficit)"], fontsize=6.8)
    ax.set_ylabel("rank")
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    ax.set_title("A  which nodes lead", loc="left", pad=4)
    ax.tick_params(axis="x", length=0)
    ax.grid(True, axis="y", zorder=0)
    ax.set_axisbelow(True)

    # ---- B: shortlist stability depends on the metric ----
    ax = axes[1]
    worst_bias = min(bias, key=lambda r: (r["risk_top6_jaccard_vs_unbiased"] +
                                          r["deficit_top6_jaccard_vs_unbiased"]))
    cases = [
        ("$k_b$ = −0.4", kb[-0.4]["risk_top6_jaccard_vs_kb_ref"],
         kb[-0.4]["deficit_top6_jaccard_vs_kb_ref"],
         kb[-0.4]["risk_spearman_vs_kb_ref"], kb[-0.4]["deficit_spearman_vs_kb_ref"]),
        ("$k_b$ = −0.6", kb[-0.6]["risk_top6_jaccard_vs_kb_ref"],
         kb[-0.6]["deficit_top6_jaccard_vs_kb_ref"],
         kb[-0.6]["risk_spearman_vs_kb_ref"], kb[-0.6]["deficit_spearman_vs_kb_ref"]),
        (f"sensor bias\nnode {worst_bias['node']} {worst_bias['offset']:+.2f}",
         worst_bias["risk_top6_jaccard_vs_unbiased"],
         worst_bias["deficit_top6_jaccard_vs_unbiased"],
         worst_bias["risk_spearman_vs_unbiased"],
         worst_bias["deficit_spearman_vs_unbiased"]),
    ]
    x = np.arange(len(cases))
    ax.bar(x - 0.19, [c[1] for c in cases], width=0.34, color=SCHEME["informal_glue"],
           zorder=2, label="$\\bar{P}$ (duration)")
    ax.bar(x + 0.19, [c[2] for c in cases], width=0.34, color=SCHEME["formal_censored"],
           zorder=2, label="E[A] (deficit)")
    for xi, c in zip(x, cases):
        r = min(c[3], c[4])
        ax.text(xi, 1.15, "$\\rho_s$ > 0.999" if r > 0.999 else f"$\\rho_s$ {r:.3f}",
                ha="center", fontsize=6.2, color=INK2)
    ax.axhline(1.0, color=INK3, linewidth=0.7, zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels([c[0] for c in cases], fontsize=6.8)
    ax.set_ylim(0, 1.30)
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_ylabel("top-6 Jaccard vs reference")
    ax.set_title("B  shortlist stability", loc="left", pad=4)
    ax.tick_params(axis="x", length=0)
    tidy(ax)
    ax.legend(loc="lower left", bbox_to_anchor=(0.03, 0.02), fontsize=6.3,
              handlelength=1.2, handletextpad=0.4, borderpad=0.2, labelspacing=0.25)

    # ---- C: duration and depth are different severities ----
    ax = axes[2]
    ax.scatter(Ed, Ea, s=13, color=INK3, alpha=0.55, linewidths=0, zorder=2)
    for n, dx, dy, ha in (("145", 5, 0, "left"), ("129", -5, 3, "right"),
                          ("131", -6, 0, "right")):
        i = int(np.where(nodes == n)[0][0])
        ax.scatter([Ed[i]], [Ea[i]], s=30, color=ZONE["old"], zorder=4,
                   edgecolors="white", linewidths=0.8)
        ax.annotate(f"node {n}", xy=(Ed[i], Ea[i]), xytext=(dx, dy),
                    textcoords="offset points", ha=ha, va="center", fontsize=6.4,
                    color=ZONE["old"], fontweight="bold")
    ax.set_xlabel("E[D], hours below 0.2 mg L$^{-1}$")
    ax.set_ylabel("E[A], deficit (mg L$^{-1}$ h)")
    ax.set_title("C  duration is not depth", loc="left", pad=4)
    n0 = int((Ed <= 1e-9).sum())
    ax.text(0.03, 0.80, f"{n0} of 92 junctions\nnever go below 0.2", transform=ax.transAxes,
            va="top", fontsize=6.3, color=INK2)
    tidy(ax, "both")

    note(fig, "Network-wide rank correlation stays high in every case; the six-node shortlist "
               "does not, and which nodes it keeps depends on the metric.", y=0.02, size=6.5)
    return finish(fig, "fig7_metric_specific_risk")


# ==================== Figure 1 — controlled study design ====================

def fig1():
    """Study system and forward baseline.

    Panel B replaces an earlier workflow schematic: the section structure already carries the
    evidence chain, whereas the baseline concentration field is a result the reader needs before
    any residual, censored observation or risk metric can be interpreted.
    """
    import wntr
    inp = os.path.join(HERE, os.pardir, "models", "net3_frozen", "Net3.inp")
    wn = wntr.network.WaterNetworkModel(inp)
    sys.path.insert(0, HERE)
    from wq_common import MATERIAL_ZONES, MONITOR_NODES

    d = npz()
    truth = d["truth_all"].astype(np.float64)
    noisy = d["noisy"].astype(np.float64)
    mon_pos = list(d["mon_pos"])
    warm = 120
    t = np.arange(warm, truth.shape[0])
    zone_of_mon = {"107": "new", "113": "new", "15": "old", "145": "old",
                   "209": "average", "231": "average"}

    fig = plt.figure(figsize=(W_FULL, 2.75))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.35, 0.85], wspace=0.30)

    # ---- A: network, zones, monitors ----
    ax = fig.add_subplot(gs[0, 0])
    xy = {n: wn.get_node(n).coordinates for n in wn.node_name_list}
    for pipe in wn.pipe_name_list:
        lk = wn.get_link(pipe)
        a, b = xy[lk.start_node_name], xy[lk.end_node_name]
        ax.plot([a[0], b[0]], [a[1], b[1]], color=ZONE[MATERIAL_ZONES[pipe]],
                linewidth=0.9, alpha=0.85, solid_capstyle="round", zorder=2)
    jx = np.array([xy[n] for n in wn.junction_name_list])
    ax.scatter(jx[:, 0], jx[:, 1], s=2.5, color=INK3, alpha=0.6, linewidths=0, zorder=3)
    for n in wn.tank_name_list:
        ax.scatter(*xy[n], s=34, marker="s", facecolor="white", edgecolor=INK,
                   linewidths=0.9, zorder=5)
    for n in wn.reservoir_name_list:
        ax.scatter(*xy[n], s=44, marker="^", facecolor="white", edgecolor=INK,
                   linewidths=0.9, zorder=5)
    mx = np.array([xy[n] for n in MONITOR_NODES])
    ax.scatter(mx[:, 0], mx[:, 1], s=26, marker="o", facecolor="white",
               edgecolor=INK, linewidths=1.1, zorder=6)
    for n, (px, py) in zip(MONITOR_NODES, mx):
        ax.annotate(n, xy=(px, py), xytext=(0, 6), textcoords="offset points",
                    ha="center", fontsize=6.0, color=INK, fontweight="bold", zorder=7)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("A  network and zones", loc="left", pad=2)
    handles = [Line2D([], [], color=ZONE[z], linewidth=2.0, label=z) for z in
               ("old", "average", "new")]
    handles += [Line2D([], [], marker="o", markerfacecolor="white", markeredgecolor=INK,
                       ls="none", markersize=4.5, label="monitor"),
                Line2D([], [], marker="^", markerfacecolor="white", markeredgecolor=INK,
                       ls="none", markersize=5, label="source"),
                Line2D([], [], marker="s", markerfacecolor="white", markeredgecolor=INK,
                       ls="none", markersize=4.5, label="tank")]
    ax.legend(handles=handles, loc="lower left", ncol=2, fontsize=5.8, handlelength=1.2,
              handletextpad=0.4, columnspacing=0.8, labelspacing=0.3, borderpad=0.2)

    # ---- B: monitored series over the assessment window ----
    ax = fig.add_subplot(gs[0, 1])
    for i, n in enumerate(MONITOR_NODES):
        c = ZONE[zone_of_mon[n]]
        ls = "-" if i % 2 == 0 else (0, (4, 1.5))
        ax.plot(t, truth[warm:, mon_pos[i]], color=c, linewidth=1.1, linestyle=ls, zorder=3)
        ax.scatter(t, noisy[warm:, i], s=3.5, color=c, alpha=0.40, linewidths=0, zorder=2)
        ax.annotate(n, xy=(t[-1], truth[-1, mon_pos[i]]), xytext=(3, 0),
                    textcoords="offset points", va="center", fontsize=5.8, color=c)
    ax.axhline(0.2, color=INK, linestyle=(0, (1, 1.8)), linewidth=1.0, zorder=4)
    ax.annotate("0.2 mg L$^{-1}$", xy=(t[0], 0.2), xytext=(2, 3),
                textcoords="offset points", fontsize=6.0, color=INK2)
    ax.set_xlim(warm, t[-1] + 6)
    ax.set_ylim(-0.03, 1.02)
    ax.set_xticks([120, 132, 144, 156, 168])
    ax.set_xlabel("time (h)")
    ax.set_ylabel("chlorine (mg L$^{-1}$)")
    ax.set_title("B  monitored series, assessment window", loc="left", pad=2)
    tidy(ax)

    # ---- C: where the network sits relative to the threshold ----
    ax = fig.add_subplot(gs[0, 2])
    wmin = truth[warm:].min(axis=0)
    below = wmin < 0.2
    ax.hist([wmin[below], wmin[~below]], bins=np.linspace(0, 1.0, 21), stacked=True,
            color=[ZONE["old"], INK3], edgecolor="white", linewidth=0.4, zorder=2)
    ax.axvline(0.2, color=INK, linestyle=(0, (1, 1.8)), linewidth=1.0, zorder=3)
    ax.set_xlabel("window minimum (mg L$^{-1}$)")
    ax.set_ylabel("junctions")
    ax.set_title("C  network minima", loc="left", pad=2)
    ax.annotate(f"{int(below.sum())} of {len(wmin)}\nbelow 0.2", xy=(0.22, 0.94),
                xycoords="axes fraction", fontsize=6.2, color=ZONE["old"],
                fontweight="bold", va="top")
    tidy(ax)

    fig.subplots_adjust(left=0.045, right=0.985, top=0.90, bottom=0.19)

    # Provenance for every forward-baseline number quoted in Section 3.1. These are properties
    # of the deterministic noise-free truth, NOT of the weighted ensemble: the count below is
    # "junctions whose window minimum drops under C_crit at least once", which is a different
    # quantity from Step 12's "junctions with P_min > 0.5" even though both happen to be 21.
    summ = {
        "source": "baseline.npz truth_all, deterministic noise-free reference trajectory",
        "window_h": [warm, int(truth.shape[0] - 1)],
        "C_crit_mgL": 0.2,
        "network_median_mgL": float(np.median(truth[warm:])),
        "network_p5_p95_mgL": [float(np.percentile(truth[warm:], 5)),
                               float(np.percentile(truth[warm:], 95))],
        "monitor_window_mean_mgL": {n: float(truth[warm:, mon_pos[i]].mean())
                                    for i, n in enumerate(MONITOR_NODES)},
        "n_junctions": int(truth.shape[1]),
        "n_junctions_min_below_C_crit": int(below.sum()),
        "definition_of_count": "#{n : min_t C_true(n,t) < C_crit over the window}",
        "not_to_be_confused_with": "step12 scenario_summary P_min_gt_0.5_nodes (ensemble-weighted)",
    }
    with open(os.path.join(CACHE, "forward_baseline.json"), "w") as fh:
        json.dump(summ, fh, indent=2)
    print("wrote baseline_cache/forward_baseline.json")
    return finish(fig, "fig1_study_design")


# ==================== Tables 1, 3–5 ====================

def tables():
    meta = load("baseline_meta.json")["summary"]
    cfg = load("baseline_meta.json")["config"]
    sd = base_sd()
    loo = load("step11_loo.json")
    _kb_rows = load("step8b_kb_sensitivity.json")["rows"]
    kb = {r["kb"]: r["by_scheme"]["formal_censored"] for r in _kb_rows}
    # the aggregate-fit statistic is a property of the candidate library at that k_b, not of a
    # weighting scheme, so it sits at row level rather than inside by_scheme
    kb_fit = {r["kb"]: r["rmse_min_over_noise_med"] for r in _kb_rows}
    bias = [r for r in load("step8c_bias_bynode.json")["rows"]
            if r["scheme"] == "formal_censored"]
    s6 = {r["sigma"]: r["by_scheme"]["formal_censored"] for r in
          load("step6_noise_sensitivity.json")["rows"]}
    ar1 = load("step7c_ar1.json")["coef"]
    zc = load("step9_zeroclip.json")
    j20 = [r for r in load("step5c_jitter_sweep.json")["rows"] if r["jitter"] == 0.2][0]
    _st_full = load("step5d_structured.json")
    st = _st_full["zones"]
    # Best achievable aggregate fit under the structured truth, against the baseline realisation's
    # noise RMSE. Same quantity the sensor-bias and k_b rows now carry, so the column is comparable.
    st_rmse = _st_full["rmse_min"] / meta["noise_rmse"]
    sb = load("step8_sensor_bias.json")
    zmap = [("old", "old"), ("avg", "average"), ("new", "new")]
    bkey = {"old": "d_old_over_sd", "avg": "d_avg_over_sd", "new": "d_new_over_sd"}
    L = []

    cfgm = load("cache_manifest.json")["config"]
    zonecfg = {"old": ("old", 14, "15, 145"), "avg": ("average", 37, "209, 231"),
          "new": ("new", 66, "107, 113")}
    L += ["## Table 1 — Baseline design and parameter ranges", "",
          "| Quantity | Old zone | Average zone | New zone | Unit |",
          "|---|---:|---:|---:|---|",
          "| Pipes | " + " | ".join(str(zonecfg[k][1]) for k in zonecfg) + " | — |",
          "| True wall decay | " + " | ".join(f"{TRUTH[k]:+.3f}" for k in zonecfg) + " | m/day |",
          "| Prior range | " + " | ".join(f"[{PRIOR_BOX[k][0]:+.3f}, {PRIOR_BOX[k][1]:+.3f}]"
                                          for k in zonecfg) + " | m/day |",
          "| Prior SD | " + " | ".join(f"{meta['prior_sd'][k]:.4f}" for k in zonecfg) + " | m/day |",
          "| Monitors | " + " | ".join(zonecfg[k][2] for k in zonecfg) + " | junction ID |",
          "",
          "| Shared setting | Value |", "|---|---|",
          f"| Bulk decay (first order, fixed) | {cfg['kb_fixed']:+.1f} 1/day |",
          f"| Source chlorine | {cfg['inlet_mgl']:.1f} mg/L |",
          f"| Tank initial chlorine | {cfgm['tank_init_mgl']:.1f} mg/L |",
          f"| Hydraulic and reporting timestep | {cfg['hydraulic_timestep_s'] // 3600} h |",
          f"| Water-quality timestep | {cfg['quality_timestep_s'] // 60} min |",
          f"| Total duration / warm-up / assessment window | {cfg['duration_h']} / "
          f"{cfg['warmup_h']} / {cfg['duration_h'] - cfg['warmup_h']} h |",
          f"| Observation error SD | {cfg['sigma_obs']:.2f} mg/L |",
          "| Operational threshold | 0.20 mg/L |",
          f"| Candidate library | {cfg['n_mc']} scrambled Sobol draws |",
          f"| Observations | {meta['n_resid']} (6 monitors x {meta['Tn']} hourly values) |",
          ""]

    L += ["## Table 3 — Baseline inference comparison", "",
          f"Same {meta['n_resid']} residuals, same {meta['n_mc']} Sobol candidates; only the "
          "weighting rule changes.", "",
          "| Coefficient | Truth | Prior SD | Scheme | Mean | SD | SD retained | ESS |",
          "|---|---:|---:|---|---:|---:|---:|---:|"]
    names = {"formal_censored": "formal censored", "formal_iid": "formal iid",
             "informal_glue": "informal GLUE 0.107",
             "informal_glue_draft_thr": "informal GLUE 0.120"}
    for k, zn in zmap:
        for i, sch in enumerate(names):
            c = meta["schemes"][sch]["coef"][k]
            ess = meta["schemes"][sch]["diagnostics"]["ess"]
            L.append(f"| {zn if i == 0 else ''} | {TRUTH[k]:+.3f} | "
                     f"{meta['prior_sd'][k]:.4f} | {names[sch]} | {c['mean']:+.4f} | "
                     f"{c['sd']:.4f} | {c['sd_retained'] * 100:.1f}% | {ess:.0f} |")
    L += ["",
          f"ESS is out of {meta['n_mc']}: the formal weights concentrate on about "
          f"{meta['schemes']['formal_censored']['diagnostics']['ess']:.0f} effective draws "
          f"({meta['schemes']['formal_censored']['diagnostics']['ess_frac'] * 100:.1f}%), the "
          f"informal score on "
          f"{meta['schemes']['informal_glue']['diagnostics']['ess']:.0f} "
          f"({meta['schemes']['informal_glue']['diagnostics']['ess_frac'] * 100:.1f}%).", ""]

    dose50 = [r for r in load("step5d_structured.json")["correlation_dose_response"]["rows"]
              if r["corr"] == 0.50][0]["by_scheme"]["formal_censored"]
    med = "/".join(f"{dose50[zn]['shift_frac_med']:.2f}" for _, zn in zmap)
    wide = [f"x{s6[0.15][zn]['sd_ret_med'] / s6[0.1][zn]['sd_ret_med']:.2f}" for _, zn in zmap]
    worst = {k: max((r[bkey[k]] for r in bias), key=abs) for k, _ in zmap}
    nb = len(bias)
    ndef = sum(1 for r in bias if r["deficit_top6_jaccard_vs_unbiased"] == 1.0)
    rho_min = min(r["risk_spearman_vs_unbiased"] for r in bias)

    L += ["## Table 4 — Robustness summary", "",
          "| Error source | old | average | new | Best fit | Risk |",
          "|---|---:|---:|---:|---:|---|"]
    L.append(f"| sigma 0.10 to 0.15 (W) | {wide[0]} | {wide[1]} | {wide[2]} | — | — |")
    L.append("| AR(1) rho = 0.4 (W) | " + " | ".join(f"x{ar1[zn]['widening']:.2f}" for _, zn in zmap)
             + " | — | — |")
    for lab, kbv, row in (("k_b = -0.4", -0.4, kb[-0.4]), ("k_b = -0.6", -0.6, kb[-0.6])):
        rk = (f"Spearman {row['risk_spearman_vs_kb_ref']:.3f}; top-6 "
              f"{row['risk_top6_jaccard_vs_kb_ref']:.2f} / "
              f"{row['deficit_top6_jaccard_vs_kb_ref']:.2f}")
        L.append(f"| {lab} (D) | " + " | ".join(f"{row['shift_over_own_sd'][k]:+.2f}"
                 for k, _ in zmap)
                 + f" | {kb_fit[kbv]:.3f} | {rk} |")
    # from `bias`, not from step8: this row's displacements are the maximum over the 24 arms of
    # the location sweep, so the fit ratio has to be measured over those same arms. Step 8 covers
    # node 15 only and reports a narrower range, which read as if the largest displacement had
    # left the residual almost untouched.
    bias_ratios = [r["rmse_min_over_noise_med"] for r in bias]
    L.append("| sensor bias, max over arms (D) | " + " | ".join(f"{worst[k]:+.2f}" for k, _ in zmap)
             + f" | {min(bias_ratios):.3f}-{max(bias_ratios):.3f} "
             # a lower bound has to be rounded DOWN. Printing this at four decimals rounded
             # 0.99928 up to 0.9993 and the table then asserted a floor the arms do not meet.
             + f"| Spearman >= {math.floor(rho_min * 1000) / 1000:.3f}; "
             + f"E[A] top-6 held in {ndef}/{nb} |")
    # "unchanged" was a literal in this cell while the artifact recorded two drift arms that drop
    # a node from the duration shortlist. A hard-coded verdict cannot go stale loudly, so the
    # counts are read from the artifact and the two metrics are reported separately, as they are
    # everywhere else in this table.
    drift = [r for r in load("step8d_sensor_drift.json")["rows"] if r["arm"] == "drift"]
    ndur = sum(1 for r in drift if r["risk_top6_jaccard_vs_unbiased"] >= 1.0)
    ndfc = sum(1 for r in drift if r["deficit_top6_jaccard_vs_unbiased"] >= 1.0)
    L.append("| sensor drift (D) | 0.89-0.99 of the mean-equivalent bias | | | — "
             + f"| E[D] top-6 held in {ndur}/{len(drift)}; E[A] in {ndfc}/{len(drift)} |")
    # Two different bases in one row, so say so: the count censuses the reference seed's
    # observation set, the displacement is a median over 30 seeds.
    L.append("| zero censoring (D) | " + " | ".join(
        f"{zc['coef'][zn]['delta_median'] / sd[k]:+.2f}" for k, zn in zmap)
        + " | — | top-3 unchanged |")
    L.append("| symmetric heterogeneity (D) | " + " | ".join(
        f"{j20['field_ensemble'][zn]['increment_mean'] / sd[k]:+.2f}" for k, zn in zmap)
        + " | — | — |")
    L.append("| structured heterogeneity (D) | " + " | ".join(
        f"{st[zn]['by_scheme']['formal_censored']['bias_net_of_baseline'] / sd[k]:+.2f}"
        for k, zn in zmap)
        + f" | {st_rmse:.3f} | — |")
    # The per-row Basis column was folded into this note to get the table down to six columns.
    # Two of its entries were caveats rather than counts and have to survive the move: the three
    # sensor-bias displacements are maxima taken separately per coefficient, so they come from
    # three different arms and are not one perturbation's vector, and the censoring count is a
    # census of the reference seed while the displacement beside it is a median over thirty.
    # Plain paragraphs, not a one-column table. As a table they repeated a "Notes" header on every
    # page they spilled onto and broke mid-sentence at the page boundary, and the splice could not
    # place them after the caption that refers to them.
    L += ["",
          "Notes. (W) widening, a factor on the interval; (D) displacement, in baseline posterior "
          "standard deviations. Risk gives the whole-network Spearman rank correlation, then the "
          "top-6 Jaccard index under $\\bar{P}$ / $E[A]$. Best fit is the "
          "smallest RMSE reachable anywhere in the candidate library under that perturbation, as "
          "a multiple of the realised observation-noise RMSE of the same draw, so a value near "
          "one means the perturbation is not visible in the aggregate residual.", "",
          f"Basis. Thirty noise realisations unless stated otherwise. AR(1) is a Fisher "
          f"calculation. Sensor bias is {nb} arms of 30, the maximum taken separately for each "
          "coefficient, so its three values come from different arms. Drift is 2 nodes x 4 "
          f"endpoints. Zero censoring is 30 realisations, but the {zc['cal_zero']} floored "
          f"readings of {meta['n_resid']} census the reference seed alone. Symmetric "
          "heterogeneity is 25 fields; structured heterogeneity is one design and a single "
          f"paired draw, whose median proxy-gap fraction over 30 draws is {med} (raw), with "
          "5-95% intervals that overlap across zones.", ""]

    L += ["## Table 5 — Held-out validation", "",
          "| Design | Withheld | Own-coefficient error | Own SD retained | Held-out RMSE | 90% coverage |",
          "|---|---|---:|---:|---:|---:|"]
    for r in loo["rows"]:
        if r["scheme"] != "formal_censored":
            continue
        L.append(f"| leave-one-monitor-out | node {r['node']} ({r['zone']}) | - | - | "
                 f"{r['pred_rmse'][0]:.4f} | {r['coverage90'][0]:.2f} |")
    for r in loo["leave_one_zone_out"]:
        if r["scheme"] != "formal_censored":
            continue
        L.append(f"| leave-one-zone-out | {r['zone_dropped']} zone "
                 f"({', '.join(r['monitors_dropped'])}) | {r['own_coef_error']:+.4f} | "
                 f"{r['own_sd_retained'][0] * 100:.0f}% | {r['pred_rmse'][0]:.4f} | "
                 f"{r['coverage90'][0]:.2f} |")
    u = loo["unmonitored_validation"]["by_scheme"]["formal_censored"]
    L += ["",
          f"At the {loo['unmonitored_validation']['n_junctions']} never-calibrated junctions "
          f"the median normalised mean absolute error, sum|error| / sum|truth| per junction, "
          f"is {u['mean_abs_rel_error'][0] * 100:.2f}% "
          f"(formal censored). Noise floor sigma = {loo['sigma']} mg/L.", ""]

    os.makedirs(OUTDIR, exist_ok=True)
    path = os.path.join(OUTDIR, "tables.md")
    with open(path, "w") as fh:
        fh.write("<!-- generated by paper_figs.py; do not hand-edit -->\n\n" + "\n".join(L))
    print("wrote", os.path.relpath(path, HERE))
    return path


FIGS = {"fig1": fig1, "fig2": fig2, "fig3": fig3, "fig4": fig4, "fig5": fig5,
        "fig6": fig6, "fig7": fig7, "tables": tables}

if __name__ == "__main__":
    want = sys.argv[1:] or ["all"]
    todo = list(FIGS) if want == ["all"] else want
    for name in todo:
        FIGS[name]()
