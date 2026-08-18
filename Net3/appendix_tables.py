"""Generate the appendix tables for the Research Paper from the cached artifacts.

Writes Figures/paper/appendix_tables.md, which 00-thesis/build.sh splices into the manuscript at
`{{A:<label>}}` markers. Labels are letters plus an index (A1, F2, ...) rather than running
numbers, so inserting an appendix does not renumber every table after it.

Nothing here recomputes a result. Every value is read from baseline_cache/, so an appendix table
cannot drift from the artifact the main text quotes.
"""
import csv
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "baseline_cache")
OUT = os.path.join(HERE, "figures", "paper", "appendix_tables.md")


def load(name):
    with open(os.path.join(CACHE, name)) as f:
        return json.load(f)


def table(rows, header):
    """A pipe table. Cells are pre-formatted strings; no rounding happens here by accident."""
    out = ["| " + " | ".join(header) + " |",
           "|" + "|".join(["---"] * len(header)) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(c) for c in r) + " |")
    return "\n".join(out)


blocks = []


def add(label, caption, body):
    blocks.append("## %s\n\n%s\n\n%s\n" % (label, caption, body))


# ---------------------------------------------------------------- A: warm-up and horizon
s0 = load("step0_warmup_convergence.json")
crit_order = ["tank_max_dLevel", "monitor_max_dC", "network_p95_dC", "tank_max_dC",
              "risk_rel_dDeficit", "age_p95_dAge"]
crit_name = {
    "tank_max_dLevel": "Tank level (m)",
    "monitor_max_dC": "Monitor chlorine (mg L-1)",
    "network_p95_dC": "Network 95th-percentile chlorine (mg L-1)",
    "tank_max_dC": "Tank chlorine (mg L-1)",
    "risk_rel_dDeficit": "Risk severity, relative change in E[A]",
    "age_p95_dAge": "Water age, 95th percentile (h)",
}
rows = []
for k in crit_order:
    c = s0["per_criterion"][k]
    w = c["worst_by_cycle"]
    fp = c.get("first_pass_warmup_h")
    if fp is not None:
        verdict = "%d h (verified)" % fp
    else:
        ex = c.get("extrapolated_warmup_h")
        verdict = "not met within 168 h; extrapolates to ~%s h (unverified)" % ex
    rows.append([crit_name[k], c["tolerance"]] +
                ["%.4f" % v if v >= 1e-4 else "%.1e" % v for v in w] + [verdict])
add("A1",
    "Table A1. Cyclostationarity diagnostics. Worst value across the truth and both corners of "
    "the prior box, for each pair of successive 24 h cycles. Tolerances were declared before the "
    "values were seen.",
    table(rows, ["Criterion", "Tolerance", "0-24", "24-48", "48-72", "72-96", "96-120",
                 "120-144", "Earliest warm-up satisfying it"]))

truth = [c["sets"]["truth"] for c in s0["cycles"]]
add("A2",
    "Table A2. Network-mean cumulative deficit for the truth, by cycle. The integrated severity "
    "falls after the start-up transient and then rises monotonically, which is why the 2% "
    "criterion is not met inside the model horizon.",
    table([["%d-%d h" % (24 * i, 24 * i + 24), "%.4f" % t["net_mean_deficit"],
            "%.4f" % t["risk_rel_dDeficit"]] for i, t in enumerate(truth)],
          ["Cycle", "Net-mean deficit (mg L-1 h)", "Relative change"]))

# ---------------------------------------------------------------- B: numerical and unit checks
s13 = load("step13_known_answer.json")
add("B1",
    "Table B1. Bulk-decay arm of the known-answer test against the analytic first-order solution "
    "at the far end of a single pipe. The relative error grows with the amount of decay, which is "
    "the signature of the quality solver's time discretisation; a wrong conversion factor would "
    "instead give a constant relative offset.",
    table([["%.2f" % r["kb_per_day"], "%.4f" % r["t_res_h"], "%.6f" % r["epanet_c"],
            "%.6f" % r["analytic_c"], "%.1e" % r["rel_error"]] for r in s13["bulk_arm"]["rows"]],
          ["k_b (1/day)", "Residence time (h)", "EPANET C", "Analytic C", "Relative error"]))

add("B2",
    "Table B2. Wall-decay arm. The exact concentration cannot be predicted from the wall "
    "coefficient alone because EPANET's first-order wall reaction is also limited by mass "
    "transfer, so what is asserted is the sign, monotonicity, and that the result stays above the "
    "mass-transfer-free limit.",
    table([["%.2f" % r["kw_per_day"], "%.6f" % r["epanet_c"], "%.6f" % r["extra_decay_vs_kw0"],
            "%.6f" % r["kf_free_lower_bound_c"],
            "yes" if r["above_bound"] else ("n/a" if r["above_bound"] is None else "NO")]
           for r in s13["wall_arm"]["rows"]],
          ["k_w (m/day)", "EPANET C", "Extra decay vs k_w=0", "Mass-transfer-free bound",
           "Above bound"]))

s15 = load("step15_unit_equivalence.json")
rows = []
for name, c in s15["comparisons"].items():
    # "legacy" is the code's own name for the superseded configuration and was never defined
    # as such in the paper, where the caption calls it superseded
    rows.append([name.replace("_", " ").replace("legacy", "superseded"),
                 "%.2e" % c["max_abs_diff_mg_L"],
                 "%.2e" % c["max_rel_diff"],
                 "%d -> %d" % (c["n_points_below_threshold_ref"],
                               c["n_points_below_threshold_other"]),
                 "%d of %d" % (c["n_cells_with_different_hours_below"], c["n_cells"])])
add("B3",
    "Table B3. Unit-equivalence test over 256 leading Sobol candidates. The corrected "
    "implementation reproduces the superseded one to floating-point precision. The units-only arm "
    "shows what correcting the concentrations while leaving the solver tolerance at the EPANET "
    "default would have cost.",
    table(rows, ["Comparison", "Max absolute difference (mg L-1)", "Max relative difference",
                 "Points below threshold", "Duration cells that differ"]))

# ---------------------------------------------------------------- C: threshold and displacement
s3 = load("step3_threshold.json")
rows = []
for r in s3["rows"]:
    rows.append(["%.3f" % r["threshold"], r["count"], "%.1f%%" % (100 * r["retention"]),
                 ] +
                ["%.4f (%.1f%%)" % (r[z]["sd"], 100 * r[z]["sd_retained"])
                 for z in ("old", "avg", "new")])
add("C1",
    "Table C1. Behavioural-threshold sweep for the informal comparator. Tightening the cut-off "
    "removes candidates but cannot restore the information factor the score omits, so the "
    "retained widths stay far above the formal values.",
    table(rows, ["Threshold (mg L-1)", "Retained", "Retention", "old SD (retained)",
                 "average SD (retained)", "new SD (retained)"]))

s4d = load("step4d_displaced_robust.json")
rows = []
for design, d in s4d["designs"].items():
    for scheme, sc in d["by_scheme"].items():
        for z in ("old", "avg", "new"):
            g = sc["groups"][z]
            rows.append([design, scheme.replace("_", " "), z,
                         "%.0f%% [%.0f-%.0f]" % (100 * g["gap_med"], 100 * g["gap_iqr"][0],
                                                 100 * g["gap_iqr"][1]),
                         "%.0f%%" % (100 * g["sd_ret_med"])])
add("C2",
    "Table C2. Displaced-prior recovery over 30 noise realisations, as the median fraction of the "
    "imposed displacement recovered with its interquartile range, and the retained prior width. "
    "DOWN displaces every coefficient towards stronger decay; OLDUP displaces the old coefficient "
    "the other way.",
    table(rows, ["Design", "Rule", "Coefficient", "Displacement recovered [25-75%]",
                 "Prior SD retained"]))

# ---------------------------------------------------------------- D: Fisher mechanics
s7 = load("step7_fisher.json")
rows = []
for case, c in s7["cases"].items():
    rows.append([case, "%.1f" % c["condition_number"]] +
                ["%.4f (%.0f%%)" % (c["coef"][z]["crlb"], 100 * c["coef"][z]["crlb_over_prior"])
                 for z in ("old", "average", "new")])
add("D1",
    "Table D1. Cramer-Rao bounds under nested nuisance models, obtained from the marginal (Schur) "
    "information. Percentages are the bound as a fraction of the prior standard deviation.",
    table(rows, ["Case", "Condition number", "old", "average", "new"]))

rows = []
for z, c in s7["fd_step_convergence"].items():
    for r in c["rows"]:
        rows.append([z, "%.4g" % r["step"], "%.1f%%" % (100 * r["step_over_truth"]),
                     "yes" if r["is_default"] else "", "%.6f" % r["crlb"]])
add("D2",
    "Table D2. Finite-difference step convergence for the Jacobian. Steps are scale-dependent "
    "rather than shared, because one absolute step would be a 2% perturbation of the old "
    "coefficient and a 40% perturbation of the new one.",
    table(rows, ["Coefficient", "Step", "As % of truth", "Selected", "Case-A bound"]))

g = s7["fisher_geometry"]
add("D3",
    "Table D3. Eigenvalues of the Case-A information matrix in raw and prior-scaled form. The "
    "raw condition number is dominated by the factor of twenty between the coefficient scales; "
    "on the dimensionless matrix the spread is far smaller.",
    table([["Raw"] + ["%.0f" % v for v in g["eigenvalues_raw"]] +
           ["%.0f" % g["condition_number_raw"]],
           ["Prior-scaled"] + ["%.1f" % v for v in g["eigenvalues_prior_scaled"]] +
           ["%.1f" % g["condition_number_prior_scaled"]]],
          ["Form", "Eigenvalue 1", "Eigenvalue 2", "Eigenvalue 3", "Condition number"]))

ar1 = load("step7c_ar1.json")
rows = []
for r in ar1["rho_sweep"]:
    rows.append(["%.1f" % r["rho"], "%.0f" % r["n_eff"]] +
                ["%.4f (%.0f%%)" % (r["coef"][z]["crlb"],
                                    100 * r["coef"][z]["crlb_over_prior"])
                 for z in ("old", "average", "new")] +
                [" / ".join("%.2f" % r["coef"][z]["widening_vs_indep"]
                            for z in ("old", "average", "new"))])
add("D4",
    "Table D4. Case-A bounds under an assumed AR(1) residual covariance, swept over the "
    "correlation parameter. The correlation is assumed rather than estimated, since the "
    "baseline observations are generated independently; the sweep therefore bounds how much "
    "the reported precision would degrade if the residuals were in fact correlated. Widening "
    "is the bound relative to the independent case, given for old, average and new. The three "
    "coefficients widen almost together up to a correlation of 0.4 and separate beyond it, "
    "the new coefficient least and the average coefficient most.",
    table(rows, ["Correlation", "Effective n", "old", "average", "new",
                 "Widening (old / avg / new)"]))

# ---------------------------------------------------------------- E: structural heterogeneity
s5c = load("step5c_jitter_sweep.json")
rows = []
for r in s5c["rows"]:
    for z in ("old", "average", "new"):
        zz = r["zones"][z]
        rows.append(["%.0f%%" % (100 * r["jitter"]), z, "%.4f" % zz["arith"],
                     "%+.4f" % zz["bias"], "%+.2f" % zz["bias_in_sd"],
                     "%+.4f" % zz["bias_increment_vs_control"]])
add("E1",
    "Table E1. Symmetric within-zone heterogeneity under the primary rule. The increment is the "
    "bias net of a paired homogeneous control run on the same noise, which is the part that can "
    "be attributed to the imposed structure.",
    table(rows, ["Jitter", "Zone", "True arithmetic mean", "Raw bias", "Raw bias / same-arm SD",
                 "Increment over control"]))

fe = {z: s5c["rows"][1]["field_ensemble"][z] for z in ("old", "average", "new")
      if "field_ensemble" in s5c["rows"][1]}
if fe:
    add("E2",
        "Table E2. The same design repeated over 25 independent heterogeneity fields at plus or "
        "minus 20%. The mean increment is a fraction of the field-to-field scatter in every zone, "
        "so no systematic bias is distinguishable from the choice of field.",
        table([[z, "%+.4f" % fe[z]["bias_mean"], "%.4f" % fe[z]["bias_sd"],
                "%+.4f" % fe[z]["increment_mean"], "%.4f" % fe[z]["increment_sd"],
                "%.2f" % fe[z]["increment_mean_over_sd"]] for z in fe],
              ["Zone", "Bias mean", "Bias SD", "Increment mean", "Increment SD",
               # escaped: an unescaped bar is a cell separator, which is how this header once parsed into
               # eight cells against a six-cell rule and left the column blank in Word
               r"\|Increment mean\| / field SD"]))

s5d = load("step5d_structured.json")
rows = []
for z, d in s5d["zones"].items():
    for scheme, sc in d["by_scheme"].items():
        rows.append([z, scheme.replace("_", " "), "%.4f" % d["arith"], "%.4f" % d["lenwt"],
                     "%.4f" % sc["mean"], "%+.4f" % sc["bias"],
                     "%+.4f" % sc["homogeneous_baseline_offset"],
                     "%.0f%%" % (100 * sc["shift_frac_of_lenwt_gap"])])
add("E3",
    "Table E3. Length-structured heterogeneity in the reference realisation. The homogeneous "
    "baseline offset is what the same rule returns with no structural error at all and must be "
    "subtracted before a shift is read as structural.",
    table(rows, ["Zone", "Rule", "Arithmetic mean", "Length-weighted proxy", "Posterior mean",
                 "Shift", "Homogeneous offset", "Fraction of gap"]))

# The amplitude sweep had no appendix table, so Section 3.4.2 was the only place the repeated
# realisation intervals appeared. Cutting them from the body without this would have left the
# zone-ordering claim resting on medians alone.
ZS = (("old", "old"), ("average", "average"), ("new", "new"))
rows = []
for r in load("step5d_structured.json")["correlation_dose_response"]["rows"]:
    cells = []
    for _, z in ZS:
        for s in ("formal_censored", "informal_glue"):
            c = r["by_scheme"][s][z]
            if c.get("shift_frac_med") is None:
                cells.append("—")
            else:
                lo, hi = c["shift_frac_5_95"]
                cells.append("%.2f [%.2f, %.2f]" % (c["shift_frac_med"], lo, hi))
    rows.append(["%.2f" % r["corr"]] + cells)
add("E4",
    "Table E4. Length-structured amplitude sweep over 30 noise realisations. Entries are the "
    "median raw displacement fraction with its 5 to 95% range, where 0 is the arithmetic zone "
    "mean and 1 the length-weighted reference. The fraction is undefined at an amplitude of zero "
    "because the two references coincide.",
    table(rows, ["Amplitude"] + ["%s %s" % (z, s) for _, z in ZS
                                 for s in ("formal", "informal")]))

# ---------------------------------------------------------------- F: per-arm sensor error
s8c = load("step8c_bias_bynode.json")
coef, rank = [], []
for r in s8c["rows"]:
    if r["scheme"] != "formal_censored":
        continue
    key = [r["node"], r["zone"], "%+.3f" % r["offset"]]
    coef.append(key + ["%+.2f" % r["d_old_over_sd"], "%+.2f" % r["d_avg_over_sd"],
                       "%+.2f" % r["d_new_over_sd"], int(r["n_censored_med"]),
                       "%.0f%s" % (r["ess_med"], " *" if r["sampling_limited"] else "")])
    rank.append(key + ["%.4f" % r["risk_spearman_vs_unbiased"],
                       "%.2f" % r["risk_top6_jaccard_vs_unbiased"],
                       "%.2f" % r["deficit_top6_jaccard_vs_unbiased"]])
# Ten columns do not fit the portrait text block: the headers alone run to 124 characters against
# about eighty, so Word broke words mid-syllable to make them fit. The data are at most seven
# characters wide, so the split is between what each half is about rather than a squeeze: the
# coefficient displacements here, the ranking statistics they propagate to next, both keyed by
# the same three columns so a reader can line up any arm across the two.
add("F1",
    "Table F1. Coefficient displacement in every arm of the sensor-bias sweep, under the primary "
    "rule, in baseline posterior standard deviations. A starred ESS is below the median-ESS "
    "criterion of Section 2.5.1, so that arm's displacement rests on a coarser library resolution "
    "than the rest. Table F2 gives what the same arms do to the risk ranking.",
    table(coef, ["Node", "Zone", "Offset (mg L-1)", "old", "average", "new",
                 "Median censored points", "Median ESS"]))

add("F2",
    "Table F2. Rank statistics for the same arms as Table F1, comparing the 92-junction risk "
    "field with the unbiased case.",
    table(rank, ["Node", "Zone", "Offset (mg L-1)", "Spearman, duration", "Jaccard, duration",
                 "Jaccard, deficit"]))

s8d = load("step8d_sensor_drift.json")
rows = []
for node, arms in s8d["equivalence"].items():
    for mag, e in arms.items():
        rows.append([node, mag, "%+.4f" % e["drift_shift"], "%+.4f" % e["const_mean_shift"],
                     "%+.4f" % e["const_end_shift"], "%.2f" % e["drift_over_const_mean"],
                     "%.2f" % e["drift_over_const_end"]])
add("F3",
    "Table F3. Sensor drift against its two constant-bias controls on identical noise. The ramp is "
    "correct at the start of the window and off by the end offset at its end; the controls hold a "
    "constant bias at half that offset and at the whole of it, which is what the last two columns "
    "divide by. A monotone drift acts, to within about a tenth, as a constant bias at its mean "
    "over the window.",
    table(rows, ["Node", "End offset", "Drift shift", "Constant at half", "Constant at end",
                 "Drift / half", "Drift / end"]))

# ---------------------------------------------------------------- G: risk-assessment mechanics
s10 = load("step10_risk_metrics.json")
a = s10["age_risk_association"]
bb = a["block_bootstrap"]
add("G1",
    "Table G1. Water-age association and the resampling design behind its interval. Whole spatial "
    "blocks are resampled rather than individual junctions, because junctions sharing pipes, flow "
    "paths and tank states are not independent samples. The interval is a descriptive width, not "
    "a significance test.",
    table([["Spearman, duration", "%.4f" % a["spearman_dur"]],
           ["Spearman, deficit", "%.4f" % a["spearman_def"]],
           ["Pearson, duration", "%.4f" % a["pearson_dur"]],
           ["Junctions", a["n_junctions"]],
           ["Blocks", bb["n_blocks"]],
           ["Block rule", bb["block_rule"]],
           ["Resamples", bb["n_resamples"]],
           ["95% interval, Spearman duration",
            "[%.3f, %.3f]" % tuple(bb["ci95_spearman_dur"])]],
          ["Quantity", "Value"]))

nav = s10["network_averages"]
G2LAB = {"E_duration_h": "Expected duration, E[D] (h)",
         "E_deficit": "Expected deficit, E[A] (mg L-1 h)",
         "min_C": "Mean member minimum (mg L-1)", "minC": "Mean member minimum (mg L-1)"}
rows = [[G2LAB.get(m, m.replace("_", " ")), "%.4f" % v["unweighted_all_junctions"],
         "%.4f" % v["consumer_only"], "%.4f" % v["demand_weighted"]]
        for m, v in nav["by_metric"].items()]
add("G2",
    "Table G2. Three network averages for each risk metric. Consumer-only is worse than "
    "unweighted while demand-weighted is better, which can only happen if the risk concentrates "
    "at small consumers.",
    table(rows, ["Metric", "Unweighted (all 92)",
                 "Consumer-only (%d)" % nav["n_consumer_junctions"], "Demand-weighted"]))

s12 = load("step12_scenarios.json")
u = s12["uncertainty_sources"]
add("G3",
    "Table G3. Scenario draw specification. One draw per retained ensemble member is reused "
    "across every scenario and dose as a common random number, so scenario differences are "
    "physical rather than Monte Carlo noise.",
    table([["Retained draws", "%d of %d (weight above %g of the maximum)"
            % (s12["n_retained_draws"], s12["n_total_draws"], s12["weight_floor_relative"])],
           ["Discarded weight mass", "%.2e" % s12["discarded_weight_mass"]],
           ["Bulk activation energy", "N(%.0f, %.0f^2) J/mol, truncated below at %.0f"
            % (u["Ea_bulk_J_per_mol"]["mean"], u["Ea_bulk_J_per_mol"]["sd"],
               u["Ea_bulk_J_per_mol"]["truncated_at"])],
           ["Realised bulk range", "%.0f to %.0f J/mol"
            % (u["Ea_bulk_J_per_mol"]["drawn_min"], u["Ea_bulk_J_per_mol"]["drawn_max"])],
           ["Wall activation energy", "N(%.0f, %.0f^2) J/mol, truncated below at %.0f"
            % (u["Ea_wall_J_per_mol"]["mean"], u["Ea_wall_J_per_mol"]["sd"],
               u["Ea_wall_J_per_mol"]["truncated_at"])],
           ["Realised wall range", "%.0f to %.0f J/mol"
            % (u["Ea_wall_J_per_mol"]["drawn_min"], u["Ea_wall_J_per_mol"]["drawn_max"])],
           ["Water temperature", u["water_temperature_C"]],
           ["Realised temperature range", "%.2f to %.2f degC"
            % tuple(u["water_temperature_span_C"])],
           ["Draws resampled at the bound", u["water_temperature_draws_resampled"]],
           ["Common random numbers", u["common_random_numbers"]]],
          ["Quantity", "Specification"]))

s16 = load("step16_window_sensitivity.json")
lo, hi = s16["risk_windows_h"][s16["paper_window"]]
rows, seen = [], None
for r in s16["rows"]:
    a, b = s16["risk_windows_h"][r["window"]]
    if r["arm"].startswith("k_b="):
        lab = "$k_b$ = " + r["arm"][4:].replace("-", "\u2212")
    else:                                     # "bias <node><signed offset>"
        node, off = r["arm"][5:-5], r["arm"][-5:]
        lab = "Sensor bias, node %s, %s mg L-1" % (node, off.replace("-", "\u2212"))
    rows.append([lab if r["arm"] != seen else "", "%d to %d" % (a, b),
                 "%.3f" % r["duration"]["spearman"], "%.2f" % r["duration"]["top6_jaccard"],
                 "%.3f" % r["deficit"]["spearman"], "%.2f" % r["deficit"]["top6_jaccard"],
                 "yes" if r["deficit"]["top6_jaccard"] >= r["duration"]["top6_jaccard"] else "no"])
    seen = r["arm"]
add("G4",
    "Table G4. Assessment-window sensitivity of the Section 3.6.2 shortlist comparisons. The "
    "calibration stays on the %d to %d h monitor series throughout; only the window over which "
    "the two risk metrics are integrated moves, to the adjacent cycle pair. Appendix A records "
    "that the deficit and water-age criteria are not met inside the horizon, so the risk field "
    "is not warm-up converged and the individual overlaps are expected to move; the question is "
    "whether the comparison between the two metrics moves with them. The %d to %d h rows "
    "reproduce the corresponding entries of Table 4. Every arm is a median over %d noise "
    "realisations, matching Steps 8b and 8c."
    % (s16["calibration_window_h"][0], s16["calibration_window_h"][1], lo, hi, s16["n_noise"]),
    table(rows, ["Perturbation", "Risk window (h)", "Duration $\\rho_s$", "Duration top-6 J",
                 "Deficit $\\rho_s$", "Deficit top-6 J", "Deficit at least as well preserved"]))

# Section 2.8.3 names two supplementary ranking diagnostics, Kendall correlation and alternative
# shortlist sizes. Both were computed and neither was shown, which left a methods statement with no
# reported result behind it. They are tabulated here from the artifacts that already held them.
s8b = load("step8b_kb_sensitivity.json")
KS = ["top3", "top5", "top6", "top10", "top15"]
rows = []
for kb in (-0.4, -0.6):
    r = next(x for x in s8b["rows"] if x["kb"] == kb)["by_scheme"][s8b["primary_weighting"]]
    for name, src_d in (("Duration, $\\bar{P}$", r["topk_jaccard_vs_kb_ref"]),
                        ("Deficit, E[A]", r["depth_deficit"]["topk_jaccard_vs_kb_ref"])):
        rows.append(["$k_b$ = \u2212%.1f" % abs(kb) if name.startswith("Duration") else "",
                     name] + ["%.2f" % src_d[k] for k in KS])
# counted rather than transcribed, because a hand-written "better at n of five" drifts silently
cmp_k = {}
for kb in (-0.4, -0.6):
    r = next(x for x in s8b["rows"] if x["kb"] == kb)["by_scheme"][s8b["primary_weighting"]]
    d, a = r["topk_jaccard_vs_kb_ref"], r["depth_deficit"]["topk_jaccard_vs_kb_ref"]
    cmp_k[kb] = (sum(a[k] > d[k] for k in KS), sum(abs(a[k] - d[k]) < 1e-12 for k in KS),
                 [k[3:] for k in KS if a[k] < d[k]])
add("G5",
    "Table G5. Shortlist-size sensitivity of the Table 4 comparisons, under a 20 per cent error in "
    "bulk decay. The main text reports k = 6. Under the underestimate the deficit shortlist is at "
    "least as well preserved at every size, better at %d of the five and equal at %d. Under the "
    "overestimate it is better at %d and worse at k = %s. The metric-specific advantage therefore "
    "belongs to the small shortlist an operator would act on rather than to ranking in general, "
    "which is the scale the section claims it at. The two metrics respond to shortlist size in "
    "opposite directions. Duration overlap is lowest at k = 3 and rises with every increase in "
    "k, whereas deficit overlap reaches its highest value by k = 5 and falls back at k = 10 "
    "and 15, where the ranking is decided by junctions whose risk is near zero."
    % (cmp_k[-0.4][0], cmp_k[-0.4][1], cmp_k[-0.6][0], " and ".join(cmp_k[-0.6][2])),
    table(rows, ["Perturbation", "Ranking metric"] + ["k = %s" % k[3:] for k in KS]))

s8 = load("step8_sensor_bias.json")
s8d = load("step8d_sensor_drift.json")


def rng(vals):
    v = sorted(vals)
    return "%.4f to %.4f" % (v[0], v[-1])


# Step 8 sweeps one monitor and includes a zero-offset row, which compares the unbiased field with
# itself: leaving it in put a 1.0000 at the top of both ranges that no arm earned. Step 8d's 24 rows
# are 8 drift ramps and their 16 constant-offset controls, not 24 drift arms. Both are split here.
bias_arms = [r for r in s8["rows"] if r["offset"] != 0.0]
drift = [r for r in s8d["rows"] if r["arm"] == "drift"]
ctrl = [r for r in s8d["rows"] if r["arm"] != "drift"]


add("G6",
    "Table G6. Kendall rank correlation beside Spearman, on the duration ranking, over every arm "
    "of the two sweeps that computed both. Step 8 sweeps the offset at node %s alone, so it is "
    "narrower than the %d-arm by-monitor sweep of Table F1, which did not compute Kendall. Its "
    "zero-offset row is excluded here because it compares the unbiased field with itself. Step "
    "8d's drift ramps are separated from the constant-offset controls they are paired against. "
    "The two coefficients agree that whole-network rank order survives sensor error. Kendall is "
    "the lower of the two throughout, because it penalises each discordant pair rather than each "
    "squared rank displacement, so reading it in place of the Spearman values quoted in Section "
    "3.6.2 would lower those figures slightly and change nothing that section concludes."
    % (s8["bias_node"], len([r for r in load("step8c_bias_bynode.json")["rows"]
                             if r["scheme"] == "formal_censored"])),
    table([["Constant offset at node %s (Step 8)" % s8["bias_node"], len(bias_arms),
            rng([r["spearman_vs_unbiased"] for r in bias_arms]),
            rng([r["kendall_vs_unbiased"] for r in bias_arms])],
           ["Linear drift at nodes 15 and 231 (Step 8d)", len(drift),
            rng([r["risk_spearman_vs_unbiased"] for r in drift]),
            rng([r["risk_kendall_vs_unbiased"] for r in drift])],
           ["Constant-offset controls for those ramps (Step 8d)", len(ctrl),
            rng([r["risk_spearman_vs_unbiased"] for r in ctrl]),
            rng([r["risk_kendall_vs_unbiased"] for r in ctrl])]],
          ["Perturbation family", "Arms", "Spearman range", "Kendall range"]))

# ---------------------------------------------------------------- H: banding, then the register
# The register's last three columns were uninterpretable on their own: the caption used to say the
# banding rules were in the main text, and they were not there. They are generated here from the
# artifact rather than transcribed, so the table cannot drift from the scale the code applied.
s12 = load("step12_scenarios.json")
d12 = s12["definitions"]
q1, q2 = s12["consequence_terciles_L_s"]
brows = []
for band, crit in d12["likelihood_bands_on_P_min"].items():
    brows.append(["Likelihood, from $P_{\\min}$" if not brows else "", band, crit])
n = len(brows)
for band, crit in d12["severity_bands_on_E_duration_h"].items():
    brows.append(["Severity, from E[D]" if len(brows) == n else "", band, crit])
n = len(brows)
for label, score in d12["consequence_scores"].items():
    name = label.split(" (")[0]
    rng = {0: "demand = 0, a non-consumer junction",
           1: "0 < demand <= %.2f L/s" % q1,
           2: "%.2f < demand <= %.2f L/s" % (q1, q2),
           3: "demand > %.2f L/s" % q2}[score]
    brows.append(["Consequence, from demand" if len(brows) == n else "", name, rng])
n = len(brows)
for band, crit in d12["risk_band_mapping"].items():
    brows.append(["Risk band, from score" if len(brows) == n else "", band, crit])
add("H1",
    "Table H1. Risk banding. The likelihood and severity axes are scored 1 to 5 on the criteria "
    "below and consequence 0 to 3; a risk score is the axis score times the consequence score, "
    "and the band follows from that score. The severity edges are pre-declared hours rather than "
    "quantiles of this network's own results, so the same scale applies across scenarios. The "
    "governing band is the higher of the breach-probability band and the severity band, since "
    "neither axis alone may reduce an action.",
    table(brows, ["Axis", "Band", "Criterion"]))

with open(os.path.join(CACHE, "step12_risk_register.csv")) as f:
    reg = list(csv.DictReader(f))
# Column names are required, not filtered. An earlier version kept only the columns that
# happened to exist, so when step12 renamed six of them the table lost them in silence and the
# sort key went with them: every row scored 0, sorted() is stable, and the register shipped in
# CSV order under a caption promising cumulative-deficit order. A generator for a published
# table must fail on a schema it does not recognise.
# Split in two rather than run fourteen columns across a page. Both halves are keyed by node
# and sorted identically, so row n of one is row n of the other.
QUANT = ["node", "P_min_current", "P_min_heatwave", "P_min_heat_ageing", "P_bar_current",
         "E_duration_current_h", "E_deficit_current_mgL_h", "demand_L_s"]
BANDS = ["node", "likelihood", "severity", "consequence", "risk_band_current",
         "risk_band_severity", "risk_band_governing"]
missing = [k for k in QUANT + BANDS if k not in reg[0]]
if missing:
    raise SystemExit("step12_risk_register.csv is missing columns the register table needs: "
                     + ", ".join(sorted(set(missing))))
ordered = sorted(reg, key=lambda x: float(x["E_deficit_current_mgL_h"]), reverse=True)
# risk_band_current is the band from the breach-probability axis; naming it "band" alone would
# sit ambiguously beside "band severity" and "band governing"
RENAME = {"risk_band_current": "band breach", "risk_band_severity": "band severity",
          "risk_band_governing": "band governing"}
RENAME.update({"P_min_current": "P_min", "P_min_heatwave": "P_min, heatwave",
               "P_min_heat_ageing": "P_min, heatwave + ageing", "P_bar_current": "P_bar",
               "E_duration_current_h": "E[D] (h)",
               "E_deficit_current_mgL_h": "E[A] (mg L-1 h)",
               "demand_L_s": "Demand (L s-1)"})
head = lambda ks: [RENAME.get(k, k.replace("_current", "").replace("_", " ")) for k in ks]
add("H2",
    "Table H2. Risk register, quantitative half, all 92 junctions ordered by expected cumulative "
    "deficit. The reference condition throughout is Scenario A, which is the 12 degree Celsius "
    "mean carrying the prescribed temperature and activation-energy uncertainty of Table G3, and "
    "not the exact-reference-temperature calibrated field of Section 3.6.1. The two differ: the "
    "network-mean expected deficit is %.4f mg L-1 h at the exact reference temperature against "
    "%.4f under Scenario A, while the count of junctions more likely than not to breach and the "
    "demand they serve are identical at %d and %.1f L s-1. P min is the probability of at least "
    "one breach in the window, under Scenario A and under the two warmer scenarios. The remaining "
    "columns are Scenario A. Table H3 carries the same junctions in the same order."
    % (s12["reference_at_T_ref_exact"]["net_mean_E_deficit"],
       s12["scenario_summary"][0]["net_mean_E_deficit"],
       s12["reference_at_T_ref_exact"]["P_min_gt_0.5_nodes"],
       s12["reference_at_T_ref_exact"]["demand_at_risk_L_s"]),
    table([[r[k] for k in QUANT] for r in ordered], head(QUANT)))
add("H3",
    "Table H3. Risk register, banded half, in the order of Table H2. Likelihood, severity and "
    "consequence are the axis scores of Table H1; the governing band is the higher of the breach "
    "and severity bands.",
    table([[r[k] for k in BANDS] for r in ordered], head(BANDS)))

# ---------------------------------------------------------------- write
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as f:
    f.write("<!-- Generated by Net3/appendix_tables.py. Do not edit by hand. -->\n\n")
    f.write("\n".join(blocks))
print("wrote %s  (%d tables)" % (os.path.relpath(OUT, HERE), len(blocks)))
