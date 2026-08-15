"""Cross-check RESULTS_LOG.md, README.md, the figures and the JSON artifacts against each other.

The prose numbers are copied in by hand, so they drift whenever a script is re-run. These checks
detect that drift until the documents are generated from the artifacts directly.

WHAT THIS SCRIPT DOES AND DOES NOT ESTABLISH — read before quoting a passing run.

  It establishes, for the numbers it covers:
    * a registered CLAIM (the list below) equals a named value at a named JSON path. This is the only
      check with real teeth, because both the sentence and the source are pinned;
    * every result artifact declares which weighting produced it, so a table cannot be read under the
      wrong inference rule;
    * no log number is a near-miss of a JSON value, which is the signature of a stale transcription;
    * no forbidden phrasing survives (superseded configurations, over-claimed wording);
    * the environment and the frozen model file match the ones that produced the cache.

  It does NOT establish:
    * that the artifacts themselves are correct. If a script computes the wrong thing, every check
      here confirms only that the prose faithfully reports a wrong number;
    * that an unregistered number in the prose came from anywhere. The general number check is a
      drift detector, not a provenance proof: with thousands of values in the pool, a two-decimal
      figure matches something by coincidence, which is why table cells are held to their own
      section's artifact from three decimals up and why the claim registry exists;
    * SEMANTIC consistency. Nothing here notices that two sections draw opposite conclusions, that a
      caveat has gone missing, or that a correct number is described with the wrong words;
    * that a figure shows what its caption says. Figure freshness is a file-timestamp comparison
      only.

A passing run therefore means "no detected drift in the covered numbers", not "the log is verified".

Usage:
    python validate_artifacts.py            # summary; exit 1 if any check fails
    python validate_artifacts.py --verbose  # also list unmatched (probably prose) numbers
    python validate_artifacts.py --release  # additionally require a clean-tree manifest at HEAD
"""
import bisect
import json
import os
import subprocess
import re
import sys

import provenance

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
CACHE = os.path.join(HERE, "baseline_cache")
LOG = os.path.join(HERE, "RESULTS_LOG.md")
DOCS = {"RESULTS_LOG.md": LOG,
        "README.md": os.path.join(REPO, "README.md"),
        "README.en.md": os.path.join(REPO, "README.en.md"),
        "REVISION_RESPONSE_MATRIX.md": os.path.join(REPO, "REVISION_RESPONSE_MATRIX.md"),
        # The submitted paper is the document that matters most and was the last to be covered.
        # Its headline numbers are now pinned to named JSON paths like any other claim, so a
        # rebuild cannot quietly carry a figure the artifacts no longer support.
        "paper.md": os.path.join(REPO, "00-thesis", "paper.md")}

# figures/ and Figures/ are the same directory on a case-insensitive filesystem; resolve whichever
# name exists so the check also works on Linux.
FIGDIR = next((os.path.join(HERE, d) for d in ("figures", "Figures")
               if os.path.isdir(os.path.join(HERE, d))), os.path.join(HERE, "figures"))

STALE_REL_TOL = 0.02      # a near-miss this close to a JSON value is treated as a stale copy
MIN_INT_TO_CHECK = 10     # below this, integers are too ambiguous to match usefully
MAX_DP = 6                # deepest rounding compared; a percent form needs two more than written

# Canonical section label -> the artifacts whose numbers that section may quote.
# The labels are also the reference for the numbering check: a "## Step X" heading in the log must
# appear here, which is what keeps section labels and script filenames aligned.
SECTIONS = {
    "Step 0": ["step0_warmup_convergence.json"],
    "Step 1": ["baseline_meta.json"],
    "Step 2": ["baseline_meta.json", "step3_threshold.json"],
    # Step 3 tabulates the formal-likelihood row from Step 1 for contrast, so it legitimately
    # quotes that artifact as well; declaring it keeps the strict table rule usable.
    "Step 3": ["step3_threshold.json", "baseline_meta.json"],
    "Step 4": ["step4_displaced_prior.json"],
    "Step 4b": ["step4b_sensitivity.json"],
    "Step 4d": ["step4d_displaced_robust.json"],
    "Step 5": ["step5_structural_error.json"],
    "Step 5a": ["step5_structural_error.json"],
    "Step 5c": ["step5c_jitter_sweep.json"],
    "Step 5c-ensemble": ["step5c_jitter_sweep.json"],
    "Step 5d": ["step5d_structured.json"],
    "Step 5d-dose": ["step5d_structured.json"],
    "Step 5e": ["step5d_structured.json"],   # grid geometry is recorded there
    "Step 6": ["step6_noise_sensitivity.json"],
    "Step 7": ["step7_fisher.json"],
    "Step 7b": ["step7b_profile.json"],
    "Step 7c": ["step7c_ar1.json", "step7c_profile_ar1.json"],
    # the bias-location sweep table lives inside the Step 8 section
    "Step 8": ["step8_sensor_bias.json", "step8c_bias_bynode.json"],
    "Step 8b": ["step8b_kb_sensitivity.json"],
    "Step 8c": ["step8c_bias_bynode.json"],
    "Step 8d": ["step8d_sensor_drift.json"],
    "Step 9": ["step9_zeroclip.json"],
    "Step 10": ["step10_risk_metrics.json"],
    "Step 11": ["step11_loo.json"],
    "Step 11b": ["step11_loo.json"],
    "Step 11c": ["step11_loo.json"],
    "Step 11d": ["step11_loo.json"],
    "Step 12": ["step12_scenarios.json"],
    "Step 13": ["step13_known_answer.json"],
    "Step 14": ["step14_repeated_noise.json"],
    "Step 15": ["step15_unit_equivalence.json", "step15_full_regression.json"],
    "Step 15b": ["step15_full_regression.json"],
}

# Every result artifact must say which weighting produced it. Without this a table of numbers is
# only implicitly attached to an inference rule, which is how an informal-GLUE result ends up quoted
# under a formal heading. Accepted keys, in order of preference.
WEIGHTING_KEYS = ("primary_weighting", "weighting")
# Artifacts that hold no weighted ensemble and so have no weighting to declare: a warm-up
# convergence test, Fisher/profile geometry, and an analytic known-answer test.
WEIGHTING_EXEMPT = {"cache_manifest.json", "step0_warmup_convergence.json",
                    # properties of the deterministic noise-free truth: no ensemble is weighted
                    "forward_baseline.json",
                    "step4b_sensitivity.json",   # single-parameter RMSE curves; no ensemble weights
                    "step7_fisher.json",         # a-priori Fisher/CRLB; ensemble SD is a side check
                    "step7b_profile.json", "step7c_ar1.json", "step7c_profile_ar1.json",
                    "step9_zeroclip.json",       # compares two formal likelihoods on one realisation
                    "step13_known_answer.json",
                    "step15_unit_equivalence.json",   # raw fields, no ensemble weighting
                    "step15_full_regression.json"}    # array + git-state comparison

FIGURE_SOURCE = {
    "step0_warmup_convergence.png": "step0_warmup_convergence.json",
    "step4b_sensitivity_curves.png": "step4b_sensitivity.json",
    "step4d_displaced_robust.png": "step4d_displaced_robust.json",
    "step5_structural_error.png": "step5_structural_error.json",
    "step5c_jitter_sweep.png": "step5c_jitter_sweep.json",
    "step5d_structured.png": "step5d_structured.json",
    "step6_noise_sensitivity.png": "step6_noise_sensitivity.json",
    "step7_fisher.png": "step7_fisher.json",
    "step7b_profile.png": "step7b_profile.json",
    "step8_sensor_bias.png": "step8_sensor_bias.json",
    "step8b_kb_sensitivity.png": "step8b_kb_sensitivity.json",
    "step8d_sensor_drift.png": "step8d_sensor_drift.json",
    "step9_zeroclip.png": "step9_zeroclip.json",
    "step10_risk_metrics.png": "step10_risk_metrics.json",
    "step11_loo.png": "step11_loo.json",
    "step12_summary.png": "step12_scenarios.json",
    "step12_scenario_maps.png": "step12_scenarios.json",
    "step12_ageing_delta.png": "step12_scenarios.json",
    "step14_repeated_noise.png": "step14_repeated_noise.json",
}
FIGURE_STALE_S = 60      # a figure written in the same script run may predate its JSON by seconds

# ---------------------------------------------------------------- claim registry
# Each entry pins ONE sentence in ONE document to ONE value at ONE JSON path. This is the check with
# real teeth: the general number scan can only notice that a figure is near something, whereas a
# registered claim fails both when the number changes AND when the sentence is reworded so that the
# anchor disappears — which is the failure mode that let a whole stale table survive before.
#
# Fields: (document, regex with exactly one capture group, artifact, dotted JSON path, scale)
# `scale` multiplies the stored value before comparison, so a stored fraction can be claimed as a
# percentage. Precision is taken from the number as written.
CLAIMS = [
    # --- headline: how much prior width each weighting retains (the project's central contrast) ---
    ("RESULTS_LOG.md", r"formal censored[^\n]{0,120}?retains\s+([0-9.]+)\s*/\s*[0-9.]+\s*/\s*[0-9.]+\s*%",
     "baseline_meta.json", "summary/schemes/formal_censored/coef/old/sd_retained", 100),
    ("RESULTS_LOG.md", r"formal censored[^\n]{0,120}?retains\s+[0-9.]+\s*/\s*([0-9.]+)\s*/\s*[0-9.]+\s*%",
     "baseline_meta.json", "summary/schemes/formal_censored/coef/avg/sd_retained", 100),
    ("RESULTS_LOG.md", r"formal censored[^\n]{0,120}?retains\s+[0-9.]+\s*/\s*[0-9.]+\s*/\s*([0-9.]+)\s*%",
     "baseline_meta.json", "summary/schemes/formal_censored/coef/new/sd_retained", 100),
    ("README.md", r"censored 正式似然只保留\s*\*\*([0-9.]+)\s*/\s*[0-9.]+\s*/\s*[0-9.]+%\*\*",
     "baseline_meta.json", "summary/schemes/formal_censored/coef/old/sd_retained", 100),
    ("README.md", r"censored 正式似然只保留\s*\*\*[0-9.]+\s*/\s*([0-9.]+)\s*/\s*[0-9.]+%\*\*",
     "baseline_meta.json", "summary/schemes/formal_censored/coef/avg/sd_retained", 100),
    ("README.md", r"censored 正式似然只保留\s*\*\*[0-9.]+\s*/\s*[0-9.]+\s*/\s*([0-9.]+)%\*\*",
     "baseline_meta.json", "summary/schemes/formal_censored/coef/new/sd_retained", 100),
    ("README.md", r"informal GLUE（草稿阈值 `0.12`）保留了先验宽度的 \*\*([0-9.]+)\s*/\s*[0-9.]+\s*/\s*[0-9.]+%\*\*",
     "baseline_meta.json", "summary/schemes/informal_glue_draft_thr/coef/old/sd_retained", 100),
    # --- the a-priori bound the posterior width is compared with ---
    ("RESULTS_LOG.md", r"Case[- ]A CRLB[^\n]{0,80}?old\s+([0-9.]+)",
     "step7_fisher.json", "cases/A: kw only/coef/old/crlb", 1),
    # --- Step 14: the repeated-sampling claims that replace "the estimator is efficient" ---
    ("RESULTS_LOG.md", r"empirical SD / CRLB[^\n]{0,80}?([0-9.]+)\s*/\s*[0-9.]+\s*/\s*[0-9.]+",
     "step14_repeated_noise.json",
     "by_scheme/formal_censored/coef/old/empirical_sd_over_crlb", 1),
    ("RESULTS_LOG.md", r"nominal 90% intervals cover[^\n]{0,80}?([0-9.]+)\s*/\s*[0-9.]+\s*/\s*[0-9.]+",
     "step14_repeated_noise.json", "by_scheme/formal_censored/coef/old/coverage/q90", 1),
    # --- the noise floor that a candidate RMSE may be compared with: window, not full record ---
    ("RESULTS_LOG.md", r"realised noise RMSE of \*?\*?([0-9.]+) mg/L\*?\*? on the same",
     "baseline_meta.json", "summary/noise_rmse", 1),
    ("RESULTS_LOG.md", r"all three at truth\) on the calibration window = \*\*([0-9.]+) mg/L",
     "step4b_sensitivity.json", "summary/noise_floor_rmse", 1),
    # --- warm-up: the residual drift that stops "converged" from being sayable ---
    ("RESULTS_LOG.md", r"residual cycle-to-cycle drift[^\n]{0,90}?([0-9.]+)\s*%",
     "step0_warmup_convergence.json", "per_criterion/risk_rel_dDeficit/worst_by_cycle[5]", 100),
    ("RESULTS_LOG.md", r"water age[^\n]{0,90}?still\s+([0-9.]+)\s*h between the last two cycles",
     "step0_warmup_convergence.json", "per_criterion/age_p95_dAge/worst_by_cycle[5]", 1),
    # --- risk: the descriptive association, and that no p-value is claimed ---
    ("RESULTS_LOG.md", r"Spearman ρ = ([0-9.]+)[^\n]{0,40}descriptive",
     "step10_risk_metrics.json", "age_risk_association/spearman_dur", 1),
    # --- sensor bias: displacement in posterior-SD units, primary rule ---
    ("RESULTS_LOG.md", r"\+0\.10 mg/L at node 15 moves[^\n]{0,60}?([0-9.]+)\s*posterior SD",
     "step8_sensor_bias.json", "rows[6]/shift_over_sd", 1),
    # --- the unit correction: is it a relabelling, and what does getting it half-right cost? ---
    ("RESULTS_LOG.md", r"unchanged to ([0-9.]+)\s*×\s*10⁻⁷ relative",
     "step15_unit_equivalence.json", "comparisons/corrected_vs_legacy/max_rel_diff", 1e7),
    ("RESULTS_LOG.md", r"introduces a ([0-9.]+)% maximum relative error",
     "step15_unit_equivalence.json", "comparisons/units_only_vs_legacy/max_rel_diff", 100),
    # --- the register's two axes: the disagreement count is the whole point of carrying both ---
    ("RESULTS_LOG.md", r"\*\*([0-9]+) in the same risk band, [0-9]+\s*\n?in a lower one",
     "step12_scenarios.json", "severity_axis/risk_band_agreement_consumers/same band", 1),
    ("RESULTS_LOG.md", r"in the same risk band,\s*\*?\*?([0-9]+)\s*\n?in a lower one",
     "step12_scenarios.json", "severity_axis/risk_band_agreement_consumers/severity band lower", 1),
    ("RESULTS_LOG.md", r"in a lower one, ([0-9]+) in a higher one",
     "step12_scenarios.json", "severity_axis/risk_band_agreement_consumers/severity band higher", 1),
    # --- Step 12's severity axis: the three statements an audit found overstated, each pinned to
    # the artifact so the wording cannot drift back. Both the "0 inversions" count and the margin
    # behind it are registered, because the count alone reads as a guarantee.
    ("RESULTS_LOG.md", r"the nearest is \*\*([0-9.]+) h\*\*",
     "step12_scenarios.json",
     "severity_axis/observed_ordering/min_hours_of_E_duration_to_first_inversion", 1),
    ("RESULTS_LOG.md", r"bound permits an inversion of ([0-9]+) band",
     "step12_scenarios.json",
     "severity_axis/ordering_is_not_guaranteed/max_inversion_permitted_by_the_bound", 1),
    ("RESULTS_LOG.md", r"depth spans a factor of ([0-9.]+)",
     "step12_scenarios.json", "severity_axis/depth_spread_over_movers/ratio", 1),
    # --- Step 12's escalation limits: the two counts an audit's follow-through produced ---
    ("RESULTS_LOG.md", r"\*\*([0-9]+) junctions \(17\.46 L/s\)\*\* can never be classified",
     "step12_scenarios.json",
     "escalation_test/two_distinct_reasons_a_band_cannot_move/n_junctions_capped_below_high", 1),
    ("RESULTS_LOG.md", r"node 243 goes to the \*\*full ([0-9]+) h window\*\*",
     "step12_scenarios.json",
     "escalation_test/two_distinct_reasons_a_band_cannot_move/worked_example/E_duration_h_D", 1),
    # --- Step 15b: the whole-repository regression. Registered because these three counts are
    # exactly the ones that were quoted for a year without an artifact behind them.
    ("RESULTS_LOG.md", r"\*\*26 artifacts, ([0-9]+) numeric fields\*\*",
     "step15_full_regression.json", "artifact_regression/n_numeric_fields_compared", 1),
    ("RESULTS_LOG.md", r"\*\*([0-9]+)\*\* fields moved by\s*\n?\s*more than 1e-4 relative",
     "step15_full_regression.json", "artifact_regression/n_fields_changed", 1),
    ("RESULTS_LOG.md", r"\*\*This log, ([0-9]+) numbers at",
     "step15_full_regression.json", "log_regression/n_numbers_at_post", 1),
    # --- sensor drift: the ratio that decides whether a drift needs its own analysis at all ---
    ("RESULTS_LOG.md", r"the ratio reaches\s+\*\*([0-9.]+)\*\*\s+at D = \+0\.10",
     "step8d_sensor_drift.json", "equivalence/231/+0.100/drift_over_const_mean", 1),
    # whitespace-tolerant: a markdown formatter pads table columns, and an anchor that depends on
    # the exact spacing silently stops matching the moment the table is reflowed
    ("RESULTS_LOG.md", r"\|\s*15 \(old\)\s*\|\s*−0\.100\s*\|\s*−([0-9.]+)\s*\|",
     "step8d_sensor_drift.json", "rows[0]/own_shift_over_sd", -1),
    # --- sensor accuracy: the answer that changed when the primary rule changed ---
    ("RESULTS_LOG.md", r"σ = 0\.10[^\n]{0,90}?retains\s+([0-9.]+)\s*/\s*[0-9.]+\s*/\s*[0-9.]+\s*% of the prior",
     "step6_noise_sensitivity.json",
     "rows[2]/by_scheme/formal_censored/old/sd_ret_med", 100),
    # --- profile likelihood: the continuous interval, which is the one that may be quoted ---
    ("RESULTS_LOG.md", r"continuous 95% interval for `k_w,old` is \[−([0-9.]+),",
     "step7b_profile.json",
     "continuous_profile/by_likelihood/censored/coef/old/lo", -1),
    # --- the submitted paper: the headline numbers of the abstract and conclusions ---
    ("paper.md", r"contracted all three coefficients to\s*\n?\s*([0-9.]+)[–-]",
     "baseline_meta.json", "summary/schemes/formal_censored/coef/old/sd_retained", 100),
    ("paper.md", r"contracted all three coefficients to\s*\n?\s*[0-9.]+[–-]([0-9.]+)%",
     "baseline_meta.json", "summary/schemes/formal_censored/coef/avg/sd_retained", 100),
    ("paper.md", r"displaced coefficients by as much as ([0-9.]+) baseline",
     "step8c_bias_bynode.json", "summary/max_own_shift_over_sd/value", -1),
    # measured over all 24 arms of the location sweep. An earlier version quoted the node-15
    # figure from step 8, which understated it: the arm producing the largest displacement is
    # not the arm producing the largest residual.
    # "never exceeded ... by more than 6.2%" read as a bound over every arm AND every realisation.
    # The stored quantity is the largest of the 24 arm-level medians, so the sentence now says so.
    ("paper.md", r"rose at most ([0-9.]+)% above the realised",
     "step8c_bias_bynode.json", "summary/max_residual_excess_pct", 1),
    # re-anchored to Section 3.2.3 when the Abstract dropped the figures for length. The numbers
    # are unchanged; only the sentence carrying them moved.
    ("paper.md", r"Cram\u00e9r\u2013Rao bound of ([0-9.]+), [0-9.]+ and [0-9.]+",
     "step14_repeated_noise.json", "by_scheme/formal_censored/coef/old/empirical_sd_over_crlb", 1),
    ("paper.md", r"Cram\u00e9r\u2013Rao bound of [0-9.]+, [0-9.]+ and ([0-9.]+)",
     "step14_repeated_noise.json", "by_scheme/formal_censored/coef/new/empirical_sd_over_crlb", 1),
    ("paper.md", r"median normalised mean absolute error of ([0-9.]+)%",
     "step11_loo.json",
     "unmonitored_validation/by_scheme/formal_censored/mean_abs_rel_error[0]", 100),
    ("paper.md", r"rank correlation stays between ([0-9.]+)\s+and",
     "step8b_kb_sensitivity.json", "rows[2]/by_scheme/formal_censored/risk_spearman_vs_kb_ref", 1),
    # --- Results section headline numbers, added after the audit found that the paper's several
    # hundred figures rested on eight anchored claims and a membership check with 14% power.
    # Each entry pins one written figure to one named JSON path; a rewording that drops the
    # anchor fails loudly rather than silently ending the check.
    # 3.1.1 forward baseline
    ("paper.md", r"network-wide median of ([0-9.]+) mg", "forward_baseline.json",
     "network_median_mgL", 1),
    # 3.2.1 the three rules on the reference realisation
    ("paper.md", r"posterior standard deviation was [0-9.]+%, ([0-9.]+)% and", "baseline_meta.json",
     "summary/schemes/formal_censored/coef/avg/sd_retained", 100),
    ("paper.md", r"posterior standard deviation was [0-9.]+%, [0-9.]+% and ([0-9.]+)%",
     "baseline_meta.json", "summary/schemes/formal_censored/coef/new/sd_retained", 100),
    ("paper.md", r"left ([0-9.]+)%, [0-9.]+% and [0-9.]+% of the prior", "baseline_meta.json",
     "summary/schemes/informal_glue_draft_thr/coef/old/sd_retained", 100),
    ("paper.md", r"left [0-9.]+%, ([0-9.]+)% and [0-9.]+% of the prior", "baseline_meta.json",
     "summary/schemes/informal_glue_draft_thr/coef/avg/sd_retained", 100),
    ("paper.md", r"returned ([0-9.]+)%, [0-9.]+% and [0-9.]+%;", "baseline_meta.json",
     "summary/schemes/formal_iid/coef/old/sd_retained", 100),
    # 3.2.2 threshold sweep and displaced prior
    ("paper.md", r"fell only to ([0-9.]+)%,", "step3_threshold.json",
     "rows[0]/old/sd_retained", 100),
    ("paper.md", r"fell only to [0-9.]+%, ([0-9.]+)%", "step3_threshold.json",
     "rows[0]/avg/sd_retained", 100),
    ("paper.md", r"7084 to ([0-9]+) candidates", "step3_threshold.json", "rows[0]/count", 1),
    ("paper.md", r"retaining\s+([0-9.]+)%, [0-9.]+% and [0-9.]+% of the displaced",
     "step4d_displaced_robust.json",
     "designs/DOWN/by_scheme/formal_censored/groups/old/sd_ret_med", 100),
    ("paper.md", r"recovered ([0-9]+)%, [0-9]+% and [0-9]+% of the imposed",
     "step4d_displaced_robust.json",
     "designs/DOWN/by_scheme/formal_censored/groups/old/gap_med", 100),
    ("paper.md", r"recovered [0-9]+%, ([0-9]+)% and [0-9]+% of the imposed",
     "step4d_displaced_robust.json",
     "designs/DOWN/by_scheme/formal_censored/groups/avg/gap_med", 100),
    # 3.2.3 Fisher, profile and repeated noise
    ("paper.md", r"equivalent to\s+([0-9.]+)%, [0-9.]+% and [0-9.]+% of the prior",
     "step7_fisher.json", "cases/A: kw only/coef/old/crlb_over_prior", 100),
    ("paper.md", r"equivalent to\s+[0-9.]+%, ([0-9.]+)% and [0-9.]+% of the prior",
     "step7_fisher.json", "cases/A: kw only/coef/average/crlb_over_prior", 100),
    ("paper.md", r"coverage of ([0-9.]+), [0-9.]+ and [0-9.]+",
     "step14_repeated_noise.json",
     "by_scheme/formal_censored/coef/old/coverage/q90", 1),
    ("paper.md", r"coverage of [0-9.]+, [0-9.]+ and ([0-9.]+)",
     "step14_repeated_noise.json",
     "by_scheme/formal_censored/coef/new/coverage/q90", 1),
    # 3.3.1 noise and assumed autocorrelation
    ("paper.md", r"leaving ([0-9.]+)%, [0-9.]+% and [0-9.]+% of the prior standard deviation\. The",
     "step6_noise_sensitivity.json",
     "rows[3]/by_scheme/formal_censored/old/sd_ret_med", 100),
    ("paper.md", r"leaving [0-9.]+%, ([0-9.]+)% and [0-9.]+% of the prior standard deviation\. The",
     "step6_noise_sensitivity.json",
     "rows[3]/by_scheme/formal_censored/average/sd_ret_med", 100),
    ("paper.md", r"widened the three bounds by ([0-9.]+),", "step7c_ar1.json",
     "coef/old/widening", 1),
    ("paper.md", r"widened the three bounds by [0-9.]+, ([0-9.]+)", "step7c_ar1.json",
     "coef/average/widening", 1),
    ("paper.md", r"times, to ([0-9.]+)%, [0-9.]+% and [0-9.]+% of the prior", "step7c_ar1.json",
     "coef/old/crlb_over_prior_ar1", 100),
    ("paper.md", r"effective sample size at about ([0-9]+) of 294", "step7c_ar1.json", "n_eff", 1),
    # 3.3.2 systematic sensor error
    ("paper.md", r"deviations for old \(node 15, \+0\.10 mg L\^-1\^\), ([0-9.]+) for average",
     "step8c_bias_bynode.json", "summary/max_own_shift_over_sd/value", -1),
    ("paper.md", r"moved the old-zone estimate by ([0-9.]+) and", "step8_sensor_bias.json",
     "rows[5]/shift_over_sd", 1),
    # 3.3.3 bulk decay compensating for wall decay
    ("paper.md", r"coefficients by ([0-9.]+) to [0-9.]+ baseline posterior",
     "step8b_kb_sensitivity.json",
     "rows[2]/by_scheme/formal_censored/shift_over_own_sd/avg", 1),
    ("paper.md", r"coefficients by [0-9.]+ to ([0-9.]+) baseline posterior",
     "step8b_kb_sensitivity.json",
     "rows[0]/by_scheme/formal_censored/shift_over_own_sd/new", -1),
    ("paper.md", r"old coefficient by only ([0-9.]+) to", "step8b_kb_sensitivity.json",
     "rows[2]/by_scheme/formal_censored/shift_over_own_sd/old", 1),
    ("paper.md", r"raises the Cramér–Rao bound from [0-9.]+%, [0-9.]+% and [0-9.]+% of the "
     r"prior standard deviation to\s*\n?\s*[0-9.]+%, ([0-9.]+)%", "step7_fisher.json",
     "cases/B: kw+kb/coef/average/crlb_over_prior", 100),
    ("paper.md", r"raises the ratios to [0-9.]+%, ([0-9.]+)%", "step7_fisher.json",
     "cases/C: kw+kb+6 offsets/coef/average/crlb_over_prior", 100),
    # 3.4 structural heterogeneity
    ("paper.md", r"floor of ([0-9.]+) mg L\^-1\^ to", "step5c_jitter_sweep.json",
     "rows[0]/struct_residual", 1),
    ("paper.md", r"RMSE remained ([0-9.]+) mg L\^-1\^ against a realised", "step5d_structured.json",
     "rmse_min", 1),
    # 3.5 held-out prediction
    ("paper.md", r"RMSEs of\s*\n?\s*([0-9.]+), [0-9.]+ and [0-9.]+ mg L\^-1\^ and coverage",
     "step11_loo.json", "leave_one_zone_out[0]/pred_rmse[0]", 1),
    ("paper.md", r"The new coefficient retained about ([0-9]+)%", "step11_loo.json",
     "leave_one_zone_out[2]/own_sd_retained[0]", 100),
    # 3.6 risk propagation
    ("paper.md", r"median window minimum of ([0-9.]+) mg", "step10_risk_metrics.json",
     "top10[3]/minC_5_50_95[1]", 1),
    ("paper.md", r"duration is ([0-9.]+) h; over the 59", "step10_risk_metrics.json",
     "network_averages/by_metric/E_duration_h/unweighted_all_junctions", 1),
    ("paper.md", r"never falls below\s*\n?\s*([0-9.]+)\.", "step8c_bias_bynode.json",
     "summary/min_risk_spearman", 1),
    # An external check of the paper against the artifacts, done by hand, found a metric defined
    # one way in the baseline and another way in the perturbation arms. Every figure it had to
    # verify by hand to reach that point is registered below, so the same audit is a script run
    # next time rather than a person reading two files side by side.
    ("paper.md", r"realised noise RMSE over the window was ([0-9.]+)", "baseline_meta.json",
     "summary/noise_rmse", 1),
    ("paper.md", r"candidate library was ([0-9.]+) mg", "baseline_meta.json",
     "summary/rmse_min", 1),
    ("paper.md", r"gives standard deviations of ([0-9.]+), [0-9.]+ and [0-9.]+ m day",
     "step14_repeated_noise.json", "crlb_case_A_from_step7/old", 1),
    ("paper.md", r"gives standard deviations of [0-9.]+, ([0-9.]+) and [0-9.]+ m day",
     "step14_repeated_noise.json", "crlb_case_A_from_step7/average", 1),
    ("paper.md", r"gives standard deviations of [0-9.]+, [0-9.]+ and ([0-9.]+) m day",
     "step14_repeated_noise.json", "crlb_case_A_from_step7/new", 1),
    ("paper.md", r"from\s*\n?\s*0.222 to 0.365 and ([0-9.]+) mg", "step12_scenarios.json",
     "scenario_summary[2]/net_mean_E_deficit", 1),
    ("paper.md", r"duration from 5.22 to ([0-9.]+) and [0-9.]+ h", "step12_scenarios.json",
     "dosing_heatwave[1]/net_mean_E_duration_h", 1),
    ("paper.md", r"duration from 5.22 to [0-9.]+ and ([0-9.]+) h", "step12_scenarios.json",
     "dosing_heatwave[2]/net_mean_E_duration_h", 1),
    ("paper.md", r"the old coefficient by only [0-9.]+ to ([0-9.]+)", "step8b_kb_sensitivity.json",
     "rows[0]/by_scheme/formal_censored/shift_over_own_sd/old", -1),
    # the two counts the duration-metric conclusion rests on, and the arm total they are out of
    ("paper.md", r"nine of the ([0-9]+) arms alter", "step8c_bias_bynode.json",
     "summary/n_arms", 1),
    ("paper.md", r"deficit-based one is unchanged in ([0-9]+)", "step8c_bias_bynode.json",
     "summary/n_arms_holding_deficit_top6", 1),
    # middles of the reported triples. The first destructive test anchored first and third and
    # left the centre value free, so a corrupted middle passed; anchoring one element of a triple
    # does not anchor the triple.
    ("paper.md", r"coverage of [0-9.]+, ([0-9.]+) and [0-9.]+",
     "step14_repeated_noise.json",
     "by_scheme/formal_censored/coef/average/coverage/q90", 1),
    ("paper.md", r"posterior standard deviation was ([0-9.]+)%,", "baseline_meta.json",
     "summary/schemes/formal_censored/coef/old/sd_retained", 100),
    ("paper.md", r"left [0-9.]+%, [0-9.]+% and ([0-9.]+)% of the prior", "baseline_meta.json",
     "summary/schemes/informal_glue_draft_thr/coef/new/sd_retained", 100),
    ("paper.md", r"fell only to [0-9.]+%, [0-9.]+%\s*\n?\s*and ([0-9.]+)%", "step3_threshold.json",
     "rows[0]/new/sd_retained", 100),
    ("paper.md", r"recovered [0-9]+%, [0-9]+% and ([0-9]+)% of the imposed",
     "step4d_displaced_robust.json",
     "designs/DOWN/by_scheme/formal_censored/groups/new/gap_med", 100),
    ("paper.md", r"leaving [0-9.]+%, [0-9.]+% and ([0-9.]+)% of the prior standard deviation\. The",
     "step6_noise_sensitivity.json",
     "rows[3]/by_scheme/formal_censored/new/sd_ret_med", 100),
    ("paper.md", r"widened the three bounds by [0-9.]+, [0-9.]+ and\s*\n?\s*([0-9.]+)",
     "step7c_ar1.json", "coef/new/widening", 1),
    ("paper.md", r"equivalent to\s+[0-9.]+%, [0-9.]+% and ([0-9.]+)% of the prior",
     "step7_fisher.json", "cases/A: kw only/coef/new/crlb_over_prior", 100),
    ("paper.md", r"RMSEs of\s*\n?\s*[0-9.]+, ([0-9.]+) and [0-9.]+ mg L\^-1\^ and coverage",
     "step11_loo.json", "leave_one_zone_out[1]/pred_rmse[0]", 1),
    ("paper.md", r"RMSEs of\s*\n?\s*[0-9.]+, [0-9.]+ and ([0-9.]+) mg L\^-1\^ and coverage",
     "step11_loo.json", "leave_one_zone_out[2]/pred_rmse[0]", 1),
    # --- scenarios: the node counts the Discussion leans on ---
    ("README.md", r"baseline ([0-9]+) nodes", "step12_scenarios.json",
     "scenario_summary[0]/P_min_gt_0.5_nodes", 1),
    ("README.md", r"heat \+ ageing ([0-9]+)", "step12_scenarios.json",
     "scenario_summary[3]/P_min_gt_0.5_nodes", 1),
]

# (document, regex, why it must not appear). These are the specific wordings a previous revision left
# behind: a superseded configuration, a claim stronger than its evidence, or a label naming the wrong
# inference rule. A number check cannot catch any of them, because each is a sentence about numbers
# that are themselves current.
FORBIDDEN = [
    (r"(?i)\b2000[- ](draw|sample|member|prior draw)", "the design is 8192 scrambled Sobol draws"),
    (r"(?i)2000 EPANET runs", "the baseline is 8192 EPANET runs"),
    (r"(?i)\b24 ?h warm-?up\b(?![^\n]{0,60}(draft|superseded|earlier|was ))",
     "the warm-up is 120 h; 24 h is the superseded draft value and must be marked as such"),
    (r"(?i)72 ?h simulation(?![^\n]{0,60}(draft|superseded|earlier))",
     "the horizon is 168 h; 72 h is the superseded draft value"),
    (r"(?i)≈\s*residence-weighted|residence-weighted[^\n]{0,20}\(i\.e\.|"
     r"recovers the residence-weighted",
     "length-weighted is an illustrative proxy, NOT the residence-weighted coefficient (Step 5d)"),
    (r"(?i)GLUE behavioural mean of|"
     r"(?<!informal )GLUE (behavioural )?(mean|posterior)(?![^\n]{0,40}comparator)",
     "name the rule: the primary results are formal censored posterior means, not GLUE means"),
    (r"(?i)p\s*[≈=]\s*1e-16|p\s*<\s*1e-1[0-9]",
     "no p-value is reported for the age-risk association; the junctions are not independent"),
    (r"(?i)(proven|demonstrably|shown to be) efficient|estimator is efficient|"
     r"(reached|reaches|attained|attains) the (Cramér–Rao|Cramer-Rao|CRLB)",
     "the formal posterior spread is LOCALLY CONSISTENT with the CRLB; efficiency needs the "
     "repeated-sampling test in Step 14"),
    (r"(?i)(complete|full) periodic steady state|water age had converged|has converged after 120",
     "120 h is a finite-horizon pragmatic choice; the deficit and water-age criteria never pass"),
    # The candidate RMSE is computed on the 49x6 assessment window, so 0.0960 (the full-record
    # realised noise RMSE) is not its comparator; 0.0973 is. Both numbers are real and live in the
    # artifacts, so a drift check cannot see this — only the pairing is wrong.
    (r"(?i)noise floor[^\n]{0,45}0\.0960|0\.0960[^\n]{0,45}noise floor|"
     r"noise floor of 0\.096\b",
     "0.0960 is the FULL-RECORD realised noise RMSE; the comparator for a window RMSE is 0.0973 "
     "(baseline_meta summary/noise_rmse). Name which of sigma / window / full-record is meant"),
    (r"(?i)all numbers are checked|every number is verified|7/7 checks passing",
     "the validator covers a subset of numbers and establishes no semantic consistency"),
    # Step 8b: two readings that were written, tested against the artifact, and refuted. Neither can
    # be caught by a number check, because both are sentences ABOUT numbers that are themselves
    # current. A line quoting one of them as superseded is exempt via the SUPERSEDED marker.
    (r"(?i)risk (product|map|ranking)[^\n]{0,30}insensitive to[^\n]{0,25}k_b|"
     r"top nodes[^\n]{0,40}unchanged at every k_b|"
     # the same claim in list form: k_b enumerated among the perturbations the ranking survives
     r"\bthreshold,\s*`?k_b`?\s*±\s*20",
     "k_b ±20% retains only four of the six hot-spots (Jaccard 0.50); what stays correlated is the "
     "full-network ordering, not the operational shortlist (Step 8b)"),
    # "baseline GLUE" named a configuration and an inference rule at once. Once the primary rule
    # became the formal likelihood the two meanings came apart, and the phrase silently asserted that
    # the headline ensemble is the informal one. Name the rule, or say "the baseline model".
    (r"(?i)\bbaseline GLUE\b",
     "ambiguous: say 'the formal likelihood-weighted ensemble', 'the informal GLUE behavioural "
     "ensemble', or 'the baseline model' if the configuration is what is meant"),
    # Narrowed after the two risk metrics were computed separately. "The turnover is a boundary
    # effect" is FALSE on P_bar (Jaccard worse at k=3 than k=6; reference rank 4 falls to 20th) and
    # TRUE on E[A] (at -0.6 the only change is a 6th/7th swap, nothing changes at -0.4). So the
    # forbidden thing is not the claim, it is making it WITHOUT NAMING THE METRIC.
    (r"(?i)(instability|turnover)[^\n]{0,40}(concentrated around|confined to)[^\n]{0,30}"
     r"(prioritisation|priority) boundary(?![^\n]{0,80}(P_bar|E\[A\]|deficit|metric))",
     "name the metric: this is true on the cumulative-deficit ranking E[A] and false on the "
     "time-averaged P_bar, where the reordering reaches rank 20 (Step 8b findings 4-5)"),
]
# Documents where FORBIDDEN is enforced. The response matrix is included because that is where the
# over-claims lived.
FORBIDDEN_DOCS = ["RESULTS_LOG.md", "README.md", "README.en.md", "REVISION_RESPONSE_MATRIX.md"]

# check_forbidden only. A sentence that REJECTS a wording necessarily contains that wording, so the
# register would flag the very lines that do the rejecting. These markers identify a rejection and
# exempt the line. Kept separate from SUPERSEDED, which also exempts numbers from the value check and
# must not be widened for a text-only purpose.
REJECTS = re.compile(r"(?i)(\bdo not write\b|\bmust not be\b|\bmust \*\*not\*\* be\b|"
                     r"\bthe earlier response said\b|\banywhere the draft says\b|"
                     r"\bthat sentence has to go\b|\bwas tested and\b|\brefuted\b|"
                     r"\bwas wrong\b|\bis false\b|\bis not used\b|"
                     r"\bnot sayable\b|\bwithout naming the metric\b)")

# Numbers that legitimately appear as prose constants anywhere in the log.
CONSTANTS = {
    0.02, 0.025, 0.05, 0.1, 0.107, 0.11, 0.12, 0.15, 0.2, 0.4, 0.5, 0.6, 0.9, 0.95, 1.0, 1.645,
    1.92, 1.96, 2.0, 6.0, 12.0, 16.0, 20.0, 24.0, 30.0, 48.0, 49.0, 72.0, 92.0, 100.0, 120.0,
    168.0, 294.0, 8192.0,
    # Algebraic consequences of N and sigma rather than outputs of any script, so no artifact holds
    # them and the drift check would otherwise match them against an unrelated value by proximity.
    17.15,      # sqrt(N) = sqrt(294), the factor by which the informal score inflates sigma
    1.71,       # sigma * sqrt(N) = the informal score's effective observation sd, mg/L
    0.0041,     # sigma / sqrt(2N), the sampling sd of the RMSE objective at the truth
}
# Quantile/percentile labels ("5-95% band", "25/75 IQR") read as bare integers.
PERCENTILES = {1, 5, 10, 25, 50, 68, 75, 90, 95, 99, 100}
# Junction names in Net3 are numeric strings, so they collide with measurements unless excluded.
NODE_IDS = {float(n) for n in provenance.B.all_junctions()}
CITATION_YEARS = range(1900, 2101)

# (regex, why it is wrong). Kept narrow on purpose: a broad unit scan is all false positives.
UNIT_RULES = [
    (re.compile(r"(?i)(coefficient|k_?w).{0,80}SD.{0,20}\(mg/L\)"),
     "wall coefficients are m/day, not mg/L"),
    (re.compile(r"(?i)k_?b.{0,40}\(mg/L\)"),
     "the bulk coefficient is per day, not mg/L"),
]
PATH_RULES = [
    (re.compile(r"/opt/anaconda3|/Users/|C:\\\\"), "machine-specific absolute path"),
]

HEADING = re.compile(r"^#{2,4}\s+(Step\s+[\w.]+?)\s*(?:—|-|–|$)", re.M)
# The second lookbehind drops numbers that are part of a hyphenated identifier (SHA-256, UTF-8,
# AR-1), which are names rather than measurements.
NUMBER = re.compile(r"(?<![\w.])(?<!\w-)(-?\d+(?:\.\d+)?)(?![\w.]*\d)")
# Lines that deliberately quote a superseded run cannot be checked: no current artifact holds those
# numbers, and that is the point of the sentence. The phrasing is the opt-out.
SUPERSEDED = re.compile(r"(?i)\b(an? earlier (run|version)|draft'?s earlier|earlier single-window|"
                        r"superseded|previously reported|before the fix)\b")


# ---------------------------------------------------------------- helpers
def walk_numbers(obj):
    """Every numeric leaf in a JSON structure, including numbers embedded in strings."""
    if isinstance(obj, bool):
        return
    if isinstance(obj, (int, float)):
        yield float(obj)
    elif isinstance(obj, dict):
        for k, v in obj.items():
            for n in walk_numbers(k):
                yield n
            for n in walk_numbers(v):
                yield n
    elif isinstance(obj, (list, tuple)):
        for v in obj:
            for n in walk_numbers(v):
                yield n
    elif isinstance(obj, str):
        for m in NUMBER.finditer(obj):
            try:
                yield float(m.group(1))
            except ValueError:
                pass


def value_set(paths):
    """Sorted absolute values of every numeric leaf in the given artifacts.

    Matching is done by interval rather than by equality of roundings: a log number written to d
    decimals is consistent with any artifact value within half a unit in the last place. Comparing
    rounded values instead misfires on ties (1755/2000 = 0.8775 prints as 87.8% but rounds to 0.877).
    Signs are dropped because the log writes negatives with U+2212, which does not parse.
    """
    vals = set()
    for p in paths:
        full = os.path.join(CACHE, p)
        if not os.path.exists(full):
            continue
        with open(full) as f:
            data = json.load(f)
        vals.update(abs(v) for v in walk_numbers(data))
    return sorted(vals)


# A short string that must appear in a document for it to be the one the claims were written
# against. README.md is the case that needs it: the published repository ships its own README at
# the same path, and without a fingerprint the six claims anchored in the project README are
# reported as broken rather than as uncovered.
DOC_FINGERPRINTS = {
    "README.md": "分组一阶余氯壁衰减",
    "README.en.md": "Uncertainty-aware calibration",
    "RESULTS_LOG.md": "## Step 1",
}


def _missing_docs():
    """Documents named in DOCS that are absent, or present but not the expected document.

    The published repository ships the code, the artifacts and RESULTS_LOG.md, but not the
    manuscript, the supervisor correspondence, or the project README. Rather than crash on the
    first absent file, or fail a claim whose document simply is not here, the run reports which
    checks it could not perform, so a passing result is never read as covering more than it did.
    """
    out = []
    for name, path in DOCS.items():
        if not os.path.exists(path):
            out.append(name)
            continue
        mark = DOC_FINGERPRINTS.get(name)
        if mark and mark not in open(path, encoding="utf8", errors="ignore").read():
            out.append(name)
    return sorted(out)


def log_sections():
    """[(label, body, first_line_no)] for each Step heading in the log, in order."""
    text = open(LOG).read()
    marks = [(m.start(), m.group(1).strip()) for m in HEADING.finditer(text)]
    out = []
    for i, (pos, label) in enumerate(marks):
        end = marks[i + 1][0] if i + 1 < len(marks) else len(text)
        line_no = text.count("\n", 0, pos) + 1
        out.append((re.sub(r"\s+", " ", label), text[pos:end], line_no))
    return out, text


def num_f(num_str):
    return float(num_str)


def is_excluded(num_str):
    """Numbers that cannot be compared against artifacts: labels, identifiers, constants."""
    v = float(num_str)
    if abs(v) in CONSTANTS:
        return True
    if "." not in num_str:                       # bare integers carry the ambiguous cases
        iv = int(v)
        return (abs(v) < MIN_INT_TO_CHECK or iv in PERCENTILES or iv in CITATION_YEARS
                or abs(v) in NODE_IDS)
    return False


def _dp(num_str):
    return len(num_str.split(".")[1]) if "." in num_str else 0


def candidates(num_str):
    """(|value|, half-width) intervals a log number may occupy relative to the stored value.

    Two rescalings are routine in the log and neither is an error: fractions are stored (0.6611) and
    quoted as percentages (66.1%), and SI quantities are stored in base units but quoted with a
    prefix (activation energies in J/mol, quoted in kJ/mol). Each form carries its own half-width,
    so a percentage is compared at the precision it implies for the fraction, not at its own.
    """
    v, dp = abs(float(num_str)), _dp(num_str)
    half = 0.5 * 10 ** -dp
    return [(v, half),
            (v / 100.0, half / 100.0), (v * 100.0, half * 100.0),
            (v / 1000.0, half / 1000.0), (v * 1000.0, half * 1000.0)]


def matches(num_str, pools):
    """True if some artifact value could have been rounded to the number as written.

    The interval is widened by a relative epsilon because rescaling introduces representation
    error: 87.8/100 is 0.8780000000000001, whose exact lower bound would exclude a stored 0.8775.
    """
    for v, half in candidates(num_str):
        half += abs(v) * 1e-9 + 1e-12
        for vals in pools:
            i = bisect.bisect_left(vals, v - half)
            if i < len(vals) and vals[i] <= v + half:
                return True
    return False


def nearest(num_str, pools):
    """Nearest artifact value to the number as written, ignoring sign (the log uses U+2212)."""
    v = abs(float(num_str))
    best = None
    for vals in pools:
        for x in vals:
            if best is None or abs(x - v) < abs(best - v):
                best = x
    return best


# ---------------------------------------------------------------- checks
def check_artifacts_exist():
    missing = sorted({p for ps in SECTIONS.values() for p in ps
                      if not os.path.exists(os.path.join(CACHE, p))})
    fig_missing = sorted(f for f in FIGURE_SOURCE if not os.path.exists(os.path.join(FIGDIR, f)))
    problems = [f"missing artifact: baseline_cache/{p}" for p in missing]
    problems += [f"missing figure: {os.path.basename(FIGDIR)}/{f}" for f in fig_missing]
    referenced = {p for ps in SECTIONS.values() for p in ps} | {"cache_manifest.json"}
    orphans = sorted(f for f in os.listdir(CACHE)
                     if f.endswith(".json") and f not in referenced
                     and not f.endswith(".npy.key.json"))   # provenance sidecars, not results
    notes = [f"artifact not mapped to any log section: {o}" for o in orphans]
    return problems, notes


def check_manifest():
    bad = provenance.check_manifest()
    return [f"manifest field {k!r}: cache={o!r} now={n!r}" for k, o, n in bad], []


def check_release_provenance(release=False):
    """Is the manifest a citable release record, or a development snapshot?

    The manifest is written from whatever tree it is run on. The working habit here has been
    edit -> provenance.py -> commit everything including the manifest, which guarantees the manifest
    records a DIRTY tree at the PREVIOUS commit and never describes the commit it ships in. That is
    fine while developing and useless for citation: `uncommitted_diff_sha256` identifies the diff
    without storing it, so nothing outside this working copy can reconstruct what ran.

    Normally this reports a note. Under --release it fails, so the release sequence has to be
    followed deliberately: commit everything -> provenance.py on a clean tree -> commit the manifest
    alone -> this check -> tag.
    """
    path = os.path.join(CACHE, "cache_manifest.json")
    MANIFEST_REL = "Net3/baseline_cache/cache_manifest.json"
    if not os.path.exists(path):
        return ["cache_manifest.json missing"], []
    with open(path) as f:
        git = json.load(f).get("git", {})
    head = subprocess.run(["git", "-C", os.path.dirname(HERE), "rev-parse", "HEAD"],
                          capture_output=True, text=True)
    head = head.stdout.strip() if head.returncode == 0 else None
    recorded = git.get("commit")
    problems, notes = [], []
    msgs = []
    if git.get("dirty"):
        n = len(git.get("dirty_files") or [])
        msgs.append(f"manifest was written on a DIRTY tree ({n} modified file(s)); the diff is "
                    f"identified by hash, not stored")
    # The manifest cannot record the commit that contains it: writing it dirties the tree, and
    # committing it advances HEAD. Requiring recorded == HEAD makes the release state unreachable.
    # The rule that is both correct and attainable: HEAD may differ from the recorded commit ONLY by
    # the manifest file itself, so every hash the manifest asserts about the code still holds.
    if recorded and head and recorded != head:
        diff = subprocess.run(["git", "-C", os.path.dirname(HERE), "diff", "--name-only",
                               recorded, head], capture_output=True, text=True)
        changed = [f for f in diff.stdout.split() if f] if diff.returncode == 0 else None
        manifest_only = changed is not None and set(changed) <= {MANIFEST_REL}
        if changed is None:
            msgs.append(f"manifest records commit {recorded[:12]}, HEAD is {head[:12]}, and the "
                        f"difference could not be read from git")
        elif not manifest_only:
            others = [f for f in changed if f != MANIFEST_REL]
            msgs.append(f"manifest records commit {recorded[:12]}, HEAD is {head[:12]}, and they "
                        f"differ in {len(others)} file(s) besides the manifest "
                        f"({', '.join(others[:3])}{'…' if len(others) > 3 else ''})")
    if not msgs:
        how = "= HEAD" if recorded == head else "HEAD differs only by the manifest itself"
        notes.append(f"release-grade: manifest records a clean tree at {(recorded or '?')[:12]} "
                     f"({how})")
        return [], notes
    msgs.append("release sequence: commit everything -> python provenance.py on a clean tree -> "
                "commit the manifest alone -> python validate_artifacts.py --release -> tag")
    return (msgs, []) if release else ([], msgs)


def check_figure_freshness():
    """A figure older than the artifact it plots is stale and must be regenerated.

    Caveat worth knowing before acting on a failure here: this compares mtimes, and **git does not
    preserve mtimes**. Checking out or merging a branch rewrites every file whose content differs and
    leaves the byte-identical ones alone, so a JSON can end up newer than its figure without either
    being stale. The remedy is still to re-run the step — never to `touch` the figure, which would
    hide a real staleness the next time.
    """
    lags = {}
    for fig, src in FIGURE_SOURCE.items():
        fp, sp = os.path.join(FIGDIR, fig), os.path.join(CACHE, src)
        if not (os.path.exists(fp) and os.path.exists(sp)):
            continue
        lags[fig] = (os.path.getmtime(sp) - os.path.getmtime(fp), src)
    late = {f: v for f, v in lags.items() if v[0] > FIGURE_STALE_S}
    if not late:
        return [], []
    # A clone or a bulk copy rewrites every artifact at once, which makes every figure look stale
    # by the same amount. Real staleness is selective: one step was re-run and its figure was not.
    # Reporting the first case as a failure would mean this check can never pass in a published
    # checkout, which is where a reader would run it.
    spread = max(v[0] for v in late.values()) - min(v[0] for v in late.values())
    if len(late) == len(lags) and spread < 600:
        return [], [f"all {len(lags)} figures trail their artifacts by a near-constant "
                    f"{min(v[0] for v in late.values()) / 60:.0f} min, which is the signature of a "
                    f"checkout or a copy rather than of staleness; mtimes are not preserved by git, "
                    f"so freshness is unverifiable here"]
    return [f"{f} is {v[0] / 60:.1f} min older than {v[1]} — regenerate it"
            for f, v in sorted(late.items())], []


# A NUMERICAL GUARD is a constant that silently changes a result when it binds: a clip, a floor, a
# cap, a truncation. Its danger is asymmetric — while inert it is invisible and harmless, and the
# moment it binds it alters the answer with no other symptom. Step 12's CLIP_LO was active for 46.8%
# of members before anyone noticed; step12's WEIGHT_FLOOR and step4d's UPPER_CAP each shipped without
# their activation in the artifact. So every guard must record HOW MUCH IT ACTED, in the artifact,
# not on stdout.
#   script -> {constant: substring that must appear in one of the artifact's keys}
GUARDS = {
    "step12_scenarios.py": {"CLIP_LO": "clipped_draws",
                            "WEIGHT_FLOOR": "discarded_weight_mass",
                            "EA_MIN": "truncated_at"},
    "step4d_displaced_robust.py": {"UPPER_CAP": "upper_cap_shift_applied"},
    "step6_noise_sensitivity.py": {"ESS_MIN": "sampling_limited"},
    "step13_known_answer.py": {"TOL_REL": "tolerance"},
}
# Constants matching this but absent from GUARDS are reported: a new guard must be registered, which
# is the only way the check can cover code that does not exist yet.
GUARD_NAME = re.compile(r"^(?!.*(?:C_MIN|THRESHOLDS|REPORT_STEPS|GRID_FLOOR|NOISE_FLOOR))"
                        r"[A-Z][A-Z0-9_]*(?:FLOOR|CLIP|CAP|GUARD|LIMIT|TOL|TOL_REL|_MIN|_MAX)\b")


def check_guards():
    """Every numerical guard must record its activation in the artifact it feeds."""
    problems, notes = [], []
    for script, entries in GUARDS.items():
        arts = [a for label, ps in SECTIONS.items() for a in ps
                if label.replace("Step ", "step").replace(" ", "") == script.split("_")[0]]
        pool = ""
        for a in set(arts):
            full = os.path.join(CACHE, a)
            if os.path.exists(full):
                pool += open(full).read()
        if not pool:
            notes.append(f"{script}: no artifact found to check its guards against")
            continue
        for const, key in entries.items():
            if key not in pool:
                problems.append(f"{script}: guard {const} does not record its activation "
                                f"(no key containing {key!r} in the artifact) — an inert guard and a "
                                f"binding one are indistinguishable from the record")
    # discovery: a guard-shaped constant that nobody registered
    for f in sorted(os.listdir(HERE)):
        if not re.fullmatch(r"step\w+\.py", f):
            continue
        for line in open(os.path.join(HERE, f)):
            m = re.match(r"([A-Z][A-Z0-9_]*)\s*=", line)
            if m and GUARD_NAME.match(m.group(1)) and m.group(1) not in GUARDS.get(f, {}):
                notes.append(f"{f}: {m.group(1)} looks like a guard but is not in GUARDS — register "
                             f"it with the artifact key that records its activation, or confirm it "
                             f"cannot bind")
    return problems, notes


def check_reproduce_list():
    """Every step script must appear in the log's run block.

    The list has drifted twice: a step is added, its row goes into the file table, and the command
    that actually runs it is forgotten — so the log documents an experiment that the stated procedure
    never performs. Checking it is one set difference.
    """
    text = read_doc("RESULTS_LOG.md")
    if text is None:
        return ["RESULTS_LOG.md missing"], []
    scripts = {f for f in os.listdir(HERE) if re.fullmatch(r"step\w+\.py", f)}
    missing = sorted(s for s in scripts if f"python {s}" not in text)
    return [f"{s} is never run by the documented procedure" for s in missing], []


def check_numbering(sections):
    problems = []
    scripts = os.listdir(HERE)
    for label, _, _ in sections:
        if label not in SECTIONS:
            problems.append(f"log heading {label!r} is not a known section "
                            f"(known: {', '.join(sorted(SECTIONS))})")
            continue
        stem = label.replace("Step ", "step").replace(" ", "")
        # sub-sections that a single script produces rather than having a file of their own
        if not any(s.startswith(stem + "_") for s in scripts) and stem not in (
                "step2", "step5a", "step5e", "step11b", "step11c",
                "step5c-ensemble", "step5d-dose", "step11d"):
            problems.append(f"log heading {label!r} has no matching {stem}_*.py script")
    return problems, []


def check_log_numbers(sections, verbose=False):
    """Compare each section's numbers with its own artifacts, then with all artifacts.

    The global fallback matters because the log legitimately cross-references numbers between
    sections (a prior SD quoted in Step 2 lives in the Step 7 artifact). Only a number that matches
    NOWHERE, yet sits within STALE_REL_TOL of some artifact value, is reported as stale.
    """
    problems, notes = [], []
    global_pool = value_set(sorted({p for ps in SECTIONS.values() for p in ps}))
    total_ok = total_un = total_xref = 0
    for label, body, first_line in sections:
        if label not in SECTIONS:
            continue
        own_pool = value_set(SECTIONS[label])
        if not own_pool:
            notes.append(f"{label}: no artifact values available, numbers unchecked")
            continue
        stale, unmatched, exempt, global_only = [], [], 0, []
        # The exemption is scoped to the markdown paragraph, not the line: prose wraps, so the
        # marker phrase and the superseded numbers it introduces often land on different lines.
        para_exempt = False
        for offset, line in enumerate(body.splitlines()):
            line_no = first_line + offset
            if not line.strip():
                para_exempt = False
                continue
            if line.lstrip().startswith("#"):     # headings carry sub-section numbers, not results
                continue
            if SUPERSEDED.search(line):
                para_exempt = True
            if para_exempt:
                exempt += 1
                continue
            # A markdown table row is where results are reported; prose is where other sections get
            # cross-referenced. Table cells are therefore held to the section's OWN artifact — the
            # global pool let a whole stale results table survive a change of weighting scheme on
            # coincidental hits. The rule applies only from three decimals up: with thousands of
            # values in the pool, a two-decimal cell matches somewhere no matter what, so demanding
            # provenance for it would be a guess dressed up as a check.
            in_table = line.lstrip().startswith("|")
            for m in NUMBER.finditer(line):
                num = m.group(1)
                if is_excluded(num):
                    continue
                if matches(num, [own_pool]):
                    total_ok += 1
                    continue
                strict = in_table and _dp(num) >= 3
                if not strict and matches(num, [global_pool]):
                    # Legal in prose, but with thousands of values a coincidental hit is possible,
                    # so these are surfaced for eyeballing rather than silently counted as verified.
                    total_ok += 1
                    global_only.append(num)
                    continue
                near = nearest(num, [own_pool, global_pool])
                if near is not None and 0 < abs(near - num_f(num)) / max(abs(near), 1e-12) \
                        <= STALE_REL_TOL:
                    stale.append((num, near, line_no))
                else:
                    unmatched.append(num)
        total_un += len(unmatched)
        for got, near, line_no in stale:
            problems.append(f"RESULTS_LOG.md:{line_no} ({label}): log says {got}; "
                            f"no artifact has it, nearest artifact value is {near:g}")
        if verbose and unmatched:
            notes.append(f"{label}: {len(unmatched)} unmatched (likely prose): "
                         f"{', '.join(sorted(set(unmatched))[:20])}")
        if verbose and exempt:
            notes.append(f"{label}: {exempt} line(s) exempt as describing a superseded run")
        if verbose and global_only:
            notes.append(f"{label}: {len(global_only)} number(s) match another section's artifact "
                         f"but not this one — check they are intended cross-references: "
                         f"{', '.join(sorted(set(global_only))[:12])}")
        total_xref += len(global_only)
    notes.append(f"numbers matched to artifacts: {total_ok} (of which {total_xref} matched only "
                 f"another section's artifact; --verbose lists them); unmatched (likely prose): "
                 f"{total_un}")
    return problems, notes


def check_text_rules(text):
    problems = []
    for i, line in enumerate(text.splitlines(), 1):
        for rx, why in UNIT_RULES + PATH_RULES:
            if rx.search(line):
                problems.append(f"RESULTS_LOG.md:{i}: {why} -> {line.strip()[:90]}")
    return problems, []


# NOT IMPLEMENTED, deliberately: a number check on README.md.
#
# The README summarises every section, so it has no per-section artifact mapping and its comparison
# pool is the union of all artifacts — thousands of values at every rounding. A check against that
# pool cannot fail: injecting three realistic regressions (SD retained 86/98/98 -> 84/97/98, a CRLB
# ratio 1.06 -> 1.11, and the warm-up 120 h -> 96 h) was caught ZERO times out of three. A check in
# this list that cannot fail is worse than no check, because it turns an unverified document into a
# green light, so it was removed rather than shipped.
#
# The workable version is an explicit claim registry — each README figure mapped to a JSON path in a
# named artifact, roughly 40 entries — which is precise but is real work. Until that exists, the
# README's numbers are verified only by whoever transcribed them, i.e. exactly the position
# RESULTS_LOG.md was in before this tooling (93 stale numbers).


_PATH_STEP = re.compile(r"([^/\[\]]+)|\[(\d+)\]")


def json_at(data, path):
    """Value at a slash-separated path with optional [i] indices.

    Slash rather than dot because several JSON keys contain a dot ("P_min_gt_0.5_nodes"), and a
    separator that appears inside key names silently mis-parses instead of failing.
    """
    cur = data
    for key, idx in _PATH_STEP.findall(path):
        cur = cur[int(idx)] if idx else cur[key]
    return cur


def read_doc(name):
    p = DOCS[name]
    if not os.path.exists(p):
        return None
    text = open(p, encoding="utf8", errors="ignore").read()
    mark = DOC_FINGERPRINTS.get(name)
    return None if (mark and mark not in text) else text


def check_claims():
    """Every registered claim: the anchor must exist and the number must equal its JSON source."""
    problems, notes = [], []
    loaded, checked = {}, 0
    uncovered = set()
    for doc, pattern, artifact, path, scale in CLAIMS:
        text = read_doc(doc)
        if text is None and doc in _missing_docs():
            uncovered.add(doc)
            continue
        if text is not None:
            # Collapse whitespace before matching. An anchor is a sentence, and a sentence wraps
            # wherever the line happens to end, so a literal space in a pattern would otherwise
            # stop matching the moment an edit reflowed the paragraph. That has now silently
            # unanchored a claim twice, once on a re-anchoring and once on a style pass.
            text = re.sub(r"\s+", " ", text)
        if text is None:
            problems.append(f"claim document missing: {doc}")
            continue
        if artifact not in loaded:
            full = os.path.join(CACHE, artifact)
            loaded[artifact] = json.load(open(full)) if os.path.exists(full) else None
        data = loaded[artifact]
        if data is None:
            problems.append(f"claim source missing: baseline_cache/{artifact}")
            continue
        try:
            want = float(json_at(data, path)) * scale
        except (KeyError, IndexError, TypeError, ValueError):
            problems.append(f"claim path not found in {artifact}: {path}")
            continue
        found = re.findall(pattern, text)
        if not found:
            problems.append(f"{doc}: registered claim has no anchor — the sentence matching "
                            f"/{pattern}/ is gone, so {artifact}:{path} is no longer checked "
                            f"anywhere. Restore the wording or update CLAIMS.")
            continue
        for got in found:
            got = got if isinstance(got, str) else got[0]
            got = got.strip().rstrip(".,;:")
            half = 0.5 * 10 ** -_dp(got) + abs(want) * 1e-9 + 1e-12
            if abs(float(got) - want) > half:
                problems.append(f"{doc}: claims {got} but {artifact}:{path} is {want:g}")
            else:
                checked += 1
    notes.append(f"{checked} registered claim(s) verified against a named JSON path "
                 f"(of {len(CLAIMS)} registered)")
    if uncovered:
        notes.append("claims in these documents were not checked, the document is not in this "
                     "checkout: " + ", ".join(sorted(uncovered)))
    return problems, notes


def check_weighting_declared():
    """Every result artifact must name the weighting that produced it."""
    problems, notes = [], []
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".json") or fn.endswith(".npy.key.json") or fn in WEIGHTING_EXEMPT:
            continue
        with open(os.path.join(CACHE, fn)) as f:
            data = json.load(f)
        if not isinstance(data, dict):
            problems.append(f"{fn}: top level is a {type(data).__name__}, so it cannot carry a "
                            f"weighting declaration — wrap it in an object")
            continue
        if not any(k in data for k in WEIGHTING_KEYS):
            problems.append(f"{fn}: no {' or '.join(WEIGHTING_KEYS)} field — a table of weighted "
                            f"numbers whose weighting is only implied by the script that wrote it")
    notes.append(f"{len(WEIGHTING_EXEMPT)} artifact(s) exempt (no weighted ensemble): "
                 f"{', '.join(sorted(WEIGHTING_EXEMPT))}")
    return problems, notes


def check_forbidden():
    """Wordings that a previous revision left behind and that no number check can catch."""
    problems = []
    for doc in FORBIDDEN_DOCS:
        text = read_doc(doc)
        if text is None:
            continue
        for i, line in enumerate(text.splitlines(), 1):
            if SUPERSEDED.search(line) or REJECTS.search(line):
                continue
            for rx, why in FORBIDDEN:
                if re.search(rx, line):
                    problems.append(f"{doc}:{i}: {why}\n           -> {line.strip()[:110]}")
    return problems, []


def check_env_claims(text):
    problems = []
    with open(os.path.join(CACHE, "cache_manifest.json")) as f:
        man = json.load(f)
    m = re.search(r"conda env `([\w.-]+)`", text)
    if m and m.group(1) != man["conda_env"]:
        problems.append(f"log claims conda env `{m.group(1)}` but the cache was produced in "
                        f"`{man['conda_env']}`")
    for pkg in ("numpy", "wntr", "scipy"):
        for got in re.findall(rf"(?i){pkg}\s+([0-9]+\.[0-9]+(?:\.[0-9]+)?)", text):
            if not man[pkg].startswith(got) and not got.startswith(man[pkg]):
                problems.append(f"log claims {pkg} {got} but the cache was produced with "
                                f"{pkg} {man[pkg]}")
    return problems, []


# ---------------------------------------------------------------- driver
def _paper_section_artifacts():
    """Paper section -> the artifacts that section may quote, read from Appendix I.

    The paper already publishes this mapping for the reader; using the same table as the check's
    binding means the two cannot silently disagree.
    """
    doc = read_doc("paper.md") or ""
    m = re.search(r"^\| Paper section \| Script \|$(.*?)(?=\n[ \t]*\n|\Z)", doc, re.M | re.S)
    if not m:
        return {}
    out = {}
    for line in m.group(1).splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) != 2 or set(cells[0]) <= set("-: "):
            continue
        arts = []
        for script in re.findall(r"`([a-z0-9_]+)\.py`", cells[1]):
            for cand in (script + ".json", EXTRA_ARTIFACTS.get(script)):
                if cand and os.path.exists(os.path.join(CACHE, cand)):
                    arts.append(cand)
        for sec in re.findall(r"\d+\.\d+(?:\.\d+)?", cells[0]):
            out.setdefault(sec, set()).update(arts)
    return out


EXTRA_ARTIFACTS = {"step1_freeze_baseline": "baseline_meta.json",
                   "step3_threshold_sensitivity": "step3_threshold.json"}


def check_paper_numbers(verbose=False):
    """Hold each paper section's numbers to the artifacts Appendix I says produced it.

    A global pool over every artifact was tried first and thrown away: with 30,000 pooled values a
    random two-decimal figure matches by coincidence 93% of the time, so it passed a deliberately
    corrupted number and would have given false confidence. Binding a section to its own artifacts
    is the same discipline check_log_numbers applies to RESULTS_LOG, and it is what gives either
    check teeth.

    Sections Appendix I does not map are not checked, and the note says how many those are. This
    still does not prove provenance; it detects drift in the numbers it covers.
    """
    doc = read_doc("paper.md")
    if doc is None:
        # not a failure in a checkout that does not ship the manuscript; say so and move on
        return [], ["paper.md is not in this checkout, so its numbers were not checked"]
    mapping = _paper_section_artifacts()
    if not mapping:
        return ["Appendix I section-to-script table not found; paper numbers unchecked"], []

    body = re.sub(r"<!--.*?-->", "", doc, flags=re.S)
    body = re.split(r"^## Appendix A\.", body, flags=re.M)[0]
    body = re.sub(r"\A---\n.*?\n---\n", "", body, flags=re.S)

    pools = {}
    def pool_for(arts):
        key = tuple(sorted(arts))
        if key in pools:
            return pools[key]
        pool = set()
        def add(v):
            # store the raw magnitude only. Pre-expanding each value over seven rounding levels
            # and three scalings was tried and abandoned: it inflated the pool until three
            # deliberately corrupted numbers all passed. Matching happens at the precision the
            # paper actually writes, which is where the discrimination is.
            try:
                f = float(v)
            except (TypeError, ValueError):
                return
            if f == f:
                pool.add(abs(f))
        def walk(o):
            if isinstance(o, dict):
                for v in o.values():
                    walk(v)
            elif isinstance(o, (list, tuple)):
                for v in o:
                    walk(v)
            elif isinstance(o, bool):
                pass
            elif isinstance(o, (int, float)):
                add(o)
            elif isinstance(o, str):
                for t in re.findall(r"-?\d+\.?\d*", o):
                    add(t)
        for a in arts:
            try:
                walk(json.load(open(os.path.join(CACHE, a))))
            except Exception:
                pass
        # The generated tables belong in every pool. Ratios such as a displacement divided by the
        # repeated-baseline SD exist in no single artifact; paper_figs.py computes them, and the
        # property worth checking is that the prose agrees with what it computed.
        # tables.md only. appendix_tables.md carries the 92-row risk register, and pooling ~8000
        # values from it let two deliberately corrupted numbers match by coincidence.
        f = os.path.join(FIGDIR, "paper", "tables.md")
        if os.path.exists(f):
            for t in re.findall(r"-?\d+\.?\d*", open(f).read()):
                add(t)
        # index by written precision so a lookup is a set hit rather than a scan of the pool
        idx = {}
        for d in range(0, 7):
            idx[d] = {round(x, d) for x in pool} | {round(x * 100, d) for x in pool}
        pools[key] = idx
        return idx

    SKIP = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "12", "20", "24",
            "30", "48", "49", "92", "100", "120", "168", "294", "0.5"}
    # Numbers the paper derives in the text rather than reading from an artifact. Listing them
    # with a reason is the point: an unexplained exemption is how a wrong number survives.
    DERIVED = {
        "0.1068": "one-sided 95% band of the RMSE objective, 0.1 * (1 + 1.645 / sqrt(2*294)), "
                  "stated in Section 2.3.4 with its inputs",
        "0.0086": "largest of the six leave-one-zone-out cross-coefficient shifts, a difference "
                  "between two step11_loo medians rather than a stored field",
    }
    problems = []
    checked = covered = probes = caught = 0
    cur = None
    for i, line in enumerate(body.split("\n"), 1):
        h = re.match(r"^#{2,4}\s+(\d+\.\d+(?:\.\d+)?)\s", line)
        if h:
            cur = h.group(1)
            continue
        if not cur or line.startswith("|") or line.lstrip().startswith("{{"):
            continue
        arts = mapping.get(cur) or mapping.get(cur.rsplit(".", 1)[0])
        if not arts:
            continue
        covered += 1
        pool = pool_for(arts)
        txt = re.sub(r"\$\$.*?\$\$", " ", line, flags=re.S)
        txt = re.sub(r"\$[^$]*\$", " ", txt)
        for m in re.finditer(r"(?<![\w.])(\d+\.\d+)(?![\w])", txt):
            v = m.group(1)
            if v in SKIP or _dp(v) < 2 or v in DERIVED:
                continue
            checked += 1
            f, d = float(v), _dp(v)
            # accept the value as written, as a percentage of a stored fraction, and one unit
            # low at the written precision, which is half-even rounding of a stored 0.3645
            at = pool.get(min(d, 6), set())
            # power: would this number still pass if it were wrong by one unit in its last place?
            # Reported below, because a membership test can match by coincidence and a reader is
            # entitled to know how often before trusting a pass.
            probes += 2
            caught += sum(1 for q in (f + 10 ** -d, f - 2 * 10 ** -d)
                          if q not in at and round(q - 10 ** -d, d) not in at)
            if f in at or round(f - 10 ** -d, d) in at:
                continue
            problems.append(f"paper.md line {i} (Section {cur}): {v} is not in "
                            f"{'/'.join(sorted(arts))} — {line.strip()[:64]}")
    power = (100.0 * caught / probes) if probes else 0.0
    notes = [f"{checked} number(s) in {len(mapping)} mapped section(s) held to their own artifacts; "
             f"{len(DERIVED)} derived value(s) exempt by name; sections Appendix I does not map "
             f"are not covered",
             f"measured power {power:.0f}%: that share of last-digit corruptions would be caught, "
             f"so a pass here is weak evidence and the CLAIMS registry remains the strong check"]
    return problems, notes


def main(verbose=False, release=False):
    absent = _missing_docs()
    if "RESULTS_LOG.md" in absent:
        print("RESULTS_LOG.md is not in this checkout; its section checks cannot run.")
        sections, text = {}, ""
    else:
        sections, text = log_sections()
    checks = [
        ("artifacts exist", lambda: check_artifacts_exist()),
        ("cache manifest", lambda: check_manifest()),
        ("release provenance", lambda: check_release_provenance(release)),
        ("weighting declared in every artifact", lambda: check_weighting_declared()),
        ("registered claims vs JSON paths", lambda: check_claims()),
        ("forbidden / superseded wording", lambda: check_forbidden()),
        ("figure freshness", lambda: check_figure_freshness()),
        ("section numbering", lambda: check_numbering(sections)),
        ("every step script is reproducible", lambda: check_reproduce_list()),
        ("numerical guards record their activation", lambda: check_guards()),
        ("log numbers vs artifacts", lambda: check_log_numbers(sections, verbose)),
        ("paper numbers vs artifacts", lambda: check_paper_numbers(verbose)),
        ("units and paths", lambda: check_text_rules(text)),
        ("environment claims", lambda: check_env_claims(text)),
    ]
    failed = 0
    for name, fn in checks:
        problems, notes = fn()
        status = "FAIL" if problems else "pass"
        print(f"[{status}] {name}" + (f" — {len(problems)} problem(s)" if problems else ""))
        for p in problems:
            print(f"         {p}")
        for n in notes:
            print(f"         note: {n}")
        failed += bool(problems)
    print(f"\n{len(checks) - failed}/{len(checks)} checks passed")
    if absent:
        print("Not covered, document not in this checkout: " + ", ".join(absent))
    print("Scope: this detects drift in the numbers it covers and the wordings it lists. It does not "
          "verify that the artifacts are correct, that unregistered numbers have a source, or that "
          "any two sections agree with each other. Read the docstring before quoting a pass.")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(verbose="--verbose" in sys.argv, release="--release" in sys.argv))
