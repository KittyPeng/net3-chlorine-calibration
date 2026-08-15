# Revision results log (post-review)

All experiments run against the frozen three-zone baseline in `wq_common.py`, whose network file is
a frozen copy at `models/net3_frozen/Net3.inp` (SHA-256 checked on import, so a WNTR upgrade cannot
change the model silently). Environment: conda env `water-supply` (python 3.13.12, numpy 2.4.2,
scipy 1.17.1, wntr 1.4.0); `environment.yml` declares it and `environment.lock.yml` holds the full
solve. Cache in `baseline_cache/`, with its provenance in `baseline_cache/cache_manifest.json`.

Convention: `σ = 0.1 mg/L` is **one standard deviation** of the Gaussian observation error.
Wall coefficients are reported in `m/day`; the bulk coefficient `k_b` in `day⁻¹`.

**Inference convention — read this before any number below.** Three weightings appear in this log
and they are not interchangeable:

- **formal censored Gaussian likelihood — PRIMARY.** Every headline parameter, risk and scenario
number uses it. It carries no behavioural threshold.
- **formal iid Gaussian likelihood** — the same thing with sensor-floor zeros treated as exact
measurements; used to isolate what the censoring correction is worth (Steps 9 and 7b).
- **informal GLUE score** `exp(−½(RMSE/σ)²)·1[RMSE<T]` — a **comparator**, retained because the
draft used it and because the contrast is itself a result. It is *not* a Gaussian likelihood: it
drops the factor `N = 294`, which makes it equivalent to assuming `σ_eff = σ√N = 1.71 mg/L`.

**Concentration unit convention (added after the Step 15 correction — read it before any mg/L
number below).** Every concentration in this log is in **mg/L**, and that is now true rather than
merely written. WNTR's Python API stores concentration internally in **kg/m³** whatever
`options.quality.inpfile_units` says, so the conversion happens at the model boundary in
`wq_common.build_model` (in) and `run_model` / `simulate_chlorine` (out), and nowhere else — which is
why no step script carries a factor of 1000. Alongside it the solver runs at a **deliberately
strict tolerance, `QUALITY_TOLERANCE = 1e-5 mg/L`** — note the unit: unlike concentrations this one
is *not* converted, EPANET receives it verbatim in the file's quality unit, and it is 1000× stricter
than EPANET's own 0.01 mg/L default. It is a chosen setting, not a converted one (§15.2). **Before Step 15 the same numbers were kg/m³ read as mg/L, i.e. 1000× too high**;
Step 15 documents what that did and did not change, and measures it.

**Naming the ensembles.** A weighting rule and the ensemble it produces are different things, and
conflating them is how a sentence ends up ambiguous between "the baseline configuration" and "the
informal rule". Two names are used throughout and only these two:

- the **formal likelihood-weighted ensemble** — the 8192 draws weighted by the censored likelihood;
  every headline parameter, risk and scenario number comes from it;
- the **informal GLUE behavioural ensemble** — the same draws under the informal score and its
  threshold; a comparator, and never the source of a headline number.

The bare phrase **"baseline GLUE" is not used**: it named a configuration (`k_b` fixed, unbiased
sensors) and an inference rule (the informal score) at once, and since the primary rule became the
formal likelihood those two meanings have come apart. Where the *configuration* is meant, write
"the baseline model"; where the *rule* is meant, name the rule.

Every artifact in `baseline_cache/` that holds a weighted ensemble carries a `primary_weighting` (or
`weighting`) field, and `validate_artifacts.py` fails if one does not, so no table in this log can be
read under the wrong rule by accident. Artifacts that hold no weighted ensemble at all — the cache
manifest, the warm-up convergence test, the single-parameter objective curves, the Fisher/profile
geometry and the analytic known-answer test — are exempt **by name**, and the validator prints that
exemption list rather than passing silently. Two sections are **deliberately and wholly informal**, because their subject is the
comparator itself: Step 3 (how much does the answer depend on the analyst's threshold?) and the
superseded first looks in Steps 4 and 5. Everywhere else the informal numbers appear beside the
formal ones and are labelled *comparator*.

**Timing convention.** 168 h simulation, **120 h warm-up** (chosen in Step 0 from a convergence test
rather than by convention — but a finite-horizon choice, not a converged one; see below), 48 h
assessment window `t = 120…168 h`, `N = 6 × 49 = 294` residuals. Numbers carried over from the
draft's 72 h / 24 h configuration are marked as superseded wherever they still appear.

**What 120 h does and does not mean.** The three chlorine-concentration criteria pass at 120 h. Two
others never pass inside the model's 168 h ceiling: with a 120 h warm-up the **residual
cycle-to-cycle drift** in the integrated deficit is still **5.1%**, and the **water age** p95 changes
by still **12.8 h between the last two cycles** available. So 120 h is a *pragmatic finite-horizon
warm-up under which the concentration field is cyclostationary to the stated tolerances*, and this
log never says the model reached periodic steady state or that water age converged. Everything
downstream inherits that limit, which is why the deficit-based severity metrics are reported with it
stated (Step 0, Step 10, Step 12).

**Sensor-error convention (Priority 2 #8).** Throughout, a stated sensor-error level `±X mg/L`
means the Gaussian **standard deviation σ = X**, *not* a 95% interval or a hard error bound. This
matters: read instead as a 95% band, "±0.1 mg/L" would imply `σ ≈ 0.1 / 1.96 ≈ 0.05 mg/L` and would
**halve every reported uncertainty interval** (±1σ covers ≈68%; ±1.96σ = ±0.196 covers ≈95%).
Draft fix — one sentence for the Methodology: *"In this study the stated sensor-error levels refer
to the Gaussian standard deviation σ, not to a 95% interval or a hard error bound."*

---



## Files and how to reproduce

Every number in this log comes from a script, but it is transcribed here by hand, so it can drift
when a script is re-run. `validate_artifacts.py` cross-checks the documents against
`baseline_cache/*.json`; run it after every change:

```
conda activate water-supply
cd Net3 && python provenance.py --check && python validate_artifacts.py
```

**What a passing run does and does not mean.** Being precise about this matters, because a green
check that cannot fail is worse than no check at all. The validator establishes:

- a **registered claim** equals a named value at a named JSON path. This is the only check with real
teeth, and it fails both when the number moves *and* when the sentence is reworded so that the
anchor disappears — the failure mode that previously let a whole stale table survive;
- every result artifact declares its weighting;
- no number in the log is a near-miss of a JSON value (the signature of a stale transcription), with
table cells held to their own section's artifact from three decimals up;
- no forbidden phrasing survives — superseded configurations (pre-Sobol draw counts, the draft
warm-up length), labels naming the wrong inference rule, and over-claims such as treating a
single-run SD–CRLB match as frequentist efficiency, or calling the warm-up fully converged;
- the environment and the frozen `.inp` match the ones that produced the cache.

It does **not** establish that the artifacts are correct, that an *unregistered* number in the prose
came from anywhere in particular, or that any two sections agree with each other. Semantic
consistency — a missing caveat, two sections drawing opposite conclusions, a correct number described
in the wrong words — is not machine-checkable here and was found by reading. Figure freshness is a
file-timestamp comparison, not a check that a figure shows what its caption says.

Numbers that describe a superseded run are exempt from the number check and must say so in the
sentence ("an earlier run…", "the draft's earlier…").


| file                             | purpose                                                                                                                                  | writes                                                                                      |
| -------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------- |
| `wq_common.py`                   | frozen baseline config + WNTR/EPANET helpers (monitors, zone assignment, seeds, priors, timing)                                          | — (imported by all step scripts)                                                            |
| `step0_warmup_convergence.py`    | is the warm-up long enough? successive-24 h-cycle convergence test with pre-declared criteria                                            | `figures/step0_warmup_convergence.png`, `step0_warmup_convergence.json`                     |
| `step1_freeze_baseline.py`       | synthetic truth + noisy observations + 8192-draw Sobol library; all three weightings; cache every prediction                             | `baseline_cache/baseline.npz`, `baseline_cache/baseline_meta.json`                          |
| `step3_threshold_sensitivity.py` | behavioural-threshold sweep of the informal comparator, from the cache (no EPANET)                                                       | `baseline_cache/step3_threshold.json`                                                       |
| `step4_displaced_prior.py`       | displaced prior, first look (informal, single realisation; superseded by 4d)                                                             | `baseline_cache/step4_displaced_prior.json`                                                 |
| `step4b_sensitivity_curves.py`   | single-parameter RMSE curves (why old is one-sided); no weighting applied                                                                | `figures/step4b_sensitivity_curves.png`, `step4b_sensitivity.json`                          |
| `step4d_displaced_robust.py`     | robust displaced prior (formal primary + informal at 2 thresholds, 30 noise, both directions)                                            | `figures/step4d_displaced_robust.png`, `step4d_displaced_robust.json`, `step4d_preds_*.npy` |
| `step5_structural_error.py`      | structural error, ±20% pipe-level jitter (informal, single field; superseded by 5c/5d)                                                   | `figures/step5_structural_error.png`, `step5_structural_error.json`                         |
| `step5c_jitter_sweep.py`         | structural error, jitter sweep ±20/35/50% + 25-field ensemble (formal primary)                                                           | `figures/step5c_jitter_sweep.png`, `step5c_jitter_sweep.json`                               |
| `step5d_structured.py`           | structural error, length-correlated within-zone (formal primary + informal comparator)                                                   | `figures/step5d_structured.png`, `step5d_structured.json`                                   |
| `step6_noise_sensitivity.py`     | sensor-accuracy sweep σ=0.02/0.05/0.10/0.15 (formal primary + informal comparator)                                                       | `figures/step6_noise_sensitivity.png`, `step6_noise_sensitivity.json`                       |
| `step7_fisher.py`                | Fisher/CRLB a-priori identifiability (k_w, +k_b, +sensor offsets); scale-dependent FD steps                                              | `figures/step7_fisher.png`, `step7_fisher.json`                                             |
| `step7b_profile.py`              | profile likelihood: CONTINUOUS 95% intervals under censored (primary) and iid likelihoods                                                | `figures/step7b_profile.png`, `step7b_profile.json`, `step7b_rmse_grid.npy`                 |
| `step7c_ar1.py`                  | AR(1) covariance Fisher/CRLB (per-coefficient widening) + ρ sensitivity sweep                                                            | `baseline_cache/step7c_ar1.json`                                                            |
| `step7c_profile_ar1.py`          | AR(1) profile likelihood (rebuilds 21³ residual grid)                                                                                    | `baseline_cache/step7c_resid_grid.npy`, `step7c_profile_ar1.json`                           |
| `step8_sensor_bias.py`           | systematic sensor bias at node 15, two-sided (formal primary; informal curvature contrast)                                               | `figures/step8_sensor_bias.png`, `step8_sensor_bias.json`                                   |
| `step8b_kb_sensitivity.py`       | calibration at k_b ± 20% — bulk–wall compensation (formal primary + comparator; Priority-2 #5)                                           | `figures/step8b_kb_sensitivity.png`, `step8b_kb_sensitivity.json`, `step8b_preds_kb*.npy`   |
| `step8c_bias_bynode.py`          | bias-location sweep across all six monitors, both signs (formal primary + comparator)                                                    | `baseline_cache/step8c_bias_bynode.json`                                                    |
| `step8d_sensor_drift.py`         | sensor DRIFT (linear ramp) against mean- and end-equivalent constant-bias controls                                                       | `figures/step8d_sensor_drift.png`, `step8d_sensor_drift.json`                               |
| `step9_zeroclip.py`              | zero-floor censoring (L=0): formal censored vs formal iid (naive exact zero)                                                             | `figures/step9_zeroclip.png`, `step9_zeroclip.json`                                         |
| `step10_risk_metrics.py`         | risk duration/depth (trapezoid, 48 h), water-age corroboration, report-step sensitivity                                                  | `figures/step10_risk_metrics.png`, `step10_risk_metrics.json`                               |
| `step11_loo.py`                  | held-out validation A–D: monitor, zone, unmonitored junctions, heterogeneous truth                                                       | `figures/step11_loo.png`, `step11_loo.json`                                                 |
| `step12_scenarios.py`            | temperature/ageing scenario projection of the formal ensemble + dosing + risk register                                                   | `figures/step12_*.png`, `step12_scenarios.json`, `step12_risk_register.csv`                 |
| `step13_known_answer.py`         | known-answer test: single pipe with an analytic first-order solution (Priority-2 #7)                                                     | `baseline_cache/step13_known_answer.json`                                                   |
| `step14_repeated_noise.py`       | repeated-noise calibration: bias, estimator SD vs CRLB, interval coverage over 100 realisations                                          | `figures/step14_repeated_noise.png`, `step14_repeated_noise.json`                           |
| `provenance.py`                  | records what produced the cache: config hash, frozen `.inp` hash, `wq_common.py` and step-script hashes, git tree hash, library versions | `baseline_cache/cache_manifest.json`                                                        |
| `validate_artifacts.py`          | registered claims, weighting declarations, number drift, forbidden wording; exits non-zero on failure                                    | — (report on stdout)                                                                        |




Run from the `Net3/` directory with the environment activated (no absolute interpreter paths, so the
commands work on any machine). Approximate timings on the development machine, where one 168 h
EPANET water-quality run takes about 35 ms:

```
conda activate water-supply
cd Net3
export MPLCONFIGDIR=../.mplcache
python step0_warmup_convergence.py      # warm-up convergence test (decides WARMUP_H)
python step1_freeze_baseline.py         # ~280 s (8192 EPANET runs, all nodes) — builds the cache
python step3_threshold_sensitivity.py   # instant (cache only)
python step4_displaced_prior.py         # ~350 s (8192 EPANET runs, monitors only)
python step4b_sensitivity_curves.py     # single-parameter RMSE curves (re-simulates a sweep)
python step4d_displaced_robust.py       # ~620 s first run, instant afterwards (prediction library cached)
python step5_structural_error.py        # heterogeneous truth re-simulated; candidates reused
python step5c_jitter_sweep.py           # jitter sweep + 25-field ensemble (one truth run per field)
python step5d_structured.py             # length-correlated truth + correlation-strength dose
python step6_noise_sensitivity.py       # ~20 s (cache only)
python step7_fisher.py                  # ~5 s   (~80 EPANET runs)
python step7_verify.py                  # Jacobian sanity check; prints only, writes no artifact
python step7b_profile.py                # ~250 s given the cached 21³ grid; ~15 min if it must rebuild
python step7c_ar1.py                    # AR(1) Fisher/CRLB (reuses the Jacobian)
python step7c_profile_ar1.py            # AR(1) profile; rebuilds the 21³ residual grid if unkeyed
python step8_sensor_bias.py             # ~20 s (cache only)
python step8b_kb_sensitivity.py         # ~620 s first run, instant afterwards
python step8c_bias_bynode.py            # ~100 s (cache only)
python step8d_sensor_drift.py           # ~12 min (cache only; 25 arms × 30 realisations)
python step9_zeroclip.py                # censored vs naive-zero comparison (cache only)
python step10_risk_metrics.py           # ~6 s   (+ report-step sensitivity runs)
python step11_loo.py                    # ~70 s
python step12_scenarios.py              # ~9 min
python step13_known_answer.py           # analytic known-answer test
python step14_repeated_noise.py         # ~10 s (cache only, 100 noise realisations)
python step15_unit_equivalence.py       # ~90 s (3 x 256 EPANET runs; no cache)
python step15b_full_regression.py       # full-library + artifact + log regression; re-runs
                                        # the 8192 legacy arm and reads two git commits
python provenance.py                    # refresh baseline_cache/cache_manifest.json
python validate_artifacts.py            # cross-check the documents against the artifacts
```

Every `step*.py` in this directory appears above; `validate_artifacts.py` fails if one does not, so
a new step cannot be added without also becoming reproducible. Timings are indicative and were
measured on the development machine; only the ones marked with a number were actually timed.

**Provenance of a release.** `provenance.py` records the git commit and whether the tree was dirty
when it ran. Running it *before* committing therefore always records `dirty: true` — the manifest
then describes a working tree that no longer exists. For a state meant to be cited, commit first,
then run `provenance.py` on the clean tree and commit the manifest by itself:

```
git commit -am "…"                                   # code, artifacts and documents
cd Net3 && python provenance.py                      # now records dirty: false
git add baseline_cache/cache_manifest.json && git commit -m "chore(net3): record clean provenance"
```

The manifest then names the commit whose code and artifacts produced the results; the manifest's own
commit is one later, which is unavoidable and is not a discrepancy.

Two step scripts must not run at the same time **from the same directory**: WNTR writes its EPANET
scratch files as `temp.inp|rpt|bin` in the working directory, so concurrent runs overwrite each
other. They *can* run in parallel from different working directories, because every path a step
script reads or writes is absolute — `cd` to a scratch directory and invoke the script by full path.

`baseline_cache/baseline.npz` (~124 MB at 8192 draws × 49 h × 92 nodes) and the `*.npy` prediction
libraries are git-ignored; rebuild them with the step that writes them. The small `*.json` summaries,
`cache_manifest.json` and this log are version-controlled.

---



## Step 0 — how long must the warm-up be? (pre-declared convergence test)

The draft discarded the first 24 h and assessed 24–72 h, but that value was never justified: the
tanks start at an assumed 0.5 mg/L, the leading-risk junctions have mean water ages of 30–70 h over
the assessed window (Step 10), and the Step 12 paired test moves the continuous severity metrics by
+8.6% to +13.0% between the draft's horizon and the current one. This step decides the question with criteria fixed *before* the numbers were
seen — and finds that **no warm-up available in this model satisfies all of them**, which is why the
answer is a defensible compromise rather than a convergence result.

**What "warmed up" means here.** Net3's demand patterns are 24 points at 1 h and pump 10 runs on a
24 h schedule, so the target state is not a constant steady state but a **cyclostationary** one: the
field over one diurnal cycle repeats in the next. Each parameter set is run once over the full
horizon and successive 24 h cycles are differenced.

**Formula**:

```
Δ(k) = max over (t in cycle k, node n) | X(t + 24, n) − X(t, n) |
```

**Properties**:

- Units follow `X` (mg/L for chlorine, h for water age, m for tank level).
- `Δ(k) → 0` means the field from hour `24k` onward is cyclostationary, so a warm-up of `24k` is
sufficient.
- Computed for the truth and for both corners of the prior box; the weak corner is the binding case
because the least reactive network forgets its initial condition most slowly.

**Hard horizon limit — this bounds what can be claimed.** Pump 10 is driven by *absolute-time*
controls (`IF SYSTEM TIME IS 01:00:00`, `25:00:00`, …) enumerated only to **159 h**, and the model
duration is 168 h. Past 168 h the pump would stay closed permanently, so a longer run is a different
system, not a longer warm-up. The test therefore has 7 cycles and cannot certify any warm-up beyond
144 h. Step 12's "long" setting (168 h / 120 h) was already sitting exactly on this ceiling.

**Result — worst value across the three parameter sets, per cycle pair:**


| criterion                            | tolerance | 0–24   | 24–48  | 48–72  | 72–96  | 96–120 | 120–144 | earliest warm-up that satisfies it                        |
| ------------------------------------ | --------- | ------ | ------ | ------ | ------ | ------ | ------- | --------------------------------------------------------- |
| tank level (m)                       | 0.05      | 0.8185 | 0.0496 | 0.0070 | 0.0009 | 0.0001 | 0.0001  | **24 h** (verified)                                       |
| monitor chlorine (mg/L)              | 0.005     | 0.9209 | 0.0437 | 0.0234 | 0.0122 | 0.0064 | 0.0034  | **120 h** (verified)                                      |
| network p95 chlorine (mg/L)          | 0.010     | 0.9511 | 0.1061 | 0.0396 | 0.0201 | 0.0103 | 0.0054  | **120 h** (verified)                                      |
| tank chlorine (mg/L)                 | 0.010     | 0.2457 | 0.0787 | 0.0416 | 0.0219 | 0.0115 | 0.0061  | **120 h** (verified)                                      |
| risk severity, rel. change in `E[A]` | 0.02      | 0.8983 | 0.6104 | 0.1206 | 0.1897 | 0.0978 | 0.0509  | not within 168 h; extrapolates to ≈168 h (**unverified**) |
| water age p95 (h)                    | 1.0       | 23.999 | 21.676 | 19.043 | 16.692 | 14.630 | 12.815  | not within 168 h; extrapolates to ≈600 h (**unverified**) |


**Three findings, in order of consequence.**

1. **The hydraulics settle almost immediately; the chemistry does not.** Tank levels are periodic to
   0.0001 m after three cycles, so the driver is not the problem. Chlorine differences decay
   geometrically by roughly a factor of two per cycle and first meet all three concentration
   tolerances at a **120 h** warm-up — five times the baseline value.
2. **The integrated risk severity has not converged even at the horizon.** Network-mean cumulative
   deficit for the truth runs 1.005 (start-up) → 0.1348 → 0.0821 → 0.092 → 0.1029 → 0.1092
   mg/L·h: it dips after the start-up transient and then climbs **monotonically**. Small
   concentration differences still move the integral, because the deficit integrates a
   threshold crossing. The relative change halves each cycle and would reach the 2% criterion at
   about cycle 6, i.e. a 144–168 h warm-up — but the model horizon cannot supply that cycle, so it
   stays an extrapolation. This is the same effect Step 12's paired test measures as +8.6% to +13.0%,
   now with its cause identified: **the 24 h window sits on the descending limb of the transient,
   not on the plateau.**
3. **Water age is horizon-dependent and cannot be converged in this model at all.** The p95
   **water age** change is still **12.8 h between the last two cycles** available and decays by only
   ~12% per cycle; extrapolation puts the crossing near 600 h, far past the 168 h ceiling. So
   `mean_age_h` in Step 10 is a property of the chosen window, not an equilibrium water age. The
   water-age corroboration must be phrased as a rank association within a fixed window; the absolute
   ages must not be reported as steady-state values. Step 10 stores all three window definitions
   (`mean_age_h`, `mean_age_last24_h`, `age_final_h`) so the choice is explicit.

**Configuration this implies for the pipeline.** Warm-up **120 h** with the existing 48 h assessment
window gives a 168 h total — *exactly* the model horizon, with nothing to spare. Conveniently the
residual count is unchanged: `6 monitors × 49 hours = 294`, because `168 − 120 = 72 − 24 = 48`, so
the informal comparator's threshold `0.107` and everything derived from `N = 294` carry over
untouched.

**What may and may not be written about this.** Three of the six criteria are met at 120 h and two
are not, so the honest formulation is bounded on both sides:

- **Sayable:** *"A 120 h warm-up was selected because the chlorine concentration field first meets
the pre-declared cyclostationarity tolerances at the monitors, network-wide and in the tanks at
that point. It is a finite-horizon choice: with a 120 h warm-up the residual cycle-to-cycle drift in the integrated deficit is still 5.1%, and water age still 12.8 h between the last two cycles within the model's 168 h ceiling."*
- **Not sayable:** that the model is fully cyclostationary after 120 h, or that water age has settled
to a horizon-independent equilibrium. Neither is true, and the second is not achievable in this
model at any warm-up available under the pump-control ceiling.

**Two routes were available and the cheaper one was taken deliberately.** The stronger route is to
extend pump 10's absolute-time controls programmatically over 15–25 diurnal cycles, verify that the
extended model reproduces the original bit-for-bit over the first 168 h, and only then run to
convergence of the severity and age criteria. That is the right thing to do if the warm-up itself is
a research question. Here it is not: it is a nuisance parameter, the concentration field — which is
what every likelihood in this log is built from — is cyclostationary to tolerance at 120 h, and
changing the warm-up would invalidate every cached result. The 5.1% deficit drift is therefore
carried as a stated limitation on the *severity* metrics (Steps 10 and 12) rather than removed.

---



## Step 1 — the frozen baseline

Config: monitors `107/113/15/145/209/231`; inlet 1.0, tank 0.5; **168 h simulation, 120 h warm-up**
(Step 0), 1 h reporting, 5 min quality step; `k_b = -0.5`; true
`(k_w,old, k_w,avg, k_w,new) = (-1.0, -0.1, -0.05) m/day`; noise seed 42, Sobol scramble seed 0;
**8192 = 2¹³ scrambled-Sobol prior draws**. Assessment window 48 h, so `N = 6 × 49 = 294`
residuals, unchanged.

This is no longer a reproduction of the draft. Three choices were changed deliberately, each with a reason recorded here rather than inherited:


| choice    | draft                    | now                        | why                                                                                       |
| --------- | ------------------------ | -------------------------- | ----------------------------------------------------------------------------------------- |
| warm-up   | 24 h                     | 120 h                      | Step 0: chlorine is not cyclostationary before then                                       |
| sampling  | pseudo-random, 2¹¹ scale | 8192 scrambled Sobol       | at the draft's sample size the formal likelihood had an effective sample size of only ~37 |
| weighting | one informal score       | three schemes side by side | the informal score is not a likelihood; see below                                         |


**Three weighting schemes, same data, same draws.** The formal censored Gaussian likelihood is the
primary analysis; the formal iid one isolates the cost of treating a sensor-floor zero as an exact
measurement; the informal GLUE score is retained as a comparator.

**Formula** (per candidate, up to an additive constant):

```
formal iid       ℓ = −(1/2) Σ_j,t ((y_jt − μ_jt)/σ)²
formal censored  ℓ = −(1/2) Σ_{y>0} ((y_jt − μ_jt)/σ)²  +  Σ_{y=0} log Φ(−μ_jt/σ)
informal GLUE    ℓ = −(1/2) (RMSE/σ)²                    · 1[RMSE < threshold]
```

**Properties**:

- The informal score is the formal iid one **divided by** `N = 294`. Equivalently it is a Gaussian
likelihood with `σ_eff = σ√N = 0.1 × 17.15 = 1.71 mg/L` — 17× the sensor noise and larger than the
1.0 mg/L inlet concentration. That is why it is nearly flat inside the behavioural set.
- The formal schemes carry **no threshold**; a hard cut-off belongs to the GLUE comparator, not to a
likelihood.
- 10 of the 294 window points are clipped at the sensor floor (43 of 1014 over the full record), so
the censored and iid schemes are not identical.


| scheme                           | ESS  | entropy (bits, of 13.00) | `k_w,old`                | `k_w,avg`                | `k_w,new`                |
| -------------------------------- | ---- | ------------------------ | ------------------------ | ------------------------ | ------------------------ |
| formal censored (primary)        | 157  | 7.97                     | -0.9876 ± 0.0938 (25.0%) | -0.1091 ± 0.0142 (30.7%) | -0.0440 ± 0.0078 (28.3%) |
| formal iid                       | 154  | 7.94                     | -0.9708 ± 0.0920 (24.5%) | -0.1087 ± 0.0141 (30.6%) | -0.0441 ± 0.0078 (28.3%) |
| informal GLUE, thr 0.107         | 4783 | 12.22                    | -1.0115 ± 0.2668 (71.1%) | -0.1170 ± 0.0424 (91.8%) | -0.0492 ± 0.0242 (88.2%) |
| informal GLUE, thr 0.120 (draft) | 7062 | 12.79                    | -0.9431 ± 0.3220 (85.8%) | -0.1196 ± 0.0455 (98.4%) | -0.0529 ± 0.0269 (98.0%) |


Bracketed figure is SD retained (= posterior SD / prior SD). The draft's configuration retained
**85.8 / 98.4 / 98.0%** of the prior width for old / average / new; the formal censored retains 25.0 / 30.7 / 28.3 % of it on exactly the same observations. Nothing about the data
changed — only how much of their information the weighting extracts. This single contrast is the most
important result in the log, and every later section is either a consequence of it or a test of how
far it survives realistic error.

**Sampling convergence.** Leading `2^k` subsets of a scrambled Sobol set are themselves balanced
designs, so these four rows are exact sub-designs rather than ad-hoc thinning:


| draws | ESS   | `k_w,old` median | `k_w,avg` median | `k_w,new` median | max median drift |
| ----- | ----- | ---------------- | ---------------- | ---------------- | ---------------- |
| 1024  | 18.8  | -0.9744          | -0.1096          | -0.0444          | —                |
| 2048  | 39.3  | -0.9752          | -0.1097          | -0.0435          | 0.033 prior SD   |
| 4096  | 78.0  | -0.9787          | -0.1094          | -0.0437          | 0.009 prior SD   |
| 8192  | 156.7 | -0.9787          | -0.1090          | -0.0437          | 0.010 prior SD   |


ESS is proportional to the number of draws (it stays at 1.9% of them), so prior sampling cannot be
made efficient by brute force — but the **quantiles have converged** anyway: the medians move by
0.01 prior SD between the last two sizes and the 5/95 endpoints by less. ESS is the conservative
diagnostic here; the Sobol design's space-filling is what buys the stability.

Minimum candidate RMSE 0.0971 mg/L, against a **realised noise RMSE of 0.0973 mg/L on the same
120–168 h calibration window** — the comparison has to be window-against-window, and the
full-record realised noise RMSE (0.0960, `noise_rmse_full_record`) is a different sample and must
not be quoted here. Behavioural counts, for the comparator:
4786/8192 at 0.107 and 7084/8192 at 0.120.

Cache `baseline.npz` holds all 8192 candidate predictions `C_all (8192 × 49 × 92)` plus both formal
log-likelihoods, so every later step recomputes weights **without re-running EPANET**.

---



## Step 2 — prior vs behavioural identifiability (corrects the overclaim)

Prior SD of a uniform range `[a, b]`:

```
prior_SD = (b − a) / √12
```


| Group   | Prior range     | Prior mid | Prior SD | True  | informal mean ± SD (retained) | formal mean ± SD (retained) |
| ------- | --------------- | --------- | -------- | ----- | ----------------------------- | --------------------------- |
| old     | [-1.5, -0.2]    | -0.850    | 0.375    | -1.00 | -0.9431 ± 0.3220 (85.8%)      | -0.9876 ± 0.0938 (25.0%)    |
| average | [-0.2, -0.04]   | -0.120    | 0.046    | -0.10 | -0.1196 ± 0.0455 (98.4%)      | -0.1091 ± 0.0142 (30.7%)    |
| new     | [-0.10, -0.005] | -0.0525   | 0.027    | -0.05 | -0.1196 → see note            | -0.0440 ± 0.0078 (28.3%)    |


The informal column uses the draft's threshold 0.120 so that it reproduces the configuration the overclaim came from; `k_w,new` there is -0.0529 ± 0.0269 (98.0% retained).

Interpretation, which has to be split by scheme because the two disagree about the *conclusion*, not just the digits:

- **Under the informal GLUE score the original criticism stands.** All three distributions retain
86–98% of the prior width and the weighted means sit near the prior midpoints. Their apparent
agreement with the true values is a consequence of the prior ranges having been centred near
those values, not of the data.
- **Under the formal likelihood it does not.** The same observations contract all three
coefficients to 25–31% of the prior width. So "the data do not identify average and new" was
never a statement about the data; it was a statement about the score used to weight them.

Deleted overclaim (was in §4.4): *"The weighted mean of every group lies close to its true value."*
The replacement claim must name the inference convention: under the informal GLUE score the means
are prior-dominated, under the formal likelihood all three are informed but still biased by up to
~0.7 posterior SD on this single noise realisation.

---



## Step 3 — behavioural-threshold sensitivity (no EPANET re-run)

Sampling SD of the RMSE statistic at the truth:

```
sd(RMSE) ≈ σ / √(2N) = 0.10 / √(2·294) = 0.0041 mg/L      (N = 6 monitors × 49 h = 294)
```

Observed minimum candidate RMSE = 0.0971 mg/L, against a realised noise RMSE of **0.0973 mg/L**
over the same 120–168 h calibration window. (The full-record realised noise RMSE is 0.0960; it
describes the observation set as a whole and is not the comparator for a window RMSE.)

**Three different quantities get called "the noise floor" in this literature, so this log names
them.** `σ = 0.10 mg/L` is the *nominal* noise level the observations were generated with;
`0.0973` is the *realised* noise RMSE on the calibration window; `0.0960` is the realised noise
RMSE over the full record. The table below is scaled on `σ`, not on either realised value.

**Scope of this step, restated.** The behavioural threshold belongs to the informal GLUE comparator
only — the formal likelihood carries no cut-off. So this is a sensitivity analysis *of the
comparator*, and the right conclusion to draw from it is about how much the comparator's answer
depends on an analyst's choice, not about how much the data say.

Retention and objective-scale table:


| Threshold | Retained | Retention | SD above σ | band node15 | band node107 | nonzero-risk nodes |
| --------- | -------- | --------- | ---------- | ----------- | ------------ | ------------------ |
| 0.107     | 4786     | 58.4%     | 1.70                 | 0.1008      | 0.0762       | 24                 |
| 0.110     | 5668     | 69.2%     | 2.42                 | 0.1121      | 0.0788       | 24                 |
| 0.120     | 7084     | 86.5%     | 4.85                 | 0.1419      | 0.0816       | 24                 |


Per-coefficient behavioural distribution — weighted `mean ± SD (m/day)` with `SD retained`
(= behavioural SD / prior SD; prior SD = 0.375 / 0.046 / 0.027 m/day for old / avg / new):


| Threshold                        | k_w,old mean ± SD (retained) | k_w,avg mean ± SD (retained) | k_w,new mean ± SD (retained) |
| -------------------------------- | ---------------------------- | ---------------------------- | ---------------------------- |
| 0.107                            | -1.0115 ± 0.2668 (71.1%)     | -0.1170 ± 0.0424 (91.8%)     | -0.0492 ± 0.0242 (88.2%)     |
| 0.110                            | -0.9956 ± 0.2844 (75.8%)     | -0.1174 ± 0.0437 (94.5%)     | -0.0512 ± 0.0254 (92.7%)     |
| 0.120                            | -0.9431 ± 0.3220 (85.8%)     | -0.1196 ± 0.0455 (98.4%)     | -0.0529 ± 0.0269 (98.0%)     |
| *formal censored (no threshold)* | *-0.9876 ± 0.0938 (25.0%)*   | *-0.1091 ± 0.0142 (30.7%)*   | *-0.0440 ± 0.0078 (28.3%)*   |


As the threshold tightens (0.120 → 0.107):

- **old sharpens a little**: SD 0.3220 → 0.2668 (SD retained 85.8% → 71.1%) and the mean moves
toward the truth (-0.9431 → -1.0115, true -1.0).
- **average / new barely change**: SD retained stays 88–98% and the means stay near the prior
midpoint.
- **But no tested threshold gets close to the formal result** — and the reason is structural, not a property of the particular values tried. A hard cut-off multiplies the score by an indicator; it cannot restore the factor `N` the informal score drops, so no choice of cut-off turns it into the Gaussian likelihood. What the sweep shows empirically is that within the range that keeps a usable acceptance rate, Tightening from 0.120 to 0.107 recovers
about a seventh of the gap for `old` and almost none for `avg`/`new`, whereas switching to the
formal likelihood cuts every SD by a factor of three. The earlier reading — "average/new cannot be
sharpened by any threshold, so their non-identifiability is a property of the monitoring array" —
is wrong: it is a property of the *score*. The monitoring array does inform them, as the formal row
and Step 7's CRLB both show.

Top-6 risk nodes are **stable** across thresholds: `131 > 129 > 141 > 127` then `15`/`139` swapping
places between 0.107 and the looser two (nonzero-risk node count fixed at 24). The identity of the
leading set is therefore insensitive to the threshold; only the order of the fifth and sixth entries
moves, and those two are separated by less than the ensemble spread on either.

Conclusions:

- The draft's 0.12 threshold sits 4.85 sd(RMSE) above the **nominal** noise level `σ = 0.10`
(`sd_above_floor = (thr − σ) / sd(RMSE)`), so a parameter set at the truth passes with
near-certainty — the filter only rejects grossly wrong sets.
- Tightening to a defensible ≈95% band (0.107) sharpens **old** and leaves **average/new** near the
prior. **This is a limit of the score, not of the monitoring array.** An earlier version of this
log drew the opposite conclusion from the same table, and that was the single most consequential
error in the draft: it turned an inefficiency of the weighting rule into a claimed property of the
sensor network. The formal row in the table above, the Case-A CRLB (Step 7) and the continuous
profile intervals (Step 7b) all show the array does inform average and new.
- The operational risk ranking is robust to the threshold choice, which is a genuinely useful
finding and is unaffected by the correction above.

---



## Step 4 — displaced-prior experiment (PRELIMINARY, single noise realisation)

Same truth, same monitoring array, same baseline observations; **only the priors move**.
Each range keeps its width but its midpoint is set to `truth − 1 prior SD`. Displacement is
downward (toward stronger decay) because moving upward would make the upper bound positive;
all displaced ranges stay non-positive and still contain the truth (the truth then sits one
prior SD above the displaced midpoint).

Question: do the observations pull the behavioural mean back toward the truth?

behavioural 3359/8192 (41.0%), min RMSE 0.0972. Weighting: informal GLUE, as in the draft.


| Group   | Orig mid | Displaced range (mid)     | True  | Beh. mean ± SD   | SD retained | Gap to truth closed |
| ------- | -------- | ------------------------- | ----- | ---------------- | ----------- | ------------------- |
| old     | -0.850   | [-2.025, -0.725] (-1.375) | -1.00 | -1.2185 ± 0.3287 | 87.6%       | 42%                 |
| average | -0.120   | [-0.226, -0.066] (-0.146) | -0.10 | -0.1272 ± 0.0407 | 88.1%       | 41%                 |
| new     | -0.0525  | [-0.125, -0.030] (-0.077) | -0.05 | -0.0556 ± 0.0175 | 63.9%       | 80%                 |


("Gap to truth closed" = distance the behavioural mean travels from the displaced midpoint
toward the truth, as a fraction of the 1 prior SD gap.)

Reading:

- The observations pull the displaced means **substantially** back toward the truth — 42%, 41% and 80% of the gap — and acceptance drops to 41%, so the objective is discriminating here.
- **This reverses what the draft's earlier run showed.** At the 24 h warm-up the same experiment
closed only 9%, 12% and 22% of the gap and retained 99.8% / 96.9% / 92.8% of the prior width, which
is what the "even the informed coefficient is not recovered" reading rested on. Under the corrected
warm-up the data are far more informative than that reading assumed.
- **Attribution is not settled.** Two things changed at once between the two runs: the warm-up
(24 h → 120 h) and the sampling design (superseded pre-Sobol draws → 8192 Sobol). The paired-window design
used in Step 0 — one long simulation, two windows, shared noise and shared draws — is what would
separate them, and it has not been run for this experiment. Until it is, the honest statement is
that the displaced-prior conclusion is warm-up-sensitive, not that the warm-up caused all of it.
- The one-sided-identifiability claim for `k_w,old` therefore needs re-testing too; it was inferred
from the 9% figure. Step 4b's objective curves and Step 4d's 30-realisation version are the
evidence that should carry it, not this table.

Status: this single-realisation table is a first look. The robust version (30 noise
realisations × two thresholds, plus an old-toward-weaker displacement) is **Step 4d** below
and supersedes these numbers for reporting.

---



## Step 4b — single-parameter objective curves (why old is one-sided, avg/new flat)

Each coefficient is swept across (and beyond) its prior while the other two are held at the
truth; RMSE is computed against the baseline noisy observations
(`step4b_sensitivity_curves.py`).
Realised noise RMSE (all three at truth) on the calibration window = **0.0973 mg/L**
(`noise_floor_rmse`); the minimum
over the sweep is 0.0972 mg/L. No weighting is applied here — the object plotted is the objective
itself, and the two thresholds are drawn only to show where the informal comparator would cut it.

Single-parameter objective: old is an asymmetric valley, avg/new are flat


| Coefficient | RMSE swing across its prior   | in units of sd(RMSE) | reading                                                                                             |
| ----------- | ----------------------------- | -------------------- | --------------------------------------------------------------------------------------------------- |
| old         | 0.0412 mg/L (0.0972 → 0.1384) | 10.0                 | steep on the weak side (weak edge exceeds the 0.12 threshold → rejectable); flat on the strong side |
| avg         | 0.0075 mg/L (0.0973 → 0.1048) | 1.8                  | whole prior below the comparator's threshold → values indistinguishable *to that rule*              |
| new         | 0.0092 mg/L (0.0972 → 0.1064) | 2.2                  | whole prior below the comparator's threshold → values indistinguishable *to that rule*              |


- **old is an asymmetric valley**: the objective climbs above the behavioural threshold only on
the weak side, so the data reject weak old but cannot distinguish among strong values (flat
floor from about -0.9 to -2.2 m/day). This is the mechanism behind one-sided identifiability,
and reconciles Steps 2–3 with Step 4: old contracts under the baseline prior (which straddles
the steep wall) but not under the displaced strong-regime prior (which sits on the flat floor).
- **avg / new**: sweeping them across their entire prior moves the objective by < 0.01 mg/L and
keeps it below the 0.107/0.12 thresholds, so every value passes the comparator's filter. But the
swing is **1.8–2.2 × sd(RMSE)**, i.e. it is *not* inside the noise: a rule that uses all `N = 294`
residuals rather than their root-mean-square can see it, which is exactly what the formal likelihood
does and why it contracts avg/new to ~30% of the prior. The correct reading is therefore "the
monitors are weakly but genuinely sensitive to avg/new, and a threshold on RMSE discards that
sensitivity", not "the monitors are insensitive to them".

---



## Step 4d — robust displaced-prior (30 noise realisations × two thresholds)

Two displacement designs (`step4d_displaced_robust.py`), **each displacing all three**
coefficients (width fixed, truth kept inside), so every gap-closed value is meaningful. The two
designs differ only in old's direction:

- **DOWN**: old / avg / new midpoints all set to `truth − 1 prior SD` (strong-decay side).
- **OLDUP**: old displaced toward the weak/steep side, nominally to `truth + 1 prior SD`, with the
upper bound capped at −0.005 to stay non-positive; avg / new displaced down as in DOWN.
**The cap binds here, and by enough to matter**: the window had to be shifted down by
**0.0303 m/day**, giving an OLDUP old prior of `[-1.305, -0.005]`. Its midpoint therefore sits about
**0.92** prior SD above the truth rather than the intended 1.0, so the displacement actually tested
is a little smaller than the design nominally specifies. The shift is recorded in the artifact
(`upper_cap_shift_applied`) so the realised design can be reconstructed from it — a guard that moves
the design and leaves no trace is not auditable.

Each design: 8192 scrambled-Sobol forward runs cached once (keyed on the design box, so editing a
design invalidates only its own library), then the gap-closing statistic recomputed over 30
independent noisy observation sets. Reported as **median [IQR]** across noise. `gap closed` =
fraction of the displaced-midpoint → truth distance covered by the posterior mean (0% = stayed at the
prior midpoint, 100% = reached the truth).

**Primary — formal censored likelihood, no threshold** (median ESS 152 DOWN / 147 OLDUP):


| coefficient | DOWN (old displaced to the strong side) | OLDUP (old displaced to the weak side) |
| ----------- | --------------------------------------- | -------------------------------------- |
| old         | **88% [67–104]** (SD ret 27%)           | **112% [95–132]** (SD ret 25%)         |
| avg         | **95% [75–110]** (SD ret 29%)           | 95% [75–111] (SD ret 29%)              |
| new         | **96% [80–122]** (SD ret 28%)           | 96% [80–122] (SD ret 29%)              |


**Comparator — informal GLUE**, at the draft's threshold 0.12 (retention 87% DOWN / 64% OLDUP) and at
the defensible ≈95% band 0.107 (retention 41% / 30%):


| coefficient | 0.12, DOWN        | 0.12, OLDUP       | 0.107, DOWN       | 0.107, OLDUP      |
| ----------- | ----------------- | ----------------- | ----------------- | ----------------- |
| old         | 9% [5–13] (99%)   | 61% [54–64] (68%) | 28% [14–55] (92%) | 85% [75–94] (52%) |
| avg         | 14% [10–19] (97%) | 15% [10–20] (97%) | 43% [27–58] (87%) | 43% [29–60] (88%) |
| new         | 23% [14–32] (92%) | 14% [9–21] (95%)  | 78% [56–82] (65%) | 66% [42–71] (69%) |


Findings, and they are not the ones the previous version of this section reported:

1. **Under the primary rule all three coefficients are pulled essentially all the way back, from
   either direction.** The medians are 88–112% with IQRs that straddle 100%, and the posterior width
   collapses to 25–29% of the prior in every cell. A displaced prior is the cleanest possible test for
   prior domination, and the formal likelihood passes it for average and new as decisively as for old.
2. **The comparator reproduces the draft's "identifiability gradient", and that gradient is an
   artefact of the comparator.** At 0.12 the informal score recovers 9–23% of the gap for two of the
   three coefficients and 61% for old-on-its-steep-side; tightening to 0.107 lifts everything.
   Reading down the four comparator columns, the apparent ordering `old ≫ new > avg` and its apparent
   threshold-dependence both vanish when the score is replaced. What survives from the old reading is
   only that the *comparator's* answer depends heavily on an analyst's threshold.
3. **old's asymmetry is real but far weaker than it looked.** Even under the primary rule old recovers
   slightly less from the strong side (88%) than from the weak side (112%), which is the residual
   trace of the asymmetric valley in Step 4b: the objective bounds how *weak* old can be more sharply
   than how *strong*. Under the informal score that same asymmetry is a 9% versus 61% chasm. So the
   one-sided-identifiability claim should be reported as a mild direction-dependence of the formal
   posterior, not as one-sided identifiability.
4. **The retention column is a diagnostic, not a result.** The primary rule's effective sample size
   is ~2% of the library in both designs. That is enough for the medians and IQRs quoted here (Step 1
   shows the quantiles stable from 4096 draws) but it is the reason no finer statistic is quoted from
   this step.

**Identifiability gradient — revised.** The controlled experiments no longer support `old ≫ new > avg`
as a statement about the data. Under the formal likelihood the three coefficients are informed to
comparable *relative* precision (SD retained 25 / 31 / 28%, Step 1) and all three recover from a
displaced prior. What differs between them is not whether the data speak but how *fragile* that
information is once realistic confounders are admitted, which is the subject of Steps 7 (Case C), 8b
and 8c. The sentence that should be reported is therefore about a gradient in **practical robustness**,
not in identifiability *per se*.

---



## Step 5 — structural-error experiment (breaks the inverse crime)



### Step 5a — mild within-zone heterogeneity (±20% pipe-level jitter)

Truth: each pipe `k_w,p = zone_mean × (1 + δ_p)`, `δ_p ~ U(-0.2, 0.2)`, seed 12345
(14 old / 37 avg / 66 new pipes). Fitted model: the same three-zone **homogeneous** model.
Same noise process (seed 42, σ=0.1, clip), same priors, **primary threshold 0.107**
(`step5_structural_error.py`; the draft's 0.12 gives the same qualitative result).

**This sub-section is a superseded first look**, retained because it is what the draft ran: one
heterogeneity field, the informal comparator, and no control. Step 5c replaces it with a 25-field
ensemble under the primary rule, and Step 5d with the structured design. Read the numbers below as the
configuration, not as the result.

- **Noise-free structural residual** (best homogeneous fit vs the heterogeneous truth) =
**0.0068 mg/L** — a small fraction of the realised window noise RMSE (0.0973). The within-zone
heterogeneity is
essentially invisible under σ = 0.1.
- Grid-search best fit (-1.067, -0.12, -0.037); 4613/8192 behavioural at thr 0.107, min RMSE 0.0973.


| zone    | true arith. mean | length-weighted | pipe range     | informal GLUE (comparator) ± SD | grid fit | bias (informal − arith) |
| ------- | ---------------- | --------------- | -------------- | ------------------------------- | -------- | ----------------------- |
| old     | -1.007           | -0.989          | [-1.17, -0.80] | -1.030 ± 0.260                  | -1.067   | -0.024                  |
| average | -0.098           | -0.096          | [-0.12, -0.08] | -0.116 ± 0.042                  | -0.120   | -0.018                  |
| new     | -0.049           | -0.047          | [-0.06, -0.04] | -0.048 ± 0.023                  | -0.037   | +0.001                  |


Risk top-6: informal `131,129,141,15,127,139`; true field `131,141,129,143,15,139` (same leading
group, minor reorder).

**Finding (null / robustness), with the caveat that one field cannot establish it.** ±20% *symmetric*
within-zone heterogeneity produces no meaningful structural error here: the effective grouped
coefficients match the true field averages to within a fraction of the posterior SD (biases −0.024 /
−0.018 / +0.001), the structural residual (0.0068) is swamped by the noise (0.0973 on the window),
and the risk map
is preserved. This is *not* the "precise but biased" pathology, because symmetric mean-zero jitter
averages out — old's length-weighted mean (−0.989) is essentially its arithmetic mean (−1.007), so the
fit has nothing to be biased toward. The claim only becomes safe in Step 5c, where the same test is
repeated over 25 independent fields and the field-to-field scatter turns out to be three times the
mean effect.

### Step 5c — heterogeneity magnitude sweep (0 control / ±20 / 35 / 50%)

Same design, jitter magnitude increased (`step5c_jitter_sweep.py`; the homogeneous candidate
predictions are truth-independent, so the 8192 cached candidates are reused and only the truth is
re-simulated). Bias = weighted mean − true per-zone arithmetic mean.

**Two corrections were needed here, and they pull in opposite directions.** The sweep now includes a
`jitter = 0` control, and it is reported under both weightings.

*Subtract the control.* At `jitter = 0` the truth is exactly homogeneous, so nothing measured there
can be structural — whatever bias appears is that realisation's noise. The raw bias column
therefore **overstates** the structural effect; `incr` (increment over the control) is the honest
part. The control also calibrates the structural residual: 0.0055 mg/L at `jitter = 0` is pure
7×7×7 grid spacing, because a homogeneous truth is exactly representable.

*Use the formal weighting.* The informal score is nearly flat inside the behavioural set, so its
weighted mean is pinned near the centre of the accepted box and cannot follow the data. It
**understates** the structural bias — by a factor of three to four in the old zone.


| jitter                         | struct. resid. | above grid floor | old bias (SD) / incr         | avg bias (SD) / incr     | new bias (SD) / incr     |
| ------------------------------ | -------------- | ---------------- | ---------------------------- | ------------------------ | ------------------------ |
| **formal censored (primary)**  |                |                  |                              |                          |                          |
| 0 (control)                    | 0.0055         | 0.0000           | +0.012 (+0.13) / —           | −0.009 (−0.64) / —       | +0.006 (+0.77) / —       |
| ±20%                           | 0.0068         | 0.0013           | −0.020 (−0.20) / **−0.0325** | −0.010 (−0.74) / −0.0012 | +0.007 (+0.93) / +0.0011 |
| ±35%                           | 0.0103         | 0.0048           | −0.033 (−0.33) / **−0.0458** | −0.011 (−0.78) / −0.0016 | +0.008 (+1.05) / +0.0020 |
| ±50%                           | 0.0145         | 0.0090           | −0.037 (−0.36) / **−0.0490** | −0.011 (−0.79) / −0.0016 | +0.009 (+1.18) / +0.0029 |
| **informal GLUE (comparator)** |                |                  |                              |                          |                          |
| 0 (control)                    | 0.0055         | 0.0000           | −0.011 (−0.04) / —           | −0.017 (−0.40) / —       | +0.001 (+0.03) / —       |
| ±20%                           | 0.0068         | 0.0013           | −0.020 (−0.07) / −0.0081     | −0.018 (−0.42) / −0.0009 | +0.001 (+0.03) / −0.0001 |
| ±35%                           | 0.0103         | 0.0048           | −0.022 (−0.09) / −0.0108     | −0.018 (−0.44) / −0.0015 | +0.001 (+0.05) / +0.0004 |
| ±50%                           | 0.0145         | 0.0090           | −0.025 (−0.10) / −0.0134     | −0.019 (−0.45) / −0.0021 | +0.002 (+0.10) / +0.0016 |


Bias vs within-zone jitter, against the homogeneous control

**What the single-field table above shows, and why it is not the finding.** Read at face value, the
formal column says symmetric heterogeneity biases the old-zone coefficient, with the increment
growing monotonically to −0.049 m/day (0.36 posterior SD) at ±50%, while the informal score puts the
same increment at −0.013 (0.10 SD). The informal-versus-formal contrast is real and is the same
inertia seen elsewhere. But every row uses **one** jitter field (seed 12345), so the increment mixes
the effect of heterogeneity with the accident of one spatial arrangement. That has to be separated
before anything is claimed.

### Step 5c-ensemble — 25 independent heterogeneity fields at ±20%

The candidate predictions are truth-independent, so repeating over many fields costs one truth
simulation each. 25 fields, formal weighting, control subtracted:

| zone    | bias mean | bias SD | structural increment mean | increment SD | |mean| / SD |
| ------- | --------- | ------- | ------------------------- | ------------ | ----------- |
| old     | +0.0231   | 0.0335  | **+0.0107**               | 0.0335       | **0.32**    |
| average | −0.0089   | 0.0046  | +0.0002                   | 0.0046       | 0.04        |
| new     | +0.0055   | 0.0021  | −0.0005                   | 0.0021       | 0.24        |

**This reverses the single-field reading, sign included.** The ensemble-mean structural increment for
`old` is **+0.0107** m/day — weaker decay — where the single field gave **−0.0325**, stronger decay.
The field-to-field scatter is 0.0335, three times the mean, and `|mean| / SD` is 0.32, 0.04 and 0.24
for the three zones. So at ±20% symmetric within-zone heterogeneity there is **no systematic
structural bias detectable above field-to-field variability**, and the single-field number is not
separable from the choice of field.

The conclusion to carry forward is therefore the original one: **symmetric, mean-zero within-zone
heterogeneity averages out** in this network at this noise level. What the corrections established
along the way still stands and still matters — the control must be subtracted, the informal score
does understate whatever effect is present, and the structural residual has a grid floor — but none
of that adds up to a precise-but-biased effect for *symmetric* heterogeneity. For that, the
heterogeneity has to be **structured** rather than mean-zero — correlated with some pipe attribute
that the grouped average does not preserve. Step 5d tests one such structure, using **pipe length**
as the correlating attribute.

A note on how this section was arrived at, because it bears on how the thesis should treat any
single-realisation result: the claim here has now been reversed twice, first by adding the `jitter = 0`
control and switching to the formal likelihood, then by replacing one field with 25. Both reversals
came from removing a hidden dependence on one arbitrary choice, and neither would have been visible
from the fit quality, which was excellent throughout.

The structural residual grows 0.0068 → 0.0145 mg/L, but only 0.0013 → 0.0090 of that is above the
grid floor, so the earlier reading "the residual is 7% of the noise floor" overstated it by about
fivefold; against the realised window noise RMSE of 0.0973 the true structural part at ±20% is
1.3%. The risk ranking
is unchanged across all magnitudes.

### Step 5d — structured (length-correlated) within-zone heterogeneity → precise-but-biased

Location-consistent design: the truth keeps the three **location** zones, but within each zone `k_w`
is **correlated with pipe length** (longer pipe → stronger decay, factor 1 + 0.5·s, s ∈ [−1,1] by
within-zone length rank). The per-zone **arithmetic** mean is held exactly at the zone mean while the
**length-weighted** mean is shifted stronger, so the two candidate targets of a three-zone homogeneous
fit are separated and can be told apart (`step5d_structured.py`; candidate library reused from cache).

Length weighting is an **illustrative proxy target**. It is *not* the residence-weighted or
hydraulically effective coefficient, and the two must not be equated — see the caveat below.

Primary rule = formal censored likelihood; the informal comparator is shown beside it. `hom. offset`
is the same rule's departure from the truth on the **homogeneous** baseline, i.e. what it does with no
structural error at all; subtract it before reading any shift as structural.


| zone    | true arith. | length-wt proxy | lenwt − arith | rule                | posterior mean ± SD  | shift       | hom. offset | frac of gap |
| ------- | ----------- | --------------- | ------------- | ------------------- | -------------------- | ----------- | ----------- | ----------- |
| old     | -1.000      | **-1.2404**     | -0.2404       | **formal censored** | **-1.1444 ± 0.1123** | **-0.1444** | +0.0124     | **60%**     |
| old     |             |                 |               | informal GLUE       | -1.0666 ± 0.2549     | -0.0666     | -0.0115     | 28%         |
| average | -0.100      | -0.1238         | -0.0238       | **formal censored** | **-0.1349 ± 0.0170** | **-0.0349** | -0.0091     | **147%**    |
| average |             |                 |               | informal GLUE       | -0.1300 ± 0.0416     | -0.0300     | -0.0170     | 126%        |
| new     | -0.050      | -0.0659         | -0.0159       | **formal censored** | **-0.0571 ± 0.0083** | **-0.0071** | +0.0060     | **45%**     |
| new     |             |                 |               | informal GLUE       | -0.0593 ± 0.0248     | -0.0093     | +0.0010     | 58%         |


Structured heterogeneity: the fit moves toward a length-weighted value

**Single-realisation reading (one noise draw).** The structural residual stays small (0.0063, so the
fit remains *precise*: RMSE 0.0968, just under the realised window noise RMSE of 0.0973), the
length-weighted proxy diverges from the
arithmetic mean (old −1.2404 vs −1.000), and under the primary rule the fitted coefficient moves off
the arithmetic mean in that direction by 60% (old), 147% (average, i.e. past the proxy) and 45% (new)
of the arithmetic → proxy gap. Two things are worth noting about the comparison with the comparator
column. The formal shift for old is **more than twice as large** (−0.1444 against −0.0666), so the
informal score understates the structural effect just as it understates sensor bias (Step 8). And the
formal rule's homogeneous-baseline offset is small (+0.012 / −0.009 / +0.006), which is what makes the
structural attribution legitimate at all: under the informal score a comparable fraction of each
"structural" shift is the score's own displacement.

These are still one noise draw. The dose-response below is what the effect looks like once the draw
is averaged over.

### Step 5d-dose — correlation strength × 30 noise realisations

Unlike Step 5c the truth here is deterministic (length rank order), so there is no heterogeneity
field to average over; the single-realisation risk is the noise draw. And if the effect is real it
must vanish at `CORR = 0`, where the two candidate targets coincide, and grow with `CORR`. Both are
tested at once. The statistic is the fraction of the arithmetic → length-weighted gap travelled: 0
means the fit sits on the arithmetic mean, 1 on the length-weighted value.


**Primary — formal censored likelihood:**

| `CORR` | gap (old) | old: median [5, 95]        | average              | new                  |
| ------ | --------- | -------------------------- | -------------------- | -------------------- |
| 0.00   | 0.0000    | control — targets coincide | control              | control              |
| 0.25   | −0.1202   | +1.10 [−0.29, +2.47]       | +1.26 [−0.95, +3.72] | +0.98 [−0.89, +2.66] |
| 0.50   | −0.2404   | +0.86 [+0.12, +1.46]       | +1.12 [−0.07, +2.43] | +0.92 [−0.03, +1.78] |
| 0.75   | −0.3606   | +0.74 [+0.23, +1.06]       | +1.07 [+0.22, +1.93] | +0.89 [+0.26, +1.46] |

**Comparator — informal GLUE** (medians only; the same runs, the same gaps):

| `CORR` | old  | average | new  |
| ------ | ---- | ------- | ---- |
| 0.25   | +0.65 | +2.00  | +0.89 |
| 0.50   | +0.43 | +1.22  | +0.78 |
| 0.75   | +0.34 | +0.97  | +0.72 |


**Finding, and it is stronger than the single draw suggested — under the primary rule.** The median
shift fraction is **0.74–1.26** at every correlation level and in every zone: the fitted coefficient
tracks the length-weighted proxy approximately fully, not partially. The single-realisation
60% / 147% / 45% sits inside these intervals, so it was a noisy draw of a quantity whose median is
near 1.

**The comparator does not reproduce this, and the difference is a factor of two for `old`.** Under
the informal score the old-zone median falls from +0.65 to +0.34 as `CORR` grows, i.e. it reads as
"the fit tracks only a third of the gap" exactly where the formal rule reads "essentially all of
it". Both rules report a FALLING fraction as `CORR` grows (formal +1.10 -> +0.86 -> +0.74; informal
+0.65 -> +0.43 -> +0.34), so the trend direction is the same, not opposite; and because the gap
itself widens (0.120 -> 0.240 -> 0.361 m/day for old) the ABSOLUTE displacement rises under both
(formal 0.132 -> 0.207 -> 0.268; informal 0.078 -> 0.104 -> 0.123). What the comparator does is
attenuate consistently, to 0.59, 0.50 and 0.46 of the formal fraction. This is the same inertia as Steps 3, 5c, 8 and 8b:
a score that is nearly flat inside its behavioural set cannot follow a systematic displacement. The
statement "the fit tracks the length-weighted proxy" is therefore a **formal-likelihood** statement
and must be quoted as one.

The intervals also behave as a real effect should. They narrow as `CORR` grows — [−0.29, +2.47] at 0.25 against [+0.23, +1.06] at 0.75 — because the gap in the denominator grows while the noise on the numerator does not. At `CORR = 0.75` all three intervals **exclude zero**, so at sufficient correlation the conclusion is resolved: the fit is inconsistent with sitting on the arithmetic mean and consistent with tracking the length-weighted value. At `CORR = 0.25` it is not resolved, which is the honest limit of what this design can detect.

Note the contrast with Step 5c, and it is the point of having both: symmetric heterogeneity produced no detectable bias across 25 fields, while length-**correlated** heterogeneity of the same magnitude produces a fully resolved one. **Structure in the heterogeneity, not heterogeneity as such, is what breaks the averaging-out.** The structure tested here is a synthetic *pipe-length* correlation, used because it cleanly separates the two candidate targets; it is **not** a flow-path, residence-time or Jacobian weighting, and this design cannot say which attribute a real network would correlate with.

**Length-weighted must not be called residence-weighted, and this bounds the claim.** The reaction weight of a pipe depends on its flow, direction, diameter, residence time and on how strongly the six monitors see it. Length is one proxy among those; it happens to separate the two candidate targets in this design, which is what makes the design useful, and nothing here shows it is the *right* target.
So the two statements available differ in strength and only the weaker one is supported:

- **Supported:** structured heterogeneity makes the homogeneous estimate move off the arithmetic mean, in the direction of a length-weighted proxy, by a fraction of the gap whose median is near 1 and whose interval excludes zero at `CORR = 0.75`.
- **Not supported:** that homogeneous fitting *recovers* the residence-weighted or effective coefficient. Establishing what the fit actually estimates needs a sensitivity- or Jacobian-weighted effective mean, which this step does not compute. That is open work, not a result.

**Honest caveat on magnitude.** Under the primary rule the old-zone offset of −0.1444 is 1.29 posterior SD (0.1123), so the structural shift is now comparable with the posterior width rather than buried inside it — a stronger and more concerning result than the 0.26 SD the informal comparator gave. What the step demonstrates is still the *mechanism*: the fit stays precise (RMSE at the noise floor) while the coefficient is pulled away from the field's simple average. Note also that the grid fit is a coarse instrument here: the old-zone grid step is 0.2170 m/day (half-step 0.1085), larger than the shift being measured, so `grid_fit` cannot resolve it and only the weighted mean can.

### Step 5e — grid-search recovery is "centred on the answer" (§3.3, second half)

The 7×7×7 grid ranges are chosen around the true values, so the nearest grid node to each truth is close by construction:


| coeff | grid step | nearest node to truth | distance | half-step |
| ----- | --------- | --------------------- | -------- | --------- |
| old   | 0.217     | -1.067                | 0.067    | 0.108     |
| avg   | 0.027     | -0.093                | 0.007    | 0.014     |
| new   | 0.016     | -0.052                | 0.002    | 0.008     |


The **noise-free** best homogeneous fit — `(-1.067, -0.093, -0.052)`, stored as
`structural_best_homogeneous` — is **exactly** these nearest nodes, so for a clean objective
"recovery to the nearest grid node" is guaranteed by the grid placement, not by the data. **The
noisy grid fits agree on `old` and disagree on the other two**, which is the more informative half
and was previously written as if every run returned the nearest-node triple: Step 5a records
`(-1.067, -0.120, -0.037)` and Step 5d records `(-1.067, -0.147, -0.052)`. Consequences:

- The deterministic grid **cannot distinguish identifiable from unidentifiable coefficients**, and it
fails in two opposite ways at once. For `old` it lands on the nearest node whatever the data say, so
a hit there is guaranteed by placement. For avg/new the fitted node instead wanders one to two steps
between runs (avg returns -0.093, -0.120 or -0.147 depending on the truth field and the noise draw),
because the objective is nearly flat across their whole prior (Step 4b) and the grid step is coarser
than anything that separates them. Either way the grid emits a triple that *reads* as a recovered
answer. This is why the draft's grid looked like successful recovery of all three, masking what GLUE
reveals.
- The draft's reported 10-realisation spread for old (±0.12 m/day) ≈ **half a grid step** (0.108),
i.e. largely quantisation, not genuine estimation uncertainty.

Fix for the draft: present the grid only as an implementation/plausibility check, and rely on GLUE
(with the displaced priors of Steps 2/4) for identifiability and uncertainty.

### Step 5 — overall

- **Uncorrelated within-zone heterogeneity (±20–50%, Step 5c)**: the grouped model is **robust**. The posterior mean stays at the arithmetic mean to within the field-to-field scatter, and the risk ranking is preserved. Across 25 independent fields at ±20% the structural increment is `+0.0107 ± 0.0335` (|mean|/sd = 0.32), i.e. undetectable.
- **Length-correlated within-zone heterogeneity (Step 5d)**: **precise but biased**. The fit stays at the realised window noise RMSE while every zone's coefficient moves toward the length-weighted proxy. Over 30 noise realisations the median fraction of the arithmetic → proxy gap travelled is **0.74–1.26 in every zone at every `CORR` level**, and at `CORR = 0.75` all three intervals exclude zero. The single-realisation 60% / 147% / 45% must **not** be quoted as the finding or used to order the zones: it is one draw of a quantity whose 5–95 interval at `CORR = 0.50` spans [0.12, 1.46] for old, [−0.07, 2.43] for average and [−0.03, 1.78] for new, and the zone that appears to pass the proxy changes between the draw and the median. On magnitude, the old-zone shift of −0.1444 is 1.29 posterior SD, i.e. comparable with the posterior width rather than buried inside it. What the fit converges to is a direction, not an identified effective coefficient.
- **Grid search (Step 5e)**: recovery to the nearest grid node is guaranteed by the grid being centred on the truth, so the deterministic "recovery of all three coefficients" is an artefact of grid placement; only a weighted ensemble can speak to identifiability.
- Directly answers §3.3 ("report how far the behavioural ensemble sits from any simple average of the
true field"): the distance is ≈0 when heterogeneity is uncorrelated and a clear, resolved shift
toward the length-weighted proxy when it is correlated. **Structure in the heterogeneity, not
heterogeneity as such, is what breaks the averaging-out** — tested here with a synthetic pipe-length
correlation, which is one structured proxy and not a demonstrated flow-path or residence weighting.
The low-residual clipping in the old zone is a further precise-but-biased source (→ Step 9).

---



## Step 6 — sensor-accuracy (noise) sensitivity  [answers the email]

Same three-zone truth; observation noise `σ = 0.02 / 0.05 / 0.10 / 0.15 mg/L` (σ = one standard
deviation; covers the supervisor's ±0.05 / ±0.10 / ±0.15). 30 noise realisations per σ, median [IQR];
reuses the baseline cache (candidates are noise-independent), no EPANET
(`step6_noise_sensitivity.py`).

**Why this section had to be re-run, and it changes the answer to the email.** "How accurate must a
sensor be?" is a question about what the data can support, so it has to be asked of an efficient
likelihood. Asked of the informal score — which is what the earlier version of this table did — it
returns the *score's* inefficiency as if it were a property of the instrument, and therefore demands a
better sensor than the measurement actually needs. Both rules are tabulated so the size of that
distortion is visible at every noise level.

**Primary — formal censored likelihood** (SD retained = posterior SD / prior SD, median [IQR]):


| σ (mg/L) | ESS med   | ESS p5 | ESS min | frac < 100 | old SD ret.     | avg SD ret.     | new SD ret.     | band @node15 |
| -------- | --------- | ------ | ------- | ---------- | --------------- | --------------- | --------------- | ------------ |
| 0.02     | 1.8       | 1.3    | 1.2     | 100%       | 4% [3–6]        | 3% [2–4]        | 3% [3–5]        | 0.006        |
| 0.05     | 20.9      | 17.2   | 16.8    | 100%       | 13% [12–13]     | 15% [14–15]     | 15% [14–15]     | 0.017        |
| **0.10** | **159.6** | 137.3  | 124.6   | 0%         | **27% [25–29]** | **30% [28–31]** | **29% [28–30]** | 0.034        |
| 0.15     | 525.7     | 405.0  | 371.8   | 0%         | 38% [37–40]     | 45% [42–48]     | 43% [42–44]     | 0.051        |


**Comparator — informal GLUE**, threshold scaling with σ as the ~95% acceptance band
(`threshold_for_sigma`, 0.107 at σ = 0.1):


| σ (mg/L) | ESS med | ESS p5 | ESS min | frac < 100 | old SD ret. | avg SD ret. | new SD ret. | band @node15 |
| -------- | ------- | ------ | ------- | ---------- | ----------- | ----------- | ----------- | ------------ |
| 0.02     | 37.0    | 3.8    | 2.0     | 97%        | 15% [12–17] | 17% [12–20] | 19% [15–21] | 0.020        |
| 0.05     | 699.1   | 94.4   | 32.0    | 7%         | 42% [32–48] | 47% [36–52] | 44% [37–52] | 0.053        |
| 0.10     | 4218.0  | 1285.7 | 1157.9  | 0%         | 65% [63–72] | 91% [83–94] | 87% [81–91] | 0.087        |
| 0.15     | 6482.6  | 5137.4 | 3743.6  | 0%         | 81% [77–88] | 98% [95–99] | 97% [95–98] | 0.126        |


Ratio of informal to formal posterior SD — how much of the "required accuracy" belongs to the score
rather than to the sensor: **×3.9 / 5.3 / 5.8** at σ = 0.02, **×3.3 / 3.2 / 3.0** at 0.05,
**×2.4 / 3.1 / 3.0** at 0.10, **×2.1 / 2.2 / 2.3** at 0.15 (old / avg / new).

Sensor-accuracy sensitivity

Findings:

1. **A σ = 0.10 sensor is already useful, and the earlier conclusion that it is not was an artefact of
   the scoring rule.** Under the primary rule a σ = 0.10 instrument retains 27 / 30 / 29 % of the prior width for old / average / new — all three well informed. The informal comparator on the same
   observations reports 65 / 91 / 87%, which is what produced the earlier "σ ≲ 0.05 is required"
   recommendation. The required accuracy is therefore milder than previously reported, and the earlier
   figure should not be quoted.
2. **Prediction uncertainty scales roughly linearly with σ under both rules**, but at different levels:
   the 5–95% band at the risk-governing old-zone node 15 grows 0.006 → 0.051 mg/L under the primary
   rule against 0.020 → 0.126 under the comparator.
3. **The two tightest rows are Monte-Carlo-limited under the primary rule, and this is the honest cost
   of a sharp likelihood on a fixed prior library.** At σ = 0.05 the median ESS is 20.9 and at σ = 0.02
   it is 1.8, i.e. essentially one member carries the posterior. Their "SD retained" figures of 13% and
   4% are therefore **not** measurements of what such a sensor would buy — with so few effective
   members the posterior width is *under*-estimated. Both rows are flagged `sampling_limited` in the
   artifact on the same ESS criterion the figure shades. Quantifying σ ≤ 0.05 properly is best done
   with a **likelihood-adapted design** (importance, sequential or adaptive sampling). More uniform
   draws would work too, but inefficiently, and it is worth being exact about which: because the ESS
   fraction is roughly constant, absolute ESS scales **linearly** with the library, so reaching the
   ESS ≈ 100 criterion needs about a **fivefold** library at σ = 0.05 and about a **fifty-fivefold**
   one at σ = 0.02. At the measured cost of the present library that is on the order of tens of
   minutes and a few hundred MB for the first, and hours and several GB for the second — expensive
   and wasteful rather than impossible. The earlier wording, that "brute force does not fix it", was
   too strong: uniform sampling does fix it, at a cost that rises with the sharpness of the
   likelihood, which is exactly the argument for adapting the sampler to it instead.
4. **The comparator's minimum ESS is the diagnostic that matters, not its median.** At σ = 0.05 its
   median ESS is 699 but its worst realisation is 32, and 7% of realisations fall below 100. A row
   summarised only by its median would look comfortably resolved while one realisation in fourteen was
   not.

**Required-accuracy conclusion (email), revised.** Under the primary formal likelihood a
**σ ≈ 0.10 mg/L** instrument — the common commercial class — already constrains all three grouped
coefficients to under a third of their prior width and supports the risk pattern; at **σ = 0.15** the
widths grow to 38–45% and calibration adds materially less. Tighter sensors would in principle help,
but *this* study cannot say by how much, because at σ ≤ 0.05 the answer is limited by the prior
sampling design rather than by the measurement. What must **not** be repeated is the earlier claim that
σ ≲ 0.05 is required for useful predictions: that was the informal score's inefficiency, measured in
sensor units. This replaces the imported foundation-notebook noise sweep flagged in §3.5.

---



## Step 7 — a-priori identifiability: Fisher information / CRLB

**Roles of the three analyses (keep distinct).** Fisher/CRLB (Step 7) is *a priori* — it uses only
the model, sensitivities and the assumed noise level, not any realised observation set, so it
*predicts* identifiability before calibration. Profile likelihood (Step 7b) uses the *observed*
data, so it is a *practical / a-posteriori* identifiability check (formal-likelihood validation),
not a prediction. **Naming, used strictly from here on:** the **formal likelihood-weighted ensemble**
is the primary empirical result and the **informal GLUE behavioural ensemble** is the comparator.
The bare phrase "baseline GLUE" is not used, because it named a configuration and an inference rule
at the same time and the two have since come apart. At the baseline model (k_b fixed, no bias) the
like-for-like chain therefore has four legs: **Fisher Case A → formal profile likelihood → formal
likelihood-weighted ensemble → informal GLUE behavioural ensemble**, the last as comparator only.
Cases B/C are expanded-model (realism) predictions.

**Two-layer reading (paper structure) — keep them separate; do NOT merge them to argue that the
informal GLUE behavioural ensemble is "realistic".**

- *Layer 1 — controlled-baseline identifiability* (identical conditions: k_b fixed, unbiased
sensors, independent Gaussian noise, only three k_w): Fisher **Case A** (a-priori prediction) ↔
formal **profile likelihood** (a-posteriori validation) ↔ the **formal likelihood-weighted
ensemble** (primary empirical result) ↔ the **informal GLUE behavioural ensemble** (comparator).
Cleanest conclusion: *the data carry information about all three; the formal
likelihood extracts it; informal GLUE is substantially broader and more prior-sensitive.*
- *Layer 2 — realism sensitivity*: Case **B** (+k_b), Case **C** (+6 monitor offsets), **AR(1)**
(Step 7c) and censored likelihood (Step 9) each *separately* relax one baseline assumption.
Conclusion: *ideal identifiability is fragile to real error sources — especially avg/new — while
old stays comparatively robust.*

These baseline assumptions are **valid by construction** for this controlled synthetic experiment
(k_b was fixed, noise was generated independent-Gaussian, no bias was injected), so **Case A is
exactly correct here** — they are only *optimistic for a real-network deployment*. So we write "valid
for the controlled baseline but optimistic for real-network application", **not** "the assumptions
are false".

**Normalization convention (so the two tables are read consistently).** Step 7 reports the **1σ**
ratio `CRLB SD / prior SD`; Step 7b reports **95%** intervals. For a uniform prior,
`prior SD = half-width/√3` and a 95% half-width `= 1.96·SD`, so the two conventions differ by
≈ `1.96/√3 ≈ 1.13`. old's `0.25` (1σ CRLB/prior-SD) and its continuous censored profile
`±0.1807 / 0.65 ≈ 0.28` (95% half-width / prior half-width, and exactly `0.25 × 1.13`)
thus describe the **same tightness expressed two ways**, not
literally the same ratio. Both the **absolute interval** and the **CRLB/prior-SD** are reported —
never only the normalized number.

Sensitivity `J[r,j] = ∂C_r/∂k_j` at the truth (central differences; r = 6 monitors × 49 h = 294),
Fisher `F = Jᵀ J / σ²`, and the CRLB (smallest achievable SD). The CRLB for the three k_w is
computed under nested models that add the reviewer's confounders, using the marginal (Schur)
information `F_pp − F_pn F_nn⁻¹ F_np` (`step7_fisher.py`, σ = 0.10; independent noise assumed —
AR(1) (#4) quantified separately in **Step 7c**):

- **A** k_w only (k_b known, no bias) — idealised best case
- **B** + k_b (k_b trades off with k_w — Priority-2 #5)
- **C** + 6 per-monitor offsets (systematic sensor bias — Priority-2 #1)

**Finite-difference steps are scale-dependent**, not shared: `H` = 0.02 / 0.005 / 0.0025 m/day for
old / average / new (`wq_common.FD_STEP`), i.e. 2 / 5 / 5% of each truth, plus 0.05 day⁻¹ for `k_b`.
A single shared step would have been a 2% perturbation of `old` and a 40% perturbation of `new`.

Direct sensitivity check (`step7_verify.py`; max |ΔC| at any monitor when one coefficient is moved
by ±H, others at truth): **old** moves only its own zone nodes 15/145, by 0.0075 over a 0.04 change;
**average** moves 209/231 by 0.0130 over 0.01; **new** (near-source) moves nodes across the whole
network by 0.0106 over 0.005. Per unit of coefficient that is **0.19 / 1.30 / 2.12 mg/L per m/day** —
new and average have the *highest* per-unit sensitivity (near-source pipes lie upstream of the whole
network), old roughly a tenth of new's.


| case                    | old        | average    | new        | condition # |
| ----------------------- | ---------- | ---------- | ---------- | ----------- |
| A — k_w only            | 0.25 ✅     | 0.29 ✅     | 0.29 ✅     | 215         |
| B — k_w + k_b           | 0.30 ✅     | 0.54 ✅     | 0.54 ✅     | 101         |
| C — k_w + k_b + offsets | **0.67 ✅** | **2.24 ❌** | **1.09 ❌** | 88          |


(values = CRLB SD / prior SD, a **1σ** ratio; **< 1** means the data could *in principle* narrow the
coefficient below its prior — a heuristic, not a strict theorem.)

**The formal ensemble's spread is locally consistent with the Case-A CRLB.** This is the check that
closes the GLUE-vs-Fisher argument as a *consistency* reading (not a frequentist-efficiency proof —
that needs the repeated-sampling test in Step 14):


| coefficient | Case-A CRLB | formal censored ensemble SD | ratio    | informal GLUE SD |
| ----------- | ----------- | --------------------------- | -------- | ---------------- |
| old         | 0.09473     | 0.09380                     | **0.99** | 0.26679          |
| average     | 0.01351     | 0.01416                     | **1.05** | 0.04240          |
| new         | 0.00796     | 0.00776                     | **0.98** | 0.02419          |


Case-A CRLB for old 0.09473 (see table). The empirical spread of the formally weighted ensemble sits
within 1–5% of that bound on this single realisation. The informal score gives 2.8–3.1× the bound on the
*same* observations. So the long-standing "GLUE is three times wider than the CRLB" gap was never
evidence that the data are uninformative, nor a bug in the Jacobian: it is the statistical
inefficiency of the score. Only the formal rows are like-for-like with a CRLB.

Fisher/CRLB with nuisance parameters

**Case A is the like-for-like benchmark** for *both* baseline ensembles — the formal
likelihood-weighted one and the informal GLUE comparator alike, since both fix k_b and assume no
sensor bias. Cases B and C are **realism sensitivity analyses**, not descriptions of either.

Findings:

1. **Idealised A — the data do contain information about all three** (CRLB/prior 0.25/0.29/0.29 < 1).
   In fact new/avg have *higher per-unit sensitivity* than old (near-source pipes lie upstream of
   the whole network; 2.12 and 1.30 against 0.19 mg/L per m/day). So the observations are
   not fundamentally devoid of information about avg/new.
2. **The A-vs-GLUE gap is entirely the score, and this is now demonstrated rather than argued.**
   The informal GLUE also fixes k_b and ignores bias, so its like-for-like benchmark is A; it gives
   0.71/0.92/0.88 against the bound's 0.25/0.29/0.29. The formal ensemble on the same data gives
   0.25/0.31/0.28 — the bound, to within 5%. The residual explanations previously offered for the
   gap (threshold, zero-clipping, non-linearity, finite sampling, old's one-sided response) are
   therefore second-order at most: dropping the factor `N` from the likelihood accounts for it, as
   Stedinger 2008 / Mantovan & Todini 2006 predict. It must **not** be attributed to k_b or bias,
   which are absent from both.
3. **Sensitivity ≠ identifiability.** new/avg respond more per unit, but their response *mimics* a
   global k_b change or a per-monitor offset (non-unique) and is easily confounded; old's response
   is smaller but *distinctive* (localised at 15/145 with a particular spatiotemporal shape), so it
   survives confounding. This is why old is the robustly identifiable one despite its lower raw
   sensitivity.
4. **Realism sensitivity (B, C).** Admitting k_b uncertainty (B) inflates the avg/new CRLB 1.8–1.9×
   against 1.2× for old; adding per-monitor bias (C) pushes **avg (2.24) and new (1.09) above
   their prior → practically unidentifiable**, while **old retains limited identifiability (0.67)**.
   Case C thus reproduces a *qualitatively similar* gradient to the informal GLUE comparator — but
   under an **expanded nuisance model neither baseline ensemble used**, so it corroborates the
   robustness of the gradient rather than "confirming" either baseline run.
   **Two quantities are easy to confuse here and they carry opposite signs** (`step7_verify.py`
   prints both side by side). The Jacobian-column cosine `F_ij/√(F_ii F_jj)` — how collinear two
   *sensitivity directions* are, a property of the design — is **+0.333 / +0.770 / +0.831** for
   old / average / new against `k_b`. The correlation of the *estimators*, normalised from
   `Cov(θ̂) ≈ σ²(JᵀJ)⁻¹`, is **−0.547 / −0.839 / −0.842**. It is the negative estimator correlation
   that carries the bulk–wall compensation story (a weaker fixed `k_b` buys a stronger `k_w`, which
   is exactly what Step 8b measures), and only it should be called a parameter correlation. An
   earlier version of this log reported the cosine as if it were the correlation.
5. **The condition number on the raw matrix is mostly a units artefact.** κ ≈ 215 for Case A
   looks like a badly ill-conditioned problem, but the three coefficients differ in scale by a
   factor of 20, so eigenvalues of the unnormalised information matrix mix magnitudes as much as
   geometry. Rescaling the parameters by their prior SD makes them dimensionless:




| form         | eigenvalues (Case A) | condition number |
| ------------ | -------------------- | ---------------- |
| raw          | 23906 / 5106 / 111   | **215**          |
| prior-scaled | 24.5 / 16.6 / 7.6    | **3.2**          |


   On the dimensionless matrix the spread is a factor of **3.2**, not 215, so the three coefficients
   are far more comparably informed than the raw number suggests. This also qualifies the framing the
   review quoted from the operational notebook — "eigenvalues spanning three orders of magnitude for
   three zonal coefficients" — which is a statement about the units, not about identifiability. The
   sloppiest direction, in prior-SD units, is `old −0.297, average −0.663, new +0.687`: essentially
   **average traded against new**, which is the same confounding every other diagnostic finds. Note
   also that the *smaller* κ for B/C does not mean better identifiability — κ is only the relative
   imbalance between directions, and adding nuisance parameters lowers the total information.
6. **The finite-difference step is not driving any of this.** A single shared `H = 0.02` would be 2%
   of the old truth but 40% of the new one, so the steps are now **scale-dependent** (0.02 / 0.005 /
   0.0025, i.e. 2 / 5 / 5% of each truth) and each coefficient was additionally re-differenced over a
   range of steps around its default:


| coefficient | default step | steps tested | step as % of truth | Case-A CRLB spread over the sweep |
| ----------- | ------------ | ------------ | ------------------ | --------------------------------- |
| old         | 0.02         | 0.005 – 0.04 | 0.5 – 4%           | **0.10%**                         |
| average     | 0.005        | 0.002 – 0.02 | 2 – 20%            | **1.00%**                         |
| new         | 0.0025       | 0.001 – 0.01 | 2 – 20%            | **1.14%**                         |


   The CRLB is flat to ~1% across a fourfold to tenfold change of step, so the derivatives are
   converged and the choice of step carries at most a ~1% error even at a 20% perturbation. The
   concern was legitimate and the answer is that it does not matter at this precision — but the
   scale-dependent steps are what make the per-coefficient Jacobian columns comparable in the first
   place, and the raw `max |ΔC|` figures quoted above must be read against each column's own step.

This delivers Priority-2 #2 (a-priori identifiability) and links to #1 (sensor bias), #5 (k_b) and
#4 (AR(1), Step 7c) as realism analyses.

**Recommended thesis statement (main identifiability line).**

> A-priori Fisher information and formal profile likelihood showed that, under the controlled
> baseline assumptions of fixed bulk decay, unbiased sensors and independent Gaussian errors, the
> six-monitor dataset contained identifiable information about all three grouped wall-reaction
> coefficients. The informal GLUE formulation produced substantially broader and more
> prior-sensitive behavioural distributions, indicating that its weighting and threshold did not
> fully exploit the information contained in the 294 residuals. Separate realism-sensitivity
> analyses showed that bulk-decay uncertainty, monitor-specific bias and temporal autocorrelation
> could substantially weaken this ideal identifiability, particularly for the average and new
> coefficients, while the old coefficient remained comparatively robust.

**What GLUE is (and is not).** The informal GLUE ensemble is a *conservative, method-dependent
behavioural envelope*, retained throughout this log as a **comparator** — the contrast with the
formal rule is itself a result. It is **not** the vehicle for risk propagation: every headline risk
and scenario number comes from the formal censored likelihood (Steps 10 and 12), and where GLUE
appears there it is reported alongside as the comparator. Nor is it a calibrated
representation of all real-world uncertainty sources (it omits k_b uncertainty, monitor offsets,
AR(1) covariance and censoring). Its breadth comes from the informal likelihood + threshold + prior

- finite sampling + non-linear parameter space, **not** from having modelled reality.

**Do not write:** "the assumptions are false" (they hold by construction here); "GLUE's width ≈
real-world identifiability" (unsupported — the nuisance terms are absent from the baseline run); or
"the formal analysis converged to the GLUE result" (Case C is a *different, expanded* model). **Do
write:** Case C produced a *qualitatively similar* identifiability gradient (old partial, avg/new
weak), corroborating the robustness of the gradient rather than confirming the informal GLUE
ensemble.

---



## Step 7b — profile likelihood (practical identifiability, a posteriori)

The third identifiability tool the reviewer asked for. **Role:** unlike Fisher/CRLB (a priori),
profile likelihood uses the *observed* data, so it is a *practical / a-posteriori* identifiability
check (formal-likelihood validation), **not** an a-priori prediction. For each coefficient, fix it
and **re-optimise the other two** (minimise SSE), giving the profile negative log-likelihood; the
95% single-parameter interval is `ΔNLL ≤ 1.92`. Formal likelihood `NLL = (N/2σ²)·RMSE²` (N = 294, σ = 0.1), read off a 21³ grid
(`step7b_profile.py`). Unlike the Step 4b *sweep* (other two fixed at the truth), the profile allows
parameter compensation.

**Two versions are reported, because the grid one is quantised.** The 21-point grid steps 0.0650
m/day in the `old` direction — coarser than several effects discussed elsewhere in this log (the
censoring shift is 0.012, structural increments ~0.01) — so a grid interval can only place its
endpoints on nodes *inside* the true interval, and is therefore biased narrow. The continuous version
minimises the two nuisance coefficients with Nelder–Mead at each target value and locates the
endpoints by Brent bisection on `ΔNLL − 1.92`; **7792** distinct EPANET evaluations, run for both
the censored and the iid likelihood.


| coefficient | grid 95%         | continuous 95%         | grid half-width | continuous half-width | change     | endpoint move   |
| ----------- | ---------------- | ---------------------- | --------------- | --------------------- | ---------- | --------------- |
| old         | [−1.110, −0.850] | **[−1.1720, −0.8106]** | ±0.1300         | ±0.1807               | **+39.0%** | 0.78 grid steps |
| average     | [−0.136, −0.088] | **[−0.1373, −0.0821]** | ±0.0240         | ±0.0276               | **+14.8%** | 0.45 grid steps |
| new         | [−0.057, −0.034] | **[−0.0599, −0.0294]** | ±0.0119         | ±0.0152               | **+28.2%** | 0.70 grid steps |


Primary continuous profile uses the **censored** likelihood; the iid continuous curve is retained as
a local-identifiability benchmark (censoring impact is isolated in Step 9). The continuous 95% interval for `k_w,old` is [−1.1720, −0.8106].

**The grid understated every interval, by ~15–39% under the primary rule.** Each endpoint moves
outward by a **sub-step** distance — 0.78 of a grid step for `old`, 0.45 for `average` and 0.70 for
`new` (the half-width change divided by that coefficient's own grid step) — and the direction is
**always outward**. Both properties are the signature of quantisation rather than of noise: a grid
cannot find a crossing that falls between nodes, so it stops at the last node inside, and the miss is
therefore bounded by one step. The fraction is **not** the same in all three, so it must not be
quoted as a single figure.

The truth −1.0 / −0.1 / −0.05 lies inside every interval under both versions. Only the continuous
intervals should be quoted, and the same caution applies to anything else read off that grid — which
is why the AR(1) widening factor for `old` in Step 7c (2.00×) should not be read as physical: its
independent-case denominator was one of these too-narrow grid intervals.

Profile likelihood: tight two-sided intervals matching the CRLB

Findings:

1. **Under the formal likelihood all three are tightly, two-sidedly identifiable** — the profile 95%
   intervals bracket the truth and **match the Fisher CRLB 95% intervals** (Case A). The two formal
   tools agree: the observations *do* identify all three (idealised, k_b fixed, no bias).
2. **Little parameter compensation**: the profile (re-optimising the other two) ≈ the sweep (others
   at truth), so the three coefficients are not strongly confounded under the formal likelihood at
   this idealised setting.
3. **old's valley is mildly asymmetric** (shallower on the strong side), the residue of its
   one-sided tendency; but the formal likelihood still bounds it two-sided. So the "one-sided"
   description from Steps 4/4b is a *GLUE/informal-likelihood* view — under the efficient formal
   likelihood old is two-sided, though its strong side is only weakly constrained.
4. **Far tighter than GLUE.** These formal intervals are much narrower than the GLUE behavioural
   distributions (which barely narrow avg/new). The gap is GLUE's informal likelihood dropping the
   N = 294 factor (→ too flat); the profile/Fisher use the full formal likelihood. This closes the
   third leg of Priority-2 #2 and confirms the Step 7 diagnosis of GLUE inefficiency.
5. As in Step 7, this is the **idealised** benchmark (Case A); admitting k_b uncertainty and sensor
   bias (Fisher Case C) is what renders avg/new *practically* unidentifiable.

---



## Step 7c — error autocorrelation AR(1) (Layer-2 realism) — Priority-2 #4

Both formal tools above assume **independent** hourly errors. Real telemetry is temporally
autocorrelated, so those intervals are a *lower bound*. Rather than mechanically multiplying every
CRLB by the scalar effective-sample-size factor, we rebuild the observation covariance with an AR(1)
structure and recompute the Fisher information (`step7c_ar1.py`, Case A = k_w only, ρ = 0.4):

**Covariance / Fisher:**

```
per-sensor 49×49 block    Σ_block[t,s] = σ² · ρ^|t−s|          (hourly AR(1))
full covariance           Σ = block_diag(6 sensor blocks)      (294×294, sensors independent)
Fisher (AR(1))            F_AR1 = Jᵀ Σ⁻¹ J        vs   independent  F = Jᵀ J / σ²
CRLB_j = √[ (F⁻¹)_jj ]
```

**Scalar reference (quick caveat):**

```
inflation ≈ √[ (1+ρ)/(1−ρ) ] = √(1.4/0.6) ≈ 1.53
N_eff     ≈ N · (1−ρ)/(1+ρ) = 294 · 0.6/1.4 ≈ 126     (of 294 hourly points)
```


| coefficient | CRLB indep | CRLB AR(1) | widening | CRLB/prior indep | CRLB/prior AR(1) |
| ----------- | ---------- | ---------- | -------- | ---------------- | ---------------- |
| old         | ±0.0947    | ±0.1406    | 1.48×    | 0.25             | 0.37             |
| average     | ±0.0135    | ±0.0201    | 1.49×    | 0.29             | 0.43             |
| new         | ±0.0080    | ±0.0114    | 1.44×    | 0.29             | 0.42             |


Findings:

1. **Each coefficient widens by its *own* factor** (1.44–1.49×), computed from the full covariance —
   *not* a single mechanical 1.53× applied to all. They land close to the scalar reference here
   because the temporal Jacobian shape is similar across the six sensors; the covariance route is
   the defensible one and would separate the factors more where sensitivities differ in time.
2. **AR(1) alone is a *modest* caveat.** Even at ρ = 0.4 the ideal identifiability survives
   (CRLB/prior-SD old 0.37, avg 0.43, new 0.42 — all still < 1). So autocorrelation *inflates* the
   idealised intervals ~1.5× but does **not** by itself destroy identifiability.
3. **The dominant realism weakening is k_b + offsets, not AR(1).** Compare: Case C (k_b + 6 offsets)
   pushes avg to 2.24 and new to 1.09 (unidentifiable; Step 7), an order of magnitude larger effect than
   AR(1)'s ×1.5. AR(1) is therefore reported as a quantified *caveat*, not the leading realism term.

**AR(1) profile likelihood — a coarse sensitivity, not a second primary interval.** Read the
intervals below as a **fixed-ρ, uncensored, grid-resolution sensitivity analysis**, not as 95%
intervals on the same footing as Step 7b's. Three things separate them: ρ = 0.4 is *assumed* and no
value of it is supported by these iid-generated observations; the objective is the plain Gaussian
`½·eᵀΣ⁻¹e`, **not** the censored likelihood that is primary everywhere else; and the endpoints come
off the 21³ grid, which Step 7b showed to be systematically narrow by 15–39%. A continuous censored
AR(1) profile has not been implemented, and until it is, only the *direction and rough size* of the
widening should be quoted. To avoid leaving the profile side "notional", the
21³ grid was rebuilt storing the full 294-residual vector `e` at every node (`step7c_profile_ar1.py`),
so the AR(1) profile could be computed exactly as `NLL(θ) = ½·e(θ)ᵀ Σ⁻¹ e(θ)` (the attached
operational notebook uses the same AR(1) observation process). The independent profile reproduces
Step 7b (sanity check), and the AR(1) profile widens by its own per-coefficient factor:


| coefficient | independent 95%  | AR(1) 95%        | half-width indep → AR(1) | widening |
| ----------- | ---------------- | ---------------- | ------------------------ | -------- |
| old         | [-1.110, -0.850] | [-1.240, -0.720] | ±0.130 → ±0.260          | 2.00×    |
| average     | [-0.136, -0.088] | [-0.144, -0.072] | ±0.024 → ±0.036          | 1.50×    |
| new         | [-0.057, -0.034] | [-0.067, -0.029] | ±0.012 → ±0.019          | 1.60×    |


The `old` factor of 2.00× should not be read as a physical result: the grid step in that direction
is 0.065 m/day, so each endpoint is quantised to about half the half-width being reported and the
ratio of two quantised numbers is coarse. `average` and `new` are on finer grids and their profile
widening (1.50–1.60×) agrees with the Fisher-CRLB widening (1.44–1.49×): both formal
tools react to AR(1) the same modest way, and the truth still lies inside every AR(1) interval — so
AR(1) inflates but does not overturn the idealised identifiability. (The scalar `√[(1+ρ)/(1−ρ)]≈1.53`
is only a reference; the per-coefficient factors above come from the full covariance.)

---



## Step 8 — systematic sensor bias (empirical, formal weighting) — Priority-2 #1

A constant offset is added to one informative old-zone monitor (node 15) and the calibration is
re-run (`step8_sensor_bias.py`; cache reused, 30 noise realisations, **formal censored weighting** —
measuring the damage with the informal score would understate it, as Step 5c shows for structural
error). The sweep is **two-sided**, because observations are censored at the sensor floor and a
negative offset is therefore not the mirror image of a positive one. Posterior SD of `k_w,old` with
no offset (the random spread) = **0.1012**.


| bias @ node 15 | old mean | old shift | shift / SD | avg mean | new mean | censored pts | ρ_Spearman | τ_Kendall | top-6 Jaccard |
| -------------- | -------- | --------- | ---------- | -------- | -------- | ------------ | ---------- | --------- | ------------- |
| −0.100         | −1.4274  | −0.3869   | −3.82      | −0.1017  | −0.0498  | 13           | 1.000      | 0.990     | 1.00          |
| −0.050         | −1.3028  | −0.2623   | −2.59      | −0.1033  | −0.0493  | 9            | 1.000      | 0.994     | 1.00          |
| −0.025         | −1.1745  | −0.1339   | −1.32      | −0.1033  | −0.0500  | 8            | 1.000      | 0.996     | 1.00          |
| 0.000          | −1.0405  | +0.0000   | +0.00      | −0.1028  | −0.0509  | 7            | 1.000      | 1.000     | 1.00          |
| +0.025         | −0.9231  | +0.1174   | +1.16      | −0.1020  | −0.0519  | 6            | 1.000      | 0.994     | 1.00          |
| +0.050         | −0.8192  | +0.2213   | +2.19      | −0.1013  | −0.0529  | 5            | 0.999      | 0.993     | 0.71          |
| +0.100         | −0.6486  | +0.3919   | +3.87      | −0.0993  | −0.0551  | 5            | 0.999      | 0.988     | 0.71          |


Rank columns compare the 92-node risk field against the unbiased case; `top-6 Jaccard` compares the
hot-spot sets.

Systematic bias pushes the calibrated coefficient

Findings:

1. **A systematic bias dominates the random spread by a factor of several.** A +0.05 mg/L offset at
   node 15 moves `k_w,old` from −1.0405 to −0.8192 — **2.19 posterior SD**, away from the truth −1.0
   — and +0.10 mg/L at node 15 moves `k_w,old` by **3.87 posterior SD**. The review's claim that "a +0.05 mg/L offset shifts the estimate by
   more than the entire random spread" is confirmed, and understated. A positive bias (sensor
   over-reads) makes the fit infer *weaker* decay.
2. **Only the locally-sensed coefficient moves.** `avg`/`new` shift by ≤0.004 across the whole ±0.10
   sweep — the bias corrupts `old`, not the others.
3. **The risk field is far more robust than the parameter.** At ±0.10, where `k_w,old` has moved
   almost 4 posterior SD, the 92-node risk ranking still has Spearman 0.999 and Kendall ≥ 0.988
   against the unbiased case, and this holds on **both** risk metrics (`P_bar` ≥ 0.9993, `E[A]`
   ≥ 0.9996). The leading set is more robust on the deficit metric Step 10 publishes: ranked by
   `E[A]` the top-6 is unchanged at every offset except +0.10 (Jaccard 0.71), whereas ranked by
   `P_bar` it already loses a node at +0.05. On `P_bar` the top-6 set is *unchanged* for every
   negative offset and loses one
   node at +0.05 and +0.10 (Jaccard 0.71). So an uncorrected sensor bias wrecks the coefficient while
   leaving the operational prioritisation nearly intact — which is reassuring operationally and a
   warning against reading the coefficient as physical.
4. **Shape — the review is right, and our first answer was wrong.** Halving the offset
   (0.05 → 0.025) retains **53%** of the shift: above the 50% that linearity would give, so the
   response is sub-linear/**concave**, as the review states, though weaker than the ~70% it cites.
   An earlier version of this section reported *convex/39%*. That figure came from the informal GLUE
   score, whose flatness inverts the curvature as well as shrinking the shift. This is the third
   place where the informal score reverses a conclusion, after Steps 3 and 5c.
5. **The response is asymmetric, and censoring is why.** `|negative| / |positive|` shift ratios are
   1.14, 1.19 and 0.99 at ±0.025, ±0.05 and ±0.10, while the number of observations sitting on the
   sensor floor runs 13 → 5 across the sweep. A negative offset pushes observations onto the floor,
   where a censored point carries different information from an unclipped one. Bias results must
   therefore always be quoted with their sign; a single ± figure is not meaningful.
6. **Empirical counterpart to Step 7 Case C.** Fisher showed per-monitor offsets *absorb* the signal
   (avg/new → unidentifiable); this run shows a real offset *biases* the most robustly informed
   coefficient. A real bias would be estimable and correctable by reference sampling — the
   operational notebook's QA step; this quantifies the **uncorrected** impact. Only a constant offset
   is modelled here, not drift.

**Why node 15, and does the location matter? (bias-location sweep,** `step8c_bias_bynode.py`**).** The
same offsets were injected at each of the six monitors in turn (two per zone: new 107/113, old
15/145, average 209/231), **two-sided** and under the **formal censored weighting**, 30 noise
realisations, risk field taken as the median over them. Baseline (unbiased) formal means are
old −1.0405 / average −0.1028 / new −0.0509 with posterior SDs 0.1012 / 0.0136 / 0.0080 m/day, and a
median of 7 observations already sitting on the sensor floor.

Own-coefficient shift in posterior SD (the coefficient of the biased monitor's own zone):

| biased node (zone) | own coef | −0.10     | −0.05 | +0.05 | +0.10 | largest cross-zone |Δ|/SD | top-6 Jaccard             |
| ------------------ | -------- | --------- | ----- | ----- | ----- | ------------------------- | ------------------------- |
| 107 (new)          | new      | −3.84     | −1.92 | +1.85 | +3.62 | 1.85 (average)            | 1.00 neg / 0.71 pos       |
| 113 (new)          | new      | **−4.44** | −2.24 | +2.08 | +3.97 | 2.02 (average)            | 1.00 neg / 0.71 pos       |
| 15 (old)           | old      | −3.82     | −2.59 | +2.19 | +3.87 | 0.51 (new)                | 1.00 neg / 0.71 pos       |
| 145 (old)          | old      | −3.31     | −1.89 | +1.59 | +2.95 | 0.80 (new)                | 1.00 neg / 0.71 pos       |
| 209 (average)      | average  | −1.89     | −0.79 | +0.60 | +1.10 | 2.66 (new)                | 1.00 except 0.71 at +0.10 |
| 231 (average)      | average  | **−5.94** | −3.25 | +2.72 | +4.22 | 0.74 (new)                | 1.00 throughout           |

Across all 24 formal rows the 92-node risk ranking stays at Spearman **≥ 0.999** against the
unbiased case on both metrics (minimum 0.9993 on `P_bar`, 0.9995 on `E[A]`), so the full-network
ordering never moves and only the leading set can.

The Jaccard column above is `P_bar`. **On `E[A]`, the metric Step 10's headline table is ranked by,
the top-6 is unchanged in 23 of the 24 arms** — the single exception is node 15 at +0.10, which
loses one node (Jaccard 0.71). So the shortlist is markedly more robust to sensor bias when ranked
by depth-weighted severity than when ranked by duration, and the four "0.71 pos" entries below are a
`P_bar` phenomenon.

Findings (location matters, and not the way the earlier informal table suggested):

1. **A bias mainly corrupts the coefficient of its own zone, but the leakage is not negligible.**
   The largest shift is always the own-zone one, yet biasing a *new*-zone monitor moves `average` by
   up to 2.0 SD and biasing an *average*-zone monitor moves `new` by up to 2.7 SD. In absolute
   m/day the leakage is small; in units of the affected coefficient's own spread it is not. This is
   the `average`↔`new` confounding that Step 7's sloppiest eigen-direction identifies, seen from the
   data side.
2. **The worst normalised corruption is not at node 15.** Measured in posterior SD the most damaging
   single sensor is **231 (average zone): −5.94 SD at −0.10 mg/L**, followed by 113 (new, −4.44 SD);
   node 15 reaches −3.82 / +3.87 SD. Node 15 remains the right *headline* case for Step 8 — `old` is
   the physically dominant and most robustly identifiable coefficient — but it is not the worst
   location, and the earlier claim that "avg/new hardly respond" was an artefact of the informal
   weighting, which flattens exactly the coefficients that are weakly informed.
3. **The two monitors within a zone differ** (113 > 107; 231 > 209; 15 > 145) — different local
   sensitivity and information content, so *which* sensor drifts, not just which zone, changes the
   magnitude by up to 3× (231 against 209).
4. **The response is sign-asymmetric**, and censoring is why: a negative offset pushes observations
   onto the sensor floor (median censored count rises to 13 at node 15, −0.10) while a positive one
   lifts them off it (down to 4–5). The negative shift is the larger one in eleven of the twelve
   node × magnitude pairs; the exception is node 15 at ±0.10 (−3.82 against +3.87), which is the same
   `|neg|/|pos| = 0.99` that Step 8's own sweep records at that magnitude. Bias results must be
   quoted with their sign.
5. **The risk product is robust in ordering but not in its top-6 membership.** Spearman never falls
   below 0.999 on either metric, and on `E[A]` only one of the 24 arms changes the leading set at
   all (node 15, +0.10). On `P_bar`: every *negative* offset leaves the top-6 set intact; the
   *positive* offsets at
   the four new/old monitors, and +0.10 at 209, each replace one hot-spot (Jaccard 0.71). Only
   node 231 leaves the set untouched at every offset. So a single drifting sensor can reorder the
   margin of the hot-spot list without disturbing the network-scale pattern.



---



## Step 8b — k_b ± 20% (empirical partner to Fisher Case B; Priority-2 #5)

Observations are generated at the true `k_b = -0.5`; the calibration is re-run with `k_b` fixed at
-0.4 / -0.5 / -0.6 (±20%) to see how the three grouped `k_w` are pushed by the bulk–wall trade-off
(`step8b_kb_sensitivity.py`; k_b=-0.5 reuses the cache, -0.4/-0.6 simulated; 30 noise realisations,
**formal censored weighting as primary**, informal GLUE retained as comparator). The risk field is
the **median over the 30 realisations**, so the ranking columns are not a single-seed result.

**Formal censored (primary).** Shift is against the `k_b = -0.5` row; `/SD` normalises it by that
row's own posterior SD (old 0.1012, average 0.0136, new 0.0080 m/day).


| k_b  | ESS | k_w,old (shift, /SD)         | k_w,avg (shift, /SD)         | k_w,new (shift, /SD)         | ρ_Spearman | top-6 Jaccard |
| ---- | --- | ---------------------------- | ---------------------------- | ---------------------------- | ---------- | ------------- |
| -0.4 | 188 | -1.1078 (-0.0672, **-0.66**) | -0.1236 (-0.0208, **-1.53**) | -0.0635 (-0.0126, **-1.58**) | 0.975      | **0.50**      |
| -0.5 | 160 | -1.0405 (ref)                | -0.1028 (ref)                | -0.0509 (ref)                | 1.000      | 1.00          |
| -0.6 | 135 | -0.9835 (+0.0570, **+0.56**) | -0.0836 (+0.0192, **+1.41**) | -0.0390 (+0.0119, **+1.49**) | 0.934      | **0.50**      |


**Informal GLUE (comparator, threshold 0.107).** Same runs, same risk field, different weighting:


| k_b  | k_w,old (/SD)   | k_w,avg (/SD)   | k_w,new (/SD)   | ρ_Spearman | top-6 Jaccard |
| ---- | --------------- | --------------- | --------------- | ---------- | ------------- |
| -0.4 | -1.0741 (-0.15) | -0.1286 (-0.28) | -0.0613 (-0.40) | 0.975      | 0.50          |
| -0.5 | -1.0368 (ref)   | -0.1170 (ref)   | -0.0517 (ref)   | 1.000      | 1.00          |
| -0.6 | -1.0086 (+0.11) | -0.1049 (+0.29) | -0.0422 (+0.40) | 0.946      | 0.50          |


k_b ±20% shifts the k_w estimates (bulk–wall compensation)

Findings:

1. **Bulk–wall compensation is confirmed empirically.** A weaker fixed k_b (−0.4) makes the fit
   infer *stronger* wall decay (k_w more negative); a stronger k_b (−0.6) makes it *weaker*. Any
   error in the fixed k_b is absorbed into the k_w estimates (Priority-2 #5).
2. **Magnitude — and this is where the weighting changes the verdict.** Under the formal likelihood
   ±20% k_b moves `average` and `new` by **1.4–1.6 posterior SD** and `old` by 0.56–0.66 SD. Under
   the informal score the same runs give 0.28–0.40 and 0.11–0.15 SD, i.e. an understatement of
   **roughly four- to fivefold** (3.9× for `new`, 4.4× for `old`, 5.5× for `average`), because a flat
   score cannot follow the data. This is the same inertia documented
   in Steps 3, 5c and 8, and it is why the primary row must be the formal one.
3. **The ordering across zones matches Fisher Case B.** Marginalising `k_b` inflates the CRLB by
   1.2× for `old` but 1.8× for `average` and `new` (Step 7), and it is `average`/`new` that move
   most here. The weakly-informed coefficients are the ones that absorb a bulk-decay error.
4. **Two risk metrics, and they do not agree — the claim has to name one.** The step8-family risk
   field is `P_bar`, the likelihood-weighted expected **fraction of the window below 0.2 mg/L**
   ("how long"). Step 10's headline top-10 is ranked by `E[A]`, the expected **cumulative deficit**
   ("how long *and* how far below"). Both are now computed here from the same cached predictions:

| `k_b` | ρ(P_bar) | top-6 J(P_bar) | ρ(E[A]) | top-6 J(E[A]) |
|---|---|---|---|---|
| −0.4 | 0.975 | **0.50** | 0.980 | **1.00** |
| −0.6 | 0.934 | **0.50** | 0.935 | **0.71** |

   The network-scale rank correlation is the same to three decimals on either metric. The **leading
   set is not**: on `P_bar` it loses two of six at both signs; on `E[A]` it is **unchanged** at
   −0.4 and loses **one** node at −0.6.

5. **On `E[A]` the change is a boundary swap, and on `P_bar` it is not.** Depth diagnostics on both:

| metric | `k_b` | top-3 | top-5 | top-6 | top-10 | top-15 | max \|Δrank\| in ref top-15 |
|---|---|---|---|---|---|---|---|
| P_bar | −0.4 | **0.20** | 0.43 | 0.50 | 0.54 | 0.67 | 16 |
| P_bar | −0.6 | **0.20** | 0.43 | 0.50 | 0.82 | 0.875 | 10 |
| E[A] | −0.4 | 0.50 | **1.00** | **1.00** | 0.67 | 0.67 | 11 |
| E[A] | −0.6 | **1.00** | **1.00** | 0.71 | 0.67 | 0.76 | 7 |

   On `E[A]` at −0.6 the *only* change to the leading set is node **143 (rank 6) leaving and node 129
   (rank 7) entering** — a pure sixth/seventh-place exchange, with top-3 and top-5 identical. At −0.4
   no node enters or leaves at all. So on the metric Step 10 publishes, the turnover **is** confined
   to the prioritisation boundary and the leading core is retained.
   On `P_bar` the same perturbation does something different: Jaccard is *worse* at k = 3 (0.20, one
   of three shared) than at k = 6, the reference 3rd node (129) falls to 8th and the 4th (127) to
   **20th**, and at −0.6 node 243 rises from 12th to 2nd. That is not a boundary effect.
   **An earlier version of this section asserted the P_bar reading as if it held generally. It does
   not; it holds for `P_bar`.**

6. **Why the two metrics differ, mechanically.** `P_bar` counts hours below the threshold and is
   therefore nearly saturated for the deepest nodes and highly sensitive to shallow, long excursions;
   `E[A]` weights by depth, so a node that dips just below 0.2 for many hours scores low. Changing
   `k_b` moves the whole concentration field slightly, which reshuffles the many nodes sitting near
   the threshold — visible in `P_bar`, largely invisible in `E[A]`, where the leaders are separated
   by their depth. Node 131 leads both metrics under every `k_b` and cannot move on `P_bar`: it is
   below the threshold for the entire window (`P_bar = 1.0000`, saturated).

7. **A cut-off plateau exists on `P_bar` and explains part of its turnover.** The reference `P_bar`
   profile is nearly tied across ranks 5–7 — `0.4585` (15), `0.4584` (143), `0.4583` (139), a spread
   of 0.0002 — so the top-6 cut-off falls inside the tie and those three are interchangeable under
   any small perturbation. It does not account for 127 falling to 20th, for 129 falling from 2nd to
   8th, or for 243 rising from 12th to 2nd.

8. **What is stable on both metrics is the bulk of the network.** The median rank change over all 92
   junctions is 1.0 place, which is why Spearman stays at 0.975 / 0.934 (`P_bar`) and 0.980 / 0.935
   (`E[A]`) while the leading sets move: the 86 low- and mid-risk junctions barely move and dominate
   the statistic.

9. **What to write.** `k_b` misspecification of ±20% shifts the average and new coefficients by
   1.4–1.6 posterior SD and leaves the **network-scale** risk ordering broadly rank-correlated on
   both metrics. Its effect on the **operational shortlist depends on the metric**: on the
   time-averaged low-chlorine probability `P_bar` the top-6 loses two of six at both signs and the
   reordering reaches well beyond the cut-off; on the cumulative-deficit ranking `E[A]` that Step 10
   publishes, the top-6 is unchanged at −0.4 and loses one node at −0.6 in a sixth/seventh-place
   exchange. Neither the claim that the risk product is insensitive to bulk-decay misspecification,
   nor the claim that bulk-decay misspecification changes the hot-spot list, is sayable without
   naming the metric it is ranked by.

10. **Caveats on this diagnostic.** Both rankings come from the median risk field over the 30
    realisations, so they inherit whatever the median smooths away; per-realisation rank stability is
    not reported. The near-ties at `P_bar` ranks 4–6 mean the *rank* of those three nodes carries
    little information in the first place — an argument for quoting the risk value alongside the rank
    whenever a shortlist is published. Steps 8, 8c and 8d now report both metrics as well, and the
    gap that opened here closes in the reassuring direction: on `E[A]` the sensor-error family barely
    touches the leading set (23 of 24 arms unchanged in Step 8c; no drift arm at all in Step 8d),
    so **`k_b` misspecification is the only perturbation tested in this log that moves the
    deficit-ranked shortlist**, and even then by one boundary exchange at one sign.


---



## Step 8d — sensor DRIFT, a time-varying offset — Priority-2 #1, second half

Steps 8 and 8c model a **constant** bias. The review's comment is headed "systematic sensor error
(bias, **drift**)", and a drifting probe — fouling, or a calibration that ages — is the more common
field failure. This step closes that half (`step8d_sensor_drift.py`; formal censored weighting,
30 noise realisations, cache reused so no EPANET is re-run).

**Ranking basis, stated because two are possible.** The top-6 set and the Spearman/Kendall columns
are both read off the **median risk field over the 30 realisations**. The alternative — the modal
top-6 *ordering* across realisations — answers a different question and is stored separately as
`modal_top6_across_realisations`. Earlier versions of this step and of Step 8 computed the rank
correlation from the median field but the Jaccard from the modal set, while declaring the median
basis for both; that is fixed, and the two turn out to agree on the **set** in all 24 arms here (and
all 7 in Step 8), differing only in the internal order of one arm. So no number in this section
changed — but a table whose two rank columns rest on different bases is not defensible even when
they happen to coincide.

**Formula** — drift is a linear ramp across the assessment window:

```
b(t) = D · (t − t₀) / 48,        t ∈ [120, 168] h
```

**Properties**:

- The sensor reads correctly at the window start and is off by `D` at its end, so its **mean** offset
over the window is `D/2`.
- A drift has both a mean and a shape. To separate them, every drift arm is run against **two
constant-bias controls on identical noise**: `const(D/2)`, the *mean-equivalent*, and `const(D)`,
the *end-equivalent*. If a drift lands on `const(D/2)` then only its mean matters and Step 8
already covers it; if it lands nearer `const(D)` the time structure matters in its own right.
- Two-sided, for the same censoring reason as Step 8. Two locations: node 15 (old, the Step 8
headline) and node 231 (average, where Step 8c found the largest normalised corruption).

Own-coefficient shift in posterior SD (baseline SDs old 0.1012, average 0.0136 m/day):


| biased node   | D      | drift (0→D) | const(D/2) | const(D) | drift / const(D/2) | censored pts (drift) |
| ------------- | ------ | ----------- | ---------- | -------- | ------------------ | -------------------- |
| 15 (old)      | −0.100 | −2.34       | −2.59      | −3.82    | **0.90**           | 10                   |
| 15 (old)      | −0.050 | −1.17       | −1.32      | −2.59    | **0.89**           | 8                    |
| 15 (old)      | +0.050 | +1.05       | +1.16      | +2.19    | **0.90**           | 6                    |
| 15 (old)      | +0.100 | +1.99       | +2.19      | +3.87    | **0.91**           | 6                    |
| 231 (average) | −0.100 | −3.19       | −3.25      | −5.94    | **0.98**           | 7                    |
| 231 (average) | −0.050 | −1.52       | −1.55      | −3.25    | **0.98**           | 7                    |
| 231 (average) | +0.050 | +1.41       | +1.44      | +2.72    | **0.98**           | 7                    |
| 231 (average) | +0.100 | +2.68       | +2.72      | +4.22    | **0.99**           | 7                    |


Sensor drift against its mean- and end-equivalent constant-bias controls

**Built-in cross-check.** The `const(D)` column reproduces Steps 8 and 8c exactly — node 15 at
−0.10 gives −3.82 SD and at +0.10 gives +3.87 SD; node 231 at −0.10 gives −5.94 SD — and
`const(D/2)` at each `D` equals `const(D)` at `D/2` (node 15: −2.59 both ways). The three step
scripts therefore agree on the shared conditions rather than merely being consistent in narrative.

Findings:

1. **A drift is, to a good approximation, its own mean.** The drift shift is **0.89–0.91** of the
   mean-equivalent constant shift at node 15 and **0.98–0.99** at node 231 (the ratio reaches
   **0.9851** at D = +0.10), and roughly half the
   end-equivalent one. So a constant-bias analysis is not a different phenomenon from a drift: to
   within 1–11% it *is* the drift analysis, evaluated at the drift's mean. Step 8's sweep therefore
   transfers, and the comment's "bias, drift" pairing is answered rather than half-answered.
2. **Where it departs from its mean, censoring is why.** Node 231's censored count is a flat 7 across
   every arm and its ratio is 0.98–0.99; node 15's runs 6 → 10 → 13 across the arms and its ratio
   drops to 0.90. A ramp spends only part of the window at the extreme offset, so it clips fewer
   observations than the constant bias of the same end value — and at an old-zone monitor, where the
   concentration is already near the sensor floor, that changes how much information the likelihood
   gets. The residual sits at 0.11–0.25 posterior SD for node 15 against 0.03–0.06 SD for node 231.
3. **The departure is signed, not symmetric.** For negative `D` the drift does *less* damage than its
   mean-equivalent (residual +0.15 to +0.25 SD toward the unbiased value); for positive `D` it does
   *more* (−0.11 to −0.19 SD). Same asymmetry, same cause, as Steps 8 and 8c.
4. **The risk product behaves exactly as under constant bias.** Spearman against the unbiased case is
   ≥ 0.999 in all 24 arms on both metrics. **On `E[A]` no drift arm changes the leading set at all** —
   the only 0.71 in the whole table is the end-equivalent *constant* control at node 15, +0.10 — so a
   drifting sensor of this size does not touch the shortlist Step 10 publishes. On `P_bar`
   the top-6 set is unchanged everywhere except node 15 at +0.10 drift
   (Jaccard 0.71) — one node replaced, the same margin effect the positive constant offsets produce.
5. **What this does and does not model.** Only a *linear* ramp within the assessment window is
   tested, with a drift that starts at zero. A drift already part-developed at the start of the
   window is the `const(D₀) + ramp` case and is not run; nor is a step change, a non-monotone
   excursion, or a fouling curve. The claim supported is the narrow one: **for a monotone drift, the
   mean offset over the calibration window is what the estimate responds to.**

---



## Step 9 — sensitivity to zero-floor censoring (zero-clipped observations, L = 0)

The draft generates observations as `C_obs = max(0, C_true + ε)` — a hard lower bound at **L = 0** —
then calibrates with unweighted RMSE that treats every clipped `0` as an *exact* measurement. The
statistically consistent treatment of a clipped point is *left-censored*: the latent reading
`Y* = μ + ε` was only observed to be `Y* ≤ 0`, so it should contribute `Φ(−μ/σ)` to the likelihood,
not `(0 − μ)²`. This is the robustness check the review asked for — performed **on the same data at
L = 0** (`step9_zeroclip.py`).

**Tobit-type censored Gaussian likelihood (L = 0)** — a censored Gaussian on the nonlinear WNTR
predictions (hence "Tobit-*type*", not a classical linear Tobit regression):

```
uncensored (obs > 0):  Gaussian    →  −½·((obs − μ)/σ)²
clipped-0  (obs = 0):  P(Y* ≤ 0)   →  log Φ(−μ/σ)          [scipy log_ndtr]
```

Since σ is fixed, minimising the Gaussian NLL orders the parameters identically to minimising
SSE/RMSE, so the "naive" arm is the formal-Gaussian counterpart of the draft's RMSE treatment.

**Clipped-zero census (the actual calibration data):**


| scope                               | clipped zeros                               | note                            |
| ----------------------------------- | ------------------------------------------- | ------------------------------- |
| full record (6 × 169 h)             | **43 / 1014**                               | the whole 168 h simulation      |
| calibration (6 × 49 post-warm-up h) | **10 / 294**                                | the points actually used to fit |
| by node                             | 15 (old): 6, 145 (old): 3, 231 (average): 1 | old 9, average 1, new 0         |


At the corrected warm-up the clipping is no longer confined to the old zone: one point now falls at
an average-zone monitor. The concentration in the average zone has dropped far enough by 120 h to
reach the sensor floor occasionally, which the 24 h window never showed.

**Naive-exact-0 vs censored-at-0 (median [IQR] over 30 noise realisations):**


| coefficient | truth | naive median [IQR]         | censored median [IQR]      | Δ median |
| ----------- | ----- | -------------------------- | -------------------------- | -------- |
| old         | −1.0  | −1.0289 [−1.1180, −0.9700] | −1.0405 [−1.1324, −0.9850] | −0.0116  |
| average     | −0.1  | −0.1027 [−0.1122, −0.0952] | −0.1028 [−0.1122, −0.0953] | −0.00003 |
| new         | −0.05 | −0.0510 [−0.0566, −0.0434] | −0.0509 [−0.0564, −0.0433] | +0.00007 |


- node-15 fraction of hours below 0.2: naive `0.4585` vs censored `0.4585` (identical);
- high-risk top-3 — ranked by expected time-fraction below 0.2 mg/L under the **formal-likelihood-
weighted ensemble** over all 8192 candidates: **131, 141, 129 — identical** under both estimators,
and now also identical to the primary risk map of Step 10, because both use the formal weighting;
- old profile 95% interval (baseline): naive `[−1.110, −0.850]` = censored `[−1.110, −0.850]`.

The profile intervals above being *exactly* equal is a **resolution artefact of the grid**, not a
null result: the 21-point grid steps 0.065 m/day in the `old` direction, five times the −0.0116
shift, so no grid-based interval can register it. The **continuous** profile of Step 7b, which is
not quantised, does resolve it — run under both likelihoods on the same data it gives
`old` = [−1.1720, −0.8106] censored against [−1.1517, −0.7973] iid, i.e. the whole interval moves
about 0.013–0.020 m/day toward stronger decay and the half-width grows from ±0.1772 to ±0.1807
(+2.0%). `average` and `new` move by well under a percent of their half-widths. So the correct
statement is: **grid-based profiles cannot resolve the censoring correction, whereas the continuous
profile shows a modest shift and a slight widening, concentrated in** `old` — consistent in sign and
size with the weighted-mean comparison below.

Zero-clipping census (10/294) and the old profile at L=0: naive vs censored overlap

Findings:

1. **10 of 294 calibration points are clipped, and they are concentrated in — but no longer confined
   to — the old zone**: 9 at the old-zone monitors (15 and 145) and **1 at average-zone monitor 231**
   (43 of 1014 over the full record). At the draft's superseded 24 h warm-up all 8 clipped points were
   in the old zone; at 120 h the average zone has depleted far enough to reach the sensor floor
   occasionally. The new zone still never clips.
2. **The censoring correction is small but measurable, and it acts on** `old` **only.** The censored
   likelihood moves the old median from −1.0289 to −1.0405, a shift of **−0.0116** m/day, and leaves
   `average` and `new` unchanged to within 0.0001. Node-15 risk (0.4585) and the top-3 risk ranking
   (131, 141, 129) are identical under both. So the draft's `max(0, ·)` + unweighted-RMSE choice did
   not bias the main conclusions — but "did not bias" is not "had no effect", and the effect it did
   have is the size of several other quantities discussed in this log.
3. **The *grid* profile intervals being exactly equal is a resolution artefact, not a null result.**
   Both estimators give `[−1.110, −0.850]` for `old` on the 21-point grid, which reads as perfect
   agreement. It is not: the grid steps 0.065 m/day in that direction, **five times** the −0.0116
   shift, so no grid-based interval could register it. The continuous profile (Step 7b) is not
   quantised and does register it — `old` = [−1.1720, −0.8106] censored against [−1.1517, −0.7973]
   iid, a shift toward stronger decay with a 2.0% wider half-width. Quote the weighted-mean
   comparison and, if a profile is quoted at all, the continuous one; never present the coincident
   *grid* intervals as evidence of equivalence.
4. **Why the effect is small here — and why that is not a general equivalence.** Two reasons specific
   to this setup: (i) clipping is rare (`10/294 ≈ 3.4%`, so 284 uncensored points dominate the
   likelihood); (ii) a clipped-0 point sits where the true concentration is already ≈ 0, so "exactly
   0" and "≤ 0" constrain a near-0 model `μ` in almost the same direction. This is an **empirical**
   result of the present experiment, **not** a theorem. A bias would appear if a larger fraction were
   censored, or if the recorded limit were *well above* the true value — i.e. a positive reporting
   limit. Note also that the censoring fraction is **warm-up dependent** (8 points at 24 h, 10 at
   120 h) and would keep growing with a longer horizon, so it is not a fixed property of the design.
5. **Reporting (≈ half a page).** State it as a robustness check with the correction quantified: the
   clipping is rare, concentrated in the low-residual old zone, and shifts the old estimate by
   −0.012 m/day without changing node-15 risk or the risk ranking. Say that the censored likelihood
   is nevertheless what the analysis uses — it is the primary weighting throughout this revision, not
   a side check — and that a censoring-aware likelihood would become necessary if a future dataset
   clipped many more points or used a positive reporting limit.

**Draft paper text (≈ half a page) — for Step 11:**

> *Method.* Because the synthetic observations were constrained to be non-negative, observations
> recorded at zero were treated as left-censored at the imposed boundary (`L = 0`). Positive
> observations contributed a Gaussian density term, whereas zero-clipped observations contributed
> the cumulative probability `P(Y* ≤ 0) = Φ(−μ(θ)/σ)`. This Tobit-type censored Gaussian likelihood
> was compared with the baseline treatment in which zero observations were entered as exact values.
>
> *Results.* Ten of the 294 post-warm-up calibration observations were clipped to zero — six at
> node 15, three at node 145 and one at node 231 (43 of 1014 over the full record), so the censoring
> is concentrated in, but not confined to, the low-residual old zone. Across 30 independent noise
> realisations the median old-group estimate changed from −1.0289 to −1.0405 m/day when the censored
> likelihood replaced the exact-zero treatment, while the average- and new-group estimates were
> unchanged; the low-chlorine summary and the risk ranking were identical.
>
> *Discussion.* The zero-floor treatment therefore had negligible influence on the principal
> calibration and risk conclusions in this synthetic experiment, reflecting the small number of
> clipped observations and their occurrence under already low-residual conditions. A censoring-aware
> likelihood would nevertheless become important if a larger fraction of observations were censored,
> or if the instrument had a positive reporting limit.



---



## Step 10 — operational risk metrics (duration + depth) and water-age corroboration

The risk map so far reduced to a single probability of being below the **selected operational
low-chlorine threshold** `C_MIN = 0.2 mg/L` (a representative operational value adopted in this
study, **not** a legal/compliance safety limit). A single probability conflates *how long* and *how
far below* the threshold a node sits. From the ensemble weighted by the formal censored likelihood
(with the informal GLUE score reported alongside as a comparator) we therefore report the original probability re-expressed as a duration,
**two** genuine severity metrics (minimum concentration, cumulative deficit) and **one**
reaction-independent hydraulic diagnostic (water age) — i.e. *one operational re-expression + two
severity metrics + one hydraulic diagnostic*, not "three new metrics" (`step10_risk_metrics.py`).

**Time axis.** The post-warm-up record is `t = 120 … 168 h` = **49 reporting points but
48 one-hour intervals**. Durations and deficits are **trapezoidally integrated over the 48 intervals**
(max duration = 48 h), not summed over 49 points.

```
duration D_i,n  = ∫ 1[C<0.2] dt          (trapezoid, h;  max 48)   frac = D̄/48 ≈ original P(C<0.2)
deficit  A_i,n  = ∫ max(0, 0.2−C) dt      (trapezoid, mg/L·h)       ← genuinely new (severity)
min C    M_i,n  = min_t C_i,n(t)          (per-member min, THEN weighted = Method A)  ← genuinely new
water age       = EPANET AGE run (reaction-independent; hydraulic diagnostic)
```

Each metric is a behavioural-weighted expectation `X̄_n = Σ_i w_i X_i,n`; 5–95 % bands are weighted
ensemble quantiles (uncertainty is **not** collapsed). Ranked by expected cumulative deficit:


| node | E[dur] h [5–95]       | E[deficit] mg/L·h [5–95] | min C med [5–95]   | mean age h |
| ---- | --------------------- | ------------------------ | ------------------ | ---------- |
| 131  | **48.0 [48.0, 48.0]** | **4.22 [4.11, 4.32]**    | 0.034 [0.03, 0.03] | 70.2       |
| 141  | 22.3 [20.0, 24.0]     | 2.20 [2.08, 2.35]        | 0.068 [0.06, 0.07] | 50.8       |
| 139  | 22.0 [22.0, 22.0]     | 1.82 [1.73, 1.91]        | 0.102 [0.10, 0.11] | 53.8       |
| 145  | 12.0 [12.0, 12.0]     | 1.68 [1.62, 1.73]        | 0.031 [0.03, 0.03] | 31.3       |
| 15   | 22.0 [22.0, 22.0]     | 1.38 [1.16, 1.59]        | 0.055 [0.05, 0.06] | 38.9       |
| 143  | 22.0 [22.0, 22.0]     | 1.16 [0.95, 1.37]        | 0.062 [0.06, 0.07] | 39.6       |
| 129  | 24.0 [24.0, 24.0]     | 0.84 [0.80, 0.88]        | 0.147 [0.15, 0.15] | 53.2       |
| 253  | 8.0 [8.0, 8.0]        | 0.75 [0.72, 0.79]        | 0.099 [0.09, 0.10] | 30.3       |


**Water age ↔ risk (all n = 92 junctions):** Spearman ρ = 0.73, descriptive association
(spatial-block bootstrap 95% `[0.455, 0.897]`) for duration and `0.72` for deficit; Pearson `0.86`.
The block bootstrap resamples 10 spatial blocks (k-means on standardised node coordinates,
2000 resamples), so its width is a **conservative descriptive width for the rank correlation, not a
significance interval** — and it is roughly twice the width an iid junction bootstrap would give,
which is the point of using it.
Spearman is the primary statistic (risk metrics are non-linear and many nodes are 0). The ordinary
p-value is deliberately not quoted: the 92 junctions share pipes, flow paths and tank states, so they
are nowhere near 92 independent samples.

**Robustness of the pattern to the inference convention.** Both weightings give effectively the same
operational answer, which is the question the review asked:


| scheme                    | ESS  | network-mean E[D] (h) | network-mean E[A] | top-10 overlap with primary |
| ------------------------- | ---- | --------------------- | ----------------- | --------------------------- |
| formal censored (primary) | 157  | 3.431                 | 0.2181            | —                           |
| informal GLUE             | 4783 | 3.417                 | 0.2207            | Jaccard 0.82                |


The aggregate severity differs by 1.2% between schemes and nine of ten hot-spot nodes are shared
(top-10 Jaccard 0.82; the two lists differ only in node 50 against node 243).
The node identities are, however, **more sensitive to the warm-up than to the weighting**: node 243,
previously ranked worst by deficit, leaves the top ten entirely at the corrected 120 h warm-up.

**Which network average?** An unweighted mean over 92 junctions counts a zero-demand node the same as
the largest consumer. 59 of the 92 junctions carry non-zero demand, totalling **690.7 L/s**, so the
choice matters and all three are reported rather than one being taken silently:


| metric          | unweighted (all 92) | consumer-only (59) | demand-weighted |
| --------------- | ------------------- | ------------------ | --------------- |
| `E[D]` (h)      | 3.4309              | **4.4610**         | **1.2311**      |
| `E[A]` (mg/L·h) | 0.2181              | 0.3053             | 0.0798          |
| min C (mg/L)    | 0.4922              | 0.4412             | 0.6348          |


**The demand weights are pattern-aware, and they have to be.** The weights come from
`wntr.metrics.hydraulic.average_expected_demand`, which applies each junction's pattern and the
global demand multiplier; the network total is **690.7 L/s**. An earlier version of this log read
the raw `base_demand` field instead, and that superseded reading is worth recording because of how
quietly it fails: four Net3 junctions (15, 35, 123, 203) encode a large demand as 1 GPM times a
large pattern, so the base field alone reports them as 0.063 L/s each instead of 16.7 / 108.4 /
75.3 / 284.6 — together **70 % of the network's actual water, counted as four of its smallest
users** — and gives a network total of 192.6 L/s. Every demand-weighted number below, and every
L/s in Step 12, changed when it was corrected.

The two adjustments move in **opposite** directions, and the reason is a result in itself:
consumer-only is *worse* than unweighted (4.46 vs 3.43 h) while demand-weighted is *better*
(1.23 h). That combination can only happen if the risk concentrates at **small** consumers rather
than large ones — network extremities and low-turnover dead-ends, which is where long residence time
is expected. So the worst chlorine conditions coincide with the least demand, and the pattern
correction **sharpens** this: once the four large patterned demands are counted at their true size,
the demand-weighted severity falls from 2.36 h to 1.23 h, i.e. the large consumers are even
better-served relative to the network mean than the base-only reading suggested.

This has two consequences for reporting. The unweighted figure used elsewhere in this log is the
**conservative** choice, not a flattering one — demand-weighting more than halves every severity
number. And a demand-weighted headline would understate exactly the customers the risk map is
supposed to find, so the consumer-only average is the one to quote alongside any demand-weighted
figure.

**Does hourly reporting miss short excursions?** The water-quality solver steps at 300 s but the risk
metrics integrate the **hourly** report, so a dip beginning and ending inside one hour is invisible,
and trapezoidal integration of a binary indicator places any crossing at the midpoint between
reports. The 12 highest-weight members were re-run at three reporting resolutions:


| report step | points in window | network-mean `E[D]` (h) | network-mean `E[A]` | top-10 Jaccard vs hourly |
| ----------- | ---------------- | ----------------------- | ------------------- | ------------------------ |
| 3600 s      | 49               | 3.436                   | 0.2180              | 1.00                     |
| 900 s       | 193              | 3.374                   | 0.2159              | 1.00                     |
| 300 s       | 577              | 3.380                   | 0.2164              | 1.00                     |


Going from hourly to the solver's own 300 s changes `E[D]` by **−1.6%** and `E[A]` by **−0.7%**, and
the top-10 set is identical. Two things are worth noting. The change is **negative**, so hourly
reporting slightly *over*-estimates duration rather than missing excursions — the smearing of
crossings by trapezoidal integration outweighs any dip lost between reports. And most of the change
appears by 900 s and does not continue to 300 s, which is what convergence looks like. **Hourly
reporting is therefore adequate here**, and it is adequate because the chlorine field is smooth on
the hour scale in this network, not because sub-hourly behaviour was assumed away.

Risk associated with water age; top nodes' expected duration with 5-95% bands

Findings:

1. **Duration and depth give different internal rankings within a broadly overlapping hotspot
   cluster.** Node **131** is worst on both counts here (below 0.2 for the full 48 h *and* the
   largest deficit at 4.22 mg/L·h), but below it the two metrics disagree: **145** is a *deep but
   short* risk (min C 0.031, only 12 h) while **129** is *long but shallow* (24 h, min C 0.147, the
   smallest deficit in the table). The cluster {131, 141, 139, 145, 15, 143, 129, 253} overlaps
   broadly across metrics while the internal order changes, which is the point of reporting duration
   and depth separately.
2. **Ensemble uncertainty is much tighter than under the informal score.** With the formal weighting
   the 5–95 % bands collapse to a fraction of their previous width (node 15: 22.0 [22.0, 22.0] h,
   where the informal ensemble gives 14–24 h). That is the same efficiency gain seen in the parameters,
   propagated into the risk metrics: the earlier wide bands were mostly the weighting's inertia, not
   irreducible parameter uncertainty.
3. **Water age gives hydraulic corroboration of the spatial pattern.** Mean water age is strongly
   *associated* with both duration and deficit (Spearman 0.73, n = 92, spatial-block bootstrap
   [0.455, 0.897]); the
   worst nodes also have the longest residence times. This is corroboration, **not** validation:
   water age comes from the *same* hydraulic model (not an independent measurement), correlation is
   not causation, and residual chlorine also depends on wall decay, source dosing and tank mixing.
4. **Water-age magnitudes are window-dependent, not steady state.** Water age remains
   **horizon-dependent within the 168 h simulation**: the field is still filling at the high-age
   nodes throughout the run. The ages in the table above are means over the fixed **120–168 h**
   assessment window and **must not be interpreted as steady-state water ages** — the same nodes give
   different numbers under the last-diurnal-cycle (144–168 h) and final-hour (`t = 168 h`)
   definitions. All three are stored in `step10_risk_metrics.json` under `age_windows`, so the choice
   is always explicit rather than implied. **Step 0 settles what this used to leave open: a longer
   run does not fix it.** The cycle-to-cycle change in the age field is still 12.8 h at a 120–144 h
   warm-up and decays by only ~12% per cycle, so the crossing extrapolates to ≈600 h — far beyond the
   168 h the model permits. Absolute water ages in Net3 are therefore horizon-dependent by
   construction, not merely under-converged here. **The rank association with chlorine risk is the
   relevant result; the absolute magnitudes remain sensitive to the selected horizon** — which is why
   Spearman, and not absolute age, is what this step reports as a finding.
5. **Interpretation, stated carefully — robustness depends on the perturbation AND on the decision
   scale.** The dominant decay term (old) is the *most robustly informed* coefficient under the
   expanded realism analyses (Steps 4/7), and water age is reaction-independent, so the leading
   hotspots are governed by residence time and the old-group decay rather than the weakly-informed
   avg/new coefficients. That buys less than the earlier version of this paragraph claimed:

   - **Robust.** The broad risk pattern is essentially unchanged between the formal and informal
     weighting schemes (network-mean `E[A]` differs by 1.2%, nine of ten hot-spots shared), and it
     stays strongly rank-correlated under every single-sensor bias tested — Spearman ≥ 0.999 across
     all 24 arms of Step 8c, with the top-6 set intact at every negative offset.
   - **Metric-dependent, and this table's metric is the milder one.** Fixing `k_b` 20% off its true
     value changes the leading set **on the time-averaged probability `P_bar`** — top-6 Jaccard 0.50
     at both signs, top-3 0.20, the reference 4th-ranked node dropping to 20th. **On the cumulative
     deficit `E[A]`, which is what the top-10 table above is ranked by, it does much less**: the
     top-6 is *unchanged* at `k_b = −0.4` and loses one node at `−0.6`, in a sixth/seventh-place
     exchange (143 out, 129 in), with top-3 and top-5 identical (Step 8b findings 4–5).

   The defensible summary is: *parameter uncertainty affects the magnitude of the node-level risk
   metrics (see the 5–95% bands); the broad spatial pattern is robust to the weighting scheme and to
   single-sensor bias; the deficit-ranked shortlist reported here is robust to `k_b` ±20% apart from
   one boundary exchange; and the duration-ranked shortlist is not.* We do **not** claim the
   uncertainty "does not propagate" — it does, into the magnitudes — and no statement about the
   stability of "the hot-spot ranking" is meaningful in this project without naming the metric it is
   ranked by.



---



## Step 11 — leave-one-monitor-out (LOO) predictive validation

The draft calibrated and evaluated on the same six monitors, with no out-of-sample check. LOO
cross-validation fills that gap: hold out one monitor, re-weight the cached ensemble on the other
five, then (a) check the three k_w stay stable and (b) *predict* the held-out sensor and measure the
out-of-sample error and predictive-band coverage (`step11_loo.py`; cache reused, 30 noise
realisations). Both weightings are run; the **formal censored likelihood is the primary** and the
informal GLUE score (σ-scaled threshold for 5 monitors) is the comparator. Full-6 reference:
formal old −1.0405 / avg −0.1028 / new −0.0509; informal old −1.0382 / avg −0.1171 / new −0.0516.

```
held-out pred RMSE  = √[ mean_t (pred_mean_m(t) − obs_m(t))² ]     compared with the noise floor σ = 0.1
90% predictive band = pred_mean ± 1.645·√(Var_ensemble + σ²)       (parameter + observation noise)
coverage            = fraction of held-out hours inside the 90% band
```


**Primary — formal censored likelihood:**

| held-out (zone) | k_old      | k_avg  | k_new  | pred RMSE @m | 90% coverage |
| --------------- | ---------- | ------ | ------ | ------------ | ------------ |
| 107 (new)       | -1.038     | -0.104 | -0.050 | 0.098        | 0.92         |
| 113 (new)       | -1.031     | -0.098 | -0.053 | 0.103        | 0.90         |
| **15 (old)**    | **-1.086** | -0.102 | -0.050 | 0.094        | 0.94         |
| 145 (old)       | -1.060     | -0.102 | -0.051 | 0.092        | 0.94         |
| 209 (average)   | -1.040     | -0.102 | -0.050 | 0.099        | 0.91         |
| 231 (average)   | -1.045     | -0.106 | -0.051 | 0.099        | 0.92         |

**Comparator — informal GLUE:**

| held-out (zone) | k_old      | k_avg  | k_new  | pred RMSE @m | 90% coverage |
| --------------- | ---------- | ------ | ------ | ------------ | ------------ |
| 107 (new)       | -1.040     | -0.116 | -0.051 | 0.097        | 0.92         |
| 113 (new)       | -1.036     | -0.116 | -0.052 | 0.102        | 0.90         |
| **15 (old)**    | **-0.985** | -0.115 | -0.053 | 0.095        | 0.94         |
| 145 (old)       | -1.025     | -0.117 | -0.051 | 0.091        | 0.94         |
| 209 (average)   | -1.037     | -0.119 | -0.053 | 0.099        | 0.92         |
| 231 (average)   | -1.035     | -0.119 | -0.052 | 0.098        | 0.92         |


(medians over 30 noise realisations; noise floor σ = 0.1 mg/L.)

LOO out-of-sample prediction error ≈ noise floor; held-out node-15 band covers the data

Findings:

1. **Out-of-sample prediction is at the noise floor.** Under the **formal primary** every held-out
   monitor is predicted with RMSE **0.092–0.103** mg/L ≈ σ = 0.1 (informal comparator 0.091–0.102) —
   i.e. the five-monitor calibration predicts the *unseen* sensor about as accurately as the
   measurement noise allows. No sign of over-fitting; the model generalises spatially. The two
   weightings agree here, which is worth noting given how often they do not: prediction is the one
   thing the informal score does about as well, because it is not being asked to extract parameter
   information.
2. **Uncertainty is close to nominal, and "conservative" overstates it.** The 90 % predictive band
   (parameter + observation-noise variance) covers **0.90–0.94** of held-out hours under both rules.
   The average sits slightly above nominal, but the minimum is **exactly** 0.90 (monitor 113), so
   the band is well-calibrated rather than systematically conservative — the genuinely conservative
   case is Step 11c's unmonitored junctions, which over-cover at 1.00.
3. **Parameter stability confirms where the information lives — but the two weightings disagree on
   the sign, and the primary is the one to quote.** Under both rules `k_old` moves perceptibly **only
   when an old-zone monitor is removed** and barely at all when a new/avg monitor is dropped, which
   is exactly what Fisher and the profile predict: old's information sits at monitors 15/145. The
   *direction* differs, though. Under the **formal primary** dropping node 15 moves `k_old` from
   −1.0405 to **−1.086** and node 145 to −1.060 — i.e. *stronger* decay, away from the truth −1.0.
   Under the informal comparator the same drop gives **−0.985** and −1.025 — *weaker* decay, toward
   the truth. The magnitude is similar (~0.05) but the sign is opposite, so a sentence like "removing
   node 15 pulls old toward −0.99" is only true of the comparator.

   **`k_avg` and `k_new` are also stable under both rules — but for different reasons, and the
   formal one must not be described as prior domination.** Under the informal comparator they are
   genuinely prior-dominated: they retain 92% and 88% of the prior width (Step 3), so there is little
   for a dropped monitor to change. Under the **formal primary they are informed** — 30.7% and 28.3%
   of the prior width retained (Step 1), and their ensemble SD sits within a few percent of the
   Case-A CRLB (Step 7). Their stability here has a different cause: **leave-one-monitor-out always
   leaves the partner monitor of the same zone in the calibration set**, so no zone ever becomes
   unobserved and each coefficient keeps its own source of information. Calling this prior
   domination would contradict the central formal result of this project. What the loss of a zone's
   information actually does is only visible in the harder leave-one-**zone**-out experiment below,
   where both of a zone's monitors are removed and the coefficient does return to its prior.
4. **But this is the easy case, and on its own it does not support the claim it was asked for.**
   Every held-out monitor has a partner in the same zone, so the zone stays observed. Two harder
   tests were added.



### Step 11b — leave-one-ZONE-out: what happens when a zone becomes unobserved

Both monitors of one zone are dropped, so the calibration has four monitors and one zone with no
sensor at all. This is the test that speaks to the actual claim.

Both weightings are run: informal for continuity with Step 11a, formal so that the heterogeneous-truth
test below differs from this one *only* in the truth.


| scheme   | zone dropped | monitors | `k_w,old`  | `k_w,avg`  | `k_w,new`  | own-coef error | own SD retained | pred RMSE | 90% cov |
| -------- | ------------ | -------- | ---------- | ---------- | ---------- | -------------- | --------------- | --------- | ------- |
| informal | old          | 15, 145  | **−0.850** | −0.115     | −0.052     | **+0.150**     | **100%**        | 0.100     | 0.96    |
| informal | average      | 209, 231 | −1.045     | **−0.120** | −0.054     | **−0.020**     | **100%**        | 0.099     | 0.93    |
| informal | new          | 107, 113 | −1.035     | −0.116     | **−0.052** | −0.002         | 95%             | 0.100     | 0.92    |
| formal   | old          | 15, 145  | **−0.851** | −0.102     | −0.051     | **+0.149**     | **100%**        | 0.100     | 0.96    |
| formal   | average      | 209, 231 | −1.038     | **−0.120** | −0.050     | **−0.020**     | **100%**        | 0.099     | 0.92    |
| formal   | new          | 107, 113 | −1.049     | −0.101     | **−0.049** | +0.001         | **57%**         | 0.101     | 0.90    |


**Finding, and it changes what the LOO result is worth.** Remove a zone's two monitors and that
zone's coefficient **reverts to its prior midpoint** — `old` to −0.850 (prior mid −0.850) with 100% of
the prior width retained, `average` to −0.120 (prior mid −0.120), also 100%. Under both weightings.
So for those two zones each coefficient is informed by its own two monitors and by nothing else.

`new` is the exception and it is instructive: under the formal likelihood it keeps **57%** of the
prior width even with both its monitors gone, so some information about it does arrive from
elsewhere. That fits Step 7's Jacobian, where `new` has the largest per-unit sensitivity at *every*
monitor (2.12 against 0.19 mg/L per m/day for `old`) because the new-zone pipes sit upstream of the whole
network. Spatial borrowing exists, but only for the zone the water passes through first.

The decisive observation is the pair of right-hand columns: **prediction quality is unchanged**.
Held-out RMSE is still ≈0.10 mg/L and coverage still 0.90–0.96 (formal; 0.92–0.96 informal), even for the zone whose coefficient
carries a +0.150 error and has learned nothing. So **predictive success is not evidence of parameter
identifiability**, and the LOO result of Step 11a cannot be read as validating the coefficients.
Prediction survives because the prior is centred near the truth — the same artefact the review
identified in §3.1 — so the prior-mean prediction is already good. §3.1's problem and §4-#6's
validation are the same problem seen twice.

### Step 11c — junctions that never enter any calibration

20 unmonitored junctions (seed 7), predicted against the **noise-free truth**. The band construction
differs from Step 11a on purpose: predicting a truth value needs the parameter spread only, with no
measurement-noise term, because there is no sensor there to add noise.

| weighting       | 90% coverage (normal approx.) | 90% coverage (weighted quantile) | mean |pred − truth| / mean|truth| |
| --------------- | ----------------------------- | -------------------------------- | --------------------------------- |
| informal GLUE   | 1.00                          | 1.00                             | 0.006                             |
| formal censored | 1.00                          | 1.00                             | 0.008                             |

Nominal coverage is 0.90, so both schemes **over-cover**: the bands are conservative rather than
calibrated. Mean prediction error at unmonitored junctions is 0.6–0.8% of the local concentration.
The weighted-quantile band is reported alongside the normal approximation because the ensemble is
skewed and censored at zero, where a symmetric band is not conservative in the direction that
matters; here the two agree, so the approximation is adequate at these junctions.

### Step 11d — the hardest case: leave-one-zone-out on a heterogeneous truth

Tests a–c all use a truth generated by the same three-zone model that is then fitted, so no structural
discrepancy can appear. Here the truth carries ±20% per-pipe heterogeneity (the Step 5c design, 8
independent fields) **and** a zone is unobserved, so extrapolation and structural error are tested
together — the combination the risk map actually relies on. Formal weighting, so the only difference
from the formal rows of Step 11b is the truth. The coefficient error is measured against each field's
own arithmetic zone mean, so it is not inflated by the heterogeneity itself.


| zone dropped | own-coef error | own SD retained | pred RMSE | 90% cov | unmonitored rel. error |
| ------------ | -------------- | --------------- | --------- | ------- | ---------------------- |
| old          | +0.139         | 100%            | 0.098     | 0.96    | 0.008                  |
| average      | −0.020         | 100%            | 0.097     | 0.94    | 0.012                  |
| new          | +0.012         | 55%             | 0.102     | 0.88    | 0.011                  |


**The structural discrepancy adds almost nothing on top of losing the sensors.** Compare row by row
with the formal rows of Step 11b: coefficient errors +0.139 / −0.020 / +0.012 against +0.149 / −0.020
/ +0.001, SD retained 100 / 100 / 55% against 100 / 100 / 57%, prediction RMSE and coverage within
0.005 and 0.02. Prediction at unmonitored junctions stays at 0.8–1.2% relative error. So in this
network, at this noise level, **the cost of an unobserved zone dominates the cost of ±20% within-zone
structural error by a wide margin** — consistent with Step 5c, where the same heterogeneity produced
no detectable bias across 25 fields.

One number does move in the right direction to be worth noting: coverage for the dropped `new` zone
falls to 0.88, just under the nominal 0.90, the only sub-nominal figure in any of these tests. It is
the case with both a structural discrepancy and the partially-informed coefficient, so it is where
the bands would be expected to lose calibration first.

**What the four tests together do and do not establish.** They establish internal consistency and
they locate its limits precisely: the calibrated ensemble reproduces unseen monitors at the noise
floor and unseen junctions to about 1%, with conservative bands, and it keeps doing so under ±20%
structural error. They do **not** establish external validity — the truth is always model-generated —
and Step 11b shows they do not even establish that the coefficients are informed, since prediction
succeeds while a coefficient sits at its prior. The honest summary is that **spatial prediction and
parameter identification are separately valid claims here, and only the first of them is supported.**

---



## Step 12 — operational temperature / ageing scenario projection (WSP application)

**Role.** After calibration (Steps 1–9), baseline risk metrics (Step 10) and LOO validation
(Step 11), this step answers the operational question the supervisor posed: *given the calibrated
ensemble, what happens to network-wide low-chlorine risk under warm-season and
heatwave conditions with an ageing-reactivity stress, and does raising the source dose
restore the baseline position?* The supervisor phrased it in terms of the GLUE ensemble, because
that is what the draft had; the ensemble propagated here is weighted by the **formal censored
likelihood**, so the word GLUE is not used for it below.

This is **not** a re-run of the enclosed homogeneous-`k_w` notebook. The scenario *method*
(Arrhenius scaling → ensemble propagation → likelihood×consequence → risk register → dosing
evaluation) is transplanted onto **this** project's six-monitor, three-zone **formally weighted**
ensemble (`k_b = −0.5` fixed; formal censored weights, ESS 157). Of the 8192 draws, 2196 are carried
forward: those whose weight exceeds `1e-6` of the maximum. That is a **numerical truncation to keep
the scenario runs affordable, not a behavioural threshold** — the formal likelihood has no
acceptance cut-off. The cost is now recorded rather than asserted: the 5996 discarded draws carry a
combined weight of **3.79 × 10⁻⁶**, so the retained ensemble holds 0.9999962 of the total, orders of
magnitude below the resolution of any number reported here (`discarded_weight_mass` in the
artifact). The supervisor's figures are not comparable with these and must never be
quoted as results of this study.

### 12.1 Two probability definitions — keep them distinct

Step 10 reported a *time-averaged* probability; Step 12's headline is a *window-breach*
probability. They answer different questions and both are reported.

**Formulas**:

```
P_min(s,n) = Σᵢ wᵢ · 1[ min over t ∈ [120, 168] h of C(s,i,n,t) < C_crit ]
P_bar(s,n) = E[D(s,n)] / 48                     (the Step-10 quantity)
E[D(s,n)]  = Σᵢ wᵢ · ∫ 1[C < C_crit] dt          (h, trapezoid over 48 intervals)
E[A(s,n)]  = Σᵢ wᵢ · ∫ max(0, C_crit − C) dt     (mg/L·h)
```

**Properties**:

- `P_min` = probability the node dips below `0.2 mg/L` **at least once** in the window;
`P_bar` = probability a randomly chosen hour is below it. `P_min ≥ P_bar` always.
- The window is **48 h** (`t = 120…168`, post warm-up), so `P_min` is a 48-hour window
minimum, **not** a 24-hour "daily minimum". For otherwise identical trajectories a
48-hour window-breach probability **cannot be lower** than the corresponding probability
on a nested 24-hour window, and the two are **equal** whenever the diurnal pattern
repeats. Window length is therefore only a secondary reason not to compare with the
supervisor's figures; the primary reasons are the different parameterisation, monitoring
array, ensemble and simulation set-up.
- "`P_min > 0.5` nodes" counts over all 92 junctions; "indeterminate" is `0.05 < P_min < 0.95`;
"demand at risk" sums base demand over junctions with `P_min > 0.5` (zero-demand nodes
contribute 0); "network mean" of `E[D]` / `E[A]` is the **unweighted arithmetic mean over
all 92 junctions**, in `h` and `mg/L·h` respectively.



### 12.2 Scenario construction

Three sources of scenario uncertainty are propagated jointly with **common random numbers** —
one `(E_a,bulk, E_a,wall, δT)` draw per behavioural member, reused across every scenario and
every dose — so scenario differences are physical, not Monte-Carlo noise.

**Formula**:

```
T(s,i)       = T_mean(s) + δTᵢ,          δTᵢ ~ N(0, 1²) °C
f(T, Ea)     = exp[ −(Ea/R) · (1/T − 1/T_ref) ],   T in K
k_b,i(s)     = (−0.5) · f(T(s,i), E_a,b,i)
k_w,g,i(s)   = k_w,g,i(T_ref) · f(T(s,i), E_a,w,i) · α_g(s),   g ∈ {old, average, new}
```

with `E_a,bulk ~ N(45, 8²)` kJ/mol and `E_a,wall ~ N(35, 10²)` kJ/mol, both **truncated below
at 5 kJ/mol** (drawn from a truncated normal rather than clipped, so no draw can be
non-physical and no probability mass piles up on the bound). The script asserts it: realised
draws span `17.59–72.05` (bulk) and `5.10–70.24` kJ/mol (wall), all above the 5 kJ/mol bound. The water-temperature
offset is now truncated too, to `±4.0 °C`, which is exactly the range that keeps every scenario mean
inside the stated `8–24 °C` validity window; realised temperatures span `8.48–23.50 °C` and no draw
needed resampling. All of these are recorded in `step12_scenarios.json` rather than only printed,
so they can be checked.


| key | label                    | mean water T | ageing multipliers `α`         |
| --- | ------------------------ | ------------ | ------------------------------ |
| A   | Baseline 12 °C           | 12 °C        | all 1.00                       |
| B   | Warm season 16 °C        | 16 °C        | all 1.00                       |
| C   | Heatwave 20 °C           | 20 °C        | all 1.00                       |
| D   | Heatwave + ageing stress | 20 °C        | new 1.00 / avg 1.35 / old 1.85 |


**Assumptions to state in the thesis**:

- Calibrated coefficients are **effective parameters at an illustrative reference temperature**
`T_ref = 12 °C`; this is inherited from the supervisor's material, **not** measured from Net3.
- Ageing multipliers are **escalation-only stress-test inputs** (`α ≥ 1`), not asset records.
Because the baseline model already distinguishes new/average/old `k_w`, `α_new < 1` is not
used — it would weaken already-weak new pipes and double-count the zone structure. Scenario
D is therefore named an *ageing-reactivity stress*, not "the effect of ageing".
- `C_crit = 0.2 mg/L` and the demand-based consequence terciles (`1.67`, `4.38 L/s`) are
illustrative operational inputs, not legal limits. Demand is a **consequence proxy**, not
population or exposure.

**Numerical guard (important).** Zonal `k_w` are clipped to `CLIP_LO = −8.0 m/day` purely as a
solver guard. An earlier run used `−3.0`, which was **active for 46.8 % of members in scenario
D** (prior lower bound `−1.5` × `f_w ≈ 1.5` × `α_old = 1.85` reaches `−5.3`) and silently
biased the reported mean. The guard is now inert: the script asserts **0 clipped draws in every
scenario**, and verifies

```
mean k_w,old(D) / mean k_w,old(C) = 1.8500 = α_old      (exact to 1e-6)
```

Any future change to the priors or to `α` must re-check this assertion.

### 12.3 Scenario results

All means below are weighted by the **formal censored likelihood** (the primary scheme from Step 1);
the informal GLUE score is not used for these headline numbers.


| Scenario                    | mean k_b | mean k_w,old | `P_min>0.5` nodes | demand at risk | % demand | high/very-high | indeterminate | net-mean `P_bar` | net-mean `E[D]` (h) | net-mean `E[A]` (mg/L·h) |
| --------------------------- | -------- | ------------ | ----------------- | -------------- | -------- | -------------- | ------------- | ---------------- | ------------------- | ------------------------ |
| A. Baseline 12 °C           | −0.504   | −0.994       | 21                | 55.4 L/s       | 8.0 %    | 10             | 0             | 0.0708           | 3.398               | 0.2216                   |
| B. Warm season 16 °C        | −0.655   | −1.222       | 28                | 64.7 L/s       | 9.4 %    | 12             | 8             | 0.0884           | 4.245               | 0.3645                   |
| C. Heatwave 20 °C           | −0.848   | −1.501       | 29                | 67.8 L/s       | 9.8 %    | 13             | 4             | 0.1087           | 5.218               | 0.5218                   |
| D. Heatwave + ageing stress | −0.848   | −2.777       | 31                | 69.4 L/s       | 10.0 %   | 13             | 4             | 0.1232           | 5.912               | 0.5790                   |


Continuity check against Steps 1–11: with the temperature held at `T_ref` exactly (cached
`C_all`, no `δT`), the baseline gives 21 nodes / 55.4 L/s and `E[A] = 0.2181`. Sampling the
stated temperature uncertainty therefore did **not materially change** the baseline
classification — the `P_min>0.5` count and demand at risk are unchanged and network-mean
`E[A]` moves by under 2 % (`0.2181 → 0.2216`).

Four-panel window-breach probability mapsAgeing increment ΔP_min = P(D) − P(C) at 20 °CScenario escalation, dosing evaluation and ageing sensitivity

**Escalation.** Seven consumer junctions change risk band between A and D, carrying 20.23 L/s
(2.9 % of network demand); **six of the seven are unmonitored**. The largest ageing-only
increment at the same temperature (D − C) is at node 249 (`ΔP_min = +0.410`), which carries **no
demand**; the largest at consumer junctions are 217 (`+0.346`), 239 (`+0.155`) and 215
(`+0.091`). Demand share is small here because the escalating nodes are network extremities —
the same "risk sits at small consumers" pattern Step 10 quantifies, so a demand-share headline
understates the operational relevance of these seven nodes.

**Which band is the escalation judged on?** On the **governing** band (§12.5.1), so the escalation
column follows the same rule as the control measures. Judged on the breach band alone it returns the
*identical* set — 7 junctions, 20.23 L/s, 6 unmonitored, and no node is found by one test and not the
other. The finding above is therefore unchanged by the governance rule.

**But a band-change test is blind to a large amount of deterioration, and the governing rule does
not fix that.** Sixteen consumer junctions (52.3 L/s) sit at `P_min ≈ 1` in *both* A and D, so their
likelihood band cannot move by construction. Their below-threshold durations move a great deal:


| node | demand (L/s) | `E[D]` A → D (h) | Δ | severity band A → D | flagged? |
| ---- | ------------ | ---------------- | ------ | ------------------- | -------- |
| 243  | 0.29         | 13.09 → **48.00**| +34.91 | medium → medium     | no       |
| 15   | **16.67**    | 21.95 → 33.51    | +11.56 | very high → very high | no     |
| 143  | 0.42         | 21.72 → 32.09    | +10.37 | medium → medium     | no       |
| 151  | **9.75**     | 7.48 → 12.00     | +4.52  | high → **very high**| **no**   |
| 131  | 2.88         | 46.55 → 48.00    | +1.45  | very high → very high | no     |


Node 243 spends the **entire 48 h window** below the threshold under scenario D and is not flagged;
node 15, the network's largest single consumer, gains 11.6 h and is not flagged; node 151 (9.75 L/s)
actually **changes its severity band** and is still flagged by neither test.

**Two different mechanisms cause this, and they were conflated in an earlier version of this log.**

1. **Saturation.** A junction already at `very high` cannot rise, so `max(breach, severity)` has no
   headroom. The governance rule decides *which axis drives an action*; it does not make a band
   comparison sensitive to change at the ceiling. Those are separate properties and only the first
   was argued for.
2. **The ceiling of the multiplicative matrix.** `score = axis score (1–5) × consequence score
   (0–3)`, so a **minor** consumer (consequence 1) tops out at score 5 → `medium`, whatever happens
   to it. **20 junctions (17.46 L/s) are capped below `high` by construction.** The register cannot
   express "this small consumer has no chlorine for the whole window" as anything worse than medium.
   That is a property of the `likelihood × consequence` form itself, not of this implementation.

So the band-change columns are **necessary but not sufficient**. For any junction whose band cannot
move, escalation has to be read on the continuous metrics (`ΔE[D]`, `ΔE[A]`), which is why they are
carried in the artifact.

**Ageing-stress sensitivity — is scenario D an artefact of** `α_old = 1.85`**?**


| ageing set | `α_avg` | `α_old` | `P_min>0.5` nodes | demand at risk | high/very-high | indeterminate | net-mean `E[D]` (h) | net-mean `E[A]` (mg/L·h) |
| ---------- | ------- | ------- | ----------------- | -------------- | -------------- | ------------- | ------------------- | ------------------------ |
| mild       | 1.15    | 1.40    | 30                | 69.4 L/s       | 13             | 4             | 5.515               | 0.550                    |
| central    | 1.35    | 1.85    | 31                | 69.4 L/s       | 13             | 4             | 5.912               | 0.579                    |
| severe     | 1.50    | 2.20    | 31                | 69.4 L/s       | 14             | 4             | 6.178               | 0.599                    |


The **headline binary metrics are nearly insensitive** to the tested multiplier range: the count
varies only from 30 to 31, demand at risk stays at 69.4 L/s and the indeterminate count is
constant. The **continuous severity metric increases monotonically** (`E[A]` 0.550 → 0.579 →
0.599 mg/L·h), and node-level probabilities do move with `α`. So the hot-spot conclusion is
robust to the multiplier choice, but the magnitude of the ageing penalty is not — it should be
reported as a stress-test range, never as "the effect of ageing".

### 12.4 Corrective dosing (control-measure evaluation, not a recommendation)

The dose scales **both** the reservoir source quality and the tank initial quality, so the whole
source regime is raised consistently. (An earlier run left tanks fixed at `0.5 mg/L`, which
confounded the result with an un-dosed boundary condition.) `inlet = 1.00` reuses scenario C.


| inlet (mg/L) | `P_min>0.5` nodes | demand at risk | % demand | net-mean `P_min` | net-mean `E[D]` (h) | net-mean `E[A]` (mg/L·h) | median-over-nodes of mean window min (mg/L) |
| ------------ | ----------------- | -------------- | -------- | ---------------- | ------------------- | ------------------------ | ------------------------------------------- |
| 1.00         | 29                | 67.8 L/s       | 9.8 %    | 0.3246           | 5.218               | 0.5218                   | 0.463                                       |
| 1.15         | 28                | 64.7 L/s       | 9.4 %    | 0.3112           | 4.628               | 0.4527                   | 0.532                                       |
| 1.30         | 28                | 64.7 L/s       | 9.4 %    | 0.2915           | 4.288               | 0.3957                   | 0.601                                       |


(The last column takes, per junction, the likelihood-weighted mean of the per-member window minimum,
then the median across the 92 junctions — it is **not** a pooled member-node median. It scales
almost exactly with the dose, `0.463 × 1.15 = 0.532`, `0.463 × 1.30 = 0.601`.)

**Linearity known-answer test.** Measured `max |C(1.3) − 1.3·C(1.0)| = 1.4e-05 mg/L`. Stated
precisely: *under fixed hydraulics and demands, first-order bulk and wall kinetics, and
proportional scaling of all source and initial chlorine concentrations, increasing the dose by a
factor* `r` *is mathematically equivalent to evaluating the unscaled concentration field against a
threshold* `C_crit / r`*.* This equivalence is a property of that idealisation — it does **not**
imply that dosing is linear in a real network with dose-dependent chemistry, DBP formation or
booster constraints.

**Paired warm-up test — the roles are now reversed.** In the draft version of this script the
168 h / 120 h horizon was the experimental arm being tested against a 72 h / 24 h baseline. Step 0
settled that question, so 168 h / 120 h is now the baseline and the short horizon is the
comparator: the rows below measure **what the draft's warm-up cost**, not whether the current one
is adequate. Identical 275-member subset (stride 8), identical weights and identical `(E_a, δT)`
draws; only the simulation horizon differs.


| inlet (mg/L) | nodes short / long | demand at risk short / long | `E[D]` short / long (h) | `E[A]` short / long (mg/L·h) | rel. change |
| ------------ | ------------------ | --------------------------- | ----------------------- | ---------------------------- | ----------- |
| 1.00         | 31 / 31            | 69.6 / 69.6 L/s             | 5.370 / 5.389           | 0.5052 / 0.5484              | +8.6 %      |
| 1.15         | 29 / 29            | 64.9 / 67.8 L/s             | 4.655 / 4.788           | 0.4306 / 0.4789              | +11.2 %     |
| 1.30         | 29 / 28            | 64.9 / 64.7 L/s             | 4.040 / 4.369           | 0.3736 / 0.4220              | +13.0 %     |


Read honestly, with one caveat that has grown: the subset is now less representative than before
(long-horizon subset `E[A] = 0.5484` against the full ensemble's `0.5218`, a 5 % gap, because the
formal weights are concentrated and a stride-8 thinning samples them unevenly). So the +8.6 % to
+13.0 % figures carry a few percent of subset error and should be read as "of order 10 %", which
is also what Step 0's cycle-to-cycle analysis independently gives. The **conclusion is unchanged**
under both horizons: node counts and demand at risk are essentially the same, every metric still
improves monotonically with dose, and demand at risk stays well above the 55.4 L/s baseline.

**Read the binary and continuous columns together.** The identical 64.7 L/s at 1.15 and 1.30 is
a **threshold artefact**, not evidence that the extra dose does nothing: the continuous metrics
improve monotonically (`E[D]` 5.22 → 4.63 → 4.29 h; `E[A]` 0.522 → 0.453 → 0.396). Even so,
+30 % dose leaves 64.7 L/s at risk against 55.4 L/s at baseline — **source dosing alone does not
restore the pre-heatwave demand-at-risk position**. Turnover measures (flushing, storage
management, rezoning) act on residence time and should be evaluated alongside dosing. DBP
formation, taste and acceptability are not modelled.

### 12.5 Risk register, product statement and review triggers

`baseline_cache/step12_risk_register.csv` (92 rows): `P_min` under current / heatwave /
heat+ageing, `P_bar`, `E[D]`, `E[A]` at baseline, demand, likelihood, consequence, risk score
and bands, escalation flag, monitor flag, sampling priority and control-measure text.

**Classification rules (put these in Methods or an appendix — the register is not reproducible
without them).**


| likelihood band | on `P_min`    | score |
| --------------- | ------------- | ----- |
| rare            | `< 0.05`      | 1     |
| unlikely        | `0.05 – 0.20` | 2     |
| possible        | `0.20 – 0.50` | 3     |
| likely          | `0.50 – 0.80` | 4     |
| almost certain  | `≥ 0.80`      | 5     |



| consequence band | on average expected demand `d` (L/s) | score |
| ---------------- | ----------------------------------- | ----- |
| non-consumer     | `d = 0`                             | 0     |
| minor            | `0 < d ≤ 1.67`                      | 1     |
| moderate         | `1.67 < d ≤ 4.38`                   | 2     |
| major            | `d > 4.38`                          | 3     |


Terciles are the `1/3` and `2/3` quantiles of demand over the **59 junctions with non-zero
demand** (WNTR parses the `.inp` and stores every demand internally in **m³/s**, so `×1000` gives
L/s; the frozen `Net3.inp` itself declares `Units GPM` — the conversion is from WNTR's internal SI
value, not from the file's units). `d` is the **pattern-aware average expected demand**
(`wntr.metrics.hydraulic.average_expected_demand`), not the `base_demand` field: four Net3
junctions encode a large demand as 1 GPM times a large pattern, and reading the base value alone
would put them in the *minor* band while they carry 70 % of the network's water. The breach risk
score is `likelihood score × consequence score` (range `0–15`), mapped as `0 → not applicable`,
`1–3 → low`, `4–6 → medium`, `7–9 → high`, `≥10 → very high`. §12.5.1 adds a second product on the
same consequence axis and a declared rule for which one drives the control measure.

**What this score is, and what it is not.** Its likelihood axis is `P_min` — the probability of *any*
breach in the 48 h window — so this product scores **whether a node breaches at all, weighted by how
much water it serves**. On its own it conflates a node marginally below the threshold for the whole
window with one far below it for two hours: a node that dips just under 0.2 mg/L once and a node
that sits at 0.03 mg/L for 48 h receive the same likelihood band. It is a **breach-probability ×
consequence** product and is named `risk_score_breach` in the register.

### 12.5.1 A second, parallel axis: severity

The register therefore carries a **second product on the same consequence axis**, scoring *how long*
rather than *whether*:


| severity band | on `E[D]` (h below 0.2 mg/L in the 48 h window) | score |
| ------------- | ---------------------------------------------- | ----- |
| negligible    | `E[D] < 1`                                     | 1     |
| brief         | `1 ≤ E[D] < 6`                                 | 2     |
| sustained     | `6 ≤ E[D] < 12`                                | 3     |
| prolonged     | `12 ≤ E[D] < 24`                               | 4     |
| persistent    | `E[D] ≥ 24` (half the window or more)          | 5     |


`risk_score_severity = severity score × consequence score`, the same `0–15` range and the same band
mapping, so the two products are directly comparable.

Three choices need stating. **The edges are absolute pre-declared hours, not quantiles of this
network's own results** — a scale taken from the data would shift with the heatwave and could
therefore never show escalation. **`E[D]` is the axis rather than `E[A]`** because hours-below has
an operational meaning that a `mg/L·h` integral does not. And because `E[D]` alone cannot separate a
long shallow excursion from a short deep one, the register also carries
`E_depth_while_below_mgL = E[A]/E[D]`, the mean depth below the threshold while below it.

**The two axes are not independent — but the constraint is weaker than it looks, and an earlier
version of this section overstated it.** Since `D_i ≤ T_window · 1[member i breaches]`, taking
weighted expectations gives

```
E[D] ≤ T_window · P_min        equivalently   P_bar ≤ P_min
```

Verified on the register: **0 of 92 junctions violate it**, and the tightest node reaches
`E[D] / (T_window · P_min) = 0.97`.

This section previously read that as "the severity score can therefore never exceed the likelihood
score, so the disagreement is necessarily one-sided". **That is wrong.** The inequality constrains
the *continuous* quantities; the two scores are banded on **different and unaligned scales**, so the
bound restricts an inversion without forbidding it. A counterexample inside the bound:

```
P_min = 0.03,  E[D] = 1.2 h        1.2 ≤ 48 × 0.03 = 1.44  ✓ satisfies the bound
                                   likelihood band 1 (P < 0.05)
                                   severity band   2 (E[D] ≥ 1 h)      -> severity is HIGHER
```

What the bound does buy is a **limit on the size** of an inversion, derived in the code from the
edges themselves rather than asserted, so it stays true if either scale is ever changed:


| likelihood band | sup `P_min` | ⇒ sup `E[D]` | max reachable severity band | max inversion |
| --------------- | ----------- | ------------ | --------------------------- | ------------- |
| 1               | 0.05        | 2.4 h        | 2                           | **+1**        |
| 2               | 0.20        | 9.6 h        | 3                           | **+1**        |
| 3               | 0.50        | 24.0 h       | 4                           | **+1**        |
| 4               | 0.80        | 38.4 h       | 5                           | **+1**        |
| 5               | 1.00        | 48.0 h       | 5                           | 0             |


So severity can exceed likelihood by **at most one band**, and cannot exceed it at all once
`P_min ≥ 0.80`. Everything below about the observed direction of disagreement is therefore an
**empirical** statement about this field, not a consequence of the inequality.

**And it is empirically close.** Measuring, for every consumer junction, how much more `E[D]` it
would need before its severity score overtook its likelihood score: the nearest is **0.939 h** away
— under an hour of additional below-threshold time in a 48 h window. Seventeen of the 59 consumers
sit at a `P_min` high enough that they can never invert, but the rest are not protected by anything
structural. "No inversions here" is a fact about this field with roughly an hour of margin, and it
should never have been written as a theorem.

**Why the duration edges are not aligned to the probability edges.** Setting them to
`2.4 / 9.6 / 24 / 38.4 h` — the probability edges times the window — would make the ordering
guaranteed. It is rejected for two reasons. It destroys the operational meaning of the severity
bands, which was the entire argument for choosing `E[D]` over `E[A]`; and it would make the ordering
a property of the *scales* rather than a finding about the network. A second axis that cannot
disagree in one direction by construction is a weaker instrument, and "severity never exceeds
breach" would then be a tautology rather than a result.

**Result (scenario A, the 59 consumer junctions).** Severity bands: negligible 42, brief 1,
sustained 6, prolonged 9, persistent 1. Against the breach product: **49 in the same risk band, 10
in a lower one, 0 in a higher one.**

**All four scenarios are banded, and this is where the absolute edges earn their keep.** The
severity scale was declared in absolute hours precisely so it could be carried across scenarios; the
portability is now exercised rather than merely claimed (`E[D]` is already computed for every
scenario, so it costs no simulation):


| scenario | negligible | brief | sustained | prolonged | **persistent** | high/v-high (breach) | high/v-high (severity) |
| -------- | ---------- | ----- | --------- | --------- | -------------- | -------------------- | ---------------------- |
| A. 12 °C            | 42 | 1 | 6 | 9 | **1** | 10 | 6 |
| B. 16 °C            | 38 | 5 | 6 | 6 | **4** | 12 | 6 |
| C. 20 °C            | 36 | 3 | 8 | 6 | **6** | 13 | 7 |
| D. 20 °C + ageing   | 35 | 4 | 6 | 6 | **8** | 13 | 7 |


**The two axes escalate differently, and the severity axis is the one that keeps moving.** The
count of *persistent* consumers — below 0.2 mg/L for at least half the window — runs
**1 → 4 → 6 → 8**, monotonically and without flattening. The breach count runs 10 → 12 → 13 → **13**
and **saturates between C and D**: at 20 °C most of the exposed set is already at `P_min > 0.8`, so
adding the ageing stress cannot move it. Demand at risk on the severity axis (sustained or worse)
goes 50.6 → 50.6 → 60.0 → 60.0 L/s.

This is the strongest argument in this log for carrying the second axis at all: **the ageing stress
in scenario D is nearly invisible to the breach product and clearly visible to the severity
product.** A register built on breach probability alone would report that ageing adds almost
nothing on top of a heatwave, which is not what the field does.

What the ten share is the *probability* and the *duration*: `P_min = 1.000` (node 125: 0.998) with
`E[D]` of 2.0–17.9 h out of 48. What they do **not** share is depth, and an earlier version of this
section wrongly flattened them into "marginally below". Sorted by depth:


| node | `P_min` | `E[D]` (h) | depth while below (mg/L) | mean C while below | `E[A]` (mg/L·h) | demand (L/s) | breach band | severity band |
| ---- | ------- | ---------- | ------------------------ | ------------------ | --------------- | ------------ | ----------- | ------------- |
| 145  | 1.000   | 12.02      | **0.1400**               | **0.0600**         | **1.6826**      | 1.86         | very high   | high          |
| 251  | 1.000   | 6.00       | 0.0954                   | 0.1046             | 0.5727          | 1.63         | medium      | low           |
| 253  | 1.000   | 8.00       | 0.0947                   | 0.1053             | 0.7577          | 3.68         | very high   | medium        |
| 247  | 1.000   | 2.00       | 0.0899                   | 0.1101             | 0.1799          | 4.75         | very high   | medium        |
| 255  | 1.000   | 8.00       | 0.0880                   | 0.1120             | 0.7039          | 2.73         | very high   | medium        |
| 147  | 1.000   | 7.45       | 0.0793                   | 0.1207             | 0.5906          | 0.58         | medium      | low           |
| 149  | 1.000   | 6.04       | 0.0487                   | 0.1513             | 0.2942          | 1.83         | very high   | medium        |
| 151  | 1.000   | 7.48       | 0.0423                   | 0.1577             | 0.3164          | 9.75         | very high   | high          |
| 153  | 1.000   | 17.89      | 0.0367                   | 0.1633             | 0.6567          | 2.98         | very high   | high          |
| 125  | 0.998   | 17.12      | 0.0254                   | 0.1746             | 0.4343          | 3.08         | very high   | high          |


**The depth spans a factor of 5.51 and must not be summarised as "shallow".** At one end node 125
averages `0.1746 mg/L` while below — genuinely marginal. At the other, **node 145 averages
`0.0600 mg/L`, less than a third of the threshold**, whenever it is under it. Calling both
"marginally below" is exactly the flattening the severity axis was introduced to stop, committed one
level down.

**And the severity band hides it.** Nodes 125, 145 and 153 are all *prolonged* (`E[D]` 12–24 h), yet
their cumulative deficits are `0.4343`, `1.6826` and `0.6567 mg/L·h` — a factor of **3.9 inside a
single band**. This is `E[D]`'s blind spot, stated in §12.5.1 as a reason for the depth column and
now demonstrated on real nodes: the band answers "for how long", and only the depth column answers
"how far below". A severity band must never be read as a statement about magnitude on its own.

So the honest summary of the ten is: they breach with near-certainty for a limited part of the
window, and their mean breach depth varies materially, from `0.025` to `0.14 mg/L`. On the breach
product they are indistinguishable from node 131, which is below the threshold for **46.6 of 48 h**;
the severity product separates them on duration, and the depth column separates them again.

**How to read the pair.** In this field the breach product never *under*-states severity and
over-states it for 10 of 59 consumer junctions. Neither ordering is the correct one: `P_min` is the
right axis for "is this node compliant", severity for "how much chlorine is missing and for how
long". This is the same lesson Step 8b reached from the other direction — a shortlist is not
meaningful without naming the metric that produced it — applied to the register itself.

**But "neither is correct" has to be a rule, not a stance.** A register exists to drive an action,
so if it emits one `control_measure` column while claiming to prefer neither axis, it has silently
made the breach axis operational. An earlier version of this register did exactly that: it sorted on
`risk_score_breach` and derived every control measure from the breach band alone. The rule is now
pre-declared:

```
governing band  =  the HIGHER of (breach band, severity band)
control_measure =  from the governing band
```

i.e. **neither axis alone may reduce an action**. The per-axis measures are emitted alongside it
(`control_measure_breach`, `control_measure_severity`) so the reason stays visible, and the register
sorts on `risk_score_governing`.

**In this baseline field the rule currently coincides with the breach axis**, because the severity
band never exceeds the breach band here — the governing band equals the breach band at all 92
junctions. That coincidence is *measured and reported* (`currently_equivalent_to_the_breach_axis`),
not assumed: the reachability table above shows the bound permits an inversion of 1 band, so the
rule is not vacuous in general and would bind on a field where it occurred.

**What it still is not.** Both products are conditioned on the calibrated ensemble and on scenario A;
neither is a measurement. And the severity axis inherits `E[D]`'s own blind spot, which is why the
depth column exists and must be read with it.

**Sampling priority is a classification-ambiguity index, not an optimal-monitoring criterion.** It
peaks where the assessment cannot classify a node either way and weights that by consequence. It
contains no Jacobian, no redundancy against the existing six monitors, no observation-noise model and
no expected posterior-variance reduction, so it must not be read as expected information gain or as
sensor placement. Name it `classification-ambiguity sampling index`. It uses the **discrete
consequence score (0–3), not raw demand**:

```
priority(n) = consequence_score(n) · P_min(n) · [1 − P_min(n)]     normalised by its maximum
priority(n) = 0  where demand = 0
```

so it peaks at `P_min = 0.5` — the junctions the assessment cannot classify either way.

**Product statement (use this wording).** The output is a **formal-likelihood calibration-conditioned
scenario projection** from the network model under stated assumptions (fixed hydraulics, first-order
kinetics, illustrative Arrhenius / ageing multipliers). It is **not** a sensor nowcast (no same-day
re-conditioning is performed), **not** a spatial measurement or geostatistical interpolation, **not**
a field forecast or dosing recommendation, and **not** a statement that water is safe. Chlorine
residual is one operational indicator among the WSP's preventive barriers.

**Review triggers (the assessment expires if):** sensor QA failure or reference-check
disagreement; hydraulic model or demand allocation revised; source/treatment regime change;
water temperature outside the 8–24 °C scenario range; a forecast heat episode; mains works
altering the assumed age profile; calibration record exceeding its approved age.

### 12.6 Findings for the Results / Discussion

1. **Temperature escalates risk materially but does not flatten the network.** `P_min>0.5`
   nodes rise 21 → 28 → 29 from 12 → 16 → 20 °C; demand at risk 8.0 % → 9.8 %; network-mean
   `E[A]` more than doubles (0.222 → 0.522 mg/L·h).
2. **Ageing at the same temperature adds a spatially selective increment** (`P_min>0.5` nodes
   29 → 31; network-mean `E[A]` 0.5218 → 0.579 mg/L·h), concentrated where old/average pipes
   control the supply path. The **`P_min`-classified set** is robust across the tested `α` range
   (30 / 31 / 31 nodes, 69.4 L/s at all three), while the severity magnitude scales with `α`
   (`E[A]` 0.550 → 0.579 → 0.599) — report it as a stress range, not as a measured effect.
   *"Set" here means the nodes classified `P_min > 0.5`, which is a breach-probability
   classification and not the `E[A]`-ranked hot-spot list of Step 10; the two are different
   products and Step 8b showed metrics of this kind can disagree.*
3. **The escalating nodes are almost all unmonitored** (**6 of the 7**, the exception being
   monitored node 231), which is an argument for information-led monitoring design rather than for
   more monitoring in general. Escalation is a change of `P_min` risk band between scenarios A
   and D, so this too is a breach-probability statement.
4. **Source dosing is an incomplete control.** +30 % inlet dose improves every continuous
   metric but does not return demand-at-risk to baseline. Under the idealisation of fixed
   hydraulics, first-order kinetics and proportional source scaling, dosing by `r` is exactly
   equivalent to testing against `C_crit / r`; it cannot buy back residence time.
5. **Scenario maps are planning products outside the calibrated regime.** The Arrhenius form,
   `T_ref` and `α_g` are assumptions; unlike the Step-11 LOO check, these projections cannot be
   verified against held-out chlorine observations. Absolute severity metrics are additionally
   horizon-dependent (paired warm-up test, §12.4).
6. **A breach-probability register and a severity register disagree, and in this field the
   disagreement runs one way.** Scoring `E[D]` on pre-declared absolute bands against the same
   consequence axis puts **10 of 59** consumer junctions in a *lower* risk band and **none** in a
   higher one; the ten all breach with certainty (`P_min = 1.000`) but only for 2.0–17.9 h and only
   a mean depth that **varies materially, 0.025–0.14 mg/L** — node 125 averages 0.175 mg/L while
   below and is genuinely marginal, node 145 averages 0.060 mg/L and is not. Within the single
   *prolonged* band the cumulative deficit still spans a factor of 3.9, so a severity band is a
   statement about duration and never about magnitude. The direction of disagreement is
   **empirical, not implied**: `E[D] ≤ T_window · P_min` holds (0 violations across 92 junctions)
   but bounds the continuous quantities, not the banded scores, whose scales are not aligned — it
   permits an inversion of one band, forbids one only once `P_min ≥ 0.80`, and the nearest consumer
   junction is **0.939 h** of `E[D]` from inverting.
7. **The ageing stress is nearly invisible to breach probability and clearly visible to severity.**
   Banding all four scenarios on the same absolute hours, the count of *persistent* consumers
   (below 0.2 mg/L for at least half the 48 h window) runs **1 → 4 → 6 → 8** from A to D, while the
   breach count runs 10 → 12 → 13 → **13** and saturates between C and D — at 20 °C most of the
   exposed set is already at `P_min > 0.8`, so the ageing stress cannot move it. A register built
   on breach probability alone would report that ageing adds almost nothing on top of a heatwave.
   This is the clearest justification in this log for carrying a second axis, and it only exists
   because the severity edges were declared in absolute hours rather than as quantiles (§12.5.1).
8. **A band-change escalation test is necessary but not sufficient, for two separate reasons.**
   Sixteen consumer junctions (52.3 L/s) sit at `P_min ≈ 1` in both A and D, so their likelihood
   band cannot move; node 243 goes to the **full 48 h window** below the threshold and node 15, the
   network's largest consumer, gains 11.6 h, and neither is flagged. `max(breach, severity)` does
   not rescue this — it **saturates**, which is a different property from the one it was adopted
   for. Separately, the multiplicative matrix **caps a minor consumer at `medium`** (score
   `5 × 1 = 5`), so **20 junctions (17.46 L/s)** can never be classified `high` however severe they
   become. That is a property of the `likelihood × consequence` form itself. Escalation at a
   saturated or capped junction has to be read on `ΔE[D]` and `ΔE[A]` (§12.3). Because neither product is the correct one, the register acts on the
   **higher** of the two bands, a pre-declared rule that currently coincides with the breach axis
   here and is reported as coinciding rather than assumed to (§12.5.1).

Outputs: `figures/step12_scenario_maps.png`, `figures/step12_ageing_delta.png`,
`figures/step12_summary.png`, `baseline_cache/step12_scenarios.json`,
`baseline_cache/step12_risk_register.csv`.

---



## Step 13 — known-answer test: are the coefficients realised as we set them? (Priority-2 #7)

Every number in this log rests on an assumption that was never tested directly: that a coefficient
written in `1/day` or `m/day`, divided by 86400 and handed to WNTR, arrives in EPANET as that
coefficient. The review found a unit error in the reported figures (§3.4), so the assumption is not
self-evidently safe. This is the single-pipe analytic check the review asked for
(`step13_known_answer.py`): one reservoir, one 1000 m pipe of 0.3 m diameter, a constant 0.05 m³/s
demand, `C0 = 1.0 mg/L`, 24 h.

**Arm 1 — pure bulk decay against the analytic solution.** With no wall reaction the concentration
at the far end once the front has passed is exact.

**Formula**:

```
C = C0 · exp(k_b · t_res),     t_res = L / v
```

**Properties**:

- Units: `k_b` in 1/s after conversion, `t_res` in s, `C` in mg/L.
- The test is a *unit* test in the literal sense: a wrong conversion factor changes `C` by a
constant relative amount at every `k_b`.


| `k_b` (1/day) | `t_res` (h) | EPANET `C` | analytic `C` | relative error |
| ------------- | ----------- | ---------- | ------------ | -------------- |
| −0.50         | 0.3927      | 0.991846   | 0.991852     | 6.7 × 10⁻⁶     |
| −1.00         | 0.3927      | 0.983744   | 0.983771     | 2.7 × 10⁻⁵     |
| −2.00         | 0.3927      | 0.967699   | 0.967805     | 1.1 × 10⁻⁴     |


Worst relative error **1.1 × 10⁻⁴** against a 10⁻³ tolerance. The error *grows with the amount of
decay* across a fourfold range of `k_b`, which is the signature of the quality solver's time
discretisation; a wrong conversion factor would instead show a constant relative offset. So the
conversion and the units are confirmed.

**Arm 2 — wall decay: monotonicity and a bound, not an exact value.** EPANET's first-order wall
reaction is limited by mass transfer as well as by `k_w`, through an overall rate constant

```
k_overall = k_b + (4 / D) · (k_w · k_f) / (k_w + k_f)
```

so `C` cannot be predicted from `k_w` alone without the mass-transfer coefficient `k_f`, which
depends on the flow regime. What *is* assertable was checked at `k_b = −0.5` 1/day over
`k_w ∈ {0, −0.1, −0.5, −1, −2}` m/day: `k_w = 0` reproduces the pure-bulk run **exactly**, `C`
falls monotonically as `k_w` strengthens, and `C` stays above the `k_f`-free limit in which the wall
term would be `(4/D)·k_w`. All three hold. That is enough to establish the wall coefficient is
applied with the intended sign and scale, and it is stated as such rather than as an exact analytic
match.

**Arm 4 — the concentration unit, against WNTR's own converter.** Arms 1–3 all test the *time*
conversion, which was documented in the inherited material and correct from the start. The
*concentration* conversion was not applied at all until Step 15, and nothing here would have caught
it: under first-order kinetics the analytic check is a ratio, so it passes on either scale. This arm
therefore pins the unit against `wntr.epanet.util.to_si` / `from_si` rather than against our belief
about them — `to_si(1.0 mg/L) = 0.001 kg/m³`, `from_si(1.0 kg/m³) = 1000 mg/L`, the helpers in
`wq_common` agree with both. Arm 5 then reads the written `.inp` back and pins the asymmetry the
description keeps drifting on: the source is written **converted** (`0.001 kg/m³ → River 1.0`) while
`TOLERANCE` is written **verbatim** (`1e-05`), so EPANET reads the tolerance in mg/L. A first-order
coefficient carries no mass unit, which is why arms 1–3 were unaffected by the error arm 4 exists to
prevent.

**Arm 3 — the conversion helper itself.**
`per_day_to_per_second(−1.0) = −0.0000115741 per second`, i.e. exactly `−1/86400`.

All checks pass. One paragraph of the thesis should state this: the toolchain realises the
bulk coefficient scale as intended and applies wall reaction with the intended sign and physically
bounded monotonic decay; the wall arm is **not** an analytic verification of the exact numerical
wall-coefficient scale under mass-transfer limitation. The calibration results are about the inverse
problem and not about a units mistake.

---



## Step 14 — repeated-noise calibration: bias, CRLB consistency, interval coverage

Every formal posterior quoted elsewhere in this log (except where a step already averages over noise)
comes from **one** noise realisation. That is enough to show the posterior is narrow and that its
width sits close to the Case-A CRLB, but it cannot establish the two properties that actually matter
for calling an estimator calibrated:

- **bias** — does the posterior mean land on the truth on average?
- **coverage** — does a nominal 90%/95% interval contain the truth about 90%/95% of the time?

Both are repeated-sampling statements. The 8192-member candidate library is noise-independent, so
this step costs no EPANET runs: only the observations are redrawn (`step14_repeated_noise.py`,
`N = 100` seeds starting at `20000`, deliberately disjoint from the `42…71` block used elsewhere).

**Shared-library caveat (bounds the whole step).** Every realisation is weighted against the *same*
Sobol prediction library, so Monte Carlo errors are not independent across realisations, and a sharp
likelihood resolves interval endpoints only through its ESS. Formal ESS median / p5 / min ≈
**158 / 122 / 113** — adequate for the reported coverages, not luxurious.

Headline calibration (old / average / new): empirical SD / CRLB ≈ 1.04 / 1.06 / 1.12; nominal 90% intervals cover 0.89 / 0.88 / 0.85.

**Primary — formal censored** (100/100 valid):


| coef    | truth | E[post. mean] | bias    | bias / emp. SD | emp. SD of mean | within-run SD | within/emp | CRLB   | emp/CRLB | cov90 | cov95 |
| ------- | ----- | ------------- | ------- | -------------- | --------------- | ------------- | ---------- | ------ | -------- | ----- | ----- |
| old     | −1.00 | −1.0055       | −0.0055 | −0.06          | 0.0981          | 0.0960        | 0.98       | 0.0947 | 1.04     | 0.89  | 0.97  |
| average | −0.10 | −0.1021       | −0.0021 | −0.15          | 0.0144          | 0.0137        | 0.95       | 0.0135 | 1.06     | 0.88  | 0.96  |
| new     | −0.05 | −0.0498       | +0.0002 | +0.02          | 0.0089          | 0.0079        | 0.89       | 0.0080 | 1.12     | 0.85  | 0.91  |


**Comparator — informal GLUE** (threshold 0.107): bias/emp. SD up to −1.76 (`average`); within-run SD
is **~2.6–2.9×** the formal width (0.2489 / 0.0393 / 0.0226 against 0.0960 / 0.0137 / 0.0079 m/day
for old / average / new); cov90 ≈ 0.98–0.99. High coverage is bought by width, not by accuracy.

Repeated-noise calibration

Findings:

1. **Under the synthetic generative model the formal posterior mean is essentially unbiased** — every
   bias is a small fraction of the estimator's own sampling SD.
2. **The formal posterior spread is locally consistent with the Case-A CRLB** (emp/CRLB ≈ 1.04–1.12;
   within-run SD / CRLB ≈ 1.00–1.01). This is a consistency check, not a proof of frequentist
   efficiency: both quantities use the same model, truth and error assumptions.
3. **Nominal intervals are close to calibrated** at 90%/95% for old and average; `new` is slightly
   under-covered at 90% (0.85), consistent with its smaller ESS footprint and the shared-library
   coarseness.
4. **The informal comparator is not a calibrated tight interval** — it over-covers by being several
   times wider, which is the same inefficiency diagnosed in Steps 1–3 and 7.

Outputs: `figures/step14_repeated_noise.png`, `baseline_cache/step14_repeated_noise.json`.
---

## Step 15 — the chlorine concentration unit was wrong by 1000×, and what that did (and did not) change

**This is not a review comment.** It was found by self-audit after the revision was otherwise
complete, and it is the largest single error in the project's history: every chlorine concentration
in the repository was on a scale 1000× too high, for the whole of its life, until this step.

**What was wrong.** WNTR's Python API stores concentration internally in **kg/m³**.
`options.quality.inpfile_units = "mg/L"` governs only how the `.inp` is read and written — it does
not change what a value assigned through the API means. Assigning `initial_quality = 1.0` intending
"1 mg/L" therefore made EPANET simulate **1000 mg/L**, and reading `run_sim()` output as mg/L
repeated the error on the way out. Measured directly, before the fix WNTR wrote this into the file
it hands EPANET:

```
[QUALITY]                       [OPTIONS]
River      1000.0               QUALITY    Chlorine mg/L
Lake       1000.0               TOLERANCE  0.01
1           500.0
```

So the source was 1000 mg/L, the tanks 500 mg/L, the sensor noise σ 100 mg/L and the "0.2 mg/L
operational threshold" 200 mg/L. After the fix the same probe shows `River 1.0`, `1 0.5`,
`TOLERANCE 1e-05`.

**Why it survived so long — three reasons, and none of them is carelessness alone.**

1. **First-order kinetics are linear in C.** The field scales exactly with the source, and σ and the
   risk threshold were expressed on the same wrong scale, so **every ratio the analysis depends on
   was preserved**: σ/C₀ = 0.1 and C_MIN/C₀ = 0.2 either way. No internal check could fail.
2. **The reaction coefficients are not affected, so the known-answer test passed.** `k_b` is
   `1/time` and first-order `k_w` is `length/time`; neither carries a mass unit, so the factor
   cannot propagate into them. Verified against the file WNTR writes: `GLOBAL BULK -0.5000` and
   `WALL -3.2808 / -0.3281 / -0.1640 ft/day`, which are exactly the intended −0.5 day⁻¹ and
   −1.0 / −0.1 / −0.05 m/day. Step 13's analytic arm therefore passed at 10⁻⁴ throughout.
3. **The factor had already been found, and mis-diagnosed.** The inherited starter notebook records,
   as a measured fact, that a **zero-order** global bulk coefficient needs `×1000`
   (`ZERO_ORDER_BULK_CAL = 1000.0`). That is the same kg/m³ fact: a zero-order coefficient carries a
   `mass/volume/time` dimension, so it shows the factor. It was recorded as a quirk of zero-order
   bulk rather than as a property of the concentration unit, so it was never applied to
   concentrations themselves — where, under first order, it has no symptom. The one place the factor
   was visible is the one place it was patched.

### 15.2 The tolerance is a second, separate change — state it as one

The fix has two parts and they are **not the same kind of change**, which an earlier version of this
section got wrong by calling the second one a unit conversion.

EPANET's water-quality `TOLERANCE` is an **absolute** concentration below which two water parcels
are treated as identical. Unlike a concentration, WNTR does **not** convert it: `write_inpfile` puts
`options.quality.tolerance` into the `.inp` verbatim and EPANET reads it in the file's declared
quality unit. Reading the written file back (Step 13, arm 5) shows both behaviours side by side:

```
[OPTIONS]  QUALITY   Chlorine mg/L
[OPTIONS]  TOLERANCE 1e-05          <- written verbatim  -> EPANET reads 1e-5 MG/L
[QUALITY]  River     1.0            <- written converted from 0.001 kg/m^3
```

So `QUALITY_TOLERANCE = 1e-5` is **1e-5 mg/L**, and calling it "1e-5 kg/m³" or "0.01 mg/L expressed
internally" is wrong on both counts. It is **1000× stricter than EPANET's 0.01 mg/L default**, and
that is an active numerical choice which has to be declared as one.

**Why this value was chosen.** The superseded configuration ran the 0.01 default against a
1000 mg/L source — a tolerance/source ratio of 1e-5. Correcting the concentrations while leaving
0.01 in place would have left 0.01 against 1 mg/L, a ratio of 1e-2: a solver three orders of
magnitude coarser than the one every earlier result was computed on. Holding the *ratio* fixed is
what lets Step 15 compare the corrected and superseded runs as the same numerical experiment and
attribute their agreement to the units alone. A stricter tolerance is very likely also *more*
accurate — it resolves the Lagrangian transport more finely — but that is not the argument being
made here, and this log does not claim it: the argument is **comparability**, and the cost is
runtime.

**Formula** — the conversion, applied at the model boundary and nowhere else:

```
internal (WNTR/EPANET)  =  reported (mg/L) / 1000        build_model, on the way in
reported (mg/L)         =  internal × 1000              run_model / simulate_chlorine, on the way out
QUALITY_TOLERANCE       =  1e-5 mg/L (chosen, not converted; EPANET default 0.01 mg/L)
```

**Properties**:

- Everything outside `wq_common` is in mg/L, so no step script changed: `σ = 0.1` and
  `C_MIN = 0.2` still read as they always did, and are now true.
- AGE runs are exempt: their "initial quality" is a time, not a concentration.

### Step 15 — measured effect on 256 leading Sobol candidates

Three arms on identical candidates (`step15_unit_equivalence.py`): **legacy** (the superseded
configuration, whose raw output is what the old cache stored and this log quoted as mg/L),
**corrected** (current), and **units_only** (units fixed, tolerance left at 0.01 — the fix done
wrong, kept as the counterfactual).

| comparison | max abs diff | max rel diff | points below 0.2 | hours-below cells that differ |
|---|---|---|---|---|
| corrected vs legacy | 1.19 × 10⁻⁷ mg/L | **2.3 × 10⁻⁷** | 79391 → 79391 | **0 of 23552** |
| units_only vs legacy | 2.09 × 10⁻² mg/L | **17.8 %** | 79391 → 79301 | **440 of 23552** |

So the correction, done properly, leaves the reported field **unchanged to 2.3 × 10⁻⁷ relative** —
floating point. Done without the tolerance it **introduces a 17.8% maximum relative error** and moves
1.9% of the risk-duration cells. The second row is why the tolerance is part of the fix and not an
afterthought.

### Step 15b — verification on the full library, and on this log

**This section used to quote numbers no artifact held.** The three bullets below were originally
produced by throwaway scripts during the correction, which breaks this repository's own rule that
every number comes from a script and a named artifact. `step15b_full_regression.py` now performs the
whole comparison and writes `step15_full_regression.json`. Running it did not merely confirm the old
bullets — **two of the three were wrong**, and the corrected figures are below.

The comparison has two halves with different reproducibility. The **full-library** half is live: the
corrected arm is read from `baseline_cache/baseline.npz` — the cache the rest of the pipeline
actually consumes, so this checks the stored artifact rather than a fresh copy of it — and only the
legacy arm is re-simulated (8192 EPANET runs). The **artifact** and **log** halves are historical,
read from git at two recorded commits: `64ec7d3` (the last state before the correction) and
`3427a9d` (the state that recorded it). They cannot be re-derived from the working tree, because
later work — the demand-pattern fix, the register's severity axis — legitimately moved numbers
afterwards. The two SHAs are stored in the artifact, and if the history is ever rewritten the script
exits with an error instead of quietly reporting nothing.

- **`baseline.npz`, all 8192 candidates**: `C_all` max absolute difference `1.2 × 10⁻⁷ mg/L`, max
  **relative difference 2.3 × 10⁻⁷**; `RMSE` 7.2 × 10⁻⁸; `loglik_censored` 1.3 × 10⁻⁷;
  `loglik_iid` 1.4 × 10⁻⁷ (full precision in the artifact). Both arms are driven by the same Sobol
  draw set, whose SHA-256 is recorded so a silently regenerated design is detectable. *(The earlier
  bullet quoted 9.2 × 10⁻⁸ and 1.7 × 10⁻⁷ for the second and third of these; neither matches what
  the tracked script measures.)*
- **26 artifacts, 4915 numeric fields** (3 excluded as runtimes and hashes): **34** fields moved by
  more than 1e-4 relative — not the five the earlier bullet claimed. The relative screen
  over-reports badly here, so they are banded by how large the quantity itself is: 4 have a
  magnitude below 1e-4 (a relative move of a number that is itself ~1e-6 is one ulp of the
  simulation), 21 move by less than 1e-4 in absolute terms, and 9 move by more. Of those 9, the
  largest is an **optimiser iteration count** (`n_epanet_evaluations` 7770 → 7792); the rest are
  finite-difference convergence diagnostics, `half_width_change_vs_grid_frac` ratios, a
  max-vs-p95 water-age diagnostic and two Kendall τ values in the fifth decimal. **No estimate,
  interval, ranking or risk number is among them.** This is the expected signature of a change in
  *numerical resolution*: what moves is the machinery that measures convergence, not the results.
- **This log, 3582 numbers at `3427a9d`** (3473 at `64ec7d3`; the sequences are aligned with
  `difflib` so inserted prose does not report every later number as changed): **three** numbers
  changed, not two, plus four insertions which are the new Step 15 text itself. The three are
  `0.991845 → 0.991846` (the EPANET value in the single-pipe analytic check, sixth decimal),
  `6.8 → 6.7 × 10⁻⁶` — which is *that same value's own relative error*, i.e. the same fact appearing
  twice in one table and counted once — and `7770 → 7792` EPANET evaluations, an iteration count.
  *(The earlier bullet's "1720 checked numbers" was the validator's count of numbers it could match
  to an artifact, not the number of numbers in the log; the two were conflated.)*

**What this means, stated carefully.** The correction changed the *meaning* of the numbers, not the
numbers. Every parameter, identifiability, structural-error, sensor-bias, `k_b`, risk and scenario
conclusion in this log stands **verbatim**, and that is measured rather than argued. What was wrong,
and is now right, is the physical interpretation: the sensor σ, the `0.2 mg/L` threshold, the dosing
levels and every `mg/L` and `mg/L·h` in the risk tables were, until this step, statements about a
1000 mg/L system. Any engineering claim that depends on the absolute concentration — required sensor
accuracy, the operational threshold, the dosing evaluation — was **not** supported before the fix and
**is** supported now.

**Guarded against return.** `step13` arm 4 pins the concentration unit against WNTR's own
`to_si`/`from_si` rather than against our belief about them, and **arm 5 reads the written `.inp`
back** and asserts the asymmetry that the prose of this section originally got wrong: the source is
written converted (`1.0` mg/L from `0.001` kg/m³) while `TOLERANCE` is written verbatim (`1e-05`,
read by EPANET as mg/L). The unit convention and the tolerance are also part of the provenance
config hash, so a cache built under a different convention can no longer be mistaken for this one.

**One label is still wrong, and correcting it is deliberately deferred.** The hashed config records
the tolerance under the key `quality_tolerance_kg_m3`. The *value* (`1e-5`) is correct and is what
every result in this log was computed with; only the unit in the key **name** is wrong — it is mg/L.

Renaming it changes `config_sha256` (`c97826ac…` → `607f2034…`). No stored result would become
wrong, because no step JSON embeds that hash — but the six `*.key.json` files would then record a
superseded config, so the keyed arrays (`step4d`, `step7b`, `step7c`, `step8b`) would silently cold-
rebuild on next use and `step9` would raise until `step7c` had been re-run. That is a re-run whose
entire output is numerically identical to what is already here. It is therefore **not** done as its
own task; it belongs with the clean-tree release re-run, which is separately required anyway. It is
recorded here so the discrepancy is not silent, which is the point of writing it down rather than
paying for it twice.

Outputs: `baseline_cache/step15_unit_equivalence.json`.
