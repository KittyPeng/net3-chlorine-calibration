# Chlorine wall-decay calibration on EPANET Net3

Code and cached results for the MSc research paper *From Parameter Identifiability to
Operational Risk: Uncertainty-Aware Calibration of Grouped Chlorine Wall Decay in EPANET Net3*
(CIVE70058, Department of Civil and Environmental Engineering, Imperial College London).

Every number in the paper is the output of a script in `Net3/`, computed from the frozen network
in `models/net3_frozen/` and stored as JSON in `Net3/baseline_cache/`. The figures and tables are
generated from those artifacts directly. Numbers quoted in the prose are transcribed from them by
hand, which is why `Net3/validate_artifacts.py` pins a registered list of those sentences to the
named JSON path each one came from. The checks below let a reader confirm this rather than take
it on trust.

## Install

The environment is pinned exactly, not loosely. This study documents behaviour that depends on
the WNTR version, so a different version can change reported digits.

```bash
conda env create -f environment.yml
conda activate water-supply
```

or, with a virtual environment instead of conda:

```bash
python -m venv .venv && source .venv/bin/activate    # Python 3.12 or newer
pip install -r requirements.txt
```

Pinned: Python 3.13.12, numpy 2.4.2, scipy 1.17.1, pandas 2.3.3, matplotlib 3.10.8, wntr 1.4.0.
The wntr wheel ships the EPANET engine, so that pin fixes **EPANET 2.2.0** as well
(`ENgetversion` reports 20200). `environment.lock.yml` holds the full transitive solve if an
exact rebuild is needed.

## Verify

Three things can be checked without running the study.

```bash
pytest -q tests/                      # the Section 2.2.4 verification checks
python Net3/provenance.py --check     # environment and frozen model against the cached results
python Net3/validate_artifacts.py     # the written numbers against the JSON they came from
```

`tests/` restates the paper's Section 2.2.4 checks as assertions: the frozen `Net3.inp` matches
the SHA-256 recorded with the results; a single pipe with bulk decay only reproduces the analytic
first-order solution; the wall reaction has the right sign, responds monotonically, and saturates
as mass transfer becomes the limiting step; and concentrations round-trip through WNTR's internal
kg/m³ representation, which is the unit error the paper documents. They build their own networks
and need no cached artifact.

`validate_artifacts.py` runs fourteen checks. The one with real teeth is a registry of 100 claims,
each pinning one written figure to one named path inside one JSON file. This repository ships
`RESULTS_LOG.md` but not the manuscript, so 30 of the 100 are checked here and the run names the
documents it could not cover rather than reporting a clean pass over a smaller set. Its docstring
sets out what a passing run does and does not establish; it is worth reading before quoting one.

## Reproduce

The scripts run in order and each writes its result to `Net3/baseline_cache/`. Step 1 rebuilds the
prediction library that the rest depend on, which is the only slow step and the only one whose
output is not tracked here.

```bash
cd Net3
python step1_freeze_baseline.py       # ~8192 forward runs; writes baseline.npz
python step3_threshold_sensitivity.py
...                                   # see Appendix I of the paper for the full mapping
python paper_figs.py                  # figures and main-text tables
python appendix_tables.py             # appendix tables
```

`Net3/baseline_cache/*.json` is tracked because the paper quotes those values directly.

`baseline.npz`, the 8192-candidate prediction library, is not in the repository tree. It is
118 MB, past GitHub's per-file limit, and twenty-four scripts read it, including `paper_figs.py`.
It is **attached to the `v1.0` release** instead, and `cache_manifest.json` records its SHA-256 as
`baseline_npz_sha256`, so a downloaded or regenerated copy can be checked against the one every
published result was computed from:

```bash
shasum -a 256 Net3/baseline_cache/baseline.npz
```

The `.npy` prediction caches for the displaced-prior and bulk-decay arms are attached to the same
release and carry their own digests in the `.npy.key.json` sidecars, and their own steps
regenerate them. Everything in `Verify` above runs without any of them. Only re-plotting the
figures and re-running the candidate-level calculations needs them.

## Layout


| Path                         | What it holds                                                         |
| ---------------------------- | --------------------------------------------------------------------- |
| `Net3/step*.py`              | one experiment each; Appendix I of the paper maps sections to scripts |
| `Net3/wq_common.py`          | shared configuration, the weighting rules, the forward model wrapper  |
| `Net3/provenance.py`         | records and checks the environment that produced the cache            |
| `Net3/validate_artifacts.py` | cross-checks the written numbers against the artifacts                |
| `Net3/paper_figs.py`         | the seven figures and the main-text table bodies                      |
| `Net3/appendix_tables.py`    | the twenty-seven appendix tables                                        |
| `Net3/baseline_cache/`       | result artifacts, JSON, quoted directly by the paper                  |
| `models/net3_frozen/`        | the frozen `Net3.inp`, hash-checked on import                         |
| `tests/`                     | the Section 2.2.4 verification checks                                 |




## Notes

Concentrations are handled in mg/L throughout and converted to WNTR's internal kg/m³ on input and
back on output. The EPANET absolute water-quality tolerance is set to 1×10⁻⁵ mg/L rather than left
at the 0.01 mg/L default; Section 2.2.4 of the paper explains why, and the tests check the
round trip.

One defect is worth knowing about before reusing this toolchain. Concentrations are converted on
file write while some reaction coefficients are not, so a **zero-order** bulk coefficient passed
through the WNTR 1.4 API reaches EPANET unscaled while the concentrations around it have been
scaled by 1000, giving reaction rates a thousand times smaller than nominal. This study uses
first-order kinetics throughout, and a first-order coefficient carries no mass unit and is
unaffected, which is why the checks in `tests/` and Appendix B of the paper pass either way.

## Use of generative AI

Generative AI (Claude Opus 5 through Claude Code, and Cursor) was used during this project for:

- manuscript restructuring, language editing and proofreading;
- revising the figure and appendix-table generation scripts, `Net3/paper_figs.py` and `Net3/appendix_tables.py`;
- code review and revision.

The experimental design and the `step*.py` computation workflow are the author's. All AI-assisted
suggestions and code changes were reviewed, tested and approved by the author.

No reported numerical result was directly generated by a language model. Every reported value in
the paper is the output of deterministic EPANET/WNTR and Python computations, stored in
`Net3/baseline_cache/`. The paper carries the full statement after its reference list.