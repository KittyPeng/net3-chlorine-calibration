"""Step 15b: the full-library and whole-repository regression behind Step 15's claims.

WHY THIS EXISTS. Step 15 measures the concentration-unit correction on 256 leading Sobol candidates
and stores that in step15_unit_equivalence.json. The log then makes three STRONGER claims — the
whole 8192-candidate library agrees to 2.3e-07, ~4000 numeric fields across every artifact moved in
five places, and two of 1720 numbers in the log changed. Those were produced by throwaway scripts
during the correction, so by this repository's own standard ("every number comes from a script and a
named artifact") the evidence chain was open. This step closes it.

THREE PARTS, and they are reproducible in two different ways.

  A. FULL CANDIDATE LIBRARY — live. The corrected arm is not re-simulated: it is read from
     baseline_cache/baseline.npz, i.e. the cache the rest of the pipeline actually consumes, which
     makes this a check on the stored artifact and not on a fresh copy of it. Only the LEGACY arm is
     re-run (8192 EPANET runs, ~5 min), and C_all / RMSE / both log-likelihoods are compared.

  B. ARTIFACT REGRESSION — historical, read from git. Every numeric leaf of every artifact is
     compared between the commit BEFORE the correction and the commit that recorded it. This cannot
     be re-derived from the working tree, because later work (the demand-pattern fix, the severity
     axis) legitimately moved numbers afterwards; the two commits are the experiment.

  C. LOG REGRESSION — historical, same two commits, same number-extraction regex the validator uses,
     aligned with difflib so an insertion does not report every later number as changed.

Parts B and C are pinned by the two commit SHAs recorded in the output. If the history is ever
rewritten they stop being reproducible, and the script says so by failing rather than by skipping.
"""
import difflib
import json
import os
import re
import subprocess
import sys
import time

import numpy as np

import wq_common as B
import provenance
import wntr

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "baseline_cache")
REPO = os.path.dirname(HERE)

# The two commits the historical half of this regression is defined against. PRE is the last state
# before the concentration-unit correction; POST is the state in which it was recorded. They are
# arguments of the experiment, not implementation details, so they are stored in the output.
PRE = "64ec7d3"
POST = "3427a9d"

REL_TOL = 1e-4          # a field must move by more than this, relatively, to be reported
# Fields that are expected to differ between any two runs and say nothing about the correction.
EXCLUDE = re.compile(r"(?i)(runtime|elapsed|generated_utc|_sha256|sha256_|git/|/git$|"
                     r"uncommitted|tree_sha|commit)")


def git_show(ref, path):
    """File contents at a commit, or None if it did not exist there."""
    r = subprocess.run(["git", "-C", REPO, "show", f"{ref}:{path}"],
                       capture_output=True, text=True)
    return r.stdout if r.returncode == 0 else None


def git_resolve(ref):
    r = subprocess.run(["git", "-C", REPO, "rev-parse", ref], capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"commit {ref} is not in this repository — the historical half of this regression "
                 f"is defined against it and cannot be reproduced without it")
    return r.stdout.strip()


def leaves(obj, prefix=""):
    """(path, value) for every numeric leaf, so two artifacts can be compared field by field."""
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from leaves(v, f"{prefix}/{k}")
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from leaves(v, f"{prefix}[{i}]")
    elif isinstance(obj, bool):
        return                                  # a flag is not a measurement
    elif isinstance(obj, (int, float)):
        yield prefix, float(obj)


def rel(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-12)


# =============================== A. full candidate library ===============================
print("=== Step 15b: full-library and whole-repository regression ===\n")
print("A. full candidate library — corrected arm read from baseline.npz, legacy arm re-run live")

npz = np.load(os.path.join(OUT, "baseline.npz"), mmap_mode="r")
C_corr = npz["C_all"]                       # (8192, 49, 92) float32, as the pipeline stores it
RMSE_corr, LLC_corr, LLI_corr = npz["RMSE"], npz["loglik_censored"], npz["loglik_iid"]
ALL_NODES = [str(n) for n in npz["all_nodes"]]
obs_glue = npz["obs_glue"]
mon_pos = npz["mon_pos"]
draws = B.prior_draws()
S = (draws["old"], draws["avg"], draws["new"])
N = C_corr.shape[0]
assert len(S[0]) == N, (len(S[0]), N)

# Both arms are driven by the SAME draws array, so "the designs agree" is true by construction and
# not worth asserting. What is worth recording is which design this was: a hash of the three draw
# vectors, so a later run that silently regenerated a different Sobol sequence is detectable.
import hashlib as _hashlib
_h = _hashlib.sha256()
for _s in S:
    _h.update(np.ascontiguousarray(np.asarray(_s, dtype=np.float64)).tobytes())
sobol_draws_sha256 = _h.hexdigest()


def legacy_field(i):
    """One candidate under the SUPERSEDED configuration: raw assignment, EPANET default tolerance."""
    wn = B.build_model(B.KB_FIXED, 0.0,
                       pre_run=B.make_kw_hook(S[0][i], S[1][i], S[2][i]))
    for r in wn.reservoir_name_list:
        wn.get_node(r).initial_quality = B.INLET_CHLORINE_MGL
    for tk in wn.tank_name_list:
        wn.get_node(tk).initial_quality = B.TANK_INIT_MGL
    wn.options.quality.tolerance = 0.01
    res = wntr.sim.EpanetSimulator(wn).run_sim(file_prefix="s15b")
    # legacy reported the raw internal array as if it were mg/L, which is the whole error
    return res.node["quality"][ALL_NODES].values[B.WARMUP_H:]


# Part A costs 8192 EPANET runs and depends only on the cache, the draws and the configuration.
# Parts B and C are re-reads of git and cost nothing, so a re-run to change how the COMPARISON is
# reported should not have to re-simulate. Reuse is keyed on config + draws, never on file existence.
OUTPATH = os.path.join(OUT, "step15_full_regression.json")
cached_a = None
if os.path.exists(OUTPATH):
    try:
        prev = json.load(open(OUTPATH))
        pa = prev.get("full_library", {})
        if (prev.get("config_sha256") == provenance.config_sha256()
                and pa.get("sobol_draws_sha256") == sobol_draws_sha256
                and pa.get("n_candidates") == int(N)):
            cached_a = pa
            print("   reusing part A from the existing artifact (same config and same draws)")
    except (json.JSONDecodeError, OSError):
        cached_a = None

max_abs_C = max_rel_C = 0.0
max_rel_RMSE = max_rel_LLC = max_rel_LLI = 0.0
t0 = time.time()
for i in range(0 if cached_a is None else N):
    leg = legacy_field(i)
    cor = np.asarray(C_corr[i], dtype=np.float64)
    d = np.abs(cor - leg)
    max_abs_C = max(max_abs_C, float(d.max()))
    max_rel_C = max(max_rel_C, float((d / np.maximum(np.abs(cor), 1e-12)).max()))

    leg_mon = leg[:, mon_pos]
    r_leg = float(np.sqrt(((leg_mon - obs_glue) ** 2).mean()))
    max_rel_RMSE = max(max_rel_RMSE, rel(float(RMSE_corr[i]), r_leg))
    max_rel_LLC = max(max_rel_LLC, rel(float(LLC_corr[i]),
                                       float(B.log_censored(leg_mon[None], obs_glue)[0])))
    max_rel_LLI = max(max_rel_LLI, rel(float(LLI_corr[i]),
                                       float(B.log_gaussian(leg_mon[None], obs_glue)[0])))
    if (i + 1) % 1024 == 0:
        print(f"   {i + 1}/{N}  ({time.time() - t0:.0f}s)  C_all max rel so far {max_rel_C:.2e}")

part_a = cached_a if cached_a is not None else {
    "n_candidates": int(N),
    "corrected_arm": "read from baseline_cache/baseline.npz — the cache the pipeline consumes",
    "legacy_arm": "re-run live: initial_quality assigned raw, EPANET default tolerance 0.01",
    "C_all_shape": list(C_corr.shape),
    "C_all_max_abs_diff_mg_L": max_abs_C,
    "C_all_max_rel_diff": max_rel_C,
    "RMSE_max_rel_diff": max_rel_RMSE,
    "loglik_censored_max_rel_diff": max_rel_LLC,
    "loglik_iid_max_rel_diff": max_rel_LLI,
    "sobol_draws_sha256": sobol_draws_sha256,
    "sobol_draws_note": "both arms are driven by this same draw set, so the comparison isolates the "
                        "configuration; the hash pins which design that was",
    "runtime_s": round(time.time() - t0, 1),
}
print(f"\n   C_all      max abs {part_a['C_all_max_abs_diff_mg_L']:.3e} mg/L, "
      f"max rel {part_a['C_all_max_rel_diff']:.3e}")
print(f"   RMSE       max rel {part_a['RMSE_max_rel_diff']:.3e}")
print(f"   loglik_cen max rel {part_a['loglik_censored_max_rel_diff']:.3e}   "
      f"loglik_iid max rel {part_a['loglik_iid_max_rel_diff']:.3e}")

# =============================== B. artifact regression ===============================
print(f"\nB. artifact regression, {PRE} (before the correction) -> {POST} (as recorded)")
pre_sha, post_sha = git_resolve(PRE), git_resolve(POST)
listing = subprocess.run(["git", "-C", REPO, "ls-tree", "--name-only", POST,
                          "Net3/baseline_cache/"], capture_output=True, text=True).stdout.split()
art_names = sorted(os.path.basename(p) for p in listing if p.endswith(".json")
                   and not p.endswith(".npy.key.json"))

per_artifact, changed_fields, n_fields_total, n_excluded = [], [], 0, 0
for name in art_names:
    path = f"Net3/baseline_cache/{name}"
    a, b = git_show(PRE, path), git_show(POST, path)
    if b is None:
        continue
    if a is None:
        per_artifact.append({"artifact": name, "status": "new at POST", "n_fields": None,
                             "n_changed": None})
        continue
    try:
        da, db = json.loads(a), json.loads(b)
    except json.JSONDecodeError:
        per_artifact.append({"artifact": name, "status": "unparseable", "n_fields": None,
                             "n_changed": None})
        continue
    la, lb = dict(leaves(da)), dict(leaves(db))
    n_here = n_ch = 0
    for k, vb in lb.items():
        if EXCLUDE.search(k):
            n_excluded += 1
            continue
        n_here += 1
        if k in la and rel(la[k], vb) > REL_TOL:
            n_ch += 1
            changed_fields.append({"artifact": name, "field": k, "pre": la[k], "post": vb,
                                   "rel_diff": rel(la[k], vb),
                                   "abs_diff": abs(la[k] - vb),
                                   "magnitude": max(abs(la[k]), abs(vb))})
    n_fields_total += n_here
    per_artifact.append({"artifact": name, "status": "compared", "n_fields": n_here,
                         "n_changed": n_ch})

changed_fields.sort(key=lambda r: -r["rel_diff"])
# A relative screen alone over-reports: most of what moves are convergence diagnostics and
# finite-difference ratios whose own magnitude is 1e-6 or smaller, where a 1% relative move is one
# ulp of the underlying simulation. Band them so the summary cannot be read as "34 results moved".
def _band(r):
    if r["magnitude"] < 1e-4:
        return "magnitude < 1e-4 (discretisation-level quantity)"
    if r["abs_diff"] < 1e-4:
        return "absolute move < 1e-4"
    return "absolute move >= 1e-4"


bands = {}
for r in changed_fields:
    bands.setdefault(_band(r), []).append(f"{r['artifact']}:{r['field']}")
print(f"   {len(art_names)} artifacts, {n_fields_total} numeric fields compared "
      f"({n_excluded} excluded as runtimes/hashes/provenance)")
print(f"   {len(changed_fields)} field(s) moved by more than {REL_TOL:g} relative")
for r in changed_fields[:10]:
    print(f"     {r['artifact']}:{r['field']}  {r['pre']:.6g} -> {r['post']:.6g} "
          f"(rel {r['rel_diff']:.2e})")

# =============================== C. log regression ===============================
print(f"\nC. log regression, same two commits")
# the validator's own number regex, so "a number in the log" means the same thing in both places
NUMBER = re.compile(r"(?<![\w.])(?<!\w-)(-?\d+(?:\.\d+)?)(?![\w.]*\d)")
log_a, log_b = git_show(PRE, "Net3/RESULTS_LOG.md"), git_show(POST, "Net3/RESULTS_LOG.md")
if log_a is None or log_b is None:
    sys.exit("RESULTS_LOG.md missing at one of the two commits")
na = NUMBER.findall(log_a)
nb = NUMBER.findall(log_b)
# align the two number sequences, so an inserted paragraph does not report every later number as
# changed; only genuine replacements of equal length are counted as a moved number
log_changes = []
sm = difflib.SequenceMatcher(a=na, b=nb, autojunk=False)
for tag, i1, i2, j1, j2 in sm.get_opcodes():
    if tag == "replace" and (i2 - i1) == (j2 - j1):
        for k in range(i2 - i1):
            log_changes.append({"pre": na[i1 + k], "post": nb[j1 + k]})
    elif tag in ("replace", "delete", "insert"):
        log_changes.append({"tag": tag, "pre": na[i1:i2][:6], "post": nb[j1:j2][:6]})
print(f"   {len(nb)} numbers in the log at POST; {len(log_changes)} change(s) after alignment")
for c in log_changes[:10]:
    print(f"     {c}")

part_b = {
    "pre_commit": pre_sha, "post_commit": post_sha,
    "pre_ref": PRE, "post_ref": POST,
    "rel_tolerance": REL_TOL,
    "excluded_field_pattern": EXCLUDE.pattern,
    "n_artifacts": len(art_names),
    "n_numeric_fields_compared": n_fields_total,
    "n_fields_excluded": n_excluded,
    "n_fields_changed": len(changed_fields),
    "changed_fields_by_magnitude_band": {k: {"n": len(v), "fields": v} for k, v in bands.items()},
    "changed_fields": changed_fields,
    "per_artifact": per_artifact,
}
part_c = {
    "pre_commit": pre_sha, "post_commit": post_sha,
    "number_regex": NUMBER.pattern,
    "n_numbers_at_post": len(nb),
    "n_numbers_at_pre": len(na),
    "n_changes_after_alignment": len(log_changes),
    "changes": log_changes,
    "alignment": "difflib.SequenceMatcher over the extracted number sequences, so inserted or "
                 "deleted prose does not report every following number as changed",
}

report = {
    "weighting": "none — this step compares stored and re-simulated arrays and two git states; no "
                 "ensemble is weighted",
    "purpose": "close the evidence chain for the three whole-repository claims Step 15's own "
               "artifact does not cover",
    "config_sha256": provenance.config_sha256(),
    "full_library": part_a,
    "artifact_regression": part_b,
    "log_regression": part_c,
    "reproducibility": {
        "part_A": "live; re-runs the legacy arm against the stored corrected cache",
        "parts_B_and_C": "historical; read from git at the two recorded commits. They are NOT "
                         "re-derivable from the working tree, because later work legitimately "
                         "moved numbers afterwards. If the history is rewritten this step fails "
                         "rather than silently reporting nothing",
    },
}
with open(os.path.join(OUT, "step15_full_regression.json"), "w") as f:
    json.dump(report, f, indent=2)
print("\nsaved step15_full_regression.json")
