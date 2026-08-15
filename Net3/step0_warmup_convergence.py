"""Step 0: is the 24 h warm-up long enough? Decide it by a pre-declared convergence test.

This is the step that SET the current configuration, so it is written from the position before that
decision. The draft's configuration -- now superseded -- was a 72 h run discarding the first 24 h and
assessing 24-72 h, and that 24 h was never justified. It had to be: the tanks start at an assumed
0.5 mg/L, several high-risk junctions have mean water ages of tens of hours, and the Step 12 paired
test already showed the continuous severity metrics moving by 10-14% when the warm-up was extended.
The answer below is 120 h, which is why wq_common now carries WARMUP_H = 120 and DURATION_H = 168
and every other step assesses 120-168 h.

Method. Demands and the pump schedule in Net3 are 24 h periodic, so the correct notion of "warmed
up" is CYCLOSTATIONARY: the field over one diurnal cycle repeats in the next. Each parameter set is
run once over the full horizon and successive 24 h cycles are differenced. The recommended warm-up
is the start of the earliest cycle whose difference from the next cycle passes every criterion
below; the criteria are fixed here, before the numbers are seen.

HORIZON LIMIT. Pump 10 is driven by ABSOLUTE-time controls enumerated only to 159 h, and the model
duration is 168 h. Beyond 168 h the pump would stay closed for good, so a longer run would be a
different system rather than a longer warm-up. The test therefore has 7 cycles to work with and
cannot certify a warm-up beyond 144 h; if convergence is not reached by then, that is the finding.

Outputs: baseline_cache/step0_warmup_convergence.json, figures/step0_warmup_convergence.png
"""
import json
import os
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import wntr

import wq_common as B

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures")
os.makedirs(FIGDIR, exist_ok=True)

HORIZON_H = 168               # the model's ceiling; see the docstring
# The warm-up this test was posed against. It is HISTORY, not configuration: the draft used 24 h,
# this test rejected it, and wq_common.WARMUP_H was then raised to the answer. Reading B.WARMUP_H
# here instead would make the step describe its own conclusion as its premise — after the config was
# updated the artifact would say "the baseline warm-up is 120 h" and the verdict would compare 120
# against 120, erasing the question the step exists to answer.
DRAFT_WARMUP_H = 24
CYCLE_H = 24                  # demand and pump-schedule period
N_CYCLES = HORIZON_H // CYCLE_H
C_MIN = 0.2

# ---- acceptance criteria, declared before the results are seen ----
TOL = {
    "monitor_max_dC": 0.005,      # mg/L, max over the 6 monitors x 24 h
    "network_p95_dC": 0.010,      # mg/L, 95th percentile over junctions of the per-node max
    "tank_max_dC": 0.010,         # mg/L, max over the 3 tanks
    "age_p95_dAge": 1.0,          # h, 95th percentile over junctions of the per-node max
    "tank_max_dLevel": 0.05,      # m, hydraulic periodicity of the driver
    "risk_rel_dDeficit": 0.02,    # relative change in the network-mean cumulative deficit
}
TOP_K = 10                        # the top-K risk set must be identical between the two cycles

# Parameter sets: the synthetic truth plus the two corners of the prior box. The weak corner is the
# binding case — the least reactive network forgets its initial condition most slowly.
PARAM_SETS = {
    "truth": (B.KW_OLD_TRUE, B.KW_AVG_TRUE, B.KW_NEW_TRUE),
    "prior_weak": (B.PRIOR["old"][1], B.PRIOR["avg"][1], B.PRIOR["new"][1]),
    "prior_strong": (B.PRIOR["old"][0], B.PRIOR["avg"][0], B.PRIOR["new"][0]),
}

ALL_NODES = B.all_junctions()
wn_ref = wntr.network.WaterNetworkModel(B.NET3_INP)
TANKS = list(wn_ref.tank_name_list)
TANK_ELEV = {t: wn_ref.get_node(t).elevation for t in TANKS}
MON_POS = [ALL_NODES.index(m) for m in B.MONITOR_NODES]


def run_long(kw, quality="CHEMICAL"):
    """One HORIZON_H run. Returns (junction field, tank field, tank levels) as (T, n) arrays."""
    # For AGE the initial "quality" is an age, so the chlorine initial conditions are meaningless
    # there; zero them so sources and tanks start at age 0, matching the convention in Step 10.
    ic = {} if quality == "CHEMICAL" else {"inlet_mgl": 0.0, "tank_mgl": 0.0}
    wn = B.build_model(B.KB_FIXED, 0.0, duration_hours=HORIZON_H, quality=quality,
                       pre_run=B.make_kw_hook(*kw), **ic)
    res, hours = B.run_model(wn)
    q = res.node["quality"]
    jun = q[ALL_NODES].values
    tank = q[TANKS].values
    head = res.node["head"][TANKS].values
    level = head - np.array([TANK_ELEV[t] for t in TANKS])[None, :]
    if quality == "AGE":
        jun, tank = jun / 3600.0, tank / 3600.0      # seconds -> hours
    return jun, tank, level, hours


def cycle(arr, k):
    return arr[k * CYCLE_H:(k + 1) * CYCLE_H]


def per_node_max(a, b):
    return np.abs(a - b).max(axis=0)


def risk_of_cycle(field, k):
    """Single-member risk severity over one diurnal cycle: hours below C_MIN and deficit.

    Uses CYCLE_H + 1 points, not CYCLE_H. Trapezoidal integration over 24 points spans 23 intervals,
    so the previous slice measured a 23 h cycle and under-reported both integrals by about 1/24. The
    cycle-to-cycle RATIO the convergence criterion uses was almost unaffected, because both cycles
    carried the same bias, but the absolute per-cycle deficit was not.
    """
    c = field[k * CYCLE_H:(k + 1) * CYCLE_H + 1]
    below = (c < C_MIN).astype(float)
    deficit = np.clip(C_MIN - c, 0, None)
    return np.trapezoid(below, dx=1.0, axis=0), np.trapezoid(deficit, dx=1.0, axis=0)


# ---------------------------------------------------------------- run
t0 = time.time()
print(f"=== Step 0: warm-up convergence over {HORIZON_H} h "
      f"({N_CYCLES} x {CYCLE_H} h cycles) ===")
print(f"horizon capped at {HORIZON_H} h: pump 10 has absolute-time controls only to 159 h\n")

age_jun, age_tank, _, hours = run_long(PARAM_SETS["truth"], quality="AGE")
fields = {}
for name, kw in PARAM_SETS.items():
    fields[name] = run_long(kw)
    print(f"  ran {name:>12}  kw = {kw}")
print(f"  ran water age (reaction-independent)\n  {time.time() - t0:.1f}s\n")

rows = []
for k in range(N_CYCLES - 1):
    rec = {"cycle": k, "window_h": [k * CYCLE_H, (k + 1) * CYCLE_H],
           "compared_with_h": [(k + 1) * CYCLE_H, (k + 2) * CYCLE_H], "sets": {}}
    for name, (jun, tank, level, _) in fields.items():
        dj = per_node_max(cycle(jun, k), cycle(jun, k + 1))
        dt = per_node_max(cycle(tank, k), cycle(tank, k + 1))
        dl = per_node_max(cycle(level, k), cycle(level, k + 1))
        dur_a, def_a = risk_of_cycle(jun, k)
        dur_b, def_b = risk_of_cycle(jun, k + 1)
        top_a = {ALL_NODES[i] for i in np.argsort(def_a)[::-1][:TOP_K]}
        top_b = {ALL_NODES[i] for i in np.argsort(def_b)[::-1][:TOP_K]}
        rel_def = (abs(def_b.mean() - def_a.mean()) / def_a.mean()) if def_a.mean() > 0 else 0.0
        rec["sets"][name] = {
            "monitor_max_dC": float(dj[MON_POS].max()),
            "network_p95_dC": float(np.percentile(dj, 95)),
            "network_max_dC": float(dj.max()),
            "tank_max_dC": float(dt.max()),
            "tank_max_dLevel": float(dl.max()),
            "risk_rel_dDeficit": float(rel_def),
            "net_mean_deficit": float(def_a.mean()),
            "net_mean_hours_below": float(dur_a.mean()),
            f"top{TOP_K}_identical": top_a == top_b,
            f"top{TOP_K}_jaccard": len(top_a & top_b) / len(top_a | top_b),
        }
    da = per_node_max(cycle(age_jun, k), cycle(age_jun, k + 1))
    rec["age_p95_dAge"] = float(np.percentile(da, 95))
    rec["age_max_dAge"] = float(da.max())
    rows.append(rec)


def passes(rec):
    """Every criterion, for every parameter set, on this cycle pair."""
    fails = []
    if rec["age_p95_dAge"] > TOL["age_p95_dAge"]:
        fails.append(f"age_p95_dAge {rec['age_p95_dAge']:.2f} > {TOL['age_p95_dAge']}")
    for name, s in rec["sets"].items():
        for key, tol in TOL.items():
            if key == "age_p95_dAge":
                continue
            if s[key] > tol:
                fails.append(f"{name}/{key} {s[key]:.4g} > {tol}")
        if not s[f"top{TOP_K}_identical"]:
            fails.append(f"{name}/top{TOP_K} set changed")
    return fails


verdicts = [(rec["cycle"], passes(rec)) for rec in rows]
first_ok = next((k for k, f in verdicts if not f), None)
recommended = None if first_ok is None else first_ok * CYCLE_H


def worst_series(key):
    """Per cycle pair, the worst value of a criterion across the parameter sets."""
    if key == "age_p95_dAge":
        return [rec[key] for rec in rows]
    return [max(s[key] for s in rec["sets"].values()) for rec in rows]


def extrapolate_cycle(series, tol, n_fit=3):
    """Cycle at which a geometrically decaying series would first reach tol.

    Reported only for criteria that never pass inside the horizon, and labelled as extrapolated:
    the model cannot be run far enough to verify it.
    """
    y = np.asarray(series[-n_fit:], dtype=float)
    if len(y) < 2 or np.any(y <= 0):
        return None
    ratios = y[1:] / y[:-1]
    r = float(np.exp(np.mean(np.log(ratios))))
    if r >= 1.0:                       # not decaying: no finite crossing
        return None
    k_last = rows[-1]["cycle"]
    extra = np.log(tol / y[-1]) / np.log(r)
    return k_last + float(np.ceil(extra))


per_criterion = {}
for key, tol in TOL.items():
    series = worst_series(key)
    ok = [rec["cycle"] for rec, v in zip(rows, series) if v <= tol]
    first = ok[0] if ok else None
    entry = {"tolerance": tol, "worst_by_cycle": [round(v, 6) for v in series],
             "first_pass_cycle": first,
             "first_pass_warmup_h": None if first is None else first * CYCLE_H}
    if first is None:
        k_ext = extrapolate_cycle(series, tol)
        entry["extrapolated_pass_cycle"] = k_ext
        entry["extrapolated_warmup_h"] = None if k_ext is None else int(k_ext * CYCLE_H)
        entry["extrapolation_note"] = ("geometric fit to the last three cycle pairs; NOT verified, "
                                       "the model horizon cannot reach this cycle")
    per_criterion[key] = entry

# ---------------------------------------------------------------- report
print(f"{'cycle':>5} {'window':>10} | {'monDC':>7} {'netP95':>7} {'netMax':>7} {'tankDC':>7} "
      f"{'tankDL':>7} {'ageP95':>7} {'relDef':>7} | verdict")
for rec, (k, fails) in zip(rows, verdicts):
    w = f"{rec['window_h'][0]}-{rec['window_h'][1]}"
    worst = {key: max(s[key] for s in rec["sets"].values())
             for key in ("monitor_max_dC", "network_p95_dC", "network_max_dC", "tank_max_dC",
                         "tank_max_dLevel", "risk_rel_dDeficit")}
    print(f"{k:>5} {w:>10} | {worst['monitor_max_dC']:7.4f} {worst['network_p95_dC']:7.4f} "
          f"{worst['network_max_dC']:7.4f} {worst['tank_max_dC']:7.4f} "
          f"{worst['tank_max_dLevel']:7.4f} {rec['age_p95_dAge']:7.3f} "
          f"{worst['risk_rel_dDeficit']:7.4f} | "
          + ("PASS" if not fails else f"fail: {'; '.join(fails[:2])}"
             + (f" (+{len(fails) - 2} more)" if len(fails) > 2 else "")))

print("\nper-criterion verdict (worst over the three parameter sets):")
print(f"{'criterion':>20} {'tol':>8} {'last value':>11} | earliest warm-up that satisfies it")
for key, e in per_criterion.items():
    last = e["worst_by_cycle"][-1]
    if e["first_pass_warmup_h"] is not None:
        verdict = f"{e['first_pass_warmup_h']} h  (verified)"
    elif e.get("extrapolated_warmup_h") is not None:
        verdict = (f"not within {HORIZON_H} h; extrapolates to ~{e['extrapolated_warmup_h']} h "
                   f"(UNVERIFIED)")
    else:
        verdict = f"not within {HORIZON_H} h and not decaying geometrically"
    print(f"{key:>20} {e['tolerance']:>8} {last:>11.4f} | {verdict}")

if recommended is None:
    print(f"\nVERDICT: no single cycle within the {HORIZON_H} h horizon satisfies EVERY criterion.")
    print("The concentration field does settle, but the integrated risk severity and the water-age")
    print("field are still drifting at the end of the longest run the model permits, so the")
    print("residual drift has to be carried as a stated limitation rather than assumed away.")
else:
    print(f"\nVERDICT: cyclostationary from hour {recommended} -> recommended warm-up "
          f"{recommended} h (the draft used {DRAFT_WARMUP_H} h; the project now runs {B.WARMUP_H} h).")
    if recommended > DRAFT_WARMUP_H:
        print(f"The draft warm-up was too short by {recommended - DRAFT_WARMUP_H} h, which is why "
              f"wq_common now carries WARMUP_H = {B.WARMUP_H}.")

# ---- what the calibration configuration should therefore be ----
CHEM_KEYS = ["monitor_max_dC", "network_p95_dC", "tank_max_dC"]
chem_h = [per_criterion[k]["first_pass_warmup_h"] for k in CHEM_KEYS]
proposed_warmup = max(h for h in chem_h if h is not None) if all(h is not None for h in chem_h) \
    else None
window_h = B.DURATION_H - B.WARMUP_H
proposed_duration = None if proposed_warmup is None else proposed_warmup + window_h
config = {
    "proposed_warmup_h": proposed_warmup,
    "proposed_duration_h": proposed_duration,
    "assessment_window_h": window_h,
    "n_resid_unchanged": proposed_warmup is not None and
                         len(B.MONITOR_NODES) * (window_h + 1) == B.N_RESID,
    "saturates_model_horizon": proposed_duration == HORIZON_H,
    "basis": "the largest warm-up that the chlorine criteria verify; the risk-severity and "
             "water-age criteria cannot be satisfied inside the model horizon",
    "residual_risk_drift_at_proposed": per_criterion["risk_rel_dDeficit"]["worst_by_cycle"][-1],
    "age_is_horizon_dependent": True,
}
print("\nconfiguration implied for the calibration pipeline:")
print(f"  warm-up {proposed_warmup} h + {window_h} h assessment window = "
      f"{proposed_duration} h total"
      + ("  (exactly the model horizon — no room to spare)" if config["saturates_model_horizon"]
         else ""))
print(f"  residual count unchanged at {B.N_RESID} "
      f"({len(B.MONITOR_NODES)} monitors x {window_h + 1} hours), so the behavioural threshold "
      f"{B.RMSE_THR} carries over")
print(f"  residual risk-severity drift at that warm-up: "
      f"{config['residual_risk_drift_at_proposed'] * 100:.1f}% per cycle — state it, do not hide it")
print("  water age stays horizon-dependent: report it as a diagnostic, not an equilibrium property")

report = {
    "horizon_h": HORIZON_H, "cycle_h": CYCLE_H, "n_cycles": N_CYCLES,
    "horizon_limit_reason": "pump 10 uses absolute-time controls enumerated only to 159 h and the "
                            "model duration is 168 h; beyond that the pump stays closed, which "
                            "would be a different system rather than a longer warm-up",
    "criterion_note": "demands and the pump schedule are 24 h periodic, so convergence means "
                      "cyclostationarity: cycle k must match cycle k+1 within every tolerance",
    "tolerances": TOL, "top_k": TOP_K,
    "param_sets": {k: list(v) for k, v in PARAM_SETS.items()},
    "kb_fixed": B.KB_FIXED, "C_MIN": C_MIN,
    # kept distinct on purpose: the value under test, and the value the project now runs with
    "draft_warmup_h": DRAFT_WARMUP_H,
    "current_warmup_h": B.WARMUP_H,
    "recommended_warmup_h": recommended,
    "implied_config": config,
    "per_criterion": per_criterion,
    "cycles": rows,
    "failures_by_cycle": {str(k): f for k, f in verdicts},
    "runtime_s": round(time.time() - t0, 1),
}
with open(os.path.join(HERE, "baseline_cache", "step0_warmup_convergence.json"), "w") as f:
    json.dump(report, f, indent=2)

# ---------------------------------------------------------------- figure
fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.4))
starts = [rec["window_h"][0] for rec in rows]
colors = {"truth": "tab:blue", "prior_weak": "tab:red", "prior_strong": "tab:green"}

ax = axes[0]
for name, c in colors.items():
    ax.semilogy(starts, [rec["sets"][name]["monitor_max_dC"] for rec in rows], "o-", color=c,
                label=f"{name} monitors")
    ax.semilogy(starts, [rec["sets"][name]["network_p95_dC"] for rec in rows], "s--", color=c,
                alpha=0.45, label=f"{name} network p95")
ax.axhline(TOL["monitor_max_dC"], color="k", ls=":", lw=1,
           label=f"monitor tol {TOL['monitor_max_dC']}")
ax.axvline(DRAFT_WARMUP_H, color="crimson", lw=1.5, label=f"draft warm-up {DRAFT_WARMUP_H} h")
ax.axvline(B.WARMUP_H, color="seagreen", lw=1.5, ls="--", label=f"adopted {B.WARMUP_H} h")
ax.set_xlabel("cycle start (h) — compared with the next 24 h cycle")
ax.set_ylabel("max |ΔC| between successive cycles (mg/L)")
ax.set_title("(a) chlorine cyclostationarity", fontsize=10)
ax.legend(fontsize=6, ncol=2)
ax.grid(alpha=0.3, which="both")

ax = axes[1]
ax.semilogy(starts, [rec["age_p95_dAge"] for rec in rows], "o-", color="tab:purple",
            label="water age, network p95")
ax.semilogy(starts, [max(r["sets"][n]["tank_max_dLevel"] for n in colors) for r in rows],
            "s-", color="tab:brown", label="tank level, max (m)")
ax.axhline(TOL["age_p95_dAge"], color="tab:purple", ls=":", lw=1, label="age tol 1.0 h")
ax.axhline(TOL["tank_max_dLevel"], color="tab:brown", ls=":", lw=1, label="level tol 0.05 m")
ax.axvline(DRAFT_WARMUP_H, color="crimson", lw=1.5)
ax.axvline(B.WARMUP_H, color="seagreen", lw=1.5, ls="--")
ax.set_xlabel("cycle start (h)")
ax.set_ylabel("difference between successive cycles")
ax.set_title("(b) water age and the hydraulic driver", fontsize=10)
ax.legend(fontsize=7)
ax.grid(alpha=0.3, which="both")

ax = axes[2]
# cycle 0 is the start-up transient (deficit ~1 mg/L.h) and would flatten everything after it
for name, c in colors.items():
    ax.plot(starts[1:], [rec["sets"][name]["net_mean_deficit"] for rec in rows[1:]], "o-",
            color=c, label=name)
ax.axvline(DRAFT_WARMUP_H, color="crimson", lw=1.5, label=f"draft warm-up {DRAFT_WARMUP_H} h")
ax.axvline(B.WARMUP_H, color="seagreen", lw=1.5, ls="--", label=f"adopted {B.WARMUP_H} h")
if proposed_warmup is not None:
    ax.axvline(proposed_warmup, color="seagreen", lw=1.5, ls="--",
               label=f"proposed {proposed_warmup} h")
ax.set_xlabel("cycle start (h)")
ax.set_ylabel("network-mean cumulative deficit (mg/L·h)")
ax.set_title("(c) the quantity that actually moves\n(start-up cycle omitted)", fontsize=10)
ax.legend(fontsize=7)
ax.grid(alpha=0.3)

fig.suptitle(f"Step 0 — warm-up convergence: chlorine, water age, tank level and risk severity "
             f"over successive 24 h cycles ({HORIZON_H} h horizon)", y=1.03)
plt.tight_layout()
figpath = os.path.join(FIGDIR, "step0_warmup_convergence.png")
plt.savefig(figpath, dpi=130, bbox_inches="tight")
print("\nfigure saved to", figpath)
print("saved step0_warmup_convergence.json")
