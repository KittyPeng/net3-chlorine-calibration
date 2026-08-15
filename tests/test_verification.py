"""The Section 2.2.4 verification checks, as assertions rather than as printed tables.

The three checks the paper describes each existed as a step script that computed a result and
printed it for a human to read. That is enough to write a paper from and not enough for a reader
to confirm, so they are restated here as tests that fail loudly.

Run with:

    pytest -q tests/

or, without pytest:

    python tests/test_verification.py

What is checked
  1. the frozen Net3 input file is the one every result was produced from, by SHA-256
  2. a single pipe with bulk decay only reproduces the analytic first-order solution
  3. the wall reaction has the expected sign, responds monotonically to its coefficient, and
     approaches the no-mass-transfer-resistance limit from the correct side
  4. concentrations round-trip through WNTR's internal kg/m^3 representation
  5. the below-threshold duration is the trapezoidal integral Section 2.8.1 defines
  6. the environment matches the one that produced the cached results

Checks 2 to 5 build their own inputs and need no cached artifact. Check 6 reads the manifest.
"""
import hashlib
import json
import os
import sys

import numpy as np
import pytest
import wntr

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "Net3"))

import wq_common as B  # noqa: E402  (needs the path above)

SECONDS_PER_DAY = 24 * 3600


# ---------------------------------------------------------------- 1. the frozen model
def test_net3_input_hash():
    """The .inp is the study's fixed point. If it moves, nothing downstream is comparable."""
    path = os.path.join(ROOT, "models", "net3_frozen", "Net3.inp")
    assert os.path.exists(path), "frozen Net3.inp is missing"
    digest = hashlib.sha256(open(path, "rb").read()).hexdigest()
    manifest = json.load(open(os.path.join(ROOT, "Net3", "baseline_cache",
                                           "cache_manifest.json")))
    assert digest == manifest["net3_inp_sha256"], (
        "Net3.inp does not match the hash recorded with the cached results.\n"
        "  file     %s\n  manifest %s" % (digest, manifest["net3_inp_sha256"]))


# ---------------------------------------------------------------- 2. analytic bulk decay
def _single_pipe(kb_per_day, kw_m_per_day=0.0, c0_mgl=1.0,
                 length_m=1000.0, diameter_m=0.3, velocity_ms=0.5, hours=12):
    """One reservoir, one pipe, one demand junction, run to steady state on quality."""
    wn = wntr.network.WaterNetworkModel()
    wn.add_reservoir("R", base_head=50.0)
    wn.add_junction("J", base_demand=velocity_ms * np.pi * (diameter_m / 2) ** 2,
                    elevation=0.0)
    wn.add_pipe("P", "R", "J", length=length_m, diameter=diameter_m, roughness=100.0)
    wn.options.quality.parameter = "CHEMICAL"
    wn.options.time.duration = hours * 3600
    wn.options.time.hydraulic_timestep = 3600
    wn.options.time.quality_timestep = 60
    wn.options.time.report_timestep = 3600
    wn.options.reaction.bulk_coeff = kb_per_day / SECONDS_PER_DAY
    wn.options.reaction.wall_coeff = kw_m_per_day / SECONDS_PER_DAY
    # NOT converted, unlike the concentrations on the next line. WNTR converts node quality
    # through its internal kg/m^3 representation but writes this option verbatim into the INP
    # file, where EPANET reads it in the file's own mg/L units; Section 2.2.4 and
    # step13_known_answer.py both check that. Dividing it by 1000 here ran these tests a
    # thousand times stricter than every production run, which is not the configuration the
    # paper describes.
    wn.options.quality.tolerance = B.QUALITY_TOLERANCE
    wn.get_node("R").initial_quality = c0_mgl / 1000.0            # mg/L -> kg/m^3
    res = wntr.sim.EpanetSimulator(wn).run_sim()
    series = res.node["quality"]["J"] * 1000.0                    # kg/m^3 -> mg/L
    return float(series.iloc[-1]), length_m / velocity_ms


@pytest.mark.parametrize("kb", [-0.5, -1.0, -2.0])
def test_bulk_decay_matches_the_analytic_solution(kb):
    """C = C0 exp(k_b t_res) once the front has passed, for bulk decay with no wall reaction."""
    c_end, t_res = _single_pipe(kb_per_day=kb)
    expected = 1.0 * np.exp(kb * t_res / SECONDS_PER_DAY)
    rel = abs(c_end - expected) / expected
    assert rel < 0.02, (
        "bulk decay at k_b = %s: simulated %.6f mg/L, analytic %.6f, relative error %.4f"
        % (kb, c_end, expected, rel))


# ---------------------------------------------------------------- 3. wall reaction
def test_wall_reaction_sign_and_monotonicity():
    """A more negative wall coefficient must remove more chlorine, and never add any."""
    baseline, _ = _single_pipe(kb_per_day=-0.5, kw_m_per_day=0.0)
    previous = baseline
    for kw in (-0.05, -0.1, -0.5, -1.0):
        c, _ = _single_pipe(kb_per_day=-0.5, kw_m_per_day=kw)
        assert c < baseline + 1e-9, "wall reaction at k_w = %s increased chlorine" % kw
        assert c <= previous + 1e-9, (
            "chlorine did not fall monotonically with k_w: %.6f at %s after %.6f"
            % (c, kw, previous))
        previous = c


def _effective_wall_rate(kw, kb=-0.5):
    """The wall contribution to the decay rate implied by the simulated concentration, 1/day.

    Comparing concentrations directly does not test saturation, because the effect is
    multiplicative: the same proportional change is a larger absolute change at a lower
    concentration. The rate is the quantity that saturates.
    """
    c, t_res = _single_pipe(kb_per_day=kb, kw_m_per_day=kw)
    return np.log(c) / (t_res / SECONDS_PER_DAY) - kb


def test_wall_reaction_is_mass_transfer_limited():
    """The wall reaction is limited by transport to the wall, so its effect saturates.

    Doubling a small coefficient nearly doubles the effective rate. Doubling an already large one
    barely moves it, because the reaction is no longer the slow step. Without mass-transfer
    resistance the response would stay proportional at every scale, which is the limiting case
    Section 2.2.4 checks against.
    """
    small = _effective_wall_rate(-0.10) / _effective_wall_rate(-0.05)
    large = _effective_wall_rate(-10.0) / _effective_wall_rate(-5.0)
    assert small > 1.7, "doubling a small k_w should nearly double the rate, got x%.2f" % small
    assert large < 1.3, "doubling a large k_w should barely move the rate, got x%.2f" % large
    assert large < small, "the response did not saturate: x%.2f then x%.2f" % (small, large)


# ---------------------------------------------------------------- 4. units
def test_concentration_round_trips_through_internal_units():
    """WNTR stores concentration in kg/m^3. A value set as mg/L must come back as mg/L.

    This is the error the study documents in Section 2.2.4: assigning 1.0 intending 1 mg/L makes
    EPANET simulate 1000 mg/L, and reading the output as mg/L repeats it on the way out.
    """
    for c0 in (0.5, 1.0, 2.0):
        c_end, t_res = _single_pipe(kb_per_day=-0.5, c0_mgl=c0)
        expected = c0 * np.exp(-0.5 * t_res / SECONDS_PER_DAY)
        assert abs(c_end - expected) / expected < 0.02, (
            "round trip at C0 = %s mg/L returned %.6f, expected about %.6f"
            % (c0, c_end, expected))


def test_first_order_kinetics_are_scale_invariant():
    """Doubling the source doubles every concentration, which is why the unit error was invisible."""
    c1, _ = _single_pipe(kb_per_day=-0.5, c0_mgl=1.0)
    c2, _ = _single_pipe(kb_per_day=-0.5, c0_mgl=2.0)
    assert abs(c2 / c1 - 2.0) < 0.01, "not linear in the source: ratio %.4f" % (c2 / c1)


# ---------------------------------------------------------------- 5. risk metric definition
def test_pbar_is_a_trapezoidal_integral_not_a_point_mean():
    """Section 2.8.1 integrates the below-threshold indicator over the T-1 one-hour intervals.

    Averaging the indicator over the T reporting points is a different quantity: it gives both
    endpoints full weight and divides by the point count. Step 10 integrated while the
    perturbation scripts averaged, so the baseline risk field and the arms it was compared
    against were not the same statistic. They agreed to three decimals on most junctions but
    ordered the top six differently, and the top six is what the risk conclusions rest on. The
    split is pinned here so that it cannot reopen unnoticed.
    """
    sys.path.insert(0, os.path.join(ROOT, "Net3"))
    import wq_common as B

    C = np.ones((1, 5, 1))
    C[0, :3, 0] = 0.0                      # below the threshold at the first three points
    got = float(B.pbar_field(C, 0.2, np.array([1.0]))[0])

    # trapezoid over the four intervals: 0.5 + 1 + 1 + 0 + 0 = 2.5 hours, over a 4 h window
    assert abs(got - 0.625) < 1e-12, "P_bar is %r, not the trapezoidal 0.625" % got
    assert abs(got - 3 / 5) > 1e-3, "P_bar has collapsed to the five-point mean 0.6"
    assert B.window_hours(C) == 4, "T reporting points bound T-1 one-hour intervals"


# ---------------------------------------------------------------- 6. environment
def test_environment_matches_the_cached_results():
    """The versions matter: this study documents version-specific WNTR behaviour."""
    manifest = json.load(open(os.path.join(ROOT, "Net3", "baseline_cache",
                                           "cache_manifest.json")))
    assert wntr.__version__ == manifest["wntr"], (
        "wntr %s, but the cached results were produced with %s"
        % (wntr.__version__, manifest["wntr"]))
    assert np.__version__ == manifest["numpy"], (
        "numpy %s, but the cached results were produced with %s"
        % (np.__version__, manifest["numpy"]))


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
