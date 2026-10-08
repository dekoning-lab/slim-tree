"""Fast checks (no SLiM) of the weak-mutation theory used by the Wright-Fisher validation suite,
and of its agreement with SLiM-Tree's own calibration code (utils/calculateFitnesses.py)."""

import os, sys
import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
import theory  # noqa: E402
from suite import validation_profiles  # noqa: E402

N, MU = 100, 2.5e-5


def test_genetic_code():
    assert len(theory.SENSE) == 61 and len(theory.STOP) == 3
    assert theory.TRANSLATE["ATG"] == "M" and theory.TRANSLATE["TGG"] == "W"
    assert {theory.CODONS[i] for i in theory.STOP} == {"TAA", "TAG", "TGA"}


def test_pfix_limits():
    assert theory.pfix(np.array([0.0]), N)[0] == pytest.approx(1 / (2 * N))
    assert theory.pfix(np.array([1e-12]), N)[0] == pytest.approx(1 / (2 * N), rel=1e-6)
    s = 0.05                                    # strongly beneficial: u ~ 1 - exp(-2s)
    assert theory.pfix(np.array([s]), N)[0] == pytest.approx(-np.expm1(-2 * s), rel=1e-6)
    assert theory.pfix(np.array([-0.05]), N)[0] < 1e-8


def test_neutral_rates():
    e = theory.expected(np.ones(20), N, MU)
    # no fitness differences: synonymous and nonsynonymous changes fix at the mutation rate
    assert e["dNdS"] == pytest.approx(1.0)
    # stops are lethal, so the sense-codon rate is mu times the fraction of mutations that stay sense
    n_sense_nbrs = sum(sum(1 for j, _ in theory.neighbours(c) if j not in theory.STOP) for c in theory.SENSE)
    assert e["K_total"] == pytest.approx(MU / 3 * n_sense_nbrs / 61)
    assert np.allclose(e["pi_codon"], 1 / 61)


@pytest.mark.parametrize("name", ["strong", "graded", "weak"])
def test_stationary_matches_closed_form(name):
    f = validation_profiles(N)[name]
    e = theory.expected(f, N, MU)
    assert np.abs(e["pi_codon"] - theory.sella_hirsh_stationary(f, N)).max() < 1e-10


@pytest.mark.parametrize("name", ["strong", "graded", "weak"])
def test_agrees_with_slimtree_calibration(name):
    """SLiM-Tree's Compute_matrix/Calc_stationary_KLD define the stationary distribution its fitness
    profiles are fitted to; ours must be the same distribution."""
    from utils import calculateFitnesses as CF
    f = validation_profiles(N)[name]
    theirs = CF.Calc_stationary_KLD(CF.Compute_matrix(list(f), N))
    ours = theory.expected(f, N, MU)["pi_codon"]
    by_name = {theory.CODONS[c]: ours[k] for k, c in enumerate(theory.SENSE)}
    diffs = [abs(by_name[CF.Get_nucleotide_from_idx(i)] - theirs[i]) for i in range(61)]
    assert max(diffs) < 1e-10


@pytest.mark.parametrize("pop", [100, 1000, 10000])
@pytest.mark.parametrize("name", ["strong", "graded", "weak"])
def test_slimtree_fitting_recovers_profile(name, pop):
    """Feeding SLiM-Tree's fitter the exact stationary distribution recovers the fitness profile
    (up to the arbitrary maximum), to within 0.001 on the population-scaled scale 4 N f, at every
    population size we simulate (the old numerical optimizer was off by up to ~10 at N = 10,000)."""
    from utils import calculateFitnesses as CF
    f = validation_profiles(pop)[name]
    pi = theory.expected(f, pop, MU)["pi_codon"]
    by_name = {theory.CODONS[c]: pi[k] for k, c in enumerate(theory.SENSE)}
    fit = np.array(CF.Find_fits([by_name[CF.Get_nucleotide_from_idx(i)] for i in range(61)], pop))
    assert 4 * pop * np.abs((fit - fit.max()) - (f - f.max())).max() < 1e-3


def test_ctmc_bootstrap_unbiased():
    """The bootstrap null simulator reproduces the stationary distribution and expected jump counts."""
    from suite import ctmc_sample, to_aa
    f = validation_profiles(N)["graded"]
    e = theory.expected(f, N, MU)
    rng = np.random.default_rng(3)
    occ, syn, non = ctmc_sample(e["Q"], e["pi_codon"], 4000, 100_000, 100, rng)
    assert 0.5 * np.abs(to_aa(occ / occ.sum()) - e["pi_aa"]).sum() < 0.03
    exp_syn, exp_non = e["K_syn"] * 4000 * 100_000, e["K_nonsyn"] * 4000 * 100_000
    assert abs(syn - exp_syn) < 4 * np.sqrt(exp_syn) + 0.02 * exp_syn
    assert abs(non - exp_non) < 4 * np.sqrt(exp_non) + 0.02 * exp_non


def test_stationary_check_passes_null_and_fails_neutral():
    """The stationary check must be able to fail: with occupancy drawn from the theoretical chain it passes,
    and with occupancy drawn from the neutral chain (what a simulation with ineffective selection produces)
    it fails. Guards against the earlier bug where too few bootstrap draws made failure impossible."""
    from suite import ctmc_sample, stationary_checks
    f = validation_profiles(N)["graded"]
    e = theory.expected(f, N, MU)
    neutral = theory.expected(np.ones(20), N, MU)
    rng = np.random.default_rng(5)
    good, _, _ = ctmc_sample(e["Q"], e["pi_codon"], 600, 50_000, 100, rng)
    bad, _, _ = ctmc_sample(neutral["Q"], neutral["pi_codon"], 600, 50_000, 100, rng)
    ok = stationary_checks("t", good, e["Q"], e["pi_codon"], 600, 50_000, 100, 300, rng, neutral["pi_codon"])
    ko = stationary_checks("t", bad, e["Q"], e["pi_codon"], 600, 50_000, 100, 300, rng, neutral["pi_codon"])
    assert all(c["pass"] for c in ok)
    assert not any(c["pass"] for c in ko)
