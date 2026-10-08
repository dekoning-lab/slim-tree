"""Exact checks of SLiM-Tree's fitness function on hand-built genotypes, for the site and gene stop models.

Each test generates a Wright-Fisher script with SLiM-Tree's real pipeline (harness.generate) and appends a probe event at
tick 2. The probe clears one haplosome, adds chosen mutations, evaluates get_genome_fitness and prints the result. The
expected values are computed here independently, from the fitness profiles SLiM-Tree passed to the writer:
  * interior coding sites use their profile; a stop's value must be the site's lowest amino-acid value minus
    S_stop / (4N - 2) (site model) or equal to it (gene model);
  * each gene's start codon has ATG optimal and its terminal stop codon has every stop optimal, and any other state there
    costs S_stop / (4N - 2) (site model);
  * non-coding sites are neutral, stops included;
  * a haplosome's selection factor is exp(sum over changed codons of f_new - f_ancestral).

Needs SLiM 5: set WFV_SLIM to the binary, or put SLiM 5 on PATH. The tests are skipped otherwise.
"""

import os, re, subprocess, sys
import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import harness, suite, theory  # noqa: E402

N = 50
PENALTY = 100.0
DELTA = PENALTY / (4 * N - 2)
NUC = "ACGT"


def slim_binary():
    slim = os.environ.get("WFV_SLIM", "slim")
    try:
        out = subprocess.run([slim, "-version"], capture_output=True, text=True).stdout
    except OSError:
        out = ""
    if not out.startswith("SLiM version 5"):
        pytest.skip("SLiM 5 not available (set WFV_SLIM)")
    return slim


def generate(tmp_path, name, extra=(), codons=40):
    csv = str(tmp_path / f"{name}_site.csv")
    suite.write_stationary_csv(csv, suite.validation_profiles(N).values(), N)
    # a placeholder observer; the probe is appended below
    script, info = harness.generate(str(tmp_path / name), N=N, mu=1e-12, codons=codons, recombination=0.0, measure_gens=N,
                                    equil_gens=N, sample_interval=N, seed=17, stationary_csv=csv, extra_args=tuple(extra),
                                    observer="")
    return script, info


def run_probe(script, cases, extra_lines=""):
    """cases: {name: [(nucleotide position, nucleotide), ...]} applied in order to a cleared haplosome."""
    lines = ["2 late() {", "\th = p1.haplosomes[0];"]
    for name, muts in cases.items():
        lines.append("\th.removeMutations();")
        for pos, nuc in muts:
            lines.append(f'\th.addNewMutation(m1, 0.0, {pos}, NULL, "{nuc}");')
        lines.append(f'\tcatn("CASE {name} " + paste(format("%.17g", get_genome_fitness(h)), sep=" "));')
    lines.append(extra_lines)
    lines.append("\tsim.simulationFinished();\n}\n")
    with open(script, "a") as fh:
        fh.write("\n".join(lines))
    res = subprocess.run([slim_binary(), script], cwd=os.path.dirname(os.path.dirname(script)), capture_output=True, text=True)
    assert res.returncode == 0, res.stderr[-2000:]
    out = {}
    for line in res.stdout.splitlines():
        if line.startswith("CASE "):
            parts = line.split()
            out[parts[1]] = [float(x) for x in parts[2:]]
        elif line.startswith("VALUE "):
            parts = line.split()
            out[parts[1]] = [float(x) for x in parts[2:]]
    return out


class Model:
    """Independent re-statement of the site values SLiM-Tree is meant to use."""

    def __init__(self, info, stop_model="site"):
        self.prof = info["fitness_profiles"]
        self.nums = info["profile_nums"]
        self.anc = [theory.CODONS[c] for c in info["ancestral_codons"]]
        self.starts = {s for s, e in info["coding_seqs"]}
        self.ends = {e for s, e in info["coding_seqs"]}
        self.coding = set()
        for s, e in info["coding_seqs"]:
            self.coding.update(range(s, e + 1))
        self.stop_model = stop_model

    @staticmethod
    def aa(codon):
        a = theory.TRANSLATE[codon]
        return "X" if a == "*" else a

    def f(self, site, codon):
        a = self.aa(codon)
        if self.stop_model == "site" and site in self.starts:
            return 1.0 if a == "M" else 1.0 - DELTA
        if self.stop_model == "site" and site in self.ends:
            return 1.0 if a == "X" else 1.0 - DELTA
        return self.prof[a][self.nums[site]]

    def log_factor(self, changes):
        """changes: {site: new codon string}"""
        return sum(self.f(s, c) - self.f(s, self.anc[s]) for s, c in changes.items())

    def interior_sites(self):
        return [s for s in sorted(self.coding) if s not in self.starts and s not in self.ends]

    def stop_neighbour(self, site):
        for j, _ in theory.neighbours(theory.CODONS.index(self.anc[site])):
            if theory.TRANSLATE[theory.CODONS[j]] == "*":
                return theory.CODONS[j]
        return None


def mutations(site, old, new):
    return [(3 * site + i, new[i]) for i in range(3) if old[i] != new[i]]


def test_site_model_stop_values_and_genotypes(tmp_path):
    script, info = generate(tmp_path, "site")
    m = Model(info)
    # stop values: lowest amino-acid value minus S_stop / (4N - 2) for every real profile; neutral profile stays 1
    aas = [a for a in m.prof if a != "X"]
    n_prof = len(m.prof["A"])
    for p in range(n_prof - 1):
        assert m.prof["X"][p] == pytest.approx(min(m.prof[a][p] for a in aas) - DELTA, abs=1e-12)
    assert m.prof["X"][n_prof - 1] == 1.0

    interior = m.interior_sites()
    stop_sites = [s for s in interior if m.stop_neighbour(s)][:2]
    s1, s2 = stop_sites
    start, end = min(m.starts), max(m.ends)
    term = m.anc[end]
    other_stop = next(theory.CODONS[j] for j, _ in theory.neighbours(theory.CODONS.index(term))
                      if theory.TRANSLATE[theory.CODONS[j]] == "*")
    read_through = next(theory.CODONS[j] for j, _ in theory.neighbours(theory.CODONS.index(term))
                        if theory.TRANSLATE[theory.CODONS[j]] != "*")
    start_loss = "ATA"
    missense = next(theory.CODONS[j] for j, _ in theory.neighbours(theory.CODONS.index(m.anc[interior[5]]))
                    if theory.TRANSLATE[theory.CODONS[j]] not in ("*", theory.TRANSLATE[m.anc[interior[5]]]))
    stop1, stop2 = m.stop_neighbour(s1), m.stop_neighbour(s2)
    cases = {
        "stop": mutations(s1, m.anc[s1], stop1),
        "two_stops": mutations(s1, m.anc[s1], stop1) + mutations(s2, m.anc[s2], stop2),
        "reversion": mutations(s1, m.anc[s1], stop1) + mutations(s1, stop1, m.anc[s1]),
        "stop_to_stop": mutations(end, term, other_stop),
        "read_through": mutations(end, term, read_through),
        "start_loss": mutations(start, "ATG", start_loss),
        "missense": mutations(interior[5], m.anc[interior[5]], missense),
    }
    expected = {
        "stop": m.log_factor({s1: stop1}),
        "two_stops": m.log_factor({s1: stop1, s2: stop2}),
        "reversion": 0.0,
        "stop_to_stop": 0.0,
        "read_through": -DELTA,
        "start_loss": -DELTA,
        "missense": m.log_factor({interior[5]: missense}),
    }
    got = run_probe(script, cases)
    for name, logf in expected.items():
        assert got[name][0] == pytest.approx(np.exp(logf), rel=1e-12), name
    # a stop is worse than the amino acid it replaces by at least the penalty
    assert -expected["stop"] >= DELTA - 1e-12
    assert expected["two_stops"] == pytest.approx(expected["stop"] + m.log_factor({s2: stop2}))


def test_site_model_profile_shift_invariance(tmp_path):
    # adding a constant to every value of a profile must not change any fitness ratio
    script, info = generate(tmp_path, "shift")
    m = Model(info)
    s1 = [s for s in m.interior_sites() if m.stop_neighbour(s)][0]
    stop1 = m.stop_neighbour(s1)
    cases = {"stop": mutations(s1, m.anc[s1], stop1)}
    text = open(script).read()
    shifted = os.path.join(os.path.dirname(script), "shifted.slim")
    n_prof = int(re.search(r'defineConstant\("N_PROF", (\d+)\)', text).group(1))
    table = re.search(r'defineConstant\("FIT_TABLE", c\(([^)]*)\)\)', text)
    values = [float(v) for v in table.group(1).split(",")]
    p = m.nums[s1]
    values = [v + 0.37 if i % n_prof == p else v for i, v in enumerate(values)]
    with open(shifted, "w") as fh:
        fh.write(text[:table.start(1)] + ", ".join(repr(v) for v in values) + text[table.end(1):])
    base = run_probe(script, cases)["stop"][0]
    moved = run_probe(shifted, cases)["stop"][0]
    assert moved == pytest.approx(base, rel=1e-12)


def test_site_model_two_genes_and_noncoding(tmp_path):
    script, info = generate(tmp_path, "genes", extra=("-G", "2", "-C", "0.8"), codons=80)
    m = Model(info)
    second_start = sorted(m.starts)[1]
    noncoding = [s for s in range(len(m.anc)) if s not in m.coding]
    nc = next(s for s in noncoding if any(theory.TRANSLATE[theory.CODONS[j]] == "*"
                                          for j, _ in theory.neighbours(theory.CODONS.index(m.anc[s]))))
    nc_stop = m.stop_neighbour(nc)
    cases = {"second_start_loss": mutations(second_start, "ATG", "ATA"),
             "noncoding_stop": mutations(nc, m.anc[nc], nc_stop)}
    got = run_probe(script, cases)
    assert got["second_start_loss"][0] == pytest.approx(np.exp(-DELTA), rel=1e-12)
    assert got["noncoding_stop"][0] == 1.0


def test_gene_model_null_alleles(tmp_path):
    script, info = generate(tmp_path, "gene", extra=("--stop_model", "gene", "--null_cost", "0.3", "--null_dominance", "0.2"))
    m = Model(info, stop_model="gene")
    aas = [a for a in m.prof if a != "X"]
    for p in range(len(m.prof["A"]) - 1):
        assert m.prof["X"][p] == pytest.approx(min(m.prof[a][p] for a in aas), abs=1e-12)
    interior = m.interior_sites()
    s1 = [s for s in interior if m.stop_neighbour(s)][0]
    stop1 = m.stop_neighbour(s1)
    # a heavily loaded copy: the most deleterious single-step missense change at each of 15 other sites
    load = {}
    for s in [x for x in interior if x != s1][:15]:
        options = [theory.CODONS[j] for j, _ in theory.neighbours(theory.CODONS.index(m.anc[s]))
                   if theory.TRANSLATE[theory.CODONS[j]] != "*"]
        load[s] = min(options, key=lambda c: m.f(s, c))
    load_muts = [x for s, c in load.items() for x in mutations(s, m.anc[s], c)]
    cases = {"stop": mutations(s1, m.anc[s1], stop1),
             "loaded": load_muts,
             "loaded_knockout": load_muts + mutations(s1, m.anc[s1], stop1)}
    got = run_probe(script, cases, extra_lines='\tcatn("VALUE null_factors " + paste(format("%.17g", NULL_FACTORS), sep=" "));')
    assert got["null_factors"] == pytest.approx([0.7, 0.94, 1.0], rel=1e-12)
    # premature stop: the gene is null, and the stop still has its own (non-beneficial) site effect
    assert got["stop"][1] == 0.0
    assert got["stop"][0] == pytest.approx(np.exp(m.log_factor({s1: stop1})), rel=1e-12)
    assert got["stop"][0] <= 1.0
    # no masking: the knocked-out copy keeps its load, so knocking out a loaded copy never raises fitness
    assert got["loaded"][1] == 1.0 and got["loaded_knockout"][1] == 0.0
    assert got["loaded_knockout"][0] == pytest.approx(got["loaded"][0] * np.exp(m.log_factor({s1: stop1})), rel=1e-9)
    w_functional = got["loaded"][0] * got["null_factors"][2]
    w_knockout = got["loaded_knockout"][0] * got["null_factors"][1]
    assert w_knockout < w_functional
