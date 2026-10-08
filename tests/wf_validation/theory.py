"""Weak-mutation theory for Wright-Fisher codon models, written independently of SLiM-Tree.

Model (the one SLiM-Tree's fitness profiles are calibrated against in utils/calculateFitnesses.py):
  * diploid Wright-Fisher population of N individuals (2N genomes)
  * Jukes-Cantor mutation: every specific nucleotide change occurs at rate mu/3, so the total
    mutation rate per nucleotide site per generation is mu (SLiM-Tree's -v; SLiM's mmJukesCantor(mu/3))
  * a new mutation changing amino-acid fitness from f_i to f_j has per-copy selection coefficient
    s = f_j - f_i (genic selection) and fixes with Kimura's probability
        u(s) = (1 - exp(-2s)) / (1 - exp(-4Ns)),   u(0) = 1/(2N)
  * stop codons are lethal: mutations to them never fix, so the chain lives on the 61 sense codons
  * weak mutation (theta = 4N mu << 1): each site is a continuous-time Markov chain on fixed states
    with rates q_ij = 2N * (mu/3) * u(f_j - f_i) between single-nucleotide neighbours.

Codons are indexed as in SLiM: nucleotides A,C,G,T = 0..3 and codon = 16*n1 + 4*n2 + n3.
Amino-acid fitness vectors are in alphabetical one-letter order, as in SLiM-Tree: ACDEFGHIKLMNPQRSTVWY.
"""

import numpy as np

NUCS = "ACGT"
AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
CODONS = [a + b + c for a in NUCS for b in NUCS for c in NUCS]          # SLiM codon index order
_CODE = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"
TRANSLATE = dict(zip(CODONS, _CODE))                                      # standard genetic code
SENSE = [i for i, c in enumerate(CODONS) if TRANSLATE[c] != "*"]          # 61 sense codons (SLiM indices)
STOP = [i for i, c in enumerate(CODONS) if TRANSLATE[c] == "*"]
AA_OF_CODON = np.array([AA_ORDER.index(TRANSLATE[c]) if TRANSLATE[c] != "*" else -1 for c in CODONS])


def neighbours(codon_idx):
    """Single-nucleotide neighbours of a codon: list of (neighbour index, codon position changed)."""
    c = CODONS[codon_idx]
    out = []
    for pos in range(3):
        for n in NUCS:
            if n != c[pos]:
                out.append((CODONS.index(c[:pos] + n + c[pos + 1:]), pos))
    return out


def pfix(s, N):
    """Kimura's fixation probability of a new mutation (initial frequency 1/(2N)) with per-copy effect s
    in a diploid Wright-Fisher population of N individuals. Numerically stable for small |s|."""
    s = np.asarray(s, dtype=float)
    out = np.empty_like(s)
    small = np.abs(4 * N * s) < 1e-8
    out[small] = 1.0 / (2 * N)
    ss = s[~small]
    out[~small] = np.expm1(-2 * ss) / np.expm1(-4 * N * ss)
    return out


def rate_matrix(aa_fitness, N, mu):
    """61x61 substitution rate matrix (per generation) over sense codons, ordered as SENSE.
    Also returns the matching mutation-only matrix (rates mu/3 between sense neighbours) and a
    synonymous-change mask."""
    f = np.asarray(aa_fitness, dtype=float)
    idx = {c: k for k, c in enumerate(SENSE)}
    Q = np.zeros((61, 61))
    M = np.zeros((61, 61))
    syn = np.zeros((61, 61), dtype=bool)
    for a, ci in enumerate(SENSE):
        fi = f[AA_OF_CODON[ci]]
        for cj, _ in neighbours(ci):
            if cj in STOP:
                continue                       # lethal: never fixes
            b = idx[cj]
            fj = f[AA_OF_CODON[cj]]
            M[a, b] = mu / 3
            Q[a, b] = 2 * N * (mu / 3) * pfix(np.array([fj - fi]), N)[0]
            syn[a, b] = AA_OF_CODON[ci] == AA_OF_CODON[cj]
        Q[a, a] = -Q[a].sum()
    return Q, M, syn


def stationary(Q):
    """Stationary distribution of a CTMC rate matrix (solves pi Q = 0, sum pi = 1)."""
    n = Q.shape[0]
    A = np.vstack([Q.T, np.ones(n)])
    b = np.zeros(n + 1)
    b[-1] = 1.0
    pi = np.linalg.lstsq(A, b, rcond=None)[0]
    pi[pi < 0] = 0.0
    return pi / pi.sum()


def sella_hirsh_stationary(aa_fitness, N):
    """Closed form for symmetric mutation with Kimura's u(s): detailed balance gives
    pi_i / pi_j = u(f_i - f_j) / u(f_j - f_i) = exp((4N - 2)(f_i - f_j)), so pi_c ∝ exp((4N-2) f_aa(c)).
    Used only to cross-check rate_matrix/stationary."""
    f = np.asarray(aa_fitness, dtype=float)
    w = np.array([np.exp((4 * N - 2) * (f[AA_OF_CODON[c]] - f.max())) for c in SENSE])
    return w / w.sum()


def expected(aa_fitness, N, mu):
    """Expected long-run quantities for one fitness profile at stationarity.

    Returns a dict with (per codon site per generation):
      K_total, K_syn, K_nonsyn      substitution rates
      dNdS                          Spielman & Wilke (2015) MutSel omega: (K_N / mutational flux_N) /
                                    (K_S / mutational flux_S), fluxes taken at stationarity
      pi_codon (length 61, SENSE order), pi_aa (length 20, AA_ORDER)
    """
    Q, M, syn = rate_matrix(aa_fitness, N, mu)
    pi = stationary(Q)
    off = ~np.eye(61, dtype=bool)
    K_syn = float(pi @ np.where(syn & off, Q, 0.0).sum(axis=1))
    K_non = float(pi @ np.where(~syn & off, Q, 0.0).sum(axis=1))
    flux_syn = float(pi @ np.where(syn, M, 0.0).sum(axis=1))
    flux_non = float(pi @ np.where(~syn, M, 0.0).sum(axis=1))
    pi_aa = np.zeros(20)
    for k, c in enumerate(SENSE):
        pi_aa[AA_OF_CODON[c]] += pi[k]
    return {
        "K_total": K_syn + K_non, "K_syn": K_syn, "K_nonsyn": K_non,
        "dNdS": (K_non / flux_non) / (K_syn / flux_syn),
        "pi_codon": pi, "pi_aa": pi_aa, "Q": Q,
    }


def neutral_rate_per_site(mu):
    """Kimura (1968): the neutral substitution rate per site equals the mutation rate."""
    return mu


def stationary_csv_rows(profiles_aa_fitness, N):
    """Rows for a SLiM-Tree stationary-distribution CSV (no header; column 0 = codon; one column per
    profile), holding the exact stationary codon distribution implied by each fitness profile.
    The mutation rate cancels out of the stationary distribution, so mu = 1."""
    pis = [expected(f, N, 1.0)["pi_codon"] for f in profiles_aa_fitness]
    rows = []
    for k, c in enumerate(SENSE):
        rows.append([CODONS[c]] + [f"{p[k]:.15g}" for p in pis])
    return rows
