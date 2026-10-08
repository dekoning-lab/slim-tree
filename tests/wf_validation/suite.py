"""Validation suite for SLiM-Tree: simulation vs weak-mutation theory, for Wright-Fisher models (default) or
SLiM-Tree's non-WF model (--model nonWF; rates are then per tick, see theory.py).

Checks
  1. neutral_mode     SLiM-Tree -N (no fitness callback), default recombination.
                      Substitution rate per nucleotide site = mu (Kimura 1968).
  2. flat_profiles    One uniform stationary distribution (all amino acids equally fit; the setup used
                      for "neutral" runs with fitness profiles). Stops are the only selected states.
                      Theory: sense-codon substitution rate, uniform stationary codon distribution,
                      no premature stop fixes.
  3. site_profiles    Three site-heterogeneous profiles (strong / graded / weak preferences), free
                      recombination so there is no background selection. For each profile:
                      synonymous and nonsynonymous substitution rates, dN/dS, stationary
                      amino-acid and codon distributions, no premature stop fixes.

Theory uses the fitness values SLiM-Tree actually wrote into the script (so these checks test the
simulation, not the fitting); a separate check compares the fitted values with the intended ones.

All theory assumes weak mutation; run with theta = 4 N mu <= 0.01 (0.001 is safer, and slower).
"""

import argparse, csv, json, math, os, shutil, sys, time
from concurrent.futures import ProcessPoolExecutor
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import theory, harness  # noqa: E402

N_DEFAULT = 100
CODONS_DEFAULT = 300
Z_MAX = 3.5          # substitution-rate checks: |z| must be below this
P_MIN = 0.005        # stationary-distribution checks: bootstrap p-value must be at least this
                     # (needs --bootstrap > 1/P_MIN - 1, otherwise a check could never fail)


def validation_profiles(N):
    """Three amino-acid fitness profiles, specified on the scaled scale S = 4 N (f - f_max)."""
    aa = theory.AA_ORDER
    strong = np.full(20, -10.0)
    for a, s in zip("LIVMF", [0.0, -1.0, -2.0, -3.0, -4.0]):
        strong[aa.index(a)] = s
    graded = -np.linspace(0, 8, 20)[np.random.default_rng(2024).permutation(20)]
    weak = -np.linspace(0, 1.5, 20)[np.random.default_rng(7).permutation(20)]
    named = {"strong": strong, "graded": graded, "weak": weak}
    return {k: 1.0 + v / (4.0 * N) for k, v in named.items()}


def write_stationary_csv(path, profiles_f, N):
    with open(path, "w", newline="") as fh:
        csv.writer(fh).writerows(theory.stationary_csv_rows(list(profiles_f), N))


# ---------------------------------------------------------------- bootstrap null (CTMC on fixed states)

def ctmc_sample(Q, pi, n_sites, gens, interval, rng):
    """Simulate the weak-mutation fixed-state CTMC for n_sites independent sites started at stationarity,
    observed at times 0, interval, ..., gens. Returns occupancy counts (61, SENSE order) and the numbers
    of synonymous and nonsynonymous substitutions in (0, gens]."""
    rate = -np.diag(Q)
    P = Q / rate[:, None]
    np.fill_diagonal(P, 0.0)
    cumP = np.cumsum(P, axis=1)
    aa = theory.AA_OF_CODON[np.array(theory.SENSE)]
    n_obs = gens // interval + 1
    state = rng.choice(61, size=n_sites, p=pi)
    t = np.zeros(n_sites)
    occ = np.zeros(61, dtype=np.int64)
    syn = non = 0
    alive = np.arange(n_sites)
    while alive.size:
        s = state[alive]
        t0 = t[alive]
        t1 = t0 + rng.exponential(1.0 / rate[s])
        first = np.ceil(t0 / interval)
        last = np.minimum(np.ceil(t1 / interval) - 1, n_obs - 1)
        np.add.at(occ, s, np.clip(last - first + 1, 0, None).astype(np.int64))
        go = t1 <= gens
        if go.any():
            u = rng.random(go.sum())
            new = (cumP[s[go]] > u[:, None]).argmax(axis=1)
            same = aa[s[go]] == aa[new]
            syn += int(same.sum())
            non += int((~same).sum())
            state[alive[go]] = new
            t[alive[go]] = t1[go]
        alive = alive[go]
    return occ, syn, non


def tvd(p, q):
    return 0.5 * float(np.abs(np.asarray(p) - np.asarray(q)).sum())


def to_aa(codon_counts_sense):
    out = np.zeros(20)
    for k, c in enumerate(theory.SENSE):
        out[theory.AA_OF_CODON[c]] += codon_counts_sense[k]
    return out


# ---------------------------------------------------------------- running

def run_config(name, jobs_spec, outdir, workers):
    """Run all replicates of one configuration in parallel; returns list of replicate results."""
    os.makedirs(outdir, exist_ok=True)
    jobs = []
    for spec in jobs_spec:
        spec = dict(spec)
        spec["workdir"] = os.path.join(outdir, f"{name}_rep{spec['seed']}")
        if os.path.exists(spec["workdir"]):
            shutil.rmtree(spec["workdir"])
        jobs.append(spec)
    with ProcessPoolExecutor(max_workers=workers) as ex:
        return list(ex.map(harness.replicate, jobs))


def load_results(name, outdir):
    """Reload the per-replicate results of an earlier run (for --evaluate-only)."""
    reps = sorted(d for d in os.listdir(outdir) if d.startswith(name + "_rep"))
    out = []
    for d in reps:
        with open(os.path.join(outdir, d, "result.json")) as fh:
            out.append(json.load(fh))
    if not out:
        raise FileNotFoundError(f"no saved {name} replicates in {outdir}")
    return out


def make_jobs(base, reps, seed0):
    return [dict(base, seed=seed0 + r) for r in range(reps)]


# ---------------------------------------------------------------- evaluation

def rate_check(label, counts, expected, note=""):
    """counts: per-replicate observed counts; expected: per-replicate expected counts."""
    counts, expected = np.asarray(counts, float), np.asarray(expected, float)
    obs, exp = counts.sum(), expected.sum()
    ratio = obs / exp
    per_rep = counts / expected
    se_rep = per_rep.std(ddof=1) / math.sqrt(len(per_rep)) if len(per_rep) > 1 else float("nan")
    se_pois = math.sqrt(exp) / exp
    se = max(se_pois, se_rep if not math.isnan(se_rep) else 0.0)
    z = (ratio - 1.0) / se
    return {"check": label, "observed": obs, "expected": exp, "ratio": ratio, "se": se, "z": z,
            "pass": abs(z) < Z_MAX, "note": note}


def zero_check(label, count, note=""):
    return {"check": label, "observed": count, "expected": 0, "pass": count == 0, "note": note}


def stationary_checks(label, occ_obs_sense, Q, pi, n_sites, gens, interval, B, rng, pi_neutral):
    """Amino-acid and codon-level checks of the stationary distribution against a parametric-bootstrap null
    (the theoretical fixed-state chain simulated with the same sites, duration and sampling)."""
    obs = np.asarray(occ_obs_sense, float)
    p_obs = obs / obs.sum()
    null = {"aa": [], "codon": []}
    for _ in range(B):
        occ, _, _ = ctmc_sample(Q, pi, n_sites, gens, interval, rng)
        p = occ / occ.sum()
        null["aa"].append(tvd(to_aa(p), to_aa(pi)))
        null["codon"].append(tvd(p, pi))
    out = []
    for level, name in (("aa", "amino-acid"), ("codon", "codon")):
        target = to_aa(pi) if level == "aa" else pi
        ref = to_aa(pi_neutral) if level == "aa" else pi_neutral
        d_obs = tvd(to_aa(p_obs) if level == "aa" else p_obs, target)
        nl = np.array(null[level])
        pval = (1 + (nl >= d_obs).sum()) / (B + 1)
        out.append({"check": f"{label}: stationary {name} distribution", "observed_tvd": d_obs,
                    "null_tvd_median": float(np.median(nl)), "null_tvd_max": float(nl.max()),
                    "null_z": float((d_obs - nl.mean()) / nl.std(ddof=1)), "p": float(pval),
                    "neutral_vs_theory_tvd": tvd(ref, target), "pass": bool(pval >= P_MIN)})
    return out


def evaluate_profiles(name, results, N, mu, intended, B, rng, checks, model="WF"):
    """Per-profile substitution and stationary checks for a fitness-profile configuration."""
    info0 = results[0]["info"]
    fit = info0["fitness_profiles"]
    scale = info0.get("selection_scale", 1.0)
    b = info0.get("offspring_per_tick", 1)
    n_prof = len(next(iter(fit.values())))
    neutral_pi = theory.expected(np.ones(20), N, mu, model, scale, b)["pi_codon"]
    names = list(intended) if intended else [f"profile{p}" for p in range(n_prof - 1)]
    for p in range(n_prof - 1):                         # last column = SLiM-Tree's neutral profile
        f = np.array([fit[a][p] for a in theory.AA_ORDER])
        if intended is not None:
            want = intended[names[p]]
            checks.append({"check": f"{name}/{names[p]}: fitted fitness vs intended (max abs diff)",
                           "observed": float(np.abs((f - f.max()) - (want - want.max())).max()),
                           "expected": 0, "pass": bool(np.abs((f - f.max()) - (want - want.max())).max() < 1e-6),
                           "note": "SLiM-Tree's fitting step; theory below uses the fitted values"})
        e = theory.expected(f, N, mu, model, scale, b)
        if model == "nonWF":
            target = theory.expected(f, N, mu)["pi_codon"]
            checks.append({"check": f"{name}/{names[p]}: non-WF theory vs Wright-Fisher target (codon TVD)",
                           "observed": tvd(e["pi_codon"], target), "expected": 0, "pass": None,
                           "note": f"selection scale {scale}; how closely the non-WF model reproduces the stationary "
                                   "distribution the profile was fitted to"})
        key = str(p)
        reps = [r for r in results if key in r["profiles"]]
        sites = [r["profiles"][key]["codon_sites"] for r in reps]
        gens = [r["generations"] for r in reps]
        syn = [r["profiles"][key]["subs"]["syn"] for r in reps]
        non = [r["profiles"][key]["subs"]["nonsyn"] for r in reps]
        stop = sum(r["profiles"][key]["subs"]["to_stop"] for r in reps)
        exp_syn = [e["K_syn"] * s * g for s, g in zip(sites, gens)]
        exp_non = [e["K_nonsyn"] * s * g for s, g in zip(sites, gens)]
        lab = f"{name}/{names[p]}"
        checks.append(rate_check(f"{lab}: synonymous substitutions", syn, exp_syn))
        checks.append(rate_check(f"{lab}: nonsynonymous substitutions", non, exp_non))
        obs_dnds = (sum(non) / sum(exp_non)) / (sum(syn) / sum(exp_syn)) * e["dNdS"]
        checks.append({"check": f"{lab}: dN/dS", "observed": obs_dnds, "expected": e["dNdS"],
                       "ratio": obs_dnds / e["dNdS"], "pass": None,
                       "note": "summary of the two rate checks above (not scored separately)"})
        checks.append(zero_check(f"{lab}: premature stop codons fixed", stop,
                                 "stop codons are excluded from SLiM-Tree's stationary distributions"))
        occ = np.sum([np.array(r["profiles"][key]["occupancy"]) for r in reps], axis=0)
        occ_sense = occ[np.array(theory.SENSE)]
        interval = reps[0]["info"]["sample_interval"]
        n_sites = int(sum(sites))
        g = int(np.median(gens))
        checks.extend(stationary_checks(lab, occ_sense, e["Q"], e["pi_codon"], n_sites, g, interval, B, rng, neutral_pi))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--theta", type=float, default=0.01, help="4 N mu (<= 0.01; weak-mutation theory)")
    ap.add_argument("--N", type=int, default=N_DEFAULT)
    ap.add_argument("--codons", type=int, default=CODONS_DEFAULT)
    ap.add_argument("--measure-gens", type=int, default=100_000)
    ap.add_argument("--equil-gens", type=int, default=None, help="default 10 N")
    ap.add_argument("--interval", type=int, default=100)
    ap.add_argument("--reps-neutral", type=int, default=10)
    ap.add_argument("--reps-flat", type=int, default=10)
    ap.add_argument("--reps-sites", type=int, default=20)
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--bootstrap", type=int, default=1000)
    ap.add_argument("--evaluate-only", action="store_true", help="re-score saved replicates in --outdir without simulating")
    ap.add_argument("--seed-offset", type=int, default=0, help="added to every replicate seed (default runs reuse seeds 1000+, 2000+, 3000+)")
    ap.add_argument("--only", nargs="*", default=None, help="subset of: neutral_mode flat_profiles site_profiles")
    ap.add_argument("--slimtree-args", default="", help="extra slim-tree arguments for every non-neutral replicate")
    ap.add_argument("--model", choices=("WF", "nonWF"), default="WF", help="Wright-Fisher (default) or SLiM-Tree's non-WF model")
    a = ap.parse_args(argv)
    if a.theta > 0.01 + 1e-12:
        sys.exit("theta must be <= 0.01: the theory assumes weak mutation")
    if 1.0 / (a.bootstrap + 1) >= P_MIN:
        sys.exit(f"--bootstrap must exceed {1 / P_MIN - 1:.0f} so that stationary checks can fail")
    N, mu = a.N, a.theta / (4 * a.N)
    equil = a.equil_gens if a.equil_gens is not None else 10 * N
    os.makedirs(a.outdir, exist_ok=True)
    rng = np.random.default_rng(12345)
    base = dict(N=N, mu=mu, codons=a.codons, measure_gens=a.measure_gens, equil_gens=equil,
                sample_interval=a.interval, extra_args=tuple(a.slimtree_args.split()), nonWF=a.model == "nonWF")
    only = set(a.only) if a.only else {"neutral_mode", "flat_profiles", "site_profiles"}
    checks, timing, all_results = [], {}, []

    if "neutral_mode" in only:
        t = time.time()
        res = load_results("neutral_mode", a.outdir) if a.evaluate_only else run_config("neutral_mode", make_jobs(dict(base, neutral=True, recombination=2.5e-8),
                                                   a.reps_neutral, 1000 + a.seed_offset), a.outdir, a.jobs)
        timing["neutral_mode"] = time.time() - t
        all_results += res
        counts = [sum(sum(p["subs"].values()) for p in r["profiles"].values()) for r in res]
        sites = [3 * sum(p["codon_sites"] for p in r["profiles"].values()) for r in res]
        b = res[0]["info"].get("offspring_per_tick", 1)
        expv = [theory.neutral_rate_per_site(mu, a.model, b) * s * r["generations"] for s, r in zip(sites, res)]
        label = "mu" if a.model == "WF" else f"mu b/(1 + b) per tick, b = {b}"
        checks.append(rate_check(f"neutral_mode: substitutions per nucleotide site = {label}", counts, expv,
                                 "SLiM-Tree -N; default recombination (linked); Kimura 1968"))

    if "flat_profiles" in only:
        flat_csv = os.path.join(a.outdir, "flat_stationary.csv")
        write_stationary_csv(flat_csv, [np.ones(20)], N)
        t = time.time()
        res = load_results("flat_profiles", a.outdir) if a.evaluate_only else run_config("flat_profiles", make_jobs(dict(base, stationary_csv=flat_csv, recombination=0.5),
                                                    a.reps_flat, 2000 + a.seed_offset), a.outdir, a.jobs)
        timing["flat_profiles"] = time.time() - t
        all_results += res
        evaluate_profiles("flat_profiles", res, N, mu, {"flat": np.ones(20)}, a.bootstrap, rng, checks, a.model)

    if "site_profiles" in only:
        prof = validation_profiles(N)
        site_csv = os.path.join(a.outdir, "site_stationary.csv")
        write_stationary_csv(site_csv, prof.values(), N)
        t = time.time()
        res = load_results("site_profiles", a.outdir) if a.evaluate_only else run_config("site_profiles", make_jobs(dict(base, stationary_csv=site_csv, recombination=0.5),
                                                    a.reps_sites, 3000 + a.seed_offset), a.outdir, a.jobs)
        timing["site_profiles"] = time.time() - t
        all_results += res
        evaluate_profiles("site_profiles", res, N, mu, prof, a.bootstrap, rng, checks, a.model)

    versions = sorted({r.get("slim_version", "unknown") for r in all_results})
    report = {"settings": {**vars(a), "mu": mu, "equil_gens": equil, "slim_version": ", ".join(versions)},
              "timing_seconds": timing, "checks": checks}
    with open(os.path.join(a.outdir, "results.json"), "w") as fh:
        json.dump(report, fh, indent=1, default=float)
    with open(os.path.join(a.outdir, "report.md"), "w") as fh:
        fh.write(render(report))
    print(render(report))
    scored = [c for c in checks if c.get("pass") is not None]
    return 0 if all(c["pass"] for c in scored) else 1


def render(report):
    s = report["settings"]
    nonwf = s.get("model", "WF") == "nonWF"
    unit = "ticks" if nonwf else "generations"
    lines = [f"# SLiM-Tree {'non-WF' if nonwf else 'Wright-Fisher'} validation\n",
             f"N = {s['N']}, theta = {s['theta']} (mu = {s['mu']:.3g}), {s['codons']} codons, "
             f"{s['measure_gens']:,} {unit} measured after {s['equil_gens']:,}, sampled every {s['interval']}.\n",
             f"Extra slim-tree arguments: `{s.get('slimtree_args') or 'none'}`. SLiM: {s.get('slim_version', 'unknown')}.\n",
             "| Check | Observed | Expected | Ratio / p | Result |", "|---|---|---|---|---|"]
    for c in report["checks"]:
        res = "—" if c.get("pass") is None else ("PASS" if c["pass"] else "**FAIL**")
        if "observed_tvd" in c:
            lines.append(f"| {c['check']} | TVD {c['observed_tvd']:.4f} | null median {c['null_tvd_median']:.4f}, max {c['null_tvd_max']:.4f} "
                         f"(fully neutral would be {c['neutral_vs_theory_tvd']:.3f}) | p = {c['p']:.3g} (z = {c['null_z']:.0f}) | {res} |")
        elif "ratio" in c and "z" in c:
            lines.append(f"| {c['check']} | {c['observed']:.0f} | {c['expected']:.1f} | "
                         f"{c['ratio']:.3f} ± {c['se']:.3f} (z = {c['z']:.1f}) | {res} |")
        elif "ratio" in c:
            lines.append(f"| {c['check']} | {c['observed']:.3f} | {c['expected']:.3f} | {c['ratio']:.3f} | {res} |")
        else:
            lines.append(f"| {c['check']} | {c['observed']:.3g} | {c['expected']} | | {res} |")
    lines.append("")
    lines.append("Timing (s): " + ", ".join(f"{k} {v:.0f}" for k, v in report["timing_seconds"].items()))
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    sys.exit(main())
