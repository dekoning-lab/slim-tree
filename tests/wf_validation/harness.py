"""Run SLiM-Tree's own script generator, add a read-only observer, run SLiM, and parse what it saw.

The script is produced by SLiM-Tree's real pipeline (readInput -> findFitness/calculateFitnesses ->
cladeReader -> writeSLiM), exactly as `slim-tree <tree> <stationary.csv> ...` would, except that SLiM
is not launched automatically. The only change to the generated script is an appended late() event
that records, every `sample_interval` generations, one genome's nucleotides and the positions of all
segregating mutations. It reads state and writes a file; it does not touch the population.

Measurements are taken in the root population p1, which runs for SLiM-Tree's burn-in
(burn_in_multiplier * N generations). The tree is "(A:1);", so the tip population only adds one
generation at the end.
"""

import hashlib, json, os, random, shutil, subprocess, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SLIM_TREE_ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
if SLIM_TREE_ROOT not in sys.path:
    sys.path.insert(0, SLIM_TREE_ROOT)
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import theory  # noqa: E402

OBSERVER = """
{t0}:{t1} late() {{
	if ((community.tick - {t0}) % {k} != 0) return;
	muts = sim.mutations;
	poly = integer(0);
	if (size(muts) > 0) {{
		f = sim.mutationFrequencies(p1, muts);
		segregating = muts[f < 1.0];
		if (size(segregating) > 0) poly = sort(unique(segregating.position));
	}}
	writeFile("{out}", paste0(community.tick, "\\t", p1.haplosomes[0].nucleotides(), "\\t", paste(poly, sep=",")), append=T);
}}
"""


def generate(workdir, *, N, mu, codons, recombination, measure_gens, equil_gens, sample_interval,
             seed, neutral=False, stationary_csv=None, extra_args=(), observer=OBSERVER):
    """Generate a SLiM-Tree Wright-Fisher script in `workdir` and append the observer.
    Returns (script_path, info) where info records what SLiM-Tree chose (profiles, ancestor, ...)."""
    workdir = os.path.abspath(workdir)
    if stationary_csv is not None:
        stationary_csv = os.path.abspath(stationary_csv)
    os.makedirs(workdir, exist_ok=False)
    total = equil_gens + measure_gens
    if total % N:
        raise ValueError("equil_gens + measure_gens must be a multiple of N (SLiM-Tree burn-in = multiplier * N)")
    with open(os.path.join(workdir, "tree.txt"), "w") as fh:
        fh.write("(A:1);\n")
    args = ["tree.txt", "null" if neutral else stationary_csv,
            "-n", str(N), "-v", repr(mu), "-g", str(codons), "-r", repr(recombination),
            "-b", str(total // N)] + (["-N"] if neutral else []) + list(extra_args)
    # WFV_SLIMTREE_ARGS adds slim-tree arguments to every run, except to neutral-mode runs
    if not neutral:
        args += os.environ.get("WFV_SLIMTREE_ARGS", "").split()

    random.seed(seed)
    np.random.seed(seed % (2 ** 32))
    old_cwd, old_argv = os.getcwd(), sys.argv
    os.chdir(workdir)
    sys.argv = ["slim-tree"] + args
    try:
        import SLiMTree
        import utils
        st = SLiMTree.SLiMTree.__new__(SLiMTree.SLiMTree)       # skip __init__, which would launch SLiM
        params = st.read_input()
        if params["nonWF"]:
            raise RuntimeError("this suite is for Wright-Fisher models only")
        redo = True
        while redo:
            params = st.process_fitness(params)
            clades = utils.cladeReader.cladeReader(params)
            clade_list = clades.get_clade_dict_list()
            redo = clades.redo_random
        st.write_slim_code(clade_list, params)
    finally:
        os.chdir(old_cwd)
        sys.argv = old_argv

    script = params["filenames"][0] + "_p1.slim"
    samples = os.path.join(workdir, "observer_samples.tsv")
    with open(script) as fh:
        text = fh.read()
    if "nonWF" in text:
        raise RuntimeError("generated script is not Wright-Fisher")
    with open(script, "a") as fh:
        fh.write(observer.format(t0=equil_gens, t1=total, k=sample_interval, out=samples))

    profiles = None
    if not neutral:
        # aa -> list of fitness values per profile (last column is SLiM-Tree's 'neutral' profile)
        profiles = {aa: [float(x) for x in vals] for aa, vals in params["fitness_profiles"].items()}
    info = {
        "N": N, "mu": mu, "codons": codons, "recombination": recombination, "neutral": neutral,
        "equil_gens": equil_gens, "measure_gens": measure_gens, "sample_interval": sample_interval,
        "seed": seed, "args": args, "script": script, "samples": samples,
        "ancestral_codons": [int(c) for c in params["ancestral_sequence"]],
        "profile_nums": [int(p) for p in params["fitness_profile_nums"]],
        "fitness_profiles": profiles,
    }
    with open(os.path.join(workdir, "run_info.json"), "w") as fh:
        json.dump(info, fh, indent=1)
    return script, info


def run_slim(script, workdir, timeout=None):
    t = time.time()
    # WFV_SLIM selects the SLiM binary (default: `slim` on PATH); its version is recorded with each run
    slim = os.environ.get("WFV_SLIM", "slim")
    check = subprocess.run([slim, "-version"], capture_output=True, text=True)
    if check.returncode != 0 or not check.stdout.strip():
        raise RuntimeError(f"SLiM binary {slim!r} could not run (`-version` failed): {check.stderr.strip()[-500:]}")
    version = check.stdout.splitlines()[0]
    binary = shutil.which(slim) or slim
    with open(binary, "rb") as fh:
        digest = hashlib.sha256(fh.read()).hexdigest()
    with open(os.path.join(workdir, "slim_version.txt"), "w") as fh:
        fh.write(binary + "\n" + version + " (sha256 " + digest[:16] + ")\n")
    res = subprocess.run([slim, script], cwd=workdir, capture_output=True, text=True, timeout=timeout)
    with open(os.path.join(workdir, "slim_stdout.txt"), "w") as fh:
        fh.write(res.stdout)
    with open(os.path.join(workdir, "slim_stderr.txt"), "w") as fh:
        fh.write(res.stderr)
    if res.returncode != 0:
        raise RuntimeError(f"SLiM failed ({res.returncode}): {res.stderr[-2000:]}")
    return time.time() - t


_NUC = {c: i for i, c in enumerate("ACGT")}
_IS_STOP = np.array([theory.TRANSLATE[c] == "*" for c in theory.CODONS])
_AA = np.array([theory.TRANSLATE[c] for c in theory.CODONS])


def parse(info):
    """Turn observer samples into substitution counts and fixed-state occupancy.

    A site's fixed state is updated whenever it is not segregating; a change in fixed state is one
    substitution. Codon context for classifying a substitution comes from the current fixed states of
    the other two positions. Codon 0 (start) and the last codon (stop) are excluded. Returns per-profile
    tallies keyed by SLiM-Tree's profile number."""
    L = info["codons"]
    anc = np.array([[_NUC[ch] for ch in theory.CODONS[c]] for c in info["ancestral_codons"]]).ravel()
    fixed = anc.copy()
    prof = np.array(info["profile_nums"])
    internal = np.zeros(L, dtype=bool)
    internal[1:L - 1] = True
    out = {}
    for p in np.unique(prof[internal]):
        out[int(p)] = {"codon_sites": int(((prof == p) & internal).sum()),
                       "subs": {"syn": 0, "nonsyn": 0, "to_stop": 0, "from_stop": 0},
                       "occupancy": np.zeros(64, dtype=np.int64)}
    n_samples = 0
    first_tick = last_tick = None
    with open(info["samples"]) as fh:
        for line in fh:
            tick, seq, poly = line.rstrip("\n").split("\t")
            tick = int(tick)
            obs = np.frombuffer(seq.encode(), dtype=np.uint8)
            obs = np.select([obs == 65, obs == 67, obs == 71, obs == 84], [0, 1, 2, 3], -1)
            mask = np.ones(3 * L, dtype=bool)
            if poly:
                mask[np.array(poly.split(","), dtype=int)] = False
            changed = np.where(mask & (obs != fixed))[0]
            if first_tick is not None:          # first sample only re-synchronises the fixed state
                for site in changed:
                    c = site // 3
                    if not internal[c]:
                        continue
                    old = fixed[3 * c:3 * c + 3].copy()
                    new = old.copy()
                    new[site % 3] = obs[site]
                    oi, ni = 16 * old[0] + 4 * old[1] + old[2], 16 * new[0] + 4 * new[1] + new[2]
                    rec = out[int(prof[c])]["subs"]
                    if _IS_STOP[ni]:
                        rec["to_stop"] += 1
                    elif _IS_STOP[oi]:
                        rec["from_stop"] += 1
                    elif _AA[oi] == _AA[ni]:
                        rec["syn"] += 1
                    else:
                        rec["nonsyn"] += 1
                    fixed[site] = obs[site]
            fixed[mask] = obs[mask]
            if first_tick is None:
                first_tick = tick
            last_tick = tick
            n_samples += 1
            codon_idx = 16 * fixed[0::3] + 4 * fixed[1::3] + fixed[2::3]
            for p, rec in out.items():
                sel = (prof == p) & internal
                rec["occupancy"] += np.bincount(codon_idx[sel], minlength=64)
    for rec in out.values():
        rec["occupancy"] = rec["occupancy"].tolist()
    return {"profiles": {str(k): v for k, v in out.items()},
            "generations": last_tick - first_tick, "n_samples": n_samples}


def replicate(job):
    """Worker entry point: generate, run and parse one replicate. `job` is a dict of generate() kwargs
    plus 'workdir'."""
    job = dict(job)
    workdir = os.path.abspath(job.pop("workdir"))
    script, info = generate(workdir, **job)
    secs = run_slim(script, workdir)
    res = parse(info)
    res["info"] = {k: v for k, v in info.items() if k not in ("ancestral_codons",)}
    res["slim_seconds"] = secs
    with open(os.path.join(workdir, "slim_version.txt")) as fh:
        res["slim_version"] = fh.read().splitlines()[1]
    with open(os.path.join(workdir, "result.json"), "w") as fh:
        json.dump(res, fh)
    return res
