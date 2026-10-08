# Wright–Fisher and non-WF validation suite

This suite checks that SLiM-Tree's Wright–Fisher simulations, and with `--model nonWF` its non-WF simulations, reproduce weak-mutation population-genetic theory. It compares long simulations with exact theoretical expectations for substitution rates and for the stationary distribution of codons at sites with a given fitness profile.

## What is tested

| Configuration | Setup | Theory compared |
|---|---|---|
| `neutral_mode` | SLiM-Tree `-N` (no fitness callback), default recombination (fully linked) | Substitutions per nucleotide site per generation = μ (Kimura 1968) |
| `flat_profiles` | One uniform stationary distribution (every amino acid equally fit) | Sense-codon substitution rate (stops lethal), dN/dS = 1, uniform stationary codon distribution, no premature stop fixes |
| `site_profiles` | Three site-heterogeneous profiles (strong, graded and weak amino-acid preferences), free recombination | Synonymous and nonsynonymous substitution rates, dN/dS, stationary amino-acid and codon distributions, no premature stop fixes |

For the profile configurations, the suite also checks that SLiM-Tree's fitting step (`utils/calculateFitnesses.py`) recovers the intended fitness values from the stationary distributions it is given.

### Theory (`theory.py`)

The theory is written independently of SLiM-Tree. It uses the same model that SLiM-Tree's fitness profiles are calibrated against:

- **Population and mutation:** a diploid Wright–Fisher population of N individuals with Jukes–Cantor mutation, so the total rate per nucleotide site is μ (SLiM-Tree `-v`).
- **Selection:** a change between amino acids with fitness values f_i and f_j has selection coefficient s = f_j − f_i. It fixes with Kimura's probability u(s) = (1 − e^(−2s)) / (1 − e^(−4Ns)).
- **Stop codons:** lethal, so each site is a Markov chain on the 61 sense codons.
- **Weak mutation:** fixed states change by single substitutions at rate q_ij = 2N (μ/3) u(f_j − f_i).
- **Outputs:** the stationary distribution, the expected synonymous and nonsynonymous substitution rates, and the mutation–selection dN/dS of Spielman & Wilke (2015, *MBE*).

`test_theory.py` (no SLiM; about 2 s) checks this module in three ways:
- against the closed form π ∝ exp((4N − 2) f) (Sella & Hirsh 2005);
- against SLiM-Tree's own calibration matrix (`Compute_matrix`, `Calc_stationary_KLD`);
- by checking that SLiM-Tree's fitter (`Find_fits`) recovers a known profile.

### How the simulations are run (`harness.py`)

- **Script generation:** each replicate uses SLiM-Tree's real pipeline: `readInput`, `findFitness`/`calculateFitnesses`, `cladeReader`, `writeSLiM`. The resulting script is the one `slim-tree` would write. The suite only stops SLiM-Tree from launching SLiM itself.
- **The observer:** the one addition to the script is an observer event. Every `--interval` generations it writes one genome's sequence and the positions of all segregating mutations. It does not modify the population.
- **Run structure:** the tree is `(A:1);`, and measurements are taken in the root population during SLiM-Tree's burn-in (`-b`). That burn-in is set to the equilibration period plus the measurement period.
- **Starting state:** SLiM-Tree draws each ancestral codon from its profile's stationary distribution, so sites start at stationarity. The first 10N generations let polymorphism reach equilibrium and are not scored.
- **Background selection:** the profile configurations use free recombination (`-r 0.5`), so there is none and the single-site theory applies.

### Statistics

- **Substitution rates:** observed / expected, with a standard error equal to the larger of the Poisson error and the between-replicate error. A check passes if |z| < 3.5.
- **Stationary distributions:** total variation distance (TVD) between the observed occupancy of fixed states and theory. The null distribution comes from simulating the theoretical Markov chain with the same number of sites, duration and sampling, starting at stationarity (parametric bootstrap). A check passes if p ≥ 0.001.
- **Scale:** each stationary check also reports how far a neutral (mutation-only) distribution would be from the target.
- **The theory uses the fitness values SLiM-Tree actually wrote into the script,** so these checks test the simulation rather than the fitting.

## Running

Weak-mutation theory requires θ = 4Nμ ≤ 0.01. θ = 0.001 is safer, but runs take about 10 times longer for the same number of substitutions.

```bash
python -m pytest tests/wf_validation/test_theory.py
python tests/wf_validation/suite.py --outdir tests/wf_validation/runs/<name> --jobs 22
python tests/wf_validation/suite.py --model nonWF --outdir tests/wf_validation/runs/<name> --jobs 22
WFV_RESULTS=tests/wf_validation/runs/<name>/results.json python -m pytest tests/wf_validation/test_wf_validation.py
```

- **Defaults:** N = 100, θ = 0.01, 300 codons, 100,000 measured generations after 1,000, sampled every 100 generations.
- **Replicates:** 10 neutral-mode, 10 flat-profile and 20 site-profile replicates.
- **Outputs:** `suite.py` writes `report.md`, `results.json` and one directory per replicate (generated script, SLiM output, observer samples, parsed counts).
- **Options:** `--only` runs a subset of configurations; `--theta 0.001 --measure-gens 1000000` is the slower, safer setting.

## Fitness model under test

Wright–Fisher scripts compute fitness multiplicatively. The fitness profiles are fitted (`utils/calculateFitnesses.py`) treating the difference between two amino acids' values as the selection coefficient of a change at one site, so each codon whose amino acid differs from the ancestor multiplies a haplosome's selection factor by exp(f_new − f_ancestral). An individual's fitness is the product of its two haplosomes' factors. Earlier versions added per-site values across the sequence and averaged the two haplosomes, so a change at one site moved fitness by only a small fraction of its intended selection coefficient. On that code the profile configurations fail: nonsynonymous substitution rates were 3–4 times theory for the strong and graded profiles, and the stationary distributions were far from their targets (total variation distance 0.53–0.59, against 0.64–0.78 for a fully neutral distribution). Non-WF models use the same fitness; see below.

- **Null alleles:** a gene is null in a haplosome if its start codon or terminal stop codon has changed, or if it carries a premature stop codon. Each gene multiplies fitness by 1 with two functional copies, 0.5 with one and 0 with none.
- **Pseudogenization:** in a haplosome where a gene is null, that gene's codons contribute nothing to the selection factor. A null protein's sequence is unconstrained.
- **Several genes:** null alleles are scored per gene, so null copies of different genes complement each other.

## Non-WF models

SLiM-Tree's non-WF scripts (`-w`) give every individual b offspring per tick (b = 1), each with a random mate. Parents and offspring then survive together with probability w·K/Σw, where w is the fitness above and K the population size. That keeps the expected size at K whatever the mean fitness.

The top-level README ("Non-WF models") and `utils/nonWFLifeCycle.py` state the life cycle's assumptions and its consequences in general form. With b = 1:
- one generation is 2 ticks;
- the effective size is 16K/11 per tick;
- SLiM-Tree multiplies selection coefficients by 11/16 (`NONWF_SELECTION_SCALE`) so that scaled selection matches the 4Kf the profiles are fitted for.

`theory.py` derives the same quantities independently for any b, including the fixation probability (a diffusion approximation). `test_theory.py` checks that the two derivations agree.

With `--model nonWF`, rates are per tick. The harness reads the selection scale and b from the generated script. Each profile also reports how far the non-WF stationary distribution is from the Wright–Fisher target, which is the distribution the profile was fitted to.

**Results** (SLiM 5.2; N = 100, θ = 0.01, 300 codons, 100,000 ticks; 10/10/20 replicates). All 25 scored checks pass:
- The neutral substitution rate is 0.986 ± 0.009 of μ/2 per tick.
- Site-profile rates are between 0.950 and 1.013 of theory. The Wright–Fisher suite shows the same small average deficit, about 2%.
- The stationary-distribution p-values range from 0.047 to 0.89.
- Non-WF theory is within codon TVD 0.0012 of the Wright–Fisher targets.

## SLiM 5

- **SLiM version.** SLiM-Tree generates SLiM 5 code: `haplosomes` instead of `genomes`, and `individual.haploidGenome1/2` instead of `genome1/2`. The suite was validated on SLiM 5.2. SLiM 5.1 and later use a new random-number generator, so seeds do not reproduce SLiM 4 runs.
- **Choosing the binary.** `WFV_SLIM=/path/to/slim` selects the SLiM binary for every test. Each run records the binary path, its version and a SHA-256 prefix in `slim_version.txt`.
- **Extra arguments.** `WFV_SLIMTREE_ARGS="..."` adds slim-tree arguments to every non-neutral run; `suite.py --slimtree-args` does the same for one suite run.
- **Fitness fitting.** `utils/calculateFitnesses.ratio_func` solves for selection coefficients in closed form, s = −ln(ratio)/(4N − 2), which is exact for Kimura's fixation probability. The numerical optimizer it replaces lost accuracy with N: errors in 4Ns reached about 0.1 at N = 1,000 and about 10 at N = 10,000. `test_theory.py` checks profile recovery at N = 100, 1,000 and 10,000.

## Scope and limits

- **Model type:** Wright–Fisher models (default) and SLiM-Tree's non-WF model (`--model nonWF`). The non-WF theory is a diffusion approximation for the life cycle described above.
- **Script writer:** it tests the local writer (`writeSLiM`). The HPC writer (`writeSLiMHPC`) inherits the same `write_initialize`, `write_fitness` and `write_repeated_commands`, so the fitness model is the same code, but HPC orchestration is not exercised.
- **Genome and mutation model:** a single gene (coding ratio 1) under Jukes–Cantor mutation. Mutation matrices and multiple genes are not covered.
- **dN/dS outputs:** SLiM-Tree's own dN/dS output (`-S`, `calculateSelectionDenominators.py`) is not checked here. The suite computes substitution counts independently from the observer.
- **Stop codons:** the theory treats stops as lethal. That matches the calibration (stationary distributions exclude stops) and the null-allele model above, under which a premature stop cannot fix in a single-gene genome.
