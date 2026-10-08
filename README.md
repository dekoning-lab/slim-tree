# SLiM-Tree

SLiM-Tree automates forward-time population genetics simulations over phylogenetic timescales under realistic models of sequence-fitness relationships. It generates and runs [SLiM](https://messerlab.org/slim/) scripts from a Newick phylogenetic tree, producing substitution histories and polymorphisms across a whole clade without the simplifying assumptions of standard phylogenetic models (infinite sites, mutation-limited evolution, no polymorphism, etc.).

## Requirements

Before installing SLiM-Tree, ensure the following are available on your system:

- **Python 3.8+**
- **SLiM 5.x** — [messerlab.org/slim](https://messerlab.org/slim/). Generated scripts use the SLiM 5 haplosome API and do not run on SLiM 4. Validated on SLiM 5.2.
- **R** with packages: `dplyr`, `BB`, `data.table`, `optparse`, `seqinr`, `doParallel`, `Rfast`
- **Java and C** — only required when computing fitness effects from protein structure

## Installation

```bash
git clone git@github.com:dekoning-lab/slim-tree.git 
cd slim-tree
git checkout main
pip install .
```

This installs the `slim-tree` command globally.

## Quick start

```bash
slim-tree my_tree.txt stationary_distributions.csv
```

- `my_tree.txt` — Newick tree with branch lengths **in generations**
  If your tree has branch lengths in **substitutions per site**, use `-s` to convert automatically: [NEED TO RUN MORE TEST TO MAKE SURE -s WORKS - WORKING ON IT]
- `stationary_distributions.csv` — codon stationary distributions (61 codons × N profiles)


```bash
slim-tree my_tree.txt stationary_distributions.csv -s -n 100 -v 2.5e-6
```

## Output

All output is written relative to the input tree file location:

- `slimScripts/` — generated `.slim` scripts (one per population)
- `nuc_FASTA/` — nucleotide FASTA sequences at each tip
- `aa_FASTA/` — amino acid FASTA sequences at each tip
- `<tree_name>_parameters.yaml` — record of all parameters used

## Options

```
positional arguments:
  input_tree                    Newick tree file; branch lengths must be in generations
                                (use -s if they are in substitutions per site)
  codon_stationary_distributions
                                CSV of codon stationary distributions (61 rows × N columns,
                                no header, codon names as row index)

optional arguments:
  -h, --help                    show this help message and exit
  -fd FILE                      CSV of amino acid fitness values; if omitted, fitnesses are
                                computed from the stationary distributions
  -s, --substitutions           convert branch lengths from substitutions per site to
                                generations using the provided -n and -v values
  -n INT                        starting population size (default: 100)
  -b INT                        burn-in multiplier × population size (default: 10)
  -v FLOAT                      per-site mutation rate (default: 2.5e-6)
  -r FLOAT                      recombination rate (default: 2.5e-8)
  -m FILE                       4×4 mutation rate matrix CSV (A,C,G,T order, no header,
                                diagonal = 0); overrides -v
  -g INT                        genome length in codons (default: 300)
  -G INT                        number of genes (default: 1)
  -C FLOAT                      coding ratio — fraction of genome that is coding (default: 1.0)
  -f FILE                       ancestral sequence FASTA (nucleotides); requires -fd
  -k INT|all|consensus          sample size per tip; 'all' = entire population,
                                'consensus' = consensus sequence (default: all)
  -sr FLOAT                     split ratio for non-WF models — fraction of parent
                                population going into first daughter (default: 0.5)
  -d FILE                       YAML file to override parameters on specific branches
                                (see Branch-specific parameters below)
  -w, --nonWF                   use a non-Wright-Fisher model
  -N, --neutral_evolution       run neutral evolution (no fitness effects)
  -c, --count_subs              count substitutions per branch (slows simulation)
  -o, --output_gens             output sequences every 100 generations
  -B, --backup                  save simulation checkpoints (increases disk/time usage)
  -P, --polymorphisms           record polymorphic sites at each branch end
  -S, --calculate_selection     compute dN/dS by counting synonymous and non-synonymous
                                fixed substitutions
  -hpc, --high_performance_computing
                                generate Slurm job scripts instead of running locally
  -p PARTITION                  Slurm partition name (required with -hpc)
  -t TIME                       maximum Slurm wall time, e.g. 12:00:00 (required with -hpc)
  -M MEMORY                     memory to request per Slurm job, e.g. 16g or 32000M
                                (optional; omit to use the cluster default)
```

### Non-WF models (`-w`)

> **Non-WF runs are calibrated to a specific life cycle.** SLiM-Tree corrects selection strength in non-WF runs with a factor derived for the life cycle below (`utils/nonWFLifeCycle.py`). If you change that life cycle (offspring numbers, age structure, which individuals face survival, how density is regulated), re-derive the factor; the values below no longer apply.

**Life cycle.** Each tick, every individual is the first parent of b offspring (b = 1), each with a random second parent. Then parents and offspring face survival together: an individual with fitness w survives with probability w·K/Σw, where K is the population size (`-n`). Fitness is multiplicative across sites, as in Wright–Fisher runs. This normalization keeps the expected population size at K whatever the mean fitness (soft selection), and leaves relative fitness unchanged.

**Assumptions** behind the values below:
- every individual faces survival every tick, with no age effects;
- there is no fecundity selection;
- mates are drawn at random;
- selection is genic (multiplicative per copy) and weak;
- K is large and roughly constant.

**Consequences, in general form.** With φ = b/(1 + b), the fraction of survivors that are newborns:

| Quantity | General form | b = 1 |
|---|---|---|
| Generation time | (1 + b)/b ticks | 2 ticks |
| Neutral substitutions per site per tick | μ·b/(1 + b) | μ/2 |
| Effective size per tick | 4K(1 + b)²/(b(4b + 7)) | 16K/11 |
| Effective size per generation | 4K(1 + b)/(4b + 7) | 8K/11 |
| Neutral diversity | 4μK·4(1 + b)/(4b + 7) | (8/11)·4Kμ |
| Selection scale applied (K / Ne per tick) | b(4b + 7)/(4(1 + b)²) | 11/16 |

**Derivation.**
- **Coalescence.** Two gene copies coalesce in one tick with probability [2φ(1 − φ) + φ²(4b − 1)/(4b)]/(2K):
  - an adult and a newborn: 1/(2K);
  - two newborns: 1/(2K), unless both copies come from their first parents, in which case they coalesce only if the newborns are siblings;
  - two adults: never.
- **Selection.** Survival selection changes an allele's frequency by s·p(1 − p) per tick, so scaled selection is 4·Ne_tick·s. The fitness profiles are fitted for 4Kf, so SLiM-Tree multiplies selection coefficients by K/Ne_tick. As b grows, the life cycle approaches Wright–Fisher and the scale tends to 1.

**Checks.**
- *Neutral, b = 1:* substitutions were 0.496 ± 0.005 μ per tick (predicted 0.5). Diversity relative to Wright–Fisher was 0.713 ± 0.007 (predicted 0.727).
- *General b:* Neutral runs with b = 1, 2 and 3 (the generated script's reproduction callback edited to give b offspring; 10 replicates each). Substitution rates were 0.497, 0.655 and 0.738 μ per tick (predicted 0.500, 0.667, 0.750). Diversity was 0.718, 0.783 and 0.827 of 4Kμ (predicted 0.727, 0.800, 0.842). The predicted dependence on b holds, but all six values are 1–2% below prediction at K = 100. For rates, the Wright–Fisher suite shows the same deficit, so that is the measurement. For diversity it is probably the small-population effects the derivation leaves out, such as census fluctuation and self-mating.
- *Site profiles:* `tests/wf_validation` (`--model nonWF`) checks that site profiles reproduce their stationary distributions. On SLiM 5.2 with b = 1, all 25 scored checks pass. Without the correction, scaled selection would be 16/11 of what the profiles intend.

**Practical consequences:**
- **Time units.** Non-WF branch lengths are in ticks. `-s` converts substitutions to Wright–Fisher generations, so a non-WF tree converted with `-s` accrues about 1/T of the intended divergence, where T is the generation time in ticks (half, for b = 1).
- **θ.** With the same `-n` and `-v`, non-WF neutral diversity is 8/11 of the Wright–Fisher value.
- **Run record.** Each run prints these values and records them in `<tree>_parameters.yaml` (`nonWF_*` keys).

## Branch-specific parameters

The `-d` flag accepts a YAML file that overrides simulation parameters on specific branches. Each top-level key must match a branch name in the Newick tree:

```yaml
A:
  n: 200
B:
  n: 600
  v: 0.000001
```

**In local mode** (no `-hpc`), only population size (`n` / `population_size`) may be changed per branch.

**In HPC mode** (`-hpc`), the following additional keys are supported:

| Key | Full name | Description |
|---|---|---|
| `n` | `population_size` | Population size |
| `v` | `mutation_rate` | Per-site mutation rate |
| `m` | `mutation_matrix` | 4×4 mutation rate matrix file |
| `r` | `recombination_rate` | Recombination rate |
| `k` | `sample_size` | Sample size at branch end |
| `sr` | `split_ratio` | Split ratio for non-WF branching |
| `ps` | `profile_shift` | Shift to a different fitness profile mid-tree |

> **Note:** Fitness profile shifts (`ps`) require `-hpc`. They are not supported in local mode.

## Examples

The `Examples/` directory contains 10 worked examples. Each folder has an input tree and a command showing the flags used. Run any example from inside its folder:

```bash
cd Examples/Ex_1_Default_Simulation
slim-tree ex1_tree.txt stationary_distributions.csv
```

## Post-processing

The `DataPostProcessing/` directory contains scripts for downstream analysis of SLiM-Tree output, including branch length estimation and dN/dS calculations.

## Testing

Unit tests for the script writers and input handling are in `tests/` (`python -m pytest tests`). `tests/wf_validation/` checks that Wright–Fisher and non-WF simulations reproduce weak-mutation theory for substitution rates and stationary codon distributions; see its README.
