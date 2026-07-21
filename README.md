# SLiM-Tree

SLiM-Tree automates forward-time population genetics simulations over phylogenetic timescales under realistic models of sequence-fitness relationships. It generates and runs [SLiM](https://messerlab.org/slim/) scripts from a Newick phylogenetic tree, producing substitution histories and polymorphisms across a whole clade without the simplifying assumptions of standard phylogenetic models (infinite sites, mutation-limited evolution, no polymorphism, etc.).

## Requirements

Before installing SLiM-Tree, ensure the following are available on your system:

- **Python 3.8+**
- **SLiM 4.x** — [messerlab.org/slim](https://messerlab.org/slim/)
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
