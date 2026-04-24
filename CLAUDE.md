# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

**Install:**
```bash
pip install .
```

**Run:**
```bash
slim-tree <input_tree.txt> <stationary_distributions.csv> [options]
# Example with substitution conversion:
slim-tree my_tree.txt stationary_distributions.csv -s -n 100 -v 2.5e-6
```

**Run tests:**
```bash
python3 -m unittest discover -s tests
```

**Run a single test file:**
```bash
python3 -m unittest tests.testFindFitness
```

**Run an example:**
```bash
cd Examples/Ex_1_Default_Simulation
slim-tree ex1_tree.txt stationary_distributions.csv
```

## Architecture

SLiM-Tree takes a Newick phylogenetic tree and generates/runs [SLiM](https://messerlab.org/slim/) forward-time simulation scripts. The pipeline runs in this order:

### Entry point
`SLiMTree.py` — the `SLiMTree` class orchestrates the full pipeline: read input → process fitness → read clades → write and execute SLiM scripts. It includes a `while` loop to regenerate the ancestral sequence if it already carries the amino acid targeted by a profile shift.

### Central data structure
`start_params` — a flat dictionary assembled by `readInput` and passed through every pipeline stage. It holds all simulation parameters (population size, mutation rate, genome length, fitness profiles, ancestral sequence, file paths, etc.).

### Pipeline modules (`utils/`)

| Module | Responsibility |
|---|---|
| `readInput.py` | Parses CLI args via `argparse`, validates them, creates output directories, saves parameters to YAML |
| `findFitness.py` | Reads codon stationary distributions CSV, computes fitness profiles via KLD optimization (or uses provided `-fd` file), generates the ancestral codon sequence |
| `calculateFitnesses.py` | Low-level fitness optimization using `scipy` — minimizes Kullback-Leibler divergence between SLiM-simulated and target stationary distributions |
| `findCoding.py` | Partitions the genome into coding/non-coding regions given genome length, number of genes, and coding ratio |
| `calculateSelectionDenominators.py` | Pre-computes dN and dS denominators for dN/dS calculations when `-S` flag is used |
| `cladeReader.py` | Parses the Newick tree with BioPython `Phylo`, recursively walks the tree depth-first, and returns a list of per-branch dicts sorted by `dist_from_start` (breadth-first execution order required by SLiM) |
| `convertTree.py` | Converts branch lengths from substitutions/site to generations using population size and mutation rate |
| `writeSLiM.py` | Base class that writes `.slim` scripts for local execution (WF and non-WF models) |
| `writeSLiMHPC.py` | Subclass that writes Slurm `.sh` job scripts instead of running locally |

### Supporting data
`fitnessDataFiles/slim_codon_nums.csv` — maps codon triplets to SLiM's internal codon index numbering; required for constructing fitness callbacks in generated scripts.

`fitnessDataFiles/table_stationary_distributions.csv` — example stationary distributions used in tests and examples.

### Output structure
All output is written relative to the input tree file's directory:
- `slimScripts/` — generated `.slim` (or `.sh`) scripts, one per branch
- `nuc_FASTA/` — nucleotide FASTA at each tip
- `aa_FASTA/` — amino acid FASTA at each tip
- `<tree_name>_parameters.yaml` — record of all run parameters

### Key design notes
- The `clade_dict_list` produced by `cladeReader` must be sorted by `dist_from_start` because SLiM scripts for each branch reference the previous branch's state and must be launched in chronological order.
- When `-hpc` is set, scripts are submitted via `sbatch` rather than executed directly; only population-size changes (and a few other parameters) can be varied per-branch in local mode, while HPC mode allows more per-branch overrides.
- Branch-specific parameter overrides (population size, mutation rate, fitness profile shifts) are specified via a YAML file passed with `-d`; `cladeReader.read_clade_data` translates short CLI abbreviations to full key names before merging into each branch's dict.
