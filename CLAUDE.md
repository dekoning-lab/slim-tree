# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Running SLiM-Tree

```bash
python3 slim-tree/ <input_tree> <codon_stationary_distributions> [options]
```

Example with common options:
```bash
python3 slim-tree/ tree.txt stationary_distributions.csv -n 100 -v 2.5e-6 -g 300
```

HPC mode (requires Slurm):
```bash
python3 slim-tree/ tree.txt stationary_distributions.csv -hpc -p <partition> -t <time>
```

## Running Tests

Run all tests from the `slim-tree/` directory:
```bash
python3 -m pytest tests/
```

Run a single test file:
```bash
python3 -m unittest tests/testReadInput.py
python3 -m unittest tests/testCladeReader.py
python3 -m unittest tests/testFindCoding.py
python3 -m unittest tests/testFindFitness.py
python3 -m unittest tests/testCalculateSelectionDenominators.py
```

Run a specific test method:
```bash
python3 -m unittest tests.testReadInput.testReadInput.test_check_arguments
```

## Architecture Overview

SLiM-Tree is a Python wrapper that generates and runs [SLiM](https://messerlab.org/slim/) simulation scripts from a Newick phylogenetic tree. The pipeline has four stages:

1. **Input parsing** (`utils/readInput.py`) — reads CLI arguments, validates them, builds a `start_params` dict, and creates output directories (`slimScripts/`, `aa_FASTA/`, `nuc_FASTA/`, and optionally `backupFiles/`, `slurmOutput/`, etc.).

2. **Fitness processing** (`utils/findFitness.py`, `utils/findCoding.py`, `utils/calculateSelectionDenominators.py`) — reads codon stationary distributions or amino acid fitness files, determines coding regions, assigns fitness profiles to codon positions in the genome, and computes dN/dS denominators when needed. The outer `while` loop in `SLiMTree.__init__` retries this stage if the ancestral sequence happens to match a profile-shift target.

3. **Tree traversal** (`utils/cladeReader.py`) — parses the Newick tree via BioPython, recursively traverses clades depth-first, and produces a list of per-clade parameter dictionaries sorted by distance from simulation start. A YAML `tree_data_file` (`-d` flag) can override population size, mutation rate/matrix, sample size, or trigger a `profile_shift` on a per-branch basis.

4. **SLiM script generation and launch** (`utils/writeSLiM.py`, `utils/writeSLiMHPC.py`) — writes one `.slim` (or `.sh` for HPC) script per clade. Scripts are chained sequentially; each script picks up where the parent population left off. The root simulation is then launched via `os.system("slim ...")` or `os.system("sbatch ...")`.

### Key data flow

`start_params` dict is the central state object threaded through all stages. It holds all CLI parameters plus computed values like `ancestral_sequence`, `fitness_profiles`, `fitness_profile_nums`, `coding_seqs`, `scaling_value`, `filenames`, and the `fitness_finder` object itself.

### Mutation model

Two modes controlled by `jukes_cantor` (set `True` when no `-m` matrix is supplied):
- **Jukes-Cantor** (`jukes_cantor=True`): uniform mutation rate `-v`
- **Matrix model** (`jukes_cantor=False`): 4×4 CSV matrix supplied via `-m`

Per-branch overrides of mutation rate/matrix are handled in `cladeReader.read_clade_data()`.

## Requirements

- Python 3, R, SLiM
- Protein-based fitness: Java and C also required
- Python packages: `BioPython`, `matplotlib`, `pandas`, `numpy`, `yaml`, `argparse`
- R packages: `dplyr`, `BB`, `data.table`, `optparse`, `seqinr`, `doParallel`, `Rfast`

## Output Structure

Generated relative to the input tree file location:
- `slimScripts/` — generated `.slim` / `.sh` scripts (one per population/clade)
- `aa_FASTA/` — amino acid FASTA outputs
- `nuc_FASTA/` — nucleotide FASTA outputs
- `backupFiles/` — simulation checkpoints (only with `-B`)
- `slurmOutput/` — Slurm stdout/stderr (only with `-hpc`)
- `selectionCalculationOutput/` — dN/dS data (only with `-S`)
- `substitutionCountingOutput/` — substitution counts (only with `-c`)
- `polymorphicSites/` — polymorphism data (only with `-P`)

## Post-Processing

The `DataPostProcessing/` folder contains standalone scripts for downstream analysis:
- `process_data.py` — Python post-processing
- `find_branch_lengths.py` — branch length analysis
- `FindPercentPolymorphic.R`, `perform_analysis_fixation_counts.R`, `perform_statistical_analysis_branch_lengths.R` — R analysis scripts
