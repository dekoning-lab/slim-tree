# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Project Is

SLiM-Tree is a Python CLI tool that automates the generation and execution of [SLiM](https://messerlab.org/slim/) (population genetics simulator) scripts from Newick phylogenetic trees. It simulates sequence evolution across multiple clades with fitness effects derived from codon stationary distributions, producing nucleotide and amino acid FASTA sequences at tree tips.

**External runtime dependencies:** SLiM must be installed and on PATH. An R installation may be required if using fitness distribution modes that invoke R optimization.

## Installation

```bash
pip install .
```

This installs the `slim-tree` CLI entry point (defined in `pyproject.toml` as `SLiMTree:SLiMTree`).

## Running Tests

Tests are plain Python unittest files in `tests/`. Run them individually or with pytest:

```bash
# Run a single test file
python -m pytest tests/testConvertTree.py

# Run all tests
python -m pytest tests/
```

No linting configuration is present in this repo.

## Architecture

### Data Flow

```
CLI args → param_dict → tree parsing → fitness computation → SLiM script generation → execution
```

1. **`utils/readInput.py`** — Argparse CLI parsing; produces `param_dict` consumed by the rest of the pipeline.
2. **`utils/cladeReader.py`** — Parses the Newick tree (via BioPython) into per-clade data structures. Handles branch-specific parameter overrides from the YAML file (`-d` flag).
3. **`utils/findFitness.py`** — Converts codon stationary distributions into per-amino-acid fitness values. Calls `utils/calculateFitnesses.py` for the underlying optimization.
4. **`utils/findCoding.py`** — Assigns coding regions and fitness profiles, generates the ancestral sequence.
5. **`utils/writeSLiM.py`** / **`utils/writeSLiMHPC.py`** — Generates `.slim` scripts (and Slurm `.sh` job scripts in HPC mode). One script per population is produced for multi-clade trees.
6. **`SLiMTree.py`** — Top-level orchestrator. Contains a `while(run_random_process)` loop that re-runs fitness/ancestral-sequence generation if the ancestral sequence accidentally contains a target amino acid during a fitness profile shift (rare but possible).

### Execution Modes

- **Local mode** (default) — runs `slim script.slim` via `os.system()` sequentially.
- **HPC mode** (`-hpc` flag) — generates Slurm submission scripts instead of running SLiM directly. Enables additional branch-level overrides (more parameters than local mode allows).
- **Neutral evolution** (`-N` flag) — skips fitness computation entirely.

### Branch-Level Parameter Overrides

The `-d <yaml>` flag accepts a YAML file mapping branch names to per-branch parameter overrides (population size, mutation rate, recombination rate, sample size, split ratio, or fitness profile shifts). HPC mode supports a superset of the parameters available in local mode.

### Output Location

All outputs (`.slim` scripts, FASTA sequences, parameter YAML) are written to the **directory containing the input tree file**, not the working directory.

### Key Files

| File | Purpose |
|---|---|
| `SLiMTree.py` | Main entry point and workflow orchestrator |
| `utils/readInput.py` | CLI argument parsing and validation |
| `utils/cladeReader.py` | Newick tree traversal and clade data structures |
| `utils/findFitness.py` | Fitness profile computation from stationary distributions |
| `utils/calculateFitnesses.py` | Fitness optimization algorithm (called by findFitness) |
| `utils/writeSLiM.py` | SLiM script generation (local mode) |
| `utils/writeSLiMHPC.py` | SLiM + Slurm script generation (HPC mode) |
| `utils/convertTree.py` | Branch length conversion (substitutions/site → generations) |
| `utils/calculateSelectionDenominators.py` | dN/dS calculation support (used with `-S` flag) |

### Examples and Post-Processing

- `Examples/` — 10 worked examples covering different tree topologies and parameter combinations.
- `DataPostProcessing/` — R/Python scripts for analyzing simulator output.
