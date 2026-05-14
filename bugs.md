# SLiM-Tree Bug Tracker

Discovered by systematic feature testing (2026-04-23). Update status as fixes are applied.

---

## Status key
- 🔴 Open — not yet fixed
- 🟡 In progress
- 🟢 Fixed — include commit or PR reference

---

## Bugs

### Bug 1 — `-m` (mutation matrix) crashes without `-fd` ❌
**Status:** 🟢 Fixed — `SLiMTree.py:70` `elif jukes_cantor` → `else`

**Symptom:**
```
AttributeError: 'findFitness' object has no attribute 'fitness_mat'
```

**Reproduction:**
```bash
slim-tree tree.txt stat_dists.csv -m mut_mat.csv -n 50 -b 2 -g 10
```

**Root cause:** `SLiMTree.py:70` has `elif (start_params["jukes_cantor"]):` before calling `find_optimal_fitnesses()`. When `-m` is supplied, `jukes_cantor=False`, so the branch is skipped and `self.fitness_mat` is never set. `process_fitness_dists()` then crashes accessing it.

**Also breaks:** Example `Ex_8_New_Mutation_Rate_Matrix` (uses `-m` without `-fd`).

**Fix:** `SLiMTree.py:70` — change `elif (start_params["jukes_cantor"]):` to `else:` so fitness optimisation runs regardless of mutation model when no `-fd` is provided.

---

### Bug 2 — Two output files written to CWD instead of output directory ❌
**Status:** 🟢 Fixed — `writeSLiM.py:123` uses `os.path.dirname(filenames[1])`; `findFitness.py:127` now also uses `os.path.join(os.path.dirname(output_dir), ...)` (previously applied incorrectly)

**Symptom:** `fitness_profile_nums.txt` and `table_fitness_dists.csv` appear in whatever directory the command is run from, not in the simulation output folder. They are silently overwritten on each run.

**Root cause:**
- `utils/writeSLiM.py:123` — `open("fitness_profile_nums.txt", "w")` — bare filename, no path
- `utils/findFitness.py:127` — `self.fitness_mat.to_csv("table_fitness_dists.csv", header=False)` — bare filename, no path

**Fix:** Use `os.path.dirname(filenames[1])` to get the output directory. Note: outputs are written to CWD (by design — allows replicate workflows where input files are shared but each replicate runs from its own directory).

---

### Bug 3 — `-s` branch-length conversion corrupts bootstrap integers ❌
**Status:** 🟢 Fixed — `convertTree.py:19,28` regex changed from `r'\d+\.?\d*'` to `r'\d+\.\d+'`

**Symptom:** A tree with bootstrap support values like `((A:0.1,B:0.2)90:0.3);` has the integer `90` treated as a branch length, scaling it to something like `600599613`.

**Root cause:** `utils/convertTree.py:19` uses the regex `r'\d+\.?\d*'` in the replacement pass, which matches **all** numbers including integers. Only the summing step correctly restricts to decimals via `r'\d+\.\d+'`; the replacement step does not.

**Fix:** `convertTree.py:19` — change replacement regex from `r'\d+\.?\d*'` to `r'\d+\.\d+'` so integers (bootstrap values, clade labels) are left untouched.

---

### Bug 4 — Example 9 parameter file is invalid YAML ❌
**Status:** 🟢 Fixed — `ex9_pop_parameters.txt` rewritten as valid YAML; `ex9_input.txt` updated to use `-hpc` flag

**Symptom:**
```
yaml.scanner.ScannerError: found character '@' that cannot start any token
  in "Examples/Ex_9_Varying_Population_Parameters/ex9_pop_parameters.txt", line 1, column 1
```

**Root cause:** `ex9_pop_parameters.txt` uses a legacy `@BranchName / -n 200` CLI-style format. The code calls `yaml.safe_load()` and expects valid YAML.

**Fix:** Rewrite `Examples/Ex_9_Varying_Population_Parameters/ex9_pop_parameters.txt` as valid YAML:
```yaml
A:
  n: 200
C:
  n: 150
B:
  n: 600
```
Also remove the `-v` (mutation rate) lines — mutation rate changes are only allowed in HPC mode.

---

### Bug 5 — Profile shifts (`ps:`) silently blocked in non-HPC mode ❌
**Status:** 🟢 Closed (by design) — profile shifts are intentionally HPC-only. Documented in README under "Branch-specific parameters".

---

### Bug 6 — Non-WF + `-S` generates a SLiM script that crashes at runtime ❌
**Status:** 🟢 Fixed — added `dN_`, `dS_`, `subs_` initialisation to child-population block in `write_subpop_nonwf()`

**Symptom:** SLiM exits with a key-not-found error when the simulation reaches a child population:
```
ERROR (SLiM): sim.getValue("dN_p2") — key not found
```

**Reproduction:**
```bash
slim-tree tree.txt stat_dists.csv -w -S -n 50 -b 2 -g 10
```

**Root cause:** `utils/writeSLiM.py` — `write_subpop_nonwf()` initialises `fixations_` and `fixations_counted_` for child populations but never initialises `dN_`, `dS_`, or `subs_`. The `write_repeated_commands()` function then tries to read those uninitialised keys.

**Fix:** Add to `write_subpop_nonwf()`, in the child-population branch, the same three `sim.setValue()` calls already present in `write_subpop()` at lines 359–361:
```python
"\n\tsim.setValue(\"dN_"+ pop_name +"\", 0);" +
"\n\tsim.setValue(\"dS_"+ pop_name +"\", 0);" +
"\n\tsim.setValue(\"subs_"+ pop_name +"\", \"\\n\\nSubstitutions:\");"
```

---

### Bug 7 — Bundled `fitnessDataFiles/table_stationary_distributions.csv` unusable as input ❌
**Status:** 🟢 Fixed — correctly formatted `stationary_distributions.csv` (61 rows, codon triplet index) added to each `Examples/Ex_*/` directory.

---

### Bug 8 — `-c`, `-S`, `-P` write output for internal (non-tip) populations
**Status:** 🟢 Closed (by design) — internal-node output is intentional. Substitution counts, dN/dS, and polymorphism files for ancestor populations are useful for studying what happened along internal branches of the tree.

---

## Code Quality Issues

### CQ1 — Dead allocation `integer(row_num*1500)` in substitution tracking
**Location:** `utils/writeSLiM.py:275`

```python
"\n\t\tmuts_mat = integer(row_num*1500);\n\t\tmuts_mat = " + pop_name + ".genomes.nucleotides(...);"
```

The first assignment creates a zero-filled array of `row_num × 1500` elements that is immediately discarded when the second line reassigns `muts_mat`. For a 1000-codon genome with 100 diploid individuals, this needlessly allocates 150,000 integers per generation tick. Remove the first line entirely.

---

### CQ2 — Unused variable `ndists` in `calculateSelectionDenominators`
**Location:** `utils/calculateSelectionDenominators.py:119`

```python
ndists = range(self.stationary_distributions.shape[1]-1)
```

This variable is defined and never referenced. Delete it.

---

### CQ3 — `fixations_p1` initialised twice per run
**Locations:** `utils/writeSLiM.py:148` (inside `setup_fitness()`) and `utils/writeSLiM.py:38` (in `1 early()` block)

`setup_fitness()` always sets `fixations_p1`, then `1 early()` calls `setup_fitness()` and immediately sets it again. The second assignment is redundant. Remove the duplicate from `setup_fitness()`.

---

### CQ4 — `dN_p*/dS_p*/subs_p*` always initialised for child populations regardless of `-S`
**Location:** `utils/writeSLiM.py:358–361` in `write_subpop()`

These three SLiM variables are written unconditionally for every child population even when `-S` (calculate_selection) is not set. Wrap lines 358–361 in:
```python
if population_parameters["count_subs"] or population_parameters["calculate_selection"]:
```

---

### Bug 9 — HPC race condition: last sibling deletes parent's chain files before slower siblings start
**Status:** 🟢 Fixed — removed `rm` cleanup block from `writeSLiMHPC.write_end_sim` (`writeSLiMHPC.py:217–220`)

**Symptom:**
```
ERROR (Species::ExecuteContextFunction_initializeAncestralNucleotides): the file at path p1.fasta could not be opened or does not exist.
Error on script line 4, character 1:
   initializeAncestralNucleotides("p1.fasta");
```

**Reproduction:**
Any HPC run on a tree with at least one internal branching node, when sibling Slurm jobs experience unequal queue wait times.

**Root cause:** `write_end_sim()` appended `system("rm p1.txt"); system("rm p1.fasta"); ...` to the `last_child_clade`'s script. All sibling jobs are submitted simultaneously at the end of the parent's simulation. On a busy cluster the last-child job can get resources immediately and finish (deleting `p1.fasta`) while sibling jobs are still queued. When a sibling eventually starts, it tries to read `p1.fasta` in its `initialize()` block and fails.

**Fix:** Remove the three `rm` lines entirely (`writeSLiMHPC.py:217–220`). The intermediate `.txt`, `.fasta`, and `_fixed_mutations.txt` files are small (KB–MB range) and safe to leave on disk until all jobs complete.

---

## No-test-coverage gaps

The 68 unit tests all pass but none of them test the full pipeline end-to-end. The following scenarios have zero test coverage and produced runtime bugs above:

| Scenario | Bug triggered |
|---|---|
| `-m` without `-fd` | Bug 1 |
| `-s` with bootstrap integers in tree | Bug 3 |
| `-w -S` combined | Bug 6 |
| `-d` YAML with `ps:` key in local mode | Bug 5 |

Add integration tests that call `SLiMTree()` directly (or via subprocess) for each of these combinations and assert a zero exit code.
