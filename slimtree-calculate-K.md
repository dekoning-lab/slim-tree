# Calculating Substitution Rate K from SLiM-Tree Output

## What K is here

K = substitutions per site per generation, measured empirically from the simulation.
It is not assumed (e.g. K = μ under neutrality) — it is counted directly from what SLiM fixed.

---

## Step 1 — Add `-c` to your slim-tree call

```r
options <- c(
  "../../tree",
  "../../table_stationary_distributions_neutral.csv",
  "-hpc",
  "-t", "05:00:00",
  "-p", "cpu2023,cpu2022,cpu2021,cpu2025,cpu2019",
  "-v", "0.00025",
  "-fd", "../../table_fitness_dists.csv",
  "-g", "3333",
  "-c"    # enables substitution counting
)
```

This creates `substitutionCountingOutput/<clade_name>_num_subs.txt` inside each replicate folder.
Each file contains a single integer: the number of nucleotide sites that fixed a new allele along that branch.

---

## Step 2 — Per-branch K

For a single branch, K is:

```
K = raw_count / genome_length_nt / branch_length_in_generations
```

where `genome_length_nt = 3333 × 3 = 9999` (your `-g 3333` codons as nucleotides),
and `branch_length_in_generations` is the length of that branch in your tree file.

```r
library(ape)

tree <- read.tree("../../tree")
genome_length_nt <- 3333 * 3

# Build a named vector of branch lengths: tip/node label -> length
# edge.length[i] is the length of the branch leading INTO tree$edge[i,2]
node_labels <- c(tree$tip.label, tree$node.label)
branch_lengths <- setNames(tree$edge.length, node_labels[tree$edge[, 2]])

# Read per-branch substitution counts from one replicate
subs_dir <- "substitutionCountingOutput"
subs_files <- list.files(subs_dir, pattern = "_num_subs.txt", full.names = TRUE)

K_per_branch <- sapply(subs_files, function(f) {
  clade_name  <- sub("_num_subs.txt$", "", basename(f))
  raw_count   <- as.integer(readLines(f))
  branch_len  <- branch_lengths[clade_name]
  raw_count / genome_length_nt / branch_len
})
names(K_per_branch) <- sub("_num_subs.txt$", "", basename(subs_files))

print(K_per_branch)
```

---

## Step 3 — Global K (whole-tree average)

Global K pools all substitutions across all branches and divides by total evolutionary time:

```
K_global = sum(raw_counts across all branches) / genome_length_nt / sum(all branch lengths)
```

This is a branch-length-weighted average. Branches that ran longer contribute proportionally more.
Counts are additive because each branch measures substitutions relative to its own starting sequence
(the parent branch's ending sequence), so there is no double-counting.

```r
library(ape)

tree <- read.tree("../../tree")
genome_length_nt <- 3333 * 3
total_generations <- sum(tree$edge.length)

subs_dir <- "substitutionCountingOutput"
subs_files <- list.files(subs_dir, pattern = "_num_subs.txt", full.names = TRUE)

total_subs <- sum(sapply(subs_files, function(f) as.integer(readLines(f))))

K_global <- total_subs / genome_length_nt / total_generations

cat("K_global =", K_global, "\n")
```

---

## Step 4 — Across replicates

```r
library(ape)

tree <- read.tree("../../tree")
genome_length_nt <- 3333 * 3
total_generations <- sum(tree$edge.length)

replicate_dirs <- list.dirs("001_005_rep_slimtree_output", recursive = FALSE)

K_global_replicates <- sapply(replicate_dirs, function(rep_dir) {
  subs_files <- list.files(file.path(rep_dir, "substitutionCountingOutput"),
                           pattern = "_num_subs.txt", full.names = TRUE)
  if (length(subs_files) == 0) return(NA)
  total_subs <- sum(sapply(subs_files, function(f) as.integer(readLines(f))))
  total_subs / genome_length_nt / total_generations
})

cat("K per replicate:\n")
print(K_global_replicates)
cat("\nMean K:", mean(K_global_replicates, na.rm = TRUE), "\n")
cat("SD K:",   sd(K_global_replicates,   na.rm = TRUE), "\n")
```

---

## Notes

- Branch lengths in your tree must be in **generations** (not substitutions/site). If you used `-s`, the tree was already converted internally by slim-tree and your original file is still in substitutions/site — use the converted lengths, which slim-tree stores in `_parameters.yaml`.
- Under neutral evolution with a uniform mutation model, K should converge to μ = 0.00025. Under selection it will be lower.
- The `-c` flag slows down the SLiM simulation slightly (one identity-by-state check per generation per branch).
