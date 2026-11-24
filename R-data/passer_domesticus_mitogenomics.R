##############################################################################
# Passer domesticus (House Sparrow) Mitogenomics Analysis in R
# Companion code for: "Complete Theoretical Framework: R Programming for
# Passer domesticus Mitogenomics"
#
# Sections match the theory document:
#   1. Get DNA sequences from NCBI
#   2. Basic sequence stats (length, GC content)
#   3. Multiple sequence alignment (MUSCLE)
#   4. Phylogenetic tree (distance + Maximum Likelihood + bootstrap)
#   5. Genetic diversity (Haplotype diversity Hd, Nucleotide diversity pi)
#   6. Tajima's D (population expansion/contraction test)
#   7. Visualization (tree, diversity bar chart, Tajima's D plot)
##############################################################################

## ---- 0. Install & load required packages ---------------------------------
required_packages <- c(
  "rentrez",     # talk to NCBI (search + download sequences)
  "ape",         # phylogenetics, distance matrices, tree I/O
  "seqinr",      # sequence stats (GC content, etc.)
  "pegas",       # population genetics (Hd, pi, Tajima's D)
  "phangorn",    # Maximum Likelihood trees + bootstrap
  "ggplot2",     # plotting
  "DECIPHER"     # Bioconductor: fast, reliable MSA (AlignSeqs = MUSCLE-like)
)

cran_pkgs <- setdiff(required_packages, "DECIPHER")
new_cran  <- cran_pkgs[!(cran_pkgs %in% installed.packages()[, "Package"])]
if (length(new_cran)) install.packages(new_cran, dependencies = TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("DECIPHER", quietly = TRUE)) BiocManager::install("DECIPHER")

library(rentrez)
library(ape)
library(seqinr)
library(pegas)
library(phangorn)
library(ggplot2)
library(DECIPHER)

set.seed(42)  # reproducibility for bootstrap / random steps


## ---- 1. Fetch DNA sequences from NCBI -------------------------------------
# Search NCBI nucleotide database for house sparrow + related Passer species
# mitochondrial complete genomes / control region sequences.

search_term <- '(Passer[Organism]) AND mitochondrion[Title] AND complete genome[Title]'

search_result <- entrez_search(db = "nuccore", term = search_term, retmax = 25)
cat("Number of records found:", search_result$count, "\n")

# Download sequences in FASTA format
fasta_data <- entrez_fetch(db = "nuccore",
                            id = search_result$ids,
                            rettype = "fasta")

# Save raw sequences to disk
dir.create("data", showWarnings = FALSE)
writeLines(fasta_data, "data/passer_mtDNA_raw.fasta")

# Read them back in with seqinr for easy handling
seqs <- read.fasta("data/passer_mtDNA_raw.fasta", as.string = TRUE, forceDNAtolower = FALSE)
cat("Sequences downloaded:", length(seqs), "\n")


## ---- 2. Basic sequence statistics -----------------------------------------
seq_lengths <- sapply(seqs, function(x) nchar(x[1]))
gc_content  <- sapply(seqs, function(x) GC(s2c(x[1])) * 100)

seq_stats <- data.frame(
  ID        = names(seqs),
  Length_bp = seq_lengths,
  GC_pct    = round(gc_content, 2)
)
print(seq_stats)
write.csv(seq_stats, "data/sequence_stats.csv", row.names = FALSE)


## ---- 3. Multiple sequence alignment (MUSCLE-style) ------------------------
# DECIPHER::AlignSeqs implements a MUSCLE-like progressive/guide-tree alignment
dna_set   <- readDNAStringSet("data/passer_mtDNA_raw.fasta")
alignment <- AlignSeqs(dna_set, anchor = NA, verbose = FALSE)

# Save alignment
writeXStringSet(alignment, filepath = "data/passer_mtDNA_aligned.fasta")

# Convert to formats needed downstream (ape "DNAbin" and pegas)
aln_dnabin <- as.DNAbin(alignment)
cat("Alignment length (bp):", ncol(as.matrix(aln_dnabin)), "\n")


## ---- 4. Phylogenetic tree building -----------------------------------------
## 4a. Distance matrix (uncorrected p-distance, as in the theory doc's
##     "Distance Matrix" example)
dist_matrix <- dist.dna(aln_dnabin, model = "raw", pairwise.deletion = TRUE)
print(as.matrix(dist_matrix))

## 4b. Neighbor-Joining starting tree
nj_tree <- nj(dist_matrix)

## 4c. Maximum Likelihood tree with phangorn (GTR+Gamma model, common for mtDNA)
phang_align <- phyDat(as.matrix(aln_dnabin), type = "DNA")
fit_start   <- pml(nj_tree, data = phang_align)
fit_gtr     <- update(fit_start, k = 4, inv = 0.2)
fit_gtr     <- optim.pml(fit_gtr, model = "GTR", optInv = TRUE, optGamma = TRUE,
                          rearrangement = "stochastic", control = pml.control(trace = 0))

## 4d. Bootstrap support (1000 pseudo-trees, as described in the theory doc)
bootstrap_trees <- bootstrap.pml(fit_gtr, bs = 1000, optNni = TRUE,
                                  control = pml.control(trace = 0))
ml_tree_boot <- plotBS(fit_gtr$tree, bootstrap_trees, p = 50, type = "phylogram")

# Save the ML tree with bootstrap values
write.tree(ml_tree_boot, file = "data/passer_ML_tree_bootstrap.nwk")


## ---- 5. Genetic diversity statistics ---------------------------------------
# Convert alignment for pegas
aln_matrix <- as.matrix(aln_dnabin)

## Haplotype diversity (Hd)
haps       <- haplotype(aln_dnabin)
hap_div    <- hap.div(aln_dnabin)          # Hd (corrected for sample size)
cat("Haplotype diversity (Hd):", round(hap_div, 3), "\n")

## Nucleotide diversity (pi)
nuc_div    <- nuc.div(aln_dnabin, pairwise.deletion = TRUE)
cat("Nucleotide diversity (pi):", round(nuc_div, 5), "\n")

diversity_summary <- data.frame(
  Statistic = c("Haplotype diversity (Hd)", "Nucleotide diversity (pi)"),
  Value     = c(round(hap_div, 3), round(nuc_div, 5))
)
print(diversity_summary)
write.csv(diversity_summary, "data/diversity_summary.csv", row.names = FALSE)


## ---- 6. Tajima's D (test for population expansion/contraction) -----------
tajima_result <- tajima.test(aln_dnabin)
print(tajima_result)

cat("\nTajima's D =", round(tajima_result$D, 3),
    " | p-value (normal approx.) =", signif(tajima_result$Pval.normal, 3), "\n")

# Interpretation helper
interpret_tajima <- function(D, p) {
  if (is.na(p) || p >= 0.05) {
    "Not statistically significant: cannot reject neutral/stable population model."
  } else if (D < 0) {
    "Significantly negative D: excess of rare variants -> signature of recent population EXPANSION (or selective sweep)."
  } else {
    "Significantly positive D: excess of intermediate-frequency variants -> signature of population DECLINE / bottleneck (or balancing selection)."
  }
}
cat(interpret_tajima(tajima_result$D, tajima_result$Pval.normal), "\n")


## ---- 7. Visualization -------------------------------------------------------
dir.create("figures", showWarnings = FALSE)

## 7a. Phylogenetic tree
png("figures/phylo_tree.png", width = 1200, height = 900, res = 150)
plot.phylo(ml_tree_boot, main = "Maximum Likelihood Tree (Passer mtDNA)\nwith Bootstrap Support")
nodelabels(ml_tree_boot$node.label, adj = c(1.1, -0.3), frame = "none", cex = 0.7)
add.scale.bar()
dev.off()

## 7b. Genetic diversity bar chart (Hd vs pi, scaled for visibility)
div_plot_df <- data.frame(
  Statistic = c("Haplotype diversity (Hd)", "Nucleotide diversity (pi x100)"),
  Value     = c(hap_div, nuc_div * 100)
)

p_div <- ggplot(div_plot_df, aes(x = Statistic, y = Value, fill = Statistic)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  labs(title = "Genetic Diversity in Passer domesticus mtDNA",
       y = "Value", x = NULL) +
  theme_minimal(base_size = 13)
ggsave("figures/diversity_barplot.png", p_div, width = 6, height = 5, dpi = 150)

## 7c. Tajima's D visual summary
tajima_df <- data.frame(Statistic = "Tajima's D", D = tajima_result$D)

p_taj <- ggplot(tajima_df, aes(x = Statistic, y = D)) +
  geom_col(fill = ifelse(tajima_result$D < 0, "steelblue", "firebrick"), width = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Tajima's D for Passer domesticus mtDNA",
       subtitle = interpret_tajima(tajima_result$D, tajima_result$Pval.normal),
       y = "D value", x = NULL) +
  theme_minimal(base_size = 12)
ggsave("figures/tajimas_D_plot.png", p_taj, width = 6, height = 5, dpi = 150)


## ---- 8. Final summary printout ---------------------------------------------
cat("\n===================== ANALYSIS SUMMARY =====================\n")
cat("Sequences analyzed :", length(seqs), "\n")
cat("Alignment length    :", ncol(aln_matrix), "bp\n")
cat("Mean GC content      :", round(mean(gc_content), 2), "%\n")
cat("Haplotype diversity  :", round(hap_div, 3), "\n")
cat("Nucleotide diversity :", round(nuc_div, 5), "\n")
cat("Tajima's D            :", round(tajima_result$D, 3),
    "(p =", signif(tajima_result$Pval.normal, 3), ")\n")
cat("Outputs saved in ./data (tables, alignment, tree) and ./figures (plots)\n")
cat("==============================================================\n")
