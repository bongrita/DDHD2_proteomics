# 03 Batch Correction De.R

# 03 Batch Correction De.R

#######################################################################
#######################################################################
# Create a DE matrix
#######################################################################
#######################################################################
# ---------------------------------------------
# Load metadata and expression matrix
# ---------------------------------------------
sample_info <- read.csv("metadata.csv")
combat_sva <- read_xlsx("Normalized_imputed_sva_batchCorrected_matrix.xlsx")

# Extract expression matrix
expr_mat <- combat_sva[, -1]
rownames(expr_mat) <- combat_sva$Protein

# Ensure metadata alignment
sample_info_ordered <- sample_info %>%
  filter(Sample_Name %in% colnames(expr_mat)) %>%
  arrange(match(Sample_Name, colnames(expr_mat)))

stopifnot(all(sample_info_ordered$Sample_Name == colnames(expr_mat)))

# Subset to conditions of interest
conditions_of_interest <- c("CTRL", "KO", "KO_Sterie", "KO_MyrPalmSterieCoA")

sub_info <- sample_info_ordered %>%
  filter(R.Condition %in% conditions_of_interest)

expr_sub <- expr_mat[, sub_info$Sample_Name]

# ---------------------------------------------
# Create design and contrast matrix
# ---------------------------------------------
condition <- factor(sub_info$R.Condition, levels = conditions_of_interest)
design <- model.matrix(~ 0 + condition)
colnames(design) <- levels(condition)

contrast_matrix <- makeContrasts(
  KO_vs_CTRL = KO - CTRL,
  single_vs_CTRL = KO_Sterie - CTRL,
  triple_vs_CTRL = KO_MyrPalmSterieCoA - CTRL,
  levels = design
)

# ---------------------------------------------
# Run DE using limma and include protein names
# ---------------------------------------------
run_de <- function(matrix, protein_names) {
  fit <- lmFit(matrix, design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  
  list(
    KO = topTable(fit2, coef = "KO_vs_CTRL", number = Inf, adjust.method = "fdr") %>%
      rownames_to_column("Protein") %>%
      mutate(Protein = protein_names[match(Protein, rownames(matrix))]),
    
    STER = topTable(fit2, coef = "single_vs_CTRL", number = Inf, adjust.method = "fdr") %>%
      rownames_to_column("Protein") %>%
      mutate(Protein = protein_names[match(Protein, rownames(matrix))]),
    
    TRIPLE = topTable(fit2, coef = "triple_vs_CTRL", number = Inf, adjust.method = "fdr") %>%
      rownames_to_column("Protein") %>%
      mutate(Protein = protein_names[match(Protein, rownames(matrix))])
  )
}

de_corrected <- run_de(expr_sub, combat_sva$Protein)

# Combine and tag with bait
de_combined <- imap_dfr(de_corrected, ~ mutate(.x, Bait = .y))

write.xlsx(de_combined, "sva_batch_corrected_DE.xlsx")
