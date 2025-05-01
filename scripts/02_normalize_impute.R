# 02 Normalize Impute.R

# 02 Normalize Impute.R

##########################################################
## Add batch info                                        #
##########################################################
od1$batch <- "d1"
od2$batch <- "d2"
##########################################################
## remove decoy                                          #
##########################################################
d1 <- od1[!od1$EG.IsDecoy, ]
d2 <- od2[!od2$EG.IsDecoy, ]
##########################################################
## Optional: Filter by Qvalue                            #
##########################################################
d1 <- d1[d1$PG.Qvalue < 0.01, ]
d2 <- d2[d2$PG.Qvalue < 0.01, ]
##########################################################
## remove contaminant                                    #
##########################################################
d1_conrev_remove <- d1[!grepl("^(CON__|REV__)", d1$PG.ProteinAccessions), ]
d2_conrev_remove <- d2[!grepl("^(CON__|REV__)", d2$PG.ProteinAccessions), ]
##########################################################
# combine two runs
##########################################################
data <- rbind(d1,d2)
##########################################################
# Log2-transform PG.Quantity & Rename Conditions
##########################################################
all_data <- data %>%
  mutate(log2_PG.Quantity = log2(PG.Quantity + 1)) %>%
  mutate(R.Condition = recode(R.Condition,
                              "CTRL"      = "CTRL",
                              "Treatment" = "KO_MyrPalmSterieCoA",
                              "Myr"       = "KO_Myr",
                              "Palm"      = "KO_Palm",
                              "Sterie"    = "KO_Sterie"
  )) %>%
  mutate(Sample_Name = paste0(R.Condition, "_", R.Replicate))  # Rename as Condition_Replicate
##########################################################
# Median normalization by sample (R.FileName)
##########################################################
all_data <- all_data %>%
  group_by(Sample_Name) %>%
  mutate(median_shift = median(log2_PG.Quantity, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(norm_log2_PG.Quantity = log2_PG.Quantity - median_shift)

# Check for missing values
table(is.na(all_data$norm_log2_PG.Quantity))

##########################################################
# Create metadata with renamed sample names
##########################################################
sample_info <- all_data %>%
  select(Sample_Name, R.Condition, R.Replicate, batch) %>%
  distinct()
write_csv(sample_info, "metadata.csv")
##########################################################
# Create wide format matrix using new sample names
##########################################################
wide_input <- all_data %>%
  group_by(PG.ProteinAccessions, Sample_Name) %>%
  summarize(norm_log2_quantity = median(norm_log2_PG.Quantity, na.rm = TRUE)) %>%
  ungroup()
# Pivot to wide format (proteins as rows, samples as columns)
wide_matrix <- wide_input %>%
  pivot_wider(names_from = Sample_Name, values_from = norm_log2_quantity) %>%
  column_to_rownames("PG.ProteinAccessions")
##########################################################
# determine the missingness of the data
# Calculate fraction of NAs per protein (row) and per sample (column)
protein_na_frac <- rowMeans(is.na(wide_matrix))
sample_na_frac <- colMeans(is.na(wide_matrix))
# Plot missingness per protein
protein_na_df <- tibble(missing_frac = protein_na_frac)
ggplot(protein_na_df, aes(x = missing_frac)) +
  geom_histogram(bins = 50, fill = "steelblue") +
  scale_y_break(c(2000, 3700)) +
  labs(title = "Missingness per Protein",
       x = "Fraction of Missing Values",
       y = "Number of Proteins") +
  theme_minimal()
ggsave("missingness_per_prot.png", width = 10, height = 5)
# Make a heatmap of missingness (1 = NA, 0 = present(blue))
missing_matrix <- is.na(wide_matrix) * 1
png("missing_matrix.png", width = 2000, height = 3000, res = 300)
pheatmap(missing_matrix, cluster_rows = TRUE, cluster_cols = TRUE,
         show_colnames = TRUE, show_rownames = FALSE)
dev.off()
# Plot missingness per sample
sample_na_df <- tibble(Sample = names(sample_na_frac),
                       missing_frac = sample_na_frac) 
sample_na_df <- sample_na_df %>%
  mutate(group = (as.integer(gl(n(), 5, n()))))

ggplot(sample_na_df, aes(x = Sample, y = missing_frac, fill = factor(group))) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Missingness per Sample",
       y = "Fraction of Missing Values") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, size = rel(2))) 
ggsave("missingness_per_sample.png", width = 12, height = 8)
## filter out protein that only exist in 1 reps
# Define your experimental groups based on sample names
groupings <- list(
  CTRL = grep("^CTRL", colnames(wide_matrix), value = TRUE),
  KO = grep("^KO_[^_]", colnames(wide_matrix), value = TRUE),  # KO_1, KO_2...
  KO_Myr = grep("^KO_Myr_", colnames(wide_matrix), value = TRUE),
  KO_Palm = grep("^KO_Palm_", colnames(wide_matrix), value = TRUE),
  KO_Sterie = grep("^KO_Sterie_", colnames(wide_matrix), value = TRUE),
  KO_MyrPalmSterieCoA = grep("^KO_MyrPalmSterieCoA_", colnames(wide_matrix), value = TRUE)
)

# Helper function to check if a row has data in ≥ 80% of group samples
keep_protein <- function(protein_row, group_samples) {
  observed <- sum(!is.na(protein_row[group_samples]))
  required <- ceiling(0.8 * length(group_samples))
  return(observed >= required)
}

# Apply the filter across all rows, keeping those that satisfy ≥ 80% in any one group
keep_rows <- apply(wide_matrix, 1, function(row) {
  any(sapply(groupings, function(samples) keep_protein(row, samples)))
})
# Filtered matrix
filtered_matrix <- wide_matrix[keep_rows, ]
# Turn into matrix 
exprs_mat <- as.matrix(filtered_matrix)
##########################################################
# impute
##########################################################
dim(exprs_mat)  # should be something like 3000 x 30
# Apply QRILC
set.seed(123)  # for reproducibility
imputed <- impute.QRILC(exprs_mat)
#Extract the imputed matrix
imputed_matrix <- as.data.frame(imputed[[1]])
imputed_matrix <- imputed_matrix %>%
  mutate(Protein = rownames(exprs_mat)) %>%
  select(Protein, everything())
write.xlsx(imputed_matrix, "imputed_matrix.xlsx")
############################################
# Normalize, impute, SVA + batch correct
# Output: combat_sva matrix with protein names
############################################
# -----------------------------
# Step 1: Load metadata and input matrix
# -----------------------------
sample_info <- read.csv("metadata.csv")
imputed_matrix <- read_xlsx("imputed_matrix.xlsx")

# Separate expression matrix and protein names
expr_mat <- imputed_matrix[,-1]
rownames(expr_mat) <- imputed_matrix[[1]]  # Assign protein names as rownames

# -----------------------------
# Step 2: Ensure column/sample order matches metadata
# -----------------------------
sample_info_ordered <- sample_info %>%
  filter(Sample_Name %in% colnames(expr_mat)) %>%
  arrange(match(Sample_Name, colnames(expr_mat)))

stopifnot(all(sample_info_ordered$Sample_Name == colnames(expr_mat)))  # Safety check

# -----------------------------
# Step 3: Subset for conditions of interest
# -----------------------------
conditions_of_interest <- c("CTRL", "KO", "KO_Sterie", "KO_MyrPalmSterieCoA")

sub_info <- sample_info_ordered %>%
  filter(R.Condition %in% conditions_of_interest)

expr_sub <- expr_mat[, sub_info$Sample_Name]

# -----------------------------
# Step 4: SVA + Batch Correction
# ----------------------------- 
mod  <- model.matrix(~ R.Condition, data = sub_info)
mod0 <- model.matrix(~ 1, data = sub_info)

# Estimate surrogate variables
svobj <- sva(as.matrix(expr_sub), mod, mod0)

# Apply batch + SVA correction
combat_sva <- removeBatchEffect(expr_sub,
                                covariates = svobj$sv,
                                design = mod)

# -----------------------------
# Step 5: Save final matrix with protein names
# -----------------------------
combat_sva_df <- as.data.frame(combat_sva) %>%
  mutate(Protein = rownames(expr_mat)) %>%
  select(Protein, everything())
write.xlsx(combat_sva_df, "Normalized_imputed_sva_batchCorrected_matrix.xlsx")
