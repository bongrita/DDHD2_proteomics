# ================================
# Load or Install All Libraries
# ================================
packages <- c(
  "readxl","conflicted", "readr", "tidyr", "tibble", "pheatmap", "ggbreak",
  "imputeLCMD", "openxlsx", "limma", "writexl", "showtext", "jsonlite", "curl",
  "ggplot2", "scales", "ggrepel", "sva", "biomaRt", "org.Mm.eg.db", 
  "clusterProfiler", "GOSemSim", "rrvgo", "enrichR", "Cairo",
  "ComplexHeatmap", "circlize", "RColorBrewer", "forcats", "dplyr"
)

install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

invisible(lapply(packages, install_if_missing))
conflict_prefer("select", "dplyr")
conflict_prefer("rename", "dplyr")
##########################################################
## Load the files
##########################################################
list.files()
od1 <- read_tsv("20240719_SABER non normalised_Report_BGS Factory Report (Normal).tsv")
od2 <- read_tsv("20241115_SAber 1258 non-ormalised_Report_BGS Factory Report (Normal).tsv")
View(d1)

nrow(od1)
nrow(od2)
unique(od1$R.Condition)
unique(od2$R.Condition)
writeClipboard(capture.output(colnames(od1)))
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
nrow(d1)
nrow(d2)
##########################################################
## Optional: Filter by Qvalue                            #
##########################################################
d1 <- d1[d1$PG.Qvalue < 0.01, ]
d2 <- d2[d2$PG.Qvalue < 0.01, ]
nrow(d1)
nrow(d2)
unique(d1$EG.IsDecoy)
##########################################################
## remove contaminant                                    #
##########################################################
d1_conrev_remove <- d1[!grepl("^(CON__|REV__)", d1$PG.ProteinAccessions), ]
d2_conrev_remove <- d2[!grepl("^(CON__|REV__)", d2$PG.ProteinAccessions), ]
nrow(d1_conrev_remove)
nrow(d2_conrev_remove)
##########################################################
# combine two runs
##########################################################
data <- rbind(d1,d2)
nrow(data)
library(dplyr)
library(tidyr)
##########################################################
# Log2-transform PG.Quantity & Rename Conditions
##########################################################
colnames(data)
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
colnames(all_data)
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
writeClipboard(capture.output(print(sample_info, n = 30)))
write_csv(sample_info, "metadata.csv")
##########################################################
# Create wide format matrix using new sample names
##########################################################
View(all_data)
wide_input <- all_data %>%
  group_by(PG.ProteinAccessions, Sample_Name) %>%
  summarize(norm_log2_quantity = median(norm_log2_PG.Quantity, na.rm = TRUE)) %>%
  ungroup()
writeClipboard(capture.output(print(wide_input, n =20)))
# Pivot to wide format (proteins as rows, samples as columns)
library(tibble)
wide_matrix <- wide_input %>%
  pivot_wider(names_from = Sample_Name, values_from = norm_log2_quantity) %>%
  column_to_rownames("PG.ProteinAccessions")
writeClipboard(capture.output(print(wide_matrix[1:20,])))
View(wide_matrix)
##########################################################
# determine the missingness of the data
library(tidyverse)
library(pheatmap)
library(ggbreak)

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
View(wide_matrix)
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
library(dplyr)

# Define your experimental groups based on sample names
groupings <- list(
  CTRL = grep("^CTRL", colnames(wide_matrix), value = TRUE),
  KO = grep("^KO_[^_]", colnames(wide_matrix), value = TRUE),  # KO_1, KO_2...
  KO_Myr = grep("^KO_Myr_", colnames(wide_matrix), value = TRUE),
  KO_Palm = grep("^KO_Palm_", colnames(wide_matrix), value = TRUE),
  KO_Sterie = grep("^KO_Sterie_", colnames(wide_matrix), value = TRUE),
  KO_MyrPalmSterieCoA = grep("^KO_MyrPalmSterieCoA_", colnames(wide_matrix), value = TRUE)
)
print(sample_info, n= 30)

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
nrow(wide_matrix)
# Filtered matrix
filtered_matrix <- wide_matrix[keep_rows, ]
nrow(filtered_matrix)
# Turn into matrix 
exprs_mat <- as.matrix(filtered_matrix)
colnames(wide_matrix)
##########################################################
# impute
##########################################################
library(imputeLCMD)
dim(exprs_mat)  # should be something like 3000 x 30
# Apply QRILC
set.seed(123)  # for reproducibility
View(exprs_mat)
imputed <- impute.QRILC(exprs_mat)
#Extract the imputed matrix
imputed_matrix <- as.data.frame(imputed[[1]])
View(imputed_matrix)
imputed_matrix <- imputed_matrix %>%
  mutate(Protein = rownames(exprs_mat)) %>%
  select(Protein, everything())
library(openxlsx)
write.xlsx(imputed_matrix, "imputed_matrix.xlsx")
############################################
# Normalize, impute, SVA + batch correct
# Output: combat_sva matrix with protein names
############################################

library(tidyverse)
library(readxl)
library(openxlsx)
library(sva)
library(limma)

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
View(combat_sva_df)
write.xlsx(combat_sva_df, "Normalized_imputed_sva_batchCorrected_matrix.xlsx")

message("✅ SVA + batch-corrected matrix saved.")

#######################################################################
#######################################################################
# Create a DE matrix
#######################################################################
#######################################################################
library(tidyverse)
library(readxl)
library(limma)
library(ggplot2)
library(openxlsx)
# ---------------------------------------------
# Step 1: Load metadata and expression matrix
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
# Step 2: Create design and contrast matrix
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
# Step 3: Run DE using limma and include protein names
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

library(purrr)

# Combine and tag with bait
de_combined <- imap_dfr(de_corrected, ~ mutate(.x, Bait = .y))

write.xlsx(de_combined, "sva_batch_corrected_DE.xlsx")


de <- read.xlsx("sva_batch_corrected_DE.xlsx")
colnames(de)
View(de)
unique(de$Bait)
proteins <- de %>%
  filter(Bait %in% "TRIPLE",
         logFC < 0 & P.Value < 0.05)
nrow(proteins)

#######################################
# Make volcano plot usig DE dataframe #
#######################################
library(tidyverse)
library(readxl)
library(ggplot2)
library(showtext)

#---------------------------------#
# Font setup
#---------------------------------#
font_add("tnr", regular = "C:/Windows/Fonts/times.ttf")
showtext_auto()

#---------------------------------#
# Load DE result with Bait column
#---------------------------------#
View(df)
de_combined <- read_xlsx("sva_batch_corrected_DE.xlsx") %>%
  mutate(Bait = factor(Bait, levels = c("KO", "STER", "TRIPLE")))

# Split into a named list
de_corrected <- split(de_combined, de_combined$Bait)

#---------------------------------#
# Volcano plot function
#---------------------------------#
make_volcano <- function(df, title, outfile_base) {
  df$custom_sig <- ifelse(df$P.Value <= 0.05, "significant", "ns")
  
  p <- ggplot(df, aes(x = logFC, y = -log10(P.Value), color = custom_sig)) +
    geom_point(alpha = 0.8, size = 2) +
    scale_color_manual(values = c("ns" = "grey", "significant" = "forestgreen")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    theme_minimal(base_family = "tnr") +
    scale_x_continuous(breaks = seq(-15, 15, 5), limits = c(-15, 15)) +
    scale_y_continuous(breaks = seq(0, 25, 5), limits = c(0, 25)) +
    labs(
      title = title,
      x = "Log2FC",
      y = expression(-log[10]~"p-value"),
      color = "Significance"
    ) +
    theme(
      axis.text = element_text(size = rel(2.5)),
      axis.title = element_text(size = rel(2.5)),
      plot.title = element_text(hjust = 0.5, size = rel(2.4), face = "bold"),
      legend.position = "none",
      axis.line.x.bottom = element_line(color = "black"),
      axis.line.y.left = element_line(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    annotate("text", x = 10, y = -log10(0.05) + 1,
             label = "p < 0.05", size = rel(8), color = "forestgreen",
             family = "tnr") +
    annotate("text", x = 10, y = -log10(0.05) - 1,
             label = "p > 0.05", size = rel(8), color = "grey",
             family = "tnr")
  
  # Save to PDF and PNG
  ggsave(paste0(outfile_base, ".pdf"), plot = p, dpi = 300, width = 8, height = 8)
  ggsave(paste0(outfile_base, ".png"), plot = p, dpi = 100, width = 8, height = 8)
}

#---------------------------------#
# Generate Plots
#---------------------------------#
make_volcano(de_corrected$KO,
             bquote("DDHD2"^"-/-" ~"vs. C57BL6/J"),
             "ctrl_vs_KO_batchcorrected")

make_volcano(de_corrected$STER,
             bquote("DDHD2"^"-/-"~ "+ Stearic-CoA vs. C57BL6/J"),
             "ctrl_vs_sterie_batchcorrected")

make_volcano(de_corrected$TRIPLE,
             bquote("DDHD2"^"-/-" ~ "+ MyrPalmStearCoA vs. C57BL6/J"),
             "ctrl_vs_triple_batchcorrected")


###############################################
###############################################
# Fig 3C Rescue plot
###############################################
###############################################


#######################################################################################
# Figure 3C Rescue plot: Batch corrected rescue plot
#######################################################################################
library(tidyverse)
library(ggplot2)
library(showtext)

# Font (if not already set up)
font_add("tnr", regular = "C:/Windows/Fonts/times.ttf")
showtext_auto()

# Use de_corrected from previously loaded split list
ko_vs_ctrl     <- de_corrected$KO
single_vs_ctrl <- de_corrected$STER
triple_vs_ctrl <- de_corrected$TRIPLE

# Merge logFC and p-values across contrasts by Protein
df <- ko_vs_ctrl %>%
  select(Protein, logFC, P.Value) %>%
  rename(Log2FC_KO_CTRL = logFC, p_KO = P.Value) %>%
  left_join(single_vs_ctrl %>% select(Protein, logFC, P.Value),
            by = "Protein", suffix = c("", "_SterCoA")) %>%
  left_join(triple_vs_ctrl %>% select(Protein, logFC, P.Value),
            by = "Protein", suffix = c("", "_MyrPalmSterie")) %>%
  rename(Log2FC_KO_SterCoA_CTRL = logFC,
         p_SterCoA = P.Value,
         Log2FC_KO_MyrCoA_CTRL = logFC_MyrPalmSterie,
         p_MyrCoA = P.Value_MyrPalmSterie)

# Identify "rescued" proteins
df$Significance_Single <- "Not Significant"
df$Significance_Single[df$p_KO <= 0.05 & df$p_SterCoA > 0.05] <- "Rescued (Red)"
df$Significance_Single[df$p_KO <= 0.05 & df$p_SterCoA <= 0.05] <- "Still Dysregulated (Black)"

# Step 1: Number of significantly dysregulated proteins in KO vs WT
total_dysregulated_ko <- sum(df$p_KO <= 0.05)

# Step 2: Number of rescued proteins (significant in KO, but non-significant after S-CoA)
rescued_single <- sum(df$Significance_Single == "Rescued (Red)")

# Step 3: Rescue percentage
rescue_percentage_single <- (rescued_single / total_dysregulated_ko) * 100

# Print it
rescue_percentage_single





# Fit regression lines
lm_red   <- lm(Log2FC_KO_SterCoA_CTRL ~ Log2FC_KO_CTRL, data = df[df$Significance_Single == "Rescued (Red)", ])
lm_black <- lm(Log2FC_KO_SterCoA_CTRL ~ Log2FC_KO_CTRL, data = df[df$Significance_Single == "Still Dysregulated (Black)", ])
max(df$Log2FC_KO_CTRL)
min(df$Log2FC_KO_CTRL)
max(df$Log2FC_KO_SterCoA_CTRL)
min(df$Log2FC_KO_SterCoA_CTRL)
# Plot
p1 <- ggplot(df, aes(x = Log2FC_KO_CTRL, y = Log2FC_KO_SterCoA_CTRL)) +
  geom_abline(intercept = coef(lm_black)[1], slope = coef(lm_black)[2], color = "black") +
  geom_abline(intercept = coef(lm_red)[1], slope = coef(lm_red)[2], color = "red") +
  geom_point(data = df[df$Significance_Single == "Not Significant", ], color = "grey", alpha = 0.5) +
  geom_point(data = df[df$Significance_Single == "Still Dysregulated (Black)", ], color = "black", alpha = 0.7) +
  geom_point(data = df[df$Significance_Single == "Rescued (Red)", ], color = "red", alpha = 0.7) +
  # Axes & style
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal(base_family = "tnr") +
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(2.5), face = "bold"),
    axis.title.x = element_text(size = rel(3)),
    axis.text.x = element_text(size = rel(3.5)),
    axis.title.y = element_text(size = rel(2.6)),
    axis.text.y = element_text(size = rel(3)),
    axis.line.x.bottom = element_line(color = "black"),
    axis.line.y.left = element_line(color = "black"),
    panel.grid = element_blank()
  ) +
  labs(
    x = expression(Log[2]~" Ratio (DDHD2"^"-/-"~"vs C57BL6/J)"),
    y = expression(Log[2]~" Ratio (DDHD2"^"-/-"~"+ Stearic-CoA vs C57BL6/J)"),
    title = "Stearic-CoA Rescue Plot"
  ) +
  scale_x_continuous(breaks = seq(-15, 15, 5), limits = c(-15, 15)) +
  scale_y_continuous(breaks = seq(-15, 15, 5), limits = c(-15, 15)) +
  annotate("text", x = 0.5, y = -8,
           label = expression("DDHD2"^"-/-" ~"+"),
           color = "red", size = rel(10), hjust = 0,
           family = "tnr") +
  annotate("text", x = 0.5, y = -10,
           label = expression("Stearic-CoA vs C57BL6/J"),
           color = "red", size = rel(10), hjust = 0,
           family = "tnr") +
  annotate("text", x = 0.5, y = -12,
           label = expression("p > 0.05"),
           color = "red", size = rel(10), hjust = 0,
           family = "tnr") +
  annotate("text", x = -10, y = 12,
           label = expression("DDHD2"^"-/-" ~" vs"),
           color = "black", size = rel(10), hjust = 0,
           family = "tnr") +
  annotate("text", x = -10, y = 10,
           label = "C57BL6/J",
           color = "black", size = rel(10), hjust = 0,
           family = "tnr") +
  annotate("text", x = -10, y = 8,
           label = "p < 0.05",
           color = "black", size = rel(10), hjust = 0,
           family = "tnr")

p1
# Save
ggsave("Rescue_batch_corrected_StearicCoA.pdf", p1, width = 10, height = 10, dpi = 300)
ggsave("Rescue_batch_corrected_StearicCoA.png", p1, width = 10, height = 10, dpi = 100)

#######################################
# Triple rescue
######################################
# Define significance categories
df$Significance_Triple <- "Not Significant"
df$Significance_Triple[df$p_KO <= 0.05 & df$p_MyrCoA > 0.05] <- "Rescued (Red)"
df$Significance_Triple[df$p_KO <= 0.05 & df$p_MyrCoA <= 0.05] <- "Still Dysregulated (Black)"


# Step 1: Number of significantly dysregulated proteins in KO vs WT
total_dysregulated_ko <- sum(df$p_KO <= 0.05)

# Step 2: Number of rescued proteins (significant in KO, but non-significant after S-CoA)
rescued_triple <- sum(df$Significance_Triple == "Rescued (Red)")

# Step 3: Rescue percentage
rescue_percentage_triple <- (rescued_triple / total_dysregulated_ko) * 100

# Print it
rescue_percentage_triple




# Fit regression lines
lm_red   <- lm(Log2FC_KO_MyrCoA_CTRL ~ Log2FC_KO_CTRL, data = df[df$Significance_Triple == "Rescued (Red)", ])
lm_black <- lm(Log2FC_KO_MyrCoA_CTRL ~ Log2FC_KO_CTRL, data = df[df$Significance_Triple == "Still Dysregulated (Black)", ])

# Plot
p2 <- ggplot(df, aes(x = Log2FC_KO_CTRL, y = Log2FC_KO_MyrCoA_CTRL)) +
  geom_abline(intercept = coef(lm_black)[1], slope = coef(lm_black)[2], color = "black") +
  geom_abline(intercept = coef(lm_red)[1], slope = coef(lm_red)[2], color = "red") +
  geom_point(data = df[df$Significance_Triple == "Not Significant", ], color = "grey", alpha = 0.5) +
  geom_point(data = df[df$Significance_Triple == "Still Dysregulated (Black)", ], color = "black", alpha = 0.7) +
  geom_point(data = df[df$Significance_Triple == "Rescued (Red)", ], color = "red", alpha = 0.7) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  
  theme_minimal(base_family = "tnr") +
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(3), face = "bold"),
    axis.title.x = element_text(size = rel(3)),
    axis.text.x = element_text(size = rel(3)),
    axis.title.y = element_text(size = rel(2.4)),
    axis.text.y = element_text(size = rel(3)),
    axis.line.x.bottom = element_line(color = "black"),
    axis.line.y.left = element_line(color = "black"),
    panel.grid = element_blank()
  ) +
  labs(
    x = expression(Log[2]~" Ratio (DDHD2"^"-/-"~"vs C57BL6/J)"),
    y = expression(Log[2]~" Ratio (DDHD2"^"-/-"~"+ MyrPalmStear-CoA vs C57BL6/J)"),
    title = "Triple-CoA Rescue Plot"
  ) +
  scale_x_continuous(breaks = seq(-15, 15, 5), limits = c(-15, 15)) +
  scale_y_continuous(breaks = seq(-15, 15, 5), limits = c(-15, 15)) +
  
  # Red annotation (rescued)
  annotate("text", x = 1, y = -5,
           label = expression("DDHD2"^"-/-" ~"+"),
           color = "red", size = rel(10), hjust = 0, family = "tnr") +
  annotate("text", x = 1, y = -7,
           label = expression("MyrPalmStear-CoA vs"),
           color = "red", size = rel(10), hjust = 0, family = "tnr") +
  annotate("text", x = 1, y = -9,
           label = expression("C57BL6/J"),
           color = "red", size = rel(10), hjust = 0, family = "tnr") +
  annotate("text", x = 1, y = -11,
           label = expression("p > 0.05"),
           color = "red", size = rel(10), hjust = 0, family = "tnr") +
  
  # Black annotation (still dysregulated)
  annotate("text", x = -10, y = 8,
           label = expression("DDHD2"^"-/-" ~" vs"),
           color = "black", size = rel(10), hjust = 0, family = "tnr") +
  annotate("text", x = -10, y = 6,
           label = "C57BL6/J",
           color = "black", size = rel(10), hjust = 0, family = "tnr") +
  annotate("text", x = -10, y = 4,
           label = "p < 0.05",
           color = "black", size = rel(10), hjust = 0, family = "tnr")

p2

# Save plot
ggsave("Rescue_batch_corrected_MyrPalmStearCoA.pdf", p2, width = 10, height = 10, dpi = 300)
ggsave("Rescue_batch_corrected_MyrPalmStearCoA.png", p2, width = 10, height = 10, dpi = 100)

#####################
# PCA ###############
#####################
sample_info <- read.csv("metadata.csv")

# Ensure metadata alignment
sample_info_ordered <- sample_info %>%
  filter(Sample_Name %in% colnames(expr_mat)) %>%
  arrange(match(Sample_Name, colnames(expr_mat)))

# Subset to conditions of interest
conditions_of_interest <- c("CTRL", "KO", "KO_Sterie", "KO_MyrPalmSterieCoA")

sub_info <- sample_info_ordered %>%
  filter(R.Condition %in% conditions_of_interest)

combat_sva <- read.xlsx("Normalized_imputed_sva_batchCorrected_matrix.xlsx")
combat <- combat_sva[,-1]
pca <- prcomp(t(combat), scale. = TRUE)

# Prepare data for plotting
pc_df <- as.data.frame(pca$x[, 1:2]) %>%
  mutate(Sample = rownames(.)) %>%
  left_join(sub_info, by = c("Sample" = "Sample_Name")) %>%
  mutate(
    Condition = factor(R.Condition,
                       levels = c("CTRL", "KO", "KO_MyrPalmSterieCoA", "KO_Sterie"),
                       labels = c("C57BL6/J", "DDHD2⁻/⁻", "DDHD2⁻/⁻ + MyrPalmStear-CoA", "DDHD2⁻/⁻ + Stearic-CoA")),
    Rep_Factor = factor(R.Replicate)
  )

pc_df$Label <- recode(pc_df$Sample,
                      "CTRL_1" = "C57BL6/J_1",
                      "CTRL_2" = "C57BL6/J_2",
                      "CTRL_3" = "C57BL6/J_3",
                      "CTRL_4" = "C57BL6/J_4",
                      "CTRL_5" = "C57BL6/J_5",
                      # 
                      "KO_1" = "DDHD2^'-/-'~'_'*1",
                      "KO_2" = "DDHD2^'-/-'~'_'*2",
                      "KO_3" = "DDHD2^'-/-'~'_'*3",
                      "KO_4" = "DDHD2^'-/-'~'_'*4",
                      "KO_5" = "DDHD2^'-/-'~'_'*5",
                      
                      "KO_MyrPalmSterieCoA_1" = "DDHD2^'-/-' + MyrPalmStear-CoA_1",
                      "KO_MyrPalmSterieCoA_2" = "DDHD2^'-/-' + MyrPalmStear-CoA_2",
                      "KO_MyrPalmSterieCoA_3" = "DDHD2^'-/-' + MyrPalmStear-CoA_3",
                      "KO_MyrPalmSterieCoA_4" = "DDHD2^'-/-' + MyrPalmStear-CoA_4",
                      "KO_MyrPalmSterieCoA_5" = "DDHD2^'-/-' + MyrPalmStear-CoA_5",
                      # 
                      "KO_Sterie_1" = "DDHD2^'-/-' + Stearic-CoA_1",
                      "KO_Sterie_2" = "DDHD2^'-/-' + Stearic-CoA_2",
                      "KO_Sterie_3" = "DDHD2^'-/-' + Stearic-CoA_3",
                      "KO_Sterie_4" = "DDHD2^'-/-' + Stearic-CoA_4",
                      "KO_Sterie_5" = "DDHD2^'-/-' + Stearic-CoA_5"
)


set.seed(123)  # for reproducibility
jitter_width <- 7
jitter_height <- 7

pc_df <- pc_df %>%
  mutate(
    jittered_PC1 = PC1 + runif(n(), -jitter_width, jitter_width),
    jittered_PC2 = PC2 + runif(n(), -jitter_height, jitter_height)
  )


p <- ggplot(pc_df, aes(x = jittered_PC1, y = jittered_PC2)) +
  # Outer square for replicate ID
  geom_point(aes(color = Rep_Factor), shape = 0, size = 2.6, stroke = 1, fill = NA) +
  
  # Inner circle for condition
  geom_point(aes(fill = Condition), shape = 21, size = 1.8, color = "black", stroke = 0.5) +
  
  
  # # Sample labels
  # ggrepel::geom_text_repel(
  #   aes(label = Label),
  #   parse = TRUE,
  #   size = 3.5,
  #   force = 2,                # Stronger repulsion
  #   force_pull = 0.5,           # Help pull away overlapping text
  #   max.overlaps = 100,
  #   box.padding = 0.4,
  #   point.padding = 0.5,
  #   segment.size = 0.1,
  #   min.segment.length = 0.05,
  #   segment.color = "grey40"
  # ) +
  # 
  labs(x = "PC1", y = "PC2", title = NULL) +
  theme_minimal(base_size = 14, base_family = "tnr") +
  scale_fill_manual(
    values = c("black", "red", "green", "blue"),
    labels = c(
      "C57BL6/J",
      expression("DDHD2"^"-/-"),
      expression("DDHD2"^"-/-"~"+ MyrPalmStear-CoA"),
      expression("DDHD2"^"-/-"~"+ Stearic-CoA")
    ),
    guide = guide_legend(override.aes = list(shape = 21, size = 5)),
    drop = FALSE
  ) +
  scale_color_manual(values = c("black", "hotpink", "orange", "purple", "darkblue"),
                     labels = c("Repetition 1",
                                "Repetition 2",
                                "Repetition 3",
                                "Repetition 4",
                                "Repetition 5")) +
  scale_x_continuous(breaks = seq(-100, 100, by = 20),)+
  scale_y_continuous(breaks = seq(-60, 60, by = 20)) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank(),
    axis.title = element_text(size = rel(1.8)),
    axis.text = element_text(size = rel(1.8)),
#    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.box = "vertical",
    legend.text = element_text(size = rel(1.8)),
    legend.background = element_blank(),
    legend.box.background = element_rect(color = "black", size = 0.2),
    legend.box.margin = margin(5, 5, 5, 5),
    plot.margin = margin(20,50,20,30), # trbl
    legend.position = c(0.65, 0.5)
  ) 
p
ggsave(plot = p, "PCA_BatchCorrected_noLabel.pdf", width = 8, height = 8, dpi = 300)
ggsave(plot = p, "PCA_BatchCorrected_noLabel.png", width = 8, height = 8, dpi = 100)

##################
# With label
##################
p <- ggplot(pc_df, aes(x = jittered_PC1, y = jittered_PC2)) +
  # Outer square for replicate ID
  geom_point(aes(color = Rep_Factor), shape = 0, size = 2.6, stroke = 1, fill = NA) +
  
  # Inner circle for condition
  geom_point(aes(fill = Condition), shape = 21, size = 1.8, color = "black", stroke = 0.5) +


  # Sample labels
  ggrepel::geom_text_repel(
    aes(label = Label),
    parse = TRUE,
    size = 3.5,
    force = 2,                # Stronger repulsion
    force_pull = 0.1,           # Help pull away overlapping text
    max.overlaps = 100,
    box.padding = 0.1,
    point.padding = 0.5,
    segment.size = 0.01,
    min.segment.length = 0.001,
    segment.color = "grey40"
  ) +

  labs(x = "PC1", y = "PC2", title = NULL) +
  theme_minimal(base_size = 14, base_family = "tnr") +
  scale_fill_manual(
    values = c("black", "red", "green", "blue"),
    labels = c(
      "C57BL6/J",
      expression("DDHD2"^"-/-"),
      expression("DDHD2"^"-/-"~"+ MyrPalmStear-CoA"),
      expression("DDHD2"^"-/-"~"+ Stearic-CoA")
    ),
    guide = guide_legend(override.aes = list(shape = 21, size = 5)),
    drop = FALSE
  ) +
  scale_color_manual(values = c("black", "hotpink", "orange", "purple", "darkblue"),
                     labels = c("Repetition 1",
                                "Repetition 2",
                                "Repetition 3",
                                "Repetition 4",
                                "Repetition 5")) +
  scale_x_continuous(breaks = seq(-100, 100, by = 20),)+
  scale_y_continuous(breaks = seq(-60, 60, by = 20)) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank(),
    axis.title = element_text(size = rel(1.5)),
    axis.text = element_text(size = rel(1.5)),
    legend.box = "vertical",
    legend.text = element_text(size = rel(0.9)),
    legend.key.size = unit(0.4, "cm"),
    legend.spacing = unit(0.2, "cm"),
    legend.background = element_blank(),
    legend.box.background = element_rect(color = "black", size = 0.2),
    legend.box.margin = margin(3, 3, 3, 3),
    plot.margin = margin(20,50,20,30),
    legend.position = c(0.8, 0.5)
  ) 
p
ggsave(plot = p, "PCA_BatchCorrected_Labeled.pdf", width = 8, height = 8, dpi = 300)
ggsave(plot = p, "PCA_BatchCorrected_Labeled.png", width = 8, height = 8, dpi = 100)


#######################################
## Hierarchical clustering         ####
#######################################

# ========================================
# Hierarchical Clustering (SVA + batch corrected)
# ========================================

library(readxl)

# Load data
combat_sva <- read_xlsx("Normalized_imputed_sva_batchCorrected_matrix.xlsx")

# Prepare expression matrix (remove protein column)
expr_sub <- combat_sva[, -1]
rownames(expr_sub) <- combat_sva$Protein  # Optional, but good to retain

# Transpose for clustering
exprs_t <- t(expr_sub)

# Distance matrix and clustering
dist_mat <- dist(exprs_t, method = "euclidean")
hc <- hclust(dist_mat, method = "ward.D2")

# Compress tree height for better layout
hc$height <- hc$height / 5

# Plot to PDF
pdf("Hierarchical_clustering_clean.pdf", width = 10, height = 10)
par(mar = c(6, 6, 4, 3))  # Bottom left top right
plot(hc,
     main = "Hierarchical Clustering",
     xlab = "",
     sub = "",
     cex = 2,
     cex.main = 3,
     cex.axis = 2,
     cex.lab = 2,
     las = 2)  # Rotate axis labels vertically
dev.off()

png("Hierarchical_clustering_clean.png", width = 2000, height = 2000, res = 100)
par(mar = c(6, 6, 6, 3))  # Bottom left top right
plot(hc,
     main = "Hierarchical Clustering",
     xlab = "",
     sub = "",
     cex = 4,
     cex.main = 8,
     cex.axis = 4,
     cex.lab = 4,
     las = 2)  # Rotate axis labels vertically
dev.off()


##############################
## GO
##############################
# Define the list of ontologies, this case run all BP, CC , MF

df <- read.xlsx("sva_batch_corrected_DE.xlsx")
nrow(df)

pval <- 0.05
fc <- 1

# Keep only the first UniProt accession before the first semicolon
df$UniProt <- sub(";.*", "", df$Protein)
df <- df %>%
  mutate(Row_ID = row_number())
unique(df$UniProt)
library(biomaRt)
library(org.Mm.eg.db)
library(clusterProfiler)
hci_bitr <- bitr(unique(df$UniProt), fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# some uniprot map to multiple entrez
multiple_mappings <- hci_bitr %>%
  group_by(UNIPROT) %>%
  filter(n() > 1) %>%
  arrange(UNIPROT)  # Sort for better readability

# Print them
print(multiple_mappings, n = nrow(multiple_mappings))
# you can only keep one 
hci_bitr_unique <- hci_bitr %>%
  group_by(UNIPROT) %>%
  dplyr::slice(1)  # Keep only the first match

# Now, join the Entrez IDs back to the original 'hci' data
enrich_mapped <- df %>%
  left_join(hci_bitr_unique, by = c("UniProt" = "UNIPROT"))
colnames(enrich_mapped)
enrich_final <- enrich_mapped %>%
  group_by(Row_ID) %>%  # Group by the original row ID
  summarise(
    Bait = dplyr::first(Bait),
    uniprot = dplyr::first(UniProt),  # Keep the original UniProt ID
    p_value = dplyr::first(P.Value),
    adj_p_value = dplyr::first(adj.P.Val),
    logFC = dplyr::first(logFC),
    ENTREZID = dplyr::first(ENTREZID),  # Keep only the first mapped ENTREZ ID
    .groups = "drop"
  ) %>%
  dplyr::select(-Row_ID) %>%
  mutate(group = case_when(p_value <= pval & logFC >= fc ~ "Upregulated",
                           p_value <= pval & logFC <= -fc ~ "Downregulated",
                           TRUE~ "Non-significant"))

# Map ENTREZID to gene symbol (SYMBOL)
gene_symbol_map <- bitr(unique(enrich_final$ENTREZID),
                        fromType = "ENTREZID",
                        toType = "SYMBOL",
                        OrgDb = org.Mm.eg.db)

# Join the SYMBOL back into the final table
enrich_final <- enrich_final %>%
  left_join(gene_symbol_map, by = "ENTREZID") %>%
  relocate(SYMBOL, .after = uniprot)  # Optional: move gene symbol next to uniprot
View(enrich_final)
library(GOSemSim)
library(rrvgo)
library(enrichR)
unique(enrich_final$group)

ontologies <- c("BP")
all_results <- data.frame()  # Initialize an empty data frame to store results
baits <- unique(enrich_final$Bait)
significance <- c("Upregulated", "Downregulated")

summary <- enrich_final %>%
  group_by(Bait, group) %>%
  summarize(count = n(), .groups = "drop")

# Loop through each ontology
for (ontology in ontologies) {
  # Loop through each bait
  for (i in 1:length(baits)) {
    # Loop through each significance group
    for (sig in significance) {
      
      bait_bp <- enrich_final %>%
        filter(Bait == baits[i],
               group == sig)
      
      # Skip if there are too few genes
      if (nrow(bait_bp) < 2) {
        cat("Too few genes for", baits[i], sig, ontology, "- skipping\n")
        next
      }
      
      # Perform GO enrichment analysis
      go_result <- enrichGO(
        gene = unique(bait_bp$ENTREZID), 
        OrgDb = org.Mm.eg.db, 
        ont = ontology, 
        pvalueCutoff = 0.05
      )
      
      go_df <- as.data.frame(go_result)
      if (nrow(go_df) < 2) {
        cat("Too few GO terms for", baits[i], sig, ontology, "- skipping\n")
        next
      }
      
      go_ids <- go_df$ID
      go_scores <- setNames(-log10(go_df$p.adjust), go_ids)
      
      # Wrap simMatrix + reduceSimMatrix in error-safe block
      tryCatch({
        sim_matrix <- calculateSimMatrix(
          go_ids,
          orgdb = 'org.Mm.eg.db',
          ont = ontology,
          method = "Rel"
        )
        
        if (nrow(sim_matrix) < 2 || ncol(sim_matrix) < 2) {
          cat("Similarity matrix too small for", baits[i], sig, ontology, "- skipping\n")
          next
        }
        
        reduced_go <- reduceSimMatrix(
          sim_matrix,
          scores = go_scores,
          threshold = 0.8,
          orgdb = 'org.Mm.eg.db'
        )
        
        go_reduced_df <- as.data.frame(reduced_go)
        go_bp <- go_df %>%
          left_join(go_reduced_df, by = c("ID" = "go"))
        
        go_bp$Bait <- baits[i]
        go_bp$Ont <- ontology
        go_bp$group <- sig
        
        all_results <- bind_rows(all_results, go_bp)
        
      }, error = function(e) {
        cat("Error during rrvgo reduction for", baits[i], sig, ontology, ":", conditionMessage(e), "\n")
      })
    }
  }
}

# Ensure parentTerm is filled
all_results_final <- all_results %>%
  dplyr::mutate(parentTerm = ifelse(is.na(parentTerm), Description, parentTerm))

# Calculate ratios
all_results_save <- all_results_final %>%
  mutate(
    GeneRatio_val = as.numeric(sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))),
    BgRatio_val = as.numeric(sapply(strsplit(BgRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))),
    ratio = GeneRatio_val / BgRatio_val
  )

# Save to CSV
write.csv(all_results_save, file = "20250403_DDHD2_BP_reduced_thres_08.csv", row.names = FALSE)



##############################################################
# recreate heatmap figures
##############################################################
# ------------------ Load data ------------------
sample_info <- read.csv("metadata.csv")
imputed_matrix <- read.xlsx("imputed_matrix.xlsx")
sva_mat <- read.xlsx("Normalized_imputed_sva_batchCorrected_matrix.xlsx")
list.files()

# ------------------ Separate expression and annotation ------------------
expr_mat <- sva_mat[,-1]               # Only expression matrix
df <- cbind(imputed_matrix[,1], expr_mat)
colnames(df)[1] <- "Protein"
# Keep only the first UniProt accession before the first semicolon
df$Protein <- sub(";.*", "", df$Protein)
df <- df %>%
  mutate(Row_ID = row_number())
hci_bitr <- bitr(unique(df$Protein), fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
# some uniprot map to multiple entrez
multiple_mappings <- hci_bitr %>%
  group_by(UNIPROT) %>%
  filter(n() > 1) %>%
  arrange(UNIPROT)  # Sort for better readability

# Print them
print(multiple_mappings, n = nrow(multiple_mappings))
# you can only keep one 
hci_bitr_unique <- hci_bitr %>%
  group_by(UNIPROT) %>%
  dplyr::slice(1)  # Keep only the first match

# Now, join the Entrez IDs back to the original 'hci' data
enrich_mapped <- df %>%
  left_join(hci_bitr_unique, by = c("Protein" = "UNIPROT"))

gene_symbol_map <- bitr(unique(enrich_mapped$ENTREZID),
                        fromType = "ENTREZID",
                        toType = "SYMBOL",
                        OrgDb = org.Mm.eg.db)

enrich_final <- enrich_mapped %>%
  left_join(gene_symbol_map, by = c("ENTREZID" = "ENTREZID"))
#####################################################################
# Now get the gene ID for presynapse
#####################################################################
cc <- read.xlsx("20250404_DDHD2_enrichGO_reduced_thres_08.xlsx")
list.files()
# presynapse
presynapse <- c("GO:0042734", "GO:0098830", "GO:0098833", "GO:0043679", "GO:0043195", "GO:0099569","GO:0099523","GO:0008021", "GO:0099143", "GO:0099182", "GO:0048786")
id <- cc %>%
  filter(ID %in% presynapse)
nrow(id)
allgeneid <- unlist(strsplit(id$geneID, split = "/"))
unique_geneid <- unique(allgeneid)

int <- enrich_final %>%
  filter(ENTREZID %in% unique_geneid)
nrow(int)

str(int)
### pheatmap 
# Required packages
labels_expr <- c(
  expression("C57BL/6J"~"_1"),
  expression("C57BL/6J"~"_2"),
  expression("C57BL/6J"~"_3"),
  expression("C57BL/6J"~"_4"),
  expression("C57BL/6J"~"_5"),
  
  
  expression("DDHD2"^"-/-"~"+ MyrPalmStearCoA_1"),
  expression("DDHD2"^"-/-"~"+ MyrPalmStearCoA_2"),
  expression("DDHD2"^"-/-"~"+ MyrPalmStearCoA_3"),
  expression("DDHD2"^"-/-"~"+ MyrPalmStearCoA_4"),
  expression("DDHD2"^"-/-"~"+ MyrPalmStearCoA_5"),
  
  expression("DDHD2"^"-/-"~"_1"),
  expression("DDHD2"^"-/-"~"_2"),
  expression("DDHD2"^"-/-"~"_3"),
  expression("DDHD2"^"-/-"~"_4"),
  expression("DDHD2"^"-/-"~"_5"),
  
  
  expression("DDHD2"^"-/-"~"+ Stearic-CoA_1"),
  expression("DDHD2"^"-/-"~"+ Stearic-CoA_2"),
  expression("DDHD2"^"-/-"~"+ Stearic-CoA_3"),
  expression("DDHD2"^"-/-"~"+ Stearic-CoA_4"),
  expression("DDHD2"^"-/-"~"+ Stearic-CoA_5")
)

colnames(int)
# Subset intensity data (remove annotation columns)
intensity_data <- int %>%
  dplyr::select(-Protein, -Row_ID, -ENTREZID, -SYMBOL)

# Convert to matrix
intensity_matrix <- as.matrix(intensity_data)

# Row-wise z-score (center and scale by row)
z_score_matrix <- t(scale(t(intensity_matrix)))

# Optional: set protein names as rownames
rownames(z_score_matrix) <- int$SYMBOL  # or use $Protein or $ENTREZID




# Add TNR from system fonts (adjust the path if needed)
font_add("tnr", regular = "C:/Windows/Fonts/times.ttf")
showtext_auto()





# Z-score color breaks

Cairo::CairoPDF("Heatmap_Presynapse_heatmap.pdf", width = 8, height = 12, family = "Times")
Cairo::CairoPNG("Heatmap_Presynapse_heatmap.png", width = 800, height = 1100, family = "Times", res = 100)
library(pheatmap)
# Plot
?pheatmap()
pheatmap(z_score_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete", 
         show_rownames = FALSE,
         legend = FALSE,
         fontsize = 20,
         cellheight = 0.8,
         cellwidth = 20,
         legend_breaks = c(-2, 0, 2),
         breaks = seq(-2, 2, length.out = 100),
         fontsize_col = 20,
         #angle_col = 90,
         main = "Presynapse",
         #border_color = "black",
         labels_col = labels_expr  # 👈 parsed labels here
)

dev.off()


#####################################################################
#####################################################################
# Now get the gene ID for ENDOPLASMIC RETICULUM
#####################################################################
#####################################################################
cc <- read.csv("20250327_DDHD2_reduced_thres_08.csv")
list.files()
# presynapse
er <- c("GO:0009511", "GO:0016529", "GO:0005788", "GO:0009510",
        "GO:0044322", "GO:0000835", "GO:0005790", "GO:0098827",
        "GO:0005791", "GO:0070971", "GO:0097038", "GO:0140534", "GO:1990007")
id <- cc %>%
  filter(ID %in% er)
nrow(id)
allgeneid <- unlist(strsplit(id$geneID, split = "/"))
unique_geneid <- unique(allgeneid)

int <- enrich_final %>%
  filter(ENTREZID %in% unique_geneid)
nrow(int)

str(int)
### pheatmap 
# Required packages
library(pheatmap)
library(dplyr)
library(pheatmap)

library(tidyverse)

labels_expr <- c(
  expression("C57BL/6J"~"_1"),
  expression("C57BL/6J"~"_2"),
  expression("C57BL/6J"~"_3"),
  expression("C57BL/6J"~"_4"),
  expression("C57BL/6J"~"_5"),
  
  
  expression("DDHD2"^"-/-"~"+ MyrPalmStearCoA_1"),
  expression("DDHD2"^"-/-"~"+ MyrPalmStearCoA_2"),
  expression("DDHD2"^"-/-"~"+ MyrPalmStearCoA_3"),
  expression("DDHD2"^"-/-"~"+ MyrPalmStearCoA_4"),
  expression("DDHD2"^"-/-"~"+ MyrPalmStearCoA_5"),
  
  expression("DDHD2"^"-/-"~"_1"),
  expression("DDHD2"^"-/-"~"_2"),
  expression("DDHD2"^"-/-"~"_3"),
  expression("DDHD2"^"-/-"~"_4"),
  expression("DDHD2"^"-/-"~"_5"),
  
  
  expression("DDHD2"^"-/-"~"+ Stearic-CoA_1"),
  expression("DDHD2"^"-/-"~"+ Stearic-CoA_2"),
  expression("DDHD2"^"-/-"~"+ Stearic-CoA_3"),
  expression("DDHD2"^"-/-"~"+ Stearic-CoA_4"),
  expression("DDHD2"^"-/-"~"+ Stearic-CoA_5")
)

colnames(int)
# Subset intensity data (remove annotation columns)
intensity_data <- int %>%
  dplyr::select(-Protein, -Row_ID, -ENTREZID, -SYMBOL)

# Convert to matrix
intensity_matrix <- as.matrix(intensity_data)

# Row-wise z-score (center and scale by row)
z_score_matrix <- t(scale(t(intensity_matrix)))

# Optional: set protein names as rownames
rownames(z_score_matrix) <- int$SYMBOL  # or use $Protein or $ENTREZID
writeClipboard(capture.output(rownames(z_score_matrix)))

library(showtext)

# Add TNR from system fonts (adjust the path if needed)
font_add("tnr", regular = "C:/Windows/Fonts/times.ttf")
showtext_auto()





# Z-score color breaks
breaks <- seq(-2, 2, length.out = 200)
Cairo::CairoPDF("Heatmap_ER_heatmap.pdf", width = 8, height = 12, family = "Times")
Cairo::CairoPNG("Heatmap_ER_heatmap.png", width = 800, height = 1200, family = "Times", res = 100)


# Plot
pheatmap(z_score_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = FALSE,
         legend = FALSE,
         fontsize = 20,
         cellheight = 1.8,
         cellwidth = 20,
         legend_breaks = c(-2, 0, 2),
         breaks = seq(-2, 2, length.out = 100),
         fontsize_col = 20,
         main = "Endoplasmic Reticulum",
         #border_color = "black",
         labels_col = labels_expr  # 👈 parsed labels here
)

dev.off()

#####################################################################
#####################################################################
# Now get the gene ID for Golgi apparatus GO:0005794
#####################################################################
#####################################################################
cc <- read.csv("20250327_DDHD2_reduced_thres_08.csv")
list.files()
# presynapse
golgi <- c("GO:0005796", "GO:1990071", "GO:0000938", "GO:1990072", "GO:0000139", "GO:0070916", "GO:0005801",
           "GO:0070931", "GO:0150051", "GO:0030130", "GO:0098791", "GO:0034044", "GO:0036063","GO:0030126", "GO:0034066", "GO:0017119")
id <- cc %>%
  filter(ID %in% golgi)
nrow(id)
allgeneid <- unlist(strsplit(id$geneID, split = "/"))
unique_geneid <- unique(allgeneid)

int <- enrich_final %>%
  filter(ENTREZID %in% unique_geneid)
nrow(int)

str(int)
### pheatmap 
# Required packages
library(pheatmap)
library(dplyr)
library(pheatmap)

library(tidyverse)

labels_expr <- c(
  expression("C57BL/6J"~"_1"),
  expression("C57BL/6J"~"_2"),
  expression("C57BL/6J"~"_3"),
  expression("C57BL/6J"~"_4"),
  expression("C57BL/6J"~"_5"),
  
  
  expression("DDHD2"^"-/-"~"+ MyrPalmStearCoA_1"),
  expression("DDHD2"^"-/-"~"+ MyrPalmStearCoA_2"),
  expression("DDHD2"^"-/-"~"+ MyrPalmStearCoA_3"),
  expression("DDHD2"^"-/-"~"+ MyrPalmStearCoA_4"),
  expression("DDHD2"^"-/-"~"+ MyrPalmStearCoA_5"),
  
  expression("DDHD2"^"-/-"~"_1"),
  expression("DDHD2"^"-/-"~"_2"),
  expression("DDHD2"^"-/-"~"_3"),
  expression("DDHD2"^"-/-"~"_4"),
  expression("DDHD2"^"-/-"~"_5"),
  
  
  expression("DDHD2"^"-/-"~"+ Stearic-CoA_1"),
  expression("DDHD2"^"-/-"~"+ Stearic-CoA_2"),
  expression("DDHD2"^"-/-"~"+ Stearic-CoA_3"),
  expression("DDHD2"^"-/-"~"+ Stearic-CoA_4"),
  expression("DDHD2"^"-/-"~"+ Stearic-CoA_5")
)

colnames(int)
# Subset intensity data (remove annotation columns)
intensity_data <- int %>%
  dplyr::select(-Protein, -Row_ID, -ENTREZID, -SYMBOL)

# Convert to matrix
intensity_matrix <- as.matrix(intensity_data)

# Row-wise z-score (center and scale by row)
z_score_matrix <- t(scale(t(intensity_matrix)))

# Optional: set protein names as rownames
rownames(z_score_matrix) <- int$SYMBOL  # or use $Protein or $ENTREZID
writeClipboard(capture.output(rownames(z_score_matrix)))

library(showtext)

# Add TNR from system fonts (adjust the path if needed)
font_add("tnr", regular = "C:/Windows/Fonts/times.ttf")
showtext_auto()





# Z-score color breaks
breaks <- seq(-2, 2, length.out = 100)
Cairo::CairoPDF("Heatmap_Golgi_heatmap.pdf", width = 8, height = 12, family = "Times")
Cairo::CairoPNG("Heatmap_Golgi_heatmap.png", width = 800, height = 1200, family = "Times", res = 100)


# Plot
pheatmap(z_score_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = FALSE,
         legend = FALSE,
         fontsize = 20,
         cellheight = 1.8,
         cellwidth = 20,
         legend_breaks = c(-2, 0, 2),
         breaks = seq(-2, 2, length.out = 100),
         fontsize_col = 20,
         main = "Golgi Apparatus",
         #border_color = "black",
         labels_col = labels_expr  # 👈 parsed labels here
)

dev.off()


#####################################################################
#####################################################################
# Now get the gene ID for Mitochondria GO:0005739
#####################################################################
#####################################################################
cc <- read.csv("20250327_DDHD2_reduced_thres_08.csv")
list.files()
# presynapse
mito <- c("GO:0019910", "GO:0005741", "GO:0098798", "GO:0043294", "GO:0000262", "GO:0017133", "GO:0044290", "GO:0002187", "GO:0031019",
           "GO:0030678", "GO:0020023", "GO:0042645", "GO:0005740", "GO:0031966", "GO:0005759", "GO:0016507", "GO:0009841", "GO:0034245")

id <- cc %>%
  filter(ID %in% mito)
nrow(id)
allgeneid <- unlist(strsplit(id$geneID, split = "/"))
unique_geneid <- unique(allgeneid)

int <- enrich_final %>%
  filter(ENTREZID %in% unique_geneid)
nrow(int)

str(int)
### pheatmap 
# Required packages
library(pheatmap)
library(dplyr)
library(pheatmap)

library(tidyverse)

labels_expr <- c(
  expression("C57BL/6J"~"_1"),
  expression("C57BL/6J"~"_2"),
  expression("C57BL/6J"~"_3"),
  expression("C57BL/6J"~"_4"),
  expression("C57BL/6J"~"_5"),
  
  
  expression("DDHD2"^"-/-"~"+ MyrPalmStearCoA_1"),
  expression("DDHD2"^"-/-"~"+ MyrPalmStearCoA_2"),
  expression("DDHD2"^"-/-"~"+ MyrPalmStearCoA_3"),
  expression("DDHD2"^"-/-"~"+ MyrPalmStearCoA_4"),
  expression("DDHD2"^"-/-"~"+ MyrPalmStearCoA_5"),
  
  expression("DDHD2"^"-/-"~"_1"),
  expression("DDHD2"^"-/-"~"_2"),
  expression("DDHD2"^"-/-"~"_3"),
  expression("DDHD2"^"-/-"~"_4"),
  expression("DDHD2"^"-/-"~"_5"),
  
  
  expression("DDHD2"^"-/-"~"+ Stearic-CoA_1"),
  expression("DDHD2"^"-/-"~"+ Stearic-CoA_2"),
  expression("DDHD2"^"-/-"~"+ Stearic-CoA_3"),
  expression("DDHD2"^"-/-"~"+ Stearic-CoA_4"),
  expression("DDHD2"^"-/-"~"+ Stearic-CoA_5")
)

colnames(int)
# Subset intensity data (remove annotation columns)
intensity_data <- int %>%
  dplyr::select(-Protein, -Row_ID, -ENTREZID, -SYMBOL)

# Convert to matrix
intensity_matrix <- as.matrix(intensity_data)

# Row-wise z-score (center and scale by row)
z_score_matrix <- t(scale(t(intensity_matrix)))

# Optional: set protein names as rownames
rownames(z_score_matrix) <- int$SYMBOL  # or use $Protein or $ENTREZID
writeClipboard(capture.output(rownames(z_score_matrix)))

library(showtext)

# Add TNR from system fonts (adjust the path if needed)
font_add("tnr", regular = "C:/Windows/Fonts/times.ttf")
showtext_auto()





# Z-score color breaks
breaks <- seq(-2, 2, length.out = 100)
Cairo::CairoPDF("Heatmap_Mito_heatmap.pdf", width = 8, height = 12, family = "Times")
Cairo::CairoPNG("Heatmap_Mito_heatmap.png", width = 800, height = 1200, family = "Times", res = 100)


# Plot
pheatmap(z_score_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = FALSE,
         legend = FALSE,
         fontsize = 20,
         cellheight = 1.1,
         cellwidth = 20,
         legend_breaks = c(-2, 0, 2),
         breaks = seq(-2, 2, length.out = 100),
         fontsize_col = 20,
         main = "Mitochondrion",
         #border_color = "black",
         labels_col = labels_expr  # 👈 parsed labels here
)

dev.off()


#####################################################################
#####################################################################
# Now get the gene ID for Postsynapse GO:0098794
#####################################################################
#####################################################################
cc <- read.csv("20250327_DDHD2_reduced_thres_08.csv")
list.files()
post <- c("GO:0098837", "GO:1990475","GO:0045211","GO:0099572","GO:0099571","GO:0098871","GO:0099189","GO:0098975","GO:0150051","GO:0098842","GO:0099160","GO:0098843","GO:0098845",
          "GO:0098845","GO:0043197","GO:0099524")

id <- cc %>%
  filter(ID %in% post) #%>%
#  filter(grepl("postsynaptic endosome", Description))
unique(id$ID)
nrow(id)
allgeneid <- unlist(strsplit(id$geneID, split = "/"))
unique_geneid <- unique(allgeneid)

int <- enrich_final %>%
  filter(ENTREZID %in% unique_geneid)
nrow(int)

str(int)
### pheatmap 
# Required packages
library(pheatmap)
library(dplyr)
library(pheatmap)

library(tidyverse)

labels_expr <- c(
  expression("C57BL/6J"~"_1"),
  expression("C57BL/6J"~"_2"),
  expression("C57BL/6J"~"_3"),
  expression("C57BL/6J"~"_4"),
  expression("C57BL/6J"~"_5"),
  
  
  expression("DDHD2"^"-/-"~"+ MyrPalmStearCoA_1"),
  expression("DDHD2"^"-/-"~"+ MyrPalmStearCoA_2"),
  expression("DDHD2"^"-/-"~"+ MyrPalmStearCoA_3"),
  expression("DDHD2"^"-/-"~"+ MyrPalmStearCoA_4"),
  expression("DDHD2"^"-/-"~"+ MyrPalmStearCoA_5"),
  
  expression("DDHD2"^"-/-"~"_1"),
  expression("DDHD2"^"-/-"~"_2"),
  expression("DDHD2"^"-/-"~"_3"),
  expression("DDHD2"^"-/-"~"_4"),
  expression("DDHD2"^"-/-"~"_5"),
  
  
  expression("DDHD2"^"-/-"~"+ Stearic-CoA_1"),
  expression("DDHD2"^"-/-"~"+ Stearic-CoA_2"),
  expression("DDHD2"^"-/-"~"+ Stearic-CoA_3"),
  expression("DDHD2"^"-/-"~"+ Stearic-CoA_4"),
  expression("DDHD2"^"-/-"~"+ Stearic-CoA_5")
)

colnames(int)
# Subset intensity data (remove annotation columns)
intensity_data <- int %>%
  dplyr::select(-Protein, -Row_ID, -ENTREZID, -SYMBOL)

# Convert to matrix
intensity_matrix <- as.matrix(intensity_data)

# Row-wise z-score (center and scale by row)
z_score_matrix <- t(scale(t(intensity_matrix)))

# Optional: set protein names as rownames
rownames(z_score_matrix) <- int$SYMBOL  # or use $Protein or $ENTREZID
writeClipboard(capture.output(rownames(z_score_matrix)))

library(showtext)

# Add TNR from system fonts (adjust the path if needed)
font_add("tnr", regular = "C:/Windows/Fonts/times.ttf")
showtext_auto()





# Z-score color breaks
breaks <- seq(-2, 2, length.out = 100)
Cairo::CairoPDF("Heatmap_Postsynapse_heatmap.pdf", width = 8, height = 12, family = "Times")
Cairo::CairoPNG("Heatmap_Postsynapse_heatmap.png", width = 800, height = 1200, family = "Times", res = 100)


# Plot
pheatmap(z_score_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = FALSE,
         legend = FALSE,
         fontsize = 20,
         cellheight = 0.9,
         cellwidth = 20,
         legend_breaks = c(-2, 0, 2),
         breaks = seq(-2, 2, length.out = 100),
         fontsize_col = 20,         main = "Postsynapse",
         #border_color = "black",
         labels_col = labels_expr  # 👈 parsed labels here
)

dev.off()


#####################################################################
#####################################################################
# Now get the gene ID for TrasportVesicle GO:0030133
#####################################################################
#####################################################################
cc <- read.csv("20250327_DDHD2_reduced_thres_08.csv")
list.files()
tv <- c("GO:0060200", "GO:0070382", "GO:0070081", "GO:0060199","GO:0030658","GO:0030142","GO:0030140","GO:1990257","GO:0030143","GO:0098566")

id <- cc %>%
filter(ID %in% tv) #%>%
  #filter(grepl("transport vesicle", Description))
unique(id$ID)
nrow(id)
allgeneid <- unlist(strsplit(id$geneID, split = "/"))
unique_geneid <- unique(allgeneid)

int <- enrich_final %>%
  filter(ENTREZID %in% unique_geneid)
nrow(int)

str(int)
### pheatmap 
# Required packages
library(pheatmap)
library(dplyr)
library(pheatmap)

library(tidyverse)

labels_expr <- c(
  expression("C57BL/6J"~"_1"),
  expression("C57BL/6J"~"_2"),
  expression("C57BL/6J"~"_3"),
  expression("C57BL/6J"~"_4"),
  expression("C57BL/6J"~"_5"),
  
  
  expression("DDHD2"^"-/-"~"+ MyrPalmStearCoA_1"),
  expression("DDHD2"^"-/-"~"+ MyrPalmStearCoA_2"),
  expression("DDHD2"^"-/-"~"+ MyrPalmStearCoA_3"),
  expression("DDHD2"^"-/-"~"+ MyrPalmStearCoA_4"),
  expression("DDHD2"^"-/-"~"+ MyrPalmStearCoA_5"),
  
  expression("DDHD2"^"-/-"~"_1"),
  expression("DDHD2"^"-/-"~"_2"),
  expression("DDHD2"^"-/-"~"_3"),
  expression("DDHD2"^"-/-"~"_4"),
  expression("DDHD2"^"-/-"~"_5"),
  
  
  expression("DDHD2"^"-/-"~"+ Stearic-CoA_1"),
  expression("DDHD2"^"-/-"~"+ Stearic-CoA_2"),
  expression("DDHD2"^"-/-"~"+ Stearic-CoA_3"),
  expression("DDHD2"^"-/-"~"+ Stearic-CoA_4"),
  expression("DDHD2"^"-/-"~"+ Stearic-CoA_5")
)

colnames(int)
# Subset intensity data (remove annotation columns)
intensity_data <- int %>%
  dplyr::select(-Protein, -Row_ID, -ENTREZID, -SYMBOL)

# Convert to matrix
intensity_matrix <- as.matrix(intensity_data)

# Row-wise z-score (center and scale by row)
z_score_matrix <- t(scale(t(intensity_matrix)))

# Optional: set protein names as rownames
rownames(z_score_matrix) <- int$SYMBOL  # or use $Protein or $ENTREZID
writeClipboard(capture.output(rownames(z_score_matrix)))

library(showtext)

# Add TNR from system fonts (adjust the path if needed)
font_add("tnr", regular = "C:/Windows/Fonts/times.ttf")
showtext_auto()





# Z-score color breaks
breaks <- seq(-2, 2, length.out = 100)
CairoPDF("Heatmap_TrasportVesicle_heatmap.pdf", width = 8, height = 12, family = "Times")
Cairo::CairoPNG("Heatmap_TrasportVesicle_heatmap.png", width = 800, height = 1200, family = "Times", res = 100)


# Plot
pheatmap(z_score_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = FALSE,
         legend = FALSE,
         fontsize = 20,
         cellheight = 1.4,
         cellwidth = 20,
         legend_breaks = c(-2, 0, 2),
         breaks = seq(-2, 2, length.out = 100),
         fontsize_col = 20,
         main = "Transport Vesicle",
         #border_color = "black",
         labels_col = labels_expr  # 👈 parsed labels here
)

dev.off()

###################################################
###################################################
# GSEA
###################################################
###################################################
colnames(cc)
list.files()
cc <- read.csv("20250402_DDHD2_CC_reduced_thres_08.csv")
colnames(cc)
unique(cc$Ont)
ccplot <- cc %>%
  filter(Ont %in% "CC") %>%
  dplyr::select(Bait, group, Description,parentTerm, ratio, p.adjust, Count)

unique(ccplot$parentTerm)
unique(ccplot$Description)
cc_sum <- ccplot %>%
  group_by(Bait, group) %>%
  summarize(Count = n())
writeClipboard(capture.output(print(cc_cleaned, n =200)))

cc_cleaned <- ccplot %>%
  group_by(Bait, group, parentTerm) %>%
  arrange(desc(ratio), p.adjust, desc(Count)) %>%
  dplyr::slice(1) %>%  # keep the top entry per group
  ungroup() %>%
  filter(Count >= 5)
View(cc_cleaned)

unique(cc_cleaned$parentTerm)

colnames(cc_cleaned)
nrow(cc_cleaned)
nrow(cc_top)

#####################################
# Heatmapdrop

#####################################
# Create enrichment score
heatmap_data <- cc_cleaned %>%
  mutate(enrichment_score = ratio) %>%
  unite("contrast", Bait, group, sep = "_") %>%
  select(parentTerm, contrast, enrichment_score) %>%
  pivot_wider(names_from = contrast, values_from = enrichment_score)

# Turn NAs to 0 or NA (as preferred)
heatmap_data[is.na(heatmap_data)] <- 0

# Convert to matrix for heatmap
rownames(heatmap_data) <- heatmap_data$parentTerm
heatmap_matrix <- as.matrix(heatmap_data[,-1])
rownames(heatmap_matrix) <- rownames(heatmap_data)
# Plot heatmap
library(pheatmap)
pheatmap(heatmap_matrix, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         fontsize_row = 8, 
         show_rownames = TRUE,
         main = "-log10(p.adjust) Enrichment Scores")
#########################################################
################################################
# mULTI CONdition dot plot
########################################
cc_cleaned$Bait <- factor(cc_cleaned$Bait, levels = c("KO","STER","TRIPLE"))

cc_cleaned <- cc_cleaned %>%
mutate(parentTerm = recode(parentTerm,
                           "endoplasmic reticulum-Golgi intermediate compartment" = "ERGIC")) %>%
  mutate(
    enrichment_score = -log10(p.adjust),
    Bait_label = recode(Bait,
                        "KO" = "DDHD2^'-/-'",
                        "STER" = "DDHD2^'-/-'~'+ Stearic-CoA'",
                        "TRIPLE" = "DDHD2^'-/-'~'+ MyrPalmStearCoA'"
    )
  )



cc_top <- cc_cleaned %>%
  group_by(Bait, group) %>%
  slice_max(order_by = -log10(p.adjust), n = 20) %>%
  mutate(parentTerm = recode(parentTerm,
                             "endoplasmic reticulum-Golgi intermediate compartment" = "ERGIC")) %>%
  mutate(
    enrichment_score = -log10(p.adjust),
    Bait_label = recode(Bait,
                        "KO" = "DDHD2^'-/-'",
                        "STER" = "DDHD2^'-/-'~'+ Stearic-CoA'",
                        "TRIPLE" = "DDHD2^'-/-'~'+ MyrPalmStearCoA'"
    )
  )

# Plot
ggplot(cc_cleaned, aes(x = ratio, y = reorder(parentTerm, ratio))) +
  geom_point(aes(size = Count, color = enrichment_score)) +
  facet_grid(
    group ~ Bait_label,
    scales = "free_y",
    space = "free",
    labeller = labeller(Bait_label = label_parsed)
  ) +
  scale_color_viridis_c(
    name = expression(-log[10]~"(adjusted p)"),
    option = "plasma"
  ) +
  scale_size(range = c(3, 10)) +
  theme_bw() +
  labs(
    x = "GeneRatio",
    y = "GO Term",
    title = "GO Enrichment across Conditions"
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = rel(4)),
    axis.text.y = element_text(size = rel(3)),
    axis.text.x = element_text(size = rel(4)),
    strip.text.y = element_text(size = rel(4), face = "bold"),
    strip.text.x = element_text(size = rel(3.5)),
    legend.title = element_text(size = rel(4)),
    legend.text = element_text(size = rel(3)),
    plot.title = element_blank()
  )




# ggsave("GSEA_all.png", width = 10, height = 10)


ggsave("GSEA_all.pdf", width = 24, height = 24, dpi = 300)


###############################################
## Bullet point 6
###############################################
# Define relevant biological process GO terms (adjust as needed)
# ===============================
# 🎨 Font Setup for Publication
# ===============================
font_add("Times", regular = "C:/Windows/Fonts/times.ttf")
showtext_auto()

# ===============================
# 📁 Load GO Enrichment & Expression Matrix
# ===============================
data <- read.xlsx("20250404_DDHD2_enrichGO_reduced_thres_08.xlsx")
df <- read.xlsx("Normalized_imputed_sva_batchCorrected_matrix_annotated.xlsx")
conflict_prefer("filter", "dplyr")
# ===============================
# 🔍 Select GO Terms of Interest (β-oxidation & Carnitine)
# ===============================
View(data)
fa <- c("GO:0006631", "GO:0015908", "GO:0009437")
beta <- data %>%
  # filter(str_detect(Description, regex("Fatty Acid|Carnitine", ignore_case = TRUE))) %>%
  filter(ID %in% fa) %>%
  separate_rows(geneID, sep = "/") %>%
  arrange(geneID, p.adjust) %>%  # <- Replace `adj_pval` with your column name
  group_by(geneID) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(geneID, Description, parentTerm, p.adjust, ID) %>%
  distinct()
View(beta)
allgeneid <- unlist(strsplit(beta$geneID, split = "/"))

# ===============================
# 🧬 Subset Expression Matrix by Selected Genes
# ===============================
df_sub <- df %>% 
  filter(ENTREZID %in% allgeneid) %>%
  left_join(beta, by = c("ENTREZID" = "geneID")) %>%
  filter(!is.na(Description))

unique(df_sub$ENTREZID)

expr_mat <- df_sub %>%
  dplyr::select(SYMBOL, starts_with("CTRL"), starts_with("KO_"), starts_with("KO_MyrPalmSterie")) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()
nrow(expr_mat)

zscore_mat <- t(scale(t(expr_mat)))

# ===============================
# 🏷️ Annotate Sample Groups
# ===============================
group_annot <- data.frame(
  Group = case_when(
    grepl("^CTRL", colnames(zscore_mat)) ~ "CTRL",
    grepl("^KO_Sterie", colnames(zscore_mat)) ~ "STER",
    grepl("^KO_MyrPalmSterie", colnames(zscore_mat)) ~ "TRIPLE",
    grepl("^KO_", colnames(zscore_mat)) ~ "KO"
  )
)
rownames(group_annot) <- colnames(zscore_mat)

# ===============================
# 🧬 Pathway Label per Protein (β-oxidation, Carnitine, or Both)
# ===============================

priority <- c("fatty acid transport", "fatty acid metabolic process", "carnitine metabolic process")
nrow(df_sub)
print(unique(df_sub %>% select(Description, ID)))
row_annot <- df_sub %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  rename(Pathway = Description) %>%
  select(SYMBOL, Pathway) %>%
  mutate(Pathway = factor(Pathway, levels = priority),
         Pathway = recode(Pathway,
                          "fatty acid metabolic process" =  "fatty acid metabolic process\nGO:0006631",
                          "fatty acid transport" = "fatty acid transport\nGO:0015908" , 
                          "carnitine metabolic process" = "carnitine metabolic process\nGO:0009437"
         )) %>%
  
  column_to_rownames("SYMBOL")

# ===============================
# 🎨 Row Annotation (Pathway)
# ===============================

color <- c(brewer.pal(3, "Set2"))
#color <- c(brewer.pal(8, "Set2"), brewer.pal(7, "Set1"))
row_ha <- rowAnnotation(
  Pathway = row_annot$Pathway,
  gp = gpar(fontsize = rel(20), fontfamily = "Times", col = "black"),
  col = list(Pathway = setNames(color , levels(row_annot$Pathway))),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = rel(20), fontfamily = "Times", fontface = "bold"), # THIS IS the word pathway 
  show_legend = FALSE
)

# ===============================
# ✅ Sample Order (manually sorted)
# ===============================
ordered_columns <- c(
  grep("^CTRL", colnames(zscore_mat), value = TRUE),
  grep("^KO_[0-9]", colnames(zscore_mat), value = TRUE),
  grep("^KO_Sterie", colnames(zscore_mat), value = TRUE),
  grep("^KO_MyrPalmSterie", colnames(zscore_mat), value = TRUE)
)
zscore_mat_ordered <- zscore_mat[, ordered_columns]
group_annot_ordered <- group_annot[ordered_columns, , drop = FALSE]

# ===============================
# 🧾 Bottom Labels as Expression
# ===============================
group_labels_expr <- list(
  "C57BL/6J_1",
  "C57BL/6J_2",
  "C57BL/6J_3",
  "C57BL/6J_4",
  "C57BL/6J_5",
  
  "DDHD2⁻/⁻_1",
  "DDHD2⁻/⁻_2",
  "DDHD2⁻/⁻_3",
  "DDHD2⁻/⁻_4",
  "DDHD2⁻/⁻_5",
  
  "DDHD2⁻/⁻ + Stearic-CoA_1",
  "DDHD2⁻/⁻ + Stearic-CoA_2",
  "DDHD2⁻/⁻ + Stearic-CoA_3",
  "DDHD2⁻/⁻ + Stearic-CoA_4",
  "DDHD2⁻/⁻ + Stearic-CoA_5",
  
  "DDHD2⁻/⁻ + MyrPalmStearCoA_1",
  "DDHD2⁻/⁻ + MyrPalmStearCoA_2",
  "DDHD2⁻/⁻ + MyrPalmStearCoA_3",
  "DDHD2⁻/⁻ + MyrPalmStearCoA_4",
  "DDHD2⁻/⁻ + MyrPalmStearCoA_5"
)

col_ha <- HeatmapAnnotation(
  label = anno_text(
    group_labels_expr,
    gp = gpar(fontsize = rel(20), fontfamily = "Times", col = "black"),
    rot = 90,
    just = "right"
  ),
  annotation_name_side = "left"
)

# ===============================
# 🎨 Heatmap Color Function
# ===============================
col_fun <- colorRamp2(
  seq(-2, 2, length.out = 11),
  brewer.pal(11, "RdYlBu")
)

# ===============================
# 📊 Plot the Heatmap
# ===============================
CairoPDF("Heatmap_FA_carnitine_heatmap.pdf", width = 16, height = 34)
CairoPDF("Heatmap_FA_carnitine_heatmap_collapse.pdf", width = 16, height = 18)
# CairoPNG("Terms_FA_carnitine_heatmap.png", width = 1600, height = 2000)
Heatmap(
  zscore_mat_ordered,
  name = "Z-score",
  col = col_fun,
  #width = unit(20, "cm"),
  #height = unit(25, "cm"),
  left_annotation = row_ha,
  bottom_annotation = col_ha,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  # show_row_names = FALSE,
  show_column_names = FALSE,
  # 🔥 Make Z-score legend text bigger
  heatmap_legend_param = list(
    legend_height = unit(8, "cm"),  # try increasing this
    title_gp = gpar(fontsize = rel(24), fontface = "bold"),
    labels_gp = gpar(fontsize = rel(24))),
  row_split = row_annot$Pathway,
  row_title_gp = gpar(fontsize = rel(24), fontfamily = "Times", fontface = "bold"),
  row_title_rot = 0,
  row_names_gp = gpar(fontsize = rel(24),fontfamily = "Times"),
  row_title_side = "left"
)

dev.off()

###############################################
###############################################
# ROS, redox etc
###############################################
# Set ROS-related parentTerm to focus on
# Define main priority terms (including parent and new "children")
colnames(data)
oxi <-  c("response to oxidative stress",
          "reactive oxygen species metabolic process")
# Create tidy data: replace oxidative stress parentTerm with its child Description
ros <- data %>%
  # filter(grepl("ROS|reactive oxygen species|oxidative stress", Description, ignore.case = FALSE))
  filter(Description %in% oxi) %>%
  select(geneID, Description, p.adjust, ID) %>%
  # rename(Description = tidy_category) %>%
  distinct()


unique(ros$Description)
priority <- c(
  "response to topologically incorrect protein",
  "response to oxidative stress",              
  "reactive oxygen species metabolic process"
)

# Expand and prioritize genes (only 1 term here)
ros_prioritized <- ros %>%
  separate_rows(geneID, sep = "/") %>%
  arrange(geneID, p.adjust) %>%  # <- Replace `adj_pval` with your column name
  group_by(geneID) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(geneID, Description, p.adjust, ID) %>%
  distinct() %>%
  rename(parentTerm = Description) %>%
  select(parentTerm, geneID)
# View(ros_prioritized)
# Get matching genes
allgeneid <- unique(ros_prioritized$geneID)
View(ros_prioritized)
# Subset intensity matrix
df_sub <- df %>%
  mutate(ENTREZID = as.character(ENTREZID)) %>%
  filter(ENTREZID %in% allgeneid) 
# Join pathway info
df_sub_prioritized <- df_sub %>%
  left_join(ros_prioritized, by = c("ENTREZID" = "geneID")) %>%
  rename(Pathway = parentTerm) %>%
  filter(!is.na(Pathway))
unique(df_sub_prioritized$Pathway)
# Prepare expression matrix
expr_mat <- df_sub_prioritized %>%
  select(SYMBOL, starts_with("CTRL"), starts_with("KO_")) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()

# Z-score transform
zscore_mat <- t(scale(t(expr_mat)))

# Sample annotations
group_labels_expr <- c(
  rep("C57BL/6J", 5),
  rep("DDHD2⁻/⁻", 5),
  rep("DDHD2⁻/⁻ + Stearic-CoA", 5),
  rep("DDHD2⁻/⁻ + MyrPalmStearCoA", 5)
)

group_annot <- data.frame(
  Group = rep(c("CTRL", "KO", "STER", "TRIPLE"), each = 5)
)
rownames(group_annot) <- colnames(zscore_mat)

# Order columns
ordered_columns <- c(
  grep("^CTRL", colnames(zscore_mat), value = TRUE),
  grep("^KO_[0-9]", colnames(zscore_mat), value = TRUE),
  grep("^KO_Sterie", colnames(zscore_mat), value = TRUE),
  grep("^KO_MyrPalmSterie", colnames(zscore_mat), value = TRUE)
)
zscore_mat_ordered <- zscore_mat[, ordered_columns]
nrow(zscore_mat_ordered)
unique(df_sub_prioritized$Pathway)
View(df_sub_prioritized)
print(unique(ros %>% select(Description,ID)))
# Row annotation
row_annot <- df_sub_prioritized %>%
  filter(SYMBOL %in% rownames(zscore_mat)) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  select(SYMBOL, Pathway) %>%
  mutate(Pathway = factor(Pathway, levels = c(
    "response to oxidative stress",
    "reactive oxygen species metabolic process")),
  Pathway = recode(Pathway, 
                   "response to oxidative stress" ="response to oxidative stress\nGO:0006979",
                   "reactive oxygen species metabolic process" = "reactive oxygen species \nmetabolic process \nGO:0072593")
  ) %>%
  column_to_rownames("SYMBOL")

color <- c(brewer.pal(3, "Set2"))
# color <- c(brewer.pal(3, "Set2"), brewer.pal(4, "Set1"))
row_ha <- rowAnnotation(
  Pathway = row_annot$Pathway,
  gp = gpar(fontsize = rel(20), fontfamily = "Times", col = "black"),
  col = list(Pathway = setNames(color[-1], levels(row_annot$Pathway))),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = rel(20), fontfamily = "Times", fontface = "bold"), # THIS IS the word pathway 
  show_legend = FALSE
)

col_ha <- HeatmapAnnotation(
  label = anno_text(group_labels_expr, rot = 90, just = "right", gp = gpar(fontsize = rel(20),
                                                                           fontfamily = "Times")),
  annotation_name_side = "left"
)

col_fun <- colorRamp2(seq(-2, 2, length.out = 11), brewer.pal(11, "RdYlBu"))

# Save heatmap
CairoPDF("Heatmap_ROS_Redox_heatmap.pdf", width = 20, height = 34)
# CairoPDF("Heatmap_ROS_Redox_heatmap_collapse.pdf", width = 16, height = 18 )
# CairoPNG("ROS_Redox_heatmap.png", width = 1800, height = 2800)
Heatmap(
  zscore_mat_ordered,
  name = "Z-score",
  col = col_fun,
  # width = unit(20, "cm"),
  # height = unit(25, "cm"),
  left_annotation = row_ha,
  bottom_annotation = col_ha,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  # show_row_names = FALSE,
  show_column_names = FALSE,
  # 🔥 Make Z-score legend text bigger
  heatmap_legend_param = list(
    legend_height = unit(8, "cm"),  # try increasing this
    title_gp = gpar(fontsize = rel(24), fontface = "bold"),
    labels_gp = gpar(fontsize = rel(24))),
  row_split = row_annot$Pathway,
  row_title_gp = gpar(fontsize = rel(28), fontfamily = "Times", fontface = "bold"),
  row_title_rot = 0,
  row_names_gp = gpar(fontsize = rel(24),fontfamily = "Times"),
  row_title_side = "left"
)
dev.off()


############################################################
############################################################
# Mitochondrial fusion / fission
############################################################
############################################################
# Set ROS-related parentTerm to focus on
# Define main priority terms (including parent and new "children")
colnames(data)
oxi <-  c("mitochondrial fusion",
          "mitochondrial fission")
# Create tidy data: replace oxidative stress parentTerm with its child Description
View(data)
ros <- data %>%
  # filter(grepl("ROS|reactive oxygen species|oxidative stress", Description, ignore.case = FALSE))
  filter(Description %in% oxi) %>%
  select(geneID, Description, p.adjust, ID) %>%
  # rename(Description = tidy_category) %>%
  distinct()


unique(ros$Description)
priority <- oxi

# Expand and prioritize genes (only 1 term here)
ros_prioritized <- ros %>%
  separate_rows(geneID, sep = "/") %>%
  arrange(geneID, p.adjust) %>%  # <- Replace `adj_pval` with your column name
  group_by(geneID) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(geneID, Description, p.adjust, ID) %>%
  distinct() %>%
  rename(parentTerm = Description) %>%
  select(parentTerm, geneID)
# View(ros_prioritized)
# Get matching genes
allgeneid <- unique(ros_prioritized$geneID)

# Subset intensity matrix
df_sub <- df %>%
  mutate(ENTREZID = as.character(ENTREZID)) %>%
  filter(ENTREZID %in% allgeneid) 
# Join pathway info
df_sub_prioritized <- df_sub %>%
  left_join(ros_prioritized, by = c("ENTREZID" = "geneID")) %>%
  rename(Pathway = parentTerm) %>%
  filter(!is.na(Pathway))
unique(df_sub_prioritized$Pathway)
# Prepare expression matrix
expr_mat <- df_sub_prioritized %>%
  select(SYMBOL, starts_with("CTRL"), starts_with("KO_")) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()

# Z-score transform
zscore_mat <- t(scale(t(expr_mat)))

# Sample annotations
group_labels_expr <- c(
  rep("C57BL/6J", 5),
  rep("DDHD2⁻/⁻", 5),
  rep("DDHD2⁻/⁻ + Stearic-CoA", 5),
  rep("DDHD2⁻/⁻ + MyrPalmStearCoA", 5)
)

group_annot <- data.frame(
  Group = rep(c("CTRL", "KO", "STER", "TRIPLE"), each = 5)
)
rownames(group_annot) <- colnames(zscore_mat)

# Order columns
ordered_columns <- c(
  grep("^CTRL", colnames(zscore_mat), value = TRUE),
  grep("^KO_[0-9]", colnames(zscore_mat), value = TRUE),
  grep("^KO_Sterie", colnames(zscore_mat), value = TRUE),
  grep("^KO_MyrPalmSterie", colnames(zscore_mat), value = TRUE)
)
zscore_mat_ordered <- zscore_mat[, ordered_columns]
nrow(zscore_mat_ordered)
unique(df_sub_prioritized$Pathway)

print(unique(ros %>% select(Description,ID)))
# Row annotation
row_annot <- df_sub_prioritized %>%
  filter(SYMBOL %in% rownames(zscore_mat)) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  select(SYMBOL, Pathway) %>%
  mutate(Pathway = factor(Pathway, levels = c(
    "mitochondrial fusion",
    "mitochondrial fission")),
    Pathway = recode(Pathway, 
                     "mitochondrial fusion" ="mitochondrial fusion\nGO:0008053",
                     "mitochondrial fission" = "mitochondrial fission\nGO:0000266")
  ) %>%
  column_to_rownames("SYMBOL")

color <- c(brewer.pal(3, "Set2"))
# color <- c(brewer.pal(3, "Set2"), brewer.pal(4, "Set1"))
row_ha <- rowAnnotation(
  Pathway = row_annot$Pathway,
  gp = gpar(fontsize = rel(20), fontfamily = "Times", col = "black"),
  col = list(Pathway = setNames(color[-1], levels(row_annot$Pathway))),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = rel(20), fontfamily = "Times", fontface = "bold"), # THIS IS the word pathway 
  show_legend = FALSE
)

col_ha <- HeatmapAnnotation(
  label = anno_text(group_labels_expr, rot = 90, just = "right", gp = gpar(fontsize = rel(20),
                                                                           fontfamily = "Times")),
  annotation_name_side = "left"
)

col_fun <- colorRamp2(seq(-2, 2, length.out = 11), brewer.pal(11, "RdYlBu"))

# Save heatmap
# CairoPDF("Heatmap_mito_heatmap.pdf", width = 20, height = 34)
CairoPDF("Heatmap_mito_heatmap_collapse.pdf", width = 16, height = 18 )
# CairoPNG("ROS_Redox_heatmap.png", width = 1800, height = 2800)
Heatmap(
  zscore_mat_ordered,
  name = "Z-score",
  col = col_fun,
  width = unit(20, "cm"),
  height = unit(25, "cm"),
  left_annotation = row_ha,
  bottom_annotation = col_ha,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  # show_row_names = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  # 🔥 Make Z-score legend text bigger
  heatmap_legend_param = list(
    legend_height = unit(8, "cm"),  # try increasing this
    title_gp = gpar(fontsize = rel(24), fontface = "bold"),
    labels_gp = gpar(fontsize = rel(24))),
  row_split = row_annot$Pathway,
  row_title_gp = gpar(fontsize = rel(28), fontfamily = "Times", fontface = "bold"),
  row_title_rot = 0,
  row_names_gp = gpar(fontsize = rel(24),fontfamily = "Times"),
  row_title_side = "left"
)
dev.off()

### START FROM HERE ###
##############################################
##############################################
# Reactome glyco
##############################################
##############################################
sig_combined_results <- read.csv("EnrichR_pathways_DE.csv")

unique(sig_combined_results$Database)
see <- sig_combined_results %>%
  filter(Database %in% "Reactome_2022") %>%
  separate(Term, into = c("Description", "ID"), sep = " R-HSA-", remove = FALSE)

ros <- c(
  "Glycolysis",
  "Glucose Metabolism",
  "Gluconeogenesis",
  "Glycogen Breakdown (Glycogenolysis)"
)

library(dplyr)
colnames(see)
beta <- see %>%
  mutate(Description = trimws(Description)) %>%
  filter(Description %in% ros) %>%
  # filter(str_detect(Description, regex("glyco|gluco", ignore_case = TRUE))) %>%
  select(Genes, Description, Adjusted.P.value, ID) %>%
  distinct()



print(unique(beta %>% select(Description, ID)))
# Expand and prioritize genes (only 1 term here)
ros_prioritized <- beta %>%
  separate_rows(Genes, sep = ";") %>%
  mutate(Genes = trimws(Genes)) %>%
  arrange(Genes, Adjusted.P.value) %>%  # <- Replace `adj_pval` with your column name
  group_by(Genes) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  rename(parentTerm = Description) %>%
  select(parentTerm, Genes)

# colnames(ros_prioritized)
allgeneid <- unique(ros_prioritized$Genes)
# print(allgeneid)
# Subset intensity matrix
df_sub <- df %>%
  filter(tolower(SYMBOL) %in% tolower(allgeneid))
# colnames(df_sub)
# Join pathway info
df_sub_prioritized <- df %>%
  mutate(SYMBOL_lower = tolower(SYMBOL)) %>%
  left_join(
    ros_prioritized %>% mutate(Genes_lower = tolower(Genes)),
    by = c("SYMBOL_lower" = "Genes_lower")
  ) %>%
  select(-SYMBOL_lower) %>%
  rename(Pathway = parentTerm) %>%
  filter(!is.na(Pathway))

# Prepare expression matrix
expr_mat <- df_sub_prioritized %>%
  select(SYMBOL, starts_with("CTRL"), starts_with("KO_")) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()
nrow(expr_mat)
# Z-score transform
zscore_mat <- t(scale(t(expr_mat)))

# Sample annotations
group_labels_expr <- c(
  rep("C57BL/6J", 5),
  rep("DDHD2⁻/⁻", 5),
  rep("DDHD2⁻/⁻ + Stearic-CoA", 5),
  rep("DDHD2⁻/⁻ + MyrPalmStearCoA", 5)
)

group_annot <- data.frame(
  Group = rep(c("CTRL", "KO", "STER", "TRIPLE"), each = 5)
)
rownames(group_annot) <- colnames(zscore_mat)

# Order columns
ordered_columns <- c(
  grep("^CTRL", colnames(zscore_mat), value = TRUE),
  grep("^KO_[0-9]", colnames(zscore_mat), value = TRUE),
  grep("^KO_Sterie", colnames(zscore_mat), value = TRUE),
  grep("^KO_MyrPalmSterie", colnames(zscore_mat), value = TRUE)
)
zscore_mat_ordered <- zscore_mat[, ordered_columns]
nrow(zscore_mat_ordered)
nrow(df_sub_prioritized)
# Row annotation
print(unique(beta %>% select(Description, ID)))
row_annot <- df_sub_prioritized %>%
  filter(SYMBOL %in% rownames(zscore_mat)) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  select(SYMBOL, Pathway) %>%
  mutate(Pathway = factor(Pathway),
         Pathway = recode(Pathway,
                          "Glucose Metabolism" = "Glucose Metabolism \nR-HSA-70326",
                          "Gluconeogenesis" = "Gluconeogenesis \nR-HSA-70263",
                          "Glycolysis" = "Glycolysis \nR-HSA-70171",
                          "Glycogen Breakdown (Glycogenolysis)" = "Glycogen Breakdown (Glycogenolysis) \nR-HSA-70221")) %>%
  column_to_rownames("SYMBOL")
nrow(row_annot)
unique(row_annot$Pathway)
# set the levels you want to show in heatmap
color <- c(brewer.pal(4, "Set2"))
#color <- c(brewer.pal(8, "Set2"), brewer.pal(7, "Set1"))
row_ha <- rowAnnotation(
  Pathway = row_annot$Pathway,
  gp = gpar(fontsize = rel(20), fontfamily = "Times", col = "black"),
  col = list(Pathway = setNames(color, levels(row_annot$Pathway))),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = rel(20), fontfamily = "Times", fontface = "bold"), # THIS IS the word pathway 
  show_legend = FALSE
)

col_ha <- HeatmapAnnotation(
  label = anno_text(group_labels_expr, rot = 90, just = "right", gp = gpar(fontsize = rel(20),
                                                                           fontfamily = "Times")),
  annotation_name_side = "left"
)

col_fun <- colorRamp2(seq(-2, 2, length.out = 11), brewer.pal(11, "RdYlBu"))
# Save heatmap
# CairoPDF("Reactome_glycolysis_heatmap.pdf", width = 14, height = 18)
CairoPDF("Reactome_glycolysis_heatmap_collapse.pdf", width = 18, height = 16)
# CairoPNG("Reactome_glycolysis_heatmap.png", width = 1400, height = 1800)
Heatmap(
  zscore_mat_ordered,
  name = "Z-score",
  col = col_fun,
  width = unit(20, "cm"),
  height = unit(25, "cm"),
  left_annotation = row_ha,
  bottom_annotation = col_ha,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  # show_row_names = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  # 🔥 Make Z-score legend text bigger
  heatmap_legend_param = list(
    legend_height = unit(8, "cm"),  # try increasing this
    title_gp = gpar(fontsize = rel(24), fontface = "bold"),
    labels_gp = gpar(fontsize = rel(24))),
  row_split = row_annot$Pathway,
  row_title_gp = gpar(fontsize = rel(28), fontfamily = "Times", fontface = "bold"),
  row_title_rot = 0,
  row_names_gp = gpar(fontsize = rel(24),fontfamily = "Times"),
  row_title_side = "left"
)
dev.off()

