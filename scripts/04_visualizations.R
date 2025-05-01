# 04 Visualizations.R

# 04 Visualizations.R

#######################################
# Make volcano plot usig DE dataframe #
#######################################
#---------------------------------#
# Font setup
#---------------------------------#
font_add("tnr", regular = "C:/Windows/Fonts/times.ttf")
showtext_auto()
#---------------------------------#
# Load DE result with Bait column
#---------------------------------
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
