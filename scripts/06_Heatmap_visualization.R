# 06_Heatmap_figures
##############################################################
# recreate heatmap figures
##############################################################
#####################################################################
# Now get the gene ID for presynapse
#####################################################################
cc <- read.xlsx("20250404_DDHD2_enrichGO_reduced_thres_08.xlsx")
list.files()
# presynapse
presynapse <- c("GO:0042734", "GO:0098830", "GO:0098833", "GO:0043679", "GO:0043195", "GO:0099569","GO:0099523","GO:0008021", "GO:0099143", "GO:0099182", "GO:0048786")
id <- cc %>%
  filter(ID %in% presynapse)
allgeneid <- unlist(strsplit(id$geneID, split = "/"))
unique_geneid <- unique(allgeneid)
int <- enrich_final %>%
  filter(ENTREZID %in% unique_geneid)
### pheatmap
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
# presynapse
er <- c("GO:0009511", "GO:0016529", "GO:0005788", "GO:0009510",
        "GO:0044322", "GO:0000835", "GO:0005790", "GO:0098827",
        "GO:0005791", "GO:0070971", "GO:0097038", "GO:0140534", "GO:1990007")
id <- cc %>%
  filter(ID %in% er)
allgeneid <- unlist(strsplit(id$geneID, split = "/"))
unique_geneid <- unique(allgeneid)
int <- enrich_final %>%
  filter(ENTREZID %in% unique_geneid)
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
allgeneid <- unlist(strsplit(id$geneID, split = "/"))
unique_geneid <- unique(allgeneid)
int <- enrich_final %>%
  filter(ENTREZID %in% unique_geneid)
### pheatmap 
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
mito <- c("GO:0019910", "GO:0005741", "GO:0098798", "GO:0043294", "GO:0000262", "GO:0017133", "GO:0044290", "GO:0002187", "GO:0031019",
          "GO:0030678", "GO:0020023", "GO:0042645", "GO:0005740", "GO:0031966", "GO:0005759", "GO:0016507", "GO:0009841", "GO:0034245")
id <- cc %>%
  filter(ID %in% mito)
allgeneid <- unlist(strsplit(id$geneID, split = "/"))
unique_geneid <- unique(allgeneid)
int <- enrich_final %>%
  filter(ENTREZID %in% unique_geneid)
### pheatmap 
# Required packages
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
post <- c("GO:0098837", "GO:1990475","GO:0045211","GO:0099572","GO:0099571","GO:0098871","GO:0099189","GO:0098975","GO:0150051","GO:0098842","GO:0099160","GO:0098843","GO:0098845",
          "GO:0098845","GO:0043197","GO:0099524")

id <- cc %>%
  filter(ID %in% post)
allgeneid <- unlist(strsplit(id$geneID, split = "/"))
unique_geneid <- unique(allgeneid)

int <- enrich_final %>%
  filter(ENTREZID %in% unique_geneid)
# Subset intensity data (remove annotation columns)
intensity_data <- int %>%
  dplyr::select(-Protein, -Row_ID, -ENTREZID, -SYMBOL)

# Convert to matrix
intensity_matrix <- as.matrix(intensity_data)

# Row-wise z-score (center and scale by row)
z_score_matrix <- t(scale(t(intensity_matrix)))

# Optional: set protein names as rownames
rownames(z_score_matrix) <- int$SYMBOL  # or use $Protein or $ENTREZID

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
tv <- c("GO:0060200", "GO:0070382", "GO:0070081", "GO:0060199","GO:0030658","GO:0030142","GO:0030140","GO:1990257","GO:0030143","GO:0098566")
id <- cc %>%
  filter(ID %in% tv) #%>%
allgeneid <- unlist(strsplit(id$geneID, split = "/"))
unique_geneid <- unique(allgeneid)

int <- enrich_final %>%
  filter(ENTREZID %in% unique_geneid)
### pheatmap 
# Subset intensity data (remove annotation columns)
intensity_data <- int %>%
  dplyr::select(-Protein, -Row_ID, -ENTREZID, -SYMBOL)

# Convert to matrix
intensity_matrix <- as.matrix(intensity_data)

# Row-wise z-score (center and scale by row)
z_score_matrix <- t(scale(t(intensity_matrix)))

# Optional: set protein names as rownames
rownames(z_score_matrix) <- int$SYMBOL  # or use $Protein or $ENTREZID

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
allgeneid <- unlist(strsplit(beta$geneID, split = "/"))

# ===============================
# 🧬 Subset Expression Matrix by Selected Genes
# ===============================
df_sub <- df %>% 
  filter(ENTREZID %in% allgeneid) %>%
  left_join(beta, by = c("ENTREZID" = "geneID")) %>%
  filter(!is.na(Description))
expr_mat <- df_sub %>%
  dplyr::select(SYMBOL, starts_with("CTRL"), starts_with("KO_"), starts_with("KO_MyrPalmSterie")) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()

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
oxi <-  c("response to oxidative stress",
          "reactive oxygen species metabolic process")
# Create tidy data: replace oxidative stress parentTerm with its child Description
ros <- data %>%
  filter(Description %in% oxi) %>%
  select(geneID, Description, p.adjust, ID) %>%
  distinct()

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
# Subset intensity matrix
df_sub <- df %>%
  mutate(ENTREZID = as.character(ENTREZID)) %>%
  filter(ENTREZID %in% allgeneid) 
# Join pathway info
df_sub_prioritized <- df_sub %>%
  left_join(ros_prioritized, by = c("ENTREZID" = "geneID")) %>%
  rename(Pathway = parentTerm) %>%
  filter(!is.na(Pathway))
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
ros <- data %>%
  # filter(grepl("ROS|reactive oxygen species|oxidative stress", Description, ignore.case = FALSE))
  filter(Description %in% oxi) %>%
  select(geneID, Description, p.adjust, ID) %>%
  # rename(Description = tidy_category) %>%
  distinct()
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
##############################################
##############################################
# Reactome glyco
##############################################
##############################################
sig_combined_results <- read.csv("EnrichR_pathways_DE.csv")
see <- sig_combined_results %>%
  filter(Database %in% "Reactome_2022") %>%
  separate(Term, into = c("Description", "ID"), sep = " R-HSA-", remove = FALSE)

ros <- c(
  "Glycolysis",
  "Glucose Metabolism",
  "Gluconeogenesis",
  "Glycogen Breakdown (Glycogenolysis)"
)

beta <- see %>%
  mutate(Description = trimws(Description)) %>%
  filter(Description %in% ros) %>%
  # filter(str_detect(Description, regex("glyco|gluco", ignore_case = TRUE))) %>%
  select(Genes, Description, Adjusted.P.value, ID) %>%
  distinct()
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

