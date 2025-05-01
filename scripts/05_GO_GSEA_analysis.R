# 05 Go Gsea Analysis.R

# 05 Go Gsea Analysis.R
##############################
## GO
##############################
# Define the list of ontologies, this case run all BP, CC , MF
df <- read.xlsx("sva_batch_corrected_DE.xlsx")
pval <- 0.05
fc <- 1
# Keep only the first UniProt accession before the first semicolon
df$UniProt <- sub(";.*", "", df$Protein)
df <- df %>%
  mutate(Row_ID = row_number())
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

