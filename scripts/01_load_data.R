# 01 Load Data.R

# 01 Load Data.R
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
