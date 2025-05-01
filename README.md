# DDHD2 Proteomics Analysis

This repository contains the full R-based analysis pipeline for label-free quantitative proteomics comparing **DDHD2 knockout (KO)** and various rescue conditions in mouse neurons. The code supports data preprocessing, normalization, imputation, batch correction, differential expression analysis, pathway enrichment, and figure generation.

---

## 📁 Repository Structure

```
DDHD2_proteomics/
├── scripts/                 # Modular R scripts
│   ├── 01_load_data.R
│   ├── 02_normalize_impute.R
│   ├── 03_batch_correction_DE.R
│   ├── 04_visualizations.R
│   ├── 05_GO_GSEA_analysis.R
│   ├── 06_Heatmap_visualization
│   └── utils.R
├── .gitignore
├── DDHD2_proteomics.Rproj
└── README.md               # You're here
```

---

## 🚀 How to Reproduce the Analysis

1. Clone or download this repository:
   ```bash
   git clone https://github.com/YOUR_USERNAME/DDHD2_proteomics.git
   ```

2. Open `DDHD2_proteomics.Rproj` in RStudio.

3. Run each script in order from the `/scripts/` folder:
   - `01_load_data.R`: Reads in raw proteomics files and sample metadata.
   - `02_normalize_impute.R`: Performs log2 transformation, normalization, missing value imputation (QRILC).
   - `03_batch_correction_DE.R`: Applies SVA + batch correction and runs limma DE analysis.
   - `04_visualizations.R`: PCA, volcano plots, rescue plots.
   - `05_GO_GSEA_analysis.R`: GO, GSEA.
   - `06_Heatmap_visualization:  heatmaps for selected pathways.

4. All plots and results are saved to `/results/` and `/figures/`.

---

## 💻 Software & Dependencies

- R v4.2.2
- Key packages:
  - `limma`, `sva`, `clusterProfiler`, `GOSemSim`, `rrvgo`
  - `ggplot2`, `pheatmap`, `ComplexHeatmap`
  - `imputeLCMD`, `enrichR`, `biomaRt`, `openxlsx`
- Full package list: see `sessionInfo.txt`

---

## 📜 License and Citation

> This repository accompanies the manuscript:  
> **DDHD2 provides a flux of saturated fatty acids for neuronal energy and functions**  
> (Author list and DOI will be added upon publication)

If using this code, please credit the original authors and link back to this repository.

---

## 📬 Contact

For questions or suggestions, contact:
- Yih Tyng Bong — [rita.bong@helsinki.fi]
- OR open an issue on this GitHub repo
