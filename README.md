# Bulk_RNA_Seq_Tutorial

This repository provides a complete **bulk RNA-sequencing analysis tutorial**, based on the code and data used for a poster presentation at the **22nd AAOT (Asian Academy of Orofacial Pain and Temporomandibular Disorder) Scientific Meeting**, held on **November 2–3, 2024, in Taipei, Taiwan**.

The study compares **Erosive Oral Lichen Planus (EOLP)** and **Non-erosive Oral Lichen Planus (NEOLP)**, as originally described in the following publication:  
📄 [PubMed ID: 37555396](https://pubmed.ncbi.nlm.nih.gov/37555396/)


---

## 📥 How to Download This Repository

You can download this tutorial repository in two easy ways:

---

### ✅ **Option 1: Download as ZIP (for those who are not familar with using git clone command)**

1. Click the green **`Code`** button near the top of this page.
2. Select **`Download ZIP`**.
3. Extract the ZIP file on your computer.
4. Open the folder and start exploring!

> ✅ No Git or command line needed.

---

### 🧪 **Option 2: Clone with Git (for those who are familar with using git clone command)**

If you're familiar with Git and want to clone this repository:

```bash
git clone https://github.com/alexjoo-kang/Bulk_RNA_Seq_Tutorial.git
```

> This will create a local copy of the project on your machine that you can update using Git.

---


## 🔍 Overview

This repository is created with **generous permission** from **Prof. Man Seok Kim** and **Prof. Hye-Ji Park** to help people who are **new to bioinformatics** become familiar with **bulk RNA-seq analysis workflows**.  

The code here was used in a real research poster and is now publicly shared to make the learning curve smoother for beginners.

It includes:

- A real-world count matrix (`final_count_matrix.tsv`)
- Tool settings for preprocessing (e.g., STAR, HTSeq)
- A unified R script to generate:
  - Volcano plot
  - GO enrichment plot
  - Lollipop plot (based on fGSEA results)
  - Heatmap (z-score scaled)

---

## 📁 Repository Structure

```
Bulk_RNA_Seq_Tutorial/
│
├── final_count_matrix.tsv                         # Count matrix used for plotting
│
├── Special_Case/                                  # NOX5-specific example scripts
│   ├── fifth_RAAS_nox5_issue_resolved.R
│   └── fifth_RAAS_nox5_issue_unresolved.R
│
├── inflammation_reference/
│   ├── for_heatmap/                               # Gene sets for z-score heatmaps
│   │   ├── Adaptive_Extracellular_Mediated_Immunity_Genes.csv
│   │   ├── ISR_Genes.csv
│   │   ├── Innate_Immune_Genes.csv
│   │   ├── Mitochondrial_Innate_Immune_Genes.csv
│   │   ├── RAAS_Genes.csv
│   │   └── UPR_Genes.csv
│   └── for_lollipop/
│       └── Inflammation_Genes_List.csv            # Used for lollipop plot grouping
│
├── preprocessing/
│   ├── ReadLength_automation.R
│   ├── bulk_rna_seq_preprocessing_using_star_and_htseq.sh
│   └── preprocessing_tools_used_for_bulk_rna_seq_olp.pdf
│
├── resulting_figures/                             # Output plots from scripts
│   ├── example_heatmap_plot.pdf
│   ├── go_enrichment_analysis_plot_custom_size.pdf
│   ├── lollipop_plot_based_on_fgsea_results_us_legal_size.pdf
│   ├── volcano_plot_customized_size.pdf
│   └── special_case_pleaes_focus_on_NADPH_Oxidase_column_NOX5_Gene/
│       ├── NOX5_Issue_Resolved_Custom_Size.pdf
│       └── NOX5_Issue_Unresolved_Custom_Size.pdf
│
└── one_unified_R_scipt_for_project_use_this.R     # Main analysis script (volcano, heatmap, etc.)
```

---

## 🧪 Data Source & Preprocessing

- **Project**: GSE213346  
- **Samples**: SRR2170343 ~ SRR21570382  
- **Genome**: GENCODE v46 (GRCh38, Primary Assembly, May 2024)  
- **Annotations**: GENCODE v46 (GTF and GFF3)  
- **Tools Used**:
  - `STAR 2.7.11b`
  - `samtools 1.20`
  - `RSeQC 5.0.2` (`infer_experiment.py`)
  - `HTSeq 2.0.8`
  - `BEDOPS 2.4.41` (`gff2bed`)

> ⚠️ Note: tool versions were confirmed after the analysis. Results may vary slightly with different versions.

---

## 📊 How to Run the Analysis

Run the following script in R:

```r
source("one_unified_R_scipt_for_project_use_this.R")
```

This script produces:
- A volcano plot  
- A GO enrichment bar plot  
- A lollipop plot (based on fGSEA)
- A z-score heatmap (based on selected gene sets)

---

## 📂 Gene Set Used for `one_unified_R_scipt_for_project_use_this.R`

The default gene set used for heatmap in the script is:
```
inflammation_reference/for_heatmap/Adaptive_Extracellular_Mediated_Immunity_Genes.csv
```

To use other sets like:
- `Innate_Immune_Genes.csv`
- `RAAS_Genes.csv`
- `ISR_Genes.csv`, etc.

Edit the `read.csv()` line in the script, and update grouping variables and heatmap titles accordingly.

---

## ❗ NOX5 Duplication Issue in RAAS Genes

- The file `RAAS_Genes.csv` includes two ENSEMBL IDs for **NOX5**.
- This can result in duplicate heatmap columns and misleading visuals.
- Use `fifth_RAAS_nox5_issue_resolved.R` to properly deduplicate the gene before plotting.

---

## ❓ Frequently Asked Questions

### Q1: How do I use other gene sets for **heatmap** analysis?

**Answer**:  
The script defaults to `Adaptive_Extracellular_Mediated_Immunity_Genes.csv` from the `for_heatmap/` folder.  
To switch:
1. Change the input file path in the script.
2. Adjust `list_ref$pathway == "..."` logic based on that file’s pathway categories.
3. Update the heatmap `name` title accordingly.

> ✅ All heatmaps use files from `inflammation_reference/for_heatmap/`.

---

### Q2: Why use 13 NEOLP samples?

**Answer**:  
There were 13 EOLP samples available. To balance group sizes, 13 NEOLP samples were selected from 27 total NEOLP samples.  
This improves z-score comparison clarity, though different NEOLP sample selections may yield slightly different results.

---

### Q3: What’s the difference between the preprint and published gene lists?

**Answer**:  
This tutorial used gene sets from the **preprint** version:  
📄 [bioRxiv 2023](https://www.biorxiv.org/content/10.1101/2023.10.08.561395v3)

The **final publication** in **PNAS 2024** added:
- *SPINK5* under **RAAS/Vascular Inflammation**
- *Syndecan-family genes* under **RAAS**

📘 Official article: [PNAS DOI](https://www.pnas.org/doi/10.1073/pnas.2401968121)

Use the updated gene list if you want full alignment with the final version.

---

## 🙏 Acknowledgements

We gratefully acknowledge the following contributors:

- 👨‍🏫 **Prof. Man Seok Kim** (Co-PI)  
- 👩‍🏫 **Prof. Hye-Ji Park** (Co-PI)  
  > They supervised the research and poster presentation and **graciously permitted public sharing** of this repository and its contents to support the broader research and learning community.

- 👨‍💻 **Jisu Jeong**  
  > Contributed significantly to the code base and generously agreed to its public release.

---

## 📬 Contact

📧 Email: [alexkang1014@naver.com](mailto:alexkang1014@naver.com)

We hope this resource helps you start your **bioinformatics and bulk RNA-seq journey** with clarity and confidence.

Feel free to reach out with feedback or questions!
