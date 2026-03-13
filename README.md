# 🧬 Mass Spectrometry Proteomics Analysis: Command-Line Workflow

[![Bash](https://img.shields.io/badge/Shell-Bash-4EAA25?style=for-the-badge&logo=gnu-bash&logoColor=white)](https://www.gnu.org/software/bash/)
[![R](https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white)](https://www.r-project.org/)
[![Conda](https://img.shields.io/badge/Conda-342B029.svg?&style=for-the-badge&logo=anaconda&logoColor=white)](https://docs.conda.io/en/latest/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=for-the-badge)](https://opensource.org/licenses/MIT)
[![GitHub stars](https://img.shields.io/github/stars/your-username/your-repo?style=for-the-badge&logo=github&logoColor=white)](https://github.com/your-username/your-repo/stargazers)

> A complete, reproducible, command-line workflow for label-free quantitative (LFQ) proteomics data analysis — from raw mass spectrometry files to statistically significant, differentially abundant proteins.

This pipeline uses best-in-class, open-source tools for each step, ensuring transparency and full control over the analysis.

**Pipeline:** `ProteoWizard` → `OpenMS` → `SearchGUI` → `PeptideShaker` → `R + limma`

---

## 📋 Table of Contents

1.  [Overview & Goal](#-overview--goal)
2.  [Prerequisites](#-prerequisites)
3.  [Installation](#-installation)
4.  [The Analysis Workflow](#-the-analysis-workflow)
5.  [Quick Start](#-quick-start)
6.  [Step-by-Step Protocol](#-step-by-step-protocol)
    *   [Step 0: Database Preparation](#-step-0-database-preparation-target--target-decoy-fasta)
    *   [Step 1: Data Conversion](#-step-1-data-conversion-raw--mzml)
    *   [Step 2: Protein Identification](#-step-2-protein-identification)
    *   [Step 3: Validation & Reporting](#-step-3-validation--reporting)
    *   [Step 4: Statistical Analysis & Visualization](#-step-4-statistical-analysis--visualization)
7.  [File Formats Explained](#-file-formats-explained)
8.  [Expected Results](#-expected-results)
9.  [Troubleshooting](#-troubleshooting)
10. [Contributing](#-contributing)
11. [License](#-license)

---

## 🎯 Overview & Goal

This protocol takes raw mass spectrometry files (e.g., `.raw`, `.wiff`) and a target protein FASTA, and delivers:

*   ✅ A validated list of identified proteins with strict **FDR control**.
*   📊 Quantitative protein abundances (LFQ/intensity).
*   📈 A list of **significantly differentially abundant proteins**.
*   📄 A **publication-quality volcano plot**.

---

## 🧰 Prerequisites

### Required Data
*   **Raw MS Files:** Thermo `.raw`, SCIEX `.wiff`, Bruker `.d`, etc.
*   **Protein Database:** Target protein database in FASTA format (e.g., UniProt proteome for your organism).

### Recommended Hardware
*   **RAM:** ≥16 GB recommended (Java tools are memory-intensive).
*   **CPU:** Multi-core processor for faster processing.

---

## 💻 Installation

Create and activate a dedicated conda environment. This is the **strongly recommended** method for managing dependencies and ensuring reproducibility.

```bash
# Create the environment with all required tools
conda create -n proteomics-cli -c bioconda -c conda-forge \
    proteowizard openms searchgui peptideshaker \
    r-base r-essentials r-limma r-ggplot2 r-readr r-dplyr

# Activate the environment
conda activate proteomics-cli
```

> **Note:** Replace `X.Y.Z` in the Java commands below with the actual version numbers installed. You can find these in your conda environment path (e.g., `~/miniconda3/envs/proteomics-cli/share/searchgui-*`).

---

## 📊 The Analysis Workflow

```mermaid
graph TD
    A[Raw Data Files<br>(.raw, .wiff)] -->|msconvert| B[Open Data Files<br>(.mzML)]
    D[Target FASTA<br>(.fasta)] -->|DecoyDatabase| E[Target-Decoy FASTA<br>(.fasta)]
    B -->|SearchGUI| C[Identification Results<br>(.mzid files)]
    E --> C
    C -->|PeptideShaker| F[Protein Report<br>(protein_reports.tsv)]
    F -->|R / limma| G[Significant Proteins<br>(significant_proteins.csv)]
    G --> H[Volcano Plot<br>(volcano_plot.png)]
```

---

## 🚀 Quick Start

For users familiar with the process. A detailed, step-by-step guide follows below.

```bash
# 1. Prepare your target-decoy database (one-time step)
DecoyDatabase -in uniprot_proteome.fasta -out database.fasta -decoy_string DECOY_ -decoy_string_position prefix

# 2. Convert all raw files to mzML
mkdir -p mzml_output
msconvert *.raw --mzML --filter "peakPicking true 1-" -o mzml_output/

# 3. Run the identification, validation, and statistical analysis
# (Follow the detailed commands in the Step-by-Step Protocol below)
```

---

## 📝 Step-by-Step Protocol

### Step 0: Database Preparation (Target → Target-Decoy FASTA)

**Purpose:** To create a database with both real (target) and fake (decoy) sequences. This is essential for accurately estimating the False Discovery Rate (FDR).

```bash
DecoyDatabase -in target_proteins.fasta -out database.fasta \
    -decoy_string DECOY_ -decoy_string_position prefix
```

### Step 1: Data Conversion (`.raw` → `.mzML`)

**Purpose:** To convert vendor-specific raw files into the open, standard `.mzML` format. Peak picking is performed to make the data suitable for search engines.

```bash
mkdir -p mzml_output
msconvert *.raw --mzML --filter "peakPicking true 1-" -o mzml_output/
# For a single file: msconvert sample.raw --mzML --filter "peakPicking true 1-" -o mzml_output/
```

### Step 2: Protein Identification

**Purpose:** To match the experimental MS/MS spectra against the theoretical spectra from your database to identify peptides.

#### 2a. Create Parameters File
This file defines the search rules for the engines.

```bash
java -cp SearchGUI-X.Y.Z.jar eu.isas.searchgui.cmd.IdentificationParametersCLI \
    -out parameters.par \
    -enzyme Trypsin \
    -enzyme_specificity C-terminal \
    -missed_cleavages 2 \
    -instrument Q_Exactive \
    -precursor_tolerance 10.0 \
    -precursor_tolerance_unit PPM \
    -fragment_tolerance 0.02 \
    -fragment_tolerance_unit DA \
    -fixed_modifications Carbamidomethyl:C \
    -variable_modifications Oxidation:M
```

#### 2b. Run Search
This command runs the search engines (e.g., MS-GF+, X!Tandem) using the parameters.

```bash
java -Xmx16G -cp SearchGUI-X.Y.Z.jar eu.isas.searchgui.cmd.SearchGUICLI \
    -spectrum_files mzml_output/*.mzML \
    -fasta_file database.fasta \
    -parameters_file parameters.par \
    -output_folder searchgui_output/
```

### Step 3: Validation & Reporting

**Purpose:** To combine all search results, validate them, and infer proteins. PeptideShaker calculates FDR and provides a final, consolidated report.

```bash
java -Xmx16G -cp PeptideShaker-X.Y.Z.jar eu.isas.peptideshaker.cmd.PeptideShakerCLI \
    -project MyProject \
    -experiment MyExperiment \
    -identification_files searchgui_output/*.mzid \
    -spectrum_files mzml_output/*.mzML \
    -fasta_file database.fasta \
    -output_folder peptideshaker_output/
```
> **Main Output:** `peptideshaker_output/protein_reports.tsv`

### Step 4: Statistical Analysis & Visualization

**Purpose:** To identify proteins that are significantly different between your experimental groups and visualize the results.

**Run the script:**

```bash
Rscript analyze.R
```

---

## 📎 File Formats Explained

| Extension | Generated By | Description |
| :--- | :--- | :--- |
| `.raw/.wiff` | Instrument | Vendor-specific, binary raw data. |
| `.mzML` | `msconvert` | Open standard MS data format (XML + binary). |
| `.fasta` | UniProt / User | Protein sequences (target + decoy). |
| `.mzid` | `SearchGUI` | Peptide-spectrum match results (XML). |
| `.tsv` | `PeptideShaker` | Protein-level report with quantification (Tab-Separated). |
| `.csv` | `R script` | Differential analysis results (Comma-Separated). |
| `.png` | `R ggplot2` | Volcano plot (Portable Network Graphics). |

---

## 🏆 Expected Results

*   `protein_reports.tsv`: A table containing **1,000–15,000 proteins** (dataset dependent).
*   `significant_proteins.csv`: A filtered table of **50–1,000 differentially abundant proteins** (using typical thresholds: FDR < 0.05, |log₂FC| > 1).
*   `volcano_plot.png`: A clear, publication-quality visualization of significant changes.

---

## 🛠️ Troubleshooting

| Issue | Likely Cause | Fix |
| :--- | :--- | :--- |
| Very few identifications | Wrong parameters / poor data | Double-check tolerances, enzyme, and modifications in `parameters.par`. QC raw data if possible. |
| PeptideShaker errors | File mismatch / no decoys | Ensure the **exact same** `.mzML` files are used for both SearchGUI and PeptideShaker. Verify `DECOY_` is in your FASTA headers. |
| "No quantitative columns in R" | Column names differ | Inspect `protein_reports.tsv` headers and update the `CONTROL_COLS`/`TREATMENT_COLS` variables in `analyze.R`. |
| Java out of memory error | Large dataset | Increase Java heap size: change `-Xmx16G` to `-Xmx32G` or higher. |

---

## 🤝 Contributing

Contributions are what make the scientific community great!

*   **Report bugs** or **suggest improvements** by opening an Issue.
*   **Submit enhancements** via a Pull Request.

---

## 📜 License

This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.
