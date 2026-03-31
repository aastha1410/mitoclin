# 🧬 Mitoclin: Clinical mtDNA Analysis Pipeline

---

## 📌 Introduction

Mitoclin is an integrated bioinformatics pipeline and web-based application designed for the analysis and clinical interpretation of mitochondrial DNA (mtDNA) sequencing data. It combines multiple computational tools into a single automated workflow, enabling efficient and accurate genomic analysis for clinical applications.

---

## ❗ Problem Statement

Mitochondrial DNA analysis is challenging due to:

- Presence of heteroplasmy (mixed variant populations)
- Complex inheritance patterns
- Lack of standardized interpretation frameworks
- Requirement of multiple independent tools and manual integration

These challenges make mtDNA analysis time-consuming, error-prone, and difficult for clinicians and forensic experts.

---

## 🎯 What Mitoclin Does

Mitoclin solves these problems by providing:

- Automated end-to-end pipeline (FASTQ → Clinical Report)
- Integration of multiple bioinformatics tools
- Heteroplasmy estimation
- Contamination detection
- Haplogroup determintion 
- Functional and clinical annotation
- User-friendly report generation via Flask interface

---

## 🧰 Tools Used

- FastQC → Quality control
- Trim Galore → Read trimming
- BWA MEM → Sequence alignment+Readgroups
- SAMtools → BAM processing
- BCFtools → Variant handling
- GATK Mutect2 → Variant calling (mitochondrial mode)
- HaploCheck → Contamination detection
- MITOMASTER → Functional annotation
- Python → Data processing & visualization
- Flask → Web interface

---

## ⚙️ Installation

### 🔹 Step 1: Create Conda Environment

```bash
conda create -n mitoclin_env python=3.10 -y
conda activate mitoclin_env
