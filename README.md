# 🧬 Mitoclin: Clinical mtDNA Analysis Pipeline

---

##  Introduction

Mitoclin is an integrated bioinformatics pipeline and web-based application designed for the analysis and clinical interpretation of mitochondrial DNA (mtDNA) sequencing data. It combines multiple computational tools into a single automated workflow, enabling efficient and accurate genomic analysis for clinical applications.

---

##  Problem Statement

Mitochondrial DNA analysis is challenging due to:

- Presence of heteroplasmy (mixed variant populations)
- Complex inheritance patterns
- Lack of standardized interpretation frameworks
- Requirement of multiple independent tools and manual integration

These challenges make mtDNA analysis time-consuming, error-prone, and difficult for clinicians and forensic experts.

---

##  What Mitoclin Does

Mitoclin solves these problems by providing:

- Automated end-to-end pipeline (FASTQ → Clinical Report)
- Integration of multiple bioinformatics tools
- Heteroplasmy estimation
- Contamination detection
- Haplogroup determintion 
- Functional and clinical annotation
- User-friendly report generation via Flask interface

---

##  Tools Used

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

## Installation & Setup 

### 1. Create Environment

```bash
conda create -n mitoclin python=3.10 -y
conda activate mitoclin
```

---

### 2. Install Required Tools

```bash
# Core bioinformatics tools
conda install -c bioconda bwa samtools bcftools fastqc trim-galore -y

# GATK (Mutect2)
conda install -c bioconda gatk4 -y

# HaploCheck (contamination analysis)
conda install -c bioconda haplocheck -y
```

---

### 3. Verify Installation

```bash
bwa
samtools
bcftools
fastqc
trim_galore
gatk
haplocheck
```

---

## Usage (Run Mitoclin Pipeline)

```bash
bash mtdna_pipeline.sh \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -r whole_genome_mtdna.fasta \
  -o output_sample
```

---

## Output

The pipeline generates:

* Quality Control Reports (FastQC)
* Trimmed Reads
* Aligned BAM files
* Variant Calls (VCF)
* Heteroplasmy Analysis
* Contamination Report (HaploCheck)
* Final Annotated Report (TSV)

---

