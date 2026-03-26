#!/bin/bash

############################################
# mtDNA Variant Calling Pipeline v2.0
# CDFD — Mitochondrial Genome Analysis
#
# Fixes vs previous version:
#  - Proper step timing + tee logging
#  - RG tags in BWA mem -R (no redundant GATK step)
#  - SAM removed after sort to save disk
#  - AF filter uses Mutect2 native AF tag, falls back to AD
#  - SNP-only filter re-added after AF≥5%
#  - MITOMAP AWK column detection robust (case-insensitive)
#  - MITOMAP ALT iteration order fixed (for i=1..n_alt)
#  - MITOMASTER: retry, timeout, empty + HTML checks
#  - MITOMASTER column mapping corrected (col10=haplogroup)
#  - HaploCheck output: multi-path search, quote stripping
#  - Mean coverage: numeric rows only
#  - Duplicate DONE block removed
#  - set -euo pipefail (catches unbound vars)
############################################

set -euo pipefail

############################################
# INPUT ARGUMENTS
############################################

if [ "$#" -ne 4 ]; then
    echo "Usage: bash mtdna_pipeline.sh <SAMPLE_NAME> <R1.fastq.gz> <R2.fastq.gz> <REFERENCE.fasta>"
    echo ""
    echo "  SAMPLE_NAME  : Unique identifier (letters, numbers, _ or -)"
    echo "  R1.fastq.gz  : Forward reads"
    echo "  R2.fastq.gz  : Reverse reads"
    echo "  REFERENCE    : Path to mtDNA reference FASTA (rCRS recommended)"
    exit 1
fi

SAMPLE=$1
R1=$2
R2=$3
REFERENCE=$4

# Validate inputs exist
for INPUT_FILE in "$R1" "$R2" "$REFERENCE"; do
    if [[ ! -f "$INPUT_FILE" ]]; then
        echo "ERROR: Input file not found: $INPUT_FILE"
        exit 1
    fi
done

if [[ ! "$SAMPLE" =~ ^[A-Za-z0-9_-]+$ ]]; then
    echo "ERROR: SAMPLE_NAME must contain only letters, numbers, _ or -"
    exit 1
fi

############################################
# DIRECTORY STRUCTURE
############################################

BASE_DIR=$(pwd)
OUT_DIR="$BASE_DIR/output_${SAMPLE}"
LOG_DIR="$OUT_DIR/logs"

mkdir -p \
    "$OUT_DIR/01_raw_fastqc" \
    "$OUT_DIR/02_trimmed" \
    "$OUT_DIR/03_mapping" \
    "$OUT_DIR/04_variant_calling" \
    "$OUT_DIR/05_bam_stats" \
    "$OUT_DIR/06_coverage" \
    "$OUT_DIR/07_haplocheck" \
    "$OUT_DIR/09_mitomaster" \
    "$OUT_DIR/09_mitomap" \
    "$OUT_DIR/10_final_report_data" \
    "$LOG_DIR"

############################################
# RUNTIME LOG + STEP WRAPPER
# tee ensures progress lines reach stdout
# (web app reads stdout for progress tracking)
############################################

RUNTIME_LOG="$LOG_DIR/pipeline_runtime.log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$RUNTIME_LOG"
}

run_command() {
    local CMD="$1"
    local LOGFILE="$2"
    local STEP_NAME="${3:-step}"
    local START_TIME=$SECONDS

    log "START: $STEP_NAME"

    if ! bash -c "$CMD" > "$LOGFILE" 2>&1; then
        log "FAILED: $STEP_NAME — see $LOGFILE"
        tail -20 "$LOGFILE" >&2
        exit 1
    fi

    local ELAPSED=$(( SECONDS - START_TIME ))
    log "DONE:  $STEP_NAME (${ELAPSED}s)"
}

log "========================================="
log "mtDNA Pipeline Started"
log "Sample    : $SAMPLE"
log "R1        : $R1"
log "R2        : $R2"
log "Reference : $REFERENCE"
log "Output    : $OUT_DIR"
log "========================================="

############################################
# REFERENCE INDEXING
############################################

log "Checking reference index files..."

if [[ ! -f "${REFERENCE}.bwt" && ! -f "${REFERENCE}.bwt.2bit.64" ]]; then
    log "Building BWA index..."
    bwa index "$REFERENCE" >> "$LOG_DIR/00_bwa_index.log" 2>&1
fi

if [[ ! -f "${REFERENCE}.fai" ]]; then
    log "Building FASTA index..."
    samtools faidx "$REFERENCE" >> "$LOG_DIR/00_faidx.log" 2>&1
fi

DICT_FILE="${REFERENCE%.*}.dict"
if [[ ! -f "$DICT_FILE" ]]; then
    log "Building sequence dictionary..."
    gatk CreateSequenceDictionary \
        -R "$REFERENCE" \
        -O "$DICT_FILE" \
        >> "$LOG_DIR/00_dict.log" 2>&1
fi

log "Reference ready."

############################################
# STEP 1: FASTQC
############################################

run_command \
    "fastqc \"$R1\" \"$R2\" -o \"$OUT_DIR/01_raw_fastqc\" -t 4" \
    "$LOG_DIR/01_fastqc.log" \
    "FastQC"

############################################
# STEP 2: TRIMMING
############################################

run_command \
    "trim_galore \
        --paired \"$R1\" \"$R2\" \
        --quality 30 \
        --length 30 \
        --retain_unpaired \
        --cores 4 \
        --output_dir \"$OUT_DIR/02_trimmed\"" \
    "$LOG_DIR/02_trimming.log" \
    "Trim Galore"

TRIM_R1=$(find "$OUT_DIR/02_trimmed" -name "*_val_1.fq*" | sort | head -n 1)
TRIM_R2=$(find "$OUT_DIR/02_trimmed" -name "*_val_2.fq*" | sort | head -n 1)

if [[ -z "$TRIM_R1" || -z "$TRIM_R2" ]]; then
    log "ERROR: Trimmed output files not found in $OUT_DIR/02_trimmed"
    ls "$OUT_DIR/02_trimmed/" >> "$RUNTIME_LOG"
    exit 1
fi

mv "$TRIM_R1" "$OUT_DIR/02_trimmed/${SAMPLE}_R1_trimmed.fastq.gz"
mv "$TRIM_R2" "$OUT_DIR/02_trimmed/${SAMPLE}_R2_trimmed.fastq.gz"
TRIM_R1="$OUT_DIR/02_trimmed/${SAMPLE}_R1_trimmed.fastq.gz"
TRIM_R2="$OUT_DIR/02_trimmed/${SAMPLE}_R2_trimmed.fastq.gz"

############################################
# STEP 3: ALIGNMENT
# RG tags embedded in BWA mem -R directly
# Eliminates the separate AddOrReplaceReadGroups step
############################################

SAM_FILE="$OUT_DIR/03_mapping/${SAMPLE}.sam"

run_command \
    "bwa mem \
        -t 16 \
        -R \"@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:lib1\tPU:unit1\" \
        \"$REFERENCE\" \"$TRIM_R1\" \"$TRIM_R2\" \
        > \"$SAM_FILE\"" \
    "$LOG_DIR/03_alignment.log" \
    "BWA MEM"

############################################
# STEP 4: SORT + INDEX BAM
############################################

SORTED_BAM="$OUT_DIR/03_mapping/${SAMPLE}_sorted.bam"

run_command \
    "samtools view -Sb \"$SAM_FILE\" \
        | samtools sort -@ 8 -o \"$SORTED_BAM\"" \
    "$LOG_DIR/04_sorting.log" \
    "Sort BAM"

run_command \
    "samtools index \"$SORTED_BAM\"" \
    "$LOG_DIR/05_index_sorted.log" \
    "Index Sorted BAM"

# Remove SAM to save disk space
rm -f "$SAM_FILE"
log "SAM removed to save disk space."

############################################
# STEP 5: MARK DUPLICATES
# (No AddOrReplaceReadGroups needed — RGs already in BAM)
############################################

DEDUP_BAM="$OUT_DIR/03_mapping/${SAMPLE}_dedup.bam"
METRICS_FILE="$OUT_DIR/03_mapping/${SAMPLE}_dup_metrics.txt"

run_command \
    "gatk MarkDuplicates \
        -I \"$SORTED_BAM\" \
        -O \"$DEDUP_BAM\" \
        -M \"$METRICS_FILE\" \
        --VALIDATION_STRINGENCY LENIENT \
        --TMP_DIR \"$OUT_DIR/03_mapping\"" \
    "$LOG_DIR/06_markduplicates.log" \
    "Mark Duplicates"

run_command \
    "samtools index \"$DEDUP_BAM\"" \
    "$LOG_DIR/07_index_dedup.log" \
    "Index Dedup BAM"

############################################
# STEP 6: BAM QC STATISTICS
############################################

run_command \
    "samtools flagstat \"$DEDUP_BAM\" \
        > \"$OUT_DIR/05_bam_stats/${SAMPLE}_flagstat.txt\"" \
    "$LOG_DIR/08_flagstat.log" \
    "Flagstat"

run_command \
    "samtools stats \"$DEDUP_BAM\" \
        > \"$OUT_DIR/05_bam_stats/${SAMPLE}_stats.txt\"" \
    "$LOG_DIR/09_stats.log" \
    "BAM Stats"

run_command \
    "samtools idxstats \"$DEDUP_BAM\" \
        > \"$OUT_DIR/05_bam_stats/${SAMPLE}_idxstats.txt\"" \
    "$LOG_DIR/10_idxstats.log" \
    "idxstats"

############################################
# STEP 7: COVERAGE DEPTH
############################################

COVERAGE_FILE="$OUT_DIR/06_coverage/${SAMPLE}_coverage.txt"
MEAN_COV_FILE="$OUT_DIR/06_coverage/${SAMPLE}_mean_coverage.txt"

run_command \
    "samtools depth -a \"$DEDUP_BAM\" > \"$COVERAGE_FILE\"" \
    "$LOG_DIR/11_coverage.log" \
    "Coverage Depth"

log "Calculating mean coverage..."
awk '$3 ~ /^[0-9]+$/ {sum+=$3; n++} END {if(n>0) printf "%.2f\n",sum/n; else print "0"}' \
    "$COVERAGE_FILE" > "$MEAN_COV_FILE"

MEAN_COV=$(cat "$MEAN_COV_FILE")
log "Mean coverage: ${MEAN_COV}x"

if awk "BEGIN{exit !($MEAN_COV < 100)}"; then
    log "WARNING: Mean coverage (${MEAN_COV}x) below recommended 100x"
fi

############################################
# STEP 8: MUTECT2 (mtDNA Mode)
############################################

UNFILTERED_VCF="$OUT_DIR/04_variant_calling/${SAMPLE}_unfiltered.vcf"

run_command \
    "gatk Mutect2 \
        -R \"$REFERENCE\" \
        -I \"$DEDUP_BAM\" \
        --mitochondria-mode \
        --annotation StrandBiasBySample \
        -O \"$UNFILTERED_VCF\"" \
    "$LOG_DIR/12_mutect2.log" \
    "Mutect2"

############################################
# STEP 9: FILTER VARIANTS
############################################

FILTERED_VCF="$OUT_DIR/04_variant_calling/${SAMPLE}_filtered.vcf"

run_command \
    "gatk FilterMutectCalls \
        -R \"$REFERENCE\" \
        --mitochondria-mode \
        -V \"$UNFILTERED_VCF\" \
        -O \"$FILTERED_VCF\"" \
    "$LOG_DIR/13_filtering.log" \
    "FilterMutectCalls"

############################################
# STEP 10: SPLIT MULTIALLELIC SITES
############################################

log "Splitting multiallelic variants..."

SPLIT_VCF="$OUT_DIR/04_variant_calling/${SAMPLE}_split.vcf"

bcftools norm \
    -f "$REFERENCE" \
    -m -any \
    "$FILTERED_VCF" \
    -o "$SPLIT_VCF" \
    2>> "$LOG_DIR/14_split.log"

############################################
# STEP 11: AF ≥ 5% FILTER
# Uses Mutect2 native AF FORMAT tag (preferred).
# Falls back to AD-based calculation if AF absent.
############################################

log "Applying AF >= 5% filter..."

AF5_VCF="$OUT_DIR/04_variant_calling/${SAMPLE}_AF5.vcf"

bcftools view -h "$SPLIT_VCF" > "$AF5_VCF"

bcftools view -H "$SPLIT_VCF" | \
awk -v OFS="\t" '
BEGIN { FS="\t" }
{
    n = split($9, fmt, ":")
    delete fmap
    for (i=1; i<=n; i++) fmap[fmt[i]] = i

    split($10, smp, ":")
    af = -1

    if ("AF" in fmap) {
        af = smp[fmap["AF"]] + 0
    } else if ("AD" in fmap) {
        split(smp[fmap["AD"]], ad, ",")
        ref_c = ad[1] + 0
        alt_c = ad[2] + 0
        total = ref_c + alt_c
        if (total > 0) af = alt_c / total
    }

    if (af >= 0.05) print
}
' >> "$AF5_VCF"

log "AF filtering done."

############################################
# STEP 12: SNP-ONLY FILTER
############################################

log "Filtering to SNPs only..."

SNP_VCF="$OUT_DIR/04_variant_calling/${SAMPLE}_AF5_snps.vcf"

bcftools view --type snps "$AF5_VCF" -o "$SNP_VCF"

if [[ ! -s "$SNP_VCF" ]]; then
    log "WARNING: No SNPs passed filters. Creating header-only VCF."
    bcftools view -h "$AF5_VCF" > "$SNP_VCF"
fi

SNP_COUNT=$(bcftools view -H "$SNP_VCF" | wc -l)
log "SNPs after filtering: $SNP_COUNT"

############################################
# STEP 13: EXTRACT HETEROPLASMY
############################################

log "Extracting heteroplasmy metrics..."

HET_FILE="$OUT_DIR/04_variant_calling/${SAMPLE}_heteroplasmy.txt"

bcftools query \
    -f '%POS\t%REF\t%ALT\t[%AD]\t[%AF]\n' \
    "$SNP_VCF" | \
awk -v OFS="\t" '
{
    pos    = $1
    ref_b  = $2
    alt_b  = $3
    split($4, ad, ",")
    ref_c  = ad[1] + 0
    alt_c  = ad[2] + 0
    af_tag = $5 + 0
    total  = ref_c + alt_c

    if (total > 0) {
        het = alt_c / total
    } else if (af_tag > 0) {
        het = af_tag; ref_c = 0; alt_c = 0
    } else {
        next
    }

    printf "%s\t%s\t%s\t%s\t%d\t%d\t%.5f\n", \
        pos, ref_b, alt_b, pos":"ref_b":"alt_b, ref_c, alt_c, het
}
' > "$HET_FILE"

log "Heteroplasmy extraction done: $(wc -l < "$HET_FILE") variants"

############################################
# STEP 14: MITOMAP ANNOTATION
############################################

log "Running MITOMAP annotation..."

MITOMAP_DB="$BASE_DIR/database/mitomap_foswiki.csv"
MITOMAP_OUT="$OUT_DIR/09_mitomap/${SAMPLE}_mitomap_annotated.tsv"

if [[ ! -f "$MITOMAP_DB" ]]; then
    log "WARNING: MITOMAP database not found: $MITOMAP_DB"
    log "         Download mitomap_foswiki.csv and place in $BASE_DIR/database/"
    log "         Creating placeholder annotation and continuing..."
    printf "POS\tREF\tALT\tAF\tPathogenicity\tDisease\n" > "$MITOMAP_OUT"
    bcftools query -f '%POS\t%REF\t%ALT\tNA\tUncertain\tNA\n' \
        "$SNP_VCF" >> "$MITOMAP_OUT"
else
    # Clean CSV → TSV
    sed 's/"//g' "$MITOMAP_DB" | tr ',' '\t' \
        > "$OUT_DIR/09_mitomap/mitomap.tsv"

    # Build position hash with robust column detection
    awk -F'\t' '
    NR==1 {
        for (i=1; i<=NF; i++) {
            col = tolower($i)
            gsub(/[^a-z]/, "", col)
            if (col=="pos"      || col=="position")   pos_col = i
            if (col=="ref"      || col=="refallele")  ref_col = i
            if (col=="alt"      || col=="altallele")  alt_col = i
            if (col=="disease"  || col=="diseases")   dis_col = i
            if (col=="status"   || col=="reportedstatus") sta_col = i
        }
        if (!pos_col || !ref_col || !alt_col) {
            print "ERROR: Required columns (pos/ref/alt) not found in MITOMAP header" \
                  > "/dev/stderr"
            exit 1
        }
        if (!dis_col) dis_col = 0
        if (!sta_col) sta_col = 0
        next
    }
    NF < 3 { next }
    tolower($pos_col) ~ /^pos/ { next }   # skip header-repeat rows
    {
        dis = dis_col ? $dis_col : "NA"
        sta = sta_col ? $sta_col : "unknown"
        key = $pos_col "_" $ref_col "_" $alt_col
        hash[key] = dis "\t" sta
    }
    END {
        for (k in hash) print k "\t" hash[k]
    }
    ' "$OUT_DIR/09_mitomap/mitomap.tsv" \
        > "$OUT_DIR/09_mitomap/mitomap_hash.tsv"

    # Annotate VCF
    awk -F'\t' '
    BEGIN {
        OFS="\t"
        print "POS","REF","ALT","AF","Pathogenicity","Disease"
    }
    FNR==NR {
        mito[$1] = $2 "\t" $3
        next
    }
    /^#/ { next }
    {
        pos  = $2
        ref  = $4
        info = $8

        af = "NA"
        n  = split(info, fields, ";")
        for (i=1; i<=n; i++) {
            if (fields[i] ~ /^AF=/) {
                split(fields[i], kv, "=")
                af = kv[2]
            }
        }

        # Process ALTs in order (i=1..n_alt, not arbitrary iteration)
        n_alt = split($5, alts, ",")
        for (i=1; i<=n_alt; i++) {
            alt = alts[i]
            key = pos "_" ref "_" alt

            if (key in mito) {
                split(mito[key], m, "\t")
                disease = m[1]
                status  = tolower(m[2])
            } else {
                disease = "NA"
                status  = "unknown"
            }

            if      (status ~ /confirmed|pathogenic/) path = "Pathogenic"
            else if (status ~ /reported/)             path = "Likely_Pathogenic"
            else if (status ~ /polymorphism|benign/)  path = "Benign"
            else                                      path = "Uncertain"

            print pos, ref, alt, af, path, disease
        }
    }
    ' "$OUT_DIR/09_mitomap/mitomap_hash.tsv" \
      "$SNP_VCF" \
      > "$MITOMAP_OUT"

    rm -f "$OUT_DIR/09_mitomap/mitomap.tsv" \
          "$OUT_DIR/09_mitomap/mitomap_hash.tsv"
fi

log "MITOMAP annotation done: $MITOMAP_OUT"

############################################
# STEP 15: HAPLOCHECK
############################################

log "Running HaploCheck..."

HAPLO_OUTDIR="$OUT_DIR/07_haplocheck"
HAPLO_PREFIX="$HAPLO_OUTDIR/${SAMPLE}"

haplocheck \
    --out "$HAPLO_PREFIX" \
    "$SNP_VCF" \
    >> "$LOG_DIR/15_haplocheck.log" 2>&1

# Locate output file — HaploCheck may write to different paths
HAPLO_RAW=""
for CANDIDATE in \
    "${HAPLO_PREFIX}" \
    "${HAPLO_PREFIX}.txt" \
    "${HAPLO_OUTDIR}/haplocheck.txt" \
    "${HAPLO_OUTDIR}/report/haplocheck.txt"; do
    if [[ -f "$CANDIDATE" ]]; then
        HAPLO_RAW="$CANDIDATE"
        break
    fi
done

if [[ -z "$HAPLO_RAW" ]]; then
    HAPLO_RAW=$(find "$HAPLO_OUTDIR" -name "*.txt" -not -name "*.log" 2>/dev/null | head -n 1)
fi

HAPLO_TSV="$HAPLO_OUTDIR/${SAMPLE}_haplocheck.tsv"

if [[ -z "$HAPLO_RAW" || ! -f "$HAPLO_RAW" ]]; then
    log "WARNING: HaploCheck output not found — using placeholder"
    printf "Sample\tHaplogroup\tContamination\tQuality\n" > "$HAPLO_TSV"
    printf "%s\tUnknown\tNA\tNA\n" "$SAMPLE"            >> "$HAPLO_TSV"
else
    # Skip comment/header lines; strip quotes; take first data row
    awk -F'\t' '
    /^#/ { next }
    /^Sample/ { next }
    NF >= 3 {
        hap  = $2; gsub(/"/, "", hap)
        cont = $3; gsub(/"/, "", cont)
        qual = (NF>=4) ? $4 : "NA"; gsub(/"/, "", qual)
        print "Sample\tHaplogroup\tContamination\tQuality"
        print $1 "\t" hap "\t" cont "\t" qual
        exit
    }
    ' "$HAPLO_RAW" > "$HAPLO_TSV"

    if [[ $(wc -l < "$HAPLO_TSV") -lt 2 ]]; then
        log "WARNING: HaploCheck TSV empty — using placeholder"
        printf "Sample\tHaplogroup\tContamination\tQuality\n" > "$HAPLO_TSV"
        printf "%s\tUnknown\tNA\tNA\n" "$SAMPLE"            >> "$HAPLO_TSV"
    fi
fi

HAPLOGROUP=$(awk -F'\t' 'NR==2{print $2}' "$HAPLO_TSV")
CONTAM=$(awk -F'\t'     'NR==2{print $3}' "$HAPLO_TSV")
log "Haplogroup: $HAPLOGROUP  |  Contamination: $CONTAM"

############################################
# STEP 16: MITOMASTER INPUT PREPARATION
############################################

MITO_DIR="$OUT_DIR/09_mitomaster"
MITO_INPUT="$MITO_DIR/${SAMPLE}_mitomaster_input.txt"

log "Preparing MITOMASTER input..."

printf "sample\tpos\tref\tvar\n" > "$MITO_INPUT"
bcftools query \
    -f "${SAMPLE}\t%POS\t%REF\t%ALT\n" \
    "$SNP_VCF" >> "$MITO_INPUT"

INPUT_VARS=$(( $(wc -l < "$MITO_INPUT") - 1 ))
log "MITOMASTER input: $INPUT_VARS variants"

############################################
# STEP 17: MITOMASTER API
############################################

MITO_RAW="$MITO_DIR/${SAMPLE}_mitomaster_raw.txt"
MITO_TSV="$MITO_DIR/${SAMPLE}_mitomaster.tsv"
MITO_FINAL="$MITO_DIR/${SAMPLE}_mitomaster_final.tsv"

log "Submitting to MITOMASTER..."

MITO_FAILED=false

HTTP_STATUS=$(curl \
    --silent --show-error \
    --retry 3 --retry-delay 15 --retry-connrefused \
    --max-time 180 --connect-timeout 30 \
    -w "%{http_code}" \
    -H "User-Agent: Mozilla/5.0 (MITOCLIN/2.0; Research)" \
    -F "file=@${MITO_INPUT}" \
    -F "fileType=snvlist" \
    -F "output=detail" \
    https://www.mitomap.org/mitomaster/websrvc.cgi \
    -o "$MITO_RAW" 2>>"$LOG_DIR/17_mitomaster.log") || true

log "MITOMASTER HTTP status: $HTTP_STATUS"

if   [[ "$HTTP_STATUS" != "200" ]];         then MITO_FAILED=true; log "WARNING: HTTP $HTTP_STATUS"
elif [[ ! -s "$MITO_RAW" ]];               then MITO_FAILED=true; log "WARNING: Empty response"
elif grep -qi "<html>" "$MITO_RAW" 2>/dev/null; then MITO_FAILED=true; log "WARNING: HTML error page returned"
fi

if [[ "$MITO_FAILED" == "true" ]]; then
    log "Creating placeholder MITOMASTER output..."
    printf "POS\tREF\tALT\tMutation\tProtein_Change\tGene\tDisease\tHaplogroup\n" \
        > "$MITO_FINAL"
    bcftools query -f '%POS\t%REF\t%ALT\n' "$SNP_VCF" | \
    awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$2,$3,"NA","NA","NA","NA","NA"}' \
        >> "$MITO_FINAL"
else
    log "Parsing MITOMASTER response..."

    sed 's/<br>/|/g' "$MITO_RAW" | sed 's/,/\t/g' > "$MITO_TSV"

    # MITOMASTER detail output (1-indexed after CSV→TSV):
    #  1=line  2=pos  3=ref  4=alt  5=mutation_name
    #  6=homoplasmy_freq  7=gene  8=codon_AA  9=codon_change
    # 10=haplogroup  11=tRNA_annot  12=disease
    awk -F'\t' 'BEGIN{OFS="\t"}
    NR==1 {
        print "POS","REF","ALT","Mutation","Protein_Change","Gene","Disease","Haplogroup"
        next
    }
    NF < 5 { next }
    $2 ~ /^[Pp]os/ { next }
    $2 == "" { next }
    {
        print $2, $3, $4, $5, $8, $7, $12, $10
    }
    ' "$MITO_TSV" > "$MITO_FINAL"

    log "MITOMASTER parsed: $(wc -l < "$MITO_FINAL") lines (including header)"
fi

############################################
# STEP FINAL: MERGE ALL FILES
############################################

log "Generating final merged report..."

FINAL_DIR="$OUT_DIR/10_final_report_data"
OUTFILE="$FINAL_DIR/${SAMPLE}_final_report.csv"

awk -F'\t' \
    -v OFS="," \
    -v sample="$SAMPLE" \
    -v cont="$CONTAM" \
    -v mean="$MEAN_COV" \
    -v hetfile="$HET_FILE" \
    -v covfile="$COVERAGE_FILE" \
    -v mitofile="$MITO_FINAL" \
    -v mitomapfile="$MITOMAP_OUT" \
'
function make_key(pos, ref, alt,    k) {
    k = pos ":" ref ":" alt
    gsub(/[[:space:]]/, "", k)
    return toupper(k)
}

FILENAME == hetfile {
    k = make_key($1, $2, $3)
    refc[k] = $5; altc[k] = $6; hetp[k] = $7
    next
}

FILENAME == covfile {
    cov_at[$2] = $3
    next
}

FILENAME == mitomapfile {
    if (FNR == 1) next
    k = make_key($1, $2, $3)
    patho[k] = $5
    next
}

FILENAME == mitofile {
    if (FNR == 1) {
        print "Sample,Position,Ref,Alt,Mutation_Type,Protein_Change,Gene," \
              "Disease,Pathogenicity,Haplogroup,Contamination,Mean_Coverage," \
              "Coverage_at_Position,Ref_Count,Alt_Count,Heteroplasmy_pct"
        next
    }

    pos = $1; ref = $2; alt = $3
    mut = $4; prot = $5; gene = $6; disease = $7; haplo = $8

    if (pos == "" || pos ~ /^[Pp]os/) next

    k    = make_key(pos, ref, alt)
    rc   = (k in refc)    ? refc[k]    : "NA"
    ac   = (k in altc)    ? altc[k]    : "NA"
    hp   = (k in hetp)    ? hetp[k]    : "NA"
    cov  = (pos in cov_at) ? cov_at[pos] : "NA"
    path = (k in patho)   ? patho[k]   : "NA"

    print sample, pos, ref, alt, mut, prot, gene, disease, \
          path, haplo, cont, mean, cov, rc, ac, hp
}
' "$HET_FILE" "$COVERAGE_FILE" "$MITOMAP_OUT" "$MITO_FINAL" \
> "$OUTFILE"

REPORT_LINES=$(( $(wc -l < "$OUTFILE") - 1 ))
log "Final report: $REPORT_LINES variants"

############################################
# PIPELINE SUMMARY
############################################

SUMMARY="$FINAL_DIR/${SAMPLE}_summary.txt"
cat > "$SUMMARY" << SUMMARY_EOF
=========================================
 mtDNA Pipeline Summary — MITOCLIN v2.0
=========================================
Sample        : $SAMPLE
Date          : $(date '+%Y-%m-%d %H:%M:%S')
Reference     : $REFERENCE
Output Dir    : $OUT_DIR
-----------------------------------------
Mean Coverage : ${MEAN_COV}x
Haplogroup    : $HAPLOGROUP
Contamination : $CONTAM
SNPs filtered : $SNP_COUNT
Report Rows   : $REPORT_LINES
-----------------------------------------
Output Files:
  FastQC      : $OUT_DIR/01_raw_fastqc/
  Trimmed     : $OUT_DIR/02_trimmed/
  BAM         : $DEDUP_BAM
  SNP VCF     : $SNP_VCF
  Heteroplasmy: $HET_FILE
  MITOMAP     : $MITOMAP_OUT
  HaploCheck  : $HAPLO_TSV
  MITOMASTER  : $MITO_FINAL
  FINAL REPORT: $OUTFILE
=========================================
SUMMARY_EOF

cat "$SUMMARY"
log "Pipeline completed successfully."
log "Final report: $OUTFILE"
