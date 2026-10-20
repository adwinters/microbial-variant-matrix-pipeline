#!/bin/bash

hostname
date

# =========================================================
# Pipeline: Single reference, multiple genomes,
# joint variant calling and per-sample variant summaries
# =========================================================

set -euo pipefail

# ---------------------------
# Check arguments
# ---------------------------
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <reference.fna> <fastq_dir> [output_dir]"
    exit 1
fi

# ---------------------------
# 0. Load required modules
# ---------------------------
module purge
module load gnu7/7.3.0
module load bwa/0.7.17
module load samtools/1.21
module load bcftools

THREADS=16

# ---------------------------
# 1. Set directories
# ---------------------------
REF_FASTA=$1
FASTQ_DIR=$2
OUT_DIR=${3:-results}

mkdir -p "$OUT_DIR"

# ---------------------------
# 2. Reference genome
# ---------------------------
echo "Using reference: $REF_FASTA"

if [[ ! -f "$REF_FASTA" ]]; then
    echo "ERROR: Reference FASTA not found."
    exit 1
fi

bwa index "$REF_FASTA"
samtools faidx "$REF_FASTA"

# ---------------------------
# 3. List reference contigs
# ---------------------------
echo "Reference contigs:"
cut -f1 "${REF_FASTA}.fai"

# ---------------------------
# 4. Process FASTQ files
# ---------------------------
BAMS=()

for FQ1 in "$FASTQ_DIR"/*_R1_001.trim.fastq.gz
do
    SAMPLE=$(basename "$FQ1" _R1_001.trim.fastq.gz)
    FQ2="$FASTQ_DIR/${SAMPLE}_R2_001.trim.fastq.gz"

    if [[ ! -f "$FQ2" ]]; then
        echo "ERROR: Paired file for $SAMPLE not found."
        exit 1
    fi

    echo "Processing sample: $SAMPLE"

    bwa mem \
        -t "$THREADS" \
        "$REF_FASTA" \
        "$FQ1" \
        "$FQ2" | \
    samtools view \
        -@ "$THREADS" \
        -Sb \
        -o "$OUT_DIR/${SAMPLE}.unsorted.bam"

    samtools addreplacerg \
        -r "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
        -o "$OUT_DIR/${SAMPLE}.bam" \
        "$OUT_DIR/${SAMPLE}.unsorted.bam"

    samtools sort \
        -@ "$THREADS" \
        "$OUT_DIR/${SAMPLE}.bam" \
        -o "$OUT_DIR/${SAMPLE}.sorted.bam"

    samtools index "$OUT_DIR/${SAMPLE}.sorted.bam"

    BAMS+=("$OUT_DIR/${SAMPLE}.sorted.bam")

    rm \
        "$OUT_DIR/${SAMPLE}.bam" \
        "$OUT_DIR/${SAMPLE}.unsorted.bam"
done

# ---------------------------
# 5. Joint variant calling
# ---------------------------
echo "Calling variants jointly..."

bcftools mpileup \
    --threads "$THREADS" \
    -Ou \
    -f "$REF_FASTA" \
    "${BAMS[@]}" | \
bcftools call \
    -mv \
    --ploidy 1 \
    -Ob \
    -o "$OUT_DIR/variants.bcf"

bcftools view \
    "$OUT_DIR/variants.bcf" \
    > "$OUT_DIR/variants.vcf"


# ---------------------------
# 6. Per-sample variant summary
# ---------------------------
OUT_SUM="$OUT_DIR/variant_summary.txt"

echo -e "Sample\tTotal_variants\tSNPs\tIndels\tMultiallelic_sites" \
    > "$OUT_SUM"

for SAMPLE in $(bcftools query -l "$OUT_DIR/variants.vcf")
do
    echo "Processing summary for $SAMPLE..."

    # Keep only non-reference calls
    # Excludes:
    #   0     (haploid reference)
    #   0/0   (diploid reference)
    #   0|0   (phased reference)
    #   ./.   (missing)
    bcftools view \
        -s "$SAMPLE" \
        "$OUT_DIR/variants.vcf" | \
    bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' | \
    awk '
        $5 != "0" &&
        $5 != "0/0" &&
        $5 != "0|0" &&
        $5 != "./."
        {print}
    ' > "$OUT_DIR/temp_${SAMPLE}.tsv"

    TOTAL=$(wc -l < "$OUT_DIR/temp_${SAMPLE}.tsv")

    SNPS=$(
        awk '
        length($3)==1 && length($4)==1 {count++}
        END {print count+0}
        ' "$OUT_DIR/temp_${SAMPLE}.tsv"
    )

    INDELS=$(
        awk '
        length($3)!=1 || length($4)!=1 {count++}
        END {print count+0}
        ' "$OUT_DIR/temp_${SAMPLE}.tsv"
    )

    MULTI=$(
        awk '
        {
            n=split($4,a,",")
            if(n>1) count++
        }
        END {print count+0}
        ' "$OUT_DIR/temp_${SAMPLE}.tsv"
    )

    echo -e \
        "${SAMPLE}\t${TOTAL}\t${SNPS}\t${INDELS}\t${MULTI}" \
        >> "$OUT_SUM"

    rm "$OUT_DIR/temp_${SAMPLE}.tsv"
done


# ---------------------------
# 7. Cleanup
# ---------------------------
rm "$OUT_DIR/variants.bcf"

# ---------------------------
# 8. Finished
# ---------------------------
echo
echo "Pipeline complete!"
echo "Joint VCF: $OUT_DIR/variants.vcf"
echo "Per-sample summary: $OUT_SUM"
echo

date
