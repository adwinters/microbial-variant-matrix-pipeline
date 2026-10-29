#!/bin/bash

hostname
date

# =========================================================
# Pipeline: Single reference, multiple genomes,
# joint variant calling and per-sample variant summaries
#
# Works for any haploid (or, with --ploidy 2, diploid) microbial species
# that has a reference genome.
#
# Inputs can be:
#   - paired-end isolate reads (auto-detected: _R1_001/_R2_001, _R1/_R2 or
#     _1/_2, with .fastq.gz or .fq.gz, optionally .trim.)
#   - assemblies, converted to error-free paired pseudo-reads
#   - or both, jointly called together
# =========================================================

set -euo pipefail

# ---------------------------
# Check arguments
# ---------------------------
usage() {
    echo "Usage: $0 -r <reference.fna> [-f <fastq_dir>] [-a <assembly_dir>] [-o <output_dir>] [-t <threads>]"
    echo "       $0 -S -o <output_dir>        (regenerate the summary from an existing run)"
    echo
    echo "  -r, --reference   reference FASTA (required unless -S)"
    echo "  -f, --fastq       directory of paired-end reads (gzipped). Recognized automatically:"
    echo "                    *_R1_001.trim.fastq.gz, *_R1_001.fastq.gz, *_R1.trim.fastq.gz,"
    echo "                    *_R1.fastq.gz, *_1.fastq.gz (and the .fq.gz / R2 / _2 equivalents)"
    echo "  -a, --assemblies  directory of assemblies (.fasta .fa .fna .fas, optionally .gz)"
    echo "                    to convert to pseudo-reads"
    echo "  -o, --output      output directory (default: results)"
    echo "  -t, --threads     number of threads (default: \$SLURM_CPUS_PER_TASK if set, else 16)"
    echo "  -c, --min-cov-depth minimum depth for a base to count as covered in the coverage"
    echo "                    files written to <output>/coverage (default: 5)"
    echo "      --r1-suffix S / --r2-suffix S"
    echo "                    use these exact read-file suffixes instead of auto-detection"
    echo "                    (e.g. --r1-suffix _R1.clean.fq.gz --r2-suffix _R2.clean.fq.gz)"
    echo "      --ploidy N    ploidy for variant calling (default: 1, haploid; use 2 for diploids)"
    echo "      --modules \"LIST\"  space-separated environment modules to load (default:"
    echo "                    \$PIPELINE_MODULES if set, else \"bwa samtools bcftools\")"
    echo "  -S, --summary-only  skip mapping and calling; rebuild <output>/variant_summary.txt"
    echo "                    from the existing <output>/variants.vcf (needs only -o)"
    echo "  -M, --no-modules  skip 'module load' (use tools already on your PATH)"
    echo "      --keep-pseudo-reads    keep the pseudo-read FASTQs in <output>/pseudo_reads"
    echo "      --delete-pseudo-reads  delete them after the run (default)"
    echo "  -h, --help        show this message"
    echo
    echo "At least one of -f or -a is required (except with -S); use both to call them jointly."
}

REF_FASTA=""
FASTQ_DIR=""
ASSEMBLY_DIR=""
OUT_DIR="results"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
USE_MODULES=1
SUMMARY_ONLY=0
MIN_DEPTH=5
PLOIDY=1
R1_SUFFIX=""
R2_SUFFIX=""
MODULES="${PIPELINE_MODULES:-bwa samtools bcftools}"
KEEP_PSEUDO=0        # set by --keep-pseudo-reads / --delete-pseudo-reads

# Error out if an option is missing its value (args: remaining-arg-count, flag)
need_value() {
    if [[ $1 -lt 2 ]]; then
        echo "ERROR: $2 requires a value."
        exit 1
    fi
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -r|--reference)   need_value $# "$1"; REF_FASTA=$2;    shift 2 ;;
        --reference=*)    REF_FASTA=${1#*=};                   shift ;;
        -f|--fastq)       need_value $# "$1"; FASTQ_DIR=$2;    shift 2 ;;
        --fastq=*)        FASTQ_DIR=${1#*=};                   shift ;;
        -a|--assemblies)  need_value $# "$1"; ASSEMBLY_DIR=$2; shift 2 ;;
        --assemblies=*)   ASSEMBLY_DIR=${1#*=};                shift ;;
        -o|--output)      need_value $# "$1"; OUT_DIR=$2;     shift 2 ;;
        --output=*)       OUT_DIR=${1#*=};                     shift ;;
        -t|--threads)     need_value $# "$1"; THREADS=$2;      shift 2 ;;
        --threads=*)      THREADS=${1#*=};                     shift ;;
        -c|--min-cov-depth) need_value $# "$1"; MIN_DEPTH=$2;    shift 2 ;;
        --min-cov-depth=*) MIN_DEPTH=${1#*=};                   shift ;;
        -S|--summary-only) SUMMARY_ONLY=1;                     shift ;;
        --ploidy)         need_value $# "$1"; PLOIDY=$2;       shift 2 ;;
        --ploidy=*)       PLOIDY=${1#*=};                      shift ;;
        --r1-suffix)      need_value $# "$1"; R1_SUFFIX=$2;    shift 2 ;;
        --r1-suffix=*)    R1_SUFFIX=${1#*=};                   shift ;;
        --r2-suffix)      need_value $# "$1"; R2_SUFFIX=$2;    shift 2 ;;
        --r2-suffix=*)    R2_SUFFIX=${1#*=};                   shift ;;
        --modules)        need_value $# "$1"; MODULES=$2;      shift 2 ;;
        --modules=*)      MODULES=${1#*=};                     shift ;;
        -M|--no-modules)  USE_MODULES=0;                       shift ;;
        --keep-pseudo-reads)   KEEP_PSEUDO=1;                  shift ;;
        --delete-pseudo-reads) KEEP_PSEUDO=0;                  shift ;;
        -h|--help)        usage; exit 0 ;;
        *)
            echo "ERROR: unrecognized argument: $1"
            echo
            usage
            exit 1
            ;;
    esac
done

if [[ -z "$REF_FASTA" && "$SUMMARY_ONLY" != "1" ]]; then
    echo "ERROR: reference FASTA (-r) is required."
    echo
    usage
    exit 1
fi

if ! [[ "$THREADS" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: threads must be a positive integer (got '$THREADS')."
    exit 1
fi

if ! [[ "$MIN_DEPTH" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: min depth must be a positive integer (got '$MIN_DEPTH')."
    exit 1
fi

if ! [[ "$PLOIDY" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: ploidy must be a positive integer (got '$PLOIDY')."
    exit 1
fi

if [[ -n "$R1_SUFFIX" || -n "$R2_SUFFIX" ]]; then
    if [[ -z "$R1_SUFFIX" || -z "$R2_SUFFIX" ]]; then
        echo "ERROR: --r1-suffix and --r2-suffix must be given together."
        exit 1
    fi
fi

# ---------------------------
# Per-sample variant summary (from $OUT_DIR/variants.vcf)
# ---------------------------
summarize_variants() {
    OUT_SUM="$OUT_DIR/variant_summary.txt"

    echo -e "Sample\tTotal_variants\tSNPs\tIndels\tMultiallelic_sites" \
        > "$OUT_SUM"

    for SAMPLE in $(bcftools query -l "$OUT_DIR/variants.vcf")
    do
        echo "Processing summary for $SAMPLE..."

        # Count only calls that carry a non-reference allele.
        # Skipped: 0 (haploid ref), 0/0, 0|0 (diploid/phased ref) and any
        # missing call (., ./., .|.), e.g. a region with no coverage.
        # SNP = REF and every ALT allele are single bases; anything else is an
        # indel. Multiallelic = more than one ALT allele at the site.
        read -r TOTAL SNPS INDELS MULTI < <(
            bcftools view \
                -s "$SAMPLE" \
                "$OUT_DIR/variants.vcf" | \
            bcftools query \
                -f '%REF\t%ALT\t[%GT]\n' | \
            awk -F'\t' '
                {
                    n = split($3, g, /[\/|]/)
                    carried = 0
                    for (i = 1; i <= n; i++) if (g[i] ~ /^[1-9]/) carried = 1
                    if (!carried) next

                    total++
                    m = split($2, a, ",")
                    if (m > 1) multi++

                    snp = (length($1) == 1)
                    for (i = 1; i <= m; i++) if (length(a[i]) != 1) snp = 0
                    if (snp) snps++; else indels++
                }
                END { print total + 0, snps + 0, indels + 0, multi + 0 }
            '
        )

        echo -e \
            "${SAMPLE}\t${TOTAL}\t${SNPS}\t${INDELS}\t${MULTI}" \
            >> "$OUT_SUM"
    done
}

# ---------------------------
# 0. Load required modules
# ---------------------------
if [[ "$USE_MODULES" == "1" ]]; then
    if type module >/dev/null 2>&1; then
        # A broken or missing module system (e.g. on some nodes) must not abort
        # here; the tool check below reports anything that did not load.
        # Cluster-specific names/versions can be set with --modules or the
        # PIPELINE_MODULES environment variable.
        for MOD in $MODULES; do
            module load "$MOD" || echo "WARNING: could not load module '$MOD'"
        done
    else
        echo "NOTE: 'module' command not available; using tools already on PATH."
    fi
fi

REQUIRED_TOOLS=(bwa samtools bcftools)
[[ "$SUMMARY_ONLY" == "1" ]] && REQUIRED_TOOLS=(bcftools)

MISSING=()
for TOOL in "${REQUIRED_TOOLS[@]}"; do
    command -v "$TOOL" >/dev/null 2>&1 || MISSING+=("$TOOL")
done
if [[ ${#MISSING[@]} -gt 0 ]]; then
    echo "ERROR: required tool(s) not found on PATH: ${MISSING[*]}"
    echo "  - Install them, or load them through your module system"
    echo "    (see --modules / PIPELINE_MODULES for cluster-specific names), then rerun."
    echo "  - If they are already on your PATH under another setup, use --no-modules."
    exit 1
fi

# Summary-only mode: rebuild the summary from an existing VCF and stop.
if [[ "$SUMMARY_ONLY" == "1" ]]; then
    if [[ ! -f "$OUT_DIR/variants.vcf" ]]; then
        echo "ERROR: $OUT_DIR/variants.vcf not found (use -o to point at an existing run)."
        exit 1
    fi
    echo "Summary-only mode: rebuilding summary from $OUT_DIR/variants.vcf"
    summarize_variants
    echo
    echo "Per-sample summary: $OUT_SUM"
    echo
    date
    exit 0
fi
# Assembly mode needs python3. If it is not on your PATH, load it here, e.g.:
# module load python3

# Pseudo-read settings (assembly mode only)
# Approx. depth per base = 2 * PSEUDO_READ_LEN / PSEUDO_STEP  (defaults: ~30x)
PSEUDO_READ_LEN=150
PSEUDO_FRAG_LEN=400
PSEUDO_STEP=10

# ---------------------------
# 1. Set directories
# ---------------------------
is_null() {
    local v="${1:-}"
    [[ -z "$v" || "${v^^}" == "NULL" || "${v^^}" == "NONE" ]]
}


if is_null "$FASTQ_DIR" && is_null "$ASSEMBLY_DIR"; then
    echo "ERROR: need a FASTQ directory, an assembly directory, or both."
    exit 1
fi

if ! is_null "$FASTQ_DIR" && [[ ! -d "$FASTQ_DIR" ]]; then
    echo "ERROR: FASTQ directory not found: $FASTQ_DIR"
    exit 1
fi

if ! is_null "$ASSEMBLY_DIR" && [[ ! -d "$ASSEMBLY_DIR" ]]; then
    echo "ERROR: Assembly directory not found: $ASSEMBLY_DIR"
    exit 1
fi

mkdir -p "$OUT_DIR"
PSEUDO_DIR="$OUT_DIR/pseudo_reads"
COVERAGE_DIR="$OUT_DIR/coverage"
mkdir -p "$COVERAGE_DIR"

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
# 4. Helper functions
# ---------------------------

# Tile every assembly in $1 into paired pseudo-reads written to $2.
# Names match the real-read pattern so the mapping loop treats them alike.
make_pseudo_reads() {
python3 - "$1" "$2" "$PSEUDO_READ_LEN" "$PSEUDO_FRAG_LEN" "$PSEUDO_STEP" <<'PY'
import gzip, os, sys

asm_dir, out_dir = sys.argv[1], sys.argv[2]
rl, frag, step = int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
EXTS = (".fasta", ".fa", ".fna", ".fas")
COMP = str.maketrans("ACGTN", "TGCAN")
QUAL = "I" * rl

if frag < rl:
    sys.exit("ERROR: PSEUDO_FRAG_LEN must be >= PSEUDO_READ_LEN")

def sample_name(f):
    n = f[:-3] if f.endswith(".gz") else f
    for e in EXTS:
        if n.lower().endswith(e):
            return n[:-len(e)]
    return None

def read_fasta(path):
    op = gzip.open(path, "rt") if path.endswith(".gz") else open(path)
    name, parts = None, []
    with op as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts).upper()
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip())
    if name is not None:
        yield name, "".join(parts).upper()

files = sorted(f for f in os.listdir(asm_dir) if sample_name(f))
if not files:
    sys.exit("ERROR: no FASTA files found in " + asm_dir)

os.makedirs(out_dir, exist_ok=True)
seen = set()
for f in files:
    s = sample_name(f)
    if s in seen:
        sys.exit("ERROR: duplicate assembly sample name: " + s)
    seen.add(s)
    pairs = skipped = 0
    pre = os.path.join(out_dir, s)
    with gzip.open(pre + "_R1_001.trim.fastq.gz", "wt") as r1, \
         gzip.open(pre + "_R2_001.trim.fastq.gz", "wt") as r2:
        for name, seq in read_fasta(os.path.join(asm_dir, f)):
            L = len(seq)
            if L < frag:
                skipped += 1
                continue
            starts = list(range(0, L - frag + 1, step))
            if starts[-1] != L - frag:
                starts.append(L - frag)
            for p in starts:
                fr = seq[p:p + frag]
                r1.write("@%s_%d/1\n%s\n+\n%s\n" % (name, p, fr[:rl], QUAL))
                r2.write("@%s_%d/2\n%s\n+\n%s\n" % (name, p, fr[-rl:].translate(COMP)[::-1], QUAL))
                pairs += 1
    msg = "  %s: %d read pairs" % (s, pairs)
    if skipped:
        msg += " (%d contigs < %d bp skipped)" % (skipped, frag)
    print(msg)
PY
}

BAMS=()
SAMPLES_SEEN=()

# Map one sample (R1, R2, name) and append its sorted BAM to BAMS
map_sample() {
    local FQ1=$1 FQ2=$2 SAMPLE=$3

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

    # Covered regions: merged intervals (BED, 0-based) where depth >= MIN_DEPTH.
    # Used downstream to tell "no variants" apart from "no coverage / absent".
    samtools depth -a "$OUT_DIR/${SAMPLE}.sorted.bam" | \
    awk -v min="$MIN_DEPTH" 'BEGIN { OFS = "\t" }
        $3 >= min {
            if (open && $1 == chr && $2 == end + 1) {
                end = $2
            } else {
                if (open) print chr, start - 1, end
                chr = $1; start = $2; end = $2; open = 1
            }
            next
        }
        { if (open) { print chr, start - 1, end; open = 0 } }
        END { if (open) print chr, start - 1, end }
    ' > "$COVERAGE_DIR/${SAMPLE}.covered.bed"

    BAMS+=("$OUT_DIR/${SAMPLE}.sorted.bam")
    SAMPLES_SEEN+=("$SAMPLE")

    rm \
        "$OUT_DIR/${SAMPLE}.bam" \
        "$OUT_DIR/${SAMPLE}.unsorted.bam"
}

# ---------------------------
# 5. Process real isolate FASTQ files
# ---------------------------
if ! is_null "$FASTQ_DIR"; then
    echo "Mapping isolate reads from: $FASTQ_DIR"

    # "R1 suffix:R2 suffix" pairs to look for (first match per sample wins)
    if [[ -n "$R1_SUFFIX" ]]; then
        READ_PATTERNS=("${R1_SUFFIX}:${R2_SUFFIX}")
    else
        READ_PATTERNS=(
            "_R1_001.trim.fastq.gz:_R2_001.trim.fastq.gz"
            "_R1_001.fastq.gz:_R2_001.fastq.gz"
            "_R1.trim.fastq.gz:_R2.trim.fastq.gz"
            "_R1.fastq.gz:_R2.fastq.gz"
            "_1.fastq.gz:_2.fastq.gz"
            "_R1_001.trim.fq.gz:_R2_001.trim.fq.gz"
            "_R1.fq.gz:_R2.fq.gz"
            "_1.fq.gz:_2.fq.gz"
        )
    fi

    for PATTERN in "${READ_PATTERNS[@]}"; do
        R1S=${PATTERN%%:*}
        R2S=${PATTERN#*:}

        for FQ1 in "$FASTQ_DIR"/*"$R1S"
        do
            [[ -e "$FQ1" ]] || continue

            SAMPLE=$(basename "$FQ1" "$R1S")
            FQ2="$FASTQ_DIR/${SAMPLE}${R2S}"

            if [[ ! -f "$FQ2" ]]; then
                echo "ERROR: Paired file for $SAMPLE not found (expected $FQ2)."
                exit 1
            fi

            # Skip a sample already mapped via a different file pattern
            # (e.g. both trimmed and untrimmed reads in the same folder).
            ALREADY=0
            for S in "${SAMPLES_SEEN[@]+"${SAMPLES_SEEN[@]}"}"; do
                [[ "$S" == "$SAMPLE" ]] && ALREADY=1
            done
            if [[ "$ALREADY" == "1" ]]; then
                echo "WARNING: $SAMPLE already mapped from another file pattern; skipping $FQ1"
                continue
            fi

            map_sample "$FQ1" "$FQ2" "$SAMPLE"
        done
    done
fi

# ---------------------------
# 6. Convert assemblies to pseudo-reads and map them
# ---------------------------
if ! is_null "$ASSEMBLY_DIR"; then
    if ! command -v python3 >/dev/null 2>&1; then
        echo "ERROR: python3 is required for assembly mode but was not found."
        exit 1
    fi

    echo "Generating pseudo-reads from assemblies in: $ASSEMBLY_DIR"
    echo "  read length ${PSEUDO_READ_LEN}, fragment ${PSEUDO_FRAG_LEN}, step ${PSEUDO_STEP}"
    make_pseudo_reads "$ASSEMBLY_DIR" "$PSEUDO_DIR"

    for FQ1 in "$PSEUDO_DIR"/*_R1_001.trim.fastq.gz
    do
        [[ -e "$FQ1" ]] || continue

        SAMPLE=$(basename "$FQ1" _R1_001.trim.fastq.gz)
        FQ2="$PSEUDO_DIR/${SAMPLE}_R2_001.trim.fastq.gz"

        # An assembly and an isolate read set must not share a sample name
        for S in "${SAMPLES_SEEN[@]+"${SAMPLES_SEEN[@]}"}"; do
            if [[ "$S" == "$SAMPLE" ]]; then
                echo "ERROR: assembly sample name '$SAMPLE' collides with an isolate read sample."
                exit 1
            fi
        done

        map_sample "$FQ1" "$FQ2" "$SAMPLE"
    done
fi

if [[ ${#BAMS[@]} -eq 0 ]]; then
    echo "ERROR: no samples were mapped. Check the read file names in the FASTQ"
    echo "       directory (see --help for recognized patterns, or use --r1-suffix/--r2-suffix)."
    exit 1
fi

# ---------------------------
# 7. Joint variant calling
# ---------------------------
echo "Calling variants jointly across ${#BAMS[@]} samples..."

bcftools mpileup \
    --threads "$THREADS" \
    -Ou \
    -f "$REF_FASTA" \
    "${BAMS[@]}" | \
bcftools call \
    -mv \
    --ploidy "$PLOIDY" \
    -Ob \
    -o "$OUT_DIR/variants.bcf"

bcftools view \
    "$OUT_DIR/variants.bcf" \
    > "$OUT_DIR/variants.vcf"


# ---------------------------
# 8. Per-sample variant summary
# ---------------------------
summarize_variants


# ---------------------------
# 9. Cleanup
# ---------------------------
rm "$OUT_DIR/variants.bcf"

if [[ -d "$PSEUDO_DIR" && "$KEEP_PSEUDO" != "1" ]]; then
    rm -r "$PSEUDO_DIR"
fi

# ---------------------------
# 10. Finished
# ---------------------------
echo
echo "Pipeline complete!"
echo "Joint VCF: $OUT_DIR/variants.vcf"
echo "Per-sample summary: $OUT_SUM"
echo "Covered regions (depth >= $MIN_DEPTH): $COVERAGE_DIR/"
echo

date
