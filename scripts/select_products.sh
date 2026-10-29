#!/bin/bash

usage() {
    echo "Usage: $0 -m <matrix.csv> -o <selected_genes.tsv> (-l <patterns.tsv> | -p <pattern>) [...]"
    echo
    echo "Select products from a variant matrix by searching their names."
    echo "Matching is case-insensitive and uses regular expressions, so a plain"
    echo "word such as 'efflux' matches any product name containing it."
    echo "Use ^ and \$ to require an exact name, e.g. '^DNA primase\$'."
    echo
    echo "  -m, --matrix    variants_per_product_matrix(.filtered).csv to search (required)"
    echo "  -o, --output    output file path, e.g. results/selected_genes.tsv (required)"
    echo "  -l, --list      file of patterns, one per line; only the first tab-separated"
    echo "                  column is used (the rest can be notes); blank lines and lines"
    echo "                  starting with # are ignored. May be given more than once."
    echo "  -p, --pattern   a single pattern. May be given more than once."
    echo "  -h, --help      show this message"
    echo
    echo "At least one -l or -p is required. A product matching several patterns is"
    echo "written once. Patterns that match nothing are reported."
}

MATRIX=""
OUTFILE=""
LISTS=()
PATTERNS=()

while [ "$#" -gt 0 ]; do
    case "$1" in
        -m|--matrix)  MATRIX="$2"; shift 2 ;;
        -o|--output)  OUTFILE="$2"; shift 2 ;;
        -l|--list)    LISTS+=("$2"); shift 2 ;;
        -p|--pattern) PATTERNS+=("$2"); shift 2 ;;
        --matrix=*)   MATRIX="${1#*=}"; shift ;;
        --output=*)   OUTFILE="${1#*=}"; shift ;;
        --list=*)     LISTS+=("${1#*=}"); shift ;;
        --pattern=*)  PATTERNS+=("${1#*=}"); shift ;;
        -h|--help)    usage; exit 0 ;;
        *) echo "ERROR: unknown argument: $1"; usage; exit 1 ;;
    esac
done

if [ -z "$MATRIX" ] || [ -z "$OUTFILE" ]; then
    echo "ERROR: both -m and -o are required."
    usage
    exit 1
fi

if [ "${#LISTS[@]}" -eq 0 ] && [ "${#PATTERNS[@]}" -eq 0 ]; then
    echo "ERROR: give at least one pattern with -l <file> or -p <pattern>."
    usage
    exit 1
fi

if [ ! -f "$MATRIX" ]; then
    echo "ERROR: matrix file not found: $MATRIX"
    exit 1
fi

for f in "${LISTS[@]}"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: pattern file not found: $f"
        exit 1
    fi
done

export MATRIX OUTFILE

# Pass lists and patterns to Python as arguments (list files first, then -p).
python3 - "${#LISTS[@]}" "${LISTS[@]}" "${PATTERNS[@]}" <<'EOF'
import os
import re
import sys
import pandas as pd

matrix_path = os.environ["MATRIX"]
outfile = os.environ["OUTFILE"]

n_lists = int(sys.argv[1])
list_files = sys.argv[2:2 + n_lists]
inline = sys.argv[2 + n_lists:]

patterns = []
for lf in list_files:
    with open(lf) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            pat = line.split("\t")[0].strip()
            if pat:
                patterns.append(pat)
patterns.extend(p for p in inline if p.strip())

compiled = []
for p in patterns:
    try:
        compiled.append((p, re.compile(p, re.IGNORECASE)))
    except re.error as e:
        sys.exit(f"ERROR: invalid pattern '{p}': {e}")

print(f"Reading: {matrix_path}")
products = list(pd.read_csv(matrix_path, index_col=0, nrows=0).columns)
print(f"Products in matrix: {len(products)}")
print(f"Patterns: {len(compiled)}")

selected = []
seen = set()
unmatched = []
broad = []

for p, rx in compiled:
    hits = [c for c in products if rx.search(c)]
    if not hits:
        unmatched.append(p)
        continue
    if len(hits) > 1:
        broad.append((p, hits))
    for h in hits:
        if h not in seen:
            seen.add(h)
            selected.append(h)

# keep matrix column order so output is stable
order = {c: i for i, c in enumerate(products)}
selected.sort(key=order.get)

os.makedirs(os.path.dirname(os.path.abspath(outfile)), exist_ok=True)
pd.Series(selected).to_csv(outfile, sep="\t", index=False, header=False)

if unmatched:
    print("\nPatterns that matched no product (check spelling, or the product "
          "may have been removed by the coverage filter):")
    for p in unmatched:
        print(f"  {p}")

if broad:
    print("\nPatterns that matched more than one product:")
    for p, hits in broad:
        shown = "; ".join(hits[:5]) + (f"; ... (+{len(hits) - 5} more)" if len(hits) > 5 else "")
        print(f"  {p}  ->  {len(hits)} products: {shown}")

print(f"\nTotal selected products: {len(selected)}")
print(f"Output file: {outfile}")
EOF
