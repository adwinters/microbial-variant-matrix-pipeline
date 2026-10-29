#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import math
import pandas as pd

from urllib.parse import unquote
from concurrent.futures import ProcessPoolExecutor, as_completed

# ------------------------------------------------------------
# Step 1: Parse command-line arguments
# ------------------------------------------------------------
import argparse

ap = argparse.ArgumentParser(
    description="Build per-product variant and coverage matrices.",
    epilog="Legacy positional form is still accepted: "
           "<vcf> <gff> [output_dir] [coverage_dir]"
)
ap.add_argument("-v", "--vcf", help="joint VCF from call_variants.sh (required)")
ap.add_argument("-g", "--gff", help="GFF annotation for the same reference (required)")
ap.add_argument("-o", "--output", dest="output_dir",
                help="output directory (default: results)")
ap.add_argument(
    "--coverage-dir", dest="coverage_dir",
    help="directory of <sample>.covered.bed files "
         "(default: <output>/coverage; skipped if absent)"
)
ap.add_argument(
    "--min-coverage", type=float, default=0.8,
    help="A sample/product cell counts as covered when at least this "
         "fraction of the product's gene bases are covered (default 0.8)"
)
ap.add_argument(
    "--min-samples-covered", type=float, default=1.0,
    help="Keep a product in the filtered matrix only if it is covered in "
         "at least this fraction of samples (default 1.0; 0 keeps all)"
)
ap.add_argument(
    "--snps-only", action="store_true",
    help="Count only single-base substitutions (ignore indels)"
)
ap.add_argument("positional", nargs="*", help=argparse.SUPPRESS)
args = ap.parse_args()

if len(args.positional) > 4:
    ap.error("too many positional arguments; use -v, -g, -o and --coverage-dir")

# Legacy positional arguments fill in whatever the flags did not set.
_legacy = list(args.positional)
for _name in ("vcf", "gff", "output_dir", "coverage_dir"):
    if getattr(args, _name) is None and _legacy:
        setattr(args, _name, _legacy.pop(0))

if not args.vcf or not args.gff:
    ap.error("both -v/--vcf and -g/--gff are required")

vcf_file = args.vcf
gff_file = args.gff
results_dir = args.output_dir or "results"
MIN_COVERAGE = args.min_coverage
MIN_SAMPLES_COVERED = args.min_samples_covered
SNPS_ONLY = args.snps_only

# Directory of <sample>.covered.bed files written by call_variants.sh.
coverage_dir = args.coverage_dir or os.path.join(results_dir, "coverage")

if not os.path.exists(vcf_file):
    sys.exit(f"Error: {vcf_file} not found.")

if not os.path.exists(gff_file):
    sys.exit(f"Error: {gff_file} not found.")

os.makedirs(results_dir, exist_ok=True)

variants_table_file = os.path.join(
    results_dir,
    "variants_table.tsv"
)

matrix_file = os.path.join(
    results_dir,
    "variants_per_product_matrix.csv"
)

coverage_file = os.path.join(
    results_dir,
    "product_coverage_matrix.csv"
)

filtered_matrix_file = os.path.join(
    results_dir,
    "variants_per_product_matrix.filtered.csv"
)

filtered_genes_file = os.path.join(
    results_dir,
    "filtered_genes.tsv"
)

no_variation_file = os.path.join(
    results_dir,
    "products_no_variation.tsv"
)

print("Starting variant analysis pipeline...\n")

# ---------------------------
# Step 1: Parse GFF file
# ---------------------------
print("Parsing GFF file to extract gene positions...")

gene_positions = []
locus_to_product = {}

gff_df = pd.read_csv(
    gff_file,
    sep="\t",
    comment="#",
    header=None,
    names=[
        "seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes"
    ]
)

print(f"Total GFF records: {len(gff_df)}")

if gff_df.empty:
    sys.exit(
        f"Error: no feature lines were read from {gff_file} "
        f"(file is empty or contains only '#' lines). "
        f"Check the path and the file size with: ls -l {gff_file}"
    )

print("\nFeature types:")
print(gff_df["type"].value_counts())

# Build locus_tag -> product lookup
for _, row in gff_df.iterrows():

    attrs = {
        kv.split("=")[0]: kv.split("=")[1]
        for kv in row["attributes"].split(";")
        if "=" in kv
    }

    locus_tag = attrs.get("locus_tag")
    product = attrs.get("product")

    if product:
        product = unquote(product)

    if locus_tag and product:
        locus_to_product[locus_tag] = product

# Parse gene entries
gff_genes = gff_df[gff_df["type"] == "gene"]

print(f"\nGene records: {len(gff_genes)}")

total_genes = len(gff_genes)
missing_product_count = 0

for count, (_, row) in enumerate(
    gff_genes.iterrows(),
    start=1
):

    attrs = {
        kv.split("=")[0]: kv.split("=")[1]
        for kv in row["attributes"].split(";")
        if "=" in kv
    }

    locus_tag = attrs.get("locus_tag")
    product = locus_to_product.get(locus_tag)

    if product is None:
        missing_product_count += 1
        product = locus_tag

    gene_positions.append(
        {
            "Chromosome": row["seqid"],
            "Start": int(row["start"]),
            "End": int(row["end"]),
            "Product": product
        }
    )

    if count % 100 == 0 or count == total_genes:
        print(
            f"  Parsed genes: {count}/{total_genes}",
            end="\r"
        )

print(f"\nGene positions loaded: {len(gene_positions)}")

if missing_product_count > 0:
    print(
        f"\nWarning: {missing_product_count} genes had no "
        f"product annotation. Using locus_tag instead."
    )

if not gene_positions:
    sys.exit(
        f"Error: no 'gene' features found in {gff_file}. "
        f"This script expects NCBI-style GFF3 with 'gene' rows."
    )

gene_df = pd.DataFrame(gene_positions)

print(
    f"\nFinished parsing GFF: "
    f"{len(gene_df)} genes loaded.\n"
)

# ---------------------------
# Step 2: Parse VCF chunks
# ---------------------------
def parse_vcf_chunk(lines_chunk, gene_df, sample_cols):

    vcf_rows_chunk = []

    for line in lines_chunk:

        fields = line.strip().split("\t")

        if len(fields) < 10:
            continue

        chrom, pos, _, ref, alt, _, _, _, fmt = fields[:9]
        samples = fields[9:]

        pos = int(pos)
        alt_alleles = alt.split(",")

        fmt_keys = fmt.split(":")

        chrom_genes = gene_df[
            gene_df["Chromosome"] == chrom
        ]

        for sample_name, sample_val in zip(sample_cols, samples):

            sample_fields = sample_val.split(":")
            sample_dict = dict(zip(fmt_keys, sample_fields))

            gt = sample_dict.get("GT", ".")

            # Skip reference and missing calls. Haploid calls with no coverage
            # (e.g. a region absent from an assembly) are written as a bare ".",
            # so test each allele instead of matching a fixed list of strings.
            gt_alleles = gt.replace("|", "/").split("/")

            # ALT allele(s) this sample actually carries (1-based indices into
            # ALT). At a multiallelic site, only the carried allele is counted
            # for this sample, not every ALT listed at the site.
            carried = sorted({
                int(a) for a in gt_alleles
                if a.isdigit() and a != "0"
            })

            if not carried:
                continue

            matched = chrom_genes[
                (chrom_genes["Start"] <= pos) &
                (chrom_genes["End"] >= pos)
            ]

            if matched.empty:
                continue

            product_name = matched.iloc[0]["Product"]

            for idx in carried:

                if idx > len(alt_alleles):
                    continue

                allele = alt_alleles[idx - 1]

                if SNPS_ONLY and (len(ref) != 1 or len(allele) != 1):
                    continue

                variant_id = (
                    f"{chrom}_{pos}_{allele}"
                )

                vcf_rows_chunk.append(
                    {
                        "Sample": sample_name,
                        "Product": product_name,
                        "VariantID": variant_id
                    }
                )

    return vcf_rows_chunk

# ---------------------------
# Step 3: Read VCF
# ---------------------------
print(
    "Reading VCF file and splitting into chunks "
    "for multiprocessing..."
)

with open(vcf_file, "r") as f:

    lines = [
        line
        for line in f
        if not line.startswith("##")
    ]

header_line = [
    line
    for line in lines
    if line.startswith("#CHROM")
][0]

sample_cols = header_line.strip().split("\t")[9:]

data_lines = [
    line
    for line in lines
    if not line.startswith("#CHROM")
]

total_variants = len(data_lines)

num_cores = os.cpu_count() or 1

base_chunk_size = max(
    10000,
    min(
        50000,
        max(1, total_variants // (num_cores * 2))
    )
)

chunk_size = base_chunk_size

chunks = [
    data_lines[i:i + chunk_size]
    for i in range(0, total_variants, chunk_size)
]

total_chunks = len(chunks)

print(f"Detected {num_cores} CPU cores")
print(f"Total variant lines: {total_variants}")
print(f"Number of chunks: {total_chunks}")
print(f"Chunk size: {chunk_size}\n")

# ---------------------------
# Step 4: Multiprocessing
# ---------------------------
print("Parsing VCF chunks...\n")

vcf_rows = []

with ProcessPoolExecutor() as executor:

    futures = {
        executor.submit(
            parse_vcf_chunk,
            chunk,
            gene_df,
            sample_cols
        ): i + 1
        for i, chunk in enumerate(chunks)
    }

    completed = 0

    for future in as_completed(futures):

        chunk_result = future.result()

        vcf_rows.extend(chunk_result)

        completed += 1

        print(
            f"  Completed chunks: "
            f"{completed}/{total_chunks}",
            end="\r"
        )

print(
    f"\nFinished VCF parsing: "
    f"{len(vcf_rows)} variant observations collected.\n"
)

# ---------------------------
# Step 5: variants_table.tsv
# ---------------------------
variants_table = pd.DataFrame(vcf_rows)

if variants_table.empty:

    raise ValueError(
        "No variants were parsed from the VCF."
    )

variants_table.to_csv(
    variants_table_file,
    sep="\t",
    index=False
)

print(
    f"Saved variants table "
    f"({len(variants_table)} rows)"
)

# ---------------------------
# Step 6: Product matrix
# ---------------------------
print(
    "\nCreating per-product variant matrix..."
)

variants_per_product_matrix = (
    variants_table
    .groupby(["Sample", "Product"])
    .size()
    .unstack(fill_value=0)
)

all_products = sorted(
    list(
        set(
            gene_df["Product"]
        )
    )
)

variants_per_product_matrix = (
    variants_per_product_matrix
    .reindex(
        columns=all_products,
        fill_value=0
    )
)

variants_per_product_matrix.to_csv(
    matrix_file,
    index=True
)

print(
    f"Saved matrix: "
    f"{variants_per_product_matrix.shape[0]} samples x "
    f"{variants_per_product_matrix.shape[1]} products"
)

print(
    f"Matrix file size: "
    f"{os.path.getsize(matrix_file):,} bytes"
)

print(f"Matrix saved to: {matrix_file}")

# ---------------------------
# Step 7: Per-product coverage matrix (optional)
#
# For each sample and product: fraction of the product's gene bases that
# are covered (depth >= the threshold used by call_variants.sh). Lets
# downstream plots distinguish "no variants" (matches reference) from
# "no coverage" (gene absent or unassembled).
# ---------------------------
if os.path.isdir(coverage_dir):

    import numpy as np

    print(f"\nBuilding per-product coverage from: {coverage_dir}")

    genes = gene_df.copy()
    genes["Length"] = genes["End"] - genes["Start"] + 1
    total_len = genes.groupby("Product")["Length"].sum()

    coverage_rows = {}

    for sample in sample_cols:

        bed_path = os.path.join(coverage_dir, f"{sample}.covered.bed")

        if not os.path.exists(bed_path):
            print(f"  Warning: no coverage file for {sample}; skipping.")
            continue

        # Merged covered intervals per chromosome (0-based, half-open)
        intervals = {}

        with open(bed_path) as bed:
            for line in bed:
                f = line.split("\t")
                if len(f) < 3:
                    continue
                intervals.setdefault(f[0], []).append(
                    (int(f[1]), int(f[2]))
                )

        covered_bases = pd.Series(0, index=genes.index, dtype="int64")

        for chrom, cg in genes.groupby("Chromosome"):

            chrom_len = int(cg["End"].max())

            if chrom in intervals:
                chrom_len = max(
                    chrom_len,
                    max(e for _, e in intervals[chrom])
                )

            diff = np.zeros(chrom_len + 2, dtype=np.int32)

            for s0, e0 in intervals.get(chrom, []):
                diff[s0] += 1
                diff[e0] -= 1

            is_covered = (np.cumsum(diff) > 0).astype(np.int64)
            prefix = np.concatenate(([0], np.cumsum(is_covered)))

            starts = cg["Start"].to_numpy() - 1
            ends = cg["End"].to_numpy()

            covered_bases.loc[cg.index] = prefix[ends] - prefix[starts]

        per_product = (
            covered_bases.groupby(genes["Product"]).sum() / total_len
        )

        coverage_rows[sample] = per_product

    if coverage_rows:

        coverage_matrix = (
            pd.DataFrame(coverage_rows).T
            .reindex(columns=all_products)
            .fillna(0)
            .round(4)
        )
        coverage_matrix.index.name = "Sample"

        coverage_matrix.to_csv(coverage_file, index=True)

        print(
            f"Saved coverage matrix: "
            f"{coverage_matrix.shape[0]} samples x "
            f"{coverage_matrix.shape[1]} products"
        )
        print(f"Coverage matrix saved to: {coverage_file}")

        # Filtered matrix: cells with < MIN_COVERAGE of the gene covered are
        # blank (unknown, not zero); products covered in fewer than
        # MIN_SAMPLES_COVERED of samples are dropped.
        covered = coverage_matrix >= MIN_COVERAGE
        covered = covered.reindex(
            index=variants_per_product_matrix.index, fill_value=True
        )
        keep = covered.mean(axis=0) >= MIN_SAMPLES_COVERED - 1e-9
        filtered = (
            variants_per_product_matrix.where(covered)
            .loc[:, keep]
        )
        filtered.to_csv(filtered_matrix_file, index=True)

        # Report of every product removed by the coverage filter.
        n_samples = covered.shape[0]
        dropped_rows = []

        for prod in keep.index[~keep.to_numpy()]:
            n_cov = int(covered[prod].sum())
            failed = [
                smp for smp in covered.index if not covered.at[smp, prod]
            ]
            dropped_rows.append({
                "Product": prod,
                "Samples_covered": f"{n_cov}/{n_samples}",
                "Reason": (
                    f"covered in {n_cov} of {n_samples} samples "
                    f"({n_cov / n_samples:.0%}), below "
                    f"--min-samples-covered {MIN_SAMPLES_COVERED:.0%}; "
                    f"<{MIN_COVERAGE:.0%} of gene bases covered in: "
                    f"{', '.join(failed)}"
                )
            })

        pd.DataFrame(
            dropped_rows, columns=["Product", "Samples_covered", "Reason"]
        ).sort_values("Product").to_csv(
            filtered_genes_file, sep="\t", index=False
        )

        # Products that survive the coverage filter but have the same value
        # in every covered sample (e.g. no variants anywhere). They carry no
        # variation, so the heatmap script would drop them.
        n_unique = filtered.nunique(axis=0, dropna=True)
        flat_products = n_unique.index[(n_unique <= 1).to_numpy()]
        flat_rows = []

        for prod in flat_products:
            vals = filtered[prod].dropna()
            blank = [
                smp for smp in filtered.index if pd.isna(filtered.at[smp, prod])
            ]

            if vals.empty:
                reason = "no covered samples"
            elif float(vals.iloc[0]) == 0:
                reason = "no variants in any covered sample"
            else:
                n_var = int(vals.iloc[0])
                reason = (
                    f"{n_var} {'variant' if n_var == 1 else 'variants'} "
                    f"in every covered sample"
                )

            if blank:
                reason += f"; not covered (blank) in: {', '.join(blank)}"

            flat_rows.append({
                "Product": prod,
                "Samples_covered": f"{len(vals)}/{filtered.shape[0]}",
                "Reason": reason
            })

        pd.DataFrame(
            flat_rows, columns=["Product", "Samples_covered", "Reason"]
        ).sort_values("Product").to_csv(
            no_variation_file, sep="\t", index=False
        )

        print(
            f"Filtered matrix: kept {int(keep.sum())} of {len(keep)} "
            f"products (>= {MIN_COVERAGE:.0%} of gene covered in "
            f">= {MIN_SAMPLES_COVERED:.0%} of samples); "
            f"{int((~covered).values.sum())} low-coverage cells blanked."
        )
        print(f"Filtered matrix saved to: {filtered_matrix_file}")
        print(
            f"Dropped products ({len(dropped_rows)}) listed in: "
            f"{filtered_genes_file}"
        )
        print(
            f"Kept but identical in every covered sample ({len(flat_rows)}) "
            f"listed in: {no_variation_file}"
        )

else:

    print(
        f"\nNote: coverage directory not found ({coverage_dir}); "
        f"skipping product coverage matrix."
    )

print("\nPipeline complete.")