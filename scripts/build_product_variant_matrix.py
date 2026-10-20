#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import math
import pandas as pd

from concurrent.futures import ProcessPoolExecutor, as_completed

# ------------------------------------------------------------
# Step 1: Parse command-line arguments
# ------------------------------------------------------------
if len(sys.argv) < 3:
    sys.exit(
        "Usage: python3 build_product_variant_matrix.py "
        "<variants.vcf> <genomic.gff> [output_dir]"
    )

vcf_file = sys.argv[1]
gff_file = sys.argv[2]
results_dir = sys.argv[3] if len(sys.argv) > 3 else "results"

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

# Build locus_tag -> product lookup
for _, row in gff_df.iterrows():

    attrs = {
        kv.split("=")[0]: kv.split("=")[1]
        for kv in row["attributes"].split(";")
        if "=" in kv
    }

    locus_tag = attrs.get("locus_tag")
    product = attrs.get("product")

    if locus_tag and product:
        locus_to_product[locus_tag] = product

# Parse gene entries
gff_genes = gff_df[gff_df["type"] == "gene"]

total_genes = len(gff_genes)
missing_product_count = 0

for count, (_, row) in enumerate(gff_genes.iterrows(), start=1):

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

if missing_product_count > 0:
    print(
        f"\nWarning: {missing_product_count} genes had no "
        f"product annotation. Using locus_tag instead."
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

            gt = sample_dict.get("GT", "./.")

            if gt in ["0", "0/0", "0|0", "./."]:
                continue

            matched = chrom_genes[
                (chrom_genes["Start"] <= pos) &
                (chrom_genes["End"] >= pos)
            ]

            if matched.empty:
                continue

            product_name = matched.iloc[0]["Product"]

            for allele in alt_alleles:

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

print("\nPipeline complete.")
