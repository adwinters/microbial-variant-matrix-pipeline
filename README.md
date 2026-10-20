# microbial-variant-matrix-pipeline

This pipeline transforms sequencing reads into product-level variant matrices for comparative genomic analysis and visualization. Starting from paired-end FASTQ files, the workflow performs read alignment, joint variant calling, variant aggregation by annotated gene product, and heatmap generation.

The pipeline performs:

1. Read alignment against a reference genome
2. Joint variant calling across isolates
3. Aggregation of variants by annotated gene product
4. Generation of sample-by-product variant matrices
5. Identification of highly variable gene products
6. Visualization of variant patterns through clustered heatmaps

The workflow can be applied to microbial isolate whole-genome sequencing datasets provided that:

- A reference FASTA file is available
- A GFF annotation file generated from the same reference assembly is available
- Paired-end genomic DNA FASTQ files are available for each isolate
  

## Software Requirements

### Variant Calling

BWA >= 0.7.17
SAMtools >= 1.21
BCFtools


### Python
Python >= 3.10
pandas
numpy
scipy
matplotlib
seaborn

Install Python dependencies:

pip install pandas numpy scipy matplotlib seaborn

## Pipeline Workflow

```text
Paired-End FASTQ Files
            │
            ▼
call_variants.sh
            │
            ▼
variants.vcf
variant_summary.txt
            │
            ▼
build_product_variant_matrix.py
            │
            ▼
variants_per_product_matrix.csv
variants_table.tsv
            │
            ├──────────────► plot_selected_products_heatmap.py
            │                            │
            │                            ▼
            │              Top-N Variable Products
            │
            ▼
select_products.sh
            │
            ▼
selected_genes.tsv
            │
            ▼
plot_selected_products_heatmap.py
            │
            ▼
Selected Product Heatmap
```

## Methods

Reads are aligned against a reference genome using BWA-MEM. Variants are jointly called using BCFtools and aggregated by annotated gene product using the supplied GFF annotation. Product-level variant counts are transformed into a sample-by-product matrix and standardized using Z-score normalization prior to heatmap generation.


# Step 1: Variant Calling

Align reads to a reference genome and perform joint variant calling.

## Usage

```text
bash call_variants.sh \
    path/to/reference/genome.fna \
    path/to/trimmed_fastq \
    results
```

## Inputs

Reference FASTA
Paired-end FASTQ files


## Outputs
```text
results/
├── variants.vcf
└── variant_summary.txt
```

### Example Variant Summary
```text
Sample     Total_variants    SNPs    Indels
Isolate1   58364             56663   1701
Isolate2   57902             56200   1702
Isolate3   57486             55792   1694
```

# Step 2: Build Product-Level Variant Matrix

Convert variant calls into a sample-by-product matrix using the reference annotation.

## Usage
```text
python3 build_product_variant_matrix.py \
    results/variants.vcf \
    reference/genome.gff \
    results
```

## Inputs
```text
variants.vcf
genome.gff
```

## Outputs
```text
results/
├── variants_table.tsv
└── variants_per_product_matrix.csv
```

The matrix structure is:
```text
Sample      Cas1   Tet(Q)   CfxA
Isolate1      2       5       0
Isolate2      1       8       3
```


# Step 3: Generate Curated Product List

Create a biologically curated list of products of interest.

Current targets include:

- CRISPR-associated proteins
- Antimicrobial resistance proteins
- Efflux transporters
- Virulence-associated proteins
- DNA repair proteins
- Replication proteins
- Transcriptional regulators

## Usage
```text
bash select_products.sh results
```

## Output
```text
results/
└── selected_genes.tsv
```

# Step 4: Generate Heatmaps

The heatmap script supports two modes:

## Option A. Most Variable Products

Plot the most variable products across all isolates.

### Default (Top 20 Products)
```text
python3 plot_selected_products_heatmap.py \
    results/variants_per_product_matrix.csv

### Custom Number of Products

python3 plot_selected_products_heatmap.py \
    results/variants_per_product_matrix.csv \
    50
```

Outputs:
```text
results/
├── variants_per_product_heatmap_top20.png
└── products_plotted_top20.tsv
```

### Example Output

![Top Variable Products](example_output/top.png)

The heatmap displays products ranked by variability across isolates.


## Option B. Curated Product Panel

Plot products listed in:
```text
selected_genes.tsv
```

Invariant products are automatically removed before plotting.

### Usage
```text
python3 plot_selected_products_heatmap.py \
    results/variants_per_product_matrix.csv \
    results/selected_genes.tsv
```

Outputs:
```text
results/
├── variants_per_product_heatmap_selected.png
└── products_plotted_selected.tsv
```

### Example Output

![Top Variable Products](example_output/selected.png)

The heatmap displays a curated panel of biologically relevant gene products selected by the user.


## Pipeline Outputs
```text
variants.vcf
variant_summary.txt
variants_table.tsv
variants_per_product_matrix.csv
selected_genes.tsv

products_plotted_selected.tsv
products_plotted_topN.tsv

variants_per_product_heatmap_selected.png
variants_per_product_heatmap_topN.png
```

