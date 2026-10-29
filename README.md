# microbial-variant-matrix-pipeline

This pipeline transforms sequencing reads and/or genome assemblies into product-level variant matrices for comparative genomic analysis and visualization. Starting from paired-end FASTQ files, assembled genomes, or both, the workflow performs alignment to a reference genome, joint variant calling, variant aggregation by annotated gene product, coverage-aware filtering, and heatmap generation.

The pipeline performs:

1. Alignment of reads (and assemblies, converted to pseudo-reads) against a reference genome
2. Joint variant calling across all isolates and assemblies
3. Per-sample coverage mapping, to distinguish "no variants" from "not covered"
4. Aggregation of variants by annotated gene product
5. Generation of sample-by-product variant matrices, with a coverage-filtered version
6. Identification of highly variable gene products
7. Visualization of variant patterns through clustered heatmaps

The workflow can be applied to microbial whole-genome sequencing datasets of any haploid species (or diploid, with `--ploidy 2`) provided that:

- A reference FASTA file is available
- A GFF annotation file corresponding to the **same** reference FASTA is available and contains `gene` features, `locus_tag` annotations, and associated `product` annotations
- Paired-end genomic DNA FASTQ files are available for each isolate, assemblies (FASTA) are available for each genome, or both

Variants are aggregated from the joint VCF and summarized at the annotated gene-product level. Coverage is handled explicitly: for every sample and product, the fraction of the product's gene bases with sufficient depth is recorded, and products or cells that are not adequately covered can be removed or shown as blank (gray in the heatmaps) instead of being interpreted as "no variants".


## Software Requirements

### Variant Calling

BWA >= 0.7.17
SAMtools >= 1.21
BCFtools

These can come from your cluster's module system (see `--modules` and `--no-modules` below) or from anywhere on your `PATH`.


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
Paired-End FASTQ Files        Assemblies (FASTA)
            │                         │
            │                tiled into pseudo-reads
            │                (error-free .fastq)
            │                         │
            └───────────┬─────────────┘
                        ▼
              call_variants.sh
                        │
                        ▼
              variants.vcf
              variant_summary.txt
              coverage/<sample>.covered.bed
                        │
                        ▼
      build_product_variant_matrix.py   (+ genome.gff)
                        │
                        ▼
      variants_per_product_matrix.csv
      variants_per_product_matrix.filtered.csv
      product_coverage_matrix.csv
      variants_table.tsv
      filtered_genes.tsv
      products_no_variation.tsv
                        │
                        ├──────────────► plot_variant_product_heatmap.py
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
      plot_variant_product_heatmap.py
                        │
                        ▼
      Selected Product Heatmap
```

## Methods

Reads are aligned against a reference genome using BWA-MEM. Assemblies are converted into error-free paired-end pseudo-reads (150 bp reads from 400 bp fragments, tiled every 10 bp, about 30x depth) and aligned the same way, so reads and assemblies are called together in one joint VCF. Variants are jointly called using BCFtools and aggregated by annotated gene product using the supplied GFF annotation. At each variant site, only the alternate allele a sample actually carries is counted for that sample. Missing genotypes are ignored.

For each sample, the bases covered at a minimum depth are written as a BED file. These are converted to the fraction of each product's gene bases that are covered. Product-level variant counts are transformed into a sample-by-product matrix. Heatmaps display `log10(variants + 1)`, which is not rescaled per product, so products can be compared directly.


# Step 1: Variant Calling

Align reads and/or assemblies to a reference genome and perform joint variant calling.

## Usage

Reads and assemblies together:

```text
bash call_variants.sh \
    -r path/to/reference/genome.fna \
    -f path/to/trimmed_fastq \
    -a path/to/assemblies \
    -o results \
    -t 8
```

Reads only (leave out `-a`) or assemblies only (leave out `-f`) also work. At least one of `-f` or `-a` is required.

## Options

```text
-r, --reference     reference FASTA (required)
-f, --fastq         directory of paired-end reads (gzipped)
-a, --assemblies    directory of assemblies (.fasta .fa .fna .fas, optionally .gz)
-o, --output        output directory (default: results)
-t, --threads       threads (default: $SLURM_CPUS_PER_TASK if set, else 16)
-c, --min-cov-depth minimum depth for a base to count as covered in the coverage
                    files (default: 5)
    --ploidy N      ploidy for variant calling (default: 1)
    --modules "LIST" environment modules to load (default: $PIPELINE_MODULES,
                    else "bwa samtools bcftools")
-M, --no-modules    skip module loading and use tools already on PATH
    --r1-suffix S / --r2-suffix S
                    exact read-file suffixes, overriding auto-detection (give both)
-S, --summary-only  rebuild variant_summary.txt from an existing variants.vcf
-h, --help          show all options
```

Read file names are detected automatically: `_R1_001.trim.fastq.gz`, `_R1_001.fastq.gz`, `_R1.trim.fastq.gz`, `_R1.fastq.gz`, `_1.fastq.gz` (and the `.fq.gz` and R2 / `_2` equivalents). The sample name is the file name before that suffix.

Assembly sample names are the file names without their extension. An assembly and a read sample must not share a name.

#### How assemblies are handled

Assemblies are not aligned directly. Each contig is cut into overlapping, error-free paired-end **pseudo-reads** (150 bp reads from 400 bp fragments, one fragment every 10 bp, about 30x depth), which are then aligned with BWA-MEM exactly like real reads. This lets assemblies and sequencing reads go through the same alignment, coverage and variant-calling steps and appear together in one joint VCF. The pseudo-reads are deleted after alignment unless `KEEP_PSEUDO=1` is set.

Because a fragment is 400 bp, contigs shorter than 400 bp cannot produce any pseudo-reads and are skipped. Genes that lie only on such contigs will show as not covered for that assembly.

## Inputs

Reference FASTA
Paired-end FASTQ files and/or assembly FASTA files

## Outputs
```text
results/
├── variants.vcf
├── variant_summary.txt
├── <sample>.sorted.bam (+ .bai)
└── coverage/
    └── <sample>.covered.bed
```

Each `.covered.bed` lists merged intervals (contig, 0-based start, end) where depth is at least `-c`.

### Example Variant Summary
```text
Sample     Total_variants    SNPs    Indels    Multiallelic_sites
Isolate1   18364             17663   701       120
Isolate2   17902             17200   702       115
Isolate3   17486             16792   694       118
```


# Step 2: Build Product-Level Variant Matrix

Convert variant calls into a sample-by-product matrix using the reference annotation, and apply coverage-based filtering.

## Usage
```text
python3 build_product_variant_matrix.py \
    -v results/variants.vcf \
    -g reference/genome.gff \
    -o results \
    --coverage-dir results/coverage \
    --min-coverage 0.8 \
    --min-samples-covered 0.8 \
    --snps-only
```

`-v` and `-g` are required. The output directory defaults to `results`, and the coverage directory defaults to `<output>/coverage`. If the coverage directory does not exist, the coverage matrix and filtered outputs are skipped. The older positional form (`<vcf> <gff> [output_dir] [coverage_dir]`) is still accepted.

## Options

```text
-v, --vcf FILE            joint VCF from call_variants.sh (required)
-g, --gff FILE            GFF annotation for the same reference (required)
-o, --output DIR          output directory (default: results)
    --coverage-dir DIR    directory of <sample>.covered.bed files
                          (default: <output>/coverage)
--min-coverage F          a sample/product cell counts as covered when at least
                          this fraction of the product's gene bases are covered
                          (default: 0.8)
--min-samples-covered F   keep a product in the filtered matrix only if it is covered
                          in at least this fraction of samples (default: 1.0;
                          0 keeps all products)
--snps-only               count single-base substitutions only (ignore indels)
```

All coverage and variant-type filtering is set here. The plotting script does not apply its own thresholds.

## Inputs
```text
variants.vcf
genome.gff
coverage/<sample>.covered.bed
```

## Outputs
```text
results/
├── variants_table.tsv
├── variants_per_product_matrix.csv
├── product_coverage_matrix.csv
├── variants_per_product_matrix.filtered.csv
├── filtered_genes.tsv
└── products_no_variation.tsv
```

- `variants_per_product_matrix.csv`: all products, unfiltered
- `product_coverage_matrix.csv`: fraction of each product's gene bases covered, per sample
- `variants_per_product_matrix.filtered.csv`: products covered in enough samples; cells that are not adequately covered are blank
- `filtered_genes.tsv`: every product removed by the coverage filter, how many samples covered it, and in which samples it was not covered
- `products_no_variation.tsv`: products kept but identical in every covered sample (for example no variants anywhere), which the heatmaps do not plot

The matrix structure is:
```text
Sample      Cas1   Tet(Q)   CfxA
Isolate1      2       5       0
Isolate2      1       8       3
Isolate3      0     blank     1
```

A blank cell means the product was not adequately covered in that sample. It is not the same as 0 (covered, matches the reference).


# Step 3: Generate Curated Product List

Create a biologically curated list of products of interest. Products must be chosen from the `product` names annotated in the reference GFF, and spelled exactly as they appear there. Only products that also survive the coverage filter (that is, are present in the filtered matrix) are kept.

The list of products is defined inside `select_products.sh`. Edit it for your organism or project.

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
bash select_products.sh \
    -m results/variants_per_product_matrix.filtered.csv \
    -o results/selected_genes.tsv
```

Both `-m` (matrix) and `-o` (output file) are required. Products on the list that are not in the filtered matrix (for example removed by the coverage filter) are left out. Check `filtered_genes.tsv` to see why.

## Output
```text
results/
└── selected_genes.tsv
```


# Step 4: Generate Heatmaps

Plot the filtered matrix (`variants_per_product_matrix.filtered.csv`) so that coverage filtering applies. The heatmap script supports two modes:

## Option A. Most Variable Products

Plot the most variable products across all isolates.

### Default (Top 20 Products)
```text
python3 plot_variant_product_heatmap.py \
    results/variants_per_product_matrix.filtered.csv
```

### Custom Number of Products
```text
python3 plot_variant_product_heatmap.py \
    results/variants_per_product_matrix.filtered.csv \
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

The heatmap displays products ranked by variance of `log10(variants + 1)` across isolates. Products with identical values in all covered samples are not plotted.


## Option B. Curated Product Panel

Plot products listed in:
```text
selected_genes.tsv
```

Products with identical values in all covered samples are automatically removed before plotting.

### Usage
```text
python3 plot_variant_product_heatmap.py \
    results/variants_per_product_matrix.filtered.csv \
    results/selected_genes.tsv
```

Outputs:
```text
results/
├── variants_per_product_heatmap_selected.png
└── products_plotted_selected.tsv
```

### Example Output

![Selected Products](example_output/selected.png)

The heatmap displays a curated panel of biologically relevant gene products selected by the user.

### Reading the heatmaps

- Color is `log10(variants + 1)`: 0 for no variants, 1 for 9 variants, 2 for 99 variants.
- Gray cells are blank cells in the filtered matrix: the product was not adequately covered in that sample. Gray does not prove the gene is absent. It can also mean the gene is too divergent to map, or is missing from an assembly.


## Notes and Limitations

- The reference and the GFF must match. Genes that are not in the reference cannot be detected.
- Variants and coverage are pooled by product name, so products with many copies or generic names (for example `hypothetical protein`) combine all genes with that name.
- Repeated and mobile-element products (transposases, integrases) can show variants caused by mis-mapped reads.
- Reads and assemblies can differ in indel counts, often in homopolymers. Use `--snps-only` for comparisons across both.
- A product with the same value in every covered sample can be conserved (0 variants) or a fixed difference from the reference (the same variants in all samples).
- Sequencing depth still matters for isolate reads. `-c` and `--min-coverage` define what counts as covered.


## Pipeline Outputs
```text
variants.vcf
variant_summary.txt
coverage/<sample>.covered.bed

variants_table.tsv
variants_per_product_matrix.csv
product_coverage_matrix.csv
variants_per_product_matrix.filtered.csv
filtered_genes.tsv
products_no_variation.tsv
selected_genes.tsv

products_plotted_selected.tsv
products_plotted_topN.tsv

variants_per_product_heatmap_selected.png
variants_per_product_heatmap_topN.png
```
