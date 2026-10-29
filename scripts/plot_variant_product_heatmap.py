#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate heatmaps from a sample-by-product variant matrix
produced by build_product_variant_matrix.py.

Modes

1. Variable-product mode

   python3 plot_variant_product_heatmap.py \
       variants_per_product_matrix.csv

   Uses TOP_N=20 by default.

2. Variable-product mode (user-defined TOP_N)

   python3 plot_variant_product_heatmap.py \
       variants_per_product_matrix.csv \
       50

   Identifies and plots the top 50 most variable products.

3. Selected-product mode

   python3 plot_variant_product_heatmap.py \
       variants_per_product_matrix.csv \
       selected_genes.tsv

   Plots all products listed in selected_genes.tsv and removes
   invariant products prior to visualization.

Color scale: log10(variants + 1), not rescaled per product.
Coverage is NOT decided here. Set --min-coverage, --min-samples-covered and
--snps-only when running build_product_variant_matrix.py, and plot its
variants_per_product_matrix.filtered.csv. Blank cells in that matrix
(gene not adequately covered in that sample) are drawn gray, so "no
variants" and "not covered" stay distinguishable. Top-N mode ranks products
by variance of the log-transformed counts over the covered samples.

Examples

   # Top 20 most variable products (default)
   python3 plot_variant_product_heatmap.py \
       results/variants_per_product_matrix.filtered.csv

   # Top 50 most variable products
   python3 plot_variant_product_heatmap.py \
       results/variants_per_product_matrix.filtered.csv \
       50

   # Curated product panel
   python3 plot_variant_product_heatmap.py \
       results/variants_per_product_matrix.filtered.csv \
       results/selected_genes.tsv
"""

import os
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

GRAY = "#c8c8c8"

import argparse

ap = argparse.ArgumentParser(
    description="Heatmap of variants per product (blank cells drawn gray)."
)
ap.add_argument("matrix")
ap.add_argument(
    "mode", nargs="?", default=None,
    help="TOP_N (integer, default 20) or a selected_genes.tsv file"
)
args = ap.parse_args()

matrix_file = args.matrix
selected_file = None
TOP_N = 20

if args.mode is not None:
    if args.mode.isdigit():
        TOP_N = int(args.mode)
    else:
        selected_file = args.mode

if not os.path.exists(matrix_file):
    sys.exit(f"Error: {matrix_file} not found.")

if selected_file and not os.path.exists(selected_file):
    sys.exit(f"Error: {selected_file} not found.")

results_dir = os.path.dirname(matrix_file) or "."

# ------------------------------------------------------------
# Output filenames
# ------------------------------------------------------------

if selected_file:

    heatmap_file = os.path.join(
        results_dir,
        "variants_per_product_heatmap_selected.png"
    )

    plotted_products_file = os.path.join(
        results_dir,
        "products_plotted_selected.tsv"
    )

else:

    heatmap_file = os.path.join(
        results_dir,
        f"variants_per_product_heatmap_top{TOP_N}.png"
    )

    plotted_products_file = os.path.join(
        results_dir,
        f"products_plotted_top{TOP_N}.tsv"
    )

print("\n--- Heatmap generation from a variants-per-product matrix ---")
print(f"Matrix file: {matrix_file}")

# ------------------------------------------------------------
# Step 2: Load the existing variant matrix
# ------------------------------------------------------------
variants_per_product_matrix = pd.read_csv(
    matrix_file,
    index_col=0
)

print(
    f"Loaded matrix: "
    f"{variants_per_product_matrix.shape[0]} samples x "
    f"{variants_per_product_matrix.shape[1]} products"
)

# ------------------------------------------------------------
# Step 3: Filter by selected products (if provided)
# ------------------------------------------------------------
if selected_file:

    if not os.path.exists(selected_file):
        sys.exit(f"Error: {selected_file} not found.")

    selected_genes = (
        pd.read_csv(selected_file, header=None, sep="\t")[0]
        .dropna()
        .unique()
        .tolist()
    )

    print(
        f"\nLoaded {len(selected_genes)} selected products "
        f"from {selected_file}"
    )

    selected_in_matrix = [
        g for g in selected_genes
        if g in variants_per_product_matrix.columns
    ]

    missing_genes = [
        g for g in selected_genes
        if g not in variants_per_product_matrix.columns
    ]

    if missing_genes:

        print(
            f"Warning: {len(missing_genes)} products "
            f"not found in matrix columns "
            f"(showing up to 10):"
        )

        for mg in missing_genes[:10]:
            print(f"  - {mg}")

    if not selected_in_matrix:
        sys.exit(
            "No matching products found in the matrix. "
            "Check your selected_genes.tsv file."
        )

    variants_per_product_matrix = (
        variants_per_product_matrix[selected_in_matrix]
    )

    print(
        f"Filtered to {len(selected_in_matrix)} "
        f"products present in the matrix."
    )

else:

    print(
        "\nNo selected_genes.tsv provided. "
        "Using all products."
    )

# ------------------------------------------------------------
# Step 4: Sanity checks
# ------------------------------------------------------------
if variants_per_product_matrix.shape[1] == 0:
    sys.exit("No products remain after filtering.")

# ------------------------------------------------------------
# Step 5: Log-transform counts
#
# Color reflects the actual number of variants: log10(count + 1).
# (0 -> 0, 9 -> 1, 99 -> 2). Values are NOT rescaled per product,
# so a 1-vs-0 difference looks small and a 40-vs-0 difference
# looks large, and products can be compared directly.
# ------------------------------------------------------------
matrix_log = np.log10(variants_per_product_matrix.astype(float) + 1)

# ------------------------------------------------------------
# Step 6: Build plotting matrix
# ------------------------------------------------------------

if selected_file:

    # Plot ALL selected products
    plotted_products = list(matrix_log.columns)

    matrix_plot = matrix_log.T

    print(
        f"\nPlotting all {len(plotted_products)} "
        f"selected products."
    )

    figure_title = (
        "Hierarchical clustering of selected products"
    )

else:

    # Plot only the most variable products
    eligible = [
        c for c in matrix_log.columns
        if matrix_log[c].nunique(dropna=True) > 1
    ]

    if not eligible:
        sys.exit("No products with variation to plot.")

    product_variance = matrix_log[eligible].var(axis=0).dropna()

    plotted_products = (
        product_variance
        .sort_values(ascending=False)
        .head(min(TOP_N, len(product_variance)))
        .index
        .tolist()
    )

    matrix_plot = matrix_log[plotted_products].T

    print(
        f"\nSelected {len(plotted_products)} most variable "
        f"products out of "
        f"{len(product_variance)} eligible products."
    )

    figure_title = (
        f"Hierarchical clustering of top "
        f"{len(plotted_products)} most variable products"
    )


# ------------------------------------------------------------
# Step 6b: Remove invariant selected products
# ------------------------------------------------------------
if selected_file:

    # A product varies if it has more than one distinct value across its
    # covered samples. (Comparing variance > 0 is unreliable: identical
    # floats can give a tiny nonzero variance, ~1e-33.)
    varies = matrix_plot.nunique(axis=1, dropna=True) > 1

    matrix_plot = matrix_plot.loc[varies]

    plotted_products = list(matrix_plot.index)

    print(
        f"Retained {len(plotted_products)} variable selected "
        f"products after removing invariant products."
    )

# ------------------------------------------------------------
# Final cleanup
# ------------------------------------------------------------
matrix_plot = matrix_plot.astype(float).replace([np.inf, -np.inf], np.nan)

# Blank cells = not adequately covered (set in build_product_variant_matrix.py).
# Remember them for the gray mask; for clustering only, fill each with its
# product's mean over the covered samples (0 if none).
blank = matrix_plot.isna()
matrix_plot = matrix_plot.apply(lambda r: r.fillna(r.mean()), axis=1).fillna(0)

print(
    f"Heatmap dimensions: "
    f"{matrix_plot.shape[0]} products x "
    f"{matrix_plot.shape[1]} samples"
)

print(
    "Contains NaN:",
    np.isnan(matrix_plot.values).any()
)

print(
    "Contains Inf:",
    np.isinf(matrix_plot.values).any()
)

print(
    "All values finite:",
    np.isfinite(matrix_plot.values).all()
)

if matrix_plot.shape[0] < 2:
    sys.exit(
        f"Need at least 2 products for clustering. "
        f"Current shape: {matrix_plot.shape}"
    )

if matrix_plot.shape[1] < 2:
    sys.exit(
        f"Need at least 2 samples for clustering. "
        f"Current shape: {matrix_plot.shape}"
    )

# ------------------------------------------------------------
# Step 6c: Gray mask for blank (not covered) cells
# ------------------------------------------------------------
mask = None

if blank.values.any():
    mask = blank
    print(
        f"Blank cells: {int(blank.values.sum())} of {blank.size} "
        f"(not covered; shown gray)."
    )

# ------------------------------------------------------------
# Step 7: Dynamic figure sizing
# ------------------------------------------------------------
num_samples = matrix_plot.shape[1]
num_products = matrix_plot.shape[0]

fig_width = max(12, num_samples * 0.3)
fig_height = max(8, num_products * 0.25)

ytick_fontsize = max(
    4,
    min(10, int(300 / num_products))
)

# ------------------------------------------------------------
# Step 8: Save plotted product list
# ------------------------------------------------------------
pd.Series(
    plotted_products
).to_csv(
    plotted_products_file,
    sep="\t",
    index=False,
    header=False
)

print(
    f"Saved list of {len(plotted_products)} products to "
    f"{plotted_products_file}"
)

# ------------------------------------------------------------
# Step 9: Plot heatmap
# ------------------------------------------------------------
sns.set(style="white")

print("\nOutput files:")
print(f"Heatmap: {heatmap_file}")
print(f"Products: {plotted_products_file}")

cg = sns.clustermap(
    matrix_plot,
    cmap="YlOrRd",
    vmin=0,
    mask=mask,
    figsize=(fig_width, fig_height),
    metric="euclidean",
    method="average",
    row_cluster=True,
    col_cluster=True,
    yticklabels=True,
    cbar_pos=(0.02, 0.80, 0.03, 0.18)
)

# Masked (low-coverage) cells show the axes background
if mask is not None:
    cg.ax_heatmap.set_facecolor(GRAY)

# Hide dendrogram panels
cg.ax_row_dendrogram.set_visible(False)
cg.ax_col_dendrogram.set_visible(False)

# Reposition heatmap upward to reclaim the
# blank space left by the hidden dendrograms
heat_pos = cg.ax_heatmap.get_position()

cg.ax_heatmap.set_position([
    heat_pos.x0,
    heat_pos.y0 + 0.05,
    heat_pos.width,
    heat_pos.height + 0.08
])

# Format labels
plt.setp(
    cg.ax_heatmap.get_yticklabels(),
    rotation=0,
    ha="right",
    fontsize=ytick_fontsize
)

plt.setp(
    cg.ax_heatmap.get_xticklabels(),
    rotation=90,
    ha="center",
    fontsize=8
)

cg.ax_heatmap.yaxis.set_label_position("left")
cg.ax_heatmap.yaxis.set_ticks_position("left")

# Title
cg.fig.suptitle(
    figure_title,
    y=0.97
)

# Move colorbar to right side of heatmap
# and slightly lower to avoid overlap with title
pos = cg.ax_heatmap.get_position()

cg.cax.set_position([
    pos.x1 + 0.02,
    pos.y0 - 0.03,
    0.02,
    pos.height
])

cg.cax.set_ylabel(
    "log10(variants + 1)",
    rotation=270,
    labelpad=15,
    fontsize=10
)

if mask is not None:
    cg.fig.legend(
        handles=[
            Patch(
                facecolor=GRAY,
                edgecolor="0.4",
                label="Not covered (blank in matrix)"
            )
        ],
        loc="lower center",
        bbox_to_anchor=(0.45, -0.04),
        frameon=False
    )

plt.savefig(
    heatmap_file,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print(f"\nHeatmap saved to: {heatmap_file}")