#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate heatmaps from a sample-by-product variant matrix
produced by build_product_variant_matrix.py.

Modes

1. Variable-product mode

   python3 plot_selected_products_heatmap.py \
       variants_per_product_matrix.csv

   Uses TOP_N=20 by default.

2. Variable-product mode (user-defined TOP_N)

   python3 plot_selected_products_heatmap.py \
       variants_per_product_matrix.csv \
       50

   Identifies and plots the top 50 most variable products.

3. Selected-product mode

   python3 plot_selected_products_heatmap.py \
       variants_per_product_matrix.csv \
       selected_genes.tsv

   Plots all products listed in selected_genes.tsv and removes
   invariant products prior to visualization.

Examples

   # Top 20 most variable products (default)
   python3 plot_selected_products_heatmap.py \
       results/variants_per_product_matrix.csv

   # Top 50 most variable products
   python3 plot_selected_products_heatmap.py \
       results/variants_per_product_matrix.csv \
       50

   # Curated product panel
   python3 plot_selected_products_heatmap.py \
       results/variants_per_product_matrix.csv \
       results/selected_genes.tsv
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.stats import zscore

# ------------------------------------------------------------
# Step 1: Parse command-line arguments
# ------------------------------------------------------------
if len(sys.argv) < 2:
    sys.exit(
        "Usage:\n"
        "  python3 plot_selected_products_heatmap.py "
        "<variants_per_product_matrix.csv>\n"
        "  python3 plot_selected_products_heatmap.py "
        "<variants_per_product_matrix.csv> <TOP_N>\n"
        "  python3 plot_selected_products_heatmap.py "
        "<variants_per_product_matrix.csv> <selected_genes.tsv>"
    )

matrix_file = sys.argv[1]

selected_file = None
TOP_N = 20

if len(sys.argv) > 2:

    # Numeric argument = variable-product mode
    if sys.argv[2].isdigit():
        TOP_N = int(sys.argv[2])

    # File argument = selected-gene mode
    else:
        selected_file = sys.argv[2]

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

print("\n--- Heatmap generation from variants_per_product_matrix.csv ---")
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
        pd.read_csv(selected_file, header=None)[0]
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
# Step 5: Z-score normalization
# ------------------------------------------------------------
warnings.filterwarnings(
    "ignore",
    message="Precision loss occurred in moment calculation*"
)

matrix_zscore = variants_per_product_matrix.apply(
    zscore,
    axis=0,
    nan_policy="omit"
)

matrix_zscore = (
    matrix_zscore
    .replace([np.inf, -np.inf], 0)
    .fillna(0)
)

# ------------------------------------------------------------
# Step 6: Build plotting matrix
# ------------------------------------------------------------

if selected_file:

    # Plot ALL selected products
    plotted_products = list(matrix_zscore.columns)

    matrix_plot = matrix_zscore.T

    print(
        f"\nPlotting all {len(plotted_products)} "
        f"selected products."
    )

    figure_title = (
        "Hierarchical clustering of selected products"
    )

else:

    # Plot only the most variable products
    product_variance = matrix_zscore.var(axis=0)

    plotted_products = (
        product_variance
        .sort_values(ascending=False)
        .head(min(TOP_N, len(product_variance)))
        .index
        .tolist()
    )

    matrix_plot = matrix_zscore[plotted_products].T

    print(
        f"\nSelected {len(plotted_products)} most variable "
        f"products out of "
        f"{matrix_zscore.shape[1]} total products."
    )

    figure_title = (
        f"Hierarchical clustering of top "
        f"{len(plotted_products)} most variable products"
    )


# ------------------------------------------------------------
# Step 6b: Remove invariant selected products
# ------------------------------------------------------------
if selected_file:

    row_var = matrix_plot.var(axis=1)

    matrix_plot = matrix_plot.loc[row_var > 0]

    plotted_products = list(matrix_plot.index)

    print(
        f"Retained {len(plotted_products)} variable selected "
        f"products after removing invariant products."
    )

# ------------------------------------------------------------
# Final cleanup
# ------------------------------------------------------------
matrix_plot = (
    matrix_plot
    .replace([np.inf, -np.inf], 0)
    .fillna(0)
    .astype(float)
)

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
    cmap="coolwarm",
    figsize=(fig_width, fig_height),
    metric="euclidean",
    method="average",
    standard_scale=0,
    row_cluster=True,
    col_cluster=True,
    yticklabels=True,
    cbar_pos=(0.02, 0.80, 0.03, 0.18)
)

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
    "Z-score",
    rotation=270,
    labelpad=15,
    fontsize=10
)

plt.savefig(
    heatmap_file,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print(f"\nHeatmap saved to: {heatmap_file}")