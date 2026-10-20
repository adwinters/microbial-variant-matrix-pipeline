#!/bin/bash

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <results_directory>"
    exit 1
fi

RESULTS_DIR=$1

python3 - <<EOF
import pandas as pd

results_dir = "${RESULTS_DIR}"

df = pd.read_csv(
    f"{results_dir}/variants_per_product_matrix.csv",
    index_col=0
)

# Curated products of biological interest.
# Modify this list for your organism or project.
selected_genes = [
    "replication initiation protein",
    "MarR family winged helix-turn-helix transcriptional regulator",
    "type II CRISPR-associated endonuclease Cas1",
    "efflux RND transporter periplasmic adaptor subunit",
    "sodium:proton antiporter NhaD",
    "FtsB family cell division protein",
    "CRISPR-associated endonuclease Cas2",
    "CfxA family broad-spectrum class A beta-lactamase",
    "DNA primase",
    "type II CRISPR RNA-guided endonuclease Cas9",
    "LrgB family protein",
    "MATE family efflux transporter",
    "penicillin-binding protein",
    "hemolysin family protein",
    "thiol-activated cytolysin family protein",
    "hydrogen peroxide-inducible genes activator",
    "TetR/AcrR family transcriptional regulator",
    "tetracycline resistance ribosomal protection protein Tet(Q)",
    "FtsK/SpoIIIE family DNA translocase",
    "DNA mismatch repair protein MutS",
    "DNA adenine methylase",
    "sialidase family protein",
    "BlaI/MecI/CopY family transcriptional regulator",
    "LuxR C-terminal-related transcriptional regulator",
    "LysE family translocator"
]

filtered = [
    c for c in df.columns
    if c in selected_genes
]

outfile = f"{results_dir}/selected_genes.tsv"

pd.Series(filtered).to_csv(
    outfile,
    index=False,
    header=False
)

print(f"Total selected genes: {len(filtered)}")
print(f"Output file: {outfile}")
EOF