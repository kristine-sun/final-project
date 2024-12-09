import subprocess
from Bio import AlignIO, Phylo
from Bio.Align import MultipleSeqAlignment
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
from typing import List

# Step 1: Install Required Tools and Dependencies
# It’s better to use a conda environment separately rather than installing via subprocess.
# Example:
# conda create -n bio_env biopython muscle fasttree
# conda activate bio_env

# Step 2: Perform Multiple Sequence Alignment and Phylogenetic Tree Construction
# Use Biopython directly for alignment instead of calling MUSCLE via subprocess.
# For MUSCLE, we use Bio.Align.Applications instead.

from Bio.Align.Applications import MuscleCommandline

muscle_cline = MuscleCommandline(input="multiplesequences.fa", out="multiplesequences.aligned.fa")
muscle_cline()

# Step 3: Construct the Phylogenetic Tree
from Bio.Phylo.Applications import FastTreeCommandline

fasttree_cline = FastTreeCommandline(input="multiplesequences.aligned.fa", out="tree.nwk")
fasttree_cline()

# Step 4: Filter Alignment for Cluster
alignment = AlignIO.read("multiplesequences.aligned.fa", "fasta")
cluster = ["CY074482", "CY073812", "CY170966", "CY167875"]

def cluster_alignment(alignment: MultipleSeqAlignment, cluster: List[str]) -> MultipleSeqAlignment:
    msa = MultipleSeqAlignment([record for record in alignment if record.id in cluster])
    return msa

cluster_align = cluster_alignment(alignment, cluster)

# Step 5: Build Cluster-Specific DataFrame
def build_cluster_dataframe(msa: MultipleSeqAlignment) -> pd.DataFrame:
    gc_prop, at_prop = [], []
    for i in range(msa.get_alignment_length()):
        alignment_column = str(msa[:, i])
        gc_tmp = sum(1 for c in alignment_column if c in 'gGcC')
        at_tmp = sum(1 for c in alignment_column if c in 'aAtT')
        gc_prop.append(gc_tmp)
        at_prop.append(at_tmp)
    
    gc_prop, at_prop = np.array(gc_prop, dtype=float), np.array(at_prop, dtype=float)
    total = gc_prop + at_prop
    gc_fraction, at_fraction = gc_prop / total, at_prop / total
    
    df = pd.DataFrame({"%GC": gc_fraction, "%AT": at_fraction})
    return df

df_cluster = build_cluster_dataframe(cluster_align)

# Step 6: Generate BED File
def create_bed_file(df: pd.DataFrame, output_file: str):
    with open(output_file, "w") as bed_file:
        for i, row in df.iterrows():
            bed_file.write(f"chr1\t{i}\t{i+1}\tGC_Content\t{row['%GC']:.2f}\n")
            bed_file.write(f"chr1\t{i}\t{i+1}\tAT_Content\t{row['%AT']:.2f}\n")

create_bed_file(df_cluster, "cluster_alignment.bed")

# Step 7: Generate WIG File
def create_wig_file(df: pd.DataFrame, output_file: str):
    with open(output_file, "w") as wig_file:
        wig_file.write('track type=wiggle_0 name="GC_Content" description="GC percentage content"\n')
        wig_file.write("fixedStep chrom=chr1 start=1 step=1\n")
        for gc_content in df["%GC"]:
            wig_file.write(f"{gc_content:.2f}\n")

create_wig_file(df_cluster, "cluster_alignment.wig")

print("BED and WIG files have been created.")
