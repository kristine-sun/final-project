!conda config --add channels bioconda
!conda install -y biopython muscle fasttree 

!muscle -super5 multiplesequences.fa -output multiplesequences.aligned.fa

!fasttree -nt < multiplesequences.aligned.fa > tree.nwk

from Bio import Phylo # Imports Bio.Phylo
from io import StringIO # Imports StringIO for reading strings
from matplotlib import pyplot as plt

# Plotting code 
fig = plt.figure()
ax=fig.add_axes([0,0,2,3])

# ToDo: One Line: Read in the tree in newick format
tree = Phylo.read("tree.nwk", "newick") # Example Newick format

# ToDo: One Line: Plot the tree
Phylo.draw(tree, axes=ax)

cluster = ["CY074482", "CY073812", "CY170966", "CY167875"]

from typing import List
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

alignment = AlignIO.read("multiplesequences.aligned.fa", "fasta")

def cluster_alignment(alignment: MultipleSeqAlignment, cluster: List[str]) -> MultipleSeqAlignment:
    """Creates a MultipleSequenceAlignment object for one cluster."""
    msa = MultipleSeqAlignment([])
    for record in alignment:
        if record.id in cluster:  # Filter by cluster membership
            msa.append(record)
       
    return msa

cluster_align = cluster_alignment(alignment, cluster)

import pandas as pd
import numpy as np

def build_cluster_dataframe(msa: MultipleSeqAlignment) -> pd.DataFrame:
    """
    Build the cluster-specific dataframe from the MSA.
    
    Returns
    -------
    A pandas DF with two columns ['%GT', '%AT'].
    Each row sums to 1.
    """
    gc_prop = []
    at_prop = []
    
    for i in range(0, msa.get_alignment_length()):
        alignment_column = str(msa[:,i])
        gc_tmp = 0
        at_tmp = 0
        
        for c in alignment_column:
            #ToDo: 2 Lines: Check if c equals a or t, and if so, increase at_tmp by 1
            if c in 'aAtT':
                at_tmp += 1
                
            #ToDo: 2 Lines: Check if c equals g or c, and if so, increase gc_tmp by 1
            if c in 'gGcC':
                gc_tmp += 1

        gc_prop.append(gc_tmp)
        at_prop.append(at_tmp)
        
        
    gc_prop = np.array(gc_prop, dtype=float)
    at_prop = np.array(at_prop, dtype=float)
    
    total = gc_prop + at_prop
    gc_fraction = gc_prop / total
    at_fraction = at_prop / total
    
    df = pd.DataFrame()
    
    df["%GC"] = gc_fraction
    df["%AT"] = at_fraction
    
    return df

df_cluster = build_cluster_dataframe(cluster_align)

df_cluster.head(50).plot(kind="bar", stacked=True, title="Cluster(First 50 Positions)")