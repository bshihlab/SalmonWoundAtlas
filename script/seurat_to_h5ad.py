# -*- coding: utf-8 -*-
"""
Created on Sun May 24 04:38:09 2026

@author: shihb

Save integrated seurat as h5ad
"""
import os
import sys
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib as mpl
import seaborn as sns
import pandas as pd
import anndata as ad

from scipy.io import mmread
from scipy.sparse import csr_matrix

working_dir="L:/ShihLab/projects/2024/20240226_rose_wounds/"
os.chdir(working_dir)

input_dir = "analysis/salmon_wound_atlas/rds/"
out_fp = input_dir + "merged_srat.h5ad"

adata_bc=pd.read_csv(input_dir+'h5ad_barcodes.csv',header=0).set_index('barcode', drop=False)
adata_features=pd.read_csv(input_dir+'gene.csv',header=0).set_index('salmon.geneid', drop=False)
sample_obs = pd.read_csv(input_dir + "h5ad_annotation_brief.csv",header=0)
umap = pd.read_csv(input_dir + "h5ad_embedding.csv",header=0).to_numpy()
RNA_counts = mmread(input_dir + "h5ad_count.mx").transpose().tocsc().astype('float32')
RNA_normalised = mmread(input_dir + "h5ad_data.mx").transpose().tocsc().astype('float32')

layers_dict = {'RNA': RNA_counts}

adata = ad.AnnData(
    X=layers_dict['RNA'], 
    obs=sample_obs, 
    var=adata_features, 
    layers=layers_dict)

adata.layers["raw"] = RNA_counts
adata.layers["X"] = RNA_normalised
adata.obsm['X_umap'] = umap

# plot umap to check 
sc.pl.umap(adata, color=["named_celltypes"])

# write it out as h5ad
adata.write_h5ad(out_fp)
