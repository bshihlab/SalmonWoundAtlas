# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 19:29:30 2026

@author: shihb
"""
import math
import os
import sys
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib as mpl
import seaborn as sns
import pandas as pd

from anndata import concat
from matplotlib import rcParams
from scipy.sparse import csr_matrix

import cell2location


working_dir="L:/ShihLab/projects/2024/20240226_rose_wounds/"
os.chdir(working_dir)

rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs

results_folder = './analysis/cell2location/'
sp_data_folder = './analysis/spaceranger/'

# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'
abundance_folder = f'{results_folder}/deconvoluted_abundance'

if not os.path.exists(abundance_folder):
    os.mkdir(abundance_folder)

### read in spaceranger data
all_samples = ["1-SI-TT-A3", "2-SI-TT-B3", "3-SI-TT-C3", "4-SI-TT-D3" ]

samples = {}
adatas = {}
slides = []
#all_samples = all_samples[0:5]
for sample_id in all_samples:
    #ignore_barcode = read.csv()
    # current visium
    current_dir = "analysis/spaceranger/" + sample_id + "/outs/"
    adata = sc.read_visium(current_dir )
    adata.obs['sample'] = sample_id
    adata.var['SYMBOL'] = adata.var_names
    adata.var.rename(columns={'gene_ids': 'ENSEMBL'}, inplace=True)
    adata.var_names = adata.var['ENSEMBL']
    adata.var.drop(columns='ENSEMBL', inplace=True)
    adata.var_names_make_unique()

    # Calculate QC metrics
    adata.X = adata.X.toarray()
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    adata.X = csr_matrix(adata.X)
    adata.var['mt'] = [gene.startswith('mt-') for gene in adata.var['SYMBOL']]
    adata.obs['mt_frac'] = adata[:, adata.var['mt'].tolist()].X.sum(1).A.squeeze()/adata.obs['total_counts']
    adata.obs['batch_key'] = sample_id

    # ignore spots if the sample id is 
    # if keep spots file are present
    current_keep_spot_fp = current_dir + "/keep_spots.csv"
    if os.path.exists(current_keep_spot_fp):
        keep_spots = pd.read_csv(current_keep_spot_fp)
        adata = adata[adata.obs_names.isin(keep_spots.Barcode)].copy()

    # add sample name to obs names
    adata.obs["sample"] = [str(i) for i in adata.obs['sample']]
    adata.obs_names = adata.obs["sample"] \
                          + '_' + adata.obs_names
    adata.obs.index.name = 'spot_id'


    # save it for later
    adatas[sample_id] = adata
    slides.append(adata)

# Combine anndata objects together
adata_vis = concat(
        slides,
        merge="unique",
        uns_merge="unique",
        label="batch",
        keys=all_samples, 
        index_unique=None)


cell2location.models.Cell2location.setup_anndata(
    adata=adata_vis,
    batch_key="batch")


#### Read in single cell reference
all_sc = []
sc_fps = os.listdir("analysis/cellbender/")
sc_fps = [val for val in sc_fps if val != "SRR24471129"] 
seurat_metadata = pd.read_csv(os.path.join("analysis/cell2location/metadata.csv"), index_col=0)
for sample in sc_fps:
    adata = sc.read_10x_h5("analysis/cellbender/" + sample + "/" + sample + ".h5")
    adata.obs_names = sample + "_" + adata.obs_names
    adata.obs['sample'] = sample_id
    adata.var['SYMBOL'] = adata.var_names
    adata.var.rename(columns={'gene_ids': 'ENSEMBL'}, inplace=True)
    adata.var_names = adata.var['ENSEMBL']
    adata.var.drop(columns='ENSEMBL', inplace=True)
    adata.var_names_make_unique()
    
    # keep only cells that are present in the seurat data
    adatas_keep_idx = [cell_id for cell_id in seurat_metadata.index if cell_id in adata.obs_names]
    adata = adata[adatas_keep_idx,:]
    # organise cell type annotation to match cell ids
    seurat_metadata_ordered = seurat_metadata.loc[adatas_keep_idx,:]
    adata.obs['celltype_Stdeconvolution'] = seurat_metadata_ordered.celltype_Stdeconvolution
    adata.obs['Sample'] = sample
    all_sc.append(adata)
    
adata_sc = concat(
        all_sc,
        merge="unique",
        label="batch",
        keys=sc_fps, 
        index_unique=None)
    

#### cell2location single cell reference setup
from cell2location.utils.filtering import filter_genes
selected = filter_genes(adata_sc, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

# filter the object
adata_ref = adata_sc[:, selected].copy()

# prepare anndata for the regression model
from cell2location.models import RegressionModel
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref, 
                        # 10X reaction / sample / batch
                        batch_key='batch', 
                        # cell type, covariate used for constructing signatures
                        labels_key='celltype_Stdeconvolution'
                       )

# create the regression model
mod = RegressionModel(adata_ref) 

# view anndata_setup as a sanity check
mod.view_anndata_setup()

mod.train(max_epochs=250, accelerator='gpu')

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'accelerator': 'gpu'}
)

# Save model
mod.save(f"{ref_run_name}", overwrite=True)

# Save anndata object with results
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref.write(adata_file)
adata_file

##Examine QC plots.
#Reconstruction accuracy to assess if there are any issues with inference. This 2D histogram plot should have most observations along a noisy diagonal.
#The estimated expression signatures are distinct from mean expression in each cluster because of batch effects. For scRNA-seq datasets which do not suffer from batch effect (this dataset does), cluster average expression can be used instead of estimating signatures with a model. When this plot is very different from a diagonal plot (e.g. very low values on Y-axis, density everywhere) it indicates problems with signature estimation.
mod.plot_QC()

## load previously generated model and output h5d
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref = sc.read_h5ad(adata_file)
mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", adata_ref)


## export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' 
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}' 
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']
inf_aver.iloc[0:5, 0:5]
inf_aver.to_csv(f"{ref_run_name}/inf_aver.csv")

#### spatial mapping
# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")

# create and train the model
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver, 
    # the expected average cell abundance: tissue-dependent 
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=30,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=20
) 
mod.view_anndata_setup()

# train cell 2 location
mod.train(max_epochs=30000, 
          # train using full data (batch_size=None)
          batch_size=None, 
          # use all data points in training because 
          # we need to estimate cell abundance at all locations
          train_size=1,
          accelerator='gpu',
         )

# plot ELBO loss history during training, removing first 100 epochs from the plot
mod.plot_history(1000)
plt.legend(labels=['full data training']);


# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'accelerator': 'gpu'}
)

# Save model
mod.save(f"{run_name}", overwrite=True)

# mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# Save anndata object with results
adata_file = f"{run_name}/sp.h5ad"
adata_vis.write(adata_file)
adata_file

# The model and output h5ad can be loaded later like this:
#adata_file = f"{run_name}/sp.h5ad"
#adata_vis = sc.read_h5ad(adata_file)
#mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

#### check QC
# Examine reconstruction accuracy to assess 
# if there are any issues with mapping. 
# The plot should be roughly diagonal, 
# strong deviations will signal problems that need to be investigated.
mod.plot_QC()

fig = mod.plot_spatial_QC_across_batches()


# add 5% quantile, representing confident cell abundance, 'at least this amount is present', 
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']

# select one slide
# print graphs for each cell type
outdir = "analysis/cell2location/figure/"
if not os.path.isdir(outdir):
    os.makedirs(outdir)


from cell2location.utils import select_slide

for slide_id in all_samples:
    slide = select_slide(adata_vis, slide_id)
    slide.obs.columns    
    # plot in spatial coordinates
    all_cell_type = seurat_metadata.celltype_Stdeconvolution.unique().tolist()
    all_cell_type = [val for val in all_cell_type if val in slide.obs]
    all_cell_type.sort()
    
    # write out the deconvoluted cell abundance so the plots can be made using R
    # This is so that each sample can be turned to the correct orientation
    slide.obs[all_cell_type].to_csv(f"{abundance_folder}/{slide_id}.csv")

    # plot 8 cell types at a time
    for i in range(math.ceil(len(all_cell_type)/8)):
        start_idx = i*8
        end_idx = start_idx + 8
        end_idx = len(all_cell_type) if len(all_cell_type) < end_idx else end_idx
        with mpl.rc_context({'axes.facecolor':  'white',
                             'figure.figsize': [4.5, 5]}):
            out_fp = working_dir + "/" + outdir + slide_id + "_" + str(i), ".png"
            #sc.set_figure_params(facecolor="grey")

            sc.pl.spatial(slide, cmap='magma',
                          # show first 8 cell types
                          color=all_cell_type[start_idx:end_idx], 
                          ncols=4, size=1.3, 
                          img_key='hires',
                          # limit color scale at 99.2% quantile of cell abundance
                          vmin=0, vmax='p99.2' 
                         )

