working_dir <- "L:/ShihLab/projects/2024/20240226_rose_wounds/"
setwd(working_dir)

####-- installation
if(FALSE){
	remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
	remotes::install_github('satijalab/azimuth', ref = 'master')
	remotes::install_github("satijalab/seurat-wrappers", "seurat5", quiet = TRUE)
	remotes::install_github('satijalab/seurat-wrappers')
	install.packages("scCustomize")
	install.packages("BiocManager")
	install.packages("devtools")
	install.packages('reticulate')
	install.packages('HGNChelper')
	install.packages('openxlsx')
	install.packages("RColorBrewer")
	BiocManager::install("scDblFinder")
	BiocManager::install("scPipe")
	BiocManager::install("biomaRt")
	BiocManager::install("slingshot")
	BiocManager::install("celda")
	BiocManager::install("scater")
	BiocManager::install("harmony")
	BiocManager::install("anndataR")
	remotes::install_github("vertesy/Seurat.utils")
	BiocManager::install("orthogene")
	reticulate::conda_install('r-reticulate', 'plotly')
	BiocManager::install("tricycle")
	devtools::install_github("haotian-zhuang/findPC")
}


####-- Load libraries and environment setting
library(Seurat)
library(findPC)
#library(SeuratDisk)
library(ggplot2)
library(scCustomize)
library(Matrix)
library(scDblFinder)
library(SingleCellExperiment)
library(tidyr)
library(viridis)
library(dplyr)
library(HGNChelper)
library(openxlsx)
library(phateR)
library(slingshot)
library(reticulate)
library(RColorBrewer)
library(tidymodels)
library(ranger)
library(rsample)
library(pbapply)
library(pheatmap)
library(scater)
library(stringr)
library(harmony)
library(SeuratWrappers)
library(anndataR)
library(scPipe)
library(celda)
library(plotly)
library(gam)
library(viridis)
library( gridExtra )
library(orthogene)
library(tricycle)
library(readxl)
library(openxlsx)


#use_condaenv("C:/Users/shihb/.conda/envs/20240429_salmon")

## Set up future module parallelisation
#plan("multicore", workers = 8)
options(future.globals.maxSize= 10000 * 1024^2)
#future_options(seed=TRUE)


####-- Organise input/output folders/files
## Input
cellbender_dir <- "analysis/cellbender"
cell_cluster_label <- "annotation/cluster/_cluster_id_to_celltype_rationale.csv"
cell_cluster_col_label <- "annotation/cluster/_clusterid_to_celltype_col.xlsx"
sample_list <- list.files(cellbender_dir)
sample_list <- sample_list[!sample_list %in% "SRR24471129"] # remove fin SAMPLE
previous_annotation_gene_list <- "analysis/_archive/ribosomal_cluster_keep_20250328/stats/all_markers_finalclustering_unique.csv"
gene_list_proliferation <- "data/cell_type_identification/proliferation_spatial_gene_matched.csv"


## Output
out_dir <- "analysis/salmon_wound_atlas"
out_rds_doublets <- paste0(out_dir, "/rds/merged_srat_doublets_20240828.rds")
out_plots <-paste0(out_dir, "/plots")
out_initial <-paste0(out_dir, "/initial_analysis")
out_tables <-paste0(out_dir, "/tables")
out_plots_qc <- paste0(out_dir , "/plots/qc")
out_tables_qc <- paste0(out_dir , "/tables/qc")
out_plots_umap <- paste0(out_dir, "/plots/umap")
out_plots_umap_gene <- paste0(out_dir, "/plots/umap/_gene")
out_plots_umap_individual <- paste0(out_dir, "/plots/umap/individual")
out_rds <- paste0(out_dir, "/rds")
out_dir_stats <- paste0(out_dir, "/stats")
out_dir_umap3d <- paste0(out_plots_umap, "/umap_3d/")
out_dir_msc <- paste0(out_plots, "/MSC/")
out_dir_immune <- paste0(out_plots, "/immune/")
out_dir_rb <- paste0(out_plots, "/ribosomal/")
out_plots_manuscript <- paste0(out_plots, "/manuscript_fig/")

#out_dir_proportion <- paste0(out_plots, "/cell_proportions/")
#out_dir_proliferation <- paste0(out_plots, "/proliferation/")

dir.create(out_dir, showWarnings = FALSE)
dir.create(out_plots, showWarnings = FALSE)
dir.create(out_initial, showWarnings = FALSE)
dir.create(out_tables, showWarnings = FALSE)
dir.create(out_plots_qc, showWarnings = FALSE)
dir.create(paste0(out_plots_qc, "/check_clusters"), showWarnings = FALSE)
dir.create(out_tables_qc, showWarnings = FALSE)
dir.create(out_plots_umap, showWarnings = FALSE)
dir.create(out_dir_stats, showWarnings = FALSE)
dir.create(out_rds, showWarnings = FALSE)
dir.create(out_plots_umap_gene, showWarnings = FALSE)
dir.create(out_plots_umap_individual, showWarnings = FALSE)
dir.create(out_dir_umap3d, showWarnings = FALSE)
dir.create(out_dir_msc, showWarnings = FALSE)
dir.create(out_dir_immune, showWarnings = FALSE)
dir.create(out_dir_rb, showWarnings = FALSE)
dir.create(out_plots_manuscript, showWarnings = FALSE)

#dir.create(out_dir_proportion, showWarnings = FALSE)
#dir.create(out_dir_proliferation, showWarnings = FALSE)


####-- Import data
## Cell annotation: sc-type
source("script/scRNAseq/_info.R")
source("script/scRNAseq/_functions.R")

## Gene annotation: ensembl 106 gene orthologue info
# Read in the full ensembl gene list from one sample in the original order and find matching gene annotation
source("script/scRNAseq/orthogene_conversion_zebrafish_human.r")
star_solo_feature_in_fp <- "annotation/gene/features.tsv"
geneid_full <- read.delim(star_solo_feature_in_fp, header=FALSE)
colnames(geneid_full) <- c("gene_id", "gene_name", "annotation")
#annotated_geneids_df <- fun.organise_gene_annotation(geneid_full)
annotated_geneids_df <- merge(geneid_full, orth_zebsalhum_df, by.x="gene_id", by.y="salmon.geneid", all.x=TRUE)
annotated_geneids_df <- annotated_geneids_df[match(geneid_full$gene_id, annotated_geneids_df$gene_id),]

# generate orthologue annotation (vs human) 
#source("script/scRNAseq/orthogene_conversion.R")


## snRNAseq data: import cellbender counts
cellbender_h5_fp <- sapply(sample_list, function(x)paste0(cellbender_dir, "/", x, "/", x, ".h5"))
cellbender_matrix_list <- sapply(cellbender_h5_fp, Read_CellBender_h5_Mat)


####-- Quality control
## QC, remove outliers and doublets
source("script/scRNAseq/qc_doublet_filter.r")


####-- Initial Merge: Merge Seurat between samples
source("script/scRNAseq/integrate_data.r")


####-- Initial merge: QC plots
source("script/scRNAseq/qc_plots_tables.r")


####-- cluster and label cell types
## Label the clusters using excel file in data folder 
source("script/scRNAseq/cluster_and_annotate.r")


####-- Draw UMAP 
source("script/scRNAseq/reduction_umap_all.r")

####-- Figure 1
source("script/scRNAseq/figure/fig1_umap_simplified.R")

####-- Figure 2
source("script/scRNAseq/figure/fig2_umap_eryth_throm.R")

####-- Figure 3
source("script/scRNAseq/figure/fig3_recruited_cells.R")

####-- Figure 5
source("script/scRNAseq/figure/fig5_myosepta.R")

####-- Figure 6
source("script/scRNAseq/figure/fig6_muscle.R")

####-- Supplementary Fig1
source("script/scRNAseq/figure/figS1_umap_full.R")

####-- Supplementary Fig2
source("script/scRNAseq/figure/figS2_dotplot_all.R")

####-- Supplementary Fig3
source("script/scRNAseq/figure/figS3_unknown_celltype.R")

####-- Supplementary Fig4
source("script/scRNAseq/figure/figS4_transitioning_cells.R")

####-- Supplementary Fig5
source("script/scRNAseq/figure/figS5_muscle_splicing.R")

####-- Table DEG
source("script/scRNAseq/table/table_deg_all_named_clusters.R")

####-- Table DEG within broad cell type groups
source("script/scRNAseq/table/table_deg_within_broadcellgroup.R")

####-- Count the number of cells in each cell type group
source("script/scRNAseq/table/table_cellcout_per_celltype.r")


####-- Save data
# Save Seurat as rds
# Save dimension reduction coordinates and metadata as csv 
source("script/scRNAseq/save_data.r")

