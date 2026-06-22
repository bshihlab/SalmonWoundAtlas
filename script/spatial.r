######## Spatial transcriptomics
####-- Load libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library("terra")
library(fs)
library('glmGamPoi')
library(cowplot)
library(grid)


if(FALSE){
  BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                         'limma', 'S4Vectors', 'SingleCellExperiment',
                         'SummarizedExperiment', 'batchelor', 'Matrix.utils', 'EBImage'))
  
  # packages from github
  devtools::install_github(repo = "kueckelj/confuns")
  devtools::install_github(repo = "theMILOlab/SPATAData")
  devtools::install_github(repo = "theMILOlab/SPATA2")
  devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
  
}


####-- Settings
working_dir <- "L:/ShihLab/projects/2024/20240226_rose_wounds/"
#gene_list_stats <- "analysis/_archive/ribosomal_cluster_keep_20250328/stats/all_markers_finalclustering_unique.csv"
##gene_list_stats <- "analysis/ribosomal_cluster_keep/stats/all_markers_finalclustering_unique.csv"
##gene_list_stats <- "analysis/ribosomal_cluster_keep/stats/all_markers_finalclustering_unique.csv"
#gene_list_stats <- "analysis/ribosomal_cluster_keep/stats/all_markers_finalclustering_unique_top20.csv"
gene_list_proliferation <- "data/cell_type_identification/proliferation_doi-10.1371-journal.pcbi.1010604_tables4.csv"
sample_interest <- c("analysis/spaceranger/1-SI-TT-A3", "analysis/spaceranger/2-SI-TT-B3" )
snRNAseq_rds <- "analysis/ribosomal_cluster_keep/rds/merged_srat.rds"

out_dir <- "analysis/spatial"
out_dir_rds <- "analysis/spatial/rds/"
out_dir_plots <- "analysis/spatial/plots/"
out_dir_avgcelltype <- "analysis/spatial/plots/averaged_celltype/"
out_dir_individualGene <- "analysis/spatial/plots/individual_gene/"
out_dir_phase <- "analysis/spatial/plots/cellcycle/"
in_dir_10x <- "analysis/spaceranger"
setwd(working_dir)
source("script/scRNAseq/_functions.R")

dir.create(dirname(out_dir), showWarnings = FALSE)
dir.create(out_dir, showWarnings = FALSE)
dir.create(out_dir_plots, showWarnings = FALSE)
dir.create(out_dir_rds, showWarnings = FALSE)
dir.create(out_dir_avgcelltype, showWarnings = FALSE)
dir.create(out_dir_individualGene, showWarnings = FALSE)
dir.create(out_dir_phase, showWarnings = FALSE)

image_rot_angle <- list(
  "1-SI-TT-A3" = 240,
  "1-SI-TT-A3_masked" = 240,
  "2-SI-TT-B3" = 215,
  "3-SI-TT-C3" = 0,
  "4-SI-TT-D3" = 90
  
)

pt_size_factor <- list(
  "1-SI-TT-A3" = 2.5,
  "1-SI-TT-A3_masked" = 2.5,
  "2-SI-TT-B3" = 2.5,
  "3-SI-TT-C3" = 2.5,
  "4-SI-TT-D3" = 2.5
)


####-- Create a annotated gene dataframe
# Read in the full ensembl gene list in the original order, repalce the gene ID with this list
source("script/scRNAseq/orthogene_conversion_zebrafish_human.r")
geneid_full <- read.delim(gzfile("analysis/spaceranger/1-SI-TT-A3/outs/raw_feature_bc_matrix/features.tsv.gz"), header=FALSE)
colnames(geneid_full) <- c("gene_id", "gene_name", "annotation")
#annotated_geneids_df <- fun.organise_gene_annotation(geneid_full)
annotated_geneids_df <- merge(geneid_full, orth_zebsalhum_df, by.x="gene_id", by.y="salmon.geneid", all.x=TRUE)
annotated_geneids_df <- annotated_geneids_df[match(geneid_full$gene_id, annotated_geneids_df$gene_id),]


## other gene lists
prolifeation_genes <- read.csv(gene_list_proliferation)
ribosomal_genes <- annotated_geneids_df[(grepl(pattern = "^RPS", x =  annotated_geneids_df$gene_name, ignore.case=TRUE)) | (grepl(pattern = "^RPS", x =  annotated_geneids_df$gene_name_base, ignore.case=TRUE)), ]
ribosomal_genes <- ribosomal_genes[ribosomal_genes$gene_id %in% geneid_full$gene_id,]
other_gene_lists <- list(Proliferation = prolifeation_genes$gene, 
                         Ribosome = ribosomal_genes$gene_name_base) 


#### Import gene lists for each cell type
# snRNAseq cluster specific genes and plot them
##cluster_specific_genes <- read.csv(gene_list_stats)
##cluster_specific_genes$named_cluster_pathsafe <- gsub(", ", "-", path_sanitize(cluster_specific_genes$cluster))



####-- Cell cycle scoring and rb genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes_salmon <- annotated_geneids_df$gene_id[grep(paste0(s.genes, collapse = "|"), annotated_geneids_df$gene_name_base , ignore.case=TRUE)]
g2m.genes_salmon <- annotated_geneids_df$gene_id[grep(paste0(s.genes, collapse = "|"), annotated_geneids_df$gene_name_base , ignore.case=TRUE)]
RB_feature <- annotated_geneids_df$gene_id[grep("^RP[SL]", annotated_geneids_df$gene.name , ignore.case=TRUE)]
#RB_feature <- c(RB_feature, annotated_geneids_df$gene_id[grep("^RP[SL]", annotated_geneids_df$gene_name  , ignore.case=TRUE)])
RB_feature <- c(RB_feature, annotated_geneids_df$gene_id[grep("^RP[SL]", annotated_geneids_df$gene_name_base , ignore.case=TRUE)])
RB_feature <- unique(RB_feature)

	
####-- Read in all spatial data
srat_spatial <- list()
sample_list <- list.dirs(in_dir_10x, recursive=FALSE)
for(sample_dir in sample_list){
#for(sample_dir in sample_interest){
  current_sample_name <- gsub(in_dir_10x, "", sample_dir)
  current_sample_name <- gsub("\\/", "", current_sample_name)
  current_image_rot_angle = image_rot_angle[[current_sample_name]]
  current_pt.size.factor <- pt_size_factor[[current_sample_name]]
  #current_wound_spots <- read.csv( paste0(sample_dir, "/outs/keep_spots.csv") ) # this is manually drawn around the wounds using Loup browser
  current_deconvolution <- read.csv(paste0("analysis/cell2location/deconvoluted_abundance/", current_sample_name, ".csv"), header=FALSE)
  annotation_names <- current_deconvolution[1,2:ncol(current_deconvolution)]
  spot_ids <- current_deconvolution[2:nrow(current_deconvolution),1]
  current_deconvolution <- current_deconvolution[2:nrow(current_deconvolution),2:ncol(current_deconvolution)]
  row.names(current_deconvolution) <- sapply(strsplit( spot_ids, "_") ,FUN=function(x)x[2])
  
  # create a vector where the annotation name is the vector, and
  # the name is for the file name/srat column name
  all_annotations <- annotation_names
  
  # remove neuron from annotation so the ordering would stay the same as before
#  all_annotations <- all_annotations[!all_annotations %in% "neuron"]
#  neuron_filename <- c("A0neuron" = "neuron")
  
  #clean file names
  all_annotations_clean <- gsub("[[:punct:]]", "_", annotation_names)
  names(all_annotations) <- paste0("A", 1:length(all_annotations), all_annotations_clean)
  
  # add neuron files
#  all_annotations <- c(neuron_filename, all_annotations)
  
  
  ## Read in high res image
  img = Read10X_Image(paste0(sample_dir, "/outs/spatial/"), image.name = "tissue_hires_image.png")
  current_srat <- Load10X_Spatial(paste0(sample_dir, "/outs"), image = img)
  #current_srat@images$slice1@scale.factors$lowres = current_srat@images$slice1@scale.factors$hires
  print(dim(current_srat))
  
  ## Replace row names with ensembl id
  geneid_full <- read.delim(gzfile(paste0(sample_dir, "/outs/filtered_feature_bc_matrix/features.tsv.gz")), header=FALSE)
  colnames(geneid_full) <- c("gene_id", "gene_name", "annotation")
  row.names(current_srat) <- geneid_full$gene_id
  
  ## Remove empty spots
  current_srat$WoundSpot <- ifelse(colnames(current_srat) %in% current_deconvolution$spot_id, "wound", "non-wound")
  current_srat<- current_srat[,(colnames(current_srat) %in% row.names(current_deconvolution))]
  
  # reorder to make sure deconvoluted spot order matches srat
  current_deconvolution <- current_deconvolution[match(colnames(current_srat), row.names(current_deconvolution)),]
  for(idx in 1:length(all_annotations)){
    current_annotation_filename <- names(all_annotations)[idx]
    current_annotation_figtitle <- all_annotations[idx]
    # add annotation
    current_srat[[current_annotation_filename]] <- as.numeric(current_deconvolution[,idx])

    # create output directories
    current_out_prefix_individual <- paste0(out_dir, "/plots/cell2location_deconvoluted/") 
    dir.create(current_out_prefix_individual, showWarnings = FALSE)
    
    # pdf
    fun.spatialPlot_individual(obj = current_srat, 
                               pt.size.factor = current_pt.size.factor,
                               features = current_annotation_filename, 
                               conversion_table = annotated_geneids_df,
                               out_prefix = current_out_prefix_individual, 
                               out_suffix = current_sample_name,
                               remove.dash=TRUE,
                               max.cutoff = "q99",
                               min.cutoff = NA,
                               file_type = "pdf",
                               legend_title = current_annotation_figtitle,
                               rot_angle = current_image_rot_angle)
    # png
    fun.spatialPlot_individual(obj = current_srat, 
                               pt.size.factor = current_pt.size.factor,
                               features = current_annotation_filename, 
                               conversion_table = annotated_geneids_df,
                               out_prefix = current_out_prefix_individual, 
                               out_suffix = current_sample_name,
                               remove.dash=TRUE,
                               max.cutoff = "q99",
                               min.cutoff = NA,
                               file_type = "png",
                               legend_title = current_annotation_figtitle,
                               rot_angle = current_image_rot_angle)
    
  }
  
  # remove spot with low count
  #current_srat <- current_srat[,unname(which(colSums(current_srat@assays$Spatial$counts)>=1000))]

  # Save this filtered file as a separate object
  current_srat_versions = list()
  exclude_barcode_fp <- paste0("data/spatial_exclude_barcodes/", basename(sample_dir), ".csv")
  if(current_sample_name %in% c("1-SI-TT-A3")){
    current_srat_versions[[paste0(current_sample_name, "_masked")]] = current_srat
  } else {
    current_srat_versions[[current_sample_name]] = current_srat
  }
  
  
  for(store_sample_id in names(current_srat_versions)){
    current_srat <- current_srat_versions[[store_sample_id]]
    current_srat[["percent.rb"]] <- PercentageFeatureSet(current_srat,  features = RB_feature)
    
    current_srat <- NormalizeData(current_srat, verbose = FALSE, assay = "Spatial")
    current_srat <- FindVariableFeatures(current_srat, verbose = FALSE)
    current_srat <- ScaleData(current_srat)
    current_srat <- CellCycleScoring(current_srat, s.features = s.genes_salmon, g2m.features = g2m.genes_salmon, set.ident = TRUE)
    
    # rerun normalisation and dimension reduction
    current_srat <- SCTransform(current_srat, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
    # also run standard log normalization for comparison
    current_srat <- NormalizeData(current_srat, verbose = FALSE, assay = "Spatial")
    current_srat <- RunPCA(current_srat, assay = "SCT", verbose = FALSE)
    current_srat <- FindNeighbors(current_srat, reduction = "pca", dims = 1:30)
    current_srat <- FindClusters(current_srat, verbose = FALSE)
    current_srat <- RunUMAP(current_srat, reduction = "pca", dims = 1:30)
    
    
    # check overall plot
    plot1 <- VlnPlot(current_srat, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
    plot2 <- SpatialFeaturePlot(current_srat, features = "nCount_Spatial", image.scale="hires") + theme(legend.position = "right") + ggtitle(store_sample_id)
    plot2 <- rotate_image(plot2, rot_angle = current_image_rot_angle)
    
    print(wrap_plots(plot1, plot2))

    #p1 <- DimPlot(current_srat, reduction = "umap", label = TRUE)
    #p2 <- SpatialDimPlot(current_srat, label = TRUE, label.size = 3, ncol=1, image.scale="hires") +ggtitle(store_sample_id)
    #p2 <- rotate_image(plot2, rot_angle = current_image_rot_angle)
    
    #print(p1 + p2)
    
    ## Find spatially variable features
    #DefaultAssay(current_srat) <- "SCT"
    #current_srat <- FindSpatiallyVariableFeatures(current_srat, assay = "SCT", layer = "scale.data", features = VariableFeatures(current_srat)[1:1000],selection.method = "moransi", x.cuts = 100, y.cuts = 100)
    #SpatialFeaturePlot(current_srat, features = head(SpatiallyVariableFeatures(current_srat, selection.method = "moransi"), 6), ncol = 3, alpha = c(0.1, 1), max.cutoff = "q95")                              
    
    ## Add scores based on other gene list
    DefaultAssay(current_srat) <- "SCT"
    for(genelist_name in names(other_gene_lists)){
      current_gene_list <- other_gene_lists[[genelist_name]]
      matched_features <- fun.find_matching_features(current_srat, current_gene_list, annotated_geneids_df, remove.dash=FALSE)
      all_matching_geneid <- matched_features[["geneid"]]
      all_matching_outname <- matched_features[["outname"]]
      current_srat <- AddModuleScore(current_srat, features=list(all_matching_geneid),name = genelist_name)
    }
    
  
    
    # store normalised spatial seurat in a list
    srat_spatial[[store_sample_id]] <- current_srat
    out_rds_fp <- paste0(out_dir_rds, store_sample_id, ".rds")
    saveRDS(current_srat, out_rds_fp)
  }
}



####-- Plot specific genes
gene_interest_specific <- c("ENSSSAG00000065665" = "rps25", 
                            "ENSSSAG00000073156" = "rpl8",
                            "ENSSSAG00000066053" = "rpl30"
                            )
if(TRUE){
  ## Loop through genes
  for(current_sample in names(srat_spatial)){
    current_srat <- srat_spatial[[current_sample]]
    current_image_rot_angle <- image_rot_angle[[current_sample]]
    current_pt.size.factor <- pt_size_factor[[current_sample]]
    for(current_gene in names(gene_interest_specific)){
      # create output directories
      current_out_prefix_individual <- paste0(out_dir_individualGene, "/specific_genes/") 
      dir.create(current_out_prefix_individual, showWarnings = FALSE)
      print(current_gene)
      fun.spatialPlot_individual(obj = current_srat, 
                                 features = current_gene, 
                                 conversion_table = annotated_geneids_df,
                                 pt.size.factor = current_pt.size.factor,
                                 out_prefix = current_out_prefix_individual, 
                                 out_suffix = current_sample,
                                 remove.dash=TRUE,
                                 rot_angle = current_image_rot_angle)
    }
  }
}
# 
# ####-- Plot average gene expression for groups of genes
# ## Get the average expression for groups of markers
# # markers from scRNAseq clustering
# genes_interest <- split(cluster_specific_genes$gene, cluster_specific_genes$named_cluster_pathsafe)
# 
# # markers for proliferation
# genes_interest <- c(genes_interest, other_gene_lists)
# 
# ##  loop through to plot each gene group
# for(current_sample in names(srat_spatial)){
#   current_srat <- srat_spatial[[current_sample]]
#   current_image_rot_angle = image_rot_angle[[current_sample]]
#   for(celltype in names(genes_interest)){
#       # print plots averaged across all markers
#       current_genes_interest <- genes_interest[[celltype]]
# 
#       # create output directories
#       current_out_prefix_avg <- paste0(out_dir_avgcelltype, "/", celltype)
#       current_out_prefix_individual <- paste0(out_dir_individualGene, "/", celltype, "/")
#       dir.create(current_out_prefix_individual, showWarnings = FALSE)
# 
#       # plot graphs
#       fun.spatialPlot_avg(obj = current_srat,
#                           features = current_genes_interest,
#                           conversion_table = annotated_geneids_df,
#                           out_prefix = current_out_prefix_avg,
#                           out_suffix = current_sample,
#                           legend_label=celltype,
#                           remove.dash=TRUE,
#                           file_type = "png",
#                           exprs_fill = FALSE,
#                           rot_angle = current_image_rot_angle)
#       fun.spatialPlot_avg(obj = current_srat,
#                           features = current_genes_interest,
#                           conversion_table = annotated_geneids_df,
#                           out_prefix = current_out_prefix_avg,
#                           out_suffix = current_sample,
#                           legend_label=celltype,
#                           remove.dash=TRUE,
#                           exprs_fill = FALSE,
#                           file_type = "pdf",
#                           rot_angle = current_image_rot_angle)
#       if(FALSE){
#       current_srat@images$slice1$centroids@coords %>%
#         bind_cols(current_srat@meta.data) %>%
#         mutate(x=x,y=y) %>%
#         ggplot(aes(x=y,y=x)) +
#         geom_point(aes(col=nCount_Spatial), shape=21) +
#         scale_colour_viridis_c(direction=-1) +
#         scale_y_reverse() +
#         coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
#         theme_bw()
#       }
# 
#       # print individual plots for each gene
#       if(TRUE){
#         fun.spatialPlot_individual(obj = current_srat,
#                                  features = current_genes_interest,
#                                  conversion_table = annotated_geneids_df,
#                                  out_prefix = current_out_prefix_individual,
#                                  out_suffix = current_sample,
#                                  remove.dash=TRUE,
#                                  rot_angle = current_image_rot_angle)
# 
#       }
#     }
# }
# 
# # dual figure
# ##  loop through to plot each gene group
# current_srat1 <- srat_spatial[["1-SI-TT-A3_masked"]]
# current_srat2 <- srat_spatial[["2-SI-TT-B3"]]
# current_image_rot_angle1 = image_rot_angle[["1-SI-TT-A3_masked"]]
# current_image_rot_angle2 = image_rot_angle[["2-SI-TT-B3"]]
# for(celltype in names(genes_interest)){
#   # print plots averaged across all markers
#   current_genes_interest <- genes_interest[[celltype]]
# 
#   # create output directories
#   current_out_prefix_avg <- paste0(out_dir_avgcelltype, "/dualplot_", celltype)
#   dir.create(current_out_prefix_individual, showWarnings = FALSE)
# 
#   # plot graphs
#   p <- fun.spatialPlot_dual(obj1 = current_srat1,
#                              obj2 = current_srat2,
#                       features = current_genes_interest,
#                       conversion_table = annotated_geneids_df,
#                       out_prefix = current_out_prefix_avg,
#                       out_suffix = current_sample,
#                       legend_label=celltype,
#                       remove.dash=TRUE,
#                       file_type = "png",
#                       exprs_fill = FALSE,
#                       rot_angle1 = current_image_rot_angle1,
#                       rot_angle2 = current_image_rot_angle2)
#   # plot graphs
#   p <- fun.spatialPlot_dual(obj1 = current_srat1,
#                             obj2 = current_srat2,
#                             features = current_genes_interest,
#                             conversion_table = annotated_geneids_df,
#                             out_prefix = current_out_prefix_avg,
#                             out_suffix = current_sample,
#                             legend_label=celltype,
#                             remove.dash=TRUE,
#                             file_type = "pdf",
#                             exprs_fill = FALSE,
#                             rot_angle1 = current_image_rot_angle1,
#                             rot_angle2 = current_image_rot_angle2)
# }
# 
# 
# 
# 
# # Add in spot deconvolution
# snRNAseq <- readRDS(snRNAseq_rds)
# DefaultAssay(snRNAseq) <- "SCT"
# 
# for(current_sample in names(srat_spatial)){
#   current_srat <- srat_spatial[[current_sample]]
#   current_image_rot_angle = image_rot_angle[[current_sample]]
#   current_srat <- SCTransform(current_srat, assay = "Spatial", verbose = FALSE) %>% RunPCA(verbose = FALSE)
#   anchors <- FindTransferAnchors(reference = snRNAseq, query = current_srat, normalization.method = "SCT")
#   predictions.assay <- TransferData(anchorset = anchors, refdata = snRNAseq$named_celltypes, prediction.assay = TRUE,
#                                     weight.reduction = current_srat[["pca"]], dims = 1:30)
#   current_srat[["predictions"]] <- predictions.assay
#   for(cell_type in unique(as.character(snRNAseq$named_celltypes))){
# 
#     # create output directories
#     current_out_prefix_avg <- paste0(out_dir_deconvolutedcelltype, "/", cell_type)
# 
#     # plot and save
#     fun.spatialPlot_description(obj = current_srat,
#                                features = cell_type,
#                                out_prefix = current_out_prefix_avg,
#                                out_suffix = current_sample,
#                                remove.dash=TRUE,
#                                rot_angle = current_image_rot_angle)
#   }
# }
# 
# 
# 
# ## Plot ribosomal score vs proliferation score for spots within the wound
# testribo <- current_srat$Ribosome1[current_srat$WoundSpot == "wound"]
# testpro <- current_srat$Proliferation1[current_srat$WoundSpot == "wound"]
# plot(testribo~testpro)
# SpatialFeaturePlot(current_srat, "Ribosome1",image.scale="hires",pt.size.factor = 2.2,max.cutoff = "q95")
# SpatialFeaturePlot(current_srat, "Proliferation1",image.scale="hires",pt.size.factor = 2.2,max.cutoff = "q95")
