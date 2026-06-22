### Figure S4
# organise output folder
current_out_dir <- paste0(out_plots_manuscript, "/figS4_transitioning_cell/")
dir.create(current_out_dir, showWarnings = FALSE)






### annotate clusters near transiting cells
annotation_transiting_cells <- read.csv("annotation/cluster/_cluster_id_to_celltype_transitioning_cells.csv")
clusterid_to_celltype_col2 <-  readxl::read_excel(cell_cluster_col_label, sheet="full")
celltype_col2 <- clusterid_to_celltype_col2$colour
names(celltype_col2) <- clusterid_to_celltype_col2$cell_type

# use the same annotation colour as figs1
celltype_col_T <- celltype_col2[names(celltype_col2) %in% annotation_transiting_cells$named_celltype_T]
celltype_col_T <- c(celltype_col_T, c("other" = "#D3D3D3"))

### plot the umap before removing high-ribosome transitioning cells
merged_srat$named_celltypes <- as.character(merged_srat$named_celltypes)
Idents(merged_srat) <- "named_celltypes"
merged_srat$named_celltype_T <- ifelse(merged_srat$named_celltypes %in% annotation_transiting_cells$named_celltype_T, 
                                       merged_srat$named_celltypes, "other")

p <- DimPlot(merged_srat, group.by="named_celltype_T", 
             reduction = "umap.integrated.harmony", label=TRUE)
p + scale_colour_manual(values=celltype_col_T) + 
  guides(colour="none")+ggtitle("")+ xlab("UMAP_1") + ylab("UMAP_2")
ggsave(paste0(current_out_dir, "/umap_named_celltypes.pdf"), width = 4, height = 5)

p <- DimPlot(merged_srat, group.by="named_celltype_T", 
             split.by="treatment_group", pt.size=0.01,
             reduction = "umap.integrated.harmony", label=FALSE)
p + scale_colour_manual(values=celltype_col_T) + 
  guides(colour="none")+ggtitle("")+ xlab("UMAP_1") + ylab("UMAP_2")
ggsave(paste0(current_out_dir, "/umap_named_celltypes_split_treat.pdf"), width = 8, height = 4.5)

### Plot ribosome content
## Other qc
for(gene in c("nFeature_RNA", "Rb.Score1", "ENSSSAG00000048424", "col1a2", "gapdh")){
  fun.FeaturePlot(merged_srat, gene, 
                  file_type = "pdf",
                  conversion_table = annotated_geneids_df, 
                  out_prefix = current_out_dir,
                  width = 6, height = 5, point_size = 0.01)
}


## work out the median nFeature_RNA

VlnPlot(merged_srat_backup, "nFeature_RNA") +  guides(colour="none")
ggsave(paste0(current_out_dir, "/vlnplot_nFeatureRNA.pdf"), width=15, height = 5)
