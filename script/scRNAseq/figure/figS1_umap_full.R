## plot by cell name (full)
# organise output folder
current_out_dir <- paste0(out_plots_manuscript, "/figS1_umap_full/")
dir.create(current_out_dir, showWarnings = FALSE)

# load colour
clusterid_to_celltype_col2 <-  readxl::read_excel(cell_cluster_col_label, sheet="full")
celltype_col2 <- clusterid_to_celltype_col2$colour
names(celltype_col2) <- clusterid_to_celltype_col2$cell_type


# plot
p <- DimPlot(merged_srat, group.by="named_celltypes", reduction = "umap.integrated.harmony", label=TRUE)
p + scale_colour_manual(values=celltype_col2) + guides(colour="none")+ggtitle("")+ xlab("UMAP_1") + ylab("UMAP_2")
ggsave(paste0(current_out_dir, "/umap_named_celltypes.pdf"), width = 7, height = 8)

# plot
p <- DimPlot(merged_srat, group.by="named_celltypes", reduction = "umap.integrated.harmony", label=FALSE)
p + scale_colour_manual(values=celltype_col2) + guides(colour="none")+ggtitle("")+ xlab("UMAP_1") + ylab("UMAP_2")
ggsave(paste0(current_out_dir, "/umap_named_celltypes_nolab.pdf"), width = 7, height = 8)


# split by treatment - vertical
p <- DimPlot(merged_srat, group.by="named_celltypes", reduction = "umap.integrated.harmony", ncol=1,split.by="treatment_group", label=FALSE, pt.size=0.05)
p + scale_colour_manual(values=celltype_col2) + guides(colour="none") +ggtitle("") + xlab("UMAP_1") + ylab("UMAP_2")
ggsave(paste0(current_out_dir, "/umap_named_celltypes_by_treat_v.pdf"), width = 4, height = 10)



# split by treatment - horizonal 
p <- DimPlot(merged_srat, group.by="named_celltypes", reduction = "umap.integrated.harmony", ncol=3,split.by="treatment_group", label=FALSE, pt.size=0.01)
p + scale_colour_manual(values=celltype_col2) + guides(colour="none") +ggtitle("")+ xlab("UMAP_1") + ylab("UMAP_2")
ggsave(paste0(current_out_dir, "/umap_named_celltypes_by_treat_h.pdf"), width = 10, height = 4)
