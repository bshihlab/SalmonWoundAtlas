### UMAP protion of Figure S1
# organise output folder
current_out_dir <- paste0(out_plots_manuscript, "/fig1_umap_simplified/")
dir.create(current_out_dir, showWarnings = FALSE)

# load colour
clusterid_to_celltype_col <-  readxl::read_excel(cell_cluster_col_label, sheet="fig1")
celltype_col <- clusterid_to_celltype_col$colour
names(celltype_col) <- clusterid_to_celltype_col$cell_type


## plot by cell name (fig1)
p <- DimPlot(merged_srat, group.by="celltype_fig1", reduction = "umap.integrated.harmony", label=TRUE)
p + scale_colour_manual(values=celltype_col)+ggtitle("")+ xlab("UMAP_1") + ylab("UMAP_2")
ggsave(paste0(current_out_dir, "/umap_celltype_fig1.pdf"), width = 9, height = 6)

# split by treatment - vertical
p <- DimPlot(merged_srat, group.by="celltype_fig1", reduction = "umap.integrated.harmony", ncol=1,split.by="treatment_group", label=FALSE, pt.size=0.05)
p + scale_colour_manual(values=celltype_col) + guides(colour="none") +ggtitle("")+ xlab("UMAP_1") + ylab("UMAP_2")
ggsave(paste0(current_out_dir, "/umap_celltype_fig1_bytreat_v.pdf"), width = 4, height = 10)

# split by treatment - horizonal 
p <- DimPlot(merged_srat, group.by="celltype_fig1", reduction = "umap.integrated.harmony", ncol=3,split.by="treatment_group", label=FALSE, pt.size=0.05)
p + scale_colour_manual(values=celltype_col) + guides(colour="none") +ggtitle("")+ xlab("UMAP_1") + ylab("UMAP_2")
ggsave(paste0(current_out_dir, "/umap_celltype_fig1_bytreat_h.pdf"), width = 10, height = 4)


