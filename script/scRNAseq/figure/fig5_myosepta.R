## Figure 3 - sub figures for cell recruitment/granulation tissue composition,
## particularly mesenchymal
# organise output folder
current_out_dir <- paste0(out_plots_manuscript, "/fig5_myosepta/")
dir.create(current_out_dir, showWarnings = FALSE)


for(gene in c("Scxa", "tnmd", "mkx")){
  fun.FeaturePlot(merged_srat, gene, 
                  file_type = "pdf",
                  conversion_table = annotated_geneids_df, 
                  out_prefix = current_out_dir,
                  width = 3, height = 2.4, point_size = 0.01)
}
merged_srat_fib <- subset(merged_srat, celltype_fig1 %in% c("fib1","fib2","fib3"))

Idents(merged_srat_fib) <- "celltype_fig1"
p <- DimPlot(merged_srat_fib, group.by="celltype_fig1", 
             reduction = "umap.integrated.harmony", 
             pt.size= 0.01, label=FALSE)
p <- p + scale_colour_manual(values=celltype_col)+ggtitle("")+ 
  xlab("UMAP_1") + ylab("UMAP_2")
p
ggsave(paste0(current_out_dir, "/umap_clusters.pdf"), width = 3, height = 2.4)


#### make dotplot for fibroblast/myosepta genes
srat_subset <- subset(merged_srat, celltype_fig1 %in% c("fib1","fib2","fib3", "pro2", "pro3", "pro4"))
markers_to_plot <- c(
  "ENSSSAG00000095221" = "scx",
  "ENSSSAG00000050783" = "mkx",
  "ENSSSAG00000085827" = "tnmd",
  "ENSSSAG00000047363" = "mmp2",
  "ENSSSAG00000020269" = "lum",
  "ENSSSAG00000044246" = "serpine1", 
  "ENSSSAG00000002530" = "pdgfra",
  "ENSSSAG00000073342" = "cd34",
  "ENSSSAG00000108961" = "hic1",
  "ENSSSAG00000096301" = "osr1",
  "ENSSSAG00000003963" = "pdgfrb",
  "ENSSSAG00000010277" = "rgs5"
  )
#markers_to_plot <- rev(markers_to_plot)
plot_geneid <- fun.find_matching_features(srat_subset, names(markers_to_plot), 
                                          conversion_table = annotated_geneids_df,
                                          remove.dash=FALSE)
Idents(srat_subset) <- "named_celltypes"
#current_celltypes <- unique(srat_subset$named_celltypes)
celltype_order <- c(
  "fibroblast, tenocyte" , "fibroblast", "fibroblast, wound-associated",  
  "progenitor, MSC", "progenitor, fibroadipogenic", 
  "progenitor, pericyte")
srat_subset$named_celltypes <- as.character(srat_subset$named_celltypes)
srat_subset$named_celltypes <- factor(srat_subset$named_celltypes, levels = rev(celltype_order))
Idents(srat_subset) <- srat_subset$named_celltypes
p <- DotPlot(object = srat_subset, features = names(markers_to_plot),  scale.max=50)
p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p <- p + scale_x_discrete(labels= markers_to_plot) + ylab("") + xlab("")
p + scale_colour_viridis(direction = -1)
ggsave(paste0(current_out_dir, "/dotplot.pdf"), width=7, height=3)

