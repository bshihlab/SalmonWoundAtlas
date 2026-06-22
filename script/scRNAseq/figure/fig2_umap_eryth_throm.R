### UMAP protion of Figure S2
# organise output folder
current_out_dir <- paste0(out_plots_manuscript, "/fig2_umap_eryth_throm/")
dir.create(current_out_dir, showWarnings = FALSE)

# load colour
clusterid_to_celltype_col2 <-  readxl::read_excel(cell_cluster_col_label, sheet="full")
celltype_col2 <- clusterid_to_celltype_col2$colour
names(celltype_col2) <- clusterid_to_celltype_col2$cell_type


# subset data
srat_subset <- subset(merged_srat, celltype_fig1 %in% c("osteoclast", "erythrocyte","macrophage", "leukocyte", "thrombocyte"))

# organise genes
markers_to_plot <- c(
                     "ENSSSAG00000100947" = "spp1", 
                     "ENSSSAG00000047020" = "csf1ra",
                     "ENSSSAG00000043023" = "jak3",
                     "ENSSSAG00000005686" = "jak1",
                     "ENSSSAG00000009589" = "vav2",
                     "ENSSSAG00000009390" = "flt3",
                     "ENSSSAG00000072535" = "csf3r", 
                     "ENSSSAG00000093862" = "hbae1",
                     "ENSSSAG00000046009" = "hbae4",
                     "ENSSSAG00000072097" = "G6fl")
markers_to_plot <- rev(markers_to_plot)
plot_geneid <- fun.find_matching_features(srat_subset, names(markers_to_plot), 
                                          conversion_table = annotated_geneids_df,
                                          remove.dash=FALSE)
Idents(srat_subset) <- "named_celltypes"
#current_celltypes <- unique(srat_subset$named_celltypes)
celltype_order <- c(
  "osteoclast" , "macrophage", "lymphocyte, T",  "lymphocyte, B", 
  "dendritic cell","neutrophil",
  "erythrocyte", "thrombocyte")
srat_subset$named_celltypes <- as.character(srat_subset$named_celltypes)
srat_subset$named_celltypes <- factor(srat_subset$named_celltypes, levels = celltype_order)
Idents(srat_subset) <- srat_subset$named_celltypes
p <- DotPlot(object = srat_subset, features = names(markers_to_plot),  scale.max=50)
p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p <- p + scale_x_discrete(labels= markers_to_plot) + ylab("") + xlab("")
p + scale_colour_viridis(direction = -1)
ggsave(paste0(current_out_dir, "/dotplot.pdf"), width=8.5, height=2.8)


# print out umap for each gene
for(gene in names(markers_to_plot)){
  fun.FeaturePlot(merged_srat, gene, 
                  conversion_table = annotated_geneids_df, 
                  out_prefix = paste0(current_out_dir, "/featureplot_"),
                  width=3, height=2.5, file_type = "pdf")
}



# # write csv on significat DEG for inflammation erythrocyte
# deg_eryth_inflam <- all_cluster_deg_named[all_cluster_deg_named$cluster %in% "erythrocyte, inflammatory",]
# deg_eryth_inflam <- deg_eryth_inflam[deg_eryth_inflam$p_val_adj < 0.0001,]
# deg_eryth_inflam <- deg_eryth_inflam[deg_eryth_inflam$avg_log2FC > 2,]
# deg_eryth_inflam <- merge(deg_eryth_inflam, annotated_geneids_df, by.x="gene", by.y="gene_id")
# deg_eryth_inflam <- data.frame(cluster="erythrocyte, inflammatory", gene=unique(deg_eryth_inflam$human.genename))
# write.csv(deg_eryth_inflam, "analysis/ribosomal_cluster_keep/cellcellcommunication/DEG_erythrocyte_inflam.csv", row.names=FALSE)