## 2D umap
DefaultAssay(merged_srat) <- "SCT"

## take PC scores, check them against percent.rb
merged_srat <- RunUMAP(merged_srat, dims = analysed_dims, 
                       reduction = "integrated.harmony", 
                       reduction.name = "umap.integrated.harmony",  
                       n.neighbors = 40L, seed.use = 42)

p <- DimPlot(merged_srat, group.by="final_integrated_clusters", reduction = "umap.integrated.harmony", label=TRUE)
#p <- DimPlot(merged_srat, group.by="integrated_clusters", reduction = "umap.integrated.harmony", label=TRUE)
tmp_colour <- c(clusters_colours, clusters_colours)
tmp_colour <- tmp_colour[1:length(unique(merged_srat$final_integrated_clusters))]
names(tmp_colour) <- unique(merged_srat$final_integrated_clusters)
p + scale_colour_manual(values=tmp_colour) 
ggsave(paste0(out_plots_umap, "/umap_original_cluster_id.png"), width = 8, height = 6)

p <- DimPlot(merged_srat, group.by="final_integrated_clusters", split.by="treatment_group", reduction = "umap.integrated.harmony", label=TRUE)
#p <- DimPlot(merged_srat, group.by="integrated_clusters", reduction = "umap.integrated.harmony", label=TRUE)
tmp_colour <- c(clusters_colours, clusters_colours)
tmp_colour <- tmp_colour[1:length(unique(merged_srat$final_integrated_clusters))]
names(tmp_colour) <- unique(merged_srat$final_integrated_clusters)
p + scale_colour_manual(values=tmp_colour)
ggsave(paste0(out_plots_umap, "/umap_original_cluster_id_bytreat.png"), width = 8, height = 6)


p <- DimPlot(merged_srat, group.by="final_integrated_clusters", split.by="orig.ident", reduction = "umap.integrated.harmony", label=TRUE)
#p <- DimPlot(merged_srat, group.by="integrated_clusters", reduction = "umap.integrated.harmony", label=TRUE)
tmp_colour <- c(clusters_colours, clusters_colours)
tmp_colour <- tmp_colour[1:length(unique(merged_srat$final_integrated_clusters))]
names(tmp_colour) <- unique(merged_srat$final_integrated_clusters)
p + scale_colour_manual(values=tmp_colour)
ggsave(paste0(out_plots_umap, "/umap_original_cluster_id_byorigid.png"), width = 8, height = 6)




for(feature in c("percent.rb")){
  fun.FeaturePlot(merged_srat, feature, 
                  conversion_table = annotated_geneids_df, 
                  split.by="orig.ident",
                  out_prefix = paste0(out_plots_umap, "/umap_qc_"))
}

for(feature in c("nFeature_RNA", "nFeature_RNA",
                 "percent.rb", "percent.krt", "scDblFinder.score")){
  fun.FeaturePlot(merged_srat, feature, 
                  conversion_table = annotated_geneids_df, 
                  #split.by="treatment_group",
                  out_prefix = paste0(out_plots_umap, "/umap_qc_"))
}

# 
# if(FALSE){
#   # plot named_celltypes
#   p <- DimPlot(merged_srat, group.by="named_celltypes", reduction = "umap.integrated.harmony", label=TRUE)
#   p + scale_colour_manual(values=celltype_col)
#   ggsave(paste0(out_plots_umap, "/umap_named_celltypes.pdf"), width = 15, height = 10)
#   
#   # Save 2D umap coordinates
#   umap2d_coord <- as.data.frame(merged_srat@reductions$umap.integrated.harmony@cell.embeddings)
#   write.csv(umap2d_coord, paste0(out_tables, "/umap2d_coord.csv"))
#   
# }