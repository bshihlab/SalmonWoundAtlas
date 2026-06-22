## Define clusters
DefaultAssay(merged_srat) <- "SCT"
merged_srat <- RunPCA(merged_srat, npcs=100)
set.seed(123)

## Elbow plot
ElbowPlot(merged_srat, ndims = 100, reduction = "pca")

# remove ribosomal cluster
# find clusters
analysed_dims <- 1:52
merged_srat <- FindNeighbors(merged_srat, reduction = "integrated.harmony",  
                             dims = analysed_dims, graph.name="clustering1")
merged_srat <- FindClusters(merged_srat, resolution = 1.1, 
                            cluster.name = "integrated_clusters", 
                            graph.name="clustering1")
merged_srat$initial_integrated_clusters <- merged_srat$integrated_clusters
merged_srat$final_integrated_clusters <- merged_srat$integrated_clusters

Idents(merged_srat) <- "integrated_clusters"
merged_srat <- FindSubCluster(merged_srat, "7", graph.name = "clustering1", resolution = 0.2)
merged_srat$final_integrated_clusters <- ifelse(as.character(merged_srat$integrated_clusters) == "7", merged_srat$sub.cluster, as.character(merged_srat$integrated_clusters))
merged_srat <- FindSubCluster(merged_srat, "13", graph.name = "clustering1", resolution = 0.1)
merged_srat$final_integrated_clusters <- ifelse(as.character(merged_srat$integrated_clusters) == "13", merged_srat$sub.cluster, as.character(merged_srat$final_integrated_clusters))
merged_srat <- FindSubCluster(merged_srat, "21", graph.name = "clustering1", resolution = 0.1)
merged_srat$final_integrated_clusters <- ifelse(as.character(merged_srat$integrated_clusters) == "21", merged_srat$sub.cluster, as.character(merged_srat$final_integrated_clusters))
merged_srat <- FindSubCluster(merged_srat, "18", graph.name = "clustering1", resolution = 0.1)
merged_srat$final_integrated_clusters <- ifelse(as.character(merged_srat$integrated_clusters) == "18", merged_srat$sub.cluster, as.character(merged_srat$final_integrated_clusters))
merged_srat <- FindSubCluster(merged_srat, "15", graph.name = "clustering1", resolution = 0.1)
merged_srat$final_integrated_clusters <- ifelse(as.character(merged_srat$integrated_clusters) == "15", merged_srat$sub.cluster, as.character(merged_srat$final_integrated_clusters))
merged_srat <- FindSubCluster(merged_srat, "17", graph.name = "clustering1", resolution = 0.1)
merged_srat$final_integrated_clusters <- ifelse(as.character(merged_srat$integrated_clusters) == "17", merged_srat$sub.cluster, as.character(merged_srat$final_integrated_clusters))
merged_srat <- FindSubCluster(merged_srat, "38", graph.name = "clustering1", resolution = 0.1)
merged_srat$final_integrated_clusters <- ifelse(as.character(merged_srat$integrated_clusters) == "38", merged_srat$sub.cluster, as.character(merged_srat$final_integrated_clusters))



Idents(merged_srat) <- "final_integrated_clusters"
#merged_srat <- FindSubCluster(merged_srat, "11_0", graph.name = "clustering1", resolution = 0.3)
#merged_srat$final_integrated_clusters <- ifelse(as.character(merged_srat$final_integrated_clusters) == "11_0", merged_srat$sub.cluster, as.character(merged_srat$final_integrated_clusters))

#### Organise files for cell2location
# write out the RNA count matrix and cell type label
cluster_annotation <- read.csv(cell_cluster_label)

# celltype for STdeconvolution
Idents(merged_srat) <- "final_integrated_clusters"
new_cluster_id <- cluster_annotation$celltype_Stdeconvolution
names(new_cluster_id) <- cluster_annotation$cluster_id
merged_srat <- RenameIdents(merged_srat, new_cluster_id)
merged_srat$celltype_Stdeconvolution <- merged_srat@active.ident

# celltype for plotting figure 1
Idents(merged_srat) <- "final_integrated_clusters"
new_cluster_id <- cluster_annotation$celltype_fig1
names(new_cluster_id) <- cluster_annotation$cluster_id
merged_srat <- RenameIdents(merged_srat, new_cluster_id)
merged_srat$celltype_fig1 <- merged_srat@active.ident

# names cell types
Idents(merged_srat) <- "final_integrated_clusters"
new_cluster_id <- cluster_annotation$celltype_full
names(new_cluster_id) <- cluster_annotation$cluster_id
merged_srat <- RenameIdents(merged_srat, new_cluster_id)
merged_srat$named_celltypes <- merged_srat@active.ident
write.csv(merged_srat@meta.data, "analysis/cell2location/metadata.csv")


# 
# #### find all markers (based on clusters)
# merged_srat<- PrepSCTFindMarkers(merged_srat)
# Idents(merged_srat) <- "final_integrated_clusters"
# all_cluster_deg <- FindAllMarkers(merged_srat, only.pos=TRUE)
# #all_cluster_deg_sig <- all_cluster_deg_sig[all_cluster_deg_sig$adj.p.]
# out_stats <- merge(all_cluster_deg, annotated_geneids_df, by.x="gene", by.y="gene_id", all.x=TRUE)
# out_stats_unique <- out_stats[order(out_stats$p_val_adj),] 
# out_stats_unique <- out_stats_unique[!duplicated(out_stats_unique$gene),]
# write.csv(out_stats_unique, paste0(out_tables_qc,"/cluster_characterisation_deg.csv"))
# 
# #### find "unique" markers
# 
# # Find markers that are unique to each cluster
# all_markers_df_sig <- out_stats[out_stats$p_val_adj < 0.0001,] # keep genes with low p value
# all_markers_df_sig <- all_markers_df_sig[all_markers_df_sig$pct.1 > 0.10,] # keep genes with at least 15% expressed in pct.1
# all_markers_df_sig <- all_markers_df_sig[all_markers_df_sig$avg_log2FC  > 1.25,] # keep genes with at least 1.25 log2FC
# all_markers_list <- split(all_markers_df_sig, all_markers_df_sig$cluster)
# 
# # keep any genes that are also present in other gene lists
# all_markers_list_unique <- lapply(names(all_markers_list), function(x){
#   out <- all_markers_list[[x]]
#   out <- out[!(out$gene %in% all_markers_df_sig$gene[!all_markers_df_sig$cluster == x]),]
#   return(out)
# }
# )
# 
# # keep genes that are in multiple cell types only 
# # if the highest differentially expressed cluster has at least 1 fold 
# # higher expression than the second cluster
# all_markers_list_repeat <- lapply(names(all_markers_list), function(x)
#   all_markers_list[[x]][all_markers_list[[x]]$gene %in% all_markers_df_sig$gene[!all_markers_df_sig$cluster == x],] 
# )
# all_markers_list_repeat <- do.call(rbind, all_markers_list_repeat)
# all_markers_list_repeat <- all_markers_list_repeat[order(all_markers_list_repeat$avg_log2FC, -all_markers_list_repeat$avg_log2FC),]
# all_markers_list_repeat <- split(all_markers_list_repeat, all_markers_list_repeat$gene)
# all_markers_list_repeat <- lapply(all_markers_list_repeat, function(x){if(x$avg_log2FC[1] - x$avg_log2FC[2] > 1){out <- x[1,]}else{out <- NULL} ; return(out)})
# 
# # combine the two lists
# all_markers_list_unique <- rbind(do.call(rbind,all_markers_list_unique), do.call(rbind,all_markers_list_repeat))
# all_markers_list_unique <- split(all_markers_list_unique, as.character(all_markers_list_unique$cluster))
# all_markers_list_unique_out <- do.call(rbind,all_markers_list_unique)
# all_markers_list_unique_out <- all_markers_list_unique_out[order(all_markers_list_unique_out$cluster, all_markers_list_unique_out$p_val_adj ), ]
# 
# # Make a table with just the top 20 unique genes
# all_markers_list_unique_top20 <- lapply(all_markers_list_unique, function(x)if(nrow(x)> 20){x[1:20,]}else{x})
# all_markers_list_unique_top20 <- do.call(rbind,all_markers_list_unique_top20)
# all_markers_list_unique_top20 <- all_markers_list_unique_top20[order(all_markers_list_unique_top20$cluster, all_markers_list_unique_top20$p_val_adj ), ]
# write.csv(all_markers_list_unique_top20, paste0(out_tables_qc,"/cluster_characterisation_deg_top20unique.csv"), row.names=FALSE)
# 
# # 
# # #### Make feature plot for each of the top unique markers
# # for(celltype in names(all_markers_list_unique)){
# #   out_celltype <- str_replace_all(celltype,  "[[:punct:]]", " ")
# #   dir.create(paste0(out_plots_umap, "/celltype/"), showWarnings = FALSE)
# #   plot_out_fp <- paste0(out_plots_umap, "/celltype/", out_celltype, ".pdf")
# #   
# #   current_genes <- all_markers_list_unique[[celltype]]$gene[1:9]
# #   FeaturePlot(merged_srat, reduction = "umap.integrated.harmony", features=current_genes, max.cutoff="q95") 
# #   ggsave(plot_out_fp, width = 15, height = 15)
# # }
# # 
# # # plot each gene from the unique marker list
# # all_markers_list_unique_top20_list <- split(all_markers_list_unique_top20, all_markers_list_unique_top20$cluster)
# # out_dir<- paste0(out_plots_umap, "/check_unique_markers_initial_clusters/")
# # dir.create(out_dir, showWarnings = FALSE)
# # for(cluster in names(all_markers_list_unique_top20_list)){
# #   current_gene_list <- all_markers_list_unique_top20_list[[cluster]]$gene
# # 
# #   for(gene in current_gene_list){
# #     fun.FeaturePlot(merged_srat, gene, 
# #                     conversion_table = annotated_geneids_df, 
# #                     out_prefix = paste0(out_dir, "/", cluster, "_"))
# #   }
# # }
# # 
# # if(FALSE){
# #   umap_cluster_col <- clusters_colours
# #   umap_cluster_col <- umap_cluster_col[1:length(unique(merged_srat$integrated_clusters))]
# #   cluster_ids <- unique(merged_srat$integrated_clusters)
# #   names(umap_cluster_col) <- cluster_ids[order(cluster_ids)]
# #   
# #   ## Annotate clusters
# #   # Data already labelled at this point; 
# #   # File cell_cluster_label:
# #   # The file is expected to have at least 2 columns: cluster and named_celltypes
# #   # File clusterid_to_celltype_col:
# #   # The file is expected to have 2 columns, named_celltypes and colour (in hex code)
# #   # If prior to knowing cell types,
# #   # change named_celltypes column content into cluster label if prior to knowing cell types
# #   # Read in cell cluster/colour label
# #   clusterid_to_celltype <- readxl::read_excel(cell_cluster_label)
# #   clusterid_to_celltype <- clusterid_to_celltype[,c("cluster", "named_celltypes")]
# #   
# #   
# #   
# #   ## Organise annotated celltype label/colour 
# #   # Label clusters
# #   cellid_cluster <- data.frame(cellid = colnames(merged_srat), cluster = as.character(merged_srat$integrated_clusters))
# #   cellid_cluster <- merge(cellid_cluster, clusterid_to_celltype, all.x=TRUE, by="cluster")
# #   cellid_cluster <- cellid_cluster[match(colnames(merged_srat), cellid_cluster$cellid),]
# #   
# #   # Add to Seurat object
# #   merged_srat$named_celltypes <- cellid_cluster$named_celltypes
# #   #merged_srat$named_celltypes2 <- cellid_cluster$named_celltypes2
# #   
# #   # reorder the celltype factor levels
# #   named_clusters_levels <- unique(as.character(merged_srat$named_celltypes))
# #   named_clusters_levels <- named_clusters_levels[order(named_clusters_levels)]
# #   merged_srat$named_celltypes <- factor(merged_srat$named_celltypes , levels=named_clusters_levels)
# # }  
