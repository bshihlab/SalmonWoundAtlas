## Save data
# Dimension reduction coordinates and cell metadata
# These are used as input for scVelo
# umap3d_coord <- as.data.frame(merged_srat@reductions$umap.integrated.harmony3D@cell.embeddings)
# umap3d_with_col <- cbind(umap3d_coord, data.frame(cell_id = colnames(merged_srat), cluster=as.character(merged_srat@meta.data$named_celltypes)))
# harmony_embedding <- merged_srat@reductions$integrated.harmony@cell.embeddings
# write.csv(umap3d_coord, paste0(out_tables, "/umap3d_coord.csv"))
# write.csv(umap3d_with_col, paste0(out_tables, "/umap3d_coord_withclusterlabel.csv"))
#write.csv(harmony_embedding, paste0(out_tables, "/harmony_embedding.csv"), row.names = TRUE)
write.csv(merged_srat@meta.data, paste0(out_tables, "/metadata.csv"), row.names = TRUE)

# Save RDS for integrated Seurat object so far
DefaultAssay(merged_srat) <- "RNA"
saveRDS(merged_srat, paste0(out_rds, "/merged_srat.rds"))
merged_srat <- readRDS(paste0(out_rds, "/merged_srat.rds"))


# Save RDS at this stage
#DefaultAssay(merged_srat_backup) <- "RNA"
#saveRDS(merged_srat_backup, paste0(out_rds, "/merged_srat_pre_ribocluster_rm.rds"))

# save seurat object as h5 AnnData
obj2 <- Convert_Assay(seurat_object = merged_srat, assay = "RNA", convert_to = "V3")
obj2$named_celltypes <- as.character(obj2$named_celltypes)
obj2$celltype_fig1 <- as.character(obj2$celltype_fig1)
obj2$celltype_Stdeconvolution <- as.character(obj2$celltype_Stdeconvolution)
harmony_reduction_embedding <- as.data.frame(obj2@reductions$umap.integrated.harmony@cell.embeddings)

writeMM(obj2@assays$RNA@counts,paste0(out_rds,"/h5ad_count.mx"))
writeMM(obj2@assays$RNA@data,paste0(out_rds,"/h5ad_data.mx"))
out_metadata <- as.data.frame(obj2@meta.data)
write.csv(out_metadata[, c("cell_id","orig.ident", "named_celltypes", "celltype_fig1", "celltype_Stdeconvolution", "treatment_group")], 
          paste0(out_rds,"/h5ad_annotation_brief.csv"), row.names=FALSE)
write.csv(data.frame("barcode"=out_metadata$cell_id), paste0(out_rds,"/h5ad_barcodes.csv"), row.names=FALSE)
write.csv(data.frame("gene"=row.names(obj2)), paste0(out_rds,"/h5ad_genes.csv"), row.names=FALSE)
write.csv(data.frame("Cell" = out_metadata$cell_id, "cell_type"=out_metadata$named_celltypes), paste0(out_rds,"/h5ad_barcode_and_celltype.csv"), row.names=FALSE)
write.csv(harmony_reduction_embedding, paste0(out_rds,"/h5ad_embedding.csv"), row.names=FALSE)
#harmony_coord <- as.data.frame(obj2@reductions$umap.integrated.harmony3D@cell.embeddings)
#write.csv(harmony_coord, "analysis/ribosomal_cluster_keep/rds/ccc_harmony_coord.csv", row.names=FALSE)

