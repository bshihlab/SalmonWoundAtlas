## Merge data and reorganise factors for metadata
# Merge doublet retained list
merged_dblKept_srat <- merge(x = srat_dblKept[[1]], y = srat_dblKept[2:6])

# Merge cells that have pasted QC filter
merged_srat <- merge(x = filtered_srat_list[[1]], y = filtered_srat_list[2:6])


## keep only genes that are expressed by more than 10 cells across all samples
count_check <- JoinLayers(merged_srat[["RNA"]])
counts_mx <- count_check@layers$counts
genes_expressed <- row.names(merged_srat)[rowSums(counts_mx) > 10]
merged_srat <- subset(merged_srat, features = genes_expressed)


# because the merge process add sample ID each time a new sample get added, redo the column names to make sure it's sampleID_cellID
merged_srat$cell_id <- colnames(merged_srat)
merged_srat$cell_id <- sapply(strsplit(merged_srat$cell_id, "_"), function(x)paste0(x[1], "_", x[length(x)]))
colnames(merged_srat) <- merged_srat$cell_id 

# reorder sample order for sample IDs so it goes from Normal skin, Early Wound to Late Wound
merged_srat$orig.ident  <- factor(as.character(merged_srat$orig.ident ), levels=c( "SRR24471130", "FHF", "W21", "W69", "W92",  "w142"))
merged_dblKept_srat$orig.ident  <- factor(as.character(merged_dblKept_srat$orig.ident ), levels=c( "SRR24471130", "FHF", "W21", "W69", "W92",  "w142"))

# Create treatment groups
merged_srat$treatment_group <- factor(as.character(merged_srat$orig.ident), levels = c("FHF", "SRR24471130", "W21", "W69", "W92", "w142"))
levels(merged_srat$treatment_group) <- c("Skin", "Skin", "WoundEarly", "WoundEarly", "WoundLate",  "WoundLate")


## Normalise data and run PCA
merged_srat <- NormalizeData(merged_srat)

# remove ribosomal genes from variable features
merged_srat <- FindVariableFeatures(merged_srat, nfeatures= 2500)
ribosomal_genes <- annotated_geneids_df[(grepl(pattern = "^RPS|^RPL", x =  annotated_geneids_df$gene_name_base, ignore.case=TRUE)) | (grepl(pattern = "^RPS|^RPL", x =  annotated_geneids_df$zebrafish.genename, ignore.case=TRUE))|  (grepl(pattern = "^RPS|^RPL", x =  annotated_geneids_df$human.genename, ignore.case=TRUE)), ]
ribosomal_genes <- fun.find_matching_features(merged_srat, ribosomal_genes$gene_id, annotated_geneids_df, remove.dash = TRUE)
ribosomal_genes_id <- ribosomal_genes$geneid
ribosomal_genes_id <- ribosomal_genes_id[ribosomal_genes_id %in% row.names(merged_srat)]

merged_srat <- AddModuleScore(
  object = merged_srat,
  features = list(ribosomal_genes_id),
  name = 'Rb.score'
)
VariableFeatures(merged_srat) <- setdiff(VariableFeatures(merged_srat), ribosomal_genes_id)

# scale data using all genes
merged_srat <- ScaleData(merged_srat, features=row.names(merged_srat))

# perform PCA using defaults (removing ribosomal genes from highly variable genes)
merged_srat <- RunPCA(merged_srat, npcs=100)

# Elbow plot
ElbowPlot(merged_srat, ndims = 100, reduction = "pca")



#merged_srat$percent.rb <- merged_srat$Rb.score1


# find keratin genes
KRT_feature <- annotated_geneids_df$gene_id[grep("^KRT", annotated_geneids_df$gene.name , ignore.case=TRUE)]
KRT_feature <- c(KRT_feature, annotated_geneids_df$gene_id[grep("^KRT", annotated_geneids_df$orthologue , ignore.case=TRUE)])
KRT_feature <- c(KRT_feature, annotated_geneids_df$gene_id[grep("^KRT", annotated_geneids_df$gene_name_base , ignore.case=TRUE)])
KRT_feature <- unique(KRT_feature)
KRT_feature <- KRT_feature[KRT_feature %in% row.names(merged_srat)]
merged_srat[["percent.krt"]] <- PercentageFeatureSet(merged_srat, features = KRT_feature)


## Normalise using SCT
merged_srat <- SCTransform(merged_srat,   vars.to.regress = c("Phase", "Rb.score1"))

## Integrate data using harmony, then join the RNA layers
merged_srat <- IntegrateLayers(object = merged_srat, normalization.method = "SCT",
                               method = HarmonyIntegration, orig.reduction = "pca",  
                               new.reduction = "integrated.harmony", verbose = TRUE, dims=1:25)
merged_srat_backup <- merged_srat

merged_srat[["RNA"]] <- JoinLayers(merged_srat[["RNA"]])

# 
# ## Remove cluster with high RPL/RPS and low nCount_RNA
# # find pc with highest correlation to ribosome, exclude it from find neighbours/ umap
analysed_dims <- 1:52
max_pc <- merged_srat@reductions$integrated.harmony@cell.embeddings
pc_corr <- apply(max_pc, 2, FUN=function(x)cor(x, merged_srat$Rb.score1))
pc_corr <- abs(pc_corr)
ribo_pc <- strsplit( names(pc_corr)[which.max(pc_corr)], "_")[[1]][2]
analysed_dims <- analysed_dims[analysed_dims != ribo_pc]

# find clusters
# #merged_srat <- merged_srat_backup
# merged_srat <- FindNeighbors(merged_srat, reduction = "integrated.harmony",
#                              dims = analysed_dims, graph.name="clustering1")
# merged_srat <- FindClusters(merged_srat, resolution = 1.2, cluster.name = "integrated_clusters", graph.name="clustering1")
# merged_srat$integrated_clusters <- as.character(merged_srat$integrated_clusters)
# Idents(merged_srat) <- "integrated_clusters"
# merged_srat <- RunUMAP(merged_srat, dims = analysed_dims,
#                        reduction = "integrated.harmony",
#                        reduction.name = "umap.integrated.harmony",
#                        n.neighbors = 30L, seed.use = 42)
# #merged_srat <- FindSubCluster(merged_srat, "16", graph.name = "clustering1", resolution = 0.2)
# #merged_srat$integrated_clusters <- ifelse(as.character(merged_srat$integrated_clusters) == "16", merged_srat$sub.cluster, as.character(merged_srat$integrated_clusters))
# 
# p <-DimPlot(merged_srat, group.by="integrated_clusters",
#             reduction = "umap.integrated.harmony", label=TRUE)
# tmp_colour <- c(clusters_colours, clusters_colours)
# tmp_colour <- tmp_colour[1:length(unique(merged_srat$integrated_clusters))]
# names(tmp_colour) <- unique(merged_srat$integrated_clusters)
# p + scale_colour_manual(values=tmp_colour)
# # 
# 
# FeaturePlot(merged_srat, feature = "Rb.score1", 
#             reduction = "umap.integrated.harmony", max.cutoff = "q95")
# FeaturePlot(merged_srat, feature = "percent.krt", 
#             reduction = "umap.integrated.harmony", max.cutoff = "q95")
# 

# 
# merged_srat$integrated_clusters <- as.character(merged_srat$integrated_clusters)
# merged_srat_sub <- subset(merged_srat, !integrated_clusters %in% c("16_0"))
# DimPlot(merged_srat_sub, group.by="integrated_clusters", 
#         reduction = "umap.integrated.harmony", label=TRUE)
# 
# 
# 
# ## Split layers and reintegrate
# DefaultAssay(merged_srat_sub) <- "RNA"
# merged_srat_sub[['SCT']] <- NULL
# merged_srat_sub[["RNA"]] <- split(merged_srat_sub[["RNA"]], f=merged_srat_sub$orig.ident)
# merged_srat_sub <- NormalizeData(merged_srat_sub)
# 
# # find variable features (needed for PCA) and remove ribosomal genes from variable features
# merged_srat <- FindVariableFeatures(merged_srat, nfeatures= 2500)
# ribosomal_genes <- annotated_geneids_df[(grepl(pattern = "^RPS|^RPL", x =  annotated_geneids_df$gene_name_base, ignore.case=TRUE)) | (grepl(pattern = "^RPS|^RPL", x =  annotated_geneids_df$zebrafish.genename, ignore.case=TRUE))|  (grepl(pattern = "^RPS|^RPL", x =  annotated_geneids_df$human.genename, ignore.case=TRUE)), ]
# ribosomal_genes <- fun.find_matching_features(merged_srat, ribosomal_genes$gene_id, annotated_geneids_df, remove.dash = TRUE)
# ribosomal_genes_id <- ribosomal_genes$geneid
# ribosomal_genes_id <- ribosomal_genes_id[ribosomal_genes_id %in% row.names(merged_srat)]
# 
# merged_srat <- AddModuleScore(
#   object = merged_srat,
#   features = list(ribosomal_genes_id),
#   name = 'Rb.score'
# )
# VariableFeatures(merged_srat) <- setdiff(VariableFeatures(merged_srat), ribosomal_genes_id)
# 
# 
# merged_srat_sub <- ScaleData(merged_srat_sub, features=row.names(merged_srat_sub))
# merged_srat_sub <- RunPCA(merged_srat_sub, npcs=100)
# merged_srat_sub <- SCTransform(merged_srat_sub)
# 
# ## reintegrate layers 
# merged_srat_sub <- IntegrateLayers(object = merged_srat_sub, normalization.method = "SCT", 
#                                    method = HarmonyIntegration, orig.reduction = "pca",  
#                                    new.reduction = "integrated.harmony", verbose = TRUE, dims=1:25)
# merged_srat_sub[["RNA"]] <- JoinLayers(merged_srat_sub[["RNA"]])
# 
# merged_srat_backup <- merged_srat
# merged_srat <- merged_srat_sub
# 
