#### Go through each of the cellbender matrix and 
# perform QC, doublet removal and decontX on each individually

filtered_srat_list <- list()
srat_dblKept <- list()
Cell_count_info <- list()


ncount_filter <- list(
  "SRR24471130" = c(180, 20000),
  "FHF" = c(180, 20000),
  "W21" = c(400, 20000),
  "W69" = c(400, 20000),
  "W92" = c(400, 20000),
  "w142" = c(400, 20000)
)

nfeature_filter <- list(
  "SRR24471130" = c(150, 8000),
  "FHF" = c(150, 8000),
  "W21" = c(300, 9000),
  "W69" = c(300, 9000),
  "W92" = c(300, 9000),
  "w142" = c(300, 9000)
)

# Loop through each sample

#### nCount_RNA and nFeature_RNA filter set by auto-filter in Loup Browser 

for(sample_id in names(cellbender_matrix_list)){
  obj <- CreateSeuratObject(counts = cellbender_matrix_list[[sample_id]],
                            project=sample_id,min.features=100)
  row.names(obj) <- geneid_full$gene_id
  
  
  # QC filter based on fixed filter
  qcIndividual_nFeature.low.filter = nfeature_filter[[sample_id]][1]
  qcIndividual_nFeature.upper.filter = nfeature_filter[[sample_id]][2]
  qcIndividual_nCount.low.filter = ncount_filter[[sample_id]][1]
  qcIndividual_nCount.upper.filter = ncount_filter[[sample_id]][2]

  # Record how many cells were there initially from cellbender output
  cell_count_initial <- ncol(obj)
  
  ## Read in passed 10x cellranger barcode, and remove any barcode that is absent from there
  cellranger_barcode_fp <- paste0("analysis/cellranger/", sample_id, "/outs/filtered_feature_bc_matrix/barcodes.tsv.gz")
  cellranger_filtered_barcode <- read.table(gzfile(cellranger_barcode_fp))
  cellranger_filtered_barcode <- strsplit(cellranger_filtered_barcode$V1, "-")
  cellranger_filtered_barcode <- sapply(cellranger_filtered_barcode, function(x)x[1])
  #cellranger_filtered_barcode <- paste0(sample_id, "_", cellranger_filtered_barcode)
  cellranger_filtered_barcode <- cellranger_filtered_barcode[cellranger_filtered_barcode %in% colnames(obj)]
  
  obj <- obj[,cellranger_filtered_barcode]
  cell_count_rm10xfilter <- ncol(obj)

  
  # Find ribosomal and mitochondrial genes, check orthologue, gene.name, gene_name_base
  # Use this to calculate the % ribosome reads
  RB_feature <- annotated_geneids_df$gene_id[grep("^RP[SL]", annotated_geneids_df$gene.name , ignore.case=TRUE)]
  RB_feature <- c(RB_feature, annotated_geneids_df$gene_id[grep("^RP[SL]", annotated_geneids_df$orthologue , ignore.case=TRUE)])
  RB_feature <- c(RB_feature, annotated_geneids_df$gene_id[grep("^RP[SL]", annotated_geneids_df$gene_name_base , ignore.case=TRUE)])
  RB_feature <- unique(RB_feature)
  obj[["percent.rb"]] <- PercentageFeatureSet(obj, features = RB_feature)
  
  ## Visualise the distribution and the threshold cutoffs for nCount and nFeature
  print("## qc_doublet_filter.r: Visualise the distribution and the threshold cutoffs for nCount and nFeature")
  metadata <- obj@meta.data
  metadata %>% 
    ggplot(aes(x=nCount_RNA)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density") +
    geom_vline(xintercept = c(qcIndividual_nCount.low.filter, qcIndividual_nCount.upper.filter))
  ggsave(paste0(out_plots_qc, "/nCount_RNA_", sample_id, ".pdf"))
  
  metadata %>% 
    ggplot(aes(x=nFeature_RNA)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density") +
    geom_vline(xintercept = c(qcIndividual_nFeature.low.filter, qcIndividual_nFeature.upper.filter)) 
  ggsave(paste0(out_plots_qc, "/nFeature_RNA_", sample_id, ".pdf"))
  
  
  print("## qc_doublet_filter.r: Subset data, removing low quality cells (low gene count/ number of genes)")
  print("## qc_doublet_filter.r: QC setting based on auto thresholds from cellranger")

  
  obj <- subset(obj, subset = nFeature_RNA > qcIndividual_nFeature.low.filter &
                  nFeature_RNA < qcIndividual_nFeature.upper.filter &
                  nCount_RNA > qcIndividual_nCount.low.filter  &
                  nCount_RNA < qcIndividual_nCount.upper.filter &
                  percent.rb < 40)
  
  
  # Count the number of cells at this stage
  cell_count_qc <- ncol(obj)
  
  
  ## Normalise and scale subsetted data
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj, features = row.names(obj))
  obj <- RunPCA(obj)
  print(dim(obj))
  
  print("## qc_doublet_filter.r: Label cell cycle scoring")
  ## Label cell cycle scoring
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  s.genes_salmon <- annotated_geneids_df$gene_id[grep(paste0(s.genes, collapse = "|"), annotated_geneids_df$gene.name , ignore.case=TRUE)]
  s.genes_salmon <- c(s.genes_salmon, annotated_geneids_df$gene_id[grep(paste0(s.genes, collapse = "|"), annotated_geneids_df$orthologue , ignore.case=TRUE)])
  s.genes_salmon <- c(s.genes_salmon, annotated_geneids_df$gene_id[grep(paste0(s.genes, collapse = "|"), annotated_geneids_df$gene_name_base , ignore.case=TRUE)])
  s.genes_salmon <- unique(s.genes_salmon)
  g2m.genes_salmon <- annotated_geneids_df$gene_id[grep(paste0(s.genes, collapse = "|"), annotated_geneids_df$gene.name , ignore.case=TRUE)]
  g2m.genes_salmon <- c(g2m.genes_salmon, annotated_geneids_df$gene_id[grep(paste0(s.genes, collapse = "|"), annotated_geneids_df$orthologue , ignore.case=TRUE)])
  g2m.genes_salmon <- c(g2m.genes_salmon, annotated_geneids_df$gene_id[grep(paste0(s.genes, collapse = "|"), annotated_geneids_df$gene_name_base , ignore.case=TRUE)])
  g2m.genes_salmon <- unique(g2m.genes_salmon)
  obj <- CellCycleScoring(obj, s.features = s.genes_salmon, g2m.features = g2m.genes_salmon, set.ident = TRUE)
  
  
  srat_dblKept[[sample_id]] <- obj
  
  print("## qc_doublet_filter.r: remove doublets using scDblFinder")
  ## remove doublets using scDblFinder
  # remove doublets
  sce <- as.SingleCellExperiment(obj)
  set.seed(123)
  sce <- scDblFinder(sce, dbr=0.15)
  sce <- sce[,!as.character(sce$scDblFinder.class) == "doublet"]
  #sce <- decontX(sce)
  
  options(Seurat.object.assay.version = "v5")
  
  ## Use decontX to remove ambient RNA. 
  # After Cellbender, there are still some ambient RNA 
  # (red blood cell genes appearing in other cell types. Improvement after decontX)
  obj <- CreateSeuratObject(counts = counts(sce), min.features=100, project=sample_id, meta.data=as.data.frame(colData(sce)))
  #obj <- CreateSeuratObject(counts = counts(sce), min.cells=10, min.features=100, project=sample_id, meta.data=as.data.frame(colData(sce)))
  #obj <- CreateSeuratObject(counts = counts(sce), project=sample_id, meta.data=as.data.frame(colData(sce)))
  obj <- RenameCells(obj, new.names = paste0(sample_id, "_", colnames(obj)))
  
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj, features = row.names(obj))
  obj <- RunPCA(obj)
  
  
  # Record the cell count after doublet removal
  cell_count_doubletrm <- ncol(obj)
  filtered_srat_list[[sample_id]] <- obj
  obj2 <- obj
  obj2 <- FindNeighbors(obj2,  dims = 1:40, graph.name="clustering0")
  obj2 <- RunUMAP(obj2,  dims = 1:30, reduction = "pca")
  obj2 <- FindClusters(obj2,  graph.name="clustering0")
  
  # add gene annotation back into obj2
  DefaultAssay(obj2) <- "RNA"
  
  print("## qc_doublet_filter.r: plot clusters")
  ## plot clusters
  p <- DimPlot(obj2, group.by = "seurat_clusters", label=TRUE)
  p + theme(legend.position = "none")
  ggsave(paste0(out_plots_umap_individual, "/", sample_id, ".pdf"), width = 4, height = 4)
  
  # Record the number of cells at each stage
  Cell_count_info[[sample_id]] <- c('sample_id' = sample_id, 
                                    'initial' = cell_count_initial, 
                                    'cellranger_filter' = cell_count_rm10xfilter,
                                    'qc' = cell_count_qc, 
                                    'doubletrm' = cell_count_doubletrm)
}
