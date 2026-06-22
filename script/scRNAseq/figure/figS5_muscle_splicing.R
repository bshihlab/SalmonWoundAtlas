## generate umap featureplot for FigS4
# organise output folder
current_out_dir <- paste0(out_plots_manuscript, "/figs5_muscle_splicing/")
dir.create(current_out_dir, showWarnings = FALSE)

### Need table from broad cell group
source("script/scRNAseq/table/table_deg_within_broadcellgroup.R")

#### check spliced vs unspliced reads
splice_genes_all <- list()
splice_genes_fg <- list()
splice_genes_nonfg <- list()
splice_cell_other <- list()
splice_cell_MFfg <- list()
splice_cell_MF <- list()
splice_cell_rbc <- list()

all_metadata <- merged_srat@meta.data
metabolic_cellid <- row.names(all_metadata)[grepl("fg", as.character(all_metadata$named_celltypes))]
other_myofibre_cellid <- row.names(all_metadata)[grepl("myofibre", as.character(all_metadata$named_celltypes))]
other_myofibre_cellid <- other_myofibre_cellid[!other_myofibre_cellid %in% metabolic_cellid]
other_rbc_cellid <- row.names(all_metadata)[grepl("erythrocyte", as.character(all_metadata$named_celltypes))]
all_other_cellid <- row.names(all_metadata)
all_other_cellid <- all_other_cellid[!all_other_cellid %in% c(other_rbc_cellid, metabolic_cellid, other_myofibre_cellid)]

gene_ids_metabolic <- metabolic_both$gene_id

### calculate proportion spliced/unspliced
for(sample_id in names(cellbender_matrix_list)){
  
  spliced_mx_fp <- paste0("analysis/star/", sample_id, "/Solo.out/Velocyto/raw/spliced.mtx")
  unspliced_mx_fp <- paste0("analysis/star/", sample_id, "/Solo.out/Velocyto/raw/unspliced.mtx")
  cells_fp <-  paste0("analysis/star/", sample_id, "/Solo.out/Velocyto/raw/barcodes.tsv")
  genes_fp <-  paste0("analysis/star/", sample_id, "/Solo.out/Velocyto/raw/features.tsv")
  
  # read matrix
  # read in velocyto unspliced
  cells_fp <-  paste0("analysis/star/", sample_id, "/Solo.out/Velocyto/raw/barcodes.tsv")
  genes_fp <-  paste0("analysis/star/", sample_id, "/Solo.out/Velocyto/raw/features.tsv")
  unspliced_mx_fp <- paste0("analysis/star/", sample_id, "/Solo.out/Velocyto/raw/unspliced.mtx")
  spliced_mx_fp <- paste0("analysis/star/", sample_id, "/Solo.out/Velocyto/raw/spliced.mtx")
  ambiguous_mx_fp <- paste0("analysis/star/", sample_id, "/Solo.out/Velocyto/raw/ambiguous.mtx")
  
  spliced_mx <- ReadMtx(mtx=spliced_mx_fp, cells=cells_fp, features=genes_fp)
  unspliced_mx <- ReadMtx(mtx=unspliced_mx_fp, cells=cells_fp, features=genes_fp)
  ambiguous_mx <- ReadMtx(mtx=ambiguous_mx_fp, cells=cells_fp, features=genes_fp)
  
  all_splice_umi <- colSums(spliced_mx)
  all_unsplice_umi <- colSums(unspliced_mx)
  all_ambiguous_umi <- colSums(ambiguous_mx)
  
  all_umi <- all_splice_umi + all_unsplice_umi + all_ambiguous_umi
  
  row.names(spliced_mx) <- geneid_full$gene_id
  row.names(unspliced_mx) <-  geneid_full$gene_id
  colnames(spliced_mx) <- paste0(sample_id, "_", colnames(spliced_mx) )
  colnames(unspliced_mx) <- paste0(sample_id, "_", colnames(unspliced_mx) )
  
  # spliced/unspliced for all genes across all other cell types
  all_spliced <- colSums(spliced_mx)
  all_unspliced <- colSums(unspliced_mx)
  percent_spliced_cell <- all_spliced/all_umi*100
  
  spliced_metabolic_genes <- colSums(spliced_mx[row.names(spliced_mx) %in% gene_ids_metabolic,])
  unspliced_metabolic_genes <- colSums(unspliced_mx[row.names(unspliced_mx) %in% gene_ids_metabolic,])
  ambiguous_metabolic_genes <- colSums(ambiguous_mx[row.names(ambiguous_mx) %in% gene_ids_metabolic,])
  percent_spliced_metabolic_genes <- spliced_metabolic_genes/(unspliced_metabolic_genes + spliced_metabolic_genes + ambiguous_metabolic_genes)*100

  spliced_nonmetabolic_genes <- colSums(spliced_mx[!row.names(spliced_mx) %in% gene_ids_metabolic,])
  unspliced_nonmetabolic_genes <- colSums(unspliced_mx[!row.names(unspliced_mx) %in% gene_ids_metabolic,])
  ambiguous_nonmetabolic_genes <- colSums(ambiguous_mx[!row.names(ambiguous_mx) %in% gene_ids_metabolic,])
  percent_spliced_nonmetabolic_genes <- spliced_nonmetabolic_genes/(unspliced_nonmetabolic_genes + spliced_nonmetabolic_genes + ambiguous_nonmetabolic_genes)*100
  
  
  # spliced/unspliced for genes upregulated in metabolic cell types
  gene_spliced_metabolic <- rowSums(spliced_mx[,colnames(spliced_mx) %in% metabolic_cellid])
  gene_unspliced_metabolic <- rowSums(unspliced_mx[,colnames(unspliced_mx) %in% metabolic_cellid])
  gene_ambiguous_metabolic <- rowSums(ambiguous_mx[,colnames(ambiguous_mx) %in% metabolic_cellid])
  percent_spliced_gene_metabolic <- gene_spliced_metabolic/(gene_unspliced_metabolic + gene_spliced_metabolic + gene_ambiguous_metabolic)*100
  
  gene_spliced_other_myofibre <- rowSums(spliced_mx[,colnames(spliced_mx) %in% other_myofibre_cellid])
  gene_unspliced_other_myofibre <- rowSums(unspliced_mx[,colnames(unspliced_mx) %in% other_myofibre_cellid])
  gene_ambiguous_other_myofibre <- rowSums(ambiguous_mx[,colnames(ambiguous_mx) %in% other_myofibre_cellid])
  percent_spliced_gene_other_myofibre <- gene_spliced_other_myofibre/(gene_spliced_other_myofibre + gene_ambiguous_other_myofibre+ gene_unspliced_other_myofibre)*100
  
  gene_spliced_rbc <- rowSums(spliced_mx[,colnames(spliced_mx) %in% other_rbc_cellid])
  gene_unspliced_rbc <- rowSums(unspliced_mx[,colnames(unspliced_mx) %in% other_rbc_cellid])
  gene_ambiguous_rbc <- rowSums(ambiguous_mx[,colnames(ambiguous_mx) %in% other_rbc_cellid])
  percent_spliced_rbc <- gene_spliced_rbc/(gene_spliced_rbc + gene_ambiguous_rbc+ gene_unspliced_rbc)*100
  
  gene_spliced_other <- rowSums(spliced_mx[,colnames(spliced_mx) %in% all_other_cellid])
  gene_unspliced_other <- rowSums(unspliced_mx[,colnames(unspliced_mx) %in% all_other_cellid])
  gene_ambiguous_other <- rowSums(ambiguous_mx[,colnames(ambiguous_mx) %in% all_other_cellid])
  percent_spliced_gene_other <- gene_spliced_other/(gene_unspliced_other + gene_ambiguous_other + gene_spliced_other)*100
  
  # record calculation
  splice_genes_all[[sample_id]]  = percent_spliced_cell
  splice_genes_fg[[sample_id]]  = percent_spliced_metabolic_genes
  splice_genes_nonfg[[sample_id]]  = percent_spliced_nonmetabolic_genes
  splice_cell_other[[sample_id]]  = percent_spliced_gene_other
  splice_cell_MF[[sample_id]]  = percent_spliced_gene_other_myofibre
  splice_cell_MFfg[[sample_id]]  = percent_spliced_gene_metabolic
  splice_cell_rbc[[sample_id]]  = percent_spliced_rbc
  
}


### annotate merged_srat with splice info
additional_splice_annotation <- list(
  splice_all_genes = splice_genes_all,
  splice_fg_genes = splice_genes_fg,
  splice_nonfg_genes = splice_genes_nonfg
)
for(current_splice_annotation in names(additional_splice_annotation)){
  current_splice <- do.call(c, additional_splice_annotation[[current_splice_annotation]])
  names(current_splice) <- sapply(strsplit( names(current_splice), "\\."), function(x)x[[2]])
  current_splice <- current_splice[names(current_splice) %in% colnames(merged_srat)]
  current_splice <- current_splice[match(colnames(merged_srat), names(current_splice))]
  names(current_splice) <- colnames(merged_srat)
  merged_srat[[current_splice_annotation]] <- current_splice
  VlnPlot(merged_srat, current_splice_annotation) +  guides(colour="none")
  ggsave(paste0(current_out_dir, "/celltype_", current_splice_annotation, ".pdf"), width=15, height = 5)
  
}
  
percent_spliced_df <- data.frame(celltype = merged_srat$named_celltypes, 
                                 percent_spliced = merged_srat$splice_all_genes, 
                                 nFeature_RNA = merged_srat$nFeature_RNA)
aggregate(percent_spliced_df, percent_spliced~celltype, mean, na.rm=TRUE)


### Compare the %splicing for metabolism associated genes 
# in metabolic muscle and other cells
splice_cell_other_df <- do.call(cbind,splice_cell_other)
splice_cell_MF_df <- do.call(cbind,splice_cell_MF)
splice_cell_MFfg_df <- do.call(cbind,splice_cell_MFfg)
splice_cell_rbc_df <- do.call(cbind,splice_cell_rbc)

avg_gene_other <- rowMeans(splice_cell_other_df, na.rm=TRUE)
avg_gene_MF <- rowMeans(splice_cell_MF_df, na.rm=TRUE)
avg_gene_MFfg <- rowMeans(splice_cell_MFfg_df, na.rm=TRUE)
avg_gene_rbc <- rowMeans(splice_cell_rbc_df, na.rm=TRUE)


### plots to make
# Correlation at gene level between muscle, metabolic (MC-M) and non-MC-M
merged_srat$MF_fg <- ifelse(colnames(merged_srat) %in% metabolic_cellid, "MFfg", "other")
all_gene_avg_exprs <- AverageExpression(merged_srat, group.by = "MF_fg")
all_gene_avg_exprs <- as.data.frame(all_gene_avg_exprs$RNA)
all_gene_avg_exprs <- all_gene_avg_exprs[names(avg_gene_other),]
# keeping only genes that are expressed in MFfg
all_gene_avg_exprs <- all_gene_avg_exprs[all_gene_avg_exprs$MFfg>0,]
genes_expressed_in_MFfg <- row.names(all_gene_avg_exprs)

# match the other averages so only genes that are exprssed in MFfg are kept
avg_gene_other <- avg_gene_other[genes_expressed_in_MFfg]
avg_gene_MF <- avg_gene_MF[genes_expressed_in_MFfg]
avg_gene_MFfg <- avg_gene_MFfg[genes_expressed_in_MFfg]
avg_gene_rbc <- avg_gene_rbc[genes_expressed_in_MFfg]

all_genes_spliced <- data.frame(gene = genes_expressed_in_MFfg, 
                                other = avg_gene_other, 
                                MF = avg_gene_MF, 
                                rbc = avg_gene_rbc,
                                MFfg = avg_gene_MFfg, 
                                MFfg_exprs = all_gene_avg_exprs$MFfg)
all_genes_spliced$gene_type <- ifelse(all_genes_spliced$gene %in% metabolic_both$gene_id,
                                      "MF_fg-upregulated genes", "other genes") 

ggplot(all_genes_spliced, aes(x=other, y=MF)) + 
  geom_point(size= 0.1, alpha = 0.6) + theme_bw() + scale_colour_viridis(direction = -1) +
  #facet_wrap(~gene_type)+ 
  geom_abline(slope=1, intercept = 0, colour="red") +
  guides(colour="none") +
  xlab("splicing in \nother cell types (%)") + ylab("splicing in \nmyofibre (%)")
ggsave(paste0(current_out_dir, "/splicing_cor_myofibre.pdf"), width=2.25, height = 2.25)

ggplot(all_genes_spliced, aes(x=other, y=MFfg)) + 
  geom_point(size= 0.1, alpha = 0.6) + theme_bw() + scale_colour_viridis(direction = -1) +
  #facet_wrap(~gene_type)+ 
  geom_abline(slope=1, intercept = 0, colour="red") +
  guides(colour="none") +  
  xlab("splicing in \nother cell types (%)") + ylab("splicing in \nmyofibre(fg) (%)")
ggsave(paste0(current_out_dir, "/splicing_cor_myofibre-fg.pdf"), width=2.25, height = 2.25)

ggplot(all_genes_spliced, aes(x=other, y=rbc)) + 
  geom_point(size= 0.1, alpha = 0.6) + theme_bw() + scale_colour_viridis(direction = -1) +
  #facet_wrap(~gene_type)+ 
  geom_abline(slope=1, intercept = 0, colour="red") +
  guides(colour="none") +  
  xlab("splicing in \nother cell types (%)") + ylab("splicing in \nerythrocyte (%)")
ggsave(paste0(current_out_dir, "/splicing_cor_erythrocyte.pdf"), width=2.25, height = 2.25)


## Plot nFeature_RNA
VlnPlot(merged_srat, "nFeature_RNA") +  guides(colour="none")


## Other qc
for(gene in c("ENSSSAG00000044737","splice_all_genes", "nFeature_RNA", "percent.krt")){
  fun.FeaturePlot(merged_srat, gene, 
                  file_type = "pdf",
                  conversion_table = annotated_geneids_df, 
                  out_prefix = current_out_dir,
                  width = 6, height = 5, point_size = 0.01)
}

## get summary on the spliced % for MF and MF-fg
splice_stats <- data.frame(celltype=merged_srat$named_celltypes, spliced=merged_srat$splice_all_genes)
out_df <- aggregate(data=splice_stats, spliced~celltype, mean)
write.csv(out_df, paste0(out_tables, "/mean_splice_per_celltype.csv"), row.names = FALSE)

out_df <- aggregate(data=splice_stats, spliced~celltype, median)
write.csv(out_df, paste0(out_tables, "/median_splice_per_celltype.csv"), row.names = FALSE)

