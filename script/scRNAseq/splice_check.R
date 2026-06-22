#### check spliced vs unspliced reads

splice_myofibre_cells <- list()
splice_myofibre_gene_other <- list()
splice_myofibre_gene_metabolic <- list()
splice_myofibre_gene_other_myofibre <- list()
splice_myofibre_gene_erythrocyte <- list()

all_metadata <- merged_srat@meta.data
metabolic_cellid <- row.names(all_metadata)[grepl("metabo", as.character(all_metadata$named_celltypes))]
other_myofibre_cellid <- row.names(all_metadata)[grepl("myofibre", as.character(all_metadata$named_celltypes))]
other_myofibre_cellid <- other_myofibre_cellid[!other_myofibre_cellid %in% metabolic_cellid]
other_rbc_cellid <- row.names(all_metadata)[grepl("erythrocyte", as.character(all_metadata$named_celltypes))]
all_other_cellid <- row.names(all_metadata)
all_other_cellid <- all_other_cellid[!all_other_cellid %in% c(other_rbc_cellid, metabolic_cellid, other_myofibre_cellid)]

### calculate proportion spliced/unspliced
for(sample_id in names(cellbender_matrix_list)){
  
  spliced_mx_fp <- paste0("analysis/star/", sample_id, "/Solo.out/Velocyto/filtered/spliced.mtx")
  unspliced_mx_fp <- paste0("analysis/star/", sample_id, "/Solo.out/Velocyto/filtered/unspliced.mtx")
  cells_fp <-  paste0("analysis/star/", sample_id, "/Solo.out/Velocyto/filtered/barcodes.tsv")
  genes_fp <-  paste0("analysis/star/", sample_id, "/Solo.out/Velocyto/filtered/features.tsv")
  
  # read matrix
  spliced_mx <- ReadMtx(mtx=spliced_mx_fp, cells=cells_fp, features=genes_fp)
  unspliced_mx <- ReadMtx(mtx=unspliced_mx_fp, cells=cells_fp, features=genes_fp)
  row.names(spliced_mx) <- geneid_full$gene_id
  row.names(unspliced_mx) <-  geneid_full$gene_id
  colnames(spliced_mx) <- paste0(sample_id, "_", colnames(spliced_mx) )
  colnames(unspliced_mx) <- paste0(sample_id, "_", colnames(unspliced_mx) )
  
  # spliced/unspliced for all genes across all other cell types
  all_spliced <- colSums(spliced_mx)
  all_unspliced <- colSums(unspliced_mx)
  percent_spliced_cell <- all_spliced/(all_unspliced + all_spliced)*100
  
  # spliced/unspliced for genes upregulated in metabolic cell types
  gene_spliced_metabolic <- rowSums(spliced_mx[,colnames(spliced_mx) %in% metabolic_cellid])
  gene_unspliced_metabolic <- rowSums(unspliced_mx[,colnames(spliced_mx) %in% metabolic_cellid])
  percent_spliced_gene_metabolic <- gene_spliced_metabolic/(gene_unspliced_metabolic + gene_spliced_metabolic)*100
  
  gene_spliced_other_myofibre <- rowSums(spliced_mx[,colnames(spliced_mx) %in% other_myofibre_cellid])
  gene_unspliced_other_myofibre <- rowSums(unspliced_mx[,colnames(spliced_mx) %in% other_myofibre_cellid])
  percent_spliced_gene_other_myofibre <- gene_spliced_other_myofibre/(gene_spliced_other_myofibre + gene_unspliced_other_myofibre)*100
  
  gene_spliced_rbc <- rowSums(spliced_mx[,colnames(spliced_mx) %in% other_rbc_cellid])
  gene_unspliced_rbc <- rowSums(unspliced_mx[,colnames(spliced_mx) %in% other_rbc_cellid])
  percent_spliced_rbc <- gene_spliced_rbc/(gene_spliced_rbc + gene_unspliced_rbc)*100
  
  gene_spliced_other <- rowSums(spliced_mx[,colnames(spliced_mx) %in% all_other_cellid])
  gene_unspliced_other <- rowSums(unspliced_mx[,colnames(spliced_mx) %in% all_other_cellid])
  percent_spliced_gene_other <- gene_spliced_other/(gene_unspliced_other + gene_spliced_other)*100
  
  # record calculation
  splice_myofibre_cells[[sample_id]]  = percent_spliced_cell
  splice_myofibre_gene_other[[sample_id]]  = percent_spliced_gene_other
  splice_myofibre_gene_other_myofibre[[sample_id]]  = percent_spliced_gene_other_myofibre
  splice_myofibre_gene_metabolic[[sample_id]]  = percent_spliced_gene_metabolic
  splice_myofibre_gene_erythrocyte[[sample_id]]  = percent_spliced_rbc
}


### annotate merged_srat with splice info
splice_myofibre_cells <- do.call(c, splice_myofibre_cells)
names(splice_myofibre_cells) <- sapply(strsplit( names(splice_myofibre_cells), "\\."), function(x)x[[2]])
splice_myofibre_cells <- splice_myofibre_cells[names(splice_myofibre_cells) %in% colnames(merged_srat)]
splice_myofibre_cells <- splice_myofibre_cells[match(colnames(merged_srat), names(splice_myofibre_cells))]
names(splice_myofibre_cells) <- colnames(merged_srat)
merged_srat$percent.spliced <- splice_myofibre_cells
percent_spliced_df <- data.frame(celltype = merged_srat$named_celltypes, 
                                 percent_spliced = merged_srat$percent.spliced, 
                                 nFeature_RNA = merged_srat$nFeature_RNA)
VlnPlot(merged_srat, "percent.spliced") +  guides(colour="none")
ggsave("")
VlnPlot(merged_srat, "nFeature_RNA") +  guides(colour="none")



### Compare the %splicing for metabolism associated genes 
# in metabolic muscle and other cells
splice_myofibre_gene_other_df <- do.call(cbind,splice_myofibre_gene_other)
splice_myofibre_gene_other_myofibre_df <- do.call(cbind,splice_myofibre_gene_other_myofibre)
splice_myofibre_gene_metabolic_df <- do.call(cbind,splice_myofibre_gene_metabolic)
splice_myofibre_gene_erythrocyte_df <- do.call(cbind,splice_myofibre_gene_erythrocyte)


avg_gene_other <- rowMeans(splice_myofibre_gene_other_df, na.rm=TRUE)
avg_gene_other_myofibre <- rowMeans(splice_myofibre_gene_other_myofibre_df, na.rm=TRUE)
avg_gene_metabolic <- rowMeans(splice_myofibre_gene_metabolic_df, na.rm=TRUE)
avg_gene_erythrocyte <- rowMeans(splice_myofibre_gene_erythrocyte_df, na.rm=TRUE)

genes_metaup_metacells <- avg_gene_metabolic[names(avg_gene_metabolic) %in% metabolic_both$gene_id]
genes_other_metacells <- avg_gene_metabolic[!names(avg_gene_metabolic) %in% metabolic_both$gene_id]

genes_metaup_othercells <- avg_gene_other[names(avg_gene_other) %in% metabolic_both$gene_id]
genes_other_othercells <- avg_gene_other[!names(avg_gene_other) %in% metabolic_both$gene_id]

genes_metaup_othermyofibrecells <- avg_gene_other_myofibre[names(avg_gene_other_myofibre) %in% metabolic_both$gene_id]
genes_other_othermyofibrecells <- avg_gene_other_myofibre[!names(avg_gene_other_myofibre) %in% metabolic_both$gene_id]

plot(genes_metaup_metacells~genes_metaup_othercells)
plot(genes_other_metacells~genes_other_othercells)
plot(genes_metaup_othermyofibrecells~genes_metaup_othercells)
plot(genes_other_othermyofibrecells~genes_other_othercells)


### plots to make
# Correlation at gene level between muscle, metabolic (MC-M) and non-MC-M
merged_srat$metabolic_muscle <- ifelse(colnames(merged_srat) %in% metabolic_cellid, "metabolic_muscle", "other")
all_gene_avg_exprs <- avg_exp_metabolic_muscle <- AverageExpression(merged_srat, group.by = "metabolic_muscle")
all_gene_avg_exprs <- as.data.frame(all_gene_avg_exprs$RNA)
all_gene_avg_exprs <- all_gene_avg_exprs[all_gene_avg_exprs$`metabolic-muscle`>0,]
all_gene_avg_exprs <- all_gene_avg_exprs[order(all_gene_avg_exprs$`metabolic-muscle`),]
all_gene_avg_exprs$`metabolic-muscle` <- log2(all_gene_avg_exprs$`metabolic-muscle`)
avg_gene_other <- avg_gene_other[row.names(all_gene_avg_exprs)]
avg_gene_metabolic <- avg_gene_metabolic[row.names(all_gene_avg_exprs)]
avg_gene_other_myofibre <- avg_gene_other_myofibre[row.names(all_gene_avg_exprs)]
avg_gene_erythrocyte <- avg_gene_erythrocyte[row.names(all_gene_avg_exprs)]

all_genes_spliced <- data.frame(gene = names(avg_gene_other), 
                                other = avg_gene_other, 
                                other_myofibre = avg_gene_other_myofibre, 
                                erythrocyte = avg_gene_erythrocyte,
                                MC_M = avg_gene_metabolic, 
                                other_exprs = all_gene_avg_exprs$other, MC_M_exprs = all_gene_avg_exprs$`metabolic-muscle`)
all_genes_spliced$gene_type <- ifelse(all_genes_spliced$gene %in% metabolic_both$gene_id,
                                      "fg-over-expressed genes", "other genes") 

ggplot(all_genes_spliced, aes(x=other, y=other_myofibre, colour = MC_M_exprs)) + 
  geom_point(size= 0.5, alpha = 0.6) + theme_bw() + scale_colour_viridis(direction = -1) +
  facet_wrap(~gene_type)+ geom_abline(slope=1, intercept = 0) +
  xlab("Other cell types") + ylab("myofibre")

ggplot(all_genes_spliced, aes(x=other, y=MC_M, colour = MC_M_exprs)) + 
  geom_point(size= 0.5, alpha = 0.6) + theme_bw() + scale_colour_viridis(direction = -1) +
  facet_wrap(~gene_type) +geom_abline(slope=1, intercept = 0)+
  xlab("Other cell types") + ylab("myofibre(fg)")

ggplot(all_genes_spliced, aes(x=other, y=erythrocyte, colour = MC_M_exprs)) + 
  geom_point(size= 0.5, alpha = 0.6) + theme_bw() + scale_colour_viridis(direction = -1) +
  facet_wrap(~gene_type) + geom_abline(slope=1, intercept = 0)+
  xlab("Other cell types") + ylab("erythrocyte")

# stacked bar chart showing % splicing across all genes, or across genes upregulated/not-upregulated in MC-M



# Splice check
plot_df <- data.frame(named_celltypes = merged_srat$named_celltypes, percent.spliced = merged_srat$percent.spliced)
plot_df <- aggregate(plot_df, percent.spliced~named_celltypes, mean, na.rm=TRUE)
ggplot(plot_df, aes(y= named_celltypes, x = percent.spliced)) + geom_point() + ylab("")


## Plot ncount_RNA
VlnPlot(merged_srat, "nFeature_RNA") +  guides(colour="none")
