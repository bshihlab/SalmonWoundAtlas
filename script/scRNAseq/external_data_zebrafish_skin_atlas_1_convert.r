### Analyse other existing data
working_dir <- "L:/ShihLab/projects/2024/20240226_rose_wounds"
setwd(working_dir)

## Zebrafish skin atlas 
library(monocle3)
#obj <- readRDS("data/other_snRNAseq_data/GSM7029635_zskin_wildtype_cds.RDS")
obj <- readRDS("data/other_snRNAseq_data/GSM7029635_zskin_all_genotypes_cds")
plot_cells(obj, color_cells_by="cell_type_sub")
plot_cells(obj, genes="ENSDARG00000035895")


# Export data
count_mx <- exprs(obj)
gene_metadata <- fData(obj)
cell_metadata <- colData(obj)

saveRDS(count_mx, "data/other_snRNAseq_data/GSM7029635_zskin_all_genotypes_cds_count_mx.RDS")
write.csv(gene_metadata, "data/other_snRNAseq_data/GSM7029635_zskin_all_genotypes_cds_gene.csv", row.names=FALSE)
write.csv(cell_metadata, "data/other_snRNAseq_data/GSM7029635_zskin_all_genotypes_cds_cellmetadata.csv", row.names=FALSE)