all_10x_passed_barcodes <- list()
for(sample in levels(merged_srat_backup$orig.ident)){
  current_10x_barcode <- read.table(gzfile(paste0("L:/ShihLab/projects/2024/20240226_rose_wounds/analysis/cellranger/", sample, "/outs/filtered_feature_bc_matrix/barcodes.tsv.gz")))
  current_10x_barcode <- strsplit(current_10x_barcode$V1, "-")
  current_10x_barcode <- sapply(current_10x_barcode, function(x)x[1])
  current_10x_barcode <- paste0(sample, "_", current_10x_barcode)
  all_10x_passed_barcodes[[sample]] = current_10x_barcode
}

all_10x_passed_barcodes <- do.call(c, all_10x_passed_barcodes)

# check against the existing filtered cells to see if any are missing from 10x passed barcode
pre_ribo_rm <- merged_srat_backup@meta.data
pre_ribo_rm_missing_10x <- pre_ribo_rm[!row.names(pre_ribo_rm) %in% all_10x_passed_barcodes,]

summary(factor(pre_ribo_rm_missing_10x$integrated_clusters))
summary(pre_ribo_rm_missing_10x$orig.ident)