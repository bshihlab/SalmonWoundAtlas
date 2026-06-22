#### Read human muscle data
working_dir <- "L:/ShihLab/projects/2024/20240226_rose_wounds"
setwd(working_dir)

# load library 
library(Seurat)
library(reticulate)
library(anndata)

data <- read_h5ad("data/other_snRNAseq_data/SKM_myonuclei2myofiber_human_2023-06-22.h5ad")
srat_h <- CreateSeuratObject(counts = t(as.matrix(data$X)), meta.data = data$obs,min.features = 500, min.cells = 30)
all_cell_type <- unique(as.character(srat_h$annotation_level2))
all_cell_type_mf <- all_cell_type[grepl("MF-|^I-|^II-", all_cell_type)]
all_cell_type_mf <- all_cell_type[grepl("MF-", all_cell_type)]
h_srat_muscle <- subset(srat_h, annotation_level2 %in% all_cell_type_mf)

# QC not performed as it is assumed that QC has already been performed
h_srat_muscle <- NormalizeData(h_srat_muscle)
h_srat_muscle <- FindVariableFeatures(h_srat_muscle, selection.method = "vst", nfeatures = 2000)
h_srat_muscle <- ScaleData(h_srat_muscle, features = row.names(h_srat_muscle))
#h_srat_muscle <- RunPCA(h_srat_muscle)
