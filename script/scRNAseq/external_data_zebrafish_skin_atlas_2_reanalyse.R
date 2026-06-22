# external data zebrafish
working_dir <- "L:/ShihLab/projects/2024/20240226_rose_wounds"
setwd(working_dir)

# load library 
library(Seurat)

# read data
orth_zebsalhum_df <- read.csv("annotation/gene/orth_zebsalhum.csv")

count_mx <- readRDS("data/other_snRNAseq_data/GSM7029635_zskin_all_genotypes_cds_count_mx.RDS")
gene <- read.csv("data/other_snRNAseq_data/GSM7029635_zskin_all_genotypes_cds_gene.csv")
cell <- read.csv("data/other_snRNAseq_data/GSM7029635_zskin_all_genotypes_cds_cellmetadata.csv")

z_srat <- CreateSeuratObject(count_mx, project = "zebrafish_skin_atlas", assay = "RNA",
                                     min.cells = 0, min.features = 0, names.field = 1,
                             names.delim = "_", 
                             meta.data = cell)
z_srat_muscle <- subset(z_srat, cell_group == "Muscle")

## normalise and scale
# QC not performed as it is assumed that QC has already been performed
z_srat_muscle <- NormalizeData(z_srat_muscle)
z_srat_muscle <- FindVariableFeatures(z_srat_muscle, selection.method = "vst", nfeatures = 2000)
z_srat_muscle <- ScaleData(z_srat_muscle, features = row.names(z_srat_muscle))
z_srat_muscle <- RunPCA(z_srat_muscle)

