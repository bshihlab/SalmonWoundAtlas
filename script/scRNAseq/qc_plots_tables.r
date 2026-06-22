# summarise cell count changes
Cell_count_info_df <- do.call(rbind,Cell_count_info)
Cell_count_info_title <- c("sample_id",	"cells passed STAR QC",
                           "cells passed Cellranger QC",
                           "QC based on nCount_RNA, nFeature_RNA and % ribosomal content (<40%)",
                           "doublet removed")
names(Cell_count_info_title) <- colnames(Cell_count_info_df)
Cell_count_info_title_df <- as.data.frame(t(as.data.frame(Cell_count_info_title)))
Cell_count_info <- rbind(Cell_count_info_title, Cell_count_info_df)

write.table(Cell_count_info, 
          sep=",",
          paste0(out_dir, "/tables/cell_count_individual.csv"), 
          col.names=FALSE)

metadata  <- merged_srat@meta.data
# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# Rename columns
metadata <- metadata %>%
        dplyr::rename(seq_folder = orig.ident,
                      nUMI = nCount_RNA,
                      nGene = nFeature_RNA)
# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  	ggplot(aes(color=seq_folder , x=nGene, fill= seq_folder )) + 
  	geom_density(alpha = 0.2) + 
  	theme_classic() +
  	scale_x_log10() + 
  	geom_vline(xintercept = 300)

# Visualize the distribution of genes detected per cell via boxplot
metadata %>% 
  	ggplot(aes(x=seq_folder, y=log10(nGene), fill=seq_folder)) + 
  	geom_boxplot() + 
  	theme_classic() +
  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  	theme(plot.title = element_text(hjust=0.5, face="bold")) +
  	ggtitle("NCells vs NGenes")


