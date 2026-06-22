## Figure 6 - muscle cells,
# organise output folder
current_out_dir <- paste0(out_plots_manuscript, "/fig6_muscle/")
dir.create(current_out_dir, showWarnings = FALSE)



## dimplot showing only muscle cell types 
all_muscle_celltypes <- unique(as.character(merged_srat$named_celltypes))
all_muscle_celltypes <- all_muscle_celltypes[grepl("myo", all_muscle_celltypes)]
all_muscle_celltypes <- c(all_muscle_celltypes, "progenitor, pax7+")

all_muscle_celltypes <- c("myofibre, white, fg",
                          "myofibre, red, fg",
                          "myofibre, white, late",
                          "myofibre, white, early",
                          "myofibre, red",
                          "myotube",
                          "myoblast",
                          "progenitor, pax7+")
srat_muscle <- subset(merged_srat, named_celltypes %in% all_muscle_celltypes)
srat_muscle$named_celltypes <- factor(as.character(srat_muscle$named_celltypes), levels =all_muscle_celltypes )


# load colour
clusterid_to_celltype_col2 <-  readxl::read_excel(cell_cluster_col_label, sheet="full")
celltype_col2 <- clusterid_to_celltype_col2$colour
names(celltype_col2) <- clusterid_to_celltype_col2$cell_type

# plot
p <- DimPlot(srat_muscle, group.by="named_celltypes", reduction = "umap.integrated.harmony", label=TRUE)
p <- p + guides(colour="none")+ggtitle("") #+ scale_colour_manual(values=celltype_col2) 
p <- p + xlim(c(-15, 0)) + ylim(c(-2, 13))+
  xlab("UMAP_1")+ ylab("UMAP_2")
p <- p + scale_colour_manual(values=celltype_col2)
p
ggsave(paste0(current_out_dir, "/umap_named_celltypes.pdf"), width = 4 ,height = 4)


## feature plot for key genes
all_muscle_celltypes <- unique(as.character(merged_srat$named_celltypes))
all_muscle_celltypes <- all_muscle_celltypes[grepl("myo", all_muscle_celltypes)]
all_muscle_celltypes <- c(all_muscle_celltypes, "progenitor, pax7+")

all_muscle_celltypes <- c(#"myofibre, white, fg",
                          #"myofibre, red, fg",
                          "myofibre, white, late",
                          "myofibre, white, early",
                          "myofibre, red",
                          "myotube",
                          "myoblast",
                          "progenitor, pax7+")
srat_muscle <- subset(merged_srat, named_celltypes %in% all_muscle_celltypes)
srat_muscle$named_celltypes <- factor(as.character(srat_muscle$named_celltypes), levels =all_muscle_celltypes )

for(gene in c("pax7", "myf5", "myod1", "myog", "mymk", 
              "smyhc1", "myom2", "ENSSSAG00000076408"
)){
  fun.FeaturePlot(srat_muscle, gene, 
                  file_type = "pdf",
                  conversion_table = annotated_geneids_df, 
                  out_prefix = current_out_dir,
                  width = 2.45, height = 2.8, 
                  point_size = 0.01,
                  xlim = c(-15, -4),
                  ylim =c(-2, 13))
}


#### Make plots for myofibres
# all_myofibre_celltypes <- c("myofibre, white, fg",
#                             "myofibre, red, fg",
#                             "myofibre, white, late",
#                             "myofibre, white, early",
#                             "myofibre, red")
# srat_myofibre <- subset(merged_srat, named_celltypes %in% all_myofibre_celltypes)
# srat_myofibre$named_celltypes <- factor(as.character(srat_myofibre$named_celltypes), levels =all_myofibre_celltypes )
# 
# # create module score with metabolic genes
# source("script/scRNAseq/table/table_deg_within_broadcellgroup.R")
# 
# srat_myofibre <- AddModuleScore(
#   object = srat_myofibre,
#   features = list(metabolic_both$gene_id),
#   name = 'MCM.score'
# )
# srat_myofibre <- AddModuleScore(
#   object = srat_myofibre,
#   features = list(structural_both$gene_id),
#   name = 'MCS.score'
# )
# 

# # organise genes
# muscle_markers_to_plot <- c(
#   "ENSSSAG00000121490" = "cmya5",
#   "ENSSSAG00000043530" =  "spegb",
#   "ENSSSAG00000053510" = "obsl1b",
#   "ENSSSAG00000067860"  = "ldb3a",
#   "ENSSSAG00000074947" = "nexn",
#   "ENSSSAG00000066053" = "rpl30",
#   "ENSSSAG00000073156" = "rpl8",
#   "ENSSSAG00000065665" = "rps25",
#   "ENSSSAG00000075922" = "ckma",
#   "ENSSSAG00000099398" = "gapdh",
#   "ENSSSAG00000081461" = "pgam2"
# )
# 
# plot_geneid <- fun.find_matching_features(srat_myofibre, names(muscle_markers_to_plot), 
#                                           conversion_table = annotated_geneids_df,
#                                           remove.dash=FALSE)
# Idents(srat_myofibre) <- "named_celltypes"
# p <- DotPlot(object = srat_myofibre, features = names(muscle_markers_to_plot), scale.max=70)
# p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# p <- p + scale_x_discrete(labels= muscle_markers_to_plot) + ylab("") + xlab("")
# p <- p + scale_colour_viridis(direction = -1)
# p
# ggsave(paste0(current_out_dir, "/dotplot_myofibre.pdf"), width=7.488, height=2.5)
# 
# 
# #### Plot the same genes in zebrafish
# source("script/scRNAseq/external_data_zebrafish_skin_atlas_2_reanalyse.R")
# muscle_markers_to_plot_zebra <- c(
#   "ENSDARG00000061379" = "cmya5",
#   "ENSDARG00000009567" =  "spegb",
#   "ENSDARG00000077388" = "obsl1b",
#   "ENSDARG00000056322"  = "ldb3a",
#   "ENSDARG00000057317" = "nexn",
#   "ENSDARG00000035871" = "rpl30",
#   "ENSDARG00000014867" = "rpl8",
#   "ENSDARG00000041811" = "rps25",
#   "ENSDARG00000035327" = "ckma",
#   "ENSDARG00000043457" = "gapdh",
#   "ENSDARG00000020007" = "pgam2")
# Idents(z_srat_muscle) <- "cluster"
# # matched <- orth_zeb2sal[orth_zeb2sal$zebrafish.genename %in% muscle_markers_to_plot,]
# # matched <- matched[matched$zebrafish.geneid %in% row.names(z_srat),]
# # matched$salmon.geneid = NULL
# # matched <- unique(matched)
# p <- DotPlot(object = z_srat_muscle, features = names(muscle_markers_to_plot_zebra), scale.max=70)
# p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# p <- p + scale_x_discrete(labels= muscle_markers_to_plot_zebra) + ylab("") + xlab("")
# p <- p + scale_colour_viridis(direction = -1)
# p
# ggsave(paste0(current_out_dir, "/dotplot_myofibre_zebrafish.pdf"), width=5.329, height=2.2)
# 
# 
# 
# #### Plot the same gene in human
# source("script/scRNAseq/external_data_human_muscle.R")
# muscle_markers_to_plot_human <- c(
#   "cmya5",
#   "speg",
#   "obsl1",
#   "ldb3",
#   "nexn",
#   "rpl30",
#   "rpl8",
#   "rps25",
#   "ckm",
#   "gapdh",
#   "pgam2")
# h_srat_muscle$annotation_level2 <- as.character(h_srat_muscle$annotation_level2)
# h_srat_muscle$annotation_level2 <- factor(h_srat_muscle$annotation_level2, 
#                                           levels = c("MF-IIsn(fg)",
#                                                      "MF-Isn(fg)",
#                                                      "MF-II",
#                                                      "MF-I"))
# Idents(h_srat_muscle) <- "annotation_level2"
# # matched <- orth_zeb2sal[orth_zeb2sal$zebrafish.genename %in% muscle_markers_to_plot,]
# # matched <- matched[matched$zebrafish.geneid %in% row.names(z_srat),]
# # matched$salmon.geneid = NULL
# # matched <- unique(matched)
# matched <- toupper(muscle_markers_to_plot_human)
# matched <- matched[matched %in% row.names(h_srat_muscle)]
# p <- DotPlot(object = h_srat_muscle, features = matched, scale.max=70)
# p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# p <- p  + ylab("") + xlab("")
# p <- p + scale_colour_viridis(direction = -1)
# p
# ggsave(paste0(current_out_dir, "/dotplot_myofibre_human.pdf"), width=6, height=2)
# 

