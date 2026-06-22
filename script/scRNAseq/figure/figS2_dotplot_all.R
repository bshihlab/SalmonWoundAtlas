### Figure S2
# organise output folder
current_out_dir <- paste0(out_plots_manuscript, "/figS2_all_celltype_dotplot/")
dir.create(current_out_dir, showWarnings = FALSE)

gene_plot_order <- read.csv("annotation/gene/supplementary_dotplot_cell_markers.csv")
cell_plot_order <- read.csv("annotation/gene/supplementary_dotplot_cell_markers_celltype_order.csv")

merged_srat_tmp <- merged_srat
DefaultAssay(merged_srat_tmp) <- "RNA"
merged_srat_tmp$named_celltypes <- as.character(merged_srat_tmp$named_celltypes)
merged_srat_tmp$named_celltypes <- factor(merged_srat_tmp$named_celltypes, levels=rev(cell_plot_order$named_celltypes))
Idents(merged_srat_tmp) <- "named_celltypes"
p <- DotPlot(object = merged_srat_tmp, features = gene_plot_order$Ensembl, scale.max=60)
p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p <- p + scale_x_discrete(labels= gene_plot_order$Key.genes) + ylab("") + xlab("")
p <- p + scale_colour_viridis(direction = -1)
p
ggsave(paste0(current_out_dir, "/dotplot_all.pdf"), width=22, height=10)

