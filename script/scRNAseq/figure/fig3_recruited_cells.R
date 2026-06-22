## Figure 3 - sub figures for cell recruitment/granulation tissue composition,
## particularly mesenchymal
# organise output folder
current_out_dir <- paste0(out_plots_manuscript, "/fig3_granulation_tissue/")
dir.create(current_out_dir, showWarnings = FALSE)


for(gene in "wnt5a"){
  fun.FeaturePlot(merged_srat, gene, 
                  file_type = "pdf",
                  conversion_table = annotated_geneids_df, 
                  out_prefix = current_out_dir)
}