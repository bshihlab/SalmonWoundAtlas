### Figure S3, undefined cell types
# organise output folder
current_out_dir <- paste0(out_plots_manuscript, "/figS3_undefined/")
dir.create(current_out_dir, showWarnings = FALSE)


### annotate clusters near transiting cells
### Plot ribosome content
## Other qc
for(gene in c("sema3d", "slit1", "aspn", "krt15", "clu", "rspo1", "rpl8","Rb.score1")){
  fun.FeaturePlot(merged_srat, gene, 
                  file_type = "pdf",
                  conversion_table = annotated_geneids_df, 
                  out_prefix = current_out_dir,
                  width = 6, height = 5, point_size = 0.01)
}
