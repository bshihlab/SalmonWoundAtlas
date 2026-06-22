#### print out cell count and proportion for fig1 and full annotation
# fig1 annotation
metadata_df <- merged_srat@meta.data
out_fp <- paste0( out_tables, "/cell_percentage.xlsx")

# Group by fig1 broad cell groups and by samples
aggregated_df <- aggregate(metadata_df,  cell_id ~ orig.ident + celltype_fig1, FUN=length)
aggregated_df <- pivot_wider(aggregated_df, names_from=orig.ident, values_from = cell_id)
celltype_lab <- aggregated_df$celltype_fig1
aggregated_df$celltype_fig1 = NULL
aggregated_df_total<- apply(aggregated_df, 2, FUN=function(x)round(x/sum(x, na.rm=TRUE)*100, 2))
aggregated_df_total_1 <- cbind(celltype = celltype_lab, as.data.frame(aggregated_df_total))
#sheet = createSheet(wb, "broad_per_sample")
#addDataFrame(aggregated_df_total, sheet=sheet, startColumn=1, row.names=FALSE)

# Group by fig1 broad cell groups and by treatment
aggregated_df <- aggregate(metadata_df,  cell_id ~ treatment_group + celltype_fig1, FUN=length)
aggregated_df <- pivot_wider(aggregated_df, names_from=treatment_group, values_from = cell_id)
celltype_lab <- aggregated_df$celltype_fig1
aggregated_df$celltype_fig1 = NULL
aggregated_df_total<- apply(aggregated_df, 2, FUN=function(x)round(x/sum(x, na.rm=TRUE)*100, 2))
aggregated_df_total_2 <- cbind(celltype = celltype_lab, as.data.frame(aggregated_df_total))
#sheet = createSheet(wb, "broad_per_treatment")
#addDataFrame(aggregated_df_total, sheet=sheet, startColumn=1, row.names=FALSE)

# Group by full cell type annotation and by orig.ident
aggregated_df <- aggregate(metadata_df,  cell_id ~ orig.ident + named_celltypes, FUN=length)
aggregated_df <- pivot_wider(aggregated_df, names_from=orig.ident, values_from = cell_id)
celltype_lab <- aggregated_df$named_celltypes
aggregated_df$named_celltypes = NULL
aggregated_df_total<- apply(aggregated_df, 2, FUN=function(x)round(x/sum(x, na.rm=TRUE)*100, 2))
aggregated_df_total_3 <- cbind(celltype = celltype_lab, as.data.frame(aggregated_df_total))
#sheet = createSheet(wb, "specific_per_sample")
#addDataFrame(aggregated_df_total, sheet=sheet, startColumn=1, row.names=FALSE)

# Group by full cell type annotation and by treatment
aggregated_df <- aggregate(metadata_df,  cell_id ~ treatment_group + named_celltypes, FUN=length)
aggregated_df <- pivot_wider(aggregated_df, names_from=treatment_group, values_from = cell_id)
celltype_lab <- aggregated_df$named_celltypes
aggregated_df$named_celltypes = NULL
aggregated_df_total<- apply(aggregated_df, 2, FUN=function(x)round(x/sum(x, na.rm=TRUE)*100, 2))
aggregated_df_total_4 <- cbind(celltype = celltype_lab, as.data.frame(aggregated_df_total))
#sheet = createSheet(wb, "specific_per_treatment")
#addDataFrame(aggregated_df_total, sheet=sheet, startColumn=1, row.names=FALSE)

# save excel
all_sheets <- list(
  broad_per_sample = aggregated_df_total_1,
  broad_per_treatment = aggregated_df_total_2,
  specific_per_sample = aggregated_df_total_3,
  specific_per_treatment = aggregated_df_total_4
)
wb = createWorkbook()
Map(function(df, tab_name){     
  
  addWorksheet(wb, tab_name)
  writeData(wb, tab_name, df)
}, 

all_sheets, names(all_sheets)
)

saveWorkbook(wb, file = out_fp, overwrite = TRUE)