# summarise cell count changes
all_cells <- data.frame(cellid = colnames(merged_srat), celltype = merged_srat$named_celltypes)
summarised_celltype <- aggregate(data= all_cells, cellid~celltype, FUN=length)

write.csv(summarised_celltype, paste0(out_tables, "/celltype_count.csv"), row.names=FALSE)

# erythocyte count per treatment
all_cells <- data.frame(cellid = colnames(merged_srat), celltype = merged_srat$named_celltypes, treatment= merged_srat$treatment_group)
summarised_treat_count <- aggregate(data= all_cells, cellid~treatment, FUN=length)
summarised_treat_celltype <- aggregate(data= all_cells, cellid~treatment + celltype, FUN=length)
summarised_treat_erythocyte <- summarised_treat_celltype[summarised_treat_celltype$celltype == "erythrocyte",]
summarised_treat_erythocyte$celltype=NULL
colnames(summarised_treat_erythocyte)[2] <- "erythrocyte"
summarised_treat_thrombocyte <- summarised_treat_celltype[summarised_treat_celltype$celltype == "thrombocyte",]
summarised_treat_thrombocyte$celltype=NULL
colnames(summarised_treat_thrombocyte)[2] <- "thrombocyte"

summarised_treat_count <- merge(summarised_treat_count, summarised_treat_erythocyte, by="treatment")
summarised_treat_count <- merge(summarised_treat_count, summarised_treat_thrombocyte, by="treatment")
summarised_treat_count$percent.erythrocyte <- summarised_treat_count$erythrocyte/summarised_treat_count$cellid*100
summarised_treat_count$percent.thrombocyte <- summarised_treat_count$thrombocyte/summarised_treat_count$cellid*100
write.csv(summarised_treat_count, paste0(out_tables, "/percent_erythrocyte_thrombocyte.csv"), row.names=FALSE)

