#### find all markers (based on annotated_cell types)
merged_srat<- PrepSCTFindMarkers(merged_srat)
Idents(merged_srat) <- "named_celltypes"
all_cluster_deg_named <- FindAllMarkers(merged_srat, only.pos=TRUE)
#all_cluster_deg_sig <- all_cluster_deg_sig[all_cluster_deg_sig$adj.p.]
out_stats_df <- merge(all_cluster_deg_named, annotated_geneids_df, by.x="gene", by.y="gene_id", all.x=TRUE)
out_stats_df <- out_stats_df[order(out_stats_df$cluster, out_stats_df$p_val_adj),] 
#out_stats_unique <- out_stats_unique[!duplicated(out_stats_unique$gene),]
write.csv(out_stats_df, paste0(out_dir_stats,"/cluster_characterisation_deg.csv"))

#### find "unique" markers

# Find markers that are unique to each cluster
all_markers_df_sig <- out_stats_df[out_stats_df$p_val_adj < 0.0001,] # keep genes with low p value
all_markers_df_sig <- all_markers_df_sig[all_markers_df_sig$pct.1 > 0.10,] # keep genes with at least 15% expressed in pct.1
all_markers_df_sig <- all_markers_df_sig[all_markers_df_sig$avg_log2FC  > 1.25,] # keep genes with at least 1.25 log2FC
all_markers_list <- split(all_markers_df_sig, all_markers_df_sig$cluster)

# keep any genes that are also present in other gene lists
all_markers_list_unique <- lapply(names(all_markers_list), function(x){
  out <- all_markers_list[[x]]
  out <- out[!(out$gene %in% all_markers_df_sig$gene[!all_markers_df_sig$cluster == x]),]
  return(out)
}
)

# keep genes that are in multiple cell types only 
# if the highest differentially expressed cluster has at least 1 fold 
# higher expression than the second cluster
all_markers_list_repeat <- lapply(names(all_markers_list), function(x)
  all_markers_list[[x]][all_markers_list[[x]]$gene %in% all_markers_df_sig$gene[!all_markers_df_sig$cluster == x],] 
)
all_markers_list_repeat <- do.call(rbind, all_markers_list_repeat)
all_markers_list_repeat <- all_markers_list_repeat[order(all_markers_list_repeat$p_val_adj, -all_markers_list_repeat$avg_log2FC),]
all_markers_list_repeat <- split(all_markers_list_repeat, all_markers_list_repeat$gene)
all_markers_list_repeat <- lapply(all_markers_list_repeat, function(x){if(x$avg_log2FC[1] - x$avg_log2FC[2] > log2(1.5)){out <- x[1,]}else{out <- NULL} ; return(out)})

# combine the two lists
all_markers_list_unique <- rbind(do.call(rbind,all_markers_list_unique), do.call(rbind,all_markers_list_repeat))
all_markers_list_unique <- split(all_markers_list_unique, as.character(all_markers_list_unique$cluster))
all_markers_list_unique_out <- do.call(rbind,all_markers_list_unique)
all_markers_list_unique_out <- all_markers_list_unique_out[order(all_markers_list_unique_out$cluster, all_markers_list_unique_out$p_val_adj ), ]

# Make a table with just the top 20 unique genes
all_markers_list_unique_top20 <- lapply(all_markers_list_unique, function(x)if(nrow(x)> 20){x[1:20,]}else{x})
all_markers_list_unique_top20 <- do.call(rbind,all_markers_list_unique_top20)
all_markers_list_unique_top20 <- all_markers_list_unique_top20[order(all_markers_list_unique_top20$cluster, all_markers_list_unique_top20$p_val_adj ), ]
write.csv(all_markers_list_unique_top20, paste0(out_dir_stats,"/cluster_characterisation_deg_top20unique.csv"), row.names=FALSE)


# organise metascape input
#### make a background gene list
metascape_bg <- annotated_geneids_df[annotated_geneids_df$gene_id %in% row.names(merged_srat),]
metascape_bg <- metascape_bg$human.genename
metascape_bg <- metascape_bg[!is.na(metascape_bg)]
metascape_bg <- unique(metascape_bg)

# filter stats df
out_stats_df$pct.diff <- out_stats_df$pct.1 - out_stats_df$pct.2
out_stats_df <- out_stats_df[out_stats_df$pct.diff > 0.2,]
out_stats_df <- out_stats_df[out_stats_df$avg_log2FC > 1,]
out_stats_df <- out_stats_df[out_stats_df$p_val_adj  <  10^-20,]
stats_list <- split(out_stats_df, out_stats_df$cluster) 

# loop through each cluster
# add genes to output 
out_metascape_df <- data.frame("bg" = c("_BACKGROUND", metascape_bg))
for(idx in 1:length(stats_list)){
  current_name <- names(stats_list)[idx]
  current_name <- gsub("[[:punct:]]", "", current_name)  # no libraries needed
  out_metascape_df[,(idx + 1)]<- rep("", (length(metascape_bg) +1))
  metascape_input <- unique(stats_list[[idx]]$human.genename)
  metascape_input <- metascape_input[!is.na(metascape_input)]
  out_metascape_df[1:(length(metascape_input)+1),(idx + 1)] <- c(current_name, metascape_input)
}

out_metascape_df_1 <- out_metascape_df[,1:22]
out_metascape_df_2 <- out_metascape_df[,c(1, 23:ncol(out_metascape_df))]

write.table(out_metascape_df_1, 
            paste0(out_dir_stats, "/Metascape_all_cell_types_part1.csv"), row.names=FALSE,
            col.names=FALSE, sep=",", quote=TRUE)
write.table(out_metascape_df_2, 
            paste0(out_dir_stats, "/Metascape_all_cell_types_part2.csv"), row.names=FALSE,
            col.names=FALSE, sep=",", quote=TRUE)
