#### make a background gene list
metascape_bg <- annotated_geneids_df[annotated_geneids_df$gene_id %in% row.names(merged_srat),]
metascape_bg <- metascape_bg$human.genename
metascape_bg <- metascape_bg[!is.na(metascape_bg)]
metascape_bg <- unique(metascape_bg)

#### Erythrocyte
# inflammatory phase
# Idents(merged_srat) <- "named_celltypes"
# erythocyte_markers_set1 <- FindMarkers(merged_srat, ident.1=c("erythrocyte, immune-associated"), 
#                                        ident.2 = c("erythrocyte"),
#                                        only.pos = TRUE)
# erythocyte_markers_set1 <- merge(erythocyte_markers_set1, annotated_geneids_df, by.x=0, by.y="gene_id")
# erythocyte_markers_set1 <- erythocyte_markers_set1[order(erythocyte_markers_set1$p_val_adj),]
# erythocyte_markers_set1 <- erythocyte_markers_set1[erythocyte_markers_set1$avg_log2FC>0,]
# erythocyte_markers_set1 <- erythocyte_markers_set1[erythocyte_markers_set1$p_val_adj < 10^-10,]
# write.csv(erythocyte_markers_set1, paste0(out_dir, "/stats/DEG_cellgroup_eryth_inflam.csv"), row.names=FALSE)
# metascape_input <- unique(erythocyte_markers_set1$human.genename)
# metascape_input <- metascape_input[!is.na(metascape_input)]
# out_df <- data.frame(  "bg" = c("_BACKGROUND", metascape_bg),
#                        "celltype" =  c("erythrocyte_inflam", rep("", length(metascape_bg))))
# out_df[2:(length(metascape_input)+1),2] <- metascape_input
# write.table(out_df, paste0(out_dir, "/stats/Metascape_eryth_inflam.csv"), 
#             row.names=FALSE, sep=",", col.names=FALSE,  quote=FALSE)
# # 
# ### Macrophage
# Idents(merged_srat) <- "named_celltypes"
# all_mac <- list(
#                 mac_m1 = "macrophage, runx3+",
#                 mac_m2 = "macrophage, arg2+",
#                 mac_spp1 = "macrophage, spp1+")
# all_mac_vect <- do.call(c, all_mac)
# out_metascape_df <- data.frame("bg" = c("_BACKGROUND", metascape_bg))
# for(idx in 1:length(names(all_mac))){
#   current_name <- names(all_mac)[idx]
#   current_cluster <- all_mac[[idx]]
#   out_deg <-  paste0(out_dir, "/stats/DEG_cellgroup_", current_name, ".csv")
#   out_metascape_df[[(idx + 1)]] =  c(current_name, rep("", length(metascape_bg)))
#   
#   # analyse deg
#   deg <- FindMarkers(merged_srat, ident.1=current_cluster,
#                      ident.2 = all_mac_vect[!all_mac_vect %in% current_cluster],
#                      only.pos = TRUE)
#   deg <- merge(deg, annotated_geneids_df, by.x=0, by.y="gene_id")
#   deg <- deg[order(deg$p_val_adj),]
#   deg <- deg[deg$avg_log2FC>0,]
#   deg <- deg[deg$p_val_adj < 10^-10,]
#   write.csv(deg, out_deg, row.names=FALSE)
#   
#   # organise metascape input
#   metascape_input <- unique(deg$human.genename)
#   metascape_input <- metascape_input[!is.na(metascape_input)]
#   out_metascape_df[2:(length(metascape_input)+1),(idx + 1)] <- metascape_input
# }
# 
# # organise metascape input
# write.table(out_metascape_df,
#             paste0(out_dir, "/stats/Metascape_mac.csv"), row.names=FALSE,
#             col.names=FALSE, sep=",", quote=FALSE)




#### Fibroblast-like cells
Idents(merged_srat) <- "named_celltypes"
all_fib_like <- c(pro4 = "progenitor, fibroadipogenic", 
                  fib1 = "fibroblast, tenocyte", 
                  fib2 = "fibroblast", 
                  fib3 = "fibroblast, wound-associated")
out_metascape_df <- data.frame("bg" = c("_BACKGROUND", metascape_bg))
for(idx in 1:length(all_fib_like)){
  current_name <- names(all_fib_like)[idx] 
  current_cluster <- all_fib_like[idx]
  out_deg <-  paste0(out_dir, "/stats/DEG_cellgroup_", current_name, ".csv")
  out_metascape_df[[idx+1]] =  c(current_name, rep("", length(metascape_bg)))
  
  # analyse deg
  deg <- FindMarkers(merged_srat, ident.1=c(current_cluster), 
                     ident.2 = all_fib_like[!all_fib_like == current_cluster],
                     only.pos = TRUE)
  deg <- merge(deg, annotated_geneids_df, by.x=0, by.y="gene_id")
  deg <- deg[order(deg$p_val_adj),]
  deg <- deg[deg$avg_log2FC>0,]
  deg <- deg[deg$p_val_adj < 10^-10,]
  write.csv(deg, out_deg, row.names=FALSE)
  
  # organise metascape input
  metascape_input <- unique(deg$human.genename)
  metascape_input <- metascape_input[!is.na(metascape_input)]
  out_metascape_df[2:(length(metascape_input)+1),(idx + 1)] <- metascape_input
}

# organise metascape input
write.table(out_metascape_df, 
            paste0(out_dir, "/stats/Metascape_fibroblast-like.csv"), row.names=FALSE,
            col.names=FALSE, sep=",", quote=FALSE)


#### Muscle (metascape enrichment - using only GO terms)
## red vs white
rvw_metabolic <- FindMarkers(merged_srat, ident.1=c("myofibre, red, fg" ), 
                             ident.2 = c("myofibre, white, fg" ))
rvw_structural <- FindMarkers(merged_srat, ident.1=c("myofibre, red"), 
                              ident.2 = c("myofibre, white, early", 
                                          "myofibre, white, late"))
# find markers that come up in both lists
metabolic_red <-  rvw_metabolic[(rvw_metabolic$avg_log2FC>0) & (rvw_metabolic$p_val_adj < 10^-10),]
metabolic_white <-  rvw_metabolic[(rvw_metabolic$avg_log2FC<0) & (rvw_metabolic$p_val_adj < 10^-10),]
structural_red <-  rvw_structural[(rvw_structural$avg_log2FC>0) & (rvw_structural$p_val_adj < 10^-10),]
structural_white <-  rvw_structural[(rvw_structural$avg_log2FC<0) & (rvw_structural$p_val_adj < 10^-10),]

colnames(metabolic_white) <- paste0("metabolic.", colnames(metabolic_white))
colnames(structural_white) <- paste0("structural.", colnames(structural_white))
colnames(metabolic_red) <- paste0("metabolic.", colnames(metabolic_red))
colnames(structural_red) <- paste0("structural.", colnames(structural_red))

# merge to find genes in both red and white comparison
red_both <- merge(metabolic_red, structural_red, by=0)
white_both <- merge(metabolic_white, structural_white, by=0)

red_both$avg_log2FC <- (red_both$structural.avg_log2FC + red_both$metabolic.avg_log2FC)/2
white_both$avg_log2FC <-(white_both$structural.avg_log2FC + white_both$metabolic.avg_log2FC)/2

red_both$p_val_adj <- (red_both$structural.p_val_adj + red_both$metabolic.p_val_adj)/2
white_both$p_val_adj <- (white_both$structural.p_val_adj + white_both$metabolic.p_val_adj)/2

red_both$pct_diff <- ((red_both$structural.pct.1-red_both$structural.pct.2) + (red_both$metabolic.pct.1-red_both$metabolic.pct.2))/2
white_both$pct_diff <- ((white_both$structural.pct.1-white_both$structural.pct.2) + (white_both$metabolic.pct.1-white_both$metabolic.pct.2))/2

colnames(red_both)[1] <- "gene_id"
colnames(white_both)[1] <- "gene_id"

# annotate gene name and reorder rows
red_both <- merge(red_both, annotated_geneids_df, by="gene_id")
red_both <- red_both[order(red_both$p_val_adj, red_both$avg_log2FC, red_both$pct_diff),]

white_both <- merge(white_both, annotated_geneids_df, by="gene_id")
white_both <- white_both[order(white_both$p_val_adj, -white_both$avg_log2FC, -white_both$pct_diff),]

write.csv(red_both, paste0(out_dir, "/stats/DEG_cellgroup_muscle_red.csv"), row.names=FALSE)
write.csv(white_both, paste0(out_dir, "/stats/DEG_cellgroup_muscle_white.csv"), row.names=FALSE)

metascape_input_white <- unique(white_both$human.genename)
metascape_input_white <- metascape_input_white[!is.na(metascape_input_white)]
metascape_input_red <- unique(red_both$human.genename)
metascape_input_red <- metascape_input_red[!is.na(metascape_input_red)]

out_df <- data.frame( "bg" = c("_BACKGROUND", metascape_bg),
                      "red" =  c("red", rep("", length(metascape_bg))),
                      "white" =  c("white", rep("", length(metascape_bg))))
out_df[2:(length(metascape_input_red)+1),2] <- metascape_input_red
out_df[2:(length(metascape_input_white)+1),3] <- metascape_input_white
write.table(out_df, paste0(out_dir, "/stats/Metascape_muscle_red_vs_white.csv"), 
            quote=TRUE,row.names=FALSE, col.names=FALSE, sep=",")


## metabolic vs structural
mvs_white <- FindMarkers(merged_srat, ident.1=c("myofibre, white, fg" ), 
                         ident.2 = c("myofibre, white, early", "myofibre, white, late"))
mvs_red <- FindMarkers(merged_srat, ident.1=c("myofibre, red, fg" ), 
                       ident.2 = c("myofibre, red"))
# find markers that come up in both lists
metabolic_white <-  mvs_white[(mvs_white$avg_log2FC>0) & (mvs_white$p_val_adj < 10^-10),]
structural_white <-  mvs_white[(mvs_white$avg_log2FC<0) & (mvs_white$p_val_adj < 10^-10),]
metabolic_red <-  mvs_red[(mvs_red$avg_log2FC>0) & (mvs_red$p_val_adj < 10^-10),]
structural_red <-  mvs_red[(mvs_red$avg_log2FC<0) & (mvs_red$p_val_adj < 10^-10),]

colnames(metabolic_white) <- paste0("white.", colnames(metabolic_white))
colnames(structural_white) <- paste0("white.", colnames(structural_white))
colnames(metabolic_red) <- paste0("red.", colnames(metabolic_red))
colnames(structural_red) <- paste0("red.", colnames(structural_red))

# merge to find genes in both red and white comparison
metabolic_both <- merge(metabolic_white, metabolic_red, by=0)
structural_both <- merge(structural_white, structural_red, by=0)

metabolic_both$avg_log2FC <- (metabolic_both$white.avg_log2FC + metabolic_both$red.avg_log2FC)/2
structural_both$avg_log2FC <-(structural_both$white.avg_log2FC + structural_both$red.avg_log2FC)/2

metabolic_both$p_val_adj <- (metabolic_both$white.p_val_adj + metabolic_both$red.p_val_adj)/2
structural_both$p_val_adj <- (structural_both$white.p_val_adj + structural_both$red.p_val_adj)/2

metabolic_both$pct_diff <- ((metabolic_both$white.pct.1-metabolic_both$white.pct.2) + (metabolic_both$red.pct.1-metabolic_both$red.pct.2))/2
structural_both$pct_diff <- ((structural_both$white.pct.1-structural_both$white.pct.2) + (structural_both$red.pct.1-structural_both$red.pct.2))/2

colnames(metabolic_both)[1] <- "gene_id"
colnames(structural_both)[1] <- "gene_id"

# annotate gene name and reorder rows
metabolic_both <- merge(metabolic_both, annotated_geneids_df, by="gene_id")
metabolic_both <- metabolic_both[order(metabolic_both$p_val_adj, metabolic_both$avg_log2FC, metabolic_both$pct_diff),]

structural_both <- merge(structural_both, annotated_geneids_df, by="gene_id")
structural_both <- structural_both[order(structural_both$p_val_adj, -structural_both$avg_log2FC, -structural_both$pct_diff),]

write.csv(metabolic_both, paste0(out_dir, "/stats/DEG_cellgroup_muscle_fg.csv"), row.names=FALSE)
write.csv(structural_both, paste0(out_dir, "/stats/DEG_cellgroup_muscle_non-fg.csv"), row.names=FALSE)

metascape_input_structural <- unique(structural_both$human.genename)
metascape_input_structural <- metascape_input_structural[!is.na(metascape_input_structural)]
metascape_input_metabolic <- unique(metabolic_both$human.genename)
metascape_input_metabolic <- metascape_input_metabolic[!is.na(metascape_input_metabolic)]

out_df <- data.frame( "bg" = c("_BACKGROUND", metascape_bg),
                      "structural" =  c("structural", rep("", length(metascape_bg))),
                      "metabolic" =  c("metabolic", rep("", length(metascape_bg))))
out_df[2:(length(metascape_input_structural)+1),2] <- metascape_input_structural
out_df[2:(length(metascape_input_metabolic)+1),3] <- metascape_input_metabolic
write.table(out_df, paste0(out_dir, "/stats/Metascape_muscle_non-fg_vs_fg.csv"), 
            quote=TRUE,row.names=FALSE, col.names=FALSE, sep=",")


#### SFC 
Idents(merged_srat) <- "named_celltypes"
all_sfc <- list(sfc = "SFC, osteoblast", 
                sfc_focus = "SFC, osteoblast-like, focus", 
                sfc_fam20c = "SFC, osteoblast-like, neurovascular")
all_sfc_vect <- do.call(c, all_sfc)
out_metascape_df <- data.frame("bg" = c("_BACKGROUND", metascape_bg))
for(idx in 1:length(names(all_sfc))){
  current_name <- names(all_sfc)[idx]
  current_cluster <- all_sfc[[idx]]
  out_deg <-  paste0(out_dir, "/stats/DEG_cellgroup_", current_name, ".csv")
  out_metascape_df[[(idx + 1)]] =  c(current_name, rep("", length(metascape_bg)))
  
  # analyse deg
  deg <- FindMarkers(merged_srat, ident.1=current_cluster, 
                     ident.2 = all_sfc_vect[!all_sfc_vect %in% current_cluster],
                     only.pos = TRUE)
  deg <- merge(deg, annotated_geneids_df, by.x=0, by.y="gene_id")
  deg <- deg[order(deg$p_val_adj),]
  deg <- deg[deg$avg_log2FC>0,]
  deg <- deg[deg$p_val_adj < 10^-10,]
  write.csv(deg, out_deg, row.names=FALSE)
  
  # organise metascape input
  metascape_input <- unique(deg$human.genename)
  metascape_input <- metascape_input[!is.na(metascape_input)]
  out_metascape_df[2:(length(metascape_input)+1),(idx + 1)] <- metascape_input
}

# organise metascape input
write.table(out_metascape_df, 
            paste0(out_dir, "/stats/Metascape_sfc.csv"), row.names=FALSE,
            col.names=FALSE, sep=",", quote=FALSE)




