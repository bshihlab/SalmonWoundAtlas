#### Functions

#### Organise gene annotation
fun.organise_gene_annotation <- function(feature_df, orthologue_fp = "data/orthologue/s_salar_ensembl_106_gene_data.csv"){
  ####-- Organise input
  orthologue_info <- read.csv(orthologue_fp)
  orthologue_info <- orthologue_info[,c("ensembl.id" , "gene.name", "external.gene.name")]
  colnames(orthologue_info) <- c("gene_id" , "gene.name", "orthologue")
  # annotate the geneid_full dataframe into a new df (geneid_full needs to stay in the same order)
  gene_ids_df <- merge(geneid_full, orthologue_info, by="gene_id", all.x=TRUE)
  gene_ids_df$gene_name_base <- sapply(strsplit(gene_ids_df$gene.name, "\\."), function(x)x[[1]][1])
  gene_ids_df$orthologue <- sapply(strsplit(gene_ids_df$orthologue, "\\."), function(x)x[[1]][1])
  return(gene_ids_df)
}

#### find matching genes given a list of genes
## organise matching features
fun.find_matching_features <- function(obj, features, conversion_table, remove.dash){
  out_list <- list(geneid=vector(),
                   outname = vector())
  # organise conversion table so it matches the order of the rows in obj
  gene_ids <- row.names(obj)
  metadata_colnames <- colnames(obj@meta.data)
  conversion_table <- conversion_table[conversion_table$gene_id %in% gene_ids,]
  conversion_table <- conversion_table[match(gene_ids, conversion_table$gene_id),]
  #print("a")
  
  # convert all to uppercase to avoid case mismatch
  features <- unique(features)
  # match each feature and determine the final output gene name
  for(feature in features){
    # if the "feature" is in metadata column names, just print out the feature as it is
    if(feature %in% metadata_colnames){
      out_list[["geneid"]] <- c(out_list[["geneid"]], feature)
      out_list[["outname"]] <- c(out_list[["outname"]], feature)
      
      # if the feature is not a metadata column
    } else {
      feature <- toupper(feature)
      
      # find matching index
      if(remove.dash){
        feature_short <- strsplit(feature, "-")[[1]][1]
        feature_short <- strsplit(feature, "\\.")[[1]][1]
      }else{
        feature_short <- strsplit(feature, "\\.")[[1]][1]
      }
      match_rowname <- which(toupper(conversion_table$gene_id) %in% feature_short)
      #match_genename <- match(feature_short, toupper(conversion_table$gene.name) %in% feature_short)
      match_genename_base <- which(toupper(conversion_table$gene_name_base) %in% feature_short)
      match_human.genename <- which(toupper(conversion_table$human.genename) %in% feature_short)
      match_zebrafish.genename <- which(toupper(conversion_table$zebrafish.genename) %in% feature_short)
      matched_idx <- unique(c(match_rowname, match_human.genename, match_zebrafish.genename, match_genename_base))
      matched_idx <- matched_idx[!is.na(matched_idx)]
      matched_idx <- unique(matched_idx)
      #print("b")
      #print(matched_idx)
      
      # print out warning on genes without match
      
      
      #print(paste0("matched_idx", matched_idx))
      # store the full list of ensembl ids/plot name for the current feature 
      current_rowname <- vector()
      current_outname <- vector()
      if(length(matched_idx) < 1){
        print(paste0("Features not found: ", feature))
      } else {
        for(idx in matched_idx){
          tmp_outname_organisation <- c(#conversion_table$orthologue[idx], 
            conversion_table$gene_id[idx],
            conversion_table$gene_name_base[idx],
            conversion_table$zebrafish.genename[idx], 
            conversion_table$human.genename[idx]
          )
          print(tmp_outname_organisation)
          # remove redundancy
          tmp_outname_organisation <- unique(tmp_outname_organisation)
          tmp_outname_organisation <- tmp_outname_organisation[!is.na(tmp_outname_organisation)]
          tmp_outname_organisation <- paste0(tmp_outname_organisation, collapse= "_")
          # store the current rowname and organised name
          #print(row.names(obj)[idx])
          current_rowname <- c(current_rowname, row.names(obj)[idx])
          current_outname <- c(current_outname, tmp_outname_organisation)
        }
      }
      
      # store all the matched feature info in the overall vector
      out_list[["geneid"]] <- c(out_list[["geneid"]], current_rowname)
      out_list[["outname"]] <- c(out_list[["outname"]], current_outname)
    }
  }
  #print(features)
  
  return(out_list)
}


#### scRNAseq: Make indivdiual feature plot given a list of genes
fun.FeaturePlot <- function(obj, features, out_prefix, conversion_table, 
                            remove.dash=FALSE, split.by="",
                            point_size = 0.2, file_type = "png", 
                            reduction = "umap.integrated.harmony", 
                            xlab = "UMAP_1", ylab = "UMAP_2",
                            max.cutoff = "q95",
                            width = 5, height = 4,
                            xlim= NULL, ylim=NULL){
  DefaultAssay(obj) <- "RNA"
  # find matching features
  matched_features <- fun.find_matching_features(obj, features, conversion_table, remove.dash=remove.dash)
  all_matching_geneid <- matched_features[["geneid"]]
  all_matching_outname <- matched_features[["outname"]]
  obj_metadata <- as.data.frame(obj@meta.data)
  
  ## Plot each matching feature. 
  
  if(length(all_matching_geneid)>0){
    DefaultAssay(obj) <- "RNA"
    
    plot_df <- as.data.frame(obj@reductions[[reduction]]@cell.embeddings)
    
    # store the original column name elsewhere and rename the colmns as Dim1 and Dim2
    original_dimnames <- colnames(plot_df)
    colnames(plot_df) <- c("Dim1", "Dim2")
    
    # Add split.by annotation if present in input parameter
    if((!split.by == "") & (split.by %in% colnames(obj_metadata))){
      plot_df[[split.by]] <- obj_metadata[[split.by]]
    }
    
    # Get ene expression data
    exprs_df <-  FetchData(object = obj, assay= "RNA", 
                           layer = "data", vars = all_matching_geneid)
    
    plot_df <- cbind(plot_df, exprs_df)
    # Add in each non-gene/metadata column
    
    obj_metadata <- as.data.frame(obj@meta.data)
    
    if(sum(features %in% colnames(obj_metadata))>0) {
      for(idx in 1:length(features)){
        
        if(features[idx] %in% colnames(obj_metadata)){
          plot_df[[features[idx]]] <- obj_metadata[,features[idx]]
          all_matching_geneid[features[idx]] <- features[idx]
        }
      }
    } 
    
    # go through each gene
    # Sort the values out so the positive dots are plotted on the top
    for(idx in 1:length(all_matching_geneid)){
      current_colname <- all_matching_geneid[idx]
      current_outname <- all_matching_outname[idx]
      if(!(split.by == "") & (split.by %in% colnames(plot_df))){
        current_plot_df <- plot_df[, c("Dim1", "Dim2", current_colname, split.by)]
        colnames(current_plot_df) <- c("Dim1", "Dim2", "exprs", "split")
      } else {
        current_plot_df <- plot_df[, c("Dim1", "Dim2", current_colname)]
        colnames(current_plot_df) <- c("Dim1", "Dim2", "exprs")
      }
      # sort so cells with expression are plotted on top
      current_plot_df <- current_plot_df[order(current_plot_df$exprs, decreasing=F),]
      # if q95 is true, turn all the values above q95 to q95
      if(max.cutoff == "q95"){
        q95_val <- quantile(current_plot_df$exprs[current_plot_df$exprs > 0], 0.95)
        current_plot_df$exprs <- ifelse(current_plot_df$exprs > q95_val, q95_val,  current_plot_df$exprs)
      }
      # Organise x and y  label
      dimention_label <- strsplit(reduction, "\\.|_")[[1]][1]
      # plot and save
      out_fp <- paste0(out_prefix, current_outname, ".", file_type)
      p <- ggplot(current_plot_df, aes(x=Dim1, y=Dim2, colour=exprs)) +
        geom_point(size=point_size) + ggtitle(current_outname) +
        xlab(original_dimnames[1]) + ylab(original_dimnames[2]) +
        scale_colour_viridis(option="magma", direction = -1) + theme_bw()
      # do facet_wrap if there is an input for split.by 
      if(!(split.by == "") & (split.by %in% colnames(plot_df))){
        p <- p + facet_wrap(~split)
      }
      p <- p + xlab(xlab) + ylab(ylab)
      # add x and y lim if entered
      if(! is.null(xlim)){
        p <- p + xlim(xlim)
      }
      if(! is.null(ylim)){
        p <- p + ylim(ylim)
      }
      print(p)
      ggsave(out_fp, width = width, height = height)
    }
  } else {
    print(paste0(paste0(features, collapse="."), " not found"))
  }
}


#### scRNAseq: Make individual Dim plot given a list of genes
fun.DimPlot <- function(obj, colour_label, colour_key = clusters_colours,
                        out_prefix, point_size = 0.2, file_type = "png", 
                        reduction = "umap.integrated.harmony", 
                        xlab = "UMAP_1", ylab = "UMAP_2",
                        label.text=TRUE, split_graph=NULL){
  # organise plot data
  obj_metadata <- as.data.frame(obj@meta.data)
  plot_df <- as.data.frame(obj@reductions[[reduction]]@cell.embeddings)
  original_dimnames <- colnames(plot_df)
  plot_df <- cbind(plot_df, obj_metadata[,c(colour_label, "cell_id", "orig.ident","treatment_group")])
  # add the split by column
  if((!is.null(split_graph))){
    if((!(split_graph %in% colnames(plot_df))) & (split_graph %in% colnames(obj_metadata))){
      plot_df <- cbind(plot_df, obj_metadata[,split_graph])
    }
  }
  colnames(plot_df)[1:3] <- c("Dim1", "Dim2", "label")
  plot_df$label <- as.character(plot_df$label)
  # calcuate midpoint for text label
  midpoint_dim1 <- aggregate(data = plot_df, Dim1~label, mean)
  midpoint_dim2 <- aggregate(data = plot_df, Dim2~label, mean)
  midpoint_df <- merge(midpoint_dim1, midpoint_dim2, by="label")
  colnames(midpoint_df) <- c( "label","midpoint_1", "midpoint_2")
  
  # edit colours
  if(sum(unique(plot_df$label) %in% names(colour_key)) == length(unique(plot_df$label))){
    colour_key <- colour_key[unique(plot_df$label)]
  } else {
    unique_labels <- unique(plot_df$label)
    colour_key <- colour_key[1:length(unique_labels)]
    names(colour_key) <- unique_labels[order(unique_labels)]
  }
  # make plots
  p <- ggplot(data = plot_df, aes(x=Dim1, y=Dim2)) + geom_point( size = point_size, aes(colour=label)) + 
    theme_bw() + xlab(original_dimnames[1]) + ylab(original_dimnames[2]) + scale_colour_manual(values=colour_key)
  p <- p + xlab(xlab) + ylab(ylab)
  # add text label
  if(label.text){
    p <- p + geom_text(data=midpoint_df, aes(x=midpoint_1 , y=midpoint_2, label=label)) 
  }
  # split graph
  if(!is.null(split_graph)){
    if(split_graph %in% colnames(plot_df))
      p <- p + facet_wrap(~.data[[split_graph]])
  }
  return(p)
}


#### rotate spatial image
rotate_image <- function(p, rot_angle) {
  gt <- ggplot_gtable(ggplot_build(p))
  panel_idx <- which(gt$layout$name == "panel")
  rot_vp <- viewport(angle = rot_angle)
  gt[["grobs"]][[panel_idx]] <- editGrob(gt[["grobs"]][[panel_idx]], vp = rot_vp)
  p_rot <- ggdraw() + draw_grob(gt)
  
  return(p_rot)
}

#### spatial: Make averaged feature plot given a list of genes (averaged Seurat data value across genes)
fun.spatialPlot_avg <- function(obj, features, 
                                conversion_table,
                                out_prefix, 
                                out_suffix = "",
                                legend_label="",
                                remove.dash = FALSE, 
                                exprs_fill = TRUE,
                                pt.size.factor = 2.2, 
                                file_type = "png", 
                                max.cutoff = "q95",
                                rot_angle=225){
  # automatically calculate pt.size.factor if not set
  #pt.size.factor <- ifelse(is.null(pt.size.factor), (0.012*(1/obj@images$slice1@spot.radius)), pt.size.factor)
  
  DefaultAssay(obj) <- "Spatial"
  # find matching features
  matched_features <- fun.find_matching_features(obj, features, conversion_table, remove.dash=remove.dash)
  all_matching_geneid <- matched_features[["geneid"]]
  all_matching_outname <- matched_features[["outname"]]
  
  # average across all found features
  obj$celltype <- rowMeans(FetchData(obj, vars=all_matching_geneid, layer= "data"))
  out_fp <- paste0(out_prefix, "_",  out_suffix, ".", file_type)
  # plot expression as fill
  if(exprs_fill){
    p <- SpatialFeaturePlot(obj, 
                            features = "celltype", 
                            pt.size.factor = pt.size.factor, 
                            max.cutoff = max.cutoff, 
                            min.cutoff = "q1",
                            image.scale="hires",
                            ncol=1)#, stroke = 0.2)
    p <- p+ scale_fill_viridis_c(option = "magma", direction = -1)
    # plot expression as outlines
  } else {
    p <- SpatialFeaturePlot(obj, 
                            features = "celltype", 
                            pt.size.factor = pt.size.factor, 
                            max.cutoff = max.cutoff, 
                            min.cutoff = "q1",
                            image.scale="hires",
                            ncol=1)#, stroke = 0.2)
  }  
  p <- p + labs(fill=legend_label)
  p <- rotate_image(p, rot_angle = rot_angle)
  print(p)
  ggsave(out_fp, width = 5, height = 5)
  return(p)
}


#### spatial: Make individual feature plot given a list of genes
fun.spatialPlot_individual <- function(obj, 
                                       features, 
                                       conversion_table,
                                       out_prefix,  
                                       out_suffix = "",
                                       legend_title = "",
                                       remove.dash = FALSE,
                                       pt.size.factor = 2.2, 
                                       file_type = "png", 
                                       max.cutoff = "q95",
                                       min.cutoff = "q5",
                                       rot_angle=225){
  
  # automatically calculate pt.size.factor if not set
  #pt.size.factor <- ifelse(is.null(pt.size.factor), (0.017*(1/obj@images$slice1@spot.radius)), pt.size.factor)
  
  DefaultAssay(obj) <- "Spatial"
  # find matching features
  matched_features <- fun.find_matching_features(obj, features, conversion_table, remove.dash=remove.dash)
  all_matching_geneid <- matched_features[["geneid"]]
  all_matching_outname <- matched_features[["outname"]]
  print(paste0("Found", paste0(all_matching_outname, collapse=",")))
  
  # Loop through all matching features and make individual plots
  for(idx in 1:length(all_matching_geneid)){
    current_rowname <- all_matching_geneid[idx]
    current_outname <- all_matching_outname[idx]
    current_outname <- str_replace_all(current_outname, "[[:punct:]]", " ")
    
    out_fp <- paste0(out_prefix, current_outname, "_", out_suffix, ".", file_type)
    p <- SpatialFeaturePlot(obj, 
                            features = current_rowname, 
                            pt.size.factor = pt.size.factor, 
                            max.cutoff = max.cutoff, 
                            min.cutoff = min.cutoff,
                            #alpha = c(0.1, 1), 
                            image.scale="hires",
                            ncol=1)#, stroke=0.1)
    current_outname <- ifelse(legend_title=="", current_outname, legend_title)
    p <- p + labs(fill=current_outname)
    p <- p +  scale_fill_viridis_c(option = "magma")
    p <- rotate_image(p, rot_angle = rot_angle)
    
    # change the background to grey - this helps visualise spots
    #p <- p + theme(plot.background = element_rect(fill = "grey"))
    
    print(p)
    ggsave(out_fp, width = 5, height = 5)
  }
}





#### spatial: Make dural feature plot given a list of genes
fun.spatialPlot_dual <- function(obj1,
                                 obj2, features, 
                                 conversion_table,
                                 out_prefix, 
                                 out_suffix = "",
                                 legend_label="",
                                 remove.dash = FALSE, 
                                 exprs_fill = TRUE,
                                 pt.size.factor = 2.2, 
                                 file_type = "png", 
                                 max.cutoff = "q95",
                                 rot_angle1=225,
                                 rot_angle2=225){
  # automatically calculate pt.size.factor if not set
  #pt.size.factor <- ifelse(is.null(pt.size.factor), (0.012*(1/obj@images$slice1@spot.radius)), pt.size.factor)
  
  DefaultAssay(obj1) <- "Spatial"
  DefaultAssay(obj2) <- "Spatial"
  # find matching features
  matched_features1 <- fun.find_matching_features(obj1, features, conversion_table, remove.dash=remove.dash)
  matched_features2 <- fun.find_matching_features(obj1, features, conversion_table, remove.dash=remove.dash)
  matched_features <- matched_features1[matched_features1 %in% matched_features2]
  all_matching_geneid <- matched_features[["geneid"]]
  all_matching_outname <- matched_features[["outname"]]
  print(paste0("Found", paste0(all_matching_outname, collapse=",")))
  
  # average across all found features
  obj1$celltype <- rowMeans(FetchData(obj1, vars=all_matching_geneid, layer= "data"))
  obj2$celltype <- rowMeans(FetchData(obj2, vars=all_matching_geneid, layer= "data"))
  out_fp <- paste0(out_prefix, "_",  out_suffix, ".", file_type)
  
  # plot fill limits
  # take the q95 as upper limit, and q1 as lower limit.
  # q95 and q1 are calculated across all samples
  all_plot_val <- c(obj1$celltype , obj2$celltype )
  all_q95 <- quantile(all_plot_val, 0.95)[[1]] 
  all_q1 <-  quantile(all_plot_val, 0.01)[[1]]
  
  obj1$celltype <- ifelse(obj1$celltype>all_q95, all_q95, ifelse(obj1$celltype < all_q1, all_q1, obj1$celltype))
  obj2$celltype <- ifelse(obj2$celltype>all_q95, all_q95, ifelse(obj2$celltype<all_q1, all_q1, obj2$celltype ))
  
  # plot expression as fill
  p1 <- SpatialFeaturePlot(obj1, 
                           features = "celltype", 
                           pt.size.factor = pt.size.factor, 
                           #max.cutoff = max.cutoff, 
                           #min.cutoff = "q1",
                           image.scale="hires",
                           ncol=1)#, stroke = 0.2)
  p2 <- SpatialFeaturePlot(obj2, 
                           features = "celltype", 
                           pt.size.factor = pt.size.factor, 
                           #max.cutoff = max.cutoff, 
                           #min.cutoff = "q1",
                           image.scale="hires",
                           ncol=1)#, stroke = 0.2)
  
  p1 <- p1+ scale_fill_distiller(palette  = "Spectral", direction = -1, limits = c(all_q1, all_q95))
  p2 <- p2+ scale_fill_distiller(palette  = "Spectral", direction = -1, limits = c(all_q1, all_q95))
  
  p1 <- p1 + labs(fill=legend_label)
  p2 <- p2 + labs(fill=legend_label)
  p1 <- rotate_image(p1, rot_angle = rot_angle1)
  p2 <- rotate_image(p2, rot_angle = rot_angle2)
  p3 <- p1 + p2
  print(p3)
  ggsave(out_fp, width = 5, height = 5)
  return(p3)
}
