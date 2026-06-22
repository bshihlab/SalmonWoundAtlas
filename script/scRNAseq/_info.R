clusters_colours <- c(
  "0" = "#8B0000",
  "1" = "#ffbf00",
  "2" = "#e31a1c",
  "3" = "#c51b8a",
  "4" = "#C4A484",
  "5" = "#FFFF00",
  "6" = "#e31a1c",
  "7" = "#fb6a4a",
  "8" = "#fcae91",
  "9" = "#cab2d6",
  "10" = "#6a3d9a",
  "11" = "#953553",
  "12" = "#A020F0",
  "13" = "#1c1c84",
  "14" = "#633200",
  "15" = "#C7ea46",
  "16" = "#FFFF00",
  "17" = "#FF69B4",
  "18" = "#969696",
  "19" = "#cccccc",
  "20" = "#525252",  
  "21" = "#fcae91",
  "22" = "#cdfbb2",
  "23" = "#2171b5",
  "24" = "#2CFF05",
  "25" = "#6baed6",
  "26" = "#a6cee3",
  "27" = "#054aaa",
  "28" = "#c6dbef",
  "29" = "#01796f",
  "30" = "#00cec8",
  "31" = "#08519c",      
  "32" = "#98fb98",
  "33" = "#633200",
  "34" = "#00FFFF",
  "35" = "#A52A2A",
  "36" = "#2CFF05",
  "37" = "#fb6a4a"
)

col_treatment <- c("Skin" = "#4575b4", 
                   "WoundEarly" = "#ff7f00", 
                   "WoundLate" = "#b2182b")

cc_col <- c("Undivided" = "#ffffbf",
            "G1" = "#2b83ba",
            "S" = "#d7191c",
            "G2M" = "#fdae61")

extra_gene_check <- c("ITGA6" = "fibroblast basement membrane", 
                      "LAMA3" = "fibroblast basement membrane", 
                      "COL6A6" = "fibroblast basement membrane", 
                      "lamc1"= "progenitor, central" ,
                      "hspg2"= "progenitor, central" ,
                      "EBF1"= "progenitor, unknown" ,
                      "pax7a"= "progenitor, unknown satellite" ,
                      "ENSSSAG00000064587"= "progenitor, unknown satellite" ,
                      "ENSSSAG00000006052"= "progenitor, unknown satellite" ,
                      "MYOG"= "muscle, differentiating" ,
                      "Ltk"= "iridophore" ,
                      "mitfa"= "iridophore" ,
                      "sox10"= "iridophore" ,
                      "TFEC"= "iridophore" ,
                      "pdgfrb"= "progenitor, pericyte" ,
                      "EBF1" = "progenitor, pericyte", 
                      "ifitm5" = "Chondrocyte", 
                      "PANX3" = "Chondrocyte", 
                      "ENSSSAG00000030457" = "Chondrocyte", 
                      "pax7" = "progenitor, muscle", 
                      "scxa" = "fasciacyte_tenocyte, tenocyte", 
                      "HAS2" = "fasciacyte_tenocyte, fasciacytes",
                      "ENSSSAG00000009429" = "fasciacyte_tenocyte", 
                      "ENSSSAG00000072358" = "oestoblast", 
                      "ENSSSAG00000054025" = "oestoblast",
                      "TAGLN" = "Vasculature",
                      "spon1b" = "Vasculature",
                      "nrk2" = "Vasculature",
                      "CLEC3b" = "Vasculature",
                      "COL5a2b"  = "Vasculature",
                      "ENSSSAG00000000109" = "Vasculature", 
                      "XLKD1"= "Endothelial, lymphatic" ,
                      "MRC1"= "Endothelial, lymphatic",
                      "EPHA4"= "Endothelial, blood",
                      "NRP1"= "Endothelial, blood",
                      "SELE"= "Endothelial, blood",
                      "Kcnj8"= "Pericyte",
                      "Abcc9"= "Pericyte",
                      "Vtn"= "Pericyte",
                      "hspg2"= "Pericyte",
                      "EBF1"= "Pericyte",
                      "pdgfrb"= "Pericyte",
                      "pdgfra" = "Progenitor"
)


# Label colours
clusterid_to_celltype_col <-  readxl::read_excel(cell_cluster_col_label, sheet="fig1")
celltype_col <- clusterid_to_celltype_col$colour
names(celltype_col) <- clusterid_to_celltype_col$cell_type
clusterid_to_celltype_col2 <-  readxl::read_excel(cell_cluster_col_label, sheet="full")
celltype_col2 <- clusterid_to_celltype_col2$colour
names(celltype_col2) <- clusterid_to_celltype_col2$cell_type

#clusterid_to_celltype_col2 <-  readxl::read_excel(gsub(".xlsx", "2.xlsx", cell_cluster_col_label))
#celltype_col2 <-  clusterid_to_celltype_col$colour
#names(celltype_col2) <- clusterid_to_celltype_col2$cell_type
