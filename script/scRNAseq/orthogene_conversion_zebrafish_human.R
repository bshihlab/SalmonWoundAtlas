### custom functions
## function for returning the most common value
fun.return_common_val <- function(x) {
  m <- which.max(table(as.character(x)))
  as.character(x)[m]
}

## function for simplifying zebrafish gene names
fun.simplify_zebrafish_genename <- function(x){
  x <- x[x!=""]
  # remove anything after .
  x <- sapply(strsplit(x, "\\."), function(x)x[[1]])
  x <- x[order(x)]
  x_unique <- unique(x)
  if(length(x_unique) == 1){
    # if there is only 1 unique value, export that value
    x_out <- x_unique
  } else {
    x_trimmed <- sapply(x_unique, FUN=function(s)substring(s,1, nchar(s)-1))
    x_trimmed_unique <- unique(x_trimmed)
    # if different values only differ by the last letter, 
    # and the last letter is a or b, trim last letter
    last_letter <- substr(x_trimmed[1], nchar(x_trimmed[1]), nchar(x_trimmed[1]))
    if(length(x_trimmed_unique) == 1 &
       last_letter %in% c("a", "b", "A", "B")
       ){
      x_out <- x_trimmed_unique
    } else {
      # otherwise return the most common value
      x_out<- fun.return_common_val(x)
    }
  }
  return(x_out)
}

## take the merged table for each gene and return summarised entries
fun.summarise_gene <- function(x){
  out_zebrafish.geneid <- paste0(unique(x$zebrafish.geneid), collapse = "|")
  out_human.geneid <- paste0(unique(x$human.geneid), collapse = "|")
  out_zebrafish.genename <- fun.simplify_zebrafish_genename(x$zebrafish.genename)
  out_zebrafish.genename <- ifelse(length(out_zebrafish.genename) == 0, NA, out_zebrafish.genename)
  out_zebrafish.genename <- ifelse(out_zebrafish.genename == "NA", NA, out_zebrafish.genename)
  out_zebrafish.genename_full <- paste0(unique(x$zebrafish.genename), collapse = "|")
  out_zebrafish.genename_full <- ifelse(out_zebrafish.genename_full == "NA", NA, out_zebrafish.genename_full)
  out_human.genename <- fun.return_common_val(x$human.genename)
  out_human.genename <- ifelse(length(out_human.genename) == 0, NA, out_human.genename)
  out_human.genename <- ifelse(out_human.genename == "NA", NA, out_human.genename)
  out_human.genename_full <- paste0(unique(x$human.genename), collapse = "|")
  out_human.genename_full <- ifelse(out_human.genename_full == "NA", NA, out_human.genename_full)
  out <- c('salmon.geneid' = x$salmon.geneid[1],
           'zebrafish.geneid' = out_zebrafish.geneid,
           'human.geneid' = out_human.geneid,
           'zebrafish.genename' = out_zebrafish.genename,
           'zebrafish.genename_full' = out_zebrafish.genename_full,
           'human.genename' = out_human.genename,
           'human.genename_full' = out_human.genename_full)
  return(out)
}


## organise a orthologue conversion table
orth_zeb2sal <- read.delim("annotation/gene/zebrafish_to_salmon_20260418.txt")
orth_zeb2hum <- read.delim("annotation/gene/zebrafish_to_human_20260418.txt")
colnames(orth_zeb2sal) <- c("zebrafish.geneid", "salmon.geneid", "zebrafish.genename")
colnames(orth_zeb2hum) <- c("human.geneid", "zebrafish.geneid", "human.genename")
orth_zebsalhum <- merge(orth_zeb2sal, orth_zeb2hum, by="zebrafish.geneid", all.x=TRUE)
# remove entries where salmon does not have geneid
orth_zebsalhum <- orth_zebsalhum[orth_zebsalhum$salmon.geneid != "",]

## take the most common non-na value for zebrafish genename and human gene name so each salmon entry only occurs once
orth_zebsalhum <- split(orth_zebsalhum, orth_zebsalhum$salmon.geneid)
# take in a list of zebrafish gene name, return the one with highest occurrence and note if there is conflict
orth_zebsalhum_summarised <- lapply(orth_zebsalhum, FUN=function(x)fun.summarise_gene(x))

orth_zebsalhum_df <- as.data.frame(do.call(rbind, orth_zebsalhum_summarised))

# create "genename_base" using human gene name if there is only 1 match, otherwise use zebrafish name.
# Human names are preferentially used due to easier gene set enrichment analysis
orth_zebsalhum_df$gene_name_base <- ifelse((!(is.na(orth_zebsalhum_df$human.genename))) & (!(grepl("\\|", orth_zebsalhum_df$human.genename_full))),
                                           orth_zebsalhum_df$human.genename, 
                                           toupper(orth_zebsalhum_df$zebrafish.genename))
orth_zebsalhum_df$gene_name_base <- tolower(orth_zebsalhum_df$gene_name_base)

write.csv(orth_zebsalhum_df, "annotation/gene/orth_zebsalhum.csv", row.names=FALSE)


## Make covnersion tables for zebrafish and human for analysing external data

