#### convert sparse matrix into 2 data frames for saving in fst
# DF1: containg the summary data for the sparse matrix
# DF2: Record the first and the last row for id in DF1 for each gene

# This way, DF2 can be used to slice out the relevant line numbers in DF1,
# And then DF1 can be used to repopulate the orignal gene expression data

library(fst)

## DF1
exprs_data <- readMM("data/h5ad_data.mx") # this is from the rds folder from the main analyss
exprs_df <- as.data.frame(summary(exprs_data))
exprs_df <- exprs_df[order(exprs_df$i),]
row.names(exprs_df) <- 1:nrow(exprs_df)
write_fst(exprs_df , "data/exprs.fst", compress = 100)

# DF2
exprs_split <- split(exprs_df, exprs_df$i)
df1_row_nums <- lapply(exprs_split, function(x)c(rowMin=min(as.numeric(row.names(x))), 
                                                rowMax=max(as.numeric(row.names(x)))))
df1_row_nums_df <- as.data.frame(do.call(rbind, df1_row_nums))
write_fst(df1_row_nums_df , "data/exprs_rownum_indicate.fst", compress = 100)
