library(arrow)
library(GenomicRanges)
library(dplyr)

args <- commandArgs(TRUE)
input = args[1]
output = args[2]

read_parquet_dir <- function(path){
  ## Define the dataset
  DS <- arrow::open_dataset(sources = path)
  ## Create a scanner
  SO <- Scanner$create(DS)
  ## Load it as n Arrow Table in memory
  AT <- SO$ToTable()
  ## Convert it to an R data frame
  DF <- as.data.frame(AT)
  remove(DS)
  remove(SO)
  remove(AT)
  return(DF)
}

contacts_to_haplotagged_pairs <- function(df){
  df$contact_idx <- table(df$read_name) %>% 
    lapply(.,function(x) seq(1,x)) %>% unlist() %>% as.numeric()
  pair_df <- data.frame(readID = paste(df$read_name,df$contact_idx,sep = "_"),
                        chr1 = df$align1_chrom,
                        pos1 = as.integer(round(0.5*(df$align1_fragment_start + df$align1_fragment_end))),
                        chr2 = df$align2_chrom,
                        pos2 = as.integer(round(0.5*(df$align2_fragment_start + df$align2_fragment_end))),
                        strand1 = factor(df$align1_strand),
                        strand2 = factor(df$align2_strand),
                        phase1 = factor(df$align1_haplotype),
                        phase2 = factor(df$align2_haplotype))
  levels(pair_df$strand1) <- c("-","+")
  levels(pair_df$strand2) <- c("-","+")
  levels(pair_df$phase1) <- c(".","0","1")
  levels(pair_df$phase2) <- c(".","0","1")
  return(pair_df)
}

contact_df <- read_parquet_dir(input)
print(head(contact_df))
pair_df <- contacts_to_haplotagged_pairs(contact_df)
print(head(pair_df))

write.table(pair_df, output,sep = "\t",quote = F,row.names = F,col.names = F)
