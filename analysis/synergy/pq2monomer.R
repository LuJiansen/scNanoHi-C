library(arrow)
library(dplyr)

args <- commandArgs(TRUE)
sample_name = args[1]
parquet.file = args[2]
file.out = args[3]
prefix = args[4]

pq <- read_parquet(parquet.file)
df <- rbind(pq[,c("align1_chrom","align1_start","align1_end","read_name","read_length","align1_haplotype")] %>%
              `colnames<-`(c("chrom","start","end","read_name","read_length","haplotype")),
            pq[,c("align2_chrom","align2_start","align2_end","read_name","read_length","align2_haplotype")] %>%
              `colnames<-`(c("chrom","start","end","read_name","read_length","haplotype")))
df <- df[!duplicated(df),]
#df <- df[df$chrom %in% paste0("chr",c(1:22,"X","Y")),]
df$read_idx <- factor(df$read_name) %>% as.numeric()
df <- droplevels(df)
df$cell <- paste(prefix,sample_name,sep = ".")

write.table(df,file = file.out,sep = "\t",row.names = F,quote = F,col.names = F)