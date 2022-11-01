library(arrow)
library(GenomicRanges)
library(dplyr)

args <- commandArgs(TRUE)
sample_name = args[1]
enzyme = args[2]
ref = args[3]
phase = args[4]
parquet.file = args[5]

data.pq.flt <- read_parquet(parquet.file)
print(head(data.pq.flt))
flt <- data.pq.flt[data.pq.flt$filter == "pass",]

flt.stat1 <- c(
    adjacent = nrow(data.pq.flt[data.pq.flt$filter == "adjacent",]),
    close = nrow(data.pq.flt[data.pq.flt$filter == "close",]),
    duplication = nrow(data.pq.flt[data.pq.flt$filter == "duplication",]),
    promiscuous = nrow(data.pq.flt[data.pq.flt$filter == "promiscuous",]),
    isolated = nrow(data.pq.flt[data.pq.flt$filter == "isolated",]),
    pass = nrow(data.pq.flt[data.pq.flt$filter == "pass",]))
print(flt.stat1)

flt_Stat2 <- function(flt) {
  trans_ratio = sum(!flt$contact_is_cis)/nrow(flt)*100
  unphased = nrow(flt[flt$align1_haplotype == -1 &
                        flt$align2_haplotype == -1,])
  phased = nrow(flt[flt$align1_haplotype > 0 &
                      flt$align2_haplotype > 0,])
  semi_phased = nrow(flt) - unphased - phased
  res = c(trans_ratio,phased/nrow(flt)*100,semi_phased/nrow(flt)*100)
  names(res) <- c("trans_ratio","full_phased_ratio","semi_phased_ratio")
  res = round(res,digits = 2)
  return(res)
}

order_stat <- function(flt){
  order.stat <- rbind(flt[,c("read_name","align1_fragment_id")] %>%
                        `colnames<-`(c("read_name","fragment_id")),
                      flt[,c("read_name","align2_fragment_id")] %>%
                        `colnames<-`(c("read_name","fragment_id")))
  order.stat <- unique(order.stat)
  order.stat <- table(order.stat$read_name) %>% table() %>% as.data.frame() %>%
    `colnames<-`(c("times","Freq"))
  order.stat$pct <- order.stat$Freq/sum(order.stat$Freq)*100
  return(order.stat)
}

mn_phased <- function(flt){
  df <- unique(rbind(flt[,c("read_idx","align1_align_idx","align1_fragment_id","align1_haplotype")] %>%
                       `colnames<-`(c("read_idx","align_idx","align_fragment_id","haplotype")),
                     flt[,c("read_idx","align1_align_idx","align2_fragment_id","align2_haplotype")] %>%
                       `colnames<-`(c("read_idx","align_idx","align_fragment_id","haplotype"))))
  dim(df)
  dim(unique(df[1:4]))
  sum.df <- table(df$haplotype) %>% as.data.frame()
  hp1_ratio = sum(sum.df[sum.df$Var1 == 1,]$Freq)/sum(sum.df$Freq)*100
  hp2_ratio = sum(sum.df[sum.df$Var1 == 2,]$Freq)/sum(sum.df$Freq)*100
  ratio = hp1_ratio + hp2_ratio
  res = setNames(c(hp1_ratio,hp2_ratio,ratio),
               c("HP1_ratio","HP2_ratio","phased_ratio"))
  res = round(res,digits = 2)
  return(res)
}

len_stat <- function(flt){
  df <- rbind(flt[,c("read_idx","align1_align_idx","align1_fragment_id","align1_start","align1_end")] %>%
                `colnames<-`(c("read_idx","align_idx","align_fragment_id","start","end")),
              flt[,c("read_idx","align2_align_idx","align2_fragment_id","align2_start","align2_end")] %>%
                `colnames<-`(c("read_idx","align_idx","align_fragment_id","start","end")))
  df <- unique(df)
  df$length <- abs(df$end - df$start)
  cat("summary of monomer length:\n")
  print(summary(df$length))
  
  df2 <- flt[,c("read_idx","read_length")] %>% unique()
  cat("summary of read length:\n")
  print(summary(df2$read_length))
  len_stat <- setNames(c(median(df$length),
                         median(df2$read_length),
                         nrow(df2)),
                       c("median_monomer_length",
                         "median_read_length",
                         "pass_reads_counts"))
  return(len_stat)
}

flt.stat2 <- flt_Stat2(flt)
print(flt.stat2)
order.stat <- order_stat(flt)
high_order = setNames(round(100-order.stat$pct[1],digits = 2),"high_order")

flt.stat3 <- mn_phased(flt)
len_stat <- len_stat(flt)

flt.stat <- c(flt.stat2,len_stat,high_order,flt.stat3,flt.stat1) %>% as.data.frame()
colnames(flt.stat) <- sample_name
flt.stat$variable <- rownames(flt.stat)
print(flt.stat)

write.table(flt.stat,paste(sample_name,enzyme,ref,phase,"contact_filtration_stats.txt",sep = "_"),
            row.names = F,quote = F,sep = "\t")
write.table(order.stat,paste(sample_name,enzyme,ref,phase,"order_stats.txt",sep = "_"),
            row.names = F,quote = F,sep = "\t")