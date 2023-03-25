library(rtracklayer)
library(MASS)
library(Matrix)
library(dplyr)
library(pbmcapply)
library(tidyverse)

prefix = "GM12878"

# high-order ratio --------------------------------------------------------

wins <- read.table("synergy/GRCh38_50k_windows.bed") 
colnames(wins) <- c("chr","start","end")
wins$ID <- 1:nrow(wins)
wins.gr <- makeGRangesFromDataFrame(wins,keep.extra.columns = T,
                                    starts.in.df.are.0based = T)

mono.gr <- readRDS("synergy/mono_all_gr.rds")
reads.sum <- table(mono.gr$read_idx) %>% as.data.frame()

mono.gr$multiway <- FALSE
mono.gr[mono.gr$read_idx %in% reads.sum[reads.sum$Freq > 2,]$Var1,]$multiway <- TRUE

ov = wins.gr %*% mono.gr[, c('read_idx','cell','multiway')] %>% as.data.table

ov.read.uniq <- ov[,c("query.id","read_idx","multiway")] %>% unique()
ov.cell.uniq <- ov[,c("query.id","cell","multiway")] %>% unique()
ov.cell.uniq <- ov.cell.uniq %>%
  group_by(query.id,cell) %>%
  summarise(multiway = as.logical(sum(multiway)))

summary_ov <- function(ov,wins.df = wins){
  c.sum <- table(ov[,c("query.id","multiway")]) %>% as.data.frame.matrix()
  c.sum$ID <- as.integer(rownames(c.sum))
  print(head(c.sum))
  
  c.sum <- left_join(wins.df,c.sum,by = "ID")
  c.sum[is.na(c.sum)] <- 0
  c.sum$total <- c.sum$`TRUE` + c.sum$`FALSE`
  c.sum$ratio <- c.sum$`TRUE`/c.sum$total
  c.sum[is.na(c.sum)] <- 0
  c.sum.gr <- makeGRangesFromDataFrame(c.sum,keep.extra.columns = T,
                                       starts.in.df.are.0based = T)
  return(c.sum.gr)
}

ov.reads.sum <- summary_ov(ov.read.uniq)
ov.cells.sum <- summary_ov(ov.cell.uniq)

saveRDS(ov.cells.sum,paste0(prefix,'_cell_high_order_ov_wins_gr.rds'))
saveRDS(ov.reads.sum,paste0(prefix,'_read_high_order_ov_wins_gr.rds'))

library(pbmcapply)
pbmclapply(1:500, function(x){
  sum(duplicated(ov.read.uniq[ov.read.uniq$query.id == x,]$read_idx))
},mc.cores = 40) %>% unlist()


ov.cells.sum <- readRDS(paste0('synergy/',prefix,'_cell_high_order_ov_wins_gr.rds'))
ov.reads.sum <- readRDS(paste0('synergy/',prefix,'_read_high_order_ov_wins_gr.rds'))

write.table(as.data.frame(ov.reads.sum)[,c("seqnames","start","end","ratio")],
            paste0('synergy/',prefix,'_reads_high_order.bedGraph'),
            row.names = F,col.names = F,sep = '\t',quote = F)
write.table(as.data.frame(ov.cells.sum)[,c("seqnames","start","end","ratio")],
            paste0('synergy/',prefix,'_cells_high_order.bedGraph'),
            row.names = F,col.names = F,sep = '\t',quote = F)