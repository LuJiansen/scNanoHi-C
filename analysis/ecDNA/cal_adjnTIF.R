library(dplyr)
library(pbmcapply)

args <- commandArgs(TRUE)
bins_file = args[1]
pair_file = args[2]
cnv_file = args[3]
resolution = args[4]
name = args[5]
threads = as.integer(args[6])

## functions ---------------------------------------------------------------

calc_TIF <- function(id,pair,bin){
  chr.ex <- bin[bin$ID == id,]$chr
  id.ex <- bin[bin$chr == chr.ex,]$ID
  message("\n chr of bin: ",chr.ex)
  pair.used <- pair %>% filter(i == id | j == id) %>%
    filter((!(i %in% id.ex)) & !(!(j %in% id.ex)))
  if (nrow(pair.used)) {
    n = sum(pair.used$x)
  } else {n = 0}
  return(n)
}

cnv_file_process <- function(file){
  df <- read.table(file,sep = '\t')
  colnames(df) <- c("chr","start","end","copy_number")
  df$chr <- factor(df$chr,levels = paste0("chr",c(1:22,"X","Y")))
  df <- df %>% arrange(chr,start)
  uni_win <- rbind(df[,c("chr","start")] %>% `colnames<-`(c("chr","pos")),
                   df[,c("chr","end")] %>% `colnames<-`(c("chr","pos"))) %>%
    unique()
  uni_win$chr <- factor(uni_win$chr,levels = paste0("chr",c(1:22,"X","Y")))
  uni_win <- uni_win %>% arrange(chr,pos)
  uni_win$ID <- paste(uni_win$chr,uni_win$pos,sep = "_")
  df$win1 <- factor(paste(df$chr,df$start,sep = "_"),levels = unique(uni_win$ID))
  df$win2 <- factor(paste(df$chr,df$end,sep = "_"),levels = unique(uni_win$ID))
  df$ID <- (as.numeric(df$win1) + as.numeric(df$win2))/2
  return(df)
}

## data --------------------------------------------------------------------

bins <- read.table(bins_file)
colnames(bins) <- c("chr","start","end")
bins <- bins[,1:3]
bins$ID <- 0:c(nrow(bins) - 1)
cnv <- cnv_file_process(cnv_file)
pair.df <- read.table(pair_file)
colnames(pair.df) <- c("i","j","x")

bins <- left_join(bins,cnv[,1:4])
bins$chr <- factor(bins$chr,levels = unique(bins$chr))

TIFbinCount <- lapply(unique(bins$chr), function(x){
  nrow(bins[bins$chr != x,])
}) %>% unlist() %>% setNames(.,unique(bins$chr))

print(head(bins))
bins$totalTIF <- pbmclapply(bins$ID, function(x){
     calc_TIF(x,pair.df,bins)
}, mc.cores = threads) %>% unlist()
bins$TIFbinCount <- TIFbinCount[bins$chr]

print(head(bins))
aveTIF = sum(bins$totalTIF)/sum(bins$TIFbinCount)
bins$nTIF <- bins$totalTIF/bins$TIFbinCount/aveTIF

nTIF.fit <- lm(nTIF~copy_number,
               data = bins[bins$copy_number < 10 &
                                 !is.na(bins$copy_number),])

bins$ExpTIF <- predict(nTIF.fit,data.frame(copy_number = bins$copy_number))
bins$adjnTIF <- bins$nTIF - bins$ExpTIF
bins$significant <- F
try(bins[bins$adjnTIF > 10 &
           !is.na(bins$adjnTIF),]$significant <- T)

saveRDS(bins,paste0(name,"_",resolution,"_adjnTIF.rds"))