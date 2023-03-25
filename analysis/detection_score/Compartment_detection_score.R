library(dplyr)
library(pbmcapply)
library(Matrix)

## identified ABA/BAB switch triplets -------------------------------------------------------------

E1.df <- read.csv("DSG_1.5mM/GM_1.5mM_all_merged_chronly_1M_E1.csv",row.names = 1)
E1.df$sign <- sign(E1.df$E1)

library(gUtils)
E1.gr <- makeGRangesFromDataFrame(E1.df,keep.extra.columns = T,starts.in.df.are.0based = T)
E1.rd.gr <- gr.reduce(E1.gr,by = "sign") %>% sort()
E1.rd.gr$win <- 1:length(E1.rd.gr)
E1.rd.gr <- E1.rd.gr[!is.na(E1.rd.gr$sign),]
E1.rd.grl <- split(E1.rd.gr,E1.rd.gr@seqnames)

require(zoo)
lapply(E1.rd.grl,print)
lapply(E1.rd.grl,length) %>% unlist()

AB.switch <- lapply(E1.rd.grl, function(x){
  head(x)
  if(length(x)>=3){
    rollapply(x$win, width = 3, by = 1,
              FUN = print, align = "left") %>% 
      as.data.frame() %>%
      filter(V2-V1==1 & V3-V2==1)
  }
}) %>% Reduce(rbind,.)

dim(AB.switch)
AB.switch.grl <- apply(AB.switch, 1, function(x){
  E1.rd.gr[E1.rd.gr$win %in% x,]
})

lapply(AB.switch.grl, function(gr){
  gr@seqnames %>% as.character() %>% unique()
}) %>% unlist() %>% table()

lapply(AB.switch.grl, function(gr){
  width(gr)/1000000
})


ABA.region <- lapply(AB.switch.grl, function(gr){
  subsetByOverlaps(E1.gr,gr)
})
lapply(ABA.region, length) %>% unlist()
saveRDS(ABA.region,"detection_score/GM_1.5mM_1Mb_AB_swith_region.rds")


# functions ---------------------------------------------------------------

pair2binmtx <- function(df, startbin, endbin){
  if (!is.null(df)) {
    df <- df[df$V1 %in% startbin:endbin &
               df$V2 %in% startbin:endbin,]
    df$V1 <- factor(df$V1,levels = startbin:endbin)
    df$V2 <- factor(df$V2,levels = startbin:endbin)
    nbin = endbin - startbin + 1
    if (nbin <= 0) {
      stop("end <= start !")
    }
    mtx <- sparseMatrix(
      i = as.numeric(df$V1),
      j = as.numeric(df$V2),
      x = 1,
      dimnames = list(levels(df$V1), levels(df$V2)),
      symmetric = T,
      dims = c(nbin,nbin)
    )
  } else {
    print("empty dataframe!")
    mtx <- NULL
  }
  return(mtx)
}

get_expect_counts <- function(gr){
  RLE <- rle(gr$sign)
  intra <- sum((RLE$lengths)^2) + 2*RLE$lengths[1]*RLE$lengths[3]
  inter <- 2*RLE$lengths[2]*(RLE$lengths[1] + RLE$lengths[3])
  return(c(inter,intra))
}
get_observed_counts <- function(mtx,gr){
  A_idx = as.character(gr[gr$sign == 1,]$X)
  B_idx = as.character(gr[gr$sign == -1,]$X)
  intra <- sum(mtx[A_idx,A_idx]) + sum(mtx[B_idx,B_idx])
  inter <- sum(mtx[A_idx,B_idx]) + sum(mtx[B_idx,A_idx])
  return(c(inter,intra))
}
detection_score <- function(df,grl,start,end){
  mtx <- pair2binmtx(df,start,end)
  DS.v <- lapply(grl, function(gr){
    Exp <- get_expect_counts(gr)
    Obs <- get_observed_counts(mtx,gr)
    DS <- diff(Obs/Exp)
    return(DS)
  }) %>% unlist()
  return(DS.v)
}


# calculate compartment DS ------------------------------------------------

chr = "all"
df.list <- list.files(path = "detection_score",
                      pattern = paste0("_GRCh38_PG_1000000_all_mtx.gz$"),full.names = T)
names(df.list) <- gsub("MboI_","",basename(df.list)) %>%
  gsub(paste0("_GRCh38_PG_1000000_all_mtx.gz$"),"",.)
df <- lapply(df.list, function(x) read.table(x))

bins <- read.table("detection_score/chronly_sort_1M_bins.txt")
colnames(bins) <- c("chrom","start","end")
bins$id <- 0:(nrow(bins)-1)

ABA.region <- readRDS("detection_score/GM_1.5mM_1Mb_AB_swith_region.rds")

ABA.region.n <- lapply(ABA.region, function(x) unique(x@seqnames)) %>% 
  unlist() %>% as.data.frame() %>% `colnames<-`("chrom") %>%
  group_by(chrom) %>% mutate(ID = paste(chrom,1:length(chrom),sep = "_"))
names(ABA.region) <- ABA.region.n$ID
saveRDS(ABA.region,"detection_score/GM_1.5mM_1Mb_AB_swith_region.rds")

DS.df <- lapply(paste0("chr",c(1:22,"X")), function(chr){
  ABA.used <- lapply(ABA.region, function(x) unique(x@seqnames) ==  chr) %>% unlist()
  ABA.region.used <- ABA.region[ABA.used]
  
  chr.range <- c(min(bins[bins$chrom == chr,]$id),
                 max(bins[bins$chrom == chr,]$id))
  
  DS <- pbmcapply::pbmclapply(names(df), function(x){
    detection_score(df[[x]],ABA.region.used,chr.range[1],chr.range[2])
  },mc.cores = 30) %>% Reduce(rbind,.) %>% 
    `row.names<-`(names(df))
  return(DS)
})
names(DS.df) <- paste0("chr",c(1:22,"X"))

saveRDS(DS.df,"GM_1.5mM_1Mb_compartment_detection_score.rds")


# annotations ----------------------------------------------------------

meta <- read.csv("scHiC_metadata.csv",row.names = 1)
meta <- meta[meta$DSG == 1.5 & meta$cell_line == "GM",]
meta$library2 <- meta$library
meta[!(meta$library %in% c("GM_24","GM_500")),]$library2 <- "GM_96"

meta$rename2 <- paste(meta$library,meta$name,sep = ".")
rownames(meta) <- meta$rename2

ann_row = meta[,c("library2","filtered_pass_count")]
ann_row$log10counts <- log10(ann_row$filtered_pass_count)
ann_row <- ann_row[,c("library2","log10counts")]
ann_color = list(library2 = setNames(brewer.pal(9,"Set3")[c(4,7,5)],
                                     unique(meta$library2)),
                 chr = setNames(c(brewer.pal(9,"Set1"),
                                  brewer.pal(12,"Set3"),
                                  brewer.pal(8,"Accent")[c(8,6,5)]),
                                paste0("chr",c(1:22,"X","Y"))))

# random compartment score ------------------------------------------------

mtx.list <- list.files(path = "detection_score",
                       pattern = "_GRCh38_PG_1000000_all_mtx.gz$",
                       full.names = T)
names(mtx.list) <- gsub("MboI_","",basename(mtx.list)) %>%
  gsub("_GRCh38_PG_1000000_all_mtx.gz","",.)
mtx.df <- lapply(mtx.list, function(x){
  print(x)
  try(read.table(x))
})

bins.used <- GRangesList(ABA.region) %>% unlist() %>% 
  as.data.table()

detection_score <- function(mtx,bin1,bin2){
  AA = nrow(mtx[mtx$V1 %in% bin1 & mtx$V2 %in% bin1,])
  BB = nrow(mtx[mtx$V1 %in% bin2 & mtx$V2 %in% bin2,])
  AB = nrow(mtx[mtx$V1 %in% bin1 & mtx$V2 %in% bin2,])
  BA = nrow(mtx[mtx$V1 %in% bin2 & mtx$V2 %in% bin1,])
  Eintra = length(bin1)^2 + length(bin2)^2
  Einter = length(bin1)*length(bin2)*2
  DS = (AA + BB)/Eintra - (AB + BA)/Einter
  return(DS)
}

run_CompDS <- function(i,grl = ABA.region,
                              mc1 = 10, mc2 = 10){
  df  <- pbmclapply(mtx.df, function(mtx){
    mclapply(grl, function(gr){
      bin.group <- split(1:length(gr),gr$sign)
      bin1 <- gr$X[bin.group[[1]]]
      bin2 <- gr$X[bin.group[[2]]]
      DS <- detection_score(mtx,bin1,bin2)
      return(DS)
    }, mc.cores = mc1) %>% unlist()
  },mc.cores = mc2) %>% Reduce(rbind,.)
  rownames(df) <- names(mtx.df)
  colnames(df) <- names(grl)
  return(df)
}

comp_DS <- run_CompDS(1)
ABA.region$chr1_1@seqnames %>% unique()

run_CompDS_random <- function(i,grl = ABA.region,
                              mc1 = 10, mc2 = 10){
  message("run ",i," ...")
  df  <- pbmclapply(mtx.df, function(mtx){
    pbmclapply(grl, function(gr){
      chr.used <- gr@seqnames %>% unique()
      bins <- unique(bins.used[bins.used$seqnames == chr.used,]$X)
      start.bin <- sample(1:(length(bins) - length(gr) + 1),1)
      rbins <- bins[start.bin:(start.bin + length(gr) - 1)]
      bin.group <- split(1:length(gr),gr$sign)
      bin1 <- rbins[bin.group[[1]]]
      bin2 <- rbins[bin.group[[2]]]
      DS <- detection_score(mtx,bin1,bin2)
      return(DS)
    }, mc.cores = mc1) %>% unlist()
  },mc.cores = mc2) %>% Reduce(rbind,.)
  rownames(df) <- names(mtx.df)
  colnames(df) <- names(grl)
  return(df)
}


comp_DS_ran[[1]] <- run_CompDS_random(1)

comp_DS_ran <- list()
for (i in 1:100) {
  comp_DS_ran[[i]] <- run_CompDS_random(i)
}
saveRDS(comp_DS_ran,"compartment_DS_random_100.rds")


# Normalized compartment DS -----------------------------------------------

comp_DS_ran <- readRDS("detection_score/compartment_DS_random_100.rds")
comp_DS_ran[[9]]
lapply(comp_DS_ran, dim)
comp_DS_ran.m <- Reduce(`+`, comp_DS_ran)
comp_DS_ran.m <- comp_DS_ran.m/100

comp_DS <- readRDS("detection_score/compartment_DS_result.rds")

comp_DS.norm <- comp_DS - comp_DS_ran.m
identical(colnames(comp_DS),colnames(comp_DS_ran.m))
colnames(comp_DS) %in% colnames(comp_DS_ran.m)
identical(rownames(comp_DS),rownames(comp_DS_ran.m))

cell.used <- meta.df[meta.df$filtered_pass_count > 10000  &
                       meta.df$cell_line == "GM12878",]$rename2
comp_DS.norm <- comp_DS.norm[cell.used,]
comp_DS.df <- data.frame(mean = rowMeans(comp_DS.norm),
                         cell = rownames(comp_DS.norm))
comp_DS.df <- cbind(comp_DS.df,ann_row[rownames(comp_DS.norm),])
comp_DS.df2 <- data.frame(mean = colMeans(comp_DS.norm),
                          pair = colnames(comp_DS.norm))

write.csv(comp_DS.df,"detection_score/cell_mean_normalized_compartment_detection_score.csv")
write.csv(comp_DS.df2,"detection_score/compSwitch_mean_normalized_compartment_detection_score.csv")

pdf("detection_score/compartment_normalized_detection_score_histogram_20221110.pdf",
    width = 8,height = 4)
ggplot(comp_DS.df,aes(x=mean)) +
  geom_histogram(color = "black",fill = "grey") +
  geom_vline(xintercept = 0,linetype = "dashed",color = "red") +
  xlab("Detection score") +
  theme_bw() + ylab("Number of cells") +
  theme(panel.grid = element_blank()) +
  ggplot(comp_DS.df2,aes(x=mean)) +
  geom_histogram(color = "black",fill = "grey") +
  geom_vline(xintercept = 0,linetype = "dashed",color = "red") +
  xlab("Detection score") +
  theme_bw() + ylab("Number of compartment switches") +
  theme(panel.grid = element_blank())
dev.off()

# summary -----------------------------------------------------------------

comp_DS.df <- read.csv("detection_score/cell_mean_normalized_compartment_detection_score.csv")
comp_DS.df2 <- read.csv("detection_score/compSwitch_mean_normalized_compartment_detection_score.csv")
comp_DS.df2$pos_ratio <- apply(comp_DS.norm, 2, function(x) sum(x>0)/length(x))


meta.df <- read.csv("scHiC_all_metadata_rename_all.csv")
rownames(meta.df) <- meta.df$rename2

cell.m <- left_join(comp_DS.df,
                    meta.df[,c("rename2","cell_name","group2")],
                    by = c("cell" = "rename2"))
write.csv(cell.m,"detection_score/cell_mean_normalized_compartment_detection_score_summary.csv")


ABA.region.df <- lapply(names(ABA.region), function(x){
  ABA.region[[x]] %>% 
    gr.reduce(by = "sign") %>%
    as.data.frame() %>% 
    arrange(start) %>%
    mutate(ID = x)
}) %>% Reduce(rbind,.)
ABA.region.df$sign <- factor(ABA.region.df$sign)
levels(ABA.region.df$sign) <- c("B","A")
ABA.region.df$range <- paste0(ABA.region.df$seqnames,":",
                              ABA.region.df$start,"-",
                              ABA.region.df$end)
ABA.region.sum <- ABA.region.df %>% 
  group_by(ID) %>%
  summarise(range = paste(range,collapse = ";"),
            width = paste(width,collapse = ";"),
            compartment = paste(sign,collapse = ""),
            CompScore = paste(format(round(E1,digits = 2),nsmall = 2),collapse = ";")) %>%
  as.data.frame() %>%
  `row.names<-`(.$ID)

library(stringr)

ABA.region.sum <- ABA.region.sum[str_sort(ABA.region.sum$ID,numeric = T),]
comp_DS.df2 <- left_join(comp_DS.df2,ABA.region.sum,
                         by = c("pair" = "ID"))
write.csv(comp_DS.df2,"detection_score/compSwitch_mean_normalized_compartment_detection_score_summary.csv")

nrow(comp_DS.df2[comp_DS.df2$mean > 0,])
nrow(comp_DS.df2)

summary(comp_DS.df$mean)
dim(comp_DS.df)
nrow(comp_DS.df[comp_DS.df$mean > 0,])

apply(comp_DS.norm, 2, function(x) sum(x>0)/length(x)) %>% 
  unlist() %>% summary()


# shuffle times test ------------------------------------------------------

comp_DS_ran <- readRDS("detection_score/compartment_DS_random_100.rds")

lapply(comp_DS_ran, dim)

comp_DS_ran.m <- Reduce(`+`, comp_DS_ran)
comp_DS_ran.m <- comp_DS_ran.m/100

comp_DS_ran.df <- lapply(1:100, function(i){
  colMeans(comp_DS_ran[[i]]) %>% 
    as.data.frame() %>%
    `colnames<-`("cell_mean") %>%
    mutate(group = i,
           ID = rownames(.))
}) %>% Reduce(rbind,.)

comp_DS_ran.df <- comp_DS_ran.df %>%
  group_by(ID) %>%
  mutate(cum = cumsum(cell_mean)) %>%
  mutate(cummean = cum/group) %>%
  mutate(dev = cummean - last(cummean)) %>%
  mutate(norm_dev = dev/last(cummean))
comp_DS_ran.df$chr <- factor(gsub("_.*","",comp_DS_ran.df$ID),
                             levels = paste0("chr",c(1:22,"X")))

pdf("figure/NM_revision/comp_shuffle_deviations_test.pdf",
    width = 5.5,height = 4)
ggplot(comp_DS_ran.df,aes(x=group,y=norm_dev*100)) +
  geom_line(aes(group = ID,color = chr)) +
  geom_hline(yintercept = c(5,1,-5,-1),
             linetype = "dashed") +
  theme_bw() +
  xlab("Number of shuffle times") +
  labs(title = "Compartment detection scores") +
  ylab("Percentage of deviations %") +
  scale_y_continuous(breaks = c(-10,-5,-1,0,1,5,10)) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,vjust = 0.5))
dev.off()
