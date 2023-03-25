library(dplyr)
library(pbmcapply)
library(Matrix)
library(ggplot2)
library(RColorBrewer)



# calculate TAD DS for each chromosome ------------------------------------

args <- commandArgs(TRUE)
chr = args[1]

bins <- read.table("detection_score/chronly_sort_50k_bins.txt")
colnames(bins) <- c("chrom","start","end")
bins$id <- 0:(nrow(bins)-1)
bins.used <- bins[bins$chrom == chr,]

insul.used <- read.table("DSG_1.5mM/GM_1.5mM_all_merged_chronly_50k_res_500k_wins_insul_score.txt",
                         header = T)
insul.used <- insul.used[insul.used$chrom == chr,]

cont.df.list <- list.files(path = "data",
                           pattern = paste0("_GRCh38_PG_50000_",chr,"_mtx.gz$"),
                           full.names = T)
names(cont.df.list) <- gsub("MboI_","",basename(cont.df.list)) %>%
  gsub(paste0("_GRCh38_PG_50000_",chr,"_mtx.gz"),"",.)
cont.df <- pbmclapply(cont.df.list, function(x){
  tryCatch(read.table(x),
           error = function(e){
             cat(x,"is empty !\n")
             return(NULL)})
}, mc.cores = 20)

pair2binmtx <- function(df, startbin, endbin){
  if (!is.null(df)) {
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

TAD_detection_score <- function(boundary,mtx,nbin=20){
  min_idx = min(as.numeric(colnames(mtx)))
  max_idx = max(as.numeric(colnames(mtx)))
  idx1 = (boundary - nbin):(boundary - 1)
  idx2 = (boundary + nbin):(boundary + 1)
  idx1 = idx1[idx1 <= max_idx & idx1 >= min_idx]
  idx2 = idx2[idx2 <= max_idx & idx2 >= min_idx]
  idx1 <- as.character(idx1)
  idx2 <- as.character(idx2)
  intra <- sum(mtx[idx1,idx1]) + sum(mtx[idx2,idx2])
  inter <- sum(mtx[idx1,idx2]) + sum(mtx[idx2,idx1])
  DS = (intra-inter)/(nbin*nbin*2)
  return(DS)
}

chr.range <- c(min(bins[bins$chrom == chr,]$id),
               max(bins[bins$chrom == chr,]$id))

cont.mtx <- pbmclapply(cont.df, function(df){
  pair2binmtx(df,chr.range[1],chr.range[2])
}, mc.cores = 20)
saveRDS(cont.mtx,paste0("GM12878_",chr,"_contact_matrix.rds"))


# merge TAD detection score -----------------------------------------------------

TAD_DS <- pbmclapply(cont.mtx, function(mtx){
  pbmclapply(insul.used[insul.used$chrom == chr,]$id,
             function(x){TAD_detection_score(x,mtx)}) %>% unlist()
},mc.cores = 20) %>% Reduce(rbind,.)
rownames(TAD_DS) <- names(cont.mtx)
write.csv(TAD_DS,paste0("GM12878_",chr,"_TAD_detection_score.csv"))

# annotations -----------------------------------------------------


meta <- read.csv("scHiC_metadata.csv",row.names = 1)
meta <- meta[meta$DSG == 1.5 & meta$cell_line == "GM",]
meta$library2 <- meta$library
meta[!(meta$library %in% c("GM_24","GM_500")),]$library2 <- "GM_96"

meta$rename2 <- paste(meta$library,meta$name,sep = ".")
rownames(meta) <- meta$rename2
ann_color = list(library2 = setNames(brewer.pal(9,"Set3")[c(4,7,5)],
                                     unique(meta$library2)),
                 chr = setNames(c(brewer.pal(9,"Set1"),
                                  brewer.pal(12,"Set3"),
                                  brewer.pal(8,"Accent")[c(8,6,5)]),
                                paste0("chr",c(1:22,"X","Y"))))

# TAD detection score matrix ----------------------------------------------

TDS.list <- list.files(path = "detection_score", full.names = T,
                       pattern = "[0-9X]_TAD_detection_score.csv")
names(TDS.list) <- gsub("GM12878_","",basename(TDS.list)) %>%
  gsub("_TAD_detection_score.csv","",.)
TDS <- lapply(names(TDS.list), function(x){
  df <- read.csv(TDS.list[[x]],row.names = 1)
  colnames(df) <- paste0(x,"_T",1:ncol(df))
  return(df)
})
names(TDS) <- names(TDS.list)
TDS <- TDS[paste0("chr",c(1:22,"X"))]

lapply(TDS, dim)

TDS.mtx <- TDS %>% Reduce(cbind,.)
for (i in colnames(TDS.mtx)) {
  if(!is.numeric(TDS.mtx[[i]])){
    print(i)
  }
}

write.csv(TDS.mtx,"detection_score/TAD_DS_mtx.csv")

# random boundary detection scores ---------------------------------------------------------

random_TDS <- function(cont.mtx,n){
  bd <- sample(bins.used$id,n)
  nov <- bd %in% insul.used[insul.used$chrom == chr,]$id
  TAD_DS <- pbmclapply(cont.mtx, function(mtx){
    lapply(bd,function(x){TAD_detection_score(x,mtx)}) %>% unlist()
  },mc.cores = 10) %>% Reduce(rbind,.)
  rownames(TAD_DS) <- names(cont.mtx)
  return(list(boundary = bd,
              overlap = nov,
              DS = TAD_DS))
}

TDS.ran <- pbmclapply(1:100, function(x){
  print(x)
  random_TDS(cont.mtx,nrow(insul.used))
}, mc.cores = 5)
saveRDS(TDS.ran,paste0("GM12878_",chr,"_random_TAD_detection_score.rds"))


TDS.ran.list <- list.files(path = "detection_score",
                           pattern = "_random_TAD_detection_score.rds$",
                           full.names = T)
names(TDS.ran.list) <- gsub("GM12878_","",basename(TDS.ran.list)) %>%
  gsub("_random_TAD_detection_score.rds","",.)

TDS.ran <- lapply(TDS.ran.list, function(x) readRDS(x))

TDS.ran.cell.df <- lapply(paste0("chr",c(1:22,"X")), function(chr){
  df <- lapply(TDS.ran[[chr]], function(x) rowMeans(x$DS)) %>% 
    Reduce(rbind,.) %>% as.data.frame() %>%
    `rownames<-`(paste0("R",1:100)) %>% 
    mutate(group = 1:100) %>%
    melt(id.vars = "group")
  df$chr <- chr
  return(df)
}) %>% Reduce(rbind,.)

TDS.ran.cell.df$chr <- factor(TDS.ran.cell.df$chr,
                              levels = paste0("chr",c(1:22,"X")))

TDS.ran.cell.df.m <- TDS.ran.cell.df %>% 
  group_by(variable,chr) %>%
  summarise(meanDS=mean(value))

TDS.ran.cell.df.m <- left_join(TDS.ran.cell.df.m,
                               meta[,c("rename2","filtered_pass_count","library2")],
                               by = c("variable" = "rename2"))


# Normalized TDS ----------------------------------------------------------

TAD_DS_norm <- lapply(paste0("chr",c(1:22,"X")), function(chr){
  mtx.used <- TDS.mtx[,grep(paste0(chr,"_"),colnames(TDS.mtx))]
  mean.used <- TDS.ran.cell.df.m[TDS.ran.cell.df.m$chr == chr,]
  print(identical(as.character(mean.used$variable),rownames(mtx.used)))
  mtx.norm <- mtx.used - mean.used$meanDS
  return(mtx.norm)
}) %>% Reduce(cbind,.)
write.csv(TAD_DS_norm,"detection_score/Normalized_TAD_DS_mtx.csv")
TAD_DS_norm <- read.csv("detection_score/Normalized_TAD_DS_mtx.csv",row.names = 1)

meta.df <- read.csv("scHiC_all_metadata_rename_all.csv")
cell.used <- meta.df[meta.df$filtered_pass_count > 10000  &
                       meta.df$cell_line == "GM12878",]$rename2
TAD_DS_norm <- TAD_DS_norm[cell.used,]

TAD_DS_norm.m <- rowMeans(TAD_DS_norm) %>% as.data.frame()
colnames(TAD_DS_norm.m) <- "meanDS"
TAD_DS_norm.m$cell <- rownames(TAD_DS_norm.m)
TAD_DS_norm.m <- left_join(TAD_DS_norm.m,
                           meta[,c("rename2","filtered_pass_count","library2")],
                           by = c("cell" = "rename2"))

TAD_DS_norm.m2 <- colMeans(TAD_DS_norm) %>% as.data.frame() %>%
  `colnames<-`("meanDS")
TAD_DS_norm.m2$boundary <- rownames(TAD_DS_norm.m2)


pdf("detection_score/TAD_detection_score_histogram_220221110.pdf",
    width = 8,height = 4)
ggplot(TAD_DS_norm.m,aes(x=meanDS)) +
  geom_histogram(color = "black",fill = "grey") +
  geom_vline(xintercept = 0,linetype = "dashed",color = "red") +
  xlab("Normalized Detection score") +
  theme_bw() + ylab("Number of cells") +
  theme(panel.grid = element_blank()) +
  ggplot(TAD_DS_norm.m2,aes(x=meanDS)) +
  geom_histogram(color = "black",fill = "grey") +
  geom_vline(xintercept = 0,linetype = "dashed",color = "red") +
  xlab("Normalized Detection score") +
  theme_bw() + ylab("Number of TADs") +
  theme(panel.grid = element_blank())
dev.off()

write.csv(TAD_DS_norm.m,"detection_score/TAD_mean_territories_detection_score.csv")
TAD_DS_norm.m <- read.csv("detection_score/TAD_mean_territories_detection_score.csv",
                          row.names = 1)

# mean norm TDS vs boundary strength --------------------------------------

TAD_DS_norm.m2$bd_strength <- insul.used$boundary_strength_500000
TAD_DS_norm.m2$insul_score <- insul.used$log2_insulation_score_500000
TAD_DS_norm.m2$sdDS <- TAD_DS_norm %>% apply(., 2, sd)
TAD_DS_norm.m2$idx <- 1:nrow(TAD_DS_norm.m2)
write.csv(TAD_DS_norm.m2,"detection_score/boundary_mean_NormDS.csv")
TAD_DS_norm.m2 <- read.csv("detection_score/boundary_mean_NormDS.csv",
                           row.names = 1)

cor.value <- round(cor(TAD_DS_norm.m2$meanDS,
                       TAD_DS_norm.m2$bd_strength,
                       method = "spearman"),
                   digits = 2)

library(ggpointdensity)
pdf("detection_score/Mean_nomalized_TAD_DS_boundary_strength_correlation_20221110.pdf",width = 6,height = 5)
ggplot(TAD_DS_norm.m2,aes(x=meanDS,y=bd_strength)) +
  geom_pointdensity(adjust = 0.05,size = 1) +
  geom_smooth(method = "lm") +
  ggplot2::annotate("text",label = paste0("sperman's r = ",cor.value),
           x = -0.01, y = 3) +
  theme_bw() + xlab("mean normalized TAD detection score") +
  ylab("ensemble boundary strength") +
  scale_color_viridis_c() +
  theme(panel.grid = element_blank())
dev.off()

TAD_DS_norm.m2$pos_ratio <- apply(TAD_DS_norm, 2, function(x) sum(x>0)/length(x))
TDS_CV <- apply(TAD_DS_norm, 2, function(x) sd(x)/mean(x))
TDS_pos <- apply(TAD_DS_norm, 2, function(x) sum(x>0)/length(x))

pdf("detection_score/positive_ratio_of_norm_TAD_DS_20221110.pdf",
    width = 5,height = 4.5)
ggplot(TAD_DS_norm.m2,aes(x=pos_ratio*100)) +
  geom_density() +
  geom_vline(xintercept = mean(TDS_pos)*100,
             linetype = "dashed",color = brewer.pal(9,"Set1")[1]) +
  geom_vline(xintercept = median(TDS_pos)*100,
             linetype = "dashed",color = brewer.pal(9,"Set1")[2]) +
  ggplot2::annotate("text",label = "mean = 59.0%",
                    x = 30, y = 0.035) +
  ggplot2::annotate("text",label = "median = 59.5%",
                    x = 30, y = 0.03) +
  theme_bw() + xlab("percentage of cells %") +
  ylab("density") + 
  labs(title = "Percentage of cells with positive TAD detection scores per TAD boundary") +
  theme(panel.grid = element_blank())
dev.off()

# summary -----------------------------------------------------------------

insul.used <- read.table("DSG_1.5mM/GM_1.5mM_all_merged_chronly_50k_res_500k_wins_insul_score.txt",
                         header = T)
TAD_DS_norm.m2 <- read.csv("detection_score/boundary_mean_NormDS.csv",
                           row.names = 1)
TAD.m <- left_join(insul.used,TAD_DS_norm.m2,by = c("tid" = "boundary"))
write.csv(TAD.m,"detection_score/boundary_mean_NormDS_summary.csv")

TAD_DS_norm.m <- read.csv("detection_score/TAD_mean_territories_detection_score.csv",
                          row.names = 1)
meta.df <- read.csv("scHiC_all_metadata_rename_all.csv")

cell.m <- left_join(TAD_DS_norm.m,
                     meta.df[,c("rename2","cell_name","group2")],
                     by = c("cell" = "rename2"))
write.csv(cell.m,"detection_score/cell_mean_TAD_NormDS_summary.csv")

# shuffle times test ------------------------------------------------------

lapply(TDS.ran, function(x){
  length(x[[1]]$boundary)
}) %>% unlist()

TDS.ran.cum <- TDS.ran.cell.df %>% 
  group_by(variable,chr) %>%
  mutate(cum = cumsum(value))

TDS.ran.cum$mean <- TDS.ran.cum$cum/TDS.ran.cum$group
TDS.ran.cum.m <- TDS.ran.cum %>%
  group_by(group,chr) %>%
  summarise(means = mean(mean))

TDS.ran.cum.m$chr <- factor(TDS.ran.cum.m$chr,
                            levels = paste0("chr",c(1:22,"X")))

TDS.ran.cum.m <- TDS.ran.cum.m %>%
  group_by(chr) %>%
  mutate(dev = means - last(means)) %>%
  mutate(norm_dev = dev/last(means))

pdf("figure/NM_revision/TDS_shuffle_deviations_test.pdf",
    width = 5.5,height = 4)
ggplot(TDS.ran.cum.m,
       aes(x = group, y = norm_dev*100)) +
  geom_line(aes(color = chr)) + theme_bw() +
  geom_hline(yintercept = c(5,1,-5,-1),
             linetype = "dashed") +
  xlab("Number of shuffle times") +
  labs(title = "TAD detection scores") +
  ylab("Percentage of deviations %") +
  scale_y_continuous(breaks = c(-10,-5,-1,0,1,5,10)) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,vjust = 0.5))
dev.off()

# TDSs vs CTCF strength ---------------------------------------------------

TAD_DS_norm.m2 <- read.csv("detection_score/boundary_mean_NormDS.csv",
                           row.names = 1)
TAD_DS_norm.m2$tid <- rownames(TAD_DS_norm.m2)
insul.used <- read.table("DSG_1.5mM/GM_1.5mM_all_merged_chronly_50k_res_500k_wins_insul_score.txt",
                         header = T)

CTCF.df <- read.table("ref/ENCODE_GM12878_RegionDB/GM12878/regions/GM12878_CTCF_ENCFF796WRU.bed")
colnames(CTCF.df)[1:3] <- c("chrom","start","end")

library(gUtils)

CTCF.gr <- dt2gr(CTCF.df)
insul.gr <- dt2gr(insul.used)

ov <- insul.gr %*% CTCF.gr

ov.df <- as.data.frame(ov) %>%
  group_by(tid) %>%
  summarise(CTCF = mean(V7))

TAD_DS_norm.m2 <- left_join(TAD_DS_norm.m2,ov.df)
TAD_DS_norm.m2[is.na(TAD_DS_norm.m2$CTCF),]$CTCF <- 0

cor(TAD_DS_norm.m2$meanDS, TAD_DS_norm.m2$CTCF)
cor(TAD_DS_norm.m2$meanDS, TAD_DS_norm.m2$CTCF,
    method = "spearman")

library(ggsignif)

pdf("figure/NM_revision/Norm_TDS_vs_CTCF_boxplot.pdf",width = 4,height = 4)
ggplot(TAD_DS_norm.m2,aes(x = factor(CTCF>0), y = meanDS)) +
  geom_boxplot(aes(fill = factor(CTCF>0)),
               outlier.size = 0.5) +
  geom_signif(
    comparisons = list(c("TRUE", "FALSE")),
    textsize = 6,
    map_signif_level = T
  ) + xlab("") + ylab("Average TAD detection scores") +
  theme_bw() + labs(title = "contains CTCF sites") +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,vjust = 0.5)) +
  scale_fill_manual("",values = brewer.pal(8,"Set2")[c(8,7)])
dev.off()
