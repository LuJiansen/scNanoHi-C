library(GenomicRanges)
library(RColorBrewer)
library(gUtils)
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)


# territories detection scores -------------------------------------------------------------

read_chr <- function(chr){
  df.list <- list.files(path = "detection_score",
                        pattern = paste0("_GRCh38_PG_1000000_",chr,"_mtx.gz"),
                        full.names = T)
  names(df.list) <- gsub("MboI_","",basename(df.list)) %>%
    gsub(paste0("_GRCh38_PG_1000000_",chr,"_mtx.gz"),"",.)
  df <- lapply(df.list, function(x) read.table(x))
  return(df)
}

terri_detection_score <- function(mtx,chr1.bin,chr2.bin){
  AA = nrow(mtx[mtx$V1 %in% chr1.bin & mtx$V2 %in% chr1.bin,])
  BB = nrow(mtx[mtx$V1 %in% chr2.bin & mtx$V2 %in% chr2.bin,])
  AB = nrow(mtx[mtx$V1 %in% chr1.bin & mtx$V2 %in% chr2.bin,])
  BA = nrow(mtx[mtx$V1 %in% chr2.bin & mtx$V2 %in% chr1.bin,])
  Eintra = length(chr1.bin)^2 + length(chr2.bin)^2
  Einter = length(chr1.bin)*length(chr2.bin)*2
  DS = (AA + BB)/Eintra - (AB + BA)/Einter
  return(DS)
}


## chr1-chr2 test ----------------------------------------------------------

inter.list <- list.files(path = "detection_score",
                         pattern = "_GRCh38_PG_1000000_chr1_chr2_mtx.gz$",
                         full.names = T)
names(inter.list) <- gsub("MboI_","",basename(inter.list)) %>%
  gsub("_GRCh38_PG_1000000_chr1_chr2_mtx.gz","",.)
inter.df <- lapply(inter.list, function(x){
  print(x)
  try(read.table(x))
})



chr1.df <- read_chr("chr1")
chr2.df <- read_chr("chr2")
length(chr2.df)

identical(names(chr1.df),names(chr2.df))
identical(names(chr1.df),names(inter.df))

ter.df <- lapply(names(inter.df), function(x){
  Reduce(rbind,list(chr1.df[[x]],chr2.df[[x]],inter.df[[x]]))
})
names(ter.df) <- names(inter.df)

bins <- read.table("detection_score/chronly_sort_1M_bins.txt")
bins$id <- 1:nrow(bins)

bins.used <- bins[bins$V1 %in% c("chr1","chr2"),]
chr1.bin <- bins.used[bins.used$V1 == "chr1",]$id
chr2.bin <- bins.used[bins.used$V1 == "chr2",]$id



ter_DS <- lapply(ter.df, function(mtx){terri_detection_score(mtx,chr1.bin,chr2.bin)})
plot(density(unlist(ter_DS)), main = "chr1-chr2 territory detection score")



## all chromosome ----------------------------------------------------------

bins <- read.table("detection_score/chronly_sort_1M_bins.txt")
bins$id <- 1:nrow(bins)


mtx.list <- list.files(path = "detection_score",
                       pattern = "_GRCh38_PG_1000000_all_mtx.gz$",
                       full.names = T)
names(mtx.list) <- gsub("MboI_","",basename(mtx.list)) %>%
  gsub("_GRCh38_PG_1000000_all_mtx.gz","",.)
mtx.df <- lapply(mtx.list, function(x){
  print(x)
  try(read.table(x))
})

chr.bin <- lapply(unique(bins$V1), function(x){
  bins[bins$V1 == x,]$id
}) %>% `names<-`(unique(bins$V1))

chr_pairs <- combn(unique(bins$V1)[c(1:22,24)],2) %>% t() %>% split(.,1:nrow(.))

lapply(chr_pairs, function(x){
  terri_detection_score(mtx.df$GM_24.P94B1,chr.bin[[x[1]]],chr.bin[[x[2]]])
}) %>% unlist()


ter_DS <- pbmclapply(mtx.df, function(mtx){
  lapply(chr_pairs, function(x){
    terri_detection_score(mtx,chr.bin[[x[1]]],chr.bin[[x[2]]])
  }) %>% unlist()
}) %>% Reduce(rbind,.)
rownames(ter_DS) <- names(mtx.df)
colnames(ter_DS) <- lapply(chr_pairs, function(x){paste0(x,collapse = "_")}) %>% unlist()

write.csv(ter_DS,"detection_score/GM12878_territories_detection_score.csv")

# random territory score --------------------------------------------------

chr.bin.len <- lapply(chr.bin, length)
bin.used <- bins[bins$V1 != "chrM",]

ter_DS_ran <- list()

run_terDS_random <- function(i,mc1 = 10, mc2 = 10){
  message("run ",i," ...")
  df  <- pbmclapply(mtx.df, function(mtx){
    mclapply(chr_pairs, function(x){
      rbin <- sample(bin.used$id, (chr.bin.len[[x[1]]] + chr.bin.len[[x[2]]]))
      bin1 <- head(rbin,chr.bin.len[[x[1]]])
      bin2 <- tail(rbin,chr.bin.len[[x[2]]])
      terri_detection_score(mtx,bin1,bin2)
    },mc.cores = mc1) %>% unlist()
  },mc.cores = mc2) %>% Reduce(rbind,.)
  rownames(df) <- names(mtx.df)
  colnames(df) <- lapply(chr_pairs,
                         function(x){paste0(x,collapse = "_")}) %>% unlist()
  return(df)
}

for (i in 1:100) {
  ter_DS_ran[[i]] <- run_terDS_random(i)
}

saveRDS(ter_DS_ran,"detection_score/territory_DS_random_100.rds")


# visualizations ----------------------------------------------------------

ter_DS_ran <- readRDS("detection_score/territory_DS_random_100.rds")

ter_DS_ran.m <- Reduce(`+`, ter_DS_ran)
ter_DS_ran.m <- ter_DS_ran.m/100

ter_DS <- read.csv("detection_score/GM12878_territories_detection_score.csv",row.names = 1)

identical(colnames(ter_DS),colnames(ter_DS_ran.m))
colnames(ter_DS) %in% colnames(ter_DS_ran.m)
identical(rownames(ter_DS),rownames(ter_DS_ran.m))

ter_DS.norm <- ter_DS - ter_DS_ran.m[,colnames(ter_DS)]

meta.df <- read.csv("scHiC_all_metadata_rename_all.csv")
cell.used <- meta.df[meta.df$filtered_pass_count > 10000  &
                       meta.df$cell_line == "GM12878",]$rename2
ter_DS.norm <- ter_DS.norm[cell.used,]

ter_DS.df <- data.frame(mean = rowMeans(ter_DS.norm),
                        cell = rownames(ter_DS.norm))
ter_DS.df <- cbind(ter_DS.df,ann_row[rownames(ter_DS.norm),])
ter_DS.df$counts <- meta[rownames(ter_DS.norm),]$filtered_pass_count
ter_DS.df2 <- data.frame(mean = colMeans(ter_DS.norm),
                         pair = colnames(ter_DS.norm))
summary(ter_DS.df$mean)

write.csv(ter_DS.df,"detection_score/cell_mean_normalized_territories_detection_score.csv")
write.csv(ter_DS.df2,"detection_score/ter_mean_normalized_territories_detection_score.csv")

pdf("detection_score/territories_detection_score_histogram_20221110.pdf",
    width = 8,height = 4)
ggplot(ter_DS.df,aes(x=mean)) +
  geom_histogram(color = "black",fill = "grey") +
  geom_vline(xintercept = 0,linetype = "dashed",color = "red") +
  xlab("Detection score") +
  theme_bw() + ylab("Number of cells") +
  theme(panel.grid = element_blank()) +
  ggplot(ter_DS.df2,aes(x=mean)) +
  geom_histogram(color = "black",fill = "grey") +
  geom_vline(xintercept = 0,linetype = "dashed",color = "red") +
  xlab("Detection score") +
  theme_bw() + ylab("Number of chromosome pairs") +
  theme(panel.grid = element_blank())
dev.off()


# summary -----------------------------------------------------------------

ter_DS.df <- read.csv("detection_score/cell_mean_normalized_territories_detection_score.csv")
ter_DS.df2 <- read.csv("detection_score/ter_mean_normalized_territories_detection_score.csv")

meta.df <- read.csv("scHiC_all_metadata_rename_all.csv")

cell.m <- left_join(ter_DS.df,
                    meta.df[,c("rename2","cell_name","group2")],
                    by = c("cell" = "rename2"))
write.csv(cell.m,"detection_score/cell_mean_normalized_territories_detection_score_summary.csv")



# shuffle times test ------------------------------------------------------

ter_DS_ran <- readRDS("detection_score/territory_DS_random_100.rds")

lapply(ter_DS_ran, dim)

ter_DS_ran.m <- Reduce(`+`, ter_DS_ran)
ter_DS_ran.m <- ter_DS_ran.m/100

ter_DS_ran.df <- lapply(1:100, function(i){
  colMeans(ter_DS_ran[[i]]) %>% 
    as.data.frame() %>%
    `colnames<-`("cell_mean") %>%
    mutate(group = i,
           ID = rownames(.))
}) %>% Reduce(rbind,.)

ter_DS_ran.df <- ter_DS_ran.df %>%
  group_by(ID) %>%
  mutate(cum = cumsum(cell_mean)) %>%
  mutate(cummean = cum/group) %>%
  mutate(dev = cummean - last(cummean)) %>%
  mutate(norm_dev = dev/last(cummean))

pdf("figure/NM_revision/territories_shuffle_deviations_test.pdf",
    width = 4,height = 4)
ggplot(ter_DS_ran.df,aes(x=group,y=norm_dev*100)) +
  geom_line(aes(group = ID),color = "grey") +
  geom_hline(yintercept = c(0.5,0.1,-0.5,-0.1),
             linetype = "dashed") +
  theme_bw() +
  xlab("Number of shuffle times") +
  labs(title = "Territories detection scores") +
  ylab("Percentage of deviations %") +
  scale_y_continuous(breaks = c(0.5,0.1,0,-0.5,-0.1)) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,vjust = 0.5))
dev.off()


