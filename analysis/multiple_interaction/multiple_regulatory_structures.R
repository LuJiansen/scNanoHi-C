library(GenomicRanges)
library(RColorBrewer)
library(gUtils)
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

prefix = "GM12878"

# load 20k bins ----------------------------------------------------------------

bins_20k <- read.table("case/chronly_sort_20k_bins.txt",header = F)
colnames(bins_20k) <- c("chrom","start","end")
bins_20k$id <- 1:(nrow(bins_20k))
bins_20k.gr <- gr0(bins_20k)

# functions ---------------------------------------------------------------

gr0 <- function(df,keep = T,zero = T){
  makeGRangesFromDataFrame(df,starts.in.df.are.0based = zero,
                           keep.extra.columns = keep)
}

multi_enh_reads <- function(gr,gene,range = gene.range,
                            pro = gene_pos.flt){
  if (gene %in% pro$gene) {
    res <- ov_table(gr,range[[gene]]) %>% 
      as.data.frame() %>% group_by(read_idx) %>%
      mutate(nbin = n_distinct(query.id))
    res$id <- range[[gene]][res$query.id]$id
    
    res <- res[res$nbin > 2,]
    pro.loc <- pro[pro$gene == gene,]$id
    read.used <- res[res$id %in% pro.loc,]$read_idx
    res <- res[res$read_idx %in% read.used,]
  } else{
    stop("Gene not in list!")
  }
  return(res)
}

multi_pro_reads <- function(gr,eid,range,enh){
  if (eid %in% unique(as.character(enh$uniq_eid))) {
    res <- ov_table(gr,range[[eid]]) %>% 
      as.data.frame() %>% group_by(read_idx) %>%
      mutate(nbin = n_distinct(query.id))
    res$id <- range[[eid]][res$query.id]$id
    res <- res[res$nbin > 2,]
    enh.loc <- unique(enh[enh$uniq_eid == eid,]$binid)
    read.used <- res[res$id %in% enh.loc,]$read_idx
    res <- res[res$read_idx %in% read.used,]
  } else{
    stop("Enhncer not in list!")
  }
  return(res)
}

add_multiway <- function(gr,filter = T){
  reads.sum <- table(gr$read_idx) %>% as.data.frame()
  if (filter) {
    cat("Numbe of singleton:",
        nrow(c),
        "\n")
    gr <- gr[gr$read_idx %in% reads.sum[reads.sum$Freq > 1,"Var1"],]
  }
  gr$multiway <- FALSE
  gr[gr$read_idx %in% reads.sum[reads.sum$Freq > 2,"Var1"],]$multiway <- TRUE
  return(gr)
}

gr2bed <- function(gr,color){
  names(gr) <- NULL
  test_bed <- as.data.frame(gr) %>%
    arrange(gene,cell,read_idx,id)
  test_bed$ID <- paste(test_bed$gene,test_bed$read_idx,sep = ":")
  df <- test_bed %>% group_by(seqnames,ID,gene,strand) %>%
    dplyr::summarise(chromStart = min(start),
                     chromEnd = max(end),
                     blockCount = length(ID),
                     blockSizes = paste(width,collapse = ","),
                     blockStarts = paste(start - min(start),collapse = ","))
  df <- df %>% arrange(seqnames,chromStart,gene)
  df$color <- factor(df$gene,levels = unique(df$gene))
  nbs <- length(levels(df$color))
  levels(df$color) <- c(rep(color,floor(nbs/length(color))),
                        color[1:(nbs%%length(color))])
  df$strand <- "."
  df <- df %>% arrange(ID)
  df_bed <- df[,c("seqnames","chromStart","chromEnd","ID",
                  "blockCount","strand","chromStart","chromEnd",
                  "color","blockCount","blockSizes","blockStarts")]
  return(df_bed)
}

# One promoter vs multiple enhancers --------------------------------------
## ABC annotations ---------------------------------------------------------

ABC_pos <- read.table(paste0("ABC_predictions/",prefix,"/",prefix,
                             ".PositivePredictions.hg38.bed"))
gene_pos <- read.table(paste0("ABC_predictions/",prefix,"/",prefix,
                              ".GeneSummary.hg38.bed"))
colnames(ABC_pos) <- c("chr","start","end","gene","ABC_score","strand")
colnames(gene_pos) <- c("chr","start","end","gene")
gene_pos <- gene_pos[gene_pos$chr %in% paste0("chr",c(1:22,"X","Y")),]

summary((ABC_pos$end - ABC_pos$start))

ABC_pos <- ABC_pos %>% 
  arrange(gene,start) %>% 
  group_by(gene) %>%
  mutate(nenh = length(gene),
         enh_dist = max(max(start) - min(end),0),
         eid = paste0(unique(gene),".e",1:length(gene)))

write.table(ABC_pos,paste0("ABC_predictions/",prefix,"/",prefix,
                           ".PositivePredictions.hg38.ann.xls"),
            row.names = F, quote = F, sep = '\t')

ABC_pos.gr <- gr0(ABC_pos)

gene_pos.gr <- gr0(gene_pos)
gene_pos.ov <- findOverlaps(bins_20k.gr,gene_pos.gr)
gene_pos.ann <- data.frame(id = gene_pos.ov@from,
                                gene = gene_pos.gr[gene_pos.ov@to]$gene,
                                name = paste0(gene_pos.gr[gene_pos.ov@to]$gene,".p"),
                                type = "promtoer")
ABC_pos.ov <- findOverlaps(bins_20k.gr,ABC_pos.gr)
ABC_pos.ann <- data.frame(id = ABC_pos.ov@from,
                               gene = ABC_pos.gr[ABC_pos.ov@to]$gene,
                               name = ABC_pos.gr[ABC_pos.ov@to]$eid,
                               type = "enhancer")
ABC_pos.ann.gr <- gr0(ABC_pos.ann[,c("id","gene","name")] %>%
                             group_by(name) %>%
                             mutate(start = min(id),
                                    end = max(id),
                                    seqnames = "bins"),
                           zero = F)
ABC_pos.ann.gr <- gr.reduce(ABC_pos.ann.gr,by = "gene")
sum(table(ABC_pos.ann.gr$gene) > 1)

# keep genes with at least two enhancer clusters
ABC_pos.flt <- as.data.frame(ABC_pos.ann.gr) %>%
  group_by(gene) %>%
  mutate(nbin = length(gene),
         name = paste0(unique(gene),".e",1:length(gene))) %>%
  filter(nbin > 1 & gene %in% gene_pos.ann$gene)

# whether enhancer in the same bin with promoter ? (at least 2 enhancer cluster excluding promoter bin)
ABC_pos.flt <- ABC_pos.flt %>% group_by(name) %>%
  mutate(promoter = gene_pos.ann[gene_pos.ann$gene == gene,]$id %in% (start:end))

gene.used <- names(which(table(ABC_pos.flt[!ABC_pos.flt$promoter,]$gene) > 1))
ABC_pos.flt <- ABC_pos.flt[ABC_pos.flt$gene %in% gene.used,]

gene_pos.flt <- gene_pos.ann[gene_pos.ann$gene %in% ABC_pos.flt$gene,]

write.table(ABC_pos.flt,paste0("case/",prefix,"_ABC_model_multi_enhancer_20k_bins_coordinates.xls"),
            row.names = F,col.names = T,sep = '\t',quote = F)
write.table(gene_pos.flt,paste0("case/",prefix,"_ABC_model_used_promoters_20k_bins_coordinates.xls"),
            row.names = F,col.names = T,sep = '\t',quote = F)

gene.range <- lapply(gene_pos.flt$gene, function(x){
  id.used <- c(as.integer(gene_pos.flt[gene_pos.flt$gene == x, "id"]),
               apply(ABC_pos.flt[ABC_pos.flt$gene == x,],
                     1, function(x) x[2]:x[3]) %>% unlist()) %>% 
    sort() %>% unique()
  print(id.used)
  gr <- bins_20k.gr[id.used]
  return(gr)
})
names(gene.range) <- gene_pos.flt$gene

saveRDS(gene.range,paste0("case/",prefix,"_ABC_used_regions.rds"))

unique(gene_pos.flt$gene) %>% length()
unique(ABC_pos.flt$gene) %>% length()


## find mutiEnh gene concatemers ------------------------------------------------------
ABC_pos.flt <- read.table(paste0("case/",prefix,"_ABC_model_multi_enhancer_20k_bins_coordinates.xls"),header = T)
gene_pos.flt <- read.table(paste0("case/",prefix,"_ABC_model_used_promoters_20k_bins_coordinates.xls"),header = T)
gene.range <- readRDS(paste0("case/",prefix,"_ABC_used_regions.rds"))

mono.gr <- readRDS(paste0("synergy/",prefix,"_all_gr.rds"))
mono.gr <- add_multiway(mono.gr)
gene.range.all <- GRangesList(gene.range) %>% 
  unlist() %>% GenomicRanges::reduce()

mono_flt.gr <- subsetByOverlaps(mono.gr,gene.range.all)
mono_flt.gr <- mono_flt.gr[mono_flt.gr$multiway,]
mono_flt.sum <- table(mono_flt.gr$read_idx) %>% as.data.frame()
mono_flt.gr <- mono_flt.gr[mono_flt.gr$read_idx %in% 
                             mono_flt.sum[mono_flt.sum$Freq > 2,"Var1"],]

multiEnh <- pbmclapply(gene_pos.flt$gene,
                            function(x) multi_enh_reads(mono_flt.gr,x,
                                                        pro = gene_pos.flt),
                            mc.cores = 30)
names(multiEnh) <- gene_pos.flt$gene
saveRDS(multiEnh,paste0("case/",prefix,"_multiEnh_concatemers.rds"))

multiEnh <- readRDS(paste0("case/",prefix,"_multiEnh_concatemers.rds"))

for (i in names(multiEnh)) {
  multiEnh[[i]]$gene <- i
}

multiEnh <- multiEnh[names(which((lapply(multiEnh, nrow) %>% unlist()) > 0))]

multiEnh.mg <- Reduce(rbind,multiEnh)
multiEnh.mg <- multiEnh.mg %>% group_by(gene) %>%
  mutate(nreads = n_distinct(read_idx))
multiEnh.mg <- multiEnh.mg %>% 
  arrange(desc(nreads),gene,cell,read_idx,id)

saveRDS(multiEnh.mg,paste0("case/",prefix,"_multiEnh_concatemers_merged_table.rds"))


## merge adjacent bins  ------------------------------------------------------------------
multiEnh.mg <- readRDS(paste0("case/",prefix,"_multiEnh_concatemers_merged_table.rds"))
multiEnh.gr <- cbind(multiEnh.mg,
                     as.data.frame(bins_20k.gr[multiEnh.mg$id])[-6]) %>%
  gr0(.,zero = F)

multiEnh.gr$ID <- paste(multiEnh.gr$gene,
                             multiEnh.gr$read_idx,sep = ":")
multiEnh.rd.gr <- gr.reduce(multiEnh.gr,by = "ID")
mcols(multiEnh.rd.gr) <- DataFrame(mcols(multiEnh.rd.gr) %>% 
                                          as.data.frame() %>% 
                                          group_by(ID) %>%
                                          mutate(nbin = length(ID)))
multiEnh.rd.gr[multiEnh.rd.gr$nbin > 2]$read_idx %>% unique()
multiEnh.rd.gr[multiEnh.rd.gr$nbin > 2]$gene %>% unique()
multiEnh.rd.gr <- multiEnh.rd.gr[multiEnh.rd.gr$nbin > 2]
mcols(multiEnh.rd.gr) <- DataFrame(mcols(multiEnh.rd.gr) %>% 
                                          as.data.frame() %>% 
                                          group_by(gene) %>%
                                          mutate(nreads = n_distinct(read_idx)))
saveRDS(multiEnh.rd.gr,paste0("case/",prefix,"_multiEnh_concatemers_reduced.rds"))
multiEnh.rd.gr <- readRDS(paste0("case/",prefix,"_multiEnh_concatemers_reduced.rds"))



## output bed files --------------------------------------------------------

multiEnh.rd.bed <- gr2bed(multiEnh.rd.gr,brewer.pal(10,"Set3")[-2])
write.table(multiEnh.rd.bed,
            paste0("case/",prefix,"_multiEnh_concatemers_20k_reduced.bed"),
            sep = "\t",quote = F,row.names = F,col.names = F)


## summary -----------------------------------------------------------------

multiEnh.rd.sum <- multiEnh.rd.gr %>% 
  as.data.frame() %>%
  group_by(gene) %>%
  summarise(nreads = n_distinct(read_idx),
            ncell = n_distinct(cell)) %>%
  arrange(desc(ncell),desc(nreads))
write.table(multiEnh.rd.sum,paste0("case/",prefix,"_multiEnh_concatemers_reduced_summary.xls"),
            row.names = F,sep = '\t',quote = F)

# One enhancer to multiple promoters --------------------------------------

## ABC annotations -----------------------------------------------------------------

ABC_pos <- read.table(paste0("ABC_predictions/",prefix,"/",prefix,
                      ".PositivePredictions.hg38.bed"))
colnames(ABC_pos) <- c("chr","start","end","gene","ABC_score","strand")
ABC_pos$chr <- factor(ABC_pos$chr,
                         levels = paste0("chr",c(1:22,"X","Y")))

ABC_pos <- ABC_pos %>% arrange(chr,start) %>%
  mutate(uniq_eid = paste(chr,start,end,sep = "_"))
ABC_pos$uniq_eid <- factor(ABC_pos$uniq_eid,
                              levels = unique(ABC_pos$uniq_eid))
levels(ABC_pos$uniq_eid) <- paste0("e",1:length(levels(ABC_pos$uniq_eid)))

ABC_pos <- ABC_pos %>% group_by(uniq_eid) %>%
  mutate(npro = n_distinct(gene))

write.csv(ABC_pos,paste0("case/",prefix,"_ABC_modelpostive_enhancer_summary.csv"))

# find enhancers regulated > 1 promoters
ABC_pos.flt <- ABC_pos[ABC_pos$npro > 1,]
unique(ABC_pos.flt$uniq_eid) %>% length()

# collapse genes into promoter clusters

gene_pos <- read.table(paste0("ABC_predictions/",prefix,"/",prefix,
                          ".GeneSummary.hg38.bed"))
colnames(gene_pos) <- c("chr","start","end","gene")
gene_pos <- gene_pos[gene_pos$chr %in% paste0("chr",c(1:22,"X","Y")),]

gene_pos.gr <- gr0(gene_pos)
gene_pos.ov <- findOverlaps(bins_20k.gr,gene_pos.gr)

gene.bin <- data.frame(id = gene_pos.ov@from,
                          gene = gene_pos.gr[gene_pos.ov@to]$gene)
gene.bin <- gene.bin %>% group_by(id) %>%
  mutate(ngenes = n_distinct(gene),
         genes = paste(gene,collapse = ","))
write.table(gene.bin,paste0("case/",prefix,"_ABC_model_20k_bins_promoter_cluters.xls"),
            row.names = F,quote = F,sep = '\t')


ABC_pos.flt <- left_join(ABC_pos.flt,gene.bin)
ABC_pos.flt <- ABC_pos.flt[!is.na(ABC_pos.flt$id),]
enh.ov <- findOverlaps(bins_20k.gr,gr0(unique(ABC_pos.flt[,c(1:3,7)])))
enh.bin <- data.frame(uniq_eid = unique(ABC_pos.flt$uniq_eid)[enh.ov@to],
                         binid = enh.ov@from)
ABC_pos.flt <- left_join(ABC_pos.flt,enh.bin)
ABC_pos.flt <- ABC_pos.flt %>% group_by(uniq_eid) %>%
  mutate(nbins = n_distinct(id))

# keep enhancers regulated at least 2 promoters bins
ABC_pos.flt$uniq_eid %>% unique() %>% length()
ABC_pos.flt[ABC_pos.flt$nbins > 1,]$uniq_eid %>% unique() %>% length()

ABC_pos.flt <- ABC_pos.flt[ABC_pos.flt$nbins > 1,]
ABC_pos.flt$uniq_eid %>% unique() %>% length()

# enhancer in promoter bins ?
ABC_pos.flt <- ABC_pos.flt %>% 
  group_by(uniq_eid) %>%
  mutate(enhInPro = sum(unique(binid) %in% id)) %>%
  mutate(noPro_nbins = nbins - enhInPro)

# keep enhancers regulated at least 2 promoters bins excluding itself bin
ABC_pos.flt[ABC_pos.flt$noPro_nbins > 1,]$uniq_eid %>% unique() %>% length()
ABC_pos.flt <- ABC_pos.flt[ABC_pos.flt$noPro_nbins > 1,]

multiPro.range <- lapply(as.character(unique(ABC_pos.flt$uniq_eid)),
                         function(x){
                           id.used <- unique(c(ABC_pos.flt[ABC_pos.flt$uniq_eid == x,]$id,
                                               ABC_pos.flt[ABC_pos.flt$uniq_eid == x,]$binid))
                           print(id.used)
                           return(id.used)
                         })
names(multiPro.range) <- as.character(unique(ABC_pos.flt$uniq_eid))

multiPro.range.rd <- lapply(names(multiPro.range), function(x){
  print(x)
  reduce(IRanges(start = multiPro.range[[x]]))
})
names(multiPro.range.rd) <- names(multiPro.range)
enh.used <- names(which((lapply(multiPro.range.rd, length) %>% unlist()) > 2))

# keep enhancers regulated at least 2 uncontinuous promoters bins
ABC_pos.flt <- ABC_pos.flt[ABC_pos.flt$uniq_eid %in% enh.used,]
write.table(ABC_pos.flt,paste0("case/",prefix,"_ABC_model_multiPro_enhancer_20k_bins_used.xls"),
            row.names = F,quote = F,sep = '\t')


multiPro.range.used <- lapply(as.character(unique(ABC_pos.flt$uniq_eid)),
                              function(x){
                                id.used <- unique(c(ABC_pos.flt[ABC_pos.flt$uniq_eid == x,]$id,
                                                    ABC_pos.flt[ABC_pos.flt$uniq_eid == x,]$binid))
                                print(id.used)
                                return(id.used)
                              })
names(multiPro.range.used) <- as.character(unique(ABC_pos.flt$uniq_eid))

multiPro.range.gr <- lapply(multiPro.range.used, function(x){
  print(x)
  gr <- bins_20k.gr[x]
  return(gr)
})

saveRDS(multiPro.range.gr,paste0("case/",prefix,"_ABC_multiPro_enh_20k_bins_range.rds"))


## find multiPro enhancer concatemers --------------------------------------------------

multiPro.range.gr <- readRDS(paste0("case/",prefix,"_ABC_multiPro_enh_20k_bins_range.rds"))
ABC_pos.flt <- read.table(paste0("case/",prefix,"_ABC_model_multiPro_enhancer_20k_bins_used.xls"),header = T)
gene.bin <- read.table(paste0("case/",prefix,"_ABC_model_20k_bins_promoter_cluters.xls"),
                          row.names = F,quote = F,sep = '\t')

multiPro.range.all <- GRangesList(multiPro.range.gr) %>% 
  unlist() %>% GenomicRanges::reduce()

flt.gr <- subsetByOverlaps(mono.gr,multiPro.range.all)
flt.gr <- flt.gr[flt.gr$multiway,]
flt.sum <- table(flt.gr$read_idx) %>% as.data.frame()
flt.gr <- flt.gr[flt.gr$read_idx %in% 
                   flt.sum[flt.sum$Freq > 2,"Var1"],]

multiPro <- pbmclapply(unique(as.character(ABC_pos.flt$uniq_eid)),
                       function(x) multi_pro_reads(flt.gr,
                                                   x,multiPro.range.gr,
                                                   ABC_pos.flt),
                       mc.cores = 30)
names(multiPro) <- unique(as.character(ABC_pos.flt$uniq_eid))
for (i in names(multiPro)) {
  multiPro[[i]]$uniq_eid <- i
}

multiPro <- multiPro[names(which((lapply(multiPro, nrow) %>% unlist()) > 0))]
saveRDS(multiPro,"case/",prefix,"_multiPro_concatemers.rds")

multiPro.mg <- Reduce(rbind,multiPro)
multiPro.mg <- multiPro.mg %>% group_by(uniq_eid) %>%
  mutate(nreads = n_distinct(read_idx))
multiPro.mg <- multiPro.mg %>% 
  arrange(desc(nreads),uniq_eid,cell,read_idx,id)

saveRDS(multiPro.mg,"case/",prefix,"_multiPro_concatemers_merged_table.rds")

## merge adjacent bins ------------------------------------------------------------------

multiPro.gr <- cbind(multiPro.mg,
                     as.data.frame(bins_20k.gr[multiPro.mg$id])[-6]) %>%
  gr0(.,zero = F)
multiPro.gr$gene <- multiPro.gr$uniq_eid
multiPro.gr$ID <- paste(multiPro.gr$gene,multiPro.gr$read_idx,sep = ":")
multiPro.rd.gr <- gr.reduce(multiPro.gr,by = "ID")
mcols(multiPro.rd.gr) <- DataFrame(mcols(multiPro.rd.gr) %>% 
                                     as.data.frame() %>% 
                                     group_by(ID) %>%
                                     mutate(nbin = length(ID)))
multiPro.rd.gr[multiPro.rd.gr$nbin > 2]$read_idx %>% unique() %>% length()
multiPro.rd.gr[multiPro.rd.gr$nbin > 2]$gene %>% unique() %>% length()
multiPro.rd.gr <- multiPro.rd.gr[multiPro.rd.gr$nbin > 2]
mcols(multiPro.rd.gr) <- DataFrame(mcols(multiPro.rd.gr) %>% 
                                     as.data.frame() %>% 
                                     group_by(gene) %>%
                                     mutate(nreads = n_distinct(read_idx)))
saveRDS(multiPro.rd.gr,"case/",prefix,"_multiPro_concatemers_reduced.rds")


## output bed files --------------------------------------------------------

multiPro.rd.bed <- gr2bed(multiPro.rd.gr,brewer.pal(10,"Set3")[-2])
write.table(multiPro.rd.bed, "case/",prefix,"_multiPro_concatemers_20k_reduced.bed",
            sep = "\t",quote = F,row.names = F,col.names = F)


## summary -----------------------------------------------------------------

multiPro.rd.sum <- multiPro.rd.gr %>% 
  as.data.frame() %>%
  group_by(uniq_eid) %>%
  summarise(nreads = n_distinct(read_idx),
            ncells = n_distinct(cell)) %>%
  arrange(desc(ncells),desc(nreads))
multiPro.rd.sum <- left_join(multiPro.rd.sum,enh.sum)
multiPro.rd.sum <- multiPro.rd.sum %>% 
  arrange(desc(ncells),desc(nreads),uniq_eid)

write.table(multiPro.rd.sum,"case/",prefix,"_multiPro_concatemers_reduced_summary.xls",
            row.names = F,sep = '\t',quote = F)
