library(dplyr)
library(gUtils)
library(rtracklayer)



# define functions --------------------------------------------------------

splitgi <- function(gi){
  gr1 <- gi@regions[gi@anchor1]
  gr2 <- gi@regions[gi@anchor2]
  meta <- mcols(gi)
  return(list(gr1 = gr1,
              gr2 = gr2,
              meta = meta))
}

unlistgi <- function(gil,add.meta = T, verbose = F){
  message("merge GInteractions")
  gr1 = gil$gr1
  gr2 = gil$gr2
  meta = gil$meta
  if (verbose) {
    print(gr1)
    print(gr2)
  }
  gi <- tryCatch(InteractionSet::GInteractions(gr1, gr2),
                 error = function(e) NULL)
  if (is.null(gi)) {
    warning("Construct GInteractions with non-identical seqlengths")
    gr1 = gr.fix(gr1, gr2)
    gr2 = gr.fix(gr2, gr1)
    gi <- InteractionSet::GInteractions(gr1, gr2)
  }
  if (verbose) {
    print(gi)
  }
  if (add.meta) {
    mcols(gi) <- meta
  }
  return(gi)
}

bedpe2gi <- function(df, loc = c(1,2,3,4,5,6),add.meta = T){
  df1 <- df[,loc[1:3]]
  df2 <- df[,loc[4:6]]
  colnames(df1) <- c("chr","start","end")
  colnames(df2) <- c("chr","start","end")
  
  gr1 <- dt2gr(df1)
  gr2 <- dt2gr(df2)
  gi <- tryCatch(InteractionSet::GInteractions(gr1, gr2), error = function(e) NULL)
  if (is.null(gi)) {
    warning("Construct GInteractions with non-identical seqlengths")
    gr1 = gr.fix(gr1, gr2)
    gr2 = gr.fix(gr2, gr1)
    gi <- InteractionSet::GInteractions(gr1, gr2)
  }
  
  if (add.meta) {
    mcols(gi) <- df[,-loc]
  }
  return(gi)
}

flankgi <- function(gi,length, reverse = F){
  gi.list <- splitgi(gi)
  if (reverse) {
    gi.list$gr1 <- gi.list$gr1 - length
    gi.list$gr2 <- gi.list$gr2 - length
  } else {
    gi.list$gr1 <- gi.list$gr1 + length
    gi.list$gr2 <- gi.list$gr2 + length
  }
  return(unlistgi(gi.list))
}

fix.ov <- function(query, subject){
  ov <- tryCatch(InteractionSet::findOverlaps(query, subject), error = function(e) NULL)
  if (is.null(ov)) {
    warning("findOverlaps applied to GInteractions with non-identical seqlengths")
    query = gr.fix(query, subject)
    subject = gr.fix(subject, query)
    ov <- InteractionSet::findOverlaps(query, subject)
  }
  return(ov)
}

addFilter <- function(df){
  df$range <- df$Start2 - df$Start1
  df$filter <- "others"
  df[df$range > 100000 &
       df$range < 1000000,]$filter <- "pass"
  print(nrow(df[df$filter == "pass",]))
  return(df)
}

liftOverbedpe <- function(gi,chain,add.meta = T,verbose = F){
  message("input GInteractions with ",length(gi)," records")
  gil <- splitgi(gi)
  gr1 <- gil$gr1
  gr1$ID <- 1:length(gr1)
  gr2 <- gil$gr2
  gr2$ID <- 1:length(gr2)
  meta <- gil$meta
  
  if (verbose) {
    print(gr1)
    print(gr2)
  }
  
  message("run liftOver ...")
  lo.test1 <- liftOver(gr1,chain = chain) %>% unlist()
  lo.test2 <- liftOver(gr2,chain = chain) %>% unlist()
  
  if (verbose) {
    print(lo.test1)
    print(lo.test2)
  }
  
  id.used <- intersect(names(which(table(lo.test1$ID) == 1)),
                       names(which(table(lo.test2$ID) == 1))) %>%
    as.numeric()
  
  message("find concordant pairs ...")
  lo1.used <- lo.test1[lo.test1$ID %in% id.used]
  lo2.used <- lo.test2[lo.test2$ID %in% id.used]
  
  if (verbose) {
    print(length(id.used))
    print(lo1.used)
    print(table(lo1.used$ID) %>% table())
    print(lo2.used)
    print(table(lo2.used$ID) %>% table())
  }
  
  message("generate GInteractions ...")
  if (identical(lo1.used$ID,lo2.used$ID)) {
    gil.out <- list(gr1 = lo1.used,
                    gr2 = lo2.used,
                    meta = meta[id.used,])
    gi.out <- unlistgi(gil.out, add.meta, verbose)
    message("output GInteractions with ",length(gi.out)," records")
    return(gi.out)
  } else {stop("length of liftOvered GRanges not identical!")}
}

resizebedpe <- function(gi,width = c(10001,10001),
                        fix = c("center","center")){
  gil <- splitgi(gi)
  gr1 <- resize(gil$gr1,width[1],fix = fix[1])
  gr2 <- resize(gil$gr2,width[1],fix = fix[1])
  gi.out <- unlistgi(list(gr1 = gr1,
                          gr2 = gr2,
                          meta = gil$meta))
  return(gi.out)
}

bedpe2lr <- function(bedpe,loc = 1:7){
  df <- bedpe[,loc]
  colnames(df) <- c("c1","s1","e1","c2","s2","e2","value")
  df <- df %>%
    mutate(t1 = paste0(c1,":",s1,"-",e1,",",value),
           t2 = paste0(c2,":",s2,"-",e2,",",value))
  lr <- rbind(df[,c("c1","s1","e1","t2")] %>%
                `colnames<-`(c("chr","start","end","target")),
              df[,c("c2","s2","e2","t1")] %>%
                `colnames<-`(c("chr","start","end","target")))
  return(lr)
}

ov_stat <- function(ov){
  df <- data.frame(nquery = ov@nLnode,
                   nsubject = ov@nRnode,
                   queryHit = length(unique(ov@from)),
                   subjectHit = length(unique(ov@to)))
  df <- df %>%
    mutate(precision = queryHit/nquery,
           recall = subjectHit/nsubject) %>%
    mutate(F1 = 2*(precision*recall)/(precision + recall))
  return(df)
}

# read and filter reference loops -----------------------------------------

ch = import.chain("F:/Tanglab/database/liftOver/hg19ToHg38.over.chain")

## GM12878 -----------------------------------------------------------------


### in situ Hi-C ------------------------------------------------------------

GM_10k.ref.df <- read.table("loop/Rao_GM_postprocessed_pixels_10000_with_motifs.bedpe",header = T)
GM_10k.ref.df$chr1 <- paste0("chr",GM_10k.ref.df$chr1)
GM_10k.ref.df$chr2 <- paste0("chr",GM_10k.ref.df$chr2)
GM_10k.ref <- bedpe2gi(GM_10k.ref.df)

### cohesin -----------------------------------------------------------------

GM.c <- list()
GM.c[["FitHiChIP"]] <- readxl::read_xlsx("ref/FitHiChIP/GM12878/Cohesin/Combined_Replicates_HIChIP_Loops/Table_GC-ALL.xlsx",
                                         sheet = 2, skip = 17)
GM.c[["MAPS"]] <- readxl::read_xlsx("ref/FitHiChIP/GM12878/Cohesin/Combined_Replicates_HIChIP_Loops/Table_GC-ALL.xlsx",
                                    sheet = 8, skip = 13)

GM.c[["FitHiChIP"]] <- addFilter(GM.c[["FitHiChIP"]])
GM.c[["MAPS"]] <- addFilter(GM.c[["MAPS"]])

GM.c.gr <- lapply(GM.c, function(x){
  bedpe2gi(x,loc = c(1,2,3,1,4,5))
})

GM.c.lo <- lapply(GM.c.gr, function(gi){liftOverbedpe(gi,ch,verbose = T)})
GM.c.lo <- lapply(GM.c.lo, resizebedpe)
GM.c.lo.flt <- lapply(GM.c.lo, function(x) x[x$filter == "pass"])
saveRDS(GM.c.lo,"ref/FitHiChIP/GM12878/Cohesin/Combined_Replicates_HIChIP_Loops/GM12878_cohesin_FitHiChIP_hg38_resize_filter.rds")


#### save bedpe (for IGV) ----------------------------------------------------

GM.c.df <- list()
GM.c.df[["FitHiChIP"]] <- as.data.frame(GM.c.lo[["FitHiChIP"]][GM.c.lo[["FitHiChIP"]]$filter == "pass"])[,c(1:3,6:8,19)]
GM.c.df[["FitHiChIP"]]$Q.Value <- -log10(GM.c.df[["FitHiChIP"]]$Q.Value)
GM.c.df[["MAPS"]] <- as.data.frame(GM.c.lo[["MAPS"]][GM.c.lo[["MAPS"]]$filter == "pass"])[,c(1:3,6:8,12)]
GM.c.df[["MAPS"]]$Q.value <- -log10(GM.c.df[["MAPS"]]$Q.value)

write.table(GM.c.df[["FitHiChIP"]],
            "ref/FitHiChIP/GM12878/Cohesin/Combined_Replicates_HIChIP_Loops/GM12878_cohesin_FitHiChIP_hg38_resize_filter.bedpe",
            row.names = F,col.names = F,sep = '\t',quote = F)
write.table(GM.c.df[["MAPS"]],
            "ref/FitHiChIP/GM12878/Cohesin/Combined_Replicates_HIChIP_Loops/GM12878_cohesin_MAPS_hg38_resize_filter.bedpe",
            row.names = F,col.names = F,sep = '\t',quote = F)


#### save longrange (for WashU) ----------------------------------------------

GM.c.lr <- list()
GM.c.lr[["FitHiChIP"]] <- bedpe2lr(GM.c.df[["FitHiChIP"]])
GM.c.lr[["MAPS"]] <- bedpe2lr(GM.c.df[["MAPS"]])
write.table(GM.c.lr[["FitHiChIP"]],
            "ref/FitHiChIP/used/GM12878_cohesin_FitHiChIP_hg38_resize_filter.lr",
            row.names = F,col.names = F,sep = '\t',quote = F)
write.table(GM.c.lr[["MAPS"]],
            "ref/FitHiChIP/used/GM12878_cohesin_MAPS_hg38_resize_filter.lr",
            row.names = F,col.names = F,sep = '\t',quote = F)
### H3K27ac -----------------------------------------------------------------

GM.h <- list()
GM.h[["FitHiChIP"]] <- readxl::read_xlsx("ref/FitHiChIP/GM12878/H3K27ac/Combined_Replicates_HIChIP_Loops/Table_GH-ALL.xlsx",
                                         sheet = 2, skip = 17)
GM.h[["MAPS"]] <- readxl::read_xlsx("ref/FitHiChIP/GM12878/H3K27ac/Combined_Replicates_HIChIP_Loops/Table_GH-ALL.xlsx",
                                    sheet = 8, skip = 13)

GM.h[["FitHiChIP"]] <- addFilter(GM.h[["FitHiChIP"]])
GM.h[["MAPS"]] <- addFilter(GM.h[["MAPS"]])

GM.h.gr <- lapply(GM.h, function(x){
  bedpe2gi(x,loc = c(1,2,3,1,4,5))
})

GM.h.lo <- lapply(GM.h.gr, function(gi){liftOverbedpe(gi,ch,verbose = T)})
GM.h.lo <- lapply(GM.h.lo, resizebedpe)
GM.h.lo.flt <- lapply(GM.h.lo, function(x) x[x$filter == "pass"])
saveRDS(GM.h.lo,"ref/FitHiChIP/GM12878/H3K27ac/Combined_Replicates_HIChIP_Loops/GM12878_H3K27ac_FitHiChIP_hg38_resize_filter.rds")

#### save bedpe (for IGV) ----------------------------------------------------

GM.h.df <- list()
GM.h.df[["FitHiChIP"]] <- as.data.frame(GM.h.lo[["FitHiChIP"]][GM.h.lo[["FitHiChIP"]]$filter == "pass"])[,c(1:3,6:8,19)]
GM.h.df[["FitHiChIP"]]$Q.Value <- -log10(GM.h.df[["FitHiChIP"]]$Q.Value)
GM.h.df[["MAPS"]] <- as.data.frame(GM.h.lo[["MAPS"]][GM.h.lo[["MAPS"]]$filter == "pass"])[,c(1:3,6:8,12)]
GM.h.df[["MAPS"]]$Q.value <- -log10(as.numeric(GM.h.df[["MAPS"]]$Q.value))

write.table(GM.h.df[["FitHiChIP"]],
            "ref/FitHiChIP/GM12878/H3K27ac/Combined_Replicates_HIChIP_Loops/GM12878_H3K27ac_FitHiChIP_hg38_resize_filter.bedpe",
            row.names = F,col.names = F,sep = '\t',quote = F)
write.table(GM.h.df[["MAPS"]],
            "ref/FitHiChIP/GM12878/H3K27ac/Combined_Replicates_HIChIP_Loops/GM12878_H3K27ac_MAPS_hg38_resize_filter.bedpe",
            row.names = F,col.names = F,sep = '\t',quote = F)

#### save longrange (for WashU) ----------------------------------------------

GM.h.lr <- list()
GM.h.lr[["FitHiChIP"]] <- bedpe2lr(GM.h.df[["FitHiChIP"]])
GM.h.lr[["MAPS"]] <- bedpe2lr(GM.h.df[["MAPS"]])
write.table(GM.h.lr[["FitHiChIP"]],
            "ref/FitHiChIP/used/GM12878_H3K27ac_FitHiChIP_hg38_resize_filter.lr",
            row.names = F,col.names = F,sep = '\t',quote = F)
write.table(GM.h.lr[["MAPS"]],
            "ref/FitHiChIP/used/GM12878_H3K27ac_MAPS_hg38_resize_filter.lr",
            row.names = F,col.names = F,sep = '\t',quote = F)


## K562 --------------------------------------------------------------------


### in situ Hi-C ------------------------------------------------------------

K562_10k.ref.df <- read.table("loop/K562_postprocessed_pixels_10000_with_motifs.bedpe",header = T)
K562_10k.ref <- bedpe2gi(K562_10k.ref.df)

### H3K27ac -----------------------------------------------------------------

K562.h <- list()
K562.h[["FitHiChIP"]] <- readxl::read_xlsx("ref/FitHiChIP/K562/H3K277ac/Combined_Replicates_HIChIP_Loops/Table_KH-ALL.xlsx",
                                           sheet = 2, skip = 20)
K562.h[["MAPS"]] <- readxl::read_xlsx("ref/FitHiChIP/K562/H3K277ac/Combined_Replicates_HIChIP_Loops/Table_KH-ALL.xlsx",
                                      sheet = 8, skip = 13)

K562.h[["FitHiChIP"]] <- addFilter(K562.h[["FitHiChIP"]])
K562.h[["MAPS"]] <- addFilter(K562.h[["MAPS"]])

K562.h.gr <- lapply(K562.h, function(x){
  bedpe2gi(x,loc = c(1,2,3,1,4,5))
})

K562.h.lo <- lapply(K562.h.gr, function(gi){liftOverbedpe(gi,ch,verbose = F)})
K562.h.lo <- lapply(K562.h.lo, resizebedpe)
K562.h.lo.flt <- lapply(K562.h.lo, function(x) x[x$filter == "pass"])
saveRDS(K562.h.lo,"ref/FitHiChIP/K562/H3K277ac/Combined_Replicates_HIChIP_Loops/K562_H3K27ac_hg38_resize_filter.rds")

#### save bedpe (for IGV) ----------------------------------------------------

K562.h.df <- list()
K562.h.df[["FitHiChIP"]] <- as.data.frame(K562.h.lo[["FitHiChIP"]][K562.h.lo[["FitHiChIP"]]$filter == "pass"])[,c(1:3,6:8,19)]
K562.h.df[["FitHiChIP"]]$Q.Value <- -log10(K562.h.df[["FitHiChIP"]]$Q.Value)
K562.h.df[["MAPS"]] <- as.data.frame(K562.h.lo[["MAPS"]][K562.h.lo[["MAPS"]]$filter == "pass"])[,c(1:3,6:8,12)]
K562.h.df[["MAPS"]]$Q.value <- -log10(as.numeric(K562.h.df[["MAPS"]]$Q.value))

write.table(K562.h.df[["FitHiChIP"]],
            "ref/FitHiChIP/K562/H3K277ac/Combined_Replicates_HIChIP_Loops/K562_H3K27ac_FitHiChIP_hg38_resize_filter.bedpe",
            row.names = F,col.names = F,sep = '\t',quote = F)
write.table(K562.h.df[["MAPS"]],
            "ref/FitHiChIP/K562/H3K277ac/Combined_Replicates_HIChIP_Loops/K562_H3K27ac_MAPS_hg38_resize_filter.bedpe",
            row.names = F,col.names = F,sep = '\t',quote = F)

#### save longrange (for WashU) ----------------------------------------------

K562.h.lr <- list()
K562.h.lr[["FitHiChIP"]] <- bedpe2lr(K562.h.df[["FitHiChIP"]])
K562.h.lr[["MAPS"]] <- bedpe2lr(K562.h.df[["MAPS"]])
write.table(K562.h.lr[["FitHiChIP"]],
            "ref/FitHiChIP/used/K562_H3K27ac_FitHiChIP_hg38_resize_filter.lr",
            row.names = F,col.names = F,sep = '\t',quote = F)
write.table(K562.h.lr[["MAPS"]],
            "ref/FitHiChIP/used/K562_H3K27ac_MAPS_hg38_resize_filter.lr",
            row.names = F,col.names = F,sep = '\t',quote = F)

### Pol2 ChIA-PET ------------------------------------------------------

K562.p <- readxl::read_xlsx("ref/FitHiChIP/K562/Reference_Loops/Table_K-R.xlsx",
                            sheet = 4, skip = 15)
K562.p <- addFilter(K562.p)
K562.p.gr <- bedpe2gi(K562.p,loc = c(1,2,3,1,4,5))
K562.p.lo <- liftOverbedpe(K562.p.gr,ch,verbose = F)
K562.p.lo <- resizebedpe(K562.p.lo)
K562.p.lo.flt <- K562.p.lo[K562.p.lo$filter == "pass"]

saveRDS(K562.p.lo,"ref/FitHiChIP/K562/Reference_Loops/K562_Pol2_ChIA-PET_hg38_resize_filter.rds")

## mESC --------------------------------------------------------------------

ch = import.chain("F:/Tanglab/database/liftOver/mm9ToMm10.over.chain")


### in situ Hi-C ------------------------------------------------------------
mESC_ref.df <- read.table("loop/mESC_postprocessed_pixels_10000_with_motifs.bedpe",header = T)
mESC_ref.df$width1 <- mESC_ref.df$x2 - mESC_ref.df$x1
mESC_ref.df$width2 <- mESC_ref.df$y2 - mESC_ref.df$y1

mESC_10k.ref <- bedpe2gi(mESC_ref.df)

### cohesin -----------------------------------------------------------------

mESC.c <- list()
mESC.c[["FitHiChIP"]] <- readxl::read_xlsx("ref/FitHiChIP/mES/Cohesin/Combined_Replicates_HIChIP_Loops/Table_MC-ALL.xlsx",
                                           sheet = 2, skip = 17)
mESC.c[["hichipper"]] <- readxl::read_xlsx("ref/FitHiChIP/mES/Cohesin/Combined_Replicates_HIChIP_Loops/Table_MC-ALL.xlsx",
                                           sheet = 6, skip = 13)

mESC.c[["FitHiChIP"]] <- addFilter(mESC.c[["FitHiChIP"]])
mESC.c[["hichipper"]] <- addFilter(mESC.c[["hichipper"]])

mESC.c.gr <- lapply(mESC.c, function(x){
  bedpe2gi(x,loc = c(1,2,3,1,4,5))
})

mESC.c.lo <- lapply(mESC.c.gr, function(gi){liftOverbedpe(gi,ch,verbose = F)})
mESC.c.lo <- lapply(mESC.c.lo, resizebedpe)
mESC.c.lo.flt <- lapply(mESC.c.lo, function(x) x[x$filter == "pass"])
saveRDS(mESC.c.lo,"ref/FitHiChIP/mES/Cohesin/Combined_Replicates_HIChIP_Loops/mESC_cohesin_mm10_resize_filter.rds")

#### save bedpe (for IGV) ----------------------------------------------------

mESC.c.df <- list()
mESC.c.df[["FitHiChIP"]] <- as.data.frame(mESC.c.lo[["FitHiChIP"]][mESC.c.lo[["FitHiChIP"]]$filter == "pass"])[,c(1:3,6:8,19)]
mESC.c.df[["FitHiChIP"]]$Q.Value <- -log10(mESC.c.df[["FitHiChIP"]]$Q.Value)
mESC.c.df[["hichipper"]] <- as.data.frame(mESC.c.lo[["hichipper"]][mESC.c.lo[["hichipper"]]$filter == "pass"])[,c(1:3,6:8,11)]
mESC.c.df[["hichipper"]]$PETCount <- log10(mESC.c.df[["hichipper"]]$PETCount)

write.table(mESC.c.df[["FitHiChIP"]],
            "ref/FitHiChIP/mES/Cohesin/Combined_Replicates_HIChIP_Loops/mESC_cohesin_FitHiChIP_mm10_resize_filter.bedpe",
            row.names = F,col.names = F,sep = '\t',quote = F)
write.table(mESC.c.df[["hichipper"]],
            "ref/FitHiChIP/mES/Cohesin/Combined_Replicates_HIChIP_Loops/mESC_cohesin_hichipper_mm10_resize_filter.bedpe",
            row.names = F,col.names = F,sep = '\t',quote = F)

#### save longrange (for WashU) ----------------------------------------------

mESC.c.lr <- list()
mESC.c.lr[["FitHiChIP"]] <- bedpe2lr(mESC.c.df[["FitHiChIP"]])
mESC.c.lr[["hichipper"]] <- bedpe2lr(mESC.c.df[["hichipper"]])
write.table(mESC.c.lr[["FitHiChIP"]],
            "ref/FitHiChIP/used/mESC_cohesin_FitHiChIP_mm10_resize_filter.lr",
            row.names = F,col.names = F,sep = '\t',quote = F)
write.table(mESC.c.lr[["hichipper"]],
            "ref/FitHiChIP/used/mESC_cohesin_hichipper_mm10_resize_filter.lr",
            row.names = F,col.names = F,sep = '\t',quote = F)


### H3K27ac -----------------------------------------------------------------

mESC.h <- list()
mESC.h[["FitHiChIP"]] <- readxl::read_xlsx("ref/FitHiChIP/mES/H3K27ac/Replicates/Table_MH-R_25M.xlsx",
                                           sheet = 2, skip = 17)
mESC.h[["hichipper"]] <- readxl::read_xlsx("ref/FitHiChIP/mES/H3K27ac/Replicates/Table_MH-R_25M.xlsx",
                                           sheet = 4, skip = 13)

mESC.h[["FitHiChIP"]] <- addFilter(mESC.h[["FitHiChIP"]])
colnames(mESC.h[["hichipper"]])[1:5] <- colnames(mESC.h[["FitHiChIP"]])[1:5]
mESC.h[["hichipper"]] <- addFilter(mESC.h[["hichipper"]])

mESC.h.gr <- lapply(mESC.h, function(x){
  bedpe2gi(x,loc = c(1,2,3,1,4,5))
})

mESC.h.lo <- lapply(mESC.h.gr, function(gi){liftOverbedpe(gi,ch,verbose = F)})
mESC.h.lo <- lapply(mESC.h.lo, resizebedpe)
mESC.h.lo.flt <- lapply(mESC.h.lo, function(x) x[x$filter == "pass"])
saveRDS(mESC.h.lo,"ref/FitHiChIP/mES/H3K27ac/Replicates/mESC_H3K27ac_mm10_resize_filter.rds")

#### save bedpe (for IGV) ----------------------------------------------------

mESC.h.df <- list()
mESC.h.df[["FitHiChIP"]] <- as.data.frame(mESC.h.lo[["FitHiChIP"]][mESC.h.lo[["FitHiChIP"]]$filter == "pass"])[,c(1:3,6:8,19)]
mESC.h.df[["FitHiChIP"]]$Q.Value <- -log10(mESC.h.df[["FitHiChIP"]]$Q.Value)
mESC.h.df[["hichipper"]] <- as.data.frame(mESC.h.lo[["hichipper"]][mESC.h.lo[["hichipper"]]$filter == "pass"])[,c(1:3,6:8,11)]
mESC.h.df[["hichipper"]]$PETcount <- log10(mESC.h.df[["hichipper"]]$PETcount)

write.table(mESC.h.df[["FitHiChIP"]],
            "ref/FitHiChIP/mES/H3K27ac/Replicates/mESC_H3K27ac_FitHiChIP_mm10_resize_filter.bedpe",
            row.names = F,col.names = F,sep = '\t',quote = F)
write.table(mESC.h.df[["hichipper"]],
            "ref/FitHiChIP/mES/H3K27ac/Replicates/mESC_H3K27ac_hichipper_mm10_resize_filter.bedpe",
            row.names = F,col.names = F,sep = '\t',quote = F)

#### save longrange (for WashU) ----------------------------------------------

mESC.h.lr <- list()
mESC.h.lr[["FitHiChIP"]] <- bedpe2lr(mESC.h.df[["FitHiChIP"]])
mESC.h.lr[["hichipper"]] <- bedpe2lr(mESC.h.df[["hichipper"]])
write.table(mESC.h.lr[["FitHiChIP"]],
            "ref/FitHiChIP/used/mESC_H3K27ac_FitHiChIP_mm10_resize_filter.lr",
            row.names = F,col.names = F,sep = '\t',quote = F)
write.table(mESC.h.lr[["hichipper"]],
            "ref/FitHiChIP/used/mESC_H3K27ac_hichipper_mm10_resize_filter.lr",
            row.names = F,col.names = F,sep = '\t',quote = F)


### SMC1 ChIA-PET ------------------------------------------------------

mESC.s <- readxl::read_xlsx("ref/FitHiChIP/mES/Reference_Loops/Table_M-R.xlsx",
                            sheet = 2, skip = 10)
mESC.s <- addFilter(mESC.s)
mESC.s.gr <- bedpe2gi(mESC.s,loc = c(1,2,3,1,4,5))
mESC.s.lo <- liftOverbedpe(mESC.s.gr,ch,verbose = F)
mESC.s.lo <- resizebedpe(mESC.s.lo)
mESC.s.lo.flt <- mESC.s.lo[mESC.s.lo$filter == "pass"]

saveRDS(mESC.s.lo,"ref/FitHiChIP/mES/Reference_Loops/mESC_SMC1_ChIA-PET_mm10_resize_filter.rds")

write.table(mESC.s.lo,"ref/FitHiChIP/mES/Reference_Loops/mESC_SMC1_ChIA-PET_mm10_resize_filter.txt",
            row.names = F,col.names = F,sep = '\t',quote = F)


# combined reference loops ------------------------------------------------
## GM12878 -----------------------------------------------------------------

GM.c.lo <- readRDS("ref/FitHiChIP/GM12878/Cohesin/Combined_Replicates_HIChIP_Loops/GM12878_cohesin_FitHiChIP_hg38_resize_filter.rds")
GM.h.lo <- readRDS("ref/FitHiChIP/GM12878/H3K27ac/Combined_Replicates_HIChIP_Loops/GM12878_H3K27ac_FitHiChIP_hg38_resize_filter.rds")
GM.c.lo.flt <- lapply(GM.c.lo, function(x) x[x$filter == "pass"])
GM.h.lo.flt <- lapply(GM.h.lo, function(x) x[x$filter == "pass"])

GM.cb.ref <- c(GM_10k.ref,GM.h.lo.flt$FitHiChIP,GM.h.lo.flt$MAPS,
               GM.c.lo.flt$FitHiChIP,GM.c.lo.flt$MAPS)
tdf <- data.frame(r1 = GM.cb.ref@anchor1,
                  r2 = GM.cb.ref@anchor2)
print(c(nrow(tdf),nrow(unique(tdf))))
GM.cb.ref <- unique(GM.cb.ref)
GM.cb.ref.f10k <- flankgi(GM.cb.ref,100001)

## K562 --------------------------------------------------------------------

K562.h.lo <- readRDS("ref/FitHiChIP/K562/H3K277ac/Combined_Replicates_HIChIP_Loops/K562_H3K27ac_hg38_resize_filter.rds")
K562.p.lo <- readRDS("ref/FitHiChIP/K562/Reference_Loops/K562_Pol2_ChIA-PET_hg38_resize_filter.rds")
K562.h.lo.flt <- lapply(K562.h.lo, function(x) x[x$filter == "pass"])
K562.h.lo.flt <- lapply(K562.h.lo, function(x) x[x$filter == "pass"])

K562.cb.ref <- c(K562_10k.ref,K562.h.lo.flt$FitHiChIP,
                 K562.h.lo.flt$MAPS,K562.p.lo.flt)
tdf <- data.frame(r1 = K562.cb.ref@anchor1,
                  r2 = K562.cb.ref@anchor2)
print(c(nrow(tdf),nrow(unique(tdf))))

K562.cb.ref <- unique(K562.cb.ref)
K562.cb.ref.f10k <- flankgi(K562.cb.ref,100001)

## mESC --------------------------------------------------------------------

mESC.c.lo <- readRDS("ref/FitHiChIP/mES/Cohesin/Combined_Replicates_HIChIP_Loops/mESC_cohesin_mm10_resize_filter.rds")
mESC.h.lo <- readRDS("ref/FitHiChIP/mES/H3K27ac/Replicates/mESC_H3K27ac_mm10_resize_filter.rds")
mESC.s.lo <- readRDS("ref/FitHiChIP/mES/Reference_Loops/mESC_SMC1_ChIA-PET_mm10_resize_filter.rds")

mESC.c.lo.flt <- lapply(mESC.c.lo, function(x) x[x$filter == "pass"])
mESC.h.lo.flt <- lapply(mESC.h.lo, function(x) x[x$filter == "pass"])
mESC.s.lo.flt <- mESC.s.lo[mESC.s.lo$filter == "pass"]

mESC.cb.ref <- c(mESC_10k.ref,mESC.h.lo.flt$FitHiChIP,
                 mESC.h.lo.flt$hichipper,
                 mESC.c.lo.flt$FitHiChIP,
                 mESC.c.lo.flt$hichipper,
                 mESC.s.lo.flt)

mESC.cb.ref <- c(mESC_10k.ref,mESC.h.lo.flt$FitHiChIP,
                 mESC.c.lo.flt$FitHiChIP,
                 mESC.s.lo.flt)

tdf <- data.frame(r1 = mESC.cb.ref@anchor1,
                  r2 = mESC.cb.ref@anchor2)
print(c(nrow(tdf),nrow(unique(tdf))))

mESC.cb.ref <- unique(mESC.cb.ref)
mESC.cb.ref.f10k <- flankgi(mESC.cb.ref,100001)


# SnapHiC loops for scNanoHi-C ------------------------------------------------------------

K562.sm.files <- read.table("loop/K562.postprocessed.summits.bedpe",header = T)
mESC.sm.files <- read.table("loop/mESC.postprocessed.summits.bedpe",header = T)
GM.sm.files <- read.table("loop/GM12878.postprocessed.summits.bedpe",header = T)
K562.sm.gi <- bedpe2gi(K562.sm.files)
GM.sm.gi <- bedpe2gi(GM.sm.files)
mESC.sm.gi <- bedpe2gi(mESC.sm.files)


# save for WashU ----------------------------------------------------------

GM.sm.df <- read.table("loop/GM12878.postprocessed.summits.bedpe",header = T)
GM.sm.df <- GM.sm.df[,c(1:6,10)]
GM.sm.df$fdr_dist <- -log10(GM.sm.df$fdr_dist)
GM.sm.lr <- bedpe2lr(GM.sm.df)
write.table(GM.sm.lr,
            "loop/GM12878.postprocessed.summits.lr",
            row.names = F,col.names = F,sep = '\t',quote = F)

mESC.sm.df <- mESC.sm.files[,c(1:6,10)]
mESC.sm.df$fdr_dist <- -log10(mESC.sm.df$fdr_dist)
mESC.sm.lr <- bedpe2lr(mESC.sm.df)
write.table(mESC.sm.lr,
            "loop/mESC.postprocessed.summits.lr",
            row.names = F,col.names = F,sep = '\t',quote = F)

# benchmarking ------------------------------------------------------------

GM.all.ov.cb <- fix.ov(GM.sm.gi, GM.cb.ref.f10k)
K562.all.ov.cb <- fix.ov(K562.sm.gi, K562.cb.ref.f10k)
mESC.all.ov.cb <- fix.ov(mESC.sm.gi, mESC.cb.ref.f10k)

# summary -----------------------------------------------------------------

ov.sum <- lapply(list(GM.all.ov.cb,
                      K562.all.ov.cb,
                      mESC.all.ov.cb), ov_stat) %>%
  Reduce(rbind,.) %>%
  `rownames<-`(c("GM12878","K562","mESC"))

write.csv(round(ov.sum,digits = 3),
          "loop/SnapHiC_combined_reference_10kloop_20kDis_benchmark_20230413.csv")

ov.sum <- read.csv("loop/SnapHiC_combined_reference_10kloop_20kDis_benchmark_20230413.csv",row.names = 1)
ov.sum$cell_line <- rownames(ov.sum)
ov.sum.m <- melt(ov.sum[,c(5:8)])

pdf("loop/SnapHiC_combined_reference_10kloop_20kDis_benchmark_20230413.pdf",
    width = 5.5,height = 4)
ggplot(ov.sum.m,aes(x=variable,y=value)) +
  geom_col(aes(fill = cell_line),
           width = 0.7,
           position = "dodge2") +
  theme_bw() + xlab("") + 
  labs(title = "Benchmarking of scNanoHi-C loops") +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,vjust = 0.5)) +
  scale_fill_manual(values = brewer.pal(10,"Set3")[c(4,7,10)])
dev.off()
