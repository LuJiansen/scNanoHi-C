library(rtracklayer)
library(MASS)
library(Matrix)
library(dplyr)
library(pbmcapply)
library(tidyverse)

prefix = "GM12878"

# chromHMM ----------------------------------------------------------------

GM.chrHMM <- read.table("ref/GM_12878_chromHMM_18_states_ENCFF338RIC.bed")
GM.chrHMM.list <- split(GM.chrHMM,GM.chrHMM$V4)
names(GM.chrHMM.list)[18] <- "ZNF_Rpts"

for (i in names(GM.chrHMM.list)) {
  write.table(GM.chrHMM.list[[i]],
              paste0("ref/GM12878_chromHMM_",i,"_ENCFF338RIC.bed"),
              row.names = F,col.names = F,sep = '\t',quote = F)
}

K562.chrHMM <- read.table("ref/ENCODE_K562_epigenome/K562_chromHMM_18state_ENCFF963KIA.bed")
K562.chrHMM.list <- split(K562.chrHMM,K562.chrHMM$V4)
names(K562.chrHMM.list)[18] <- "ZNF_Rpts"

for (i in names(K562.chrHMM.list)) {
  write.table(K562.chrHMM.list[[i]],paste0("ref/K562_chromHMM_",i,"_ENCFF963KIA.bed"),
              row.names = F,col.names = F,sep = '\t',quote = F)
}


# load annotation database ------------------------------------------------

library(LOLA)

GM.ann <- loadRegionDB('ref/ENCODE_GM12878_RegionDB')
K562.ann <- loadRegionDB('ref/ENCODE_K562_RegionDB')

chrHMM.col <- setNames(list(c(255,0,0),c(255,69,0),c(255,69,0),c(255,69,0),
                            c(0,128,0),c(0,100,0),c(194,225,5),c(194,225,5),
                            c(255,195,77),c(255,195,77),c(255,255,0),c(102,205,170),
                            c(138,145,208),c(205,92,92),c(189,183,107),c(128,128,128),
                            c(192,192,192),c(255,255,255)),
                       c("TssA","TssFlnk","TssFlnkU","TssFlnkD",
                         "Tx","TxWk","EnhG1","EnhG2","EnhA1",
                         "EnhA2","EnhWk","ZNF_Rpts","Het",
                         "TssBiv","EnhBiv","ReprPC","ReprPCWk",
                         "Quies")) %>%
  lapply(., function(x) rgb(x[1],x[2],x[3],maxColorValue = 255))


epiAnno <- setNames(c("#5AAE61","#9970AB","#8DD3C7","#FDB462","#FFD92F","#FB8072","#B3DE69","#80B1D3","#BEBADA",
                      brewer.pal(9,"Spectral")[4:1],"#E41A1C"),
                    c("A compartment","B compartment","DNase-seq",
                      "H3K27ac","H3K4me1","H3K4me3","H3K36me3","H3K27me3","H3K9me3",
                      "Unexpressed genes","Low expressed genes","High expressed genes",
                      "Very high expressed genes","SuperEnhancer"))

TFAnno <- setNames(brewer.pal(8,"Set2")[c(6,rep(2,10),rep(1,3),rep(7,3))],
                   c("EP300","BATF","EBF1","ETS1",
                     "IKZF1","MAX","MYC","PAX5","RAD21","RUNX3","SPI1",
                     "RNAPolIIA","RNAPolIIA-pS2","RNAPolIIA-pS5",
                     "CTCF","SMC3","YY1"))

lolatest <- function(lola,mc.cores = 5){
  lola <- as.data.frame(lola)
  rownames(lola) <- lola$description
  lola.split <- split(lola,lola$description)
  test.res <- pbmclapply(lola.split,function(x){
    fisher.test(matrix(as.integer(x[,c("support","b","c","d")]),nrow = 2))
  },mc.cores = mc.cores)
  
  res <- lapply(test.res, function(x){
    c(x$p.value,x$estimate,x$conf.int[1],x$conf.int[2])
  }) %>% Reduce(rbind,.) %>%
    as.data.frame() %>%
    `colnames<-`(c("pval","oddsRatio","confIntLow","confIntHigh"))
  res <- cbind(lola[names(lola.split),
                    c("collection","description","support",
                      "b","c","d","filename","size")],res)
  res$LogPadj <- -log10(p.adjust(res$pval,method = "BH"))
  return(res)
}

# merged monomer enrichment --------------------------------------------------------------------

mono.ann <- loadRegionDB(paste0('ref/ENCODE_',prefix,'_RegionDB'))
mono.gr <- readRDS(paste0(prefix,"_all_gr.rds"))

mono.lola <- runLOLA(userSets = mono.gr[mono.gr$multiway == TRUE,],
                     userUniverse = mono.gr,
                     mono.ann, cores=30)
mono.lola <- mono.lola[order(support, decreasing=TRUE),]
saveRDS(mono.lola,paste0(prefix,'_all_highorder_LOLA_enrich_20230309.rds'))

mono.lola2 <- lolatest(mono.lola,30)
saveRDS(mono.lola2,paste0(prefix,'_mono_all_highorder_fisher_enrich_20230309.rds'))

mono.lola2$description <- factor(mono.lola2$description,
                                 levels = c(names(epiAnno),names(TFAnno),names(chrHMM.col)))

pdf(paste0("synergy/",prefix,"_all_high_order_monomer_fisher_test_new.pdf"),
    width = 8,height = 5)
ggplot(mono.lola2[names(epiAnno)[-3],],
       aes(x=description,y=log2(oddsRatio))) +
  geom_col(aes(fill = description),
           width = 0.7,color = "black") +
  geom_errorbar(aes(ymin = log2(confIntLow),
                    ymax = log2(confIntHigh)),
                width = 0.3) +
  geom_hline(yintercept = 0) + ylab("") +
  theme_bw() + labs(title = paste0(prefix," Enrichment of high-order monomers")) +
  scale_fill_manual(values = epiAnno[-3]) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        panel.grid = element_blank())
dev.off()

pdf(paste0("synergy/",prefix,"_all_high_order_monomer_chromHMM_fisher_test_new.pdf"),
    width = 8,height = 5)
ggplot(mono.lola2[names(chrHMM.col),],
       aes(x=description,y=log2(oddsRatio))) +
  geom_col(aes(fill = description),
           width = 0.7,color = "black") +
  geom_errorbar(aes(ymin = log2(confIntLow),
                    ymax = log2(confIntHigh)),
                width = 0.3) +
  geom_hline(yintercept = 0) +
  theme_bw() + labs(title = paste0(prefix," enrichment of high-order monomers in 18-state chromHMM")) +
  scale_fill_manual(values = chrHMM.col) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        panel.grid = element_blank())
dev.off()


# single-cell enrichment --------------------------------------------------

mono.grl <- split(mono.gr,mono.gr$cell)

mono.cell.res <- pbmclapply(mono.grl, function(gr){
  locResults <- runLOLA(userSets = gr[gr$multiway == TRUE,],
                        userUniverse = gr,
                        mono.ann, cores=10)
  locResults <- locResults[order(support, decreasing=TRUE),]
  lola2 <- lolatest(locResults,10)
  return(lola2)
}, mc.cores = 5)

for (i in names(mono.cell.res)) {
  mono.cell.res[[i]]$cell <- i
}
saveRDS(mono.cell.res,paste0(prefix,"_mono_all_cell_highorder_ftest_enrich_20230309.rds"))

mono.cell.res <- readRDS(paste0("synergy/",prefix,"_mono_all_cell_highorder_ftest_enrich_20230309.rds"))
mono.cell.res.df <- Reduce(rbind,mono.cell.res)
mono.cell.res.df.m <- mono.cell.res.df %>% 
  group_by(description) %>%
  summarise(mean_logpval = mean(-log10(pval)))
mono.cell.res.df <- left_join(mono.cell.res.df,mono.cell.res.df.m)


mono.cell.res.df$description <- factor(mono.cell.res.df$description,
                                       levels = c(names(epiAnno),names(TFAnno),names(chrHMM.col)))

sum(!(mono.cell.res.df$cell %in% meta$rename2))

mono.cell.res.df <- left_join(mono.cell.res.df,
                              meta[,c("rename2","filtered_pass_count",
                                      "filtered_trans_ratio",
                                      "filtered_high_order")],
                              by = c("cell" = "rename2"))

pdf(paste0("synergy/",prefix,"_filter_singlecell_high_order_monomer_fisher_test_new.pdf"),
    width = 8,height = 5)
ggplot(mono.cell.res.df[mono.cell.res.df$description %in% names(epiAnno)[-3] &
                          mono.cell.res.df$filtered_pass_count > 10000 &
                          mono.cell.res.df$filtered_high_order > 30,],
       aes(x=description,y=log2(oddsRatio))) +
  geom_boxplot(aes(fill = description),
               width = 0.7,color = "black",
               outlier.shape = NA) + ylim(-1,1) +
  geom_hline(yintercept = 0) + xlab("") +
  theme_bw() + labs(title = "",prefix," Enrichment of high-order monomers in single-cells") +
  scale_fill_manual(values = epiAnno[-3]) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        panel.grid = element_blank())
dev.off()

pdf(paste0("synergy/",prefix,"_filter_singlecell_high_order_monomer_chromHMM_fisher_test_new.pdf"),
    width = 8,height = 5)
ggplot(mono.cell.res.df[mono.cell.res.df$description %in% names(chrHMM.col) &
                          mono.cell.res.df$filtered_pass_count > 10000 &
                          mono.cell.res.df$filtered_high_order > 30,],
       aes(x=description,y=log2(oddsRatio))) +
  geom_boxplot(aes(fill = description),
               width = 0.7,color = "black",
               outlier.shape = NA) +
  geom_hline(yintercept = 0) +
  theme_bw() + labs(title = paste0(prefix," enrichment of high-order monomers in 18-state chromHMM")) +
  scale_fill_manual(values = chrHMM.col) + 
  scale_y_continuous(limits = c(-1,1),oob = oob_squish_any) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        panel.grid = element_blank())
dev.off()
