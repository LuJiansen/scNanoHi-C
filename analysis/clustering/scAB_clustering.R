library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(tidyverse)
library(pheatmap)

cpg_rank <- function(cpg){
  cpg <- cpg[apply(cpg, 1, function(x) sum(is.na(x)) == 0),]
  cpg.rank <- apply(cpg, 2, function(x) (rank(x)-1)/c(nrow(cpg)-1))
  return(cpg.rank)
}

meta.df <- read.csv("scHiC_all_used_metadata.csv",row.names = 1)
meta.used <- meta.df[meta.df$filtered_trans_ratio < 30 &
                       meta.df$filtered_pass_count > 10000 &
                       meta.df$cell_line %in% c("GM","HG002","K562"),]

cpg.mg <- read.csv("cpg_score/all_used_raw_cpg_value_matrix.csv",row.names = 1)
auto.chr <- rownames(cpg.mg) %>% gsub("_.*","",.) %in% paste0("chr",1:22)
cpg.mg.used <- cpg.mg[auto.chr,rownames(meta.used)]

### filtered PCA ------------------------------------------------------------

cell.used <- meta.used[meta.used$filtered_pass_count > 50000,]$rename
cpg.flt <- cpg_rank(cpg.mg.used[,cell.used])

cpg.flt.pca <- prcomp(t(cpg.flt)) 

plot_pca_var <- function(pca,n){
  sum <- summary(pca)
  sum$importance[2,1:n] %>% 
    as.data.frame() %>% `colnames<-`("Var") %>%
    mutate(name = factor(rownames(.),levels = rownames(.))) %>%
    ggplot(.,aes(x=name,y=Var)) +
    geom_point() +
    geom_line(group=1) +
    theme_bw() +
    xlab("") +ylab("Explained variance %") +
    labs(title = "percentage of explained variance of PC") +
    theme(panel.grid = element_blank()) -> p
  return(list(sum = sum,
              plot = p))
}

cpg.flt.pca.var <- plot_pca_var(cpg.flt.pca,10)
pdf("comp_pca/CpG_2d_PCA_Var_rank_flt_20221209.pdf",width = 4,height = 4)
cpg.flt.pca.var$plot
dev.off()
flt.pca.label <- paste0(colnames(cpg.flt.pca.var$sum$importance)," (",
                        round(cpg.flt.pca.var$sum$importance[2,]*100,
                              digits = 1),"%)")

cpg.flt.pca.df <- cpg.flt.pca$x[,1:10] %>% 
  as.data.frame() %>%
  `colnames<-`(paste("rank_flt",colnames(.),sep = "_")) %>%
  mutate(cell = rownames(.))

sum(!(cpg.flt.pca.df$cell %in% meta.flt$rename))

cpg.flt.pca.df <- left_join(cpg.flt.pca.df,meta.df,
                            by = c("cell" = "rename"))

library(patchwork)
pdf("comp_pca/CpG_2d_PCA_rank_normalize_filtered_20221209.pdf",width = 9,height = 4.3)
ggplot(cpg.flt.pca.df,aes(x=rank_flt_PC1,y=rank_flt_PC2)) +
  geom_point(aes(color = cell_line)) +
  theme_bw() +
  xlab(flt.pca.label[1]) +
  ylab(flt.pca.label[2]) +
  labs(title = "rank NA removed PC1-PC2") +
  scale_color_manual(values = brewer.pal(9,"Set1")) +
  theme(legend.position = "none") +
  ggplot(cpg.flt.pca.df,aes(x=rank_flt_PC3,y=rank_flt_PC2)) +
  geom_point(aes(color = cell_line)) +
  theme_bw() + 
  ylab(flt.pca.label[2]) +
  xlab(flt.pca.label[3]) +
  labs(title = "rank NA removed PC2-PC3") +
  scale_color_manual(values = brewer.pal(9,"Set1")) &
  plot_annotation(title = 'PCA of Normalized nearby CpG') &
  theme(plot.title = element_text(hjust = 0.5,vjust = 0.5),
        panel.grid = element_blank())
dev.off()

pdf("comp_pca/CpG_2d_PCA_rank_normalize_by_group_filtered_20221209.pdf",width = 9,height = 4.3)
ggplot(cpg.flt.pca.df,aes(x=rank_flt_PC1,y=rank_flt_PC2)) +
  geom_point(aes(color = group3)) +
  theme_bw() +
  xlab(flt.pca.label[1]) +
  ylab(flt.pca.label[2]) +
  labs(title = "rank NA removed PC1-PC2") +
  scale_color_manual(values = group3.col) +
  theme(legend.position = "none") +
  ggplot(cpg.flt.pca.df,aes(x=rank_flt_PC3,y=rank_flt_PC2)) +
  geom_point(aes(color = group3)) +
  theme_bw() + 
  ylab(flt.pca.label[2]) +
  xlab(flt.pca.label[3]) +
  labs(title = "rank NA removed PC2-PC3") +
  scale_color_manual(values = group3.col) &
  plot_annotation(title = 'PCA of Normalized nearby CpG') &
  theme(plot.title = element_text(hjust = 0.5,vjust = 0.5),
        panel.grid = element_blank())
dev.off()


# UMAP under different cut-off --------------------------------------------

cpg.cov <- apply(cpg.mg.used, 2, function(x) sum(!is.na(x))) %>%
  as.data.frame() %>%
  `colnames<-`("nbins")
cpg.cov$ratio <- cpg.cov$nbins/nrow(cpg.mg.used)*100

summary(cpg.cov$nbins)
summary(cpg.cov$ratio)

set.seed(42)
run_cpg_umap <- function(df = cpg.mg.used,cutoff){
  cell.used <- meta.used[meta.used$filtered_pass_count > cutoff,]$rename
  df.used <- cpg_rank(df[,cell.used])
  print(dim(df.used))
  
  cpg.flt.umap <- uwot::umap(t(df.used),metric = "cosine",
                             min_dist = 0.3,verbose = F,
                             pca = 5,n_neighbors = 10)
  cpg.flt.umap.df <- as.data.frame(cpg.flt.umap) %>%
    `colnames<-`(c("UMAP_1","UMAP_2")) %>%
    mutate(cell = rownames(.))
  
  cpg.flt.umap.df <- left_join(cpg.flt.umap.df,meta.used,
                               by = c("cell" = "rename"))
  return(cpg.flt.umap.df)
}


scAB.umap <- lapply(c(10000,20000,30000,40000,50000),
                    function(x){
                      run_cpg_umap(cutoff = x)
                    })
names(scAB.umap) <- c(10000,20000,30000,40000,50000)

nstats <- lapply(c(10000,20000,30000,40000,50000),function(x){
  cell.used <- meta.used[meta.used$filtered_pass_count > x,]$rename
  df.used <- cpg_rank(cpg.mg.used[,cell.used])
  print(dim(df.used))})
names(nstats) <- c(10000,20000,30000,40000,50000)

lapply(nstats, function(x) x[1]) %>% unlist()
lapply(nstats, function(x) round(x[1]/2930*100,digits = 1)) %>% unlist()

lapply(nstats, function(x) x[2]) %>% unlist()
lapply(nstats, function(x) round(x[2]/1057*100,digits = 1)) %>% unlist()

scAB.umap.p <- lapply(names(scAB.umap), function(x){
  ggplot(scAB.umap[[x]],aes(x=UMAP_1,y=UMAP_2)) +
    geom_point(aes(color = group3),size = 1) + 
    scale_color_manual(values = group3.col) +
    labs(title = paste0("cutoff = ",x,
                        "; C = ",nstats[[x]][2],
                        " (",round(nstats[[x]][2]/1057*100,digits = 1),"%)",
                        "; W = ",nstats[[x]][1],
                        " (",round(nstats[[x]][1]/2930*100,digits = 1),"%)")) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5,vjust = 0.5))
})
names(scAB.umap.p) <- c(10000,20000,30000,40000,50000)

pdf("scAB_cluster_under_cutoff.pdf",width = 10,height = 8)
scAB.umap.p$`10000` + scAB.umap.p$`20000` +
  scAB.umap.p$`30000` + scAB.umap.p$`40000`
dev.off()