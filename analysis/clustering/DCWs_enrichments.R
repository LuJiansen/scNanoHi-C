library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(tidyverse)
library(pheatmap)

# diffrential windows -----------------------------------------------------

cpg.mg <- read.csv("cpg_score/all_used_raw_cpg_value_matrix.csv",row.names = 1)
auto.chr <- rownames(cpg.mg) %>% gsub("_.*","",.) %in% paste0("chr",1:22)
cpg.mg.used <- cpg.mg[auto.chr,rownames(meta.used)]

meta.df <- read.csv("scHiC_all_used_metadata.csv",row.names = 1)
meta.used <- meta.df[meta.df$filtered_trans_ratio < 30 &
                       meta.df$filtered_pass_count > 10000 &
                       meta.df$cell_line %in% c("GM","HG002","K562"),]

cell.used <- meta.used[meta.used$filtered_pass_count > 50000,]$rename
cpg.flt <- cpg_rank(cpg.mg.used[,cell.used])

df <- cpg.flt %>% as.data.frame() %>%
  mutate(windows = rownames(.)) %>% melt() %>%
  `colnames<-`(c("win","cell","value")) %>%
  mutate(cellline = gsub("_.*","",cell))

res <- df[df$cellline %in% c("GM","K562"),] %>% 
  group_by(win) %>% 
  do(w = wilcox.test(value~cellline, data=., paired=FALSE),
     t = t.test(value~cellline, data=., paired=FALSE)) %>% 
  summarise(win, w_pval = w$p.value,
            t_pval = t$p.value,
            mean1 = t$estimate[1],
            mean2 = t$estimate[2],
            log2FC = log2(t$estimate[1]/t$estimate[2]))
res <- res %>% mutate(w_padj = p.adjust(w_pval,method = "BH"),
                      t_padj = p.adjust(t_pval,method = "BH"),
                      mean_diff = mean1 - mean2)

res.sig <- res %>% filter(w_padj < 0.01 & abs(mean_diff) > 0.1)
res.sig$chr <- factor(gsub("_.*","",res.sig$win),
                      paste0("chr",c(1:22,"X")))
res.sig$cell_line <- "GM"
res.sig[res.sig$mean_diff < 0,]$cell_line <- "K562"
res.sig <- as.data.frame(res.sig)
rownames(res.sig) <- res.sig$win
res.sig <- res.sig %>% arrange(desc(log2FC))

pdf("comp_pca/DCWs_chromosome_distribution_barplot.pdf",width = 6,height = 4)
ggplot(res.sig,aes(x=chr)) +
  geom_bar(aes(fill = chr),width = 0.7) +
  scale_fill_manual(values = chr_col) +
  facet_wrap(vars(cell_line),scales = "free_y",ncol = 1) +
  theme_bw() +labs(title = "Distribution of differential compartmentalized windows") +
  scale_x_discrete("",breaks = levels(res.sig$chr),
                   label = c(1:22,"X")) +ylab("Number of windows") +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,vjust = 0.5))
dev.off()

library(pheatmap)
library(cowplot)
chr_col <- setNames(c(brewer.pal(9,"Set1"),
                      brewer.pal(12,"Set3"),
                      brewer.pal(8,"Accent")[c(8,6,5)]),
                    paste0("chr",c(1:22,"X","Y")))
meta.df <- read.csv("scHiC_all_used_metadata.csv",row.names = 1)
cell.ann = meta.df[,c("cell_line","group3","filtered_pass_count")]
cell.ann$log10counts <- log10(cell.ann$filtered_pass_count)
cell.ann <- cell.ann[colnames(cpg.flt),c("cell_line","group3","log10counts")]

sum(meta.flt[meta.flt$cell_line %in% c("GM","K562"),]$rename %in% colnames(cpg.flt))
hm.cell.used <- meta.flt[meta.flt$cell_line %in% c("GM","K562"),]$rename

p <- pheatmap(cpg.flt[res.sig$win,hm.cell.used],cluster_rows = T,
              annotation_col = cell.ann,show_rownames = F,
              annotation_row = res.sig[,c("chr","cell_line")],
              show_colnames = F,scale = "row",
              clustering_method = "ward.D2",
              main = "Differential scA/B windows (1Mb)",
              annotation_colors = list(cell_line = setNames(brewer.pal(3,"Set1")[1:3],
                                                            c("GM","HG002","K562")),
                                       chr = chr_col,group3 = group3.col),
              breaks = seq(-2,2,length.out=100),
              color = c(colorRampPalette(brewer.pal(11,"PRGn"))(100)))

save_plot("comp_pca/GM_vs_K562_Differential_2d_scAB_1M_windows_row_scale.pdf",
          p,base_width = 8,base_height = 6)

table(res.sig$cell_line)[1]
res.sig$type <- "up"
res.sig[res.sig$cell_line == "GM",]$type <- "down"

p <- pheatmap(cpg.flt[rev(res.sig$win),hm.cell.used],cluster_rows = F,
              annotation_col = cell.ann,show_rownames = F,
              annotation_row = res.sig[,c("chr","type")],
              show_colnames = F,#gaps_row = table(res.sig$cell_line)[1],
              clustering_method = "ward.D2",
              main = "Differential scA/B windows (1Mb)",
              annotation_colors = list(cell_line = setNames(brewer.pal(3,"Set1")[1:3],
                                                            c("GM","HG002","K562")),
                                       chr = chr_col,group3 = group3.col,
                                       type = setNames(brewer.pal(9,"Set3")[4:5],
                                                       c("up","down"))),
              color = c(colorRampPalette(brewer.pal(9,"PRGn"))(100)))
save_plot("comp_pca/GM_vs_K562_Differential_2d_scAB_1M_windows.pdf",
          p,base_width = 8,base_height = 6)


# enrichment analysis -----------------------------------------------------

library(LOLA)

GM.ann <- loadRegionDB('ref/ENCODE_GM12878_RegionDB',useCache = F)
K562.ann <- loadRegionDB('ref/ENCODE_K562_RegionDB',useCache = F)

flt.df <- flt.pca.ld.df[,c("chr","pos")] %>% unique()
flt.df$start <- pmax(as.numeric(flt.df$pos)-500000,0)
flt.df$end <- as.numeric(flt.df$pos)+500000

library(GenomicRanges)
flt.gr <- makeGRangesFromDataFrame(flt.df,seqnames.field = "chr",
                                   start.field = "start",
                                   end.field = "end",
                                   starts.in.df.are.0based = T)

win_enrich <- function(win.used,ann = K562.ann,
                       univ = flt.gr,title = NULL){
  res = runLOLA(userSets = univ[win.used],
                userUniverse = univ, ann, cores=4)
  res <- res[order(support, decreasing=TRUE),]
  
  p <- ggplot(res,aes(x=description,y=log2(oddsRatio))) +
    labs(title = title) +
    geom_col(aes(fill = pValueLog)) +
    theme_bw() + coord_flip() +
    scale_fill_viridis_c() +
    theme(plot.title = element_text(hjust = 0.5,vjust = 0.5))
  return(list(res = res,
              plot = p))
}

enrich_align <- function(pos1,pos2,neg1,neg2,
                         title1=NULL,title2=NULL){
  var.used <- c("description","oddsRatio","pValueLog","collection")
  df.pos <- rbind(as.data.frame(pos1)[,var.used],
                  as.data.frame(pos2)[,var.used])
  df.neg <- rbind(as.data.frame(neg1)[,var.used],
                  as.data.frame(neg2)[,var.used])
  
  p <- (ggplot(df.pos,aes(x=description,y=log2(oddsRatio))) +
          labs(title = title1) +
          facet_wrap(vars(collection))) /
    (ggplot(df.neg,aes(x=description,y=log2(oddsRatio))) +
       labs(title = title2) +
       facet_wrap(vars(collection))) &
    geom_col(aes(fill = pValueLog),width = 0.7) &
    geom_hline(yintercept = 0) &
    scale_y_continuous(expand = c(0.1,0.1)) &
    theme_bw() & coord_flip() &
    scale_fill_viridis_c() &
    plot_annotation(title = "") &
    theme(plot.title = element_text(hjust = 0.5,vjust = 0.5),
          panel.grid = element_blank())
  return(list(df.pos = df.pos,
              df.neg = df.neg,
              plot = p))
}
ftest <- function(df1,df2){
  mg <- full_join(df1[,c("support","c","description")],
                  df2[,c("support","c","description")],
                  by = "description") %>%
    dplyr::select(description,support.x,c.x,support.y,c.y) %>%
    column_to_rownames("description")
  
  test.res <- apply(mg, 1, function(x) fisher.test(matrix(x,nrow = 2)))
  
  mg$description <- rownames(mg)
  mg$pval <- lapply(test.res,function(x) x$p.value) %>% unlist()
  mg$oddsRatio <- lapply(test.res,function(x) x$estimate) %>% unlist()
  mg$LogPadj <- -log10(p.adjust(mg$pval,method = "BH"))
  
  return(mg)
}

dif.enrich <- list()
dif.enrich[["up"]][["K562"]] <- win_enrich(res.sig[res.sig$log2FC > 0,]$win,title = "GM12878 up-regulated (n=136)")
dif.enrich[["down"]][["K562"]] <- win_enrich(res.sig[res.sig$log2FC < 0,]$win,title = "GM12878 down-regulated (n=103)")
dif.enrich[["up"]][["GM12878"]] <- win_enrich(res.sig[res.sig$log2FC > 0,]$win,title = "GM12878 up-regulated (n=136)",ann = GM.ann)
dif.enrich[["down"]][["GM12878"]] <- win_enrich(res.sig[res.sig$log2FC < 0,]$win,title = "GM12878 down-regulated (n=103)",ann = GM.ann)

ann.used <- c("A compartment","B compartment","DNase-seq filtered",
              "H3K27ac","H3K27me3","H3K36me3","H3K4me1",
              "H3K4me3","H3K9me3","SuperEnhancer")

dif.mg.pos <- ftest(dif.enrich[["up"]][["GM12878"]]$res %>% filter(description %in% ann.used),
                    dif.enrich[["up"]][["K562"]]$res %>% filter(description %in% ann.used))
dif.mg.neg <- ftest(dif.enrich[["down"]][["GM12878"]]$res %>% filter(description %in% ann.used),
                    dif.enrich[["down"]][["K562"]]$res %>% filter(description %in% ann.used))

pdf("comp_pca/CpG_scAB_2d_flt_diff_windows_merged_enrichment_new.pdf",width = 10,height = 4)
ggplot(dif.mg.neg,aes(x=description,y=-log(pval))) +
  labs(title = paste("K562 up-regulated windows")) +
  ggplot(dif.mg.pos,aes(x=description,y=-log(pval))) +
  labs(title = paste("K562 down-regulated windows")) &
  geom_col(aes(fill = log2(oddsRatio)),width = 0.7) &
  geom_hline(yintercept = 2,linetype = "dashed") &
  theme_bw() & coord_flip() &
  scale_y_continuous("-Log10(p-value)",
                     expand = expansion(mult = c(0, .05))) &
  scale_fill_gradientn(colors = c(brewer.pal(11,"RdBu")),
                       oob = oob_squish_any,limits = c(-4,4)) &
  plot_annotation(title = "") &
  theme(plot.title = element_text(hjust = 0.5,vjust = 0.5),
        panel.grid = element_blank())
dev.off()

pdf("CpG_scAB_2d_flt_diff_windows_enrichment.pdf",width = 10,height = 4)
dif.enrich$up$K562$plot +
  dif.enrich$down$K562$plot &
  plot_annotation("K562 annotations") &
  theme(plot.title = element_text(hjust = 0.5,vjust = 0.5))

dif.enrich$up$GM12878$plot +
  dif.enrich$down$GM12878$plot &
  plot_annotation("GM annotations") &
  theme(plot.title = element_text(hjust = 0.5,vjust = 0.5))
dev.off()

pdf("comp_pca/CpG_scAB_2d_flt_diff_windows_enrichment_align_new.pdf",width = 10,height = 8)
enrich_align(dif.enrich$down$GM12878$res %>% filter(description %in% ann.used),
             dif.enrich$down$K562$res %>% filter(description %in% ann.used),
             dif.enrich$up$GM12878$res %>% filter(description %in% ann.used),
             dif.enrich$up$K562$res %>% filter(description %in% ann.used),
             "K562 up-regulated windows",
             "K562 down-regulated windows")
dev.off()


p <- pheatmap(win_expr[rev(rownames(res.sig)),"log2fc"],scale = "none",
              cluster_cols = F,cluster_rows = F,show_rownames = F,
              breaks = seq(-4,4,length.out = 100),
              color = colorRampPalette(c(brewer.pal(9,"RdBu")))(100))
save_plot("comp_pca/windows_log2fc_RNA_expr_of_DCWs_heatmap.pdf",
          p,base_width = 1,base_height = 6)
