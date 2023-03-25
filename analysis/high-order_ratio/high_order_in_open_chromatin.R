library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(ggrastr)

prefix = "GM12878"
prefix = "K562"
ratio.df <- read.csv(paste0("synergy/",prefix,"_high_order_ratio_new.csv"))

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

ratio.df$total_density <- get_density(ratio.df$total_cells,
                                      ratio.df$total_reads, n = 500)

summary(ratio.df$total_reads)
summary(ratio.df$total_cells)

ratio.df$filter <- "filtered"

ratio.flt2 <- ratio.df %>% 
  filter(total_reads > 0 & total_cells > 0) %>%
  filter(total_reads > quantile(total_reads,0.05) &
           total_reads < quantile(total_reads,0.95),
         total_cells > quantile(total_cells,0.05),
         total_cells < quantile(total_cells,0.95))
ratio.flt2$filter <- "pass"
ratio.df[ratio.df$ID %in% ratio.flt2$ID,]$filter <- "pass"
identical(ratio.flt2$ID,ratio.df[ratio.df$ID %in% ratio.flt2$ID,]$ID)
identical(ratio.flt2$comp,ratio.df[ratio.df$ID %in% ratio.flt2$ID,]$comp)

summary(ratio.flt2$total_reads)
summary(ratio.flt2$total_cells)

dim(ratio.df)
dim(ratio.flt2)

write.csv(ratio.df,paste0("synergy/",prefix,"_high_order_ratio_new.csv"))

library(ggrastr)
library(ggside)
pdf(paste0("synergy/",prefix,"_reads_cells_number_in_50k_windows_20230315.pdf"),
    width = 6, height = 4)
ggplot(ratio.df,aes(x=total_reads,y=total_cells)) +
  geom_point_rast(aes(color = total_density),size = 0.1) +
  geom_xsidedensity(aes(y = after_stat(density)), position = "stack") +
  geom_ysidedensity(aes(x = after_stat(density)), position = "stack") +
  geom_hline(yintercept = c(min(ratio.flt2$total_cells),
                            max(ratio.flt2$total_cells)),
             linetype = "dashed") +
  geom_vline(xintercept = c(min(ratio.flt2$total_reads),
                            max(ratio.flt2$total_reads)),
             linetype = "dashed") +
  scale_color_viridis_c() +
  theme_bw() +
  theme(panel.grid = element_blank())
dev.off()


ov.c.len <- read.table(paste0("synergy/",prefix,"_all_concatemer_length_by_50k_wins.txt"),
                       header = T)
ov.len2 <- read.table(paste0("synergy/",prefix,"_all_monomer_length_by_50k_wins2.txt"),header = T)

sum(ratio.flt2$ID %in% ov.c.len$query.id)
sum(ratio.flt2$ID %in% ov.len2$query.id)

colnames(ov.c.len) <- c("ID","mean_conc_len","median_conc_len")
colnames(ov.len2) <- c("ID","mean_mono_len","median_mono_len")

ratio.flt2 <- left_join(ratio.flt2,ov.c.len)
ratio.flt2 <- left_join(ratio.flt2,ov.len2)



# high-order ratio vs compartment score -----------------------------------------------------------------------
ratio.flt2$d_rr_comp <- get_density(ratio.flt2$read_ratio,ratio.flt2$comp, n = 500)
ratio.flt2$d_cr_comp <- get_density(ratio.flt2$cell_ratio,ratio.flt2$comp, n = 500)

# high-order ratio vs monomer length --------------------------------------
ratio.flt2$d_rr_mean_mono_len <- get_density(ratio.flt2$read_ratio,
                                             ratio.flt2$mean_mono_len, n = 500)
ratio.flt2$d_cr_mean_mono_len <- get_density(ratio.flt2$cell_ratio,
                                             ratio.flt2$mean_mono_len, n = 500)

ratio.flt2$d_rr_med_mono_len <- get_density(ratio.flt2$read_ratio,
                                             ratio.flt2$median_mono_len, n = 500)
ratio.flt2$d_cr_med_mono_len <- get_density(ratio.flt2$cell_ratio,
                                             ratio.flt2$median_mono_len, n = 500)

# high-order ratio vs concatemer length --------------------------------------
ratio.flt2$d_rr_mean_conc_len <- get_density(ratio.flt2$read_ratio,
                                             ratio.flt2$mean_conc_len, n = 500)
ratio.flt2$d_cr_mean_conc_len <- get_density(ratio.flt2$cell_ratio,
                                             ratio.flt2$mean_conc_len, n = 500)

ratio.flt2$d_rr_med_conc_len <- get_density(ratio.flt2$read_ratio,
                                            ratio.flt2$median_conc_len, n = 500)
ratio.flt2$d_cr_med_conc_len <- get_density(ratio.flt2$cell_ratio,
                                            ratio.flt2$median_conc_len, n = 500)

# high-order ratio vs coverage --------------------------------------
ratio.flt2$d_rr_rcov <- get_density(ratio.flt2$read_ratio,
                                    ratio.flt2$total_reads, n = 500)
ratio.flt2$d_cr_ccov <- get_density(ratio.flt2$cell_ratio,
                                    ratio.flt2$total_cells, n = 500)


# compartment score vs monomer length -------------------------------------------------------
ratio.flt2$d_comp_mean_mono_len <- get_density(ratio.flt2$comp,
                                               ratio.flt2$mean_mono_len, n = 500)
ratio.flt2$d_comp_median_mono_len <- get_density(ratio.flt2$comp,
                                                 ratio.flt2$median_mono_len, n = 500)

# compartment score vs concatemer length -------------------------------------------------------
ratio.flt2$d_comp_mean_conc_len <- get_density(ratio.flt2$comp,
                                               ratio.flt2$mean_conc_len, n = 500)
ratio.flt2$d_comp_median_conc_len <- get_density(ratio.flt2$comp,
                                                 ratio.flt2$median_conc_len, n = 500)

write.csv(ratio.flt2,paste0("synergy/",prefix,"_high_order_ratio_filter2.csv"),
          row.names = F)

# density scatterplot -----------------------------------------------------

ratio.flt2 <- read.csv(paste0("synergy/",prefix,"_high_order_ratio_filter2.csv"))

plot_cor <- function(df,x,y,d,xlab=NULL,ylab=NULL,
                     cor.method = "spearman",
                     sm.method = NULL,
                     zoom = NULL){
  p <- ggplot(df,aes_string(x=x,y=y)) +
    geom_point_rast(aes_string(color = d),size = 0.5) +
    xlab(xlab) + ylab(ylab) + 
    labs(title = paste0(paste0(cor.method,"'s r = "),
                        round(cor(df[,x],df[,y],
                                  method = cor.method),
                              digits = 2))) +
    geom_smooth(method = sm.method) + theme_bw() + 
    scale_color_gradientn("density",colors = rev(brewer.pal(11,"Spectral"))) &
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5,vjust = 0.5))
  if (!is.null(zoom)) {
    p <- p + coord_cartesian(xlim = quantile(df[,x],c(zoom,1-zoom)),
                             ylim = quantile(df[,y],c(zoom,1-zoom)))
  }
  return(p)
}


## high-order vs CS --------------------------------------------------------


pdf(paste0("synergy/",prefix,"_high_order_ratio_compartment_score_scatterplot_filter2.pdf"),
    width = 10,height = 4)
plot_cor(ratio.flt2,x = "read_ratio", y = "comp",
         d = "d_rr_comp",
         xlab = "Ratio of high-order concatemers",
         ylab = "Compartment Score") +
plot_cor(ratio.flt2,x = "cell_ratio", y = "comp",
         d = "d_cr_comp",
         xlab = "Ratio of high-order cells",
         ylab = "Compartment Score")
dev.off()

pdf(paste0("synergy/",prefix,"_high_order_ratio_compartment_score_scatterplot_filter2_zoom.pdf"),
    width = 10,height = 4)
plot_cor(ratio.flt2,x = "read_ratio", y = "comp",
         d = "d_rr_comp",
         xlab = "Ratio of high-order concatemers",
         ylab = "Compartment Score",
         zoom = 0.05) +
  plot_cor(ratio.flt2,x = "cell_ratio", y = "comp",
           d = "d_cr_comp",
           xlab = "Ratio of high-order cells",
           ylab = "Compartment Score",
           zoom = 0.05)
dev.off()


## high-order vs mono_len --------------------------------------------------

pdf(paste0("synergy/",prefix,"_high_order_ratio_mono_length_scatterplot_filter2.pdf"),
    width = 10,height = 8)
plot_cor(ratio.flt2,x = "mean_mono_len", y = "read_ratio",
         d = "d_rr_mean_mono_len",
         xlab = "Mean length of monomers",
         ylab = "Ratio of high-order concatemers") +
  plot_cor(ratio.flt2,x = "mean_mono_len", y = "cell_ratio",
           d = "d_cr_mean_mono_len",
           xlab = "Mean length of monomers",
           ylab = "Ratio of high-order cells") +
  plot_cor(ratio.flt2,x = "median_mono_len", y = "read_ratio",
           d = "d_rr_med_mono_len",
           xlab = "Median length of monomers",
           ylab = "Ratio of high-order concatemers") +
  plot_cor(ratio.flt2,x = "median_mono_len", y = "cell_ratio",
           d = "d_cr_med_mono_len",
           xlab = "Median length of monomers",
           ylab = "Ratio of high-order cells")
dev.off()


pdf(paste0("synergy/",prefix,"_high_order_ratio_mono_length_scatterplot_filter2_zoom.pdf"),
    width = 10,height = 8)
plot_cor(ratio.flt2,x = "mean_mono_len", y = "read_ratio",
         d = "d_rr_mean_mono_len",
         xlab = "Mean length of monomers",
         ylab = "Ratio of high-order concatemers",
         zoom = 0.05) +
  plot_cor(ratio.flt2,x = "mean_mono_len", y = "cell_ratio",
           d = "d_cr_mean_mono_len",
           xlab = "Mean length of monomers",
           ylab = "Ratio of high-order cells",
           zoom = 0.05) +
  plot_cor(ratio.flt2,x = "median_mono_len", y = "read_ratio",
           d = "d_rr_med_mono_len",
           xlab = "Median length of monomers",
           ylab = "Ratio of high-order concatemers",
           zoom = 0.05) +
  plot_cor(ratio.flt2,x = "median_mono_len", y = "cell_ratio",
           d = "d_cr_med_mono_len",
           xlab = "Median length of monomers",
           ylab = "Ratio of high-order cells",
           zoom = 0.05)
dev.off()

## high-order vs conc_len --------------------------------------------------


pdf(paste0("synergy/",prefix,"_high_order_ratio_concatemer_length_scatterplot_filter2.pdf"),
    width = 10,height = 8)
plot_cor(ratio.flt2,x = "mean_conc_len", y = "read_ratio",
         d = "d_rr_mean_conc_len",
         xlab = "Mean length of concatemers",
         ylab = "Ratio of high-order concatemers") +
  plot_cor(ratio.flt2,x = "mean_conc_len", y = "cell_ratio",
           d = "d_cr_mean_conc_len",
           xlab = "Mean length of concatemers",
           ylab = "Ratio of high-order cells") +
  plot_cor(ratio.flt2,x = "median_conc_len", y = "read_ratio",
           d = "d_rr_med_conc_len",
           xlab = "Median length of concatemers",
           ylab = "Ratio of high-order concatemers") +
  plot_cor(ratio.flt2,x = "median_conc_len", y = "cell_ratio",
           d = "d_cr_med_conc_len",
           xlab = "Median length of concatemers",
           ylab = "Ratio of high-order cells")
dev.off()

pdf(paste0("synergy/",prefix,"_high_order_ratio_concatemer_length_scatterplot_filter2_zoom.pdf"),
    width = 10,height = 8)
plot_cor(ratio.flt2,x = "mean_conc_len", y = "read_ratio",
         d = "d_rr_mean_conc_len",
         xlab = "Mean length of concatemers",
         ylab = "Ratio of high-order concatemers",
         zoom = 0.05) +
  plot_cor(ratio.flt2,x = "mean_conc_len", y = "cell_ratio",
           d = "d_cr_mean_conc_len",
           xlab = "Mean length of concatemers",
           ylab = "Ratio of high-order cells",
           zoom = 0.05) +
  plot_cor(ratio.flt2,x = "median_conc_len", y = "read_ratio",
           d = "d_rr_med_conc_len",
           xlab = "Median length of concatemers",
           ylab = "Ratio of high-order concatemers",
           zoom = 0.05) +
  plot_cor(ratio.flt2,x = "median_conc_len", y = "cell_ratio",
           d = "d_cr_med_conc_len",
           xlab = "Median length of concatemers",
           ylab = "Ratio of high-order cells",
           zoom = 0.05)
dev.off()


## high-order vs cov -------------------------------------------------------

pdf(paste0("synergy/",prefix,"_high_order_ratio_coverage_scatterplot_filter2.pdf"),
    width = 10,height = 4)
plot_cor(ratio.flt2,x = "total_reads", y = "read_ratio",
         d = "d_rr_rcov",
         xlab = "coverage of concatemers",
         ylab = "Ratio of high-order concatemers") +
  plot_cor(ratio.flt2,x = "total_cells", y = "cell_ratio",
           d = "d_cr_ccov",
           xlab = "coverage of cells",
           ylab = "Ratio of high-order cells")
dev.off()

## CS vs mono_len ----------------------------------------------------------

pdf(paste0("synergy/",prefix,"_comp_ratio_mono_length_scatterplot_filter2.pdf"),
    width = 10,height = 4)
plot_cor(ratio.flt2,x = "mean_mono_len", y = "comp",
         d = "d_comp_mean_mono_len",
         xlab = "Mean length of monomers",
         ylab = "Compartment scores") +
  plot_cor(ratio.flt2,x = "median_mono_len", y = "comp",
           d = "d_comp_median_mono_len",
           xlab = "Median length of monomers",
           ylab = "Compartment scores")
dev.off()

pdf(paste0("synergy/",prefix,"_comp_ratio_mono_length_scatterplot_filter2_zoom.pdf"),
    width = 10,height = 4)
plot_cor(ratio.flt2,x = "mean_mono_len", y = "comp",
         d = "d_comp_mean_mono_len",
         xlab = "Mean length of monomers",
         ylab = "Compartment scores",
         zoom = 0.05) +
  plot_cor(ratio.flt2,x = "median_mono_len", y = "comp",
           d = "d_comp_median_mono_len",
           xlab = "Median length of monomers",
           ylab = "Compartment scores",
           zoom = 0.05)
dev.off()

## CS vs conc_len ----------------------------------------------------------


pdf(paste0("synergy/",prefix,"_comp_ratio_conc_length_scatterplot_filter2.pdf"),
    width = 10,height = 4)
plot_cor(ratio.flt2,x = "mean_conc_len", y = "comp",
         d = "d_comp_mean_conc_len",
         xlab = "Mean length of concatemers",
         ylab = "Compartment scores") +
  plot_cor(ratio.flt2,x = "median_conc_len", y = "comp",
           d = "d_comp_median_conc_len",
           xlab = "Median length of concatemers",
           ylab = "Compartment scores")
dev.off()

pdf(paste0("synergy/",prefix,"_comp_ratio_conc_length_scatterplot_filter2_zoom.pdf"),
    width = 10,height = 4)
plot_cor(ratio.flt2,x = "mean_conc_len", y = "comp",
         d = "d_comp_mean_conc_len",
         xlab = "Mean length of concatemers",
         ylab = "Compartment scores",
         zoom = 0.05) +
  plot_cor(ratio.flt2,x = "median_conc_len", y = "comp",
           d = "d_comp_median_conc_len",
           xlab = "Median length of concatemers",
           ylab = "Compartment scores",
           zoom = 0.05)
dev.off()


# Partial correlation analysis --------------------------------------------


## CS vs monomer/concatemer length -----------------------------------------

library(ppcor)
library(pheatmap)
library(cowplot)

comp_len.pcor1 <- pcor(ratio.flt2[,c("comp","mean_mono_len",
                                    "mean_conc_len")],
                      method = "spearman")
comp_len.pcor2 <- pcor(ratio.flt2[,c("comp","median_mono_len",
                                    "median_conc_len")],
                      method = "spearman")

diag(comp_len.pcor1$estimate) <- 0
p1 <- pheatmap(comp_len.pcor1$estimate,cluster_rows = T,cluster_cols = T,
               breaks = seq(-max(abs(comp_len.pcor1$estimate)),
                            max(abs(comp_len.pcor1$estimate)),length.out = 50),
               color = rev(colorRampPalette(brewer.pal(11,"RdBu"))(50)),
               display_numbers = T,main = paste0(prefix," partial spearman's correlations"))

diag(comp_len.pcor2$estimate) <- 0
p2 <- pheatmap(comp_len.pcor2$estimate,cluster_rows = T,cluster_cols = T,
               breaks = seq(-max(abs(comp_len.pcor2$estimate)),
                            max(abs(comp_len.pcor2$estimate)),length.out = 50),
               color = rev(colorRampPalette(brewer.pal(11,"RdBu"))(50)),
               display_numbers = T,main = paste0(prefix," partial spearman's correlations"))

save_plot(paste0("synergy/",prefix,"_comp_mean_length_pcor_heatmap.pdf"),p1,
          base_height = 4.5,base_width = 5)
save_plot(paste0("synergy/",prefix,"_comp_median_length_pcor_heatmap.pdf"),p2,
          base_height = 4.5,base_width = 5)

comp_len.lm1 <- lm(mean_mono_len ~ mean_conc_len,data = ratio.flt2)
comp_len.lm2 <- lm(comp ~ mean_conc_len,data = ratio.flt2)

comp_len.lm3 <- lm(mean_conc_len ~ mean_mono_len,data = ratio.flt2)
comp_len.lm4 <- lm(comp ~ mean_mono_len,data = ratio.flt2)

cor(comp_len.lm1$residuals,comp_len.lm2$residuals,method = "pearson")
cor(comp_len.lm3$residuals,comp_len.lm4$residuals,method = "pearson")

comp_len.lm5 <- lm(median_mono_len ~ median_conc_len,data = ratio.flt2)
comp_len.lm6 <- lm(comp ~ median_conc_len,data = ratio.flt2)

comp_len.lm7 <- lm(median_conc_len ~ median_mono_len,data = ratio.flt2)
comp_len.lm8 <- lm(comp ~ median_mono_len,data = ratio.flt2)

cor(comp_len.lm5$residuals,comp_len.lm6$residuals,method = "pearson")
cor(comp_len.lm7$residuals,comp_len.lm8$residuals,method = "pearson")


pcor.df <- data.frame(resid_mean_mono_len = resid(comp_len.lm1) + 
                        mean(ratio.flt2$mean_mono_len),
                      resid_comp1 = resid(comp_len.lm2) + 
                        mean(ratio.flt2$comp),
                      resid_mean_conc_len = resid(comp_len.lm3) + 
                        mean(ratio.flt2$mean_conc_len),
                      resid_comp2 = resid(comp_len.lm4) + 
                        mean(ratio.flt2$comp),
                      resid_median_mono_len = resid(comp_len.lm5) + 
                        mean(ratio.flt2$median_mono_len),
                      resid_comp3 = resid(comp_len.lm6) + 
                        mean(ratio.flt2$comp),
                      resid_median_conc_len = resid(comp_len.lm7) + 
                        mean(ratio.flt2$median_conc_len),
                      resid_comp4 = resid(comp_len.lm8) + 
                        mean(ratio.flt2$comp))

pcor.df$density1 <- get_density(pcor.df$resid_mean_mono_len,
                                pcor.df$resid_comp1, n = 500)
pcor.df$density2 <- get_density(pcor.df$resid_mean_conc_len,
                                pcor.df$resid_comp2, n = 500)
pcor.df$density3 <- get_density(pcor.df$resid_median_mono_len,
                                pcor.df$resid_comp3, n = 500)
pcor.df$density4 <- get_density(pcor.df$resid_median_conc_len,
                                pcor.df$resid_comp4, n = 500)


pdf(paste0("synergy/",prefix,"_comp_length_pcor_scatterplot_filter2.pdf"),
    width = 10,height = 8)
plot_cor(pcor.df,x = "resid_mean_mono_len",
         y = "resid_comp1", d = "density1",
         xlab = "Mean length of monomers",
         ylab = "Compartment scores",
         cor.method = "pearson") +
  plot_cor(pcor.df,x = "resid_mean_conc_len", 
           y = "resid_comp2", d = "density2",
           xlab = "Mean length of concatemers",
           ylab = "Compartment scores",
           cor.method = "pearson") +
  plot_cor(pcor.df,x = "resid_median_mono_len", 
           y = "resid_comp3", d = "density3",
           xlab = "Median length of monomers",
           ylab = "Compartment scores",
           cor.method = "pearson") +
  plot_cor(pcor.df,x = "resid_median_conc_len", 
           y = "resid_comp4", d = "density4",
           xlab = "Median length of concatemers",
           ylab = "Compartment scores",
           cor.method = "pearson")
dev.off()

pdf(paste0("synergy/",prefix,"_comp_length_pcor_scatterplot_filter2_zoom.pdf"),
    width = 10,height = 8)
plot_cor(pcor.df,x = "resid_mean_mono_len",
         y = "resid_comp1", d = "density1",
         xlab = "Mean length of monomers",
         ylab = "Compartment scores",
         zoom = 0.05,
         cor.method = "pearson") +
  plot_cor(pcor.df,x = "resid_mean_conc_len", 
           y = "resid_comp2", d = "density2",
           xlab = "Mean length of concatemers",
           ylab = "Compartment scores",
           zoom = 0.05,
           cor.method = "pearson") +
  plot_cor(pcor.df,x = "resid_median_mono_len", 
           y = "resid_comp3", d = "density3",
           xlab = "Median length of monomers",
           ylab = "Compartment scores",
           zoom = 0.05,
           cor.method = "pearson") +
  plot_cor(pcor.df,x = "resid_median_conc_len", 
           y = "resid_comp4", d = "density4",
           xlab = "Median length of concatemers",
           ylab = "Compartment scores",
           zoom = 0.05,
           cor.method = "pearson")
dev.off()


# high-order ratio vs others ----------------------------------------------

ho.r.pcor1 <- pcor(ratio.flt2[,c("read_ratio","mean_mono_len",
                                 "mean_conc_len","total_reads","comp")],
                   method = "spearman")
ho.r.pcor2 <- pcor(ratio.flt2[,c("read_ratio","median_mono_len",
                                 "median_conc_len","total_reads","comp")],
                   method = "spearman")

ho.c.pcor1 <- pcor(ratio.flt2[,c("cell_ratio","mean_mono_len",
                                 "mean_conc_len","total_cells","comp")],
                   method = "spearman")
ho.c.pcor2 <- pcor(ratio.flt2[,c("cell_ratio","median_mono_len",
                                 "median_conc_len","total_cells","comp")],
                   method = "spearman")


diag(ho.r.pcor1$estimate) <- 0
p1 <- pheatmap(ho.r.pcor1$estimate,cluster_rows = T,cluster_cols = T,
               breaks = seq(-max(abs(ho.r.pcor1$estimate)),
                            max(abs(ho.r.pcor1$estimate)),length.out = 50),
               color = rev(colorRampPalette(brewer.pal(11,"RdBu"))(50)),
               display_numbers = T,main = paste0(prefix," partial spearman's correlations"))

diag(ho.r.pcor2$estimate) <- 0
p2 <- pheatmap(ho.r.pcor2$estimate,cluster_rows = T,cluster_cols = T,
               breaks = seq(-max(abs(ho.r.pcor2$estimate)),
                            max(abs(ho.r.pcor2$estimate)),length.out = 50),
               color = rev(colorRampPalette(brewer.pal(11,"RdBu"))(50)),
               display_numbers = T,main = paste0(prefix," partial spearman's correlations"))

diag(ho.c.pcor1$estimate) <- 0
p3 <- pheatmap(ho.c.pcor1$estimate,cluster_rows = T,cluster_cols = T,
               breaks = seq(-max(abs(ho.c.pcor1$estimate)),
                            max(abs(ho.c.pcor1$estimate)),length.out = 50),
               color = rev(colorRampPalette(brewer.pal(11,"RdBu"))(50)),
               display_numbers = T,main = paste0(prefix," partial spearman's correlations"))

diag(ho.c.pcor2$estimate) <- 0
p4 <- pheatmap(ho.c.pcor2$estimate,cluster_rows = T,cluster_cols = T,
               breaks = seq(-max(abs(ho.c.pcor2$estimate)),
                            max(abs(ho.c.pcor2$estimate)),length.out = 50),
               color = rev(colorRampPalette(brewer.pal(11,"RdBu"))(50)),
               display_numbers = T,main = paste0(prefix," partial spearman's correlations"))


save_plot(paste0("synergy/",prefix,"_read_highorder_vs_others1_pcor_heatmap.pdf"),p1,
          base_height = 4.5,base_width = 5)
save_plot(paste0("synergy/",prefix,"_read_highorder_vs_others2_pcor_heatmap.pdf"),p2,
          base_height = 4.5,base_width = 5)
save_plot(paste0("synergy/",prefix,"_cell_highorder_vs_others1_pcor_heatmap.pdf"),p3,
          base_height = 4.5,base_width = 5)
save_plot(paste0("synergy/",prefix,"_cell_highorder_vs_others2_pcor_heatmap.pdf"),p4,
          base_height = 4.5,base_width = 5)



pcor_plot <- function(df,x,y,cov,xlab=NULL,ylab=NULL,zoom=0.05){
  df <- df[,c(x,y,cov)]
  colnames(df) <- c("x","y",paste0("cov",1:length(cov)))
  lm1 <-  lm(formula(paste("x ~ ",paste0("cov",1:length(cov),collapse = "+"))),data = df)
  lm2 <-  lm(formula(paste("y ~ ",paste0("cov",1:length(cov),collapse = "+"))),data = df)
  
  print(cor(lm1$residuals,lm2$residuals,method = "pearson"))
  
  pcor_df <- data.frame(resid_x = resid(lm1) + mean(df$x),
                        resid_y = resid(lm2) + mean(df$y))
  pcor_df$density <- get_density(pcor_df$resid_x,
                                 pcor_df$resid_y, n = 500)
  p1 <- plot_cor(pcor_df,x = "resid_x",
                 y = "resid_y", d = "density",
                 xlab = xlab,
                 ylab = ylab,
                 cor.method = "pearson")
  p2 <- plot_cor(pcor_df,x = "resid_x",
                 y = "resid_y", d = "density",
                 xlab = xlab,
                 ylab = ylab,
                 cor.method = "pearson",
                 zoom = zoom)
  p <- p1 + p2
  return(list(df = pcor_df,
              p1 = p))
}

pcor.p <- list()
pcor.p[["rr_mono"]] <- pcor_plot(ratio.flt2,"mean_mono_len","read_ratio",
                                 c("mean_conc_len","total_reads","comp"),
                                 "Mean length of monomers",
                                 "Ratio of high-order concatemers")

pcor.p[["rr_comp"]] <- pcor_plot(ratio.flt2,"comp","read_ratio",
                                 c("mean_conc_len","total_reads","mean_mono_len"),
                                 "Compartment scores",
                                 "Ratio of high-order concatemers")

pcor.p[["cr_mono"]] <- pcor_plot(ratio.flt2,"mean_mono_len","cell_ratio",
                                 c("mean_conc_len","total_reads","comp"),
                                 "Mean length of monomers",
                                 "Ratio of high-order cells")

pcor.p[["cr_comp"]] <- pcor_plot(ratio.flt2,"comp","cell_ratio",
                                 c("mean_conc_len","total_reads","mean_mono_len"),
                                 "Compartment scores",
                                 "Ratio of high-order cells")

save_plot(paste0("synergy/",prefix,"_read_highorder_vs_mono_pcor_scatterplot.pdf"),
          pcor.p[["rr_mono"]]$p,base_height = 4,base_width = 10)
save_plot(paste0("synergy/",prefix,"_read_highorder_vs_comp_pcor_scatterplot.pdf"),
          pcor.p[["rr_comp"]]$p,base_height = 4,base_width = 10)
save_plot(paste0("synergy/",prefix,"_cell_highorder_vs_mono_pcor_scatterplot.pdf"),
          pcor.p[["cr_mono"]]$p,base_height = 4,base_width = 10)
save_plot(paste0("synergy/",prefix,"_cell_highorder_vs_comp_pcor_scatterplot.pdf"),
          pcor.p[["cr_comp"]]$p,base_height = 4,base_width = 10)





