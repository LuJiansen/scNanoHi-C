source('synergy_funcs.R')

prefix="all"
mc.cores.used = 10
hg38.ref <- readRDS('hg38_seqinfo.rds')
bins <- readRDS('scCh_20k_bins.rds')
gc_frag_cov <- readRDS('gc_frag_cov.rds')


all.df <- read.table(paste0("GM_",prefix,"_merged_monomer.txt.gz"))
colnames(all.df) <- c("chr", "start", "end", "read_name", "read_idx", "cell")

all.gr <- makeGRangesFromDataFrame(all.df, keep.extra.columns = T,
                                    starts.in.df.are.0based = T)
seqlengths(all.gr) <- seqlengths(hg38.ref)[names(seqlengths(all.gr))]
saveRDS(all.gr, paste0("mono_", prefix, "_gr.rds"))


## get chromunity ----------------------------------------------------------

all.cc <- get_chromunity(all.gr, bin_obj = bins, mc.cores = mc.cores.used)
saveRDS(all.cc, paste0("mono_", prefix, "_cc.rds"))

all.cc$read_idx <- all.cc$cid
all.cc$cell_idx <- as.numeric(factor(all.cc$cell))


## get candidate binsets ---------------------------------------------------

all_bs_ov <- find_ov_binsets(all.cc, unique(all.cc$chid),
                            bins = bins$bins, mc.cores = mc.cores.used)
saveRDS(all_bs_ov, paste0("mono_",prefix,"_bs_ov.rds"))

unique(all_bs_ov$binsets$chid) %>% length()

all_bs_ov.ann <- sc_annotate(all_bs_ov$binsets, concatemers = all.cc,
                           covariates = gc_frag_cov,
                           resolution = 2e4, mc.cores = mc.cores.used)

saveRDS(all_bs_ov.ann, paste0("mono_",prefix,"_bs_ov_cc_scann.rds"))

sum(all_bs_ov.ann$ncells > all_bs_ov.ann$nreads)


## generate background binsets ---------------------------------------------

all_back_gr <- make_background(all_bs_ov$binsets)
all_back_gr.ann <- sc_annotate(all_back_gr,concatemers = all.cc,
                                covariates = gc_frag_cov,
                                resolution = 2e4,mc.cores = mc.cores.used) 
saveRDS(all_back_gr, paste0("mono_",prefix,"_background.rds"))
saveRDS(all_back_gr.ann, paste0("mono_",prefix,"_background_ann.rds"))

## train model -------------------------------------------------------------

all_reads_model = fit(all_back_gr.ann, "nreads")
all_cells_model = fit(all_back_gr.ann, "ncells")

all_bs_ov.score <- sscore(all_bs_ov.ann, all_reads_model, all_cells_model) 


## synergy test ------------------------------------------------------------

all_synergy = synergy(all_bs_ov.score, all_reads_model, all_cells_model)

plot(density(all_synergy$nreads$fdr - all_synergy$ncells$fdr))

all.support.sum <- all_synergy$binsets[cardinality > 2,
                                  .(max_reads = max(nreads),
                                    max_cells = max(ncells)),
                                  by = .(bid)]
all.support.flt <- all.support.sum[max_reads > 2 & max_cells > 2,]

all.sng <- all_bs_ov$binsets[all_bs_ov$binsets$bid %in% all.support.flt$bid,]
unique(all.sng$bid) %>% length()

saveRDS(all_synergy, paste0("mono_",prefix,"_synergy_result.rds"))
saveRDS(all.sng, paste0("mono_",prefix,"_sng.rds"))

## export bed --------------------------------------------------------------

library(RColorBrewer)
all.sng.bed <- bs2bed(all.sng,brewer.pal(10,"Set3")[-2])

write.table(all.sng.bed, paste0("mono_",prefix,"_20k_synergies.bed"),
            sep = "\t",quote = F,row.names = F,col.names = F)