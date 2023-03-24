source('synergy_funcs.R')

prefix="all"
mc.cores.used = 10
hg38.ref <- readRDS('hg38_seqinfo.rds')
bins <- readRDS('scCh_20k_bins.rds')
gc_frag_cov <- readRDS('gc_frag_cov.rds')
targets_EP <- readRDS("GM12878_chromHMM_20k_RE_targets.rds")
EP.bins <- make_window_bins(windows = targets_EP)

all.gr <- readRDS('mono_all_gr.rds')

all.EP = all.gr[gr.in(all.gr,targets_EP)]
all.EP$cid = all.EP$read_idx

re.cc <- get_chromunity(concatemers = all.EP, EP.bins,
                        piecewise = FALSE, shave = TRUE,
                        seed = 42, verbose = TRUE,
                        cthresh = 3, bthresh = 2,
                        mc.cores = mc.cores.used)
re.cc$read_idx <- re.cc$cid
re.cc$cell_idx <- as.numeric(factor(re.cc$cell))

re.bs <- find_ov_binsets(re.cc, unique(re.cc$chid),
                        bins = EP.bins$bins,
                        mc.cores = mc.cores.used)
unique(re.bs$binsets$bid) %>% length()

re.bs.scann <- sc_annotate(re.bs$binsets,concatemers = all.EP,
                           covariates = gc_frag_cov,resolution = 2e4,
                           mc.cores = mc.cores.used)

re_back_gr <- make_background(re.bs$binsets,resolution = 2e4)
re_back_gr.ann <- sc_annotate(re_back_gr,concatemers = all.EP,
                               covariates = gc_frag_cov,
                               resolution = 2e4,
                               mc.cores = mc.cores.used)

re_reads_model = fit(re_back_gr.ann, "nreads")
re_cells_model = fit(re_back_gr.ann, "ncells")

re_bs_ov.score <- sscore(re.bs.scann, re_reads_model, re_cells_model) 
re_synergy = synergy(re_bs_ov.score, re_reads_model, re_cells_model)
unique(re_synergy$binsets$bid) %>% length()

saveRDS(re.cc,"RE_cc.rds")
saveRDS(re.bs,"RE_bs.rds")
saveRDS(re.bs.scann,"RE_bs_scann.rds")
saveRDS(re_back_gr,"RE_back.rds")
saveRDS(re_back_gr.ann,"RE_back_ann.rds")
saveRDS(re_synergy,"RE_sng.rds")