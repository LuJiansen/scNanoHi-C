library(rtracklayer)
library(skitools)
library(chromunity)
library(MASS)
library(arrow)
library(Matrix)
library(dplyr)
library(parallel)
library(pbmcapply)
library(tidyverse)

make_window_bins = function(region = NULL,resolution = 5e4,
                            windows = NULL,window.size = 2e6,
                            stride = window.size/2,verbose = TRUE){
  if (is.null(windows))
    windows = gr.start(gr.tile(region, stride))+window.size/2
  
  windows = dt2gr(gr2dt(windows)[, start := ifelse(start < 0, 1, start)])
  
  if (inherits(windows, 'GRanges'))
    windows = split(windows, 1:length(windows))
  
  if (!is.null(resolution)){
    bins = gr.tile(GenomicRanges::reduce(gr.stripstrand(unlist(windows))), resolution)[, c()]
  } else {
    bins = windows %>% unlist %>% gr.stripstrand %>% disjoin
  }
  
  if (verbose) cmessage('Generated ', length(bins), ' bins across ', length(windows), ' windows')
  
  ## match window ids and bins 
  binmap = bins %*% grl.unlist(windows)[, c('grl.ix')] %>% 
    as.data.table %>% setnames('query.id', 'binid') %>% 
    setnames('grl.ix', 'winid') %>% setkeyv('winid')
  
  return(list(windows = windows,
              bins = bins,
              binmap = binmap))
}

get_chromunity = function(concatemers, bin_obj, max.size = 2^31-1,
                          subsample.frac = NULL, max.slice = 1e6,
                          min.reads = 3,min.cells = 3, mc.cores = 5,
                          k.knn = 100, k.min = 0,piecewise = TRUE,
                          shave = FALSE, bthresh = 2,cthresh = 3,
                          seed = 42, verbose = TRUE) {
  if (is.null(bin_obj))
  {
    stop("bins object is required")
  } else {
    windows = bin_obj$windows
    bins = bin_obj$bins
    binmap = bin_obj$binmap
  }
  
  if (is.null(concatemers$cid))
  {
    if ('read_idx' %in% names(values(concatemers)))
      names(values(concatemers))[match('read_idx', names(values(concatemers)))] = 'cid'
    else
      stop("concatemer GRanges must have metadata column $cid or $read_idx specifying the concatemer id")
  }
  
  params = data.table(k.knn = k.knn, k.min = k.min, seed = seed)
  
  if (verbose) cmessage('Matching concatemers with bins, and bins with windows using gr.match with max.slice ', max.slice, ' and ', mc.cores, ' cores')
  
  ## (batch) match up concatemers with binids
  concatemers$binid = gr.match(concatemers, bins, max.slice = max.slice,
                               mc.cores =  mc.cores, verbose = verbose)
  
  ## maybe NA need to be removed
  concatemers = concatemers %Q% (!is.na(binid))
  
  ## cycle through (possibly complex) windows call cluster_concatemers and convert to gr.sums
  ## winids = unique(binmap$winid)
  winids = unique(binmap[binid %in% unique(concatemers$binid)]$winid)
  
  if (verbose) cmessage('Starting concatemer community detection across ',
                        length(winids), ' windows')
  if(piecewise){
    ##cc = mclapply(winids, mc.cores = mc.cores, mc.preschedule = TRUE, function(win)
    cc = pbmclapply(winids, mc.cores = mc.cores, mc.preschedule = TRUE, function(win)
    {
      print(win)
      suppressWarnings({
        these.bins = binmap[.(win), ]
        cc = sc_concatemer_communities(concatemers %Q% (binid %in% these.bins$binid),
                                       k.knn = k.knn, max.size = max.size,
                                       k.min = k.min, seed = seed, verbose = verbose>1)
        if (length(cc))
        {
          cc = cc[cc$nreads >= min.reads & cc$ncell >= min.cells]
          if (length(cc)) {
            cc$winid = win
          }
        }
      })
      cc
    })
    
    cc = cc %>% do.call(grbind, .)
    cc = dt2gr(gr2dt(cc)[, chid := .GRP, by = .(chid, winid)])
  } else {
    if (shave)
    {
      if (verbose)
        cmessage('Shaving concatemers with bthresh = ', bthresh, ' and cthresh = ', cthresh)
      concatemers = shave_concatemers(concatemers, bthresh = bthresh,
                                      cthresh = cthresh, verbose = verbose)
    }
    
    ncat = concatemers$cid %>% unique %>% length
    nbin = concatemers$binid %>% unique %>% length
    
    if (verbose)
      cmessage(sprintf('Running concatemer communities with %s concatemers and %s bins', ncat, nbin))
    
    cc = sc_concatemer_communities(concatemers,
                                   k.knn = k.knn, k.min = k.min,
                                   seed = seed, max.size = max.size,
                                   verbose = verbose,
                                   subsample.frac = subsample.frac)
    
    if (length(cc))
      cc = cc[cc$nreads >= min.reads & cc$ncell >= min.cells]
      cc$winid <- 1
  }
  return(cc)
}

shave_concatemers = function(concatemers, cthresh = 3, bthresh = 2,
                             verbose = TRUE)
{
  
  .shave = function(concatemers, bthresh = 2, cthresh = 2)
  {
    dt = unique(gr2dt(concatemers), by = c('cid', 'binid'))
    ccount = dt[, .N, by = cid] # number of bins in each concatemer
    bcount = dt[, .N, by = binid] # number of concatemers in each bin
    coolc = ccount[N>=bthresh, cid]
    coolb = bcount[N>=cthresh, binid]
    concatemers %Q% (cid %in% coolc) %Q% (binid %in% coolb)
  }
  old = concatemers
  new = .shave(concatemers, bthresh = bthresh, cthresh = cthresh)
  while (length(new)<length(old))
  {
    if (verbose)
      cmessage('shaving --> Diff: ', length(old)-length(new), ', Old: ', length(old), ', New:', length(new))
    old = new;
    new = .shave(old, bthresh = bthresh, cthresh = cthresh)
  }
  new
}

sc_concatemer_communities = function (concatemers, k.knn = 100, k.min = 0,
                                      drop.small = FALSE, small = 1e4, 
                                      max.size = 2^31-1, subsample.frac = NULL, 
                                      seed = 42, verbose = TRUE, debug = FALSE)  
{
  reads = concatemers
  
  if (is.null(reads$cid))
  {
    if ('read_idx' %in% names(values(reads)))
      names(values(reads))[match('read_idx', names(values(reads)))] = 'cid'
    else
      stop("concatemer GRanges must have metadata column $cid or $read_idx specifying the concatemer id")
  }
  
  if (drop.small) {
    if (verbose) cmessage(paste0("Filtering out reads < ", small))
    reads = gr2dt(reads)
    setkeyv(reads, c("seqnames", "start"))
    reads[, `:=`(max.local.dist, end[.N] - start[1]), by = cid]
    reads = reads[max.local.dist > small]
    reads = dt2gr(reads)
  }
  
  if (debug)
    browser()
  
  if (verbose) cmessage("Matching reads to tiles")
  reads = as.data.table(reads)[, `:=`(count, .N), by = cid]
  reads2 = reads[count > 2, ] 
  
  if (!nrow(reads2))
  {
    warning('no high order concatemers, returning empty result')
    return(reads[, chid := NA][c(), ])
  }
  
  reads2$cid = factor(reads2$cid)
  ucid = levels(reads2$cid)
  
  if (!is.null(subsample.frac))
  {
    if (verbose) cmessage("Using fraction subsampling")
    setkey(reads2, cid)
    reads2 = reads2[.(sample(ucid, length(ucid)*subsample.frac)), ]
    reads2$cid = factor(reads2$cid)
    ucid = levels(reads2$cid)
  }
  
  
  if (verbose) cmessage("Matrices made")
  ## gc()
  
  ## remove bins hit only by one concatemer
  reads2$binid = factor(reads2$binid)
  setkey(reads2, binid)
  ubid = reads2[, .N, by = binid][N>1 , binid]
  
  if (!length(ubid))
  {
    warning('No bins hit by two concatemers, returning empty result')
    return(reads[, chid := NA][c(), ])
  }
  
  reads2 = reads2[.(ubid), ]
  ## added for subsamp
  reads2$binid = factor(reads2$binid) 
  
  ## refactor may be necessary
  reads2$cid = factor(reads2$cid)
  ucid = levels(reads2$cid)
  
  ## size is the number of pairs which can't exceed the max.size (or integer max)
  size = data.table(cid = reads2$cid, binid = reads2$binid)[, choose(.N,2), by = binid][, sum(as.numeric(V1))]
  
  ncat = reads2$cid %>% unique %>% length
  nbin = reads2$binid %>% unique %>% length
  
  if (size > max.size)
    stop(sprintf('size %s of the problem %s with %s concatemers and %s bins  exceeds max.size %s, please considering subsampling concatemers, using fewer windows or bins, or shaving concatemers with more aggressive parameters (bthresh, cthresh)', size, ifelse(shave, 'after shaving', ''), ncat, nbin, max.size))
  
  
  ## all pairs of concatemers that share a bin
  reads2[, cidi := as.integer(cid)]
  
  if (verbose) cmessage("Making Pairs object") 
  pairs = reads2[, t(combn(cidi, 2)) %>% as.data.table, by = binid]
  if (verbose) cmessage("Pairs object made") 
  
  setkey(reads2, cid)
  concatm = reads2[.(ucid),  sparseMatrix(factor(cid, ucid) %>% as.integer, binid %>% as.integer, x = 1, dimnames = list(ucid, levels(binid)))]
  
  p1 = concatm[pairs[, V1], -1, drop = FALSE]
  p2 = concatm[pairs[, V2], -1, drop = FALSE]
  matching = Matrix::rowSums(p1 & p2)
  total = Matrix::rowSums(p1 | p2)
  dt = data.table(bx1 = pairs[, V1], bx2 = pairs[, V2], mat = matching, 
                  tot = total)[, `:=`(frac, mat/tot)]
  dt2 = copy(dt)
  dt2$bx2 = dt$bx1
  dt2$bx1 = dt$bx2
  dt3 = rbind(dt, dt2)
  dt3$nmat = dt3$mat
  dt3$nfrac = dt3$frac
  setkeyv(dt3, c("nfrac", "nmat"))
  dt3 = unique(dt3)
  dt3.2 = dt3[order(nfrac, nmat, decreasing = T)]
  if (verbose) cmessage("Pairs made")
  gc()
  k = k.knn
  knn.dt = dt3.2[mat > 2 & tot > 2, .(knn = bx2[1:k]), by = bx1][!is.na(knn), ]
  if (!nrow(knn.dt))
  {
    warning('no concatemers with neighbors, perhaps data too sparse? returning empty result')
    return(reads[, chid := NA][c(), ])
  }
  setkey(knn.dt)
  knn = sparseMatrix(knn.dt$bx1, knn.dt$knn, x = 1)
  knn.shared = knn %*% knn
  if (verbose) cmessage("KNN done")
  KMIN = k.min
  A = knn.shared * sign(knn.shared > KMIN)
  A[cbind(1:nrow(A), 1:nrow(A))] = 0
  #  A <- as(A, "sparseMatrix")
  A = A + t(A)
  G = graph.adjacency(A, weighted = TRUE, mode = "undirected")
  cl.l = cluster_fast_greedy(G)
  cl = cl.l$membership
  ## rename so chid has most support
  cls = 1:max(cl)
  names(cls) = cl %>% table %>% sort %>% rev %>% names
  cl = cls[as.character(cl)]
  if (verbose) cmessage("Communities made")
  memb.dt = data.table(cid = ucid[1:nrow(A)], chid = cl)
  reads[, cid := as.character(cid)]
  reads = merge(reads, memb.dt, by = "cid")
  reads[, c("nreads","ncell") := list(uniqueN(cid),uniqueN(cell)), by = chid]
  reads = dt2gr(reads)
  return(reads)
}

find_ov_binsets <- function(cc,uchid,bins,
                            mc.cores,pad=0,
                            min.reads=3,
                            min.cells=3,
                            min.bins=3,
                            reduce_binset = T){
  if (is.null(cc$read_idx) | is.null(cc$cell_idx))
  {
      stop("concatemer GRanges must have metadata column $read_idx and $cell_idx")
  }
  
  bins$binid <- 1:length(bins)
  result = pbmclapply(uchid, mc.cores = mc.cores, function(this.chid)
  {
    print(this.chid)
    {
      this.cc = cc %Q% (chid == this.chid)
      cc.pad <- this.cc + pad
      cc.pad@seqinfo <- bins@seqinfo
      cc.pad@seqnames <- Rle(factor(cc.pad@seqnames,
                                    levels = seqlevels(seqinfo(cc.pad))))
      
      ov <- findOverlaps(bins,cc.pad)
      ov.cc <- cc.pad[ov@to,]
      ov.cc$binid <- ov@from
      bin.sum <- as.data.frame(ov.cc) %>% group_by(binid) %>%
        summarise(nreads = n_distinct(read_idx),
                  ncell = n_distinct(cell_idx)) %>%
        filter(nreads >= min.reads & ncell >= min.cells)
      
      # stop when available bins < min.bins
      if (nrow(bin.sum)>=min.bins)
      {
        bins.used <- bins[bin.sum$binid,]
        bins.used$ncell <- bin.sum$ncell
        bins.used$nreads <- bin.sum$nreads
        
        if (length(bins.used)) {
          cc.used <- as.data.frame(ov.cc[ov.cc$binid %in% bins.used$binid])
          keep.col = c("seqnames","start","end","width","strand","cid",
                       "cell","binid","chid","winid")
          cc.used <- merge(cc.used[,keep.col],
                           as.data.frame(mcols(bins.used)),
                           by = "binid") %>% dt2gr()
          cc.used$bid <- this.chid
          binset = bins.used
        } else {
          binset = NULL
        }
        
        if (length(binset))
        {
          binset$chid = this.chid
          binset$bid = this.chid
          binset$winid = this.cc$winid[1]
          if (reduce_binset)
          {
            binset = gr.reduce(binset)
          }
        }
        if (!is.null(binset) & length(binset) >= min.bins) {
          return(list(binset = binset,
                      cc = cc.used))
        } else {
          return(NULL)
        }
      } else {
        return(NULL)
      }
    }
  })
  binsets <- lapply(result, function(x) x$binset) %>% do.call(grbind, .)
  binsets <- as.data.frame(binsets) %>%
    mutate(binid = as.numeric(as.factor(paste(seqnames,start,end,sep = "_")))) %>%
    dt2gr()
  
  # remove binset duplicates
  bid.list <- split(binsets$binid,binsets$bid)
  bid.uniq <- names(bid.list[!duplicated(bid.list)])
  
  binsets <- binsets[binsets$bid %in% bid.uniq,]
  cc.used <- lapply(result, function(x) x$cc) %>% do.call(grbind, .)
  cc.used <- cc.used[cc.used$bid %in% bid.uniq,]
  return(list(binsets = binsets,
              cc = cc.used))
}

sc_annotate = function(binsets, concatemers, covariates = NULL,
                       k = 5, interchromosomal.dist = 1e8,
                       interchromosomal.table = NULL,
                       gg = NULL, mc.cores = 5,
                       numchunks = 200*mc.cores-1,
                       seed = 42, verbose = TRUE,
                       unique.per.setid = TRUE, resolution = 5e4)
{
  set.seed(seed)
  
  if (!inherits(binsets, 'GRanges'))
    binsets = dt2gr(binsets)
  
  binsets <- gr2dt(binsets) %>% 
    mutate(binid = as.numeric(as.factor(paste(seqnames,start,end,sep = "_")))) %>% 
    dt2gr()
  binsets = gr.stripstrand(binsets)
  binsets$bid = as.factor(binsets$bid)
  ## each row of binsets is a bin
  #binsets$binid = 1:length(binsets)
  #####
  ## frontload / batch some useful calculations
  #####
  ## bin vs concatemer overlaps
  if (verbose) smessage('Overlapping ', length(binsets),
                        ' bins with ', length(concatemers), ' monomers')
  ov = binsets %*% concatemers[, c('cid','cell')] %>% as.data.table
  ##
  if (verbose) smessage('Computing bin by bin pairwise distance')
  ## bin x bin pairwise distance within each set
  if (is.null(gg)){
    bindist = gr2dt(binsets)[, as.data.table(expand.grid(i = binid, j = binid))[i<j, ], by = bid] %>% setkeyv(c('i', 'j'))
    bindist[, distance := GenomicRanges::distance(binsets[i], binsets[j])]
  } else {
    if (verbose) smessage('Using graph distance')
    ######
    bindist = gr2dt(binsets)[, as.data.table(expand.grid(i = binid, j = binid))[i<j, ], by = bid] %>% setkeyv(c('i', 'j'))
    bindist[, distance := GenomicRanges::distance(binsets[i], binsets[j])]
    bindist.intra = bindist[!is.na(distance)]    
    ######
    gg = gg$copy$disjoin(disjoin(binsets))
    binsetd = data.table(binid = binsets$binid, bid = binsets$bid) %>% setkey('binid')
    binsetd[, sun.bin.id := 1:.N, by = bid]
    
    ## Better approach to get graph distances, need further testing
    bindist.g = pbmclapply(unique(binsetd$bid), function(this.bid){
      this.dist = tryCatch(gg$dist(gr.nochr(binsets[binsetd[bid == this.bid, binid]])) %>% melt %>% as.data.table %>% setnames(., c('gi', 'gj', 'distance')), error = function(e) NULL)
      if(!is.null(this.dist)){
        this.dist[, bid := as.factor(this.bid)]
        return(this.dist)
      }
    }, mc.cores = mc.cores) %>% rbindlist
    bindist = merge(bindist.g, binsetd, by.x = c("gi", "bid"), by.y = c("sun.bin.id", "bid"), allow.cartesian=TRUE)
    setnames(bindist, "binid", "i")
    bindist = merge(bindist, binsetd, by.x = c("gj", "bid"), by.y = c("sun.bin.id", "bid"), allow.cartesian=TRUE)
    setnames(bindist, "binid", "j")
    bindist = unique(bindist[, .(bid, i, j, distance)][i<j])
    bindist = merge(bindist, bindist.intra, by = c("bid", "i", "j"), all.x = T)
    bindist[, distance := ifelse(distance.x <= 1, distance.y, distance.x)]  
    bindist[, distance := round_any(distance, resolution)]
    bindist = bindist[, .(bid, i , j, distance)]  
    ## bindist[, distance := ifelse(distance < resolution, resolution, distance)]
    setkeyv(bindist, c('i', 'j'))
  }
  if (!is.null(interchromosomal.table)){
    interchromosomal.table[, dist := round_any(dist, resolution)]
    bindist = merge(bindist, gr2dt(binsets)[, .(seqnames, binid)], by.x = "i", by.y = "binid")
    bindist = merge(bindist, gr2dt(binsets)[, .(seqnames, binid)], by.x = "j", by.y = "binid")
    setnames(bindist, c("seqnames.x", "seqnames.y"), c("V1", "V2"))
    bindist = merge(bindist, interchromosomal.table, by = c("V1", "V2"), all.x = T)
    bindist[, distance := ifelse(is.na(distance), dist, distance)]
    bindist = bindist[, .(j,   i,  bid,    distance)]
    bindist[, distance := ifelse(is.na(distance), max(interchromosomal.table$dist), distance)]
    setkeyv(bindist, c('i', 'j'))
  } else {
    bindist[is.infinite(distance), distance := interchromosomal.dist]
    bindist[is.na(distance), distance := interchromosomal.dist]
  }
  #####
  ## now start annotating sub.binsets aka sets
  #####
  if (verbose) smessage('Making sub binsets')
  ## make all sub power sets up to k for all bids
  ## each row is now a binid with a setid and bid
  ## adding ones for sum.cov to be removed later
  sub.binsets = gr2dt(binsets)[, powerset(binid, 1, k), by = bid] %>% setnames(c('bid', 'setid', 'binid')) %>% setkey(bid)
  sub.binsets[, ":="(iid = 1:.N, tot = .N), by = .(setid, bid)] ## label each item in each sub-binset, and total count will be useful below
  
  if (verbose) smessage('Made ', nrow(sub.binsets), ' sub-binsets')
  
  ## first use ov to count how many concatemers fully overlap all the bins in the subbinset
  if (verbose) smessage('Counting concatemers across sub-binsets across ', mc.cores, ' cores')
  ref.counts = unique(sub.binsets[, .(bid, setid)])[, c("nreads","ncells") := list(0,0)]
  if (nrow(ov)){
    ubid = unique(sub.binsets$bid) ## split up to lists to leverage pbmclapply
    ubidl = split(ubid, ceiling(runif(length(ubid))*numchunks)) ## randomly chop up ubid into twice the number of mc.coreso
    counts = pbmclapply(ubidl, mc.cores = mc.cores, function(bids)
    {
      out = tryCatch(merge(sub.binsets[.(as.factor(bids)), ],
                           ov[bid %in% bids], by = c('binid', 'bid'),
                           allow.cartesian = TRUE),
                     error = function(e) NULL)
      if (!is.null(out) & nrow(out) > 0){
        if (unique.per.setid){
          out = out[, .(binid, bid, setid, iid, tot, cid, cell)]
          
          # conut number of fully overlapping concatemers
          cid.tmp = out[, .(hit = all(1:tot[1] %in% iid)),
                        by = .(cid, cell, setid, bid, tot)][hit == TRUE]
          setkeyv(cid.tmp, c("bid", "cid", "cell", "tot"))
          
          # for each bid, each concatemer only count once for the longest sub-binset it overlapping 
          cid.tmp = cid.tmp[, tail(.SD, 1), by = .(cid, bid)]
          cid.counts = cid.tmp[, .(nreads = sum(hit, na.rm = TRUE)),
                               by = .(setid, bid)]
          
          # count number of unique fully overlapping cells
          cell.counts = cid.tmp[, .(ncells = uniqueN(cell)),
                                by = .(setid, bid)]
          
          # cell.tmp = out[, .(hit = all(1:tot[1] %in% iid)),
          #                by = .(cell, setid, bid, tot)][hit == TRUE]
          # setkeyv(cell.tmp, c("bid", "cell","tot"))
          # cell.tmp = cell.tmp[, tail(.SD, 1), by = .(cell, bid)]
          # cell.counts = cell.tmp[, .(ncells = sum(hit, na.rm = TRUE)),
          #                        by = .(setid, bid)]
          
          this.counts = merge(ref.counts[bid %in% bids],
                              cid.counts, by = c("bid", "setid"),
                              all.x = T)
          this.counts = merge(this.counts[bid %in% bids],
                              cell.counts, by = c("bid", "setid"),
                              all.x = T)
          this.counts[, c("nreads","ncells") := list(sum(nreads.x, nreads.y, na.rm = T),
                                                     sum(ncells.x, ncells.y, na.rm = T)),
                      by = .(bid, setid)][, .(bid, setid, nreads, ncells)]
        } else {
          out[, .(hit = all(1:tot[1] %in% iid)),
              by = .(cid, cell, setid, bid)] %>%
            .[(hit)] %>%
            .[, .(nreads = sum(hit, na.rm = TRUE),
                  ncells = uniqueN(cell)),
              by = .(setid, bid)]
        }
      } else {
        NULL
      }
    })  %>% rbindlist
  }
  
  ## other goodies
  ## changing to mean instead of median
  if (verbose) smessage('Computing min median max distances per setid')
  ## dists = sub.binsets[, bindist[as.data.table(expand.grid(i = binid, j = binid))[i<j, ], .(dist = c('min.dist', 'mean.dist', 'max.dist'), value = as.numeric(summary(distance+1, na.rm = T)[c(1, 4, 6)]))], by = .(setid, bid)] %>% dcast(bid + setid ~ dist, value.var = 'value')
  
  sub.binsets = sub.binsets[tot > 1]
  ubid = unique(sub.binsets$bid) ## split up to lists to leverage pbmclapply
  ubidl = split(ubid, ceiling(runif(length(ubid))*numchunks)) ## randomly chop up ubid into twice the number of mc.coreso
  dists = pbmclapply(ubidl, mc.cores = mc.cores, function(bids)
  {
    this.sub.binsets = sub.binsets[bid %in% bids]
    this.dists = this.sub.binsets[, bindist[as.data.table(expand.grid(i = binid, j = binid))[i<j, ], .(dist = c('min.dist',  'max.dist'), value = quantile(distance+1,  c(0, 1)))], by = .(setid, bid)] %>% dcast(bid + setid ~ dist, value.var = 'value')
  })  %>% rbindlist
  
  
  if (verbose) smessage('Computing marginal sum per bid')
  margs = counts[, .(sum.nreads = sum(nreads),
                     sum.ncells = sum(ncells)),
                 by = bid]
  
  if (verbose) smessage('Computing total width and cardinality per setid')
  ubid = unique(sub.binsets$bid) ## split up to lists to leverage pbmclapply
  ubidl = split(ubid, ceiling(runif(length(ubid))*numchunks)) ## randomly chop up ubid into twice the number of mc.coreso
  widths = pbmclapply(ubidl, mc.cores = mc.cores, function(bids){
    this.sub.binsets = sub.binsets[bid %in% bids]
    this.widths = this.sub.binsets[, .(width = sum(width(binsets)[binid])+1, cardinality = .N), by = .(bid, setid)]
  })  %>% rbindlist
  
  
  if (verbose) smessage('Merging counts, distance, and width')
  annotated.binsets = unique(sub.binsets[, .(bid, setid)]) %>%
    merge(counts, by = c('bid', 'setid'),allow.cartesian=TRUE) %>%
    merge(dists, all.x = TRUE, by = c('bid', 'setid'),allow.cartesian=TRUE) %>%
    merge(widths, all.x = TRUE, by = c('bid', 'setid'),allow.cartesian=TRUE) %>%
    merge(margs, all.x = TRUE, by = c('bid'),allow.cartesian=TRUE)  
  
  annotated.binsets = annotated.binsets[cardinality > 1]
  
  ## now cycle through interval / numeric covariates
  if (!is.null(covariates))
  {
    setkeyv(sub.binsets, c('bid', 'binid'))
    if (verbose) smessage('Adding covariates')
    for (i in 1:length(covariates))
    {
      if (covariates$type[i] == 'numeric')
      {
        dat = covariates$data[[i]][, covariates$field[[i]]]
        names(values(dat)) = 'val'
      }
      else
        dat = covariates$data[[i]][, c(0)]
      
      if (verbose) smessage('Overlapping ', length(dat), ' ranges from covariate ', covariates$names[i])
      ov = binsets %*% dat %>% as.data.table
      
      if (verbose) smessage('Tallying values from ', covariates$type[i], ' covariate ', covariates$names[i])
      covdata = data.table(val = numeric(), bid = factor(), setid = numeric()) %>% setnames('val', covariates$names[i])
      
      if (nrow(ov))
      {
        sov = merge(sub.binsets, ov, by = c('bid', 'binid'), allow.cartesian = TRUE)
        
        if (covariates[i]$type == 'interval') ## count vs ?sum width TODO: add flag on covariates specifying whether to count
        {
          covdata = sov[ , .(val = .N), by = .(bid, setid)]
          ## covdata = ov[sub.binsets, , allow.cartesian = TRUE][, .(val = sum(width, na.rm = TRUE)), by = .(bid, setid)]          
        }
        else if (covariates[i]$type == 'numeric') ## average value
          covdata = sov[, .(val = sum(width*val, na.rm = TRUE)/sum(width + 0*val, na.rm = TRUE)), by = .(bid, setid)]
        
        ## eps is small amount to add to covariate to make it loggable
        eps = range(covdata$val, na.rm = TRUE) %>% diff * 1e-4
        covdata$val = covdata$val + eps
        
        if (any(is.na(covdata$val)) || any(covdata$val<=0))
          warning(paste('Covariate ', covariates$names[[i]], ' has NA or <= values, check or fix before fitting '))
        setnames(covdata, 'val', covariates$names[i])
      }
      
      if (verbose) smessage('Merging data from ', covariates$type[i], ' covariate ', covariates$names[i])
      annotated.binsets = merge(annotated.binsets, covdata, all.x = TRUE, by = c('bid', 'setid'))
    }
  }
  return(annotated.binsets)
}

make_background <- function(binsets,n=NULL,resolution=2e4,
                            genome = "BSgenome.Hsapiens.UCSC.hg38::Hsapiens"){
  if(is.null(n)){
    n = signif(length(unique(binsets$bid)) + 500,1)
  }
  message("Genome: ",genome)
  message("resolution: ",resolution)
  message("Generating ",n," background binsets ...")
  back_gr = gr2dt(dt2gr(background(binsets = binsets, n = n,
                                   resolution = resolution)))
  
  upper.bound = as.data.table(hg_seqlengths(genome = genome), keep.rownames = T) 
  setkeyv(back_gr, c("seqnames", "start"))
  back_gr = back_gr[!bid %in% back_gr[width < (resolution-1)]$bid]
  back_gr = gr2dt(gr.reduce(dt2gr(back_gr), by = "bid"))
  back_gr$bid <- as.factor(back_gr$bid)
  back_gr = merge(back_gr, upper.bound, by.x = "seqnames", by.y = "V1",
                  all.x = T, allow.cartesian = T)
  back_gr = back_gr[end < V2][start < V2]
  back_gr[, overall.cardinality := .N, by = bid]
  back_gr = back_gr[overall.cardinality > 1]
  back_gr = dt2gr(back_gr)
  return(back_gr)
}

fit = function(annotated.binsets,response = "count",
               covariates = NULL, nb = TRUE,
               return.model = TRUE, verbose = TRUE, maxit = 50)
{
  if (!nb) stop('not yet supported')
  ## added sumcounts as cov and width as only offset
  if (is.null(covariates)) {
    covariates = setdiff(names(annotated.binsets), c('bid', 'setid', 'mean.dist', 'ncells', "nreads"))
  }
  fmstring = paste(response,'~', paste(paste0('log(', covariates, '+ 1)', collapse = ' + ')))
  ## fmstring = paste0(fmstring, " + ", "offset(log(width))")
  print(fmstring)
  fm = formula(fmstring)
  ##
  model = tryCatch(glm.nb(formula = fm, data = annotated.binsets,
                          control = glm.control(maxit = maxit)),
                   error = function(e) NULL)
  
  return(list(model = model, covariates = covariates))
}

sscore = function(annotated.binsets, reads_model, cells_model, verbose = TRUE)
{
  if (is.null(annotated.binsets$nreads) | is.null(annotated.binsets$ncells))
    stop('annotsynated.binsets need to have $nreads and $ncells column, did you forget to annotate?')
  
  if (length(missing1 <- setdiff(reads_model$covariates, names(annotated.binsets))) |
      length(missing2 <- setdiff(cells_model$covariates, names(annotated.binsets))))
    stop('annotated.binsets missing covariates: ', paste(missing1, missing2, collapse = ', '))
  
  annotated.binsets$nreads.predicted = (predict(reads_model$model, type = "response", newdata = annotated.binsets))
  annotated.binsets$ncells.predicted = (predict(cells_model$model, type = "response", newdata = annotated.binsets))
  
  ## randomized p value procedure for negative binomial
  r.pval = annotated.binsets[, pnbinom(nreads - 1, mu = nreads.predicted, size = reads_model$model$theta, lower.tail = F)]
  r.pval.right = annotated.binsets[, pnbinom(nreads, mu = nreads.predicted, size = reads_model$model$theta, lower.tail = F)]
  r.pval.right = ifelse(is.na(r.pval.right), 1, r.pval.right)
  r.pval = ifelse(is.na(r.pval), 1, r.pval)
  
  c.pval = annotated.binsets[, pnbinom(nreads - 1, mu = nreads.predicted, size = reads_model$model$theta, lower.tail = F)]
  c.pval.right = annotated.binsets[, pnbinom(nreads, mu = nreads.predicted, size = reads_model$model$theta, lower.tail = F)]
  c.pval.right = ifelse(is.na(c.pval.right), 1, c.pval.right)
  c.pval = ifelse(is.na(c.pval), 1, c.pval)
  
  annotated.binsets[, c("nreads_enrichment","ncells_enrichment") := list(nreads / nreads.predicted, ncells / ncells.predicted)]
  annotated.binsets$nreads_pval = runif(nrow(annotated.binsets), min = r.pval.right, max = r.pval)
  annotated.binsets$ncells_pval = runif(nrow(annotated.binsets), min = c.pval.right, max = c.pval)
  return(annotated.binsets)
}

glm.nb.fh = function (formula, data, weights, subset, na.action, start = NULL,
                      etastart, mustart, control = glm.control(...), method = "glm.fit",
                      model = TRUE, x = FALSE, y = TRUE, contrasts = NULL, ...,  theta = NULL,
                      init.theta, link = log)
{
  loglik <- function(n, th, mu, y, w) sum(w * (lgamma(th +
                                                        y) - lgamma(th) - lgamma(y + 1) + th * log(th) + y *
                                                 log(mu + (y == 0)) - (th + y) * log(th + mu)))
  link <- substitute(link)
  fam0 <- if (missing(init.theta))
    do.call("poisson", list(link = link))
  else do.call("negative.binomial", list(theta = init.theta,
                                         link = link))
  mf <- Call <- match.call()
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "etastart", "mustart", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval.parent(mf)
  Terms <- attr(mf, "terms")
  if (method == "model.frame")
    return(mf)
  Y <- model.response(mf, "numeric")
  X <- if (!is.empty.model(Terms))
    model.matrix(Terms, mf, contrasts)
  else matrix(, NROW(Y), 0)
  w <- model.weights(mf)
  if (!length(w))
    w <- rep(1, nrow(mf))
  else if (any(w < 0))
    stop("negative weights not allowed")
  offset <- model.offset(mf)
  mustart <- model.extract(mf, "mustart")
  etastart <- model.extract(mf, "etastart")
  n <- length(Y)
  if (!missing(method)) {
    if (!exists(method, mode = "function"))
      stop(gettextf("unimplemented method: %s", sQuote(method)),
           domain = NA)
    glm.fitter <- get(method)
  }
  else {
    method <- "glm.fit"
    glm.fitter <- stats::glm.fit
  }
  if (control$trace > 1)
    message("Initial fit:")
  
  fit <- glm.fitter(x = X, y = Y, w = w, start = start, etastart = etastart,
                    mustart = mustart, offset = offset, family = fam0, control = list(maxit = control$maxit,
                                                                                      epsilon = control$epsilon, trace = control$trace >
                                                                                        1), intercept = attr(Terms, "intercept") > 0)
  class(fit) <- c("glm", "lm")
  mu <- fit$fitted.values
  if (is.null(theta))
  {
    th <- as.vector(theta.ml(Y, mu, sum(w), w, limit = control$maxit,
                             trace = control$trace > 2))
  }
  else
  {
    th = theta
    
    if (control$trace > 1)
      message(gettextf("Fixing theta value to 'theta': %f", signif(th)),
              domain = NA)
  }
  
  if (control$trace > 1)
    message(gettextf("Initial value for 'theta': %f", signif(th)),
            domain = NA)
  fam <- do.call("negative.binomial", list(theta = th[1], link = link))
  iter <- 0
  d1 <- sqrt(2 * max(1, fit$df.residual))
  d2 <- del <- 1
  g <- fam$linkfun
  Lm <- loglik(n, th, mu, Y, w)
  Lm0 <- Lm + 2 * d1
  while (
    ((iter <- iter + 1) <= control$maxit) &
    ((abs(Lm0 - Lm)/d1 + abs(del)/d2) > control$epsilon)
  ){
    eta <- g(mu)
    fit <- glm.fitter(x = X, y = Y, w = w, etastart = eta,
                      offset = offset, family = fam, control = list(maxit = control$maxit,
                                                                    epsilon = control$epsilon, trace = control$trace >
                                                                      1), intercept = attr(Terms, "intercept") >
                        0)
    t0 <- th
    if (is.null(theta))
    {
      th <- theta.ml(Y, mu, sum(w), w, limit = control$maxit,
                     trace = control$trace > 2)
    } else
    {
      th = theta
    }
    
    fam <- do.call("negative.binomial", list(theta = th[1],  ## we don't need all the thetas here if theta is vectorized
                                             link = link))
    
    mu <- fit$fitted.values
    del <- t0 - th ## this is where the vectorized theta makes a difference
    Lm0 <- Lm
    Lm <- loglik(n, th, mu, Y, w) ## and here - log likelihood computation
    if (control$trace) {
      Ls <- loglik(n, th, Y, Y, w)
      Dev <- 2 * (Ls - Lm)
      message(sprintf("Theta(%d) = %f, 2(Ls - Lm) = %f",
                      iter, signif(th), signif(Dev)), domain = NA)
    }
  }
  if (!is.null(attr(th, "warn")))
    fit$th.warn <- attr(th, "warn")
  if (iter > control$maxit) {
    warning("alternation limit reached")
    fit$th.warn <- gettext("alternation limit reached")
  }
  if (length(offset) && attr(Terms, "intercept")) {
    null.deviance <- if (length(Terms))
      glm.fitter(X[, "(Intercept)", drop = FALSE], Y, w,
                 offset = offset, family = fam, control = list(maxit = control$maxit,
                                                               epsilon = control$epsilon, trace = control$trace >
                                                                 1), intercept = TRUE)$deviance
    else fit$deviance
    fit$null.deviance <- null.deviance
  }
  class(fit) <- c("negbin", "glm", "lm")
  fit$terms <- Terms
  fit$formula <- as.vector(attr(Terms, "formula"))
  Call$init.theta <- signif(as.vector(th), 10)
  Call$link <- link
  fit$call <- Call
  if (model)
    fit$model <- mf
  fit$na.action <- attr(mf, "na.action")
  if (x)
    fit$x <- X
  if (!y)
    fit$y <- NULL
  fit$theta <- as.vector(th)
  fit$SE.theta <- attr(th, "SE")
  fit$twologlik <- as.vector(2 * Lm)
  fit$aic <- -fit$twologlik + 2 * fit$rank + 2
  fit$contrasts <- attr(X, "contrasts")
  fit$xlevels <- .getXlevels(Terms, mf)
  fit$method <- method
  fit$control <- control
  fit$offset <- offset
  fit
}

synergy = function(scored.binsets, reads_model, cells_model,
                   mc.cores = 5, verbose = TRUE,maxit = 50,
                   p_threshold = 0.05,
                   min.cells = 3, min.reads = 3,
                   p_adjust_methods = "BH")
{
  if (verbose) smessage('Scoring binsets')
  
  scored.binsets[, multiway := cardinality > 2]
  setkey(scored.binsets, bid)
  ubid = unique(scored.binsets$bid)
  
  r.res = pbmclapply(ubid, function(this.bid){
    muffle(
      dflm(glm.nb.fh(data = scored.binsets[.(this.bid),],
                     nreads ~ multiway + offset(log(nreads.predicted)),
                     theta = reads_model$model$theta,
                     control = glm.control(maxit = maxit)))[2, ][, name := this.bid]
      )
  }, mc.cores = mc.cores) %>% rbindlist
  #setnames(r.res, 'name', 'bid')
  
  c.res = pbmclapply(ubid, function(this.bid){
    muffle(
      dflm(glm.nb.fh(data = scored.binsets[.(this.bid),],
                     ncells ~ multiway + offset(log(ncells.predicted)),
                     theta = cells_model$model$theta,
                     control = glm.control(maxit = maxit)))[2, ][, name := this.bid]
    )
  }, mc.cores = mc.cores) %>% rbindlist
  #setnames(c.res, 'name', 'bid')
  
  r.res$fdr = signif(p.adjust(r.res$p, p_adjust_methods), 2)
  c.res$fdr = signif(p.adjust(c.res$p, p_adjust_methods), 2)
  
  synergy.name <- intersect(r.res[r.res$fdr < p_threshold,]$name,
                            c.res[c.res$fdr < p_threshold,]$name)
  
  binsets.sig <- scored.binsets[bid %in% synergy.name,]
  return(list(nreads=r.res,
              ncells=c.res,
              binsets = binsets.sig))
}

bs2bed <- function(bs,color){
  names(bs) <- NULL
  bs.bed <- as.data.frame(bs) %>%
    arrange(as.numeric(bid),as.numeric(binid),as.numeric(winid))
  bs.bed$ID <- paste(bs.bed$winid,bs.bed$bid,sep = "_")
  #bs.bed$start <- bs.bed$start - 1
  bs.bed$width <- bs.bed$end - bs.bed$start
  df <- bs.bed %>% group_by(seqnames,ID,strand) %>% 
    arrange(start) %>%
    summarise(chromStart = min(start),
              chromEnd = max(end),
              blockCount = length(binid),
              blockSizes = paste(width,collapse = ","),
              blockStarts = paste(start - min(start),collapse = ","))
  df$color <- factor(df$ID,levels = unique(df$ID))
  nbs <- length(levels(df$color))
  levels(df$color) <- c(rep(color,floor(nbs/length(color))),
                        color[1:(nbs%%length(color))])
  df$strand <- "."
  df_bed <- df[,c("seqnames","chromStart","chromEnd","ID",
                  "blockCount","strand","chromStart","chromEnd",
                  "color","blockCount","blockSizes","blockStarts")]
  return(df_bed)
}


run_logit_enrich <- function(binsets,back,annDB){
  require(LOLA)
  require(epiDisplay)
  binsets$group <- 1
  back$group <- 0
  sng.c <- rbind(gr2dt(back)[,.(seqnames,start,end,bid,group)],
                 gr2dt(binsets)[,.(seqnames,start,end,bid,group)]) %>%
    dt2gr()
  
  message("Overlapping reference peaks with binsets ...")
  ann.ov <- pbmclapply(annDB, function(x){
    ov <- findOverlaps(sng.c,x)
    df <- as.data.frame(sng.c)
    df$hit <- 0
    df[unique(ov@from),]$hit <- 1
    df <- df %>% group_by(bid,group) %>%
      summarise(logwidth = log(sum(width)),
                hits = sign(sum(hit)))
    return(df)
  })
  
  message("Fit glm model and retrieve odd ratios ...")
  logit.res <- pbmclapply(ann.ov, function(df){
    fits <- glm(df$hits~df$group + df$logwidth,
                family = binomial(link = "logit"),
                data = df)
    res <- logistic.display(fits)
    return(res$table[1,])
  })
  
  all.or.df <- Reduce(rbind,logit.res) %>%
    as.data.frame() %>%
    `rownames<-`(NULL)
  all.or.df$group <- names(logit.res)
  all.or.df$log2OR <- log2(all.or.df$OR)
  all.or.df$log2lower95ci <- log2(all.or.df$lower95ci)
  all.or.df$log2upper95ci <- log2(all.or.df$upper95ci)
  all.or.df$log10P <- log2(all.or.df$`Pr(>|Z|)`)
  
  p <- ggplot(all.or.df,aes(x=group,y=log2OR)) +
    geom_col(aes(fill = -log10P)) +
    geom_errorbar(aes(ymin = log2lower95ci,
                      ymax = log2upper95ci),
                  width = 0.3) +
    theme_bw() + xlab("") +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) +
    scale_fill_viridis_c()
  
  return(list(data = all.or.df,
              plot = p))
}
