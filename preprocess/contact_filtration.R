library(arrow)
library(GenomicRanges)
library(dplyr)

args <- commandArgs(TRUE)
sample_name = args[1]
parquet.file = args[2]
parquet.out = args[3]

read_parquet_dir <- function(path){
  ## Define the dataset
  DS <- arrow::open_dataset(sources = path)
  ## Create a scanner
  SO <- Scanner$create(DS)
  ## Load it as n Arrow Table in memory
  AT <- SO$ToTable()
  ## Convert it to an R data frame
  DF <- as.data.frame(AT)
  remove(DS)
  remove(SO)
  remove(AT)
  return(DF)
}

contact_filter <- function(cont.df){
  require(GenomicRanges)
  cont.df$idx <- seq(1:nrow(cont.df))
  cont.df$filter <- "pass"
  cont.df.flt1 <- cont.df[!cont.df$contact_fragment_adjacent,]
  cont.df[!(cont.df$idx %in% cont.df.flt1$idx),]$filter <- "adjacent"

  cont.df.flt2 <- cont.df.flt1[!(cont.df.flt1$contact_is_cis &
                                   abs(cont.df.flt1$contact_genome_distance) < 1000),]
  cont.df[cont.df$filter == "pass" & !(cont.df$idx %in% cont.df.flt2$idx),]$filter <- "close"

  read_stat <- table(cont.df.flt2$read_name) %>%
    as.data.frame() %>% `colnames<-`(c("read_name","n_contacts"))

  cont.df.flt2 <- full_join(cont.df.flt2,read_stat)
  cont.df.flt2 <- cont.df.flt2 %>% arrange(desc(n_contacts))
  cont.df.flt3 <- cont.df.flt2[!duplicated(cont.df.flt2[,c("align1_fragment_id","align2_fragment_id")]
),]
  cont.df[cont.df$filter == "pass" & !(cont.df$idx %in% cont.df.flt3$idx),]$filter <- "duplication"

  frag_df <- table(c(cont.df.flt3$align2_fragment_id,
                     cont.df.flt3$align1_fragment_id)) %>%
    as.data.frame() %>% arrange(desc(Freq))

  prom_fg <- frag_df[frag_df$Freq > 10,]$Var1
  cont.df.flt4 <- cont.df.flt3[!(cont.df.flt3$align1_fragment_id %in% prom_fg |
                                   cont.df.flt3$align2_fragment_id %in% prom_fg),]
  cont.df[cont.df$filter == "pass" & !(cont.df$idx %in% cont.df.flt4$idx),]$filter <- "promiscuous"

  frag.gr.df <- rbind(cont.df.flt4[,c("align1_fragment_id","align1_chrom",
                                   "align1_fragment_start","align1_fragment_end")] %>%
                     `colnames<-`(c("id","chr","start","end")),
                   cont.df.flt4[,c("align2_fragment_id","align2_chrom",
                                   "align2_fragment_start","align2_fragment_end")] %>%
                     `colnames<-`(c("id","chr","start","end"))) %>% unique()
  #frag.gr.df <- frag.gr.df %>% arrange(chr,start)
  frag.gr <- makeGRangesFromDataFrame(frag.gr.df,keep.extra.columns = T)
  dis <- distanceToNearest(frag.gr)
  iso_frag = frag.gr$id[which(dis@elementMetadata$distance > 10000000)]

  cont.df.flt5 <- cont.df.flt4[!(cont.df.flt4$align1_fragment_id %in% iso_frag |
                                   cont.df.flt4$align2_fragment_id %in% iso_frag),]
  cont.df[cont.df$filter == "pass" & !(cont.df$idx %in% cont.df.flt5$idx),]$filter <- "isolated"
  remove(cont.df.flt1)
  remove(cont.df.flt2)
  remove(cont.df.flt3)
  remove(cont.df.flt4)
  remove(cont.df.flt5)
  return(cont.df)
}

data.pq <- read_parquet_dir(parquet.file)
data.pq.flt <- contact_filter(data.pq)

write_parquet(data.pq.flt, parquet.out)