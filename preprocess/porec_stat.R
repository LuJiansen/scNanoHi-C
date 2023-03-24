library(dplyr)
library(tibble)

args <- commandArgs(TRUE)
sample = args[1]
enzyme = args[2]
ref = args[3]
phase = args[4]
con_file = args[5] #.concatemer_summary.csv

# concatamer summary
con_df = read.csv(con_file)
con_df <- split(con_df,f=con_df$section)
con_df.v <- lapply(con_df, function(x){
  x$variable <- paste(x$section,x$level_0,x$level_1,x$level_2,sep = "_")
  x <- x[,c("value","variable")]
  y <- x$value
  names(y) <- x$variable
  return(y)
})
con_df.v <- con_df.v[-7] # remove user_metadata

con_df.m <- lapply(con_df.v, as.data.frame) %>%
  Reduce(rbind,.) %>% as.data.frame() %>%
  `colnames<-`(sample)
con_df.m[[sample]] <- as.numeric(con_df.m[[sample]]) %>% round(digits = 2)
print(head(con_df.m))

# haplotag summary
if (phase != "unphased") {
  haptag.list <- list.files(path = "../mapping",
                            pattern = paste0(enzyme,"_",sample,"_batch.*_",ref,"_",
                                             phase,".haplotagged.txt$"),
                            full.names = T)
  names(haptag.list) <- gsub(".haplotagged.txt","",basename(haptag.list))
  haptag <- lapply(haptag.list, function(x) read.table(x)) %>% Reduce(rbind,.)
  colnames(haptag) <- c("align","haplotype","phased_set","chromosome")
  haptag.df <- c(nrow(haptag[haptag$haplotype != "none",])/nrow(haptag)*100,
                 nrow(haptag[haptag$haplotype == "H1",])/nrow(haptag)*100,
                 nrow(haptag[haptag$haplotype == "H2",])/nrow(haptag)*100) %>%
    as.data.frame() %>% `rownames<-`(c("haptag_ratio","haptag_h1_ratio","haptag_h2_ratio")) %>%
    `colnames<-`(sample)
  print(head(haptag.df))
} else {
  haptag.df <- rep(0,3) %>% as.data.frame() %>%
    `rownames<-`(c("haptag_ratio","haptag_h1_ratio","haptag_h2_ratio")) %>%
    `colnames<-`(sample)
}

# alignment summary
library(data.table)
align.list <- list.files(path = "../align_table",
                          pattern = paste0(enzyme,"_",sample,"_batch.*_",ref,"_",
                                           phase,".at.pore_c.parquet.log$"),
                         full.names = T)
names(align.list) <- gsub(".at.pore_c.parquet.log","",basename(align.list))
align.df <- lapply(align.list, function(x){
  fread(x,skip = grep("Summary", readLines(x))) %>% as.data.frame()
}) %>% Reduce(function(x,y) full_join(x,y,by="V1"),.) %>%
  column_to_rownames("V1") %>%
  rowSums() %>% as.data.frame() %>%
  `colnames<-`(sample)
rownames(align.df) <- paste0("align_",rownames(align.df))
print(head(align.df))

sum_df <- list(con_df.m,haptag.df,align.df) %>% Reduce(rbind,.)
sum_df$variable <- rownames(sum_df)

write.table(sum_df,paste(sample,enzyme,ref,phase,"pore_c_stats.txt",sep = "_"),
            row.names = F,quote = F,sep = "\t")