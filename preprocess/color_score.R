library(dplyr)
library(reshape2)
library(tidyverse)

args <- commandArgs(TRUE)
enzyme = args[1]
indir = args[2]
file_pattern1 = args[3]
file_pattern2 = args[4]
output = args[5]

color.files <- list.files(path = indir,
                        pattern = paste0(file_pattern1,
                                         ".*",file_pattern2),
                        full.names = T)
names(color.files) <- gsub(paste0(enzyme,"_"),"",basename(color.files)) %>%
  gsub(file_pattern1,"",.) %>% gsub(file_pattern2,"",.)

color <- lapply(color.files, function(x){
  df <- read.table(x,header = F)
  name <- gsub(paste0(enzyme,"_"),"",basename(x)) %>%
    gsub(file_pattern1,"",.) %>% 
    gsub(file_pattern2,"",.)
  colnames(df) <- c("chr","pos",name)
  df$id <- paste(df$chr,df$pos,sep = "_")
  df <- df %>% select(-c(chr,pos))
  return(df)
}) %>% Reduce(full_join,.) %>%
  column_to_rownames("id")

write.csv(color,output)