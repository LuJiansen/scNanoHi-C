library(arrow)
library(stringr)
library(dplyr)

args <- commandArgs(TRUE)
parquet.file = args[1]
regex_str = args[2]
parquet.out = args[3]

data.pq <- read_parquet(parquet.file)
print(head(data.pq))

bool_value = str_detect(data.pq$chrom,regex_str)
data.pq.flt = data.pq[bool_value,] %>% droplevels()

write_parquet(data.pq.flt, parquet.out)