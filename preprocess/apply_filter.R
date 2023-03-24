library(arrow)
library(GenomicRanges)
library(dplyr)

args <- commandArgs(TRUE)
parquet.file = args[1]
parquet.out = args[2]

data.pq <- read_parquet(parquet.file)
print(head(data.pq))

data.pq.flt <- data.pq[data.pq$filter == "pass",]
write_parquet(data.pq.flt, parquet.out)