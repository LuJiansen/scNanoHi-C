library(arrow)
library(stringr)

args <- commandArgs(TRUE)
parquet.file = args[1]
regex_str = args[2]
parquet.out = args[3]

data.pq <- read_parquet(parquet.file)
print(head(data.pq))

bool_value = str_detect(data.pq$align1_chrom,regex_str) & str_detect(data.pq$align2_chrom,regex_str)
data.pq.flt = data.pq[bool_value,]

write_parquet(data.pq.flt, parquet.out)