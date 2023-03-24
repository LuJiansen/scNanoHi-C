library(arrow)
library(GenomicRanges)
library(dplyr)

args <- commandArgs(TRUE)
parquet.file = args[1]
parquet.out = args[2]

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

if(file_test('-f',parquet.file)){
  data.pq <- read_parquet(parquet.file)
} else if (file_test('-d',parquet.file)) {
  data.pq <- read_parquet_dir(parquet.file)
} else{
  stop("Wrong input!")
}

print(head(data.pq))

order.stat <- table(data.pq$read_name)

data.pq.high <- data.pq[data.pq$read_name %in% names(order.stat[order.stat>1]),]
sum(table(data.pq.high$read_name) == 1)

write_parquet(data.pq.high, parquet.out)