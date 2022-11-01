library(dplyr)
library(tibble)

args <- commandArgs(TRUE)
enzyme <- args[1]
ref <- args[2]
phase <- args[3]

stat1.list <- list.files(pattern = paste0("_",enzyme,"_",ref,"_", phase, "_pore_c_stats.txt"),
  full.names = TRUE)
names(stat1.list) <- gsub(paste0("_",enzyme,"_",ref,"_",phase,"_pore_c_stats.txt"),"",basename(stat1.list))
stat1 <- lapply(stat1.list, function(x) read.table(x,header = T)) %>%
  Reduce(full_join,.) %>% column_to_rownames("variable") %>% 
  t() %>% as.data.frame()
stat1$cell <- rownames(stat1)
print(head(stat1))

stat2.list <- list.files(pattern = paste0("_",enzyme,"_",ref,"_",phase,"_contact_filtration_stats.txt"),full.names = T)
names(stat2.list) <- gsub(paste0("_",enzyme,"_",ref,"_",phase,"_contact_filtration_stats.txt"),"",basename(stat2.list))
stat2 <- lapply(stat2.list, function(x) read.table(x,header = T)) %>%
  Reduce(full_join,.) %>% column_to_rownames("variable") %>% 
  t() %>% as.data.frame()
colnames(stat2) <- paste("filtered",colnames(stat2),sep = "_")
stat2$cell <- rownames(stat2)
print(head(stat2))

stat.df <- full_join(stat1,stat2,by = "cell") %>% column_to_rownames("cell")
stat.df$mapping_rate <- (1 - stat.df$align_unmapped/stat.df$reads__count_all)*100
stat.df$filtered_contacts_per_gb <- stat.df$filtered_pass/stat.df$bases__bases_all*10^9
stat.df$all_high_order_ratio <- 100 - stat.df$concatemer_order__perc_2
stat.df$filtered_adjacent_ratio <- (stat.df$filtered_adjacent/stat.df$contacts_total_count_contacts)*100
stat.df$filtered_close_ratio <- (stat.df$filtered_close/stat.df$contacts_total_count_contacts)*100
stat.df$filtered_duplication_ratio <- (stat.df$filtered_duplication/stat.df$contacts_total_count_contacts)*100
stat.df$filtered_promiscuous_ratio <- (stat.df$filtered_promiscuous/stat.df$contacts_total_count_contacts)*100
stat.df$filtered_isolated_ratio <- (stat.df$filtered_isolated/stat.df$contacts_total_count_contacts)*100
stat.df$filtered_pass_ratio <- (stat.df$filtered_pass/stat.df$contacts_total_count_contacts)*100

stat_select <- stat.df[,c("bases__bases_all","reads__count_all","all_high_order_ratio",
                          "contacts_total_count_contacts","contacts_total_perc_trans_contacts",
                          "density__contacts_per_gb_all","reads__perc_concatemer",
                          "haptag_ratio","mapping_rate",
                          "filtered_trans_ratio","filtered_high_order",
                          "filtered_full_phased_ratio","filtered_semi_phased_ratio",
                          "filtered_median_monomer_length","filtered_median_read_length",
                          "filtered_pass_reads_counts",
                          "filtered_HP1_ratio","filtered_HP2_ratio","filtered_phased_ratio",
                          "filtered_adjacent_ratio","filtered_close_ratio",
                          "filtered_duplication_ratio","filtered_promiscuous_ratio",
                          "filtered_isolated_ratio","filtered_pass_ratio",
                          "filtered_pass","filtered_contacts_per_gb")]
colnames(stat_select) <- c("all_total_base","all_reads_count","all_high_order",
                           "all_contacts_count","all_trans_ratio",
                           "all_contacts_per_gb","all_concatemaer_perc",
                           "all_haptag_ratio","all_mapping_rate",
                           "filtered_trans_ratio","filtered_high_order",
                           "filtered_full_phased_ratio","filtered_semi_phased_ratio",
                           "filtered_median_monomer_length","filtered_median_read_length",
                           "filtered_pass_reads_counts",
                           "filtered_HP1_ratio","filtered_HP2_ratio","filtered_monomer_phased_ratio",
                           "filtered_adjacent_ratio","filtered_close_ratio",
                           "filtered_duplication_ratio","filtered_promiscuous_ratio",
                           "filtered_isolated_ratio","filtered_pass_ratio",
                           "filtered_pass_count","filtered_contacts_per_gb")

stat_sum <- apply(stat_select, 2, function(x) {
  min = min(x)
  q1 = as.numeric(quantile(x,0.25))
  mean = mean(x)
  median = median(x)
  q3 = as.numeric(quantile(x,0.75))
  max = max(x)
  return(c(min = min,
           q1 = q1,
           mean = mean,
           median = median,
           q3 = q3,
           max = max))
}) %>% t() %>%
  as.data.frame()

write.csv(stat.df, paste("merged", enzyme, ref, phase, "stats_all.csv", sep = "_"), quote = F)
write.csv(stat_select, paste("merged", enzyme, ref, phase, "stats_selected.csv", sep = "_"), quote = F)
write.csv(stat_sum,paste("merged", enzyme, ref, phase, "stats_summary.csv", sep = "_"), quote = F)

# high-order summary
ho.list <- list.files(pattern = paste0("_",enzyme,"_",ref,"_",phase,"_order_stats.txt"),
                      full.names = T)
names(ho.list) <- gsub(pattern = paste0("_",enzyme,"_",ref,"_",phase,"_order_stats.txt"),
                       "",basename(ho.list))
ho <- lapply(ho.list, function(x) read.table(x,header = T) %>% dplyr::select(times,pct)) %>%
  Reduce(function(x,y) full_join(x,y,by = "times"),.) %>%
  `colnames<-`(c("order",names(ho.list)))
ho[is.na(ho)] <- 0
write.csv(ho,paste("merged", enzyme, ref, phase, "high_order_summary.csv", sep = "_"),
          quote = F,row.names = F)