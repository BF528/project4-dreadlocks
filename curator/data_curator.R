setwd("/projectnb2/bf528/users/dreadlocks/project_4")

library(tidyverse)

# set file path before loading
file_path <- list.files("curator/result/barcode_count", pattern = "^SRR",
                        full.names = F)
meta_data <- read.csv("curator/data/meta_data.txt")

# reading tables
file_list <- vector("list", length = length(file_path))
for(i in seq_along(file_path)){
  
  file_list[[i]] <- read.table(paste0("curator/result/barcode_count/",
                                    file_path[i]),
                               row.names = NULL,
                               col.names = "barcode")
}
count_barcode <- do.call(rbind, file_list)

b_n <- sapply(file_list, nrow)
run_name <- lapply(seq_along(b_n), 
                   function(x) rep(unlist(strsplit(file_path[x], 
                                                   split = "_"))[1],
                                   b_n[x])) %>% unlist()
length(run_name) == nrow(count_barcode)

sample_name <- meta_data[match(run_name,meta_data$Run), "Sample.Name"]
length(sample_name) == nrow(count_barcode)

count_barcode$run_name <- run_name
count_barcode$sample_name <- sample_name

colnames(count_barcode)[1] <- "count"
count_barcode$count <- count_barcode$count %>% as.numeric()

head(count_barcode)
str(count_barcode)

count_barcode <- count_barcode %>% filter(nchar(barcode) == 19) 
count_barcode %>% head()

count_barcode %>% filter(sample_name == "GSM2230758") %>%
  ggplot(aes(log(count))) + 
  stat_ecdf(geom = "step")+
  labs(title="Empirical Cumulative \n Density Function",
       y = "Density", x="log(count) in barcode")+
  facet_wrap(~ run_name) +
  theme_bw()

write.csv(count_barcode, file = "curator/result/barcode_count/count_barcode",
          row.names = F)

