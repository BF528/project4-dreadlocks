setwd("/projectnb2/bf528/users/dreadlocks/project_4")

library(tidyverse)

count_barcode <- read.csv("curator/result/barcode_count/count_barcode")
mean_barcode <- count_barcode %>% group_by(barcode, sample_name) %>% 
  summarise(mean_count = mean(count))

cut_off <- 10000
GSM2230757_wl <- mean_barcode %>% 
  filter(mean_count > cut_off & sample_name == "GSM2230757") %>% 
  ungroup() %>% select(barcode)

GSM2230758_wl <- mean_barcode %>% 
  filter(mean_count > cut_off & sample_name == "GSM2230758") %>% 
  ungroup() %>% select(barcode)

uf_count <- mean_barcode %>% 
  filter(sample_name == "GSM2230758") %>% 
  ungroup() %>% select(mean_count)

f_count <- mean_barcode %>% 
  filter(mean_count > cut_off & sample_name == "GSM2230758") %>% 
  ungroup() %>% select(mean_count)

par(mfrow = c(1,2))
hist(log(uf_count$mean_count), probability = F,
     main = "", xlab = "log(count) in barcode")
title("A", adj = 0)
hist(log(f_count$mean_count), probability = F,
     main = "", xlab = "log(count) in barcode")
title("B", adj = 0)

GSM2230759_wl <- mean_barcode %>% 
  filter(mean_count > cut_off & sample_name == "GSM2230759") %>% 
  ungroup() %>% select(barcode)

GSM2230760_wl <- mean_barcode %>% 
  filter(mean_count > cut_off & sample_name == "GSM2230760") %>% 
  ungroup() %>% select(barcode)

write.table(GSM2230757_wl, file = "curator/result/barcode_count/GSM2230757_wl.csv", 
            sep=",",  col.names=FALSE,
            row.names = FALSE, quote = F)
write.table(GSM2230758_wl, file = "curator/result/barcode_count/GSM2230758_wl.csv", 
            sep=",",  col.names=FALSE,
            row.names = FALSE, quote = F)
write.table(GSM2230759_wl, file = "curator/result/barcode_count/GSM2230759_wl.csv", 
            sep=",",  col.names=FALSE,
            row.names = FALSE, quote = F)
write.table(GSM2230760_wl, file = "curator/result/barcode_count/GSM2230760_wl.csv", 
            sep=",",  col.names=FALSE,
            row.names = FALSE, quote = F)
