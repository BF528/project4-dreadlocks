setwd("/projectnb2/bf528/users/dreadlocks/project_4/biologist")

dat <- read.csv('/projectnb/bf528/users/dreadlocks/project_4/analyst/sample_markers.csv')
outdir <- '/projectnb/bf528/users/dreadlocks/project_4/biologist/output/'
library(enrichR)
library(tidyverse)


summary(dat)
dat$cluster <- factor(dat$cluster)

# select databases
dbs <- listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", 
         "GO_Biological_Process_2018", "KEGG_2019_Human")


# split data into clusters 
clusters <- split(dat$gene, f=dat$cluster)
# run enrichr for each cluster
enriched <- lapply(clusters, enrichr, databases=dbs)

# copied from follow_up.R, sequence is the same for clusters 
cell_types <- c("alpha_1", "alpha_2", "beta", 
                     "ductal", "acinar", "delta",
                     "endothelial_1", "endothelial_2")

# for each cluster, for each annotation db
for (l in 1:length(enriched)){
  for (i in 1:4){
    enriched[[l]][[i]] <- enriched[[l]][[i]][,c(1,4,8)] # subset columns 
    names(enriched[[l]][[i]]) <- c('Term', 'P_adj', 'Combined.Score') # change column names 
    enriched[[l]][[i]] <- cbind(enriched[[l]][[i]], AnnotationDB = dbs[i]) # add additional column for names of annotation databases  
    enriched[[l]][[i]] <- enriched[[l]][[i]] %>% filter(P_adj < 0.1) %>% # apply P-adjusted cutoff of 0.1
      slice_max(Combined.Score, n = 5, with_ties = FALSE) # select top 10 ranked by combined score 
  }
}


# combine 4 dbs into one dataframe for each cluster
for (l in 1:length(enriched)){ # for each cluster, combine 4 df
  enriched[[l]] <- bind_rows(enriched[[l]][[1]], enriched[[l]][[2]], enriched[[l]][[3]], enriched[[l]][[4]])
  # format digits
  enriched[[l]]['P_adj'] <- format(enriched[[l]]['P_adj'], scientific = TRUE, digits = 3)
  enriched[[l]]['Combined.Score'] <- format(round(enriched[[l]]['Combined.Score'], 1), nsmall = 1)
}


# write output- 1 df for each cluster 
for (i in 1:length(enriched)) {
  filename <- paste0(outdir, "cluster", i-1, cell_types[i], ".csv", sep = "")
  write.csv(x = enriched[[i]],file = filename, row.names = T)
}





