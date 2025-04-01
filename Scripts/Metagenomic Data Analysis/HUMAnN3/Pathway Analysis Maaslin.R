setwd("~/Dropbox/Ongoing Analysis/StomaMerged_November_2018/2021/Humann3 Output/PathwayAnalysis")

library(Maaslin2)
library(pheatmap)
library(viridis)

input_metadata = read.table("metadata_pathways.txt", header=TRUE,row.names="sample")
input_metadata = read.table("metadata_pathways_0_2h.txt", header=TRUE,row.names="sample")
input_metadata = read.table("metadata_pathways_0_6h.txt", header=TRUE,row.names="sample")
input_metadata = read.table("metadata_pathways_0_8h.txt", header=TRUE,row.names="sample")
input_data = read.table("pathways.txt", header=TRUE,row.names="sample")



# Analyzing the statistical differences between feeding and non-feeding time for metabolic functions


fit_data <- Maaslin2(
  input_data, input_metadata, 'Pathway_Oh_6h_mp0.3_ms_0.25_ma0.00000001', transform = "LOG", min_prevalence=0.3, min_abundance=0.00000001, max_significance=0.25, analysis_method="LM",normalization="NONE",
  fixed_effects = c('TimePoint', 'FeedingTime'),reference = c("TimePoint,Zero"),
  random_effects = c('Patient'),
  standardize = FALSE)

fit_data <- Maaslin2(
  input_data, input_metadata, 'Pathway_Oh_48h_mp0.3_ms_0.2', transform = "LOG", min_prevalence=0.3, max_significance=0.2, analysis_method="LM",normalization="NONE",
  fixed_effects = c('TimePoint', 'FeedingTime'),reference = c("TimePoint,Zero"),
  random_effects = c('Patient'),
  standardize = FALSE)



# Heatmap
annotation_table= read.table("annotation_heatmap.txt")
pathways<-read.table("feedingmetabolicfunctions_0h_6h_average10e-5.txt", row.names=NULL)
row.names(pathways) <- make.unique(as.character(pathways$row.names))
pathways$row.names <- NULL

ann_colors = list(
  TimePoint = c(Zero ="#fff5f0", Two="#fcbba1", Four="#fc9272", Six="#ef3b2c", Eight="#cb181d", Ten = "#a50f15", Fortyeight="#67000d"),
  Patient = c(P1="#B5BBE3", P2="#aaf2ff", P3="#3399FF",P4="#efb94c", P6="#cc6666", P7="#f4fcfc"),
  FeedingTime =c(Yes="red", No ="green"))


rwbcols <- c( "#781d1a","#E07B91", "#4A6FE3", "white")

 
feedingmetabolites <- pheatmap(log10(pathways+1),
                               fontsize_row = 10, fontsize_col = 12, 
                               cluster_cols = FALSE, cluster_rows = FALSE,
                               color = turbo(10000),
                               cellwidth = 40, cellheight = 12,# cutree_rows = 4,  
                               annotation_col=annotation_table, annotation_colors = ann_colors, border_color =  "grey60")
feedingmetabolites
ggsave(feedingmetabolites, file="~/Dropbox/Ongoing Analysis/StomaMerged_November_2018/2021/Humann3 Output/PathwayAnalysis/feedingmetabolicfunctions_0h_6h_average10e-5.pdf", width=30, height=30, useDingbats=FALSE,limitsize = FALSE)



