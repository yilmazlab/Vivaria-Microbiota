library(pheatmap)
library(RColorBrewer)
setwd("~/Dropbox/Ongoing Analysis/Wild-Type Mice vs Wild Mice Analysis/2022/Humann3/Heatmap")
annotation_table= read.table("pathabundancesubset_anno.txt")
annotation_row= read.table("annotation_row.txt")
Taxa<-read.table("pathabundancesubset.txt")
Taxa[1:233] <- lapply(Taxa[1:233], scale) 

Taxa<-read.table("pathway_coverage_for_main_pathways.txt")
pheatmap(Taxa, 
         fontsize_row = 15, fontsize_col = 5, cluster_cols = TRUE, cluster_rows = TRUE,
         color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu"))),
         cellwidth = 8, cellheight = 20, cutree_cols = 6, cutree_rows = 11, 
         annotation_col=annotation_table, annotation_colors = ann_colors, border_color = "grey60")



ann_colors = list(
  Gender = c(Female="#cc6666", Male="#B5BBE3", Unknown="#f4fcfc"),
  Content = c(Cecum="#e8a188", Feces="brown"),
  Group = c(WildMice="#38af60", SPF="#3399FF", Human="#efb94c"),
  Classes =c(Biosynthesis="green", Generation_of_Precursor_Metabolites_and_Energy="#41b6c4" , Degradation_Utilization_Assimilation="orange"),
  Child_Classes=c(Amino_Acid_Biosynthesis="#abf5be", Amino_Acid_Degradation="#8c705b", Aromatic_Compound_Biosynthesis="#63c97d", Carbohydrate_Biosynthesis="#24853c", Cell_Structure_Biosynthesis="#023b10", Cofactor_Carrier_and_Vitamin_Biosynthesis="blue", Fatty_Acid_and_Lipid_Biosynthesis="yellow", Fermentation="orange", Glycolysis="#aaf2ff",Nucleoside_and_Nucleotide_Biosynthesis="red", Nucleoside_and_Nucleotide_Degradation="purple",Pentose_Phosphate_Pathways="#023b10", Secondary_Metabolite_Biosynthesis="#E6AFB9"))

rwbcols <- c("#D33F6A","#E07B91","#E6AFB9","#f4fcfc","#B5BBE3", "#8595E1", "#4A6FE3")

paletteLength <- 1000 #chose as you wish
myBreaks <- c(seq(min(Taxa[1:83,]), 0, length.out=ceiling(paletteLength/2) + 1), #Taxa[22:83,] in this code are the columns with values that will be displayed in the heatmap
              seq(max(Taxa[1:83,])/paletteLength, max(Taxa[1:83,]), length.out=floor(paletteLength/2)))



pheatmap(Taxa, 
         fontsize_row = 15, fontsize_col = 5, cluster_cols = TRUE, cluster_rows = TRUE,
         color = colorRampPalette(rev(brewer.pal(n = 11, name =
                                                   "RdBu")))(paletteLength),
         breaks = myBreaks,
         cellwidth = 8, cellheight = 20, cutree_cols = 6, cutree_rows = 11, 
         annotation_col=annotation_table, annotation_row=annotation_row,annotation_colors = ann_colors, border_color = "grey60")


# Pathway coverage
Taxa<-read.table("pathway_coverage_for_main_pathways.txt")

ann_colors = list(
  Gender = c(Female="#cc6666", Male="#B5BBE3", Unknown="#f4fcfc"),
  Content = c(Cecum="#e8a188", Feces="brown"),
  Group = c(WildMice="#38af60", SPF="#3399FF", Human="#efb94c"))

pheatmap(Taxa, 
         fontsize_row = 5, fontsize_col = 10, cluster_cols = TRUE, cluster_rows = TRUE,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="Spectral")))(1000), # Generates a palette of 100 colors,
         cellwidth = 15, cellheight = 5,  cutree_cols = 6, cutree_rows = 10, 
         annotation_col=annotation_table, annotation_colors = ann_colors, border_color = "grey60")

# Diversity to Humann3 pipeline -------------------------------------------
source("metaphlanToPhyloseq.R")
setwd("/Users/bahti/Dropbox/Ongoing Analysis/Wild-Type Mice vs Wild Mice Analysis/2022/Humann3/Diversity to Human")
data = read.delim("pathabundance_merged_ID.txt",sep = "")
data = read.delim("pathcoverage_merged-cpm_2.txt",sep = "")
metadata = read.delim("pathabundancesubset_anno.txt",sep = "")

data[, 2:ncol(data)] <- lapply(data[, 2:ncol(data)], as.numeric)
str(data)
My_Theme = theme(
  axis.title.x = element_text(size = 18),
  axis.text.x = element_text(size = 18),
  axis.title.y = element_text(size = 18),
  axis.text.y = element_text(size = 18))

Humann3 = metaphlanToPhyloseq(data, metadata, simplenames = TRUE, 
                           roundtointeger = FALSE)
rank_names(Humann3)

Humann3.int = Humann3
otu_table(Humann3.int) = round(otu_table(Humann3)*1e5)
alpha_meas = c("Shannon", "Simpson")
p <- plot_richness(Humann3.int, "Group",  measures=alpha_meas, color="Group") + geom_boxplot(alpha = 0.8) + theme_bw() 
p 
ggsave(p, file="alpha_diversity_for all groups.pdf", width=16.69, height=8.27,useDingbats=FALSE)

alpha_diversity_Humann3 <-estimate_richness(Humann3.int, split = TRUE, measures = NULL)
data_cbind <- cbind(sample_data(Humann3.int), alpha_diversity_Humann3)
write.csv(data_cbind, "alphadiversity_Humann3_ALLunintegratedareincluded.csv")

ord= ordinate(Humann3, method="PCoA", distance="jaccard")
plot<-plot_ordination(Humann3, ord, color="Group", shape="Content") + geom_point(size=5, alpha=1) + theme_bw()  + My_Theme  + stat_ellipse(aes(group=Group)) + geom_text_repel(size=5,  aes(label = Location), max.overlaps = Inf) # + stat_ellipse(aes(group=Patient))
plot  
ggsave(plot, file="Beta Diversity of Humann3 for Metagenomic Data.pdf", width=14, height=9,useDingbats=FALSE)

