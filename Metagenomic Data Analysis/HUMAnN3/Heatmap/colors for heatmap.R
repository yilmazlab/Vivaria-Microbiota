##Transform all cytokine values so that mean = 0 and SD = 1
heatmap_data_cytokines[4:59] <- lapply(heatmap_data_cytokines[4:59], scale) #heatmap_data_cytokines[4:59] in this code are the columns with values that will be displayed in the heatmap

### then code for preparing the correct matrix format for the heatmap function followed

#set the heat color breaks (you want 1000 color gradient steps starting at the lowest normalised value and ending at the highest normalised value of the dataset)
paletteLength <- 1000 #chose as you wish
myBreaks <- c(seq(min(heatmap_matrix_flipped_num[1:56,]), 0, length.out=ceiling(paletteLength/2) + 1), #heatmap_matrix_flipped_num[1:56,] in this code are the columns with values that will be displayed in the heatmap
              seq(max(heatmap_matrix_flipped_num[1:56,])/paletteLength, max(heatmap_matrix_flipped_num[1:56,]), length.out=floor(paletteLength/2)))


#produce heatmap
pheatmap(#put the data and parameters as you want and use the color settings as follows
         color = colorRampPalette(rev(brewer.pal(n = 11, name =
                                                   "RdBu")))(paletteLength),
         breaks = myBreaks)