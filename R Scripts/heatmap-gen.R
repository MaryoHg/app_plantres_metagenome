### heatmap drawing for each database:
rm (list = ls()); dev.off()
library(RColorBrewer); library(pheatmap)
library(plyr); library(dplyr)
set.seed(999)
options(scipen = 10000, digits = 2)


	########### genus (family for virus) taxonomic level ##

	# Loading the realtive abundance of each database
	# Un-comment the properly-to-use database
		#data <- read.delim('/Users/mariohg/Dropbox/Norma/data/relab/bact/bact_rare_L6.txt', skip = 1)[,-1]
		data <- read.delim('/Users/mariohg/Dropbox/Norma/data/relab/arch/arch_rare_L6.txt', skip = 1)[,-1]
		#data <- read.delim('/Users/mariohg/Dropbox/Norma/data/relab/fungi/fungi_rare_L6.txt', skip = 1)[,-1]
		#data <- read.delim('/Users/mariohg/Dropbox/Norma/data/relab/protozoa/protozoa_rare_L6.txt', skip = 1)[,-1]
		#data <- read.delim('/Users/mariohg/Dropbox/Norma/data/relab/viral/viral_rare_L6.txt', skip = 1)[,-1]


		# format data for later simple use: properly change for viruses and protozoa and fungi (less than 50 genus)
		data <- data[!(data$Level_6 %in% c("g__", "Other")),] #ignore unassgined and other
		data[,6:9] <- data[,6:9] * 100
		data <- data %>% dplyr::arrange(desc(apply(X =data[,6:9], MARGIN = 1, FUN = mean)))
		data <- data[1:50,] # select the 50 most abundant; ignore this for FUNGI AND PROTOZOA
		
		# average and sd
		data$UNS <- apply(X = data[,6:7], MARGIN = 1, FUN = mean)
		data$MAS <- apply(data[,8:9], MARGIN = 1, FUN = mean)
		data$UNSsd <- apply(X = data[,6:7], MARGIN = 1, FUN = sd) 
		data$MASsd <- apply(X = data[,8:9], MARGIN = 1, FUN = sd)
		
		# add factor column for later heatmap for UNS and MAS treatment:
		data <- within(data, {
			UNSf <- 0
			UNSf [UNS <=0.01]<- 0
			UNSf [UNS > 0.01 & UNS <= 0.05]<- 1
			UNSf [UNS > 0.05 & UNS <= 0.10]<- 2
			UNSf [UNS > 0.10 & UNS <= 0.20]<- 3
			UNSf [UNS > 0.20 & UNS  <= 0.50]<- 4
			UNSf [UNS > 0.50 & UNS <= 1.00]<- 5
			UNSf [UNS > 1.00 & UNS <= 2.00]<- 6
			UNSf [UNS > 2.00 & UNS <= 4.00]<- 7
			UNSf [UNS > 4.00 & UNS <= 8.00]<- 8
			UNSf [UNS > 8.00 & UNS <= 12.00]<- 9
			UNSf [UNS > 12.00 & UNS <= 16.00]<- 10
			UNSf [UNS > 16.00 & UNS <= 20.00]<- 11
			UNSf [UNS > 20.00]<- 12})
		
		data <- within(data, {
			MASf <- 0
			MASf [MAS <=0.01]<- 0
			MASf [MAS > 0.01 & MAS <= 0.05]<- 1
			MASf [MAS > 0.05 & MAS <= 0.10]<- 2
			MASf [MAS > 0.10 & MAS <= 0.20]<- 3
			MASf [MAS > 0.20 & MAS  <= 0.50]<- 4
			MASf [MAS > 0.50 & MAS <= 1.00]<- 5
			MASf [MAS > 1.00 & MAS <= 2.00]<- 6
			MASf [MAS > 2.00 & MAS <= 4.00]<- 7
			MASf [MAS > 4.00 & MAS <= 8.00]<- 8
			MASf [MAS > 8.00 & MAS <= 12.00]<- 9
			MASf [MAS > 12.00 & MAS <= 16.00]<- 10
			MASf [MAS > 16.00 & MAS <= 20.00]<- 11
			MASf [MAS > 20.00]<- 12})
		
		# before heatmap: editing names to delete "patterns"
		# CHANGE properly for VIRUSES (change to "Level_4": order)
		data$Level_2 <- as.character(gsub(pattern = "p__", replacement = "", data$Level_2))
		data$Level_6 <- as.character(gsub(pattern = "g__", replacement = "", data$Level_6))
		
		# metadata for colors per Phyla and Treatments
		samp.names <- data.frame ('Sample'=c('UNSf', 'MASf'),
					  'Treatment'=c('UNS', 'MAS'), row.names = 1); samp.names
		
		# CHANGE properly only for VIRUS database (Order & Level_4)
		row.data <- data.frame("Phylum" = data$Level_2, data$Level_6, row.names = 2); head(row.data)
		
		# color for heatmapping:
		# Phyla first, then the list complete
		Phylum <- RColorBrewer::brewer.pal(n = length(levels(factor(data$Level_2))), name = "Set3"); Phylum
		names(Phylum) <- levels(factor(data$Level_2)); Phylum
		
		heat_colors <- list (Treatment = c('UNS' = 'blue', 'MAS' = 'red'),
				     Phylum = Phylum); heat_colors
		
		# Ordering from low to high to plot better
		data_mat <- data %>% dplyr::arrange(desc(apply(X =data[,6:9], MARGIN = 1, FUN = mean)))
		rownames(data_mat) <- data$Level_6
		
		# italic rownames and breaks
		genusnames <- lapply(rownames (data_mat), function(x) bquote(italic(.(x))))
		breaks <- base::seq(from = 0, to = 12, by = 0.01)
		
		# Drawing the heatmap
		#pdf(width = 8.5, height = 11, paper = "letter", )
		pheatmap::pheatmap(mat = data_mat[,14:15],
				   # colors
			 annotation = samp.names,
			 annotation_colors = heat_colors,
			 annotation_row = row.data,
			 #settings
			 #clustering_method = "complete", #clustering_distance_rows = "euclidean",
			 fontsize = 8, show_colnames = F, labels_row = as.expression(genusnames),
			 cellwidth = 20, cellheight = 10, treeheight_row = 15,
			 cluster_rows = F, cluster_cols = T, cutree_cols = 2, treeheight_col = 15,
			 color = colorRampPalette(c('white', 'yellow', 'red'))(length(breaks)),
			 breaks = breaks, legend_breaks = 0:12, legend_labels = c('<=0.01',
			 				   '> 0.01 & <= 0.05',
			 				   '> 0.05 & <= 0.10',
			 				   '> 0.10 & <= 0.20',
			 				   '> 0.20 &  <= 0.50',
			 				   '> 0.50 & <= 1.00',
			 				   '> 1.00 & <= 2.00',
			 				   '> 2.00 & <= 4.00',
			 				   '> 4.00 & <= 8.00',
			 				   '> 8.00 & <= 12.00',
			 				   '> 12.00 & <= 16.00',
			 				   '> 16.00 & <= 20.00',
			 				   '> 20.00'))
		