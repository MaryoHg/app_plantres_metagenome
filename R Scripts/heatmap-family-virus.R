					
	# loading the family virus level
	data <- read.delim('/Users/mariohg/Dropbox/Norma/data/relab/viral/viral_rare_L5.txt', skip = 1)[,-1]
	
	
	# format data for later simple use: properly change for viruses and protozoa and fungi (less than 50 genus)
	data <- data[!(data$Level_5 %in% c("f__", "Other")),] #ignore unassgined and other
	data[,5:8] <- data[,5:8] * 100
	data <- data %>% dplyr::arrange(desc(apply(X =data[,5:8], MARGIN = 1, FUN = mean)))
	data <- data[1:50,] # select the 50 most abundant; ignore this for FUNGI AND PROTOZOA
	# average and sd
	data$UNS <- rowMeans(data[,5:6]);
	data$MAS <- rowMeans(data[,7:8])
	data$UNSsd <- apply(data[,5:6], 1, sd); 
	data$MASsd <- apply(data[,7:8],1, sd)
	# add factor column for later heatmap per treatment: UNS and MAS
	data <- within(data, {
		UNSf <- 0
		UNSf [UNS <=  0.01]   <- 0
		UNSf [UNS >  0.01 & UNS <= 0.05]<- 1
		UNSf [UNS >  0.05 & UNS <= 0.10]<- 2
		UNSf [UNS >  0.10 & UNS <= 0.20]<- 3
		UNSf [UNS >  0.2 & UNS  <=  0.50]<- 4
		UNSf [UNS >  0.50 & UNS <=  1.00]<- 5
		UNSf [UNS >  1.00 & UNS <= 2.00]<- 6
		UNSf [UNS > 2.00 & UNS <= 5.00]<- 7
		UNSf [UNS > 5.00 & UNS <= 10.00]<- 8
		UNSf [UNS > 10.00 & UNS <= 25.00]<- 9
		UNSf [UNS > 25.00 & UNS <= 50.00]<- 10
		UNSf [UNS > 50.00 & UNS <= 75.00]<- 11
		UNSf [UNS > 75.00] <- 12})
	
	data <- within(data, {
		MASf <- 0
		MASf [MAS <=  0.01]   <- 0
		MASf [MAS >  0.01 & MAS <= 0.05]<- 1
		MASf [MAS >  0.05 & MAS <= 0.10]<- 2
		MASf [MAS >  0.10 & MAS <= 0.20]<- 3
		MASf [MAS >  0.2 & MAS  <=  0.50]<- 4
		MASf [MAS >  0.50 & MAS <=  1.00]<- 5
		MASf [MAS >  1.00 & MAS <= 2.00]<- 6
		MASf [MAS > 2.00 & MAS <= 5.00]<- 7
		MASf [MAS > 5.00 & MAS <= 10.00]<- 8
		MASf [MAS > 10.00 & MAS <= 25.00]<- 9
		MASf [MAS > 25.00 & MAS <= 50.00]<- 10
		MASf [MAS > 50.00 & MAS <= 75.00]<- 11
		MASf [MAS > 75.00] <- 12})
	# before heatmap: editing names to delete "patterns"
	# CHANGE properly for VIRUSES (change to "Level_4": order)
	data$Level_4 <- as.character(gsub(pattern = "o__", replacement = "", data$Level_4))
	data$Level_5 <- as.character(gsub(pattern = "f__", replacement = "", data$Level_5))
	
	# metadata for colors per Phyla and Treatments
	samp.names <- data.frame ('Sample'=c('UNSf', 'MASf'),
				  'Treatment'=c('UNS', 'MAS'), row.names = 1); samp.names
	# CHANGE propoerly only for VIRUS database (Order & Level_4)
	row.data <- data.frame("Order" = data$Level_4, data$Level_5, row.names = 2); head(row.data)
	
	# color for heatmapping:
	# Phyla first, then the list complete
	Phylum <- RColorBrewer::brewer.pal(n = length(levels(factor(data$Level_4))), name = "Set3"); Phylum
	names(Phylum) <- levels(factor(data$Level_4)); Phylum
	
	heat_colors <- list (Treatment = c('UNS' = 'blue', 'MAS' = 'red'),
			     Order = Phylum); heat_colors
	
	# Ordering from low to high to plot better
	data_mat <- data %>% dplyr::arrange(desc(apply(X =data[,5:8], 1, mean)))
	rownames(data_mat) <- data$Level_5
	
	# italic rownames and breaks
	genusnames <- lapply(rownames (data_mat), function(x) bquote(italic(.(x))))
	breaks <- base::seq(from = 0, to = 12, by = 0.01)
	
	# Drawing the heatmap
	#pdf(width = 8.5, height = 11, paper = "letter", )
	pheatmap(mat = data_mat[,13:14],
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
		 breaks = breaks, legend_breaks = 0:12, legend_labels = c('<= 0.01',
		 							 '> 0.01 & <= 0.10',
		 							 '> 0.05 & <= 0.10',
		 							 '> 0.10 & <= 0.20',
		 							 '> 0.20 & <= 0.50',
		 							 '> 0.50 & <= 1.00',
		 							 '> 1.00 & <= 2.00',
		 							 '> 2.00 & <= 5.00',
		 							 '> 5.00 & <= 10.0',
		 							 '> 10.0 & <= 25.0',
		 							 '> 25.0 & <= 50.0',
		 							 '> 50.0 & <= 75.0',
		 							 '> 75.0'))
	