# 1. Loading data at the genus taxonomic level: COUNTS VALUES
	# b) Bacteria at the genus taxonomic level (Fig. 3b)
		# Editing and formating the data for proper analyses

	data <- read.delim('/Users/mariohg/Dropbox/Norma/data/counts/bact/bact_rare_L6.txt', comment.char = '#')[,-1]
	rownames(data) <- paste ("Genus", 1:nrow(data), sep = "_")
	taxonomy <- cbind ("TaxaID" = rownames(data), data[,1:5]) #change to 2 for virus family
	head(taxonomy)

	# take the matrix of data only
	data <- data[,-1:-5]; head(data)
	
	# Loading the metadata
	mapping <- read.delim('/Users/mariohg/Dropbox/Norma/data/mapping.txt'); mapping
	
	# 2. DESeq2 Matrix procedure
	countData <- as.matrix(data+1); head (countData)
	dds <- DESeq2::DESeqDataSetFromMatrix(countData = countData, colData = mapping, design = as.formula(~ Treatment))
	res <- DESeq2::DESeq(object = dds, parallel = TRUE)
	
	# 3. Comparison between treatments: something like ratio but with statistics
	# Comparison between treatments: MAS vs UNS
	results <- DESeq2::results(object = res, contrast = c("Treatment", "MAS", "UNS")) #original
	#results <- DESeq2::results(object = res, contrast = c("Treatment", "UNS", "MAS")) #edited
	
	# convert to dataframe
	output <- data.frame(results)
	nrow(output)
	
	# filter those highly significant and log2FoldChange higher than abs (1)
	outfiltered <- subset(output, padj <= 0.05 & abs(log2FoldChange) > 1)
	nrow (outfiltered)
	
	# merge with their corresponding taxonomic levels
	outfiltered$TaxaID <- rownames(outfiltered)
	outfiltered <- base::merge (x = outfiltered, y = taxonomy, by = "TaxaID")
	
	# editing the results for plotting
		# delete the g__ pattern,
		# add a factor for color plotting
		# Ordering by logFoldChange
	library(dplyr)
	outfiltered$Level_6 <- as.character(gsub(pattern = "g__", replacement = "", outfiltered$Level_6))
	outfiltered$Level_2 <- as.character(gsub(pattern = "p__", replacement = "", outfiltered$Level_2))
	#outfiltered <- outfiltered %>% dplyr::mutate (mycolor = ifelse(test = log2FoldChange > 0, yes = "type1", no = "type2"))
	#outfiltered <- outfiltered %>% dplyr::arrange(desc(outfiltered$log2FoldChange))
	
	# 5. Plotting
	colnames(outfiltered)
	#pdf(file = '/Users/mariohg/Renv/Fig.3b-legends.pdf', width = 4, height = 6, paper = "letter")
	ggplot(data = outfiltered, aes(x= Level_6, y=log2FoldChange)) +
		geom_point(aes(color = Level_2), size = 2, shape=16, stroke=2) +
		theme_light() + labs (x = "", y = "Log2FoldChange") +
		#scale_color_manual(values = RColorBrewer::brewer.pal(5, "Accent")) +
		scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")) +
		scale_x_discrete(limits=c(as.character(outfiltered$Level_6))) +
		geom_hline(yintercept = -2, linetype="dashed", color = "gray90", size= 0.75) +
		geom_hline(yintercept = 2, linetype="dashed", color = "gray90", size= 0.75) +
		theme(legend.position = "none", #"left",
		      panel.border = element_rect(colour = "black", size = 1), 
		      axis.text.y = element_text(face = "italic", color ="black", size = 12),
		      axis.text.x = element_text(size = 12, color = "black"),
		      panel.grid = element_blank()) +
		scale_y_continuous(limits = c(-3, 3), breaks = seq (from = -4, to = 3, by = 1)) +
		geom_hline(yintercept = 0, linetype = "solid", color = "gray90", size = 0.5) +
		coord_flip()
	dev.off()
	
	
	