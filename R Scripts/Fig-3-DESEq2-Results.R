## DESeq2 on metagenome:
#######################
rm (list = ls()); dev.off()
library(ggplot2)
set.seed(999)
library(DESeq2); BiocParallel::register(MulticoreParam(4))
library(pheatmap)

# 1. Loading data at the genus taxonomic level: COUNTS VALUES
	# Uncomment the dataset you want to test:
	#a) Bacteria at the phylum taxonomic level (Fig. 3a)
	data <- read.delim('/Users/mariohg/Dropbox/Norma/data/counts/bact/bact_rare_L2.txt', comment.char = '#')[,-1]
	rownames(data) <- data$Level_2; data <- data[,-1]
	head(data)
	
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
	outfiltered$TaxaID <- rownames(outfiltered); outfiltered
	
	
	# editing the results for plotting
		# delete the g__ pattern,
		# add a factor for color plotting
		# Ordering by logFoldChange
	library(dplyr)
	outfiltered$TaxaID <- as.character(gsub(pattern = "p__", replacement = "", outfiltered$TaxaID))
	
	# 5. Plotting
	outfiltered
	pdf(file = '/Users/mariohg/Renv/Fig.3a.pdf', width = 4, height = 2, paper = "letter")
	ggplot(data = outfiltered, aes(x= TaxaID, y=log2FoldChange)) +
		geom_point(size = 2, color= "black", shape=16, stroke=2) +
		theme_light() + labs (x = "", y = "Log2FoldChange") +
		scale_x_discrete(limits=c(as.character(outfiltered$TaxaID))) +
		geom_hline(yintercept = -2, linetype="dashed", color = "gray90", size= 0.75) +
		geom_hline(yintercept = 2, linetype="dashed", color = "gray90", size= 0.75) +
		theme(legend.position = "none",
		      panel.border = element_rect(colour = "black", size = 1), 
		      axis.text.y = element_text(color ="black", size = 12),
		      axis.text.x = element_text(size = 12, color = "black"),
		      panel.grid = element_blank()) +
		scale_y_continuous(limits = c(-3, 3), breaks = seq (from = -4, to = 3, by = 1)) +
		geom_hline(yintercept = 0, linetype = "solid", color = "gray90", size = 0.5) +
		coord_flip()
	dev.off()
