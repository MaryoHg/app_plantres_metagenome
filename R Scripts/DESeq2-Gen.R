## DESeq2 on metagenome:
#######################
rm (list = ls()); dev.off()
library(ggplot2)
set.seed(999)
library(DESeq2); BiocParallel::register(MulticoreParam(4))
library(pheatmap)

# 1. Loading data at the genus taxonomic level: COUNTS VALUES
	# LUncomment the dataset you want to test:
	data <- read.delim('/Users/mariohg/Dropbox/Norma/data/counts/bact/bact_rare_L6.txt', comment.char = '#')[,-1]
	#data <- read.delim('/Users/mariohg/Dropbox/Norma/data/counts/arch/arch_rare_L6.txt', comment.char = '#')[,-1]
	#data <- read.delim('/Users/mariohg/Dropbox/Norma/data/counts/fungi/fungi_rare_L6.txt', comment.char = '#')[,-1]
	#data <- read.delim('/Users/mariohg/Dropbox/Norma/data/counts/protozoa/protozoa_rare_L6.txt', comment.char = '#')[,-1]
	#data <- read.delim('/Users/mariohg/Dropbox/Norma/data/counts/virus/viral_rare_L6.txt', comment.char = '#')[,-1]
	#data <- read.delim('/Users/mariohg/Dropbox/Norma/data/counts/virus/viral_rare_L5.txt', comment.char = '#')[,c(-1:-3)]

	
	
	# Editing and formating the data for proper analyses
	rownames(data) <- paste ("Genus", 1:nrow(data), sep = "_")
	taxonomy <- cbind ("GenusID" = rownames(data), data[,1:5]) #change to 2 for virus family
	head(taxonomy)
	
	# take the matrix of data only
	data <- data[,-1:-5]; head(data) # chabge to 2 for virus families

	# Loading the metadata
	mapping <- read.delim('/Users/mariohg/Dropbox/Norma/data/mapping.txt'); mapping
	
# 2. DESeq2 Matrix procedure
	countData <- as.matrix(data+1); head (countData)
	dds <- DESeq2::DESeqDataSetFromMatrix(countData = countData, colData = mapping, design = as.formula(~ Treatment))
	res <- DESeq2::DESeq(object = dds, parallel = TRUE)
	
# 3. Comparison between treatments: something like ratio but with statistics
	# Comparison between treatments: MAS vs UNS
	results <- DESeq2::results(object = res, contrast = c("Treatment", "MAS", "UNS"))
	
	# convert to dataframe
	output <- data.frame(results); nrow(output)
	
	# filter those highly significant and log2FoldChange higher than abs (1)
	outfiltered <- subset(output, padj <= 0.05 & abs(log2FoldChange) > 1)
	nrow (outfiltered)
	
	# merge with their corresponding taxonomic levels
	outfiltered$GenusID <- rownames(outfiltered)
	outfiltered <- base::merge (x = outfiltered, y = taxonomy, by = "GenusID")
	
	# editting the results for plotting
	# delete the g__ patter,
	# add a factor for color plotting
	# Ordering by logFoldChange
	outfiltered$Level_6 <- as.character(gsub(pattern = "g__", replacement = "", outfiltered$Level_6))
	outfiltered <- outfiltered %>% dplyr::mutate (mycolor = ifelse(test = log2FoldChange > 0, yes = "type1", no = "type2"))
	outfiltered <- outfiltered %>% dplyr::arrange(desc(outfiltered$log2FoldChange))
	
	#outfiltered$Genera <- paste (outfiltered$Level_6, paste("(P = ", round(outfiltered$padj, digits = 3), ")", sep = ""), sep = " ")
	
# 5. Plotting
	colnames(outfiltered)
	#virus_gen <- 
	ggplot(data = outfiltered, aes(x= Level_6, y=log2FoldChange)) +
		geom_point(size = 2, color= "black", shape=16, stroke=2) +
		theme_light() + labs (x = "Virus: Genus level", y = "Log2FoldChange") +
		scale_color_manual(values = c("black", "gray50")) +
		scale_x_discrete(limits=c(as.character(outfiltered$Level_6))) +
		theme(legend.position = "none",
		      #panel.border = element_blank(), axis.text.x = element_text(size = 11),
		      axis.text.y.left = element_text(face = "italic", color ="black"),
		      axis.text.y = element_text(size = 12)) +
		scale_y_continuous(limits = c(-4, 3), breaks = seq (from = -4, to = 3, by = 1)) +
		geom_hline(yintercept = 0, linetype="solid", color = "black", size = 0.3) +
		coord_flip()
	
	#cowplot::plot_grid(bacteria, virus, labels = c("a)", "b)"), label_size = 12, nrow = 2)
	
	# UNÍ BACTERIA Y VIRUS. LOS DEMÁS NO TUVIERON CAMBIOS.
	
	