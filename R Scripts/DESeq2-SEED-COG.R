## DESeq2 on metagenome: NINGUNA FUNCIÓN CAMBIO SIGNIFICATIVAMENTE
#######################
rm (list = ls()); dev.off()
library(ggplot2)
set.seed(999)
library(DESeq2); BiocParallel::register(MulticoreParam(4))
library(pheatmap)

# 1. Loading data at the genus taxonomic level: COUNTS VALUES
	# LUncomment the dataset you want to test:
	#data <- read.delim('/Users/mariohg/Dropbox/Norma/data/Metagenome/SEED/SEED_L1_COUNT.txt')
	data <- read.delim('/Users/mariohg/Dropbox/Norma/data/Metagenome/COGs/COG_COUNT_L2.txt')

	
	# Editing and formating the data for proper analyses
	rownames(data) <- data[,1]; View(data)
	data <- data[,-1]; head(data)

	# Loading the metadata
	mapping <- data.frame ("SampleID" = colnames(data), 
			       "Treatment" = rep (c("MAS", "UNS"), each = 2)); mapping
	
# 2. DESeq2 Matrix procedure
	countData <- as.matrix(data+1); head(countData)
	dds <- DESeq2::DESeqDataSetFromMatrix(countData = countData, colData = mapping, design = as.formula(~ Treatment))
	res <- DESeq2::DESeq(object = dds, parallel = TRUE)
	
# 3. Comparison between treatments: something like ratio but with statistics
	# Comparison between treatments: MAS vs UNS
	results <- DESeq2::results(object = res, contrast = c("Treatment", "MAS", "UNS"))
	
	# convert to dataframe
	output <- data.frame(results); nrow(output)
	
	# filter those highly significant and log2FoldChange higher than abs (1)
	outfiltered <- subset(output, padj <= 0.05) #& abs(log2FoldChange) > 1)
	nrow (outfiltered)
	
	# merge with their corresponding taxonomic levels
	outfiltered$Phylum <- rownames(outfiltered)

	# editting the results for plotting
	# delete the g__ patter,
	# add a factor for color plotting
	# Ordering by logFoldChange
	outfiltered$Phylum <- as.character(gsub(pattern = "p__", replacement = "", outfiltered$Phylum))
	outfiltered <- outfiltered %>% dplyr::mutate (mycolor = ifelse(test = log2FoldChange > 0, yes = "type1", no = "type2"))
	outfiltered <- outfiltered %>% dplyr::arrange(desc(outfiltered$log2FoldChange))
	
# 5. Plotting
	colnames(outfiltered)
	#virus <- 
	ggplot(data = outfiltered, aes(x= Phylum, y=log2FoldChange)) +
		geom_point(size = 2, color= "black", shape=16, stroke=2) +
		theme_light() + labs (x = "Biomarkers (Phylum level)", y = "LFC (Log2FoldChange") +
		scale_color_manual(values = c("black", "gray50")) +
		scale_x_discrete(limits=c(as.character(outfiltered$Phylum))) +
		scale_y_continuous(limits = c(-3, 3), breaks = seq (from = -3, to = 3, by = 1)) +
		geom_hline(yintercept = 0, linetype="solid", 
			   color = "black", size = 0.3) +
		coord_flip() + theme(legend.position = "none",
		      panel.border = element_blank(), axis.text.x = element_text(size = 11),
		      axis.text.y.left = element_text(color ="black"),
		      axis.text.y = element_text(size = 11))
	
	#cowplot::plot_grid(bacteria, virus, labels = c("a)", "b)"), label_size = 12, nrow = 2)
	
	# UNÍ BACTERIA Y VIRUS. LOS DEMÁS NO TUVIERON CAMBIOS.
	
	