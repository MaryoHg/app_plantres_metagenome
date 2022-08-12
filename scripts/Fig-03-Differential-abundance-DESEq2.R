		#############################################################
		## Figure 3a : Differential abundance analysis ob bacteria ##
		#############################################################
rm (list = ls()); dev.off()
library(ggplot2)
set.seed(999)
library(DESeq2); BiocParallel::register(BiocParallel::MulticoreParam(4)); library(dplyr)
metadata <- read.delim('../data/metadata.txt', check.names = F, stringsAsFactors = T)[,c(7,3)]

## 1. Loading data at the genus taxonomic level: COUNTS VALUES
	#a) Bacteria at the phylum taxonomic level (Fig. 3a)
	data <- read.delim('../data/taxonomy_counts/bact_rare_L2.txt', comment.char = '#')[,-1] %>% 
	tibble::column_to_rownames('Level_2')
	head(data)
	
## 2. DESeq2 Matrix procedure: pseudo-count of +1
	countData <- as.matrix(data+1); head (countData)
	dds <- DESeq2::DESeqDataSetFromMatrix(countData = countData, colData = metadata, design = as.formula(~ Treatment))
	res <- DESeq2::DESeq(object = dds, parallel = TRUE)
	
# 3. Comparison between treatments:
	# a. Comparison between treatments: MAS vs UNS
	# b. convert deseq2 output to data.frame()
	# c. filter those highly significantly different between treatments and log2FoldChange higher than abs (1)
	# d. delete the p__ pattern,
	
	results <- DESeq2::results(object = res, contrast = c("Treatment", "MAS", "UNS")) %>% 
		as.data.frame() %>% 
		filter(padj <= 0.05 & abs(log2FoldChange) > 1) %>% 
		print() %>%
		tibble::rownames_to_column(., var = "TaxaID") %>% # use rownames and add it as column
		mutate(TaxaID = gsub(pattern = "p__", replacement = "", TaxaID)) # delete the p__ patter


	# 5. Plotting
	results
	pdf(file = '../figures/Fig-3a-Phyla-diff-abundance.pdf', width = 4, height = 2, paper = "letter")
	ggplot(data = results, aes(x= TaxaID, y=log2FoldChange)) +
		geom_point(size = 2, color= "black", shape=16, stroke=2) +
		theme_light() + labs (x = "", y = "Log2FoldChange") +
		scale_x_discrete(limits=c(as.character(results$TaxaID))) +
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
	
	## save the .RData file
	save.image(file = "../environment_Rdata/Fig-03-Differential-abundance-DESEq2.RData")
	
	####################################################################
	## Figure 3b: Differential abundance analysis on bacterial genera ##
	####################################################################
	
	## 1. Loading data at the genus taxonomic level: COUNTS VALUES
	# Editing and formating the data for proper analyses
	# keep the matrix of samples, i.e., only 4 columns
	datagen <- read.delim('../data/taxonomy_counts/bact_rare_L6.txt', check.names = F, comment.char = '#')[,-1] %>%
		magrittr::set_rownames(paste ("Genus", 1:nrow(.), sep = "_")) #%>% # . indicate lazytable
	# select(ends_with(c("A", "B"))) %>%
	# as.matrix()
	
	## taxonomy levels information: loading data and keep the taxonomic levels of bacterial genera
	taxonomy <- cbind("GenusID"=rownames(datagen), datagen[,1:5])
	
	# taxonomy <- read.delim('../data/taxonomy_counts/bact_rare_L6.txt', check.names = F, comment.char = '#')[,-1] %>%
	# 	magrittr::set_rownames(paste ("Genus", 1:nrow(.), sep = "_")) %>%
	# 	select(starts_with("Level")) %>% 
	# 	tibble::rownames_to_column(., var = "TaxaID")
	
	# Loading the metadata
	mapping <- read.delim('../data/metadata.txt')[,c(7,3)]; mapping
	
	## 2. DESeq2 Matrix procedure
	## adding +1 to each count value: pseudovalue to avoid zeros
	## run de deseq2 function
	countDatag <- as.matrix(datagen[,6:9]+1); head (countDatag)
	ddsg <- DESeq2::DESeqDataSetFromMatrix(countData = countDatag, colData = mapping, design = as.formula(~ Treatment))
	resg <- DESeq2::DESeq(object = ddsg, parallel = TRUE)
	
	## 3. Comparison between treatments: 
	# Comparison between treatments: MAS vs UNS
	# convert to dataframe
	# filter those highly significant and log2FoldChange higher than abs (1)
	# merge with their corresponding taxonomic levels
	# add a factor for color plotting
	# Ordering by logFoldChange
	resultsgen <- DESeq2::results(object = resg, contrast = c("Treatment", "MAS", "UNS")) %>% 
		as.data.frame() %>% 
		filter(padj <= 0.05 & abs(log2FoldChange) > 1) %>% 
		tibble::rownames_to_column(., var = "GenusID") %>% # use rownames and add it as column
		inner_join(x = taxonomy[,c(1,2,6)], by = "GenusID") %>%
		mutate(Level_6 = gsub(pattern = "g__", replacement = "", Level_6)) %>%  # delete the g__ pattern
		mutate(Level_2 = gsub(pattern = "p__", replacement = "", Level_2)) %>% 
		mutate (mycolor = ifelse(test = log2FoldChange > 0, yes = "type1", no = "type2")) %>% # add a factor column for ggplot coloring
		#arrange(desc(log2FoldChange))
		arrange(desc(Level_2))
	
	
	# 5. Plotting
	colnames(resultsgen)
	pdf(file = '../figures/Fig-3b-Phyla-diff-abundance.pdf', width = 5.5, height = 6, paper = "letter")
	ggplot(data = resultsgen, aes(x= Level_6, y=log2FoldChange)) +
		geom_point(aes(color = Level_2), size = 2, shape=16, stroke=2) +
		theme_light() + labs (x = "", y = "Log2FoldChange") +
		#scale_color_manual(values = RColorBrewer::brewer.pal(5, "Accent")) +
		scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")) +
		scale_x_discrete(limits=c(as.character(resultsgen$Level_6))) +
		geom_hline(yintercept = -2, linetype="dashed", color = "gray90", size= 0.75) +
		geom_hline(yintercept = 2, linetype="dashed", color = "gray90", size= 0.75) +
		theme(#legend.position = "none", #"left",
			panel.border = element_rect(colour = "black", size = 1), 
			axis.text.y = element_text(face = "italic", color ="black", size = 12),
			axis.text.x = element_text(size = 12, color = "black"),
			panel.grid = element_blank()) +
		scale_y_continuous(limits = c(-3, 3), breaks = seq (from = -4, to = 3, by = 1)) +
		geom_hline(yintercept = 0, linetype = "solid", color = "gray90", size = 0.5) +
		coord_flip()
	dev.off()
	
	## save the .RData file
	save.image(file = "../environment_Rdata/Fig-03-Differential-abundance-DESEq2.RData")
	