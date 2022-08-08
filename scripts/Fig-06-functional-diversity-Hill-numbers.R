## Figure 6: Functional diversity profile based on Hill numbers: COG and SEED databases at level 3
##################################################################################################
options(scipen = 10000, digits = 3)
set.seed(999)
rm(list = ls()); dev.off()
library(ggplot2); library(Rmisc)
metadata <- read.delim('../data/metadata.txt', check.names = F); metadata

## 1 Loading the raw SEED level 3 counts to transpose and determine the Hill numbers with Ma et al.
	data <- read.delim('../data/functionality/COG_COUNT_L3.txt'); View(data)
	rownames(data) <- paste ("Funct", 1:nrow(data), sep = ""); View(data)
	
	# before transposing we check NA's:
	which(!rowSums(!is.na(data)))
	
	data_t <- data.frame (t(data[,-1]))
	data_t <- cbind("SampleID" = rownames(data_t), data_t)
	
	write.table(x = data_t, sep = '\t', row.names = T, quote = F,
		    file = '/Users/mariohg/Dropbox/Norma/data/functional-diversity/COG_L3_transpose.txt')
	
	
## 2. Hill number per database per metagenome: uncomment the properly database
	## SEED level 3 functional diversity based on hill numbers
	data <- read.csv('../data/functionality/alpha-diversity-hill-numbers-COGL3_SEEDL3.csv', header = T)[,-6:-7]
	head(data)
	colnames(data)[1] <- "name"; head(data)
	data <- base::merge(x = metadata[,c(3,5)], y = data, by = "name"); data
	
	## make long to summarize per treatment
	data.long <- tidyr::gather(data = data[,-1],
				   key = "q.order",
				   value = "diversity",
				   colnames(data)[4:6]); data.long
	
	div.avg <- Rmisc::summarySE(data = data.long, 
				    measurevar = "diversity", 
				    groupvars = c("database", "Treatment", "q.order"), 
				    na.rm = T); div.avg
	
	# Plotting the results
	SEED <- 
	ggplot(data = subset(div.avg, database=="SEEDL3"), aes (x = q.order, y = diversity, group = Treatment)) +
		geom_errorbar (aes (ymin = diversity-se, ymax = diversity+se), colour="black", width=0.1) +
		geom_line (size=0.5) +
		geom_point (aes (shape = Treatment), size=4, fill="black") +
		scale_shape_manual(values = c(16,1)) +
		scale_x_discrete(name = "",
				 labels = c("q = 0","q = 1","q = 2")) +
		scale_y_continuous(name = "Effective number of functions", expand = c(0,0), limits = c(0,4000), breaks = seq(from=0, to=4000, by = 1000)) +
		theme_bw() + theme(legend.position = c(0.85,0.85), 
				   legend.title = element_text(colour="black", size=12, face="bold"),
				   legend.background =  element_rect(fill= NA, colour = NA),
				   legend.text = element_text(colour="black", size=12, face = "bold"),
				   axis.text.x=element_text(colour="black", size = 10),
				   axis.text.y=element_text(colour="black", size = 10)); SEED
	
		
		
		theme(legend.position = c(0.8,0.8),
				   axis.text.x=element_text(colour="black", size = 10),
				   axis.text.y=element_text(colour="black", size = 10)); SEED
	
	COG <- 
		ggplot(data = subset(div.avg, database=="COGL3"), aes (x = q.order, y = diversity, group = Treatment)) +
		geom_errorbar (aes (ymin = diversity-se, ymax = diversity+se), colour="black", width=0.1) +
		geom_line (size=0.5) +
		geom_point (aes (shape = Treatment), size=4, fill="black") +
		scale_shape_manual(values = c(16,1)) +
		#labs (x = "Diversity order (q values)", y = "Effective number of species") +
		scale_x_discrete(name = "Diversity order (q values)",
				 labels = c("q = 0","q = 1","q = 2")) +
		scale_y_continuous(name = "Effective number of functions", expand = c(0,0), limits = c(0,4000), breaks = seq(from=0, to=4000, by = 1000)) +
		theme_bw() + theme(legend.position = "none",
				   axis.text.x=element_text(colour="black", size = 10),
				   axis.text.y=element_text(colour="black", size = 10)); COG

	pdf(file ='../figures/Fig-06-fucntional-diversity.pdf', paper = "letter", width = 5, height = 7)
	cowplot::plot_grid(SEED, COG, labels = c("a)", "b)"), label_size = 12, nrow = 2,
					align = "v")
	dev.off()
	