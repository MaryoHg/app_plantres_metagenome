## Hill numbers: taxonomy at species count level dataset
	
rm(list = ls()); dev.off()
set.seed(999)
library(ggplot2);
source('/Users/mariohg/Dropbox/R_functions/summarySE.R')
options(scipen = 10000, digits = 3)

	# 1 Loading the raw SEED 3 count to transpose and determine the Hill numbers with Ma et al.
	data <- read.delim('/Users/mariohg/Dropbox/Norma/data/functional-diversity/COG_COUNT_L3.txt')
	rownames(data) <- paste ("Funct", 1:nrow(data), sep = ""); head(data)
	
	# before transposing we check NA's:
	which(!rowSums(!is.na(data)))
	
	data_t <- data.frame (t(data[,-1]))
	data_t <- cbind("SampleID" = rownames(data_t), data_t)
	
	write.table(x = data_t, sep = '\t', row.names = T, quote = F,
		    file = '/Users/mariohg/Dropbox/Norma/data/functional-diversity/COG_L3_transpose.txt')
	
	
	# Hill number per database per metagenome: uncomment the properly database
	data <- read.csv('/Users/mariohg/Dropbox/Norma/data/functional-diversity/alpha-SEED3.csv', header = T)[,-5:-6]
	#data <- read.csv('/Users/mariohg/Dropbox/Norma/data/functional-diversity/alpha-COG3.csv', header = T)[,-5:-6]
	data$Treatment <- rep (c("MAS", "UNS"), each = 2); data
	data.long <- tidyr::gather(data = data,
				   key = "q.order",
				   value = "diversity",
				   colnames(data)[2:4]); data.long
	
	div.avg <- summarySE(data = data.long, measurevar = "diversity", groupvars = c("Treatment", "q.order"), na.rm = T)
	
	# Plotting the results
	SEED <- 
	ggplot(data = div.avg, aes (x = q.order, y = diversity, group = Treatment)) +
		geom_errorbar (aes (ymin = diversity-se, ymax = diversity+se), colour="black", width=0.1) +
		geom_line (size=0.5) +
		geom_point (aes (shape = Treatment), size=4, fill="black") +
		scale_shape_manual(values = c(16,1)) +
		#labs (x = "Diversity order (q values)", y = "Effective number of species") +
		scale_x_discrete(name = "",
				 labels = c("q = 0","q = 1","q = 2")) +
		scale_y_continuous(name = "Effective number of functions", expand = c(0,0), limits = c(0,4000), breaks = seq(from=0, to=4000, by = 1000)) +
		theme_bw() + theme(legend.position = c(0.8,0.8),
				   axis.text.x=element_text(colour="black", size = 10),
				   axis.text.y=element_text(colour="black", size = 10)); SEED
	
	COG <- 
		ggplot(data = div.avg, aes (x = q.order, y = diversity, group = Treatment)) +
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

	pdf(file = "/Users/mariohg/Renv/fig6_plot.pdf", paper = "letter", width = 5, height = 7)
	cowplot::plot_grid(SEED, COG, labels = c("a)", "b)"), label_size = 12, nrow = 2,
					align = "v")
	dev.off()
	