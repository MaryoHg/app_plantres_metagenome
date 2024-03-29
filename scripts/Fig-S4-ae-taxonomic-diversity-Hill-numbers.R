## Hill numbers: taxonomic diversity at the species count level datasets
########################################################################

rm(list = ls()); dev.off()
set.seed(999)
library(ggplot2); library(Rmisc)
#source('/Users/mariohg/Dropbox/R_functions/summarySE.R')
options(scipen = 10000, digits = 3)
metadata <- read.delim('../data/metadata.txt', check.names = F); head(metadata)


## 1. Loading the alpha taxonomic diversity values for all databases
	# Hill number per database per metagenome
	# merge with the metadata
	diversity <- read.delim('../data/alpha-diversity-metagenome.csv', sep = ',',
			   header = T, comment.char = "#")[,c(-6,-7)]; head(diversity)
	diversity$treatment <- as.character(gsub(pattern = "Bracken", replacement = "", diversity$treatment))
	colnames(diversity)[2] <- "nameid1"
	head (diversity)
			   #/data.tsv, 
	div.long <- tidyr::gather(data = diversity,
				   key = "q.order",
				   value = "diversity",
				   colnames(diversity)[3:5], factor_key = TRUE)
	div.long <- base::merge(x = metadata[,c(3,6)], y = div.long, by = "nameid1")
	
## 2. Average per treatment and taxonomic database
	div.avg <- Rmisc::summarySE(data = div.long, 
				    measurevar = "diversity", 
				    groupvars = c("database", "Treatment", "q.order"), na.rm = T)
	
	# Bacteria
	bacteria <- 
	ggplot(data = subset(div.avg, database=="Bacteria"), aes (x = q.order, y = diversity, group = Treatment)) +
		geom_errorbar (aes (ymin = diversity-se, ymax = diversity+se), colour="black", width=0.1) +
		geom_line (size=0.5) +
		geom_point (aes (shape = Treatment), size=4, fill="black") +
		scale_shape_manual(values = c(16,1)) +
		#labs (x = "Diversity order (q values)", y = "Effective number of species") +
		scale_x_discrete(name = "", labels = c("q = 0","q = 1","q = 2")) +
		scale_y_continuous(name = "Effective number of species", expand = c(0,0), limits = c(0,4000), breaks = seq(from=0, to=4000, by = 1000)) +
		theme_bw() + theme(legend.position = c(0.85,0.85), 
				   legend.title = element_text(colour="black", size=12, face="bold"),
				   legend.background =  element_rect(fill= NA, colour = NA),
				   legend.text = element_text(colour="black", size=12, face = "bold")); bacteria
	
	# Archaea
	archaea <- 
	ggplot(data = subset(div.avg, database=="Archaea"), aes (x = q.order, y = diversity, group = Treatment)) +
		geom_errorbar (aes (ymin = diversity-se, ymax = diversity+se), colour="black", width=0.1, position = position_dodge(width = 0.1)) +
		geom_line (size=0.5, position = position_dodge(width = 0.1)) +
		geom_point (aes (shape = Treatment), size=4, fill="black", position = position_dodge(width = 0.1)) +
		scale_shape_manual(values = c(16,1)) +
		#labs (x = "Diversity order (q values)", y = "Effective number of species") +
		scale_x_discrete(name = "",
				 labels = c("q = 0","q = 1","q = 2")) +
		scale_y_continuous(name = "", expand = c(0,0), limits = c(0,250), breaks = seq(from=0, to=250, by = 50)) +
		theme_bw() + theme(legend.position = "none"); archaea
	
	# Fungi
	fungi <- 
	ggplot(data = subset(div.avg, database=="Fungi"), aes (x = q.order, y = diversity, group = Treatment)) +
		geom_errorbar (aes (ymin = diversity-se, ymax = diversity+se), colour="black", width=0.1, position = position_dodge(width = 0.1)) +
		geom_line (size=0.5, position = position_dodge(width = 0.1)) +
		geom_point (aes (shape = Treatment), size=4, fill="black", position = position_dodge(width = 0.1)) +
		scale_shape_manual(values = c(16,1)) +
		#labs (x = "Diversity order (q values)", y = "Effective number of species") +
		scale_x_discrete(name = "",
				 labels = c("q = 0","q = 1","q = 2")) +
		scale_y_continuous(name = "Effective number of species", expand = c(0,0), limits = c(0,60), breaks = seq(from=0, to=60, by = 10)) +
		theme_bw() + theme(legend.position = "none"); fungi
	
	# Protozoa
	protozoa <- 
	ggplot(data = subset(div.avg, database=="Protozoa"), aes (x = q.order, y = diversity, group = Treatment)) +
		geom_errorbar (aes (ymin = diversity-se, ymax = diversity+se), colour="black", width=0.1, position = position_dodge(width = 0.1)) +
		geom_line (size=0.5, position = position_dodge(width = 0.1)) +
		geom_point (aes (shape = Treatment), size=4, fill="black", position = position_dodge(width = 0.1)) +
		scale_shape_manual(values = c(16,1)) +
		#labs (x = "Diversity order (q values)", y = "Effective number of species") +
		scale_x_discrete(name = "",
				 labels = c("q = 0","q = 1","q = 2")) +
		scale_y_continuous(name = "", expand = c(0,0), limits = c(10,40), breaks = seq(from=10, to=40, by = 5)) +
		theme_bw() + theme(legend.position = "none"); protozoa
	
	# Virus
	virus <- 
	ggplot(data = subset(div.avg, database=="Virus"), aes (x = q.order, y = diversity, group = Treatment)) +
		geom_errorbar (aes (ymin = diversity-se, ymax = diversity+se), colour="black", width=0.1, position = position_dodge(width = 0.1)) +
		geom_line (size=0.5, position = position_dodge(width = 0.1)) +
		geom_point (aes (shape = Treatment), size=4, fill="black", position = position_dodge(width = 0.1)) +
		scale_shape_manual(values = c(16,1)) +
		#labs (x = "Diversity order (q values)", y = "Effective number of species") +
		scale_x_discrete(name = "Diversity order (q values)",
				 labels = c("q = 0","q = 1","q = 2")) +
		scale_y_continuous(name = "Effective number of species", expand = c(0,0), limits = c(0,5500), breaks = seq(from=0, to=5500, by = 1000)) +
		theme_bw() + theme(legend.position = "none", plot.margin = unit(c(0,0,0,0), "cm")); virus
	
	# Joining all the plots:
	cowplot::plot_grid(bacteria, archaea, fungi, protozoa, virus,
			   labels = c("a)", "b)", "c)", "d)", "e)"), label_size = 12,
			   nrow = 3, ncol = 2, align = "v")
	