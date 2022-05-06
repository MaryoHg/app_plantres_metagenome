#####Plotting the ratio at the genus taxonomic level
set.seed(999)
options(scipen = 1000, digits = 3)
library(ggplot2); library(dplyr)

	# read the raw relative abundance table
	# change to percent, i.e. multiply by 100
	# ordering: from the dominant to the low dominant on average across samples
	
# Bacteria
	data <- read.delim('/Users/mariohg/Dropbox/Norma/data/relab/bact/bact_rare_L6.txt', comment.char = "#")
	
# Archaea
	data <- read.delim('/Users/mariohg/Dropbox/Norma/data/relab/arch/arch_rare_L6.txt', comment.char = "#")

# Fungi
	data <- read.delim('/Users/mariohg/Dropbox/Norma/data/relab/fungi/fungi_rare_L6.txt', comment.char = "#")

# Protists
	data <- read.delim('/Users/mariohg/Dropbox/Norma/data/relab/protozoa/protozoa_rare_L6.txt', comment.char = "#")
	
# Virus
	data <- read.delim('/Users/mariohg/Dropbox/Norma/data/relab/viral/viral_rare_L6.txt', comment.char = "#")
	
	
	# edit the table for downstream calculation
	data <- data[!(data$Level_6 %in% c("g__", "Other")),]
	data$Level_6 <- as.character(gsub(pattern = "g__", replacement = "", x = data$Level_6))
	
	# change to percent and order from high to low abundance on average
	data[,7:10] <- data[,7:10]*100
	data <- data %>% dplyr::arrange(desc(apply(X = data[,7:10], MARGIN = 1, FUN = mean)))
	
	# determine the mean±SD for all treatments:
	data$AVG <- apply(X = data[,7:10], MARGIN = 1, FUN = mean)
	data$SD <- apply(X = data[,7:10], MARGIN = 1, FUN = sd)
	
	# determine the MAS±sd and UNS±sd average treatments
	colnames(data)
	data$MAS <- apply(X = data[,9:10], MARGIN = 1, FUN = mean)
	data$MASsd <- apply(X = data[,9:10], MARGIN = 1, FUN = sd)
	
	data$UNS <- apply(X = data[,7:8], MARGIN = 1, FUN = mean)
	data$UNSsd <- apply(X = data[,7:8], MARGIN = 1, FUN = sd)
	
	# determine the ratio based on both average vales
	data$ratio <- (data$MAS - data$UNS) / data$UNS
	
	# determine the X-times of increased: negative of positive taken into account with if-else()
	data <- data %>% 
		dplyr::mutate(times_of_increased = ifelse(test = data$MAS > data$UNS, yes = data$MAS/data$UNS, no = data$UNS/data$MAS))
	
	# add the "positive" or "negative" change as a factor to plot with colors:
	data <- data %>% 
		dplyr::mutate(ratio_type = ifelse(test = ratio > 0, yes = "positive", no = "negative"))

	# filter the 50 more dominant genera
	# ordering by lower to higher ratio change
	# plot the results
	data1 <- data[1:50,] 
	data1 <- data1 %>% dplyr::arrange(-desc(AVG))

	ggplot2::ggplot(data1, aes(x=Level_6, y=ratio)) +
		geom_segment(aes(x=Level_6, xend=Level_6, y=0, yend=ratio, color=ratio_type), size=1.3, alpha=0.9) +
		geom_point(size=0.1, fill = "black", color="black", shape=21, stroke=2) +
		scale_x_discrete(limits=c(as.character(data1$Level_6))) +
		#scale_y_continuous(limits = c(-1, 3), breaks = seq(-1, 3, by= 0.5)) + #defining ratio axis values
		#scale_y_continuous(limits = c(-0.8, 0.8), breaks = seq(from = -0.8, to = 0.8, by = 0.2)) +
		theme_bw() + #theme_light() +
		theme(axis.text.y = element_text(face = "italic", hjust = 1, size = 11, color = "black"),#, colour = "black"), #tick labels
		      axis.text.x= element_text(colour ="black", size = 11, color = "black"), #tick labels
		      legend.position = "none") + #panel.border = element_blank()
		xlab("") +
		ylab("Ratio change") + coord_flip() +
		geom_hline(yintercept=0, linetype="solid", color = "black", alpha=0.5) +
		theme(panel.grid.minor.x = element_blank())
	
	