#####Plotting the ratio at the genus taxonomic level
set.seed(999)
options(scipen = 1000, digits = 3)
library(ggplot2)
library(dplyr)

	# loading the raw entire dataset for phylum and family (only virus) relative abundance
	data <- read.delim('/Users/mariohg/Dropbox/Norma/data/relabavg/ratio_L6_databases.txt', comment.char = "#")
	data$Level_6 <- as.character(gsub(pattern = "g__", replacement = "", data$Level_6))
	colnames(data)[8] <- "taxonomy"
	
	# determine ratio
	data$ratio <- (data$MAS-data$UNS)/data$UNS
	
	## check if ratio was positive or negative: add color by type
	data <- data %>% 
		dplyr::mutate(mycolor = ifelse(test = ratio > 0, yes = "type1", no = "type2"))
	
# bacteria
	data1 <- subset (data, database=="bacteria")
	data1 <- data1 %>% dplyr::arrange(-desc(ratio))

	bact <- 
	ggplot(data1, aes(x=taxonomy, y=ratio)) +
		geom_segment( aes(x=taxonomy, xend=taxonomy, y=0, yend=ratio, color=mycolor), size=1.3, alpha=0.9) +
		geom_point( size=0.1, fill = "black", color="black", shape=21, stroke=2) +
		scale_x_discrete(limits=c(as.character(data1$taxonomy))) +
		scale_y_continuous(limits = c(-1, 3), breaks = seq(-1, 3, by= 1)) + #defining ratio axis values
		theme_bw() + #theme_light() +
		theme(axis.text.y = element_text(face = "italic", hjust = 1, size = 11, color = "black"),#, colour = "black"), #tick labels
		      axis.text.x= element_text(colour ="black", size = 11, color = "black"), #tick labels
		      legend.position = "none") + #panel.border = element_blank()
		xlab("") +
		ylab("Ratio change") + coord_flip() +
		geom_hline(yintercept=0, linetype="solid", color = "black", alpha=0.5); bact
	
	
# arch
	data1 <- subset (data, database=="archaea")
	data1 <- data1 %>% dplyr::arrange(-desc(ratio))
	
	arch <- 
		ggplot(data1, aes(x=taxonomy, y=ratio)) +
		geom_segment( aes(x=taxonomy, xend=taxonomy, y=0, yend=ratio, color=mycolor), size=1.3, alpha=0.9) +
		geom_point( size=0.1, fill = "black", color="black", shape=21, stroke=2) +
		scale_x_discrete(limits=c(as.character(data1$taxonomy))) +
		scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by= 0.5)) + #defining ratio axis values
		theme_bw() + #theme_light() +
		theme(axis.text.y = element_text(face= "italic", hjust = 1, size = 11, color = "black"),#, colour = "black"), #tick labels
		      axis.text.x= element_text(colour ="black", size = 11, color = "black"), #tick labels
		      legend.position = "none") + #panel.border = element_blank()
		xlab("") +
		ylab("Ratio change") + coord_flip() +
		geom_hline(yintercept=0, linetype="solid", color = "black", alpha=0.5); arch
	
	
# fungi
	data1 <- subset (data, database=="fungi")
	data1 <- data1 %>% dplyr::arrange(-desc(ratio))
	
	fungi <- 
		ggplot(data1, aes(x=taxonomy, y=ratio)) +
		geom_segment( aes(x=taxonomy, xend=taxonomy, y=0, yend=ratio, color=mycolor), size=1.3, alpha=0.9) +
		geom_point( size=0.1, fill = "black", color="black", shape=21, stroke=2) +
		scale_x_discrete(limits=c(as.character(data1$taxonomy))) +
		#scale_y_continuous(limits = c(-1, 0.6), breaks = seq(-1, 0.6, by= 0.2)) + #defining ratio axis values
		theme_bw() + #theme_light() +
		theme(axis.text.y = element_text(face = "italic", hjust = 1, size = 11, color = "black"),#, colour = "black"), #tick labels
		      axis.text.x= element_text(colour ="black", size = 11, color = "black"), #tick labels
		      legend.position = "none") + #panel.border = element_blank()
		xlab("") +
		ylab("Ratio changes") + coord_flip() +
		geom_hline(yintercept=0, linetype="solid", color = "black", alpha=0.5); fungi
	
	
# protoz
	data1 <- subset (data, database=="protozoa")
	data1 <- data1 %>% dplyr::arrange(-desc(ratio))
	
	protoz <- 
		ggplot(data1, aes(x=taxonomy, y=ratio)) +
		geom_segment( aes(x=taxonomy, xend=taxonomy, y=0, yend=ratio, color=mycolor), size=1.3, alpha=0.9) +
		geom_point( size=0.1, fill = "black", color="black", shape=21, stroke=2) +
		scale_x_discrete(limits=c(as.character(data1$taxonomy))) +
		scale_y_continuous(limits = c(-0.2, .60), breaks = seq(-0.2, .6, by= 0.2)) + #defining ratio axis values
		theme_bw() + #theme_light() +
		theme(axis.text.y = element_text(face = "italic", hjust = 1, size = 11, color = "black"),#, colour = "black"), #tick labels
		      axis.text.x= element_text(colour ="black", size = 11, color = "black"), #tick labels
		      legend.position = "none") + #panel.border = element_blank()
		xlab("") +
		ylab("Ratio changes") + coord_flip() +
		geom_hline(yintercept=0, linetype="solid", color = "black", alpha=0.5); protoz
	
	
# virus: the 50 most abundant genera across treatments
	data1 <- subset (data, database=="virus")
	data1 <- data1 %>% dplyr::arrange(-desc(ratio))
	
	virus <- 
		ggplot(data1, aes(x=taxonomy, y=ratio)) +
		geom_segment( aes(x=taxonomy, xend=taxonomy, y=0, yend=ratio, color=mycolor), size=1.3, alpha=0.9) +
		geom_point( size=0.1, fill = "black", color="black", shape=21, stroke=2) +
		scale_x_discrete(limits=c(as.character(data1$taxonomy))) +
		scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, by = .25)) + #defining ratio axis values
		theme_bw() + #theme_light() +
		theme(axis.text.y = element_text(face = "italic", hjust = 1, size = 11, color = "black"),#, colour = "black"), #tick labels
		      axis.text.x= element_text(colour ="black", size = 11, color = "black"), #tick labels
		      legend.position = "none") + #panel.border = element_blank()
		xlab("") +
		ylab("Ratio change") + coord_flip() +
		geom_hline(yintercept=0, linetype="solid", color = "black", alpha=0.5); virus
	
## Joining the plots
	
	# Joining the three small plots: arch, fungi and protozoa
	cowplot::plot_grid(bact, arch, fungi, protoz, virus, labels = c("a)", "b)", "c)", "d)", "e)"), label_size = 12, ncol = 5, align = "h")

	