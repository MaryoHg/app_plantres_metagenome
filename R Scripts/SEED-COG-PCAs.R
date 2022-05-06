### PCA computing for each database: Functionality or Metagenome

rm (list = ls()); dev.off()
library(ggplot2); library(FactoMineR)
set.seed(999)
options(scipen = 10000, digits = 3)
source('/Users/mariohg/Dropbox/R_functions/mytheme-PCA.R')

	#############
	############# Metagenome PCAs

	# Reading the clr-database for each database at higher taxonomic-level: uncomment the needed
	# SEED L1: phylum
		txmedata <- read.delim('/Users/mariohg/Dropbox/Norma/results/clr-database/SEED-L1-clr.tsv'); rownames(txmedata) <- txmedata$SampleID
		taxonomy <- read.delim('/Users/mariohg/Dropbox/Norma/results/clr-database/SEED-L1-ID.tsv')
		
	# COG L2
		txmedata <- read.delim('/Users/mariohg/Dropbox/Norma/results/clr-database/COG-L2-clr.tsv'); rownames(txmedata) <- txmedata$SampleID
		taxonomy <- read.delim('/Users/mariohg/Dropbox/Norma/results/clr-database/COG-L2-ID.tsv')
	
		
	# Compute the PCA
	pca1 <- FactoMineR::PCA(txmedata[,3:ncol(txmedata)])
	
	# Get Samples coords: joining with their metadata for samples and complete taxonomic levels.
	data.coords <- data.frame (pca1$ind$coord);
	data.coords$SampleID <- rownames(data.coords);
	data.coords <- merge (x = txmedata[,1:2], data.coords, by = "SampleID")
	
	# Get taxonomy information
	taxas <- data.frame(pca1$var$coord);
	taxas$FunctionID <- rownames(taxas)
	taxas <- base::merge(x = taxas, y = taxonomy[,1:2], by = "FunctionID")
	write.table(x = taxas, sep = '\t', quote = F, file = '/Users/mariohg/taxas.tsv', row.names = T)

	## Plotting the PCA plot
	head (data.coords)
	ggplot(data = data.coords, aes(x=Dim.1, y=Dim.2, color = Treatment, shape = Treatment)) + 
		geom_point(size = 3) +
		scale_x_continuous(limits=c(round(min (data.coords$Dim.1)-1), 
					    round(max (data.coords$Dim.1)+1)), 
				   breaks=seq(from = round(min (data.coords$Dim.1)-1), 
				   	   to = round(max (data.coords$Dim.1)+1), by = 2)) + 
		scale_y_continuous(limits=c(round(min (data.coords$Dim.2)-1), 
					    round(max (data.coords$Dim.2)+1)), 
				   breaks=seq(from = round(min (data.coords$Dim.2)-1), 
				   	   to = round(max (data.coords$Dim.2)+1), by = 2)) +
		mytheme +
		scale_shape_manual(values=c(16,16)) +
		scale_color_manual(values=c("red", "blue")) +
		geom_hline(yintercept=0, linetype="solid", color = "black") + 
		geom_vline (xintercept = 0, linetype="solid", color = "black") + 
		labs (x = paste ("PC 1 - Percent variation explained",
				 paste(round (pca1$eig[1,2], digits = 2), "%", sep = "")),
				 y = paste ("PC 2 - Percent variation explained",
				 	   paste(round (pca1$eig[2,2], digits = 2), "%", sep = "")))
	