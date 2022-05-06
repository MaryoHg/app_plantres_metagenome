### PCA computing for each database:

rm (list = ls()); dev.off()
library(ggplot2); library(FactoMineR)
set.seed(999)
options(scipen = 10000, digits = 3)
source('/Users/mariohg/Dropbox/R_functions/mytheme-PCA.R')
mapping <- read.delim('/Users/mariohg/mapping.txt')

	#############
	############# phylum taxonomic level

	# Reading the clr-database for each database at higher taxonomic-level: uncomment the needed
	# archaea: phylum
		#txmedata <- read.delim('/Users/mariohg/Dropbox/Norma/results/clr-database/arch-phyl-clr.tsv'); rownames(txmedata) <- txmedata$SampleID
	# bacteria: phylum
		#txmedata <- read.delim('/Users/mariohg/Dropbox/Norma/results/clr-database/bact-phyl-clr.tsv'); rownames(txmedata) <- txmedata$SampleID
	# Protozoa
		#txmedata <- read.delim('/Users/mariohg/Dropbox/Norma/results/clr-database/protoza-phyl-clr.tsv'); rownames(txmedata) <- txmedata$SampleID
	# Fungi
		#txmedata <- read.delim('/Users/mariohg/Dropbox/Norma/results/clr-database/fungi-phyl-clr.tsv'); rownames(txmedata) <- txmedata$SampleID
	# Virus
		txmedata <- read.delim('/Users/mariohg/Dropbox/Norma/results/clr-database/virus-family-clr.tsv'); rownames(txmedata) <- txmedata$SampleID
		
				
	# Compute the PCA
	pca1 <- FactoMineR::PCA(txmedata[,4:ncol(txmedata)])
	
	# Get Samples coords: joining with their metadata for samples and complete taxonomic levels.
	data.coords <- data.frame (pca1$ind$coord);
	data.coords$SampleID <- rownames(data.coords);
	data.coords <- merge (mapping, data.coords, by = "SampleID")
	#write.table(x = data.coords, file = '/Users/mariohg/site.coords.tsv', sep = '\t', row.names = F)
	
	# Get taxonomy information
	taxas <- data.frame(pca1$var$coord);
	taxas$GenusID <- rownames(taxas)
	#write.table(x = taxas, file = '/Users/mariohg/taxas.tsv', sep = '\t', row.names = T)
	
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
	
	
	
	#####################################
	########### genus/family (virus) taxonomic level

	# Loadint the crl-database at the family (virus) level
	# also load the ttest dataset to obtain the full taxonomic path
	# Bacteria
		txmedata <- read.delim('/Users/mariohg/Dropbox/Norma/results/clr-database/bact-genus-clr.tsv'); rownames(txmedata) <- txmedata$SampleID
		taxonomy <- read.delim('/Users/mariohg/Dropbox/Norma/results/clr-database/bact-genus-ID.tsv')
		
		data <- read.delim('/Users/mariohg/Dropbox/Norma/data/relab/bact/bact_rare_L6.txt', skip = 1)[,-1]
		data <- data[!(data$Level_6 %in% c("g__", "Other")),]
		data <- data %>% dplyr::arrange(desc(apply(X =data[,6:9], MARGIN = 1, FUN = mean))); data <- data[1:50,]
	
	# Archaea
		txmedata <- read.delim('/Users/mariohg/Dropbox/Norma/results/clr-database/arch-genus-clr.tsv'); rownames(txmedata) <- txmedata$SampleID
		taxonomy <- read.delim('/Users/mariohg/Dropbox/Norma/results/clr-database/arch-genus-ID.tsv')
		
		data <- read.delim('/Users/mariohg/Dropbox/Norma/data/relab/arch/arch_rare_L6.txt', skip = 1)[,-1]
		data <- data[!(data$Level_6 %in% c("g__", "Other")),]
		data <- data %>% dplyr::arrange(desc(apply(X =data[,6:9], MARGIN = 1, FUN = mean))); data <- data[1:50,]
	

	# Virus genus
		txmedata <- read.delim('/Users/mariohg/Dropbox/Norma/results/clr-database/virus-genus-clr.tsv'); rownames(txmedata) <- txmedata$SampleID
		taxonomy <- read.delim('/Users/mariohg/Dropbox/Norma/results/clr-database/virus-genus-ID.tsv')
		
		data <- read.delim('/Users/mariohg/Dropbox/Norma/data/relab/viral/viral_rare_L6.txt', skip = 1)[,c(-1:-3)]
		data <- data[!(data$Level_6 %in% c("g__", "Other")),]
		data <- data %>% dplyr::arrange(desc(apply(X =data[,4:7], MARGIN = 1, FUN = mean))); data <- data[1:50,]
		
	# Virus FAMILIES
		txmedata <- read.delim('/Users/mariohg/Dropbox/Norma/results/clr-database/virus-family-clr.tsv'); rownames(txmedata) <- txmedata$SampleID
		taxonomy <- read.delim('/Users/mariohg/Dropbox/Norma/results/clr-database/virus-family-ID.tsv')
		
		data <- read.delim('/Users/mariohg/Dropbox/Norma/data/relab/viral/viral_rare_L5.txt', skip = 1)[,c(-1:-3)]
		data <- data[!(data$Level_5 %in% c("f__", "Other")),]
		data <- data %>% dplyr::arrange(desc(apply(X =data[,3:6], MARGIN = 1, FUN = mean))); data <- data[1:50,]
		
		
	# Fungi: less than 50 genera
		txmedata <- read.delim('/Users/mariohg/Dropbox/Norma/results/clr-database/fungi-genus-clr.tsv'); rownames(txmedata) <- txmedata$SampleID
		taxonomy <- read.delim('/Users/mariohg/Dropbox/Norma/results/clr-database/fungi-genus-ID.tsv')
		
		#data <- read.delim('/Users/mariohg/Dropbox/Norma/data/relab/fungi/fungi_rare_L6.txt', skip = 1)[,-1]
		#data <- data[!(data$Level_6 %in% c("g__", "Other")),]
		#data <- data %>% dplyr::arrange(desc(apply(X =data[,6:9], MARGIN = 1, FUN = mean))); data <- data[1:50,]
		
	# Protozoa: less than 50 genera
		txmedata <- read.delim('/Users/mariohg/Dropbox/Norma/results/clr-database/protozoa-genus-clr.tsv'); rownames(txmedata) <- txmedata$SampleID
		taxonomy <- read.delim('/Users/mariohg/Dropbox/Norma/results/clr-database/protoza-genus-ID.tsv')
		
		#data <- read.delim('/Users/mariohg/Dropbox/Norma/data/relab/arch/arch_rare_L6.txt', skip = 1)[,-1]
		#data <- data[!(data$Level_6 %in% c("g__", "Other")),]
		#data <- data %>% dplyr::arrange(desc(apply(X =data[,6:9], MARGIN = 1, FUN = mean))); data <- data[1:50,]
		
		
		
	# Compute the PCA
	pca1 <- FactoMineR::PCA(txmedata[,4:ncol(txmedata)])
	
	# Get Samples coords: joining with their metadata for samples and complete taxonomic levels.
	# extract the 50+ abundant base on "data" dataframe
	data.coords <- data.frame (pca1$ind$coord);
	data.coords$SampleID <- rownames(data.coords);
	data.coords <- merge (mapping, data.coords, by = "SampleID")

	# Get taxonomy information
	taxas <- data.frame(pca1$var$coord);
	taxas$GenusID <- rownames(taxas) #change to FamilyID for virus family: same below
	taxas <- base::merge(x = taxas, y = taxonomy, by = "GenusID"); head (taxas)
	write.table(x = taxas, file = '/Users/mariohg/taxas.tsv', sep = '\t', row.names = F)
	# the 50+ abundant & save it:
	taxas50 <- taxas[taxas$Level_6 %in% data$Level_6,] #change to level_5 for virus family
	write.table(x = taxas50, file = '/Users/mariohg/taxas.tsv', sep = '\t', row.names = F)
	
	
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
	