## Source of commands: https://github.com/anttonalberdi/hilldiv/blob/master/documentation/div_test.md
options(scipen = 1000)
set.seed(999)
library(ggplot2)
library(hilldiv)

## analyzing the data from the SEED functional data
	## Reading the OTU table: rows are functions and columns are samples
	## Uncomment the functional level you need to test:
	#data <- read.delim('/Users/mariohg/Dropbox/Norma/data/Metagenome/SEED/SEED_L1_COUNT.txt')
	data <- read.delim('../data/functionality/SEED_L2_COUNT.txt', row.names = 1)[,1:4]; head(data)

	## Hierarchy or mapping file
	metadata <- read.delim('../data/metadata.txt', check.names = F)[,c(5,3)]; head(metadata)

	## Testing different q orders: uncomment the needed
	## first column of metadata must be colnmaes of "data"
	hilldiv::div_test(countable = data, qvalue = 0, hierarchy = metadata)
	hilldiv::div_test(countable = data, qvalue = 1, hierarchy = metadata)
	hilldiv::div_test(countable = data, qvalue = 2, hierarchy = metadata)
	

## analyzing the data from the SEED functional data
	## Reading the OTU table: rows are functions and columns are samples
	## Uncomment the functional level you need to test:
	data <- read.delim('/Users/mariohg/Dropbox/Norma/data/Metagenome/COGs/COG_COUNT_L2.txt')
	data <- read.delim('/Users/mariohg/Dropbox/Norma/data/Metagenome/COGs/COG_COUNT_L3.txt')
	rownames(data) <- data[,1]; data[,1]<-NULL
	
	## Hierarchy or mapping file
	metadata <- read.delim('/Users/mariohg/Dropbox/Norma/data/mapping.txt')[,-2]; metadata
	
	## Testing different q orders: uncomment the needed
	hilldiv::div_test(countable = data, qvalue = 0, hierarchy = metadata) # not significantly different
	hilldiv::div_test(countable = data, qvalue = 1, hierarchy = metadata) # not significantly different
	hilldiv::div_test(countable = data, qvalue = 2, hierarchy = metadata) # not significantly different
	
		
	