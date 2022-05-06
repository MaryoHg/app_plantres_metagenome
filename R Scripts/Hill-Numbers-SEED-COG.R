## Source of commands: https://github.com/anttonalberdi/hilldiv/blob/master/documentation/div_test.md
options(scipen = 1000)
set.seed(999)
library(ggplot2)
library(hilldiv)


## tutorial
data(bat.diet.otutable); str (bat.diet.otutable)
data(bat.diet.tree); str (bat.diet.tree)
data(bat.diet.hierarchy); str (bat.diet.hierarchy)

hilldiv::div_test(countable = bat.diet.otutable, qvalue = 0, hierarchy = bat.diet.hierarchy)

hilldiv::div_test(countable = bat.diet.otutable,
		  qvalue = 1,
		  hierarchy = bat.diet.hierarchy,
		  posthoc = TRUE)

## analyzing the data from the SEED functional data
	## Reading the OTU table: rows are functions and columns are samples
	## Uncomment the functional level you need to test:
	#data <- read.delim('/Users/mariohg/Dropbox/Norma/data/Metagenome/SEED/SEED_L1_COUNT.txt')
	data <- read.delim('/Users/mariohg/Dropbox/Norma/data/Metagenome/SEED/SEED_L2_COUNT.txt')
	rownames(data) <- data$Function; data$Function<-NULL
	
	## Hierarchy or mapping file
	mapping <- read.delim('/Users/mariohg/Dropbox/Norma/data/mapping.txt')[,-2]; mapping
	
	## Testing different q orders: uncomment the needed
	hilldiv::div_test(countable = data, qvalue = 0, hierarchy = mapping)
	hilldiv::div_test(countable = data, qvalue = 1, hierarchy = mapping)
	hilldiv::div_test(countable = data, qvalue = 2, hierarchy = mapping)
	

## analyzing the data from the SEED functional data
	## Reading the OTU table: rows are functions and columns are samples
	## Uncomment the functional level you need to test:
	data <- read.delim('/Users/mariohg/Dropbox/Norma/data/Metagenome/COGs/COG_COUNT_L2.txt')
	data <- read.delim('/Users/mariohg/Dropbox/Norma/data/Metagenome/COGs/COG_COUNT_L3.txt')
	rownames(data) <- data[,1]; data[,1]<-NULL
	
	## Hierarchy or mapping file
	mapping <- read.delim('/Users/mariohg/Dropbox/Norma/data/mapping.txt')[,-2]; mapping
	
	## Testing different q orders: uncomment the needed
	hilldiv::div_test(countable = data, qvalue = 0, hierarchy = mapping) # not significantly different
	hilldiv::div_test(countable = data, qvalue = 1, hierarchy = mapping) # not significantly different
	hilldiv::div_test(countable = data, qvalue = 2, hierarchy = mapping) # not significantly different
	
		
	