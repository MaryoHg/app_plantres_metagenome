
## This function takes a dataframe in wide format with the first column as factor and 
##	do Welch's t-test (assumes unequal variances)
##
##	the data frame must have the following format
##
##	Treatment	Gene1		gene2		gene3		gene_n
##	A		count		count		count		count
##	A		count		count		count		count
##	A		count		count		count		count
##	B		count		count		count		count
##	B		count		count		count		count
##	.		.		.		.		.
##	.		.		.		.		.
##	.		.		.		.		.
##	Z		count_n		count_n		count_n		count_n		


# As function (with loop inside)
library ("WRS2")

loop.for.kruskal.wallis.test <- function(input.tab) {
        
        input.table <- input.tab
        
        taxa.names <- c(names(input.table[2:(ncol(input.table))]))
        outp.table <- data.frame(taxa.names, (rep(0, length(taxa.names)) ))
        names(outp.table) <-  c("Taxon", "p_value")
        
        for (i in 2 : (length(taxa.names)+1) ) {  
                
                
                res.perm <- kruskal.test(input.table[,i] ~ input.table[,1])$p.value
                p_value <- res.perm
                p_value <- round(p_value, 5)
                
                outp.table[(i-1),2] <- p_value
        }
        
        ;  outp.table
}
