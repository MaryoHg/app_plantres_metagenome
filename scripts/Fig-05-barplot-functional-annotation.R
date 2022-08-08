################################################
## Dando formato a SEED functional annotation ##
################################################
rm (list = ls())
library(ggplot2)
library (tidyr)
library (reshape2)
source('/home/mario/Dropbox/R_functions/summarySE.R')
source('/home/mario/Dropbox/R_functions/tranpose.any.table.R')
source('/home/mario/Dropbox/R_functions/loop_for_Welch_test.R')

##1. Loading data and rename headers; create the factor for this data
        seed.L1.raw <- read.delim('/home/mario/Documents/Doctorado/PROYECTOS/norma-jimenez/shotgun_sonora/megahit_sonora/prokka_vs_nr/MEGAN_ANALYSIS/SEED/SEED_L1_COUNT.txt')
        colnames(seed.L1.raw)[1] <- "Function";
        View (seed.L1.raw)
        Functions <- seed.L1.raw[,1]
        
        ## Creating factors for statistics
        treatments <- factor(rep (c("MAS", "UNS"), each =2))
        treatments
        
        
        ## a) Formatting for barplot: wide to long
        seed.L1.long <- melt (seed.L1.raw,
                              id.vars="Function",
                              variable.name='Replicates',
                              value.name="reads_assigned")

        ## Adding TREATMENT (UNS and MAS):
        seed.L1.long$TREATMENT <- rep (c("MAS", "UNS"), each =90)
        View(seed.L1.long)
        
        ## b) AVERAGE and statistics: summarySE()
        ## Filtering samples with >=50 reads assigned:
        seed.L1.long.50.most <- subset(seed.L1.long, reads_assigned>=33) #46*
        
        seed.L1.long.50.most.summary <- summarySE(data = seed.L1.long.50.most, measurevar = "reads_assigned",
                                                    groupvars = c("Function", "TREATMENT"), na.rm = T)
        View (seed.L1.long.50.most.summary)
        

## 2.0 Ploting: BAR PLOT
        SEED_L1 <- ggplot(seed.L1.long.50.most.summary[seed.L1.long.50.most.summary$Function!="Plant cell walls and outer surfaces",],
                          aes(x=Function, y=reads_assigned,fill=TREATMENT)) +
                geom_errorbar (aes (ymin=reads_assigned-se, ymax=reads_assigned+se), 
                               colour="black", width=0.3, position=position_dodge(0.9)) +
                geom_bar(position=position_dodge(), stat = 'identity', width = 0.9) +
                scale_fill_manual(values = c("#FF0000","#0000FF"), # red and blue # 
                                  name="Treatment", 
                                  labels=c("MAS","UNS")) +
                xlab("") + ylab("Number of reads") +
                theme_bw() + scale_y_continuous(expand = c(0, 0), limits = c(0,8000)) + ## must be up here
                theme (axis.text.x = element_text(angle = 70, hjust = 1)) ## x-axis names angle.
        SEED_L1
        
        SEED_L1 + coord_flip()
        
        ## Saving as pdf file:
        # a4r is a landscape plot.
        pdf(file="/home/mario/Documents/Doctorado/PROYECTOS/norma-jimenez/shotgun_sonora/megahit_sonora/prokka_vs_nr/MEGAN_ANALYSIS/SEED/SEED.L1.pdf",
            width=11, height=8.5, paper='a4r')
        
        SEED_L1 + theme( 
                # remove the vertical grid lines
                panel.grid.major.x = element_blank() ,
                # explicitly set the horizontal lines (or they will disappear too)
                panel.grid.major.y = element_line( size=0.1, color="black" ),
                # delete the border)
                panel.border = element_blank(),
                # add only x/y axes only
                axis.line = element_line(colour = "black"),
                # x axis labels
                axis.text.x = element_text(color = "black", size = 11),
                # y axis labels
                axis.text.y = element_text(color = "black", size = 11),
                # add a margin to the graph: margin(up,right,down,left)
                panel.background = element_rect(fill = "white"),
                plot.margin = margin(1.5, 1, 0, 1, "cm"),
                plot.background = element_rect(fill = "white", colour = "black", size = 1)     )
        dev.off()
                ## NOTE: THE BARPLOT WAS SAVED AS - SEED_L1.pdf
        
        
## 3.0 Batch Welch's two sample t-test:
        ## transpose the raw table 
        seed.wide <- transpose.any.table(seed.L1.raw)
        seed.L1.wide <- cbind("Treatment"=treatments, seed.wide)
        View (seed.L1.wide)
        
        Pvalues <- data.frame (loop.for.welch.test(seed.L1.wide))
        subset(Pvalues, p_value <=0.05)
        #View(Pvalues)
        
## 4.0 PCA Construction
        library(FactoMineR)

        seed.pca <- PCA (seed.wide)
        print (seed.pca$var$coord)
        print (seed.pca$ind$coord)        
        
        write.table(seed.pca$var$coord, '/media/mario/DTSE3/SEED_Funct_OK.tsv', sep = '\t', row.names = T)
        write.table(seed.pca$ind$coord, '/media/mario/DTSE3/SEED_Treat_OK.tsv', sep = '\t', row.names = T)
        