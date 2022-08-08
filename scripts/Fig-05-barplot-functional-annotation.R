#######################################################################
## Figure 05 - Functional annotation based on SEED and COG databases ##
#######################################################################
rm (list = ls()); dev.off()
set.seed(999)
library(ggplot2); library (tidyr); library (reshape2); library(Rmisc); library(dplyr)
source('../scripts/loop_for_Welch_test.R.')
metadata <- read.delim('../data/metadata.txt', check.names = F)

		### PLOTTING FIGURE 5A: SEED FUNCTIONAL ANNOTATION BARPLOT ###
		##############################################################

##1. Loading SEED data
	## load the raw dataframe
	## change from wide to long format
	## merge with their metadata
	## filter functions with counts >=33 contigs (46 functions are left)
	data <- read.delim('../data/functionality/SEED_L1_COUNT.txt') %>% 
        	tidyr::pivot_longer(cols = !Function, names_to = "nameid2", values_to = "count") %>% 
        	inner_join(x = metadata[,c(3,7)], by = "nameid2") %>% 
		dplyr::filter(count >= 33)

	## summarizing data before plotting
	data.avg <- Rmisc::summarySE(data = data, measurevar = "count",
				     groupvars = c("Treatment", "Function"), na.rm = T); head(data.avg)

## 2.0 Ploting: BAR PLOT
        SEED_L1 <- 
        	ggplot(data.avg[data.avg$Function!="Plant cell walls and outer surfaces",],
                          aes(x = Function, y = count,fill= Treatment)) +
                geom_errorbar (aes (ymin=count-se, ymax=count+se), 
                               colour="black", width=0.3, position=position_dodge(0.9)) +
                geom_bar(position=position_dodge(), stat = 'identity', width = 0.9) +
                scale_fill_manual(values = c("#FF0000","#0000FF"), # red and blue # 
                                  name="Treatment", 
                                  labels=c("MAS","UNS")) +
                xlab("") + ylab("Number of annotated contigs") +
                theme_bw() + scale_y_continuous(expand = c(0, 0), limits = c(0,8000)) + ## must be up here
                theme (axis.text.x = element_text(angle = 70, hjust = 1)) ## x-axis names angle.
        SEED_L1
        
        SEED_L1 + coord_flip()
        
        ## 2.1. Saving as pdf file:
        # a4r is a landscape plot.
        pdf(file="../figures/Figure-5a-Functional-annotation-barplot-SEEDL1.pdf",
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

       
	        #################################################
	        ## Dando formato a COG-L2 functional annotation #
	        #################################################

## 1. Loading COG Level 2 data
        
        ## load the raw dataframe
        ## change from wide to long format
        ## merge with their metadata
        ## filter functions with counts >=33 contigs (46 functions are left)
        datacog <- read.delim('../data/functionality/COG_COUNT_L2.txt') %>%
        	rename("Function" = 1) %>%
        	tidyr::pivot_longer(cols = !Function, names_to = "nameid2", values_to = "count") %>% 
        	inner_join(x = metadata[,c(3,7)], by = "nameid2")
        
## 2. Summarizing data before plotting: AVERAGE and SD with Rmisc::summarySE()
        datacog.avg <- Rmisc::summarySE(data = datacog, measurevar = "count",
        				 groupvars = c("Function", "Treatment"), na.rm = T); 
        head(datacog.avg)
        
        ## Ploting: averaged BAR PLOT
        COG_L2 <- ggplot(datacog.avg,
        		  aes(x=Function, y=count,fill=Treatment)) +
        	geom_errorbar (aes (ymin=count-se, ymax=count+se),
        		       colour="black", width=0.3, position=position_dodge(0.9)) +
        	geom_bar(position=position_dodge(), stat = 'identity', width = 0.9) +
        	scale_fill_manual(values = c("#FF0000","#0000FF"), # red and blue # 
        			  name="Treatment", 
        			  labels=c("MAS","UNS")) +
        	xlab("") + ylab("Number of annotated contigs") +
        	theme_bw() + scale_y_continuous(expand = c(0, 0), limits = c(0,2000)) + ## must be up here
        	theme (axis.text.x = element_text(angle = 70, hjust = 1)) ## x-axis names angle.
        
        COG_L2
        
        ## 3.0 ¡SAVIN THE PLOT AS A PDF! 
        ### DO NOT RUN TO THE SAME PLACE IF YOU WANT SOMETHING DIFFERENT!
        # a4r is a landscape plot.
        pdf(file = "../figures/Figure-5b-Functional-annotation-barplot-COGL2.pdf",
            width=11, height=8.5, paper='a4r') ## also work with: postscript(..., paper="a4")
        
        COG_L2 + theme( 
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

## save the environment to a file:
	save.image("../environment_Rdata/Fig.05-barplot-functional-annotation.RData")
        