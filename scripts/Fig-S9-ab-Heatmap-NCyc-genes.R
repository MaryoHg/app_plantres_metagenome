### Figure S9-A,B: Heatmap of the Nitrogen Cycle related genes
##############################################################
options(scipen = 100000, digits = 3); rm (list = ls())
set.seed(999)
	zscorefunct <- function(x){
	        (x - mean (x)) / sd (x)
	        }
	
	relabpercfunt <- function (x) {
	        (x/sum(x)) * 100}

metadata <- read.table('../data/metadata.txt',  header = T, stringsAsFactors = T, check.names = F)

# 1. Loading the data-frame: counts per each paired end metagenome (n = 8)
        data <- read.delim('../data/functionality/nitro_genes.tsv', 
        		   comment.char = "#", row.names = 1)
        colnames(data) <- as.character(gsub(pattern = "X", replacement = "", colnames(data)))
         
        # summarizing data by replicate and treatment
        data_t <- data.frame(t(data))
        meta_table <- data.frame(row.names=rownames(data_t), 
                                 t(as.data.frame(strsplit(rownames(data_t),"_"))))
        data_N <- stats::aggregate(data_t, by = list(meta_table$X2, meta_table$X1), FUN = mean)
        rownames(data_N) <- paste(data_N$Group.1, data_N$Group.2, sep = "_")
        #data_N1 <- aggregate(data_t, by = list(meta_table$X2, meta_table$X3), mean)
        
## 2. Nomamalizing the data according to: 
        ## z-score (by gene)
        data_z <- data.frame(t(apply(X = data_N[,-1:-2], MARGIN = 2, FUN = zscorefunct)))
        data_z <- data_z %>%
        	dplyr::filter(!is.na(Control_15))
        
        ## relative abundance per sample (SampleID)
        data_perc <- data.frame(apply(X = data_N[,-1:-2], MARGIN = 1, FUN = relabpercfunt))
        colSums(data_perc)
        
        
## 3. Annotation metadata for heatmap
        
        ## metadata for samples and their respective treatment
        col_anot <- data.frame ("Treatment" = colnames(data_z),
                                "Metagenome" = rep (c("UNS", "MAS"), 
                                		    times =2), row.names = 1); col_anot
        
        # Row annotation data: 
        ## load annotation for N pathways classification of genes
        ## set colors for each pathway
        annotation <- read.delim('../data/functionality/annot_genes.tsv', row.names = 2)
        annotation$Annotation <- NULL

        pathways <- RColorBrewer::brewer.pal(n = 9, name = "Set3"); pathways
        names(pathways) <- levels(factor(annotation$Pathways)); pathways
        
        # annotations colors:
        heatmap_colors <- list (Metagenome = c("UNS" = "blue", "MAS" = "red"),
                                Pathways = pathways); heatmap_colors
        
## 4. Drawing the heatmap: z-score heatmap
        # italic gene names:
        genes <- lapply(rownames(data_z), function (x) bquote(italic(.(x))))
        
        # breaks fpr heatmap legend
        breaks = base::seq(from = -2, to = 2, by = 0.01)
        
        pheatmap::pheatmap(mat = data_z, 
        		   labels_row = as.expression(genes),
                           cellwidth = 30, cellheight = 9,
                           annotation_col = col_anot,
                           annotation_row = annotation,
                           annotation_colors = heatmap_colors,
                           show_colnames = F, 
        		   treeheight_col = 25, treeheight_row = 25,
                           cluster_rows = T, cluster_cols = T, 
        		   cutree_cols = 2,
                           color = colorRampPalette(c("white", "yellow", "red"))(length(breaks)),
                           breaks = breaks, clustering_method = "average")
        
# 4. Drawing heatmap: relative abundance heatmap
        # Adding the defines factor to plot the colors with:
        colnames(data_perc)
        data_perc <- within(data_perc, {
                UNS1 <- 0
                UNS1 [Control_15 <=0.01]<- 0
                UNS1 [Control_15 > 0.01 & Control_15 <= 0.05]<- 1
                UNS1 [Control_15 > 0.05 & Control_15 <= 0.10]<- 2
                UNS1 [Control_15 > 0.10 & Control_15 <= 0.20]<- 3
                UNS1 [Control_15 > 0.20 & Control_15  <= 0.50]<- 4
                UNS1 [Control_15 > 0.50 & Control_15 <= 1.00]<- 5
                UNS1 [Control_15 > 1.00 & Control_15 <= 2.00]<- 6
                UNS1 [Control_15 > 2.00 & Control_15 <= 4.00]<- 7
                UNS1 [Control_15 > 4.00 & Control_15 <= 8.00]<- 8
                UNS1 [Control_15 > 8.00 & Control_15 <= 12.00]<- 9
                UNS1 [Control_15 > 12.00 & Control_15 <= 16.00]<- 10
                UNS1 [Control_15 > 16.00 & Control_15 <= 20.00]<- 11
                UNS1 [Control_15 > 20.00]<- 12})
        
        data_perc <- within(data_perc, {
                UNS2 <- 0
                UNS2 [Control_25 <=0.01]<- 0
                UNS2 [Control_25 > 0.01 & Control_25 <= 0.05]<- 1
                UNS2 [Control_25 > 0.05 & Control_25 <= 0.10]<- 2
                UNS2 [Control_25 > 0.10 & Control_25 <= 0.20]<- 3
                UNS2 [Control_25 > 0.20 & Control_25  <= 0.50]<- 4
                UNS2 [Control_25 > 0.50 & Control_25 <= 1.00]<- 5
                UNS2 [Control_25 > 1.00 & Control_25 <= 2.00]<- 6
                UNS2 [Control_25 > 2.00 & Control_25 <= 4.00]<- 7
                UNS2 [Control_25 > 4.00 & Control_25 <= 8.00]<- 8
                UNS2 [Control_25 > 8.00 & Control_25 <= 12.00]<- 9
                UNS2 [Control_25 > 12.00 & Control_25 <= 16.00]<- 10
                UNS2 [Control_25 > 16.00 & Control_25 <= 20.00]<- 11
                UNS2 [Control_25 > 20.00]<- 12})
        
        
        data_perc <- within(data_perc, {
                MAS1 <- 0
                MAS1 [Maiz_25 <=0.01]<- 0
                MAS1 [Maiz_25 > 0.01 & Maiz_25 <= 0.05]<- 1
                MAS1 [Maiz_25 > 0.05 & Maiz_25 <= 0.10]<- 2
                MAS1 [Maiz_25 > 0.10 & Maiz_25 <= 0.20]<- 3
                MAS1 [Maiz_25 > 0.20 & Maiz_25  <= 0.50]<- 4
                MAS1 [Maiz_25 > 0.50 & Maiz_25 <= 1.00]<- 5
                MAS1 [Maiz_25 > 1.00 & Maiz_25 <= 2.00]<- 6
                MAS1 [Maiz_25 > 2.00 & Maiz_25 <= 4.00]<- 7
                MAS1 [Maiz_25 > 4.00 & Maiz_25 <= 8.00]<- 8
                MAS1 [Maiz_25 > 8.00 & Maiz_25 <= 12.00]<- 9
                MAS1 [Maiz_25 > 12.00 & Maiz_25 <= 16.00]<- 10
                MAS1 [Maiz_25 > 16.00 & Maiz_25 <= 20.00]<- 11
                MAS1 [Maiz_25 > 20.00]<- 12})
        
        data_perc <- within(data_perc, {
                MAS2 <- 0
                MAS2 [Maiz_35 <=0.01]<- 0
                MAS2 [Maiz_35 > 0.01 & Maiz_35 <= 0.05]<- 1
                MAS2 [Maiz_35 > 0.05 & Maiz_35 <= 0.10]<- 2
                MAS2 [Maiz_35 > 0.10 & Maiz_35 <= 0.20]<- 3
                MAS2 [Maiz_35 > 0.20 & Maiz_35  <= 0.50]<- 4
                MAS2 [Maiz_35 > 0.50 & Maiz_35 <= 1.00]<- 5
                MAS2 [Maiz_35 > 1.00 & Maiz_35 <= 2.00]<- 6
                MAS2 [Maiz_35 > 2.00 & Maiz_35 <= 4.00]<- 7
                MAS2 [Maiz_35 > 4.00 & Maiz_35 <= 8.00]<- 8
                MAS2 [Maiz_35 > 8.00 & Maiz_35 <= 12.00]<- 9
                MAS2 [Maiz_35 > 12.00 & Maiz_35 <= 16.00]<- 10
                MAS2 [Maiz_35 > 16.00 & Maiz_35 <= 20.00]<- 11
                MAS2 [Maiz_35 > 20.00]<- 12})
        
        # annotation:
        # Colums (samples) data percent
        col_anot2 <- data.frame ("Treatment" = colnames(data_perc[,5:8]),
                                "Metagenome" = rep (c("UNS", "MAS"), each =2), row.names = 1); col_anot2
        # breaks
        breaks2 <- base::seq(from = 0, to = 12, by = 0.01)
        
        # heatmap
        pheatmap::pheatmap(mat = data_perc[,5:8], labels_row = as.expression(genes),
                           cellwidth = 30, cellheight = 9,
                           annotation_col = col_anot2,
                           annotation_row = annotation,
                           annotation_colors = heatmap_colors,
                           show_colnames = F, treeheight_col = 25, treeheight_row = 25,
                           cluster_rows = T, cluster_cols = T, cutree_cols = 2,
                           color = colorRampPalette(c("white", "yellow", "red"))(length(breaks2)),
                           breaks = breaks2, clustering_method = "average",
                           legend_breaks = 0:12, legend_labels = c('<=0.01',
                                                                   '> 0.01 & <= 0.05',
                                                                   '> 0.05 & <= 0.10',
                                                                   '> 0.10 & <= 0.20',
                                                                   '> 0.20 &  <= 0.50',
                                                                   '> 0.50 & <= 1.00',
                                                                   '> 1.00 & <= 2.00',
                                                                   '> 2.00 & <= 4.00',
                                                                   '> 4.00 & <= 8.00',
                                                                   '> 8.00 & <= 12.00',
                                                                   '> 12.00 & <= 16.00',
                                                                   '> 16.00 & <= 20.00',
                                                                   '> 20.00'))
