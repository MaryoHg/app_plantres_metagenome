### Heatmap of the Nitrogen Cycle related genes
###############################################
rm (list = ls())
set.seed(999)
cal_z_score <- function(x){
        (x - mean (x)) / sd (x)
        }

cal_percent <- function (x) {
        (x/sum(x)) * 100
        }

# 1. Loadint the dataframe: counts per each paired end metagenome (n = 8)
        # Loading the raw data counts:
        data <- read.delim('/Users/mariohg/Dropbox/Norma/data/NCyc-genes/nitro_genes.tsv', comment.char = "#", row.names = 1)
        data_t <- data.frame(t(data))
        
        meta_table <- data.frame(row.names=rownames(data_t), 
                                 t(as.data.frame(strsplit(rownames(data_t),"_"))))
        
        data_N <- aggregate(data_t, by = list(meta_table$X2, meta_table$X3), mean)
        
        
        # Nomrmalizing the data according to:
        # Determining the z-score centered dataframe by gene (colum):
        data_z <- data.frame(t(apply(X = data_N[,-1:-2], MARGIN = 2, FUN = cal_z_score)))
        
        # Total percent by row (SampleID)
        data_perc <- data.frame(apply(X = data_N[,-1:-2], MARGIN = 1, FUN = cal_percent))
        colSums(data_perc)
        
        
# 2. Annotation metadata for heatmap
        
        # Colums (samples) data
        col_anot <- data.frame ("Treatment" = colnames(data_z),
                                "Metagenome" = rep (c("UNS", "MAS"), times =2), row.names = 1); col_anot
        # Row annotation data
        annotation <- read.delim('/Users/mariohg/Dropbox/Norma/data/NCyc-genes/annot_genes.tsv', row.names = 2)
        annotation$Annotation <- NULL

        pathways <- RColorBrewer::brewer.pal(n = 9, name = "Set3"); pathways
        names(pathways) <- levels(factor(annotation$Pathways)); pathways
        
        # annotations colors:
        heatmap_colors <- list (Metagenome = c("UNS" = "blue", "MAS" = "red"),
                                Pathways = pathways); heatmap_colors
        
# 3. Drawing the heatmap: z-score heatmap
        # italic gene names:
        genes <- lapply(rownames(data_z), function (x) bquote(italic(.(x))))
        
        # breaks of heatmap legend
        breaks = base::seq(from = -2, to = 2, by = 0.01)
        
        pheatmap::pheatmap(mat = data_z, labels_row = as.expression(genes),
                           cellwidth = 30, cellheight = 9,
                           annotation_col = col_anot,
                           annotation_row = annotation,
                           annotation_colors = heatmap_colors,
                           show_colnames = F, treeheight_col = 25, treeheight_row = 25,
                           cluster_rows = T, cluster_cols = T, cutree_cols = 2,
                           color = colorRampPalette(c("white", "yellow", "red"))(length(breaks)),
                           breaks = breaks, clustering_method = "average")
        
# 4. Drawing heatmap: relative abundance heatmap
        # Adding the defines factor to plot the colors with:
        colnames(data_perc)
        data_perc <- within(data_perc, {
                UNS1 <- 0
                UNS1 [Control1 <=0.01]<- 0
                UNS1 [Control1 > 0.01 & Control1 <= 0.05]<- 1
                UNS1 [Control1 > 0.05 & Control1 <= 0.10]<- 2
                UNS1 [Control1 > 0.10 & Control1 <= 0.20]<- 3
                UNS1 [Control1 > 0.20 & Control1  <= 0.50]<- 4
                UNS1 [Control1 > 0.50 & Control1 <= 1.00]<- 5
                UNS1 [Control1 > 1.00 & Control1 <= 2.00]<- 6
                UNS1 [Control1 > 2.00 & Control1 <= 4.00]<- 7
                UNS1 [Control1 > 4.00 & Control1 <= 8.00]<- 8
                UNS1 [Control1 > 8.00 & Control1 <= 12.00]<- 9
                UNS1 [Control1 > 12.00 & Control1 <= 16.00]<- 10
                UNS1 [Control1 > 16.00 & Control1 <= 20.00]<- 11
                UNS1 [Control1 > 20.00]<- 12})
        
        data_perc <- within(data_perc, {
                UNS2 <- 0
                UNS2 [Control2 <=0.01]<- 0
                UNS2 [Control2 > 0.01 & Control2 <= 0.05]<- 1
                UNS2 [Control2 > 0.05 & Control2 <= 0.10]<- 2
                UNS2 [Control2 > 0.10 & Control2 <= 0.20]<- 3
                UNS2 [Control2 > 0.20 & Control2  <= 0.50]<- 4
                UNS2 [Control2 > 0.50 & Control2 <= 1.00]<- 5
                UNS2 [Control2 > 1.00 & Control2 <= 2.00]<- 6
                UNS2 [Control2 > 2.00 & Control2 <= 4.00]<- 7
                UNS2 [Control2 > 4.00 & Control2 <= 8.00]<- 8
                UNS2 [Control2 > 8.00 & Control2 <= 12.00]<- 9
                UNS2 [Control2 > 12.00 & Control2 <= 16.00]<- 10
                UNS2 [Control2 > 16.00 & Control2 <= 20.00]<- 11
                UNS2 [Control2 > 20.00]<- 12})
        
        
        data_perc <- within(data_perc, {
                MAS1 <- 0
                MAS1 [Maiz1 <=0.01]<- 0
                MAS1 [Maiz1 > 0.01 & Maiz1 <= 0.05]<- 1
                MAS1 [Maiz1 > 0.05 & Maiz1 <= 0.10]<- 2
                MAS1 [Maiz1 > 0.10 & Maiz1 <= 0.20]<- 3
                MAS1 [Maiz1 > 0.20 & Maiz1  <= 0.50]<- 4
                MAS1 [Maiz1 > 0.50 & Maiz1 <= 1.00]<- 5
                MAS1 [Maiz1 > 1.00 & Maiz1 <= 2.00]<- 6
                MAS1 [Maiz1 > 2.00 & Maiz1 <= 4.00]<- 7
                MAS1 [Maiz1 > 4.00 & Maiz1 <= 8.00]<- 8
                MAS1 [Maiz1 > 8.00 & Maiz1 <= 12.00]<- 9
                MAS1 [Maiz1 > 12.00 & Maiz1 <= 16.00]<- 10
                MAS1 [Maiz1 > 16.00 & Maiz1 <= 20.00]<- 11
                MAS1 [Maiz1 > 20.00]<- 12})
        
        data_perc <- within(data_perc, {
                MAS2 <- 0
                MAS2 [Maiz2 <=0.01]<- 0
                MAS2 [Maiz2 > 0.01 & Maiz2 <= 0.05]<- 1
                MAS2 [Maiz2 > 0.05 & Maiz2 <= 0.10]<- 2
                MAS2 [Maiz2 > 0.10 & Maiz2 <= 0.20]<- 3
                MAS2 [Maiz2 > 0.20 & Maiz2  <= 0.50]<- 4
                MAS2 [Maiz2 > 0.50 & Maiz2 <= 1.00]<- 5
                MAS2 [Maiz2 > 1.00 & Maiz2 <= 2.00]<- 6
                MAS2 [Maiz2 > 2.00 & Maiz2 <= 4.00]<- 7
                MAS2 [Maiz2 > 4.00 & Maiz2 <= 8.00]<- 8
                MAS2 [Maiz2 > 8.00 & Maiz2 <= 12.00]<- 9
                MAS2 [Maiz2 > 12.00 & Maiz2 <= 16.00]<- 10
                MAS2 [Maiz2 > 16.00 & Maiz2 <= 20.00]<- 11
                MAS2 [Maiz2 > 20.00]<- 12})
        
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
