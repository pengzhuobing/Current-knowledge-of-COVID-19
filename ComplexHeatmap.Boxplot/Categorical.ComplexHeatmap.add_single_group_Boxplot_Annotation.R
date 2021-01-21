library ("ComplexHeatmap")
library(circlize)
Args <- commandArgs(TRUE)
corr <- read.table(Args[1],header=T,check.names=F)

colors = structure(c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462"),
           names = c("C1", "C2", "C3","C4","L1","L2"))

aG <- read.table(Args[2],header=T,check.names=F)
aGroup <- log2(as.matrix(scale(aG, center=FALSE, scale=colSums(aG))))
bG <- read.table(Args[3],header=T,check.names=F)
bGroup <-log2(as.matrix(scale(bG, center=FALSE, scale=colSums(bG))))

row_ha = rowAnnotation(A = anno_boxplot(aGroup, height = unit(4, "cm"), width = unit(5, "cm"), box_width = 0.5, outline = FALSE,gp = gpar(fill = "#8BB1D3")),
					   B = anno_boxplot(bGroup, height = unit(4, "cm"), width = unit(5, "cm"), ,box_width = 0.5, outline = FALSE,gp = gpar(fill = "#FFC080"))
)


pdf(Args[4],width=10,height=10)
sample_order <- Heatmap(as.matrix(corr),col=colors,
                        heatmap_legend_param = list(title = "class",labels_gp=gpar(fontsize=8)),
						#show_row_names = FALSE,
						show_column_names = FALSE,
						row_names_side = "left", show_row_dend = FALSE,
                        row_names_gp = gpar(fontsize = 8),
                        column_names_gp = gpar(fontsize = 8),
                        width = unit(0.15*ncol(corr), "cm"), height = unit(0.25*nrow(corr), "cm"),
                        left_annotation = row_ha
						)
sample_order
dev.off()
