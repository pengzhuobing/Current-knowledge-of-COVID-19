library ("ComplexHeatmap")
library(circlize)
library(dtwclust)


Args <- commandArgs(TRUE)

aG <- read.table(Args[2],header=T,check.names=F)
aGroup <- log2(as.matrix(scale(aG, center=FALSE, scale=colSums(aG))))
ZaGroup <- zscore(as.matrix(scale(aG, center=FALSE, scale=colSums(aG))))
bG <- read.table(Args[3],header=T,check.names=F)
bGroup <- log2(as.matrix(scale(bG, center=FALSE, scale=colSums(bG))))
ZbGroup <- zscore(as.matrix(scale(bG, center=FALSE, scale=colSums(bG))))

m1 <- aGroup
m2 <- bGroup

corr <- read.table(Args[1],header=T,check.names=F)
colors = structure(c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462"),
           names = c("C1", "C2", "C3","C4","L1","L2"))

nr = nrow(m1)
rg = range(c(m1, m2))
rg[1] = rg[1] - (rg[2] - rg[1])* 0.02
rg[2] = rg[2] + (rg[2] - rg[1])* 0.02
anno_multiple_boxplot = function(index) {
    pushViewport(viewport(xscale = rg, yscale = c(0.5, nr + 0.5)))
    for(i in seq_along(index)) {
        grid.rect(y = nr-i+1, height = 1, default.units = "native")
        grid.boxplot(m1[ index[i], ], pos = nr-i+1 + 0.2, box_width = 0.3, 
            gp = gpar(fill = "#8BB1D3"), direction = "horizontal",outline = FALSE,pch = 20, size = unit(1, "mm"))
        grid.boxplot(m2[ index[i], ], pos = nr-i+1 - 0.2, box_width = 0.3, 
            gp = gpar(fill = "#FFC080"), direction = "horizontal",pch = 20, size = unit(1, "mm"))
    }
    grid.xaxis()
    popViewport()
}

boxplotAnno = rowAnnotation(boxplot = anno_multiple_boxplot, width = unit(5, "cm"),show_annotation_name = FALSE,border = FALSE,gp = gpar(col = "red"))


ht_list = Heatmap(as.matrix(corr),col=colors,heatmap_legend_param = list(title = "class",labels_gp=gpar(fontsize=8)),show_column_names = FALSE,show_row_names = TRUE,width = unit(0.2*ncol(corr), "cm"), height = unit(0.4*nrow(corr), "cm"),left_annotation = boxplotAnno,row_names_gp = gpar(fontsize = 8),column_names_gp = gpar(fontsize = 8),row_names_side = "left", show_row_dend = FALSE)

#lgd = Legend(title = "boxplots",legend_gp = gpar(fill = c("#8BB1D3","#FFC080")))

pdf(Args[4],width=10,height=10)
draw(ht_list, padding = unit(c(20, 2, 2, 2), "mm"))
dev.off()

