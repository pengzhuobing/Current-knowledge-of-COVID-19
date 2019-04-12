##绘制 zscore 的热图，并且用星星表示显著性差异，用三角形表示丰度变化趋势，涉及的文件放到同级的 Data 的目录下；
library('ComplexHeatmap',lib.loc  = "/home/pengzhuobing/R/x86_64-pc-linux-gnu-library/3.4") #这个包需要下载github最新版本
library(png)
#将丰度文件zscore标准化
profile <- read.table("add_HC_CHF_AHF_genus_qvalue_0.05.xls.profile", header=TRUE, sep="\t", comment.char = "", check.names = F, row.names=1)
zscore <- t(profile)
zscore <- scale(zscore, center = T, scale = TRUE)
write.table(zscore,"genus_z_score.txt",quote = F,sep = "\t",col.names = NA) #设置col.names = NA将第一行第一个设为\t；
dt = read.table("genus_z_score.txt", header=TRUE, sep="\t", comment.char = "", check.names = F, row.names=1)
#读取注释文件，这个文件也是需要提前处理好，用来表示显著性差异的
temp = read.table("HC_CHF_AHF_genus_annotation_qvalue.txt", sep = "\t", header = TRUE, row.names = 1, comment.char = "", check.names = F)
#这个文件也是需要提前处理好，代表pch的数字代码，用来表示丰度趋势的
trend <- read.table("HC_CHF_AHF_genus_annotation_profile1.txt", sep = "\t", header = TRUE, row.names = 1, comment.char = "", check.names = F)

#首先画热图，需要得到聚类的顺序
p = Heatmap(dt,rect_gp = gpar(col = "grey50", lwd = 1), name = "", heatmap_legend_param = list(border = "grey60"), cluster_rows = F, cluster_columns = T,show_row_names = T)
temp_anno = temp[colnames(dt)[column_order(p)],]
trend_anno <- trend[colnames(dt)[column_order(p)],]
aa <- dt[,colnames(dt)[column_order(p)]]
a <- as.matrix(aa)
#注释文件的图片路径，这个脚本实际上是把星星的图片，绘制在了热图下边，所以需要提供图片文件，在 Data 目录下。类推的话，如果想画其他图片也是可以的，
p_png = "data/p.png"
lessp_png = "data/lessp.png"
no_sig_png = "data/nosig.png"
#注释丰度趋势
#down_down = "data/down_down.png"
#up_up = "data/up_up.png"
#up_down = "data/up_down.png"
#down_up = "data/down_up.png"

#将注释文件的字符串转变为变量
HC_vs_CHF = as.character(lapply(as.character(temp_anno$`HC vs CHF`), function(x) eval(parse(text = x))))
HC_vs_AHF = as.character(lapply(as.character(temp_anno$`HC vs AHF`), function(x) eval(parse(text = x))))
CHF_vs_AHF = as.character(lapply(as.character(temp_anno$`CHF vs AHF`), function(x) eval(parse(text = x))))

#生成热图的注释
value = rep(0.5,22)
ha = HeatmapAnnotation(
#这个会有底色   "HC vs CHF vs AHF" = anno_simple(value, pch = trend_anno,pt_gp = gpar(fill = "black")),
 	"HC vs CHF vs AHF" = anno_points(value,axis=F,pch = trend_anno,size = unit(5, "mm"),axis_param = NULL,border=FALSE,show_annotation_name = FALSE,ylim=c(0,1),gp = gpar(fill = "black")),
	"HC vs CHF" = anno_image(HC_vs_CHF,border = F),
  "HC vs AHF" = anno_image(HC_vs_AHF, border = F),
	"CHF vs AHF" = anno_image(CHF_vs_AHF,border = F),
	annotation_name_gp = gpar(fontsize = 20),
	  gap = unit(c(0.3,0.15,0.15), "cm")
    ## 上边是门水平的脚本，如果画属或者种的话，需要调节图片的参数
    #"NCA vs sCAD" = anno_image(NCA_vs_sCAD, border = F, height = unit(4,"mm"), space = unit(4, "mm")),  #需要通过调节参数控制图片大小
    #"NCA vs AMI" = anno_image(NCA_vs_AMI,border = F,height = unit(4,"mm"),space = unit(4, "mm")), 
    #"sCAD vs AMI" = anno_image(sCAD_vs_AMI,border = F,height = unit(4,"mm"), space = unit(4, "mm"))
)


pdf("genus_heatmap.pdf", height = 8, width =25)
Heatmap(a, color = c("navy", "white", "firebrick3"), rect_gp = gpar(col = "grey50", lwd = 1), name = "", heatmap_legend_param = list(border = "grey50"), cluster_rows = F, cluster_columns = T,show_row_names = T, bottom_annotation = ha, width = unit(2*ncol(dt), "cm"), height = unit(2*nrow(dt), "cm"), row_names_gp = gpar(fontsize = 20), column_names_gp = gpar(fontsize = 20))
#Heatmap(dt,row_dend_width = unit(2, "cm"),color = c("navy", "white", "firebrick3"), rect_gp = gpar(col = "grey50", lwd = 1), name = "", heatmap_legend_param = list(border = "grey50"), cluster_rows= F, cluster_columns = T,show_row_names = T, bottom_annotation = ha, width = unit(2*ncol(dt), "cm"), height = unit(2*nrow(dt), "cm"), row_names_gp = gpar(fontsize = 20), column_names_gp = gpar(fontsize = 20))
dev.off()
