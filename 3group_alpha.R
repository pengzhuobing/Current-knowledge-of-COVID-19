library(ggplot2)
library(ggthemes)
library(ggpubr)
Args <- commandArgs(TRUE)
#aa<-read.table("Healthy_HTN_genus.number.txt.new.xls" ,head=T)
aa<-read.table("Disease_genus.alpha.num.xls",head=T)
aa$Group<-factor(aa$Group,levels = c("HTN","OB","T2D"),order=TRUE)
my_comparisons <- list(c("HTN","OB"),c("HTN","T2D"),c("OB","T2D"))
print=rbind(aggregate(aa$alpha, by=list(aa$Group), FUN=function(x)c(median=median(x),mean=mean(x),max=max(x),min=min(x))))
diff=rbind(compare_means(alpha~Group, data=aa))
#write.table(print,file="Disease_genus.alpha.shannon.log",sep = '\t',quote=F,append=F)
#write.table(diff,file="Disease_genus.alpha.shannon.diff.log",sep = '\t',quote=F,append=F)
pdf("Disease_genus.alpha.num.pdf",width=2.5,height=2.5)
p <- ggplot(aa,aes(x=Group, y=alpha))
p +geom_violin(aes(color=Group),show.legend=F) + 
	geom_boxplot(aes(color=Group),width=0.2,show.legend=F,size=0.4)+
	geom_hline(yintercept = mean(aa$alpha), linetype=2,color="#737373",size=0.3)+
	stat_compare_means(comparisons=my_comparisons,label = "p.signif",label.y = c(150,170,190),method = "wilcox.test",hide.ns=FALSE,size=0.3)+
	theme_wsj()+
	theme(axis.text.x = element_text(color="black",size = 10,angle=45,hjust=1),
		axis.text.y = element_text(color="black",size=6),
		axis.line = element_line(color="black"),
		axis.ticks.x = element_line(color="black"),
		axis.ticks.y = element_line(color="black"),
		panel.background = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.border = element_blank(),
		plot.margin = unit(c(0.3, 0.1, 0.2, 0.2), "in"),
		legend.title = element_text(size=10),
		legend.text = element_text(size=10),
		axis.title=element_text(size=10),
		axis.title.y=element_text(size=8),
		legend.key.size=unit(0.4,"cm")) + 
		labs(y="Number of genus",x="",color="")+
	coord_cartesian(ylim = c(min(aa$alpha),max(aa$alpha)+max(aa$alpha)*0.4))+
	scale_fill_manual(values = c("#d95f02","#a555b5","#025dac"),labels=c("HTN","OB","T2D"))+ 
	scale_color_manual(values = c("#d95f02","#a555b5","#025dac"),labels=c("HTN","OB","T2D")) 
dev.off()
