library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggthemes)
aa<-read.table("HTN-OB-T2D_genus.alpha.num.xls" ,head=T)
pdf("HTN-OB-T2D_genus.alpha.num.pdf",width=2.5,height=2.5)
aa$Group<-factor(aa$Group,levels = c("Healthy-HTN","Healthy-OB","Healthy-T2D"),order=TRUE)
aa$type<-factor(aa$type,levels = c("Healthy","Disease"),order=TRUE)
dodge <- position_dodge(width = 1)
#colour <- c("#1b9e77", "#f781bf")
p2 <-ggplot(data = aa, aes(x = Group, y = alpha)) +
	#geom_violin(aes(color=colour),fill=colour,position = dodge,show.legend=F)+
	geom_violin(aes(color=type),position = dodge,show.legend=F)+
	#geom_violin(position = dodge,show.legend=F)+
	geom_boxplot(aes(color=type),width=0.2,outlier.colour=NA, position = dodge,size=0.4,show.legend=T) +
	#geom_boxplot(width=0.2,outlier.colour=NA, position = dodge,size=0.4,show.legend=T) +
	#geom_boxplot(color = "black",width=0.2,outlier.colour=NA, position = dodge,size=0.4,show.legend=T) +
	stat_compare_means(aes(group=type),label = "p.signif",label.y = 150,hide.ns=TRUE,method="wilcox.test")+
	scale_color_manual(values=c("#1b9e77", "#d95f02"),labels=c("Healthy","Disease"))+
	theme_wsj()
p2+theme(axis.text.x = element_text(color="black",angle=45,hjust=1,size=10),
	axis.text.y = element_text(color="black",size=8),
	axis.line = element_line(color="black"),
	axis.ticks.x = element_line(color="black"),
	axis.ticks.y = element_line(color="black"),
	panel.background = element_blank(),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	panel.border = element_blank(),
	plot.margin = unit(c(0.05, 0.1, 0.05, 0.1), "in"),
	legend.text = element_text(size=6),
	axis.title=element_text(size=10),
	axis.title.y=element_text(size=10),
	legend.position = "top",
	legend.key.size=unit(0.4,"cm")) + 
	labs(y="Number of genus",color="",x="")+
coord_cartesian(ylim = c(50,170))
dev.off()
