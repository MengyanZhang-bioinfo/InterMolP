#BRCA####
 
clusters<-as.data.frame(clusters)
clusters_tide<-data.frame(Patient=rownames(clusters),group=clusters$clusters)

#1.clustertide
 
tide_final<-left_join(tide,clusters_tide,by="Patient")
T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
p_tide_zanshi<-T_result$p.value 
library(ggplot2)
library(gghalves)
library(ggprism)   
library(ggpubr)    
tide_plot<-data.frame(Tide=as.numeric(tide_final$TIDE),Cluster=as.character(tide_final$group))
p<-ggplot(tide_plot,aes(x=Cluster,y=Tide,fill=Cluster,color=Cluster))+
  scale_fill_manual(values=c("#80D1C8","#FFD4A9"))+
  scale_colour_manual(values=c("#80D1C8","#FFD4A9"))
p
 
p1 <- p + geom_half_violin(
  position = position_nudge(x = 0.1, y = 0),
  side = 'R',
  adjust = 1.2,
  trim = FALSE,
  color = NA,
  alpha = 0.6
)
p1

p2 <- p1 + geom_half_point(
  position = position_nudge(x = -0.35, y = 0),
  size = 3,
  shape = 19,
  range_scale = 0.5,
  alpha = 0.6
)
p2 
 
p3 <- p2 + geom_boxplot(
  outlier.shape = NA,   
  width = 0.1,
  alpha = 0.8
)
p3
 
p4 <- p3 + stat_compare_means(
 
  method = "wilcox.test",   
  label = "p.format",   
  label.x = 1.5,   
  size = 5   
)
p4
 
# windowsFonts(A = windowsFont("TimesNewRoman"),
#              B = windowsFont("Arial"))

p5 <- p4 + theme_light() + theme(panel.grid = element_blank()) +
  theme(
    axis.ticks.length = unit(-0.25, "cm"),
    axis.text.x = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    axis.text.y = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    panel.grid.major = element_line(linetype = "dotted", color = "lightgray", size = 0.5),  #  
    text = element_text(size = 20)
  )
p5

p5 <- p4 + theme_light() + theme(panel.grid = element_blank()) +
  theme(
    axis.ticks.length = unit(-0.25, "cm"),
    axis.text.x = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    axis.text.y = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    panel.grid.major = element_line(linetype = "dotted", color = "lightgray", size = 0.5),
    text = element_text(size = 20)
  ) +
  coord_cartesian(ylim = c(-4, 4))  #  


print(p5)

 
pdf("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\BRCA_cluster_tide_plot1.pdf",width = 7,height = 7)
p5
dev.off()
#ggsave("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\BRCA_cluster_tide_plot.pdf",p5,width = 10,height = 10)

 
label<-read.csv("G:\\pan_result\\gpu_result_tcga\\BRCA\\label.csv")

clusters_tide<-data.frame(Patient=label$X,group=label$label)
library(tidyverse)
 
tide_final<-left_join(tide,clusters_tide,by="Patient")
T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
p_tide_zanshi<-T_result$p.value 
library(ggplot2)
library(gghalves)
library(ggprism)   
library(ggpubr)    
tide_plot<-data.frame(Tide=as.numeric(tide_final$TIDE),Stage=as.character(tide_final$group))
p<-ggplot(tide_plot,aes(x=Stage,y=Tide,fill=Stage,color=Stage))+
  scale_fill_manual(values=c("#D55E00","#56B4E9"))+
  scale_colour_manual(values=c("#D55E00","#56B4E9"))

p1 <- p + geom_half_violin(
  position = position_nudge(x = 0.1, y = 0),
  side = 'R',
  adjust = 1.2,
  trim = FALSE,
  color = NA,
  alpha = 0.6
)

 
p2 <- p1 + geom_half_point(
  position = position_nudge(x = -0.35, y = 0),
  size = 3,
  shape = 19,
  range_scale = 0.5,
  alpha = 0.6
)

 
p3 <- p2 + geom_boxplot(
  outlier.shape = NA,   
  width = 0.1,
  alpha = 0.8
)

 
p4 <- p3 + stat_compare_means(

  method = "wilcox.test",  
  label = "p.format",   
  label.x = 1.5,   
  size = 5   
)

 
windowsFonts(A = windowsFont("TimesNewRoman"),
             B = windowsFont("Arial"))

p5 <- p4 + theme_light() + theme(panel.grid = element_blank()) +
  theme(
    axis.ticks.length = unit(-0.25, "cm"),
    axis.text.x = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    axis.text.y = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    panel.grid.major = element_line(linetype = "dotted", color = "lightgray", size = 0.5),  # 
    text = element_text(size = 20)
  )


print(p5)
 

ggsave("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\BRCA_satge_tide_plot.pdf",p5,width = 10,height = 10)


#COAD####

 
clusters<-as.data.frame(clusters)
clusters_tide<-data.frame(Patient=rownames(clusters),group=clusters$clusters)

#1.clustertide
 
tide_final<-left_join(tide,clusters_tide,by="Patient")
T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
p_tide_zanshi<-T_result$p.value 
library(ggplot2)
library(gghalves)
library(ggprism)   
library(ggpubr)    
tide_plot<-data.frame(Tide=as.numeric(tide_final$TIDE),Cluster=as.character(tide_final$group))
p<-ggplot(tide_plot,aes(x=Cluster,y=Tide,fill=Cluster,color=Cluster))+
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8"))+
  scale_colour_manual(values=c("#FFD4A9", "#80D1C8"))

p1 <- p + geom_half_violin(
  position = position_nudge(x = 0.1, y = 0),
  side = 'R',
  adjust = 1.2,
  trim = FALSE,
  color = NA,
  alpha = 0.6
)

 
p2 <- p1 + geom_half_point(
  position = position_nudge(x = -0.35, y = 0),
  size = 3,
  shape = 19,
  range_scale = 0.5,
  alpha = 0.6
)

 
p3 <- p2 + geom_boxplot(
  outlier.shape = NA,   
  width = 0.1,
  alpha = 0.8
)

 
p4 <- p3 + stat_compare_means(
 
  method = "wilcox.test",   
  label = "p.format",   
  label.x = 1.5,   
  size = 5   
)

 
windowsFonts(A = windowsFont("TimesNewRoman"),
             B = windowsFont("Arial"))

p5 <- p4 + theme_light() + theme(panel.grid = element_blank()) +
  theme(
    axis.ticks.length = unit(-0.25, "cm"),
    axis.text.x = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    axis.text.y = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    panel.grid.major = element_line(linetype = "dotted", color = "lightgray", size = 0.5),  #  
    text = element_text(size = 20)
  )+
  coord_cartesian(ylim = c(-4, 4))  # 


 
print(p5)
 

ggsave("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\COAD_cluster_tide_plot1.pdf",p5,width = 10,height = 10)


 
label<-read.csv("label.csv")

clusters_tide<-data.frame(Patient=label$X,group=label$label)
library(tidyverse)
 
tide_final<-left_join(tide,clusters_tide,by="Patient")
T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
p_tide_zanshi<-T_result$p.value 
library(ggplot2)
library(ggprism)  
library(ggpubr)    
tide_plot<-data.frame(Tide=as.numeric(tide_final$TIDE),Stage=as.character(tide_final$group))
p<-ggplot(tide_plot,aes(x=Stage,y=Tide,fill=Stage,color=Stage))+
  scale_fill_manual(values=c("#D55E00","#56B4E9"))+
  scale_colour_manual(values=c("#D55E00","#56B4E9"))
 
p1 <- p + geom_half_violin(
  position = position_nudge(x = 0.1, y = 0),
  side = 'R',
  adjust = 1.2,
  trim = FALSE,
  color = NA,
  alpha = 0.6
)

 
p2 <- p1 + geom_half_point(
  position = position_nudge(x = -0.35, y = 0),
  size = 3,
  shape = 19,
  range_scale = 0.5,
  alpha = 0.6
)

 
p3 <- p2 + geom_boxplot(
  outlier.shape = NA,   
  width = 0.1,
  alpha = 0.8
)


p4 <- p3 + stat_compare_means(
 
  method = "wilcox.test",  
  label = "p.format",   
  label.x = 1.5,   
  size = 5   
)


windowsFonts(A = windowsFont("TimesNewRoman"),
             B = windowsFont("Arial"))

p5 <- p4 + theme_light() + theme(panel.grid = element_blank()) +
  theme(
    axis.ticks.length = unit(-0.25, "cm"),
    axis.text.x = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    axis.text.y = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    panel.grid.major = element_line(linetype = "dotted", color = "lightgray", size = 0.5),  #  
    text = element_text(size = 20)
  )

 
print(p5)
 

ggsave("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\COAD_satge_tide_plot.pdf",p5,width = 10,height = 10)



#ESCA####

 
clusters<-as.data.frame(clusters)
clusters_tide<-data.frame(Patient=rownames(clusters),group=clusters$clusters)

#1.clustertide
 
tide_final<-left_join(tide,clusters_tide,by="Patient")
T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
p_tide_zanshi<-T_result$p.value 
library(ggplot2)
library(gghalves)
library(ggprism)  
library(ggpubr)    
tide_plot<-data.frame(Tide=as.numeric(tide_final$TIDE),Cluster=as.character(tide_final$group))
p<-ggplot(tide_plot,aes(x=Cluster,y=Tide,fill=Cluster,color=Cluster))+
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8"))+
  scale_colour_manual(values=c("#FFD4A9", "#80D1C8"))
p
 
p1 <- p + geom_half_violin(
  position = position_nudge(x = 0.1, y = 0),
  side = 'R',
  adjust = 1.2,
  trim = FALSE,
  color = NA,
  alpha = 0.6
)
p1
 
p2 <- p1 + geom_half_point(
  position = position_nudge(x = -0.35, y = 0),
  size = 3,
  shape = 19,
  range_scale = 0.5,
  alpha = 0.6
)
p2 
 
p3 <- p2 + geom_boxplot(
  outlier.shape = NA,   
  width = 0.1,
  alpha = 0.8
)
p3
 
p4 <- p3 + stat_compare_means(
 
  method = "wilcox.test",   
  label = "p.format",   
  label.x = 1.5,   
  size = 5   
)
p4

# windowsFonts(A = windowsFont("TimesNewRoman"),
#              B = windowsFont("Arial"))

p5 <- p4 + theme_light() + theme(panel.grid = element_blank()) +
  theme(
    axis.ticks.length = unit(-0.25, "cm"),
    axis.text.x = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    axis.text.y = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    panel.grid.major = element_line(linetype = "dotted", color = "lightgray", size = 0.5),  #  
    text = element_text(size = 20)
  )+
  coord_cartesian(ylim = c(-4, 4))  #  
p5
 
print(p5)
 
pdf("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\ESCA_cluster_tide_plot1.pdf",width = 7,height = 7)
p5
dev.off()
#ggsave("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\ESCA_cluster_tide_plot.pdf",p5,width = 10,height = 10)

 
label<-read.csv("G:\\pan_result\\gpu_result_tcga\\ESCA\\label.csv")

clusters_tide<-data.frame(Patient=label$X,group=label$label)
library(tidyverse)
 
tide_final<-left_join(tide,clusters_tide,by="Patient")
T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
p_tide_zanshi<-T_result$p.value 
library(ggplot2)
library(gghalves)
library(ggprism)   
library(ggpubr)   
tide_plot<-data.frame(Tide=as.numeric(tide_final$TIDE),Stage=as.character(tide_final$group))
p<-ggplot(tide_plot,aes(x=Stage,y=Tide,fill=Stage,color=Stage))+
  scale_fill_manual(values=c("#D55E00","#56B4E9"))+
  scale_colour_manual(values=c("#D55E00","#56B4E9"))
 
p1 <- p + geom_half_violin(
  position = position_nudge(x = 0.1, y = 0),
  side = 'R',
  adjust = 1.2,
  trim = FALSE,
  color = NA,
  alpha = 0.6
)

 
p2 <- p1 + geom_half_point(
  position = position_nudge(x = -0.35, y = 0),
  size = 3,
  shape = 19,
  range_scale = 0.5,
  alpha = 0.6
)

 
p3 <- p2 + geom_boxplot(
  outlier.shape = NA,  
  width = 0.1,
  alpha = 0.8
)

 
p4 <- p3 + stat_compare_means(
 
  method = "wilcox.test",  
  label = "p.format",   
  label.x = 1.5,   
  size = 5  
)

 
windowsFonts(A = windowsFont("TimesNewRoman"),
             B = windowsFont("Arial"))

p5 <- p4 + theme_light() + theme(panel.grid = element_blank()) +
  theme(
    axis.ticks.length = unit(-0.25, "cm"),
    axis.text.x = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    axis.text.y = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    panel.grid.major = element_line(linetype = "dotted", color = "lightgray", size = 0.5),  #  
    text = element_text(size = 20)
  )

 
print(p5)
 

ggsave("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\ESCA_satge_tide_plot.pdf",p5,width = 10,height = 10)




#HNSC####

 
clusters<-as.data.frame(clusters)
clusters_tide<-data.frame(Patient=rownames(clusters),group=clusters$clusters)

#1.clustertide
 
tide_final<-left_join(tide,clusters_tide,by="Patient")
T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
p_tide_zanshi<-T_result$p.value 
library(ggplot2)
library(gghalves)
library(ggprism)   
library(ggpubr)   
tide_plot<-data.frame(Tide=as.numeric(tide_final$TIDE),Cluster=as.character(tide_final$group))
p<-ggplot(tide_plot,aes(x=Cluster,y=Tide,fill=Cluster,color=Cluster))+
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8"))+
  scale_colour_manual(values=c("#FFD4A9", "#80D1C8"))
p
 
p1 <- p + geom_half_violin(
  position = position_nudge(x = 0.1, y = 0),
  side = 'R',
  adjust = 1.2,
  trim = FALSE,
  color = NA,
  alpha = 0.6
)
p1
 
p2 <- p1 + geom_half_point(
  position = position_nudge(x = -0.35, y = 0),
  size = 3,
  shape = 19,
  range_scale = 0.5,
  alpha = 0.6
)
p2 
 
p3 <- p2 + geom_boxplot(
  outlier.shape = NA,   
  width = 0.1,
  alpha = 0.8
)
p3
 
p4 <- p3 + stat_compare_means(
 
  method = "wilcox.test",   
  label = "p.format",   
  label.x = 1.5,   
  size = 5  
)
p4
 
# windowsFonts(A = windowsFont("TimesNewRoman"),
#              B = windowsFont("Arial"))

p5 <- p4 + theme_light() + theme(panel.grid = element_blank()) +
  theme(
    axis.ticks.length = unit(-0.25, "cm"),
    axis.text.x = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    axis.text.y = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    panel.grid.major = element_line(linetype = "dotted", color = "lightgray", size = 0.5),  #  
    text = element_text(size = 20)
  )+
  coord_cartesian(ylim = c(-4, 4))  #  
p5
 
print(p5)
 
pdf("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_cluster_tide_plot1.pdf",width = 7,height = 7)
p5
dev.off()
#ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_cluster_tide_plot.pdf",p5,width = 10,height = 10)

 
label<-read.csv("G:\\pan_result\\gpu_result_tcga\\HNSC\\label.csv")

clusters_tide<-data.frame(Patient=label$X,group=label$label)
library(tidyverse)
 
tide_final<-left_join(tide,clusters_tide,by="Patient")
T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
p_tide_zanshi<-T_result$p.value
library(ggplot2)
library(gghalves)
library(ggprism)   
library(ggpubr)    
tide_plot<-data.frame(Tide=as.numeric(tide_final$TIDE),Stage=as.character(tide_final$group))
p<-ggplot(tide_plot,aes(x=Stage,y=Tide,fill=Stage,color=Stage))+
  scale_fill_manual(values=c("#D55E00","#56B4E9"))+
  scale_colour_manual(values=c("#D55E00","#56B4E9"))
 
p1 <- p + geom_half_violin(
  position = position_nudge(x = 0.1, y = 0),
  side = 'R',
  adjust = 1.2,
  trim = FALSE,
  color = NA,
  alpha = 0.6
)

 
p2 <- p1 + geom_half_point(
  position = position_nudge(x = -0.35, y = 0),
  size = 3,
  shape = 19,
  range_scale = 0.5,
  alpha = 0.6
)


p3 <- p2 + geom_boxplot(
  outlier.shape = NA,   
  width = 0.1,
  alpha = 0.8
)

 
p4 <- p3 + stat_compare_means(
 
  method = "wilcox.test",  
  label = "p.format",   
  label.x = 1.5,  
  size = 5   
)

 
windowsFonts(A = windowsFont("TimesNewRoman"),
             B = windowsFont("Arial"))

p5 <- p4 + theme_light() + theme(panel.grid = element_blank()) +
  theme(
    axis.ticks.length = unit(-0.25, "cm"),
    axis.text.x = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    axis.text.y = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    panel.grid.major = element_line(linetype = "dotted", color = "lightgray", size = 0.5),  #  
    text = element_text(size = 20)
  )

 
print(p5)
 

ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_satge_tide_plot.pdf",p5,width = 10,height = 10)


#KIRC####

 
clusters<-as.data.frame(clusters)
clusters_tide<-data.frame(Patient=rownames(clusters),group=clusters$clusters)

#1.clustertide
 
tide_final<-left_join(tide,clusters_tide,by="Patient")
T_result<-kruskal.test(as.numeric(tide_final$TIDE)~tide_final$group)
p_tide_zanshi<-T_result$p.value
library(ggplot2)
library(gghalves)
library(ggprism)   
library(ggpubr)    
tide_plot<-data.frame(Tide=as.numeric(tide_final$TIDE),Cluster=as.character(tide_final$group))
p<-ggplot(tide_plot,aes(x=Cluster,y=Tide,fill=Cluster,color=Cluster))+
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+
  scale_colour_manual(values=c("#FFD4A9", "#80D1C8","pink"))
p
 
p1 <- p + geom_half_violin(
  position = position_nudge(x = 0.1, y = 0),
  side = 'R',
  adjust = 1.2,
  trim = FALSE,
  color = NA,
  alpha = 0.6
)
p1
 
p2 <- p1 + geom_half_point(
  position = position_nudge(x = -0.35, y = 0),
  size = 3,
  shape = 19,
  range_scale = 0.5,
  alpha = 0.6
)
p2 
 
p3 <- p2 + geom_boxplot(
  outlier.shape = NA,   
  width = 0.1,
  alpha = 0.8
)
p3
 
p4 <- p3 + stat_compare_means(
 
  method = "kruskal.test",   
  label = "p.format",   
  label.x = 1.5,  
  size = 5   
)
p4
 
# windowsFonts(A = windowsFont("TimesNewRoman"),
#              B = windowsFont("Arial"))

p5 <- p4 + theme_light() + theme(panel.grid = element_blank()) +
  theme(
    axis.ticks.length = unit(-0.25, "cm"),
    axis.text.x = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    axis.text.y = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    panel.grid.major = element_line(linetype = "dotted", color = "lightgray", size = 0.5),  #  
    text = element_text(size = 20)
  )+
  coord_cartesian(ylim = c(-4, 4))  #  
p5
 
print(p5)
 
pdf("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\KIRC_cluster_tide_plot1.pdf",width = 7,height = 7)
p5
dev.off()
#ggsave("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\KIRC_cluster_tide_plot.pdf",p5,width = 10,height = 10)

 
label<-read.csv("G:\\pan_result\\gpu_result_tcga\\KIRC\\label.csv")

clusters_tide<-data.frame(Patient=label$X,group=label$label)
library(tidyverse)

tide_final<-left_join(tide,clusters_tide,by="Patient")
T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
p_tide_zanshi<-T_result$p.value 
library(ggplot2)
library(gghalves)
library(ggprism)   
library(ggpubr)    
tide_plot<-data.frame(Tide=as.numeric(tide_final$TIDE),Stage=as.character(tide_final$group))
p<-ggplot(tide_plot,aes(x=Stage,y=Tide,fill=Stage,color=Stage))+
  scale_fill_manual(values=c("#D55E00","#56B4E9"))+
  scale_colour_manual(values=c("#D55E00","#56B4E9"))
 
p1 <- p + geom_half_violin(
  position = position_nudge(x = 0.1, y = 0),
  side = 'R',
  adjust = 1.2,
  trim = FALSE,
  color = NA,
  alpha = 0.6
)


p2 <- p1 + geom_half_point(
  position = position_nudge(x = -0.35, y = 0),
  size = 3,
  shape = 19,
  range_scale = 0.5,
  alpha = 0.6
)

 
p3 <- p2 + geom_boxplot(
  outlier.shape = NA,   
  width = 0.1,
  alpha = 0.8
)

 
p4 <- p3 + stat_compare_means(
 
  method = "wilcox.test",  
  label = "p.format",  
  label.x = 1.5,   
  size = 5   
)

 
windowsFonts(A = windowsFont("TimesNewRoman"),
             B = windowsFont("Arial"))

p5 <- p4 + theme_light() + theme(panel.grid = element_blank()) +
  theme(
    axis.ticks.length = unit(-0.25, "cm"),
    axis.text.x = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    axis.text.y = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    panel.grid.major = element_line(linetype = "dotted", color = "lightgray", size = 0.5),  #  
    text = element_text(size = 20)
  )

 
print(p5)
 

ggsave("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\KIRC_satge_tide_plot.pdf",p5,width = 10,height = 10)



#KIRP####

 
clusters<-as.data.frame(clusters)
clusters_tide<-data.frame(Patient=rownames(clusters),group=clusters$clusters)

#1.clustertide

tide_final<-left_join(tide,clusters_tide,by="Patient")
T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
p_tide_zanshi<-T_result$p.value 
library(ggplot2)
library(gghalves)
library(ggprism)   
library(ggpubr)    


lower_percentile <- 0.075  
upper_percentile <- 0.925


lower_limit <- quantile(tide_plot$Tide, lower_percentile)
upper_limit <- quantile(tide_plot$Tide, upper_percentile)


tide_plot_clean <- tide_plot[tide_plot$Tide >= lower_limit & tide_plot$Tide <= upper_limit, ]




T_result<-wilcox.test(as.numeric(tide_plot_clean$Tide)~tide_plot_clean$Cluster)
p_tide_zanshi<-T_result$p.value 


tide_plot<-data.frame(Tide=as.numeric(tide_final$TIDE),Cluster=as.character(tide_final$group))
p<-ggplot(tide_plot,aes(x=Cluster,y=Tide,fill=Cluster,color=Cluster))+
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8"))+
  scale_colour_manual(values=c("#FFD4A9", "#80D1C8"))
 
p1 <- p + geom_half_violin(
  position = position_nudge(x = 0.1, y = 0),
  side = 'R',
  adjust = 1.2,
  trim = FALSE,
  color = NA,
  alpha = 0.6
)

 
p2 <- p1 + geom_half_point(
  position = position_nudge(x = -0.35, y = 0),
  size = 3,
  shape = 19,
  range_scale = 0.5,
  alpha = 0.6
)

 
p3 <- p2 + geom_boxplot(
  outlier.shape = NA,   
  width = 0.1,
  alpha = 0.8
)

 
p4 <- p3 + stat_compare_means(
 
  method = "wilcox.test",   
  label = "p.format",  
  label.x = 1.5,   
  size = 5   
)

 
windowsFonts(A = windowsFont("TimesNewRoman"),
             B = windowsFont("Arial"))

p5 <- p4 + theme_light() + theme(panel.grid = element_blank()) +
  theme(
    axis.ticks.length = unit(-0.25, "cm"),
    axis.text.x = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    axis.text.y = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    panel.grid.major = element_line(linetype = "dotted", color = "lightgray", size = 0.5),  #  
    text = element_text(size = 20)
  )+
  coord_cartesian(ylim = c(-4, 4))  #  

 
print(p5)
 

ggsave("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\KIRP_cluster_tide_plot1.pdf",p5,width = 10,height = 10)


label<-read.csv("G:\\pan_result\\gpu_result_tcga\\KIRP\\label.csv")

clusters_tide<-data.frame(Patient=label$X,group=label$label)
library(tidyverse)
 
tide_final<-left_join(tide,clusters_tide,by="Patient")
T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
p_tide_zanshi<-T_result$p.value 
library(ggplot2)
library(ggprism)   
library(ggpubr)    
tide_plot<-data.frame(Tide=as.numeric(tide_final$TIDE),Stage=as.character(tide_final$group))
p<-ggplot(tide_plot,aes(x=Stage,y=Tide,fill=Stage,color=Stage))+
  scale_fill_manual(values=c("#D55E00","#56B4E9"))+
  scale_colour_manual(values=c("#D55E00","#56B4E9"))
 
p1 <- p + geom_half_violin(
  position = position_nudge(x = 0.1, y = 0),
  side = 'R',
  adjust = 1.2,
  trim = FALSE,
  color = NA,
  alpha = 0.6
)

 
p2 <- p1 + geom_half_point(
  position = position_nudge(x = -0.35, y = 0),
  size = 3,
  shape = 19,
  range_scale = 0.5,
  alpha = 0.6
)

 
p3 <- p2 + geom_boxplot(
  outlier.shape = NA,   
  width = 0.1,
  alpha = 0.8
)

 
p4 <- p3 + stat_compare_means(
 
  method = "wilcox.test", 
  label = "p.format",   
  label.x = 1.5,   
  size = 5   
)

 
windowsFonts(A = windowsFont("TimesNewRoman"),
             B = windowsFont("Arial"))

p5 <- p4 + theme_light() + theme(panel.grid = element_blank()) +
  theme(
    axis.ticks.length = unit(-0.25, "cm"),
    axis.text.x = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    axis.text.y = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    panel.grid.major = element_line(linetype = "dotted", color = "lightgray", size = 0.5),  #  
    text = element_text(size = 20)
  )

 
print(p5)
 

ggsave("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\KIRP_satge_tide_plot.pdf",p5,width = 10,height = 10)


#LIHC####


clusters<-as.data.frame(clusters)
clusters_tide<-data.frame(Patient=rownames(clusters),group=clusters$clusters)

#1.clustertide
 
tide_final<-left_join(tide,clusters_tide,by="Patient")
T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
p_tide_zanshi<-T_result$p.value 
library(ggplot2)
library(gghalves)
library(ggprism)   
library(ggpubr)    
tide_plot<-data.frame(Tide=as.numeric(tide_final$TIDE),Cluster=as.character(tide_final$group))
p<-ggplot(tide_plot,aes(x=Cluster,y=Tide,fill=Cluster,color=Cluster))+
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8"))+
  scale_colour_manual(values=c("#FFD4A9", "#80D1C8"))
p
 
p1 <- p + geom_half_violin(
  position = position_nudge(x = 0.1, y = 0),
  side = 'R',
  adjust = 1.2,
  trim = FALSE,
  color = NA,
  alpha = 0.6
)
p1
 
p2 <- p1 + geom_half_point(
  position = position_nudge(x = -0.35, y = 0),
  size = 3,
  shape = 19,
  range_scale = 0.5,
  alpha = 0.6
)
p2 
 
p3 <- p2 + geom_boxplot(
  outlier.shape = NA,   
  width = 0.1,
  alpha = 0.8
)
p3

p4 <- p3 + stat_compare_means(
 
  method = "wilcox.test",   
  label = "p.format",   
  label.x = 1.5,   
  size = 5   
)
p4
 
# windowsFonts(A = windowsFont("TimesNewRoman"),
#              B = windowsFont("Arial"))

p5 <- p4 + theme_light() + theme(panel.grid = element_blank()) +
  theme(
    axis.ticks.length = unit(-0.25, "cm"),
    axis.text.x = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    axis.text.y = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    panel.grid.major = element_line(linetype = "dotted", color = "lightgray", size = 0.5),  # 
    text = element_text(size = 20)
  )+
  coord_cartesian(ylim = c(-4, 4))  #  
p5
 
print(p5)
 
pdf("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\LIHC_cluster_tide_plot1.pdf",width = 7,height = 7)
p5
dev.off()
#ggsave("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\LIHC_cluster_tide_plot.pdf",p5,width = 10,height = 10)

 
label<-read.csv("G:\\pan_result\\gpu_result_tcga\\LIHC\\label.csv")

clusters_tide<-data.frame(Patient=label$X,group=label$label)
library(tidyverse)
 
tide_final<-left_join(tide,clusters_tide,by="Patient")
T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
p_tide_zanshi<-T_result$p.value 
library(ggplot2)
library(gghalves)
library(ggprism)   
library(ggpubr)    
tide_plot<-data.frame(Tide=as.numeric(tide_final$TIDE),Stage=as.character(tide_final$group))
p<-ggplot(tide_plot,aes(x=Stage,y=Tide,fill=Stage,color=Stage))+
  scale_fill_manual(values=c("#D55E00","#56B4E9"))+
  scale_colour_manual(values=c("#D55E00","#56B4E9"))
 
p1 <- p + geom_half_violin(
  position = position_nudge(x = 0.1, y = 0),
  side = 'R',
  adjust = 1.2,
  trim = FALSE,
  color = NA,
  alpha = 0.6
)

 
p2 <- p1 + geom_half_point(
  position = position_nudge(x = -0.35, y = 0),
  size = 3,
  shape = 19,
  range_scale = 0.5,
  alpha = 0.6
)

 
p3 <- p2 + geom_boxplot(
  outlier.shape = NA,   
  width = 0.1,
  alpha = 0.8
)

 
p4 <- p3 + stat_compare_means(
 
  method = "wilcox.test",  
  label = "p.format",   
  label.x = 1.5,   
  size = 5   
)

 
windowsFonts(A = windowsFont("TimesNewRoman"),
             B = windowsFont("Arial"))

p5 <- p4 + theme_light() + theme(panel.grid = element_blank()) +
  theme(
    axis.ticks.length = unit(-0.25, "cm"),
    axis.text.x = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    axis.text.y = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    panel.grid.major = element_line(linetype = "dotted", color = "lightgray", size = 0.5),  #  
    text = element_text(size = 20)
  )

 
print(p5)
 

ggsave("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\LIHC_satge_tide_plot.pdf",p5,width = 10,height = 10)




#LUNG####

 
clusters<-as.data.frame(clusters)
clusters_tide<-data.frame(Patient=rownames(clusters),group=clusters$clusters)

#1.clustertide
 
tide_final<-left_join(tide,clusters_tide,by="Patient")
T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
p_tide_zanshi<-T_result$p.value 
library(ggplot2)
library(gghalves)
library(ggprism)   
library(ggpubr)    
tide_plot<-data.frame(Tide=as.numeric(tide_final$TIDE),Cluster=as.character(tide_final$group))
p<-ggplot(tide_plot,aes(x=Cluster,y=Tide,fill=Cluster,color=Cluster))+
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8"))+
  scale_colour_manual(values=c("#FFD4A9", "#80D1C8"))
p
 
p1 <- p + geom_half_violin(
  position = position_nudge(x = 0.1, y = 0),
  side = 'R',
  adjust = 1.2,
  trim = FALSE,
  color = NA,
  alpha = 0.6
)
p1
 
p2 <- p1 + geom_half_point(
  position = position_nudge(x = -0.35, y = 0),
  size = 3,
  shape = 19,
  range_scale = 0.5,
  alpha = 0.6
)
p2 
 
p3 <- p2 + geom_boxplot(
  outlier.shape = NA,   
  width = 0.1,
  alpha = 0.8
)
p3
 
p4 <- p3 + stat_compare_means(

  method = "wilcox.test",   
  label = "p.format",   
  label.x = 1.5,   
  size = 5   
)
p4
 
# windowsFonts(A = windowsFont("TimesNewRoman"),
#              B = windowsFont("Arial"))

p5 <- p4 + theme_light() + theme(panel.grid = element_blank()) +
  theme(
    axis.ticks.length = unit(-0.25, "cm"),
    axis.text.x = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    axis.text.y = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    panel.grid.major = element_line(linetype = "dotted", color = "lightgray", size = 0.5),  #  
    text = element_text(size = 20)
  )+
  coord_cartesian(ylim = c(-4, 4))  #  
p5
 
print(p5)
 
pdf("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\NSCLC_cluster_tide_plot1.pdf",width = 7,height = 7)
p5
dev.off()
#ggsave("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\NSCLC_cluster_tide_plot.pdf",p5,width = 10,height = 10)
getwd()
#####
write.table(clusters,"G:/pan_result/gpu_result_tcga/NSCLC/450k/NSCLC_cluster.txt")

 
label<-read.csv("G:\\pan_result\\gpu_result_tcga\\NSCLC\\label.csv")

clusters_tide<-data.frame(Patient=label$X,group=label$label)
library(tidyverse)
 
tide_final<-left_join(tide,clusters_tide,by="Patient")
T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
p_tide_zanshi<-T_result$p.value 
library(ggplot2)
library(gghalves)
library(ggprism)   
library(ggpubr)    
tide_plot<-data.frame(Tide=as.numeric(tide_final$TIDE),Stage=as.character(tide_final$group))
p<-ggplot(tide_plot,aes(x=Stage,y=Tide,fill=Stage,color=Stage))+
  scale_fill_manual(values=c("#D55E00","#56B4E9"))+
  scale_colour_manual(values=c("#D55E00","#56B4E9"))
 
p1 <- p + geom_half_violin(
  position = position_nudge(x = 0.1, y = 0),
  side = 'R',
  adjust = 1.2,
  trim = FALSE,
  color = NA,
  alpha = 0.6
)

 
p2 <- p1 + geom_half_point(
  position = position_nudge(x = -0.35, y = 0),
  size = 3,
  shape = 19,
  range_scale = 0.5,
  alpha = 0.6
)

 
p3 <- p2 + geom_boxplot(
  outlier.shape = NA,   
  width = 0.1,
  alpha = 0.8
)

 
p4 <- p3 + stat_compare_means(
 
  method = "wilcox.test",  
  label = "p.format",   
  label.x = 1.5,   
  size = 5   
)

 
windowsFonts(A = windowsFont("TimesNewRoman"),
             B = windowsFont("Arial"))

p5 <- p4 + theme_light() + theme(panel.grid = element_blank()) +
  theme(
    axis.ticks.length = unit(-0.25, "cm"),
    axis.text.x = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    axis.text.y = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    panel.grid.major = element_line(linetype = "dotted", color = "lightgray", size = 0.5),  #  
    text = element_text(size = 20)
  )

 
print(p5)
 

ggsave("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\NSCLC_satge_tide_plot.pdf",p5,width = 10,height = 10)

#PAAD####

 
clusters<-as.data.frame(clusters)
clusters_tide<-data.frame(Patient=rownames(clusters),group=clusters$clusters)

#1.clustertide
 
tide_final<-left_join(tide,clusters_tide,by="Patient")
T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
p_tide_zanshi<-T_result$p.value 
library(ggplot2)
library(gghalves)
library(ggprism)   
library(ggpubr)    
tide_plot<-data.frame(Tide=as.numeric(tide_final$TIDE),Cluster=as.character(tide_final$group))
p<-ggplot(tide_plot,aes(x=Cluster,y=Tide,fill=Cluster,color=Cluster))+
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8"))+
  scale_colour_manual(values=c("#FFD4A9", "#80D1C8"))
p
 
p1 <- p + geom_half_violin(
  position = position_nudge(x = 0.1, y = 0),
  side = 'R',
  adjust = 1.2,
  trim = FALSE,
  color = NA,
  alpha = 0.6
)
p1
 
p2 <- p1 + geom_half_point(
  position = position_nudge(x = -0.35, y = 0),
  size = 3,
  shape = 19,
  range_scale = 0.5,
  alpha = 0.6
)
p2 
 
p3 <- p2 + geom_boxplot(
  outlier.shape = NA,   
  width = 0.1,
  alpha = 0.8
)
p3
 
p4 <- p3 + stat_compare_means(
 
  method = "wilcox.test",  
  label = "p.format",   
  label.x = 1.5,   
  size = 5   
)
p4
 
# windowsFonts(A = windowsFont("TimesNewRoman"),
#              B = windowsFont("Arial"))

p5 <- p4 + theme_light() + theme(panel.grid = element_blank()) +
  theme(
    axis.ticks.length = unit(-0.25, "cm"),
    axis.text.x = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    axis.text.y = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    panel.grid.major = element_line(linetype = "dotted", color = "lightgray", size = 0.5),  #  
    text = element_text(size = 20)
  )+
  coord_cartesian(ylim = c(-4, 4))  #  
p5
 
print(p5)
 
pdf("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\PAAD_cluster_tide_plot1.pdf",width = 7,height = 7)
p5
dev.off()
#ggsave("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\PAAD_cluster_tide_plot.pdf",p5,width = 10,height = 10)


 
label<-read.csv("G:\\pan_result\\gpu_result_tcga\\PAAD\\label.csv")

clusters_tide<-data.frame(Patient=label$X,group=label$label)
library(tidyverse)
 
tide_final<-left_join(tide,clusters_tide,by="Patient")
T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
p_tide_zanshi<-T_result$p.value 
library(ggplot2)
library(gghalves)
library(ggprism)   
library(ggpubr)    
tide_plot<-data.frame(Tide=as.numeric(tide_final$TIDE),Stage=as.character(tide_final$group))
p<-ggplot(tide_plot,aes(x=Stage,y=Tide,fill=Stage,color=Stage))+
  scale_fill_manual(values=c("#D55E00","#56B4E9"))+
  scale_colour_manual(values=c("#D55E00","#56B4E9"))
 
p1 <- p + geom_half_violin(
  position = position_nudge(x = 0.1, y = 0),
  side = 'R',
  adjust = 1.2,
  trim = FALSE,
  color = NA,
  alpha = 0.6
)

 
p2 <- p1 + geom_half_point(
  position = position_nudge(x = -0.35, y = 0),
  size = 3,
  shape = 19,
  range_scale = 0.5,
  alpha = 0.6
)

 
p3 <- p2 + geom_boxplot(
  outlier.shape = NA,   
  width = 0.1,
  alpha = 0.8
)

 
p4 <- p3 + stat_compare_means(
 
  method = "wilcox.test",  
  label = "p.format",   
  label.x = 1.5,   
  size = 5   
)

 
windowsFonts(A = windowsFont("TimesNewRoman"),
             B = windowsFont("Arial"))

p5 <- p4 + theme_light() + theme(panel.grid = element_blank()) +
  theme(
    axis.ticks.length = unit(-0.25, "cm"),
    axis.text.x = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    axis.text.y = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    panel.grid.major = element_line(linetype = "dotted", color = "lightgray", size = 0.5),  #  
    text = element_text(size = 20)
  )

 
print(p5)
 

ggsave("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\PAAD_satge_tide_plot.pdf",p5,width = 10,height = 10)


#THCA####

 
clusters<-as.data.frame(clusters)
clusters_tide<-data.frame(Patient=rownames(clusters),group=clusters$clusters)

#1.clustertide
 
tide_final<-left_join(tide,clusters_tide,by="Patient")
T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
p_tide_zanshi<-T_result$p.value 
library(ggplot2)
library(gghalves)
library(ggprism)   
library(ggpubr)    
tide_plot<-data.frame(Tide=as.numeric(tide_final$TIDE),Cluster=as.character(tide_final$group))
p<-ggplot(tide_plot,aes(x=Cluster,y=Tide,fill=Cluster,color=Cluster))+
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8"))+
  scale_colour_manual(values=c("#FFD4A9", "#80D1C8"))
p
 
p1 <- p + geom_half_violin(
  position = position_nudge(x = 0.1, y = 0),
  side = 'R',
  adjust = 1.2,
  trim = FALSE,
  color = NA,
  alpha = 0.6
)
p1
 
p2 <- p1 + geom_half_point(
  position = position_nudge(x = -0.35, y = 0),
  size = 3,
  shape = 19,
  range_scale = 0.5,
  alpha = 0.6
)
p2 
 
p3 <- p2 + geom_boxplot(
  outlier.shape = NA,  
  width = 0.1,
  alpha = 0.8
)
p3
 
p4 <- p3 + stat_compare_means(
 
  method = "wilcox.test",   
  label = "p.format",   
  label.x = 1.5,   
  size = 5   
)
p4
 
# windowsFonts(A = windowsFont("TimesNewRoman"),
#              B = windowsFont("Arial"))

p5 <- p4 + theme_light() + theme(panel.grid = element_blank()) +
  theme(
    axis.ticks.length = unit(-0.25, "cm"),
    axis.text.x = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    axis.text.y = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    panel.grid.major = element_line(linetype = "dotted", color = "lightgray", size = 0.5),  #  
    text = element_text(size = 20)
  )+
  coord_cartesian(ylim = c(-4, 4))  #  
p5
 
print(p5)
 
pdf("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\THCA_cluster_tide_plot1.pdf",width = 7,height = 7)
p5
dev.off()

 
label<-read.csv("G:\\pan_result\\gpu_result_tcga\\THCA\\label.csv")

clusters_tide<-data.frame(Patient=label$X,group=label$label)
library(tidyverse)
 
tide_final<-left_join(tide,clusters_tide,by="Patient")
T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
p_tide_zanshi<-T_result$p.value 
library(ggplot2)
library(gghalves)
library(ggprism)   
library(ggpubr)    
tide_plot<-data.frame(Tide=as.numeric(tide_final$TIDE),Stage=as.character(tide_final$group))
p<-ggplot(tide_plot,aes(x=Stage,y=Tide,fill=Stage,color=Stage))+
  scale_fill_manual(values=c("#D55E00","#56B4E9"))+
  scale_colour_manual(values=c("#D55E00","#56B4E9"))
 
p1 <- p + geom_half_violin(
  position = position_nudge(x = 0.1, y = 0),
  side = 'R',
  adjust = 1.2,
  trim = FALSE,
  color = NA,
  alpha = 0.6
)

 
p2 <- p1 + geom_half_point(
  position = position_nudge(x = -0.35, y = 0),
  size = 3,
  shape = 19,
  range_scale = 0.5,
  alpha = 0.6
)

 
p3 <- p2 + geom_boxplot(
  outlier.shape = NA,   
  width = 0.1,
  alpha = 0.8
)

 
p4 <- p3 + stat_compare_means(
 
  method = "wilcox.test",  
  label = "p.format",   
  label.x = 1.5,   
  size = 5   
)

 
windowsFonts(A = windowsFont("TimesNewRoman"),
             B = windowsFont("Arial"))

p5 <- p4 + theme_light() + theme(panel.grid = element_blank()) +
  theme(
    axis.ticks.length = unit(-0.25, "cm"),
    axis.text.x = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    axis.text.y = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    panel.grid.major = element_line(linetype = "dotted", color = "lightgray", size = 0.5),  #  
    text = element_text(size = 20)
  )

 
print(p5)
 

ggsave("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\THCA_satge_tide_plot.pdf",p5,width = 10,height = 10)


#UCEC####

 
clusters<-as.data.frame(clusters)
clusters_tide<-data.frame(Patient=rownames(clusters),group=clusters$clusters)

#1.clustertide
 
tide_final<-left_join(tide,clusters_tide,by="Patient")
T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
p_tide_zanshi<-T_result$p.value 
library(ggplot2)
library(gghalves)
library(ggprism)   
library(ggpubr)    
tide_plot<-data.frame(Tide=as.numeric(tide_final$TIDE),Cluster=as.character(tide_final$group))
p<-ggplot(tide_plot,aes(x=Cluster,y=Tide,fill=Cluster,color=Cluster))+
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8"))+
  scale_colour_manual(values=c("#FFD4A9", "#80D1C8"))
p
 
p1 <- p + geom_half_violin(
  position = position_nudge(x = 0.1, y = 0),
  side = 'R',
  adjust = 1.2,
  trim = FALSE,
  color = NA,
  alpha = 0.6
)
p1
 
p2 <- p1 + geom_half_point(
  position = position_nudge(x = -0.35, y = 0),
  size = 3,
  shape = 19,
  range_scale = 0.5,
  alpha = 0.6
)
p2 

p3 <- p2 + geom_boxplot(
  outlier.shape = NA,   
  width = 0.1,
  alpha = 0.8
)
p3
 
p4 <- p3 + stat_compare_means(
 
  method = "wilcox.test",   
  label = "p.format",   
  label.x = 1.5,   
  size = 5   
)
p4
 
# windowsFonts(A = windowsFont("TimesNewRoman"),
#              B = windowsFont("Arial"))

p5 <- p4 + theme_light() + theme(panel.grid = element_blank()) +
  theme(
    axis.ticks.length = unit(-0.25, "cm"),
    axis.text.x = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    axis.text.y = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    panel.grid.major = element_line(linetype = "dotted", color = "lightgray", size = 0.5),  #  
    text = element_text(size = 20)
  )+
  coord_cartesian(ylim = c(-4, 4))  #  
p5
 
print(p5)

pdf("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\UCEC_cluster_tide_plot1.pdf",width = 7,height = 7)
p5
dev.off()
#ggsave("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\UCEC_cluster_tide_plot.pdf",p5,width = 10,height = 10)

 
label<-read.csv("G:\\pan_result\\gpu_result_tcga\\UCEC\\label.csv")

clusters_tide<-data.frame(Patient=label$X,group=label$label)
library(tidyverse)
 
tide_final<-left_join(tide,clusters_tide,by="Patient")
T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
p_tide_zanshi<-T_result$p.value 
library(ggplot2)
library(gghalves)
library(ggprism)   
library(ggpubr)    
tide_plot<-data.frame(Tide=as.numeric(tide_final$TIDE),Stage=as.character(tide_final$group))
p<-ggplot(tide_plot,aes(x=Stage,y=Tide,fill=Stage,color=Stage))+
  scale_fill_manual(values=c("#D55E00","#56B4E9"))+
  scale_colour_manual(values=c("#D55E00","#56B4E9"))
 
p1 <- p + geom_half_violin(
  position = position_nudge(x = 0.1, y = 0),
  side = 'R',
  adjust = 1.2,
  trim = FALSE,
  color = NA,
  alpha = 0.6
)

 
p2 <- p1 + geom_half_point(
  position = position_nudge(x = -0.35, y = 0),
  size = 3,
  shape = 19,
  range_scale = 0.5,
  alpha = 0.6
)

 
p3 <- p2 + geom_boxplot(
  outlier.shape = NA,   
  width = 0.1,
  alpha = 0.8
)

 
p4 <- p3 + stat_compare_means(
 
  method = "wilcox.test",  
  label = "p.format",   
  label.x = 1.5,   
  size = 5   
)

 
windowsFonts(A = windowsFont("TimesNewRoman"),
             B = windowsFont("Arial"))

p5 <- p4 + theme_light() + theme(panel.grid = element_blank()) +
  theme(
    axis.ticks.length = unit(-0.25, "cm"),
    axis.text.x = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    axis.text.y = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
    panel.grid.major = element_line(linetype = "dotted", color = "lightgray", size = 0.5),  #  
    text = element_text(size = 20)
  )

 
print(p5)


ggsave("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\UCEC_satge_tide_plot.pdf",p5,width = 10,height = 10)

