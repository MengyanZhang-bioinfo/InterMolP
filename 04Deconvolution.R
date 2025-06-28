#BRCA####
#1.3.2  ####
# 
library(ggplot2)
library(ggthemes)
class<-rep(c("cluster1","cluster2"),c(nrow(cluster1_ciber),nrow(cluster2_ciber)))
mianyi_huatu<-c()
for (i in 1:nrow(ciber_zong)) {
  b<-as.numeric(ciber_zong[i,])
  tonglu1<-rep(rownames(ciber_zong)[i],ncol(ciber_zong))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)

class(mianyi_huatu)
table(mianyi_huatu$type)
mianyi_huatu$scores<-as.numeric(mianyi_huatu$scores)
mianyi_huatu$class<-as.factor(mianyi_huatu$class)

p_zanshi2<-data.frame(p_zanshi)
rownames(p_zanshi2)<-rownames(ciber_zong)
# p_zanshi2<-data.frame(p=p_zanshi2[name,])
# rownames(p_zanshi2)<-name
# label<-p_zanshi2$p_zanshi
# label1<-p_zanshi2$p_zanshi
# label[which(label1>=0.05)]<-"ns"
# label[which(label1<0.001)]<-"***"
# label[which(label1<0.01&label1>=0.001)]<-"**"
# label[which(label1<0.05&label1>=0.01)]<-"*"
library(ggpubr)
bbb<-ggplot(mianyi_huatu,aes(y=scores,x=factor(type),fill=class))+
  geom_boxplot(outlier.size=1,outlier.shape = NA)+
  theme_classic()+
  xlab("")+ylab("")+guides()+ylim(0,0.5)+
  scale_fill_manual(values=c("#80D1C8","#FFD4A9", "pink"))+stat_compare_means(method = "wilcox.test",label = "p.signif")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  #  
# stat_compare_means(aes(group = class),method = "wilcox.test",label = "p.signif")
# geom_text(aes(x =1, y=0.5, label = label[1], group = type), size = 4) +
#   geom_text(aes(x =2, y=0.5, label = label[2], group = type), size = 4)+
#   geom_text(aes(x =3, y=0.5, label = label[3], group = type), size = 4)+
#   geom_text(aes(x =4, y=0.5, label = label[4], group = type), size = 4)+
#   geom_text(aes(x =5, y=0.5, label = label[5], group = type), size = 4)+
#   geom_text(aes(x =6, y=0.5, label = label[6], group = type), size = 4)+
#   geom_text(aes(x =7, y=0.5, label = label[7], group = type), size = 4)
bbb
ggsave("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\BRCA_cluster_scrna_plot.pdf",bbb,width = 10,height = 10)
#COAD####
#1.3.2 ####
# 
library(ggplot2)
library(ggthemes)
class<-rep(c("cluster1","cluster2"),c(nrow(cluster1_ciber),nrow(cluster2_ciber)))
mianyi_huatu<-c()
for (i in 1:nrow(ciber_zong)) {
  b<-as.numeric(ciber_zong[i,])
  tonglu1<-rep(rownames(ciber_zong)[i],ncol(ciber_zong))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)

class(mianyi_huatu)
table(mianyi_huatu$type)
mianyi_huatu$scores<-as.numeric(mianyi_huatu$scores)
mianyi_huatu$class<-as.factor(mianyi_huatu$class)

p_zanshi2<-data.frame(p_zanshi)
rownames(p_zanshi2)<-rownames(ciber_zong)
# p_zanshi2<-data.frame(p=p_zanshi2[name,])
# rownames(p_zanshi2)<-name
# label<-p_zanshi2$p_zanshi
# label1<-p_zanshi2$p_zanshi
# label[which(label1>=0.05)]<-"ns"
# label[which(label1<0.001)]<-"***"
# label[which(label1<0.01&label1>=0.001)]<-"**"
# label[which(label1<0.05&label1>=0.01)]<-"*"
library(ggpubr)
bbb<-ggplot(mianyi_huatu,aes(y=scores,x=factor(type),fill=class))+
  geom_boxplot(outlier.size=1,outlier.shape = NA)+
  theme_classic()+
  xlab("")+ylab("")+guides()+ylim(0,0.5)+
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+stat_compare_means(method = "wilcox.test",label = "p.signif")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  #  
# stat_compare_means(aes(group = class),method = "wilcox.test",label = "p.signif")
# geom_text(aes(x =1, y=0.5, label = label[1], group = type), size = 4) +
#   geom_text(aes(x =2, y=0.5, label = label[2], group = type), size = 4)+
#   geom_text(aes(x =3, y=0.5, label = label[3], group = type), size = 4)+
#   geom_text(aes(x =4, y=0.5, label = label[4], group = type), size = 4)+
#   geom_text(aes(x =5, y=0.5, label = label[5], group = type), size = 4)+
#   geom_text(aes(x =6, y=0.5, label = label[6], group = type), size = 4)+
#   geom_text(aes(x =7, y=0.5, label = label[7], group = type), size = 4)
bbb
ggsave("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\COAD_cluster_scrna_plot.pdf",bbb,width = 10,height = 10)
#COAD####
#1.3.2  ####
#
library(ggplot2)
library(ggthemes)
class<-rep(c("cluster1","cluster2"),c(nrow(cluster1_ciber),nrow(cluster2_ciber)))
mianyi_huatu<-c()
for (i in 1:nrow(ciber_zong)) {
  b<-as.numeric(ciber_zong[i,])
  tonglu1<-rep(rownames(ciber_zong)[i],ncol(ciber_zong))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)

class(mianyi_huatu)
table(mianyi_huatu$type)
mianyi_huatu$scores<-as.numeric(mianyi_huatu$scores)
mianyi_huatu$class<-as.factor(mianyi_huatu$class)

p_zanshi2<-data.frame(p_zanshi)
rownames(p_zanshi2)<-rownames(ciber_zong)
# p_zanshi2<-data.frame(p=p_zanshi2[name,])
# rownames(p_zanshi2)<-name
# label<-p_zanshi2$p_zanshi
# label1<-p_zanshi2$p_zanshi
# label[which(label1>=0.05)]<-"ns"
# label[which(label1<0.001)]<-"***"
# label[which(label1<0.01&label1>=0.001)]<-"**"
# label[which(label1<0.05&label1>=0.01)]<-"*"
library(ggpubr)
bbb<-ggplot(mianyi_huatu,aes(y=scores,x=factor(type),fill=class))+
  geom_boxplot(outlier.size=1,outlier.shape = NA)+
  theme_classic()+
  xlab("")+ylab("")+guides()+ylim(0,0.5)+
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+stat_compare_means(method = "wilcox.test",label = "p.signif")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  #  
# stat_compare_means(aes(group = class),method = "wilcox.test",label = "p.signif")
# geom_text(aes(x =1, y=0.5, label = label[1], group = type), size = 4) +
#   geom_text(aes(x =2, y=0.5, label = label[2], group = type), size = 4)+
#   geom_text(aes(x =3, y=0.5, label = label[3], group = type), size = 4)+
#   geom_text(aes(x =4, y=0.5, label = label[4], group = type), size = 4)+
#   geom_text(aes(x =5, y=0.5, label = label[5], group = type), size = 4)+
#   geom_text(aes(x =6, y=0.5, label = label[6], group = type), size = 4)+
#   geom_text(aes(x =7, y=0.5, label = label[7], group = type), size = 4)
bbb
ggsave("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\COAD_cluster_scrna_plot.pdf",bbb,width = 10,height = 10)
#ESCA####

#1.3.2  ####
# 
library(ggplot2)
library(ggthemes)
class<-rep(c("cluster1","cluster2"),c(nrow(cluster1_ciber),nrow(cluster2_ciber)))
mianyi_huatu<-c()
for (i in 1:nrow(ciber_zong)) {
  b<-as.numeric(ciber_zong[i,])
  tonglu1<-rep(rownames(ciber_zong)[i],ncol(ciber_zong))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)

class(mianyi_huatu)
table(mianyi_huatu$type)
mianyi_huatu$scores<-as.numeric(mianyi_huatu$scores)
mianyi_huatu$class<-as.factor(mianyi_huatu$class)

p_zanshi2<-data.frame(p_zanshi)
rownames(p_zanshi2)<-rownames(ciber_zong)
# p_zanshi2<-data.frame(p=p_zanshi2[name,])
# rownames(p_zanshi2)<-name
# label<-p_zanshi2$p_zanshi
# label1<-p_zanshi2$p_zanshi
# label[which(label1>=0.05)]<-"ns"
# label[which(label1<0.001)]<-"***"
# label[which(label1<0.01&label1>=0.001)]<-"**"
# label[which(label1<0.05&label1>=0.01)]<-"*"
library(ggpubr)
bbb<-ggplot(mianyi_huatu,aes(y=scores,x=factor(type),fill=class))+
  geom_boxplot(outlier.size=1,outlier.shape = NA)+
  theme_classic()+
  xlab("")+ylab("")+guides()+#ylim(0,0.5)+
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+stat_compare_means(method = "wilcox.test",label = "p.signif")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  #  
# stat_compare_means(aes(group = class),method = "wilcox.test",label = "p.signif")
# geom_text(aes(x =1, y=0.5, label = label[1], group = type), size = 4) +
#   geom_text(aes(x =2, y=0.5, label = label[2], group = type), size = 4)+
#   geom_text(aes(x =3, y=0.5, label = label[3], group = type), size = 4)+
#   geom_text(aes(x =4, y=0.5, label = label[4], group = type), size = 4)+
#   geom_text(aes(x =5, y=0.5, label = label[5], group = type), size = 4)+
#   geom_text(aes(x =6, y=0.5, label = label[6], group = type), size = 4)+
#   geom_text(aes(x =7, y=0.5, label = label[7], group = type), size = 4)
bbb
ggsave("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\ESCA_cluster_scrna_plot.pdf",bbb,width = 10,height = 10)

#HNSC####

#1.3.2  ####
# 
library(ggplot2)
library(ggthemes)
class<-rep(c("cluster1","cluster2"),c(nrow(cluster1_ciber),nrow(cluster2_ciber)))
mianyi_huatu<-c()
for (i in 1:nrow(ciber_zong)) {
  b<-as.numeric(ciber_zong[i,])
  tonglu1<-rep(rownames(ciber_zong)[i],ncol(ciber_zong))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)

class(mianyi_huatu)
table(mianyi_huatu$type)
mianyi_huatu$scores<-as.numeric(mianyi_huatu$scores)
mianyi_huatu$class<-as.factor(mianyi_huatu$class)

p_zanshi2<-data.frame(p_zanshi)
rownames(p_zanshi2)<-rownames(ciber_zong)
# p_zanshi2<-data.frame(p=p_zanshi2[name,])
# rownames(p_zanshi2)<-name
# label<-p_zanshi2$p_zanshi
# label1<-p_zanshi2$p_zanshi
# label[which(label1>=0.05)]<-"ns"
# label[which(label1<0.001)]<-"***"
# label[which(label1<0.01&label1>=0.001)]<-"**"
# label[which(label1<0.05&label1>=0.01)]<-"*"
library(ggpubr)
bbb<-ggplot(mianyi_huatu,aes(y=scores,x=factor(type),fill=class))+
  geom_boxplot(outlier.size=1,outlier.shape = NA)+
  theme_classic()+
  xlab("")+ylab("")+guides()+ylim(0,0.5)+
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+stat_compare_means(method = "wilcox.test",label = "p.signif")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  #  
# stat_compare_means(aes(group = class),method = "wilcox.test",label = "p.signif")
# geom_text(aes(x =1, y=0.5, label = label[1], group = type), size = 4) +
#   geom_text(aes(x =2, y=0.5, label = label[2], group = type), size = 4)+
#   geom_text(aes(x =3, y=0.5, label = label[3], group = type), size = 4)+
#   geom_text(aes(x =4, y=0.5, label = label[4], group = type), size = 4)+
#   geom_text(aes(x =5, y=0.5, label = label[5], group = type), size = 4)+
#   geom_text(aes(x =6, y=0.5, label = label[6], group = type), size = 4)+
#   geom_text(aes(x =7, y=0.5, label = label[7], group = type), size = 4)
bbb
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_cluster_scrna_plot.pdf",bbb,width = 10,height = 10)

#KIRC####

#1.3.2  ####
# 
library(ggplot2)
library(ggthemes)
class<-rep(c("cluster1","cluster2","cluster3"),c(nrow(cluster1_ciber),nrow(cluster2_ciber),nrow(cluster3_ciber)))
mianyi_huatu<-c()
for (i in 1:nrow(ciber_zong)) {
  b<-as.numeric(ciber_zong[i,])
  tonglu1<-rep(rownames(ciber_zong)[i],ncol(ciber_zong))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)

class(mianyi_huatu)
table(mianyi_huatu$type)
mianyi_huatu$scores<-as.numeric(mianyi_huatu$scores)
mianyi_huatu$class<-as.factor(mianyi_huatu$class)

p_zanshi2<-data.frame(p_zanshi)
rownames(p_zanshi2)<-rownames(ciber_zong)
# p_zanshi2<-data.frame(p=p_zanshi2[name,])
# rownames(p_zanshi2)<-name
# label<-p_zanshi2$p_zanshi
# label1<-p_zanshi2$p_zanshi
# label[which(label1>=0.05)]<-"ns"
# label[which(label1<0.001)]<-"***"
# label[which(label1<0.01&label1>=0.001)]<-"**"
# label[which(label1<0.05&label1>=0.01)]<-"*"
library(ggpubr)
bbb<-ggplot(mianyi_huatu,aes(y=scores,x=factor(type),fill=class))+
  geom_boxplot(outlier.size=1,outlier.shape = NA)+
  theme_classic()+
  xlab("")+ylab("")+guides()+#ylim(0,0.5)+
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+stat_compare_means(method = "kruskal.test",label = "p.signif")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  #  

bbb
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\KIRC_cluster_scrna_plot.pdf",bbb,width = 10,height = 10)
#KIRP####

#1.3.2  ####
# 
library(ggplot2)
library(ggthemes)
class<-rep(c("cluster1","cluster2"),c(nrow(cluster1_ciber),nrow(cluster2_ciber)))
mianyi_huatu<-c()
for (i in 1:nrow(ciber_zong)) {
  b<-as.numeric(ciber_zong[i,])
  tonglu1<-rep(rownames(ciber_zong)[i],ncol(ciber_zong))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)

class(mianyi_huatu)
table(mianyi_huatu$type)
mianyi_huatu$scores<-as.numeric(mianyi_huatu$scores)
mianyi_huatu$class<-as.factor(mianyi_huatu$class)

p_zanshi2<-data.frame(p_zanshi)
rownames(p_zanshi2)<-rownames(ciber_zong)
# p_zanshi2<-data.frame(p=p_zanshi2[name,])
# rownames(p_zanshi2)<-name
# label<-p_zanshi2$p_zanshi
# label1<-p_zanshi2$p_zanshi
# label[which(label1>=0.05)]<-"ns"
# label[which(label1<0.001)]<-"***"
# label[which(label1<0.01&label1>=0.001)]<-"**"
# label[which(label1<0.05&label1>=0.01)]<-"*"
library(ggpubr)
bbb<-ggplot(mianyi_huatu,aes(y=scores,x=factor(type),fill=class))+
  geom_boxplot(outlier.size=1,outlier.shape = NA)+
  theme_classic()+
  xlab("")+ylab("")+guides()+ylim(0,0.5)+
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+stat_compare_means(method = "wilcox.test",label = "p.signif")
# stat_compare_means(aes(group = class),method = "wilcox.test",label = "p.signif")
# geom_text(aes(x =1, y=0.5, label = label[1], group = type), size = 4) +
#   geom_text(aes(x =2, y=0.5, label = label[2], group = type), size = 4)+
#   geom_text(aes(x =3, y=0.5, label = label[3], group = type), size = 4)+
#   geom_text(aes(x =4, y=0.5, label = label[4], group = type), size = 4)+
#   geom_text(aes(x =5, y=0.5, label = label[5], group = type), size = 4)+
#   geom_text(aes(x =6, y=0.5, label = label[6], group = type), size = 4)+
#   geom_text(aes(x =7, y=0.5, label = label[7], group = type), size = 4)
bbb
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\KIRP_cluster_scrna_plot.pdf",bbb,width = 10,height = 10)
#LIHC####

#1.3.2  ####
# 
library(ggplot2)
library(ggthemes)
class<-rep(c("cluster1","cluster2"),c(nrow(cluster1_ciber),nrow(cluster2_ciber)))
mianyi_huatu<-c()
for (i in 1:nrow(ciber_zong)) {
  b<-as.numeric(ciber_zong[i,])
  tonglu1<-rep(rownames(ciber_zong)[i],ncol(ciber_zong))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)

class(mianyi_huatu)
table(mianyi_huatu$type)
mianyi_huatu$scores<-as.numeric(mianyi_huatu$scores)
mianyi_huatu$class<-as.factor(mianyi_huatu$class)

p_zanshi2<-data.frame(p_zanshi)
rownames(p_zanshi2)<-rownames(ciber_zong)
mianyi_huatu <- mianyi_huatu %>%
  filter(mianyi_huatu$scores < quantile(mianyi_huatu$scores, 0.8)) # 95%
# p_zanshi2<-data.frame(p=p_zanshi2[name,])
# rownames(p_zanshi2)<-name
# label<-p_zanshi2$p_zanshi
# label1<-p_zanshi2$p_zanshi
# label[which(label1>=0.05)]<-"ns"
# label[which(label1<0.001)]<-"***"
# label[which(label1<0.01&label1>=0.001)]<-"**"
# label[which(label1<0.05&label1>=0.01)]<-"*"
library(ggpubr)
bbb<-ggplot(mianyi_huatu,aes(y=scores,x=factor(type),fill=class))+
  geom_boxplot(outlier.size=1,outlier.shape = NA)+
  theme_classic()+
  xlab("")+ylab("")+guides()+ylim(0,0.5)+
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+stat_compare_means(method = "wilcox.test",label = "p.signif")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  #  
# stat_compare_means(aes(group = class),method = "wilcox.test",label = "p.signif")
# geom_text(aes(x =1, y=0.5, label = label[1], group = type), size = 4) +
#   geom_text(aes(x =2, y=0.5, label = label[2], group = type), size = 4)+
#   geom_text(aes(x =3, y=0.5, label = label[3], group = type), size = 4)+
#   geom_text(aes(x =4, y=0.5, label = label[4], group = type), size = 4)+
#   geom_text(aes(x =5, y=0.5, label = label[5], group = type), size = 4)+
#   geom_text(aes(x =6, y=0.5, label = label[6], group = type), size = 4)+
#   geom_text(aes(x =7, y=0.5, label = label[7], group = type), size = 4)
bbb
ggsave("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\LIHC_cluster_scrna_plot.pdf",bbb,width = 10,height = 10)
#LUNG####

#1.3.2  ####
# 
library(ggplot2)
library(ggthemes)
class<-rep(c("cluster1","cluster2"),c(nrow(cluster1_ciber),nrow(cluster2_ciber)))
mianyi_huatu<-c()
for (i in 1:nrow(ciber_zong)) {
  b<-as.numeric(ciber_zong[i,])
  tonglu1<-rep(rownames(ciber_zong)[i],ncol(ciber_zong))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)

class(mianyi_huatu)
table(mianyi_huatu$type)
mianyi_huatu$scores<-as.numeric(mianyi_huatu$scores)
mianyi_huatu$class<-as.factor(mianyi_huatu$class)

p_zanshi2<-data.frame(p_zanshi)
rownames(p_zanshi2)<-rownames(ciber_zong)
# p_zanshi2<-data.frame(p=p_zanshi2[name,])
# rownames(p_zanshi2)<-name
# label<-p_zanshi2$p_zanshi
# label1<-p_zanshi2$p_zanshi
# label[which(label1>=0.05)]<-"ns"
# label[which(label1<0.001)]<-"***"
# label[which(label1<0.01&label1>=0.001)]<-"**"
# label[which(label1<0.05&label1>=0.01)]<-"*"
library(ggpubr)
bbb<-ggplot(mianyi_huatu,aes(y=scores,x=factor(type),fill=class))+
  geom_boxplot(outlier.size=1,outlier.shape = NA)+
  theme_classic()+
  xlab("")+ylab("")+guides()+ylim(0,0.5)+
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+stat_compare_means(method = "wilcox.test",label = "p.signif")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  #  
# stat_compare_means(aes(group = class),method = "wilcox.test",label = "p.signif")
# geom_text(aes(x =1, y=0.5, label = label[1], group = type), size = 4) +
#   geom_text(aes(x =2, y=0.5, label = label[2], group = type), size = 4)+
#   geom_text(aes(x =3, y=0.5, label = label[3], group = type), size = 4)+
#   geom_text(aes(x =4, y=0.5, label = label[4], group = type), size = 4)+
#   geom_text(aes(x =5, y=0.5, label = label[5], group = type), size = 4)+
#   geom_text(aes(x =6, y=0.5, label = label[6], group = type), size = 4)+
#   geom_text(aes(x =7, y=0.5, label = label[7], group = type), size = 4)
bbb
ggsave("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\NSCLC_cluster_scrna_plot.pdf",bbb,width = 10,height = 10)

#PAAD####

#1.3.2  ####
# 
library(ggplot2)
library(ggthemes)
class<-rep(c("cluster1","cluster2"),c(nrow(cluster1_ciber),nrow(cluster2_ciber)))
mianyi_huatu<-c()
for (i in 1:nrow(ciber_zong)) {
  b<-as.numeric(ciber_zong[i,])
  tonglu1<-rep(rownames(ciber_zong)[i],ncol(ciber_zong))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)

class(mianyi_huatu)
table(mianyi_huatu$type)
mianyi_huatu$scores<-as.numeric(mianyi_huatu$scores)
mianyi_huatu$class<-as.factor(mianyi_huatu$class)

p_zanshi2<-data.frame(p_zanshi)
rownames(p_zanshi2)<-rownames(ciber_zong)
# p_zanshi2<-data.frame(p=p_zanshi2[name,])
# rownames(p_zanshi2)<-name
# label<-p_zanshi2$p_zanshi
# label1<-p_zanshi2$p_zanshi
# label[which(label1>=0.05)]<-"ns"
# label[which(label1<0.001)]<-"***"
# label[which(label1<0.01&label1>=0.001)]<-"**"
# label[which(label1<0.05&label1>=0.01)]<-"*"
library(ggpubr)
bbb<-ggplot(mianyi_huatu,aes(y=scores,x=factor(type),fill=class))+
  geom_boxplot(outlier.size=1,outlier.shape = NA)+
  theme_classic()+
  xlab("")+ylab("")+guides()+ylim(0,0.5)+
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+stat_compare_means(method = "wilcox.test",label = "p.signif")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  #  
# stat_compare_means(aes(group = class),method = "wilcox.test",label = "p.signif")
# geom_text(aes(x =1, y=0.5, label = label[1], group = type), size = 4) +
#   geom_text(aes(x =2, y=0.5, label = label[2], group = type), size = 4)+
#   geom_text(aes(x =3, y=0.5, label = label[3], group = type), size = 4)+
#   geom_text(aes(x =4, y=0.5, label = label[4], group = type), size = 4)+
#   geom_text(aes(x =5, y=0.5, label = label[5], group = type), size = 4)+
#   geom_text(aes(x =6, y=0.5, label = label[6], group = type), size = 4)+
#   geom_text(aes(x =7, y=0.5, label = label[7], group = type), size = 4)
bbb
ggsave("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\PAAD_cluster_scrna_plot.pdf",bbb,width = 10,height = 10)
#THCA####

#1.3.2  ####
# 
library(ggplot2)
library(ggthemes)
class<-rep(c("cluster1","cluster2"),c(nrow(cluster1_ciber),nrow(cluster2_ciber)))
mianyi_huatu<-c()
for (i in 1:nrow(ciber_zong)) {
  b<-as.numeric(ciber_zong[i,])
  tonglu1<-rep(rownames(ciber_zong)[i],ncol(ciber_zong))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)

class(mianyi_huatu)
table(mianyi_huatu$type)
mianyi_huatu$scores<-as.numeric(mianyi_huatu$scores)
mianyi_huatu$class<-as.factor(mianyi_huatu$class)

p_zanshi2<-data.frame(p_zanshi)
rownames(p_zanshi2)<-rownames(ciber_zong)
# p_zanshi2<-data.frame(p=p_zanshi2[name,])
# rownames(p_zanshi2)<-name
# label<-p_zanshi2$p_zanshi
# label1<-p_zanshi2$p_zanshi
# label[which(label1>=0.05)]<-"ns"
# label[which(label1<0.001)]<-"***"
# label[which(label1<0.01&label1>=0.001)]<-"**"
# label[which(label1<0.05&label1>=0.01)]<-"*"
library(ggpubr)
bbb<-ggplot(mianyi_huatu,aes(y=scores,x=factor(type),fill=class))+
  geom_boxplot(outlier.size=1,outlier.shape = NA)+
  theme_classic()+
  xlab("")+ylab("")+guides()+ylim(0,0.5)+
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+stat_compare_means(method = "wilcox.test",label = "p.signif")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  #  
# stat_compare_means(aes(group = class),method = "wilcox.test",label = "p.signif")
# geom_text(aes(x =1, y=0.5, label = label[1], group = type), size = 4) +
#   geom_text(aes(x =2, y=0.5, label = label[2], group = type), size = 4)+
#   geom_text(aes(x =3, y=0.5, label = label[3], group = type), size = 4)+
#   geom_text(aes(x =4, y=0.5, label = label[4], group = type), size = 4)+
#   geom_text(aes(x =5, y=0.5, label = label[5], group = type), size = 4)+
#   geom_text(aes(x =6, y=0.5, label = label[6], group = type), size = 4)+
#   geom_text(aes(x =7, y=0.5, label = label[7], group = type), size = 4)
bbb
ggsave("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\THCA_cluster_scrna_plot.pdf",bbb,width = 10,height = 10)

#UCEC####

#1.3.2  ####
# 
library(ggplot2)
library(ggthemes)
class<-rep(c("cluster1","cluster2"),c(nrow(cluster1_ciber),nrow(cluster2_ciber)))
mianyi_huatu<-c()
for (i in 1:nrow(ciber_zong)) {
  b<-as.numeric(ciber_zong[i,])
  tonglu1<-rep(rownames(ciber_zong)[i],ncol(ciber_zong))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)

class(mianyi_huatu)
table(mianyi_huatu$type)
mianyi_huatu$scores<-as.numeric(mianyi_huatu$scores)
mianyi_huatu$class<-as.factor(mianyi_huatu$class)

p_zanshi2<-data.frame(p_zanshi)
rownames(p_zanshi2)<-rownames(ciber_zong)
# p_zanshi2<-data.frame(p=p_zanshi2[name,])
# rownames(p_zanshi2)<-name
# label<-p_zanshi2$p_zanshi
# label1<-p_zanshi2$p_zanshi
# label[which(label1>=0.05)]<-"ns"
# label[which(label1<0.001)]<-"***"
# label[which(label1<0.01&label1>=0.001)]<-"**"
# label[which(label1<0.05&label1>=0.01)]<-"*"
library(ggpubr)
bbb<-ggplot(mianyi_huatu,aes(y=scores,x=factor(type),fill=class))+
  geom_boxplot(outlier.size=1,outlier.shape = NA)+
  theme_classic()+
  xlab("")+ylab("")+guides()+ylim(0,0.5)+
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+stat_compare_means(method = "wilcox.test",label = "p.signif")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # 
# stat_compare_means(aes(group = class),method = "wilcox.test",label = "p.signif")
# geom_text(aes(x =1, y=0.5, label = label[1], group = type), size = 4) +
#   geom_text(aes(x =2, y=0.5, label = label[2], group = type), size = 4)+
#   geom_text(aes(x =3, y=0.5, label = label[3], group = type), size = 4)+
#   geom_text(aes(x =4, y=0.5, label = label[4], group = type), size = 4)+
#   geom_text(aes(x =5, y=0.5, label = label[5], group = type), size = 4)+
#   geom_text(aes(x =6, y=0.5, label = label[6], group = type), size = 4)+
#   geom_text(aes(x =7, y=0.5, label = label[7], group = type), size = 4)
bbb
ggsave("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\UCEC_cluster_scrna_plot.pdf",bbb,width = 10,height = 10)

