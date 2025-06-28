#BRCA####

#1.8 emt####
#1.8.1####
library(genefilter)
#BiocManager::install("GSVA")
library(GSVA)
library(Biobase)
library(stringr)
library(openxlsx)
gene_set<- read.xlsx('E:\\1-s2.0-S2211124716317090-mmc3.xlsx')##
colnames(gene_set)<-c("Metagene","Cell.type","Immunity")
colnames(gene_set)<-gene_set[1,]
gene_set<-gene_set[-1,]
gene_set<-gene_set[, 1:2]#
head(gene_set)
#
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
#Estimates GSVA enrichment scores.
cluster1_exp<-expr[,rownames(clusters)[which(clusters$clusters=="1")]]
cluster2_exp<-expr[,intersect(rownames(clusters)[which(clusters$clusters=="2")],colnames(expr))]
str(cluster1_exp)
#cluster12
gsva_cluster1<- gsva(as.matrix(cluster1_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_cluster2<- gsva(as.matrix(cluster2_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
# cluster3_exp<-cluster3_exp[,-1]
# gsva_cluster3<- gsva(as.matrix(cluster3_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

###
gsva_zong<-cbind( gsva_cluster2,gsva_cluster1)
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))
for (j in 1:28) {
  T_result<-kruskal.test(as.numeric(gsva_zong[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_mianyi_cluster<-p_zanshi
#fdr_zanshi<-p.adjust(p_zanshi,method="fdr")
de_mianyi<-as.data.frame(gsva_zong[which(p_zanshi<0.05),])#22
#

class<-rep(c("cluster2","cluster1"),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))

mianyi_huatu<-c()
for (i in 1:nrow(de_mianyi)) {
  b<-as.numeric(de_mianyi[i,])
  tonglu1<-rep(rownames(de_mianyi)[i],ncol(de_mianyi))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)

library(ggplot2)
library(ggthemes)
class(mianyi_huatu)
table(mianyi_huatu$type)
mianyi_huatu$scores<-as.numeric(mianyi_huatu$scores)
mianyi_huatu$class<-as.factor(mianyi_huatu$class)
table(mianyi_huatu$type)

p_zanshi2<-data.frame(p=p_zanshi[which(p_zanshi<0.05)])
rownames(p_zanshi2)<-rownames(de_mianyi)

comlabel<-p_zanshi2$p
comlabel1<-p_zanshi2$p
comlabel[which(comlabel1<0.001)]<-"***"
comlabel[which(comlabel1<0.01&comlabel1>=0.001)]<-"**"
comlabel[which(comlabel1<0.05&comlabel1>=0.01)]<-"*"

comlabel<-rev(comlabel)
# 
comparisons <- combn(unique(mianyi_huatu$class), 2, simplify = FALSE)
aaa<-ggplot(mianyi_huatu,aes(x=scores,y=factor(type),fill=class))+
  geom_boxplot(outlier.size=1)+theme_classic()+scale_y_discrete(position = "right",limits = rev(levels(factor(mianyi_huatu$type))))+xlab("")+ylab("")+guides()+
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+#stat_compare_means(method = "kruskal.test",label = "p.signif",comparisons = comparisons)
  geom_text(aes(y =1, x=-0.6, label = comlabel[1], group = type), size = 4) +
  geom_text(aes(y =2, x=-0.6, label = comlabel[2], group = type), size = 4)+
  geom_text(aes(y =3, x=-0.6, label = comlabel[3], group = type), size = 4)+
  geom_text(aes(y =4, x=-0.6, label = comlabel[4], group = type), size = 4)+
  geom_text(aes(y =5, x=-0.6, label = comlabel[5], group = type), size = 4)+
  geom_text(aes(y =6, x=-0.6, label = comlabel[6], group = type), size = 4)+
  geom_text(aes(y =7, x=-0.6, label = comlabel[7], group = type), size = 4)+
  geom_text(aes(y =8, x=-0.6, label = comlabel[8], group = type), size = 4)+
  geom_text(aes(y =9, x=-0.6, label = comlabel[9], group = type), size = 4)+
  geom_text(aes(y =10, x=-0.6, label = comlabel[10], group = type), size = 4)+
  geom_text(aes(y =11, x=-0.6, label = comlabel[11], group = type), size = 4)+
  geom_text(aes(y =12, x=-0.6, label = comlabel[12], group = type), size = 4)+
  geom_text(aes(y =13, x=-0.6, label = comlabel[13], group = type), size = 4)+
  geom_text(aes(y =14, x=-0.6, label = comlabel[14], group = type), size = 4)+
  geom_text(aes(y =15, x=-0.6, label = comlabel[15], group = type), size = 4)+
  geom_text(aes(y =16, x=-0.6, label = comlabel[16], group = type), size = 4)+
  geom_text(aes(y =17, x=-0.6, label = comlabel[17], group = type), size = 4)+
  geom_text(aes(y =18, x=-0.6, label = comlabel[18], group = type), size = 4)+
  geom_text(aes(y =19, x=-0.6, label = comlabel[19], group = type), size = 4)+
  geom_text(aes(y =20, x=-0.6, label = comlabel[20], group = type), size = 4)+
  geom_text(aes(y =21, x=-0.6, label = comlabel[21], group = type), size = 4)+
  geom_text(aes(y =22, x=-0.6, label = comlabel[22], group = type), size = 4)+
  geom_text(aes(y =23, x=-0.6, label = comlabel[23], group = type), size = 4)+
  geom_text(aes(y =24, x=-0.6, label = comlabel[24], group = type), size = 4)+
  geom_text(aes(y =25, x=-0.6, label = comlabel[25], group = type), size = 4)+
  geom_text(aes(y =26, x=-0.6, label = comlabel[26], group = type), size = 4)+
  geom_text(aes(y =27, x=-0.6, label = comlabel[27], group = type), size = 4)+
  geom_text(aes(y =28, x=-0.6, label = comlabel[28], group = type), size = 4)
# +
# geom_text(aes(y =23, x=0, label = label[23], group = type), size = 4)+
# geom_text(aes(y =24, x=0, label = label[24], group = type), size = 4)+
# geom_text(aes(y =25, x=0, label = label[25], group = type), size = 4)+
# geom_text(aes(y =26, x=0, label = label[26], group = type), size = 4)+
# geom_text(aes(y =27, x=0, label = label[27], group = type), size = 4)+
# geom_text(aes(y =28, x=0, label = label[28], group = type), size = 4)
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\BRCA_ggssea_plot.pdf",aaa,width = 10,height = 10)



library(ggsci)
library(tidyr)
library(ggpubr)

#1.8.2 emt####
#BiocManager::install("KEGGREST")
library(KEGGREST)
#hsa01212:
# KEGG
# path <- keggGet("hsa01212")
# # 
# gene.info <- path[[1]]$GENE
# # gene symbolEntrez ID
# genes <- unlist(lapply(gene.info,function(x) strsplit(x,";")))#
# gene.symbol <- genes[1:length(genes)%%3 == 2]#
#gene.id <- genes[1:length(genes)%%3 == 1]
#
gene_set<- read.xlsx("E:\\pathway2.xlsx")
gene_set<-gene_set[-which(gene_set$pathway=="fatty acid metabolism"),]
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
#Estimates GSVA enrichment scores.
cluster1_exp<-cluster1_exp[,-1]
cluster2_exp<-cluster2_exp[,-1]

gsva_cluster1<- gsva(as.matrix(cluster1_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_cluster2<- gsva(as.matrix(cluster2_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
#gsva_cluster3<- gsva(as.matrix(cluster3_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

###
gsva_zong<-cbind(gsva_cluster2,gsva_cluster1)
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))

for (j in 1:nrow(gsva_zong)) {
  T_result<-kruskal.test(as.numeric(gsva_zong[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_zanshi
#fdr_zanshi<-p.adjust(p_zanshi,method="fdr")
de_mianyi<-as.data.frame(gsva_zong[which(p_zanshi<0.05),])#1
dim(de_mianyi)
#

class<-rep(c("cluster2","cluster1"),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))

mianyi_huatu<-c()
for (i in 1:nrow(de_mianyi)) {
  b<-as.numeric(de_mianyi[i,])
  tonglu1<-rep(rownames(de_mianyi)[i],ncol(de_mianyi))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}

#1
mianyi_huatu<-data.frame(scores=de_mianyi$`gsva_zong[which(p_zanshi < 0.05), ]`,type="Cell migration",class=class)


colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)
library(ggsignif)

library(ggsci)
library(tidyr)
library(ggpubr)
class(mianyi_huatu)
table(mianyi_huatu$type)
mianyi_huatu$scores<-as.numeric(mianyi_huatu$scores)
mianyi_huatu$class<-as.factor(mianyi_huatu$class)
###
p <- ggviolin(mianyi_huatu, x = "class", # xsurstat
              y = "scores", # y
              fill = "class", #surstat
              facet.by = "type",# variable
              alpha = 1,width = 0.5,legend = "right",legend.title = "Dead",
              ylab="Normalized Expression",  xlab=FALSE,
              font.y = 15,x.text.angle = 45, y.text.angle = 90,
              font.tickslab = c(15,"plain","black"),
              add = "boxplot", 
              add.params = list(fill = "white", width = 0.1,linetype = 1)) 
p <- p + stat_compare_means(method = "wilcox.test", label = "p.format",
                            label.x.npc ="left", size = 5) +
  scale_fill_manual(values=c("#80D1C8","#FFD4A9", "pink"))
p




# aaa<-ggplot(mianyi_huatu,aes(y=scores,x=factor(type),fill=class))+geom_violin()+
#   geom_boxplot(width=0.1,position=position_dodge(0.9),outlier.colour = NA)+
#   theme_classic()+xlab("")+ylab("Scores")+guides()+
#   scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+ stat_compare_means(aes(group = class), label = "p.signif")
# aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\BRCA_zhuanyi_plot.pdf",p,width = 10,height = 10)

#COAD####

#1.8 emt####
#1.8.1####
library(genefilter)
#BiocManager::install("GSVA")
library(GSVA)
library(Biobase)
library(stringr)
library(openxlsx)
gene_set<- read.xlsx('E:\\1-s2.0-S2211124716317090-mmc3.xlsx')##
colnames(gene_set)<-c("Metagene","Cell.type","Immunity")
colnames(gene_set)<-gene_set[1,]
gene_set<-gene_set[-1,]
gene_set<-gene_set[, 1:2]#
head(gene_set)
#
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
#Estimates GSVA enrichment scores.
cluster1_exp<-expr[,rownames(clusters)[which(clusters$clusters=="1")]]
cluster2_exp<-expr[,intersect(rownames(clusters)[which(clusters$clusters=="2")],colnames(expr))]
str(cluster1_exp)
#cluster12
gsva_cluster1<- gsva(as.matrix(cluster1_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_cluster2<- gsva(as.matrix(cluster2_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
# cluster3_exp<-cluster3_exp[,-1]
# gsva_cluster3<- gsva(as.matrix(cluster3_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

###
gsva_zong<-cbind( gsva_cluster2,gsva_cluster1)
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))
for (j in 1:28) {
  T_result<-kruskal.test(as.numeric(gsva_zong[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_mianyi_cluster<-p_zanshi
#fdr_zanshi<-p.adjust(p_zanshi,method="fdr")
de_mianyi<-as.data.frame(gsva_zong[which(p_zanshi<0.05),])#21
#

class<-rep(c("cluster2","cluster1"),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))

mianyi_huatu<-c()
for (i in 1:nrow(de_mianyi)) {
  b<-as.numeric(de_mianyi[i,])
  tonglu1<-rep(rownames(de_mianyi)[i],ncol(de_mianyi))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)

library(ggplot2)
library(ggthemes)
class(mianyi_huatu)
table(mianyi_huatu$type)
mianyi_huatu$scores<-as.numeric(mianyi_huatu$scores)
mianyi_huatu$class<-as.factor(mianyi_huatu$class)
table(mianyi_huatu$type)

p_zanshi2<-data.frame(p=p_zanshi[which(p_zanshi<0.05)])
rownames(p_zanshi2)<-rownames(de_mianyi)

comlabel<-p_zanshi2$p
comlabel1<-p_zanshi2$p
comlabel[which(comlabel1<0.001)]<-"***"
comlabel[which(comlabel1<0.01&comlabel1>=0.001)]<-"**"
comlabel[which(comlabel1<0.05&comlabel1>=0.01)]<-"*"

comlabel<-rev(comlabel)
# 
comparisons <- combn(unique(mianyi_huatu$class), 2, simplify = FALSE)
aaa<-ggplot(mianyi_huatu,aes(x=scores,y=factor(type),fill=class))+
  geom_boxplot(outlier.size=1)+theme_classic()+scale_y_discrete(position = "right",limits = rev(levels(factor(mianyi_huatu$type))))+xlab("")+ylab("")+guides()+
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+#stat_compare_means(method = "kruskal.test",label = "p.signif",comparisons = comparisons)
  geom_text(aes(y =1, x=-0.6, label = comlabel[1], group = type), size = 4) +
  geom_text(aes(y =2, x=-0.6, label = comlabel[2], group = type), size = 4)+
  geom_text(aes(y =3, x=-0.6, label = comlabel[3], group = type), size = 4)+
  geom_text(aes(y =4, x=-0.6, label = comlabel[4], group = type), size = 4)+
  geom_text(aes(y =5, x=-0.6, label = comlabel[5], group = type), size = 4)+
  geom_text(aes(y =6, x=-0.6, label = comlabel[6], group = type), size = 4)+
  geom_text(aes(y =7, x=-0.6, label = comlabel[7], group = type), size = 4)+
  geom_text(aes(y =8, x=-0.6, label = comlabel[8], group = type), size = 4)+
  geom_text(aes(y =9, x=-0.6, label = comlabel[9], group = type), size = 4)+
  geom_text(aes(y =10, x=-0.6, label = comlabel[10], group = type), size = 4)+
  geom_text(aes(y =11, x=-0.6, label = comlabel[11], group = type), size = 4)+
  geom_text(aes(y =12, x=-0.6, label = comlabel[12], group = type), size = 4)+
  geom_text(aes(y =13, x=-0.6, label = comlabel[13], group = type), size = 4)+
  geom_text(aes(y =14, x=-0.6, label = comlabel[14], group = type), size = 4)+
  geom_text(aes(y =15, x=-0.6, label = comlabel[15], group = type), size = 4)+
  geom_text(aes(y =16, x=-0.6, label = comlabel[16], group = type), size = 4)+
  geom_text(aes(y =17, x=-0.6, label = comlabel[16], group = type), size = 4)+
  geom_text(aes(y =18, x=-0.6, label = comlabel[16], group = type), size = 4)+
  geom_text(aes(y =19, x=-0.6, label = comlabel[16], group = type), size = 4)+
  geom_text(aes(y =20, x=-0.6, label = comlabel[16], group = type), size = 4)+
  geom_text(aes(y =21, x=-0.6, label = comlabel[16], group = type), size = 4)
# geom_text(aes(y =17, x=0, label = label[17], group = type), size = 4)+
# geom_text(aes(y =18, x=0, label = label[18], group = type), size = 4)+
# geom_text(aes(y =19, x=0, label = label[19], group = type), size = 4)+
# geom_text(aes(y =20, x=0, label = label[20], group = type), size = 4)+
# geom_text(aes(y =21, x=0, label = label[21], group = type), size = 4)+
# geom_text(aes(y =22, x=0, label = label[22], group = type), size = 4)+
# geom_text(aes(y =23, x=0, label = label[23], group = type), size = 4)+
# geom_text(aes(y =24, x=0, label = label[24], group = type), size = 4)+
# geom_text(aes(y =25, x=0, label = label[25], group = type), size = 4)+
# geom_text(aes(y =26, x=0, label = label[26], group = type), size = 4)+
# geom_text(aes(y =27, x=0, label = label[27], group = type), size = 4)+
# geom_text(aes(y =28, x=0, label = label[28], group = type), size = 4)
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\COAD_ggssea_plot.pdf",aaa,width = 10,height = 10)



library(ggsci)
library(tidyr)
library(ggpubr)

#1.8.2 emt####
#BiocManager::install("KEGGREST")
library(KEGGREST)
#hsa01212:
# KEGG
# path <- keggGet("hsa01212")
# # 
# gene.info <- path[[1]]$GENE
# # gene symbolEntrez ID
# genes <- unlist(lapply(gene.info,function(x) strsplit(x,";")))#
# gene.symbol <- genes[1:length(genes)%%3 == 2]#
#gene.id <- genes[1:length(genes)%%3 == 1]
#
gene_set<- read.xlsx("E:\\pathway2.xlsx")
gene_set<-gene_set[-which(gene_set$pathway=="fatty acid metabolism"),]
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
#Estimates GSVA enrichment scores.

gsva_cluster1<- gsva(as.matrix(cluster1_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_cluster2<- gsva(as.matrix(cluster2_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
#gsva_cluster3<- gsva(as.matrix(cluster3_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

###
gsva_zong<-cbind(gsva_cluster2,gsva_cluster1)
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))

for (j in 1:nrow(gsva_zong)) {
  T_result<-kruskal.test(as.numeric(gsva_zong[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_zanshi
#fdr_zanshi<-p.adjust(p_zanshi,method="fdr")
de_mianyi<-as.data.frame(gsva_zong[which(p_zanshi<0.05),])#2
#

class<-rep(c("cluster2","cluster1"),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))

mianyi_huatu<-c()
for (i in 1:nrow(de_mianyi)) {
  b<-as.numeric(de_mianyi[i,])
  tonglu1<-rep(rownames(de_mianyi)[i],ncol(de_mianyi))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)
library(ggsignif)

library(ggsci)
library(tidyr)
library(ggpubr)
class(mianyi_huatu)
table(mianyi_huatu$type)
mianyi_huatu$scores<-as.numeric(mianyi_huatu$scores)
mianyi_huatu$class<-as.factor(mianyi_huatu$class)
###
p <- ggviolin(mianyi_huatu, x = "class", # xsurstat
              y = "scores", # y
              fill = "class", #surstat
              facet.by = "type",# variable
              alpha = 1,width = 0.5,legend = "right",legend.title = "Dead",
              ylab="Normalized Expression",  xlab=FALSE,
              font.y = 15,x.text.angle = 45, y.text.angle = 90,
              font.tickslab = c(15,"plain","black"),
              add = "boxplot", 
              add.params = list(fill = "white", width = 0.1,linetype = 1)) 
p <- p + stat_compare_means(method = "wilcox.test", label = "p.format",
                            label.x.npc ="left", size = 5) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))
p




# aaa<-ggplot(mianyi_huatu,aes(y=scores,x=factor(type),fill=class))+geom_violin()+
#   geom_boxplot(width=0.1,position=position_dodge(0.9),outlier.colour = NA)+
#   theme_classic()+xlab("")+ylab("Scores")+guides()+
#   scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+ stat_compare_means(aes(group = class), label = "p.signif")
# aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\COAD_zhuanyi_plot.pdf",p,width = 10,height = 10)

#ESCA####

#1.8 emt####
#1.8.1####
library(genefilter)
#BiocManager::install("GSVA")
library(GSVA)
library(Biobase)
library(stringr)
library(openxlsx)
gene_set<- read.xlsx('E:\\1-s2.0-S2211124716317090-mmc3.xlsx')##
colnames(gene_set)<-c("Metagene","Cell.type","Immunity")
colnames(gene_set)<-gene_set[1,]
gene_set<-gene_set[-1,]
gene_set<-gene_set[, 1:2]#
head(gene_set)
#
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
#Estimates GSVA enrichment scores.
cluster1_exp<-expr[,rownames(clusters)[which(clusters$clusters=="1")]]
cluster2_exp<-expr[,intersect(rownames(clusters)[which(clusters$clusters=="2")],colnames(expr))]
str(cluster1_exp)
#cluster12
gsva_cluster1<- gsva(as.matrix(cluster1_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_cluster2<- gsva(as.matrix(cluster2_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
# cluster3_exp<-cluster3_exp[,-1]
# gsva_cluster3<- gsva(as.matrix(cluster3_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

###
gsva_zong<-cbind( gsva_cluster2,gsva_cluster1)
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))
for (j in 1:28) {
  T_result<-kruskal.test(as.numeric(gsva_zong[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_mianyi_cluster<-p_zanshi
#fdr_zanshi<-p.adjust(p_zanshi,method="fdr")
de_mianyi<-as.data.frame(gsva_zong[which(p_zanshi<0.05),])#21
#

class<-rep(c("cluster2","cluster1"),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))

mianyi_huatu<-c()
for (i in 1:nrow(de_mianyi)) {
  b<-as.numeric(de_mianyi[i,])
  tonglu1<-rep(rownames(de_mianyi)[i],ncol(de_mianyi))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)

library(ggplot2)
library(ggthemes)
class(mianyi_huatu)
table(mianyi_huatu$type)
mianyi_huatu$scores<-as.numeric(mianyi_huatu$scores)
mianyi_huatu$class<-as.factor(mianyi_huatu$class)
table(mianyi_huatu$type)

p_zanshi2<-data.frame(p=p_zanshi[which(p_zanshi<0.05)])
rownames(p_zanshi2)<-rownames(de_mianyi)

comlabel<-p_zanshi2$p
comlabel1<-p_zanshi2$p
comlabel[which(comlabel1<0.001)]<-"***"
comlabel[which(comlabel1<0.01&comlabel1>=0.001)]<-"**"
comlabel[which(comlabel1<0.05&comlabel1>=0.01)]<-"*"

comlabel<-rev(comlabel)
# 
comparisons <- combn(unique(mianyi_huatu$class), 2, simplify = FALSE)
aaa<-ggplot(mianyi_huatu,aes(x=scores,y=factor(type),fill=class))+
  geom_boxplot(outlier.size=1)+theme_classic()+scale_y_discrete(position = "right",limits = rev(levels(factor(mianyi_huatu$type))))+xlab("")+ylab("")+guides()+
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+#stat_compare_means(method = "kruskal.test",label = "p.signif",comparisons = comparisons)
  geom_text(aes(y =1, x=-0.6, label = comlabel[1], group = type), size = 4) +
  geom_text(aes(y =2, x=-0.6, label = comlabel[2], group = type), size = 4)+
  geom_text(aes(y =3, x=-0.6, label = comlabel[3], group = type), size = 4)+
  geom_text(aes(y =4, x=-0.6, label = comlabel[4], group = type), size = 4)+
  geom_text(aes(y =5, x=-0.6, label = comlabel[5], group = type), size = 4)+
  geom_text(aes(y =6, x=-0.6, label = comlabel[6], group = type), size = 4)+
  geom_text(aes(y =7, x=-0.6, label = comlabel[7], group = type), size = 4)+
  geom_text(aes(y =8, x=-0.6, label = comlabel[8], group = type), size = 4)+
  geom_text(aes(y =9, x=-0.6, label = comlabel[9], group = type), size = 4)+
  geom_text(aes(y =10, x=-0.6, label = comlabel[10], group = type), size = 4)+
  geom_text(aes(y =11, x=-0.6, label = comlabel[11], group = type), size = 4)+
  geom_text(aes(y =12, x=-0.6, label = comlabel[12], group = type), size = 4)+
  geom_text(aes(y =13, x=-0.6, label = comlabel[13], group = type), size = 4)+
  geom_text(aes(y =14, x=-0.6, label = comlabel[14], group = type), size = 4)+
  geom_text(aes(y =15, x=-0.6, label = comlabel[15], group = type), size = 4)+
  geom_text(aes(y =16, x=-0.6, label = comlabel[16], group = type), size = 4)+
  geom_text(aes(y =17, x=-0.6, label = comlabel[17], group = type), size = 4)+
  geom_text(aes(y =18, x=-0.6, label = comlabel[18], group = type), size = 4)
# geom_text(aes(y =17, x=0, label = label[17], group = type), size = 4)+
# geom_text(aes(y =18, x=0, label = label[18], group = type), size = 4)+
# geom_text(aes(y =19, x=0, label = label[19], group = type), size = 4)+
# geom_text(aes(y =20, x=0, label = label[20], group = type), size = 4)+
# geom_text(aes(y =21, x=0, label = label[21], group = type), size = 4)+
# geom_text(aes(y =22, x=0, label = label[22], group = type), size = 4)+
# geom_text(aes(y =23, x=0, label = label[23], group = type), size = 4)+
# geom_text(aes(y =24, x=0, label = label[24], group = type), size = 4)+
# geom_text(aes(y =25, x=0, label = label[25], group = type), size = 4)+
# geom_text(aes(y =26, x=0, label = label[26], group = type), size = 4)+
# geom_text(aes(y =27, x=0, label = label[27], group = type), size = 4)+
# geom_text(aes(y =28, x=0, label = label[28], group = type), size = 4)
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\ESCA_ggssea_plot.pdf",aaa,width = 10,height = 10)



library(ggsci)
library(tidyr)
library(ggpubr)

#1.8.2 emt####
#BiocManager::install("KEGGREST")
library(KEGGREST)
#hsa01212:
# KEGG
# path <- keggGet("hsa01212")
# # 
# gene.info <- path[[1]]$GENE
# # gene symbolEntrez ID
# genes <- unlist(lapply(gene.info,function(x) strsplit(x,";")))#
# gene.symbol <- genes[1:length(genes)%%3 == 2]#
#gene.id <- genes[1:length(genes)%%3 == 1]
#
gene_set<- read.xlsx("E:\\pathway2.xlsx")
gene_set<-gene_set[-which(gene_set$pathway=="fatty acid metabolism"),]
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
#Estimates GSVA enrichment scores.

gsva_cluster1<- gsva(as.matrix(cluster1_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_cluster2<- gsva(as.matrix(cluster2_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
#gsva_cluster3<- gsva(as.matrix(cluster3_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

###
gsva_zong<-cbind(gsva_cluster2,gsva_cluster1)
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))

for (j in 1:nrow(gsva_zong)) {
  T_result<-kruskal.test(as.numeric(gsva_zong[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_zanshi
#fdr_zanshi<-p.adjust(p_zanshi,method="fdr")
de_mianyi<-as.data.frame(gsva_zong)#2
#

class<-rep(c("cluster2","cluster1"),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))

mianyi_huatu<-c()
for (i in 1:nrow(de_mianyi)) {
  b<-as.numeric(de_mianyi[i,])
  tonglu1<-rep(rownames(de_mianyi)[i],ncol(de_mianyi))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)
library(ggsignif)

library(ggsci)
library(tidyr)
library(ggpubr)
class(mianyi_huatu)
table(mianyi_huatu$type)
mianyi_huatu$scores<-as.numeric(mianyi_huatu$scores)
mianyi_huatu$class<-as.factor(mianyi_huatu$class)
###
p <- ggviolin(mianyi_huatu, x = "class", # xsurstat
              y = "scores", # y
              fill = "class", #surstat
              facet.by = "type",# variable
              alpha = 1,width = 0.5,legend = "right",legend.title = "Dead",
              ylab="Normalized Expression",  xlab=FALSE,
              font.y = 15,x.text.angle = 45, y.text.angle = 90,
              font.tickslab = c(15,"plain","black"),
              add = "boxplot", 
              add.params = list(fill = "white", width = 0.1,linetype = 1)) 
p <- p + stat_compare_means(method = "wilcox.test", label = "p.format",
                            label.x.npc ="left", size = 5) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))
p




# aaa<-ggplot(mianyi_huatu,aes(y=scores,x=factor(type),fill=class))+geom_violin()+
#   geom_boxplot(width=0.1,position=position_dodge(0.9),outlier.colour = NA)+
#   theme_classic()+xlab("")+ylab("Scores")+guides()+
#   scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+ stat_compare_means(aes(group = class), label = "p.signif")
# aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\ESCA_zhuanyi_plot.pdf",p,width = 10,height = 10)

#HNSC####
#1.8 emt####
#1.8.1####
library(genefilter)
#BiocManager::install("GSVA")
library(GSVA)
library(Biobase)
library(stringr)
library(openxlsx)
gene_set<- read.xlsx('E:\\1-s2.0-S2211124716317090-mmc3.xlsx')##
colnames(gene_set)<-c("Metagene","Cell.type","Immunity")
colnames(gene_set)<-gene_set[1,]
gene_set<-gene_set[-1,]
gene_set<-gene_set[, 1:2]#
head(gene_set)
#
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
#Estimates GSVA enrichment scores.
cluster1_exp<-expr[,rownames(clusters)[which(clusters$clusters=="1")]]
cluster2_exp<-expr[,intersect(rownames(clusters)[which(clusters$clusters=="2")],colnames(expr))]
str(cluster1_exp)
#cluster12
gsva_cluster1<- gsva(as.matrix(cluster1_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_cluster2<- gsva(as.matrix(cluster2_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
# cluster3_exp<-cluster3_exp[,-1]
# gsva_cluster3<- gsva(as.matrix(cluster3_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

###
gsva_zong<-cbind( gsva_cluster2,gsva_cluster1)
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))
for (j in 1:28) {
  T_result<-kruskal.test(as.numeric(gsva_zong[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_mianyi_cluster<-p_zanshi
#fdr_zanshi<-p.adjust(p_zanshi,method="fdr")
de_mianyi<-as.data.frame(gsva_zong[which(p_zanshi<0.05),])#22
#

class<-rep(c("cluster2","cluster1"),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))

mianyi_huatu<-c()
for (i in 1:nrow(de_mianyi)) {
  b<-as.numeric(de_mianyi[i,])
  tonglu1<-rep(rownames(de_mianyi)[i],ncol(de_mianyi))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)

library(ggplot2)
library(ggthemes)
class(mianyi_huatu)
table(mianyi_huatu$type)
mianyi_huatu$scores<-as.numeric(mianyi_huatu$scores)
mianyi_huatu$class<-as.factor(mianyi_huatu$class)
table(mianyi_huatu$type)

p_zanshi2<-data.frame(p=p_zanshi[which(p_zanshi<0.05)])
rownames(p_zanshi2)<-rownames(de_mianyi)

comlabel<-p_zanshi2$p
comlabel1<-p_zanshi2$p
comlabel[which(comlabel1<0.001)]<-"***"
comlabel[which(comlabel1<0.01&comlabel1>=0.001)]<-"**"
comlabel[which(comlabel1<0.05&comlabel1>=0.01)]<-"*"

comlabel<-rev(comlabel)
# 
comparisons <- combn(unique(mianyi_huatu$class), 2, simplify = FALSE)
aaa<-ggplot(mianyi_huatu,aes(x=scores,y=factor(type),fill=class))+
  geom_boxplot(outlier.size=1)+theme_classic()+scale_y_discrete(position = "right",limits = rev(levels(factor(mianyi_huatu$type))))+xlab("")+ylab("")+guides()+
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+#stat_compare_means(method = "kruskal.test",label = "p.signif",comparisons = comparisons)
  geom_text(aes(y =1, x=-0.6, label = comlabel[1], group = type), size = 4) +
  geom_text(aes(y =2, x=-0.6, label = comlabel[2], group = type), size = 4)+
  geom_text(aes(y =3, x=-0.6, label = comlabel[3], group = type), size = 4)+
  geom_text(aes(y =4, x=-0.6, label = comlabel[4], group = type), size = 4)+
  geom_text(aes(y =5, x=-0.6, label = comlabel[5], group = type), size = 4)+
  geom_text(aes(y =6, x=-0.6, label = comlabel[6], group = type), size = 4)+
  geom_text(aes(y =7, x=-0.6, label = comlabel[7], group = type), size = 4)+
  geom_text(aes(y =8, x=-0.6, label = comlabel[8], group = type), size = 4)+
  geom_text(aes(y =9, x=-0.6, label = comlabel[9], group = type), size = 4)+
  geom_text(aes(y =10, x=-0.6, label = comlabel[10], group = type), size = 4)+
  geom_text(aes(y =11, x=-0.6, label = comlabel[11], group = type), size = 4)+
  geom_text(aes(y =12, x=-0.6, label = comlabel[12], group = type), size = 4)
# +
#   geom_text(aes(y =13, x=-0.6, label = comlabel[13], group = type), size = 4)+
#   geom_text(aes(y =14, x=-0.6, label = comlabel[14], group = type), size = 4)+
#   geom_text(aes(y =15, x=-0.6, label = comlabel[15], group = type), size = 4)+
#   geom_text(aes(y =16, x=-0.6, label = comlabel[16], group = type), size = 4)+
#   geom_text(aes(y =17, x=-0.6, label = comlabel[17], group = type), size = 4)+
#   geom_text(aes(y =18, x=-0.6, label = comlabel[18], group = type), size = 4)+
#   geom_text(aes(y =19, x=-0.6, label = comlabel[19], group = type), size = 4)+
#   geom_text(aes(y =20, x=-0.6, label = comlabel[20], group = type), size = 4)
# +
#   geom_text(aes(y =21, x=-0.6, label = comlabel[21], group = type), size = 4)+
#   geom_text(aes(y =22, x=-0.6, label = comlabel[22], group = type), size = 4)
# geom_text(aes(y =17, x=0, label = label[17], group = type), size = 4)+
# geom_text(aes(y =18, x=0, label = label[18], group = type), size = 4)+
# geom_text(aes(y =19, x=0, label = label[19], group = type), size = 4)+
# geom_text(aes(y =20, x=0, label = label[20], group = type), size = 4)+
# geom_text(aes(y =21, x=0, label = label[21], group = type), size = 4)+
# geom_text(aes(y =22, x=0, label = label[22], group = type), size = 4)+
# geom_text(aes(y =23, x=0, label = label[23], group = type), size = 4)+
# geom_text(aes(y =24, x=0, label = label[24], group = type), size = 4)+
# geom_text(aes(y =25, x=0, label = label[25], group = type), size = 4)+
# geom_text(aes(y =26, x=0, label = label[26], group = type), size = 4)+
# geom_text(aes(y =27, x=0, label = label[27], group = type), size = 4)+
# geom_text(aes(y =28, x=0, label = label[28], group = type), size = 4)
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_ggssea_plot.pdf",aaa,width = 10,height = 10)



library(ggsci)
library(tidyr)
library(ggpubr)

#1.8.2 emt####
#BiocManager::install("KEGGREST")
library(KEGGREST)
#hsa01212:
# KEGG
# path <- keggGet("hsa01212")
# # 
# gene.info <- path[[1]]$GENE
# # gene symbolEntrez ID
# genes <- unlist(lapply(gene.info,function(x) strsplit(x,";")))#
# gene.symbol <- genes[1:length(genes)%%3 == 2]#
#gene.id <- genes[1:length(genes)%%3 == 1]
#
gene_set<- read.xlsx("E:\\pathway2.xlsx")
gene_set<-gene_set[-which(gene_set$pathway=="fatty acid metabolism"),]
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
#Estimates GSVA enrichment scores.

gsva_cluster1<- gsva(as.matrix(cluster1_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_cluster2<- gsva(as.matrix(cluster2_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
#gsva_cluster3<- gsva(as.matrix(cluster3_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

###
gsva_zong<-cbind(gsva_cluster2,gsva_cluster1)
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))

for (j in 1:nrow(gsva_zong)) {
  T_result<-kruskal.test(as.numeric(gsva_zong[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_zanshi
#fdr_zanshi<-p.adjust(p_zanshi,method="fdr")
de_mianyi<-as.data.frame(gsva_zong)#1
#

class<-rep(c("cluster2","cluster1"),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))

mianyi_huatu<-c()
for (i in 1:nrow(de_mianyi)) {
  b<-as.numeric(de_mianyi[i,])
  tonglu1<-rep(rownames(de_mianyi)[i],ncol(de_mianyi))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}

#1
mianyi_huatu<-data.frame(scores=de_mianyi$`gsva_zong[which(p_zanshi < 0.05), ]`,type="Cell migration",class=class)
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)
library(ggsignif)

library(ggsci)
library(tidyr)
library(ggpubr)
class(mianyi_huatu)
table(mianyi_huatu$type)
mianyi_huatu$scores<-as.numeric(mianyi_huatu$scores)
mianyi_huatu$class<-as.factor(mianyi_huatu$class)
###
p <- ggviolin(mianyi_huatu, x = "class", # xsurstat
              y = "scores", # y
              fill = "class", #surstat
              facet.by = "type",# variable
              alpha = 1,width = 0.5,legend = "right",legend.title = "Dead",
              ylab="Normalized Expression",  xlab=FALSE,
              font.y = 15,x.text.angle = 45, y.text.angle = 90,
              font.tickslab = c(15,"plain","black"),
              add = "boxplot", 
              add.params = list(fill = "white", width = 0.1,linetype = 1)) 
p <- p + stat_compare_means(method = "wilcox.test", label = "p.format",
                            label.x.npc ="left", size = 5) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))
p




# aaa<-ggplot(mianyi_huatu,aes(y=scores,x=factor(type),fill=class))+geom_violin()+
#   geom_boxplot(width=0.1,position=position_dodge(0.9),outlier.colour = NA)+
#   theme_classic()+xlab("")+ylab("Scores")+guides()+
#   scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+ stat_compare_means(aes(group = class), label = "p.signif")
# aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_zhuanyi_plot.pdf",p,width = 10,height = 10)

#KIRC####

#1.8 emt####
#1.8.1####
library(genefilter)
#BiocManager::install("GSVA")
library(GSVA)
library(Biobase)
library(stringr)
library(openxlsx)
gene_set<- read.xlsx('E:\\1-s2.0-S2211124716317090-mmc3.xlsx')##
colnames(gene_set)<-c("Metagene","Cell.type","Immunity")
colnames(gene_set)<-gene_set[1,]
gene_set<-gene_set[-1,]
gene_set<-gene_set[, 1:2]#
head(gene_set)
#
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
#Estimates GSVA enrichment scores.
cluster1_exp<-expr[,rownames(clusters)[which(clusters$clusters=="1")]]
cluster2_exp<-expr[,intersect(rownames(clusters)[which(clusters$clusters=="2")],colnames(expr))]
str(cluster1_exp)
#cluster12
gsva_cluster1<- gsva(as.matrix(cluster1_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_cluster2<- gsva(as.matrix(cluster2_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
cluster3_exp<-expr[,rownames(clusters)[which(clusters$clusters=="3")]]
gsva_cluster3<- gsva(as.matrix(cluster3_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

###
gsva_zong<-cbind(gsva_cluster3, gsva_cluster2,gsva_cluster1)
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(3,2,1),c(ncol(gsva_cluster3),ncol(gsva_cluster2),ncol(gsva_cluster1)))
for (j in 1:28) {
  T_result<-kruskal.test(as.numeric(gsva_zong[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_mianyi_cluster<-p_zanshi
#fdr_zanshi<-p.adjust(p_zanshi,method="fdr")
de_mianyi<-as.data.frame(gsva_zong[which(p_zanshi<0.05),])#21
#

class<-rep(c("cluster3","cluster2","cluster1"),c(ncol(gsva_cluster3),ncol(gsva_cluster2),ncol(gsva_cluster1)))

mianyi_huatu<-c()
for (i in 1:nrow(de_mianyi)) {
  b<-as.numeric(de_mianyi[i,])
  tonglu1<-rep(rownames(de_mianyi)[i],ncol(de_mianyi))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)

library(ggplot2)
library(ggthemes)
class(mianyi_huatu)
table(mianyi_huatu$type)
mianyi_huatu$scores<-as.numeric(mianyi_huatu$scores)
mianyi_huatu$class<-as.factor(mianyi_huatu$class)
table(mianyi_huatu$type)

p_zanshi2<-data.frame(p=p_zanshi[which(p_zanshi<0.05)])
rownames(p_zanshi2)<-rownames(de_mianyi)

comlabel<-p_zanshi2$p
comlabel1<-p_zanshi2$p
comlabel[which(comlabel1<0.001)]<-"***"
comlabel[which(comlabel1<0.01&comlabel1>=0.001)]<-"**"
comlabel[which(comlabel1<0.05&comlabel1>=0.01)]<-"*"

comlabel<-rev(comlabel)
# 
comparisons <- combn(unique(mianyi_huatu$class), 2, simplify = FALSE)
aaa<-ggplot(mianyi_huatu,aes(x=scores,y=factor(type),fill=class))+
  geom_boxplot(outlier.size=1)+theme_classic()+scale_y_discrete(position = "right",limits = rev(levels(factor(mianyi_huatu$type))))+xlab("")+ylab("")+guides()+
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+#stat_compare_means(method = "kruskal.test",label = "p.signif",comparisons = comparisons)
  geom_text(aes(y =1, x=-0.6, label = comlabel[1], group = type), size = 4) +
  geom_text(aes(y =2, x=-0.6, label = comlabel[2], group = type), size = 4)+
  geom_text(aes(y =3, x=-0.6, label = comlabel[3], group = type), size = 4)+
  geom_text(aes(y =4, x=-0.6, label = comlabel[4], group = type), size = 4)+
  geom_text(aes(y =5, x=-0.6, label = comlabel[5], group = type), size = 4)+
  geom_text(aes(y =6, x=-0.6, label = comlabel[6], group = type), size = 4)+
  geom_text(aes(y =7, x=-0.6, label = comlabel[7], group = type), size = 4)+
  geom_text(aes(y =8, x=-0.6, label = comlabel[8], group = type), size = 4)+
  geom_text(aes(y =9, x=-0.6, label = comlabel[9], group = type), size = 4)+
  geom_text(aes(y =10, x=-0.6, label = comlabel[10], group = type), size = 4)+
  geom_text(aes(y =11, x=-0.6, label = comlabel[11], group = type), size = 4)+
  geom_text(aes(y =12, x=-0.6, label = comlabel[12], group = type), size = 4)+
  geom_text(aes(y =13, x=-0.6, label = comlabel[13], group = type), size = 4)+
  geom_text(aes(y =14, x=-0.6, label = comlabel[14], group = type), size = 4)+
  geom_text(aes(y =15, x=-0.6, label = comlabel[15], group = type), size = 4)+
  geom_text(aes(y =16, x=-0.6, label = comlabel[16], group = type), size = 4)+
  geom_text(aes(y =17, x=-0.6, label = comlabel[16], group = type), size = 4)+
  geom_text(aes(y =18, x=-0.6, label = comlabel[16], group = type), size = 4)+
  geom_text(aes(y =19, x=-0.6, label = comlabel[16], group = type), size = 4)+
  geom_text(aes(y =20, x=-0.6, label = comlabel[16], group = type), size = 4)+
  geom_text(aes(y =21, x=-0.6, label = comlabel[16], group = type), size = 4)
# geom_text(aes(y =17, x=0, label = label[17], group = type), size = 4)+
# geom_text(aes(y =18, x=0, label = label[18], group = type), size = 4)+
# geom_text(aes(y =19, x=0, label = label[19], group = type), size = 4)+
# geom_text(aes(y =20, x=0, label = label[20], group = type), size = 4)+
# geom_text(aes(y =21, x=0, label = label[21], group = type), size = 4)+
# geom_text(aes(y =22, x=0, label = label[22], group = type), size = 4)+
# geom_text(aes(y =23, x=0, label = label[23], group = type), size = 4)+
# geom_text(aes(y =24, x=0, label = label[24], group = type), size = 4)+
# geom_text(aes(y =25, x=0, label = label[25], group = type), size = 4)+
# geom_text(aes(y =26, x=0, label = label[26], group = type), size = 4)+
# geom_text(aes(y =27, x=0, label = label[27], group = type), size = 4)+
# geom_text(aes(y =28, x=0, label = label[28], group = type), size = 4)
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\KIRC_ggssea_plot.pdf",aaa,width = 10,height = 10)



library(ggsci)
library(tidyr)
library(ggpubr)

#1.8.2 emt####
#BiocManager::install("KEGGREST")
library(KEGGREST)
#hsa01212:
# KEGG
# path <- keggGet("hsa01212")
# # 
# gene.info <- path[[1]]$GENE
# # gene symbolEntrez ID
# genes <- unlist(lapply(gene.info,function(x) strsplit(x,";")))#
# gene.symbol <- genes[1:length(genes)%%3 == 2]#
#gene.id <- genes[1:length(genes)%%3 == 1]
#
gene_set<- read.xlsx("E:\\pathway2.xlsx")
gene_set<-gene_set[-which(gene_set$pathway=="fatty acid metabolism"),]
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
#Estimates GSVA enrichment scores.
cluster1_exp<-cluster1_exp[,-1]
cluster2_exp<-cluster2_exp[,-1]
cluster3_exp<-cluster3_exp[,-1]
gsva_cluster1<- gsva(as.matrix(cluster1_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_cluster2<- gsva(as.matrix(cluster2_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_cluster3<- gsva(as.matrix(cluster3_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

###
gsva_zong<-cbind(gsva_cluster3,gsva_cluster2,gsva_cluster1)
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(3,2,1),c(ncol(gsva_cluster3),ncol(gsva_cluster2),ncol(gsva_cluster1)))

for (j in 1:nrow(gsva_zong)) {
  T_result<-kruskal.test(as.numeric(gsva_zong[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_zanshi
#fdr_zanshi<-p.adjust(p_zanshi,method="fdr")
de_mianyi<-as.data.frame(gsva_zong[which(p_zanshi<0.05),])#2
#

class<-rep(c("cluster3","cluster2","cluster1"),c(ncol(gsva_cluster3),ncol(gsva_cluster2),ncol(gsva_cluster1)))

mianyi_huatu<-c()
for (i in 1:nrow(de_mianyi)) {
  b<-as.numeric(de_mianyi[i,])
  tonglu1<-rep(rownames(de_mianyi)[i],ncol(de_mianyi))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)
library(ggsignif)

library(ggsci)
library(tidyr)
library(ggpubr)
class(mianyi_huatu)
table(mianyi_huatu$type)
mianyi_huatu$scores<-as.numeric(mianyi_huatu$scores)
mianyi_huatu$class<-as.factor(mianyi_huatu$class)
###
p <- ggviolin(mianyi_huatu, x = "class", # xsurstat
              y = "scores", # y
              fill = "class", #surstat
              facet.by = "type",# variable
              alpha = 1,width = 0.5,legend = "right",legend.title = "Dead",
              ylab="Normalized Expression",  xlab=FALSE,
              font.y = 15,x.text.angle = 45, y.text.angle = 90,
              font.tickslab = c(15,"plain","black"),
              add = "boxplot", 
              add.params = list(fill = "white", width = 0.1,linetype = 1)) 
p <- p + stat_compare_means(method = "kruskal.test", label = "p.format",
                            label.x.npc ="left", size = 5) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))
p




# aaa<-ggplot(mianyi_huatu,aes(y=scores,x=factor(type),fill=class))+geom_violin()+
#   geom_boxplot(width=0.1,position=position_dodge(0.9),outlier.colour = NA)+
#   theme_classic()+xlab("")+ylab("Scores")+guides()+
#   scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+ stat_compare_means(aes(group = class), label = "p.signif")
# aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\KIRC_zhuanyi_plot.pdf",p,width = 10,height = 10)

#KIRP####

#1.8 emt####
#1.8.1####
library(genefilter)
#BiocManager::install("GSVA")
library(GSVA)
library(Biobase)
library(stringr)
library(openxlsx)
gene_set<- read.xlsx('E:\\1-s2.0-S2211124716317090-mmc3.xlsx')##
colnames(gene_set)<-c("Metagene","Cell.type","Immunity")
colnames(gene_set)<-gene_set[1,]
gene_set<-gene_set[-1,]
gene_set<-gene_set[, 1:2]#
head(gene_set)
#
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
#Estimates GSVA enrichment scores.
cluster1_exp<-expr[,rownames(clusters)[which(clusters$clusters=="1")]]
cluster2_exp<-expr[,intersect(rownames(clusters)[which(clusters$clusters=="2")],colnames(expr))]
str(cluster1_exp)
#cluster12
gsva_cluster1<- gsva(as.matrix(cluster1_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_cluster2<- gsva(as.matrix(cluster2_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
# cluster3_exp<-cluster3_exp[,-1]
# gsva_cluster3<- gsva(as.matrix(cluster3_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

###
gsva_zong<-cbind( gsva_cluster2,gsva_cluster1)
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))
for (j in 1:28) {
  T_result<-kruskal.test(as.numeric(gsva_zong[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_mianyi_cluster<-p_zanshi
#fdr_zanshi<-p.adjust(p_zanshi,method="fdr")
de_mianyi<-as.data.frame(gsva_zong[which(p_zanshi<0.05),])#7
#

class<-rep(c("cluster2","cluster1"),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))

mianyi_huatu<-c()
for (i in 1:nrow(de_mianyi)) {
  b<-as.numeric(de_mianyi[i,])
  tonglu1<-rep(rownames(de_mianyi)[i],ncol(de_mianyi))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)

library(ggplot2)
library(ggthemes)
class(mianyi_huatu)
table(mianyi_huatu$type)
mianyi_huatu$scores<-as.numeric(mianyi_huatu$scores)
mianyi_huatu$class<-as.factor(mianyi_huatu$class)
table(mianyi_huatu$type)

p_zanshi2<-data.frame(p=p_zanshi[which(p_zanshi<0.05)])
rownames(p_zanshi2)<-rownames(de_mianyi)

comlabel<-p_zanshi2$p
comlabel1<-p_zanshi2$p
comlabel[which(comlabel1<0.001)]<-"***"
comlabel[which(comlabel1<0.01&comlabel1>=0.001)]<-"**"
comlabel[which(comlabel1<0.05&comlabel1>=0.01)]<-"*"

comlabel<-rev(comlabel)
# 
comparisons <- combn(unique(mianyi_huatu$class), 2, simplify = FALSE)
aaa<-ggplot(mianyi_huatu,aes(x=scores,y=factor(type),fill=class))+
  geom_boxplot(outlier.size=1)+theme_classic()+scale_y_discrete(position = "right",limits = rev(levels(factor(mianyi_huatu$type))))+xlab("")+ylab("")+guides()+
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+#stat_compare_means(method = "kruskal.test",label = "p.signif",comparisons = comparisons)
  geom_text(aes(y =1, x=-0.6, label = comlabel[1], group = type), size = 4) +
  geom_text(aes(y =2, x=-0.6, label = comlabel[2], group = type), size = 4)+
  geom_text(aes(y =3, x=-0.6, label = comlabel[3], group = type), size = 4)+
  geom_text(aes(y =4, x=-0.6, label = comlabel[4], group = type), size = 4)+
  geom_text(aes(y =5, x=-0.6, label = comlabel[5], group = type), size = 4)+
  geom_text(aes(y =6, x=-0.6, label = comlabel[6], group = type), size = 4)+
  geom_text(aes(y =7, x=-0.6, label = comlabel[7], group = type), size = 4)
# geom_text(aes(y =17, x=0, label = label[17], group = type), size = 4)+
# geom_text(aes(y =18, x=0, label = label[18], group = type), size = 4)+
# geom_text(aes(y =19, x=0, label = label[19], group = type), size = 4)+
# geom_text(aes(y =20, x=0, label = label[20], group = type), size = 4)+
# geom_text(aes(y =21, x=0, label = label[21], group = type), size = 4)+
# geom_text(aes(y =22, x=0, label = label[22], group = type), size = 4)+
# geom_text(aes(y =23, x=0, label = label[23], group = type), size = 4)+
# geom_text(aes(y =24, x=0, label = label[24], group = type), size = 4)+
# geom_text(aes(y =25, x=0, label = label[25], group = type), size = 4)+
# geom_text(aes(y =26, x=0, label = label[26], group = type), size = 4)+
# geom_text(aes(y =27, x=0, label = label[27], group = type), size = 4)+
# geom_text(aes(y =28, x=0, label = label[28], group = type), size = 4)
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\KIRP_ggssea_plot.pdf",aaa,width = 10,height = 10)



library(ggsci)
library(tidyr)
library(ggpubr)

#1.8.2 emt####
#BiocManager::install("KEGGREST")
library(KEGGREST)
#hsa01212:
# KEGG
# path <- keggGet("hsa01212")
# # 
# gene.info <- path[[1]]$GENE
# # gene symbolEntrez ID
# genes <- unlist(lapply(gene.info,function(x) strsplit(x,";")))#
# gene.symbol <- genes[1:length(genes)%%3 == 2]#
#gene.id <- genes[1:length(genes)%%3 == 1]
#
gene_set<- read.xlsx("E:\\pathway2.xlsx")
gene_set<-gene_set[-which(gene_set$pathway=="fatty acid metabolism"),]
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
#Estimates GSVA enrichment scores.

gsva_cluster1<- gsva(as.matrix(cluster1_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_cluster2<- gsva(as.matrix(cluster2_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
#gsva_cluster3<- gsva(as.matrix(cluster3_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

###
gsva_zong<-cbind(gsva_cluster2,gsva_cluster1)
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))

for (j in 1:nrow(gsva_zong)) {
  T_result<-kruskal.test(as.numeric(gsva_zong[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_zanshi
#fdr_zanshi<-p.adjust(p_zanshi,method="fdr")
# de_mianyi<-as.data.frame(gsva_zong[which(p_zanshi<0.05),])#2

de_mianyi<-as.data.frame(gsva_zong)
#

class<-rep(c("cluster2","cluster1"),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))

mianyi_huatu<-c()
for (i in 1:nrow(de_mianyi)) {
  b<-as.numeric(de_mianyi[i,])
  tonglu1<-rep(rownames(de_mianyi)[i],nrow(de_mianyi))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)
# mianyi_huatu$type="Cell migration"

library(ggsignif)

library(ggsci)
library(tidyr)
library(ggpubr)
class(mianyi_huatu)
table(mianyi_huatu$type)
mianyi_huatu$scores<-as.numeric(mianyi_huatu$scores)
mianyi_huatu$class<-as.factor(mianyi_huatu$class)
###
p <- ggviolin(mianyi_huatu, x = "class", # xsurstat
              y = "scores", # y
              fill = "class", #surstat
              facet.by = "type",# variable
              alpha = 1,width = 0.5,legend = "right",legend.title = "Dead",
              ylab="Normalized Expression",  xlab=FALSE,
              font.y = 15,x.text.angle = 45, y.text.angle = 90,
              font.tickslab = c(15,"plain","black"),
              add = "boxplot", 
              add.params = list(fill = "white", width = 0.1,linetype = 1)) 
p <- p + stat_compare_means(method = "wilcox.test", label = "p.format",
                            label.x.npc ="left", size = 5) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))
p




# aaa<-ggplot(mianyi_huatu,aes(y=scores,x=factor(type),fill=class))+geom_violin()+
#   geom_boxplot(width=0.1,position=position_dodge(0.9),outlier.colour = NA)+
#   theme_classic()+xlab("")+ylab("Scores")+guides()+
#   scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+ stat_compare_means(aes(group = class), label = "p.signif")
# aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\KIRP_zhuanyi_plot.pdf",p,width = 10,height = 10)
#LIHC####

#1.8 emt####
#1.8.1####
library(genefilter)
#BiocManager::install("GSVA")
library(GSVA)
library(Biobase)
library(stringr)
library(openxlsx)
gene_set<- read.xlsx('E:\\1-s2.0-S2211124716317090-mmc3.xlsx')##
colnames(gene_set)<-c("Metagene","Cell.type","Immunity")
colnames(gene_set)<-gene_set[1,]
gene_set<-gene_set[-1,]
gene_set<-gene_set[, 1:2]#
head(gene_set)
#
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
#Estimates GSVA enrichment scores.
cluster1_exp<-expr[,rownames(clusters)[which(clusters$clusters=="1")]]
cluster2_exp<-expr[,intersect(rownames(clusters)[which(clusters$clusters=="2")],colnames(expr))]
str(cluster1_exp)
#cluster12
gsva_cluster1<- gsva(as.matrix(cluster1_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_cluster2<- gsva(as.matrix(cluster2_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
# cluster3_exp<-cluster3_exp[,-1]
# gsva_cluster3<- gsva(as.matrix(cluster3_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

###
gsva_zong<-cbind( gsva_cluster2,gsva_cluster1)
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))
for (j in 1:28) {
  T_result<-kruskal.test(as.numeric(gsva_zong[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_mianyi_cluster<-p_zanshi
#fdr_zanshi<-p.adjust(p_zanshi,method="fdr")
de_mianyi<-as.data.frame(gsva_zong[which(p_zanshi<0.05),])#22
#

class<-rep(c("cluster2","cluster1"),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))

mianyi_huatu<-c()
for (i in 1:nrow(de_mianyi)) {
  b<-as.numeric(de_mianyi[i,])
  tonglu1<-rep(rownames(de_mianyi)[i],ncol(de_mianyi))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)

library(ggplot2)
library(ggthemes)
class(mianyi_huatu)
table(mianyi_huatu$type)
mianyi_huatu$scores<-as.numeric(mianyi_huatu$scores)
mianyi_huatu$class<-as.factor(mianyi_huatu$class)
table(mianyi_huatu$type)

p_zanshi2<-data.frame(p=p_zanshi[which(p_zanshi<0.05)])
rownames(p_zanshi2)<-rownames(de_mianyi)

comlabel<-p_zanshi2$p
comlabel1<-p_zanshi2$p
comlabel[which(comlabel1<0.001)]<-"***"
comlabel[which(comlabel1<0.01&comlabel1>=0.001)]<-"**"
comlabel[which(comlabel1<0.05&comlabel1>=0.01)]<-"*"

comlabel<-rev(comlabel)
# 
comparisons <- combn(unique(mianyi_huatu$class), 2, simplify = FALSE)
aaa<-ggplot(mianyi_huatu,aes(x=scores,y=factor(type),fill=class))+
  geom_boxplot(outlier.size=1)+theme_classic()+scale_y_discrete(position = "right",limits = rev(levels(factor(mianyi_huatu$type))))+xlab("")+ylab("")+guides()+
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+#stat_compare_means(method = "kruskal.test",label = "p.signif",comparisons = comparisons)
  geom_text(aes(y =1, x=-0.6, label = comlabel[1], group = type), size = 4) +
  geom_text(aes(y =2, x=-0.6, label = comlabel[2], group = type), size = 4)+
  geom_text(aes(y =3, x=-0.6, label = comlabel[3], group = type), size = 4)+
  geom_text(aes(y =4, x=-0.6, label = comlabel[4], group = type), size = 4)+
  geom_text(aes(y =5, x=-0.6, label = comlabel[5], group = type), size = 4)+
  geom_text(aes(y =6, x=-0.6, label = comlabel[6], group = type), size = 4)+
  geom_text(aes(y =7, x=-0.6, label = comlabel[7], group = type), size = 4)+
  geom_text(aes(y =8, x=-0.6, label = comlabel[8], group = type), size = 4)+
  geom_text(aes(y =9, x=-0.6, label = comlabel[9], group = type), size = 4)+
  geom_text(aes(y =10, x=-0.6, label = comlabel[10], group = type), size = 4)+
  geom_text(aes(y =11, x=-0.6, label = comlabel[11], group = type), size = 4)+
  geom_text(aes(y =12, x=-0.6, label = comlabel[12], group = type), size = 4)+
  geom_text(aes(y =13, x=-0.6, label = comlabel[13], group = type), size = 4)+
  geom_text(aes(y =14, x=-0.6, label = comlabel[14], group = type), size = 4)+
  geom_text(aes(y =15, x=-0.6, label = comlabel[15], group = type), size = 4)+
  geom_text(aes(y =16, x=-0.6, label = comlabel[16], group = type), size = 4)+
  geom_text(aes(y =17, x=-0.6, label = comlabel[17], group = type), size = 4)+
  geom_text(aes(y =18, x=-0.6, label = comlabel[18], group = type), size = 4)+
  geom_text(aes(y =19, x=-0.6, label = comlabel[19], group = type), size = 4)+
  geom_text(aes(y =20, x=-0.6, label = comlabel[20], group = type), size = 4)+
  geom_text(aes(y =21, x=-0.6, label = comlabel[21], group = type), size = 4)+
  geom_text(aes(y =22, x=-0.6, label = comlabel[22], group = type), size = 4)
# geom_text(aes(y =17, x=0, label = label[17], group = type), size = 4)+
# geom_text(aes(y =18, x=0, label = label[18], group = type), size = 4)+
# geom_text(aes(y =19, x=0, label = label[19], group = type), size = 4)+
# geom_text(aes(y =20, x=0, label = label[20], group = type), size = 4)+
# geom_text(aes(y =21, x=0, label = label[21], group = type), size = 4)+
# geom_text(aes(y =22, x=0, label = label[22], group = type), size = 4)+
# geom_text(aes(y =23, x=0, label = label[23], group = type), size = 4)+
# geom_text(aes(y =24, x=0, label = label[24], group = type), size = 4)+
# geom_text(aes(y =25, x=0, label = label[25], group = type), size = 4)+
# geom_text(aes(y =26, x=0, label = label[26], group = type), size = 4)+
# geom_text(aes(y =27, x=0, label = label[27], group = type), size = 4)+
# geom_text(aes(y =28, x=0, label = label[28], group = type), size = 4)
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\LIHC_ggssea_plot.pdf",aaa,width = 10,height = 10)



library(ggsci)
library(tidyr)
library(ggpubr)

#1.8.2 emt####
#BiocManager::install("KEGGREST")
library(KEGGREST)
#hsa01212:
# KEGG
# path <- keggGet("hsa01212")
# # 
# gene.info <- path[[1]]$GENE
# # gene symbolEntrez ID
# genes <- unlist(lapply(gene.info,function(x) strsplit(x,";")))#
# gene.symbol <- genes[1:length(genes)%%3 == 2]#
#gene.id <- genes[1:length(genes)%%3 == 1]
#
gene_set<- read.xlsx("E:\\pathway2.xlsx")
gene_set<-gene_set[-which(gene_set$pathway=="fatty acid metabolism"),]
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
#Estimates GSVA enrichment scores.
cluster1_exp<-cluster1_exp[,-1]
cluster2_exp<-cluster2_exp[,-1]
gsva_cluster1<- gsva(as.matrix(cluster1_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_cluster2<- gsva(as.matrix(cluster2_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
#gsva_cluster3<- gsva(as.matrix(cluster3_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

###
gsva_zong<-cbind(gsva_cluster2,gsva_cluster1)
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))

for (j in 1:nrow(gsva_zong)) {
  T_result<-kruskal.test(as.numeric(gsva_zong[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_zanshi
#fdr_zanshi<-p.adjust(p_zanshi,method="fdr")
de_mianyi<-as.data.frame(gsva_zong[which(p_zanshi<0.05),])#2
#

class<-rep(c("cluster2","cluster1"),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))

mianyi_huatu<-c()
for (i in 1:nrow(de_mianyi)) {
  b<-as.numeric(de_mianyi[i,])
  tonglu1<-rep(rownames(de_mianyi)[i],ncol(de_mianyi))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)
library(ggsignif)

library(ggsci)
library(tidyr)
library(ggpubr)
class(mianyi_huatu)
table(mianyi_huatu$type)
mianyi_huatu$scores<-as.numeric(mianyi_huatu$scores)
mianyi_huatu$class<-as.factor(mianyi_huatu$class)
###
p <- ggviolin(mianyi_huatu, x = "class", # xsurstat
              y = "scores", # y
              fill = "class", #surstat
              facet.by = "type",# variable
              alpha = 1,width = 0.5,legend = "right",legend.title = "Dead",
              ylab="Normalized Expression",  xlab=FALSE,
              font.y = 15,x.text.angle = 45, y.text.angle = 90,
              font.tickslab = c(15,"plain","black"),
              add = "boxplot", 
              add.params = list(fill = "white", width = 0.1,linetype = 1)) 
p <- p + stat_compare_means(method = "wilcox.test", label = "p.format",
                            label.x.npc ="left", size = 5) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))
p




# aaa<-ggplot(mianyi_huatu,aes(y=scores,x=factor(type),fill=class))+geom_violin()+
#   geom_boxplot(width=0.1,position=position_dodge(0.9),outlier.colour = NA)+
#   theme_classic()+xlab("")+ylab("Scores")+guides()+
#   scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+ stat_compare_means(aes(group = class), label = "p.signif")
# aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\LIHC_zhuanyi_plot.pdf",p,width = 10,height = 10)


#LUNG####

#1.8 emt####
#1.8.1####
library(genefilter)
#BiocManager::install("GSVA")
library(GSVA)
library(Biobase)
library(stringr)
library(openxlsx)
gene_set<- read.xlsx('E:\\1-s2.0-S2211124716317090-mmc3.xlsx')##
colnames(gene_set)<-c("Metagene","Cell.type","Immunity")
colnames(gene_set)<-gene_set[1,]
gene_set<-gene_set[-1,]
gene_set<-gene_set[, 1:2]#
head(gene_set)
#
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
#Estimates GSVA enrichment scores.
cluster1_exp<-expr[,rownames(clusters)[which(clusters$clusters=="1")]]
cluster2_exp<-expr[,intersect(rownames(clusters)[which(clusters$clusters=="2")],colnames(expr))]
str(cluster1_exp)
#cluster12
gsva_cluster1<- gsva(as.matrix(cluster1_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_cluster2<- gsva(as.matrix(cluster2_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
# cluster3_exp<-cluster3_exp[,-1]
# gsva_cluster3<- gsva(as.matrix(cluster3_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

###
gsva_zong<-cbind( gsva_cluster2,gsva_cluster1)
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))
for (j in 1:28) {
  T_result<-kruskal.test(as.numeric(gsva_zong[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_mianyi_cluster<-p_zanshi
#fdr_zanshi<-p.adjust(p_zanshi,method="fdr")
de_mianyi<-as.data.frame(gsva_zong[which(p_zanshi<0.05),])#22
#

class<-rep(c("cluster2","cluster1"),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))

mianyi_huatu<-c()
for (i in 1:nrow(de_mianyi)) {
  b<-as.numeric(de_mianyi[i,])
  tonglu1<-rep(rownames(de_mianyi)[i],ncol(de_mianyi))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)

library(ggplot2)
library(ggthemes)
class(mianyi_huatu)
table(mianyi_huatu$type)
mianyi_huatu$scores<-as.numeric(mianyi_huatu$scores)
mianyi_huatu$class<-as.factor(mianyi_huatu$class)
table(mianyi_huatu$type)

p_zanshi2<-data.frame(p=p_zanshi[which(p_zanshi<0.05)])
rownames(p_zanshi2)<-rownames(de_mianyi)

comlabel<-p_zanshi2$p
comlabel1<-p_zanshi2$p
comlabel[which(comlabel1<0.001)]<-"***"
comlabel[which(comlabel1<0.01&comlabel1>=0.001)]<-"**"
comlabel[which(comlabel1<0.05&comlabel1>=0.01)]<-"*"

comlabel<-rev(comlabel)
# 
comparisons <- combn(unique(mianyi_huatu$class), 2, simplify = FALSE)
aaa<-ggplot(mianyi_huatu,aes(x=scores,y=factor(type),fill=class))+
  geom_boxplot(outlier.size=1)+theme_classic()+scale_y_discrete(position = "right",limits = rev(levels(factor(mianyi_huatu$type))))+xlab("")+ylab("")+guides()+
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+#stat_compare_means(method = "kruskal.test",label = "p.signif",comparisons = comparisons)
  geom_text(aes(y =1, x=-0.6, label = comlabel[1], group = type), size = 4) +
  geom_text(aes(y =2, x=-0.6, label = comlabel[2], group = type), size = 4)+
  geom_text(aes(y =3, x=-0.6, label = comlabel[3], group = type), size = 4)+
  geom_text(aes(y =4, x=-0.6, label = comlabel[4], group = type), size = 4)+
  geom_text(aes(y =5, x=-0.6, label = comlabel[5], group = type), size = 4)+
  geom_text(aes(y =6, x=-0.6, label = comlabel[6], group = type), size = 4)+
  geom_text(aes(y =7, x=-0.6, label = comlabel[7], group = type), size = 4)+
  geom_text(aes(y =8, x=-0.6, label = comlabel[8], group = type), size = 4)+
  geom_text(aes(y =9, x=-0.6, label = comlabel[9], group = type), size = 4)+
  geom_text(aes(y =10, x=-0.6, label = comlabel[10], group = type), size = 4)+
  geom_text(aes(y =11, x=-0.6, label = comlabel[11], group = type), size = 4)+
  geom_text(aes(y =12, x=-0.6, label = comlabel[12], group = type), size = 4)+
  geom_text(aes(y =13, x=-0.6, label = comlabel[13], group = type), size = 4)+
  geom_text(aes(y =14, x=-0.6, label = comlabel[14], group = type), size = 4)+
  geom_text(aes(y =15, x=-0.6, label = comlabel[15], group = type), size = 4)+
  geom_text(aes(y =16, x=-0.6, label = comlabel[16], group = type), size = 4)+
  geom_text(aes(y =17, x=-0.6, label = comlabel[17], group = type), size = 4)
# +
#   geom_text(aes(y =18, x=-0.6, label = comlabel[18], group = type), size = 4)+
#   geom_text(aes(y =19, x=-0.6, label = comlabel[19], group = type), size = 4)+
#   geom_text(aes(y =20, x=-0.6, label = comlabel[20], group = type), size = 4)+
#   geom_text(aes(y =21, x=-0.6, label = comlabel[21], group = type), size = 4)+
#   geom_text(aes(y =22, x=-0.6, label = comlabel[22], group = type), size = 4)
# geom_text(aes(y =17, x=0, label = label[17], group = type), size = 4)+
# geom_text(aes(y =18, x=0, label = label[18], group = type), size = 4)+
# geom_text(aes(y =19, x=0, label = label[19], group = type), size = 4)+
# geom_text(aes(y =20, x=0, label = label[20], group = type), size = 4)+
# geom_text(aes(y =21, x=0, label = label[21], group = type), size = 4)+
# geom_text(aes(y =22, x=0, label = label[22], group = type), size = 4)+
# geom_text(aes(y =23, x=0, label = label[23], group = type), size = 4)+
# geom_text(aes(y =24, x=0, label = label[24], group = type), size = 4)+
# geom_text(aes(y =25, x=0, label = label[25], group = type), size = 4)+
# geom_text(aes(y =26, x=0, label = label[26], group = type), size = 4)+
# geom_text(aes(y =27, x=0, label = label[27], group = type), size = 4)+
# geom_text(aes(y =28, x=0, label = label[28], group = type), size = 4)
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\NSCLC_ggssea_plot.pdf",aaa,width = 10,height = 10)



library(ggsci)
library(tidyr)
library(ggpubr)

#1.8.2 emt####
#BiocManager::install("KEGGREST")
library(KEGGREST)
#hsa01212:
# KEGG
# path <- keggGet("hsa01212")
# # 
# gene.info <- path[[1]]$GENE
# # gene symbolEntrez ID
# genes <- unlist(lapply(gene.info,function(x) strsplit(x,";")))#
# gene.symbol <- genes[1:length(genes)%%3 == 2]#
#gene.id <- genes[1:length(genes)%%3 == 1]
#
gene_set<- read.xlsx("E:\\pathway2.xlsx")
gene_set<-gene_set[-which(gene_set$pathway=="fatty acid metabolism"),]
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
#Estimates GSVA enrichment scores.

gsva_cluster1<- gsva(as.matrix(cluster1_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_cluster2<- gsva(as.matrix(cluster2_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
#gsva_cluster3<- gsva(as.matrix(cluster3_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

###
gsva_zong<-cbind(gsva_cluster2,gsva_cluster1)
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))

for (j in 1:nrow(gsva_zong)) {
  T_result<-kruskal.test(as.numeric(gsva_zong[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_zanshi
#fdr_zanshi<-p.adjust(p_zanshi,method="fdr")
de_mianyi<-as.data.frame(gsva_zong[which(p_zanshi<0.05),])#2
#

class<-rep(c("cluster2","cluster1"),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))

mianyi_huatu<-c()
for (i in 1:nrow(de_mianyi)) {
  b<-as.numeric(de_mianyi[i,])
  tonglu1<-rep(rownames(de_mianyi)[i],ncol(de_mianyi))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)
library(ggsignif)

library(ggsci)
library(tidyr)
library(ggpubr)
class(mianyi_huatu)
table(mianyi_huatu$type)
mianyi_huatu$scores<-as.numeric(mianyi_huatu$scores)
mianyi_huatu$class<-as.factor(mianyi_huatu$class)
###
p <- ggviolin(mianyi_huatu, x = "class", # xsurstat
              y = "scores", # y
              fill = "class", #surstat
              facet.by = "type",# variable
              alpha = 1,width = 0.5,legend = "right",legend.title = "Dead",
              ylab="Normalized Expression",  xlab=FALSE,
              font.y = 15,x.text.angle = 45, y.text.angle = 90,
              font.tickslab = c(15,"plain","black"),
              add = "boxplot", 
              add.params = list(fill = "white", width = 0.1,linetype = 1)) 
p <- p + stat_compare_means(method = "wilcox.test", label = "p.format",
                            label.x.npc ="left", size = 5) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))
p




# aaa<-ggplot(mianyi_huatu,aes(y=scores,x=factor(type),fill=class))+geom_violin()+
#   geom_boxplot(width=0.1,position=position_dodge(0.9),outlier.colour = NA)+
#   theme_classic()+xlab("")+ylab("Scores")+guides()+
#   scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+ stat_compare_means(aes(group = class), label = "p.signif")
# aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\NSCLC_zhuanyi_plot.pdf",p,width = 10,height = 10)

#PAAD####

#1.8 emt####
#1.8.1####
library(genefilter)
#BiocManager::install("GSVA")
library(GSVA)
library(Biobase)
library(stringr)
library(openxlsx)
gene_set<- read.xlsx('E:\\1-s2.0-S2211124716317090-mmc3.xlsx')##
colnames(gene_set)<-c("Metagene","Cell.type","Immunity")
colnames(gene_set)<-gene_set[1,]
gene_set<-gene_set[-1,]
gene_set<-gene_set[, 1:2]#
head(gene_set)
#
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
#Estimates GSVA enrichment scores.
cluster1_exp<-expr[,rownames(clusters)[which(clusters$clusters=="1")]]
cluster2_exp<-expr[,intersect(rownames(clusters)[which(clusters$clusters=="2")],colnames(expr))]
str(cluster1_exp)
#cluster12
gsva_cluster1<- gsva(as.matrix(cluster1_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_cluster2<- gsva(as.matrix(cluster2_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
# cluster3_exp<-cluster3_exp[,-1]
# gsva_cluster3<- gsva(as.matrix(cluster3_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

###
gsva_zong<-cbind( gsva_cluster2,gsva_cluster1)
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))
for (j in 1:28) {
  T_result<-kruskal.test(as.numeric(gsva_zong[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_mianyi_cluster<-p_zanshi
#fdr_zanshi<-p.adjust(p_zanshi,method="fdr")
de_mianyi<-as.data.frame(gsva_zong[which(p_zanshi<0.05),])#22
#

class<-rep(c("cluster2","cluster1"),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))

mianyi_huatu<-c()
for (i in 1:nrow(de_mianyi)) {
  b<-as.numeric(de_mianyi[i,])
  tonglu1<-rep(rownames(de_mianyi)[i],ncol(de_mianyi))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)

library(ggplot2)
library(ggthemes)
class(mianyi_huatu)
table(mianyi_huatu$type)
mianyi_huatu$scores<-as.numeric(mianyi_huatu$scores)
mianyi_huatu$class<-as.factor(mianyi_huatu$class)
table(mianyi_huatu$type)

p_zanshi2<-data.frame(p=p_zanshi[which(p_zanshi<0.05)])
rownames(p_zanshi2)<-rownames(de_mianyi)

comlabel<-p_zanshi2$p
comlabel1<-p_zanshi2$p
comlabel[which(comlabel1<0.001)]<-"***"
comlabel[which(comlabel1<0.01&comlabel1>=0.001)]<-"**"
comlabel[which(comlabel1<0.05&comlabel1>=0.01)]<-"*"

comlabel<-rev(comlabel)
# 
comparisons <- combn(unique(mianyi_huatu$class), 2, simplify = FALSE)
aaa<-ggplot(mianyi_huatu,aes(x=scores,y=factor(type),fill=class))+
  geom_boxplot(outlier.size=1)+theme_classic()+scale_y_discrete(position = "right",limits = rev(levels(factor(mianyi_huatu$type))))+xlab("")+ylab("")+guides()+
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+#stat_compare_means(method = "kruskal.test",label = "p.signif",comparisons = comparisons)
  geom_text(aes(y =1, x=-0.6, label = comlabel[1], group = type), size = 4) +
  geom_text(aes(y =2, x=-0.6, label = comlabel[2], group = type), size = 4)+
  geom_text(aes(y =3, x=-0.6, label = comlabel[3], group = type), size = 4)+
  geom_text(aes(y =4, x=-0.6, label = comlabel[4], group = type), size = 4)+
  geom_text(aes(y =5, x=-0.6, label = comlabel[5], group = type), size = 4)+
  geom_text(aes(y =6, x=-0.6, label = comlabel[6], group = type), size = 4)+
  geom_text(aes(y =7, x=-0.6, label = comlabel[7], group = type), size = 4)+
  geom_text(aes(y =8, x=-0.6, label = comlabel[8], group = type), size = 4)+
  geom_text(aes(y =9, x=-0.6, label = comlabel[9], group = type), size = 4)+
  geom_text(aes(y =10, x=-0.6, label = comlabel[10], group = type), size = 4)+
  geom_text(aes(y =11, x=-0.6, label = comlabel[11], group = type), size = 4)+
  geom_text(aes(y =12, x=-0.6, label = comlabel[12], group = type), size = 4)+
  geom_text(aes(y =13, x=-0.6, label = comlabel[13], group = type), size = 4)+
  geom_text(aes(y =14, x=-0.6, label = comlabel[14], group = type), size = 4)+
  geom_text(aes(y =15, x=-0.6, label = comlabel[15], group = type), size = 4)+
  geom_text(aes(y =16, x=-0.6, label = comlabel[16], group = type), size = 4)+
  geom_text(aes(y =17, x=-0.6, label = comlabel[17], group = type), size = 4)+
  geom_text(aes(y =18, x=-0.6, label = comlabel[18], group = type), size = 4)+
  geom_text(aes(y =19, x=-0.6, label = comlabel[19], group = type), size = 4)+
  geom_text(aes(y =20, x=-0.6, label = comlabel[20], group = type), size = 4)
# +
#   geom_text(aes(y =21, x=-0.6, label = comlabel[21], group = type), size = 4)+
#   geom_text(aes(y =22, x=-0.6, label = comlabel[22], group = type), size = 4)
# geom_text(aes(y =17, x=0, label = label[17], group = type), size = 4)+
# geom_text(aes(y =18, x=0, label = label[18], group = type), size = 4)+
# geom_text(aes(y =19, x=0, label = label[19], group = type), size = 4)+
# geom_text(aes(y =20, x=0, label = label[20], group = type), size = 4)+
# geom_text(aes(y =21, x=0, label = label[21], group = type), size = 4)+
# geom_text(aes(y =22, x=0, label = label[22], group = type), size = 4)+
# geom_text(aes(y =23, x=0, label = label[23], group = type), size = 4)+
# geom_text(aes(y =24, x=0, label = label[24], group = type), size = 4)+
# geom_text(aes(y =25, x=0, label = label[25], group = type), size = 4)+
# geom_text(aes(y =26, x=0, label = label[26], group = type), size = 4)+
# geom_text(aes(y =27, x=0, label = label[27], group = type), size = 4)+
# geom_text(aes(y =28, x=0, label = label[28], group = type), size = 4)
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\PAAD_ggssea_plot.pdf",aaa,width = 10,height = 10)



library(ggsci)
library(tidyr)
library(ggpubr)

#1.8.2 emt####
#BiocManager::install("KEGGREST")
library(KEGGREST)
#hsa01212:
# KEGG
# path <- keggGet("hsa01212")
# # 
# gene.info <- path[[1]]$GENE
# # gene symbolEntrez ID
# genes <- unlist(lapply(gene.info,function(x) strsplit(x,";")))#
# gene.symbol <- genes[1:length(genes)%%3 == 2]#
#gene.id <- genes[1:length(genes)%%3 == 1]
#
gene_set<- read.xlsx("E:\\pathway2.xlsx")
gene_set<-gene_set[-which(gene_set$pathway=="fatty acid metabolism"),]
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
#Estimates GSVA enrichment scores.

gsva_cluster1<- gsva(as.matrix(cluster1_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_cluster2<- gsva(as.matrix(cluster2_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
#gsva_cluster3<- gsva(as.matrix(cluster3_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

###
gsva_zong<-cbind(gsva_cluster2,gsva_cluster1)
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))

for (j in 1:nrow(gsva_zong)) {
  T_result<-kruskal.test(as.numeric(gsva_zong[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_zanshi
#fdr_zanshi<-p.adjust(p_zanshi,method="fdr")
de_mianyi<-as.data.frame(gsva_zong)#2
#

class<-rep(c("cluster2","cluster1"),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))

mianyi_huatu<-c()
for (i in 1:nrow(de_mianyi)) {
  b<-as.numeric(de_mianyi[i,])
  tonglu1<-rep(rownames(de_mianyi)[i],ncol(de_mianyi))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)
library(ggsignif)

library(ggsci)
library(tidyr)
library(ggpubr)
class(mianyi_huatu)
table(mianyi_huatu$type)
mianyi_huatu$scores<-as.numeric(mianyi_huatu$scores)
mianyi_huatu$class<-as.factor(mianyi_huatu$class)
###
p <- ggviolin(mianyi_huatu, x = "class", # xsurstat
              y = "scores", # y
              fill = "class", #surstat
              facet.by = "type",# variable
              alpha = 1,width = 0.5,legend = "right",legend.title = "Dead",
              ylab="Normalized Expression",  xlab=FALSE,
              font.y = 15,x.text.angle = 45, y.text.angle = 90,
              font.tickslab = c(15,"plain","black"),
              add = "boxplot", 
              add.params = list(fill = "white", width = 0.1,linetype = 1)) 
p <- p + stat_compare_means(method = "wilcox.test", label = "p.format",
                            label.x.npc ="left", size = 5) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))
p




# aaa<-ggplot(mianyi_huatu,aes(y=scores,x=factor(type),fill=class))+geom_violin()+
#   geom_boxplot(width=0.1,position=position_dodge(0.9),outlier.colour = NA)+
#   theme_classic()+xlab("")+ylab("Scores")+guides()+
#   scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+ stat_compare_means(aes(group = class), label = "p.signif")
# aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\PAAD_zhuanyi_plot.pdf",p,width = 10,height = 10)
#THCA####

#1.8 emt####
#1.8.1####
library(genefilter)
#BiocManager::install("GSVA")
library(GSVA)
library(Biobase)
library(stringr)
library(openxlsx)
gene_set<- read.xlsx('E:\\1-s2.0-S2211124716317090-mmc3.xlsx')##
colnames(gene_set)<-c("Metagene","Cell.type","Immunity")
colnames(gene_set)<-gene_set[1,]
gene_set<-gene_set[-1,]
gene_set<-gene_set[, 1:2]#
head(gene_set)
#
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
#Estimates GSVA enrichment scores.
cluster1_exp<-expr[,rownames(clusters)[which(clusters$clusters=="1")]]
cluster2_exp<-expr[,intersect(rownames(clusters)[which(clusters$clusters=="2")],colnames(expr))]
str(cluster1_exp)
#cluster12
gsva_cluster1<- gsva(as.matrix(cluster1_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_cluster2<- gsva(as.matrix(cluster2_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
# cluster3_exp<-cluster3_exp[,-1]
# gsva_cluster3<- gsva(as.matrix(cluster3_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

###
gsva_zong<-cbind( gsva_cluster2,gsva_cluster1)
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))
for (j in 1:28) {
  T_result<-kruskal.test(as.numeric(gsva_zong[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_mianyi_cluster<-p_zanshi
#fdr_zanshi<-p.adjust(p_zanshi,method="fdr")
de_mianyi<-as.data.frame(gsva_zong[which(p_zanshi<0.05),])#22
#

class<-rep(c("cluster2","cluster1"),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))

mianyi_huatu<-c()
for (i in 1:nrow(de_mianyi)) {
  b<-as.numeric(de_mianyi[i,])
  tonglu1<-rep(rownames(de_mianyi)[i],ncol(de_mianyi))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)

library(ggplot2)
library(ggthemes)
class(mianyi_huatu)
table(mianyi_huatu$type)
mianyi_huatu$scores<-as.numeric(mianyi_huatu$scores)
mianyi_huatu$class<-as.factor(mianyi_huatu$class)
table(mianyi_huatu$type)

p_zanshi2<-data.frame(p=p_zanshi[which(p_zanshi<0.05)])
rownames(p_zanshi2)<-rownames(de_mianyi)

comlabel<-p_zanshi2$p
comlabel1<-p_zanshi2$p
comlabel[which(comlabel1<0.001)]<-"***"
comlabel[which(comlabel1<0.01&comlabel1>=0.001)]<-"**"
comlabel[which(comlabel1<0.05&comlabel1>=0.01)]<-"*"

comlabel<-rev(comlabel)
# 
comparisons <- combn(unique(mianyi_huatu$class), 2, simplify = FALSE)
aaa<-ggplot(mianyi_huatu,aes(x=scores,y=factor(type),fill=class))+
  geom_boxplot(outlier.size=1)+theme_classic()+scale_y_discrete(position = "right",limits = rev(levels(factor(mianyi_huatu$type))))+xlab("")+ylab("")+guides()+
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+#stat_compare_means(method = "kruskal.test",label = "p.signif",comparisons = comparisons)
  geom_text(aes(y =1, x=-0.6, label = comlabel[1], group = type), size = 4) +
  geom_text(aes(y =2, x=-0.6, label = comlabel[2], group = type), size = 4)+
  geom_text(aes(y =3, x=-0.6, label = comlabel[3], group = type), size = 4)+
  geom_text(aes(y =4, x=-0.6, label = comlabel[4], group = type), size = 4)+
  geom_text(aes(y =5, x=-0.6, label = comlabel[5], group = type), size = 4)+
  geom_text(aes(y =6, x=-0.6, label = comlabel[6], group = type), size = 4)+
  geom_text(aes(y =7, x=-0.6, label = comlabel[7], group = type), size = 4)+
  geom_text(aes(y =8, x=-0.6, label = comlabel[8], group = type), size = 4)+
  geom_text(aes(y =9, x=-0.6, label = comlabel[9], group = type), size = 4)+
  geom_text(aes(y =10, x=-0.6, label = comlabel[10], group = type), size = 4)+
  geom_text(aes(y =11, x=-0.6, label = comlabel[11], group = type), size = 4)+
  geom_text(aes(y =12, x=-0.6, label = comlabel[12], group = type), size = 4)+
  geom_text(aes(y =13, x=-0.6, label = comlabel[13], group = type), size = 4)+
  geom_text(aes(y =14, x=-0.6, label = comlabel[14], group = type), size = 4)+
  geom_text(aes(y =15, x=-0.6, label = comlabel[15], group = type), size = 4)+
  geom_text(aes(y =16, x=-0.6, label = comlabel[16], group = type), size = 4)+
  geom_text(aes(y =17, x=-0.6, label = comlabel[17], group = type), size = 4)+
  geom_text(aes(y =18, x=-0.6, label = comlabel[18], group = type), size = 4)+
  geom_text(aes(y =19, x=-0.6, label = comlabel[19], group = type), size = 4)+
  geom_text(aes(y =20, x=-0.6, label = comlabel[20], group = type), size = 4)+
  geom_text(aes(y =21, x=-0.6, label = comlabel[21], group = type), size = 4)+
  geom_text(aes(y =22, x=-0.6, label = comlabel[22], group = type), size = 4)
# geom_text(aes(y =17, x=0, label = label[17], group = type), size = 4)+
# geom_text(aes(y =18, x=0, label = label[18], group = type), size = 4)+
# geom_text(aes(y =19, x=0, label = label[19], group = type), size = 4)+
# geom_text(aes(y =20, x=0, label = label[20], group = type), size = 4)+
# geom_text(aes(y =21, x=0, label = label[21], group = type), size = 4)+
# geom_text(aes(y =22, x=0, label = label[22], group = type), size = 4)+
# geom_text(aes(y =23, x=0, label = label[23], group = type), size = 4)+
# geom_text(aes(y =24, x=0, label = label[24], group = type), size = 4)+
# geom_text(aes(y =25, x=0, label = label[25], group = type), size = 4)+
# geom_text(aes(y =26, x=0, label = label[26], group = type), size = 4)+
# geom_text(aes(y =27, x=0, label = label[27], group = type), size = 4)+
# geom_text(aes(y =28, x=0, label = label[28], group = type), size = 4)
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\THCA_ggssea_plot.pdf",aaa,width = 10,height = 10)



library(ggsci)
library(tidyr)
library(ggpubr)

#1.8.2 emt####
#BiocManager::install("KEGGREST")
library(KEGGREST)
#hsa01212:
# KEGG
# path <- keggGet("hsa01212")
# # 
# gene.info <- path[[1]]$GENE
# # gene symbolEntrez ID
# genes <- unlist(lapply(gene.info,function(x) strsplit(x,";")))#
# gene.symbol <- genes[1:length(genes)%%3 == 2]#
#gene.id <- genes[1:length(genes)%%3 == 1]
#
gene_set<- read.xlsx("E:\\pathway2.xlsx")
gene_set<-gene_set[-which(gene_set$pathway=="fatty acid metabolism"),]
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
#Estimates GSVA enrichment scores.

gsva_cluster1<- gsva(as.matrix(cluster1_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_cluster2<- gsva(as.matrix(cluster2_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
#gsva_cluster3<- gsva(as.matrix(cluster3_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

###
gsva_zong<-cbind(gsva_cluster2,gsva_cluster1)
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))

for (j in 1:nrow(gsva_zong)) {
  T_result<-kruskal.test(as.numeric(gsva_zong[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_zanshi
#fdr_zanshi<-p.adjust(p_zanshi,method="fdr")
de_mianyi<-as.data.frame(gsva_zong[which(p_zanshi<0.05),])#2
#

class<-rep(c("cluster2","cluster1"),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))

mianyi_huatu<-c()
for (i in 1:nrow(de_mianyi)) {
  b<-as.numeric(de_mianyi[i,])
  tonglu1<-rep(rownames(de_mianyi)[i],ncol(de_mianyi))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)
library(ggsignif)

library(ggsci)
library(tidyr)
library(ggpubr)
class(mianyi_huatu)
table(mianyi_huatu$type)
mianyi_huatu$scores<-as.numeric(mianyi_huatu$scores)
mianyi_huatu$class<-as.factor(mianyi_huatu$class)
###
p <- ggviolin(mianyi_huatu, x = "class", # xsurstat
              y = "scores", # y
              fill = "class", #surstat
              facet.by = "type",# variable
              alpha = 1,width = 0.5,legend = "right",legend.title = "Dead",
              ylab="Normalized Expression",  xlab=FALSE,
              font.y = 15,x.text.angle = 45, y.text.angle = 90,
              font.tickslab = c(15,"plain","black"),
              add = "boxplot", 
              add.params = list(fill = "white", width = 0.1,linetype = 1)) 
p <- p + stat_compare_means(method = "wilcox.test", label = "p.format",
                            label.x.npc ="left", size = 5) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))
p




# aaa<-ggplot(mianyi_huatu,aes(y=scores,x=factor(type),fill=class))+geom_violin()+
#   geom_boxplot(width=0.1,position=position_dodge(0.9),outlier.colour = NA)+
#   theme_classic()+xlab("")+ylab("Scores")+guides()+
#   scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+ stat_compare_means(aes(group = class), label = "p.signif")
# aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\THCA_zhuanyi_plot.pdf",p,width = 10,height = 10)
#UCEC####

#1.8 emt####
#1.8.1####
library(genefilter)
#BiocManager::install("GSVA")
library(GSVA)
library(Biobase)
library(stringr)
library(openxlsx)
gene_set<- read.xlsx('E:\\1-s2.0-S2211124716317090-mmc3.xlsx')##
colnames(gene_set)<-c("Metagene","Cell.type","Immunity")
colnames(gene_set)<-gene_set[1,]
gene_set<-gene_set[-1,]
gene_set<-gene_set[, 1:2]#
head(gene_set)
#
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
#Estimates GSVA enrichment scores.
cluster1_exp<-expr[,rownames(clusters)[which(clusters$clusters=="1")]]
cluster2_exp<-expr[,intersect(rownames(clusters)[which(clusters$clusters=="2")],colnames(expr))]
str(cluster1_exp)
#cluster12
gsva_cluster1<- gsva(as.matrix(cluster1_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_cluster2<- gsva(as.matrix(cluster2_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
# cluster3_exp<-cluster3_exp[,-1]
# gsva_cluster3<- gsva(as.matrix(cluster3_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

###
gsva_zong<-cbind( gsva_cluster2,gsva_cluster1)
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))
for (j in 1:28) {
  T_result<-kruskal.test(as.numeric(gsva_zong[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_mianyi_cluster<-p_zanshi
#fdr_zanshi<-p.adjust(p_zanshi,method="fdr")
de_mianyi<-as.data.frame(gsva_zong[which(p_zanshi<0.05),])#16
#

class<-rep(c("cluster2","cluster1"),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))

mianyi_huatu<-c()
for (i in 1:nrow(de_mianyi)) {
  b<-as.numeric(de_mianyi[i,])
  tonglu1<-rep(rownames(de_mianyi)[i],ncol(de_mianyi))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)

library(ggplot2)
library(ggthemes)
class(mianyi_huatu)
table(mianyi_huatu$type)
mianyi_huatu$scores<-as.numeric(mianyi_huatu$scores)
mianyi_huatu$class<-as.factor(mianyi_huatu$class)
table(mianyi_huatu$type)

p_zanshi2<-data.frame(p=p_zanshi[which(p_zanshi<0.05)])
rownames(p_zanshi2)<-rownames(de_mianyi)

comlabel<-p_zanshi2$p
comlabel1<-p_zanshi2$p
comlabel[which(comlabel1<0.001)]<-"***"
comlabel[which(comlabel1<0.01&comlabel1>=0.001)]<-"**"
comlabel[which(comlabel1<0.05&comlabel1>=0.01)]<-"*"

comlabel<-rev(comlabel)
# 
comparisons <- combn(unique(mianyi_huatu$class), 2, simplify = FALSE)
aaa<-ggplot(mianyi_huatu,aes(x=scores,y=factor(type),fill=class))+
  geom_boxplot(outlier.size=1)+theme_classic()+scale_y_discrete(position = "right",limits = rev(levels(factor(mianyi_huatu$type))))+xlab("")+ylab("")+guides()+
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+#stat_compare_means(method = "kruskal.test",label = "p.signif",comparisons = comparisons)
  geom_text(aes(y =1, x=-0.6, label = comlabel[1], group = type), size = 4) +
  geom_text(aes(y =2, x=-0.6, label = comlabel[2], group = type), size = 4)+
  geom_text(aes(y =3, x=-0.6, label = comlabel[3], group = type), size = 4)+
  geom_text(aes(y =4, x=-0.6, label = comlabel[4], group = type), size = 4)+
  geom_text(aes(y =5, x=-0.6, label = comlabel[5], group = type), size = 4)+
  geom_text(aes(y =6, x=-0.6, label = comlabel[6], group = type), size = 4)+
  geom_text(aes(y =7, x=-0.6, label = comlabel[7], group = type), size = 4)+
  geom_text(aes(y =8, x=-0.6, label = comlabel[8], group = type), size = 4)+
  geom_text(aes(y =9, x=-0.6, label = comlabel[9], group = type), size = 4)+
  geom_text(aes(y =10, x=-0.6, label = comlabel[10], group = type), size = 4)+
  geom_text(aes(y =11, x=-0.6, label = comlabel[11], group = type), size = 4)+
  geom_text(aes(y =12, x=-0.6, label = comlabel[12], group = type), size = 4)+
  geom_text(aes(y =13, x=-0.6, label = comlabel[13], group = type), size = 4)+
  geom_text(aes(y =14, x=-0.6, label = comlabel[14], group = type), size = 4)+
  geom_text(aes(y =15, x=-0.6, label = comlabel[15], group = type), size = 4)+
  geom_text(aes(y =16, x=-0.6, label = comlabel[16], group = type), size = 4)
# +
#   geom_text(aes(y =17, x=-0.6, label = comlabel[17], group = type), size = 4)+
#   geom_text(aes(y =18, x=-0.6, label = comlabel[18], group = type), size = 4)+
#   geom_text(aes(y =19, x=-0.6, label = comlabel[19], group = type), size = 4)+
#   geom_text(aes(y =20, x=-0.6, label = comlabel[20], group = type), size = 4)+
#   geom_text(aes(y =21, x=-0.6, label = comlabel[21], group = type), size = 4)+
#   geom_text(aes(y =22, x=-0.6, label = comlabel[22], group = type), size = 4)
# geom_text(aes(y =17, x=0, label = label[17], group = type), size = 4)+
# geom_text(aes(y =18, x=0, label = label[18], group = type), size = 4)+
# geom_text(aes(y =19, x=0, label = label[19], group = type), size = 4)+
# geom_text(aes(y =20, x=0, label = label[20], group = type), size = 4)+
# geom_text(aes(y =21, x=0, label = label[21], group = type), size = 4)+
# geom_text(aes(y =22, x=0, label = label[22], group = type), size = 4)+
# geom_text(aes(y =23, x=0, label = label[23], group = type), size = 4)+
# geom_text(aes(y =24, x=0, label = label[24], group = type), size = 4)+
# geom_text(aes(y =25, x=0, label = label[25], group = type), size = 4)+
# geom_text(aes(y =26, x=0, label = label[26], group = type), size = 4)+
# geom_text(aes(y =27, x=0, label = label[27], group = type), size = 4)+
# geom_text(aes(y =28, x=0, label = label[28], group = type), size = 4)
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\UCEC_ggssea_plot.pdf",aaa,width = 10,height = 10)



library(ggsci)
library(tidyr)
library(ggpubr)

#1.8.2 emt####
#BiocManager::install("KEGGREST")
library(KEGGREST)
#hsa01212:
# KEGG
# path <- keggGet("hsa01212")
# # 
# gene.info <- path[[1]]$GENE
# # gene symbolEntrez ID
# genes <- unlist(lapply(gene.info,function(x) strsplit(x,";")))#
# gene.symbol <- genes[1:length(genes)%%3 == 2]#
#gene.id <- genes[1:length(genes)%%3 == 1]
#
gene_set<- read.xlsx("E:\\pathway2.xlsx")
gene_set<-gene_set[-which(gene_set$pathway=="fatty acid metabolism"),]
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
#Estimates GSVA enrichment scores.

gsva_cluster1<- gsva(as.matrix(cluster1_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_cluster2<- gsva(as.matrix(cluster2_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
#gsva_cluster3<- gsva(as.matrix(cluster3_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

###
gsva_zong<-cbind(gsva_cluster2,gsva_cluster1)
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))

for (j in 1:nrow(gsva_zong)) {
  T_result<-kruskal.test(as.numeric(gsva_zong[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_zanshi
#fdr_zanshi<-p.adjust(p_zanshi,method="fdr")
de_mianyi<-as.data.frame(gsva_zong[which(p_zanshi<0.05),])#2
#

class<-rep(c("cluster2","cluster1"),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))

mianyi_huatu<-c()
for (i in 1:nrow(de_mianyi)) {
  b<-as.numeric(de_mianyi[i,])
  tonglu1<-rep(rownames(de_mianyi)[i],ncol(de_mianyi))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)
library(ggsignif)

library(ggsci)
library(tidyr)
library(ggpubr)
class(mianyi_huatu)
table(mianyi_huatu$type)
mianyi_huatu$scores<-as.numeric(mianyi_huatu$scores)
mianyi_huatu$class<-as.factor(mianyi_huatu$class)
###
p <- ggviolin(mianyi_huatu, x = "class", # xsurstat
              y = "scores", # y
              fill = "class", #surstat
              facet.by = "type",# variable
              alpha = 1,width = 0.5,legend = "right",legend.title = "Dead",
              ylab="Normalized Expression",  xlab=FALSE,
              font.y = 15,x.text.angle = 45, y.text.angle = 90,
              font.tickslab = c(15,"plain","black"),
              add = "boxplot", 
              add.params = list(fill = "white", width = 0.1,linetype = 1)) 
p <- p + stat_compare_means(method = "wilcox.test", label = "p.format",
                            label.x.npc ="left", size = 5) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))
p




# aaa<-ggplot(mianyi_huatu,aes(y=scores,x=factor(type),fill=class))+geom_violin()+
#   geom_boxplot(width=0.1,position=position_dodge(0.9),outlier.colour = NA)+
#   theme_classic()+xlab("")+ylab("Scores")+guides()+
#   scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+ stat_compare_means(aes(group = class), label = "p.signif")
# aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\UCEC_zhuanyi_plot.pdf",p,width = 10,height = 10)

