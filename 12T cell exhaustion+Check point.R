#BRCA####
#T cell exhaustion
# devtools::install_github("GuoBioinfoLab/TCellSI")
library(TCellSI)
ResultScores <- TCellSI::TCSS_Calculate(expr) 
#cluster12
gsva_cluster1<- ResultScores[,cluster1]
gsva_cluster2<- ResultScores[,intersect(cluster2,colnames(gsva_cluster2))]
# intersect(cluster2,colnames(gsva_cluster2))
# cluster3_exp<-cluster3_exp[,-1]
# gsva_cluster3<- gsva(as.matrix(cluster3_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

###
gsva_zong<-cbind( gsva_cluster2,gsva_cluster1)
gsva_zong<-gsva_zong[c("Progenitor_exhaustion","Terminal_exhaustion"),]
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))
for (j in 1:2) {
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
for (i in 1:nrow(gsva_zong)) {
  b<-as.numeric(gsva_zong[i,])
  tonglu1<-rep(rownames(gsva_zong)[i],ncol(gsva_zong))
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
mianyi_huatu <- mianyi_huatu %>%
  filter(mianyi_huatu$scores < quantile(mianyi_huatu$scores, 0.95)) # 95%
aaa<-ggplot(mianyi_huatu, aes(x = factor(type), y = scores, fill = class)) +
  geom_violin(width = 0.8,position = position_dodge(0.8), trim = TRUE, alpha = 0.6) +  # alpha 
  geom_boxplot(width = 0.2, position = position_dodge(0.8), color = "black", 
               outlier.shape = NA, alpha = 0.5) +  # 
  scale_fill_manual(values=c("#80D1C8","#FFD4A9")) +
  labs(x = "", y = 'Score') +
  stat_compare_means(method = "kruskal.test", label = "p.signif") +  # Kruskal-Wallis 
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  #  X 
    axis.text.y = element_text(size = 12),
    legend.title = element_blank(),  # 
    legend.text = element_text(size = 12)
  )

aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\BRCA_T_exhau_plot.pdf",aaa,width = 10,height = 10)

#CHECK POINT

checkpoint<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")
check_exp<-expr[checkpoint$Checkpoint.gene[which(checkpoint$Activator.A..Inhibitor.I.=="I")],]
check_exp<-expr[checkpoint$Checkpoint.gene,]

check_exp<-na.omit(check_exp)
cluster1_check<-check_exp[,cluster1]
cluster2_check<-check_exp[,intersect(cluster2,colnames(check_exp))]
cluster_check<-cbind(cluster1_check,cluster2_check)
p_zanshi<-c()
mianyi_class<-rep(c(1,2),c(ncol(cluster1_check),ncol(cluster2_check)))
for (j in 1:nrow(cluster_check)) {
  T_result<-kruskal.test(as.numeric(cluster_check[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
which(p_zanshi<0.05)

checkgene_sig<-rownames(cluster_check)[which(p_zanshi<0.05)]
checkgene_sig_p<-cbind(checkgene_sig,p_zanshi[which(p_zanshi<0.05)])
write.table(checkgene_sig_p,"G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\BRCA_checkgene_sig_p.txt")

checkgene_sig2<-checkpoint[checkpoint$Checkpoint.gene %in%checkgene_sig,]
checkgene_sig2<-checkgene_sig2[order(checkgene_sig2$Activator.A..Inhibitor.I.), ]
cluster1_check<-check_exp[checkgene_sig2$Checkpoint.gene,cluster1]
cluster2_check<-check_exp[checkgene_sig2$Checkpoint.gene,intersect(cluster2,colnames(check_exp))]

cluster1_check<-apply(cluster1_check, 1, mean)
cluster2_check<-apply(cluster2_check, 1, mean)
library(fmsb)
leida<-rbind(cluster1_check,cluster2_check)
rownames(leida)<-c("cluster1","cluster2")
library(ggradar)
classleida<-c("cluster1","cluster2")
leida<-cbind(classleida,leida)
leida<-as.data.frame(leida)
str(leida)
leida[,2:ncol(leida)]<-apply(leida[,2:ncol(leida)],1,as.numeric)
max(leida[1,2:ncol(leida)])
min(leida[1,2:ncol(leida)])
max(leida[2,2:ncol(leida)])
min(leida[2,2:ncol(leida)])

plot <-ggradar(plot.data =leida,
               grid.min = 0,
               grid.mid = 45,
               grid.max = 90, 
               values.radar = c(0,45,90),group.colours = c("#80D1C8","#FFD4A9"),axis.label.size=4)
plot
colnames(leida)

ggsave("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\BRCA_CHECK_plot.pdf",plot,width = 10,height = 10)
write.table(checkgene_sig2,"G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\BRCA_checkgene_sig2.txt")
#COAD####

#T cell exhaustion
# devtools::install_github("GuoBioinfoLab/TCellSI")
library(TCellSI)
ResultScores <- TCellSI::TCSS_Calculate(expr) 
#cluster12
gsva_cluster1<- ResultScores[,cluster1]
gsva_cluster2<- ResultScores[,cluster2]
# cluster3_exp<-cluster3_exp[,-1]
# gsva_cluster3<- gsva(as.matrix(cluster3_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

###
gsva_zong<-cbind( gsva_cluster2,gsva_cluster1)
gsva_zong<-gsva_zong[c("Progenitor_exhaustion","Terminal_exhaustion"),]
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))
for (j in 1:2) {
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
for (i in 1:nrow(gsva_zong)) {
  b<-as.numeric(gsva_zong[i,])
  tonglu1<-rep(rownames(gsva_zong)[i],ncol(gsva_zong))
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

aaa<-ggplot(mianyi_huatu, aes(x = factor(type), y = scores, fill = class)) +
  geom_violin(position = position_dodge(0.8), trim = TRUE, alpha = 0.6) +  # alpha 
  geom_boxplot(width = 0.2, position = position_dodge(0.8), color = "black", 
               outlier.shape = NA, alpha = 0.5) +  # 
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8")) +
  labs(x = "", y = 'Score') +
  stat_compare_means(method = "kruskal.test", label = "p.signif") +  # Kruskal-Wallis 
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  #  X 
    axis.text.y = element_text(size = 12),
    legend.title = element_blank(),  # 
    legend.text = element_text(size = 12)
  )

aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\COAD_T_exhau_plot.pdf",aaa,width = 10,height = 10)

#CHECK POINT

checkpoint<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")
check_exp<-expr[checkpoint$Checkpoint.gene[which(checkpoint$Activator.A..Inhibitor.I.=="I")],]
check_exp<-expr[checkpoint$Checkpoint.gene,]

check_exp<-na.omit(check_exp)
cluster1_check<-check_exp[,cluster1]
cluster2_check<-check_exp[,cluster2]
cluster_check<-cbind(cluster1_check,cluster2_check)
p_zanshi<-c()
mianyi_class<-rep(c(1,2),c(ncol(cluster1_check),ncol(cluster2_check)))
for (j in 1:nrow(cluster_check)) {
  T_result<-kruskal.test(as.numeric(cluster_check[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
which(p_zanshi<0.05)

checkgene_sig<-rownames(cluster_check)[which(p_zanshi<0.05)]
checkgene_sig_p<-cbind(checkgene_sig,p_zanshi[which(p_zanshi<0.05)])
write.table(checkgene_sig_p,"G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\coad_checkgene_sig_p.txt")

checkgene_sig2<-checkpoint[checkpoint$Checkpoint.gene %in%checkgene_sig,]
checkgene_sig2<-checkgene_sig2[order(checkgene_sig2$Activator.A..Inhibitor.I.), ]
cluster1_check<-check_exp[checkgene_sig2$Checkpoint.gene,cluster1]
cluster2_check<-check_exp[checkgene_sig2$Checkpoint.gene,cluster2]

cluster1_check<-apply(cluster1_check, 1, mean)
cluster2_check<-apply(cluster2_check, 1, mean)
library(fmsb)
leida<-rbind(cluster1_check,cluster2_check)
rownames(leida)<-c("cluster1","cluster2")
library(ggradar)
classleida<-c("cluster1","cluster2")
leida<-cbind(classleida,leida)
leida<-as.data.frame(leida)
str(leida)
leida[,2:ncol(leida)]<-apply(leida[,2:ncol(leida)],1,as.numeric)
max(leida[1,2:ncol(leida)])
min(leida[1,2:ncol(leida)])
max(leida[2,2:ncol(leida)])
min(leida[2,2:ncol(leida)])

plot <-ggradar(plot.data =leida,
               grid.min = 0,
               grid.mid = 0.21,
               grid.max = 0.42, 
               values.radar = c(0,0.21,0.42),group.colours = c("#FFD4A9", "#80D1C8"),axis.label.size=4)
plot
colnames(leida)

ggsave("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\COAD_CHECK_plot.pdf",plot,width = 10,height = 10)
write.table(checkgene_sig2,"G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\COAD_checkgene_sig2.txt")


#ESCA####\

#T cell exhaustion
# devtools::install_github("GuoBioinfoLab/TCellSI")
library(TCellSI)
ResultScores <- TCellSI::TCSS_Calculate(expr) 
#cluster12
gsva_cluster1<- ResultScores[,cluster1]
gsva_cluster2<- ResultScores[,cluster2]
# cluster3_exp<-cluster3_exp[,-1]
# gsva_cluster3<- gsva(as.matrix(cluster3_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

###
gsva_zong<-cbind( gsva_cluster2,gsva_cluster1)
gsva_zong<-gsva_zong[c("Progenitor_exhaustion","Terminal_exhaustion"),]
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))
for (j in 1:2) {
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
for (i in 1:nrow(gsva_zong)) {
  b<-as.numeric(gsva_zong[i,])
  tonglu1<-rep(rownames(gsva_zong)[i],ncol(gsva_zong))
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
mianyi_huatu <- mianyi_huatu %>%
  filter(mianyi_huatu$scores < quantile(mianyi_huatu$scores, 0.95)) # 95%


aaa<-ggplot(mianyi_huatu, aes(x = factor(type), y = scores, fill = class)) +
  geom_violin(position = position_dodge(0.8), trim = TRUE, alpha = 0.6) +  # alpha 
  geom_boxplot(width = 0.2, position = position_dodge(0.8), color = "black", 
               outlier.shape = NA, alpha = 0.5) +  # 
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8")) +
  labs(x = "", y = 'Score') +
  stat_compare_means(method = "kruskal.test", label = "p.signif") +  # Kruskal-Wallis 
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  #  X 
    axis.text.y = element_text(size = 12),
    legend.title = element_blank(),  # 
    legend.text = element_text(size = 12)
  )

aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\ESCA_T_exhau_plot.pdf",aaa,width = 10,height = 10)

#CHECK POINT

checkpoint<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")
check_exp<-expr[checkpoint$Checkpoint.gene[which(checkpoint$Activator.A..Inhibitor.I.=="I")],]
check_exp<-expr[checkpoint$Checkpoint.gene,]

check_exp<-na.omit(check_exp)
cluster1_check<-check_exp[,cluster1]
cluster2_check<-check_exp[,cluster2]
cluster_check<-cbind(cluster1_check,cluster2_check)
p_zanshi<-c()
mianyi_class<-rep(c(1,2),c(ncol(cluster1_check),ncol(cluster2_check)))
for (j in 1:nrow(cluster_check)) {
  T_result<-kruskal.test(as.numeric(cluster_check[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
which(p_zanshi<0.05)

checkgene_sig<-rownames(cluster_check)[which(p_zanshi<0.05)]
checkgene_sig_p<-cbind(checkgene_sig,p_zanshi[which(p_zanshi<0.05)])
write.table(checkgene_sig_p,"ESCA_checkgene_sig_p.txt")

checkgene_sig2<-checkpoint[checkpoint$Checkpoint.gene %in%checkgene_sig,]
checkgene_sig2<-checkgene_sig2[order(checkgene_sig2$Activator.A..Inhibitor.I.), ]
cluster1_check<-check_exp[checkgene_sig2$Checkpoint.gene,cluster1]
cluster2_check<-check_exp[checkgene_sig2$Checkpoint.gene,cluster2]

cluster1_check<-apply(cluster1_check, 1, mean)
cluster2_check<-apply(cluster2_check, 1, mean)
library(fmsb)
leida<-rbind(cluster1_check,cluster2_check)
rownames(leida)<-c("cluster1","cluster2")
library(ggradar)
classleida<-c("cluster1","cluster2")
leida<-cbind(classleida,leida)
leida<-as.data.frame(leida)
str(leida)
leida[,2:ncol(leida)]<-apply(leida[,2:ncol(leida)],1,as.numeric)
max(leida[1,2:ncol(leida)])
min(leida[1,2:ncol(leida)])
max(leida[2,2:ncol(leida)])
min(leida[2,2:ncol(leida)])

plot <-ggradar(plot.data =leida,
               grid.min = 0,
               grid.mid = 0.21,
               grid.max = 0.42, 
               values.radar = c(0,0.21,0.42),group.colours = c("#FFD4A9", "#80D1C8"),axis.label.size=4)

colnames(leida)
write.table(checkgene_sig2,"ESCA_checkgene_sig2.txt")
ggsave("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\ESCA_CHECK_plot.pdf",plot,width = 10,height = 10)
#HNSC####

#T cell exhaustion
# devtools::install_github("GuoBioinfoLab/TCellSI")
library(TCellSI)
ResultScores <- TCellSI::TCSS_Calculate(expr) 
#cluster12
gsva_cluster1<- ResultScores[,cluster1]
gsva_cluster2<- ResultScores[,cluster2]
# cluster3_exp<-cluster3_exp[,-1]
# gsva_cluster3<- gsva(as.matrix(cluster3_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

###
gsva_zong<-cbind( gsva_cluster2,gsva_cluster1)
gsva_zong<-gsva_zong[c("Progenitor_exhaustion","Terminal_exhaustion"),]
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))
for (j in 1:2) {
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
for (i in 1:nrow(gsva_zong)) {
  b<-as.numeric(gsva_zong[i,])
  tonglu1<-rep(rownames(gsva_zong)[i],ncol(gsva_zong))
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

aaa<-ggplot(mianyi_huatu, aes(x = factor(type), y = scores, fill = class)) +
  geom_violin(position = position_dodge(0.8), trim = TRUE, alpha = 0.6) +  # alpha 
  geom_boxplot(width = 0.2, position = position_dodge(0.8), color = "black", 
               outlier.shape = NA, alpha = 0.5) +  # 
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8")) +
  labs(x = "", y = 'Score') +
  stat_compare_means(method = "kruskal.test", label = "p.signif") +  # Kruskal-Wallis 
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  #  X 
    axis.text.y = element_text(size = 12),
    legend.title = element_blank(),  # 
    legend.text = element_text(size = 12)
  )

aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_T_exhau_plot.pdf",aaa,width = 10,height = 10)

#CHECK POINT

checkpoint<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")
check_exp<-expr[checkpoint$Checkpoint.gene[which(checkpoint$Activator.A..Inhibitor.I.=="I")],]
check_exp<-expr[checkpoint$Checkpoint.gene,]

check_exp<-na.omit(check_exp)
cluster1_check<-check_exp[,cluster1]
cluster2_check<-check_exp[,cluster2]
cluster_check<-cbind(cluster1_check,cluster2_check)
p_zanshi<-c()
mianyi_class<-rep(c(1,2),c(ncol(cluster1_check),ncol(cluster2_check)))
for (j in 1:nrow(cluster_check)) {
  T_result<-kruskal.test(as.numeric(cluster_check[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
which(p_zanshi<0.05)

checkgene_sig<-rownames(cluster_check)[which(p_zanshi<0.05)]
checkgene_sig_p<-cbind(checkgene_sig,p_zanshi[which(p_zanshi<0.05)])
write.table(checkgene_sig_p,"HNSC_checkgene_sig_p.txt")

checkgene_sig2<-checkpoint[checkpoint$Checkpoint.gene %in%checkgene_sig,]
checkgene_sig2<-checkgene_sig2[order(checkgene_sig2$Activator.A..Inhibitor.I.), ]
cluster1_check<-check_exp[checkgene_sig2$Checkpoint.gene,cluster1]
cluster2_check<-check_exp[checkgene_sig2$Checkpoint.gene,cluster2]

cluster1_check<-apply(cluster1_check, 1, mean)
cluster2_check<-apply(cluster2_check, 1, mean)
library(fmsb)
leida<-rbind(cluster1_check,cluster2_check)
rownames(leida)<-c("cluster1","cluster2")
library(ggradar)
classleida<-c("cluster1","cluster2")
leida<-cbind(classleida,leida)
leida<-as.data.frame(leida)
str(leida)
leida[,2:ncol(leida)]<-apply(leida[,2:ncol(leida)],1,as.numeric)
max(leida[1,2:ncol(leida)])
min(leida[1,2:ncol(leida)])
max(leida[2,2:ncol(leida)])
min(leida[2,2:ncol(leida)])

plot <-ggradar(plot.data =leida,
               grid.min = 0,
               grid.mid = 0.19,
               grid.max = 0.38, 
               values.radar = c(0,0.19,0.38),group.colours = c("#FFD4A9", "#80D1C8"),axis.label.size=4)
plot
colnames(leida)

ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_CHECK_plot.pdf",plot,width = 10,height = 10)
write.table(checkgene_sig2,"HNSC_checkgene_sig2.txt")



#KIRC####

#T cell exhaustion
# devtools::install_github("GuoBioinfoLab/TCellSI")
library(TCellSI)
ResultScores <- TCellSI::TCSS_Calculate(expr) 
#cluster12
gsva_cluster1<- ResultScores[,cluster2]
gsva_cluster2<- ResultScores[,cluster3]
gsva_cluster3<- ResultScores[,cluster1]

# cluster3_exp<-cluster3_exp[,-1]
# gsva_cluster3<- gsva(as.matrix(cluster3_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

###
gsva_zong<-cbind(gsva_cluster3, gsva_cluster2,gsva_cluster1)
gsva_zong<-gsva_zong[c("Progenitor_exhaustion","Terminal_exhaustion"),]
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(3,2,1),c(ncol(gsva_cluster3),ncol(gsva_cluster2),ncol(gsva_cluster1)))
for (j in 1:2) {
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
for (i in 1:nrow(gsva_zong)) {
  b<-as.numeric(gsva_zong[i,])
  tonglu1<-rep(rownames(gsva_zong)[i],ncol(gsva_zong))
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
mianyi_huatu <- mianyi_huatu %>%
  filter(mianyi_huatu$scores < quantile(mianyi_huatu$scores, 0.95)) # 95%
aaa<-ggplot(mianyi_huatu, aes(x = factor(type), y = scores, fill = class)) +
  geom_violin(position = position_dodge(0.8), trim = TRUE, alpha = 0.6) +  # alpha 
  geom_boxplot(width = 0.2, position = position_dodge(0.8), color = "black", 
               outlier.shape = NA, alpha = 0.5) +  # 
  scale_fill_manual(values=c("#FFD4A9","#80D1C8","pink")) +
  labs(x = "", y = 'Score') +
  stat_compare_means(method = "kruskal.test", label = "p.signif") +  # Kruskal-Wallis 
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  #  X 
    axis.text.y = element_text(size = 12),
    legend.title = element_blank(),  # 
    legend.text = element_text(size = 12)
  )

aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\KIRC_T_exhau_plot.pdf",aaa,width = 10,height = 10)

#CHECK POINT

checkpoint<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")
check_exp<-expr[checkpoint$Checkpoint.gene[which(checkpoint$Activator.A..Inhibitor.I.=="I")],]
check_exp<-expr[checkpoint$Checkpoint.gene,]

check_exp<-na.omit(check_exp)
cluster1_check<-check_exp[,cluster1]
cluster2_check<-check_exp[,cluster2]
cluster3_check<-check_exp[,cluster3]
cluster_check<-cbind(cluster1_check,cluster2_check,cluster3_check)
p_zanshi<-c()
mianyi_class<-rep(c(1,2,3),c(ncol(cluster1_check),ncol(cluster2_check),ncol(cluster3_check)))
for (j in 1:nrow(cluster_check)) {
  T_result<-kruskal.test(as.numeric(cluster_check[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
which(p_zanshi<0.05)

checkgene_sig<-rownames(cluster_check)[which(p_zanshi<0.05)]
checkgene_sig_p<-cbind(checkgene_sig,p_zanshi[which(p_zanshi<0.05)])
write.table(checkgene_sig_p,"KIRC_checkgene_sig_p.txt")

checkgene_sig2<-checkpoint[checkpoint$Checkpoint.gene %in%checkgene_sig,]
checkgene_sig2<-checkgene_sig2[order(checkgene_sig2$Activator.A..Inhibitor.I.), ]
cluster1_check<-check_exp[checkgene_sig2$Checkpoint.gene,cluster1]
cluster2_check<-check_exp[checkgene_sig2$Checkpoint.gene,cluster2]
cluster3_check<-check_exp[checkgene_sig2$Checkpoint.gene,cluster3]


cluster1_check<-apply(cluster1_check, 1, mean)
cluster2_check<-apply(cluster2_check, 1, mean)
cluster3_check<-apply(cluster3_check, 1, mean)

library(fmsb)
leida<-rbind(cluster1_check,cluster2_check,cluster3_check)
rownames(leida)<-c("cluster1","cluster2","cluster3")
library(ggradar)
classleida<-c("cluster1","cluster2","cluster3")
leida<-cbind(classleida,leida)
leida<-as.data.frame(leida)
str(leida)
leida[,2:ncol(leida)]<-apply(leida[,2:ncol(leida)],1,as.numeric)
max(leida[1,2:ncol(leida)])
min(leida[1,2:ncol(leida)])
max(leida[2,2:ncol(leida)])
min(leida[2,2:ncol(leida)])

plot <-ggradar(plot.data =leida,
               grid.min = 0,
               grid.mid = 0.19,
               grid.max = 0.38, 
               values.radar = c(0,0.19,0.38),group.colours = c("#FFD4A9", "#80D1C8","pink"),axis.label.size=4)
plot
colnames(leida)

ggsave("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\KIRC_CHECK_plot.pdf",plot,width = 10,height = 10)
write.table(checkgene_sig2,"KIRC_checkgene_sig2.txt")
#KIRP####

#T cell exhaustion
# devtools::install_github("GuoBioinfoLab/TCellSI")
library(TCellSI)
ResultScores <- TCellSI::TCSS_Calculate(expr) 
#cluster12
gsva_cluster1<- ResultScores[,cluster1]
gsva_cluster2<- ResultScores[,cluster2]
# cluster3_exp<-cluster3_exp[,-1]
# gsva_cluster3<- gsva(as.matrix(cluster3_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

###
gsva_zong<-cbind( gsva_cluster2,gsva_cluster1)
gsva_zong<-gsva_zong[c("Progenitor_exhaustion","Terminal_exhaustion"),]
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))
for (j in 1:2) {
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
for (i in 1:nrow(gsva_zong)) {
  b<-as.numeric(gsva_zong[i,])
  tonglu1<-rep(rownames(gsva_zong)[i],ncol(gsva_zong))
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
mianyi_huatu <- mianyi_huatu %>%
  filter(mianyi_huatu$scores < quantile(mianyi_huatu$scores, 0.95)) # 95%
aaa<-ggplot(mianyi_huatu, aes(x = factor(type), y = scores, fill = class)) +
  geom_violin(position = position_dodge(0.8), trim = TRUE, alpha = 0.6) +  # alpha 
  geom_boxplot(width = 0.2, position = position_dodge(0.8), color = "black", 
               outlier.shape = NA, alpha = 0.5) +  # 
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8")) +
  labs(x = "", y = 'Score') +
  stat_compare_means(method = "kruskal.test", label = "p.signif") +  # Kruskal-Wallis 
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  #  X 
    axis.text.y = element_text(size = 12),
    legend.title = element_blank(),  # 
    legend.text = element_text(size = 12)
  )

aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\KIRP_T_exhau_plot.pdf",aaa,width = 10,height = 10)

#CHECK POINT

checkpoint<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")
check_exp<-expr[checkpoint$Checkpoint.gene[which(checkpoint$Activator.A..Inhibitor.I.=="I")],]
check_exp<-expr[checkpoint$Checkpoint.gene,]

check_exp<-na.omit(check_exp)
cluster1_check<-check_exp[,cluster1]
cluster2_check<-check_exp[,cluster2]
cluster_check<-cbind(cluster1_check,cluster2_check)
p_zanshi<-c()
mianyi_class<-rep(c(1,2),c(ncol(cluster1_check),ncol(cluster2_check)))
for (j in 1:nrow(cluster_check)) {
  T_result<-kruskal.test(as.numeric(cluster_check[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
which(p_zanshi<0.05)

checkgene_sig<-rownames(cluster_check)[which(p_zanshi<0.05)]
checkgene_sig_p<-cbind(checkgene_sig,p_zanshi[which(p_zanshi<0.05)])
write.table(checkgene_sig_p,"KIRP_checkgene_sig_p.txt")

checkgene_sig2<-checkpoint[checkpoint$Checkpoint.gene %in%checkgene_sig,]
checkgene_sig2<-checkgene_sig2[order(checkgene_sig2$Activator.A..Inhibitor.I.), ]
cluster1_check<-check_exp[checkgene_sig2$Checkpoint.gene,cluster1]
cluster2_check<-check_exp[checkgene_sig2$Checkpoint.gene,cluster2]

cluster1_check<-apply(cluster1_check, 1, mean)
cluster2_check<-apply(cluster2_check, 1, mean)
library(fmsb)
leida<-rbind(cluster1_check,cluster2_check)
rownames(leida)<-c("cluster1","cluster2")
library(ggradar)
classleida<-c("cluster1","cluster2")
leida<-cbind(classleida,leida)
leida<-as.data.frame(leida)
str(leida)
leida[,2:ncol(leida)]<-apply(leida[,2:ncol(leida)],1,as.numeric)
max(leida[1,2:ncol(leida)])
min(leida[1,2:ncol(leida)])
max(leida[2,2:ncol(leida)])
min(leida[2,2:ncol(leida)])

plot <-ggradar(plot.data =leida,
               grid.min = 0,
               grid.mid = 0.17,
               grid.max = 0.34, 
               values.radar = c(0,0.17,0.34),group.colours = c("#FFD4A9", "#80D1C8"),axis.label.size=4)
plot
colnames(leida)

ggsave("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\KIRP_CHECK_plot.pdf",plot,width = 10,height = 10)
write.table(checkgene_sig2,"KIRP_checkgene_sig2.txt")

#LIHC####
#T cell exhaustion
# devtools::install_github("GuoBioinfoLab/TCellSI")
library(TCellSI)
ResultScores <- TCellSI::TCSS_Calculate(expr) 
#cluster12
gsva_cluster1<- ResultScores[,cluster1]
gsva_cluster2<- ResultScores[,cluster2]
# cluster3_exp<-cluster3_exp[,-1]
# gsva_cluster3<- gsva(as.matrix(cluster3_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

###
gsva_zong<-cbind( gsva_cluster2,gsva_cluster1)
gsva_zong<-gsva_zong[c("Progenitor_exhaustion","Terminal_exhaustion"),]
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))
for (j in 1:2) {
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
for (i in 1:nrow(gsva_zong)) {
  b<-as.numeric(gsva_zong[i,])
  tonglu1<-rep(rownames(gsva_zong)[i],ncol(gsva_zong))
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
mianyi_huatu <- mianyi_huatu %>%
  filter(mianyi_huatu$scores < quantile(mianyi_huatu$scores, 0.95)) # 95%
aaa<-ggplot(mianyi_huatu, aes(x = factor(type), y = scores, fill = class)) +
  geom_violin(position = position_dodge(0.8), trim = TRUE, alpha = 0.6) +  # alpha 
  geom_boxplot(width = 0.2, position = position_dodge(0.8), color = "black", 
               outlier.shape = NA, alpha = 0.5) +  # 
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8")) +
  labs(x = "", y = 'Score') +
  stat_compare_means(method = "kruskal.test", label = "p.signif") +  # Kruskal-Wallis 
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  #  X 
    axis.text.y = element_text(size = 12),
    legend.title = element_blank(),  # 
    legend.text = element_text(size = 12)
  )

aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\LIHC_T_exhau_plot.pdf",aaa,width = 10,height = 10)

#CHECK POINT

checkpoint<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")
check_exp<-expr[checkpoint$Checkpoint.gene[which(checkpoint$Activator.A..Inhibitor.I.=="I")],]
check_exp<-expr[checkpoint$Checkpoint.gene,]

check_exp<-na.omit(check_exp)
cluster1_check<-check_exp[,cluster1]
cluster2_check<-check_exp[,cluster2]
cluster_check<-cbind(cluster1_check,cluster2_check)
p_zanshi<-c()
mianyi_class<-rep(c(1,2),c(ncol(cluster1_check),ncol(cluster2_check)))
for (j in 1:nrow(cluster_check)) {
  T_result<-kruskal.test(as.numeric(cluster_check[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
which(p_zanshi<0.06)

checkgene_sig<-rownames(cluster_check)[which(p_zanshi<0.06)]
checkgene_sig_p<-cbind(checkgene_sig,p_zanshi[which(p_zanshi<0.06)])
write.table(checkgene_sig_p,"LIHC_checkgene_sig_p.txt")

checkgene_sig2<-checkpoint[checkpoint$Checkpoint.gene %in%checkgene_sig,]
checkgene_sig2<-checkgene_sig2[order(checkgene_sig2$Activator.A..Inhibitor.I.), ]
cluster1_check<-check_exp[checkgene_sig2$Checkpoint.gene,cluster1]
cluster2_check<-check_exp[checkgene_sig2$Checkpoint.gene,cluster2]

cluster1_check<-apply(cluster1_check, 1, mean)
cluster2_check<-apply(cluster2_check, 1, mean)
library(fmsb)
leida<-rbind(cluster1_check,cluster2_check)
rownames(leida)<-c("cluster1","cluster2")
library(ggradar)
classleida<-c("cluster1","cluster2")
leida<-cbind(classleida,leida)
leida<-as.data.frame(leida)
str(leida)
leida[,2:ncol(leida)]<-apply(leida[,2:ncol(leida)],1,as.numeric)
max(leida[1,2:ncol(leida)])
min(leida[1,2:ncol(leida)])
max(leida[2,2:ncol(leida)])
min(leida[2,2:ncol(leida)])

plot <-ggradar(plot.data =leida,
               grid.min = 0,
               grid.mid = 0.14,
               grid.max = 0.28, 
               values.radar = c(0,0.14,0.28),group.colours = c("#FFD4A9", "#80D1C8"),axis.label.size=4)
plot
colnames(leida)

ggsave("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\LIHC_CHECK_plot.pdf",plot,width = 10,height = 10)
write.table(checkgene_sig2,"LIHC_checkgene_sig2.txt")

#LUNG####
#T cell exhaustion
# devtools::install_github("GuoBioinfoLab/TCellSI")
library(TCellSI)
ResultScores <- TCellSI::TCSS_Calculate(expr) 
#cluster12
gsva_cluster1<- ResultScores[,cluster1]
gsva_cluster2<- ResultScores[,cluster2]
# cluster3_exp<-cluster3_exp[,-1]
# gsva_cluster3<- gsva(as.matrix(cluster3_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

###
gsva_zong<-cbind( gsva_cluster2,gsva_cluster1)
gsva_zong<-gsva_zong[c("Progenitor_exhaustion","Terminal_exhaustion"),]
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))
for (j in 1:2) {
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
for (i in 1:nrow(gsva_zong)) {
  b<-as.numeric(gsva_zong[i,])
  tonglu1<-rep(rownames(gsva_zong)[i],ncol(gsva_zong))
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
mianyi_huatu <- mianyi_huatu %>%
  filter(mianyi_huatu$scores < quantile(mianyi_huatu$scores, 0.95)) # 95%
aaa<-ggplot(mianyi_huatu, aes(x = factor(type), y = scores, fill = class)) +
  geom_violin(position = position_dodge(0.8), trim = TRUE, alpha = 0.6) +  # alpha 
  geom_boxplot(width = 0.2, position = position_dodge(0.8), color = "black", 
               outlier.shape = NA, alpha = 0.5) +  # 
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8")) +
  labs(x = "", y = 'Score') +
  stat_compare_means(method = "kruskal.test", label = "p.signif") +  # Kruskal-Wallis 
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  #  X 
    axis.text.y = element_text(size = 12),
    legend.title = element_blank(),  # 
    legend.text = element_text(size = 12)
  )

aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\NSCLC_T_exhau_plot.pdf",aaa,width = 10,height = 10)

#CHECK POINT

checkpoint<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")
check_exp<-expr[checkpoint$Checkpoint.gene[which(checkpoint$Activator.A..Inhibitor.I.=="I")],]
check_exp<-expr[checkpoint$Checkpoint.gene,]

check_exp<-na.omit(check_exp)
cluster1_check<-check_exp[,cluster1]
cluster2_check<-check_exp[,cluster2]
cluster_check<-cbind(cluster1_check,cluster2_check)
p_zanshi<-c()
mianyi_class<-rep(c(1,2),c(ncol(cluster1_check),ncol(cluster2_check)))
for (j in 1:nrow(cluster_check)) {
  T_result<-kruskal.test(as.numeric(cluster_check[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
which(p_zanshi<0.05)

checkgene_sig<-rownames(cluster_check)[which(p_zanshi<0.05)]
checkgene_sig_p<-cbind(checkgene_sig,p_zanshi[which(p_zanshi<0.05)])
write.table(checkgene_sig_p,"NSCLC_checkgene_sig_p.txt")

checkgene_sig2<-checkpoint[checkpoint$Checkpoint.gene %in%checkgene_sig,]
checkgene_sig2<-checkgene_sig2[order(checkgene_sig2$Activator.A..Inhibitor.I.), ]
cluster1_check<-check_exp[checkgene_sig2$Checkpoint.gene,cluster1]
cluster2_check<-check_exp[checkgene_sig2$Checkpoint.gene,cluster2]

cluster1_check<-apply(cluster1_check, 1, mean)
cluster2_check<-apply(cluster2_check, 1, mean)
library(fmsb)
leida<-rbind(cluster1_check,cluster2_check)
rownames(leida)<-c("cluster1","cluster2")
library(ggradar)
classleida<-c("cluster1","cluster2")
leida<-cbind(classleida,leida)
leida<-as.data.frame(leida)
str(leida)
leida[,2:ncol(leida)]<-apply(leida[,2:ncol(leida)],1,as.numeric)
max(leida[1,2:ncol(leida)])
min(leida[1,2:ncol(leida)])
max(leida[2,2:ncol(leida)])
min(leida[2,2:ncol(leida)])

plot <-ggradar(plot.data =leida,
               grid.min = 0,
               grid.mid = 0.21,
               grid.max = 0.42, 
               values.radar = c(0,0.21,0.42),group.colours = c("#FFD4A9", "#80D1C8"),axis.label.size=4)
plot
colnames(leida)

ggsave("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\NSCLC_CHECK_plot.pdf",plot,width = 10,height = 10)
write.table(checkgene_sig2,"NSCLC_checkgene_sig2.txt")
#PAAD####
#T cell exhaustion
# devtools::install_github("GuoBioinfoLab/TCellSI")
library(TCellSI)
ResultScores <- TCellSI::TCSS_Calculate(expr) 
#cluster12
gsva_cluster1<- ResultScores[,cluster1]
gsva_cluster2<- ResultScores[,cluster2]
# cluster3_exp<-cluster3_exp[,-1]
# gsva_cluster3<- gsva(as.matrix(cluster3_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

###
gsva_zong<-cbind( gsva_cluster2,gsva_cluster1)
gsva_zong<-gsva_zong[c("Progenitor_exhaustion","Terminal_exhaustion"),]
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))
for (j in 1:2) {
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
for (i in 1:nrow(gsva_zong)) {
  b<-as.numeric(gsva_zong[i,])
  tonglu1<-rep(rownames(gsva_zong)[i],ncol(gsva_zong))
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
mianyi_huatu <- mianyi_huatu %>%
  filter(mianyi_huatu$scores < quantile(mianyi_huatu$scores, 0.95)) # 95%
aaa<-ggplot(mianyi_huatu, aes(x = factor(type), y = scores, fill = class)) +
  geom_violin(position = position_dodge(0.8), trim = TRUE, alpha = 0.6) +  # alpha 
  geom_boxplot(width = 0.2, position = position_dodge(0.8), color = "black", 
               outlier.shape = NA, alpha = 0.5) +  # 
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8")) +
  labs(x = "", y = 'Score') +
  stat_compare_means(method = "kruskal.test", label = "p.signif") +  # Kruskal-Wallis 
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  #  X 
    axis.text.y = element_text(size = 12),
    legend.title = element_blank(),  # 
    legend.text = element_text(size = 12)
  )

aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\PAAD_T_exhau_plot.pdf",aaa,width = 10,height = 10)

#CHECK POINT

checkpoint<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")
check_exp<-expr[checkpoint$Checkpoint.gene[which(checkpoint$Activator.A..Inhibitor.I.=="I")],]
check_exp<-expr[checkpoint$Checkpoint.gene,]

check_exp<-na.omit(check_exp)
cluster1_check<-check_exp[,cluster1]
cluster2_check<-check_exp[,cluster2]
cluster_check<-cbind(cluster1_check,cluster2_check)
p_zanshi<-c()
mianyi_class<-rep(c(1,2),c(ncol(cluster1_check),ncol(cluster2_check)))
for (j in 1:nrow(cluster_check)) {
  T_result<-kruskal.test(as.numeric(cluster_check[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
which(p_zanshi<0.05)

checkgene_sig<-rownames(cluster_check)[which(p_zanshi<0.05)]
checkgene_sig_p<-cbind(checkgene_sig,p_zanshi[which(p_zanshi<0.05)])
write.table(checkgene_sig_p,"PAAD_checkgene_sig_p.txt")

checkgene_sig2<-checkpoint[checkpoint$Checkpoint.gene %in%checkgene_sig,]
checkgene_sig2<-checkgene_sig2[order(checkgene_sig2$Activator.A..Inhibitor.I.), ]
cluster1_check<-check_exp[checkgene_sig2$Checkpoint.gene,cluster1]
cluster2_check<-check_exp[checkgene_sig2$Checkpoint.gene,cluster2]

cluster1_check<-apply(cluster1_check, 1, mean)
cluster2_check<-apply(cluster2_check, 1, mean)
library(fmsb)
leida<-rbind(cluster1_check,cluster2_check)
rownames(leida)<-c("cluster1","cluster2")
library(ggradar)
classleida<-c("cluster1","cluster2")
leida<-cbind(classleida,leida)
leida<-as.data.frame(leida)
str(leida)
leida[,2:ncol(leida)]<-apply(leida[,2:ncol(leida)],1,as.numeric)
max(leida[1,2:ncol(leida)])
min(leida[1,2:ncol(leida)])
max(leida[2,2:ncol(leida)])
min(leida[2,2:ncol(leida)])

plot <-ggradar(plot.data =leida,
               grid.min = 0,
               grid.mid = 0.18,
               grid.max = 0.36, 
               values.radar = c(0,0.18,0.36),group.colours = c("#FFD4A9", "#80D1C8"),axis.label.size=4)
plot
colnames(leida)

ggsave("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\PAAD_CHECK_plot.pdf",plot,width = 10,height = 10)
write.table(checkgene_sig2,"PAAD_checkgene_sig2.txt")

#THCA####


# devtools::install_github("GuoBioinfoLab/TCellSI")
library(TCellSI)
ResultScores <- TCellSI::TCSS_Calculate(expr) 
#cluster12
gsva_cluster1<- ResultScores[,cluster1]
gsva_cluster2<- ResultScores[,cluster2]
# cluster3_exp<-cluster3_exp[,-1]
# gsva_cluster3<- gsva(as.matrix(cluster3_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

###
gsva_zong<-cbind( gsva_cluster2,gsva_cluster1)
gsva_zong<-gsva_zong[c("Progenitor_exhaustion","Terminal_exhaustion"),]
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))
for (j in 1:2) {
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
for (i in 1:nrow(gsva_zong)) {
  b<-as.numeric(gsva_zong[i,])
  tonglu1<-rep(rownames(gsva_zong)[i],ncol(gsva_zong))
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
mianyi_huatu <- mianyi_huatu %>%
  filter(mianyi_huatu$scores < quantile(mianyi_huatu$scores, 0.95)) # 95%
aaa<-ggplot(mianyi_huatu, aes(x = factor(type), y = scores, fill = class)) +
  geom_violin(position = position_dodge(0.8), trim = TRUE, alpha = 0.6) +  # alpha 
  geom_boxplot(width = 0.2, position = position_dodge(0.8), color = "black", 
               outlier.shape = NA, alpha = 0.5) +  # 
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8")) +
  labs(x = "", y = 'Score') +
  stat_compare_means(method = "kruskal.test", label = "p.signif") +  # Kruskal-Wallis 
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  #  X 
    axis.text.y = element_text(size = 12),
    legend.title = element_blank(),  # 
    legend.text = element_text(size = 12)
  )

aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\THCA_T_exhau_plot.pdf",aaa,width = 10,height = 10)

#CHECK POINT

checkpoint<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")
check_exp<-expr[checkpoint$Checkpoint.gene[which(checkpoint$Activator.A..Inhibitor.I.=="I")],]
check_exp<-expr[checkpoint$Checkpoint.gene,]

check_exp<-na.omit(check_exp)
cluster1_check<-check_exp[,cluster1]
cluster2_check<-check_exp[,cluster2]
cluster_check<-cbind(cluster1_check,cluster2_check)
p_zanshi<-c()
mianyi_class<-rep(c(1,2),c(ncol(cluster1_check),ncol(cluster2_check)))
for (j in 1:nrow(cluster_check)) {
  T_result<-kruskal.test(as.numeric(cluster_check[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
which(p_zanshi<0.05)

checkgene_sig<-rownames(cluster_check)[which(p_zanshi<0.05)]
checkgene_sig_p<-cbind(checkgene_sig,p_zanshi[which(p_zanshi<0.05)])
write.table(checkgene_sig_p,"THCA_checkgene_sig_p.txt")

checkgene_sig2<-checkpoint[checkpoint$Checkpoint.gene %in%checkgene_sig,]
checkgene_sig2<-checkgene_sig2[order(checkgene_sig2$Activator.A..Inhibitor.I.), ]
cluster1_check<-check_exp[checkgene_sig2$Checkpoint.gene,cluster1]
cluster2_check<-check_exp[checkgene_sig2$Checkpoint.gene,cluster2]

cluster1_check<-apply(cluster1_check, 1, mean)
cluster2_check<-apply(cluster2_check, 1, mean)
library(fmsb)
leida<-rbind(cluster1_check,cluster2_check)
rownames(leida)<-c("cluster1","cluster2")
library(ggradar)
classleida<-c("cluster1","cluster2")
leida<-cbind(classleida,leida)
leida<-as.data.frame(leida)
str(leida)
leida[,2:ncol(leida)]<-apply(leida[,2:ncol(leida)],1,as.numeric)
max(leida[1,2:ncol(leida)])
min(leida[1,2:ncol(leida)])
max(leida[2,2:ncol(leida)])
min(leida[2,2:ncol(leida)])

plot <-ggradar(plot.data =leida,
               grid.min = 0,
               grid.mid = 0.15,
               grid.max = 0.3, 
               values.radar = c(0,0.15,0.3),group.colours = c("#FFD4A9", "#80D1C8"),axis.label.size=4)
plot
colnames(leida)

ggsave("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\THCA_CHECK_plot.pdf",plot,width = 10,height = 10)
write.table(checkgene_sig2,"THCA_checkgene_sig2.txt")

#UCEC####

#T cell exhaustion
# devtools::install_github("GuoBioinfoLab/TCellSI")
library(TCellSI)
ResultScores <- TCellSI::TCSS_Calculate(expr) 
#cluster12
gsva_cluster1<- ResultScores[,cluster1]
gsva_cluster2<- ResultScores[,cluster2]
# cluster3_exp<-cluster3_exp[,-1]
# gsva_cluster3<- gsva(as.matrix(cluster3_exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

###
gsva_zong<-cbind( gsva_cluster2,gsva_cluster1)
gsva_zong<-gsva_zong[c("Progenitor_exhaustion","Terminal_exhaustion"),]
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_cluster2),ncol(gsva_cluster1)))
for (j in 1:2) {
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
for (i in 1:nrow(gsva_zong)) {
  b<-as.numeric(gsva_zong[i,])
  tonglu1<-rep(rownames(gsva_zong)[i],ncol(gsva_zong))
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
mianyi_huatu <- mianyi_huatu %>%
  filter(mianyi_huatu$scores < quantile(mianyi_huatu$scores, 0.95)) # 95%
aaa<-ggplot(mianyi_huatu, aes(x = factor(type), y = scores, fill = class)) +
  geom_violin(position = position_dodge(0.8), trim = TRUE, alpha = 0.6) +  # alpha 
  geom_boxplot(width = 0.2, position = position_dodge(0.8), color = "black", 
               outlier.shape = NA, alpha = 0.5) +  # 
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8")) +
  labs(x = "", y = 'Score') +
  stat_compare_means(method = "kruskal.test", label = "p.signif") +  # Kruskal-Wallis 
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  #  X 
    axis.text.y = element_text(size = 12),
    legend.title = element_blank(),  # 
    legend.text = element_text(size = 12)
  )

aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\UCEC_T_exhau_plot.pdf",aaa,width = 10,height = 10)

#CHECK POINT

checkpoint<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")
check_exp<-expr[checkpoint$Checkpoint.gene[which(checkpoint$Activator.A..Inhibitor.I.=="I")],]
check_exp<-expr[checkpoint$Checkpoint.gene,]

check_exp<-na.omit(check_exp)
cluster1_check<-check_exp[,cluster1]
cluster2_check<-check_exp[,cluster2]
cluster_check<-cbind(cluster1_check,cluster2_check)
p_zanshi<-c()
mianyi_class<-rep(c(1,2),c(ncol(cluster1_check),ncol(cluster2_check)))
for (j in 1:nrow(cluster_check)) {
  T_result<-kruskal.test(as.numeric(cluster_check[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
which(p_zanshi<0.05)

checkgene_sig<-rownames(cluster_check)[which(p_zanshi<0.05)]
checkgene_sig_p<-cbind(checkgene_sig,p_zanshi[which(p_zanshi<0.05)])
write.table(checkgene_sig_p,"G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\UCEC_checkgene_sig_p.txt")

checkgene_sig2<-checkpoint[checkpoint$Checkpoint.gene %in%checkgene_sig,]
checkgene_sig2<-checkgene_sig2[order(checkgene_sig2$Activator.A..Inhibitor.I.), ]
cluster1_check<-check_exp[checkgene_sig2$Checkpoint.gene,cluster1]
cluster2_check<-check_exp[checkgene_sig2$Checkpoint.gene,cluster2]

cluster1_check<-apply(cluster1_check, 1, mean)
cluster2_check<-apply(cluster2_check, 1, mean)
library(fmsb)
leida<-rbind(cluster1_check,cluster2_check)
rownames(leida)<-c("cluster1","cluster2")
library(ggradar)
classleida<-c("cluster1","cluster2")
leida<-cbind(classleida,leida)
leida<-as.data.frame(leida)
str(leida)
leida[,2:ncol(leida)]<-apply(leida[,2:ncol(leida)],1,as.numeric)
max(leida[1,2:ncol(leida)])
min(leida[1,2:ncol(leida)])
max(leida[2,2:ncol(leida)])
min(leida[2,2:ncol(leida)])

plot <-ggradar(plot.data =leida,
               grid.min = 0,
               grid.mid = 0.24,
               grid.max = 0.48, 
               values.radar = c(0,0.24,0.48),group.colours = c("#FFD4A9", "#80D1C8"),axis.label.size=4)
plot
colnames(leida)

ggsave("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\UCEC_CHECK_plot.pdf",plot,width = 10,height = 10)
write.table(checkgene_sig2,"G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\UCEC_checkgene_sig2.txt")


