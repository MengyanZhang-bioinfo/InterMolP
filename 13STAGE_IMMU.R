#BRCA####
#stage####
#T####
label<-read.csv("G:\\pan_result\\gpu_result_tcga\\BRCA\\label.csv")
gsva_stage1<- ResultScores[,label[which(label$label=="1"),1]]
gsva_stage2<- ResultScores[,intersect(label[which(label$label=="2"),1],colnames(ResultScores))]

###
gsva_zong<-cbind( gsva_stage2,gsva_stage1)
gsva_zong<-gsva_zong[c("Progenitor_exhaustion","Terminal_exhaustion"),]
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_stage2),ncol(gsva_stage1)))
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

class<-rep(c("stage2","stage1"),c(ncol(gsva_stage2),ncol(gsva_stage1)))

mianyi_huatu<-c()
for (i in 1:nrow(gsva_zong)) {
  b<-as.numeric(gsva_zong[i,])
  tonglu1<-rep(rownames(gsva_zong)[i],ncol(gsva_zong))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)
library(ggpubr)
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
  scale_fill_manual(values=c("#D55E00","#56B4E9")) +
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
ggsave("G:\\pan_result\\gpu_result_tcga\\all_stage_texh\\BRCA_T_exhau_plot_stage.pdf",aaa,width = 10,height = 10)
save(mianyi_huatu,file = "G:\\pan_result\\gpu_result_tcga\\all_stage_texh\\BRCA_T.Rdata")

#(4)####

checkpoint<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")
check_exp<-expr[checkpoint$Checkpoint.gene[which(checkpoint$Activator.A..Inhibitor.I.=="I")],]
check_exp<-expr[checkpoint$Checkpoint.gene,]

stage1<-label[which(label$label=="1"),1]
stage2<-label[which(label$label=="2"),1]
check_exp<-na.omit(check_exp)
stage1_check<-check_exp[,stage1]
stage2_check<-check_exp[,intersect(stage2,colnames(check_exp))]
stage_check<-cbind(stage1_check,stage2_check)
p_zanshi<-c()
mianyi_class<-rep(c(1,2),c(ncol(stage1_check),ncol(stage2_check)))
for (j in 1:nrow(stage_check)) {
  T_result<-kruskal.test(as.numeric(stage_check[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
which(p_zanshi<0.05)

checkgene_sig<-rownames(stage_check)[which(p_zanshi<0.05)]
checkgene_sig_p<-cbind(checkgene_sig,p_zanshi[which(p_zanshi<0.05)])
write.table(checkgene_sig_p,"G:\\pan_result\\gpu_result_tcga\\all_stage_checkpoint\\BRCAstage_checkgene_sig_p.txt")

checkgene_sig2<-checkpoint[checkpoint$Checkpoint.gene %in%checkgene_sig,]
checkgene_sig2<-checkgene_sig2[order(checkgene_sig2$Activator.A..Inhibitor.I.), ]
stage1_check<-check_exp[checkgene_sig2$Checkpoint.gene,stage1]
stage2_check<-check_exp[checkgene_sig2$Checkpoint.gene,intersect(stage2,colnames(check_exp))]

stage1_check<-apply(stage1_check, 1, mean)
stage2_check<-apply(stage2_check, 1, mean)
library(fmsb)
leida<-rbind(stage1_check,stage2_check)
rownames(leida)<-c("stage1","stage2")
library(ggradar)
classleida<-c("stage1","stage2")
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
               grid.mid = 128,
               grid.max = 256, 
               values.radar = c(0,128,256),group.colours = c("#D55E00","#56B4E9"),axis.label.size=4)
plot
colnames(leida)

ggsave("G:\\pan_result\\gpu_result_tcga\\all_stage_checkpoint\\BRCA_CHECK_plot_stage.pdf",plot,width = 10,height = 10)
save(leida,file = "G:\\pan_result\\gpu_result_tcga\\all_stage_checkpoint\\BRCA_leida.Rdata")
#COAD####
#stage####
#T####
label<-read.csv("G:\\pan_result\\gpu_result_tcga\\COAD\\label.csv")
gsva_stage1<- ResultScores[,label[which(label$label=="1"),1]]
gsva_stage2<- ResultScores[,intersect(label[which(label$label=="2"),1],colnames(ResultScores))]

###
gsva_zong<-cbind( gsva_stage2,gsva_stage1)
gsva_zong<-gsva_zong[c("Progenitor_exhaustion","Terminal_exhaustion"),]
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_stage2),ncol(gsva_stage1)))
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

class<-rep(c("stage2","stage1"),c(ncol(gsva_stage2),ncol(gsva_stage1)))

mianyi_huatu<-c()
for (i in 1:nrow(gsva_zong)) {
  b<-as.numeric(gsva_zong[i,])
  tonglu1<-rep(rownames(gsva_zong)[i],ncol(gsva_zong))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)
library(ggpubr)
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
  scale_fill_manual(values=c("#D55E00","#56B4E9")) +
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
ggsave("G:\\pan_result\\gpu_result_tcga\\all_stage_texh\\COAD_T_exhau_plot_stage.pdf",aaa,width = 10,height = 10)
save(mianyi_huatu,file = "G:\\pan_result\\gpu_result_tcga\\all_stage_texh\\COAD_T.Rdata")

#(4)####

checkpoint<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")
check_exp<-expr[checkpoint$Checkpoint.gene[which(checkpoint$Activator.A..Inhibitor.I.=="I")],]
check_exp<-expr[checkpoint$Checkpoint.gene,]

stage1<-label[which(label$label=="1"),1]
stage2<-label[which(label$label=="2"),1]
check_exp<-na.omit(check_exp)
stage1_check<-check_exp[,stage1]
stage2_check<-check_exp[,intersect(stage2,colnames(check_exp))]
stage_check<-cbind(stage1_check,stage2_check)
p_zanshi<-c()
mianyi_class<-rep(c(1,2),c(ncol(stage1_check),ncol(stage2_check)))
for (j in 1:nrow(stage_check)) {
  T_result<-kruskal.test(as.numeric(stage_check[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
which(p_zanshi<0.05)

checkgene_sig<-rownames(stage_check)[which(p_zanshi<0.05)]
checkgene_sig_p<-cbind(checkgene_sig,p_zanshi[which(p_zanshi<0.05)])
write.table(checkgene_sig_p,"G:\\pan_result\\gpu_result_tcga\\all_stage_checkpoint\\COADstage_checkgene_sig_p.txt")

checkgene_sig2<-checkpoint[checkpoint$Checkpoint.gene %in%checkgene_sig,]
checkgene_sig2<-checkgene_sig2[order(checkgene_sig2$Activator.A..Inhibitor.I.), ]
stage1_check<-check_exp[checkgene_sig2$Checkpoint.gene,stage1]
stage2_check<-check_exp[checkgene_sig2$Checkpoint.gene,intersect(stage2,colnames(check_exp))]

stage1_check<-apply(stage1_check, 1, mean)
stage2_check<-apply(stage2_check, 1, mean)
library(fmsb)
leida<-rbind(stage1_check,stage2_check)
rownames(leida)<-c("stage1","stage2")
library(ggradar)
classleida<-c("stage1","stage2")
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
               grid.mid = 128,
               grid.max = 256, 
               values.radar = c(0,128,256),group.colours = c("#D55E00","#56B4E9"),axis.label.size=4)
plot
colnames(leida)

ggsave("G:\\pan_result\\gpu_result_tcga\\all_stage_checkpoint\\COAD_CHECK_plot_stage.pdf",plot,width = 10,height = 10)
save(leida,file = "G:\\pan_result\\gpu_result_tcga\\all_stage_checkpoint\\COAD_leida.Rdata")

#ESCA####

#stage####
#T####
label<-read.csv("G:\\pan_result\\gpu_result_tcga\\ESCA\\label.csv")
gsva_stage1<- ResultScores[,label[which(label$label=="1"),1]]
gsva_stage2<- ResultScores[,intersect(label[which(label$label=="2"),1],colnames(ResultScores))]

###
gsva_zong<-cbind( gsva_stage2,gsva_stage1)
gsva_zong<-gsva_zong[c("Progenitor_exhaustion","Terminal_exhaustion"),]
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_stage2),ncol(gsva_stage1)))
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

class<-rep(c("stage2","stage1"),c(ncol(gsva_stage2),ncol(gsva_stage1)))

mianyi_huatu<-c()
for (i in 1:nrow(gsva_zong)) {
  b<-as.numeric(gsva_zong[i,])
  tonglu1<-rep(rownames(gsva_zong)[i],ncol(gsva_zong))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)
library(ggpubr)
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
  scale_fill_manual(values=c("#D55E00","#56B4E9")) +
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
ggsave("G:\\pan_result\\gpu_result_tcga\\all_stage_texh\\ESCA_T_exhau_plot_stage.pdf",aaa,width = 10,height = 10)
save(mianyi_huatu,file = "G:\\pan_result\\gpu_result_tcga\\all_stage_texh\\ESCA_T.Rdata")

#(4)####

checkpoint<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")
check_exp<-expr[checkpoint$Checkpoint.gene[which(checkpoint$Activator.A..Inhibitor.I.=="I")],]
check_exp<-expr[checkpoint$Checkpoint.gene,]

stage1<-label[which(label$label=="1"),1]
stage2<-label[which(label$label=="2"),1]
check_exp<-na.omit(check_exp)
stage1_check<-check_exp[,stage1]
stage2_check<-check_exp[,intersect(stage2,colnames(check_exp))]
stage_check<-cbind(stage1_check,stage2_check)
p_zanshi<-c()
mianyi_class<-rep(c(1,2),c(ncol(stage1_check),ncol(stage2_check)))
for (j in 1:nrow(stage_check)) {
  T_result<-kruskal.test(as.numeric(stage_check[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
which(p_zanshi<0.05)

checkgene_sig<-rownames(stage_check)[which(p_zanshi<0.05)]
checkgene_sig_p<-cbind(checkgene_sig,p_zanshi[which(p_zanshi<0.05)])
write.table(checkgene_sig_p,"G:\\pan_result\\gpu_result_tcga\\all_stage_checkpoint\\ESCAstage_checkgene_sig_p.txt")

checkgene_sig2<-checkpoint[checkpoint$Checkpoint.gene %in%checkgene_sig,]
checkgene_sig2<-checkgene_sig2[order(checkgene_sig2$Activator.A..Inhibitor.I.), ]
stage1_check<-check_exp[checkgene_sig2$Checkpoint.gene,stage1]
stage2_check<-check_exp[checkgene_sig2$Checkpoint.gene,intersect(stage2,colnames(check_exp))]

stage1_check<-apply(stage1_check, 1, mean)
stage2_check<-apply(stage2_check, 1, mean)
library(fmsb)
leida<-rbind(stage1_check,stage2_check)
rownames(leida)<-c("stage1","stage2")
library(ggradar)
classleida<-c("stage1","stage2")
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
               grid.mid = 0.2,
               grid.max = 0.4, 
               values.radar = c(0,0.2,0.4),group.colours = c("#D55E00","#56B4E9"),axis.label.size=4)
plot
colnames(leida)

ggsave("G:\\pan_result\\gpu_result_tcga\\all_stage_checkpoint\\ESCA_CHECK_plot_stage.pdf",plot,width = 10,height = 10)
save(leida,file = "G:\\pan_result\\gpu_result_tcga\\all_stage_checkpoint\\ESCA_leida.Rdata")

#HNSC####
#stage####
#T####
label<-read.csv("G:\\pan_result\\gpu_result_tcga\\HNSC\\label.csv")
gsva_stage1<- ResultScores[,label[which(label$label=="1"),1]]
gsva_stage2<- ResultScores[,intersect(label[which(label$label=="2"),1],colnames(ResultScores))]

###
gsva_zong<-cbind( gsva_stage2,gsva_stage1)
gsva_zong<-gsva_zong[c("Progenitor_exhaustion","Terminal_exhaustion"),]
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_stage2),ncol(gsva_stage1)))
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

class<-rep(c("stage2","stage1"),c(ncol(gsva_stage2),ncol(gsva_stage1)))

mianyi_huatu<-c()
for (i in 1:nrow(gsva_zong)) {
  b<-as.numeric(gsva_zong[i,])
  tonglu1<-rep(rownames(gsva_zong)[i],ncol(gsva_zong))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)
library(ggpubr)
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
  scale_fill_manual(values=c("#D55E00","#56B4E9")) +
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
ggsave("G:\\pan_result\\gpu_result_tcga\\all_stage_texh\\HNSC_T_exhau_plot_stage.pdf",aaa,width = 10,height = 10)
save(mianyi_huatu,file = "G:\\pan_result\\gpu_result_tcga\\all_stage_texh\\HNSC_T.Rdata")

#(4)####

checkpoint<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")
check_exp<-expr[checkpoint$Checkpoint.gene[which(checkpoint$Activator.A..Inhibitor.I.=="I")],]
check_exp<-expr[checkpoint$Checkpoint.gene,]

stage1<-label[which(label$label=="1"),1]
stage2<-label[which(label$label=="2"),1]
check_exp<-na.omit(check_exp)
stage1_check<-check_exp[,stage1]
stage2_check<-check_exp[,intersect(stage2,colnames(check_exp))]
stage_check<-cbind(stage1_check,stage2_check)
p_zanshi<-c()
mianyi_class<-rep(c(1,2),c(ncol(stage1_check),ncol(stage2_check)))
for (j in 1:nrow(stage_check)) {
  T_result<-kruskal.test(as.numeric(stage_check[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
which(p_zanshi<0.05)

checkgene_sig<-rownames(stage_check)[which(p_zanshi<0.05)]
checkgene_sig_p<-cbind(checkgene_sig,p_zanshi[which(p_zanshi<0.05)])
write.table(checkgene_sig_p,"G:\\pan_result\\gpu_result_tcga\\all_stage_checkpoint\\HNSCstage_checkgene_sig_p.txt")

checkgene_sig2<-checkpoint[checkpoint$Checkpoint.gene %in%checkgene_sig,]
checkgene_sig2<-checkgene_sig2[order(checkgene_sig2$Activator.A..Inhibitor.I.), ]
stage1_check<-check_exp[checkgene_sig2$Checkpoint.gene,stage1]
stage2_check<-check_exp[checkgene_sig2$Checkpoint.gene,intersect(stage2,colnames(check_exp))]

stage1_check<-apply(stage1_check, 1, mean)
stage2_check<-apply(stage2_check, 1, mean)
library(fmsb)
leida<-rbind(stage1_check,stage2_check)
rownames(leida)<-c("stage1","stage2")
library(ggradar)
classleida<-c("stage1","stage2")
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
               values.radar = c(0,0.18,0.36),group.colours = c("#D55E00","#56B4E9"),axis.label.size=4)
plot
colnames(leida)

ggsave("G:\\pan_result\\gpu_result_tcga\\all_stage_checkpoint\\HNSC_CHECK_plot_stage.pdf",plot,width = 10,height = 10)
save(leida,file = "G:\\pan_result\\gpu_result_tcga\\all_stage_checkpoint\\HNSC_leida.Rdata")

#KIRC####
#stage####
#T####
label<-read.csv("G:\\pan_result\\gpu_result_tcga\\KIRC\\label.csv")
gsva_stage1<- ResultScores[,label[which(label$label=="1"),1]]
gsva_stage2<- ResultScores[,intersect(label[which(label$label=="2"),1],colnames(ResultScores))]

###
gsva_zong<-cbind( gsva_stage2,gsva_stage1)
gsva_zong<-gsva_zong[c("Progenitor_exhaustion","Terminal_exhaustion"),]
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_stage2),ncol(gsva_stage1)))
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

class<-rep(c("stage2","stage1"),c(ncol(gsva_stage2),ncol(gsva_stage1)))

mianyi_huatu<-c()
for (i in 1:nrow(gsva_zong)) {
  b<-as.numeric(gsva_zong[i,])
  tonglu1<-rep(rownames(gsva_zong)[i],ncol(gsva_zong))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)
library(ggpubr)
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
  scale_fill_manual(values=c("#D55E00","#56B4E9")) +
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
ggsave("G:\\pan_result\\gpu_result_tcga\\all_stage_texh\\KIRC_T_exhau_plot_stage.pdf",aaa,width = 10,height = 10)
save(mianyi_huatu,file = "G:\\pan_result\\gpu_result_tcga\\all_stage_texh\\KIRC_T.Rdata")

#(4)####

checkpoint<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")
check_exp<-expr[checkpoint$Checkpoint.gene[which(checkpoint$Activator.A..Inhibitor.I.=="I")],]
check_exp<-expr[checkpoint$Checkpoint.gene,]

stage1<-label[which(label$label=="1"),1]
stage2<-label[which(label$label=="2"),1]
check_exp<-na.omit(check_exp)
stage1_check<-check_exp[,stage1]
stage2_check<-check_exp[,intersect(stage2,colnames(check_exp))]
stage_check<-cbind(stage1_check,stage2_check)
p_zanshi<-c()
mianyi_class<-rep(c(1,2),c(ncol(stage1_check),ncol(stage2_check)))
for (j in 1:nrow(stage_check)) {
  T_result<-kruskal.test(as.numeric(stage_check[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
which(p_zanshi<0.05)

checkgene_sig<-rownames(stage_check)[which(p_zanshi<0.05)]
checkgene_sig_p<-cbind(checkgene_sig,p_zanshi[which(p_zanshi<0.05)])
write.table(checkgene_sig_p,"G:\\pan_result\\gpu_result_tcga\\all_stage_checkpoint\\KIRCstage_checkgene_sig_p.txt")

checkgene_sig2<-checkpoint[checkpoint$Checkpoint.gene %in%checkgene_sig,]
checkgene_sig2<-checkgene_sig2[order(checkgene_sig2$Activator.A..Inhibitor.I.), ]
stage1_check<-check_exp[checkgene_sig2$Checkpoint.gene,stage1]
stage2_check<-check_exp[checkgene_sig2$Checkpoint.gene,intersect(stage2,colnames(check_exp))]

stage1_check<-apply(stage1_check, 1, mean)
stage2_check<-apply(stage2_check, 1, mean)
library(fmsb)
leida<-rbind(stage1_check,stage2_check)
rownames(leida)<-c("stage1","stage2")
library(ggradar)
classleida<-c("stage1","stage2")
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
               values.radar = c(0,0.15,0.3),group.colours = c("#D55E00","#56B4E9"),axis.label.size=4)
plot
colnames(leida)

ggsave("G:\\pan_result\\gpu_result_tcga\\all_stage_checkpoint\\KIRC_CHECK_plot_stage.pdf",plot,width = 10,height = 10)
save(leida,file = "G:\\pan_result\\gpu_result_tcga\\all_stage_checkpoint\\KIRC_leida.Rdata")

#KIRP####

#stage####
#T####
label<-read.csv("G:\\pan_result\\gpu_result_tcga\\KIRP\\label.csv")
gsva_stage1<- ResultScores[,label[which(label$label=="1"),1]]
gsva_stage2<- ResultScores[,intersect(label[which(label$label=="2"),1],colnames(ResultScores))]

###
gsva_zong<-cbind( gsva_stage2,gsva_stage1)
gsva_zong<-gsva_zong[c("Progenitor_exhaustion","Terminal_exhaustion"),]
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_stage2),ncol(gsva_stage1)))
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

class<-rep(c("stage2","stage1"),c(ncol(gsva_stage2),ncol(gsva_stage1)))

mianyi_huatu<-c()
for (i in 1:nrow(gsva_zong)) {
  b<-as.numeric(gsva_zong[i,])
  tonglu1<-rep(rownames(gsva_zong)[i],ncol(gsva_zong))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)
library(ggpubr)
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
  scale_fill_manual(values=c("#D55E00","#56B4E9")) +
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
ggsave("G:\\pan_result\\gpu_result_tcga\\all_stage_texh\\KIRP_T_exhau_plot_stage.pdf",aaa,width = 10,height = 10)
save(mianyi_huatu,file = "G:\\pan_result\\gpu_result_tcga\\all_stage_texh\\KIRP_T.Rdata")

#(4)####

checkpoint<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")
check_exp<-expr[checkpoint$Checkpoint.gene[which(checkpoint$Activator.A..Inhibitor.I.=="I")],]
check_exp<-expr[checkpoint$Checkpoint.gene,]

stage1<-label[which(label$label=="1"),1]
stage2<-label[which(label$label=="2"),1]
check_exp<-na.omit(check_exp)
stage1_check<-check_exp[,stage1]
stage2_check<-check_exp[,intersect(stage2,colnames(check_exp))]
stage_check<-cbind(stage1_check,stage2_check)
p_zanshi<-c()
mianyi_class<-rep(c(1,2),c(ncol(stage1_check),ncol(stage2_check)))
for (j in 1:nrow(stage_check)) {
  T_result<-kruskal.test(as.numeric(stage_check[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
which(p_zanshi<0.05)

checkgene_sig<-rownames(stage_check)[which(p_zanshi<0.05)]
checkgene_sig_p<-cbind(checkgene_sig,p_zanshi[which(p_zanshi<0.05)])
write.table(checkgene_sig_p,"G:\\pan_result\\gpu_result_tcga\\all_stage_checkpoint\\KIRPstage_checkgene_sig_p.txt")

checkgene_sig2<-checkpoint[checkpoint$Checkpoint.gene %in%checkgene_sig,]
checkgene_sig2<-checkgene_sig2[order(checkgene_sig2$Activator.A..Inhibitor.I.), ]
stage1_check<-check_exp[checkgene_sig2$Checkpoint.gene,stage1]
stage2_check<-check_exp[checkgene_sig2$Checkpoint.gene,intersect(stage2,colnames(check_exp))]

stage1_check<-apply(stage1_check, 1, mean)
stage2_check<-apply(stage2_check, 1, mean)
library(fmsb)
leida<-rbind(stage1_check,stage2_check)
rownames(leida)<-c("stage1","stage2")
library(ggradar)
classleida<-c("stage1","stage2")
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
               values.radar = c(0,0.15,0.3),group.colours = c("#D55E00","#56B4E9"),axis.label.size=4)
plot
colnames(leida)

ggsave("G:\\pan_result\\gpu_result_tcga\\all_stage_checkpoint\\KIRP_CHECK_plot_stage.pdf",plot,width = 10,height = 10)
save(leida,file = "G:\\pan_result\\gpu_result_tcga\\all_stage_checkpoint\\KIRP_leida.Rdata")

#LIHC####
#stage####
#T####
label<-read.csv("G:\\pan_result\\gpu_result_tcga\\LIHC\\label.csv")
gsva_stage1<- ResultScores[,label[which(label$label=="1"),1]]
gsva_stage2<- ResultScores[,intersect(label[which(label$label=="2"),1],colnames(ResultScores))]

###
gsva_zong<-cbind( gsva_stage2,gsva_stage1)
gsva_zong<-gsva_zong[c("Progenitor_exhaustion","Terminal_exhaustion"),]
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_stage2),ncol(gsva_stage1)))
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

class<-rep(c("stage2","stage1"),c(ncol(gsva_stage2),ncol(gsva_stage1)))

mianyi_huatu<-c()
for (i in 1:nrow(gsva_zong)) {
  b<-as.numeric(gsva_zong[i,])
  tonglu1<-rep(rownames(gsva_zong)[i],ncol(gsva_zong))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)
library(ggpubr)
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
  scale_fill_manual(values=c("#D55E00","#56B4E9")) +
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
ggsave("G:\\pan_result\\gpu_result_tcga\\all_stage_texh\\LIHC_T_exhau_plot_stage.pdf",aaa,width = 10,height = 10)
save(mianyi_huatu,file = "G:\\pan_result\\gpu_result_tcga\\all_stage_texh\\LIHC_T.Rdata")

#(4)####

checkpoint<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")
check_exp<-expr[checkpoint$Checkpoint.gene[which(checkpoint$Activator.A..Inhibitor.I.=="I")],]
check_exp<-expr[checkpoint$Checkpoint.gene,]

stage1<-label[which(label$label=="1"),1]
stage2<-label[which(label$label=="2"),1]
check_exp<-na.omit(check_exp)
stage1_check<-check_exp[,stage1]
stage2_check<-check_exp[,intersect(stage2,colnames(check_exp))]
stage_check<-cbind(stage1_check,stage2_check)
p_zanshi<-c()
mianyi_class<-rep(c(1,2),c(ncol(stage1_check),ncol(stage2_check)))
for (j in 1:nrow(stage_check)) {
  T_result<-kruskal.test(as.numeric(stage_check[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
which(p_zanshi<0.05)

checkgene_sig<-rownames(stage_check)[which(p_zanshi<0.05)]
checkgene_sig_p<-cbind(checkgene_sig,p_zanshi[which(p_zanshi<0.05)])
write.table(checkgene_sig_p,"G:\\pan_result\\gpu_result_tcga\\all_stage_checkpoint\\LIHCstage_checkgene_sig_p.txt")

checkgene_sig2<-checkpoint[checkpoint$Checkpoint.gene %in%checkgene_sig,]
checkgene_sig2<-checkgene_sig2[order(checkgene_sig2$Activator.A..Inhibitor.I.), ]
stage1_check<-check_exp[checkgene_sig2$Checkpoint.gene,stage1]
stage2_check<-check_exp[checkgene_sig2$Checkpoint.gene,intersect(stage2,colnames(check_exp))]

stage1_check<-apply(stage1_check, 1, mean)
stage2_check<-apply(stage2_check, 1, mean)
library(fmsb)
leida<-rbind(stage1_check,stage2_check)
rownames(leida)<-c("stage1","stage2")
library(ggradar)
classleida<-c("stage1","stage2")
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
               values.radar = c(0,0.15,0.3),group.colours = c("#D55E00","#56B4E9"),axis.label.size=4)
plot
colnames(leida)

ggsave("G:\\pan_result\\gpu_result_tcga\\all_stage_checkpoint\\LIHC_CHECK_plot_stage.pdf",plot,width = 10,height = 10)
save(leida,file = "G:\\pan_result\\gpu_result_tcga\\all_stage_checkpoint\\LIHC_leida.Rdata")

#LUNG####
#stage####
#T####
label<-read.csv("G:\\pan_result\\gpu_result_tcga\\NSCLC\\label.csv")
gsva_stage1<- ResultScores[,label[which(label$label=="1"),1]]
gsva_stage2<- ResultScores[,intersect(label[which(label$label=="2"),1],colnames(ResultScores))]

###
gsva_zong<-cbind( gsva_stage2,gsva_stage1)
gsva_zong<-gsva_zong[c("Progenitor_exhaustion","Terminal_exhaustion"),]
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_stage2),ncol(gsva_stage1)))
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

class<-rep(c("stage2","stage1"),c(ncol(gsva_stage2),ncol(gsva_stage1)))

mianyi_huatu<-c()
for (i in 1:nrow(gsva_zong)) {
  b<-as.numeric(gsva_zong[i,])
  tonglu1<-rep(rownames(gsva_zong)[i],ncol(gsva_zong))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)
library(ggpubr)
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
  scale_fill_manual(values=c("#D55E00","#56B4E9")) +
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
ggsave("G:\\pan_result\\gpu_result_tcga\\all_stage_texh\\NSCLC_T_exhau_plot_stage.pdf",aaa,width = 10,height = 10)
save(mianyi_huatu,file = "G:\\pan_result\\gpu_result_tcga\\all_stage_texh\\NSCLC_T.Rdata")

#(4)####

checkpoint<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")
check_exp<-expr[checkpoint$Checkpoint.gene[which(checkpoint$Activator.A..Inhibitor.I.=="I")],]
check_exp<-expr[checkpoint$Checkpoint.gene,]

stage1<-label[which(label$label=="1"),1]
stage2<-label[which(label$label=="2"),1]
check_exp<-na.omit(check_exp)
stage1_check<-check_exp[,stage1]
stage2_check<-check_exp[,intersect(stage2,colnames(check_exp))]
stage_check<-cbind(stage1_check,stage2_check)
p_zanshi<-c()
mianyi_class<-rep(c(1,2),c(ncol(stage1_check),ncol(stage2_check)))
for (j in 1:nrow(stage_check)) {
  T_result<-kruskal.test(as.numeric(stage_check[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
which(p_zanshi<0.05)

checkgene_sig<-rownames(stage_check)[which(p_zanshi<0.05)]
checkgene_sig_p<-cbind(checkgene_sig,p_zanshi[which(p_zanshi<0.05)])
write.table(checkgene_sig_p,"G:\\pan_result\\gpu_result_tcga\\all_stage_checkpoint\\NSCLCstage_checkgene_sig_p.txt")

checkgene_sig2<-checkpoint[checkpoint$Checkpoint.gene %in%checkgene_sig,]
checkgene_sig2<-checkgene_sig2[order(checkgene_sig2$Activator.A..Inhibitor.I.), ]
stage1_check<-check_exp[checkgene_sig2$Checkpoint.gene,stage1]
stage2_check<-check_exp[checkgene_sig2$Checkpoint.gene,intersect(stage2,colnames(check_exp))]

stage1_check<-apply(stage1_check, 1, mean)
stage2_check<-apply(stage2_check, 1, mean)
library(fmsb)
leida<-rbind(stage1_check,stage2_check)
rownames(leida)<-c("stage1","stage2")
library(ggradar)
classleida<-c("stage1","stage2")
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
               values.radar = c(0,0.21,0.42),group.colours = c("#D55E00","#56B4E9"),axis.label.size=4)
plot
colnames(leida)

ggsave("G:\\pan_result\\gpu_result_tcga\\all_stage_checkpoint\\NSCLC_CHECK_plot_stage.pdf",plot,width = 10,height = 10)
save(leida,file = "G:\\pan_result\\gpu_result_tcga\\all_stage_checkpoint\\NSCLC_leida.Rdata")


#PAAD####

#stage####
#T####
label<-read.csv("G:\\pan_result\\gpu_result_tcga\\PAAD\\label.csv")
gsva_stage1<- ResultScores[,label[which(label$label=="1"),1]]
gsva_stage2<- ResultScores[,intersect(label[which(label$label=="2"),1],colnames(ResultScores))]

###
gsva_zong<-cbind( gsva_stage2,gsva_stage1)
gsva_zong<-gsva_zong[c("Progenitor_exhaustion","Terminal_exhaustion"),]
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_stage2),ncol(gsva_stage1)))
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

class<-rep(c("stage2","stage1"),c(ncol(gsva_stage2),ncol(gsva_stage1)))

mianyi_huatu<-c()
for (i in 1:nrow(gsva_zong)) {
  b<-as.numeric(gsva_zong[i,])
  tonglu1<-rep(rownames(gsva_zong)[i],ncol(gsva_zong))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)
library(ggpubr)
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
  scale_fill_manual(values=c("#D55E00","#56B4E9")) +
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
ggsave("G:\\pan_result\\gpu_result_tcga\\all_stage_texh\\PAAD_T_exhau_plot_stage.pdf",aaa,width = 10,height = 10)
save(mianyi_huatu,file = "G:\\pan_result\\gpu_result_tcga\\all_stage_texh\\PAAD_T.Rdata")

#(4)####

checkpoint<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")
check_exp<-expr[checkpoint$Checkpoint.gene[which(checkpoint$Activator.A..Inhibitor.I.=="I")],]
check_exp<-expr[checkpoint$Checkpoint.gene,]

stage1<-label[which(label$label=="1"),1]
stage2<-label[which(label$label=="2"),1]
check_exp<-na.omit(check_exp)
stage1_check<-check_exp[,stage1]
stage2_check<-check_exp[,intersect(stage2,colnames(check_exp))]
stage_check<-cbind(stage1_check,stage2_check)
p_zanshi<-c()
mianyi_class<-rep(c(1,2),c(ncol(stage1_check),ncol(stage2_check)))
for (j in 1:nrow(stage_check)) {
  T_result<-kruskal.test(as.numeric(stage_check[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
which(p_zanshi<0.05)

checkgene_sig<-rownames(stage_check)[which(p_zanshi<0.05)]
checkgene_sig_p<-cbind(checkgene_sig,p_zanshi[which(p_zanshi<0.05)])
write.table(checkgene_sig_p,"G:\\pan_result\\gpu_result_tcga\\all_stage_checkpoint\\PAADstage_checkgene_sig_p.txt")

checkgene_sig2<-checkpoint[checkpoint$Checkpoint.gene %in%checkgene_sig,]
checkgene_sig2<-checkgene_sig2[order(checkgene_sig2$Activator.A..Inhibitor.I.), ]
stage1_check<-check_exp[checkgene_sig2$Checkpoint.gene,stage1]
stage2_check<-check_exp[checkgene_sig2$Checkpoint.gene,intersect(stage2,colnames(check_exp))]

stage1_check<-apply(stage1_check, 1, mean)
stage2_check<-apply(stage2_check, 1, mean)
library(fmsb)
leida<-rbind(stage1_check,stage2_check)
rownames(leida)<-c("stage1","stage2")
library(ggradar)
classleida<-c("stage1","stage2")
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
               grid.mid = 0.22,
               grid.max = 0.44, 
               values.radar = c(0,0.22,0.44),group.colours = c("#D55E00","#56B4E9"),axis.label.size=4)
plot
colnames(leida)

ggsave("G:\\pan_result\\gpu_result_tcga\\all_stage_checkpoint\\PAAD_CHECK_plot_stage.pdf",plot,width = 10,height = 10)
save(leida,file = "G:\\pan_result\\gpu_result_tcga\\all_stage_checkpoint\\PAAD_leida.Rdata")

#THCA####
#stage####
#T####
label<-read.csv("G:\\pan_result\\gpu_result_tcga\\THCA\\label.csv")
gsva_stage1<- ResultScores[,label[which(label$label=="1"),1]]
gsva_stage2<- ResultScores[,intersect(label[which(label$label=="2"),1],colnames(ResultScores))]

###
gsva_zong<-cbind( gsva_stage2,gsva_stage1)
gsva_zong<-gsva_zong[c("Progenitor_exhaustion","Terminal_exhaustion"),]
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_stage2),ncol(gsva_stage1)))
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

class<-rep(c("stage2","stage1"),c(ncol(gsva_stage2),ncol(gsva_stage1)))

mianyi_huatu<-c()
for (i in 1:nrow(gsva_zong)) {
  b<-as.numeric(gsva_zong[i,])
  tonglu1<-rep(rownames(gsva_zong)[i],ncol(gsva_zong))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)
library(ggpubr)
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
  scale_fill_manual(values=c("#D55E00","#56B4E9")) +
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
ggsave("G:\\pan_result\\gpu_result_tcga\\all_stage_texh\\THCA_T_exhau_plot_stage.pdf",aaa,width = 10,height = 10)
save(mianyi_huatu,file = "G:\\pan_result\\gpu_result_tcga\\all_stage_texh\\THCA_T.Rdata")

#(4)####

checkpoint<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")
check_exp<-expr[checkpoint$Checkpoint.gene[which(checkpoint$Activator.A..Inhibitor.I.=="I")],]
check_exp<-expr[checkpoint$Checkpoint.gene,]

stage1<-label[which(label$label=="1"),1]
stage2<-label[which(label$label=="2"),1]
check_exp<-na.omit(check_exp)
stage1_check<-check_exp[,stage1]
stage2_check<-check_exp[,intersect(stage2,colnames(check_exp))]
stage_check<-cbind(stage1_check,stage2_check)
p_zanshi<-c()
mianyi_class<-rep(c(1,2),c(ncol(stage1_check),ncol(stage2_check)))
for (j in 1:nrow(stage_check)) {
  T_result<-kruskal.test(as.numeric(stage_check[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
which(p_zanshi<0.05)

checkgene_sig<-rownames(stage_check)[which(p_zanshi<0.05)]
checkgene_sig_p<-cbind(checkgene_sig,p_zanshi[which(p_zanshi<0.05)])
write.table(checkgene_sig_p,"G:\\pan_result\\gpu_result_tcga\\all_stage_checkpoint\\THCAstage_checkgene_sig_p.txt")

checkgene_sig2<-checkpoint[checkpoint$Checkpoint.gene %in%checkgene_sig,]
checkgene_sig2<-checkgene_sig2[order(checkgene_sig2$Activator.A..Inhibitor.I.), ]
stage1_check<-check_exp[checkgene_sig2$Checkpoint.gene,stage1]
stage2_check<-check_exp[checkgene_sig2$Checkpoint.gene,intersect(stage2,colnames(check_exp))]

stage1_check<-apply(stage1_check, 1, mean)
stage2_check<-apply(stage2_check, 1, mean)
library(fmsb)
leida<-rbind(stage1_check,stage2_check)
rownames(leida)<-c("stage1","stage2")
library(ggradar)
classleida<-c("stage1","stage2")
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
               grid.mid = 0.16,
               grid.max = 0.32, 
               values.radar = c(0,0.16,0.32),group.colours = c("#D55E00","#56B4E9"),axis.label.size=4)
plot
colnames(leida)

ggsave("G:\\pan_result\\gpu_result_tcga\\all_stage_checkpoint\\THCA_CHECK_plot_stage.pdf",plot,width = 10,height = 10)
save(leida,file = "G:\\pan_result\\gpu_result_tcga\\all_stage_checkpoint\\THCA_leida.Rdata")


#UCEC####
#stage####
#T####
label<-read.csv("G:\\pan_result\\gpu_result_tcga\\UCEC\\label.csv")
gsva_stage1<- ResultScores[,label[which(label$label=="1"),1]]
gsva_stage2<- ResultScores[,intersect(label[which(label$label=="2"),1],colnames(ResultScores))]

###
gsva_zong<-cbind( gsva_stage2,gsva_stage1)
gsva_zong<-gsva_zong[c("Progenitor_exhaustion","Terminal_exhaustion"),]
str(gsva_zong)
p_zanshi<-c()
mianyi_class<-rep(c(2,1),c(ncol(gsva_stage2),ncol(gsva_stage1)))
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

class<-rep(c("stage2","stage1"),c(ncol(gsva_stage2),ncol(gsva_stage1)))

mianyi_huatu<-c()
for (i in 1:nrow(gsva_zong)) {
  b<-as.numeric(gsva_zong[i,])
  tonglu1<-rep(rownames(gsva_zong)[i],ncol(gsva_zong))
  zanshi<-cbind(b,tonglu1,class)
  mianyi_huatu<-rbind(mianyi_huatu,zanshi)
}
colnames(mianyi_huatu)<-c("scores","type","class")
mianyi_huatu<-as.data.frame(mianyi_huatu)
library(ggpubr)
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
  scale_fill_manual(values=c("#D55E00","#56B4E9")) +
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
ggsave("G:\\pan_result\\gpu_result_tcga\\all_stage_texh\\UCEC_T_exhau_plot_stage.pdf",aaa,width = 10,height = 10)
save(mianyi_huatu,file = "G:\\pan_result\\gpu_result_tcga\\all_stage_texh\\UCEC_T.Rdata")

#(4)####

checkpoint<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")
check_exp<-expr[checkpoint$Checkpoint.gene[which(checkpoint$Activator.A..Inhibitor.I.=="I")],]
check_exp<-expr[checkpoint$Checkpoint.gene,]

stage1<-label[which(label$label=="1"),1]
stage2<-label[which(label$label=="2"),1]
check_exp<-na.omit(check_exp)
stage1_check<-check_exp[,stage1]
stage2_check<-check_exp[,intersect(stage2,colnames(check_exp))]
stage_check<-cbind(stage1_check,stage2_check)
p_zanshi<-c()
mianyi_class<-rep(c(1,2),c(ncol(stage1_check),ncol(stage2_check)))
for (j in 1:nrow(stage_check)) {
  T_result<-kruskal.test(as.numeric(stage_check[j,])~mianyi_class)
  
  p<-T_result$p.value#p
  
  # t1<-T_result[[1]]#
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
which(p_zanshi<0.05)

checkgene_sig<-rownames(stage_check)[which(p_zanshi<0.05)]
checkgene_sig_p<-cbind(checkgene_sig,p_zanshi[which(p_zanshi<0.05)])
write.table(checkgene_sig_p,"G:\\pan_result\\gpu_result_tcga\\all_stage_checkpoint\\UCECstage_checkgene_sig_p.txt")

checkgene_sig2<-checkpoint[checkpoint$Checkpoint.gene %in%checkgene_sig,]
checkgene_sig2<-checkgene_sig2[order(checkgene_sig2$Activator.A..Inhibitor.I.), ]
stage1_check<-check_exp[checkgene_sig2$Checkpoint.gene,stage1]
stage2_check<-check_exp[checkgene_sig2$Checkpoint.gene,intersect(stage2,colnames(check_exp))]

stage1_check<-apply(stage1_check, 1, mean)
stage2_check<-apply(stage2_check, 1, mean)
library(fmsb)
leida<-rbind(stage1_check,stage2_check)
rownames(leida)<-c("stage1","stage2")
library(ggradar)
classleida<-c("stage1","stage2")
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
               grid.mid = 0.11,
               grid.max = 0.22, 
               values.radar = c(0,0.11,0.22),group.colours = c("#D55E00","#56B4E9"),axis.label.size=4)
plot
colnames(leida)

ggsave("G:\\pan_result\\gpu_result_tcga\\all_stage_checkpoint\\UCEC_CHECK_plot_stage.pdf",plot,width = 10,height = 10)
save(leida,file = "G:\\pan_result\\gpu_result_tcga\\all_stage_checkpoint\\UCEC_leida.Rdata")

















