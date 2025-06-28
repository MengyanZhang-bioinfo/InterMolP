#BRCA####

setwd("G:\\pan_result\\gpu_result_tcga\\BRCA")
gene<-read.csv("validategene.csv",row.names = 1)
gene<-gene$X0
gene<-gsub("-",".",gene)
exp<-read.csv("sigmoid.csv",row.names = 1)
exp<-exp[,gene]
# setdiff(gene,colnames(exp))
# #exp<-read.csv("297featuresweightme.csv",row.names = 1)
# exp[, "HLA.DRA"]
label<-read.csv("label.csv")
exp11<-exp[rownames(exp)%in%(label$X[which(label$label=="1")]),]
exp22<-exp[rownames(exp)%in%(label$X[which(label$label=="2")]),]
exp<-rbind(exp11,exp22)
# rownames(exp)<-exp$X
# exp<-exp[,-1]
expscale<-scale(exp)

library(factoextra)
dir.create("450k")
setwd("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k")
set.seed(123)
fviz_nbclust<-fviz_nbclust(expscale, kmeans, k.max = 15)#2
fviz_nbclust
ggsave("BRCA_fviz_nbclust.pdf",fviz_nbclust,width = 10,height = 10)
ggsave("BRCA_fviz_nbclust.pdf",fviz_nbclust,width = 10,height = 10)
save(expscale,file = "G:\\pan_result\\gpu_result_tcga\\all_exp\\BRCA_EXPSCALE.Rdata")



# expr<-read.table("tcga_brca_tpm.txt")
# expr<-t(expr)
# expr<-as.data.frame(expr)
# expr$name<-substr(rownames(expr),1,15)
# library(limma)
# expr<-avereps(expr,ID=expr$name)
# ?avereps
# expr<-as.data.frame(expr)
# rownames(expr)<-expr$name
# expr<-expr[,-ncol(expr)]
# expr<-as.data.frame(t(expr))
# write.table(expr,"brca_tpm_quchong.txt")
expr1<-read.table("brca_tpm_stage1.txt",row.names = 1,header = T)
expr2<-read.table("brca_tpm_stage2.txt",row.names = 1,header = T)
expr<-cbind(expr1,expr2)
# expr1<-expr[,colnames(expr)%in%exp11$X]
# expr2<-expr[,colnames(expr)%in%exp22$X]
str(expr1)
# expr1<-apply(expr1, 2, as.numeric)
# rownames(expr1)<-rownames(expr2)
# expr2<-apply(expr2, 2, as.numeric)
# rownames(expr2)<-rownames(expr1)


library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(ggsignif)


exp_tide<-t(apply(expr,1,function(x){x-(mean(x))}))
write.table(data.frame(ID=rownames(exp_tide),exp_tide,check.names=F),'BRCA_tide.txt',sep = "\t",quote=F,row.names = FALSE)


getwd()



new_folder <- "G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\stage1"
dir.create(new_folder)

source_folder <- "G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\stage1"
r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/BRCA/LM4_mean.txt"
file.copy(source_file, new_folder)
write.table(expr1,"stage1\\BRCA_tpm_stage1.txt",sep = "\t")



new_folder <- "G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\stage2"
dir.create(new_folder)

source_folder <- "G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\stage2"
r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/BRCA/LM4_mean.txt"
file.copy(source_file, new_folder)
write.table(expr2,"stage2\\BRCA_tpm_stage2.txt",sep = "\t")




setwd("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\stage1")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","BRCA_tpm_stage1.txt",perm = 100,QN = T,"stage1_CIBERSORT-Results_LM4_mean.txt")

setwd("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\stage2")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","BRCA_tpm_stage2.txt",perm = 100,QN = T,"stage2_CIBERSORT-Results_LM4_mean.txt")


stage1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\stage1\\stage1_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")
stage2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\stage2\\stage2_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")


ciber_zong<-cbind(t(stage1_ciber[,1:(ncol(stage1_ciber)-3)]),(t(stage2_ciber[,1:(ncol(stage1_ciber)-3)])))
colnames(ciber_zong)<-ciber_zong[1,]
ciber_zong<-ciber_zong[-1,]
ciber_zong<-as.data.frame(ciber_zong)
ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],2,as.numeric)
str(ciber_zong)
p_zanshi<-c()
for (j in 1:nrow(ciber_zong)) {
  T_result<-wilcox.test(as.numeric(ciber_zong[j,1:nrow(stage1_ciber)]),as.numeric(ciber_zong[j,(nrow(stage1_ciber)+1):ncol(ciber_zong)]),paired = F,var.equal = T)
  
  p<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_zanshi
p_stage_yuanshi<-p_zanshi
#[1] 0.23361657        NaN 0.06673382 0.04491467 0.45374285 0.69362902 0.86320080 0.72045094 0.09366834
p_stage_yuanshi
# 
library(ggplot2)
library(ggthemes)
class<-rep(c("zaoqi","jinzhan"),c(nrow(stage1_ciber),nrow(stage2_ciber)))
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
library(dplyr)
library(ggpubr)
aaa<-ggplot(mianyi_huatu,aes(y=scores,x=factor(type),fill=class))+
  geom_boxplot(outlier.size=1)+
  theme_classic()+
  xlab("")+ylab("")+guides()+
  scale_fill_manual(values=c("#eed65d", "#5d75ee"))+stat_compare_means(method = "wilcox.test",label = "p.signif")
# geom_text(aes(x =1, y=1, label = label[1], group = type), size = 4) +
# geom_text(aes(x =2, y=1, label = label[2], group = type), size = 4)+
# geom_text(aes(x =3, y=1, label = label[3], group = type), size = 4)+
# geom_text(aes(x =4, y=1, label = label[4], group = type), size = 4)+
# geom_text(aes(x =5, y=1, label = label[5], group = type), size = 4)+
# geom_text(aes(x =6, y=1, label = label[6], group = type), size = 4)+
# geom_text(aes(x =7, y=1, label = label[7], group = type), size = 4)
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\BRCA_stage_scrna_plot.pdf",aaa,width = 10,height = 10)



# cluster1 
new_folder <- "G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\cluster1"
dir.create(new_folder)

r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/BRCA/LM4_mean.txt"
file.copy(source_file, new_folder)


# cluster2 
new_folder <- "G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\cluster2"
dir.create(new_folder)

r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/BRCA/LM4_mean.txt"
file.copy(source_file, new_folder)

# 
# 
# # cluster3 
# new_folder <- "G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\cluster3"
# dir.create(new_folder)
# 
# r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)
# 
# file.copy(r_files, new_folder)
# 
# source_file <- "G:/pan_result/gpu_result_tcga/BRCA/scrna/LM4_mean.txt"
# file.copy(source_file, new_folder)
# 

exp2<-cbind(expr1,expr2)
#1.2  ####
#rm(exp2)

a<-c("euclidean","maximum","manhattan","canberra","binary","minkowski")
b<-c("ward.D","ward.D2","single","complete","average","mcquitty","median", "centroid")
zuhe<-c()
for (i in 1:6) {
  zuhe1<-combn(c(a[i],b),2)[,1:8]
  zuhe<-cbind(zuhe,zuhe1)
}
?hclust
tide<-read.csv("G:\\pan_result\\tide\\BRCA_tide_result.csv",header = T)

library(survival)
library(openxlsx)
pan_pfs<-read.xlsx("E:\\450k\\TCGAPFS.xlsx")
brca_os<-read.table("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\UCSC+BRCA_OS.txt",sep = "\t",header = T)


p_shengcun_pfs<-c()
p_shengcun_os<-c()
p_tide<-c()

i=1
for (i in 1:ncol(zuhe)) {
  hc<-hclust(dist(expscale,method = zuhe[1,i]),method = zuhe[2,i])
  clusters<-cutree(hc,k=2)
  h_class<-data.frame(sample=rownames(expscale),class=clusters)

  cluster1<-h_class[which(h_class$class=="1"),1]
  cluster2<-h_class[which(h_class$class=="2"),1]
  #cluster3<-h_class[which(h_class$class=="3"),1]
  
  cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
  cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
  #cluster3_exp<-exp2[,colnames(exp2) %in%cluster3]
  
  cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
  cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
  #cluster3_exp<-data.frame(gene=rownames(exp2),cluster3_exp)
  
  #PFS
  clusters<-as.data.frame(clusters)

  clusters_tide<-data.frame(Patient=rownames(clusters),group=clusters$clusters)

  samples<-substr(gsub('[.]','-',rownames(clusters)),1,15)
  clusters2<-data.frame(sample=samples,clusters)
  #pfs
  brca_samp<-gsub('[.]','-',substr(rownames(clusters),1,12))
  brca_pfs<-pan_pfs[pan_pfs$bcr_patient_barcode%in%brca_samp,]
  clusters<-data.frame(bcr_patient_barcode=brca_samp,clusters)
  brca_pfs<-merge(brca_pfs,clusters,by="bcr_patient_barcode")
  survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)
  p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  p_shengcun_pfs<-c(p_shengcun_pfs,p1)
  
  #os
  brca_os2<-merge(brca_os,clusters2,by="sample")
  str(brca_os2)
  fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  p2 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  # p1
  p_shengcun_os<-c(p_shengcun_os,p2)
  
  #tide

  tide_final<-left_join(tide,clusters_tide,by="Patient")
  T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
  p_tide_zanshi<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_tide<-c(p_tide,p_tide_zanshi)
  
  
  print(i)
}
hcluster_result<-cbind(p_shengcun_pfs,p_shengcun_os,p_tide)


save(hcluster_result,file = "hcluster_result.Rdata")
# save(p_shengcun_os,file = "405k_p_shengcun_os_hclust.Rdata")
# save(p_shengcun_pfs,file = "405k_p_shengcun_pfs_hclust.Rdata")




which(p_ciber$V1<0.05&p_shengcun_pfs<0.05)
#[1]  1  2  9 10 17 18 20 25 26 28 30 41 42
#p_ciber_hclust<-p_ciber


library(ConsensusClusterPlus)
STAD_tumor = sweep(t(exp), 1, apply(t(exp), 1, median, na.rm = T))
a<-c("hc","pam","km")
b<-c("pearson","spearman","euclidean","binary","maximum","canberra","minkowski")
zuhe<-c()
for (i in 1:3) {
  zuhe1<-combn(c(a[i],b),2)[,1:7]
  zuhe<-cbind(zuhe,zuhe1)
}


library(survival)
library(openxlsx)
pan_pfs<-read.xlsx("E:\\450k\\TCGAPFS.xlsx")
brca_os<-read.table("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\UCSC+BRCA_OS.txt",sep = "\t",header = T)
if (is.null(STAD_tumor) || nrow(STAD_tumor) == 0) {
  stop("STAD_tumor is null or empty.")
}

if (is.null(zuhe) || ncol(zuhe) < 2 || i > ncol(zuhe)) {
  stop("Check zuhe matrix and index i.")
}
p_shengcun_pfs<-c()
p_shengcun_os<-c()
p_tide<-c()

#p_ciber_hclust<-p_ciber
i=2
for (i in 1:ncol(zuhe)) {
  results = ConsensusClusterPlus(STAD_tumor, maxK = 2, reps = 1000, pItem = 0.8, pFeature = 1, title = "untitled_consensus_cluster", clusterAlg = zuhe[1,i],
                                 distance = zuhe[2,i], seed = 123, tmyPal=NULL, writeTable=FALSE, plot = "pdf")
  clusters<-results[[2]][["consensusClass"]]
  h_class<-data.frame(sample=rownames(expscale),class=clusters)

  cluster1<-h_class[which(h_class$class=="1"),1]
  cluster2<-h_class[which(h_class$class=="2"),1]
  #cluster3<-h_class[which(h_class$class=="3"),1]
  
  cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
  cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
  #cluster3_exp<-exp2[,colnames(exp2) %in%cluster3]
  
  cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
  cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
  #cluster3_exp<-data.frame(gene=rownames(exp2),cluster3_exp)
  
  #PFS
  clusters<-as.data.frame(clusters)

  clusters_tide<-data.frame(Patient=rownames(clusters),group=clusters$clusters)

  samples<-substr(gsub('[.]','-',rownames(clusters)),1,15)
  clusters2<-data.frame(sample=samples,clusters)
  #pfs
  brca_samp<-gsub('[.]','-',substr(rownames(clusters),1,12))
  brca_pfs<-pan_pfs[pan_pfs$bcr_patient_barcode%in%brca_samp,]
  clusters<-data.frame(bcr_patient_barcode=brca_samp,clusters)
  brca_pfs<-merge(brca_pfs,clusters,by="bcr_patient_barcode")
  survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)
  p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  p_shengcun_pfs<-c(p_shengcun_pfs,p1)
  
  #os
  brca_os2<-merge(brca_os,clusters2,by="sample")
  str(brca_os2)
  fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  p2 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  # p1
  p_shengcun_os<-c(p_shengcun_os,p2)
  
  #tide

  tide_final<-left_join(tide,clusters_tide,by="Patient")
  T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
  p_tide_zanshi<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_tide<-c(p_tide,p_tide_zanshi)
  
  print(i)
  
  
}
p_ciber<-as.data.frame(p_ciber)

save(p_ciber,file="45kp_ciber_consencus.Rdata")
save(p_shengcun_os,file="45kp_p_shengcun_os_consencus.Rdata")
save(p_shengcun_pfs,file="45kp_p_shengcun_pfs_consencus.Rdata")


load(file = "cluster2\\45kp_ciber_consencus.Rdata")
zuhe[,15]#[1] "km"      "pearson"
zuhe[,16]
zuhe[,17]
zuhe[,18]
zuhe[,19]
zuhe[,20]
zuhe[,21]
 
i=1
results = ConsensusClusterPlus(STAD_tumor, maxK = 3, reps = 1000, pItem = 0.8, pFeature = 1, title = "untitled_consensus_cluster", clusterAlg = zuhe[1,i],
                               distance = zuhe[2,i], seed = 123, tmyPal=NULL, writeTable=FALSE, plot = "pdf")
clusters<-results[[2]][["consensusClass"]]
h_class<-data.frame(sample=rownames(expscale),class=clusters)
 
cluster1<-h_class[which(h_class$class=="1"),1]
cluster2<-h_class[which(h_class$class=="2"),1]
cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
write.table(cluster1_exp,"G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\cluster1\\h_cluster1_exp.txt",sep="\t",row.names = F)
write.table(cluster2_exp,"G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\cluster2\\h_cluster2_exp.txt",sep="\t",row.names = F)
setwd("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\cluster1")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster1_exp.txt",perm = 100,QN = T,"cluster1_CIBERSORT-Results_LM4_mean.txt")
setwd("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\cluster2")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster2_exp.txt",perm = 100,QN = T,"cluster2_CIBERSORT-Results_LM4_mean.txt")

cluster1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\cluster1\\cluster1_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")
cluster2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\cluster2\\cluster2_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")

ciber_zong<-cbind(t(cluster1_ciber[,1:(ncol(cluster1_ciber)-3)]),(t(cluster2_ciber[,1:(ncol(cluster1_ciber)-3)])))
colnames(ciber_zong)<-ciber_zong[1,]
ciber_zong<-ciber_zong[-1,]
ciber_zong<-as.data.frame(ciber_zong)
ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],1,as.numeric)
str(ciber_zong)
p_zanshi<-c()
for (j in 1:nrow(ciber_zong)) {
  T_result<-wilcox.test(as.numeric(ciber_zong[j,1:nrow(cluster1_ciber)]),as.numeric(ciber_zong[j,(nrow(cluster1_ciber)+1):ncol(ciber_zong)]),paired = F,var.equal = T)
  
  p<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}


 
i=1
hc<-hclust(dist(expscale,method = zuhe[1,i]),method = zuhe[2,i])
clusters<-cutree(hc,k=2)

h_class<-data.frame(sample=rownames(expscale),class=clusters)
table(h_class$class)
 
 
cluster1<-h_class[which(h_class$class=="1"),1]
cluster2<-h_class[which(h_class$class=="2"),1]
#cluster3<-h_class[which(h_class$class=="3"),1]

cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
#cluster3_exp<-exp2[,colnames(exp2) %in%cluster3]

cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
#cluster3_exp<-data.frame(gene=rownames(exp2),cluster3_exp)

write.table(cluster1_exp,"G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\cluster1\\h_cluster1_exp.txt",sep="\t",row.names = F)
write.table(cluster2_exp,"G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\cluster2\\h_cluster2_exp.txt",sep="\t",row.names = F)
#write.table(cluster3_exp,"G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\cluster3\\h_cluster3_exp.txt",sep="\t",row.names = F)

setwd("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\cluster1")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster1_exp.txt",perm = 100,QN = T,"cluster1_CIBERSORT-Results_LM4_mean.txt")
setwd("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\cluster2")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster2_exp.txt",perm = 100,QN = T,"cluster2_CIBERSORT-Results_LM4_mean.txt")
# setwd("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\cluster3")
# source("cibersort.R")
# CIBERSORT("LM4_mean.txt","h_cluster3_exp.txt",perm = 100,QN = T,"cluster3_CIBERSORT-Results_LM4_mean.txt")


cluster1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\cluster1\\cluster1_CIBERSORT-Results_LM4_mean.txt",header=T)
cluster2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\cluster2\\cluster2_CIBERSORT-Results_LM4_mean.txt",header=T)
#cluster3_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\cluster3\\cluster3_CIBERSORT-Results_LM4_mean.txt",header=T)


ciber_zong<-cbind(t(cluster1_ciber[,1:(ncol(cluster1_ciber)-3)]),t(cluster2_ciber[,1:(ncol(cluster1_ciber)-3)]))
colnames(ciber_zong)<-ciber_zong[1,]
ciber_zong<-ciber_zong[-1,]
ciber_zong<-as.data.frame(ciber_zong)
ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],1,as.numeric)
str(ciber_zong)
p_zanshi<-c()

ciber_class<-rep(c(1,2),c(nrow(cluster1_ciber),nrow(cluster2_ciber)))
for (j in 1:nrow(ciber_zong)) {
  T_result<-wilcox.test(as.numeric(ciber_zong[j,])~ciber_class)
  
  p<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_cluster_new<-p_zanshi
p_cluster_new



#COAD####


setwd("G:\\pan_result\\gpu_result_tcga\\COAD")
gene<-read.csv("424featuresweightme.csv",row.names = 1)
gene<-colnames(gene)
#exp<-read.csv("sigmoid.csv")
exp<-read.csv("424featuresweightme.csv",row.names = 1)

label<-read.csv("label.csv")
exp11<-exp[rownames(exp)%in%(label$X[which(label$label=="1")]),]
exp22<-exp[rownames(exp)%in%(label$X[which(label$label=="2")]),]
exp<-rbind(exp11,exp22)
# rownames(exp)<-exp$X
# exp<-exp[,-1]
expscale<-scale(exp)

library(factoextra)
setwd("G:\\pan_result\\gpu_result_tcga\\COAD\\450k")
set.seed(123)
fviz_nbclust<-fviz_nbclust(expscale, kmeans, k.max = 15)#2
fviz_nbclust
ggsave("COAD_fviz_nbclust.pdf",fviz_nbclust,width = 10,height = 10)
 
# ggsave("COAD_kmeans_plot.pdf",kmeans_plot,width = 10,height = 10)


# expr<-read.table("tcga_brca_tpm.txt")
# expr<-t(expr)
# expr<-as.data.frame(expr)
# expr$name<-substr(rownames(expr),1,15)
# library(limma)
# expr<-avereps(expr,ID=expr$name)
# ?avereps
# expr<-as.data.frame(expr)
# rownames(expr)<-expr$name
# expr<-expr[,-ncol(expr)]
# expr<-as.data.frame(t(expr))
# write.table(expr,"brca_tpm_quchong.txt")
expr1<-read.csv("G:\\dnameth\\COAD\\exp1.csv",row.names = 1)
expr2<-read.csv("G:\\dnameth\\COAD\\exp2.csv",row.names = 1)
expr<-cbind(expr1,expr2)
# expr1<-expr[,colnames(expr)%in%exp11$X]
# expr2<-expr[,colnames(expr)%in%exp22$X]
str(expr1)
# expr1<-apply(expr1, 2, as.numeric)
# rownames(expr1)<-rownames(expr2)
# expr2<-apply(expr2, 2, as.numeric)
# rownames(expr2)<-rownames(expr1)


library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(ggsignif)


exp_tide<-t(apply(expr,1,function(x){x-(mean(x))}))
write.table(data.frame(ID=rownames(exp_tide),exp_tide,check.names=F),'coad_tide.txt',sep = "\t",quote=F,row.names = FALSE)


getwd()



new_folder <- "G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\stage1"
dir.create(new_folder)

source_folder <- "G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\stage1"
r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/COAD/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)
write.table(expr1,"stage1\\COAD_tpm_stage1.txt",sep = "\t")



new_folder <- "G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\stage2"
dir.create(new_folder)

source_folder <- "G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\stage2"
r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/COAD/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)
write.table(expr2,"stage2\\COAD_tpm_stage2.txt",sep = "\t")




setwd("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\stage1")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","COAD_tpm_stage1.txt",perm = 100,QN = T,"stage1_CIBERSORT-Results_LM4_mean.txt")

setwd("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\stage2")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","COAD_tpm_stage2.txt",perm = 100,QN = T,"stage2_CIBERSORT-Results_LM4_mean.txt")


stage1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\stage1\\stage1_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")
stage2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\stage2\\stage2_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")


ciber_zong<-cbind(t(stage1_ciber[,1:(ncol(stage1_ciber)-3)]),(t(stage2_ciber[,1:(ncol(stage1_ciber)-3)])))
colnames(ciber_zong)<-ciber_zong[1,]
ciber_zong<-ciber_zong[-1,]
ciber_zong<-as.data.frame(ciber_zong)
ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],2,as.numeric)
str(ciber_zong)
p_zanshi<-c()
for (j in 1:nrow(ciber_zong)) {
  T_result<-wilcox.test(as.numeric(ciber_zong[j,1:nrow(stage1_ciber)]),as.numeric(ciber_zong[j,(nrow(stage1_ciber)+1):ncol(ciber_zong)]),paired = F,var.equal = T)
  
  p<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_zanshi
p_stage_yuanshi<-p_zanshi
#[1] 0.23361657        NaN 0.06673382 0.04491467 0.45374285 0.69362902 0.86320080 0.72045094 0.09366834

# 
library(ggplot2)
library(ggthemes)
class<-rep(c("zaoqi","jinzhan"),c(nrow(stage1_ciber),nrow(stage2_ciber)))
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
library(dplyr)
library(ggpubr)
aaa<-ggplot(mianyi_huatu,aes(y=scores,x=factor(type),fill=class))+
  geom_boxplot(outlier.size=1)+
  theme_classic()+
  xlab("")+ylab("")+guides()+
  scale_fill_manual(values=c("#eed65d", "#5d75ee"))+stat_compare_means(method = "wilcox.test",label = "p.signif")
# geom_text(aes(x =1, y=1, label = label[1], group = type), size = 4) +
# geom_text(aes(x =2, y=1, label = label[2], group = type), size = 4)+
# geom_text(aes(x =3, y=1, label = label[3], group = type), size = 4)+
# geom_text(aes(x =4, y=1, label = label[4], group = type), size = 4)+
# geom_text(aes(x =5, y=1, label = label[5], group = type), size = 4)+
# geom_text(aes(x =6, y=1, label = label[6], group = type), size = 4)+
# geom_text(aes(x =7, y=1, label = label[7], group = type), size = 4)
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\COAD_stage_scrna_plot.pdf",aaa,width = 10,height = 10)



# cluster1 
new_folder <- "G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\cluster1"
dir.create(new_folder)

r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/COAD/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)


# cluster2 
new_folder <- "G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\cluster2"
dir.create(new_folder)

r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/COAD/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)

# 
# 
# # cluster3 
# new_folder <- "G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\cluster3"
# dir.create(new_folder)
# 
# r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)
# 
# file.copy(r_files, new_folder)
# 
# source_file <- "G:/pan_result/gpu_result_tcga/COAD/scrna/LM4_mean.txt"
# file.copy(source_file, new_folder)
# 

exp2<-cbind(expr1,expr2)
#1.2  ####
#rm(exp2)

a<-c("euclidean","maximum","manhattan","canberra","binary","minkowski")
b<-c("ward.D","ward.D2","single","complete","average","mcquitty","median", "centroid")
zuhe<-c()
for (i in 1:6) {
  zuhe1<-combn(c(a[i],b),2)[,1:8]
  zuhe<-cbind(zuhe,zuhe1)
}
?hclust
tide<-read.csv("G:\\pan_result\\tide\\coad_tide_result.csv",header = T)

library(survival)
library(openxlsx)
pan_pfs<-read.xlsx("E:\\450k\\.xlsx")
brca_os<-read.table("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\UCSC+COAD_OS.txt",sep = "\t",header = T)


p_shengcun_pfs<-c()
p_shengcun_os<-c()
p_tide<-c()

i=1
for (i in 1:ncol(zuhe)) {
  hc<-hclust(dist(expscale,method = zuhe[1,i]),method = zuhe[2,i])
  clusters<-cutree(hc,k=2)
  h_class<-data.frame(sample=rownames(expscale),class=clusters)

  cluster1<-h_class[which(h_class$class=="1"),1]
  cluster2<-h_class[which(h_class$class=="2"),1]
  #cluster3<-h_class[which(h_class$class=="3"),1]
  
  cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
  cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
  #cluster3_exp<-exp2[,colnames(exp2) %in%cluster3]
  
  cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
  cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
  #cluster3_exp<-data.frame(gene=rownames(exp2),cluster3_exp)
  
  #PFS
  clusters<-as.data.frame(clusters)

  clusters_tide<-data.frame(Patient=rownames(clusters),group=clusters$clusters)

  samples<-substr(gsub('[.]','-',rownames(clusters)),1,15)
  clusters2<-data.frame(sample=samples,clusters)
  #pfs
  brca_samp<-gsub('[.]','-',substr(rownames(clusters),1,12))
  brca_pfs<-pan_pfs[pan_pfs$bcr_patient_barcode%in%brca_samp,]
  clusters<-data.frame(bcr_patient_barcode=brca_samp,clusters)
  brca_pfs<-merge(brca_pfs,clusters,by="bcr_patient_barcode")
  survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)
  p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  p_shengcun_pfs<-c(p_shengcun_pfs,p1)
  
  #os
  brca_os2<-merge(brca_os,clusters2,by="sample")
  str(brca_os2)
  fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  p2 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  # p1
  p_shengcun_os<-c(p_shengcun_os,p2)
  
  #tide

  tide_final<-left_join(tide,clusters_tide,by="Patient")
  T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
  p_tide_zanshi<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_tide<-c(p_tide,p_tide_zanshi)
  
  
  print(i)
}
hcluster_result<-cbind(p_shengcun_pfs,p_shengcun_os,p_tide)


save(hcluster_result,file = "hcluster_result.Rdata")
# save(p_shengcun_os,file = "405k_p_shengcun_os_hclust.Rdata")
# save(p_shengcun_pfs,file = "405k_p_shengcun_pfs_hclust.Rdata")




which(p_ciber$V1<0.05&p_shengcun_pfs<0.05)
#[1]  1  2  9 10 17 18 20 25 26 28 30 41 42
#p_ciber_hclust<-p_ciber


library(ConsensusClusterPlus)
STAD_tumor = sweep(t(exp), 1, apply(t(exp), 1, median, na.rm = T))
a<-c("hc","pam","km")
b<-c("pearson","spearman","euclidean","binary","maximum","canberra","minkowski")
zuhe<-c()
for (i in 1:3) {
  zuhe1<-combn(c(a[i],b),2)[,1:7]
  zuhe<-cbind(zuhe,zuhe1)
}
library(survival)
library(openxlsx)
pan_pfs<-read.xlsx("E:\\450k\\TCGAPFS.xlsx")
brca_os<-read.table("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\UCSC+COAD_OS.txt",sep = "\t",header = T)


p_ciber<-c()
p_shengcun_pfs<-c()
p_shengcun_os<-c()


#p_ciber_hclust<-p_ciber
i=1
for (i in 1:ncol(zuhe)) {
  results = ConsensusClusterPlus(STAD_tumor, maxK = 3, reps = 1000, pItem = 0.8, pFeature = 1, title = "untitled_consensus_cluster", clusterAlg = zuhe[1,i],
                                 distance = zuhe[2,i], seed = 123, tmyPal=NULL, writeTable=FALSE, plot = "pdf")
  clusters<-results[[3]][["consensusClass"]]
  if(min(table(clusters))>15){
    h_class<-data.frame(sample=rownames(expscale),class=clusters)
  
    cluster1<-h_class[which(h_class$class=="1"),1]
    cluster2<-h_class[which(h_class$class=="2"),1]
    cluster3<-h_class[which(h_class$class=="3"),1]
    
    cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
    cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
    cluster3_exp<-exp2[,colnames(exp2) %in%cluster3]
    
    cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
    cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
    cluster3_exp<-data.frame(gene=rownames(exp2),cluster3_exp)
    
    write.table(cluster1_exp,"G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\cluster1\\h_cluster1_exp.txt",sep="\t",row.names = F)
    write.table(cluster2_exp,"G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\cluster2\\h_cluster2_exp.txt",sep="\t",row.names = F)
    write.table(cluster3_exp,"G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\cluster3\\h_cluster3_exp.txt",sep="\t",row.names = F)
    
    setwd("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\cluster1")
    source("cibersort.R")
    CIBERSORT("LM4_mean.txt","h_cluster1_exp.txt",perm = 100,QN = T,"cluster1_CIBERSORT-Results_LM4_mean.txt")
    setwd("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\cluster2")
    source("cibersort.R")
    CIBERSORT("LM4_mean.txt","h_cluster2_exp.txt",perm = 100,QN = T,"cluster2_CIBERSORT-Results_LM4_mean.txt")
    setwd("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\cluster3")
    source("cibersort.R")
    CIBERSORT("LM4_mean.txt","h_cluster3_exp.txt",perm = 100,QN = T,"cluster3_CIBERSORT-Results_LM4_mean.txt")
    
    
    cluster1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\cluster1\\cluster1_CIBERSORT-Results_LM4_mean.txt",header=T)
    cluster2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\cluster2\\cluster2_CIBERSORT-Results_LM4_mean.txt",header=T)
    cluster3_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\cluster3\\cluster3_CIBERSORT-Results_LM4_mean.txt",header=T)
    
    
    ciber_zong<-cbind(t(cluster1_ciber[,1:(ncol(cluster1_ciber)-3)]),t(cluster2_ciber[,1:(ncol(cluster1_ciber)-3)]),t(cluster3_ciber[,1:(ncol(cluster1_ciber)-3)]))
    colnames(ciber_zong)<-ciber_zong[1,]
    ciber_zong<-ciber_zong[-1,]
    ciber_zong<-as.data.frame(ciber_zong)
    ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],1,as.numeric)
    str(ciber_zong)
    p_zanshi<-c()
    
    ciber_class<-rep(c(1,2,3),c(nrow(cluster1_ciber),nrow(cluster2_ciber),nrow(cluster3_ciber)))
    for (j in 1:nrow(ciber_zong)) {
      T_result<-kruskal.test(as.numeric(ciber_zong[j,])~ciber_class)
      
      p<-T_result$p.value
      
      # t1<-T_result[[1]]
      p_zanshi<-c(p_zanshi,p)
      #t_zanshi<-c(t_zanshi,t1)
    }
  }else{
    p_zanshi<-rep(1,nrow(ciber_zong))
  }
  p_ciber<-rbind(p_ciber,p_zanshi)
  #PFS
  clusters<-as.data.frame(clusters)

  samples<-substr(gsub('[.]','-',rownames(clusters)),1,15)
  clusters2<-data.frame(sample=samples,clusters)
  #pfs
  brca_samp<-gsub('[.]','-',substr(rownames(clusters),1,12))
  brca_pfs<-pan_pfs[pan_pfs$bcr_patient_barcode%in%brca_samp,]
  clusters<-data.frame(bcr_patient_barcode=brca_samp,clusters)
  brca_pfs<-merge(brca_pfs,clusters,by="bcr_patient_barcode")
  survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)
  p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  p_shengcun_pfs<-c(p_shengcun_pfs,p1)
  
  #os
  brca_os2<-merge(brca_os,clusters2,by="sample")
  str(brca_os2)
  fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  p2 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  # p1
  p_shengcun_os<-c(p_shengcun_os,p2)
  
  
  print(i)
  
  
}
p_ciber<-as.data.frame(p_ciber)

save(p_ciber,file="45kp_ciber_consencus.Rdata")
save(p_shengcun_os,file="45kp_p_shengcun_os_consencus.Rdata")
save(p_shengcun_pfs,file="45kp_p_shengcun_pfs_consencus.Rdata")


load(file = "cluster2\\45kp_ciber_consencus.Rdata")
zuhe[,15]#[1] "km"      "pearson"
zuhe[,16]
zuhe[,17]
zuhe[,18]
zuhe[,19]
zuhe[,20]
zuhe[,21]

i=1
results = ConsensusClusterPlus(STAD_tumor, maxK = 3, reps = 1000, pItem = 0.8, pFeature = 1, title = "untitled_consensus_cluster", clusterAlg = zuhe[1,i],
                               distance = zuhe[2,i], seed = 123, tmyPal=NULL, writeTable=FALSE, plot = "pdf")
clusters<-results[[2]][["consensusClass"]]
h_class<-data.frame(sample=rownames(expscale),class=clusters)
 
cluster1<-h_class[which(h_class$class=="1"),1]
cluster2<-h_class[which(h_class$class=="2"),1]
cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
write.table(cluster1_exp,"G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\cluster1\\h_cluster1_exp.txt",sep="\t",row.names = F)
write.table(cluster2_exp,"G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\cluster2\\h_cluster2_exp.txt",sep="\t",row.names = F)
setwd("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\cluster1")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster1_exp.txt",perm = 100,QN = T,"cluster1_CIBERSORT-Results_LM4_mean.txt")
setwd("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\cluster2")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster2_exp.txt",perm = 100,QN = T,"cluster2_CIBERSORT-Results_LM4_mean.txt")

cluster1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\cluster1\\cluster1_CIBERSORT-Results_LM4_mean.txt",header=T)
cluster2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\cluster2\\cluster2_CIBERSORT-Results_LM4_mean.txt",header=T)

ciber_zong<-cbind(t(cluster1_ciber[,1:(ncol(cluster1_ciber)-3)]),(t(cluster2_ciber[,1:(ncol(cluster1_ciber)-3)])))
colnames(ciber_zong)<-ciber_zong[1,]
ciber_zong<-ciber_zong[-1,]
ciber_zong<-as.data.frame(ciber_zong)
ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],1,as.numeric)
str(ciber_zong)
p_zanshi<-c()
for (j in 1:nrow(ciber_zong)) {
  T_result<-wilcox.test(as.numeric(ciber_zong[j,1:nrow(cluster1_ciber)]),as.numeric(ciber_zong[j,(nrow(cluster1_ciber)+1):ncol(ciber_zong)]),paired = F,var.equal = T)
  
  p<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}


 
i=26
hc<-hclust(dist(expscale,method = zuhe[1,i]),method = zuhe[2,i])
clusters<-cutree(hc,k=2)

h_class<-data.frame(sample=rownames(expscale),class=clusters)
 
 
cluster1<-h_class[which(h_class$class=="1"),1]
cluster2<-h_class[which(h_class$class=="2"),1]
#cluster3<-h_class[which(h_class$class=="3"),1]

cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
#cluster3_exp<-exp2[,colnames(exp2) %in%cluster3]

cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
#cluster3_exp<-data.frame(gene=rownames(exp2),cluster3_exp)

write.table(cluster1_exp,"G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\cluster1\\h_cluster1_exp.txt",sep="\t",row.names = F)
write.table(cluster2_exp,"G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\cluster2\\h_cluster2_exp.txt",sep="\t",row.names = F)
#write.table(cluster3_exp,"G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\cluster3\\h_cluster3_exp.txt",sep="\t",row.names = F)

setwd("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\cluster1")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster1_exp.txt",perm = 100,QN = T,"cluster1_CIBERSORT-Results_LM4_mean.txt")
setwd("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\cluster2")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster2_exp.txt",perm = 100,QN = T,"cluster2_CIBERSORT-Results_LM4_mean.txt")
# setwd("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\cluster3")
# source("cibersort.R")
# CIBERSORT("LM4_mean.txt","h_cluster3_exp.txt",perm = 100,QN = T,"cluster3_CIBERSORT-Results_LM4_mean.txt")


cluster1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\cluster1\\cluster1_CIBERSORT-Results_LM4_mean.txt",header=T)
cluster2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\cluster2\\cluster2_CIBERSORT-Results_LM4_mean.txt",header=T)
#cluster3_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\cluster3\\cluster3_CIBERSORT-Results_LM4_mean.txt",header=T)


ciber_zong<-cbind(t(cluster1_ciber[,1:(ncol(cluster1_ciber)-3)]),t(cluster2_ciber[,1:(ncol(cluster1_ciber)-3)]))
colnames(ciber_zong)<-ciber_zong[1,]
ciber_zong<-ciber_zong[-1,]
ciber_zong<-as.data.frame(ciber_zong)
ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],1,as.numeric)
str(ciber_zong)
p_zanshi<-c()

ciber_class<-rep(c(1,2),c(nrow(cluster1_ciber),nrow(cluster2_ciber)))
for (j in 1:nrow(ciber_zong)) {
  T_result<-wilcox.test(as.numeric(ciber_zong[j,])~ciber_class)
  
  p<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_cluster_new<-p_zanshi

#ESCA####

setwd("G:\\pan_result\\gpu_result_tcga\\ESCA")
gene<-read.csv("202featuresweightme.csv",row.names = 1)
gene<-colnames(gene)
#exp<-read.csv("sigmoid.csv")
exp<-read.csv("202featuresweightme.csv",row.names = 1)

label<-read.csv("label.csv")
exp11<-exp[rownames(exp)%in%(label$X[which(label$label=="1")]),]
exp22<-exp[rownames(exp)%in%(label$X[which(label$label=="2")]),]
exp<-rbind(exp11,exp22)
# rownames(exp)<-exp$X
# exp<-exp[,-1]
expscale<-scale(exp)

library(factoextra)
setwd("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k")
set.seed(123)
fviz_nbclust<-fviz_nbclust(expscale, kmeans, k.max = 15)#2
fviz_nbclust
ggsave("ESCA_fviz_nbclust.pdf",fviz_nbclust,width = 10,height = 10)
 
# ggsave("ESCA_kmeans_plot.pdf",kmeans_plot,width = 10,height = 10)


# expr<-read.table("tcga_brca_tpm.txt")
# expr<-t(expr)
# expr<-as.data.frame(expr)
# expr$name<-substr(rownames(expr),1,15)
# library(limma)
# expr<-avereps(expr,ID=expr$name)
# ?avereps
# expr<-as.data.frame(expr)
# rownames(expr)<-expr$name
# expr<-expr[,-ncol(expr)]
# expr<-as.data.frame(t(expr))
# write.table(expr,"brca_tpm_quchong.txt")
expr1<-read.csv("G:\\dnameth\\ESCA\\exp1.csv",row.names = 1)
expr2<-read.csv("G:\\dnameth\\ESCA\\exp2.csv",row.names = 1)
expr<-cbind(expr1,expr2)
# expr1<-expr[,colnames(expr)%in%exp11$X]
# expr2<-expr[,colnames(expr)%in%exp22$X]
str(expr1)
# expr1<-apply(expr1, 2, as.numeric)
# rownames(expr1)<-rownames(expr2)
# expr2<-apply(expr2, 2, as.numeric)
# rownames(expr2)<-rownames(expr1)


library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(ggsignif)


exp_tide<-t(apply(expr,1,function(x){x-(mean(x))}))
write.table(data.frame(ID=rownames(exp_tide),exp_tide,check.names=F),'ESCA_tide.txt',sep = "\t",quote=F,row.names = FALSE)


getwd()



new_folder <- "G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\stage1"
dir.create(new_folder)

source_folder <- "G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\stage1"
r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/ESCA/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)
write.table(expr1,"stage1\\ESCA_tpm_stage1.txt",sep = "\t")



new_folder <- "G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\stage2"
dir.create(new_folder)

source_folder <- "G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\stage2"
r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/ESCA/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)
write.table(expr2,"stage2\\ESCA_tpm_stage2.txt",sep = "\t")




setwd("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\stage1")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","ESCA_tpm_stage1.txt",perm = 100,QN = T,"stage1_CIBERSORT-Results_LM4_mean.txt")

setwd("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\stage2")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","ESCA_tpm_stage2.txt",perm = 100,QN = T,"stage2_CIBERSORT-Results_LM4_mean.txt")


stage1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\stage1\\stage1_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")
stage2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\stage2\\stage2_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")


ciber_zong<-cbind(t(stage1_ciber[,1:(ncol(stage1_ciber)-3)]),(t(stage2_ciber[,1:(ncol(stage1_ciber)-3)])))
colnames(ciber_zong)<-ciber_zong[1,]
ciber_zong<-ciber_zong[-1,]
ciber_zong<-as.data.frame(ciber_zong)
ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],2,as.numeric)
str(ciber_zong)
p_zanshi<-c()
for (j in 1:nrow(ciber_zong)) {
  T_result<-wilcox.test(as.numeric(ciber_zong[j,1:nrow(stage1_ciber)]),as.numeric(ciber_zong[j,(nrow(stage1_ciber)+1):ncol(ciber_zong)]),paired = F,var.equal = T)
  
  p<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_zanshi
p_stage_yuanshi<-p_zanshi
#[1] 0.23361657        NaN 0.06673382 0.04491467 0.45374285 0.69362902 0.86320080 0.72045094 0.09366834

# 
library(ggplot2)
library(ggthemes)
class<-rep(c("zaoqi","jinzhan"),c(nrow(stage1_ciber),nrow(stage2_ciber)))
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
library(dplyr)
library(ggpubr)
aaa<-ggplot(mianyi_huatu,aes(y=scores,x=factor(type),fill=class))+
  geom_boxplot(outlier.size=1)+
  theme_classic()+
  xlab("")+ylab("")+guides()+
  scale_fill_manual(values=c("#eed65d", "#5d75ee"))+stat_compare_means(method = "wilcox.test",label = "p.signif")
# geom_text(aes(x =1, y=1, label = label[1], group = type), size = 4) +
# geom_text(aes(x =2, y=1, label = label[2], group = type), size = 4)+
# geom_text(aes(x =3, y=1, label = label[3], group = type), size = 4)+
# geom_text(aes(x =4, y=1, label = label[4], group = type), size = 4)+
# geom_text(aes(x =5, y=1, label = label[5], group = type), size = 4)+
# geom_text(aes(x =6, y=1, label = label[6], group = type), size = 4)+
# geom_text(aes(x =7, y=1, label = label[7], group = type), size = 4)
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\ESCA_stage_scrna_plot.pdf",aaa,width = 10,height = 10)



# cluster1 
new_folder <- "G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\cluster1"
dir.create(new_folder)

r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/ESCA/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)


# cluster2 
new_folder <- "G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\cluster2"
dir.create(new_folder)

r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/ESCA/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)

# 
# 
# # cluster3 
# new_folder <- "G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\cluster3"
# dir.create(new_folder)
# 
# r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)
# 
# file.copy(r_files, new_folder)
# 
# source_file <- "G:/pan_result/gpu_result_tcga/ESCA/scrna/LM4_mean.txt"
# file.copy(source_file, new_folder)
# 

exp2<-cbind(expr1,expr2)
#1.2  ####
#rm(exp2)

a<-c("euclidean","maximum","manhattan","canberra","binary","minkowski")
b<-c("ward.D","ward.D2","single","complete","average","mcquitty","median", "centroid")
zuhe<-c()
for (i in 1:6) {
  zuhe1<-combn(c(a[i],b),2)[,1:8]
  zuhe<-cbind(zuhe,zuhe1)
}
?hclust
tide<-read.csv("G:\\pan_result\\tide\\ESCA_tide_result.csv",header = T)

library(survival)
library(openxlsx)
pan_pfs<-read.xlsx("E:\\450k\\TCGAPFS.xlsx")
brca_os<-read.table("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\UCSC+ESCA_OS.txt",sep = "\t",header = T)


p_shengcun_pfs<-c()
p_shengcun_os<-c()
p_tide<-c()

i=1
for (i in 1:ncol(zuhe)) {
  hc<-hclust(dist(expscale,method = zuhe[1,i]),method = zuhe[2,i])
  clusters<-cutree(hc,k=2)
  h_class<-data.frame(sample=rownames(expscale),class=clusters)

  cluster1<-h_class[which(h_class$class=="1"),1]
  cluster2<-h_class[which(h_class$class=="2"),1]
  #cluster3<-h_class[which(h_class$class=="3"),1]
  
  cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
  cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
  #cluster3_exp<-exp2[,colnames(exp2) %in%cluster3]
  
  cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
  cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
  #cluster3_exp<-data.frame(gene=rownames(exp2),cluster3_exp)
  
  #PFS
  clusters<-as.data.frame(clusters)

  clusters_tide<-data.frame(Patient=rownames(clusters),group=clusters$clusters)

  samples<-substr(gsub('[.]','-',rownames(clusters)),1,15)
  clusters2<-data.frame(sample=samples,clusters)
  #pfs
  brca_samp<-gsub('[.]','-',substr(rownames(clusters),1,12))
  brca_pfs<-pan_pfs[pan_pfs$bcr_patient_barcode%in%brca_samp,]
  clusters<-data.frame(bcr_patient_barcode=brca_samp,clusters)
  brca_pfs<-merge(brca_pfs,clusters,by="bcr_patient_barcode")
  survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)
  p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  p_shengcun_pfs<-c(p_shengcun_pfs,p1)
  
  #os
  brca_os2<-merge(brca_os,clusters2,by="sample")
  str(brca_os2)
  fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  p2 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  # p1
  p_shengcun_os<-c(p_shengcun_os,p2)
  
  #tide

  tide_final<-left_join(tide,clusters_tide,by="Patient")
  T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
  p_tide_zanshi<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_tide<-c(p_tide,p_tide_zanshi)
  
  
  print(i)
}
hcluster_result<-cbind(p_shengcun_pfs,p_shengcun_os,p_tide)


save(hcluster_result,file = "hcluster_result.Rdata")
# save(p_shengcun_os,file = "405k_p_shengcun_os_hclust.Rdata")
# save(p_shengcun_pfs,file = "405k_p_shengcun_pfs_hclust.Rdata")




which(p_ciber$V1<0.05&p_shengcun_pfs<0.05)
#[1]  1  2  9 10 17 18 20 25 26 28 30 41 42
#p_ciber_hclust<-p_ciber


library(ConsensusClusterPlus)
STAD_tumor = sweep(t(exp), 1, apply(t(exp), 1, median, na.rm = T))
a<-c("hc","pam","km")
b<-c("pearson","spearman","euclidean","binary","maximum","canberra","minkowski")
zuhe<-c()
for (i in 1:3) {
  zuhe1<-combn(c(a[i],b),2)[,1:7]
  zuhe<-cbind(zuhe,zuhe1)
}


library(survival)
library(openxlsx)
pan_pfs<-read.xlsx("E:\\450k\\TCGAPFS.xlsx")
brca_os<-read.table("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\UCSC+ESCA_OS.txt",sep = "\t",header = T)
if (is.null(STAD_tumor) || nrow(STAD_tumor) == 0) {
  stop("STAD_tumor is null or empty.")
}

if (is.null(zuhe) || ncol(zuhe) < 2 || i > ncol(zuhe)) {
  stop("Check zuhe matrix and index i.")
}
p_shengcun_pfs<-c()
p_shengcun_os<-c()
p_tide<-c()

#p_ciber_hclust<-p_ciber
i=2
for (i in 1:ncol(zuhe)) {
  results = ConsensusClusterPlus(STAD_tumor, maxK = 5, reps = 1000, pItem = 0.8, pFeature = 1, title = "untitled_consensus_cluster", clusterAlg = zuhe[1,i],
                                 distance = zuhe[2,i], seed = 123, tmyPal=NULL, writeTable=FALSE, plot = "pdf")
  clusters<-results[[2]][["consensusClass"]]
  h_class<-data.frame(sample=rownames(expscale),class=clusters)

  cluster1<-h_class[which(h_class$class=="1"),1]
  cluster2<-h_class[which(h_class$class=="2"),1]
  #cluster3<-h_class[which(h_class$class=="3"),1]
  
  cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
  cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
  #cluster3_exp<-exp2[,colnames(exp2) %in%cluster3]
  
  cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
  cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
  #cluster3_exp<-data.frame(gene=rownames(exp2),cluster3_exp)
  
  #PFS
  clusters<-as.data.frame(clusters)

  clusters_tide<-data.frame(Patient=rownames(clusters),group=clusters$clusters)

  samples<-substr(gsub('[.]','-',rownames(clusters)),1,15)
  clusters2<-data.frame(sample=samples,clusters)
  #pfs
  brca_samp<-gsub('[.]','-',substr(rownames(clusters),1,12))
  brca_pfs<-pan_pfs[pan_pfs$bcr_patient_barcode%in%brca_samp,]
  clusters<-data.frame(bcr_patient_barcode=brca_samp,clusters)
  brca_pfs<-merge(brca_pfs,clusters,by="bcr_patient_barcode")
  survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)
  p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  p_shengcun_pfs<-c(p_shengcun_pfs,p1)
  
  #os
  brca_os2<-merge(brca_os,clusters2,by="sample")
  str(brca_os2)
  fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  p2 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  # p1
  p_shengcun_os<-c(p_shengcun_os,p2)
  
  #tide

  tide_final<-left_join(tide,clusters_tide,by="Patient")
  T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
  p_tide_zanshi<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_tide<-c(p_tide,p_tide_zanshi)
  
  print(i)
  
  
}
concluster_result<-cbind(p_shengcun_pfs,p_shengcun_os,p_tide)


save(concluster_result,file = "concluster_result.Rdata")

p_ciber<-as.data.frame(p_ciber)

save(p_ciber,file="45kp_ciber_consencus.Rdata")
save(p_shengcun_os,file="45kp_p_shengcun_os_consencus.Rdata")
save(p_shengcun_pfs,file="45kp_p_shengcun_pfs_consencus.Rdata")


load(file = "cluster2\\45kp_ciber_consencus.Rdata")
zuhe[,15]#[1] "km"      "pearson"
zuhe[,16]
zuhe[,17]
zuhe[,18]
zuhe[,19]
zuhe[,20]
zuhe[,21]
 
i=1
results = ConsensusClusterPlus(STAD_tumor, maxK = 3, reps = 1000, pItem = 0.8, pFeature = 1, title = "untitled_consensus_cluster", clusterAlg = zuhe[1,i],
                               distance = zuhe[2,i], seed = 123, tmyPal=NULL, writeTable=FALSE, plot = "pdf")
clusters<-results[[2]][["consensusClass"]]
h_class<-data.frame(sample=rownames(expscale),class=clusters)

cluster1<-h_class[which(h_class$class=="1"),1]
cluster2<-h_class[which(h_class$class=="2"),1]
cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
write.table(cluster1_exp,"G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\cluster1\\h_cluster1_exp.txt",sep="\t",row.names = F)
write.table(cluster2_exp,"G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\cluster2\\h_cluster2_exp.txt",sep="\t",row.names = F)
setwd("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\cluster1")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster1_exp.txt",perm = 100,QN = T,"cluster1_CIBERSORT-Results_LM4_mean.txt")
setwd("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\cluster2")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster2_exp.txt",perm = 100,QN = T,"cluster2_CIBERSORT-Results_LM4_mean.txt")

cluster1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\cluster1\\cluster1_CIBERSORT-Results_LM4_mean.txt",header=T)
cluster2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\cluster2\\cluster2_CIBERSORT-Results_LM4_mean.txt",header=T)

ciber_zong<-cbind(t(cluster1_ciber[,1:(ncol(cluster1_ciber)-3)]),(t(cluster2_ciber[,1:(ncol(cluster1_ciber)-3)])))
colnames(ciber_zong)<-ciber_zong[1,]
ciber_zong<-ciber_zong[-1,]
ciber_zong<-as.data.frame(ciber_zong)
ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],1,as.numeric)
str(ciber_zong)
p_zanshi<-c()
for (j in 1:nrow(ciber_zong)) {
  T_result<-wilcox.test(as.numeric(ciber_zong[j,1:nrow(cluster1_ciber)]),as.numeric(ciber_zong[j,(nrow(cluster1_ciber)+1):ncol(ciber_zong)]),paired = F,var.equal = T)
  
  p<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}



i=9
hc<-hclust(dist(expscale,method = zuhe[1,i]),method = zuhe[2,i])
clusters<-cutree(hc,k=2)
table(clusters)
h_class<-data.frame(sample=rownames(expscale),class=clusters)
 
 
cluster1<-h_class[which(h_class$class=="1"),1]
cluster2<-h_class[which(h_class$class=="2"),1]
#cluster3<-h_class[which(h_class$class=="3"),1]

cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
#cluster3_exp<-exp2[,colnames(exp2) %in%cluster3]

cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
#cluster3_exp<-data.frame(gene=rownames(exp2),cluster3_exp)

write.table(cluster1_exp,"G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\cluster1\\h_cluster1_exp.txt",sep="\t",row.names = F)
write.table(cluster2_exp,"G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\cluster2\\h_cluster2_exp.txt",sep="\t",row.names = F)
#write.table(cluster3_exp,"G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\cluster3\\h_cluster3_exp.txt",sep="\t",row.names = F)

setwd("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\cluster1")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster1_exp.txt",perm = 100,QN = T,"cluster1_CIBERSORT-Results_LM4_mean.txt")
setwd("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\cluster2")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster2_exp.txt",perm = 100,QN = T,"cluster2_CIBERSORT-Results_LM4_mean.txt")
# setwd("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\cluster3")
# source("cibersort.R")
# CIBERSORT("LM4_mean.txt","h_cluster3_exp.txt",perm = 100,QN = T,"cluster3_CIBERSORT-Results_LM4_mean.txt")


cluster1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\cluster1\\cluster1_CIBERSORT-Results_LM4_mean.txt",header=T)
cluster2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\cluster2\\cluster2_CIBERSORT-Results_LM4_mean.txt",header=T)
#cluster3_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\cluster3\\cluster3_CIBERSORT-Results_LM4_mean.txt",header=T)


ciber_zong<-cbind(t(cluster1_ciber[,1:(ncol(cluster1_ciber)-3)]),t(cluster2_ciber[,1:(ncol(cluster1_ciber)-3)]))
colnames(ciber_zong)<-ciber_zong[1,]
ciber_zong<-ciber_zong[-1,]
ciber_zong<-as.data.frame(ciber_zong)
ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],1,as.numeric)
str(ciber_zong)
p_zanshi<-c()

ciber_class<-rep(c(1,2),c(nrow(cluster1_ciber),nrow(cluster2_ciber)))
for (j in 1:nrow(ciber_zong)) {
  T_result<-wilcox.test(as.numeric(ciber_zong[j,])~ciber_class)
  
  p<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_cluster_new<-p_zanshi

#HNSC####

setwd("G:\\pan_result\\gpu_result_tcga\\HNSC")
gene<-read.csv("297featuresweightme.csv",row.names = 1)
gene<-colnames(gene)
#exp<-read.csv("sigmoid.csv")
exp<-read.csv("297featuresweightme.csv",row.names = 1)

label<-read.csv("label.csv")
exp11<-exp[rownames(exp)%in%(label$X[which(label$label=="1")]),]
exp22<-exp[rownames(exp)%in%(label$X[which(label$label=="2")]),]
exp<-rbind(exp11,exp22)
# rownames(exp)<-exp$X
# exp<-exp[,-1]
expscale<-scale(exp)

library(factoextra)
dir.create("450k")
setwd("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k")
set.seed(123)
fviz_nbclust<-fviz_nbclust(expscale, kmeans, k.max = 15)#2
fviz_nbclust
ggsave("HNSC_fviz_nbclust.pdf",fviz_nbclust,width = 10,height = 10)
 
# ggsave("HNSC_kmeans_plot.pdf",kmeans_plot,width = 10,height = 10)


# expr<-read.table("tcga_brca_tpm.txt")
# expr<-t(expr)
# expr<-as.data.frame(expr)
# expr$name<-substr(rownames(expr),1,15)
# library(limma)
# expr<-avereps(expr,ID=expr$name)
# ?avereps
# expr<-as.data.frame(expr)
# rownames(expr)<-expr$name
# expr<-expr[,-ncol(expr)]
# expr<-as.data.frame(t(expr))
# write.table(expr,"brca_tpm_quchong.txt")
expr1<-read.csv("G:\\dnameth\\HNSC\\exp1.csv",row.names = 1)
expr2<-read.csv("G:\\dnameth\\HNSC\\exp2.csv",row.names = 1)
expr<-cbind(expr1,expr2)
# expr1<-expr[,colnames(expr)%in%exp11$X]
# expr2<-expr[,colnames(expr)%in%exp22$X]
str(expr1)
# expr1<-apply(expr1, 2, as.numeric)
# rownames(expr1)<-rownames(expr2)
# expr2<-apply(expr2, 2, as.numeric)
# rownames(expr2)<-rownames(expr1)


library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(ggsignif)


exp_tide<-t(apply(expr,1,function(x){x-(mean(x))}))
write.table(data.frame(ID=rownames(exp_tide),exp_tide,check.names=F),'HNSC_tide.txt',sep = "\t",quote=F,row.names = FALSE)


getwd()



new_folder <- "G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\stage1"
dir.create(new_folder)

source_folder <- "G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\stage1"
r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/HNSC/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)
write.table(expr1,"stage1\\HNSC_tpm_stage1.txt",sep = "\t")



new_folder <- "G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\stage2"
dir.create(new_folder)

source_folder <- "G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\stage2"
r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/HNSC/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)
write.table(expr2,"stage2\\HNSC_tpm_stage2.txt",sep = "\t")




setwd("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\stage1")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","HNSC_tpm_stage1.txt",perm = 100,QN = T,"stage1_CIBERSORT-Results_LM4_mean.txt")

setwd("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\stage2")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","HNSC_tpm_stage2.txt",perm = 100,QN = T,"stage2_CIBERSORT-Results_LM4_mean.txt")


stage1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\stage1\\stage1_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")
stage2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\stage2\\stage2_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")


ciber_zong<-cbind(t(stage1_ciber[,1:(ncol(stage1_ciber)-3)]),(t(stage2_ciber[,1:(ncol(stage1_ciber)-3)])))
colnames(ciber_zong)<-ciber_zong[1,]
ciber_zong<-ciber_zong[-1,]
ciber_zong<-as.data.frame(ciber_zong)
ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],2,as.numeric)
str(ciber_zong)
p_zanshi<-c()
for (j in 1:nrow(ciber_zong)) {
  T_result<-wilcox.test(as.numeric(ciber_zong[j,1:nrow(stage1_ciber)]),as.numeric(ciber_zong[j,(nrow(stage1_ciber)+1):ncol(ciber_zong)]),paired = F,var.equal = T)
  
  p<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_zanshi
p_stage_yuanshi<-p_zanshi
#[1] 0.23361657        NaN 0.06673382 0.04491467 0.45374285 0.69362902 0.86320080 0.72045094 0.09366834
p_stage_yuanshi
# 
library(ggplot2)
library(ggthemes)
class<-rep(c("zaoqi","jinzhan"),c(nrow(stage1_ciber),nrow(stage2_ciber)))
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
library(dplyr)
library(ggpubr)
aaa<-ggplot(mianyi_huatu,aes(y=scores,x=factor(type),fill=class))+
  geom_boxplot(outlier.size=1)+
  theme_classic()+
  xlab("")+ylab("")+guides()+
  scale_fill_manual(values=c("#eed65d", "#5d75ee"))+stat_compare_means(method = "wilcox.test",label = "p.signif")
# geom_text(aes(x =1, y=1, label = label[1], group = type), size = 4) +
# geom_text(aes(x =2, y=1, label = label[2], group = type), size = 4)+
# geom_text(aes(x =3, y=1, label = label[3], group = type), size = 4)+
# geom_text(aes(x =4, y=1, label = label[4], group = type), size = 4)+
# geom_text(aes(x =5, y=1, label = label[5], group = type), size = 4)+
# geom_text(aes(x =6, y=1, label = label[6], group = type), size = 4)+
# geom_text(aes(x =7, y=1, label = label[7], group = type), size = 4)
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_stage_scrna_plot.pdf",aaa,width = 10,height = 10)



# cluster1 
new_folder <- "G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\cluster1"
dir.create(new_folder)

r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/HNSC/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)


# cluster2 
new_folder <- "G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\cluster2"
dir.create(new_folder)

r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/HNSC/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)

# 
# 
# # cluster3 
# new_folder <- "G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\cluster3"
# dir.create(new_folder)
# 
# r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)
# 
# file.copy(r_files, new_folder)
# 
# source_file <- "G:/pan_result/gpu_result_tcga/HNSC/scrna/LM4_mean.txt"
# file.copy(source_file, new_folder)
# 

exp2<-cbind(expr1,expr2)
#1.2  ####
#rm(exp2)

a<-c("euclidean","maximum","manhattan","canberra","binary","minkowski")
b<-c("ward.D","ward.D2","single","complete","average","mcquitty","median", "centroid")
zuhe<-c()
for (i in 1:6) {
  zuhe1<-combn(c(a[i],b),2)[,1:8]
  zuhe<-cbind(zuhe,zuhe1)
}
?hclust
tide<-read.csv("G:\\pan_result\\tide\\HNSC_tide_result.csv",header = T)

library(survival)
library(openxlsx)
pan_pfs<-read.xlsx("E:\\450k\\TCGAPFS.xlsx")
brca_os<-read.table("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\UCSC+HNSC_OS.txt",sep = "\t",header = T)


p_shengcun_pfs<-c()
p_shengcun_os<-c()
p_tide<-c()

i=1
for (i in 1:ncol(zuhe)) {
  hc<-hclust(dist(expscale,method = zuhe[1,i]),method = zuhe[2,i])
  clusters<-cutree(hc,k=2)
  h_class<-data.frame(sample=rownames(expscale),class=clusters)

  cluster1<-h_class[which(h_class$class=="1"),1]
  cluster2<-h_class[which(h_class$class=="2"),1]
  #cluster3<-h_class[which(h_class$class=="3"),1]
  
  cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
  cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
  #cluster3_exp<-exp2[,colnames(exp2) %in%cluster3]
  
  cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
  cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
  #cluster3_exp<-data.frame(gene=rownames(exp2),cluster3_exp)
  
  #PFS
  clusters<-as.data.frame(clusters)

  clusters_tide<-data.frame(Patient=rownames(clusters),group=clusters$clusters)

  samples<-substr(gsub('[.]','-',rownames(clusters)),1,15)
  clusters2<-data.frame(sample=samples,clusters)
  #pfs
  brca_samp<-gsub('[.]','-',substr(rownames(clusters),1,12))
  brca_pfs<-pan_pfs[pan_pfs$bcr_patient_barcode%in%brca_samp,]
  clusters<-data.frame(bcr_patient_barcode=brca_samp,clusters)
  brca_pfs<-merge(brca_pfs,clusters,by="bcr_patient_barcode")
  survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)
  p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  p_shengcun_pfs<-c(p_shengcun_pfs,p1)
  
  #os
  brca_os2<-merge(brca_os,clusters2,by="sample")
  str(brca_os2)
  fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  p2 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  # p1
  p_shengcun_os<-c(p_shengcun_os,p2)
  
  #tide

  tide_final<-left_join(tide,clusters_tide,by="Patient")
  T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
  p_tide_zanshi<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_tide<-c(p_tide,p_tide_zanshi)
  
  
  print(i)
}
hcluster_result<-cbind(p_shengcun_pfs,p_shengcun_os,p_tide)


save(hcluster_result,file = "hcluster_result.Rdata")
# save(p_shengcun_os,file = "405k_p_shengcun_os_hclust.Rdata")
# save(p_shengcun_pfs,file = "405k_p_shengcun_pfs_hclust.Rdata")




which(p_ciber$V1<0.05&p_shengcun_pfs<0.05)
#[1]  1  2  9 10 17 18 20 25 26 28 30 41 42
#p_ciber_hclust<-p_ciber


library(ConsensusClusterPlus)
STAD_tumor = sweep(t(exp), 1, apply(t(exp), 1, median, na.rm = T))
a<-c("hc","pam","km")
b<-c("pearson","spearman","euclidean","binary","maximum","canberra","minkowski")
zuhe<-c()
for (i in 1:3) {
  zuhe1<-combn(c(a[i],b),2)[,1:7]
  zuhe<-cbind(zuhe,zuhe1)
}


library(survival)
library(openxlsx)
pan_pfs<-read.xlsx("E:\\450k\\TCGAPFS.xlsx")
brca_os<-read.table("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\UCSC+HNSC_OS.txt",sep = "\t",header = T)
if (is.null(STAD_tumor) || nrow(STAD_tumor) == 0) {
  stop("STAD_tumor is null or empty.")
}

if (is.null(zuhe) || ncol(zuhe) < 2 || i > ncol(zuhe)) {
  stop("Check zuhe matrix and index i.")
}
p_shengcun_pfs<-c()
p_shengcun_os<-c()
p_tide<-c()

#p_ciber_hclust<-p_ciber
i=2
for (i in 1:ncol(zuhe)) {
  results = ConsensusClusterPlus(STAD_tumor, maxK = 2, reps = 1000, pItem = 0.8, pFeature = 1, title = "untitled_consensus_cluster", clusterAlg = zuhe[1,i],
                                 distance = zuhe[2,i], seed = 123, tmyPal=NULL, writeTable=FALSE, plot = "pdf")
  clusters<-results[[2]][["consensusClass"]]
  h_class<-data.frame(sample=rownames(expscale),class=clusters)

  cluster1<-h_class[which(h_class$class=="1"),1]
  cluster2<-h_class[which(h_class$class=="2"),1]
  #cluster3<-h_class[which(h_class$class=="3"),1]
  
  cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
  cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
  #cluster3_exp<-exp2[,colnames(exp2) %in%cluster3]
  
  cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
  cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
  #cluster3_exp<-data.frame(gene=rownames(exp2),cluster3_exp)
  
  #PFS
  clusters<-as.data.frame(clusters)

  clusters_tide<-data.frame(Patient=rownames(clusters),group=clusters$clusters)

  samples<-substr(gsub('[.]','-',rownames(clusters)),1,15)
  clusters2<-data.frame(sample=samples,clusters)
  #pfs
  brca_samp<-gsub('[.]','-',substr(rownames(clusters),1,12))
  brca_pfs<-pan_pfs[pan_pfs$bcr_patient_barcode%in%brca_samp,]
  clusters<-data.frame(bcr_patient_barcode=brca_samp,clusters)
  brca_pfs<-merge(brca_pfs,clusters,by="bcr_patient_barcode")
  survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)
  p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  p_shengcun_pfs<-c(p_shengcun_pfs,p1)
  
  #os
  brca_os2<-merge(brca_os,clusters2,by="sample")
  str(brca_os2)
  fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  p2 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  # p1
  p_shengcun_os<-c(p_shengcun_os,p2)
  
  #tide

  tide_final<-left_join(tide,clusters_tide,by="Patient")
  T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
  p_tide_zanshi<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_tide<-c(p_tide,p_tide_zanshi)
  
  print(i)
  
  
}
p_ciber<-as.data.frame(p_ciber)

save(p_ciber,file="45kp_ciber_consencus.Rdata")
save(p_shengcun_os,file="45kp_p_shengcun_os_consencus.Rdata")
save(p_shengcun_pfs,file="45kp_p_shengcun_pfs_consencus.Rdata")


load(file = "cluster2\\45kp_ciber_consencus.Rdata")
zuhe[,15]#[1] "km"      "pearson"
zuhe[,16]
zuhe[,17]
zuhe[,18]
zuhe[,19]
zuhe[,20]
zuhe[,21]
 
i=1
results = ConsensusClusterPlus(STAD_tumor, maxK = 3, reps = 1000, pItem = 0.8, pFeature = 1, title = "untitled_consensus_cluster", clusterAlg = zuhe[1,i],
                               distance = zuhe[2,i], seed = 123, tmyPal=NULL, writeTable=FALSE, plot = "pdf")
clusters<-results[[2]][["consensusClass"]]
h_class<-data.frame(sample=rownames(expscale),class=clusters)
 
cluster1<-h_class[which(h_class$class=="1"),1]
cluster2<-h_class[which(h_class$class=="2"),1]
cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
write.table(cluster1_exp,"G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\cluster1\\h_cluster1_exp.txt",sep="\t",row.names = F)
write.table(cluster2_exp,"G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\cluster2\\h_cluster2_exp.txt",sep="\t",row.names = F)
setwd("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\cluster1")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster1_exp.txt",perm = 100,QN = T,"cluster1_CIBERSORT-Results_LM4_mean.txt")
setwd("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\cluster2")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster2_exp.txt",perm = 100,QN = T,"cluster2_CIBERSORT-Results_LM4_mean.txt")

cluster1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\cluster1\\cluster1_CIBERSORT-Results_LM4_mean.txt",header=T)
cluster2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\cluster2\\cluster2_CIBERSORT-Results_LM4_mean.txt",header=T)

ciber_zong<-cbind(t(cluster1_ciber[,1:(ncol(cluster1_ciber)-3)]),(t(cluster2_ciber[,1:(ncol(cluster1_ciber)-3)])))
colnames(ciber_zong)<-ciber_zong[1,]
ciber_zong<-ciber_zong[-1,]
ciber_zong<-as.data.frame(ciber_zong)
ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],1,as.numeric)
str(ciber_zong)
p_zanshi<-c()
for (j in 1:nrow(ciber_zong)) {
  T_result<-wilcox.test(as.numeric(ciber_zong[j,1:nrow(cluster1_ciber)]),as.numeric(ciber_zong[j,(nrow(cluster1_ciber)+1):ncol(ciber_zong)]),paired = F,var.equal = T)
  
  p<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}


 
i=29
hc<-hclust(dist(expscale,method = zuhe[1,i]),method = zuhe[2,i])
clusters<-cutree(hc,k=2)

h_class<-data.frame(sample=rownames(expscale),class=clusters)
table(h_class$class)
 
 
cluster1<-h_class[which(h_class$class=="1"),1]
cluster2<-h_class[which(h_class$class=="2"),1]
#cluster3<-h_class[which(h_class$class=="3"),1]

cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
#cluster3_exp<-exp2[,colnames(exp2) %in%cluster3]

cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
#cluster3_exp<-data.frame(gene=rownames(exp2),cluster3_exp)

write.table(cluster1_exp,"G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\cluster1\\h_cluster1_exp.txt",sep="\t",row.names = F)
write.table(cluster2_exp,"G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\cluster2\\h_cluster2_exp.txt",sep="\t",row.names = F)
#write.table(cluster3_exp,"G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\cluster3\\h_cluster3_exp.txt",sep="\t",row.names = F)

setwd("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\cluster1")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster1_exp.txt",perm = 100,QN = T,"cluster1_CIBERSORT-Results_LM4_mean.txt")
setwd("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\cluster2")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster2_exp.txt",perm = 100,QN = T,"cluster2_CIBERSORT-Results_LM4_mean.txt")
# setwd("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\cluster3")
# source("cibersort.R")
# CIBERSORT("LM4_mean.txt","h_cluster3_exp.txt",perm = 100,QN = T,"cluster3_CIBERSORT-Results_LM4_mean.txt")


cluster1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\cluster1\\cluster1_CIBERSORT-Results_LM4_mean.txt",header=T)
cluster2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\cluster2\\cluster2_CIBERSORT-Results_LM4_mean.txt",header=T)
#cluster3_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\cluster3\\cluster3_CIBERSORT-Results_LM4_mean.txt",header=T)


ciber_zong<-cbind(t(cluster1_ciber[,1:(ncol(cluster1_ciber)-3)]),t(cluster2_ciber[,1:(ncol(cluster1_ciber)-3)]))
colnames(ciber_zong)<-ciber_zong[1,]
ciber_zong<-ciber_zong[-1,]
ciber_zong<-as.data.frame(ciber_zong)
ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],1,as.numeric)
str(ciber_zong)
p_zanshi<-c()

ciber_class<-rep(c(1,2),c(nrow(cluster1_ciber),nrow(cluster2_ciber)))
for (j in 1:nrow(ciber_zong)) {
  T_result<-wilcox.test(as.numeric(ciber_zong[j,])~ciber_class)
  
  p<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_cluster_new<-p_zanshi
p_cluster_new

#KIRC####

setwd("G:\\pan_result\\gpu_result_tcga\\KIRC")
gene<-read.csv("563featuresweightme.csv",row.names = 1)
gene<-colnames(gene)
#exp<-read.csv("sigmoid.csv")
exp<-read.csv("563featuresweightme.csv",row.names = 1)

label<-read.csv("label.csv")
exp11<-exp[rownames(exp)%in%(label$X[which(label$label=="1")]),]
exp22<-exp[rownames(exp)%in%(label$X[which(label$label=="2")]),]
exp<-rbind(exp11,exp22)
# rownames(exp)<-exp$X
# exp<-exp[,-1]
expscale<-scale(exp)

library(factoextra)
setwd("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k")
set.seed(123)
fviz_nbclust<-fviz_nbclust(expscale, kmeans, k.max = 15)#3
fviz_nbclust
ggsave("KIRC_fviz_nbclust.pdf",fviz_nbclust,width = 10,height = 10)

# ggsave("KIRC_kmeans_plot.pdf",kmeans_plot,width = 10,height = 10)


# expr<-read.table("tcga_brca_tpm.txt")
# expr<-t(expr)
# expr<-as.data.frame(expr)
# expr$name<-substr(rownames(expr),1,15)
# library(limma)
# expr<-avereps(expr,ID=expr$name)
# ?avereps
# expr<-as.data.frame(expr)
# rownames(expr)<-expr$name
# expr<-expr[,-ncol(expr)]
# expr<-as.data.frame(t(expr))
# write.table(expr,"brca_tpm_quchong.txt")
expr1<-read.csv("G:\\dnameth\\KIRC\\exp1.csv",row.names = 1)
expr2<-read.csv("G:\\dnameth\\KIRC\\exp2.csv",row.names = 1)
expr<-cbind(expr1,expr2)
# expr1<-expr[,colnames(expr)%in%exp11$X]
# expr2<-expr[,colnames(expr)%in%exp22$X]
str(expr1)
# expr1<-apply(expr1, 2, as.numeric)
# rownames(expr1)<-rownames(expr2)
# expr2<-apply(expr2, 2, as.numeric)
# rownames(expr2)<-rownames(expr1)


library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(ggsignif)


exp_tide<-t(apply(expr,1,function(x){x-(mean(x))}))
write.table(data.frame(ID=rownames(exp_tide),exp_tide,check.names=F),'KIRC_tide.txt',sep = "\t",quote=F,row.names = FALSE)


getwd()



new_folder <- "G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\stage1"
dir.create(new_folder)

source_folder <- "G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\stage1"
r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/KIRC/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)
write.table(expr1,"stage1\\KIRC_tpm_stage1.txt",sep = "\t")



new_folder <- "G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\stage2"
dir.create(new_folder)

source_folder <- "G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\stage2"
r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/KIRC/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)
write.table(expr2,"stage2\\KIRC_tpm_stage2.txt",sep = "\t")




setwd("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\stage1")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","KIRC_tpm_stage1.txt",perm = 100,QN = T,"stage1_CIBERSORT-Results_LM4_mean.txt")

setwd("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\stage2")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","KIRC_tpm_stage2.txt",perm = 100,QN = T,"stage2_CIBERSORT-Results_LM4_mean.txt")


stage1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\stage1\\stage1_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")
stage2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\stage2\\stage2_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")


ciber_zong<-cbind(t(stage1_ciber[,1:(ncol(stage1_ciber)-3)]),(t(stage2_ciber[,1:(ncol(stage1_ciber)-3)])))
colnames(ciber_zong)<-ciber_zong[1,]
ciber_zong<-ciber_zong[-1,]
ciber_zong<-as.data.frame(ciber_zong)
ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],2,as.numeric)
str(ciber_zong)
p_zanshi<-c()
for (j in 1:nrow(ciber_zong)) {
  T_result<-wilcox.test(as.numeric(ciber_zong[j,1:nrow(stage1_ciber)]),as.numeric(ciber_zong[j,(nrow(stage1_ciber)+1):ncol(ciber_zong)]),paired = F,var.equal = T)
  
  p<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_zanshi
p_stage_yuanshi<-p_zanshi
#[1] 0.23361657        NaN 0.06673382 0.04491467 0.45374285 0.69362902 0.86320080 0.72045094 0.09366834

# 
library(ggplot2)
library(ggthemes)
class<-rep(c("zaoqi","jinzhan"),c(nrow(stage1_ciber),nrow(stage2_ciber)))
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
library(dplyr)
library(ggpubr)
aaa<-ggplot(mianyi_huatu,aes(y=scores,x=factor(type),fill=class))+
  geom_boxplot(outlier.size=1)+
  theme_classic()+
  xlab("")+ylab("")+guides()+
  scale_fill_manual(values=c("#eed65d", "#5d75ee"))+stat_compare_means(method = "wilcox.test",label = "p.signif")
# geom_text(aes(x =1, y=1, label = label[1], group = type), size = 4) +
# geom_text(aes(x =2, y=1, label = label[2], group = type), size = 4)+
# geom_text(aes(x =3, y=1, label = label[3], group = type), size = 4)+
# geom_text(aes(x =4, y=1, label = label[4], group = type), size = 4)+
# geom_text(aes(x =5, y=1, label = label[5], group = type), size = 4)+
# geom_text(aes(x =6, y=1, label = label[6], group = type), size = 4)+
# geom_text(aes(x =7, y=1, label = label[7], group = type), size = 4)
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\KIRC_stage_scrna_plot.pdf",aaa,width = 10,height = 10)



# cluster1 
new_folder <- "G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\cluster1"
dir.create(new_folder)

r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/KIRC/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)


# cluster2 
new_folder <- "G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\cluster2"
dir.create(new_folder)

r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/KIRC/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)

# 

# cluster3 
new_folder <- "G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\cluster3"
dir.create(new_folder)

r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/KIRC/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)


exp2<-cbind(expr1,expr2)
#1.2  ####
#rm(exp2)

a<-c("euclidean","maximum","manhattan","canberra","binary","minkowski")
b<-c("ward.D","ward.D2","single","complete","average","mcquitty","median", "centroid")
zuhe<-c()
for (i in 1:6) {
  zuhe1<-combn(c(a[i],b),2)[,1:8]
  zuhe<-cbind(zuhe,zuhe1)
}
?hclust
tide<-read.csv("G:\\pan_result\\tide\\KIRC_tide_result.csv",header = T)

library(survival)
library(openxlsx)
pan_pfs<-read.xlsx("E:\\450k\\TCGAPFS.xlsx")
brca_os<-read.table("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\UCSC+KIRC_OS.txt",sep = "\t",header = T)


p_shengcun_pfs<-c()
p_shengcun_os<-c()
p_tide<-c()

i=1
for (i in 1:ncol(zuhe)) {
  hc<-hclust(dist(expscale,method = zuhe[1,i]),method = zuhe[2,i])
  clusters<-cutree(hc,k=3)
  h_class<-data.frame(sample=rownames(expscale),class=clusters)

  cluster1<-h_class[which(h_class$class=="1"),1]
  cluster2<-h_class[which(h_class$class=="2"),1]
  cluster3<-h_class[which(h_class$class=="3"),1]
  
  cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
  cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
  cluster3_exp<-exp2[,colnames(exp2) %in%cluster3]
  
  cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
  cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
  cluster3_exp<-data.frame(gene=rownames(exp2),cluster3_exp)
  
  #PFS
  clusters<-as.data.frame(clusters)

  clusters_tide<-data.frame(Patient=rownames(clusters),group=clusters$clusters)

  samples<-substr(gsub('[.]','-',rownames(clusters)),1,15)
  clusters2<-data.frame(sample=samples,clusters)
  #pfs
  brca_samp<-gsub('[.]','-',substr(rownames(clusters),1,12))
  brca_pfs<-pan_pfs[pan_pfs$bcr_patient_barcode%in%brca_samp,]
  clusters<-data.frame(bcr_patient_barcode=brca_samp,clusters)
  brca_pfs<-merge(brca_pfs,clusters,by="bcr_patient_barcode")
  survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)
  p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  p_shengcun_pfs<-c(p_shengcun_pfs,p1)
  
  #os
  brca_os2<-merge(brca_os,clusters2,by="sample")
  str(brca_os2)
  fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  p2 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  # p1
  p_shengcun_os<-c(p_shengcun_os,p2)
  
  #tide

  tide_final<-left_join(tide,clusters_tide,by="Patient")
  T_result<-kruskal.test(as.numeric(tide_final$TIDE)~tide_final$group)
  p_tide_zanshi<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_tide<-c(p_tide,p_tide_zanshi)
  
  
  print(i)
}
hcluster_result<-cbind(p_shengcun_pfs,p_shengcun_os,p_tide)


save(hcluster_result,file = "hcluster_result.Rdata")
# save(p_shengcun_os,file = "405k_p_shengcun_os_hclust.Rdata")
# save(p_shengcun_pfs,file = "405k_p_shengcun_pfs_hclust.Rdata")




which(p_ciber$V1<0.05&p_shengcun_pfs<0.05)
#[1]  1  2  9 10 17 18 20 25 26 28 30 41 42
#p_ciber_hclust<-p_ciber


library(ConsensusClusterPlus)
STAD_tumor = sweep(t(exp), 1, apply(t(exp), 1, median, na.rm = T))
a<-c("hc","pam","km")
b<-c("pearson","spearman","euclidean","binary","maximum","canberra","minkowski")
zuhe<-c()
for (i in 1:3) {
  zuhe1<-combn(c(a[i],b),2)[,1:7]
  zuhe<-cbind(zuhe,zuhe1)
}


library(survival)
library(openxlsx)
pan_pfs<-read.xlsx("E:\\450k\\TCGAPFS.xlsx")
brca_os<-read.table("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\UCSC+KIRC_OS.txt",sep = "\t",header = T)
if (is.null(STAD_tumor) || nrow(STAD_tumor) == 0) {
  stop("STAD_tumor is null or empty.")
}

if (is.null(zuhe) || ncol(zuhe) < 2 || i > ncol(zuhe)) {
  stop("Check zuhe matrix and index i.")
}
p_shengcun_pfs<-c()
p_shengcun_os<-c()
p_tide<-c()

#p_ciber_hclust<-p_ciber
i=2
for (i in 1:ncol(zuhe)) {
  results = ConsensusClusterPlus(STAD_tumor, maxK = 2, reps = 1000, pItem = 0.8, pFeature = 1, title = "untitled_consensus_cluster", clusterAlg = zuhe[1,i],
                                 distance = zuhe[2,i], seed = 123, tmyPal=NULL, writeTable=FALSE, plot = "pdf")
  clusters<-results[[2]][["consensusClass"]]
  h_class<-data.frame(sample=rownames(expscale),class=clusters)

  cluster1<-h_class[which(h_class$class=="1"),1]
  cluster2<-h_class[which(h_class$class=="2"),1]
  #cluster3<-h_class[which(h_class$class=="3"),1]
  
  cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
  cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
  #cluster3_exp<-exp2[,colnames(exp2) %in%cluster3]
  
  cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
  cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
  #cluster3_exp<-data.frame(gene=rownames(exp2),cluster3_exp)
  
  #PFS
  clusters<-as.data.frame(clusters)

  clusters_tide<-data.frame(Patient=rownames(clusters),group=clusters$clusters)

  samples<-substr(gsub('[.]','-',rownames(clusters)),1,15)
  clusters2<-data.frame(sample=samples,clusters)
  #pfs
  brca_samp<-gsub('[.]','-',substr(rownames(clusters),1,12))
  brca_pfs<-pan_pfs[pan_pfs$bcr_patient_barcode%in%brca_samp,]
  clusters<-data.frame(bcr_patient_barcode=brca_samp,clusters)
  brca_pfs<-merge(brca_pfs,clusters,by="bcr_patient_barcode")
  survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)
  p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  p_shengcun_pfs<-c(p_shengcun_pfs,p1)
  
  #os
  brca_os2<-merge(brca_os,clusters2,by="sample")
  str(brca_os2)
  fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  p2 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  # p1
  p_shengcun_os<-c(p_shengcun_os,p2)
  
  #tide

  tide_final<-left_join(tide,clusters_tide,by="Patient")
  T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
  p_tide_zanshi<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_tide<-c(p_tide,p_tide_zanshi)
  
  print(i)
  
  
}
p_ciber<-as.data.frame(p_ciber)

save(p_ciber,file="45kp_ciber_consencus.Rdata")
save(p_shengcun_os,file="45kp_p_shengcun_os_consencus.Rdata")
save(p_shengcun_pfs,file="45kp_p_shengcun_pfs_consencus.Rdata")


load(file = "cluster2\\45kp_ciber_consencus.Rdata")
zuhe[,15]#[1] "km"      "pearson"
zuhe[,16]
zuhe[,17]
zuhe[,18]
zuhe[,19]
zuhe[,20]
zuhe[,21]
 
i=1
results = ConsensusClusterPlus(STAD_tumor, maxK = 3, reps = 1000, pItem = 0.8, pFeature = 1, title = "untitled_consensus_cluster", clusterAlg = zuhe[1,i],
                               distance = zuhe[2,i], seed = 123, tmyPal=NULL, writeTable=FALSE, plot = "pdf")
clusters<-results[[2]][["consensusClass"]]
h_class<-data.frame(sample=rownames(expscale),class=clusters)
 
cluster1<-h_class[which(h_class$class=="1"),1]
cluster2<-h_class[which(h_class$class=="2"),1]
cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
write.table(cluster1_exp,"G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\cluster1\\h_cluster1_exp.txt",sep="\t",row.names = F)
write.table(cluster2_exp,"G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\cluster2\\h_cluster2_exp.txt",sep="\t",row.names = F)
setwd("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\cluster1")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster1_exp.txt",perm = 100,QN = T,"cluster1_CIBERSORT-Results_LM4_mean.txt")
setwd("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\cluster2")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster2_exp.txt",perm = 100,QN = T,"cluster2_CIBERSORT-Results_LM4_mean.txt")

cluster1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\cluster1\\cluster1_CIBERSORT-Results_LM4_mean.txt",header=T)
cluster2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\cluster2\\cluster2_CIBERSORT-Results_LM4_mean.txt",header=T)

ciber_zong<-cbind(t(cluster1_ciber[,1:(ncol(cluster1_ciber)-3)]),(t(cluster2_ciber[,1:(ncol(cluster1_ciber)-3)])))
colnames(ciber_zong)<-ciber_zong[1,]
ciber_zong<-ciber_zong[-1,]
ciber_zong<-as.data.frame(ciber_zong)
ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],1,as.numeric)
str(ciber_zong)
p_zanshi<-c()
for (j in 1:nrow(ciber_zong)) {
  T_result<-wilcox.test(as.numeric(ciber_zong[j,1:nrow(cluster1_ciber)]),as.numeric(ciber_zong[j,(nrow(cluster1_ciber)+1):ncol(ciber_zong)]),paired = F,var.equal = T)
  
  p<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}


 
i=26
hc<-hclust(dist(expscale,method = zuhe[1,i]),method = zuhe[2,i])
clusters<-cutree(hc,k=3)

h_class<-data.frame(sample=rownames(expscale),class=clusters)
table(h_class$class)
 
 
cluster1<-h_class[which(h_class$class=="1"),1]
cluster2<-h_class[which(h_class$class=="2"),1]
cluster3<-h_class[which(h_class$class=="3"),1]
table(h_class$class)
cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
cluster3_exp<-exp2[,colnames(exp2) %in%cluster3]

cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
cluster3_exp<-data.frame(gene=rownames(exp2),cluster3_exp)

write.table(cluster1_exp,"G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\cluster1\\h_cluster1_exp.txt",sep="\t",row.names = F)
write.table(cluster2_exp,"G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\cluster2\\h_cluster2_exp.txt",sep="\t",row.names = F)
write.table(cluster3_exp,"G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\cluster3\\h_cluster3_exp.txt",sep="\t",row.names = F)


setwd("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\cluster1")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster1_exp.txt",perm = 100,QN = T,"cluster1_CIBERSORT-Results_LM4_mean.txt")
setwd("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\cluster2")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster2_exp.txt",perm = 100,QN = T,"cluster2_CIBERSORT-Results_LM4_mean.txt")
setwd("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\cluster3")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster3_exp.txt",perm = 100,QN = T,"cluster3_CIBERSORT-Results_LM4_mean.txt")


cluster1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\cluster1\\cluster1_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")
cluster2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\cluster2\\cluster2_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")
cluster3_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\cluster3\\cluster3_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")


ciber_zong<-cbind(t(cluster1_ciber[,1:(ncol(cluster1_ciber)-3)]),t(cluster2_ciber[,1:(ncol(cluster1_ciber)-3)]),t(cluster3_ciber[,1:(ncol(cluster1_ciber)-3)]))
colnames(ciber_zong)<-ciber_zong[1,]
ciber_zong<-ciber_zong[-1,]
ciber_zong<-as.data.frame(ciber_zong)
ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],1,as.numeric)
str(ciber_zong)
p_zanshi<-c()

ciber_class<-rep(c(1,2,3),c(nrow(cluster1_ciber),nrow(cluster2_ciber),nrow(cluster3_ciber)))
for (j in 1:nrow(ciber_zong)) {
  T_result<-kruskal.test(as.numeric(ciber_zong[j,])~ciber_class)
  
  p<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_cluster_new<-p_zanshi
p_cluster_new

#KIRP####

setwd("G:\\pan_result\\gpu_result_tcga\\KIRP")
gene<-read.csv("444featuresweightme.csv",row.names = 1)
gene<-colnames(gene)
#exp<-read.csv("sigmoid.csv")
exp<-read.csv("444featuresweightme.csv",row.names = 1)

label<-read.csv("label.csv")
exp11<-exp[rownames(exp)%in%(label$X[which(label$label=="1")]),]
exp22<-exp[rownames(exp)%in%(label$X[which(label$label=="2")]),]
exp<-rbind(exp11,exp22)
# rownames(exp)<-exp$X
# exp<-exp[,-1]
expscale<-scale(exp)

library(factoextra)
setwd("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k")
set.seed(123)
fviz_nbclust<-fviz_nbclust(expscale, kmeans, k.max = 15)#2
fviz_nbclust
ggsave("KIRP_fviz_nbclust.pdf",fviz_nbclust,width = 10,height = 10)
 
# ggsave("KIRP_kmeans_plot.pdf",kmeans_plot,width = 10,height = 10)


# expr<-read.table("tcga_brca_tpm.txt")
# expr<-t(expr)
# expr<-as.data.frame(expr)
# expr$name<-substr(rownames(expr),1,15)
# library(limma)
# expr<-avereps(expr,ID=expr$name)
# ?avereps
# expr<-as.data.frame(expr)
# rownames(expr)<-expr$name
# expr<-expr[,-ncol(expr)]
# expr<-as.data.frame(t(expr))
# write.table(expr,"brca_tpm_quchong.txt")
expr1<-read.csv("G:\\dnameth\\KIRP\\exp1.csv",row.names = 1)
expr2<-read.csv("G:\\dnameth\\KIRP\\exp2.csv",row.names = 1)
expr<-cbind(expr1,expr2)
# expr1<-expr[,colnames(expr)%in%exp11$X]
# expr2<-expr[,colnames(expr)%in%exp22$X]
str(expr1)
# expr1<-apply(expr1, 2, as.numeric)
# rownames(expr1)<-rownames(expr2)
# expr2<-apply(expr2, 2, as.numeric)
# rownames(expr2)<-rownames(expr1)


library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(ggsignif)


exp_tide<-t(apply(expr,1,function(x){x-(mean(x))}))
write.table(data.frame(ID=rownames(exp_tide),exp_tide,check.names=F),'KIRP_tide.txt',sep = "\t",quote=F,row.names = FALSE)


getwd()


library(utils)
# rforge <- "http://r-forge.r-project.org"#
# if(!"estimate" %in% installed.packages())

library(estimate)

dat = log2(edgeR::cpm(expr)+1)
write.table(dat, file = "KIRP_estimate_input.txt", sep = '\t', quote = F)
exp.file = "KIRP_estimate_input.txt"
in.gct.file = "ESTIMATE_KIRP_input.gct"
 
outputGCT(exp.file, in.gct.file)
rt <- read.table("ESTIMATE_KIRP_input.gct", 
                 skip = 2,  
                 header = TRUE, 
                 sep = "\t")

 
filterCommonGenes(input.f = exp.file, 
                  output.f = "ESTIMATE_KIRP_gene.gct",  
                  id = "GeneSymbol"  
) 

#2esitimate 

out.score.file = "ESTIMATE_score.gct"
 
estimateScore(in.gct.file, 
              out.score.file, 
              platform = "illumina")  

 
ESTIMATE_score = read.table(out.score.file,
                            skip = 2,
                            header = T,
                            row.names = 1)
ESTIMATE_score = as.data.frame(t(ESTIMATE_score[, 2:ncol(ESTIMATE_score)]))  
ESTIMATE_score[1:6, 1:3]
write.table(ESTIMATE_score, file = "KIRP_estimate_score.txt", sep = '\t', quote = F)




new_folder <- "G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\stage1"
dir.create(new_folder)

source_folder <- "G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\stage1"
r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/KIRC/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)
write.table(expr1,"stage1\\KIRP_tpm_stage1.txt",sep = "\t")



new_folder <- "G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\stage2"
dir.create(new_folder)

source_folder <- "G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\stage2"
r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/KIRC/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)
write.table(expr2,"stage2\\KIRP_tpm_stage2.txt",sep = "\t")




setwd("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\stage1")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","KIRP_tpm_stage1.txt",perm = 100,QN = T,"stage1_CIBERSORT-Results_LM4_mean.txt")

setwd("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\stage2")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","KIRP_tpm_stage2.txt",perm = 100,QN = T,"stage2_CIBERSORT-Results_LM4_mean.txt")


stage1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\stage1\\stage1_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")
stage2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\stage2\\stage2_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")


ciber_zong<-cbind(t(stage1_ciber[,1:(ncol(stage1_ciber)-3)]),(t(stage2_ciber[,1:(ncol(stage1_ciber)-3)])))
colnames(ciber_zong)<-ciber_zong[1,]
ciber_zong<-ciber_zong[-1,]
ciber_zong<-as.data.frame(ciber_zong)
ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],2,as.numeric)
str(ciber_zong)
p_zanshi<-c()
for (j in 1:nrow(ciber_zong)) {
  T_result<-wilcox.test(as.numeric(ciber_zong[j,1:nrow(stage1_ciber)]),as.numeric(ciber_zong[j,(nrow(stage1_ciber)+1):ncol(ciber_zong)]),paired = F,var.equal = T)
  
  p<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_zanshi
p_stage_yuanshi<-p_zanshi
#[1] 0.23361657        NaN 0.06673382 0.04491467 0.45374285 0.69362902 0.86320080 0.72045094 0.09366834

#
library(ggplot2)
library(ggthemes)
class<-rep(c("zaoqi","jinzhan"),c(nrow(stage1_ciber),nrow(stage2_ciber)))
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
library(dplyr)
library(ggpubr)
aaa<-ggplot(mianyi_huatu,aes(y=scores,x=factor(type),fill=class))+
  geom_boxplot(outlier.size=1)+
  theme_classic()+
  xlab("")+ylab("")+guides()+
  scale_fill_manual(values=c("#eed65d", "#5d75ee"))+stat_compare_means(method = "wilcox.test",label = "p.signif")
# geom_text(aes(x =1, y=1, label = label[1], group = type), size = 4) +
# geom_text(aes(x =2, y=1, label = label[2], group = type), size = 4)+
# geom_text(aes(x =3, y=1, label = label[3], group = type), size = 4)+
# geom_text(aes(x =4, y=1, label = label[4], group = type), size = 4)+
# geom_text(aes(x =5, y=1, label = label[5], group = type), size = 4)+
# geom_text(aes(x =6, y=1, label = label[6], group = type), size = 4)+
# geom_text(aes(x =7, y=1, label = label[7], group = type), size = 4)
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\KIRP_stage_scrna_plot.pdf",aaa,width = 10,height = 10)



# cluster1 
new_folder <- "G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\cluster1"
dir.create(new_folder)

r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/KIRC/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)


# cluster2 
new_folder <- "G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\cluster2"
dir.create(new_folder)

r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/KIRC/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)

# 
# 
# # cluster3 
# new_folder <- "G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\cluster3"
# dir.create(new_folder)
# 
# r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)
# 
# file.copy(r_files, new_folder)
# 
# source_file <- "G:/pan_result/gpu_result_tcga/KIRP/scrna/LM4_mean.txt"
# file.copy(source_file, new_folder)
# 

exp2<-cbind(expr1,expr2)
#1.2  ####
#rm(exp2)

a<-c("euclidean","maximum","manhattan","canberra","binary","minkowski")
b<-c("ward.D","ward.D2","single","complete","average","mcquitty","median", "centroid")
zuhe<-c()
for (i in 1:6) {
  zuhe1<-combn(c(a[i],b),2)[,1:8]
  zuhe<-cbind(zuhe,zuhe1)
}
?hclust
tide<-read.csv("G:\\pan_result\\tide\\KIRP_tide_result.csv",header = T)

library(survival)
library(openxlsx)
pan_pfs<-read.xlsx("E:\\450k\\TCGAPFS.xlsx")
brca_os<-read.table("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\UCSC+KIRP_OS.txt",sep = "\t",header = T)


p_shengcun_pfs<-c()
p_shengcun_os<-c()
p_tide<-c()

i=1
for (i in 1:ncol(zuhe)) {
  hc<-hclust(dist(expscale,method = zuhe[1,i]),method = zuhe[2,i])
  clusters<-cutree(hc,k=2)
  h_class<-data.frame(sample=rownames(expscale),class=clusters)

  cluster1<-h_class[which(h_class$class=="1"),1]
  cluster2<-h_class[which(h_class$class=="2"),1]
  #cluster3<-h_class[which(h_class$class=="3"),1]
  
  cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
  cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
  #cluster3_exp<-exp2[,colnames(exp2) %in%cluster3]
  
  cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
  cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
  #cluster3_exp<-data.frame(gene=rownames(exp2),cluster3_exp)
  
  #PFS
  clusters<-as.data.frame(clusters)

  clusters_tide<-data.frame(Patient=rownames(clusters),group=clusters$clusters)

  samples<-substr(gsub('[.]','-',rownames(clusters)),1,15)
  clusters2<-data.frame(sample=samples,clusters)
  #pfs
  brca_samp<-gsub('[.]','-',substr(rownames(clusters),1,12))
  brca_pfs<-pan_pfs[pan_pfs$bcr_patient_barcode%in%brca_samp,]
  clusters<-data.frame(bcr_patient_barcode=brca_samp,clusters)
  brca_pfs<-merge(brca_pfs,clusters,by="bcr_patient_barcode")
  survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)
  p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  p_shengcun_pfs<-c(p_shengcun_pfs,p1)
  
  #os
  brca_os2<-merge(brca_os,clusters2,by="sample")
  str(brca_os2)
  fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  p2 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  # p1
  p_shengcun_os<-c(p_shengcun_os,p2)
  
  #tide

  tide_final<-left_join(tide,clusters_tide,by="Patient")
  T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
  p_tide_zanshi<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_tide<-c(p_tide,p_tide_zanshi)
  
  
  print(i)
}
hcluster_result<-cbind(p_shengcun_pfs,p_shengcun_os,p_tide)


save(hcluster_result,file = "hcluster_result.Rdata")
# save(p_shengcun_os,file = "405k_p_shengcun_os_hclust.Rdata")
# save(p_shengcun_pfs,file = "405k_p_shengcun_pfs_hclust.Rdata")




which(p_ciber$V1<0.05&p_shengcun_pfs<0.05)
#[1]  1  2  9 10 17 18 20 25 26 28 30 41 42
#p_ciber_hclust<-p_ciber


library(ConsensusClusterPlus)
STAD_tumor = sweep(t(exp), 1, apply(t(exp), 1, median, na.rm = T))
a<-c("hc","pam","km")
b<-c("pearson","spearman","euclidean","binary","maximum","canberra","minkowski")
zuhe<-c()
for (i in 1:3) {
  zuhe1<-combn(c(a[i],b),2)[,1:7]
  zuhe<-cbind(zuhe,zuhe1)
}
library(survival)
library(openxlsx)
pan_pfs<-read.xlsx("E:\\450k\\TCGAPFS.xlsx")
brca_os<-read.table("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\UCSC+KIRP_OS.txt",sep = "\t",header = T)


p_ciber<-c()
p_shengcun_pfs<-c()
p_shengcun_os<-c()


#p_ciber_hclust<-p_ciber
i=1
for (i in 1:ncol(zuhe)) {
  results = ConsensusClusterPlus(STAD_tumor, maxK = 3, reps = 1000, pItem = 0.8, pFeature = 1, title = "untitled_consensus_cluster", clusterAlg = zuhe[1,i],
                                 distance = zuhe[2,i], seed = 123, tmyPal=NULL, writeTable=FALSE, plot = "pdf")
  clusters<-results[[3]][["consensusClass"]]
  if(min(table(clusters))>15){
    h_class<-data.frame(sample=rownames(expscale),class=clusters)
  
    cluster1<-h_class[which(h_class$class=="1"),1]
    cluster2<-h_class[which(h_class$class=="2"),1]
    cluster3<-h_class[which(h_class$class=="3"),1]
    
    cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
    cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
    cluster3_exp<-exp2[,colnames(exp2) %in%cluster3]
    
    cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
    cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
    cluster3_exp<-data.frame(gene=rownames(exp2),cluster3_exp)
    
    write.table(cluster1_exp,"G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\cluster1\\h_cluster1_exp.txt",sep="\t",row.names = F)
    write.table(cluster2_exp,"G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\cluster2\\h_cluster2_exp.txt",sep="\t",row.names = F)
    write.table(cluster3_exp,"G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\cluster3\\h_cluster3_exp.txt",sep="\t",row.names = F)
    
    setwd("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\cluster1")
    source("cibersort.R")
    CIBERSORT("LM4_mean.txt","h_cluster1_exp.txt",perm = 100,QN = T,"cluster1_CIBERSORT-Results_LM4_mean.txt")
    setwd("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\cluster2")
    source("cibersort.R")
    CIBERSORT("LM4_mean.txt","h_cluster2_exp.txt",perm = 100,QN = T,"cluster2_CIBERSORT-Results_LM4_mean.txt")
    setwd("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\cluster3")
    source("cibersort.R")
    CIBERSORT("LM4_mean.txt","h_cluster3_exp.txt",perm = 100,QN = T,"cluster3_CIBERSORT-Results_LM4_mean.txt")
    
    
    cluster1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\cluster1\\cluster1_CIBERSORT-Results_LM4_mean.txt",header=T)
    cluster2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\cluster2\\cluster2_CIBERSORT-Results_LM4_mean.txt",header=T)
    cluster3_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\cluster3\\cluster3_CIBERSORT-Results_LM4_mean.txt",header=T)
    
    
    ciber_zong<-cbind(t(cluster1_ciber[,1:(ncol(cluster1_ciber)-3)]),t(cluster2_ciber[,1:(ncol(cluster1_ciber)-3)]),t(cluster3_ciber[,1:(ncol(cluster1_ciber)-3)]))
    colnames(ciber_zong)<-ciber_zong[1,]
    ciber_zong<-ciber_zong[-1,]
    ciber_zong<-as.data.frame(ciber_zong)
    ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],1,as.numeric)
    str(ciber_zong)
    p_zanshi<-c()
    
    ciber_class<-rep(c(1,2,3),c(nrow(cluster1_ciber),nrow(cluster2_ciber),nrow(cluster3_ciber)))
    for (j in 1:nrow(ciber_zong)) {
      T_result<-kruskal.test(as.numeric(ciber_zong[j,])~ciber_class)
      
      p<-T_result$p.value
      
      # t1<-T_result[[1]]
      p_zanshi<-c(p_zanshi,p)
      #t_zanshi<-c(t_zanshi,t1)
    }
  }else{
    p_zanshi<-rep(1,nrow(ciber_zong))
  }
  p_ciber<-rbind(p_ciber,p_zanshi)
  #PFS
  clusters<-as.data.frame(clusters)

  samples<-substr(gsub('[.]','-',rownames(clusters)),1,15)
  clusters2<-data.frame(sample=samples,clusters)
  #pfs
  brca_samp<-gsub('[.]','-',substr(rownames(clusters),1,12))
  brca_pfs<-pan_pfs[pan_pfs$bcr_patient_barcode%in%brca_samp,]
  clusters<-data.frame(bcr_patient_barcode=brca_samp,clusters)
  brca_pfs<-merge(brca_pfs,clusters,by="bcr_patient_barcode")
  survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)
  p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  p_shengcun_pfs<-c(p_shengcun_pfs,p1)
  
  #os
  brca_os2<-merge(brca_os,clusters2,by="sample")
  str(brca_os2)
  fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  p2 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  # p1
  p_shengcun_os<-c(p_shengcun_os,p2)
  
  
  print(i)
  
  
}
p_ciber<-as.data.frame(p_ciber)

save(p_ciber,file="45kp_ciber_consencus.Rdata")
save(p_shengcun_os,file="45kp_p_shengcun_os_consencus.Rdata")
save(p_shengcun_pfs,file="45kp_p_shengcun_pfs_consencus.Rdata")


load(file = "cluster2\\45kp_ciber_consencus.Rdata")
zuhe[,15]#[1] "km"      "pearson"
zuhe[,16]
zuhe[,17]
zuhe[,18]
zuhe[,19]
zuhe[,20]
zuhe[,21]
 
i=1
results = ConsensusClusterPlus(STAD_tumor, maxK = 2, reps = 1000, pItem = 0.8, pFeature = 1, title = "untitled_consensus_cluster", clusterAlg = zuhe[1,i],
                               distance = zuhe[2,i], seed = 123, tmyPal=NULL, writeTable=FALSE, plot = "pdf")
clusters<-results[[2]][["consensusClass"]]
h_class<-data.frame(sample=rownames(expscale),class=clusters)
 
cluster1<-h_class[which(h_class$class=="1"),1]
cluster2<-h_class[which(h_class$class=="2"),1]
cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
write.table(cluster1_exp,"G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\cluster1\\h_cluster1_exp.txt",sep="\t",row.names = F)
write.table(cluster2_exp,"G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\cluster2\\h_cluster2_exp.txt",sep="\t",row.names = F)
setwd("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\cluster1")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster1_exp.txt",perm = 100,QN = T,"cluster1_CIBERSORT-Results_LM4_mean.txt")
setwd("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\cluster2")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster2_exp.txt",perm = 100,QN = T,"cluster2_CIBERSORT-Results_LM4_mean.txt")

cluster1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\cluster1\\cluster1_CIBERSORT-Results_LM4_mean.txt",header=T)
cluster2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\cluster2\\cluster2_CIBERSORT-Results_LM4_mean.txt",header=T)

ciber_zong<-cbind(t(cluster1_ciber[,1:(ncol(cluster1_ciber)-3)]),(t(cluster2_ciber[,1:(ncol(cluster1_ciber)-3)])))
colnames(ciber_zong)<-ciber_zong[1,]
ciber_zong<-ciber_zong[-1,]
ciber_zong<-as.data.frame(ciber_zong)
ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],1,as.numeric)
str(ciber_zong)
p_zanshi<-c()
for (j in 1:nrow(ciber_zong)) {
  T_result<-wilcox.test(as.numeric(ciber_zong[j,1:nrow(cluster1_ciber)]),as.numeric(ciber_zong[j,(nrow(cluster1_ciber)+1):ncol(ciber_zong)]),paired = F,var.equal = T)
  
  p<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}


 
i=9


hc<-hclust(dist(expscale,method = zuhe[1,i]),method = zuhe[2,i])
clusters<-cutree(hc,k=2)

h_class<-data.frame(sample=rownames(expscale),class=clusters)
 
 
cluster1<-h_class[which(h_class$class=="1"),1]
cluster2<-h_class[which(h_class$class=="2"),1]
#cluster3<-h_class[which(h_class$class=="3"),1]

cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
#cluster3_exp<-exp2[,colnames(exp2) %in%cluster3]

cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
#cluster3_exp<-data.frame(gene=rownames(exp2),cluster3_exp)

write.table(cluster1_exp,"G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\cluster1\\h_cluster1_exp.txt",sep="\t",row.names = F)
write.table(cluster2_exp,"G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\cluster2\\h_cluster2_exp.txt",sep="\t",row.names = F)
#write.table(cluster3_exp,"G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\cluster3\\h_cluster3_exp.txt",sep="\t",row.names = F)

setwd("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\cluster1")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster1_exp.txt",perm = 100,QN = T,"cluster1_CIBERSORT-Results_LM4_mean.txt")
setwd("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\cluster2")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster2_exp.txt",perm = 100,QN = T,"cluster2_CIBERSORT-Results_LM4_mean.txt")
# setwd("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\cluster3")
# source("cibersort.R")
# CIBERSORT("LM4_mean.txt","h_cluster3_exp.txt",perm = 100,QN = T,"cluster3_CIBERSORT-Results_LM4_mean.txt")


cluster1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\cluster1\\cluster1_CIBERSORT-Results_LM4_mean.txt",header=T)
cluster2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\cluster2\\cluster2_CIBERSORT-Results_LM4_mean.txt",header=T)
#cluster3_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\cluster3\\cluster3_CIBERSORT-Results_LM4_mean.txt",header=T)


ciber_zong<-cbind(t(cluster1_ciber[,1:(ncol(cluster1_ciber)-3)]),t(cluster2_ciber[,1:(ncol(cluster1_ciber)-3)]))
colnames(ciber_zong)<-ciber_zong[1,]
ciber_zong<-ciber_zong[-1,]
ciber_zong<-as.data.frame(ciber_zong)
ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],1,as.numeric)
str(ciber_zong)
p_zanshi<-c()

ciber_class<-rep(c(1,2),c(nrow(cluster1_ciber),nrow(cluster2_ciber)))
for (j in 1:nrow(ciber_zong)) {
  T_result<-wilcox.test(as.numeric(ciber_zong[j,])~ciber_class)
  
  p<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_cluster_new<-p_zanshi
p_cluster_new

#LIHC####

setwd("G:\\pan_result\\gpu_result_tcga\\LIHC")
gene<-read.csv("642featuresweightme.csv",row.names = 1)
gene<-colnames(gene)
#exp<-read.csv("sigmoid.csv")
exp<-read.csv("642featuresweightme.csv",row.names = 1)

label<-read.csv("label.csv")
exp11<-exp[rownames(exp)%in%(label$X[which(label$label=="1")]),]
exp22<-exp[rownames(exp)%in%(label$X[which(label$label=="2")]),]
exp<-rbind(exp11,exp22)
# rownames(exp)<-exp$X
# exp<-exp[,-1]
expscale<-scale(exp)

library(factoextra)
setwd("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k")
set.seed(123)
fviz_nbclust<-fviz_nbclust(expscale, kmeans, k.max = 15)#2
fviz_nbclust
ggsave("LIHC_fviz_nbclust.pdf",fviz_nbclust,width = 10,height = 10)
 
# ggsave("LIHC_kmeans_plot.pdf",kmeans_plot,width = 10,height = 10)


# expr<-read.table("tcga_brca_tpm.txt")
# expr<-t(expr)
# expr<-as.data.frame(expr)
# expr$name<-substr(rownames(expr),1,15)
# library(limma)
# expr<-avereps(expr,ID=expr$name)
# ?avereps
# expr<-as.data.frame(expr)
# rownames(expr)<-expr$name
# expr<-expr[,-ncol(expr)]
# expr<-as.data.frame(t(expr))
# write.table(expr,"brca_tpm_quchong.txt")
expr1<-read.csv("G:\\dnameth\\LIHC\\exp1.csv",row.names = 1)
expr2<-read.csv("G:\\dnameth\\LIHC\\exp2.csv",row.names = 1)
expr<-cbind(expr1,expr2)
# expr1<-expr[,colnames(expr)%in%exp11$X]
# expr2<-expr[,colnames(expr)%in%exp22$X]
str(expr1)
# expr1<-apply(expr1, 2, as.numeric)
# rownames(expr1)<-rownames(expr2)
# expr2<-apply(expr2, 2, as.numeric)
# rownames(expr2)<-rownames(expr1)


library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(ggsignif)


exp_tide<-t(apply(expr,1,function(x){x-(mean(x))}))
write.table(data.frame(ID=rownames(exp_tide),exp_tide,check.names=F),'LIHC_tide.txt',sep = "\t",quote=F,row.names = FALSE)


library(estimate)
 
dat = log2(edgeR::cpm(expr)+1)
write.table(dat, file = "LIHC_estimate_input.txt", sep = '\t', quote = F)
exp.file = "LIHC_estimate_input.txt"
in.gct.file = "ESTIMATE_LIHC_input.gct"
 
outputGCT(exp.file, in.gct.file)
rt <- read.table("ESTIMATE_LIHC_input.gct", 
                 skip = 2,  
                 header = TRUE, 
                 sep = "\t")

 
filterCommonGenes(input.f = exp.file,  
                  output.f = "ESTIMATE_LIHC_gene.gct", 
                  id = "GeneSymbol"  
) 

#2esitimate 

out.score.file = "ESTIMATE_score.gct"
 
estimateScore(in.gct.file,  
              out.score.file,
              platform = "illumina")  

 
ESTIMATE_score = read.table(out.score.file,
                            skip = 2,
                            header = T,
                            row.names = 1)
ESTIMATE_score = as.data.frame(t(ESTIMATE_score[, 2:ncol(ESTIMATE_score)]))  
ESTIMATE_score[1:6, 1:3]
write.table(ESTIMATE_score, file = "LIHC_estimate_score.txt", sep = '\t', quote = F)






getwd()



new_folder <- "G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\stage1"
dir.create(new_folder)

source_folder <- "G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\stage1"
r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/LIHC/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)
write.table(expr1,"stage1\\LIHC_tpm_stage1.txt",sep = "\t")



new_folder <- "G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\stage2"
dir.create(new_folder)

source_folder <- "G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\stage2"
r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/LIHC/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)
write.table(expr2,"stage2\\LIHC_tpm_stage2.txt",sep = "\t")




setwd("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\stage1")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","LIHC_tpm_stage1.txt",perm = 100,QN = T,"stage1_CIBERSORT-Results_LM4_mean.txt")

setwd("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\stage2")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","LIHC_tpm_stage2.txt",perm = 100,QN = T,"stage2_CIBERSORT-Results_LM4_mean.txt")


stage1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\stage1\\stage1_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")
stage2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\stage2\\stage2_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")


ciber_zong<-cbind(t(stage1_ciber[,1:(ncol(stage1_ciber)-3)]),(t(stage2_ciber[,1:(ncol(stage1_ciber)-3)])))
colnames(ciber_zong)<-ciber_zong[1,]
ciber_zong<-ciber_zong[-1,]
ciber_zong<-as.data.frame(ciber_zong)
ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],2,as.numeric)
str(ciber_zong)
p_zanshi<-c()
for (j in 1:nrow(ciber_zong)) {
  T_result<-wilcox.test(as.numeric(ciber_zong[j,1:nrow(stage1_ciber)]),as.numeric(ciber_zong[j,(nrow(stage1_ciber)+1):ncol(ciber_zong)]),paired = F,var.equal = T)
  
  p<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_zanshi
p_stage_yuanshi<-p_zanshi
#[1] 0.23361657        NaN 0.06673382 0.04491467 0.45374285 0.69362902 0.86320080 0.72045094 0.09366834
p_stage_yuanshi
# 
library(ggplot2)
library(ggthemes)
class<-rep(c("zaoqi","jinzhan"),c(nrow(stage1_ciber),nrow(stage2_ciber)))
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
library(dplyr)
library(ggpubr)
aaa<-ggplot(mianyi_huatu,aes(y=scores,x=factor(type),fill=class))+
  geom_boxplot(outlier.size=1)+
  theme_classic()+
  xlab("")+ylab("")+guides()+
  scale_fill_manual(values=c("#eed65d", "#5d75ee"))+stat_compare_means(method = "wilcox.test",label = "p.signif")
# geom_text(aes(x =1, y=1, label = label[1], group = type), size = 4) +
# geom_text(aes(x =2, y=1, label = label[2], group = type), size = 4)+
# geom_text(aes(x =3, y=1, label = label[3], group = type), size = 4)+
# geom_text(aes(x =4, y=1, label = label[4], group = type), size = 4)+
# geom_text(aes(x =5, y=1, label = label[5], group = type), size = 4)+
# geom_text(aes(x =6, y=1, label = label[6], group = type), size = 4)+
# geom_text(aes(x =7, y=1, label = label[7], group = type), size = 4)
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\LIHC_stage_scrna_plot.pdf",aaa,width = 10,height = 10)



# cluster1
new_folder <- "G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\cluster1"
dir.create(new_folder)

r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/LIHC/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)


# cluster2 
new_folder <- "G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\cluster2"
dir.create(new_folder)

r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/LIHC/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)

# 
# 
# # cluster3 
# new_folder <- "G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\cluster3"
# dir.create(new_folder)
# 
# r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)
# 
# file.copy(r_files, new_folder)
# 
# source_file <- "G:/pan_result/gpu_result_tcga/LIHC/scrna/LM4_mean.txt"
# file.copy(source_file, new_folder)
# 

exp2<-cbind(expr1,expr2)
#1.2 ####
#rm(exp2)

a<-c("euclidean","maximum","manhattan","canberra","binary","minkowski")
b<-c("ward.D","ward.D2","single","complete","average","mcquitty","median", "centroid")
zuhe<-c()
for (i in 1:6) {
  zuhe1<-combn(c(a[i],b),2)[,1:8]
  zuhe<-cbind(zuhe,zuhe1)
}
?hclust
tide<-read.csv("G:\\pan_result\\tide\\LIHC_tide_result.csv",header = T)

library(survival)
library(openxlsx)
pan_pfs<-read.xlsx("E:\\450k\\TCGAPFS.xlsx")
brca_os<-read.table("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\UCSC+LIHC_OS.txt",sep = "\t",header = T)


p_shengcun_pfs<-c()
p_shengcun_os<-c()
p_tide<-c()

i=1
for (i in 1:ncol(zuhe)) {
  hc<-hclust(dist(expscale,method = zuhe[1,i]),method = zuhe[2,i])
  clusters<-cutree(hc,k=2)
  h_class<-data.frame(sample=rownames(expscale),class=clusters)

  cluster1<-h_class[which(h_class$class=="1"),1]
  cluster2<-h_class[which(h_class$class=="2"),1]
  #cluster3<-h_class[which(h_class$class=="3"),1]
  
  cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
  cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
  #cluster3_exp<-exp2[,colnames(exp2) %in%cluster3]
  
  cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
  cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
  #cluster3_exp<-data.frame(gene=rownames(exp2),cluster3_exp)
  
  #PFS
  clusters<-as.data.frame(clusters)

  clusters_tide<-data.frame(Patient=rownames(clusters),group=clusters$clusters)

  samples<-substr(gsub('[.]','-',rownames(clusters)),1,15)
  clusters2<-data.frame(sample=samples,clusters)
  #pfs
  brca_samp<-gsub('[.]','-',substr(rownames(clusters),1,12))
  brca_pfs<-pan_pfs[pan_pfs$bcr_patient_barcode%in%brca_samp,]
  clusters<-data.frame(bcr_patient_barcode=brca_samp,clusters)
  brca_pfs<-merge(brca_pfs,clusters,by="bcr_patient_barcode")
  survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)
  p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  p_shengcun_pfs<-c(p_shengcun_pfs,p1)
  
  #os
  brca_os2<-merge(brca_os,clusters2,by="sample")
  str(brca_os2)
  fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  p2 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  # p1
  p_shengcun_os<-c(p_shengcun_os,p2)
  
  #tide

  tide_final<-left_join(tide,clusters_tide,by="Patient")
  T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
  p_tide_zanshi<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_tide<-c(p_tide,p_tide_zanshi)
  
  
  print(i)
}
hcluster_result<-cbind(p_shengcun_pfs,p_shengcun_os,p_tide)


save(hcluster_result,file = "hcluster_result.Rdata")
# save(p_shengcun_os,file = "405k_p_shengcun_os_hclust.Rdata")
# save(p_shengcun_pfs,file = "405k_p_shengcun_pfs_hclust.Rdata")




which(p_ciber$V1<0.05&p_shengcun_pfs<0.05)
#[1]  1  2  9 10 17 18 20 25 26 28 30 41 42
#p_ciber_hclust<-p_ciber


library(ConsensusClusterPlus)
STAD_tumor = sweep(t(exp), 1, apply(t(exp), 1, median, na.rm = T))
a<-c("hc","pam","km")
b<-c("pearson","spearman","euclidean","binary","maximum","canberra","minkowski")
zuhe<-c()
for (i in 1:3) {
  zuhe1<-combn(c(a[i],b),2)[,1:7]
  zuhe<-cbind(zuhe,zuhe1)
}


library(survival)
library(openxlsx)
pan_pfs<-read.xlsx("E:\\450k\\TCGAPFS.xlsx")
brca_os<-read.table("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\UCSC+LIHC_OS.txt",sep = "\t",header = T)
if (is.null(STAD_tumor) || nrow(STAD_tumor) == 0) {
  stop("STAD_tumor is null or empty.")
}

if (is.null(zuhe) || ncol(zuhe) < 2 || i > ncol(zuhe)) {
  stop("Check zuhe matrix and index i.")
}
p_shengcun_pfs<-c()
p_shengcun_os<-c()
p_tide<-c()

#p_ciber_hclust<-p_ciber
i=2
for (i in 1:ncol(zuhe)) {
  results = ConsensusClusterPlus(STAD_tumor, maxK = 3, reps = 1000, pItem = 0.8, pFeature = 1, title = "untitled_consensus_cluster", clusterAlg = zuhe[1,i],
                                 distance = zuhe[2,i], seed = 123, tmyPal=NULL, writeTable=FALSE, plot = "pdf")
  clusters<-results[[2]][["consensusClass"]]
  h_class<-data.frame(sample=rownames(expscale),class=clusters)

  cluster1<-h_class[which(h_class$class=="1"),1]
  cluster2<-h_class[which(h_class$class=="2"),1]
  #cluster3<-h_class[which(h_class$class=="3"),1]
  
  cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
  cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
  #cluster3_exp<-exp2[,colnames(exp2) %in%cluster3]
  
  cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
  cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
  #cluster3_exp<-data.frame(gene=rownames(exp2),cluster3_exp)
  
  #PFS
  clusters<-as.data.frame(clusters)

  clusters_tide<-data.frame(Patient=rownames(clusters),group=clusters$clusters)

  samples<-substr(gsub('[.]','-',rownames(clusters)),1,15)
  clusters2<-data.frame(sample=samples,clusters)
  #pfs
  brca_samp<-gsub('[.]','-',substr(rownames(clusters),1,12))
  brca_pfs<-pan_pfs[pan_pfs$bcr_patient_barcode%in%brca_samp,]
  clusters<-data.frame(bcr_patient_barcode=brca_samp,clusters)
  brca_pfs<-merge(brca_pfs,clusters,by="bcr_patient_barcode")
  survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)
  p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  p_shengcun_pfs<-c(p_shengcun_pfs,p1)
  
  #os
  brca_os2<-merge(brca_os,clusters2,by="sample")
  str(brca_os2)
  fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  p2 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  # p1
  p_shengcun_os<-c(p_shengcun_os,p2)
  
  #tide

  tide_final<-left_join(tide,clusters_tide,by="Patient")
  T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
  p_tide_zanshi<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_tide<-c(p_tide,p_tide_zanshi)
  
  print(i)
  
  
}
concluster_result<-cbind(p_shengcun_pfs,p_shengcun_os,p_tide)


save(concluster_result,file = "concluster_result.Rdata")
p_ciber<-as.data.frame(p_ciber)

# save(p_ciber,file="45kp_ciber_consencus.Rdata")
# save(p_shengcun_os,file="45kp_p_shengcun_os_consencus.Rdata")
# save(p_shengcun_pfs,file="45kp_p_shengcun_pfs_consencus.Rdata")

setwd("G:/pan_result/gpu_result_tcga/LIHC/450k/")
load(file = "cluster2\\45kp_ciber_consencus.Rdata")
zuhe[,15]#[1] "km"      "pearson"
zuhe[,16]
zuhe[,17]
zuhe[,18]
zuhe[,19]
zuhe[,20]
zuhe[,21]
 
i=10
results = ConsensusClusterPlus(STAD_tumor, maxK = 3, reps = 1000, pItem = 0.8, pFeature = 1, title = "untitled_consensus_cluster", clusterAlg = zuhe[1,i],
                               distance = zuhe[2,i], seed = 123, tmyPal=NULL, writeTable=FALSE, plot = "pdf")
clusters<-results[[2]][["consensusClass"]]
table(clusters)
clusters2<-clusters
clusters2[which(clusters=="1")]=2
clusters2[which(clusters=="2")]=1
clusters<-clusters2
table(clusters)
h_class<-data.frame(sample=rownames(expscale),class=clusters)
 
cluster1<-h_class[which(h_class$class=="1"),1]
cluster2<-h_class[which(h_class$class=="2"),1]
cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
write.table(cluster1_exp,"G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\cluster1\\h_cluster1_exp.txt",sep="\t",row.names = F)
write.table(cluster2_exp,"G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\cluster2\\h_cluster2_exp.txt",sep="\t",row.names = F)
setwd("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\cluster1")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster1_exp.txt",perm = 100,QN = T,"cluster1_CIBERSORT-Results_LM4_mean.txt")
setwd("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\cluster2")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster2_exp.txt",perm = 100,QN = T,"cluster2_CIBERSORT-Results_LM4_mean.txt")

cluster1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\cluster1\\cluster1_CIBERSORT-Results_LM4_mean.txt",header=T)
cluster2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\cluster2\\cluster2_CIBERSORT-Results_LM4_mean.txt",header=T)

ciber_zong<-cbind(t(cluster1_ciber[,1:(ncol(cluster1_ciber)-3)]),(t(cluster2_ciber[,1:(ncol(cluster1_ciber)-3)])))
colnames(ciber_zong)<-ciber_zong[1,]
ciber_zong<-ciber_zong[-1,]
ciber_zong<-as.data.frame(ciber_zong)
ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],1,as.numeric)
str(ciber_zong)
p_zanshi<-c()
for (j in 1:nrow(ciber_zong)) {
  T_result<-wilcox.test(as.numeric(ciber_zong[j,1:nrow(cluster1_ciber)]),as.numeric(ciber_zong[j,(nrow(cluster1_ciber)+1):ncol(ciber_zong)]),paired = F,var.equal = T)
  
  p<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}


 
i=2

hc<-hclust(dist(expscale,method = zuhe[1,i]),method = zuhe[2,i])
clusters<-cutree(hc,k=2)
table(clusters)
h_class<-data.frame(sample=rownames(expscale),class=clusters)
table(h_class$class)
 
 
cluster1<-h_class[which(h_class$class=="1"),1]
cluster2<-h_class[which(h_class$class=="2"),1]
#cluster3<-h_class[which(h_class$class=="3"),1]

cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
#cluster3_exp<-exp2[,colnames(exp2) %in%cluster3]

cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
#cluster3_exp<-data.frame(gene=rownames(exp2),cluster3_exp)

write.table(cluster1_exp,"G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\cluster1\\h_cluster1_exp.txt",sep="\t",row.names = F)
write.table(cluster2_exp,"G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\cluster2\\h_cluster2_exp.txt",sep="\t",row.names = F)
#write.table(cluster3_exp,"G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\cluster3\\h_cluster3_exp.txt",sep="\t",row.names = F)

setwd("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\cluster1")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster1_exp.txt",perm = 100,QN = T,"cluster1_CIBERSORT-Results_LM4_mean.txt")
setwd("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\cluster2")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster2_exp.txt",perm = 100,QN = T,"cluster2_CIBERSORT-Results_LM4_mean.txt")
# setwd("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\cluster3")
# source("cibersort.R")
# CIBERSORT("LM4_mean.txt","h_cluster3_exp.txt",perm = 100,QN = T,"cluster3_CIBERSORT-Results_LM4_mean.txt")


cluster1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\cluster1\\cluster1_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")
cluster2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\cluster2\\cluster2_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")
#cluster3_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\cluster3\\cluster3_CIBERSORT-Results_LM4_mean.txt",header=T)


ciber_zong<-cbind(t(cluster1_ciber[,1:(ncol(cluster1_ciber)-3)]),t(cluster2_ciber[,1:(ncol(cluster1_ciber)-3)]))
colnames(ciber_zong)<-ciber_zong[1,]
ciber_zong<-ciber_zong[-1,]
ciber_zong<-as.data.frame(ciber_zong)
ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],1,as.numeric)
str(ciber_zong)
p_zanshi<-c()

ciber_class<-rep(c(1,2),c(nrow(cluster1_ciber),nrow(cluster2_ciber)))
for (j in 1:nrow(ciber_zong)) {
  T_result<-wilcox.test(as.numeric(ciber_zong[j,])~ciber_class)
  
  p<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_cluster_new<-p_zanshi
p_cluster_new

#LUNG####

setwd("G:\\pan_result\\gpu_result_tcga\\NSCLC")
gene<-read.csv("256_64_789features_weightme.csv",row.names = 1)
gene<-colnames(gene)
#exp<-read.csv("G:\\pan_result\\gpu_result_tcga\\NSCLC\\sigmoid.csv")
exp<-read.csv("256_64_789features_weightme.csv",row.names = 1)

label<-read.csv("label.csv")
exp11<-exp[rownames(exp)%in%(label$X[which(label$label=="1")]),]
exp22<-exp[rownames(exp)%in%(label$X[which(label$label=="2")]),]
exp<-rbind(exp11,exp22)
# rownames(exp)<-exp$X
# exp<-exp[,-1]
expscale<-scale(exp)

library(factoextra)
setwd("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k")
set.seed(1)
fviz_nbclust<-fviz_nbclust(exp, kmeans, k.max = 15)#2
fviz_nbclust
ggsave("NSCLC_fviz_nbclust.pdf",fviz_nbclust,width = 10,height = 10)
 
# ggsave("NSCLC_kmeans_plot.pdf",kmeans_plot,width = 10,height = 10)


# expr<-read.table("tcga_brca_tpm.txt")
# expr<-t(expr)
# expr<-as.data.frame(expr)
# expr$name<-substr(rownames(expr),1,15)
# library(limma)
# expr<-avereps(expr,ID=expr$name)
# ?avereps
# expr<-as.data.frame(expr)
# rownames(expr)<-expr$name
# expr<-expr[,-ncol(expr)]
# expr<-as.data.frame(t(expr))
# write.table(expr,"brca_tpm_quchong.txt")
expr1<-read.csv("G:\\pan_result\\gpu_result_tcga\\NSCLC\\exp1.csv",row.names = 1)
expr2<-read.csv("G:\\pan_result\\gpu_result_tcga\\NSCLC\\exp2.csv",row.names = 1)
expr<-cbind(expr1,expr2)
# expr1<-expr[,colnames(expr)%in%exp11$X]
# expr2<-expr[,colnames(expr)%in%exp22$X]
str(expr1)
# expr1<-apply(expr1, 2, as.numeric)
# rownames(expr1)<-rownames(expr2)
# expr2<-apply(expr2, 2, as.numeric)
# rownames(expr2)<-rownames(expr1)


library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(ggsignif)


exp_tide<-t(apply(expr,1,function(x){x-(mean(x))}))
write.table(data.frame(ID=rownames(exp_tide),exp_tide,check.names=F),'NSCLC_tide.txt',sep = "\t",quote=F,row.names = FALSE)

 
library(estimate)
 
dat = log2(edgeR::cpm(expr)+1)
write.table(dat, file = "NSCLC_estimate_input.txt", sep = '\t', quote = F)
exp.file = "NSCLC_estimate_input.txt"
in.gct.file = "ESTIMATE_NSCLC_input.gct"

outputGCT(exp.file, in.gct.file)
rt <- read.table("ESTIMATE_NSCLC_input.gct", 
                 skip = 2,  
                 header = TRUE, 
                 sep = "\t")

 
filterCommonGenes(input.f = exp.file,  
                  output.f = "ESTIMATE_NSCLC_gene.gct",  
                  id = "GeneSymbol" 
) 

#2esitimate 

out.score.file = "ESTIMATE_score.gct"
 
estimateScore(in.gct.file,  
              out.score.file, 
              platform = "illumina") 

 
ESTIMATE_score = read.table(out.score.file,
                            skip = 2,
                            header = T,
                            row.names = 1)
ESTIMATE_score = as.data.frame(t(ESTIMATE_score[, 2:ncol(ESTIMATE_score)]))  
ESTIMATE_score[1:6, 1:3]
write.table(ESTIMATE_score, file = "NSCLC_estimate_score.txt", sep = '\t', quote = F)











getwd()



new_folder <- "G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\stage1"
dir.create(new_folder)

source_folder <- "G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\stage1"
r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/NSCLC/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)
write.table(expr1,"stage1\\NSCLC_tpm_stage1.txt",sep = "\t")



new_folder <- "G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\stage2"
dir.create(new_folder)

source_folder <- "G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\stage2"
r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/NSCLC/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)
write.table(expr2,"stage2\\NSCLC_tpm_stage2.txt",sep = "\t")




setwd("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\stage1")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","NSCLC_tpm_stage1.txt",perm = 100,QN = T,"stage1_CIBERSORT-Results_LM4_mean.txt")

setwd("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\stage2")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","NSCLC_tpm_stage2.txt",perm = 100,QN = T,"stage2_CIBERSORT-Results_LM4_mean.txt")


stage1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\stage1\\stage1_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")
stage2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\stage2\\stage2_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")


ciber_zong<-cbind(t(stage1_ciber[,1:(ncol(stage1_ciber)-3)]),(t(stage2_ciber[,1:(ncol(stage1_ciber)-3)])))
colnames(ciber_zong)<-ciber_zong[1,]
ciber_zong<-ciber_zong[-1,]
ciber_zong<-as.data.frame(ciber_zong)
ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],2,as.numeric)
str(ciber_zong)
p_zanshi<-c()
for (j in 1:nrow(ciber_zong)) {
  T_result<-wilcox.test(as.numeric(ciber_zong[j,1:nrow(stage1_ciber)]),as.numeric(ciber_zong[j,(nrow(stage1_ciber)+1):ncol(ciber_zong)]),paired = F,var.equal = T)
  
  p<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_zanshi
p_stage_yuanshi<-p_zanshi
#[1] 0.23361657        NaN 0.06673382 0.04491467 0.45374285 0.69362902 0.86320080 0.72045094 0.09366834
p_stage_yuanshi
# 
library(ggplot2)
library(ggthemes)
class<-rep(c("zaoqi","jinzhan"),c(nrow(stage1_ciber),nrow(stage2_ciber)))
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
library(dplyr)
library(ggpubr)
aaa<-ggplot(mianyi_huatu,aes(y=scores,x=factor(type),fill=class))+
  geom_boxplot(outlier.size=1)+
  theme_classic()+
  xlab("")+ylab("")+guides()+
  scale_fill_manual(values=c("#eed65d", "#5d75ee"))+stat_compare_means(method = "wilcox.test",label = "p.signif")
# geom_text(aes(x =1, y=1, label = label[1], group = type), size = 4) +
# geom_text(aes(x =2, y=1, label = label[2], group = type), size = 4)+
# geom_text(aes(x =3, y=1, label = label[3], group = type), size = 4)+
# geom_text(aes(x =4, y=1, label = label[4], group = type), size = 4)+
# geom_text(aes(x =5, y=1, label = label[5], group = type), size = 4)+
# geom_text(aes(x =6, y=1, label = label[6], group = type), size = 4)+
# geom_text(aes(x =7, y=1, label = label[7], group = type), size = 4)
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\NSCLC_stage_scrna_plot.pdf",aaa,width = 10,height = 10)



# cluster1 
new_folder <- "G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\cluster1"
dir.create(new_folder)

r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/NSCLC/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)


# cluster2 
new_folder <- "G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\cluster2"
dir.create(new_folder)

r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/NSCLC/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)

# 
# 
# # cluster3 
# new_folder <- "G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\cluster3"
# dir.create(new_folder)
# 
# r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)
# 
# file.copy(r_files, new_folder)
# 
# source_file <- "G:/pan_result/gpu_result_tcga/NSCLC/scrna/LM4_mean.txt"
# file.copy(source_file, new_folder)
# 

exp2<-cbind(expr1,expr2)
#1.2  ####
#rm(exp2)

a<-c("euclidean","maximum","manhattan","canberra","binary","minkowski")
b<-c("ward.D","ward.D2","single","complete","average","mcquitty","median", "centroid")
zuhe<-c()
for (i in 1:6) {
  zuhe1<-combn(c(a[i],b),2)[,1:8]
  zuhe<-cbind(zuhe,zuhe1)
}
?hclust
tide<-read.csv("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\NSCLC_tide_result.csv",header = T)

library(survival)
library(openxlsx)
pan_pfs<-read.xlsx("E:\\450k\\TCGAPFS.xlsx")
luad_os<-read_tsv("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\TCGA-LUAD.survival.tsv")
lusc_os<-read_tsv("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\TCGA-LUSC.survival.tsv")
brca_os<-rbind(luad_os,lusc_os)


p_shengcun_pfs<-c()
p_shengcun_os<-c()
p_tide<-c()

i=1
for (i in 1:ncol(zuhe)) {
  hc<-hclust(dist(exp,method = zuhe[1,i]),method = zuhe[2,i])
  clusters<-cutree(hc,k=2)
  h_class<-data.frame(sample=rownames(exp),class=clusters)

  cluster1<-h_class[which(h_class$class=="1"),1]
  cluster2<-h_class[which(h_class$class=="2"),1]
  #cluster3<-h_class[which(h_class$class=="3"),1]
  
  cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
  cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
  #cluster3_exp<-exp2[,colnames(exp2) %in%cluster3]
  
  cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
  cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
  #cluster3_exp<-data.frame(gene=rownames(exp2),cluster3_exp)
  
  #PFS
  clusters<-as.data.frame(clusters)

  clusters_tide<-data.frame(Patient=rownames(clusters),group=clusters$clusters)

  #samples<-substr(gsub('[.]','-',rownames(clusters)),1,15)
  clusters2<-data.frame(sample=gsub('[.]','-',rownames(clusters)),clusters)
  #pfs
  brca_samp<-gsub('[.]','-',substr(rownames(clusters),1,12))
  brca_pfs<-pan_pfs[pan_pfs$bcr_patient_barcode%in%brca_samp,]
  clusters<-data.frame(bcr_patient_barcode=brca_samp,clusters)
  brca_pfs<-merge(brca_pfs,clusters,by="bcr_patient_barcode")
  survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)
  p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  p_shengcun_pfs<-c(p_shengcun_pfs,p1)
  
  #os
  brca_os2<-merge(brca_os,clusters2,by="sample")
  str(brca_os2)
  fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  p2 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  # p1
  p_shengcun_os<-c(p_shengcun_os,p2)
  
  #tide

  tide_final<-left_join(tide,clusters_tide,by="Patient")
  T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
  p_tide_zanshi<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_tide<-c(p_tide,p_tide_zanshi)
  
  
  print(i)
}
hcluster_result<-cbind(p_shengcun_pfs,p_shengcun_os,p_tide)


save(hcluster_result,file = "hcluster_result.Rdata")
# save(p_shengcun_os,file = "405k_p_shengcun_os_hclust.Rdata")
# save(p_shengcun_pfs,file = "405k_p_shengcun_pfs_hclust.Rdata")




which(p_ciber$V1<0.05&p_shengcun_pfs<0.05)
#[1]  1  2  9 10 17 18 20 25 26 28 30 41 42
#p_ciber_hclust<-p_ciber


library(ConsensusClusterPlus)
STAD_tumor = sweep(t(exp), 1, apply(t(exp), 1, median, na.rm = T))
a<-c("hc","pam","km")
b<-c("pearson","spearman","euclidean","binary","maximum","canberra","minkowski")
zuhe<-c()
for (i in 1:3) {
  zuhe1<-combn(c(a[i],b),2)[,1:7]
  zuhe<-cbind(zuhe,zuhe1)
}


library(survival)
library(openxlsx)
# pan_pfs<-read.xlsx("E:\\450k\\TCGAPFS.xlsx")
# brca_os<-read.table("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\UCSC+NSCLC_OS.txt",sep = "\t",header = T)
if (is.null(STAD_tumor) || nrow(STAD_tumor) == 0) {
  stop("STAD_tumor is null or empty.")
}

if (is.null(zuhe) || ncol(zuhe) < 2 || i > ncol(zuhe)) {
  stop("Check zuhe matrix and index i.")
}
p_shengcun_pfs<-c()
p_shengcun_os<-c()
p_tide<-c()

#p_ciber_hclust<-p_ciber
i=2
for (i in 1:ncol(zuhe)) {
  results = ConsensusClusterPlus(STAD_tumor, maxK = 3, reps = 1000, pItem = 0.8, pFeature = 1, title = "untitled_consensus_cluster", clusterAlg = zuhe[1,i],
                                 distance = zuhe[2,i], seed = 123, tmyPal=NULL, writeTable=FALSE, plot = "pdf")
  clusters<-results[[2]][["consensusClass"]]
  h_class<-data.frame(sample=rownames(expscale),class=clusters)

  cluster1<-h_class[which(h_class$class=="1"),1]
  cluster2<-h_class[which(h_class$class=="2"),1]
  #cluster3<-h_class[which(h_class$class=="3"),1]
  
  cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
  cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
  #cluster3_exp<-exp2[,colnames(exp2) %in%cluster3]
  
  cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
  cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
  #cluster3_exp<-data.frame(gene=rownames(exp2),cluster3_exp)
  
  #PFS
  clusters<-as.data.frame(clusters)

  clusters_tide<-data.frame(Patient=rownames(clusters),group=clusters$clusters)

  #samples<-substr(gsub('[.]','-',rownames(clusters)),1,15)
  clusters2<-data.frame(sample=gsub('[.]','-',rownames(clusters)),clusters)
  #pfs
  brca_samp<-gsub('[.]','-',substr(rownames(clusters),1,12))
  brca_pfs<-pan_pfs[pan_pfs$bcr_patient_barcode%in%brca_samp,]
  clusters<-data.frame(bcr_patient_barcode=brca_samp,clusters)
  brca_pfs<-merge(brca_pfs,clusters,by="bcr_patient_barcode")
  survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)
  p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  p_shengcun_pfs<-c(p_shengcun_pfs,p1)
  
  #os
  brca_os2<-merge(brca_os,clusters2,by="sample")
  str(brca_os2)
  fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  p2 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  # p1
  p_shengcun_os<-c(p_shengcun_os,p2)
  
  #tide

  tide_final<-left_join(tide,clusters_tide,by="Patient")
  T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
  p_tide_zanshi<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_tide<-c(p_tide,p_tide_zanshi)
  
  print(i)
  
  
}
p_ciber<-as.data.frame(p_ciber)

save(p_tide,file="45kp_tide_consencus.Rdata")
save(p_shengcun_os,file="45kp_p_shengcun_os_consencus.Rdata")
save(p_shengcun_pfs,file="45kp_p_shengcun_pfs_consencus.Rdata")
concluster_result<-cbind(p_shengcun_pfs,p_shengcun_os,p_tide)

load(file = "cluster2\\45kp_ciber_consencus.Rdata")
zuhe[,15]#[1] "km"      "pearson"
zuhe[,16]
zuhe[,17]
zuhe[,18]
zuhe[,19]
zuhe[,20]
zuhe[,21]
 
i=1
results = ConsensusClusterPlus(STAD_tumor, maxK = 3, reps = 1000, pItem = 0.8, pFeature = 1, title = "untitled_consensus_cluster", clusterAlg = zuhe[1,i],
                               distance = zuhe[2,i], seed = 123, tmyPal=NULL, writeTable=FALSE, plot = "pdf")
clusters<-results[[2]][["consensusClass"]]
h_class<-data.frame(sample=rownames(expscale),class=clusters)
 
cluster1<-h_class[which(h_class$class=="1"),1]
cluster2<-h_class[which(h_class$class=="2"),1]
cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
write.table(cluster1_exp,"G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\cluster1\\h_cluster1_exp.txt",sep="\t",row.names = F)
write.table(cluster2_exp,"G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\cluster2\\h_cluster2_exp.txt",sep="\t",row.names = F)
setwd("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\cluster1")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster1_exp.txt",perm = 100,QN = T,"cluster1_CIBERSORT-Results_LM4_mean.txt")
setwd("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\cluster2")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster2_exp.txt",perm = 100,QN = T,"cluster2_CIBERSORT-Results_LM4_mean.txt")

cluster1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\cluster1\\cluster1_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")
cluster2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\cluster2\\cluster2_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")

ciber_zong<-cbind(t(cluster1_ciber[,1:(ncol(cluster1_ciber)-3)]),(t(cluster2_ciber[,1:(ncol(cluster1_ciber)-3)])))
colnames(ciber_zong)<-ciber_zong[1,]
ciber_zong<-ciber_zong[-1,]
ciber_zong<-as.data.frame(ciber_zong)
ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],1,as.numeric)
str(ciber_zong)
p_zanshi<-c()
for (j in 1:nrow(ciber_zong)) {
  T_result<-wilcox.test(as.numeric(ciber_zong[j,1:nrow(cluster1_ciber)]),as.numeric(ciber_zong[j,(nrow(cluster1_ciber)+1):ncol(ciber_zong)]),paired = F,var.equal = T)
  
  p<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}


 
i=9
hc<-hclust(dist(expscale,method = zuhe[1,i]),method = zuhe[2,i])
clusters<-cutree(hc,k=2)

h_class<-data.frame(sample=rownames(expscale),class=clusters)
table(h_class$class)
 
 
cluster1<-h_class[which(h_class$class=="1"),1]
cluster2<-h_class[which(h_class$class=="2"),1]
#cluster3<-h_class[which(h_class$class=="3"),1]

cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
#cluster3_exp<-exp2[,colnames(exp2) %in%cluster3]

cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
#cluster3_exp<-data.frame(gene=rownames(exp2),cluster3_exp)

write.table(cluster1_exp,"G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\cluster1\\h_cluster1_exp.txt",sep="\t",row.names = F)
write.table(cluster2_exp,"G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\cluster2\\h_cluster2_exp.txt",sep="\t",row.names = F)
#write.table(cluster3_exp,"G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\cluster3\\h_cluster3_exp.txt",sep="\t",row.names = F)

setwd("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\cluster1")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster1_exp.txt",perm = 100,QN = T,"cluster1_CIBERSORT-Results_LM4_mean.txt")
setwd("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\cluster2")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster2_exp.txt",perm = 100,QN = T,"cluster2_CIBERSORT-Results_LM4_mean.txt")
# setwd("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\cluster3")
# source("cibersort.R")
# CIBERSORT("LM4_mean.txt","h_cluster3_exp.txt",perm = 100,QN = T,"cluster3_CIBERSORT-Results_LM4_mean.txt")


cluster1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\cluster1\\cluster1_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")
cluster2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\cluster2\\cluster2_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")
#cluster3_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\cluster3\\cluster3_CIBERSORT-Results_LM4_mean.txt",header=T)


ciber_zong<-cbind(t(cluster1_ciber[,1:(ncol(cluster1_ciber)-3)]),t(cluster2_ciber[,1:(ncol(cluster1_ciber)-3)]))
colnames(ciber_zong)<-ciber_zong[1,]
ciber_zong<-ciber_zong[-1,]
ciber_zong<-as.data.frame(ciber_zong)
ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],1,as.numeric)
str(ciber_zong)
p_zanshi<-c()

ciber_class<-rep(c(1,2),c(nrow(cluster1_ciber),nrow(cluster2_ciber)))
for (j in 1:nrow(ciber_zong)) {
  T_result<-wilcox.test(as.numeric(ciber_zong[j,])~ciber_class)
  
  p<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_cluster_new<-p_zanshi
p_cluster_new

#PAAD####

setwd("G:\\pan_result\\gpu_result_tcga\\PAAD")
gene<-read.csv("67featuresweightme.csv",row.names = 1)
gene<-colnames(gene)
#exp<-read.csv("sigmoid.csv")
exp<-read.csv("67featuresweightme.csv",row.names = 1)

label<-read.csv("label.csv")
exp11<-exp[rownames(exp)%in%(label$X[which(label$label=="1")]),]
exp22<-exp[rownames(exp)%in%(label$X[which(label$label=="2")]),]
exp<-rbind(exp11,exp22)
# rownames(exp)<-exp$X
# exp<-exp[,-1]
expscale<-scale(exp)

library(factoextra)
dir.create("450k")
setwd("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k")
set.seed(123)
fviz_nbclust<-fviz_nbclust(expscale, kmeans, k.max = 15)#2
fviz_nbclust
ggsave("PAAD_fviz_nbclust.pdf",fviz_nbclust,width = 10,height = 10)
 
# ggsave("PAAD_kmeans_plot.pdf",kmeans_plot,width = 10,height = 10)


# expr<-read.table("tcga_brca_tpm.txt")
# expr<-t(expr)
# expr<-as.data.frame(expr)
# expr$name<-substr(rownames(expr),1,15)
# library(limma)
# expr<-avereps(expr,ID=expr$name)
# ?avereps
# expr<-as.data.frame(expr)
# rownames(expr)<-expr$name
# expr<-expr[,-ncol(expr)]
# expr<-as.data.frame(t(expr))
# write.table(expr,"brca_tpm_quchong.txt")
expr1<-read.csv("G:\\dnameth\\PAAD\\exp1.csv",row.names = 1)
expr2<-read.csv("G:\\dnameth\\PAAD\\exp2.csv",row.names = 1)
expr<-cbind(expr1,expr2)
# expr1<-expr[,colnames(expr)%in%exp11$X]
# expr2<-expr[,colnames(expr)%in%exp22$X]
str(expr1)
# expr1<-apply(expr1, 2, as.numeric)
# rownames(expr1)<-rownames(expr2)
# expr2<-apply(expr2, 2, as.numeric)
# rownames(expr2)<-rownames(expr1)


library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(ggsignif)


exp_tide<-t(apply(expr,1,function(x){x-(mean(x))}))
write.table(data.frame(ID=rownames(exp_tide),exp_tide,check.names=F),'PAAD_tide.txt',sep = "\t",quote=F,row.names = FALSE)


getwd()



new_folder <- "G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\stage1"
dir.create(new_folder)

source_folder <- "G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\stage1"
r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/PAAD/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)
write.table(expr1,"stage1\\PAAD_tpm_stage1.txt",sep = "\t")



new_folder <- "G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\stage2"
dir.create(new_folder)

source_folder <- "G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\stage2"
r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/PAAD/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)
write.table(expr2,"stage2\\PAAD_tpm_stage2.txt",sep = "\t")




setwd("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\stage1")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","PAAD_tpm_stage1.txt",perm = 100,QN = T,"stage1_CIBERSORT-Results_LM4_mean.txt")

setwd("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\stage2")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","PAAD_tpm_stage2.txt",perm = 100,QN = T,"stage2_CIBERSORT-Results_LM4_mean.txt")


stage1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\stage1\\stage1_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")
stage2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\stage2\\stage2_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")


ciber_zong<-cbind(t(stage1_ciber[,1:(ncol(stage1_ciber)-3)]),(t(stage2_ciber[,1:(ncol(stage1_ciber)-3)])))
colnames(ciber_zong)<-ciber_zong[1,]
ciber_zong<-ciber_zong[-1,]
ciber_zong<-as.data.frame(ciber_zong)
ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],2,as.numeric)
str(ciber_zong)
p_zanshi<-c()
for (j in 1:nrow(ciber_zong)) {
  T_result<-wilcox.test(as.numeric(ciber_zong[j,1:nrow(stage1_ciber)]),as.numeric(ciber_zong[j,(nrow(stage1_ciber)+1):ncol(ciber_zong)]),paired = F,var.equal = T)
  
  p<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_zanshi
p_stage_yuanshi<-p_zanshi
#[1] 0.23361657        NaN 0.06673382 0.04491467 0.45374285 0.69362902 0.86320080 0.72045094 0.09366834
p_stage_yuanshi
# 
library(ggplot2)
library(ggthemes)
class<-rep(c("zaoqi","jinzhan"),c(nrow(stage1_ciber),nrow(stage2_ciber)))
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
library(dplyr)
library(ggpubr)
aaa<-ggplot(mianyi_huatu,aes(y=scores,x=factor(type),fill=class))+
  geom_boxplot(outlier.size=1)+
  theme_classic()+
  xlab("")+ylab("")+guides()+
  scale_fill_manual(values=c("#eed65d", "#5d75ee"))+stat_compare_means(method = "wilcox.test",label = "p.signif")
# geom_text(aes(x =1, y=1, label = label[1], group = type), size = 4) +
# geom_text(aes(x =2, y=1, label = label[2], group = type), size = 4)+
# geom_text(aes(x =3, y=1, label = label[3], group = type), size = 4)+
# geom_text(aes(x =4, y=1, label = label[4], group = type), size = 4)+
# geom_text(aes(x =5, y=1, label = label[5], group = type), size = 4)+
# geom_text(aes(x =6, y=1, label = label[6], group = type), size = 4)+
# geom_text(aes(x =7, y=1, label = label[7], group = type), size = 4)
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\PAAD_stage_scrna_plot.pdf",aaa,width = 10,height = 10)



# cluster1 
new_folder <- "G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\cluster1"
dir.create(new_folder)

r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/PAAD/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)


# cluster2 
new_folder <- "G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\cluster2"
dir.create(new_folder)

r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/PAAD/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)

# 
# 
# # cluster3 
# new_folder <- "G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\cluster3"
# dir.create(new_folder)
# 
# r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)
# 
# file.copy(r_files, new_folder)
# 
# source_file <- "G:/pan_result/gpu_result_tcga/PAAD/scrna/LM4_mean.txt"
# file.copy(source_file, new_folder)
# 

exp2<-cbind(expr1,expr2)
#1.2  ####
#rm(exp2)

a<-c("euclidean","maximum","manhattan","canberra","binary","minkowski")
b<-c("ward.D","ward.D2","single","complete","average","mcquitty","median", "centroid")
zuhe<-c()
for (i in 1:6) {
  zuhe1<-combn(c(a[i],b),2)[,1:8]
  zuhe<-cbind(zuhe,zuhe1)
}
?hclust
tide<-read.csv("G:\\pan_result\\tide\\PAAD_tide_result.csv",header = T)

library(survival)
library(openxlsx)
pan_pfs<-read.xlsx("E:\\450k\\TCGAPFS.xlsx")
brca_os<-read.table("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\UCSC+PAAD_OS.txt",sep = "\t",header = T)


p_shengcun_pfs<-c()
p_shengcun_os<-c()
p_tide<-c()

i=1
for (i in 1:ncol(zuhe)) {
  hc<-hclust(dist(expscale,method = zuhe[1,i]),method = zuhe[2,i])
  clusters<-cutree(hc,k=2)
  h_class<-data.frame(sample=rownames(expscale),class=clusters)

  cluster1<-h_class[which(h_class$class=="1"),1]
  cluster2<-h_class[which(h_class$class=="2"),1]
  #cluster3<-h_class[which(h_class$class=="3"),1]
  
  cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
  cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
  #cluster3_exp<-exp2[,colnames(exp2) %in%cluster3]
  
  cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
  cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
  #cluster3_exp<-data.frame(gene=rownames(exp2),cluster3_exp)
  
  #PFS
  clusters<-as.data.frame(clusters)

  clusters_tide<-data.frame(Patient=rownames(clusters),group=clusters$clusters)

  samples<-substr(gsub('[.]','-',rownames(clusters)),1,15)
  clusters2<-data.frame(sample=samples,clusters)
  #pfs
  brca_samp<-gsub('[.]','-',substr(rownames(clusters),1,12))
  brca_pfs<-pan_pfs[pan_pfs$bcr_patient_barcode%in%brca_samp,]
  clusters<-data.frame(bcr_patient_barcode=brca_samp,clusters)
  brca_pfs<-merge(brca_pfs,clusters,by="bcr_patient_barcode")
  survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)
  p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  p_shengcun_pfs<-c(p_shengcun_pfs,p1)
  
  #os
  brca_os2<-merge(brca_os,clusters2,by="sample")
  str(brca_os2)
  fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  p2 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  # p1
  p_shengcun_os<-c(p_shengcun_os,p2)
  
  #tide

  tide_final<-left_join(tide,clusters_tide,by="Patient")
  T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
  p_tide_zanshi<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_tide<-c(p_tide,p_tide_zanshi)
  
  
  print(i)
}
hcluster_result<-cbind(p_shengcun_pfs,p_shengcun_os,p_tide)


save(hcluster_result,file = "hcluster_result.Rdata")
# save(p_shengcun_os,file = "405k_p_shengcun_os_hclust.Rdata")
# save(p_shengcun_pfs,file = "405k_p_shengcun_pfs_hclust.Rdata")




which(p_ciber$V1<0.05&p_shengcun_pfs<0.05)
#[1]  1  2  9 10 17 18 20 25 26 28 30 41 42
#p_ciber_hclust<-p_ciber


library(ConsensusClusterPlus)
STAD_tumor = sweep(t(exp), 1, apply(t(exp), 1, median, na.rm = T))
a<-c("hc","pam","km")
b<-c("pearson","spearman","euclidean","binary","maximum","canberra","minkowski")
zuhe<-c()
for (i in 1:3) {
  zuhe1<-combn(c(a[i],b),2)[,1:7]
  zuhe<-cbind(zuhe,zuhe1)
}


library(survival)
library(openxlsx)
pan_pfs<-read.xlsx("E:\\450k\\TCGAPFS.xlsx")
brca_os<-read.table("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\UCSC+PAAD_OS.txt",sep = "\t",header = T)
if (is.null(STAD_tumor) || nrow(STAD_tumor) == 0) {
  stop("STAD_tumor is null or empty.")
}

if (is.null(zuhe) || ncol(zuhe) < 2 || i > ncol(zuhe)) {
  stop("Check zuhe matrix and index i.")
}
p_shengcun_pfs<-c()
p_shengcun_os<-c()
p_tide<-c()

#p_ciber_hclust<-p_ciber
i=2
for (i in 1:ncol(zuhe)) {
  results = ConsensusClusterPlus(STAD_tumor, maxK = 2, reps = 1000, pItem = 0.8, pFeature = 1, title = "untitled_consensus_cluster", clusterAlg = zuhe[1,i],
                                 distance = zuhe[2,i], seed = 123, tmyPal=NULL, writeTable=FALSE, plot = "pdf")
  clusters<-results[[2]][["consensusClass"]]
  h_class<-data.frame(sample=rownames(expscale),class=clusters)

  cluster1<-h_class[which(h_class$class=="1"),1]
  cluster2<-h_class[which(h_class$class=="2"),1]
  #cluster3<-h_class[which(h_class$class=="3"),1]
  
  cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
  cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
  #cluster3_exp<-exp2[,colnames(exp2) %in%cluster3]
  
  cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
  cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
  #cluster3_exp<-data.frame(gene=rownames(exp2),cluster3_exp)
  
  #PFS
  clusters<-as.data.frame(clusters)

  clusters_tide<-data.frame(Patient=rownames(clusters),group=clusters$clusters)

  samples<-substr(gsub('[.]','-',rownames(clusters)),1,15)
  clusters2<-data.frame(sample=samples,clusters)
  #pfs
  brca_samp<-gsub('[.]','-',substr(rownames(clusters),1,12))
  brca_pfs<-pan_pfs[pan_pfs$bcr_patient_barcode%in%brca_samp,]
  clusters<-data.frame(bcr_patient_barcode=brca_samp,clusters)
  brca_pfs<-merge(brca_pfs,clusters,by="bcr_patient_barcode")
  survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)
  p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  p_shengcun_pfs<-c(p_shengcun_pfs,p1)
  
  #os
  brca_os2<-merge(brca_os,clusters2,by="sample")
  str(brca_os2)
  fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  p2 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  # p1
  p_shengcun_os<-c(p_shengcun_os,p2)
  
  #tide

  tide_final<-left_join(tide,clusters_tide,by="Patient")
  T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
  p_tide_zanshi<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_tide<-c(p_tide,p_tide_zanshi)
  
  print(i)
  
  
}
p_ciber<-as.data.frame(p_ciber)

save(p_ciber,file="45kp_ciber_consencus.Rdata")
save(p_shengcun_os,file="45kp_p_shengcun_os_consencus.Rdata")
save(p_shengcun_pfs,file="45kp_p_shengcun_pfs_consencus.Rdata")


load(file = "cluster2\\45kp_ciber_consencus.Rdata")
zuhe[,15]#[1] "km"      "pearson"
zuhe[,16]
zuhe[,17]
zuhe[,18]
zuhe[,19]
zuhe[,20]
zuhe[,21]
 
i=1
results = ConsensusClusterPlus(STAD_tumor, maxK = 3, reps = 1000, pItem = 0.8, pFeature = 1, title = "untitled_consensus_cluster", clusterAlg = zuhe[1,i],
                               distance = zuhe[2,i], seed = 123, tmyPal=NULL, writeTable=FALSE, plot = "pdf")
clusters<-results[[2]][["consensusClass"]]
h_class<-data.frame(sample=rownames(expscale),class=clusters)
 
cluster1<-h_class[which(h_class$class=="1"),1]
cluster2<-h_class[which(h_class$class=="2"),1]
cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
write.table(cluster1_exp,"G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\cluster1\\h_cluster1_exp.txt",sep="\t",row.names = F)
write.table(cluster2_exp,"G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\cluster2\\h_cluster2_exp.txt",sep="\t",row.names = F)
setwd("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\cluster1")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster1_exp.txt",perm = 100,QN = T,"cluster1_CIBERSORT-Results_LM4_mean.txt")
setwd("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\cluster2")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster2_exp.txt",perm = 100,QN = T,"cluster2_CIBERSORT-Results_LM4_mean.txt")

cluster1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\cluster1\\cluster1_CIBERSORT-Results_LM4_mean.txt",header=T)
cluster2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\cluster2\\cluster2_CIBERSORT-Results_LM4_mean.txt",header=T)

ciber_zong<-cbind(t(cluster1_ciber[,1:(ncol(cluster1_ciber)-3)]),(t(cluster2_ciber[,1:(ncol(cluster1_ciber)-3)])))
colnames(ciber_zong)<-ciber_zong[1,]
ciber_zong<-ciber_zong[-1,]
ciber_zong<-as.data.frame(ciber_zong)
ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],1,as.numeric)
str(ciber_zong)
p_zanshi<-c()
for (j in 1:nrow(ciber_zong)) {
  T_result<-wilcox.test(as.numeric(ciber_zong[j,1:nrow(cluster1_ciber)]),as.numeric(ciber_zong[j,(nrow(cluster1_ciber)+1):ncol(ciber_zong)]),paired = F,var.equal = T)
  
  p<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}


 
i=2
hc<-hclust(dist(expscale,method = zuhe[1,i]),method = zuhe[2,i])
clusters<-cutree(hc,k=2)

h_class<-data.frame(sample=rownames(expscale),class=clusters)
table(h_class$class)
 
 
cluster1<-h_class[which(h_class$class=="1"),1]
cluster2<-h_class[which(h_class$class=="2"),1]
#cluster3<-h_class[which(h_class$class=="3"),1]

cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
#cluster3_exp<-exp2[,colnames(exp2) %in%cluster3]

cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
#cluster3_exp<-data.frame(gene=rownames(exp2),cluster3_exp)

write.table(cluster1_exp,"G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\cluster1\\h_cluster1_exp.txt",sep="\t",row.names = F)
write.table(cluster2_exp,"G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\cluster2\\h_cluster2_exp.txt",sep="\t",row.names = F)
#write.table(cluster3_exp,"G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\cluster3\\h_cluster3_exp.txt",sep="\t",row.names = F)

setwd("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\cluster1")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster1_exp.txt",perm = 100,QN = T,"cluster1_CIBERSORT-Results_LM4_mean.txt")
setwd("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\cluster2")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster2_exp.txt",perm = 100,QN = T,"cluster2_CIBERSORT-Results_LM4_mean.txt")
# setwd("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\cluster3")
# source("cibersort.R")
# CIBERSORT("LM4_mean.txt","h_cluster3_exp.txt",perm = 100,QN = T,"cluster3_CIBERSORT-Results_LM4_mean.txt")


cluster1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\cluster1\\cluster1_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")
cluster2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\cluster2\\cluster2_CIBERSORT-Results_LM4_mean.txt",header=T,sep="\t")
#cluster3_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\cluster3\\cluster3_CIBERSORT-Results_LM4_mean.txt",header=T)


ciber_zong<-cbind(t(cluster1_ciber[,1:(ncol(cluster1_ciber)-3)]),t(cluster2_ciber[,1:(ncol(cluster1_ciber)-3)]))
colnames(ciber_zong)<-ciber_zong[1,]
ciber_zong<-ciber_zong[-1,]
ciber_zong<-as.data.frame(ciber_zong)
ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],1,as.numeric)
str(ciber_zong)
p_zanshi<-c()

ciber_class<-rep(c(1,2),c(nrow(cluster1_ciber),nrow(cluster2_ciber)))
for (j in 1:nrow(ciber_zong)) {
  T_result<-wilcox.test(as.numeric(ciber_zong[j,])~ciber_class)
  
  p<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_cluster_new<-p_zanshi
p_cluster_new

#THCA####

setwd("G:\\pan_result\\gpu_result_tcga\\THCA")
gene<-read.csv("448featuresweightme.csv",row.names = 1)
gene<-colnames(gene)
#exp<-read.csv("sigmoid.csv")
exp<-read.csv("448featuresweightme.csv",row.names = 1)

label<-read.csv("label.csv")
exp11<-exp[rownames(exp)%in%(label$X[which(label$label=="1")]),]
exp22<-exp[rownames(exp)%in%(label$X[which(label$label=="2")]),]
exp<-rbind(exp11,exp22)
# rownames(exp)<-exp$X
# exp<-exp[,-1]
expscale<-scale(exp)

library(factoextra)
setwd("G:\\pan_result\\gpu_result_tcga\\THCA\\450k")
set.seed(123)
fviz_nbclust<-fviz_nbclust(expscale, kmeans, k.max = 15)#2
fviz_nbclust
ggsave("THCA_fviz_nbclust.pdf",fviz_nbclust,width = 10,height = 10)
 
# ggsave("THCA_kmeans_plot.pdf",kmeans_plot,width = 10,height = 10)


# expr<-read.table("tcga_brca_tpm.txt")
# expr<-t(expr)
# expr<-as.data.frame(expr)
# expr$name<-substr(rownames(expr),1,15)
# library(limma)
# expr<-avereps(expr,ID=expr$name)
# ?avereps
# expr<-as.data.frame(expr)
# rownames(expr)<-expr$name
# expr<-expr[,-ncol(expr)]
# expr<-as.data.frame(t(expr))
# write.table(expr,"brca_tpm_quchong.txt")
expr1<-read.csv("G:\\dnameth\\THCA\\exp1.csv",row.names = 1)
expr2<-read.csv("G:\\dnameth\\THCA\\exp2.csv",row.names = 1)
expr<-cbind(expr1,expr2)
# expr1<-expr[,colnames(expr)%in%exp11$X]
# expr2<-expr[,colnames(expr)%in%exp22$X]
str(expr1)
# expr1<-apply(expr1, 2, as.numeric)
# rownames(expr1)<-rownames(expr2)
# expr2<-apply(expr2, 2, as.numeric)
# rownames(expr2)<-rownames(expr1)


library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(ggsignif)


exp_tide<-t(apply(expr,1,function(x){x-(mean(x))}))
write.table(data.frame(ID=rownames(exp_tide),exp_tide,check.names=F),'THCA_tide.txt',sep = "\t",quote=F,row.names = FALSE)


getwd()


write.table(clusters,"450k\\THCA_cluster.txt")




#estimate####
library(estimate)
 
dat = log2(edgeR::cpm(expr)+1)
write.table(dat, file = "THCA_estimate_input.txt", sep = '\t', quote = F)
exp.file = "THCA_estimate_input.txt"
in.gct.file = "ESTIMATE_THCA_input.gct"
 
outputGCT(exp.file, in.gct.file)
rt <- read.table("ESTIMATE_THCA_input.gct", 
                 skip = 2, 
                 header = TRUE, 
                 sep = "\t")

 
filterCommonGenes(input.f = exp.file,  
                  output.f = "ESTIMATE_THCA_gene.gct",  
                  id = "GeneSymbol"  
) 

#2esitimate

out.score.file = "ESTIMATE_score.gct"
 
estimateScore(in.gct.file,  
              out.score.file, 
              platform = "illumina")  


ESTIMATE_score = read.table(out.score.file,
                            skip = 2,
                            header = T,
                            row.names = 1)
ESTIMATE_score = as.data.frame(t(ESTIMATE_score[, 2:ncol(ESTIMATE_score)]))  
ESTIMATE_score[1:6, 1:3]
write.table(ESTIMATE_score, file = "THCA_estimate_score.txt", sep = '\t', quote = F)




new_folder <- "G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\stage1"
dir.create(new_folder)

source_folder <- "G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\stage1"
r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/THCA/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)
write.table(expr1,"stage1\\THCA_tpm_stage1.txt",sep = "\t")



new_folder <- "G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\stage2"
dir.create(new_folder)

source_folder <- "G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\stage2"
r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/THCA/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)
write.table(expr2,"stage2\\THCA_tpm_stage2.txt",sep = "\t")




setwd("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\stage1")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","THCA_tpm_stage1.txt",perm = 100,QN = T,"stage1_CIBERSORT-Results_LM4_mean.txt")

setwd("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\stage2")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","THCA_tpm_stage2.txt",perm = 100,QN = T,"stage2_CIBERSORT-Results_LM4_mean.txt")


stage1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\stage1\\stage1_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")
stage2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\stage2\\stage2_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")


ciber_zong<-cbind(t(stage1_ciber[,1:(ncol(stage1_ciber)-3)]),(t(stage2_ciber[,1:(ncol(stage1_ciber)-3)])))
colnames(ciber_zong)<-ciber_zong[1,]
ciber_zong<-ciber_zong[-1,]
ciber_zong<-as.data.frame(ciber_zong)
ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],2,as.numeric)
str(ciber_zong)
p_zanshi<-c()
for (j in 1:nrow(ciber_zong)) {
  T_result<-wilcox.test(as.numeric(ciber_zong[j,1:nrow(stage1_ciber)]),as.numeric(ciber_zong[j,(nrow(stage1_ciber)+1):ncol(ciber_zong)]),paired = F,var.equal = T)
  
  p<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_zanshi
p_stage_yuanshi<-p_zanshi
#[1] 0.23361657        NaN 0.06673382 0.04491467 0.45374285 0.69362902 0.86320080 0.72045094 0.09366834
p_stage_yuanshi
# 
library(ggplot2)
library(ggthemes)
class<-rep(c("zaoqi","jinzhan"),c(nrow(stage1_ciber),nrow(stage2_ciber)))
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
library(dplyr)
library(ggpubr)
aaa<-ggplot(mianyi_huatu,aes(y=scores,x=factor(type),fill=class))+
  geom_boxplot(outlier.size=1)+
  theme_classic()+
  xlab("")+ylab("")+guides()+
  scale_fill_manual(values=c("#eed65d", "#5d75ee"))+stat_compare_means(method = "wilcox.test",label = "p.signif")
# geom_text(aes(x =1, y=1, label = label[1], group = type), size = 4) +
# geom_text(aes(x =2, y=1, label = label[2], group = type), size = 4)+
# geom_text(aes(x =3, y=1, label = label[3], group = type), size = 4)+
# geom_text(aes(x =4, y=1, label = label[4], group = type), size = 4)+
# geom_text(aes(x =5, y=1, label = label[5], group = type), size = 4)+
# geom_text(aes(x =6, y=1, label = label[6], group = type), size = 4)+
# geom_text(aes(x =7, y=1, label = label[7], group = type), size = 4)
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\THCA_stage_scrna_plot.pdf",aaa,width = 10,height = 10)



# cluster1 
new_folder <- "G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\cluster1"
dir.create(new_folder)

r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/THCA/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)


# cluster2 
new_folder <- "G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\cluster2"
dir.create(new_folder)

r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/THCA/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)

# 
# 
# # cluster3 
# new_folder <- "G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\cluster3"
# dir.create(new_folder)
# 
# r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)
# 
# file.copy(r_files, new_folder)
# 
# source_file <- "G:/pan_result/gpu_result_tcga/THCA/scrna/LM4_mean.txt"
# file.copy(source_file, new_folder)
# 

exp2<-cbind(expr1,expr2)
#1.2  ####
#rm(exp2)

a<-c("euclidean","maximum","manhattan","canberra","binary","minkowski")
b<-c("ward.D","ward.D2","single","complete","average","mcquitty","median", "centroid")
zuhe<-c()
for (i in 1:6) {
  zuhe1<-combn(c(a[i],b),2)[,1:8]
  zuhe<-cbind(zuhe,zuhe1)
}
?hclust
tide<-read.csv("G:\\pan_result\\tide\\THCA_tide_result.csv",header = T)

library(survival)
library(openxlsx)
pan_pfs<-read.xlsx("E:\\450k\\TCGAPFS.xlsx")
brca_os<-read.table("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\UCSC+THCA_OS.txt",sep = "\t",header = T)


p_shengcun_pfs<-c()
p_shengcun_os<-c()
p_tide<-c()

i=1
for (i in 1:ncol(zuhe)) {
  hc<-hclust(dist(expscale,method = zuhe[1,i]),method = zuhe[2,i])
  clusters<-cutree(hc,k=2)
  h_class<-data.frame(sample=rownames(expscale),class=clusters)

  cluster1<-h_class[which(h_class$class=="1"),1]
  cluster2<-h_class[which(h_class$class=="2"),1]
  #cluster3<-h_class[which(h_class$class=="3"),1]
  
  cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
  cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
  #cluster3_exp<-exp2[,colnames(exp2) %in%cluster3]
  
  cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
  cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
  #cluster3_exp<-data.frame(gene=rownames(exp2),cluster3_exp)
  
  #PFS
  clusters<-as.data.frame(clusters)

  clusters_tide<-data.frame(Patient=rownames(clusters),group=clusters$clusters)

  samples<-substr(gsub('[.]','-',rownames(clusters)),1,15)
  clusters2<-data.frame(sample=samples,clusters)
  #pfs
  brca_samp<-gsub('[.]','-',substr(rownames(clusters),1,12))
  brca_pfs<-pan_pfs[pan_pfs$bcr_patient_barcode%in%brca_samp,]
  clusters<-data.frame(bcr_patient_barcode=brca_samp,clusters)
  brca_pfs<-merge(brca_pfs,clusters,by="bcr_patient_barcode")
  survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)
  p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  p_shengcun_pfs<-c(p_shengcun_pfs,p1)
  
  #os
  brca_os2<-merge(brca_os,clusters2,by="sample")
  str(brca_os2)
  fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  p2 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  # p1
  p_shengcun_os<-c(p_shengcun_os,p2)
  
  #tide

  tide_final<-left_join(tide,clusters_tide,by="Patient")
  T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
  p_tide_zanshi<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_tide<-c(p_tide,p_tide_zanshi)
  
  
  print(i)
}
hcluster_result<-cbind(p_shengcun_pfs,p_shengcun_os,p_tide)


save(hcluster_result,file = "hcluster_result.Rdata")
# save(p_shengcun_os,file = "405k_p_shengcun_os_hclust.Rdata")
# save(p_shengcun_pfs,file = "405k_p_shengcun_pfs_hclust.Rdata")




which(p_ciber$V1<0.05&p_shengcun_pfs<0.05)
#[1]  1  2  9 10 17 18 20 25 26 28 30 41 42
#p_ciber_hclust<-p_ciber


library(ConsensusClusterPlus)
STAD_tumor = sweep(t(exp), 1, apply(t(exp), 1, median, na.rm = T))
a<-c("hc","pam","km")
b<-c("pearson","spearman","euclidean","binary","maximum","canberra","minkowski")
zuhe<-c()
for (i in 1:3) {
  zuhe1<-combn(c(a[i],b),2)[,1:7]
  zuhe<-cbind(zuhe,zuhe1)
}


library(survival)
library(openxlsx)
pan_pfs<-read.xlsx("E:\\450k\\TCGAPFS.xlsx")
brca_os<-read.table("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\UCSC+THCA_OS.txt",sep = "\t",header = T)
if (is.null(STAD_tumor) || nrow(STAD_tumor) == 0) {
  stop("STAD_tumor is null or empty.")
}

if (is.null(zuhe) || ncol(zuhe) < 2 || i > ncol(zuhe)) {
  stop("Check zuhe matrix and index i.")
}
p_shengcun_pfs<-c()
p_shengcun_os<-c()
p_tide<-c()

#p_ciber_hclust<-p_ciber
i=2
for (i in 1:ncol(zuhe)) {
  results = ConsensusClusterPlus(STAD_tumor, maxK = 2, reps = 1000, pItem = 0.8, pFeature = 1, title = "untitled_consensus_cluster", clusterAlg = zuhe[1,i],
                                 distance = zuhe[2,i], seed = 123, tmyPal=NULL, writeTable=FALSE, plot = "pdf")
  clusters<-results[[2]][["consensusClass"]]
  h_class<-data.frame(sample=rownames(expscale),class=clusters)

  cluster1<-h_class[which(h_class$class=="1"),1]
  cluster2<-h_class[which(h_class$class=="2"),1]
  #cluster3<-h_class[which(h_class$class=="3"),1]
  
  cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
  cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
  #cluster3_exp<-exp2[,colnames(exp2) %in%cluster3]
  
  cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
  cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
  #cluster3_exp<-data.frame(gene=rownames(exp2),cluster3_exp)
  
  #PFS
  clusters<-as.data.frame(clusters)

  clusters_tide<-data.frame(Patient=rownames(clusters),group=clusters$clusters)

  samples<-substr(gsub('[.]','-',rownames(clusters)),1,15)
  clusters2<-data.frame(sample=samples,clusters)
  #pfs
  brca_samp<-gsub('[.]','-',substr(rownames(clusters),1,12))
  brca_pfs<-pan_pfs[pan_pfs$bcr_patient_barcode%in%brca_samp,]
  clusters<-data.frame(bcr_patient_barcode=brca_samp,clusters)
  brca_pfs<-merge(brca_pfs,clusters,by="bcr_patient_barcode")
  survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)
  p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  p_shengcun_pfs<-c(p_shengcun_pfs,p1)
  
  #os
  brca_os2<-merge(brca_os,clusters2,by="sample")
  str(brca_os2)
  fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  p2 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  # p1
  p_shengcun_os<-c(p_shengcun_os,p2)
  
  #tide

  tide_final<-left_join(tide,clusters_tide,by="Patient")
  T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
  p_tide_zanshi<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_tide<-c(p_tide,p_tide_zanshi)
  
  print(i)
  
  
}
p_ciber<-as.data.frame(p_ciber)

save(p_ciber,file="45kp_ciber_consencus.Rdata")
save(p_shengcun_os,file="45kp_p_shengcun_os_consencus.Rdata")
save(p_shengcun_pfs,file="45kp_p_shengcun_pfs_consencus.Rdata")


load(file = "cluster2\\45kp_ciber_consencus.Rdata")
zuhe[,15]#[1] "km"      "pearson"
zuhe[,16]
zuhe[,17]
zuhe[,18]
zuhe[,19]
zuhe[,20]
zuhe[,21]
 
i=1
results = ConsensusClusterPlus(STAD_tumor, maxK = 3, reps = 1000, pItem = 0.8, pFeature = 1, title = "untitled_consensus_cluster", clusterAlg = zuhe[1,i],
                               distance = zuhe[2,i], seed = 123, tmyPal=NULL, writeTable=FALSE, plot = "pdf")
clusters<-results[[2]][["consensusClass"]]
h_class<-data.frame(sample=rownames(expscale),class=clusters)
 
cluster1<-h_class[which(h_class$class=="1"),1]
cluster2<-h_class[which(h_class$class=="2"),1]
cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
write.table(cluster1_exp,"G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\cluster1\\h_cluster1_exp.txt",sep="\t",row.names = F)
write.table(cluster2_exp,"G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\cluster2\\h_cluster2_exp.txt",sep="\t",row.names = F)
setwd("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\cluster1")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster1_exp.txt",perm = 100,QN = T,"cluster1_CIBERSORT-Results_LM4_mean.txt")
setwd("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\cluster2")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster2_exp.txt",perm = 100,QN = T,"cluster2_CIBERSORT-Results_LM4_mean.txt")

cluster1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\cluster1\\cluster1_CIBERSORT-Results_LM4_mean.txt",header=T)
cluster2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\cluster2\\cluster2_CIBERSORT-Results_LM4_mean.txt",header=T)

ciber_zong<-cbind(t(cluster1_ciber[,1:(ncol(cluster1_ciber)-3)]),(t(cluster2_ciber[,1:(ncol(cluster1_ciber)-3)])))
colnames(ciber_zong)<-ciber_zong[1,]
ciber_zong<-ciber_zong[-1,]
ciber_zong<-as.data.frame(ciber_zong)
ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],1,as.numeric)
str(ciber_zong)
p_zanshi<-c()
for (j in 1:nrow(ciber_zong)) {
  T_result<-wilcox.test(as.numeric(ciber_zong[j,1:nrow(cluster1_ciber)]),as.numeric(ciber_zong[j,(nrow(cluster1_ciber)+1):ncol(ciber_zong)]),paired = F,var.equal = T)
  
  p<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}


 
i=10
hc<-hclust(dist(expscale,method = zuhe[1,i]),method = zuhe[2,i])
clusters<-cutree(hc,k=2)

h_class<-data.frame(sample=rownames(expscale),class=clusters)
table(h_class$class)
 
 
cluster1<-h_class[which(h_class$class=="1"),1]
cluster2<-h_class[which(h_class$class=="2"),1]
#cluster3<-h_class[which(h_class$class=="3"),1]

cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
#cluster3_exp<-exp2[,colnames(exp2) %in%cluster3]

cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
#cluster3_exp<-data.frame(gene=rownames(exp2),cluster3_exp)

write.table(cluster1_exp,"G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\cluster1\\h_cluster1_exp.txt",sep="\t",row.names = F)
write.table(cluster2_exp,"G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\cluster2\\h_cluster2_exp.txt",sep="\t",row.names = F)
#write.table(cluster3_exp,"G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\cluster3\\h_cluster3_exp.txt",sep="\t",row.names = F)

setwd("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\cluster1")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster1_exp.txt",perm = 100,QN = T,"cluster1_CIBERSORT-Results_LM4_mean.txt")
setwd("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\cluster2")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster2_exp.txt",perm = 100,QN = T,"cluster2_CIBERSORT-Results_LM4_mean.txt")
# setwd("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\cluster3")
# source("cibersort.R")
# CIBERSORT("LM4_mean.txt","h_cluster3_exp.txt",perm = 100,QN = T,"cluster3_CIBERSORT-Results_LM4_mean.txt")


cluster1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\cluster1\\cluster1_CIBERSORT-Results_LM4_mean.txt",header=T,sep="\t")
cluster2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\cluster2\\cluster2_CIBERSORT-Results_LM4_mean.txt",header=T,sep="\t")
#cluster3_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\cluster3\\cluster3_CIBERSORT-Results_LM4_mean.txt",header=T)


ciber_zong<-cbind(t(cluster1_ciber[,1:(ncol(cluster1_ciber)-3)]),t(cluster2_ciber[,1:(ncol(cluster1_ciber)-3)]))
colnames(ciber_zong)<-ciber_zong[1,]
ciber_zong<-ciber_zong[-1,]
ciber_zong<-as.data.frame(ciber_zong)
ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],1,as.numeric)
str(ciber_zong)
p_zanshi<-c()

ciber_class<-rep(c(1,2),c(nrow(cluster1_ciber),nrow(cluster2_ciber)))
for (j in 1:nrow(ciber_zong)) {
  T_result<-wilcox.test(as.numeric(ciber_zong[j,])~ciber_class)
  
  p<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_cluster_new<-p_zanshi
p_cluster_new

#UCEC####

setwd("G:\\pan_result\\gpu_result_tcga\\UCEC")
gene<-read.csv("541featuresweightme.csv",row.names = 1)
gene<-colnames(gene)
#exp<-read.csv("sigmoid.csv")
exp<-read.csv("541featuresweightme.csv",row.names = 1)

label<-read.csv("label.csv")
exp11<-exp[rownames(exp)%in%(label$X[which(label$label=="1")]),]
exp22<-exp[rownames(exp)%in%(label$X[which(label$label=="2")]),]
exp<-rbind(exp11,exp22)
# rownames(exp)<-exp$X
# exp<-exp[,-1]
expscale<-scale(exp)

library(factoextra)
setwd("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k")
set.seed(123)
fviz_nbclust<-fviz_nbclust(expscale, kmeans, k.max = 15)#2
fviz_nbclust
ggsave("UCEC_fviz_nbclust.pdf",fviz_nbclust,width = 10,height = 10)
 
# ggsave("UCEC_kmeans_plot.pdf",kmeans_plot,width = 10,height = 10)


# expr<-read.table("tcga_brca_tpm.txt")
# expr<-t(expr)
# expr<-as.data.frame(expr)
# expr$name<-substr(rownames(expr),1,15)
# library(limma)
# expr<-avereps(expr,ID=expr$name)
# ?avereps
# expr<-as.data.frame(expr)
# rownames(expr)<-expr$name
# expr<-expr[,-ncol(expr)]
# expr<-as.data.frame(t(expr))
# write.table(expr,"brca_tpm_quchong.txt")
expr1<-read.csv("G:\\dnameth\\UCEC\\exp1.csv",row.names = 1)
expr2<-read.csv("G:\\dnameth\\UCEC\\exp2.csv",row.names = 1)
expr<-cbind(expr1,expr2)
# expr1<-expr[,colnames(expr)%in%exp11$X]
# expr2<-expr[,colnames(expr)%in%exp22$X]
str(expr1)
# expr1<-apply(expr1, 2, as.numeric)
# rownames(expr1)<-rownames(expr2)
# expr2<-apply(expr2, 2, as.numeric)
# rownames(expr2)<-rownames(expr1)


library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(ggsignif)


exp_tide<-t(apply(expr,1,function(x){x-(mean(x))}))
write.table(data.frame(ID=rownames(exp_tide),exp_tide,check.names=F),'UCEC_tide.txt',sep = "\t",quote=F,row.names = FALSE)


getwd()



new_folder <- "G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\stage1"
dir.create(new_folder)

source_folder <- "G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\stage1"
r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/UCEC/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)
write.table(expr1,"stage1\\UCEC_tpm_stage1.txt",sep = "\t")



new_folder <- "G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\stage2"
dir.create(new_folder)

source_folder <- "G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\stage2"
r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/COAD/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)
write.table(expr2,"stage2\\UCEC_tpm_stage2.txt",sep = "\t")


 
library(estimate)
 
dat = log2(edgeR::cpm(expr)+1)
write.table(dat, file = "UCEC_estimate_input.txt", sep = '\t', quote = F)
exp.file = "UCEC_estimate_input.txt"
in.gct.file = "ESTIMATE_UCEC_input.gct"
 
outputGCT(exp.file, in.gct.file)
rt <- read.table("ESTIMATE_UCEC_input.gct", 
                 skip = 2,  
                 header = TRUE, 
                 sep = "\t")


filterCommonGenes(input.f = exp.file,  
                  output.f = "ESTIMATE_UCEC_gene.gct",  
                  id = "GeneSymbol"  
) 

#2esitimate 

out.score.file = "ESTIMATE_score.gct"

estimateScore(in.gct.file,  
              out.score.file, 
              platform = "illumina")  

 
ESTIMATE_score = read.table(out.score.file,
                            skip = 2,
                            header = T,
                            row.names = 1)
ESTIMATE_score = as.data.frame(t(ESTIMATE_score[, 2:ncol(ESTIMATE_score)])) 
ESTIMATE_score[1:6, 1:3]
write.table(ESTIMATE_score, file = "UCEC_estimate_score.txt", sep = '\t', quote = F)





setwd("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\stage1")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","UCEC_tpm_stage1.txt",perm = 100,QN = T,"stage1_CIBERSORT-Results_LM4_mean.txt")

setwd("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\stage2")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","UCEC_tpm_stage2.txt",perm = 100,QN = T,"stage2_CIBERSORT-Results_LM4_mean.txt")


stage1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\stage1\\stage1_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")
stage2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\stage2\\stage2_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")


ciber_zong<-cbind(t(stage1_ciber[,1:(ncol(stage1_ciber)-3)]),(t(stage2_ciber[,1:(ncol(stage1_ciber)-3)])))
colnames(ciber_zong)<-ciber_zong[1,]
ciber_zong<-ciber_zong[-1,]
ciber_zong<-as.data.frame(ciber_zong)
ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],2,as.numeric)
str(ciber_zong)
p_zanshi<-c()
for (j in 1:nrow(ciber_zong)) {
  T_result<-wilcox.test(as.numeric(ciber_zong[j,1:nrow(stage1_ciber)]),as.numeric(ciber_zong[j,(nrow(stage1_ciber)+1):ncol(ciber_zong)]),paired = F,var.equal = T)
  
  p<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_zanshi
p_stage_yuanshi<-p_zanshi
#[1] 0.23361657        NaN 0.06673382 0.04491467 0.45374285 0.69362902 0.86320080 0.72045094 0.09366834
p_stage_yuanshi
# 
library(ggplot2)
library(ggthemes)
class<-rep(c("zaoqi","jinzhan"),c(nrow(stage1_ciber),nrow(stage2_ciber)))
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
library(dplyr)
library(ggpubr)
aaa<-ggplot(mianyi_huatu,aes(y=scores,x=factor(type),fill=class))+
  geom_boxplot(outlier.size=1)+
  theme_classic()+
  xlab("")+ylab("")+guides()+
  scale_fill_manual(values=c("#eed65d", "#5d75ee"))+stat_compare_means(method = "wilcox.test",label = "p.signif")
# geom_text(aes(x =1, y=1, label = label[1], group = type), size = 4) +
# geom_text(aes(x =2, y=1, label = label[2], group = type), size = 4)+
# geom_text(aes(x =3, y=1, label = label[3], group = type), size = 4)+
# geom_text(aes(x =4, y=1, label = label[4], group = type), size = 4)+
# geom_text(aes(x =5, y=1, label = label[5], group = type), size = 4)+
# geom_text(aes(x =6, y=1, label = label[6], group = type), size = 4)+
# geom_text(aes(x =7, y=1, label = label[7], group = type), size = 4)
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\UCEC_stage_scrna_plot.pdf",aaa,width = 10,height = 10)



# cluster1 
new_folder <- "G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\cluster1"
dir.create(new_folder)

r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/UCEC/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)


# cluster2 
new_folder <- "G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\cluster2"
dir.create(new_folder)

r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)

file.copy(r_files, new_folder)

source_file <- "G:/pan_result/gpu_result_tcga/UCEC/scrna/LM4_mean.txt"
file.copy(source_file, new_folder)

# 
# 
# # cluster3 
# new_folder <- "G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\cluster3"
# dir.create(new_folder)
# 
# r_files <- list.files("G:/pan_result/gpu_result_tcga/ciber", pattern = "\\.R$", full.names = TRUE)
# 
# file.copy(r_files, new_folder)
# 
# source_file <- "G:/pan_result/gpu_result_tcga/UCEC/scrna/LM4_mean.txt"
# file.copy(source_file, new_folder)
# 

exp2<-cbind(expr1,expr2)
#1.2  ####
#rm(exp2)

a<-c("euclidean","maximum","manhattan","canberra","binary","minkowski")
b<-c("ward.D","ward.D2","single","complete","average","mcquitty","median", "centroid")
zuhe<-c()
for (i in 1:6) {
  zuhe1<-combn(c(a[i],b),2)[,1:8]
  zuhe<-cbind(zuhe,zuhe1)
}
?hclust
tide<-read.csv("G:\\pan_result\\tide\\UCEC_tide_result.csv",header = T)

library(survival)
library(openxlsx)
pan_pfs<-read.xlsx("E:\\450k\\TCGAPFS.xlsx")
brca_os<-read.table("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\UCSC+UCEC_OS.txt",sep = "\t",header = T)


p_shengcun_pfs<-c()
p_shengcun_os<-c()
p_tide<-c()

i=1
for (i in 1:ncol(zuhe)) {
  hc<-hclust(dist(expscale,method = zuhe[1,i]),method = zuhe[2,i])
  clusters<-cutree(hc,k=2)
  h_class<-data.frame(sample=rownames(expscale),class=clusters)

  cluster1<-h_class[which(h_class$class=="1"),1]
  cluster2<-h_class[which(h_class$class=="2"),1]
  #cluster3<-h_class[which(h_class$class=="3"),1]
  
  cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
  cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
  #cluster3_exp<-exp2[,colnames(exp2) %in%cluster3]
  
  cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
  cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
  #cluster3_exp<-data.frame(gene=rownames(exp2),cluster3_exp)
  
  #PFS
  clusters<-as.data.frame(clusters)

  clusters_tide<-data.frame(Patient=rownames(clusters),group=clusters$clusters)

  samples<-substr(gsub('[.]','-',rownames(clusters)),1,15)
  clusters2<-data.frame(sample=samples,clusters)
  #pfs
  brca_samp<-gsub('[.]','-',substr(rownames(clusters),1,12))
  brca_pfs<-pan_pfs[pan_pfs$bcr_patient_barcode%in%brca_samp,]
  clusters<-data.frame(bcr_patient_barcode=brca_samp,clusters)
  brca_pfs<-merge(brca_pfs,clusters,by="bcr_patient_barcode")
  survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)
  p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  p_shengcun_pfs<-c(p_shengcun_pfs,p1)
  
  #os
  brca_os2<-merge(brca_os,clusters2,by="sample")
  str(brca_os2)
  fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  p2 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  # p1
  p_shengcun_os<-c(p_shengcun_os,p2)
  
  #tide

  tide_final<-left_join(tide,clusters_tide,by="Patient")
  T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
  p_tide_zanshi<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_tide<-c(p_tide,p_tide_zanshi)
  
  
  print(i)
}
hcluster_result<-cbind(p_shengcun_pfs,p_shengcun_os,p_tide)


save(hcluster_result,file = "hcluster_result.Rdata")
# save(p_shengcun_os,file = "405k_p_shengcun_os_hclust.Rdata")
# save(p_shengcun_pfs,file = "405k_p_shengcun_pfs_hclust.Rdata")




which(p_ciber$V1<0.05&p_shengcun_pfs<0.05)
#[1]  1  2  9 10 17 18 20 25 26 28 30 41 42
#p_ciber_hclust<-p_ciber


library(ConsensusClusterPlus)
STAD_tumor = sweep(t(exp), 1, apply(t(exp), 1, median, na.rm = T))
a<-c("hc","pam","km")
b<-c("pearson","spearman","euclidean","binary","maximum","canberra","minkowski")
zuhe<-c()
for (i in 1:3) {
  zuhe1<-combn(c(a[i],b),2)[,1:7]
  zuhe<-cbind(zuhe,zuhe1)
}


library(survival)
library(openxlsx)
pan_pfs<-read.xlsx("E:\\450k\\TCGAPFS.xlsx")
brca_os<-read.table("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\UCSC+UCEC_OS.txt",sep = "\t",header = T)
if (is.null(STAD_tumor) || nrow(STAD_tumor) == 0) {
  stop("STAD_tumor is null or empty.")
}

if (is.null(zuhe) || ncol(zuhe) < 2 || i > ncol(zuhe)) {
  stop("Check zuhe matrix and index i.")
}
p_shengcun_pfs<-c()
p_shengcun_os<-c()
p_tide<-c()

#p_ciber_hclust<-p_ciber
i=2
for (i in 1:ncol(zuhe)) {
  results = ConsensusClusterPlus(STAD_tumor, maxK = 2, reps = 1000, pItem = 0.8, pFeature = 1, title = "untitled_consensus_cluster", clusterAlg = zuhe[1,i],
                                 distance = zuhe[2,i], seed = 123, tmyPal=NULL, writeTable=FALSE, plot = "pdf")
  clusters<-results[[2]][["consensusClass"]]
  h_class<-data.frame(sample=rownames(expscale),class=clusters)

  cluster1<-h_class[which(h_class$class=="1"),1]
  cluster2<-h_class[which(h_class$class=="2"),1]
  #cluster3<-h_class[which(h_class$class=="3"),1]
  
  cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
  cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
  #cluster3_exp<-exp2[,colnames(exp2) %in%cluster3]
  
  cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
  cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
  #cluster3_exp<-data.frame(gene=rownames(exp2),cluster3_exp)
  
  #PFS
  clusters<-as.data.frame(clusters)

  clusters_tide<-data.frame(Patient=rownames(clusters),group=clusters$clusters)

  samples<-substr(gsub('[.]','-',rownames(clusters)),1,15)
  clusters2<-data.frame(sample=samples,clusters)
  #pfs
  brca_samp<-gsub('[.]','-',substr(rownames(clusters),1,12))
  brca_pfs<-pan_pfs[pan_pfs$bcr_patient_barcode%in%brca_samp,]
  clusters<-data.frame(bcr_patient_barcode=brca_samp,clusters)
  brca_pfs<-merge(brca_pfs,clusters,by="bcr_patient_barcode")
  survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)
  p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  p_shengcun_pfs<-c(p_shengcun_pfs,p1)
  
  #os
  brca_os2<-merge(brca_os,clusters2,by="sample")
  str(brca_os2)
  fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
  p2 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
  # p1
  p_shengcun_os<-c(p_shengcun_os,p2)
  
  #tide

  tide_final<-left_join(tide,clusters_tide,by="Patient")
  T_result<-wilcox.test(as.numeric(tide_final$TIDE)~tide_final$group)
  p_tide_zanshi<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_tide<-c(p_tide,p_tide_zanshi)
  
  print(i)
  
  
}
p_ciber<-as.data.frame(p_ciber)

save(p_ciber,file="45kp_ciber_consencus.Rdata")
save(p_shengcun_os,file="45kp_p_shengcun_os_consencus.Rdata")
save(p_shengcun_pfs,file="45kp_p_shengcun_pfs_consencus.Rdata")


load(file = "cluster2\\45kp_ciber_consencus.Rdata")
zuhe[,15]#[1] "km"      "pearson"
zuhe[,16]
zuhe[,17]
zuhe[,18]
zuhe[,19]
zuhe[,20]
zuhe[,21]
 
i=1
results = ConsensusClusterPlus(STAD_tumor, maxK = 3, reps = 1000, pItem = 0.8, pFeature = 1, title = "untitled_consensus_cluster", clusterAlg = zuhe[1,i],
                               distance = zuhe[2,i], seed = 123, tmyPal=NULL, writeTable=FALSE, plot = "pdf")
clusters<-results[[2]][["consensusClass"]]
h_class<-data.frame(sample=rownames(expscale),class=clusters)
 
cluster1<-h_class[which(h_class$class=="1"),1]
cluster2<-h_class[which(h_class$class=="2"),1]
cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
write.table(cluster1_exp,"G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\cluster1\\h_cluster1_exp.txt",sep="\t",row.names = F)
write.table(cluster2_exp,"G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\cluster2\\h_cluster2_exp.txt",sep="\t",row.names = F)
setwd("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\cluster1")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster1_exp.txt",perm = 100,QN = T,"cluster1_CIBERSORT-Results_LM4_mean.txt")
setwd("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\cluster2")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster2_exp.txt",perm = 100,QN = T,"cluster2_CIBERSORT-Results_LM4_mean.txt")

cluster1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\cluster1\\cluster1_CIBERSORT-Results_LM4_mean.txt",header=T)
cluster2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\cluster2\\cluster2_CIBERSORT-Results_LM4_mean.txt",header=T)

ciber_zong<-cbind(t(cluster1_ciber[,1:(ncol(cluster1_ciber)-3)]),(t(cluster2_ciber[,1:(ncol(cluster1_ciber)-3)])))
colnames(ciber_zong)<-ciber_zong[1,]
ciber_zong<-ciber_zong[-1,]
ciber_zong<-as.data.frame(ciber_zong)
ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],1,as.numeric)
str(ciber_zong)
p_zanshi<-c()
for (j in 1:nrow(ciber_zong)) {
  T_result<-wilcox.test(as.numeric(ciber_zong[j,1:nrow(cluster1_ciber)]),as.numeric(ciber_zong[j,(nrow(cluster1_ciber)+1):ncol(ciber_zong)]),paired = F,var.equal = T)
  
  p<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}


 
i=30
hc<-hclust(dist(expscale,method = zuhe[1,i]),method = zuhe[2,i])
clusters<-cutree(hc,k=2)

h_class<-data.frame(sample=rownames(expscale),class=clusters)
table(h_class$class)
 
 
cluster1<-h_class[which(h_class$class=="1"),1]
cluster2<-h_class[which(h_class$class=="2"),1]
#cluster3<-h_class[which(h_class$class=="3"),1]

cluster1_exp<-exp2[,colnames(exp2) %in%cluster1]
cluster2_exp<-exp2[,colnames(exp2) %in%cluster2]
#cluster3_exp<-exp2[,colnames(exp2) %in%cluster3]

cluster1_exp<-data.frame(gene=rownames(exp2),cluster1_exp)
cluster2_exp<-data.frame(gene=rownames(exp2),cluster2_exp)
#cluster3_exp<-data.frame(gene=rownames(exp2),cluster3_exp)

write.table(cluster1_exp,"G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\cluster1\\h_cluster1_exp.txt",sep="\t",row.names = F)
write.table(cluster2_exp,"G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\cluster2\\h_cluster2_exp.txt",sep="\t",row.names = F)
#write.table(cluster3_exp,"G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\cluster3\\h_cluster3_exp.txt",sep="\t",row.names = F)

setwd("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\cluster1")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster1_exp.txt",perm = 100,QN = T,"cluster1_CIBERSORT-Results_LM4_mean.txt")
setwd("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\cluster2")
source("cibersort.R")
CIBERSORT("LM4_mean.txt","h_cluster2_exp.txt",perm = 100,QN = T,"cluster2_CIBERSORT-Results_LM4_mean.txt")
# setwd("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\cluster3")
# source("cibersort.R")
# CIBERSORT("LM4_mean.txt","h_cluster3_exp.txt",perm = 100,QN = T,"cluster3_CIBERSORT-Results_LM4_mean.txt")


cluster1_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\cluster1\\cluster1_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")
cluster2_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\cluster2\\cluster2_CIBERSORT-Results_LM4_mean.txt",header=T,sep = "\t")
#cluster3_ciber<-read.table("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\cluster3\\cluster3_CIBERSORT-Results_LM4_mean.txt",header=T)


ciber_zong<-cbind(t(cluster1_ciber[,1:(ncol(cluster1_ciber)-3)]),t(cluster2_ciber[,1:(ncol(cluster1_ciber)-3)]))
colnames(ciber_zong)<-ciber_zong[1,]
ciber_zong<-ciber_zong[-1,]
ciber_zong<-as.data.frame(ciber_zong)
ciber_zong[,1:ncol(ciber_zong)]<-apply(ciber_zong[,1:ncol(ciber_zong)],1,as.numeric)
str(ciber_zong)
p_zanshi<-c()

ciber_class<-rep(c(1,2),c(nrow(cluster1_ciber),nrow(cluster2_ciber)))
for (j in 1:nrow(ciber_zong)) {
  T_result<-wilcox.test(as.numeric(ciber_zong[j,])~ciber_class)
  
  p<-T_result$p.value
  
  # t1<-T_result[[1]]
  p_zanshi<-c(p_zanshi,p)
  #t_zanshi<-c(t_zanshi,t1)
}
p_cluster_new<-p_zanshi
p_cluster_new


