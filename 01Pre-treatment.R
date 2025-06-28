#PAAD,ESCA,KIRP,COAD,KIRC,LIHC,UCEC,SKCM,PRAD,THCA,HNSC,BRCA
library(DMwR2)
library(data.table)
library(parallel)
setwd("G:\\dnameth")
n_hebing1<-c()
extracted<-c("PAAD","ESCA","KIRP","COAD","KIRC","LIHC","BLCA","UCEC","SKCM","PRAD","THCA","HNSC","BRCA")
i=11


  LUADcan<-read.table(paste0(extracted[i],"\\",extracted[i],"cancerKNN_exp.txt"),header=TRUE,sep='\t')
  genemap=read.table(paste0(extracted[i],"\\gencode.v22.annotation.gene.probeMap"),header=TRUE,sep='\t')
  genemap<-genemap[,1:2]
  LUADcan$id=rownames(LUADcan)
  LUADcan=merge(LUADcan,genemap,by='id')
  LUADcan=aggregate(LUADcan[,2:(ncol(LUADcan)-1)],by=list(gene=LUADcan[,ncol(LUADcan)]),FUN=mean)

    if( file.exists((paste0(extracted[i],"\\",extracted[i],"normalKNN_exp.txt")))){
      LUADnormal_exp=read.table(paste0(extracted[i],"\\",extracted[i],"normalKNN_exp.txt"),header=TRUE,row.names=1)

  }else{
    LUADnormal_exp=read.table(paste0(extracted[i],"\\",extracted[i],"normalnoKNN_exp.txt"),header=TRUE,row.names=1)

  }
  
 # LUADnormal_exp<-read.table(paste0(extracted[i],"\\",extracted[i],"normalKNN_exp.txt"),header=TRUE,sep='\t')
  LUADnormal_exp$id=rownames(LUADnormal_exp)
  LUADnormal_exp=merge(LUADnormal_exp,genemap,by='id')
  LUADnormal_exp=aggregate(LUADnormal_exp[,2:(ncol(LUADnormal_exp)-1)],by=list(gene=LUADnormal_exp[,ncol(LUADnormal_exp)]),FUN=mean)
  
  n_expcanc_yuan<-ncol(LUADcan)
  n_expnorm_yuan<-ncol(LUADnormal_exp)
  n_expcanc_yuan
  n_expnorm_yuan
 
  rownames(LUADcan)=LUADcan[,1]
  LUADcan=LUADcan[,-1]
  test <- LUADcan
  normalize <- function(x) {
    return((x - min(x)) / (max(x) - min(x)))
  }
  
  test_n <- as.data.frame(lapply(test, normalize))
  head(test_n)
  rownames(test_n)=rownames(LUADcan)
  LUADcan1=test_n 
  
  rownames(LUADnormal_exp)=LUADnormal_exp[,1]
  LUADnormal_exp=LUADnormal_exp[,-1]
  test1 <- LUADnormal_exp
  test_n1 <- as.data.frame(lapply(test1, normalize))
  head(test_n1)
  rownames(test_n1)=rownames(LUADnormal_exp)
  LUADnormal_exp1=test_n1 
  

  LUAD=read.table(paste0(extracted[i],"\\",extracted[i],"cancerKNN.txt"),header=TRUE,row.names=1)
  file.exists((paste0(extracted[i],"\\",extracted[i],"normalKNN.txt")))
  LUADnormal=read.table(paste0(extracted[i],"\\",extracted[i],"normalKNN.txt"),header=TRUE,row.names=1)
 
  n_menorm_yuan<-ncol(LUADnormal)
  n_mecanc_yuan<-ncol(LUAD)
  n_menorm_yuan
  n_mecanc_yuan

  LUADsurvival=fread(paste0(extracted[i],"\\","TCGA-",extracted[i],".GDC_phenotype.tsv.gz"),header=TRUE)
  LUADsurvival<-as.data.frame(LUADsurvival)
  rownames(LUADsurvival)<-LUADsurvival[,1]
  LUADsurvival<-LUADsurvival[,-1]
  rownames(LUADsurvival)=gsub("-","\\.",rownames(LUADsurvival))

  unique(LUADsurvival$clinical_stage)
  unique(LUADsurvival$tumor_stage)
 # unique(LUADsurvival$clinical_stage[matches_stage_iii])
  # matches_stage_i <-grepl("\\bstage i(a|b)?\\b", LUADsurvival$tumor_stage)
  # matches_stage_ii <-grepl("\\bstage ii(a|b|c)?\\b", LUADsurvival$tumor_stage)
  # matches_stage_iii <-grepl("\\bstage iii(a|b|c)?\\b", LUADsurvival$tumor_stage)
  # LUADsurvival$tumor_stage[matches_stage_i]<-"1"
  # LUADsurvival$tumor_stage[matches_stage_ii]<-"2"
  # LUADsurvival$tumor_stage[matches_stage_iii]<-"3"
  # LUADstage1=LUADsurvival[LUADsurvival$tumor_stage==1,]
  # LUADstage2=LUADsurvival[LUADsurvival$tumor_stage==2,]
  # LUADstage3=LUADsurvival[LUADsurvival$tumor_stage==3,]
  matches_stage_i <-grepl("\\bStage I(A\\d*|B\\d*|C\\d*)?\\b", LUADsurvival$clinical_stage)
  matches_stage_ii <-grepl("\\bStage II(A\\d*|B\\d*|C\\d*)?\\b", LUADsurvival$clinical_stage)
  matches_stage_iii <-grepl("\\bStage III(A\\d*|B\\d*|C\\d*)?\\b", LUADsurvival$clinical_stage)
  LUADsurvival$clinical_stage[matches_stage_i]<-"1"
  LUADsurvival$clinical_stage[matches_stage_ii]<-"2"
  LUADsurvival$clinical_stage[matches_stage_iii]<-"3"

  LUADstage1=LUADsurvival[LUADsurvival$clinical_stage==1,]
  LUADstage2=LUADsurvival[LUADsurvival$clinical_stage==2,]
  LUADstage3=LUADsurvival[LUADsurvival$clinical_stage==3,]

  #colnames(LUAD)=gsub("-","\\.",colnames(LUAD))
  LUADstage1_cg=LUAD[,intersect(rownames(LUADstage1),colnames(LUAD))]
  LUADstage2_cg=LUAD[,intersect(rownames(LUADstage2),colnames(LUAD))]
  LUADstage3_cg=LUAD[,intersect(rownames(LUADstage3),colnames(LUAD))]
  LUAD=cbind(LUADstage1_cg,LUADstage2_cg,LUADstage3_cg)
  
  
  patient=intersect(colnames(LUAD),colnames(LUADcan1))
  LUAD=LUAD[,patient]
  LUADcan1=LUADcan1[,patient]
  normalsample=intersect(colnames(LUADnormal),colnames(LUADnormal_exp1))
  LUADnormal=LUADnormal[,normalsample]
  LUADnormal_exp1=LUADnormal_exp1[,normalsample] 
  

  n_norm_pei<-ncol(LUADnormal)
  n_canc_pei<-ncol(LUAD)
  n_norm_pei
  n_canc_pei
  
  cgtest<-read.table("/cgtest.txt",sep = "\t")
  LUADcg=intersect(rownames(LUAD),cgtest[,1])
  rownames(cgtest)=cgtest[,1]
  LUADmethy=cbind(cgtest[LUADcg,],LUAD) 
  LUADmethy=aggregate(LUADmethy[,6:ncol(LUADmethy)],by=list(gene=LUADmethy[,3]),FUN=mean)
  
  
  #LUADnormal=read.table("F:/zmy/F/Pancancer/LUNG/LUADnormalKNN.txt",header=TRUE,sep='\t',row.names=1,check.names=FALSE)
  LUADnormalcg=intersect(rownames(LUADnormal),cgtest[,1])
  #rownames(cgtest)=cgtest[,1]
  LUADnormalmethy=cbind(cgtest[LUADnormalcg,],LUADnormal)
  LUADnormalmethy=aggregate(LUADnormalmethy[,6:ncol(LUADnormalmethy)],by=list(gene=LUADnormalmethy[,3]),FUN=mean)
  
  #write.csv(LUSCnormalmethy,"F:/zmy/F/Pancancer/LUNG/LUSCnormalmethylationgene.csv")
  
  rownames(LUADmethy)<-LUADmethy[,1]
  LUADmethy<-LUADmethy[,-1]
  rownames(LUADnormalmethy)<-LUADnormalmethy[,1]
  LUADnormalmethy<-LUADnormalmethy[,-1]

  methygene=intersect(rownames(LUADmethy),rownames(LUADnormalmethy)) 
  expgene=intersect(rownames(LUADcan1),rownames(LUADnormal_exp1)) 
  gene=intersect(methygene,expgene) 
  

  n_megenes<-length(methygene)
  n_expgenes<-length(expgene)
  n_genes_pei<-length(gene)
  
   
  
  LUADmethy=LUADmethy[gene,]
  LUADcan=LUADcan1[gene,]
  LUADnormalmethy=LUADnormalmethy[gene,]
  LUADnormal_exp=LUADnormal_exp1[gene,]
  
  
  #LUADmethy=read.csv("F:/zmy/F/Pancancer/LUNG/LUADmethylationgene.csv",header=TRUE,sep=',',row.names=1,check.names=FALSE)
  LUADstage1_methy=LUADmethy[,intersect(rownames(LUADstage1),colnames(LUADmethy))]
  LUADstage2_methy=LUADmethy[,intersect(rownames(LUADstage2),colnames(LUADmethy))]
  LUADstage3_methy=LUADmethy[,intersect(rownames(LUADstage3),colnames(LUADmethy))]
  LUADstage1_exp=LUADcan[,colnames(LUADstage1_methy)]
  LUADstage2_exp=LUADcan[,colnames(LUADstage2_methy)]
  LUADstage3_exp=LUADcan[,colnames(LUADstage3_methy)]
  
  me1=LUADstage1_methy
  me2=cbind(LUADstage2_methy,LUADstage3_methy)
  mn=LUADnormalmethy
  exp1=LUADstage1_exp
  exp2=cbind(LUADstage2_exp,LUADstage3_exp)
  en=LUADnormal_exp

  n_stage1<-ncol(LUADstage1_methy)
  n_stage2<-ncol(LUADstage2_methy)
  n_stage3<-ncol(LUADstage3_methy)
  
  n_stage11<-ncol(me1)
  n_stage22<-ncol(me2)
  n_mn<-ncol(mn)
  n_megenes<-0
  n_expgenes<-0
  n_genes_pei<-0
  n_stage1<-0
  n_stage2<-0
  n_stage3<-0
  n_stage11<-0
  n_stage22<-0
  n_mn<-0
  n_hebing<-c(n_expcanc_yuan,n_expnorm_yuan,n_menorm_yuan,n_mecanc_yuan,
              n_norm_pei,n_canc_pei,n_megenes,n_expgenes,n_genes_pei,n_stage1,n_stage2,n_stage3,
              n_stage11,n_stage22,n_mn)
  
  
  n_hebing1<-rbind(n_hebing1,n_hebing)
  colnames(n_hebing1)<-c("n_expcanc_yuan","n_expnorm_yuan","n_menorm_yuan","n_mecanc_yuan",
                         "n_norm_pei","n_canc_pei","n_megenes","n_expgenes","n_genes_pei","n_stage1","n_stage2","n_stage3",
                         "n_stage11","n_stage22","n_mn")
  rownames(n_hebing1)<-c(extracted[1:11])
  write.csv(n_hebing1,"n_hebing1_11.csv")

  sample_sum<-data.table(sample=c(colnames(me1),colnames(me2),colnames(mn)),label=rep(c(1,2,0),c(ncol(me1),ncol(me2),ncol(mn))))
  write.csv(sample_sum,paste0(extracted[i],"\\sample_sum.csv"))
  
  write.csv(me1,paste0(extracted[i],"\\me1.csv"))
  write.csv(me2,paste0(extracted[i],"\\me2.csv"))
  write.csv(mn,paste0(extracted[i],"\\mn.csv"))
  write.csv(exp1,paste0(extracted[i],"\\exp1.csv"))
  write.csv(exp2,paste0(extracted[i],"\\exp2.csv"))
  write.csv(en,paste0(extracted[i],"\\en.csv"))
  
  #me1=cbind(LUADstage1_methy,LUSCstage1_methy)
  #me2=cbind(LUADstage2_methy,LUADstage3_methy,LUSCstage2_methy,LUSCstage3_methy)
  #mn=cbind(LUADnormalmethy,LUSCnormalmethy)
  

  #####cor
  
  expstage=cbind(exp1,exp2)
  mestage=cbind(me1,me2)
  ncol(expstage)
  ncol(mestage)
  nrow(expstage)
  nrow(mestage)
  data=data.frame()
  for (j in 1:nrow(mestage)){
    
    cor=cor(as.numeric(expstage[j,]),as.numeric(mestage[j,]))
    data=rbind(data,cor)
    
  }
  rownames(data)=rownames(expstage)
  
  write.csv(data,paste0(extracted[i],"\\genecor.csv"))
  print(i)
}


#####LUAD###
cgtest=read.table("F:/zmy/F/Pancancer/cgtest.txt",sep='\t',header=FALSE)
GPL3=read.table("F:/zmy/F/breast_stage/data/GPL3.txt",sep='\t',header=TRUE)
cg=intersect(cgtest[,1],GPL3[,1])
LUAD=read.table("F:\\zmy\\F\\Pancancer\\LUAD\\TCGA-LUAD.methylation450.tsv",header=TRUE,sep='\t')
rownames(LUAD)=LUAD[,1]
cancer=LUAD[,grep("\\.01",colnames(LUAD))]
cancer1=LUAD[,grep("\\.06",colnames(LUAD))]
cancer=cbind(cancer,cancer1)
normal=LUAD[,grep("\\.11",colnames(LUAD))]
cancerknn=cancer[cg,]
normalknn=normal[cg,]

cancerknn=cancerknn[which(rowMeans(!is.na(cancerknn)) > 0.5), ]
normalknn=normalknn[which(rowMeans(!is.na(cancerknn)) > 0.5), ]

##knn
install.packages("devtools")
library(devtools)  # You need to install this package!
install_github("ltorgo/DMwR2",ref="master")
library(DMwR2)##
cancerKNN=knnImputation(cancerknn, k = 10, scale = F, meth = "weighAvg", distData = NULL)
write.table(cancerKNN,"F:/Pancancer/LUNG/cancerKNN.txt",sep='\t')
LUADnormal=knnImputation(normalknn, k = 10, scale = F, meth = "weighAvg", distData = NULL)
write.table(normal,"F:/Pancancer/LUNG/normalKNN.txt",sep='\t')

LUSC=read.table("F:/Pancancer/LUSC/TCGA-LUSC.methylation450.tsv",header=TRUE,sep='\t')
rownames(LUSC)=LUSC[,1]
LUSCcancer=LUSC[,grep("\\.01",colnames(LUSC))]
LUSCcancer1=LUSC[,grep("\\.06",colnames(LUSC))]
LUSCcancer=cbind(LUSCcancer,LUSCcancer1)
LUSCnormal=LUSC[,grep("\\.11",colnames(LUSC))]
LUSCcancerknn=LUSCcancer[cg,]
LUSCnormalknn=LUSCnormal[cg,]
write.table(LUSCcancerknn,"F:/Pancancer/LUNG/LUSCcancerpreknn.txt",sep='\t')
write.table(LUSCnormalknn,"F:/Pancancer/LUNG/LUSCnormalpreknn.txt",sep='\t')
LUSCcancerknn=LUSCcancerknn[which(rowMeans(!is.na(LUSCcancerknn)) > 0.5), ]
LUSCnormalknn=LUSCnormalknn[which(rowMeans(!is.na(LUSCnormalknn)) > 0.5), ]

##knn
install.packages("devtools")
library(devtools)  # You need to install this package!
install_github("ltorgo/DMwR2",ref="master")
library(DMwR2)##
LUSCcancerKNN=knnImputation(LUSCcancerknn, k = 10, scale = F, meth = "weighAvg", distData = NULL)
write.table(LUSCcancerKNN,"F:/Pancancer/LUNG/LUSCcancerKNN.txt",sep='\t')
LUSCnormal=knnImputation(LUSCnormalknn, k = 10, scale = F, meth = "weighAvg", distData = NULL)
write.table(LUSCnormal,"F:/Pancancer/LUNG/LUSCnormalKNN.txt",sep='\t')

####exp KNN#############
LUADexp=read.table("F:\\zmy\\F\\Pancancer\\LUAD\\TCGA-LUAD.htseq_fpkm.tsv",header=TRUE,sep='\t',row.names=1)
LUADnormal_exp=LUADexp[,grep("\\.11",colnames(LUADexp))]
LUADexp=LUADexp[,grep("\\.01",colnames(LUADexp))]
LUADexp1=LUADexp[,grep("\\.06",colnames(LUADexp))]
LUADexp=cbind(LUADexp,LUADexp1)

n0 <- apply(LUADexp == 0, 1, sum)
i0 <- which(n0 > 367)
LUADexp=LUADexp[-i0, ]
n0 <- apply(LUADnormal_exp == 0, 1, sum)
i0 <- which(n0 > 41)
LUADnormal_exp=LUADnormal_exp[-i0, ]
gene=intersect(rownames(LUADnormal_exp),rownames(LUADexp))
LUADnormal_exp=LUADnormal_exp[gene,]
LUADexp=LUADexp[gene,]

LUADexp=knnImputation(LUADexp, k = 10, scale = F, meth = "weighAvg", distData = NULL) 
LUADnormal_exp=knnImputation(LUADnormal_exp, k = 10, scale = F, meth = "weighAvg", distData = NULL)
genemap=read.table("F:\\zmy\\F\\Pancancer\\LUAD\\gencode.v22.annotation.gene.probeMap",header=TRUE,sep='\t')
LUADexp$id=rownames(LUADexp)
LUADexp=merge(LUADexp,genemap,by='id')
LUADexp=aggregate(LUADexp[,2:525],by=list(gene=LUADexp[,526]),FUN=mean)
LUADnormal_exp$id=rownames(LUADnormal_exp)
LUADnormal_exp=merge(LUADnormal_exp,genemap,by='id')
LUADnormal_exp=aggregate(LUADnormal_exp[,2:60],by=list(gene=LUADnormal_exp[,61]),FUN=mean)

####0-1 Standardization
rownames(LUADexp)=LUADexp[,1]
LUADexp=LUADexp[,-1]
test <- LUADexp
normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

test_n <- as.data.frame(lapply(test, normalize))
head(test_n)
rownames(test_n)=rownames(LUADexp)
LUADexp1=test_n 
rownames(LUADnormal_exp)=LUADnormal_exp[,1]
LUADnormal_exp=LUADnormal_exp[,-1]
test1 <- LUADnormal_exp
test_n1 <- as.data.frame(lapply(test1, normalize))
head(test_n1)
rownames(test_n1)=rownames(LUADnormal_exp)
LUADnormal_exp1=test_n1 

#####stage
LUAD=read.table("F:/zmy/F/Pancancer/LUNG/LUADcancerKNN.txt",header=TRUE,sep='\t',row.names=1,check.names=FALSE)
LUADsurvival=read.table("F:/zmy/F/Pancancer/LUAD/TCGA-LUAD.GDC_phenotype.tsv",header=TRUE,sep='\t',row.names=1)
rownames(LUADsurvival)=gsub("-","\\.",rownames(LUADsurvival))
LUADstage1=LUADsurvival[LUADsurvival$tumor_stage==1,]
LUADstage2=LUADsurvival[LUADsurvival$tumor_stage==2,]
LUADstage3=LUADsurvival[LUADsurvival$tumor_stage==3,]
LUADstage1_cg=LUAD[,intersect(rownames(LUADstage1),colnames(LUAD))]
LUADstage2_cg=LUAD[,intersect(rownames(LUADstage2),colnames(LUAD))]
LUADstage3_cg=LUAD[,intersect(rownames(LUADstage3),colnames(LUAD))]
LUAD=cbind(LUADstage1_cg,LUADstage2_cg,LUADstage3_cg)
patient=intersect(colnames(LUAD),colnames(LUADexp1))
LUAD=LUAD[,patient]
LUADexp1=LUADexp1[,patient]
normalsample=intersect(colnames(LUADnormal),colnames(LUADnormal_exp1))
LUADnormal=LUADnormal[,normalsample]
LUADnormal_exp1=LUADnormal_exp1[,normalsample]

LUSC=read.table("F:/zmy/F/Pancancer/LUNG/LUSCcancerKNN.txt",header=TRUE,sep='\t',row.names=1,check.names=FALSE)
LUSCsurvival=read.table("F:/zmy/F/Pancancer/LUSC/TCGA-LUSC.GDC_phenotype.tsv",header=TRUE,sep='\t',row.names=1)
rownames(LUSCsurvival)=gsub("-","\\.",rownames(LUSCsurvival))
LUSCstage1=LUSCsurvival[LUSCsurvival$tumor_stage.diagnoses==1,]
LUSCstage2=LUSCsurvival[LUSCsurvival$tumor_stage.diagnoses==2,]
LUSCstage3=LUSCsurvival[LUSCsurvival$tumor_stage.diagnoses==3,]
LUSCstage1_cg=LUSC[,intersect(rownames(LUSCstage1),colnames(LUSC))]
LUSCstage2_cg=LUSC[,intersect(rownames(LUSCstage2),colnames(LUSC))]
LUSCstage3_cg=LUSC[,intersect(rownames(LUSCstage3),colnames(LUSC))]
LUSC=cbind(LUSCstage1_cg,LUSCstage2_cg,LUSCstage3_cg)
patient=intersect(colnames(LUSC),colnames(LUSCexp))
LUSC=LUSC[,patient]
LUSCexp=LUSCexp[,patient]
normalsample1=intersect(colnames(LUSCnormal),colnames(LUSCnormal_exp1))
LUSCnormal=LUSCnormal[,normalsample1]
LUSCnormal_exp1=LUSCnormal_exp1[,normalsample1] ###


LUADcg=intersect(rownames(LUAD),cgtest[,1])
rownames(cgtest)=cgtest[,1]
LUADmethy=cbind(cgtest[LUADcg,],LUAD)
LUADmethy=aggregate(LUADmethy[,6:443],by=list(gene=LUADmethy[,3]),FUN=mean)

LUSCcg=intersect(rownames(LUSC),cgtest[,1])
LUSCmethy=cbind(cgtest[LUSCcg,],LUSC)
LUSCmethy=aggregate(LUSCmethy[,6:368],by=list(gene=LUSCmethy[,3]),FUN=mean)

#LUADnormal=read.table("F:/zmy/F/Pancancer/LUNG/LUADnormalKNN.txt",header=TRUE,sep='\t',row.names=1,check.names=FALSE)
LUADnormalcg=intersect(rownames(LUADnormal),cgtest[,1])
rownames(cgtest)=cgtest[,1]
LUADnormalmethy=cbind(cgtest[LUADnormalcg,],LUADnormal)
LUADnormalmethy=aggregate(LUADnormalmethy[,6:26],by=list(gene=LUADnormalmethy[,3]),FUN=mean)

LUSCnormal=read.table("F:/zmy/F/Pancancer/LUNG/LUSCnormalKNN.txt",header=TRUE,sep='\t',row.names=1,check.names=FALSE)
LUSCnormalcg=intersect(rownames(LUSCnormal),cgtest[,1])
rownames(cgtest)=cgtest[,1]
LUSCnormalmethy=cbind(cgtest[LUSCnormalcg,],LUSCnormal)
LUSCnormalmethy=aggregate(LUSCnormalmethy[,6:47],by=list(gene=LUSCnormalmethy[,3]),FUN=mean)
#write.csv(LUSCnormalmethy,"F:/zmy/F/Pancancer/LUNG/LUSCnormalmethylationgene.csv")



methygene=intersect(rownames(LUADmethy),rownames(LUADnormalmethy)) 
expgene=intersect(rownames(LUADexp1),rownames(LUADnormal_exp1)) 
gene=intersect(methygene,expgene) 

LUADmethy=LUADmethy[gene,]
LUADexp=LUADexp1[gene,]
LUADnormalmethy=LUADnormalmethy[gene,]
LUADnormal_exp=LUADnormal_exp1[gene,]


#LUADmethy=read.csv("F:/zmy/F/Pancancer/LUNG/LUADmethylationgene.csv",header=TRUE,sep=',',row.names=1,check.names=FALSE)
LUADstage1_methy=LUADmethy[,intersect(rownames(LUADstage1),colnames(LUADmethy))]
LUADstage2_methy=LUADmethy[,intersect(rownames(LUADstage2),colnames(LUADmethy))]
LUADstage3_methy=LUADmethy[,intersect(rownames(LUADstage3),colnames(LUADmethy))]
LUADstage1_exp=LUADexp[,colnames(LUADstage1_methy)]
LUADstage2_exp=LUADexp[,colnames(LUADstage2_methy)]
LUADstage3_exp=LUADexp[,colnames(LUADstage3_methy)]

LUSCmethy=read.csv("F:/zmy/F/Pancancer/LUNG/LUSCmethylationgene.csv",header=TRUE,sep=',',row.names=1,check.names=FALSE)
LUSCstage1_methy=LUSCmethy[,intersect(rownames(LUSCstage1),colnames(LUSCmethy))]
LUSCstage2_methy=LUSCmethy[,intersect(rownames(LUSCstage2),colnames(LUSCmethy))]
LUSCstage3_methy=LUSCmethy[,intersect(rownames(LUSCstage3),colnames(LUSCmethy))]

me1=LUADstage1_methy
me2=cbind(LUADstage2_methy,LUADstage3_methy)
mn=LUADnormalmethy
exp1=LUADstage1_exp
exp2=cbind(LUADstage2_exp,LUADstage3_exp)
en=LUADnormalmethy
write.csv(me1,"E:\\Pancancer\\LUAD\\me1.csv")
write.csv(me2,"E:\\Pancancer\\LUAD\\me2.csv")
write.csv(mn,"E:\\Pancancer\\LUAD\\mn.csv")
write.csv(exp1,"E:\\Pancancer\\LUAD\\exp1.csv")
write.csv(exp2,"E:\\Pancancer\\LUAD\\exp2.csv")
write.csv(en,"E:\\Pancancer\\LUAD\\en.csv")

#me1=cbind(LUADstage1_methy,LUSCstage1_methy)
#me2=cbind(LUADstage2_methy,LUADstage3_methy,LUSCstage2_methy,LUSCstage3_methy)
#mn=cbind(LUADnormalmethy,LUSCnormalmethy)

exp1=cbind(LUADstage1_exp,LUSCstage1_exp)
exp2=cbind(LUADstage2_exp,LUADstage3_exp,LUSCstage2_exp,LUSCstage3_exp)
mn=cbind(LUADnormalmethy,LUSCnormalmethy)



write.csv(expstage1,"./model/expstage11.csv")
write.csv(expstage2,"./model/expstage22.csv")
write.csv(mestage1,"./model/mestage11.csv")
write.csv(mestage2,"./model/mestage22.csv")
write.csv(normalexp,"./model/normalexp.csv")
write.csv(normalme,"./model/normalme.csv")


expstage=cbind(expstage1,expstage2)
mestage=cbind(mestage1,mestage2)
data=data.frame()
for (i in 1:14355){
  
  cor=cor(as.numeric(expstage[i,]),as.numeric(mestage[i,]))
  data=rbind(data,cor)
  
}
rownames(data)=rownames(expstage)

write.csv(data,"./model/genecor.csv")




