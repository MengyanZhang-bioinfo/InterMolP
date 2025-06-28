#BRCA####
 
library(survival)
library(openxlsx)
 
pan_pfs<-read.xlsx("E:\\ PFS.xlsx")
brca_os<-read.table("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\UCSC+BRCA_OS.txt",sep = "\t",header = T)


# i=9
# results = ConsensusClusterPlus(STAD_tumor, maxK = 3, reps = 1000, pItem = 0.8, pFeature = 1, title = "untitled_consensus_cluster", clusterAlg = zuhe[1,i],
#                                distance = zuhe[2,i], seed = 123, tmyPal=NULL, writeTable=FALSE, plot = "pdf")
# clusters<-results[[2]][["consensusClass"]]
clusters<-as.data.frame(clusters)
samples<-substr(gsub('[.]','-',rownames(clusters)),1,15)
clusters<-data.frame(sample=samples,clusters)
brca_os2<-merge(brca_os,clusters,by="sample")
fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
p1
library(survminer)
ccc<-ggsurvplot(fit2,
                pval = TRUE, #conf.int = TRUE, 
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#80D1C8","#FFD4A9","pink"))
ccc

ggsave("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\BRCA_cluster_os_plot.pdf",ccc$plot,width = 10,height = 10)

#pfs
brca_samp<-gsub('[.]','-',substr(rownames(clusters),1,12))
brca_pfs<-pan_pfs[pan_pfs$bcr_patient_barcode%in%brca_samp,]
clusters<-data.frame(bcr_patient_barcode=brca_samp,clusters)
brca_pfs<-merge(brca_pfs,clusters,by="bcr_patient_barcode")
fit2 <- survfit(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)

survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)
p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
p1
ccc<-ggsurvplot(fit2,
                pval = TRUE, #conf.int = TRUE, 
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#80D1C8","#FFD4A9","pink"))
ccc

ggsave("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\BRCA_cluster_pfs_plot.pdf",ccc$plot,width = 10,height = 10)



 
label2<-label[-which(label$label=="0"),]
samples<-substr(gsub('[.]','-',label2$X),1,15)
clusters11<-data.frame(sample=samples,stage=label2$label)
brca_os2<-merge(brca_os,clusters11,by="sample")
fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ stage,data  = brca_os2)
survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ stage,data  = brca_os2)
p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
p1
library(survminer)
ccc<-ggsurvplot(fit2,
                pval = TRUE, #conf.int = TRUE, 
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#D55E00","#56B4E9"))
ccc

ggsave("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\BRCA_stage_os_plot.pdf",ccc$plot,width = 10,height = 10)



clusters22<-data.frame(bcr_patient_barcode=substr(samples,1,12),stage=label2$label)
brca_pfs22<-merge(brca_pfs,clusters22,by="bcr_patient_barcode")
fit2 <- survfit(Surv(as.numeric(PFS.time), PFS) ~ stage,data  = brca_pfs22)

survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ stage,data  = brca_pfs22)
p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
p1
ccc<-ggsurvplot(fit2,
                pval = TRUE, #conf.int = TRUE, 
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#D55E00","#56B4E9"))
ccc
ggsave("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\BRCA_stage_pfs_plot.pdf",ccc$plot,width = 10,height = 10)

#COAD####

 
library(survival)
library(openxlsx)

pan_pfs<-read.xlsx("E:\\ PFS.xlsx")
brca_os<-read.table("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\UCSC+COAD_OS.txt",sep = "\t",header = T)
 

# i=9
# results = ConsensusClusterPlus(STAD_tumor, maxK = 3, reps = 1000, pItem = 0.8, pFeature = 1, title = "untitled_consensus_cluster", clusterAlg = zuhe[1,i],
#                                distance = zuhe[2,i], seed = 123, tmyPal=NULL, writeTable=FALSE, plot = "pdf")
# clusters<-results[[2]][["consensusClass"]]
clusters<-as.data.frame(clusters)
samples<-substr(gsub('[.]','-',rownames(clusters)),1,15)
clusters<-data.frame(sample=samples,clusters)
brca_os2<-merge(brca_os,clusters,by="sample")
fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
p1
library(survminer)
ccc<-ggsurvplot(fit2,
                pval = TRUE, #conf.int = TRUE, 
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#FFD4A9", "#80D1C8","pink"))
ccc

ggsave("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\COAD_cluster_os_plot.pdf",ccc$plot,width = 10,height = 10)


 
label2<-label[-which(label$label=="0"),]
samples<-substr(gsub('[.]','-',label2$X),1,15)
clusters11<-data.frame(sample=samples,stage=label2$label)
brca_os2<-merge(brca_os,clusters11,by="sample")
fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ stage,data  = brca_os2)
survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ stage,data  = brca_os2)
p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
p1
library(survminer)
ccc<-ggsurvplot(fit2,
                pval = TRUE, #conf.int = TRUE, 
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#D55E00","#56B4E9"))
ccc

ggsave("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\COAD_stage_os_plot.pdf",ccc$plot,width = 10,height = 10)

#ESCA####

library(survival)
library(openxlsx)
 
pan_pfs<-read.xlsx("E:\\PFS.xlsx")
brca_os<-read.table("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\UCSC+ESCA_OS.txt",sep = "\t",header = T)
 

# i=9
# results = ConsensusClusterPlus(STAD_tumor, maxK = 3, reps = 1000, pItem = 0.8, pFeature = 1, title = "untitled_consensus_cluster", clusterAlg = zuhe[1,i],
#                                distance = zuhe[2,i], seed = 123, tmyPal=NULL, writeTable=FALSE, plot = "pdf")
# clusters<-results[[2]][["consensusClass"]]
clusters<-as.data.frame(clusters)
samples<-substr(gsub('[.]','-',rownames(clusters)),1,15)
clusters<-data.frame(sample=samples,clusters)
brca_os2<-merge(brca_os,clusters,by="sample")
fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
p1
library(survminer)
ccc<-ggsurvplot(fit2,
                pval = TRUE, #conf.int = TRUE, 
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#FFD4A9", "#80D1C8","pink"))
ccc

ggsave("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\ESCA_cluster_os_plot.pdf",ccc$plot,width = 10,height = 10)



label2<-label[-which(label$label=="0"),]
samples<-substr(gsub('[.]','-',label2$X),1,15)
clusters11<-data.frame(sample=samples,stage=label2$label)
brca_os2<-merge(brca_os,clusters11,by="sample")
fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ stage,data  = brca_os2)
survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ stage,data  = brca_os2)
p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
p1
library(survminer)
ccc<-ggsurvplot(fit2,
                pval = TRUE, #conf.int = TRUE, 
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#D55E00","#56B4E9"))
ccc

ggsave("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\ESCA_stage_os_plot.pdf",ccc$plot,width = 10,height = 10)

#HNSC####

 
library(survival)
library(openxlsx)
 
pan_pfs<-read.xlsx("E:\\ PFS.xlsx")
brca_os<-read.table("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\UCSC+HNSC_OS.txt",sep = "\t",header = T)
 

# i=9
# results = ConsensusClusterPlus(STAD_tumor, maxK = 3, reps = 1000, pItem = 0.8, pFeature = 1, title = "untitled_consensus_cluster", clusterAlg = zuhe[1,i],
#                                distance = zuhe[2,i], seed = 123, tmyPal=NULL, writeTable=FALSE, plot = "pdf")
# clusters<-results[[2]][["consensusClass"]]
clusters<-as.data.frame(clusters)
samples<-substr(gsub('[.]','-',rownames(clusters)),1,15)
clusters<-data.frame(sample=samples,clusters)
brca_os2<-merge(brca_os,clusters,by="sample")
fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
p1
library(survminer)
ccc<-ggsurvplot(fit2,
                pval = TRUE, #conf.int = TRUE, 
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#FFD4A9", "#80D1C8","pink"))
ccc

ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_cluster_os_plot.pdf",ccc$plot,width = 10,height = 10)

#pfs
brca_samp<-gsub('[.]','-',substr(rownames(clusters),1,12))
brca_pfs<-pan_pfs[pan_pfs$bcr_patient_barcode%in%brca_samp,]
clusters<-data.frame(bcr_patient_barcode=brca_samp,clusters)
brca_pfs<-merge(brca_pfs,clusters,by="bcr_patient_barcode")
fit2 <- survfit(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)

survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)
p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
p1
ccc<-ggsurvplot(fit2,
                pval = TRUE, #conf.int = TRUE, 
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#FFD4A9", "#80D1C8","pink"))
ccc

ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_cluster_pfs_plot.pdf",ccc$plot,width = 10,height = 10)



 
label2<-label[-which(label$label=="0"),]
samples<-substr(gsub('[.]','-',label2$X),1,15)
clusters11<-data.frame(sample=samples,stage=label2$label)
brca_os2<-merge(brca_os,clusters11,by="sample")
fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ stage,data  = brca_os2)
survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ stage,data  = brca_os2)
p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
p1
library(survminer)
ccc<-ggsurvplot(fit2,
                pval = TRUE, #conf.int = TRUE, 
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#D55E00","#56B4E9"))
ccc

ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_stage_os_plot.pdf",ccc$plot,width = 10,height = 10)



clusters22<-data.frame(bcr_patient_barcode=substr(samples,1,12),stage=label2$label)
brca_pfs22<-merge(brca_pfs,clusters22,by="bcr_patient_barcode")
fit2 <- survfit(Surv(as.numeric(PFS.time), PFS) ~ stage,data  = brca_pfs22)

survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ stage,data  = brca_pfs22)
p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
p1
ccc<-ggsurvplot(fit2,
                pval = TRUE, #conf.int = TRUE, 
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#D55E00","#56B4E9"))
ccc
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_stage_pfs_plot.pdf",ccc$plot,width = 10,height = 10)

#KIRC####
 
library(survival)
library(openxlsx)
 
pan_pfs<-read.xlsx("E:\\ PFS.xlsx")
brca_os<-read.table("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\UCSC+KIRC_OS.txt",sep = "\t",header = T)
 

# i=9
# results = ConsensusClusterPlus(STAD_tumor, maxK = 3, reps = 1000, pItem = 0.8, pFeature = 1, title = "untitled_consensus_cluster", clusterAlg = zuhe[1,i],
#                                distance = zuhe[2,i], seed = 123, tmyPal=NULL, writeTable=FALSE, plot = "pdf")
# clusters<-results[[2]][["consensusClass"]]
clusters<-as.data.frame(clusters)
samples<-substr(gsub('[.]','-',rownames(clusters)),1,15)
clusters<-data.frame(sample=samples,clusters)
brca_os2<-merge(brca_os,clusters,by="sample")
fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
p1
library(survminer)
ccc<-ggsurvplot(fit2,
                pval = TRUE, #conf.int = TRUE, 
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#FFD4A9", "#80D1C8","pink"))
ccc

ggsave("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\KIRC_cluster_os_plot.pdf",ccc$plot,width = 10,height = 10)

#pfs
brca_samp<-gsub('[.]','-',substr(rownames(clusters),1,12))
brca_pfs<-pan_pfs[pan_pfs$bcr_patient_barcode%in%brca_samp,]
clusters<-data.frame(bcr_patient_barcode=brca_samp,clusters)
brca_pfs<-merge(brca_pfs,clusters,by="bcr_patient_barcode")
fit2 <- survfit(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)

survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)
p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
p1
ccc<-ggsurvplot(fit2,
                pval = TRUE, #conf.int = TRUE, 
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#FFD4A9", "#80D1C8","pink"))
ccc

ggsave("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\KIRC_cluster_pfs_plot.pdf",ccc$plot,width = 10,height = 10)



 
label2<-label[-which(label$label=="0"),]
samples<-substr(gsub('[.]','-',label2$X),1,15)
clusters11<-data.frame(sample=samples,stage=label2$label)
brca_os2<-merge(brca_os,clusters11,by="sample")
fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ stage,data  = brca_os2)
survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ stage,data  = brca_os2)
p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
p1
library(survminer)
ccc<-ggsurvplot(fit2,
                pval = TRUE, #conf.int = TRUE, 
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#D55E00","#56B4E9","#009E73"))
ccc

ggsave("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\KIRC_stage_os_plot.pdf",ccc$plot,width = 10,height = 10)



clusters22<-data.frame(bcr_patient_barcode=substr(samples,1,12),stage=label2$label)
brca_pfs22<-merge(brca_pfs,clusters22,by="bcr_patient_barcode")
fit2 <- survfit(Surv(as.numeric(PFS.time), PFS) ~ stage,data  = brca_pfs22)

survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ stage,data  = brca_pfs22)
p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
p1
ccc<-ggsurvplot(fit2,
                pval = TRUE, #conf.int = TRUE, 
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#D55E00","#56B4E9"))
ccc
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\KIRC_stage_pfs_plot.pdf",ccc$plot,width = 10,height = 10)

#KIRP####

 
library(survival)
library(openxlsx)
 
pan_pfs<-read.xlsx("E:\\ PFS.xlsx")
brca_os<-read.table("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\UCSC+KIRP_OS.txt",sep = "\t",header = T)
 

# i=9
# results = ConsensusClusterPlus(STAD_tumor, maxK = 3, reps = 1000, pItem = 0.8, pFeature = 1, title = "untitled_consensus_cluster", clusterAlg = zuhe[1,i],
#                                distance = zuhe[2,i], seed = 123, tmyPal=NULL, writeTable=FALSE, plot = "pdf")
# clusters<-results[[2]][["consensusClass"]]
clusters<-as.data.frame(clusters)
samples<-substr(gsub('[.]','-',rownames(clusters)),1,15)
clusters<-data.frame(sample=samples,clusters)
brca_os2<-merge(brca_os,clusters,by="sample")
fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
p1
library(survminer)
ccc<-ggsurvplot(fit2,
                pval = TRUE, #conf.int = TRUE, 
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#FFD4A9", "#80D1C8","pink"))
ccc

ggsave("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\KIRP_cluster_os_plot.pdf",ccc$plot,width = 10,height = 10)


 
label2<-label[-which(label$label=="0"),]
samples<-substr(gsub('[.]','-',label2$X),1,15)
clusters11<-data.frame(sample=samples,stage=label2$label)
brca_os2<-merge(brca_os,clusters11,by="sample")
fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ stage,data  = brca_os2)
survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ stage,data  = brca_os2)
p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
p1
library(survminer)
ccc<-ggsurvplot(fit2,
                pval = TRUE, #conf.int = TRUE, 
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#D55E00","#56B4E9"))
ccc

ggsave("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\KIRP_stage_os_plot.pdf",ccc$plot,width = 10,height = 10)

#LIHC####

 
library(survival)
library(openxlsx)
 
pan_pfs<-read.xlsx("E:\\ PFS.xlsx")
brca_os<-read.table("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\UCSC+KIRP_OS.txt",sep = "\t",header = T)
 

# i=9
# results = ConsensusClusterPlus(STAD_tumor, maxK = 3, reps = 1000, pItem = 0.8, pFeature = 1, title = "untitled_consensus_cluster", clusterAlg = zuhe[1,i],
#                                distance = zuhe[2,i], seed = 123, tmyPal=NULL, writeTable=FALSE, plot = "pdf")
# clusters<-results[[2]][["consensusClass"]]
clusters<-as.data.frame(clusters)
samples<-substr(gsub('[.]','-',rownames(clusters)),1,15)
clusters<-data.frame(sample=samples,clusters)
brca_os2<-merge(brca_os,clusters,by="sample")
fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
p1
library(survminer)
ccc<-ggsurvplot(fit2,
                pval = TRUE, #conf.int = TRUE, 
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#FFD4A9", "#80D1C8","pink"))
ccc

ggsave("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\KIRP_cluster_os_plot.pdf",ccc$plot,width = 10,height = 10)


 
label2<-label[-which(label$label=="0"),]
samples<-substr(gsub('[.]','-',label2$X),1,15)
clusters11<-data.frame(sample=samples,stage=label2$label)
brca_os2<-merge(brca_os,clusters11,by="sample")
fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ stage,data  = brca_os2)
survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ stage,data  = brca_os2)
p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
p1
library(survminer)
ccc<-ggsurvplot(fit2,
                pval = TRUE, #conf.int = TRUE, 
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#D55E00","#56B4E9"))
ccc

ggsave("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\KIRP_stage_os_plot.pdf",ccc$plot,width = 10,height = 10)

#LUNG####
 
library(survival)
library(openxlsx)
 
pan_pfs<-read.xlsx("E:\\ PFS.xlsx")
brca_os<-read.table("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\UCSC+NSCLC_OS.txt",sep = "\t",header = T)
 

# i=9
# results = ConsensusClusterPlus(STAD_tumor, maxK = 3, reps = 1000, pItem = 0.8, pFeature = 1, title = "untitled_consensus_cluster", clusterAlg = zuhe[1,i],
#                                distance = zuhe[2,i], seed = 123, tmyPal=NULL, writeTable=FALSE, plot = "pdf")
# clusters<-results[[2]][["consensusClass"]]
clusters<-as.data.frame(clusters)
samples<-substr(gsub('[.]','-',rownames(clusters)),1,15)
clusters<-data.frame(sample=samples,clusters)
brca_os2<-merge(brca_os,clusters,by="sample")
fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
p1
library(survminer)
ccc<-ggsurvplot(fit2,
                pval = TRUE, #conf.int = TRUE, 
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#FFD4A9", "#80D1C8","pink"))
ccc

ggsave("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\NSCLC_cluster_os_plot.pdf",ccc$plot,width = 10,height = 10)

#pfs
brca_samp<-gsub('[.]','-',substr(rownames(clusters),1,12))
brca_pfs<-pan_pfs[pan_pfs$bcr_patient_barcode%in%brca_samp,]
clusters<-data.frame(bcr_patient_barcode=brca_samp,clusters)
brca_pfs<-merge(brca_pfs,clusters,by="bcr_patient_barcode")
fit2 <- survfit(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)

survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)
p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
p1
ccc<-ggsurvplot(fit2,
                pval = TRUE, #conf.int = TRUE, 
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#FFD4A9", "#80D1C8","pink"))
ccc

ggsave("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\NSCLC_cluster_pfs_plot.pdf",ccc$plot,width = 10,height = 10)

#PAAD####
 
library(survival)
library(openxlsx)
 
pan_pfs<-read.xlsx("E:\\ PFS.xlsx")
brca_os<-read.table("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\UCSC+PAAD_OS.txt",sep = "\t",header = T)
 

# i=9
# results = ConsensusClusterPlus(STAD_tumor, maxK = 3, reps = 1000, pItem = 0.8, pFeature = 1, title = "untitled_consensus_cluster", clusterAlg = zuhe[1,i],
#                                distance = zuhe[2,i], seed = 123, tmyPal=NULL, writeTable=FALSE, plot = "pdf")
# clusters<-results[[2]][["consensusClass"]]
clusters<-as.data.frame(clusters)
samples<-substr(gsub('[.]','-',rownames(clusters)),1,15)
clusters<-data.frame(sample=samples,clusters)
brca_os2<-merge(brca_os,clusters,by="sample")
fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
p1
library(survminer)
ccc<-ggsurvplot(fit2,
                pval = TRUE, #conf.int = TRUE, 
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#FFD4A9", "#80D1C8","pink"))
ccc

ggsave("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\PAAD_cluster_os_plot.pdf",ccc$plot,width = 10,height = 10)

#pfs
brca_samp<-gsub('[.]','-',substr(rownames(clusters),1,12))
brca_pfs<-pan_pfs[pan_pfs$bcr_patient_barcode%in%brca_samp,]
clusters<-data.frame(bcr_patient_barcode=brca_samp,clusters)
brca_pfs<-merge(brca_pfs,clusters,by="bcr_patient_barcode")
fit2 <- survfit(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)

survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)
p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
p1
ccc<-ggsurvplot(fit2,
                pval = TRUE, #conf.int = TRUE, 
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#FFD4A9", "#80D1C8","pink"))
ccc

ggsave("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\PAAD_cluster_pfs_plot.pdf",ccc$plot,width = 10,height = 10)



 
label2<-label[-which(label$label=="0"),]
samples<-substr(gsub('[.]','-',label2$X),1,15)
clusters11<-data.frame(sample=samples,stage=label2$label)
brca_os2<-merge(brca_os,clusters11,by="sample")
fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ stage,data  = brca_os2)
survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ stage,data  = brca_os2)
p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
p1
library(survminer)
ccc<-ggsurvplot(fit2,
                pval = TRUE, #conf.int = TRUE, 
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#D55E00","#56B4E9"))
ccc

ggsave("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\PAAD_stage_os_plot.pdf",ccc$plot,width = 10,height = 10)



clusters22<-data.frame(bcr_patient_barcode=substr(samples,1,12),stage=label2$label)
brca_pfs22<-merge(brca_pfs,clusters22,by="bcr_patient_barcode")
fit2 <- survfit(Surv(as.numeric(PFS.time), PFS) ~ stage,data  = brca_pfs22)

survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ stage,data  = brca_pfs22)
p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
p1
ccc<-ggsurvplot(fit2,
                pval = TRUE, #conf.int = TRUE, 
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#D55E00","#56B4E9"))
ccc
ggsave("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\PAAD_stage_pfs_plot.pdf",ccc$plot,width = 10,height = 10)
#THCA####
 
library(survival)
library(openxlsx)
 
pan_pfs<-read.xlsx("E:\\ PFS.xlsx")
brca_os<-read.table("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\UCSC+THCA_OS.txt",sep = "\t",header = T)
 

# i=9
# results = ConsensusClusterPlus(STAD_tumor, maxK = 3, reps = 1000, pItem = 0.8, pFeature = 1, title = "untitled_consensus_cluster", clusterAlg = zuhe[1,i],
#                                distance = zuhe[2,i], seed = 123, tmyPal=NULL, writeTable=FALSE, plot = "pdf")
# clusters<-results[[2]][["consensusClass"]]
clusters<-as.data.frame(clusters)
samples<-substr(gsub('[.]','-',rownames(clusters)),1,15)
clusters<-data.frame(sample=samples,clusters)
brca_os2<-merge(brca_os,clusters,by="sample")
fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
p1
library(survminer)
ccc<-ggsurvplot(fit2,
                pval = TRUE, #conf.int = TRUE, 
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#FFD4A9", "#80D1C8","pink"))
ccc

ggsave("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\THCA_cluster_os_plot.pdf",ccc$plot,width = 10,height = 10)

#pfs
brca_samp<-gsub('[.]','-',substr(rownames(clusters),1,12))
brca_pfs<-pan_pfs[pan_pfs$bcr_patient_barcode%in%brca_samp,]
clusters<-data.frame(bcr_patient_barcode=brca_samp,clusters)
brca_pfs<-merge(brca_pfs,clusters,by="bcr_patient_barcode")
fit2 <- survfit(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)

survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)
p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
p1
ccc<-ggsurvplot(fit2,
                pval = TRUE, #conf.int = TRUE, 
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#FFD4A9", "#80D1C8","pink"))
ccc

ggsave("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\THCA_cluster_pfs_plot.pdf",ccc$plot,width = 10,height = 10)



 
label2<-label[-which(label$label=="0"),]
samples<-substr(gsub('[.]','-',label2$X),1,15)
clusters11<-data.frame(sample=samples,stage=label2$label)
brca_os2<-merge(brca_os,clusters11,by="sample")
fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ stage,data  = brca_os2)
survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ stage,data  = brca_os2)
p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
p1
library(survminer)
ccc<-ggsurvplot(fit2,
                pval = TRUE, #conf.int = TRUE, 
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#D55E00","#56B4E9"))
ccc

ggsave("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\THCA_stage_os_plot.pdf",ccc$plot,width = 10,height = 10)



clusters22<-data.frame(bcr_patient_barcode=substr(samples,1,12),stage=label2$label)
brca_pfs22<-merge(brca_pfs,clusters22,by="bcr_patient_barcode")
fit2 <- survfit(Surv(as.numeric(PFS.time), PFS) ~ stage,data  = brca_pfs22)

survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ stage,data  = brca_pfs22)
p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
p1
ccc<-ggsurvplot(fit2,
                pval = TRUE, #conf.int = TRUE, 
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#D55E00","#56B4E9"))
ccc
ggsave("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\THCA_stage_pfs_plot.pdf",ccc$plot,width = 10,height = 10)

#UCEC####

 
library(survival)
library(openxlsx)
 
pan_pfs<-read.xlsx("E:\\ PFS.xlsx")
brca_os<-read.table("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\UCSC+UCEC_OS.txt",sep = "\t",header = T)
 

# i=9
# results = ConsensusClusterPlus(STAD_tumor, maxK = 3, reps = 1000, pItem = 0.8, pFeature = 1, title = "untitled_consensus_cluster", clusterAlg = zuhe[1,i],
#                                distance = zuhe[2,i], seed = 123, tmyPal=NULL, writeTable=FALSE, plot = "pdf")
# clusters<-results[[2]][["consensusClass"]]
clusters<-as.data.frame(clusters)
samples<-substr(gsub('[.]','-',rownames(clusters)),1,15)
clusters<-data.frame(sample=samples,clusters)
brca_os2<-merge(brca_os,clusters,by="sample")
fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ clusters,data  = brca_os2)
p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
p1
library(survminer)
ccc<-ggsurvplot(fit2,
                pval = TRUE, #conf.int = TRUE, 
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#FFD4A9", "#80D1C8","pink"))
ccc

ggsave("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\UCEC_cluster_os_plot.pdf",ccc$plot,width = 10,height = 10)

#pfs
brca_samp<-gsub('[.]','-',substr(rownames(clusters),1,12))
brca_pfs<-pan_pfs[pan_pfs$bcr_patient_barcode%in%brca_samp,]
clusters<-data.frame(bcr_patient_barcode=brca_samp,clusters)
brca_pfs<-merge(brca_pfs,clusters,by="bcr_patient_barcode")
fit2 <- survfit(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)

survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ clusters,data  = brca_pfs)
p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
p1
ccc<-ggsurvplot(fit2,
                pval = TRUE, #conf.int = TRUE, 
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#FFD4A9", "#80D1C8","pink"))
ccc

ggsave("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\UCEC_cluster_pfs_plot.pdf",ccc$plot,width = 10,height = 10)



 
label2<-label[-which(label$label=="0"),]
samples<-substr(gsub('[.]','-',label2$X),1,15)
clusters11<-data.frame(sample=samples,stage=label2$label)
brca_os2<-merge(brca_os,clusters11,by="sample")
fit2 <- survfit(Surv(as.numeric(OS.time), OS) ~ stage,data  = brca_os2)
survtest <- survdiff(Surv(as.numeric(OS.time), OS) ~ stage,data  = brca_os2)
p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
p1
library(survminer)
ccc<-ggsurvplot(fit2,
                pval = TRUE, #conf.int = TRUE, 
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#D55E00","#56B4E9"))
ccc

ggsave("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\UCEC_stage_os_plot.pdf",ccc$plot,width = 10,height = 10)



clusters22<-data.frame(bcr_patient_barcode=substr(samples,1,12),stage=label2$label)
brca_pfs22<-merge(brca_pfs,clusters22,by="bcr_patient_barcode")
fit2 <- survfit(Surv(as.numeric(PFS.time), PFS) ~ stage,data  = brca_pfs22)

survtest <- survdiff(Surv(as.numeric(PFS.time), PFS) ~ stage,data  = brca_pfs22)
p1 <- 1 - pchisq(survtest$chisq,  length(survtest$n) - 1)
p1
ccc<-ggsurvplot(fit2,
                pval = TRUE, #conf.int = TRUE,
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = "strata", # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                ggtheme = theme_bw(), # Change ggplot2 theme
                palette = c("#D55E00","#56B4E9"))
ccc
ggsave("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\UCEC_stage_pfs_plot.pdf",ccc$plot,width = 10,height = 10)
