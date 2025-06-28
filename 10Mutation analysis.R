#BRCA####
#1.12####
# R
# setwd("G:\\pan_result\\maftools")
# library(TCGAbiolinks)
# 
# query <- GDCquery(
#   project = "TCGA-BRCA",
#   data.category = "Simple Nucleotide Variation",
#   data.type = "Masked Somatic Mutation",
#   access = "open"
# )
# 
# GDCdownload(query)
# 
# GDCprepare(query, save = T,save.filename = "TCGA-BRCA_SNP.Rdata") # Rdatamaftools
# setwd("G:\\pan_result\\maftools")
# 
# cancer_name<-toupper(c("brca","BRCA","BRCA","BRCA","kirc","kirp","BRCA","lung","BRCA","BRCA","ucec","LUAD","LUSC","BRCA","BRCA","ucec"))
# # cancer_name<-toupper(c("LUAD","LUSC","BRCA","BRCA","ucec"))
# 
# # name<-cancer_name[1]
# for (name in cancer_name) {
#   query <- GDCquery(
#     project = paste0("TCGA-",name),
#     data.category = "Simple Nucleotide Variation",
#     data.type = "Masked Somatic Mutation",
#     access = "open"
#   )
#   
#   GDCdownload(query)
#   
#  # GDCprepare(query, save = T,save.filename =paste0("TCGA-",name,"_SNP.Rdata")) # Rdatamaftools
#   
# }
# 
# sessionInfo()

setwd("G:\\pan_result\\maftools\\GDCdata\\TCGA-BRCA\\Simple_Nucleotide_Variation")
library(tidyverse)
# devtools::install_github("PoisonAlien/maftools")
library(maftools)
#dir"maf_data/"
maffilepath=dir(path = "Masked_Somatic_Mutation\\", pattern="masked.maf.gz$", full.names = T, recursive=T)#dir(path = "maf_data/") ,"masked.maf.gz$"masked.maf.gz
head(maffilepath)
#mafdata  MafFrame  MAF  NA
mafdata <- lapply(maffilepath, function(x) {  
  tryCatch({    read.maf(x, isTCGA = TRUE)  
  }, error = function(e) {    
    message(paste("Error reading file:", x, "\n", e$message))  
    return(NA)  })
})
#MAF
mafdata_clean1 <- mafdata[!sapply(mafdata, is.na)]
#maf
merger_maf=merge_mafs(mafdata_clean1)
save(merger_maf,file=c("TCGA-BRCA_SNP.Rdata"))#rda

class(merger_maf)
slotNames(merger_maf)


maf_cluster<-data.frame(sample=substr(gsub('[.]','-',rownames(clusters)),1,12),cluster=clusters$clusters)



# MAF
plotmafSummary(maf = merger_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
merger_maf_1 <- subsetMaf(maf = merger_maf, tsb = maf_cluster$sample[which(maf_cluster$cluster == "1")])
merger_maf_2 <- subsetMaf(maf = merger_maf, tsb = maf_cluster$sample[which(maf_cluster$cluster == "2")])
pdf("G:\\pan_result\\maftools\\BRCA_cluster1.pdf")
plotmafSummary(maf = merger_maf_1, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
dev.off()
pdf("G:\\pan_result\\maftools\\BRCA_cluster2.pdf")
plotmafSummary(maf = merger_maf_2, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
dev.off()
class(merger_maf)
showMethods(classes = class(merger_maf))
merger_maf@variants.per.sample

#samplevariants
variants_per_sample <- merger_maf@variants.per.sample
variants_per_sample1<- merger_maf_1@variants.per.sample
variants_per_sample2<- merger_maf_2@variants.per.sample
variants_per_sample<-as.data.frame(variants_per_sample)
t.test(variants_per_sample1$Variants,variants_per_sample2$Variants)

#tmb
BRCA.tmb1 <- tmb(merger_maf_1, captureSize = 38, logScale = T)
BRCA.tmb2 <- tmb(merger_maf_2, captureSize = 38, logScale = T)
wilcox.test(BRCA.tmb1$total_perMB,BRCA.tmb2$total_perMB)
BRCA.tmb1$cluster<-"1"
BRCA.tmb2$cluster<-"2"
BRCA.tmb<-rbind(BRCA.tmb1,BRCA.tmb2)
colnames(BRCA.tmb)[3]<-"TMB"
str(BRCA.tmb) # 
BRCA.tmb <- BRCA.tmb %>%
  filter(BRCA.tmb$TMB < quantile(BRCA.tmb$TMB, 0.95)) # 95%



p <- ggplot(BRCA.tmb, aes(x = cluster, y = TMB, color = cluster)) +
  geom_boxplot(width=0.5,size=2, alpha = 0.2, notch = TRUE, notchwidth = 0.5) + # 
  scale_color_manual(values = c("#80D1C8","#FFD4A9",  "pink")) + # 
  theme_bw() + # 
  theme(
    panel.grid = element_blank(), # 
    panel.border = element_blank(), # 
    axis.line = element_line(colour = "black") # 
  ) +
  ylab("TMB score") + # y
  xlab(NULL) + # x
  guides(scale = "none") + # 
  geom_signif(comparisons = list(c("1", "2")),
              map_signif_level = TRUE,
              textsize = 8,color="black") + # 
  
  geom_jitter(
    width = 0.1, # 
    size = 1.5, # 
    aes(color = cluster), #  fill
    shape = 21, # 
    alpha = 0.6 # 
  )

print(p)
ggsave("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\BRCA_TMB_plot.pdf",p,width = 10,height = 10)

#cluster1 cluster2####
length(which(DEG_all$FC>=1))
length(which(DEG_all$FC<1))

#COAD####

#1.12####
# R
# library(TCGAbiolinks)
# 
# query <- GDCquery(
#   project = "TCGA-COAD", 
#   data.category = "Simple Nucleotide Variation",
#   data.type = "Masked Somatic Mutation",
#   access = "open"
# )
# 
# GDCdownload(query)
# 
# GDCprepare(query, save = T,save.filename = "TCGA-COAD_SNP.Rdata") # Rdatamaftools
# setwd("G:\\pan_result\\maftools")
# 
# cancer_name<-toupper(c("brca","coad","esca","hnsc","kirc","kirp","lihc","lung","paad","thca","ucec","LUAD","LUSC","paad","thca","ucec"))
# # cancer_name<-toupper(c("LUAD","LUSC","paad","thca","ucec"))
# 
# # name<-cancer_name[1]
# for (name in cancer_name) {
#   query <- GDCquery(
#     project = paste0("TCGA-",name),
#     data.category = "Simple Nucleotide Variation",
#     data.type = "Masked Somatic Mutation",
#     access = "open"
#   )
#   
#   GDCdownload(query)
#   
#  # GDCprepare(query, save = T,save.filename =paste0("TCGA-",name,"_SNP.Rdata")) # Rdatamaftools
#   
# }
# 
# sessionInfo()

setwd("G:\\pan_result\\maftools\\GDCdata\\TCGA-COAD\\Simple_Nucleotide_Variation")
library(tidyverse)
# devtools::install_github("PoisonAlien/maftools")
library(maftools)
#dir"maf_data/"
maffilepath=dir(path = "Masked_Somatic_Mutation\\", pattern="masked.maf.gz$", full.names = T, recursive=T)#dir(path = "maf_data/") ,"masked.maf.gz$"masked.maf.gz
head(maffilepath)
#mafdata  MafFrame  MAF  NA
mafdata <- lapply(maffilepath, function(x) {  
  tryCatch({    read.maf(x, isTCGA = TRUE)  
  }, error = function(e) {    
    message(paste("Error reading file:", x, "\n", e$message))  
    return(NA)  })
})
#MAF
mafdata_clean1 <- mafdata[!sapply(mafdata, is.na)]
#maf
merger_maf=merge_mafs(mafdata_clean1)
save(merger_maf,file=c("TCGA-COAD_SNP.Rdata"))#rda

class(merger_maf)
slotNames(merger_maf)
maf_cluster<-data.frame(sample=substr(gsub('[.]','-',rownames(clusters)),1,12),cluster=clusters$clusters)



# MAF
plotmafSummary(maf = merger_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
merger_maf_1 <- subsetMaf(maf = merger_maf, tsb = maf_cluster$sample[which(maf_cluster$cluster == "1")])
merger_maf_2 <- subsetMaf(maf = merger_maf, tsb = maf_cluster$sample[which(maf_cluster$cluster == "2")])
pdf("G:\\pan_result\\maftools\\coad_cluster1.pdf")
plotmafSummary(maf = merger_maf_1, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
dev.off()
pdf("G:\\pan_result\\maftools\\coad_cluster2.pdf")
plotmafSummary(maf = merger_maf_2, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
dev.off()
class(merger_maf)
showMethods(classes = class(merger_maf))
merger_maf@variants.per.sample

#samplevariants
variants_per_sample <- merger_maf@variants.per.sample
variants_per_sample1<- merger_maf_1@variants.per.sample
variants_per_sample2<- merger_maf_2@variants.per.sample
variants_per_sample<-as.data.frame()
t.test(variants_per_sample1$Variants,variants_per_sample2$Variants)

#tmb
COAD.tmb1 <- tmb(merger_maf_1, captureSize = 38, logScale = T)
COAD.tmb2 <- tmb(merger_maf_2, captureSize = 38, logScale = T)
wilcox.test(COAD.tmb1$total_perMB,COAD.tmb2$total_perMB)
COAD.tmb1$cluster<-"1"
COAD.tmb2$cluster<-"2"
COAD.tmb<-rbind(COAD.tmb1,COAD.tmb2)
colnames(COAD.tmb)[3]<-"TMB"
COAD.tmb <- COAD.tmb %>%
  filter(COAD.tmb$TMB < quantile(COAD.tmb$TMB, 0.95)) # 95%



p <- ggplot(COAD.tmb, aes(x = cluster, y = TMB, color = cluster)) +
  geom_boxplot(width=0.5,size=2, alpha = 0.2, notch = TRUE, notchwidth = 0.5) + # 
  scale_color_manual(values = c("#FFD4A9", "#80D1C8", "pink")) + # 
  theme_bw() + # 
  theme(
    panel.grid = element_blank(), # 
    panel.border = element_blank(), # 
    axis.line = element_line(colour = "black") # 
  ) +
  ylab("HRD score") + # y
  xlab(NULL) + # x
  guides(scale = "none") + # 
  geom_signif(comparisons = list(c("1", "2")),
              map_signif_level = TRUE,
              textsize = 8,color="black") + # 
  
  geom_jitter(
    width = 0.1, # 
    size = 1.5, # 
    aes(fill = cluster), # 
    shape = 21, # 
    alpha = 0.6 # 
  )

print(p)
ggsave("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\COAD_TMB_plot.pdf",p,width = 10,height = 10)
#ESCA####
#1.12####
# R
# library(TCGAbiolinks)
# 
# query <- GDCquery(
#   project = "TCGA-ESCA", 
#   data.category = "Simple Nucleotide Variation",
#   data.type = "Masked Somatic Mutation",
#   access = "open"
# )
# 
# GDCdownload(query)
# 
# GDCprepare(query, save = T,save.filename = "TCGA-ESCA_SNP.Rdata") # Rdatamaftools
# setwd("G:\\pan_result\\maftools")
# 
# cancer_name<-toupper(c("brca","ESCA","esca","hnsc","kirc","kirp","lihc","lung","paad","thca","ucec","LUAD","LUSC","paad","thca","ucec"))
# # cancer_name<-toupper(c("LUAD","LUSC","paad","thca","ucec"))
# 
# # name<-cancer_name[1]
# for (name in cancer_name) {
#   query <- GDCquery(
#     project = paste0("TCGA-",name),
#     data.category = "Simple Nucleotide Variation",
#     data.type = "Masked Somatic Mutation",
#     access = "open"
#   )
#   
#   GDCdownload(query)
#   
#  # GDCprepare(query, save = T,save.filename =paste0("TCGA-",name,"_SNP.Rdata")) # Rdatamaftools
#   
# }
# 
# sessionInfo()

setwd("G:\\pan_result\\maftools\\GDCdata\\TCGA-ESCA\\Simple_Nucleotide_Variation")
library(tidyverse)
# devtools::install_github("PoisonAlien/maftools")
library(maftools)
#dir"maf_data/"
maffilepath=dir(path = "Masked_Somatic_Mutation\\", pattern="masked.maf.gz$", full.names = T, recursive=T)#dir(path = "maf_data/") ,"masked.maf.gz$"masked.maf.gz
head(maffilepath)
#mafdata  MafFrame  MAF  NA
mafdata <- lapply(maffilepath, function(x) {  
  tryCatch({    read.maf(x, isTCGA = TRUE)  
  }, error = function(e) {    
    message(paste("Error reading file:", x, "\n", e$message))  
    return(NA)  })
})
#MAF
mafdata_clean1 <- mafdata[!sapply(mafdata, is.na)]
#maf
merger_maf=merge_mafs(mafdata_clean1)
save(merger_maf,file=c("TCGA-ESCA_SNP.Rdata"))#rda

class(merger_maf)
slotNames(merger_maf)
maf_cluster<-data.frame(sample=substr(gsub('[.]','-',rownames(clusters)),1,12),cluster=clusters$clusters)



# MAF
plotmafSummary(maf = merger_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
merger_maf_1 <- subsetMaf(maf = merger_maf, tsb = maf_cluster$sample[which(maf_cluster$cluster == "1")])
merger_maf_2 <- subsetMaf(maf = merger_maf, tsb = maf_cluster$sample[which(maf_cluster$cluster == "2")])
pdf("G:\\pan_result\\maftools\\ESCA_cluster1.pdf")
plotmafSummary(maf = merger_maf_1, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
dev.off()
pdf("G:\\pan_result\\maftools\\ESCA_cluster2.pdf")
plotmafSummary(maf = merger_maf_2, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
dev.off()
class(merger_maf)
showMethods(classes = class(merger_maf))
merger_maf@variants.per.sample

#samplevariants
variants_per_sample <- merger_maf@variants.per.sample
variants_per_sample1<- merger_maf_1@variants.per.sample
variants_per_sample2<- merger_maf_2@variants.per.sample
variants_per_sample<-as.data.frame()
t.test(variants_per_sample1$Variants,variants_per_sample2$Variants)

#tmb
ESCA.tmb1 <- tmb(merger_maf_1, captureSize = 38, logScale = T)
ESCA.tmb2 <- tmb(merger_maf_2, captureSize = 38, logScale = T)
wilcox.test(ESCA.tmb1$total_perMB,ESCA.tmb2$total_perMB)
ESCA.tmb1$cluster<-"1"
ESCA.tmb2$cluster<-"2"
ESCA.tmb<-rbind(ESCA.tmb1,ESCA.tmb2)
colnames(ESCA.tmb)[3]<-"TMB"

str(ESCA.tmb) # 
ESCA.tmb <- ESCA.tmb %>%
  filter(ESCA.tmb$TMB < quantile(ESCA.tmb$TMB, 0.95)) # 95%


p <- ggplot(ESCA.tmb, aes(x = cluster, y = TMB, color = cluster)) +
  geom_boxplot(width=0.5,size=2, alpha = 0.2, notch = TRUE, notchwidth = 0.5) + # 
  scale_color_manual(values = c("#FFD4A9", "#80D1C8", "pink")) + # 
  theme_bw() + # 
  theme(
    panel.grid = element_blank(), # 
    panel.border = element_blank(), # 
    axis.line = element_line(colour = "black") # 
  ) +
  ylab("HRD score") + # y
  xlab(NULL) + # x
  guides(scale = "none") + # 
  geom_signif(comparisons = list(c("1", "2")),
              map_signif_level = TRUE,
              textsize = 8,color="black") + # 
  
  geom_jitter(
    width = 0.1, # 
    size = 1.5, # 
    aes(fill = cluster), # 
    shape = 21, # 
    alpha = 0.6 # 
  )

print(p)
ggsave("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\ESCA_TMB_plot.pdf",p,width = 10,height = 10)

#HNSC####
#1.12####
# R
# setwd("G:\\pan_result\\maftools")
# library(TCGAbiolinks)
# 
# query <- GDCquery(
#   project = "TCGA-HNSC",
#   data.category = "Simple Nucleotide Variation",
#   data.type = "Masked Somatic Mutation",
#   access = "open"
# )
# 
# GDCdownload(query)
# 
# GDCprepare(query, save = T,save.filename = "TCGA-HNSC_SNP.Rdata") # Rdatamaftools
# setwd("G:\\pan_result\\maftools")
# 
# cancer_name<-toupper(c("brca","HNSC","HNSC","hnsc","kirc","kirp","HNSC","lung","HNSC","HNSC","ucec","LUAD","LUSC","HNSC","HNSC","ucec"))
# # cancer_name<-toupper(c("LUAD","LUSC","HNSC","HNSC","ucec"))
# 
# # name<-cancer_name[1]
# for (name in cancer_name) {
#   query <- GDCquery(
#     project = paste0("TCGA-",name),
#     data.category = "Simple Nucleotide Variation",
#     data.type = "Masked Somatic Mutation",
#     access = "open"
#   )
#   
#   GDCdownload(query)
#   
#  # GDCprepare(query, save = T,save.filename =paste0("TCGA-",name,"_SNP.Rdata")) # Rdatamaftools
#   
# }
# 
# sessionInfo()

setwd("G:\\pan_result\\maftools\\GDCdata\\TCGA-HNSC\\Simple_Nucleotide_Variation")
library(tidyverse)
# devtools::install_github("PoisonAlien/maftools")
library(maftools)
#dir"maf_data/"
maffilepath=dir(path = "Masked_Somatic_Mutation\\", pattern="masked.maf.gz$", full.names = T, recursive=T)#dir(path = "maf_data/") ,"masked.maf.gz$"masked.maf.gz
head(maffilepath)
#mafdata  MafFrame  MAF  NA
mafdata <- lapply(maffilepath, function(x) {  
  tryCatch({    read.maf(x, isTCGA = TRUE)  
  }, error = function(e) {    
    message(paste("Error reading file:", x, "\n", e$message))  
    return(NA)  })
})
#MAF
mafdata_clean1 <- mafdata[!sapply(mafdata, is.na)]
#maf
merger_maf=merge_mafs(mafdata_clean1)
save(merger_maf,file=c("TCGA-HNSC_SNP.Rdata"))#rda

class(merger_maf)
slotNames(merger_maf)
maf_cluster<-data.frame(sample=substr(gsub('[.]','-',rownames(clusters)),1,12),cluster=clusters$clusters)



# MAF
plotmafSummary(maf = merger_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
merger_maf_1 <- subsetMaf(maf = merger_maf, tsb = maf_cluster$sample[which(maf_cluster$cluster == "1")])
merger_maf_2 <- subsetMaf(maf = merger_maf, tsb = maf_cluster$sample[which(maf_cluster$cluster == "2")])
pdf("G:\\pan_result\\maftools\\HNSC_cluster1.pdf")
plotmafSummary(maf = merger_maf_1, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
dev.off()
pdf("G:\\pan_result\\maftools\\HNSC_cluster2.pdf")
plotmafSummary(maf = merger_maf_2, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
dev.off()
class(merger_maf)
showMethods(classes = class(merger_maf))
merger_maf@variants.per.sample

#samplevariants
variants_per_sample <- merger_maf@variants.per.sample
variants_per_sample1<- merger_maf_1@variants.per.sample
variants_per_sample2<- merger_maf_2@variants.per.sample
variants_per_sample<-as.data.frame()
t.test(variants_per_sample1$Variants,variants_per_sample2$Variants)



#tmb
HNSC.tmb1 <- tmb(merger_maf_1, captureSize = 38, logScale = T)
HNSC.tmb2 <- tmb(merger_maf_2, captureSize = 38, logScale = T)
wilcox.test(HNSC.tmb1$total_perMB,HNSC.tmb2$total_perMB)
HNSC.tmb1$cluster<-"1"
HNSC.tmb2$cluster<-"2"
HNSC.tmb<-rbind(HNSC.tmb1,HNSC.tmb2)
colnames(HNSC.tmb)[3]<-"TMB"
HNSC.tmb <- HNSC.tmb %>%
  filter(HNSC.tmb$TMB < quantile(HNSC.tmb$TMB, 0.95)) # 95%


p <- ggplot(HNSC.tmb, aes(x = cluster, y = TMB, color = cluster)) +
  geom_boxplot(width=0.5,size=2, alpha = 0.2, notch = TRUE, notchwidth = 0.5) + # 
  scale_color_manual(values = c("#FFD4A9", "#80D1C8", "pink")) + # 
  theme_bw() + # 
  theme(
    panel.grid = element_blank(), # 
    panel.border = element_blank(), # 
    axis.line = element_line(colour = "black") # 
  ) +
  ylab("HRD score") + # y
  xlab(NULL) + # x
  guides(scale = "none") + # 
  geom_signif(comparisons = list(c("1", "2")),
              map_signif_level = TRUE,
              textsize = 8,color="black") + # 
  
  geom_jitter(
    width = 0.1, # 
    size = 1.5, # 
    aes(fill = cluster), # 
    shape = 21, # 
    alpha = 0.6 # 
  )

print(p)
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_TMB_plot.pdf",p,width = 10,height = 10)

#cluster1 cluster2####
length(which(DEG_all$FC>=1))
length(which(DEG_all$FC<1))


#KIRC####
#1.12####
# R
# library(TCGAbiolinks)
# 
# query <- GDCquery(
#   project = "TCGA-KIRC", 
#   data.category = "Simple Nucleotide Variation",
#   data.type = "Masked Somatic Mutation",
#   access = "open"
# )
# 
# GDCdownload(query)
# 
# GDCprepare(query, save = T,save.filename = "TCGA-KIRC_SNP.Rdata") # Rdatamaftools
# setwd("G:\\pan_result\\maftools")
# 
# cancer_name<-toupper(c("brca","KIRC","KIRC","hnsc","kirc","kirp","lihc","lung","paad","thca","ucec","LUAD","LUSC","paad","thca","ucec"))
# # cancer_name<-toupper(c("LUAD","LUSC","paad","thca","ucec"))
# 
# # name<-cancer_name[1]
# for (name in cancer_name) {
#   query <- GDCquery(
#     project = paste0("TCGA-",name),
#     data.category = "Simple Nucleotide Variation",
#     data.type = "Masked Somatic Mutation",
#     access = "open"
#   )
#   
#   GDCdownload(query)
#   
#  # GDCprepare(query, save = T,save.filename =paste0("TCGA-",name,"_SNP.Rdata")) # Rdatamaftools
#   
# }
# 
# sessionInfo()

setwd("G:\\pan_result\\maftools\\GDCdata\\TCGA-KIRC\\Simple_Nucleotide_Variation")
library(tidyverse)
# devtools::install_github("PoisonAlien/maftools")
library(maftools)
#dir"maf_data/"
maffilepath=dir(path = "Masked_Somatic_Mutation\\", pattern="masked.maf.gz$", full.names = T, recursive=T)#dir(path = "maf_data/") ,"masked.maf.gz$"masked.maf.gz
head(maffilepath)
#mafdata  MafFrame  MAF  NA
mafdata <- lapply(maffilepath, function(x) {  
  tryCatch({    read.maf(x, isTCGA = TRUE)  
  }, error = function(e) {    
    message(paste("Error reading file:", x, "\n", e$message))  
    return(NA)  })
})
#MAF
mafdata_clean1 <- mafdata[!sapply(mafdata, is.na)]
#maf
merger_maf=merge_mafs(mafdata_clean1)
save(merger_maf,file=c("TCGA-KIRC_SNP.Rdata"))#rda

class(merger_maf)
slotNames(merger_maf)
maf_cluster<-data.frame(sample=substr(gsub('[.]','-',rownames(clusters)),1,12),cluster=clusters$clusters)



# MAF
plotmafSummary(maf = merger_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
merger_maf_1 <- subsetMaf(maf = merger_maf, tsb = maf_cluster$sample[which(maf_cluster$cluster == "1")])
merger_maf_2 <- subsetMaf(maf = merger_maf, tsb = maf_cluster$sample[which(maf_cluster$cluster == "2")])
merger_maf_3 <- subsetMaf(maf = merger_maf, tsb = maf_cluster$sample[which(maf_cluster$cluster == "3")])
pdf("G:\\pan_result\\maftools\\KIRC_cluster1.pdf")
plotmafSummary(maf = merger_maf_1, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
dev.off()
pdf("G:\\pan_result\\maftools\\KIRC_cluster2.pdf")
plotmafSummary(maf = merger_maf_2, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
dev.off()
pdf("G:\\pan_result\\maftools\\KIRC_cluster3.pdf")
plotmafSummary(maf = merger_maf_3, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
dev.off()

class(merger_maf)
showMethods(classes = class(merger_maf))
merger_maf@variants.per.sample

#samplevariants
variants_per_sample <- merger_maf@variants.per.sample
variants_per_sample1<- merger_maf_1@variants.per.sample
variants_per_sample2<- merger_maf_2@variants.per.sample
variants_per_sample<-as.data.frame()
t.test(variants_per_sample1$Variants,variants_per_sample2$Variants)



#tmb####
KIRC.tmb1 <- tmb(merger_maf_1, captureSize = 38, logScale = T)
KIRC.tmb2 <- tmb(merger_maf_2, captureSize = 38, logScale = T)
KIRC.tmb3 <- tmb(merger_maf_3, captureSize = 38, logScale = T)
cluster_tmb<-rbind(KIRC.tmb1,KIRC.tmb2,KIRC.tmb3)
mianyi_class<-rep(c(1,2,3),c(nrow(KIRC.tmb1),nrow(KIRC.tmb2),nrow(KIRC.tmb3)))
kruskal.test(as.numeric(cluster_tmb$total_perMB)~mianyi_class)

KIRC.tmb1$cluster<-"1"
KIRC.tmb2$cluster<-"2"
KIRC.tmb3$cluster<-"3"
KIRC.tmb<-rbind(KIRC.tmb1,KIRC.tmb2,KIRC.tmb3)
colnames(KIRC.tmb)[3]<-"TMB"
KIRC.tmb$cluster<-as.factor(KIRC.tmb$cluster)
str(KIRC.tmb$cluster)
class(KIRC.tmb$cluster)

kruskal.test(as.numeric(KIRC.tmb$TMB) ~ KIRC.tmb$cluster)

# p <- ggplot(KIRC.tmb, aes(x = cluster, y = TMB, color = cluster)) +
#   geom_boxplot(width=0.5,size=2, alpha = 0.2, notch = TRUE, notchwidth = 0.5) + # 
#   scale_color_manual(values = c("#FFD4A9", "#80D1C8", "pink")) + # 
#   theme_bw() + # 
#   theme(
#     panel.grid = element_blank(), # 
#     panel.border = element_blank(), # 
#     axis.line = element_line(colour = "black") # 
#   ) +
#   ylab("HRD score") + # y
#   xlab(NULL) + # x
#   guides(scale = "none") + # 
#   geom_signif(comparisons = list(c("1", "2","3")),
#               map_signif_level = TRUE,
#               textsize = 8,color="black") + # 
#   
#   geom_jitter(
#     width = 0.1, # 
#     size = 1.5, # 
#     aes(fill = cluster), # 
#     shape = 21, # 
#     alpha = 0.6 # 
#   )

#  Kruskal-Wallis 
kruskal_test_result <- kruskal.test(as.numeric(KIRC.tmb$total_perMB) ~ KIRC.tmb$cluster)


#  p 
p_value <- kruskal_test_result$p.value
p <- KIRC.tmb <- KIRC.tmb %>%
  filter(KIRC.tmb$TMB < quantile(KIRC.tmb$TMB, 0.95)) # 95%
# 
ggplot(KIRC.tmb, aes(x = cluster, y = TMB, color = cluster)) +
  geom_boxplot(width = 0.5, size = 2, alpha = 0.2, notch = TRUE, notchwidth = 0.5) +  # 
  scale_color_manual(values = c("#FFD4A9", "#80D1C8", "pink")) +  # 
  theme_bw() +  # 
  theme(
    panel.grid = element_blank(),  # 
    panel.border = element_blank(),  # 
    axis.line = element_line(colour = "black")  # 
  ) +
  ylab("HRD score") +  # y
  xlab(NULL) +  # x
  guides(scale = "none") +  # 
  geom_signif(comparisons = list(c("1", "2", "3")), 
              annotations = paste("p =", format.pval(p_value)),  #  p 
              map_signif_level = TRUE,
              textsize = 8, color = "black") +  #  p 
  geom_jitter(
    width = 0.1,  # 
    size = 1.5,  # 
    aes(fill = cluster),  # 
    shape = 21,  # 
    alpha = 0.6  # 
  )

print(p)
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\KIRC_TMB_plot.pdf",p,width = 10,height = 10)
#####

data<-cbind(t(gene1),t(gene2),t(gene3))
class<-rep(c(1,2,3),c(nrow(gene1),nrow(gene2),nrow(gene3)))
T_p<-c()
FC_result<-c()
for(i in 1:nrow(data)){
  #T_result<-kruskal.test(as.numeric(data[i,])~class)
  #####################################t############
  #p1<-T_result$p.value#p
  # t1<-T_result[[1]]#T
  # T_p<-c(T_p,p1)
  #T_value<-c(T_value,t1)
  FC<-mean(as.numeric(data[i,(nrow(gene1)+nrow(gene2)+1):ncol(data)]))/mean(as.numeric(data[i,1:(nrow(gene1)+nrow(gene2))]))
  ##fold change########
  FC_result<-c(FC_result,FC)
}
length(which(FC_result<1))
#KIRP####
#1.12####
# R
# library(TCGAbiolinks)
# 
# query <- GDCquery(
#   project = "TCGA-KIRP", 
#   data.category = "Simple Nucleotide Variation",
#   data.type = "Masked Somatic Mutation",
#   access = "open"
# )
# 
# GDCdownload(query)
# 
# GDCprepare(query, save = T,save.filename = "TCGA-KIRP_SNP.Rdata") # Rdatamaftools
# setwd("G:\\pan_result\\maftools")
# 
# cancer_name<-toupper(c("brca","KIRP","esca","hnsc","kirc","kirp","lihc","lung","paad","thca","ucec","LUAD","LUSC","paad","thca","ucec"))
# # cancer_name<-toupper(c("LUAD","LUSC","paad","thca","ucec"))
# 
# # name<-cancer_name[1]
# for (name in cancer_name) {
#   query <- GDCquery(
#     project = paste0("TCGA-",name),
#     data.category = "Simple Nucleotide Variation",
#     data.type = "Masked Somatic Mutation",
#     access = "open"
#   )
#   
#   GDCdownload(query)
#   
#  # GDCprepare(query, save = T,save.filename =paste0("TCGA-",name,"_SNP.Rdata")) # Rdatamaftools
#   
# }
# 
# sessionInfo()

setwd("G:\\pan_result\\maftools\\GDCdata\\TCGA-KIRP\\Simple_Nucleotide_Variation")
library(tidyverse)
# devtools::install_github("PoisonAlien/maftools")
library(maftools)
#dir"maf_data/"
maffilepath=dir(path = "Masked_Somatic_Mutation\\", pattern="masked.maf.gz$", full.names = T, recursive=T)#dir(path = "maf_data/") ,"masked.maf.gz$"masked.maf.gz
head(maffilepath)
#mafdata  MafFrame  MAF  NA
mafdata <- lapply(maffilepath, function(x) {  
  tryCatch({    read.maf(x, isTCGA = TRUE)  
  }, error = function(e) {    
    message(paste("Error reading file:", x, "\n", e$message))  
    return(NA)  })
})
#MAF
mafdata_clean1 <- mafdata[!sapply(mafdata, is.na)]
#maf
merger_maf=merge_mafs(mafdata_clean1)
save(merger_maf,file=c("TCGA-KIRP_SNP.Rdata"))#rda

class(merger_maf)
slotNames(merger_maf)
maf_cluster<-data.frame(sample=substr(gsub('[.]','-',rownames(clusters)),1,12),cluster=clusters$clusters)



# MAF
plotmafSummary(maf = merger_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
merger_maf_1 <- subsetMaf(maf = merger_maf, tsb = maf_cluster$sample[which(maf_cluster$cluster == "1")])
merger_maf_2 <- subsetMaf(maf = merger_maf, tsb = maf_cluster$sample[which(maf_cluster$cluster == "2")])
pdf("G:\\pan_result\\maftools\\KIRP_cluster1.pdf")
plotmafSummary(maf = merger_maf_1, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
dev.off()
pdf("G:\\pan_result\\maftools\\KIRP_cluster2.pdf")
plotmafSummary(maf = merger_maf_2, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
dev.off()


class(merger_maf)
showMethods(classes = class(merger_maf))
merger_maf@variants.per.sample

#samplevariants
# variants_per_sample <- merger_maf@variants.per.sample
# variants_per_sample1<- merger_maf_1@variants.per.sample
# variants_per_sample2<- merger_maf_2@variants.per.sample
# variants_per_sample<-as.data.frame()
# t.test(variants_per_sample1$Variants,variants_per_sample2$Variants)
#tmb
KIRP.tmb1 <- tmb(merger_maf_1, captureSize = 38, logScale = T)
KIRP.tmb2 <- tmb(merger_maf_2, captureSize = 38, logScale = T)
wilcox.test(KIRP.tmb1$total_perMB,KIRP.tmb2$total_perMB)
KIRP.tmb1$cluster<-"1"
KIRP.tmb2$cluster<-"2"
KIRP.tmb<-rbind(KIRP.tmb1,KIRP.tmb2)
colnames(KIRP.tmb)[3]<-"TMB"

p <- ggplot(KIRP.tmb, aes(x = cluster, y = TMB, color = cluster)) +
  geom_boxplot(width=0.5,size=2, alpha = 0.2, notch = TRUE, notchwidth = 0.5) + # 
  scale_color_manual(values = c("#FFD4A9", "#80D1C8", "pink")) + # 
  theme_bw() + # 
  theme(
    panel.grid = element_blank(), # 
    panel.border = element_blank(), # 
    axis.line = element_line(colour = "black") # 
  ) +
  ylab("HRD score") + # y
  xlab(NULL) + # x
  guides(scale = "none") + # 
  geom_signif(comparisons = list(c("1", "2")),
              map_signif_level = TRUE,
              textsize = 8,color="black") + # 
  
  geom_jitter(
    width = 0.1, # 
    size = 1.5, # 
    aes(fill = cluster), # 
    shape = 21, # 
    alpha = 0.6 # 
  )

print(p)
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\KIRP_TMB_plot.pdf",p,width = 10,height = 10)

#LIHC####
#1.12####
# R
# setwd("G:\\pan_result\\maftools")
# library(TCGAbiolinks)
# 
# query <- GDCquery(
#   project = "TCGA-LIHC",
#   data.category = "Simple Nucleotide Variation",
#   data.type = "Masked Somatic Mutation",
#   access = "open"
# )
# 
# GDCdownload(query)
# 
# GDCprepare(query, save = T,save.filename = "TCGA-LIHC_SNP.Rdata") # Rdatamaftools
# setwd("G:\\pan_result\\maftools")
# 
# cancer_name<-toupper(c("brca","LIHC","LIHC","hnsc","kirc","kirp","lihc","lung","paad","thca","ucec","LUAD","LUSC","paad","thca","ucec"))
# # cancer_name<-toupper(c("LUAD","LUSC","paad","thca","ucec"))
# 
# # name<-cancer_name[1]
# for (name in cancer_name) {
#   query <- GDCquery(
#     project = paste0("TCGA-",name),
#     data.category = "Simple Nucleotide Variation",
#     data.type = "Masked Somatic Mutation",
#     access = "open"
#   )
#   
#   GDCdownload(query)
#   
#  # GDCprepare(query, save = T,save.filename =paste0("TCGA-",name,"_SNP.Rdata")) # Rdatamaftools
#   
# }
# 
# sessionInfo()

setwd("G:\\pan_result\\maftools\\GDCdata\\TCGA-LIHC\\Simple_Nucleotide_Variation")
library(tidyverse)
# devtools::install_github("PoisonAlien/maftools")
library(maftools)
#dir"maf_data/"
maffilepath=dir(path = "Masked_Somatic_Mutation\\", pattern="masked.maf.gz$", full.names = T, recursive=T)#dir(path = "maf_data/") ,"masked.maf.gz$"masked.maf.gz
head(maffilepath)
#mafdata  MafFrame  MAF  NA
mafdata <- lapply(maffilepath, function(x) {  
  tryCatch({    read.maf(x, isTCGA = TRUE)  
  }, error = function(e) {    
    message(paste("Error reading file:", x, "\n", e$message))  
    return(NA)  })
})
#MAF
mafdata_clean1 <- mafdata[!sapply(mafdata, is.na)]
#maf
merger_maf=merge_mafs(mafdata_clean1)
save(merger_maf,file=c("TCGA-LIHC_SNP.Rdata"))#rda

class(merger_maf)
slotNames(merger_maf)
maf_cluster<-data.frame(sample=substr(gsub('[.]','-',rownames(clusters)),1,12),cluster=clusters$clusters)



# MAF
plotmafSummary(maf = merger_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
merger_maf_1 <- subsetMaf(maf = merger_maf, tsb = maf_cluster$sample[which(maf_cluster$cluster == "1")])
merger_maf_2 <- subsetMaf(maf = merger_maf, tsb = maf_cluster$sample[which(maf_cluster$cluster == "2")])
pdf("G:\\pan_result\\maftools\\LIHC_cluster1.pdf")
plotmafSummary(maf = merger_maf_1, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
dev.off()
pdf("G:\\pan_result\\maftools\\LIHC_cluster2.pdf")
plotmafSummary(maf = merger_maf_2, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
dev.off()
class(merger_maf)
showMethods(classes = class(merger_maf))
merger_maf@variants.per.sample

#samplevariants
variants_per_sample <- merger_maf@variants.per.sample
variants_per_sample1<- merger_maf_1@variants.per.sample
variants_per_sample2<- merger_maf_2@variants.per.sample
variants_per_sample<-as.data.frame()
t.test(variants_per_sample1$Variants,variants_per_sample2$Variants)




#tmb
LIHC.tmb1 <- tmb(merger_maf_1, captureSize = 38, logScale = T)
LIHC.tmb2 <- tmb(merger_maf_2, captureSize = 38, logScale = T)
wilcox.test(LIHC.tmb1$total_perMB,LIHC.tmb2$total_perMB)
LIHC.tmb1$cluster<-"1"
LIHC.tmb2$cluster<-"2"
LIHC.tmb<-rbind(LIHC.tmb1,LIHC.tmb2)
colnames(LIHC.tmb)[3]<-"TMB"
LIHC.tmb <- LIHC.tmb %>%
  filter(LIHC.tmb$TMB < quantile(LIHC.tmb$TMB, 0.95)) # 95%
p <- ggplot(LIHC.tmb, aes(x = cluster, y = TMB, color = cluster)) +
  geom_boxplot(width=0.5,size=2, alpha = 0.2, notch = TRUE, notchwidth = 0.5) + # 
  scale_color_manual(values = c("#FFD4A9", "#80D1C8", "pink")) + # 
  theme_bw() + # 
  theme(
    panel.grid = element_blank(), # 
    panel.border = element_blank(), # 
    axis.line = element_line(colour = "black") # 
  ) +
  ylab("HRD score") + # y
  xlab(NULL) + # x
  guides(scale = "none") + # 
  geom_signif(comparisons = list(c("1", "2")),
              map_signif_level = TRUE,
              textsize = 8,color="black") + # 
  
  geom_jitter(
    width = 0.1, # 
    size = 1.5, # 
    aes(fill = cluster), # 
    shape = 21, # 
    alpha = 0.6 # 
  )

print(p)
ggsave("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\LIHC_TMB_plot.pdf",p,width = 10,height = 10)

#LUNG####
#1.12####
# R
# setwd("G:\\pan_result\\maftools")
# library(TCGAbiolinks)
# 
# query <- GDCquery(
#   project = "TCGA-NSCLC",
#   data.category = "Simple Nucleotide Variation",
#   data.type = "Masked Somatic Mutation",
#   access = "open"
# )
# 
# GDCdownload(query)
# 
# GDCprepare(query, save = T,save.filename = "TCGA-NSCLC_SNP.Rdata") # Rdatamaftools
# setwd("G:\\pan_result\\maftools")
# 
# cancer_name<-toupper(c("brca","NSCLC","NSCLC","hnsc","kirc","kirp","NSCLC","NSCLC","paad","NSCLC","ucec","LUAD","LUSC","paad","NSCLC","ucec"))
# # cancer_name<-toupper(c("LUAD","LUSC","paad","NSCLC","ucec"))
# 
# # name<-cancer_name[1]
# for (name in cancer_name) {
# query <- GDCquery(
#   project = paste0("TCGA-",NSCLC),
#   data.category = "Simple Nucleotide Variation",
#   data.type = "Masked Somatic Mutation",
#   access = "open"
# )
#   
#   GDCdownload(query)
#   
#  # GDCprepare(query, save = T,save.filename =paste0("TCGA-",name,"_SNP.Rdata")) # Rdatamaftools
#   
# }
# 
# sessionInfo()

setwd("G:\\pan_result\\maftools\\GDCdata\\TCGA-LUSC\\Simple_Nucleotide_Variation")
library(tidyverse)
# devtools::install_github("PoisonAlien/maftools")
library(maftools)
#dir"maf_data/"
maffilepath=dir(path = "Masked_Somatic_Mutation\\", pattern="masked.maf.gz$", full.names = T, recursive=T)#dir(path = "maf_data/") ,"masked.maf.gz$"masked.maf.gz
head(maffilepath)
#mafdata  MafFrame  MAF  NA
mafdata <- lapply(maffilepath, function(x) {  
  tryCatch({    read.maf(x, isTCGA = TRUE)  
  }, error = function(e) {    
    message(paste("Error reading file:", x, "\n", e$message))  
    return(NA)  })
})
#MAF
mafdata_clean1 <- mafdata[!sapply(mafdata, is.na)]
#maf
merger_maf=merge_mafs(mafdata_clean1)
save(merger_maf,file=c("TCGA-LUSC_SNP.Rdata"))#rda

class(merger_maf)
slotNames(merger_maf)
maf_cluster<-data.frame(sample=substr(gsub('[.]','-',rownames(clusters)),1,12),cluster=clusters$clusters)


load("G:\\pan_result\\maftools\\GDCdata\\TCGA-LUAD\\Simple_Nucleotide_Variation\\TCGA-LUAD_SNP.Rdata")
maf1<-merger_maf
load("G:\\pan_result\\maftools\\GDCdata\\TCGA-LUSC\\Simple_Nucleotide_Variation\\TCGA-LUSC_SNP.Rdata")
maf2<-merger_maf
# MAF
maf_cluster <- merge_mafs(maf1, maf2)
#  maf1  maf2  MAF 
merged_data <- rbind(maf1@data, maf2@data)

#  MAF 
maf_cluster <- MAF(merged_data)



library(maftools)
library(dplyr)
# 
# maf1maf2MAF
# 
# maf1maf2
mutations_maf1 <- maf1@data
mutations_maf2 <- maf2@data

# 
head(mutations_maf1)
head(mutations_maf2)

# 
combined_mutations <- bind_rows(mutations_maf1, mutations_maf2)

# MAF
# ：combined_mutationsMAF
merged_maf <- read.maf(data = combined_mutations)
merged_maf <-maf_cluster
# MAF
print(merged_maf)

# MAF
plotmafSummary(maf = merger_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
merger_maf_1 <- subsetMaf(maf = merged_maf, tsb = maf_cluster$sample[which(maf_cluster$cluster == "1")])
merger_maf_11 <- subsetMaf(maf = merger_maf, tsb = maf_cluster$sample[which(maf_cluster$cluster == "1")])



merger_maf_2 <- subsetMaf(maf = maf1, tsb = maf_cluster$sample[which(maf_cluster$cluster == "2")])

#xin
merger_maf_1 <- subsetMaf(maf = merged_maf, tsb = clusters$bcr_patient_barcode[which(clusters$clusters == "1")])





pdf("G:\\pan_result\\maftools\\NSCLC_cluster1.pdf")
plotmafSummary(maf = merger_maf_1, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
dev.off()
pdf("G:\\pan_result\\maftools\\NSCLC_cluster2.pdf")
plotmafSummary(maf = merger_maf_2, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
dev.off()
class(merger_maf)
showMethods(classes = class(merger_maf))
merger_maf@variants.per.sample

#samplevariants
variants_per_sample <- merger_maf@variants.per.sample
variants_per_sample1<- maf1@variants.per.sample
variants_per_sample2<- maf2@variants.per.sample
variants_per_sample<-as.data.frame()
t.test(variants_per_sample1$Variants,variants_per_sample2$Variants)

#tmb
NSCLC.tmb1 <- tmb(maf1, captureSize = 38, logScale = T)
NSCLC.tmb2 <- tmb(maf2, captureSize = 38, logScale = T)
NSCLC.tmbhebing<-rbind(NSCLC.tmb1,NSCLC.tmb2)
NSCLC.tmb1<-NSCLC.tmbhebing[NSCLC.tmbhebing$Tumor_Sample_Barcode %in% c(clusters$bcr_patient_barcode[which(clusters$clusters=="1")]),]

NSCLC.tmb2<-NSCLC.tmbhebing[NSCLC.tmbhebing$Tumor_Sample_Barcode %in% c(clusters$bcr_patient_barcode[which(clusters$clusters=="2")]),]

wilcox.test(NSCLC.tmb1$total_perMB,NSCLC.tmb2$total_perMB)
NSCLC.tmb1$cluster<-"1"
NSCLC.tmb2$cluster<-"2"
NSCLC.tmb<-rbind(NSCLC.tmb1,NSCLC.tmb2)
colnames(NSCLC.tmb)[3]<-"TMB"

p <- ggplot(NSCLC.tmb, aes(x = cluster, y = TMB, color = cluster)) +
  geom_boxplot(width=0.5,size=2, alpha = 0.2, notch = TRUE, notchwidth = 0.5) + # 
  scale_color_manual(values = c("#FFD4A9", "#80D1C8", "pink")) + # 
  theme_bw() + # 
  theme(
    panel.grid = element_blank(), # 
    panel.border = element_blank(), # 
    axis.line = element_line(colour = "black") # 
  ) +
  ylab("HRD score") + # y
  xlab(NULL) + # x
  guides(scale = "none") + # 
  geom_signif(comparisons = list(c("1", "2")),
              map_signif_level = TRUE,
              textsize = 8,color="black") + # 
  
  geom_jitter(
    width = 0.1, # 
    size = 1.5, # 
    aes(fill = cluster), # 
    shape = 21, # 
    alpha = 0.6 # 
  )

print(p)
ggsave("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\NSCLC_TMB_plot.pdf",p,width = 10,height = 10)

#PAAD####
#1.12####
# R
# setwd("G:\\pan_result\\maftools")
# library(TCGAbiolinks)
# 
# query <- GDCquery(
#   project = "TCGA-PAAD",
#   data.category = "Simple Nucleotide Variation",
#   data.type = "Masked Somatic Mutation",
#   access = "open"
# )
# 
# GDCdownload(query)
# 
# GDCprepare(query, save = T,save.filename = "TCGA-PAAD_SNP.Rdata") # Rdatamaftools
# setwd("G:\\pan_result\\maftools")
# 
# cancer_name<-toupper(c("brca","PAAD","PAAD","hnsc","kirc","kirp","PAAD","lung","paad","PAAD","ucec","LUAD","LUSC","paad","PAAD","ucec"))
# # cancer_name<-toupper(c("LUAD","LUSC","paad","PAAD","ucec"))
# 
# # name<-cancer_name[1]
# for (name in cancer_name) {
#   query <- GDCquery(
#     project = paste0("TCGA-",name),
#     data.category = "Simple Nucleotide Variation",
#     data.type = "Masked Somatic Mutation",
#     access = "open"
#   )
#   
#   GDCdownload(query)
#   
#  # GDCprepare(query, save = T,save.filename =paste0("TCGA-",name,"_SNP.Rdata")) # Rdatamaftools
#   
# }
# 
# sessionInfo()

setwd("G:\\pan_result\\maftools\\GDCdata\\TCGA-PAAD\\Simple_Nucleotide_Variation")
library(tidyverse)
# devtools::install_github("PoisonAlien/maftools")
library(maftools)
#dir"maf_data/"
maffilepath=dir(path = "Masked_Somatic_Mutation\\", pattern="masked.maf.gz$", full.names = T, recursive=T)#dir(path = "maf_data/") ,"masked.maf.gz$"masked.maf.gz
head(maffilepath)
#mafdata  MafFrame  MAF  NA
mafdata <- lapply(maffilepath, function(x) {  
  tryCatch({    read.maf(x, isTCGA = TRUE)  
  }, error = function(e) {    
    message(paste("Error reading file:", x, "\n", e$message))  
    return(NA)  })
})
#MAF
mafdata_clean1 <- mafdata[!sapply(mafdata, is.na)]
#maf
merger_maf=merge_mafs(mafdata_clean1)
save(merger_maf,file=c("TCGA-PAAD_SNP.Rdata"))#rda

class(merger_maf)
slotNames(merger_maf)
maf_cluster<-data.frame(sample=substr(gsub('[.]','-',rownames(clusters)),1,12),cluster=clusters$clusters)



# MAF
plotmafSummary(maf = merger_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
merger_maf_1 <- subsetMaf(maf = merger_maf, tsb = maf_cluster$sample[which(maf_cluster$cluster == "1")])
merger_maf_2 <- subsetMaf(maf = merger_maf, tsb = maf_cluster$sample[which(maf_cluster$cluster == "2")])
pdf("G:\\pan_result\\maftools\\PAAD_cluster1.pdf")
plotmafSummary(maf = merger_maf_1, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
dev.off()
pdf("G:\\pan_result\\maftools\\PAAD_cluster2.pdf")
plotmafSummary(maf = merger_maf_2, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
dev.off()
class(merger_maf)
showMethods(classes = class(merger_maf))
merger_maf@variants.per.sample

#samplevariants
variants_per_sample <- merger_maf@variants.per.sample
variants_per_sample1<- merger_maf_1@variants.per.sample
variants_per_sample2<- merger_maf_2@variants.per.sample
variants_per_sample<-as.data.frame()
t.test(variants_per_sample1$Variants,variants_per_sample2$Variants)
#tmb
PAAD.tmb1 <- tmb(merger_maf_1, captureSize = 38, logScale = T)
PAAD.tmb2 <- tmb(merger_maf_2, captureSize = 38, logScale = T)
wilcox.test(PAAD.tmb1$total_perMB,PAAD.tmb2$total_perMB)
PAAD.tmb1$cluster<-"1"
PAAD.tmb2$cluster<-"2"
PAAD.tmb<-rbind(PAAD.tmb1,PAAD.tmb2)
colnames(PAAD.tmb)[3]<-"TMB"
str(PAAD.tmb) # 
PAAD.tmb <- PAAD.tmb %>%
  filter(PAAD.tmb$TMB < quantile(PAAD.tmb$TMB, 0.95)) # 95%


p <- ggplot(PAAD.tmb, aes(x = cluster, y = TMB, color = cluster)) +
  geom_boxplot(width=0.5,size=2, alpha = 0.2, notch = TRUE, notchwidth = 0.5) + # 
  scale_color_manual(values = c("#FFD4A9", "#80D1C8", "pink")) + # 
  theme_bw() + # 
  theme(
    panel.grid = element_blank(), # 
    panel.border = element_blank(), # 
    axis.line = element_line(colour = "black") # 
  ) +
  ylab("HRD score") + # y
  xlab(NULL) + # x
  guides(scale = "none") + # 
  geom_signif(comparisons = list(c("1", "2")),
              map_signif_level = TRUE,
              textsize = 8,color="black") + # 
  
  geom_jitter(
    size = 1.5, # 
    aes(fill = cluster), # 
    shape = 21, # 
    alpha = 0.6 # 
  )

print(p)
ggsave("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\PAAD_TMB_plot.pdf",p,width = 10,height = 10)


#cluster1 cluster2####
length(which(DEG_all$FC>=1))
length(which(DEG_all$FC<1))
#THCA####
#1.12####
# R
# setwd("G:\\pan_result\\maftools")
# library(TCGAbiolinks)
# 
# query <- GDCquery(
#   project = "TCGA-THCA",
#   data.category = "Simple Nucleotide Variation",
#   data.type = "Masked Somatic Mutation",
#   access = "open"
# )
# 
# GDCdownload(query)
# 
# GDCprepare(query, save = T,save.filename = "TCGA-THCA_SNP.Rdata") # Rdatamaftools
# setwd("G:\\pan_result\\maftools")
# 
# cancer_name<-toupper(c("brca","THCA","THCA","hnsc","kirc","kirp","THCA","THCA","paad","thca","ucec","LUAD","LUSC","paad","thca","ucec"))
# # cancer_name<-toupper(c("LUAD","LUSC","paad","thca","ucec"))
# 
# # name<-cancer_name[1]
# for (name in cancer_name) {
#   query <- GDCquery(
#     project = paste0("TCGA-",name),
#     data.category = "Simple Nucleotide Variation",
#     data.type = "Masked Somatic Mutation",
#     access = "open"
#   )
#   
#   GDCdownload(query)
#   
#  # GDCprepare(query, save = T,save.filename =paste0("TCGA-",name,"_SNP.Rdata")) # Rdatamaftools
#   
# }
# 
# sessionInfo()

setwd("G:\\pan_result\\maftools\\GDCdata\\TCGA-THCA\\Simple_Nucleotide_Variation")
library(tidyverse)
# devtools::install_github("PoisonAlien/maftools")
library(maftools)
#dir"maf_data/"
maffilepath=dir(path = "Masked_Somatic_Mutation\\", pattern="masked.maf.gz$", full.names = T, recursive=T)#dir(path = "maf_data/") ,"masked.maf.gz$"masked.maf.gz
head(maffilepath)
#mafdata  MafFrame  MAF  NA
mafdata <- lapply(maffilepath, function(x) {  
  tryCatch({    read.maf(x, isTCGA = TRUE)  
  }, error = function(e) {    
    message(paste("Error reading file:", x, "\n", e$message))  
    return(NA)  })
})
#MAF
mafdata_clean1 <- mafdata[!sapply(mafdata, is.na)]
#maf
merger_maf=merge_mafs(mafdata_clean1)
save(merger_maf,file=c("TCGA-THCA_SNP.Rdata"))#rda

class(merger_maf)
slotNames(merger_maf)
maf_cluster<-data.frame(sample=substr(gsub('[.]','-',rownames(clusters)),1,12),cluster=clusters$clusters)



# MAF
plotmafSummary(maf = merger_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
merger_maf_1 <- subsetMaf(maf = merger_maf, tsb = maf_cluster$sample[which(maf_cluster$cluster == "1")])
merger_maf_2 <- subsetMaf(maf = merger_maf, tsb = maf_cluster$sample[which(maf_cluster$cluster == "2")])
pdf("G:\\pan_result\\maftools\\THCA_cluster1.pdf")
plotmafSummary(maf = merger_maf_1, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
dev.off()
pdf("G:\\pan_result\\maftools\\THCA_cluster2.pdf")
plotmafSummary(maf = merger_maf_2, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
dev.off()
class(merger_maf)
showMethods(classes = class(merger_maf))
merger_maf@variants.per.sample

#samplevariants
variants_per_sample <- merger_maf@variants.per.sample
variants_per_sample1<- merger_maf_1@variants.per.sample
variants_per_sample2<- merger_maf_2@variants.per.sample
variants_per_sample<-as.data.frame()
t.test(variants_per_sample1$Variants,variants_per_sample2$Variants)


#tmb
THCA.tmb1 <- tmb(merger_maf_1, captureSize = 38, logScale = T)
THCA.tmb2 <- tmb(merger_maf_2, captureSize = 38, logScale = T)
wilcox.test(THCA.tmb1$total_perMB,THCA.tmb2$total_perMB)
THCA.tmb1$cluster<-"1"
THCA.tmb2$cluster<-"2"
THCA.tmb<-rbind(THCA.tmb1,THCA.tmb2)
colnames(THCA.tmb)[3]<-"TMB"
THCA.tmb <- THCA.tmb %>%
  filter(THCA.tmb$TMB < quantile(THCA.tmb$TMB, 0.95)) # 95%
p <- ggplot(THCA.tmb, aes(x = cluster, y = TMB, color = cluster)) +
  geom_boxplot(width=0.5,size=2, alpha = 0.2, notch = TRUE, notchwidth = 0.5) + # 
  scale_color_manual(values = c("#FFD4A9", "#80D1C8", "pink")) + # 
  theme_bw() + # 
  theme(
    panel.grid = element_blank(), # 
    panel.border = element_blank(), # 
    axis.line = element_line(colour = "black") # 
  ) +
  ylab("HRD score") + # y
  xlab(NULL) + # x
  guides(scale = "none") + # 
  geom_signif(comparisons = list(c("1", "2")),
              map_signif_level = TRUE,
              textsize = 8,color="black") + # 
  
  geom_jitter(
    width = 0.1, # 
    size = 1.5, # 
    aes(fill = cluster), # 
    shape = 21, # 
    alpha = 0.6 # 
  )

print(p)
ggsave("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\THCA_TMB_plot.pdf",p,width = 10,height = 10)

#UCEC####
#1.12####
# R
# setwd("G:\\pan_result\\maftools")
# library(TCGAbiolinks)
# 
# query <- GDCquery(
#   project = "TCGA-UCEC",
#   data.category = "Simple Nucleotide Variation",
#   data.type = "Masked Somatic Mutation",
#   access = "open"
# )
# 
# GDCdownload(query)
# 
# GDCprepare(query, save = T,save.filename = "TCGA-UCEC_SNP.Rdata") # Rdatamaftools
# setwd("G:\\pan_result\\maftools")
# 
# cancer_name<-toupper(c("brca","UCEC","UCEC","hnsc","kirc","kirp","UCEC","lung","paad","UCEC","ucec","LUAD","LUSC","paad","UCEC","ucec"))
# # cancer_name<-toupper(c("LUAD","LUSC","paad","UCEC","ucec"))
# 
# # name<-cancer_name[1]
# for (name in cancer_name) {
#   query <- GDCquery(
#     project = paste0("TCGA-",name),
#     data.category = "Simple Nucleotide Variation",
#     data.type = "Masked Somatic Mutation",
#     access = "open"
#   )
#   
#   GDCdownload(query)
#   
#  # GDCprepare(query, save = T,save.filename =paste0("TCGA-",name,"_SNP.Rdata")) # Rdatamaftools
#   
# }
# 
# sessionInfo()

setwd("G:\\pan_result\\maftools\\GDCdata\\TCGA-UCEC\\Simple_Nucleotide_Variation")
library(tidyverse)
# devtools::install_github("PoisonAlien/maftools")
library(maftools)
#dir"maf_data/"
maffilepath=dir(path = "Masked_Somatic_Mutation\\", pattern="masked.maf.gz$", full.names = T, recursive=T)#dir(path = "maf_data/") ,"masked.maf.gz$"masked.maf.gz
head(maffilepath)
#mafdata  MafFrame  MAF  NA
mafdata <- lapply(maffilepath, function(x) {  
  tryCatch({    read.maf(x, isTCGA = TRUE)  
  }, error = function(e) {    
    message(paste("Error reading file:", x, "\n", e$message))  
    return(NA)  })
})
#MAF
mafdata_clean1 <- mafdata[!sapply(mafdata, is.na)]
#maf
merger_maf=merge_mafs(mafdata_clean1)
save(merger_maf,file=c("TCGA-UCEC_SNP.Rdata"))#rda

class(merger_maf)
slotNames(merger_maf)
maf_cluster<-data.frame(sample=substr(gsub('[.]','-',rownames(clusters)),1,12),cluster=clusters$clusters)



# MAF
plotmafSummary(maf = merger_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
merger_maf_1 <- subsetMaf(maf = merger_maf, tsb = maf_cluster$sample[which(maf_cluster$cluster == "1")])
merger_maf_2 <- subsetMaf(maf = merger_maf, tsb = maf_cluster$sample[which(maf_cluster$cluster == "2")])
pdf("G:\\pan_result\\maftools\\UCEC_cluster1.pdf")
plotmafSummary(maf = merger_maf_1, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
dev.off()
pdf("G:\\pan_result\\maftools\\UCEC_cluster2.pdf")
plotmafSummary(maf = merger_maf_2, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
dev.off()
class(merger_maf)
showMethods(classes = class(merger_maf))
merger_maf@variants.per.sample

#samplevariants
variants_per_sample <- merger_maf@variants.per.sample
variants_per_sample1<- merger_maf_1@variants.per.sample
variants_per_sample2<- merger_maf_2@variants.per.sample
variants_per_sample<-as.data.frame()
t.test(variants_per_sample1$Variants,variants_per_sample2$Variants)

#tmb
UCEC.tmb1 <- tmb(merger_maf_1, captureSize = 38, logScale = T)
UCEC.tmb2 <- tmb(merger_maf_2, captureSize = 38, logScale = T)
wilcox.test(UCEC.tmb1$total_perMB,UCEC.tmb2$total_perMB)
UCEC.tmb1$cluster<-"1"
UCEC.tmb2$cluster<-"2"
UCEC.tmb<-rbind(UCEC.tmb1,UCEC.tmb2)
colnames(UCEC.tmb)[3]<-"TMB"
UCEC.tmb <- UCEC.tmb %>%
  filter(UCEC.tmb$TMB < quantile(UCEC.tmb$TMB, 0.9)) # 95%
p <- ggplot(UCEC.tmb, aes(x = cluster, y = TMB, color = cluster)) +
  geom_boxplot(width=0.5,size=2, alpha = 0.2, notch = TRUE, notchwidth = 0.5) + # 
  scale_color_manual(values = c("#FFD4A9", "#80D1C8", "pink")) + # 
  theme_bw() + # 
  theme(
    panel.grid = element_blank(), # 
    panel.border = element_blank(), # 
    axis.line = element_line(colour = "black") # 
  ) +
  ylab("HRD score") + # y
  xlab(NULL) + # x
  guides(scale = "none") + # 
  geom_signif(comparisons = list(c("1", "2")),
              map_signif_level = TRUE,
              textsize = 8,color="black") + # 
  
  geom_jitter(
    width = 0.1, # 
    size = 1.5, # 
    aes(fill = cluster), # 
    shape = 21, # 
    alpha = 0.6 # 
  )

print(p)
ggsave("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\UCEC_TMB_plot.pdf",p,width = 10,height = 10)


