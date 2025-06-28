#BRCA####
#BiocManager::install("monocle3")
#install.packages("leidenbase")
#devtools::install_github("r-lib/rlang")
library(rlang)
#packageVersion("dplyr")
#install.packages("https://www.bioconductor.org/packages/3.8/bioc/src/contrib/BiocGenerics_0.28.0.tar.gz",repos=NULL,type="source")
library(monocle)
library(monocle)
library(Seurat)
pmbc<-readRDS("sample.obj.rds")
GetAssay(pmbc,assay = "RNA")
unique(pmbc@meta.data[["orig.ident"]])

Idents(pmbc) <- "orig.ident"
ob1 <- subset(x = pmbc, downsample = 300)

expr_matrix<-as(as.matrix(ob1@assays$RNA@counts),'sparseMatrix')
p_data<-ob1@meta.data
p_data$celltype<-ob1@active.ident
f_data<-data.frame(gene_short_name=row.names(ob1),row.names = row.names(ob1))
dim(expr_matrix)

pd<-new('AnnotatedDataFrame',data=p_data)
fd<-new('AnnotatedDataFrame',data=f_data)

cds<-newCellDataSet(expr_matrix,
                    phenoData=pd,
                    featureData=fd,
                    lowerDetectionLimit=0.5,
                    expressionFamily=negbinomial.size())

cds<-estimateSizeFactors(cds)
#1gc()
# memory.limit(123456)
# rm(expr_matrix)
# rm(f_data)
# rm(p_data)
cds<-estimateDispersions(cds)

cds<-detectGenes(cds,min_expr = 0.1)

print(head(fData(cds)))


expressed_genes = read.table("validategene.csv",header = T,sep = ",")
expressed_genes = expressed_genes$X0

cds <- setOrderingFilter(cds, expressed_genes)
pdf("train.ordergenes.pdf")
plot_ordering_genes(cds)
dev.off()

#cds2<-cds
cds<-cds2
cds <- reduceDimension(cds, max_components = 5,
                       method = 'DDRTree')##method
cds <- orderCells(cds)
pdf("train.monocle.state5.pdf",width = 7,height = 7)
plot_cell_trajectory(cds, color_by = "Pseudotime",size=1,show_backbone=TRUE)
dev.off()

save(cds,file = "cds.RData")

pdf("train.monocle.celltype5.pdf",width = 7,height = 7)
plot_cell_trajectory(cds,color_by = "cell_type",size=1,show_backbone = TRUE)+
  geom_point(aes(size = Pseudotime , color = cell_type)) + 
  scale_color_manual(breaks = c("X", "Y", "Z"), values=c("darkblue", "darkred", "orange","#4F7E57")) 
plot_cell_trajectory(cds,color_by = "cell_type",size=1,show_backbone = TRUE) + scale_color_manual( values=c("darkblue", "darkred", "orange","#4F7E57")) + theme(legend.position = "right")
dev.off()
?plot_cell_trajectory
library(ggpubr)
library(tidyr)
packageVersion("tidyr")

#install.packages("tidyr")

load("cds.RData")
df <- pData(cds) 
#View(df)
pdf("train.monocle.pseudotime.pdf",width = 7,height = 7)
ggplot(df, aes(Pseudotime, colour = cell_type, fill=cell_type)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2( )
dev.off()
write.csv(df,file = "pData.csv")

library(SingleCellExperiment)
library(Seurat)



# devtools::install_github("alanocallaghan/scater")
library(scater)
# data.combined.new@meta.data[["Pseudotime"]]<-cds@phenoData@data[["Pseudotime"]]
# sce_gene <- as.SingleCellExperiment(data.combined.new)
exp1<-read.csv("G:\\pan_result\\gpu_result_tcga\\BRCA\\sigmoid.csv")
setdiff(expressed_genes,colnames(exp1))
expressed_genes<-gsub("-",".",expressed_genes)
exp<-exp1[,expressed_genes]
# exp<-as.data.frame(t(exp))
gene<-exp

label<-read.csv("G:\\pan_result\\gpu_result_tcga\\BRCA\\label.csv")
label$X=substr(label$X,1,15)
gene1<-gene[label$X[which(label$label=="1")],]
gene2<-gene[label$X[which(label$label=="2")],]
data<-cbind(t(gene1),t(gene2))
T_p<-c()
FC_result<-c()
for(i in 1:nrow(data)){
  T_result<-wilcox.test(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]),as.numeric(data[i,1:nrow(gene1)]),paired = F,var.equal = T)
  p1<-T_result$p.value

  T_p<-c(T_p,p1)

  FC<-mean(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]))/mean(as.numeric(data[i,1:nrow(gene1)]))

  FC_result<-c(FC_result,FC)
}

FDR_p<-p.adjust(T_p,method="fdr",length(T_p))
DEG_all<-data.frame(rownames(data),T_p,FC_result,stringsAsFactors = FALSE)
colnames(DEG_all)<-c("Gene","pval","FC") 


DEG<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC>2),]#1295

expression_matrix <- exprs(cds)


gene_expression <- expression_matrix["TH", ]#2.340444e-09

table(gene_expression)
# gene_expression <- FetchData(cds, vars = "NPTX2", slot = "data")

cds@phenoData@data[["TH"]]<-gene_expression

getwd()
pdf("monocle.pseudotime_gene2.pdf",width = 7,height = 7)

plot_cell_trajectory(cds, color_by = "TH") +
  scale_color_gradient(low = "blue", high = "red") +   
  geom_point(aes(size = Pseudotime, color = TH), alpha = 0.6) +   
  theme_minimal() +
  theme(legend.position = "right")   

dev.off()
rownames(gene)<-exp1$X
#cluster####
label<-read.table("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\BRCA_cluster.txt")
gene1<-gene[rownames(label)[which(label$clusters=="1")],]
gene2<-gene[rownames(label)[which(label$clusters=="2")],]
data<-cbind(t(gene1),t(gene2))
T_p<-c()
FC_result<-c()
for(i in 1:nrow(data)){
  T_result<-wilcox.test(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]),as.numeric(data[i,1:nrow(gene1)]),paired = F,var.equal = T)
  
  p1<-T_result$p.value

  T_p<-c(T_p,p1)
  #T_value<-c(T_value,t1)
  FC<-mean(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]))/mean(as.numeric(data[i,1:nrow(gene1)]))

  FC_result<-c(FC_result,FC)
}


DEG_all<-data.frame(rownames(data),T_p,FC_result,stringsAsFactors = FALSE)

colnames(DEG_all)<-c("Gene","pval","FC") 

write.csv(DEG_all,file = "DEG_all.csv")
DEG<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC>1.3),]#1295
DEG_down<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC<0.5),]#1295

expression_matrix <- exprs(cds)


gene_expression <- expression_matrix["TYMP", ]

# gene_expression <- FetchData(cds, vars = "ZNF487", slot = "data")

cds@phenoData@data[["TYMP"]]<-gene_expression
gene_expression2<-0

DEG<-intersect(DEG$Gene,expression_matrix@Dimnames[[1]])
for (i in 1:length(DEG)) {
  gene_expression <- expression_matrix[DEG[i], ]
  gene_expression2<-gene_expression2+gene_expression
}

cds@phenoData@data[["upgene"]]<-gene_expression2


getwd()
pdf("monocle.pseudotime_upgene.pdf",width = 7,height = 7)

plot_cell_trajectory(cds, color_by = "upgene") +
  scale_color_gradient(low = "blue", high = "red") +   
  geom_point(aes(size = Pseudotime, color = upgene), alpha = 0.6) +   
  theme_minimal() +
  theme(legend.position = "right")   

dev.off()

#downgene####
gene_expression2<-0

DEG_down<-intersect(DEG_down$Gene,expression_matrix@Dimnames[[1]])
for (i in 1:length(DEG_down)) {
  gene_expression <- expression_matrix[DEG_down[i], ]
  gene_expression2<-gene_expression2+gene_expression
}

cds@phenoData@data[["downgene"]]<-gene_expression2




getwd()
pdf("monocle.pseudotime_downgene.pdf",width = 7,height = 7)

plot_cell_trajectory(cds, color_by = "downgene") +
  scale_color_gradient(low = "blue", high = "red") +   
  geom_point(aes(size = Pseudotime, color = downgene), alpha = 0.6) +   
  theme_minimal() +
  theme(legend.position = "right")   

dev.off()

#COAD####
#COAD
library(rlang)
#packageVersion("dplyr")
#install.packages("https://www.bioconductor.org/packages/3.8/bioc/src/contrib/BiocGenerics_0.28.0.tar.gz",repos=NULL,type="source")
library(monocle)
library(monocle)
library(Seurat)
library(SeuratObject)
#load("death_result_2\\GSE135337\\HP\\GSE146771_geneplotneed.Rdata")
#load("result\\GSE146771\\HP\\GSE146771_geneplotneed.Rdata")coad
#load("H:\\database\\SCDTD\\data\\chenlong\\Esophageal Cancer\\disease\\GSE196756_geneplotneed.Rdata")
#load("H:\\database\\SCDTD\\data\\chenlong\\Colon Cancer\\disease\\GSE140288_geneplotneed.Rdata")
load("H:\\database\\SCDTD\\data\\chenlong\\Colon Cancer\\disease\\GSE222300_geneplotneed.Rdata")
pmbc<-data.combined.new


Idents(pmbc) <- "orig.ident"
ob1 <- subset(x = pmbc, downsample = 300)

expr_matrix<-as(as.matrix(pmbc@assays[["RNA"]]@layers[["counts"]]),'sparseMatrix')
expr_matrix<-as(as.matrix(pmbc@assays[["SCT"]]@counts),'sparseMatrix')
p_data<-pmbc@meta.data
p_data$celltype<-pmbc@active.ident
f_data<-data.frame(gene_short_name=row.names(pmbc),row.names = row.names(pmbc))
dim(expr_matrix)

pd<-new('AnnotatedDataFrame',data=p_data)
fd<-new('AnnotatedDataFrame',data=f_data)





cds<-newCellDataSet(expr_matrix,
                    phenoData=pd,
                    featureData=fd,
                    lowerDetectionLimit=0.1,
                    expressionFamily=negbinomial.size())



cds<-estimateSizeFactors(cds)
#1gc()
# memory.limit(123456)
# rm(expr_matrix)
# rm(f_data)
# rm(p_data)
cds<-estimateDispersions(cds)

cds<-detectGenes(cds,min_expr = 0.1)

print(head(fData(cds)))


gene<-read.csv("G:\\pan_result\\gpu_result_tcga\\COAD\\424featuresweightme.csv",row.names = 1)
expressed_genes<-colnames(gene)

# expressed_genes = read.table("validategene.csv",header = T,sep = ",")
# expressed_genes = expressed_genes$X0


setwd("G:\\pan_result\\gpu_result_tcga\\COAD\\scrna")
cds <- setOrderingFilter(cds, expressed_genes)
pdf("train.ordergenes_allcell.pdf")
plot_ordering_genes(cds)
dev.off()



cds <- reduceDimension(cds, max_components = 5,
                       method = 'DDRTree')##method
cds <- orderCells(cds)
setwd("G:\\pan_result\\gpu_result_tcga\\COAD\\scrna")
pdf("train.monocle.state5_allcell.pdf",width = 7,height = 7)
plot_cell_trajectory(cds, color_by = "Pseudotime",size=1,show_backbone=TRUE)
dev.off()

cds_subset<-cds['SAMD11',]
plot_genes_in_pseudotime(cds_subset, color_by = "cell_type")

save(cds,file = "cds.RData")

pdf("train.monocle.celltype5_allcell.pdf",width = 7,height = 7)
plot_cell_trajectory(cds,color_by = "cell_type",size=1,show_backbone = TRUE)+
  geom_point(aes(size = Pseudotime , color = cell_type)) + 
  scale_color_manual(breaks = c("X", "Y", "Z"), values=c("darkblue", "darkred", "orange","grey","yellow","green","skyblue","pink","purple")) 
plot_cell_trajectory(cds,color_by = "cell_type",size=1,show_backbone = TRUE) + scale_color_manual( values=c("darkblue", "darkred", "orange","grey","yellow","green","skyblue","pink","purple")) + theme(legend.position = "right")
dev.off()

library(ggpubr)
library(tidyr)
packageVersion("tidyr")

#install.packages("tidyr")

#load("cds.RData")
df <- pData(cds) 

pdf("train.monocle.pseudotime_allcell.pdf",width = 7,height = 7)
ggplot(df, aes(Pseudotime, colour = cell_type, fill=cell_type)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2( )
dev.off()
write.csv(df,file = "pData.csv")



# 
library(SingleCellExperiment)
library(Seurat)




# devtools::install_github("alanocallaghan/scater")
library(scater)
data.combined.new@meta.data[["Pseudotime"]]<-cds@phenoData@data[["Pseudotime"]]
sce_gene <- as.SingleCellExperiment(data.combined.new)


label<-read.csv("G:\\pan_result\\gpu_result_tcga\\COAD\\label.csv")
gene1<-gene[label$X[which(label$label=="1")],]
gene2<-gene[label$X[which(label$label=="2")],]
data<-cbind(t(gene1),t(gene2))
T_p<-c()
FC_result<-c()
for(i in 1:nrow(data)){
  T_result<-wilcox.test(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]),as.numeric(data[i,1:nrow(gene1)]),paired = F,var.equal = T)
  
  p1<-T_result$p.value

  T_p<-c(T_p,p1)
  #T_value<-c(T_value,t1)
  FC<-mean(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]))/mean(as.numeric(data[i,1:nrow(gene1)]))

  FC_result<-c(FC_result,FC)
}

FDR_p<-p.adjust(T_p,method="fdr",length(T_p))#
DEG_all<-data.frame(rownames(data),T_p,FC_result,stringsAsFactors = FALSE)
colnames(DEG_all)<-c("Gene","pval","FC") 


DEG<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC>2),]#1295
plotExpression(sce_gene, features="RGS16", color_by = "Pseudotime")#p=0.002117181

pdf("monocle.pseudotime_allcell_gene.pdf",width = 7,height = 7)
plotExpression(sce_gene, features="RGS16", color_by = "Pseudotime")#p=0.002117181

dev.off()
write.csv(df,file = "pData.csv")




p1 <- plot_cell_trajectory(cds, color_by = "cell_type", size=1, show_backbone = TRUE) +
  scale_color_manual(values = c("darkblue", "darkred", "orange", "grey", "yellow", 
                                "green", "skyblue", "pink", "purple")) +
  theme(legend.position = "right")


p2 <- plotExpression(sce_gene, features="RGS16", color_by = "Pseudotime") + 
  theme_minimal()
library(gridExtra)
library(ggplot2)
library(monocle)

grid.arrange(p1, p2, ncol = 2)




library(BayesPrism)
# sc.dat <- pmbc@assays[["RNA"]]
# 
# dim(sc.dat)#33538 17135
# sc.dat<-t(sc.dat)
# cell.type.labels <- pmbc@meta.data[["cell_type"]]
# sc.dat<-as.matrix(sc.dat)

exing<-read.table("G:\\pan_result\\gpu_result_tcga\\COAD\\scrna\\COAD_EPI_cells_E.txt")
cell_types <- rep("Malignant epithelial cells",nrow(exing))
cell_names<-exing$x

named_cell_types <- setNames(cell_types, cell_names)


pmbc<-data.combined.new

pmbc@meta.data$cell_type <- as.character(pmbc@meta.data$cell_type)


match_index <- match(rownames(pmbc@meta.data), names(named_cell_types))

pmbc@meta.data$cell_type[!is.na(match_index)] <- named_cell_types[match_index[!is.na(match_index)]]


pmbc@meta.data$cell_type <- factor(pmbc@meta.data$cell_type)


sc.dat <- pmbc[["SCT"]]@counts
dim(sc.dat)#20932 12480
sc.dat<-t(sc.dat)
cell.type.labels <- pmbc@meta.data[["cell_type"]]
sc.dat<-as.matrix(sc.dat)

nrow(sc.dat)*0.3#3744
sc.dat.filtered <- cleanup.genes (input=sc.dat,
                                  input.type="count.matrix",
                                  species="hs",
                                  gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY","hb","act"),
                                  exp.cells=nrow(sc.dat)*0.3)
dim(sc.dat.filtered)#12480  1041

save(sc.dat.filtered,file="sc.dat.filtered_70.Rdata")


rownames(sc.dat.filtered)[1:4]
sc.dat.filtered<-as.data.frame(sc.dat.filtered)
sc.dat.filtered$cellname<-as.character(cell.type.labels)
qiuhe<-aggregate(sc.dat.filtered[,-(ncol(sc.dat.filtered))], by = sc.dat.filtered[ncol(sc.dat.filtered)], FUN = sum)
qiujunzhi<-aggregate(sc.dat.filtered[,-(ncol(sc.dat.filtered))], by = sc.dat.filtered[ncol(sc.dat.filtered)], FUN = mean)
qiuhe<-t(qiuhe)
qiujunzhi<-t(qiujunzhi)
write.table(qiuhe,"LM4_sum.txt",col.names = F,sep = "\t")
write.table(qiujunzhi,"LM4_mean.txt",col.names = F,sep = "\t")
# a<-t(sc.dat.filtered)
# max(a[3,])
# #b<-as.matrix(pmbc@active.ident)
# colnames(a)<-cell.type.labels
# class(a)
# a<-as.data.frame(a)
# b<-sapply(split.default(a, names(a)), rowMeans)
# b<-sapply(split.default(a, names(a)), rowMeans)
# dim(b)#21178     4
# head(b)
# max(b[,1])
 


#ESCA####
# 
library(rlang)
#packageVersion("dplyr")
#install.packages("https://www.bioconductor.org/packages/3.8/bioc/src/contrib/BiocGenerics_0.28.0.tar.gz",repos=NULL,type="source")
library(monocle)
library(monocle)
library(Seurat)
library(SeuratObject)
#load("death_result_2\\GSE135337\\HP\\GSE146771_geneplotneed.Rdata")
#load("result\\GSE146771\\HP\\GSE146771_geneplotneed.Rdata")
load("H:\\database\\SCDTD\\data\\chenlong\\Esophageal Cancer\\disease\\GSE196756_geneplotneed.Rdata")
pmbc<-data.combined.new
GetAssay(data.combined.new,assay = "RNA")
unique(data.combined.new@meta.data[["orig.ident"]])

Idents(pmbc) <- "orig.ident"
ob1 <- subset(x = pmbc, downsample = 300)

expr_matrix<-as(as.matrix(ob1@assays[["RNA"]]@layers[["counts"]]),'sparseMatrix')
expr_matrix<-as(as.matrix(ob1@assays[["SCT"]]@counts),'sparseMatrix')
p_data<-ob1@meta.data
p_data$celltype<-ob1@active.ident
f_data<-data.frame(gene_short_name=row.names(ob1),row.names = row.names(ob1))
dim(expr_matrix)

pd<-new('AnnotatedDataFrame',data=p_data)
fd<-new('AnnotatedDataFrame',data=f_data)


library(igraph)



cds<-newCellDataSet(expr_matrix,
                    phenoData=pd,
                    featureData=fd,
                    lowerDetectionLimit=0.1,
                    expressionFamily=negbinomial.size())



cds<-estimateSizeFactors(cds)
#1gc()
# memory.limit(123456)
# rm(expr_matrix)
# rm(f_data)
# rm(p_data)
cds<-estimateDispersions(cds)

cds<-detectGenes(cds,min_expr = 0.1)

print(head(fData(cds)))




gene<-read.csv("G:\\pan_result\\gpu_result_tcga\\ESCA\\202featuresweightme.csv",row.names = 1)
expressed_genes<-colnames(gene)

# expressed_genes = read.table("validategene.csv",header = T,sep = ",")
# expressed_genes = expressed_genes$X0


setwd("G:\\pan_result\\gpu_result_tcga\\ESCA\\scrna")
cds <- setOrderingFilter(cds, expressed_genes)
pdf("train.ordergenes.pdf")
plot_ordering_genes(cds)
dev.off()



cds <- reduceDimension(cds, max_components = 5,
                       method = 'DDRTree')##method
cds <- orderCells(cds)
setwd("G:\\pan_result\\gpu_result_tcga\\ESCA\\scrna")
pdf("train.monocle.state5.pdf",width = 7,height = 7)
plot_cell_trajectory(cds, color_by = "Pseudotime",size=1,show_backbone=TRUE)
dev.off()

save(cds,file = "cds.RData")

pdf("train.monocle.celltype5.pdf",width = 7,height = 7)
plot_cell_trajectory(cds,color_by = "cell_type",size=1,show_backbone = TRUE)+
  geom_point(aes(size = Pseudotime , color = cell_type)) + 
  scale_color_manual(breaks = c("X", "Y", "Z"), values=c("darkblue", "darkred", "orange","grey","yellow","green","skyblue","pink","purple")) 
plot_cell_trajectory(cds,color_by = "cell_type",size=1,show_backbone = TRUE) + scale_color_manual( values=c("darkblue", "darkred", "orange","grey","yellow","green","skyblue","pink","purple")) + theme(legend.position = "right")
dev.off()
?plot_cell_trajectory

library(ggpubr)
library(tidyr)
packageVersion("tidyr")

#install.packages("tidyr")

#load("cds.RData")
df <- pData(cds) 

pdf("train.monocle.pseudotime.pdf",width = 7,height = 7)
ggplot(df, aes(Pseudotime, colour = cell_type, fill=cell_type)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2( )
dev.off()
write.csv(df,file = "pData.csv")


#####
# 
library(SingleCellExperiment)
library(Seurat)




# devtools::install_github("alanocallaghan/scater")
library(scater)
# data.combined.new@meta.data[["Pseudotime"]]<-cds@phenoData@data[["Pseudotime"]]
# sce_gene <- as.SingleCellExperiment(data.combined.new)


label<-read.csv("G:\\pan_result\\gpu_result_tcga\\ESCA\\label.csv")
gene1<-gene[label$X[which(label$label=="1")],]
gene2<-gene[label$X[which(label$label=="2")],]
data<-cbind(t(gene1),t(gene2))
T_p<-c()
FC_result<-c()
for(i in 1:nrow(data)){
  T_result<-wilcox.test(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]),as.numeric(data[i,1:nrow(gene1)]),paired = F,var.equal = T)
  
  p1<-T_result$p.value

  T_p<-c(T_p,p1)
  #T_value<-c(T_value,t1)
  FC<-mean(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]))/mean(as.numeric(data[i,1:nrow(gene1)]))

  FC_result<-c(FC_result,FC)
}

FDR_p<-p.adjust(T_p,method="fdr",length(T_p))#
DEG_all<-data.frame(rownames(data),T_p,FC_result,stringsAsFactors = FALSE)
colnames(DEG_all)<-c("Gene","pval","FC") 


DEG<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC>2),]#1295

getwd()
# pdf("monocle.pseudotime_allcell_gene2.pdf",width = 7,height = 7)
# 
# plot_cell_trajectory(cds, color_by = "NPTX2") +
#   scale_color_gradient(low = "blue", high = "red") +   
#   geom_point(aes(size = Pseudotime, color = NPTX2), alpha = 0.6) +   
#   theme_minimal() +
#   theme(legend.position = "right")   
# 
# dev.off()



expression_matrix <- exprs(cds)


gene_expression <- expression_matrix["OSM", ]#2.340444e-09

# gene_expression <- FetchData(cds, vars = "NPTX2", slot = "data")

cds@phenoData@data[["OSM"]]<-gene_expression

getwd()
pdf("monocle.pseudotime_gene2.pdf",width = 7,height = 7)

plot_cell_trajectory(cds, color_by = "OSM") +
  scale_color_gradient(low = "blue", high = "red") +   
  geom_point(aes(size = Pseudotime, color = OSM), alpha = 0.6) +   
  theme_minimal() +
  theme(legend.position = "right")   

dev.off()



write.csv(DEG_all,file = "DEG_all.csv")
DEG<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC>2),]#1295
DEG_down<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC<0.5),]#1295


expression_matrix <- exprs(cds)


# gene_expression <- expression_matrix["TYMP", ]
# 
# # gene_expression <- FetchData(cds, vars = "ZNF487", slot = "data")
# 
# cds@phenoData@data[["TYMP"]]<-gene_expression
gene_expression2<-0

DEG<-intersect(DEG$Gene,expression_matrix@Dimnames[[1]])
for (i in 1:length(DEG)) {
  gene_expression <- expression_matrix[DEG[i], ]
  gene_expression2<-gene_expression2+gene_expression
}


cds@phenoData@data[["upgene"]]<-gene_expression2




getwd()
pdf("monocle.pseudotime_upgene.pdf",width = 7,height = 7)

plot_cell_trajectory(cds, color_by = "upgene") +
  scale_color_gradient(low = "blue", high = "red") +   
  geom_point(aes(size = Pseudotime, color = upgene), alpha = 0.6) +   
  theme_minimal() +
  theme(legend.position = "right")   

dev.off()

#downgene####
gene_expression2<-0

DEG_down<-intersect(DEG_down$Gene,expression_matrix@Dimnames[[1]])
for (i in 1:length(DEG_down)) {
  gene_expression <- expression_matrix[DEG_down[i], ]
  gene_expression2<-gene_expression2+gene_expression
}

cds@phenoData@data[["downgene"]]<-gene_expression2


#coad

getwd()
pdf("monocle.pseudotime_downgene.pdf",width = 7,height = 7)

plot_cell_trajectory(cds, color_by = "downgene") +
  scale_color_gradient(low = "blue", high = "red") +   
  geom_point(aes(size = Pseudotime, color = downgene), alpha = 0.6) +   
  theme_minimal() +
  theme(legend.position = "right")   

dev.off()



library(BayesPrism)
sc.dat <- pmbc[["RNA"]]@counts

dim(sc.dat)#33538 17135
sc.dat<-t(sc.dat)
cell.type.labels <- pmbc@meta.data[["cell_type"]]
sc.dat<-as.matrix(sc.dat)


sc.dat <- pmbc[["SCT"]]@counts
dim(sc.dat)#20932 12480
sc.dat<-t(sc.dat)
cell.type.labels <- pmbc@meta.data[["cell_type"]]
sc.dat<-as.matrix(sc.dat)

nrow(sc.dat)*0.3#3744
sc.dat.filtered <- cleanup.genes (input=sc.dat,
                                  input.type="count.matrix",
                                  species="hs",
                                  gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY","hb","act"),
                                  exp.cells=nrow(sc.dat)*0.3)
dim(sc.dat.filtered)#12480  1041
save(sc.dat.filtered,file="sc.dat.filtered_70.Rdata")

rownames(sc.dat.filtered)[1:4]
sc.dat.filtered<-as.data.frame(sc.dat.filtered)
sc.dat.filtered$cellname<-as.character(cell.type.labels)
qiuhe<-aggregate(sc.dat.filtered[,-(ncol(sc.dat.filtered))], by = sc.dat.filtered[ncol(sc.dat.filtered)], FUN = sum)
qiujunzhi<-aggregate(sc.dat.filtered[,-(ncol(sc.dat.filtered))], by = sc.dat.filtered[ncol(sc.dat.filtered)], FUN = mean)
qiuhe<-t(qiuhe)
qiujunzhi<-t(qiujunzhi)
write.table(qiuhe,"LM4_sum.txt",col.names = F,sep = "\t")
write.table(qiujunzhi,"LM4_mean.txt",col.names = F,sep = "\t")
# a<-t(sc.dat.filtered)
# max(a[3,])
# #b<-as.matrix(pmbc@active.ident)
# colnames(a)<-cell.type.labels
# class(a)
# a<-as.data.frame(a)
# b<-sapply(split.default(a, names(a)), rowMeans)
# b<-sapply(split.default(a, names(a)), rowMeans)
# dim(b)#21178     4
# head(b)
# max(b[,1])
 


#HNSC####
# HNSC
library(rlang)
#packageVersion("dplyr")
#install.packages("https://www.bioconductor.org/packages/3.8/bioc/src/contrib/BiocGenerics_0.28.0.tar.gz",repos=NULL,type="source")
library(monocle)
library(monocle)
library(Seurat)
library(SeuratObject)
#load("death_result_2\\GSE135337\\HP\\GSE146771_geneplotneed.Rdata")
#load("result\\GSE146771\\HP\\GSE146771_geneplotneed.Rdata")
#load("H:\\database\\SCDTD\\data\\chenlong\\Esophageal Cancer\\disease\\GSE196756_geneplotneed.Rdata")
#load("H:\\database\\SCDTD\\data\\chenlong\\Colon Cancer\\disease\\GSE140288_geneplotneed.Rdata")
#load("H:\\database\\SCDTD\\data\\chenlong\\Cervical Cancer\\disease\\GSE168652_geneplotneed.Rdata")
#load("H:\\database\\SCDTD\\data\\chenlong\\Liver Cancer\\Liver Cancer\\GSE125449_geneplotneed.Rdata")
# load("H:\\database\\SCDTD\\data\\result\\GSE263365_HP\\GSE263365_geneplotneed.Rdata")
load("H:\\database\\HNSC\\geo\\GSE243359_seruat.Rdata")
GetAssay(data.combined.new,assay = "RNA")
unique(data.combined.new@meta.data[["orig.ident"]])
data.combined.new@assays[["RNA"]]@layers[["counts"]]@Dimnames[[1]]<-rownames(data.combined.new)
data.combined.new@assays[["RNA"]]@layers[["counts"]]@Dimnames[[2]]<-colnames(data.combined.new)

pmbc<-data.combined.new


Idents(pmbc) <- "orig.ident"
ob1 <- subset(x = pmbc, downsample = 300)

expr_matrix<-as(as.matrix(ob1@assays[["RNA"]]@layers[["counts"]]),'sparseMatrix')
expr_matrix<-as(as.matrix(ob1@assays[["SCT"]]@counts),'sparseMatrix')
p_data<-ob1@meta.data
p_data$celltype<-ob1@active.ident
f_data<-data.frame(gene_short_name=row.names(ob1),row.names = row.names(ob1))
dim(expr_matrix)

pd<-new('AnnotatedDataFrame',data=p_data)
fd<-new('AnnotatedDataFrame',data=f_data)





cds<-newCellDataSet(expr_matrix,
                    phenoData=pd,
                    featureData=fd,
                    lowerDetectionLimit=0.1,
                    expressionFamily=negbinomial.size())



cds<-estimateSizeFactors(cds)
#1gc()
# memory.limit(123456)
# rm(expr_matrix)
# rm(f_data)
# rm(p_data)
cds<-estimateDispersions(cds)

cds<-detectGenes(cds,min_expr = 0.1)

print(head(fData(cds)))


gene<-read.csv("G:\\pan_result\\gpu_result_tcga\\HNSC\\297featuresweightme.csv",row.names = 1)
expressed_genes<-colnames(gene)

# expressed_genes = read.table("validategene.csv",header = T,sep = ",")
# expressed_genes = expressed_genes$X0


setwd("G:\\pan_result\\gpu_result_tcga\\HNSC\\scrna")
cds <- setOrderingFilter(cds, expressed_genes)
pdf("train.ordergenes.pdf")
plot_ordering_genes(cds)
dev.off()



cds <- reduceDimension(cds, max_components = 5,
                       method = 'DDRTree')##method
cds <- orderCells(cds)
dir.create("G:\\pan_result\\gpu_result_tcga\\HNSC\\scrna")
setwd("G:\\pan_result\\gpu_result_tcga\\HNSC\\scrna")
pdf("train.monocle.state5.pdf",width = 7,height = 7)
plot_cell_trajectory(cds, color_by = "Pseudotime",size=1,show_backbone=TRUE)
dev.off()

save(cds,file = "cds.RData")

pdf("train.monocle.celltype5.pdf",width = 7,height = 7)
plot_cell_trajectory(cds,color_by = "cell_type",size=1,show_backbone = TRUE)+
  geom_point(aes(size = Pseudotime , color = cell_type)) + 
  scale_color_manual(breaks = c("X", "Y", "Z"), values=c("darkblue", "darkred", "orange","grey","yellow","green","skyblue","pink","purple")) 
plot_cell_trajectory(cds,color_by = "cell_type",size=1,show_backbone = TRUE) + scale_color_manual( values=c("darkblue", "darkred", "orange","grey","yellow","green","skyblue","pink","purple")) + theme(legend.position = "right")
dev.off()

library(ggpubr)
library(tidyr)
packageVersion("tidyr")

#install.packages("tidyr")

#load("cds.RData")
df <- pData(cds) 

pdf("train.monocle.pseudotime.pdf",width = 7,height = 7)
ggplot(df, aes(Pseudotime, colour = cell_type, fill=cell_type)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2( )
dev.off()
write.csv(df,file = "pData.csv")

#####
# 
library(SingleCellExperiment)
library(Seurat)




# devtools::install_github("alanocallaghan/scater")
library(scater)
# data.combined.new@meta.data[["Pseudotime"]]<-cds@phenoData@data[["Pseudotime"]]
# sce_gene <- as.SingleCellExperiment(data.combined.new)


label<-read.csv("G:\\pan_result\\gpu_result_tcga\\HNSC\\label.csv")
gene1<-gene[label$X[which(label$label=="1")],]
gene2<-gene[label$X[which(label$label=="2")],]
data<-cbind(t(gene1),t(gene2))
T_p<-c()
FC_result<-c()
for(i in 1:nrow(data)){
  T_result<-wilcox.test(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]),as.numeric(data[i,1:nrow(gene1)]),paired = F,var.equal = T)
  
  p1<-T_result$p.value

  T_p<-c(T_p,p1)
  #T_value<-c(T_value,t1)
  FC<-mean(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]))/mean(as.numeric(data[i,1:nrow(gene1)]))

  FC_result<-c(FC_result,FC)
}

FDR_p<-p.adjust(T_p,method="fdr",length(T_p))
DEG_all<-data.frame(rownames(data),T_p,FC_result,stringsAsFactors = FALSE)
colnames(DEG_all)<-c("Gene","pval","FC") 


DEG<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC>2),]#1295

getwd()
# pdf("monocle.pseudotime_allcell_gene2.pdf",width = 7,height = 7)
# 
# plot_cell_trajectory(cds, color_by = "NPTX2") +
#   scale_color_gradient(low = "blue", high = "red") +   
#   geom_point(aes(size = Pseudotime, color = NPTX2), alpha = 0.6) +   
#   theme_minimal() +
#   theme(legend.position = "right")   
# 
# dev.off()



expression_matrix <- exprs(cds)


gene_expression <- expression_matrix["PPP1R16A", ]#2.340444e-09

# gene_expression <- FetchData(cds, vars = "NPTX2", slot = "data")

cds@phenoData@data[["PPP1R16A"]]<-gene_expression

getwd()
pdf("monocle.pseudotime_gene2.pdf",width = 7,height = 7)

plot_cell_trajectory(cds, color_by = "PPP1R16A") +
  scale_color_gradient(low = "blue", high = "red") +   
  geom_point(aes(size = Pseudotime, color = PPP1R16A), alpha = 0.6) +   
  theme_minimal() +
  theme(legend.position = "right")   

dev.off()


#cluster####
label<-read.table("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_cluster.txt")
gene1<-gene[rownames(label)[which(label$clusters=="1")],]
gene2<-gene[rownames(label)[which(label$clusters=="2")],]
data<-cbind(t(gene1),t(gene2))
T_p<-c()
FC_result<-c()
for(i in 1:nrow(data)){
  T_result<-wilcox.test(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]),as.numeric(data[i,1:nrow(gene1)]),paired = F,var.equal = T)
  p1<-T_result$p.value

  T_p<-c(T_p,p1)
  #T_value<-c(T_value,t1)
  FC<-mean(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]))/mean(as.numeric(data[i,1:nrow(gene1)]))

  FC_result<-c(FC_result,FC)
}


DEG_all<-data.frame(rownames(data),T_p,FC_result,stringsAsFactors = FALSE)

colnames(DEG_all)<-c("Gene","pval","FC") 

write.csv(DEG_all,file = "DEG_all.csv")
DEG<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC>2),]#1295
DEG_down<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC<0.5),]#1295


expression_matrix <- exprs(cds)


gene_expression <- expression_matrix["TYMP", ]

# gene_expression <- FetchData(cds, vars = "ZNF487", slot = "data")

cds@phenoData@data[["TYMP"]]<-gene_expression
gene_expression2<-0

DEG<-intersect(DEG$Gene,expression_matrix@Dimnames[[1]])
for (i in 1:length(DEG)) {
  gene_expression <- expression_matrix[DEG[i], ]
  gene_expression2<-gene_expression2+gene_expression
}


cds@phenoData@data[["upgene"]]<-gene_expression2




getwd()
pdf("monocle.pseudotime_upgene.pdf",width = 7,height = 7)

plot_cell_trajectory(cds, color_by = "upgene") +
  scale_color_gradient(low = "blue", high = "red") +   
  geom_point(aes(size = Pseudotime, color = upgene), alpha = 0.6) +   
  theme_minimal() +
  theme(legend.position = "right")   

dev.off()

#downgene####
gene_expression2<-0

DEG_down<-intersect(DEG_down$Gene,expression_matrix@Dimnames[[1]])
for (i in 1:length(DEG_down)) {
  gene_expression <- expression_matrix[DEG_down[i], ]
  gene_expression2<-gene_expression2+gene_expression
}

cds@phenoData@data[["downgene"]]<-gene_expression2




getwd()
pdf("monocle.pseudotime_downgene.pdf",width = 7,height = 7)

plot_cell_trajectory(cds, color_by = "downgene") +
  scale_color_gradient(low = "blue", high = "red") +   
  geom_point(aes(size = Pseudotime, color = downgene), alpha = 0.6) +   
  theme_minimal() +
  theme(legend.position = "right")   

dev.off()




library(BayesPrism)
# sc.dat <- pmbc@assays[["RNA"]]
# 
# dim(sc.dat)#33538 17135
# sc.dat<-t(sc.dat)
# cell.type.labels <- pmbc@meta.data[["cell_type"]]
# sc.dat<-as.matrix(sc.dat)

exing<-read.table("G:\\pan_result\\gpu_result_tcga\\HNSC\\scrna\\HNSC_EPI_cells_E.txt")
cell_types <- rep("Malignant epithelial cells",nrow(exing))
cell_names<-exing$x

named_cell_types <- setNames(cell_types, cell_names)


pmbc<-data.combined.new

pmbc@meta.data$cell_type <- as.character(pmbc@meta.data$cell_type)


match_index <- match(rownames(pmbc@meta.data), names(named_cell_types))

pmbc@meta.data$cell_type[!is.na(match_index)] <- named_cell_types[match_index[!is.na(match_index)]]


pmbc@meta.data$cell_type <- factor(pmbc@meta.data$cell_type)


sc.dat <- pmbc[["SCT"]]@counts
sc.dat <- pmbc@assays[["RNA"]]@layers[["counts"]]
dim(sc.dat)#22393  1849
sc.dat<-t(sc.dat)
cell.type.labels <- pmbc@meta.data[["cell_type"]]
sc.dat<-as.matrix(sc.dat)

nrow(sc.dat)*0.3#3744
sc.dat.filtered <- cleanup.genes (input=sc.dat,
                                  input.type="count.matrix",
                                  species="hs",
                                  gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY","hb","act"),
                                  exp.cells=nrow(sc.dat)*0.3)
dim(sc.dat.filtered)#12480  1041

save(sc.dat.filtered,file="sc.dat.filtered_70.Rdata")


rownames(sc.dat.filtered)[1:4]
sc.dat.filtered<-as.data.frame(sc.dat.filtered) 
sc.dat.filtered$cellname<-as.character(cell.type.labels)
qiuhe<-aggregate(sc.dat.filtered[,-(ncol(sc.dat.filtered))], by = sc.dat.filtered[ncol(sc.dat.filtered)], FUN = sum)
qiujunzhi<-aggregate(sc.dat.filtered[,-(ncol(sc.dat.filtered))], by = sc.dat.filtered[ncol(sc.dat.filtered)], FUN = mean)
qiuhe<-t(qiuhe)
qiujunzhi<-t(qiujunzhi)
write.table(qiuhe,"LM4_sum.txt",col.names = F,sep = "\t")
write.table(qiujunzhi,"LM4_mean.txt",col.names = F,sep = "\t")
# a<-t(sc.dat.filtered)
# max(a[3,])
# #b<-as.matrix(pmbc@active.ident)
# colnames(a)<-cell.type.labels
# class(a)
# a<-as.data.frame(a)
# b<-sapply(split.default(a, names(a)), rowMeans)
# b<-sapply(split.default(a, names(a)), rowMeans)
# dim(b)#21178     4
# head(b)
# max(b[,1])
 


#KIRC####
# KIRC
library(rlang)
#packageVersion("dplyr")
#install.packages("https://www.bioconductor.org/packages/3.8/bioc/src/contrib/BiocGenerics_0.28.0.tar.gz",repos=NULL,type="source")
library(monocle)
library(monocle)
library(Seurat)
library(SeuratObject)
#load("death_result_2\\GSE135337\\HP\\GSE146771_geneplotneed.Rdata")
#load("result\\GSE146771\\HP\\GSE146771_geneplotneed.Rdata")KIRC
#load("H:\\database\\SCDTD\\data\\chenlong\\Esophageal Cancer\\disease\\GSE196756_geneplotneed.Rdata")
#load("H:\\database\\SCDTD\\data\\chenlong\\Colon Cancer\\disease\\GSE140288_geneplotneed.Rdata")
load("H:\\database\\SCDTD\\data\\chenlong\\Kidney Chromophobe\\disease\\HP\\GSE159115_geneplotneed.Rdata")
GetAssay(data.combined.new,assay = "RNA")
unique(data.combined.new@meta.data[["orig.ident"]])
pmbc<-data.combined.new


Idents(pmbc) <- "orig.ident"
ob1 <- subset(x = pmbc, downsample = 300)

expr_matrix<-as(as.matrix(ob1@assays[["RNA"]]@layers[["counts"]]),'sparseMatrix')
expr_matrix<-as(as.matrix(ob1@assays[["SCT"]]@counts),'sparseMatrix')
p_data<-ob1@meta.data
p_data$celltype<-ob1@active.ident
f_data<-data.frame(gene_short_name=row.names(ob1),row.names = row.names(ob1))
dim(expr_matrix)

pd<-new('AnnotatedDataFrame',data=p_data)
fd<-new('AnnotatedDataFrame',data=f_data)





cds<-newCellDataSet(expr_matrix,
                    phenoData=pd,
                    featureData=fd,
                    lowerDetectionLimit=0.1,
                    expressionFamily=negbinomial.size())



cds<-estimateSizeFactors(cds)
#1gc()
memory.limit(123456)
rm(expr_matrix)
rm(f_data)
rm(p_data)
cds<-estimateDispersions(cds)

cds<-detectGenes(cds,min_expr = 0.1)

print(head(fData(cds)))


gene<-read.csv("G:\\pan_result\\gpu_result_tcga\\KIRC\\563featuresweightme.csv",row.names = 1)
expressed_genes<-colnames(gene)

# expressed_genes = read.table("validategene.csv",header = T,sep = ",")
# expressed_genes = expressed_genes$X0


setwd("G:\\pan_result\\gpu_result_tcga\\KIRC\\scrna")
cds <- setOrderingFilter(cds, expressed_genes)
pdf("train.ordergenes.pdf")
plot_ordering_genes(cds)
dev.off()



cds <- reduceDimension(cds, max_components = 5,
                       method = 'DDRTree')##method
cds <- orderCells(cds)
setwd("G:\\pan_result\\gpu_result_tcga\\KIRC\\scrna")
pdf("train.monocle.state5.pdf",width = 7,height = 7)
plot_cell_trajectory(cds, color_by = "Pseudotime",size=1,show_backbone=TRUE)
dev.off()

save(cds,file = "cdsl.RData")

pdf("train.monocle.celltype5.pdf",width = 7,height = 7)
plot_cell_trajectory(cds,color_by = "cell_type",size=1,show_backbone = TRUE)+
  geom_point(aes(size = Pseudotime , color = cell_type)) + 
  scale_color_manual(breaks = c("X", "Y", "Z"), values=c("darkblue", "darkred", "orange","grey","yellow","green","skyblue","pink","purple","#458B74","darkgoldenrod","#CD3333")) 
plot_cell_trajectory(cds,color_by = "cell_type",size=1,show_backbone = TRUE) + scale_color_manual( values=c("darkblue", "darkred", "orange","grey","yellow","green","skyblue","pink","purple","#458B74","darkgoldenrod","#CD3333")) + theme(legend.position = "right")
dev.off()

library(ggpubr)
library(tidyr)
packageVersion("tidyr")

#install.packages("tidyr")

#load("cds.RData")
df <- pData(cds) 

pdf("train.monocle.pseudotime.pdf",width = 7,height = 7)
ggplot(df, aes(Pseudotime, colour = cell_type, fill=cell_type)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2( )
dev.off()
write.csv(df,file = "pData.csv")


#####
# 
library(SingleCellExperiment)
library(Seurat)




# devtools::install_github("alanocallaghan/scater")
library(scater)
# data.combined.new@meta.data[["Pseudotime"]]<-cds@phenoData@data[["Pseudotime"]]
# sce_gene <- as.SingleCellExperiment(data.combined.new)


label<-read.csv("G:\\pan_result\\gpu_result_tcga\\KIRC\\label.csv")
gene1<-gene[label$X[which(label$label=="1")],]
gene2<-gene[label$X[which(label$label=="2")],]
data<-cbind(t(gene1),t(gene2))
T_p<-c()
FC_result<-c()
for(i in 1:nrow(data)){
  T_result<-wilcox.test(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]),as.numeric(data[i,1:nrow(gene1)]),paired = F,var.equal = T)
  
  p1<-T_result$p.value

  T_p<-c(T_p,p1)
  #T_value<-c(T_value,t1)
  FC<-mean(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]))/mean(as.numeric(data[i,1:nrow(gene1)]))

  FC_result<-c(FC_result,FC)
}

FDR_p<-p.adjust(T_p,method="fdr",length(T_p))#
DEG_all<-data.frame(rownames(data),T_p,FC_result,stringsAsFactors = FALSE)
colnames(DEG_all)<-c("Gene","pval","FC") 


DEG<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC>2),]#1295

getwd()
# pdf("monocle.pseudotime_allcell_gene2.pdf",width = 7,height = 7)
# 
# plot_cell_trajectory(cds, color_by = "NPTX2") +
#   scale_color_gradient(low = "blue", high = "red") +   
#   geom_point(aes(size = Pseudotime, color = NPTX2), alpha = 0.6) +   
#   theme_minimal() +
#   theme(legend.position = "right")   
# 
# dev.off()



expression_matrix <- exprs(cds)


gene_expression <- expression_matrix["NPTX2", ]#2.340444e-09

# gene_expression <- FetchData(cds, vars = "NPTX2", slot = "data")

cds@phenoData@data[["NPTX2"]]<-gene_expression

getwd()
pdf("monocle.pseudotime_gene2.pdf",width = 7,height = 7)

plot_cell_trajectory(cds, color_by = "NPTX2") +
  scale_color_gradient(low = "blue", high = "red") +   
  geom_point(aes(size = Pseudotime, color = NPTX2), alpha = 0.6) +   
  theme_minimal() +
  theme(legend.position = "right")   

dev.off()
#cluster####
label<-read.table("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\KIRC_cluster.txt")
gene1<-gene[rownames(label)[which(label$clusters=="1")],]
gene2<-gene[rownames(label)[which(label$clusters=="2")],]
data<-cbind(t(gene1),t(gene2))
T_p<-c()
FC_result<-c()
for(i in 1:nrow(data)){
  T_result<-wilcox.test(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]),as.numeric(data[i,1:nrow(gene1)]),paired = F,var.equal = T)
  
  p1<-T_result$p.value

  T_p<-c(T_p,p1)
  #T_value<-c(T_value,t1)
  FC<-mean(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]))/mean(as.numeric(data[i,1:nrow(gene1)]))

  FC_result<-c(FC_result,FC)
}


DEG_all<-data.frame(rownames(data),T_p,FC_result,stringsAsFactors = FALSE)

colnames(DEG_all)<-c("Gene","pval","FC") 

write.csv(DEG_all,file = "DEG_all.csv")
DEG<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC>2),]#1295
DEG_down<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC<1),]#1295


expression_matrix <- exprs(cds)


gene_expression <- expression_matrix["TYMP", ]

# gene_expression <- FetchData(cds, vars = "ZNF487", slot = "data")

cds@phenoData@data[["TYMP"]]<-gene_expression
gene_expression2<-0

DEG<-intersect(DEG$Gene,expression_matrix@Dimnames[[1]])
for (i in 1:length(DEG)) {
  gene_expression <- expression_matrix[DEG[i], ]
  gene_expression2<-gene_expression2+gene_expression
}


cds@phenoData@data[["upgene"]]<-gene_expression2




getwd()
pdf("monocle.pseudotime_upgene.pdf",width = 7,height = 7)

plot_cell_trajectory(cds, color_by = "upgene") +
  scale_color_gradient(low = "blue", high = "red") +   
  geom_point(aes(size = Pseudotime, color = upgene), alpha = 0.6) +   
  theme_minimal() +
  theme(legend.position = "right")   

dev.off()

#downgene####
gene_expression2<-0

DEG_down<-intersect(DEG_down$Gene,expression_matrix@Dimnames[[1]])
for (i in 1:length(DEG_down)) {
  gene_expression <- expression_matrix[DEG_down[i], ]
  gene_expression2<-gene_expression2+gene_expression
}

cds@phenoData@data[["downgene"]]<-gene_expression2




getwd()
pdf("monocle.pseudotime_downgene.pdf",width = 7,height = 7)

plot_cell_trajectory(cds, color_by = "downgene") +
  scale_color_gradient(low = "blue", high = "red") +   
  geom_point(aes(size = Pseudotime, color = downgene), alpha = 0.6) +   
  theme_minimal() +
  theme(legend.position = "right")   

dev.off()


library(BayesPrism)
# sc.dat <- pmbc@assays[["RNA"]]
# 
# dim(sc.dat)#33538 17135
# sc.dat<-t(sc.dat)
# cell.type.labels <- pmbc@meta.data[["cell_type"]]
# sc.dat<-as.matrix(sc.dat)


exing<-read.table("G:\\pan_result\\gpu_result_tcga\\KIRC\\scrna\\KIRC_EPI_cells_E.txt")
cell_types <- rep("Malignant epithelial cells",nrow(exing))
cell_names<-exing$x

named_cell_types <- setNames(cell_types, cell_names)


pmbc<-data.combined.new

pmbc@meta.data$cell_type <- as.character(pmbc@meta.data$cell_type)


match_index <- match(rownames(pmbc@meta.data), names(named_cell_types))

pmbc@meta.data$cell_type[!is.na(match_index)] <- named_cell_types[match_index[!is.na(match_index)]]


pmbc@meta.data$cell_type <- factor(pmbc@meta.data$cell_type)



sc.dat <- pmbc[["SCT"]]@counts
dim(sc.dat)#20932 12480
sc.dat<-t(sc.dat)
cell.type.labels <- pmbc@meta.data[["cell_type"]]
sc.dat<-as.matrix(sc.dat)

nrow(sc.dat)*0.3#3744
sc.dat.filtered <- cleanup.genes (input=sc.dat,
                                  input.type="count.matrix",
                                  species="hs",
                                  gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY","hb","act"),
                                  exp.cells=nrow(sc.dat)*0.3)
dim(sc.dat.filtered)#12480  1041

save(sc.dat.filtered,file="sc.dat.filtered_70.Rdata")


rownames(sc.dat.filtered)[1:4]
sc.dat.filtered<-as.data.frame(sc.dat.filtered)
sc.dat.filtered$cellname<-as.character(cell.type.labels)
qiuhe<-aggregate(sc.dat.filtered[,-(ncol(sc.dat.filtered))], by = sc.dat.filtered[ncol(sc.dat.filtered)], FUN = sum)
qiujunzhi<-aggregate(sc.dat.filtered[,-(ncol(sc.dat.filtered))], by = sc.dat.filtered[ncol(sc.dat.filtered)], FUN = mean)
qiuhe<-t(qiuhe)
qiujunzhi<-t(qiujunzhi)
write.table(qiuhe,"LM4_sum.txt",col.names = F,sep = "\t")
write.table(qiujunzhi,"LM4_mean.txt",col.names = F,sep = "\t")
# a<-t(sc.dat.filtered)
# max(a[3,])
# #b<-as.matrix(pmbc@active.ident)
# colnames(a)<-cell.type.labels
# class(a)
# a<-as.data.frame(a)
# b<-sapply(split.default(a, names(a)), rowMeans)
# b<-sapply(split.default(a, names(a)), rowMeans)
# dim(b)#21178     4
# head(b)
# max(b[,1])
 


#KIRP####
# KIRP
library(rlang)
#packageVersion("dplyr")
#install.packages("https://www.bioconductor.org/packages/3.8/bioc/src/contrib/BiocGenerics_0.28.0.tar.gz",repos=NULL,type="source")
library(monocle)
library(monocle)
library(Seurat)
library(SeuratObject)
#load("death_result_2\\GSE135337\\HP\\GSE146771_geneplotneed.Rdata")
#load("result\\GSE146771\\HP\\GSE146771_geneplotneed.Rdata")KIRP
#load("H:\\database\\SCDTD\\data\\chenlong\\Esophageal Cancer\\disease\\GSE196756_geneplotneed.Rdata")
#load("H:\\database\\SCDTD\\data\\chenlong\\Colon Cancer\\disease\\GSE140288_geneplotneed.Rdata")
load("H:\\database\\SCDTD\\data\\chenlong\\Kidney Chromophobe\\disease\\HP\\GSE152938_geneplotneed.Rdata")
GetAssay(data.combined.new,assay = "RNA")
unique(data.combined.new@meta.data[["orig.ident"]])
pmbc<-data.combined.new


Idents(pmbc) <- "orig.ident"
ob1 <- subset(x = pmbc, downsample = 300)

expr_matrix<-as(as.matrix(ob1@assays[["RNA"]]@counts),'sparseMatrix')
expr_matrix<-as(as.matrix(ob1@assays[["SCT"]]@counts),'sparseMatrix')
#p_data<-ob1@meta.data
p_data<-ob1@meta.data
p_data$celltype<-ob1@active.ident
f_data<-data.frame(gene_short_name=row.names(ob1),row.names = row.names(ob1))
dim(expr_matrix)

pd<-new('AnnotatedDataFrame',data=p_data)
fd<-new('AnnotatedDataFrame',data=f_data)





cds<-newCellDataSet(expr_matrix,
                    phenoData=pd,
                    featureData=fd,
                    lowerDetectionLimit=0.1,
                    expressionFamily=negbinomial.size())



cds<-estimateSizeFactors(cds)
#1gc()
# memory.limit(123456)
# rm(expr_matrix)
# rm(f_data)
# rm(p_data)
print(memory.limit())
memory.limit(size=25000)  # 

cds<-estimateDispersions(cds)

cds<-detectGenes(cds,min_expr = 0.1)

print(head(fData(cds)))


gene<-read.csv("G:\\pan_result\\gpu_result_tcga\\KIRP\\444featuresweightme.csv",row.names = 1)
expressed_genes<-colnames(gene)

# expressed_genes = read.table("validategene.csv",header = T,sep = ",")
# expressed_genes = expressed_genes$X0


setwd("G:\\pan_result\\gpu_result_tcga\\KIRP\\scrna")
cds <- setOrderingFilter(cds, expressed_genes)
pdf("train.ordergenes_allcell.pdf")
plot_ordering_genes(cds)
dev.off()



cds <- reduceDimension(cds, max_components = 5,
                       method = 'DDRTree')##method
cds <- orderCells(cds)
setwd("G:\\pan_result\\gpu_result_tcga\\KIRP\\scrna")
pdf("train.monocle.state5_allcell.pdf",width = 7,height = 7)
plot_cell_trajectory(cds, color_by = "Pseudotime",size=1,show_backbone=TRUE)
dev.off()

save(cds,file = "cds_allcell.RData")

pdf("train.monocle.celltype5_allcell.pdf",width = 7,height = 7)
plot_cell_trajectory(cds,color_by = "cell_type",size=1,show_backbone = TRUE)+
  geom_point(aes(size = Pseudotime , color = cell_type)) + 
  scale_color_manual(breaks = c("X", "Y", "Z"), values=c("darkblue", "darkred", "orange","grey","yellow","green","skyblue","pink","purple","#458B74","darkgoldenrod","#CD3333")) 
plot_cell_trajectory(cds,color_by = "cell_type",size=1,show_backbone = TRUE) + scale_color_manual( values=c("darkblue", "darkred", "orange","grey","yellow","green","skyblue","pink","purple","#458B74","darkgoldenrod","#CD3333")) + theme(legend.position = "right")
dev.off()

library(ggpubr)
library(tidyr)
packageVersion("tidyr")

#install.packages("tidyr")

#load("cds.RData")
df <- pData(cds) 

pdf("train.monocle.pseudotime.pdf",width = 7,height = 7)
ggplot(df, aes(Pseudotime, colour = cell_type, fill=cell_type)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2( )
dev.off()
write.csv(df,file = "pData.csv")

library(SingleCellExperiment)
library(Seurat)




# devtools::install_github("alanocallaghan/scater")
library(scater)
data.combined.new@meta.data[["Pseudotime"]]<-cds@phenoData@data[["Pseudotime"]]
sce_gene <- as.SingleCellExperiment(data.combined.new)


label<-read.csv("G:\\pan_result\\gpu_result_tcga\\KIRP\\label.csv")
gene1<-gene[label$X[which(label$label=="1")],]
gene2<-gene[label$X[which(label$label=="2")],]
data<-cbind(t(gene1),t(gene2))
T_p<-c()
FC_result<-c()
for(i in 1:nrow(data)){
  T_result<-wilcox.test(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]),as.numeric(data[i,1:nrow(gene1)]),paired = F,var.equal = T)
  
  p1<-T_result$p.value

  T_p<-c(T_p,p1)
  #T_value<-c(T_value,t1)
  FC<-mean(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]))/mean(as.numeric(data[i,1:nrow(gene1)]))

  FC_result<-c(FC_result,FC)
}

FDR_p<-p.adjust(T_p,method="fdr",length(T_p))#
DEG_all<-data.frame(rownames(data),T_p,FC_result,stringsAsFactors = FALSE)
colnames(DEG_all)<-c("Gene","pval","FC") 


DEG<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC>2),]#1295
# plotExpression(sce_gene, features="NPTX2", color_by = "Pseudotime")#p=2.340444e-09
# 
# pdf("monocle.pseudotime_allcell_gene.pdf",width = 7,height = 7)
# plotExpression(sce_gene, features="NPTX2", color_by = "Pseudotime")#p=2.340444e-09
# 
# dev.off()





expression_matrix <- exprs(cds)


gene_expression <- expression_matrix["PKMYT1", ]#1.397517e-11
# gene_expression <- FetchData(cds, vars = "NPTX2", slot = "data")

cds@phenoData@data[["PKMYT1"]]<-gene_expression

getwd()
pdf("monocle.pseudotime_gene2.pdf",width = 7,height = 7)

plot_cell_trajectory(cds, color_by = "PKMYT1") +
  scale_color_gradient(low = "blue", high = "red") +   
  geom_point(aes(size = Pseudotime, color = PKMYT1), alpha = 0.6) +   
  theme_minimal() +
  theme(legend.position = "right")   

dev.off()

#cluster####
label<-read.table("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\KIRP_cluster.txt")
gene1<-gene[rownames(label)[which(label$clusters=="1")],]
gene2<-gene[rownames(label)[which(label$clusters=="2")],]
data<-cbind(t(gene1),t(gene2))
T_p<-c()
FC_result<-c()
for(i in 1:nrow(data)){
  T_result<-wilcox.test(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]),as.numeric(data[i,1:nrow(gene1)]),paired = F,var.equal = T)
  
  p1<-T_result$p.value

  T_p<-c(T_p,p1)
  #T_value<-c(T_value,t1)
  FC<-mean(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]))/mean(as.numeric(data[i,1:nrow(gene1)]))

  FC_result<-c(FC_result,FC)
}


DEG_all<-data.frame(rownames(data),T_p,FC_result,stringsAsFactors = FALSE)

colnames(DEG_all)<-c("Gene","pval","FC") 

write.csv(DEG_all,file = "DEG_all.csv")
DEG<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC>2),]#1295
DEG_down<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC<0.5),]#1295


expression_matrix <- exprs(cds)


gene_expression <- expression_matrix["TYMP", ]

# gene_expression <- FetchData(cds, vars = "ZNF487", slot = "data")

cds@phenoData@data[["TYMP"]]<-gene_expression
gene_expression2<-0

DEG<-intersect(DEG$Gene,expression_matrix@Dimnames[[1]])
for (i in 1:length(DEG)) {
  gene_expression <- expression_matrix[DEG[i], ]
  gene_expression2<-gene_expression2+gene_expression
}


cds@phenoData@data[["upgene"]]<-gene_expression2




getwd()
pdf("monocle.pseudotime_upgene.pdf",width = 7,height = 7)

plot_cell_trajectory(cds, color_by = "upgene") +
  scale_color_gradient(low = "blue", high = "red") +   
  geom_point(aes(size = Pseudotime, color = upgene), alpha = 0.6) +   
  theme_minimal() +
  theme(legend.position = "right")   

dev.off()

#downgene####
gene_expression2<-0

DEG_down<-intersect(DEG_down$Gene,expression_matrix@Dimnames[[1]])
for (i in 1:length(DEG_down)) {
  gene_expression <- expression_matrix[DEG_down[i], ]
  gene_expression2<-gene_expression2+gene_expression
}

cds@phenoData@data[["downgene"]]<-gene_expression2




getwd()
pdf("monocle.pseudotime_downgene.pdf",width = 7,height = 7)

plot_cell_trajectory(cds, color_by = "downgene") +
  scale_color_gradient(low = "blue", high = "red") +   
  geom_point(aes(size = Pseudotime, color = downgene), alpha = 0.6) +   
  theme_minimal() +
  theme(legend.position = "right")   

dev.off()

library(BayesPrism)
# sc.dat <- pmbc@assays[["RNA"]]
# 
# dim(sc.dat)#33538 17135
# sc.dat<-t(sc.dat)
# cell.type.labels <- pmbc@meta.data[["cell_type"]]
# sc.dat<-as.matrix(sc.dat)


sc.dat <- pmbc[["SCT"]]@counts
dim(sc.dat)#20932 12480
sc.dat<-t(sc.dat)
cell.type.labels <- pmbc@meta.data[["cell_type"]]
sc.dat<-as.matrix(sc.dat)

nrow(sc.dat)*0.3#3744
sc.dat.filtered <- cleanup.genes (input=sc.dat,
                                  input.type="count.matrix",
                                  species="hs",
                                  gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY","hb","act"),
                                  exp.cells=nrow(sc.dat)*0.3)
dim(sc.dat.filtered)#12480  1041

save(sc.dat.filtered,file="sc.dat.filtered_70.Rdata")


rownames(sc.dat.filtered)[1:4]
sc.dat.filtered<-as.data.frame(sc.dat.filtered)
sc.dat.filtered$cellname<-as.character(cell.type.labels)
qiuhe<-aggregate(sc.dat.filtered[,-(ncol(sc.dat.filtered))], by = sc.dat.filtered[ncol(sc.dat.filtered)], FUN = sum)
qiujunzhi<-aggregate(sc.dat.filtered[,-(ncol(sc.dat.filtered))], by = sc.dat.filtered[ncol(sc.dat.filtered)], FUN = mean)
qiuhe<-t(qiuhe)
qiujunzhi<-t(qiujunzhi)
write.table(qiuhe,"LM4_sum.txt",col.names = F,sep = "\t")
write.table(qiujunzhi,"LM4_mean.txt",col.names = F,sep = "\t")
# a<-t(sc.dat.filtered)
# max(a[3,])
# #b<-as.matrix(pmbc@active.ident)
# colnames(a)<-cell.type.labels
# class(a)
# a<-as.data.frame(a)
# b<-sapply(split.default(a, names(a)), rowMeans)
# b<-sapply(split.default(a, names(a)), rowMeans)
# dim(b)#21178     4
# head(b)
# max(b[,1])
 

#LIHC####
# LIHC
library(rlang)
#packageVersion("dplyr")
#install.packages("https://www.bioconductor.org/packages/3.8/bioc/src/contrib/BiocGenerics_0.28.0.tar.gz",repos=NULL,type="source")
library(monocle)
library(monocle)
library(Seurat)
library(SeuratObject)
#load("death_result_2\\GSE135337\\HP\\GSE146771_geneplotneed.Rdata")
#load("result\\GSE146771\\HP\\GSE146771_geneplotneed.Rdata")LIHC
#load("H:\\database\\SCDTD\\data\\chenlong\\Esophageal Cancer\\disease\\GSE196756_geneplotneed.Rdata")
#load("H:\\database\\SCDTD\\data\\chenlong\\Colon Cancer\\disease\\GSE140288_geneplotneed.Rdata")
#load("H:\\database\\SCDTD\\data\\chenlong\\Cervical Cancer\\disease\\GSE168652_geneplotneed.Rdata")
load("H:\\database\\SCDTD\\data\\chenlong\\Liver Cancer\\Liver Cancer\\GSE125449_geneplotneed.Rdata")
GetAssay(data.combined.new,assay = "RNA")
pmbc<-data.combined.new


Idents(pmbc) <- "orig.ident"
ob1 <- subset(x = pmbc, downsample = 300)

expr_matrix<-as(as.matrix(ob1@assays[["RNA"]]@layers[["counts"]]),'sparseMatrix')
expr_matrix<-as(as.matrix(ob1@assays[["SCT"]]@counts),'sparseMatrix')
p_data<-ob1@meta.data
p_data$celltype<-ob1@active.ident
f_data<-data.frame(gene_short_name=row.names(ob1),row.names = row.names(ob1))
dim(expr_matrix)

pd<-new('AnnotatedDataFrame',data=p_data)
fd<-new('AnnotatedDataFrame',data=f_data)





cds<-newCellDataSet(expr_matrix,
                    phenoData=pd,
                    featureData=fd,
                    lowerDetectionLimit=0.1,
                    expressionFamily=negbinomial.size())



cds<-estimateSizeFactors(cds)
#1gc()
# memory.limit(123456)
# rm(expr_matrix)
# rm(f_data)
# rm(p_data)
cds<-estimateDispersions(cds)

cds<-detectGenes(cds,min_expr = 0.1)

print(head(fData(cds)))


gene<-read.csv("G:\\pan_result\\gpu_result_tcga\\LIHC\\642featuresweightme.csv",row.names = 1)
expressed_genes<-colnames(gene)

# expressed_genes = read.table("validategene.csv",header = T,sep = ",")
# expressed_genes = expressed_genes$X0


setwd("G:\\pan_result\\gpu_result_tcga\\LIHC\\scrna")
cds <- setOrderingFilter(cds, expressed_genes)
pdf("train.ordergenes.pdf")
plot_ordering_genes(cds)
dev.off()



cds <- reduceDimension(cds, max_components = 5,
                       method = 'DDRTree')##method
cds <- orderCells(cds)
setwd("G:\\pan_result\\gpu_result_tcga\\LIHC\\scrna")
pdf("train.monocle.state5.pdf",width = 7,height = 7)
plot_cell_trajectory(cds, color_by = "Pseudotime",size=1,show_backbone=TRUE)
dev.off()

save(cds,file = "cds.RData")

pdf("train.monocle.celltype5.pdf",width = 7,height = 7)
plot_cell_trajectory(cds,color_by = "cell_type",size=1,show_backbone = TRUE)+
  geom_point(aes(size = Pseudotime , color = cell_type)) + 
  scale_color_manual(breaks = c("X", "Y", "Z"), values=c("darkblue", "darkred", "orange","grey","yellow","green","skyblue","pink","purple")) 
plot_cell_trajectory(cds,color_by = "cell_type",size=1,show_backbone = TRUE) + scale_color_manual( values=c("darkblue", "darkred", "orange","grey","yellow","green","skyblue","pink","purple")) + theme(legend.position = "right")
dev.off()

library(ggpubr)
library(tidyr)
packageVersion("tidyr")

#install.packages("tidyr")

#load("cds.RData")
df <- pData(cds) 

pdf("train.monocle.pseudotime.pdf",width = 7,height = 7)
ggplot(df, aes(Pseudotime, colour = cell_type, fill=cell_type)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2( )
dev.off()
write.csv(df,file = "pData.csv")


library(SingleCellExperiment)
library(Seurat)




# devtools::install_github("alanocallaghan/scater")
library(scater)
# data.combined.new@meta.data[["Pseudotime"]]<-cds@phenoData@data[["Pseudotime"]]
# sce_gene <- as.SingleCellExperiment(data.combined.new)


label<-read.csv("G:\\pan_result\\gpu_result_tcga\\LIHC\\label.csv")
gene1<-gene[label$X[which(label$label=="1")],]
gene2<-gene[label$X[which(label$label=="2")],]
data<-cbind(t(gene1),t(gene2))
T_p<-c()
FC_result<-c()
for(i in 1:nrow(data)){
  T_result<-wilcox.test(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]),as.numeric(data[i,1:nrow(gene1)]),paired = F,var.equal = T)
  
  p1<-T_result$p.value

  T_p<-c(T_p,p1)
  #T_value<-c(T_value,t1)
  FC<-mean(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]))/mean(as.numeric(data[i,1:nrow(gene1)]))

  FC_result<-c(FC_result,FC)
}


DEG_all<-data.frame(rownames(data),T_p,FC_result,stringsAsFactors = FALSE)

colnames(DEG_all)<-c("Gene","pval","FC") 

write.csv(DEG_all,file = "DEG_all.csv")
#DEG<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC>2),]#1295

expression_matrix <- exprs(cds)


gene_expression <- expression_matrix["GNMT", ]

# gene_expression <- FetchData(cds, vars = "ZNF487", slot = "data")

cds@phenoData@data[["GNMT"]]<-gene_expression


getwd()
pdf("monocle.pseudotime_gene2.pdf",width = 7,height = 7)

plot_cell_trajectory(cds, color_by = "GNMT") +
  scale_color_gradient(low = "blue", high = "red") +   
  geom_point(aes(size = Pseudotime, color = GNMT), alpha = 0.6) +   
  theme_minimal() +
  theme(legend.position = "right")   

dev.off()


#cluster####
label<-read.table("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\LIHC_cluster.txt")
gene1<-gene[rownames(label)[which(label$clusters=="1")],]
gene2<-gene[rownames(label)[which(label$clusters=="2")],]
data<-cbind(t(gene1),t(gene2))
T_p<-c()
FC_result<-c()
for(i in 1:nrow(data)){
  T_result<-wilcox.test(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]),as.numeric(data[i,1:nrow(gene1)]),paired = F,var.equal = T)
  
  p1<-T_result$p.value

  T_p<-c(T_p,p1)
  #T_value<-c(T_value,t1)
  FC<-mean(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]))/mean(as.numeric(data[i,1:nrow(gene1)]))

  FC_result<-c(FC_result,FC)
}


DEG_all<-data.frame(rownames(data),T_p,FC_result,stringsAsFactors = FALSE)

colnames(DEG_all)<-c("Gene","pval","FC") 

write.csv(DEG_all,file = "DEG_all.csv")
DEG<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC>2),]#1295
DEG_down<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC<0.5),]#1295


expression_matrix <- exprs(cds)


gene_expression <- expression_matrix["TYMP", ]

# gene_expression <- FetchData(cds, vars = "ZNF487", slot = "data")

cds@phenoData@data[["TYMP"]]<-gene_expression
gene_expression2<-0

DEG<-intersect(DEG$Gene,expression_matrix@Dimnames[[1]])
for (i in 1:length(DEG)) {
  gene_expression <- expression_matrix[DEG[i], ]
  gene_expression2<-gene_expression2+gene_expression
}


cds@phenoData@data[["upgene"]]<-gene_expression2




getwd()
pdf("monocle.pseudotime_upgene.pdf",width = 7,height = 7)

plot_cell_trajectory(cds, color_by = "upgene") +
  scale_color_gradient(low = "blue", high = "red") +   
  geom_point(aes(size = Pseudotime, color = upgene), alpha = 0.6) +   
  theme_minimal() +
  theme(legend.position = "right")   

dev.off()

#downgene####
gene_expression2<-0

DEG_down<-intersect(DEG_down$Gene,expression_matrix@Dimnames[[1]])
for (i in 1:length(DEG_down)) {
  gene_expression <- expression_matrix[DEG_down[i], ]
  gene_expression2<-gene_expression2+gene_expression
}

cds@phenoData@data[["downgene"]]<-gene_expression2




getwd()
pdf("monocle.pseudotime_downgene.pdf",width = 7,height = 7)

plot_cell_trajectory(cds, color_by = "downgene") +
  scale_color_gradient(low = "blue", high = "red") +   
  geom_point(aes(size = Pseudotime, color = downgene), alpha = 0.6) +   
  theme_minimal() +
  theme(legend.position = "right")   

dev.off()




getwd()



library(BayesPrism)
# sc.dat <- pmbc@assays[["RNA"]]
# 
# dim(sc.dat)#33538 17135
# sc.dat<-t(sc.dat)
# cell.type.labels <- pmbc@meta.data[["cell_type"]]
# sc.dat<-as.matrix(sc.dat)


exing<-read.table("G:\\pan_result\\gpu_result_tcga\\LIHC\\scrna\\LIHC_EPI_cells_E.txt")
cell_types <- rep("Malignant epithelial cells",nrow(exing))
cell_names<-exing$x

named_cell_types <- setNames(cell_types, cell_names)


pmbc<-data.combined.new

pmbc@meta.data$cell_type <- as.character(pmbc@meta.data$cell_type)


match_index <- match(rownames(pmbc@meta.data), names(named_cell_types))

pmbc@meta.data$cell_type[!is.na(match_index)] <- named_cell_types[match_index[!is.na(match_index)]]


pmbc@meta.data$cell_type <- factor(pmbc@meta.data$cell_type)


sc.dat <- pmbc[["SCT"]]@counts
dim(sc.dat)#20932 12480
sc.dat<-t(sc.dat)
cell.type.labels <- pmbc@meta.data[["cell_type"]]
sc.dat<-as.matrix(sc.dat)

nrow(sc.dat)*0.3#3744
sc.dat.filtered <- cleanup.genes (input=sc.dat,
                                  input.type="count.matrix",
                                  species="hs",
                                  gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY","hb","act"),
                                  exp.cells=nrow(sc.dat)*0.3)
dim(sc.dat.filtered)#12480  1041

save(sc.dat.filtered,file="sc.dat.filtered_70.Rdata")


rownames(sc.dat.filtered)[1:4]
sc.dat.filtered<-as.data.frame(sc.dat.filtered)
sc.dat.filtered$cellname<-as.character(cell.type.labels)
qiuhe<-aggregate(sc.dat.filtered[,-(ncol(sc.dat.filtered))], by = sc.dat.filtered[ncol(sc.dat.filtered)], FUN = sum)
qiujunzhi<-aggregate(sc.dat.filtered[,-(ncol(sc.dat.filtered))], by = sc.dat.filtered[ncol(sc.dat.filtered)], FUN = mean)
qiuhe<-t(qiuhe)
qiujunzhi<-t(qiujunzhi)
write.table(qiuhe,"LM4_sum.txt",col.names = F,sep = "\t")
write.table(qiujunzhi,"LM4_mean.txt",col.names = F,sep = "\t")
# a<-t(sc.dat.filtered)
# max(a[3,])
# #b<-as.matrix(pmbc@active.ident)
# colnames(a)<-cell.type.labels
# class(a)
# a<-as.data.frame(a)
# b<-sapply(split.default(a, names(a)), rowMeans)
# b<-sapply(split.default(a, names(a)), rowMeans)
# dim(b)#21178     4
# head(b)
# max(b[,1])
 

#LUNG####
# LUNG
library(rlang)
#packageVersion("dplyr")
#install.packages("https://www.bioconductor.org/packages/3.8/bioc/src/contrib/BiocGenerics_0.28.0.tar.gz",repos=NULL,type="source")
library(monocle)
library(monocle)
library(Seurat)
library(SeuratObject)
#load("death_result_2\\GSE135337\\HP\\GSE146771_geneplotneed.Rdata")
#load("result\\GSE146771\\HP\\GSE146771_geneplotneed.Rdata")LUNG
#load("H:\\database\\SCDTD\\data\\chenlong\\Esophageal Cancer\\disease\\GSE196756_geneplotneed.Rdata")
#load("H:\\database\\SCDTD\\data\\chenlong\\Colon Cancer\\disease\\GSE140288_geneplotneed.Rdata")
#load("H:\\database\\SCDTD\\data\\chenlong\\Cervical Cancer\\disease\\GSE168652_geneplotneed.Rdata")
#load("H:\\database\\SCDTD\\data\\chenlong\\Liver Cancer\\Liver Cancer\\GSE125449_geneplotneed.Rdata")
#load("H:\\database\\SCDTD\\data\\chenlong\\Thyroid Cancer\\disease\\GSE191288_geneplotneed.Rdata")
#load("death_result_1\\GSE117570\\HP\\GSE117570_geneplotneed.Rdata")
load("H:\\database\\SCDTD\\data\\chenlong\\Liver Cancer\\COVID-19\\GSE148881_geneplotneed.Rdata")
GetAssay(data.combined.new,assay = "RNA")
unique(data.combined.new@meta.data[["orig.ident"]])
pmbc<-data.combined.new


Idents(pmbc) <- "orig.ident"
ob1 <- subset(x = pmbc, downsample = 300)

expr_matrix<-as(as.matrix(ob1@assays[["RNA"]]@layers[["counts"]]),'sparseMatrix')
expr_matrix<-as(as.matrix(ob1@assays[["SCT"]]@counts),'sparseMatrix')
p_data<-ob1@meta.data
p_data$celltype<-ob1@active.ident
f_data<-data.frame(gene_short_name=row.names(ob1),row.names = row.names(ob1))
dim(expr_matrix)

pd<-new('AnnotatedDataFrame',data=p_data)
fd<-new('AnnotatedDataFrame',data=f_data)





cds<-newCellDataSet(expr_matrix,
                    phenoData=pd,
                    featureData=fd,
                    lowerDetectionLimit=0.1,
                    expressionFamily=negbinomial.size())



cds<-estimateSizeFactors(cds)
#1gc()
# memory.limit(123456)
# rm(expr_matrix)
# rm(f_data)
# rm(p_data)
cds<-estimateDispersions(cds)

cds<-detectGenes(cds,min_expr = 0.1)

print(head(fData(cds)))


gene<-read.csv("G:\\pan_result\\gpu_result_tcga\\LUNG\\957featuresweightme.csv",row.names = 1)
expressed_genes<-colnames(gene)

# expressed_genes = read.table("validategene.csv",header = T,sep = ",")
# expressed_genes = expressed_genes$X0


setwd("G:\\pan_result\\gpu_result_tcga\\LUNG\\scrna")
cds <- setOrderingFilter(cds, expressed_genes)
pdf("train.ordergenes.pdf")
plot_ordering_genes(cds)
dev.off()



cds <- reduceDimension(cds, max_components = 5,
                       method = 'DDRTree')##method
cds <- orderCells(cds)
setwd("G:\\pan_result\\gpu_result_tcga\\LUNG\\scrna")
pdf("train.monocle.state5.pdf",width = 7,height = 7)
plot_cell_trajectory(cds, color_by = "Pseudotime",size=1,show_backbone=TRUE)
dev.off()

save(cds,file = "cds.RData")

pdf("train.monocle.celltype5.pdf",width = 7,height = 7)
plot_cell_trajectory(cds,color_by = "cell_type",size=1,show_backbone = TRUE)+
  geom_point(aes(size = Pseudotime , color = cell_type)) + 
  scale_color_manual(breaks = c("X", "Y", "Z"), values=c("darkblue", "darkred", "orange","grey","yellow","green","skyblue","pink","purple")) 
plot_cell_trajectory(cds,color_by = "cell_type",size=1,show_backbone = TRUE) + scale_color_manual( values=c("darkblue", "darkred", "orange","grey","yellow","green","skyblue","pink","purple")) + theme(legend.position = "right")
dev.off()

library(ggpubr)
library(tidyr)
packageVersion("tidyr")

#install.packages("tidyr")

#load("cds.RData")
df <- pData(cds) 

pdf("train.monocle.pseudotime.pdf",width = 7,height = 7)
ggplot(df, aes(Pseudotime, colour = cell_type, fill=cell_type)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2( )
dev.off()
write.csv(df,file = "pData.csv")


library(SingleCellExperiment)
library(Seurat)




# devtools::install_github("alanocallaghan/scater")
library(scater)
data.combined.new@meta.data[["Pseudotime"]]<-cds@phenoData@data[["Pseudotime"]]
sce_gene <- as.SingleCellExperiment(data.combined.new)


label<-read.csv("G:\\pan_result\\gpu_result_tcga\\LUNG\\label.csv")
gene1<-gene[label$X[which(label$label=="1")],]
gene2<-gene[label$X[which(label$label=="2")],]
data<-cbind(t(gene1),t(gene2))
T_p<-c()
FC_result<-c()
for(i in 1:nrow(data)){
  T_result<-wilcox.test(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]),as.numeric(data[i,1:nrow(gene1)]),paired = F,var.equal = T)
  
  p1<-T_result$p.value

  T_p<-c(T_p,p1)
  #T_value<-c(T_value,t1)
  FC<-mean(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]))/mean(as.numeric(data[i,1:nrow(gene1)]))

  FC_result<-c(FC_result,FC)
}

FDR_p<-p.adjust(T_p,method="fdr",length(T_p))#
DEG_all<-data.frame(rownames(data),T_p,FC_result,stringsAsFactors = FALSE)
colnames(DEG_all)<-c("Gene","pval","FC") 


DEG<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC>2),]#1295
# plotExpression(sce_gene, features="NPTX2", color_by = "Pseudotime")#p=2.340444e-09
# 
# pdf("monocle.pseudotime_allcell_gene.pdf",width = 7,height = 7)
# plotExpression(sce_gene, features="NPTX2", color_by = "Pseudotime")#p=2.340444e-09
# 
# dev.off()





expression_matrix <- exprs(cds)


gene_expression <- expression_matrix["MIF", ]#1.397517e-11
# gene_expression <- FetchData(cds, vars = "NPTX2", slot = "data")

cds@phenoData@data[["MIF"]]<-gene_expression

getwd()
pdf("monocle.pseudotime_gene2.pdf",width = 7,height = 7)

plot_cell_trajectory(cds, color_by = "MIF") +
  scale_color_gradient(low = "blue", high = "red") +   
  geom_point(aes(size = Pseudotime, color = MIF), alpha = 0.6) +   
  theme_minimal() +
  theme(legend.position = "right")   

dev.off()


#cluster####
label<-read.table("G:\\pan_result\\gpu_result_tcga\\LUNG\\450k\\LUNG_cluster.txt")
gene1<-gene[rownames(label)[which(label$clusters=="1")],]
gene2<-gene[rownames(label)[which(label$clusters=="2")],]
data<-cbind(t(gene1),t(gene2))
T_p<-c()
FC_result<-c()
for(i in 1:nrow(data)){
  T_result<-wilcox.test(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]),as.numeric(data[i,1:nrow(gene1)]),paired = F,var.equal = T)
  
  p1<-T_result$p.value

  T_p<-c(T_p,p1)
  #T_value<-c(T_value,t1)
  FC<-mean(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]))/mean(as.numeric(data[i,1:nrow(gene1)]))

  FC_result<-c(FC_result,FC)
}


DEG_all<-data.frame(rownames(data),T_p,FC_result,stringsAsFactors = FALSE)

colnames(DEG_all)<-c("Gene","pval","FC") 

write.csv(DEG_all,file = "DEG_all.csv")
DEG<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC>2),]#1295
DEG_down<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC<0.5),]#1295


expression_matrix <- exprs(cds)


gene_expression <- expression_matrix["TYMP", ]

# gene_expression <- FetchData(cds, vars = "ZNF487", slot = "data")

cds@phenoData@data[["TYMP"]]<-gene_expression
gene_expression2<-0

DEG<-intersect(DEG$Gene,expression_matrix@Dimnames[[1]])
for (i in 1:length(DEG)) {
  gene_expression <- expression_matrix[DEG[i], ]
  gene_expression2<-gene_expression2+gene_expression
}


cds@phenoData@data[["upgene"]]<-gene_expression2




getwd()
pdf("monocle.pseudotime_upgene.pdf",width = 7,height = 7)

plot_cell_trajectory(cds, color_by = "upgene") +
  scale_color_gradient(low = "blue", high = "red") +   
  geom_point(aes(size = Pseudotime, color = upgene), alpha = 0.6) +   
  theme_minimal() +
  theme(legend.position = "right")   

dev.off()

#downgene####
gene_expression2<-0

DEG_down<-intersect(DEG_down$Gene,expression_matrix@Dimnames[[1]])
for (i in 1:length(DEG_down)) {
  gene_expression <- expression_matrix[DEG_down[i], ]
  gene_expression2<-gene_expression2+gene_expression
}

cds@phenoData@data[["downgene"]]<-gene_expression2




getwd()
pdf("monocle.pseudotime_downgene.pdf",width = 7,height = 7)

plot_cell_trajectory(cds, color_by = "downgene") +
  scale_color_gradient(low = "blue", high = "red") +   
  geom_point(aes(size = Pseudotime, color = downgene), alpha = 0.6) +   
  theme_minimal() +
  theme(legend.position = "right")   

dev.off()





library(BayesPrism)
# sc.dat <- pmbc@assays[["RNA"]]
# 
# dim(sc.dat)#33538 17135
# sc.dat<-t(sc.dat)
# cell.type.labels <- pmbc@meta.data[["cell_type"]]
# sc.dat<-as.matrix(sc.dat)

exing<-read.table("G:\\pan_result\\gpu_result_tcga\\LUNG\\scrna\\LUNG_EPI_cells_E.txt")
cell_types <- rep("Malignant epithelial cells",nrow(exing))
cell_names<-exing$x

named_cell_types <- setNames(cell_types, cell_names)


pmbc<-data.combined.new

pmbc@meta.data$cell_type <- as.character(pmbc@meta.data$cell_type)


match_index <- match(rownames(pmbc@meta.data), names(named_cell_types))

pmbc@meta.data$cell_type[!is.na(match_index)] <- named_cell_types[match_index[!is.na(match_index)]]


pmbc@meta.data$cell_type <- factor(pmbc@meta.data$cell_type)


sc.dat <- pmbc[["SCT"]]@counts
dim(sc.dat)#20932 12480
sc.dat<-t(sc.dat)
cell.type.labels <- pmbc@meta.data[["cell_type"]]
sc.dat<-as.matrix(sc.dat)

nrow(sc.dat)*0.3#3744
sc.dat.filtered <- cleanup.genes (input=sc.dat,
                                  input.type="count.matrix",
                                  species="hs",
                                  gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY","hb","act"),
                                  exp.cells=nrow(sc.dat)*0.3)
dim(sc.dat.filtered)#12480  1041

save(sc.dat.filtered,file="sc.dat.filtered_70.Rdata")


rownames(sc.dat.filtered)[1:4]
sc.dat.filtered<-as.data.frame(sc.dat.filtered) 
sc.dat.filtered$cellname<-as.character(cell.type.labels)
qiuhe<-aggregate(sc.dat.filtered[,-(ncol(sc.dat.filtered))], by = sc.dat.filtered[ncol(sc.dat.filtered)], FUN = sum)
qiujunzhi<-aggregate(sc.dat.filtered[,-(ncol(sc.dat.filtered))], by = sc.dat.filtered[ncol(sc.dat.filtered)], FUN = mean)
qiuhe<-t(qiuhe)
qiujunzhi<-t(qiujunzhi)
write.table(qiuhe,"LM4_sum.txt",col.names = F,sep = "\t")
write.table(qiujunzhi,"LM4_mean.txt",col.names = F,sep = "\t")
# a<-t(sc.dat.filtered)
# max(a[3,])
# #b<-as.matrix(pmbc@active.ident)
# colnames(a)<-cell.type.labels
# class(a)
# a<-as.data.frame(a)
# b<-sapply(split.default(a, names(a)), rowMeans)
# b<-sapply(split.default(a, names(a)), rowMeans)
# dim(b)#21178     4
# head(b)
# max(b[,1])
 


#PAAD####
# PAAD
library(rlang)
#packageVersion("dplyr")
#install.packages("https://www.bioconductor.org/packages/3.8/bioc/src/contrib/BiocGenerics_0.28.0.tar.gz",repos=NULL,type="source")
library(monocle)
library(monocle)
library(Seurat)
library(SeuratObject)
#load("death_result_2\\GSE135337\\HP\\GSE146771_geneplotneed.Rdata")
#load("result\\GSE146771\\HP\\GSE146771_geneplotneed.Rdata")PAAD
#load("H:\\database\\SCDTD\\data\\chenlong\\Esophageal Cancer\\disease\\GSE196756_geneplotneed.Rdata")
#load("H:\\database\\SCDTD\\data\\chenlong\\Colon Cancer\\disease\\GSE140288_geneplotneed.Rdata")
#load("H:\\database\\SCDTD\\data\\chenlong\\Cervical Cancer\\disease\\GSE168652_geneplotneed.Rdata")
#load("H:\\database\\SCDTD\\data\\chenlong\\Liver Cancer\\Liver Cancer\\GSE125449_geneplotneed.Rdata")
# load("H:\\database\\SCDTD\\data\\result\\GSE263365_HP\\GSE263365_geneplotneed.Rdata")
load("G:\\pan_result\\gpu_result_tcga\\ESCA\\scrna\\harmony_0.4_20\\singleR\\singleR\\GSM5032772_geneplotneed.Rdata")
GetAssay(data.combined.new,assay = "RNA")
unique(data.combined.new@meta.data[["orig.ident"]])
pmbc<-data.combined.new


Idents(pmbc) <- "orig.ident"
ob1 <- subset(x = pmbc, downsample = 300)

expr_matrix<-as(as.matrix(ob1@assays[["RNA"]]@layers[["counts"]]),'sparseMatrix')
expr_matrix<-as(as.matrix(ob1@assays[["SCT"]]@counts),'sparseMatrix')
p_data<-ob1@meta.data
p_data$celltype<-ob1@active.ident
f_data<-data.frame(gene_short_name=row.names(ob1),row.names = row.names(ob1))
dim(expr_matrix)

pd<-new('AnnotatedDataFrame',data=p_data)
fd<-new('AnnotatedDataFrame',data=f_data)





cds<-newCellDataSet(expr_matrix,
                    phenoData=pd,
                    featureData=fd,
                    lowerDetectionLimit=0.1,
                    expressionFamily=negbinomial.size())



cds<-estimateSizeFactors(cds)
#1gc()
# memory.limit(123456)
# rm(expr_matrix)
# rm(f_data)
# rm(p_data)
cds<-estimateDispersions(cds)

cds<-detectGenes(cds,min_expr = 0.1)

print(head(fData(cds)))


gene<-read.csv("G:\\pan_result\\gpu_result_tcga\\PAAD\\67featuresweightme.csv",row.names = 1)
expressed_genes<-colnames(gene)

# expressed_genes = read.table("validategene.csv",header = T,sep = ",")
# expressed_genes = expressed_genes$X0


setwd("G:\\pan_result\\gpu_result_tcga\\PAAD\\scrna")
cds <- setOrderingFilter(cds, expressed_genes)
pdf("train.ordergenes.pdf")
plot_ordering_genes(cds)
dev.off()



cds <- reduceDimension(cds, max_components = 5,
                       method = 'DDRTree')##method
cds <- orderCells(cds)
setwd("G:\\pan_result\\gpu_result_tcga\\PAAD\\scrna")
pdf("train.monocle.state5.pdf",width = 7,height = 7)
plot_cell_trajectory(cds, color_by = "Pseudotime",size=1,show_backbone=TRUE)
dev.off()

save(cds,file = "cds.RData")

pdf("train.monocle.celltype5.pdf",width = 7,height = 7)
plot_cell_trajectory(cds,color_by = "cell_type",size=1,show_backbone = TRUE)+
  geom_point(aes(size = Pseudotime , color = cell_type)) + 
  scale_color_manual(breaks = c("X", "Y", "Z"), values=c("darkblue", "darkred", "orange","grey","yellow","green","skyblue","pink","purple")) 
plot_cell_trajectory(cds,color_by = "cell_type",size=1,show_backbone = TRUE) + scale_color_manual( values=c("darkblue", "darkred", "orange","grey","yellow","green","skyblue","pink","purple")) + theme(legend.position = "right")
dev.off()

library(ggpubr)
library(tidyr)
packageVersion("tidyr")

#install.packages("tidyr")

#load("cds.RData")
df <- pData(cds) 

pdf("train.monocle.pseudotime.pdf",width = 7,height = 7)
ggplot(df, aes(Pseudotime, colour = cell_type, fill=cell_type)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2( )
dev.off()
write.csv(df,file = "pData.csv")

#####
# 
library(SingleCellExperiment)
library(Seurat)




# devtools::install_github("alanocallaghan/scater")
library(scater)
# data.combined.new@meta.data[["Pseudotime"]]<-cds@phenoData@data[["Pseudotime"]]
# sce_gene <- as.SingleCellExperiment(data.combined.new)


label<-read.csv("G:\\pan_result\\gpu_result_tcga\\PAAD\\label.csv")
gene1<-gene[label$X[which(label$label=="1")],]
gene2<-gene[label$X[which(label$label=="2")],]
data<-cbind(t(gene1),t(gene2))
T_p<-c()
FC_result<-c()
for(i in 1:nrow(data)){
  T_result<-wilcox.test(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]),as.numeric(data[i,1:nrow(gene1)]),paired = F,var.equal = T)
  
  p1<-T_result$p.value

  T_p<-c(T_p,p1)
  #T_value<-c(T_value,t1)
  FC<-mean(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]))/mean(as.numeric(data[i,1:nrow(gene1)]))

  FC_result<-c(FC_result,FC)
}

FDR_p<-p.adjust(T_p,method="fdr",length(T_p))#
DEG_all<-data.frame(rownames(data),T_p,FC_result,stringsAsFactors = FALSE)
colnames(DEG_all)<-c("Gene","pval","FC") 


DEG<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC>2),]#1295

getwd()
# pdf("monocle.pseudotime_allcell_gene2.pdf",width = 7,height = 7)
# 
# plot_cell_trajectory(cds, color_by = "NPTX2") +
#   scale_color_gradient(low = "blue", high = "red") +   
#   geom_point(aes(size = Pseudotime, color = NPTX2), alpha = 0.6) +   
#   theme_minimal() +
#   theme(legend.position = "right")   
# 
# dev.off()



expression_matrix <- exprs(cds)


gene_expression <- expression_matrix["ACTR3B", ]#2.340444e-09

# gene_expression <- FetchData(cds, vars = "NPTX2", slot = "data")

cds@phenoData@data[["ACTR3B"]]<-gene_expression

getwd()
pdf("monocle.pseudotime_gene2.pdf",width = 7,height = 7)

plot_cell_trajectory(cds, color_by = "ACTR3B") +
  scale_color_gradient(low = "blue", high = "red") +   
  geom_point(aes(size = Pseudotime, color = ACTR3B), alpha = 0.6) +   
  theme_minimal() +
  theme(legend.position = "right")   

dev.off()

#cluster####
label<-read.table("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\PAAD_cluster.txt")
gene1<-gene[rownames(label)[which(label$clusters=="1")],]
gene2<-gene[rownames(label)[which(label$clusters=="2")],]
data<-cbind(t(gene1),t(gene2))
T_p<-c()
FC_result<-c()
for(i in 1:nrow(data)){
  T_result<-wilcox.test(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]),as.numeric(data[i,1:nrow(gene1)]),paired = F,var.equal = T)
  
  p1<-T_result$p.value

  T_p<-c(T_p,p1)
  #T_value<-c(T_value,t1)
  FC<-mean(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]))/mean(as.numeric(data[i,1:nrow(gene1)]))

  FC_result<-c(FC_result,FC)
}


DEG_all<-data.frame(rownames(data),T_p,FC_result,stringsAsFactors = FALSE)

colnames(DEG_all)<-c("Gene","pval","FC") 

write.csv(DEG_all,file = "DEG_all.csv")
DEG<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC>1),]#1295
DEG_down<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC<1),]#1295


expression_matrix <- exprs(cds)


gene_expression <- expression_matrix["TYMP", ]

# gene_expression <- FetchData(cds, vars = "ZNF487", slot = "data")

cds@phenoData@data[["TYMP"]]<-gene_expression
gene_expression2<-0

DEG<-intersect(DEG$Gene,expression_matrix@Dimnames[[1]])
for (i in 1:length(DEG)) {
  gene_expression <- expression_matrix[DEG[i], ]
  gene_expression2<-gene_expression2+gene_expression
}


cds@phenoData@data[["upgene"]]<-gene_expression2




getwd()
pdf("monocle.pseudotime_upgene.pdf",width = 7,height = 7)

plot_cell_trajectory(cds, color_by = "upgene") +
  scale_color_gradient(low = "blue", high = "red") +   
  geom_point(aes(size = Pseudotime, color = upgene), alpha = 0.6) +   
  theme_minimal() +
  theme(legend.position = "right")   

dev.off()

#downgene####
gene_expression2<-0

DEG_down<-intersect(DEG_down$Gene,expression_matrix@Dimnames[[1]])
for (i in 1:length(DEG_down)) {
  gene_expression <- expression_matrix[DEG_down[i], ]
  gene_expression2<-gene_expression2+gene_expression
}

cds@phenoData@data[["downgene"]]<-gene_expression2




getwd()
pdf("monocle.pseudotime_downgene.pdf",width = 7,height = 7)

plot_cell_trajectory(cds, color_by = "downgene") +
  scale_color_gradient(low = "blue", high = "red") +   
  geom_point(aes(size = Pseudotime, color = downgene), alpha = 0.6) +   
  theme_minimal() +
  theme(legend.position = "right")   

dev.off()





library(BayesPrism)
# sc.dat <- pmbc@assays[["RNA"]]
# 
# dim(sc.dat)#33538 17135
# sc.dat<-t(sc.dat)
# cell.type.labels <- pmbc@meta.data[["cell_type"]]
# sc.dat<-as.matrix(sc.dat)

exing<-read.table("G:\\pan_result\\gpu_result_tcga\\PAAD\\scrna\\PAAD_EPI_cells_E.txt")
cell_types <- rep("Malignant epithelial cells",nrow(exing))
cell_names<-exing$x

named_cell_types <- setNames(cell_types, cell_names)


pmbc<-data.combined.new

pmbc@meta.data$cell_type <- as.character(pmbc@meta.data$cell_type)


match_index <- match(rownames(pmbc@meta.data), names(named_cell_types))

pmbc@meta.data$cell_type[!is.na(match_index)] <- named_cell_types[match_index[!is.na(match_index)]]


pmbc@meta.data$cell_type <- factor(pmbc@meta.data$cell_type)



sc.dat <- pmbc[["SCT"]]@counts
sc.dat <- pmbc@assays[["RNA"]]@layers[["counts"]]
dim(sc.dat)#22393  1849
sc.dat<-t(sc.dat)
cell.type.labels <- pmbc@meta.data[["cell_type"]]
sc.dat<-as.matrix(sc.dat)

nrow(sc.dat)*0.3#3744
sc.dat.filtered <- cleanup.genes (input=sc.dat,
                                  input.type="count.matrix",
                                  species="hs",
                                  gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY","hb","act"),
                                  exp.cells=nrow(sc.dat)*0.3)
dim(sc.dat.filtered)#12480  1041

save(sc.dat.filtered,file="sc.dat.filtered_70.Rdata")


rownames(sc.dat.filtered)[1:4]
sc.dat.filtered<-as.data.frame(sc.dat.filtered) 
sc.dat.filtered$cellname<-as.character(cell.type.labels)
qiuhe<-aggregate(sc.dat.filtered[,-(ncol(sc.dat.filtered))], by = sc.dat.filtered[ncol(sc.dat.filtered)], FUN = sum)
qiujunzhi<-aggregate(sc.dat.filtered[,-(ncol(sc.dat.filtered))], by = sc.dat.filtered[ncol(sc.dat.filtered)], FUN = mean)
qiuhe<-t(qiuhe)
qiujunzhi<-t(qiujunzhi)
write.table(qiuhe,"LM4_sum.txt",col.names = F,sep = "\t")
write.table(qiujunzhi,"LM4_mean.txt",col.names = F,sep = "\t")
# a<-t(sc.dat.filtered)
# max(a[3,])
# #b<-as.matrix(pmbc@active.ident)
# colnames(a)<-cell.type.labels
# class(a)
# a<-as.data.frame(a)
# b<-sapply(split.default(a, names(a)), rowMeans)
# b<-sapply(split.default(a, names(a)), rowMeans)
# dim(b)#21178     4
# head(b)
# max(b[,1])
 

#THCA####
# THCA
library(rlang)
#packageVersion("dplyr")
#install.packages("https://www.bioconductor.org/packages/3.8/bioc/src/contrib/BiocGenerics_0.28.0.tar.gz",repos=NULL,type="source")
library(monocle)
library(monocle)
library(Seurat)
library(SeuratObject)
#load("death_result_2\\GSE135337\\HP\\GSE146771_geneplotneed.Rdata")
#load("result\\GSE146771\\HP\\GSE146771_geneplotneed.Rdata")THCA
#load("H:\\database\\SCDTD\\data\\chenlong\\Esophageal Cancer\\disease\\GSE196756_geneplotneed.Rdata")
#load("H:\\database\\SCDTD\\data\\chenlong\\Colon Cancer\\disease\\GSE140288_geneplotneed.Rdata")
#load("H:\\database\\SCDTD\\data\\chenlong\\Cervical Cancer\\disease\\GSE168652_geneplotneed.Rdata")
#load("H:\\database\\SCDTD\\data\\chenlong\\Liver Cancer\\Liver Cancer\\GSE125449_geneplotneed.Rdata")
load("H:\\database\\SCDTD\\data\\chenlong\\Thyroid Cancer\\disease\\GSE191288_geneplotneed.Rdata")
GetAssay(data.combined.new,assay = "RNA")
unique(data.combined.new@meta.data[["orig.ident"]])
pmbc<-data.combined.new


Idents(pmbc) <- "orig.ident"
ob1 <- subset(x = pmbc, downsample = 300)

expr_matrix<-as(as.matrix(ob1@assays[["RNA"]]@layers[["counts"]]),'sparseMatrix')
expr_matrix<-as(as.matrix(ob1@assays[["SCT"]]@counts),'sparseMatrix')
p_data<-ob1@meta.data
p_data$celltype<-ob1@active.ident
f_data<-data.frame(gene_short_name=row.names(ob1),row.names = row.names(ob1))
dim(expr_matrix)

pd<-new('AnnotatedDataFrame',data=p_data)
fd<-new('AnnotatedDataFrame',data=f_data)





cds<-newCellDataSet(expr_matrix,
                    phenoData=pd,
                    featureData=fd,
                    lowerDetectionLimit=0.1,
                    expressionFamily=negbinomial.size())



cds<-estimateSizeFactors(cds)
#1gc()
# memory.limit(123456)
# rm(expr_matrix)
# rm(f_data)
# rm(p_data)
cds<-estimateDispersions(cds)

cds<-detectGenes(cds,min_expr = 0.1)

print(head(fData(cds)))


gene<-read.csv("G:\\pan_result\\gpu_result_tcga\\THCA\\448featuresweightme.csv",row.names = 1)
expressed_genes<-colnames(gene)

# expressed_genes = read.table("validategene.csv",header = T,sep = ",")
# expressed_genes = expressed_genes$X0


setwd("G:\\pan_result\\gpu_result_tcga\\THCA\\scrna")
cds <- setOrderingFilter(cds, expressed_genes)
pdf("train.ordergenes.pdf")
plot_ordering_genes(cds)
dev.off()



cds <- reduceDimension(cds, max_components = 5,
                       method = 'DDRTree')##method
cds <- orderCells(cds)
setwd("G:\\pan_result\\gpu_result_tcga\\THCA\\scrna")
pdf("train.monocle.state5.pdf",width = 7,height = 7)
plot_cell_trajectory(cds, color_by = "Pseudotime",size=1,show_backbone=TRUE)
dev.off()

save(cds,file = "cds.RData")

pdf("train.monocle.celltype5.pdf",width = 7,height = 7)
plot_cell_trajectory(cds,color_by = "cell_type",size=1,show_backbone = TRUE)+
  geom_point(aes(size = Pseudotime , color = cell_type)) + 
  scale_color_manual(breaks = c("X", "Y", "Z"), values=c("darkblue", "darkred", "orange","grey","yellow","green","skyblue","pink","purple")) 
plot_cell_trajectory(cds,color_by = "cell_type",size=1,show_backbone = TRUE) + scale_color_manual( values=c("darkblue", "darkred", "orange","grey","yellow","green","skyblue","pink","purple")) + theme(legend.position = "right")
dev.off()

library(ggpubr)
library(tidyr)
packageVersion("tidyr")

#install.packages("tidyr")

#load("cds.RData")
df <- pData(cds) 

pdf("train.monocle.pseudotime.pdf",width = 7,height = 7)
ggplot(df, aes(Pseudotime, colour = cell_type, fill=cell_type)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2( )
dev.off()
write.csv(df,file = "pData.csv")



library(SingleCellExperiment)
library(Seurat)




# devtools::install_github("alanocallaghan/scater")
library(scater)
# data.combined.new@meta.data[["Pseudotime"]]<-cds@phenoData@data[["Pseudotime"]]
# sce_gene <- as.SingleCellExperiment(data.combined.new)


label<-read.table("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\THCA_cluster.txt")
gene1<-gene[rownames(label)[which(label$clusters=="1")],]
gene2<-gene[rownames(label)[which(label$clusters=="2")],]
data<-cbind(t(gene1),t(gene2))
T_p<-c()
FC_result<-c()
for(i in 1:nrow(data)){
  T_result<-wilcox.test(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]),as.numeric(data[i,1:nrow(gene1)]),paired = F,var.equal = T)
  
  p1<-T_result$p.value

  T_p<-c(T_p,p1)
  #T_value<-c(T_value,t1)
  FC<-mean(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]))/mean(as.numeric(data[i,1:nrow(gene1)]))

  FC_result<-c(FC_result,FC)
}


DEG_all<-data.frame(rownames(data),T_p,FC_result,stringsAsFactors = FALSE)

colnames(DEG_all)<-c("Gene","pval","FC") 

write.csv(DEG_all,file = "DEG_all.csv")
DEG<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC>1),]#1295
DEG_down<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC<1),]#1295


expression_matrix <- exprs(cds)


gene_expression <- expression_matrix["TYMP", ]

# gene_expression <- FetchData(cds, vars = "ZNF487", slot = "data")

cds@phenoData@data[["TYMP"]]<-gene_expression
gene_expression2<-0
gene_expression2<-gene_expression2+gene_expression

for (i in 1:nrow(DEG)) {
  gene_expression <- expression_matrix[DEG$Gene[i], ]
  gene_expression2<-gene_expression2+gene_expression
}

for (i in 31:nrow(DEG)) {
  gene_expression <- expression_matrix[DEG$Gene[i], ]
  gene_expression2<-gene_expression2+gene_expression
}
cds@phenoData@data[["upgene"]]<-gene_expression2




getwd()
pdf("monocle.pseudotime_upgene.pdf",width = 7,height = 7)

plot_cell_trajectory(cds, color_by = "upgene") +
  scale_color_gradient(low = "blue", high = "red") +   
  geom_point(aes(size = Pseudotime, color = upgene), alpha = 0.6) +   
  theme_minimal() +
  theme(legend.position = "right")   

dev.off()

#downgene####
gene_expression2<-0

for (i in 1:nrow(DEG_down)) {
  gene_expression <- expression_matrix[DEG_down$Gene[i], ]
  gene_expression2<-gene_expression2+gene_expression
}

for (i in 192:nrow(DEG_down)) {
  gene_expression <- expression_matrix[DEG_down$Gene[i], ]
  gene_expression2<-gene_expression2+gene_expression
}
cds@phenoData@data[["downgene"]]<-gene_expression2




getwd()
pdf("monocle.pseudotime_downgene.pdf",width = 7,height = 7)

plot_cell_trajectory(cds, color_by = "downgene") +
  scale_color_gradient(low = "blue", high = "red") +   
  geom_point(aes(size = Pseudotime, color = downgene), alpha = 0.6) +   
  theme_minimal() +
  theme(legend.position = "right")   

dev.off()







library(BayesPrism)
# sc.dat <- pmbc@assays[["RNA"]]
# 
# dim(sc.dat)#33538 17135
# sc.dat<-t(sc.dat)
# cell.type.labels <- pmbc@meta.data[["cell_type"]]
# sc.dat<-as.matrix(sc.dat)

exing<-read.table("G:\\pan_result\\gpu_result_tcga\\THCA\\scrna\\THCA_EPI_cells_E.txt")
cell_types <- rep("Malignant epithelial cells",nrow(exing))
cell_names<-exing$x

named_cell_types <- setNames(cell_types, cell_names)


pmbc<-data.combined.new

pmbc@meta.data$cell_type <- as.character(pmbc@meta.data$cell_type)


match_index <- match(rownames(pmbc@meta.data), names(named_cell_types))

pmbc@meta.data$cell_type[!is.na(match_index)] <- named_cell_types[match_index[!is.na(match_index)]]


pmbc@meta.data$cell_type <- factor(pmbc@meta.data$cell_type)


sc.dat <- pmbc[["SCT"]]@counts
dim(sc.dat)#20932 12480
sc.dat<-t(sc.dat)
cell.type.labels <- pmbc@meta.data[["cell_type"]]
sc.dat<-as.matrix(sc.dat)

nrow(sc.dat)*0.3#3744
sc.dat.filtered <- cleanup.genes (input=sc.dat,
                                  input.type="count.matrix",
                                  species="hs",
                                  gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY","hb","act"),
                                  exp.cells=nrow(sc.dat)*0.3)
dim(sc.dat.filtered)#12480  1041

save(sc.dat.filtered,file="sc.dat.filtered_70.Rdata")


rownames(sc.dat.filtered)[1:4]
sc.dat.filtered<-as.data.frame(sc.dat.filtered) 
sc.dat.filtered$cellname<-as.character(cell.type.labels)
qiuhe<-aggregate(sc.dat.filtered[,-(ncol(sc.dat.filtered))], by = sc.dat.filtered[ncol(sc.dat.filtered)], FUN = sum)
qiujunzhi<-aggregate(sc.dat.filtered[,-(ncol(sc.dat.filtered))], by = sc.dat.filtered[ncol(sc.dat.filtered)], FUN = mean)
qiuhe<-t(qiuhe)
qiujunzhi<-t(qiujunzhi)
write.table(qiuhe,"LM4_sum.txt",col.names = F,sep = "\t")
write.table(qiujunzhi,"LM4_mean.txt",col.names = F,sep = "\t")
# a<-t(sc.dat.filtered)
# max(a[3,])
# #b<-as.matrix(pmbc@active.ident)
# colnames(a)<-cell.type.labels
# class(a)
# a<-as.data.frame(a)
# b<-sapply(split.default(a, names(a)), rowMeans)
# b<-sapply(split.default(a, names(a)), rowMeans)
# dim(b)#21178     4
# head(b)
# max(b[,1])
 


#UCEC####
# THCA
library(rlang)
#packageVersion("dplyr")
#install.packages("https://www.bioconductor.org/packages/3.8/bioc/src/contrib/BiocGenerics_0.28.0.tar.gz",repos=NULL,type="source")
library(monocle)
library(monocle)
library(Seurat)
library(SeuratObject)
#load("death_result_2\\GSE135337\\HP\\GSE146771_geneplotneed.Rdata")
#load("result\\GSE146771\\HP\\GSE146771_geneplotneed.Rdata")THCA
#load("H:\\database\\SCDTD\\data\\chenlong\\Esophageal Cancer\\disease\\GSE196756_geneplotneed.Rdata")
#load("H:\\database\\SCDTD\\data\\chenlong\\Colon Cancer\\disease\\GSE140288_geneplotneed.Rdata")
#load("H:\\database\\SCDTD\\data\\chenlong\\Cervical Cancer\\disease\\GSE168652_geneplotneed.Rdata")
#load("H:\\database\\SCDTD\\data\\chenlong\\Liver Cancer\\Liver Cancer\\GSE125449_geneplotneed.Rdata")
load("H:\\database\\SCDTD\\data\\chenlong\\Thyroid Cancer\\disease\\GSE191288_geneplotneed.Rdata")
GetAssay(data.combined.new,assay = "RNA")
unique(data.combined.new@meta.data[["orig.ident"]])
pmbc<-data.combined.new


Idents(pmbc) <- "orig.ident"
ob1 <- subset(x = pmbc, downsample = 300)

expr_matrix<-as(as.matrix(ob1@assays[["RNA"]]@layers[["counts"]]),'sparseMatrix')
expr_matrix<-as(as.matrix(ob1@assays[["SCT"]]@counts),'sparseMatrix')
p_data<-ob1@meta.data
p_data$celltype<-ob1@active.ident
f_data<-data.frame(gene_short_name=row.names(ob1),row.names = row.names(ob1))
dim(expr_matrix)

pd<-new('AnnotatedDataFrame',data=p_data)
fd<-new('AnnotatedDataFrame',data=f_data)





cds<-newCellDataSet(expr_matrix,
                    phenoData=pd,
                    featureData=fd,
                    lowerDetectionLimit=0.1,
                    expressionFamily=negbinomial.size())



cds<-estimateSizeFactors(cds)
#1gc()
# memory.limit(123456)
# rm(expr_matrix)
# rm(f_data)
# rm(p_data)
cds<-estimateDispersions(cds)

cds<-detectGenes(cds,min_expr = 0.1)

print(head(fData(cds)))


gene<-read.csv("G:\\pan_result\\gpu_result_tcga\\THCA\\448featuresweightme.csv",row.names = 1)
expressed_genes<-colnames(gene)

# expressed_genes = read.table("validategene.csv",header = T,sep = ",")
# expressed_genes = expressed_genes$X0


setwd("G:\\pan_result\\gpu_result_tcga\\THCA\\scrna")
cds <- setOrderingFilter(cds, expressed_genes)
pdf("train.ordergenes.pdf")
plot_ordering_genes(cds)
dev.off()



cds <- reduceDimension(cds, max_components = 5,
                       method = 'DDRTree')##method
cds <- orderCells(cds)
setwd("G:\\pan_result\\gpu_result_tcga\\THCA\\scrna")
pdf("train.monocle.state5.pdf",width = 7,height = 7)
plot_cell_trajectory(cds, color_by = "Pseudotime",size=1,show_backbone=TRUE)
dev.off()

save(cds,file = "cds.RData")

pdf("train.monocle.celltype5.pdf",width = 7,height = 7)
plot_cell_trajectory(cds,color_by = "cell_type",size=1,show_backbone = TRUE)+
  geom_point(aes(size = Pseudotime , color = cell_type)) + 
  scale_color_manual(breaks = c("X", "Y", "Z"), values=c("darkblue", "darkred", "orange","grey","yellow","green","skyblue","pink","purple")) 
plot_cell_trajectory(cds,color_by = "cell_type",size=1,show_backbone = TRUE) + scale_color_manual( values=c("darkblue", "darkred", "orange","grey","yellow","green","skyblue","pink","purple")) + theme(legend.position = "right")
dev.off()

library(ggpubr)
library(tidyr)
packageVersion("tidyr")

#install.packages("tidyr")

#load("cds.RData")
df <- pData(cds) 

pdf("train.monocle.pseudotime.pdf",width = 7,height = 7)
ggplot(df, aes(Pseudotime, colour = cell_type, fill=cell_type)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2( )
dev.off()
write.csv(df,file = "pData.csv")



library(SingleCellExperiment)
library(Seurat)




# devtools::install_github("alanocallaghan/scater")
library(scater)
# data.combined.new@meta.data[["Pseudotime"]]<-cds@phenoData@data[["Pseudotime"]]
# sce_gene <- as.SingleCellExperiment(data.combined.new)


label<-read.table("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\THCA_cluster.txt")
gene1<-gene[rownames(label)[which(label$clusters=="1")],]
gene2<-gene[rownames(label)[which(label$clusters=="2")],]
data<-cbind(t(gene1),t(gene2))
T_p<-c()
FC_result<-c()
for(i in 1:nrow(data)){
  T_result<-wilcox.test(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]),as.numeric(data[i,1:nrow(gene1)]),paired = F,var.equal = T)
  
  p1<-T_result$p.value

  T_p<-c(T_p,p1)
  #T_value<-c(T_value,t1)
  FC<-mean(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]))/mean(as.numeric(data[i,1:nrow(gene1)]))

  FC_result<-c(FC_result,FC)
}


DEG_all<-data.frame(rownames(data),T_p,FC_result,stringsAsFactors = FALSE)

colnames(DEG_all)<-c("Gene","pval","FC") 

write.csv(DEG_all,file = "DEG_all.csv")
DEG<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC>1),]#1295
DEG_down<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC<1),]#1295


expression_matrix <- exprs(cds)


gene_expression <- expression_matrix["TYMP", ]

# gene_expression <- FetchData(cds, vars = "ZNF487", slot = "data")

cds@phenoData@data[["TYMP"]]<-gene_expression
gene_expression2<-0
gene_expression2<-gene_expression2+gene_expression

for (i in 1:nrow(DEG)) {
  gene_expression <- expression_matrix[DEG$Gene[i], ]
  gene_expression2<-gene_expression2+gene_expression
}

for (i in 31:nrow(DEG)) {
  gene_expression <- expression_matrix[DEG$Gene[i], ]
  gene_expression2<-gene_expression2+gene_expression
}
cds@phenoData@data[["upgene"]]<-gene_expression2




getwd()
pdf("monocle.pseudotime_upgene.pdf",width = 7,height = 7)

plot_cell_trajectory(cds, color_by = "upgene") +
  scale_color_gradient(low = "blue", high = "red") +   
  geom_point(aes(size = Pseudotime, color = upgene), alpha = 0.6) +   
  theme_minimal() +
  theme(legend.position = "right")   

dev.off()

#downgene####
gene_expression2<-0

for (i in 1:nrow(DEG_down)) {
  gene_expression <- expression_matrix[DEG_down$Gene[i], ]
  gene_expression2<-gene_expression2+gene_expression
}

for (i in 192:nrow(DEG_down)) {
  gene_expression <- expression_matrix[DEG_down$Gene[i], ]
  gene_expression2<-gene_expression2+gene_expression
}
cds@phenoData@data[["downgene"]]<-gene_expression2




getwd()
pdf("monocle.pseudotime_downgene.pdf",width = 7,height = 7)

plot_cell_trajectory(cds, color_by = "downgene") +
  scale_color_gradient(low = "blue", high = "red") +   
  geom_point(aes(size = Pseudotime, color = downgene), alpha = 0.6) +   
  theme_minimal() +
  theme(legend.position = "right")   

dev.off()







library(BayesPrism)
# sc.dat <- pmbc@assays[["RNA"]]
# 
# dim(sc.dat)#33538 17135
# sc.dat<-t(sc.dat)
# cell.type.labels <- pmbc@meta.data[["cell_type"]]
# sc.dat<-as.matrix(sc.dat)

exing<-read.table("G:\\pan_result\\gpu_result_tcga\\THCA\\scrna\\THCA_EPI_cells_E.txt")
cell_types <- rep("Malignant epithelial cells",nrow(exing))
cell_names<-exing$x

named_cell_types <- setNames(cell_types, cell_names)


pmbc<-data.combined.new

pmbc@meta.data$cell_type <- as.character(pmbc@meta.data$cell_type)


match_index <- match(rownames(pmbc@meta.data), names(named_cell_types))

pmbc@meta.data$cell_type[!is.na(match_index)] <- named_cell_types[match_index[!is.na(match_index)]]


pmbc@meta.data$cell_type <- factor(pmbc@meta.data$cell_type)


sc.dat <- pmbc[["SCT"]]@counts
dim(sc.dat)#20932 12480
sc.dat<-t(sc.dat)
cell.type.labels <- pmbc@meta.data[["cell_type"]]
sc.dat<-as.matrix(sc.dat)

nrow(sc.dat)*0.3#3744
sc.dat.filtered <- cleanup.genes (input=sc.dat,
                                  input.type="count.matrix",
                                  species="hs",
                                  gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY","hb","act"),
                                  exp.cells=nrow(sc.dat)*0.3)
dim(sc.dat.filtered)#12480  1041

save(sc.dat.filtered,file="sc.dat.filtered_70.Rdata")


rownames(sc.dat.filtered)[1:4]
sc.dat.filtered<-as.data.frame(sc.dat.filtered) 
sc.dat.filtered$cellname<-as.character(cell.type.labels)
qiuhe<-aggregate(sc.dat.filtered[,-(ncol(sc.dat.filtered))], by = sc.dat.filtered[ncol(sc.dat.filtered)], FUN = sum)
qiujunzhi<-aggregate(sc.dat.filtered[,-(ncol(sc.dat.filtered))], by = sc.dat.filtered[ncol(sc.dat.filtered)], FUN = mean)
qiuhe<-t(qiuhe)
qiujunzhi<-t(qiujunzhi)
write.table(qiuhe,"LM4_sum.txt",col.names = F,sep = "\t")
write.table(qiujunzhi,"LM4_mean.txt",col.names = F,sep = "\t")
 





