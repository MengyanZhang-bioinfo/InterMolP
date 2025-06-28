#BRCA####
# devtools::install_github("broadinstitute/infercnv")

# install.packages("rjags")
library(rjags)
library(infercnv)
library(remotes)
library(AnnoProbe)
library(Seurat)
library(ggplot2)
# https://git-scm.com/downloads 
# git!!!
# url='https://gitee.com/jmzeng/annoprobe.git'
# install_git(url)
# BiocManager::install("GEOquery")
# library(GEOquery)
load("G:\\pan_result\\gpu_result_tcga\\BRCA\\scRNA_new\\harmony_0.4_20\\shoudong\\GSE161529.seu.obj.Rdata")
setwd("G:\\pan_result\\gpu_result_tcga\\BRCA\\scRNA_new")
options("Seurat.object.assay.version" = "v3")
#1####
#
seu.obj$cell_type <- seu.obj@active.ident
# seu.obj@assays[["RNA"]]@layers[["counts"]]@Dimnames[[1]]<-rownames(seu.obj)
# seu.obj@assays[["RNA"]]@layers[["counts"]]@Dimnames[[2]]<-colnames(seu.obj)
epi.cells <- subset(seu.obj, idents=c('Epithelium'))
# epi.cells  <- row.names(seu.obj@meta.data)[which(seu.obj@meta.data[["cell_type"]]=='Epithelial_cells')]
epi_sce <- subset(seu.obj, cells = epi.cells@assays[["RNA"]]@counts@Dimnames[[2]])
epi_sce@meta.data
length(epi.cells)
#Raw Counts Matrix for Genes x Cells
epiMat=as.data.frame(GetAssayData(subset(seu.obj, cells = epi.cells@assays[["RNA"]]@counts@Dimnames[[2]])))
# CBtab
groupinfo=data.frame(v1=colnames(epiMat),
                     v2=epi_sce@meta.data$cell_type)
#
geneInfor=annoGene(rownames(epiMat),"SYMBOL",'human')
colnames(geneInfor)
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)

## 
# 
dat=epiMat[rownames(epiMat) %in% geneInfor[,1],]
dat=dat[match( geneInfor[,1], rownames(dat) ),] 
dim(dat)
expFile='expFile.txt'
write.table(dat,file = expFile,sep = '\t',quote = F)
groupFiles='groupFiles.txt'
head(groupinfo)
write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
head(geneInfor)
geneFile='geneFile.txt'
write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)
library(Matrix)
dat_matrix <- as.matrix(dat)

# dat_sparse <- as(dat_matrix, "dgCMatrix")
#2inferCNV####
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=dat_matrix,
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= geneFile,
                                    ref_group_names=NULL)  ## 
# 
out_dir <- "path_to_save_inferCNV_result"
# infercnv_obj <- NormalizeData(infercnv_obj) # 
# infercnv_obj <- ScaleData(infercnv_obj) # 
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=out_dir, 
                             cluster_by_groups=T, 
                             denoise=TRUE,
                             HMM=T,
                             num_threads=8)
#3####
# CNV
cnv_regions <- read.table("G:\\pan_result\\gpu_result_tcga\\BRCA\\scRNA_new\\path_to_save_inferCNV_result\\17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_regions.dat", header = TRUE)
cnv_genes <- read.table("path_to_save_inferCNV_result/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_genes.dat", header = TRUE)

# 
cell_groupings <- read.table("path_to_save_inferCNV_result/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings", header = TRUE)

# CNV
malignant_cells <- cell_groupings[cell_groupings$cell_group_name %in% cnv_regions$cell_group, ]

# 
malignant_cells_need <- malignant_cells$cell
malignant_cells_need
write.table(malignant_cells_need,"BRCA_EPI_cells_E.txt")

#COAD####\
# devtools::install_github("broadinstitute/infercnv")

# install.packages("rjags")
library(rjags)
library(infercnv)
library(remotes)
library(AnnoProbe)
library(Seurat)
library(ggplot2)
# https://git-scm.com/downloads 
# git!!!
# url='https://gitee.com/jmzeng/annoprobe.git'
# install_git(url)
# BiocManager::install("GEOquery")
# library(GEOquery)
load("G:\\pan_result\\gpu_result_tcga\\COAD\\scRNA_new\\harmony_0.4_20\\shoudong\\GSE188711.seu.obj.Rdata")
setwd("G:\\pan_result\\gpu_result_tcga\\COAD\\scRNA_new")
options("Seurat.object.assay.version" = "v3")
#1####
#
seu.obj$cell_type <- seu.obj@active.ident
seu.obj@assays[["RNA"]]@layers[["counts"]]@Dimnames[[1]]<-rownames(seu.obj)
seu.obj@assays[["RNA"]]@layers[["counts"]]@Dimnames[[2]]<-colnames(seu.obj)
epi.cells <- subset(seu.obj, idents=c('Epithelium'))
epi.cells@assays[["RNA"]]@layers[["counts"]]@Dimnames[[1]]<-rownames(epi.cells)
epi.cells@assays[["RNA"]]@layers[["counts"]]@Dimnames[[2]]<-colnames(epi.cells)
# epi.cells  <- row.names(seu.obj@meta.data)[which(seu.obj@meta.data[["cell_type"]]=='Epithelial_cells')]
epi_sce <- subset(seu.obj, cells = epi.cells@assays[["RNA"]]@layers[["counts"]]@Dimnames[[2]])
epi_sce@meta.data
length(epi.cells)
#Raw Counts Matrix for Genes x Cells
epiMat=as.data.frame(GetAssayData(subset(seu.obj, cells=epi.cells@assays[["RNA"]]@layers[["counts"]]@Dimnames[[2]])))
# CBtab
groupinfo=data.frame(v1=colnames(epiMat),
                     v2=epi_sce@meta.data$cell_type)
#
geneInfor=annoGene(rownames(epiMat),"SYMBOL",'human')
colnames(geneInfor)
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)

## 
# 
dat=epiMat[rownames(epiMat) %in% geneInfor[,1],]
dat=dat[match( geneInfor[,1], rownames(dat) ),] 
dim(dat)
expFile='expFile.txt'
write.table(dat,file = expFile,sep = '\t',quote = F)
groupFiles='groupFiles.txt'
head(groupinfo)
write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
head(geneInfor)
geneFile='geneFile.txt'
write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)
library(Matrix)
dat_matrix <- as.matrix(dat)

#dat_sparse <- as(dat_matrix, "dgCMatrix")
#2inferCNV####
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=dat_matrix,
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= geneFile,
                                    ref_group_names=NULL)  ## 
# 
out_dir <- "path_to_save_inferCNV_result"
# infercnv_obj <- NormalizeData(infercnv_obj) # 
# infercnv_obj <- ScaleData(infercnv_obj) # 
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=out_dir, 
                             cluster_by_groups=T, 
                             denoise=TRUE,
                             HMM=T,
                             num_threads=8)
#3####
# CNV
cnv_regions <- read.table("G:\\pan_result\\gpu_result_tcga\\COAD\\scRNA_new\\path_to_save_inferCNV_result\\17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_regions.dat", header = TRUE)
cnv_genes <- read.table("path_to_save_inferCNV_result/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_genes.dat", header = TRUE)

# 
cell_groupings <- read.table("path_to_save_inferCNV_result/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings", header = TRUE)

# CNV
malignant_cells <- cell_groupings[cell_groupings$cell_group_name %in% cnv_regions$cell_group, ]

# 
malignant_cells_need <- malignant_cells$cell
malignant_cells_need
write.table(malignant_cells_need,"COAD_EPI_cells_E.txt")

#ESCA####
# devtools::install_github("broadinstitute/infercnv")

# install.packages("rjags")
library(rjags)
library(infercnv)
library(remotes)
library(AnnoProbe)
library(Seurat)
library(ggplot2)
# https://git-scm.com/downloads 
# git!!!
# url='https://gitee.com/jmzeng/annoprobe.git'
# install_git(url)
# BiocManager::install("GEOquery")
# library(GEOquery)
load("G:\\pan_result\\gpu_result_tcga\\ESCA\\scRNA_new\\harmony_0.4_20\\shoudong\\GSE196756.seu.obj.Rdata")
setwd("G:\\pan_result\\gpu_result_tcga\\ESCA\\scRNA_new")
options("Seurat.object.assay.version" = "v3")
#1####
#
seu.obj$cell_type <- seu.obj@active.ident
# seu.obj@assays[["RNA"]]@layers[["counts"]]@Dimnames[[1]]<-rownames(seu.obj)
# seu.obj@assays[["RNA"]]@layers[["counts"]]@Dimnames[[2]]<-colnames(seu.obj)
epi.cells <- subset(seu.obj, idents=c('Epithelium'))
# epi.cells  <- row.names(seu.obj@meta.data)[which(seu.obj@meta.data[["cell_type"]]=='Epithelial_cells')]
epi_sce <- subset(seu.obj, cells = epi.cells@assays[["SCT"]]@counts@Dimnames[[2]])
epi_sce@meta.data
length(epi.cells)
#Raw Counts Matrix for Genes x Cells
epiMat=as.data.frame(GetAssayData(subset(seu.obj, cells = epi.cells@assays[["SCT"]]@counts@Dimnames[[2]])))
# CBtab
groupinfo=data.frame(v1=colnames(epiMat),
                     v2=epi_sce@meta.data$cell_type)
#
geneInfor=annoGene(rownames(epiMat),"SYMBOL",'human')
colnames(geneInfor)
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)

## 
# 
dat=epiMat[rownames(epiMat) %in% geneInfor[,1],]
dat=dat[match( geneInfor[,1], rownames(dat) ),] 
dim(dat)
expFile='expFile.txt'
write.table(dat,file = expFile,sep = '\t',quote = F)
groupFiles='groupFiles.txt'
head(groupinfo)
write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
head(geneInfor)
geneFile='geneFile.txt'
write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)
library(Matrix)
dat_matrix <- as.matrix(dat)

dat_sparse <- as(dat_matrix, "dgCMatrix")
#2inferCNV####
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=dat_matrix,
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= geneFile,
                                    ref_group_names=NULL)  ## 
# 
out_dir <- "path_to_save_inferCNV_result"
# infercnv_obj <- NormalizeData(infercnv_obj) # 
# infercnv_obj <- ScaleData(infercnv_obj) # 
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=out_dir, 
                             cluster_by_groups=T, 
                             denoise=TRUE,
                             HMM=T,
                             num_threads=8)
#3####
# CNV
cnv_regions <- read.table("G:\\pan_result\\gpu_result_tcga\\ESCA\\scRNA_new\\path_to_save_inferCNV_result\\17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_regions.dat", header = TRUE)
cnv_genes <- read.table("path_to_save_inferCNV_result/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_genes.dat", header = TRUE)

# 
cell_groupings <- read.table("path_to_save_inferCNV_result/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings", header = TRUE)

# CNV
malignant_cells <- cell_groupings[cell_groupings$cell_group_name %in% cnv_regions$cell_group, ]

# 
malignant_cells_need <- malignant_cells$cell
malignant_cells_need
write.table(malignant_cells_need,"ESCA_EPI_cells_E.txt")

#HNSC####
# devtools::install_github("broadinstitute/infercnv")

# install.packages("rjags")
library(rjags)
library(infercnv)
library(remotes)
library(AnnoProbe)
library(Seurat)
library(ggplot2)
# https://git-scm.com/downloads 
# git!!!
# url='https://gitee.com/jmzeng/annoprobe.git'
# install_git(url)
# BiocManager::install("GEOquery")
# library(GEOquery)
load("H:\\\\database\\HNSC\\geo\\GSE243359_seruat.Rdata")
data.combined.new@assays[["RNA"]]@layers[["counts"]]@Dimnames[[1]]<-rownames(data.combined.new)
data.combined.new@assays[["RNA"]]@layers[["counts"]]@Dimnames[[2]]<-colnames(data.combined.new)
setwd("G:\\pan_result\\gpu_result_tcga\\HNSC\\scrna")
options("Seurat.object.assay.version" = "v3")
#1####
#
epi.cells <- subset(data.combined.new, idents=c('Epithelial_cells'))
epi.cells@assays[["RNA"]]@layers[["counts"]]@Dimnames[[1]]<-rownames(epi.cells)
epi.cells@assays[["RNA"]]@layers[["counts"]]@Dimnames[[2]]<-colnames(epi.cells)
row.names(epi.cells )
# epi.cells  <- row.names(data.combined.new@meta.data)[which(data.combined.new@meta.data[["cell_type"]]=='Epithelial_cells')]
epi_sce <- subset(data.combined.new, cells = epi.cells@assays[["RNA"]]@layers[["counts"]]@Dimnames[[2]])
epi_sce@meta.data
length(epi.cells)
#Raw Counts Matrix for Genes x Cells
epiMat=as.data.frame(GetAssayData(subset(data.combined.new, cells=epi.cells@assays[["RNA"]]@layers[["counts"]]@Dimnames[[2]])))
# CBtab
groupinfo=data.frame(v1=colnames(epiMat),
                     v2=epi_sce@meta.data$cell_type)
#
geneInfor=annoGene(rownames(epiMat),"SYMBOL",'human')
colnames(geneInfor)
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)

## 
# 
dat=epiMat[rownames(epiMat) %in% geneInfor[,1],]
dat=dat[match( geneInfor[,1], rownames(dat) ),] 
dim(dat)
expFile='expFile.txt'
write.table(dat,file = expFile,sep = '\t',quote = F)
groupFiles='groupFiles.txt'
head(groupinfo)
write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
head(geneInfor)
geneFile='geneFile.txt'
write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)
library(Matrix)
dat_matrix <- as.matrix(dat)

dat_sparse <- as(dat_matrix, "dgCMatrix")
#2inferCNV####
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=dat_matrix,
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= geneFile,
                                    ref_group_names=NULL)  ## 
# 
out_dir <- "path_to_save_inferCNV_result"
# infercnv_obj <- NormalizeData(infercnv_obj) # 
# infercnv_obj <- ScaleData(infercnv_obj) # 
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=out_dir, 
                             cluster_by_groups=T, 
                             denoise=TRUE,
                             HMM=T,
                             num_threads=8)
#3####
# CNV
cnv_regions <- read.table("G:\\pan_result\\gpu_result_tcga\\HNSC\\scrna\\path_to_save_inferCNV_result\\17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_regions.dat", header = TRUE)
cnv_genes <- read.table("path_to_save_inferCNV_result/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_genes.dat", header = TRUE)

# 
cell_groupings <- read.table("path_to_save_inferCNV_result/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings", header = TRUE)

# CNV
malignant_cells <- cell_groupings[cell_groupings$cell_group_name %in% cnv_regions$cell_group, ]

# 
malignant_cells_need <- malignant_cells$cell
malignant_cells_need
write.table(malignant_cells_need,"HNSC_EPI_cells_E.txt")

#KIRC####
# devtools::install_github("broadinstitute/infercnv")

# install.packages("rjags")
library(rjags)
library(infercnv)
library(remotes)
library(AnnoProbe)
library(Seurat)
library(ggplot2)
# https://git-scm.com/downloads 
# git!!!
# url='https://gitee.com/jmzeng/annoprobe.git'
# install_git(url)
# BiocManager::install("GEOquery")
# library(GEOquery)
load("G:\\pan_result\\gpu_result_tcga\\KIRC\\scRNA_new\\harmony_0.4_20\\shoudong\\GSE159115.seu.obj.Rdata")
setwd("G:\\pan_result\\gpu_result_tcga\\KIRC\\scRNA_new")
options("Seurat.object.assay.version" = "v3")
#1####
#
seu.obj$cell_type <- seu.obj@active.ident
# seu.obj@assays[["RNA"]]@layers[["counts"]]@Dimnames[[1]]<-rownames(seu.obj)
# seu.obj@assays[["RNA"]]@layers[["counts"]]@Dimnames[[2]]<-colnames(seu.obj)
epi.cells <- subset(seu.obj, idents=c('Epithelium'))
# epi.cells  <- row.names(seu.obj@meta.data)[which(seu.obj@meta.data[["cell_type"]]=='Epithelial_cells')]
epi_sce <- subset(seu.obj, cells = epi.cells@assays[["SCT"]]@counts@Dimnames[[2]])
epi_sce@meta.data
length(epi.cells)
#Raw Counts Matrix for Genes x Cells
epiMat=as.data.frame(GetAssayData(subset(seu.obj, cells = epi.cells@assays[["SCT"]]@counts@Dimnames[[2]])))
# CBtab
groupinfo=data.frame(v1=colnames(epiMat),
                     v2=epi_sce@meta.data$cell_type)
#
geneInfor=annoGene(rownames(epiMat),"SYMBOL",'human')
colnames(geneInfor)
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)

## 
# 
dat=epiMat[rownames(epiMat) %in% geneInfor[,1],]
dat=dat[match( geneInfor[,1], rownames(dat) ),] 
dim(dat)
expFile='expFile.txt'
write.table(dat,file = expFile,sep = '\t',quote = F)
groupFiles='groupFiles.txt'
head(groupinfo)
write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
head(geneInfor)
geneFile='geneFile.txt'
write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)
library(Matrix)
dat_matrix <- as.matrix(dat)

dat_sparse <- as(dat_matrix, "dgCMatrix")
#2inferCNV####
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=dat_matrix,
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= geneFile,
                                    ref_group_names=NULL)  ## 
# 
out_dir <- "path_to_save_inferCNV_result"
# infercnv_obj <- NormalizeData(infercnv_obj) # 
# infercnv_obj <- ScaleData(infercnv_obj) # 
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=out_dir, 
                             cluster_by_groups=T, 
                             denoise=TRUE,
                             HMM=T,
                             num_threads=8)
#3####
# CNV
cnv_regions <- read.table("G:\\pan_result\\gpu_result_tcga\\KIRC\\scRNA_new\\path_to_save_inferCNV_result\\17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_regions.dat", header = TRUE)
cnv_genes <- read.table("path_to_save_inferCNV_result/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_genes.dat", header = TRUE)

# 
cell_groupings <- read.table("path_to_save_inferCNV_result/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings", header = TRUE)

# CNV
malignant_cells <- cell_groupings[cell_groupings$cell_group_name %in% cnv_regions$cell_group, ]

# 
malignant_cells_need <- malignant_cells$cell
malignant_cells_need
write.table(malignant_cells_need,"KIRC_EPI_cells_E.txt")


#KIRP####
# devtools::install_github("broadinstitute/infercnv")

# install.packages("rjags")
library(rjags)
library(infercnv)
library(remotes)
library(AnnoProbe)
library(Seurat)
library(ggplot2)
# https://git-scm.com/downloads 
# git!!!
# url='https://gitee.com/jmzeng/annoprobe.git'
# install_git(url)
# BiocManager::install("GEOquery")
# library(GEOquery)
load("G:\\pan_result\\gpu_result_tcga\\KIRP\\scRNA_new\\harmony_0.4_20\\shoudong\\GSE152938.seu.obj.Rdata")
setwd("G:\\pan_result\\gpu_result_tcga\\KIRP\\scRNA_new")
options("Seurat.object.assay.version" = "v3")
#1####
#
seu.obj$cell_type <- seu.obj@active.ident
# seu.obj@assays[["RNA"]]@layers[["counts"]]@Dimnames[[1]]<-rownames(seu.obj)
# seu.obj@assays[["RNA"]]@layers[["counts"]]@Dimnames[[2]]<-colnames(seu.obj)
epi.cells <- subset(seu.obj, idents=c('Epithelium'))
# epi.cells  <- row.names(seu.obj@meta.data)[which(seu.obj@meta.data[["cell_type"]]=='Epithelial_cells')]
epi_sce <- subset(seu.obj, cells = epi.cells@assays[["RNA"]]@counts@Dimnames[[2]])
epi_sce@meta.data
length(epi.cells)
#Raw Counts Matrix for Genes x Cells
epiMat=as.data.frame(GetAssayData(subset(seu.obj, cells = epi.cells@assays[["RNA"]]@counts@Dimnames[[2]])))
# CBtab
groupinfo=data.frame(v1=colnames(epiMat),
                     v2=epi_sce@meta.data$cell_type)
#
geneInfor=annoGene(rownames(epiMat),"SYMBOL",'human')
colnames(geneInfor)
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)

## 
# 
dat=epiMat[rownames(epiMat) %in% geneInfor[,1],]
dat=dat[match( geneInfor[,1], rownames(dat) ),] 
dim(dat)
expFile='expFile.txt'
write.table(dat,file = expFile,sep = '\t',quote = F)
groupFiles='groupFiles.txt'
head(groupinfo)
write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
head(geneInfor)
geneFile='geneFile.txt'
write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)
library(Matrix)
dat_matrix <- as.matrix(dat)

dat_sparse <- as(dat_matrix, "dgCMatrix")
#2inferCNV####
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=dat_matrix,
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= geneFile,
                                    ref_group_names=NULL)  ## 
# 
out_dir <- "path_to_save_inferCNV_result"
# infercnv_obj <- NormalizeData(infercnv_obj) # 
# infercnv_obj <- ScaleData(infercnv_obj) # 
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=out_dir, 
                             cluster_by_groups=T, 
                             denoise=TRUE,
                             HMM=T,
                             num_threads=8)
#3####
# CNV
cnv_regions <- read.table("G:\\pan_result\\gpu_result_tcga\\KIRP\\scRNA_new\\path_to_save_inferCNV_result\\17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_regions.dat", header = TRUE)
cnv_genes <- read.table("path_to_save_inferCNV_result/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_genes.dat", header = TRUE)

# 
cell_groupings <- read.table("path_to_save_inferCNV_result/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings", header = TRUE)

# CNV
malignant_cells <- cell_groupings[cell_groupings$cell_group_name %in% cnv_regions$cell_group, ]

# 
malignant_cells_need <- malignant_cells$cell
malignant_cells_need
write.table(malignant_cells_need,"KIRP_EPI_cells_E.txt")


#LIHC####
# devtools::install_github("broadinstitute/infercnv")

# install.packages("rjags")
library(rjags)
library(infercnv)
library(remotes)
library(AnnoProbe)
library(Seurat)
library(ggplot2)
# https://git-scm.com/downloads 
# git!!!
# url='https://gitee.com/jmzeng/annoprobe.git'
# install_git(url)
# BiocManager::install("GEOquery")
# library(GEOquery)
load("G:\\pan_result\\gpu_result_tcga\\LIHC\\scRNA_new\\harmony_0.4_20\\shoudong\\GSE112271.seu.obj.Rdata")
setwd("G:\\pan_result\\gpu_result_tcga\\LIHC\\scRNA_new")
options("Seurat.object.assay.version" = "v3")
#1####
#
seu.obj$cell_type <- seu.obj@active.ident
seu.obj@assays[["RNA"]]@layers[["counts"]]@Dimnames[[1]]<-rownames(seu.obj)
seu.obj@assays[["RNA"]]@layers[["counts"]]@Dimnames[[2]]<-colnames(seu.obj)
epi.cells <- subset(seu.obj, idents=c('Epithelium'))
epi.cells@assays[["RNA"]]@layers[["counts"]]@Dimnames[[1]]<-rownames(epi.cells)
epi.cells@assays[["RNA"]]@layers[["counts"]]@Dimnames[[2]]<-colnames(epi.cells)
# epi.cells  <- row.names(seu.obj@meta.data)[which(seu.obj@meta.data[["cell_type"]]=='Epithelial_cells')]
epi_sce <- subset(seu.obj, cells = epi.cells@assays[["RNA"]]@layers[["counts"]]@Dimnames[[2]])
epi_sce@meta.data
length(epi.cells)
#Raw Counts Matrix for Genes x Cells
epiMat=as.data.frame(GetAssayData(subset(seu.obj, cells=epi.cells@assays[["RNA"]]@layers[["counts"]]@Dimnames[[2]])))
# CBtab
groupinfo=data.frame(v1=colnames(epiMat),
                     v2=epi_sce@meta.data$cell_type)
#
geneInfor=annoGene(rownames(epiMat),"SYMBOL",'human')
colnames(geneInfor)
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)

## 
# 
dat=epiMat[rownames(epiMat) %in% geneInfor[,1],]
dat=dat[match( geneInfor[,1], rownames(dat) ),] 
dim(dat)
expFile='expFile.txt'
write.table(dat,file = expFile,sep = '\t',quote = F)
groupFiles='groupFiles.txt'
head(groupinfo)
write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
head(geneInfor)
geneFile='geneFile.txt'
write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)
library(Matrix)
dat_matrix <- as.matrix(dat)

# dat_sparse <- as(dat_matrix, "dgCMatrix")
rm(seu.obj)
#2inferCNV####
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=dat_matrix,
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= geneFile,
                                    ref_group_names=NULL)  ## 
# 
out_dir <- "path_to_save_inferCNV_result"
# infercnv_obj <- NormalizeData(infercnv_obj) # 
# infercnv_obj <- ScaleData(infercnv_obj) # 
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=out_dir, 
                             cluster_by_groups=T, 
                             denoise=TRUE,
                             HMM=T,
                             num_threads=8)
#3####
# CNV
cnv_regions <- read.table("G:\\pan_result\\gpu_result_tcga\\LIHC\\scRNA_new\\path_to_save_inferCNV_result\\17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_regions.dat", header = TRUE)
cnv_genes <- read.table("path_to_save_inferCNV_result/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_genes.dat", header = TRUE)

# 
cell_groupings <- read.table("path_to_save_inferCNV_result/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings", header = TRUE)

# CNV
malignant_cells <- cell_groupings[cell_groupings$cell_group_name %in% cnv_regions$cell_group, ]

# 
malignant_cells_need <- malignant_cells$cell
malignant_cells_need
write.table(malignant_cells_need,"LIHC_EPI_cells_E.txt")


#LUNG####
# devtools::install_github("broadinstitute/infercnv")

# install.packages("rjags")
library(rjags)
library(infercnv)
library(remotes)
library(AnnoProbe)
library(Seurat)
library(ggplot2)
# https://git-scm.com/downloads 
# git!!!
# url='https://gitee.com/jmzeng/annoprobe.git'
# install_git(url)
# BiocManager::install("GEOquery")
# library(GEOquery)
load("G:\\pan_result\\gpu_result_tcga\\LUNG\\scRNA_new\\harmony_0.4_20\\shoudong\\GSE131907.seu.obj.Rdata")
setwd("G:\\pan_result\\gpu_result_tcga\\LUNG\\scRNA_new")
options("Seurat.object.assay.version" = "v3")
#1####
#
seu.obj$cell_type <- seu.obj@active.ident
seu.obj@assays[["RNA"]]@layers[["counts"]]@Dimnames[[1]]<-rownames(seu.obj)
seu.obj@assays[["RNA"]]@layers[["counts"]]@Dimnames[[2]]<-colnames(seu.obj)
epi.cells <- subset(seu.obj, idents=c('Epithelium'))
epi.cells@assays[["RNA"]]@layers[["counts"]]@Dimnames[[1]]<-rownames(epi.cells)
epi.cells@assays[["RNA"]]@layers[["counts"]]@Dimnames[[2]]<-colnames(epi.cells)
# epi.cells  <- row.names(seu.obj@meta.data)[which(seu.obj@meta.data[["cell_type"]]=='Epithelial_cells')]
epi_sce <- subset(seu.obj, cells = epi.cells@assays[["RNA"]]@layers[["counts"]]@Dimnames[[2]])
epi_sce@meta.data
length(epi.cells)
#Raw Counts Matrix for Genes x Cells
epiMat=as.data.frame(GetAssayData(subset(seu.obj, cells=epi.cells@assays[["RNA"]]@layers[["counts"]]@Dimnames[[2]])))
rm(seu.obj)

# CBtab
groupinfo=data.frame(v1=colnames(epiMat),
                     v2=epi_sce@meta.data$cell_type)
#
geneInfor=annoGene(rownames(epiMat),"SYMBOL",'human')
colnames(geneInfor)
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)

## 
# 
dat=epiMat[rownames(epiMat) %in% geneInfor[,1],]
dat=dat[match( geneInfor[,1], rownames(dat) ),] 
dim(dat)
expFile='expFile.txt'
write.table(dat,file = expFile,sep = '\t',quote = F)
groupFiles='groupFiles.txt'
head(groupinfo)
write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
head(geneInfor)
geneFile='geneFile.txt'
write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)
library(Matrix)
dat_matrix <- as.matrix(dat)

dat_sparse <- as(dat_matrix, "dgCMatrix")
#2inferCNV####
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=dat_matrix,
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= geneFile,
                                    ref_group_names=NULL)  ## 
# 
out_dir <- "path_to_save_inferCNV_result"
# infercnv_obj <- NormalizeData(infercnv_obj) # 
# infercnv_obj <- ScaleData(infercnv_obj) # 
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=out_dir, 
                             cluster_by_groups=T, 
                             denoise=TRUE,
                             HMM=T,
                             num_threads=8)
#3####
# CNV
cnv_regions <- read.table("G:\\pan_result\\gpu_result_tcga\\LUNG\\scRNA_new\\path_to_save_inferCNV_result\\17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_regions.dat", header = TRUE)
cnv_genes <- read.table("path_to_save_inferCNV_result/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_genes.dat", header = TRUE)

# 
cell_groupings <- read.table("path_to_save_inferCNV_result/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings", header = TRUE)

# CNV
malignant_cells <- cell_groupings[cell_groupings$cell_group_name %in% cnv_regions$cell_group, ]

# 
malignant_cells_need <- malignant_cells$cell
malignant_cells_need
write.table(malignant_cells_need,"LUNG_EPI_cells_E.txt")


#PAAD####

# install.packages("rjags")
library(rjags)
library(infercnv)
library(remotes)
library(AnnoProbe)
library(Seurat)
library(ggplot2)
# https://git-scm.com/downloads 
# git!!!
# url='https://gitee.com/jmzeng/annoprobe.git'
# install_git(url)
# BiocManager::install("GEOquery")
# library(GEOquery)
load("G:\\pan_result\\gpu_result_tcga\\PAAD\\scRNA_new\\harmony_0.4_20\\shoudong\\GSE155698.seu.obj.Rdata")
setwd("G:\\pan_result\\gpu_result_tcga\\PAAD\\scRNA_new")
options("Seurat.object.assay.version" = "v3")
#1####
#
seu.obj$cell_type <- seu.obj@active.ident
seu.obj@assays[["RNA"]]@layers[["counts"]]@Dimnames[[1]]<-rownames(seu.obj)
seu.obj@assays[["RNA"]]@layers[["counts"]]@Dimnames[[2]]<-colnames(seu.obj)
epi.cells <- subset(seu.obj, idents=c('Epithelium'))
epi.cells@assays[["RNA"]]@layers[["counts"]]@Dimnames[[1]]<-rownames(epi.cells)
epi.cells@assays[["RNA"]]@layers[["counts"]]@Dimnames[[2]]<-colnames(epi.cells)
# epi.cells  <- row.names(seu.obj@meta.data)[which(seu.obj@meta.data[["cell_type"]]=='Epithelial_cells')]
epi_sce <- subset(seu.obj, cells = epi.cells@assays[["RNA"]]@layers[["counts"]]@Dimnames[[2]])
epi_sce@meta.data
length(epi.cells)
#Raw Counts Matrix for Genes x Cells
epiMat=as.data.frame(GetAssayData(subset(seu.obj, cells=epi.cells@assays[["RNA"]]@layers[["counts"]]@Dimnames[[2]])))
# CBtab
groupinfo=data.frame(v1=colnames(epiMat),
                     v2=epi_sce@meta.data$cell_type)
#
geneInfor=annoGene(rownames(epiMat),"SYMBOL",'human')
colnames(geneInfor)
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)

## 
# 
dat=epiMat[rownames(epiMat) %in% geneInfor[,1],]
dat=dat[match( geneInfor[,1], rownames(dat) ),] 
dim(dat)
expFile='expFile.txt'
write.table(dat,file = expFile,sep = '\t',quote = F)
groupFiles='groupFiles.txt'
head(groupinfo)
write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
head(geneInfor)
geneFile='geneFile.txt'
write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)
library(Matrix)
dat_matrix <- as.matrix(dat)

#dat_sparse <- as(dat_matrix, "dgCMatrix")
#2inferCNV####
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=dat_matrix,
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= geneFile,
                                    ref_group_names=NULL)  ## 
# 
out_dir <- "path_to_save_inferCNV_result"
# infercnv_obj <- NormalizeData(infercnv_obj) # 
# infercnv_obj <- ScaleData(infercnv_obj) # 
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=out_dir, 
                             cluster_by_groups=T, 
                             denoise=TRUE,
                             HMM=T,
                             num_threads=8)
#3####
# CNV
cnv_regions <- read.table("G:\\pan_result\\gpu_result_tcga\\PAAD\\scRNA_new\\path_to_save_inferCNV_result\\17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_regions.dat", header = TRUE)
cnv_genes <- read.table("path_to_save_inferCNV_result/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_genes.dat", header = TRUE)

# 
cell_groupings <- read.table("path_to_save_inferCNV_result/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings", header = TRUE)

# CNV
malignant_cells <- cell_groupings[cell_groupings$cell_group_name %in% cnv_regions$cell_group, ]

# 
malignant_cells_need <- malignant_cells$cell
malignant_cells_need
write.table(malignant_cells_need,"PAAD_EPI_cells_E.txt")


#THCA####
# devtools::install_github("broadinstitute/infercnv")

# install.packages("rjags")
library(rjags)
library(infercnv)
library(remotes)
library(AnnoProbe)
library(Seurat)
library(ggplot2)
# https://git-scm.com/downloads 
# git!!!
# url='https://gitee.com/jmzeng/annoprobe.git'
# install_git(url)
# BiocManager::install("GEOquery")
# library(GEOquery)
load("G:\\pan_result\\gpu_result_tcga\\THCA\\scRNA_new\\harmony_0.4_20\\shoudong\\GSE191288.seu.obj.Rdata")
setwd("G:\\pan_result\\gpu_result_tcga\\THCA\\scRNA_new")
options("Seurat.object.assay.version" = "v3")
#1####
#
seu.obj$cell_type <- seu.obj@active.ident
# seu.obj@assays[["RNA"]]@layers[["counts"]]@Dimnames[[1]]<-rownames(seu.obj)
# seu.obj@assays[["RNA"]]@layers[["counts"]]@Dimnames[[2]]<-colnames(seu.obj)
epi.cells <- subset(seu.obj, idents=c('Epithelium'))
# epi.cells  <- row.names(seu.obj@meta.data)[which(seu.obj@meta.data[["cell_type"]]=='Epithelial_cells')]
epi_sce <- subset(seu.obj, cells = epi.cells@assays[["SCT"]]@counts@Dimnames[[2]])
epi_sce@meta.data
length(epi.cells)
#Raw Counts Matrix for Genes x Cells
epiMat=as.data.frame(GetAssayData(subset(seu.obj, cells = epi.cells@assays[["SCT"]]@counts@Dimnames[[2]])))
# CBtab
groupinfo=data.frame(v1=colnames(epiMat),
                     v2=epi_sce@meta.data$cell_type)
#
geneInfor=annoGene(rownames(epiMat),"SYMBOL",'human')
colnames(geneInfor)
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)

## 
# 
dat=epiMat[rownames(epiMat) %in% geneInfor[,1],]
dat=dat[match( geneInfor[,1], rownames(dat) ),] 
dim(dat)
expFile='expFile.txt'
write.table(dat,file = expFile,sep = '\t',quote = F)
groupFiles='groupFiles.txt'
head(groupinfo)
write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
head(geneInfor)
geneFile='geneFile.txt'
write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)
library(Matrix)
dat_matrix <- as.matrix(dat)

dat_sparse <- as(dat_matrix, "dgCMatrix")
#2inferCNV####
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=dat_matrix,
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= geneFile,
                                    ref_group_names=NULL)  ## 
# 
out_dir <- "path_to_save_inferCNV_result"
# infercnv_obj <- NormalizeData(infercnv_obj) # 
# infercnv_obj <- ScaleData(infercnv_obj) # 
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=out_dir, 
                             cluster_by_groups=T, 
                             denoise=TRUE,
                             HMM=T,
                             num_threads=8)
#3####
# CNV
cnv_regions <- read.table("G:\\pan_result\\gpu_result_tcga\\THCA\\scRNA_new\\path_to_save_inferCNV_result\\17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_regions.dat", header = TRUE)
cnv_genes <- read.table("path_to_save_inferCNV_result/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_genes.dat", header = TRUE)

# 
cell_groupings <- read.table("path_to_save_inferCNV_result/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings", header = TRUE)

# CNV
malignant_cells <- cell_groupings[cell_groupings$cell_group_name %in% cnv_regions$cell_group, ]

# 
malignant_cells_need <- malignant_cells$cell
malignant_cells_need
write.table(malignant_cells_need,"THCA_EPI_cells_E.txt")


#UCEC####
# devtools::install_github("broadinstitute/infercnv")

# install.packages("rjags")
library(rjags)
library(infercnv)
library(remotes)
library(AnnoProbe)
library(Seurat)
library(ggplot2)
# https://git-scm.com/downloads 
# git!!!
# url='https://gitee.com/jmzeng/annoprobe.git'
# install_git(url)
# BiocManager::install("GEOquery")
# library(GEOquery)
load("G:\\pan_result\\gpu_result_tcga\\UCEC\\scRNA_new\\harmony_0.4_20\\shoudong\\GSE173682.seu.obj.Rdata")
setwd("G:\\pan_result\\gpu_result_tcga\\UCEC\\scRNA_new")
options("Seurat.object.assay.version" = "v3")
#1####
#
epi.cells <- subset(seu.obj, idents=c('Epithelium'))
# epi.cells  <- row.names(seu.obj@meta.data)[which(seu.obj@meta.data[["cell_type"]]=='Epithelium')]
epi_sce <- subset(seu.obj, cells = epi.cells@assays[["SCT"]]@counts@Dimnames[[2]])
epi_sce@meta.data
length(epi.cells)
#Raw Counts Matrix for Genes x Cells
epiMat=as.data.frame(GetAssayData(subset(seu.obj, cells=epi.cells@assays[["SCT"]]@counts@Dimnames[[2]])))
rm(seu.obj)
# CBtab
groupinfo=data.frame(v1=colnames(epiMat),
                     v2=epi_sce@meta.data$cell_type)
#
geneInfor=annoGene(rownames(epiMat),"SYMBOL",'human')
colnames(geneInfor)
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)

## 
# 
dat=epiMat[rownames(epiMat) %in% geneInfor[,1],]
dat=dat[match( geneInfor[,1], rownames(dat) ),] 
dim(dat)
expFile='expFile.txt'
write.table(dat,file = expFile,sep = '\t',quote = F)
groupFiles='groupFiles.txt'
head(groupinfo)
write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
head(geneInfor)
geneFile='geneFile.txt'
write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)
library(Matrix)
dat_matrix <- as.matrix(dat)

# dat_sparse <- as(dat_matrix, "dgCMatrix")
#2inferCNV####
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=dat_matrix,
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= geneFile,
                                    ref_group_names=NULL)  ## 
# 
out_dir <- "path_to_save_inferCNV_result"
# infercnv_obj <- NormalizeData(infercnv_obj) # 
# infercnv_obj <- ScaleData(infercnv_obj) # 
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=out_dir, 
                             cluster_by_groups=T, 
                             denoise=TRUE,
                             HMM=T,
                             num_threads=8)
#3####
# CNV
cnv_regions <- read.table("G:\\pan_result\\gpu_result_tcga\\UCEC\\scRNA_new\\path_to_save_inferCNV_result\\17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_regions.dat", header = TRUE)
cnv_genes <- read.table("path_to_save_inferCNV_result/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_genes.dat", header = TRUE)

# 
cell_groupings <- read.table("path_to_save_inferCNV_result/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings", header = TRUE)

# CNV
malignant_cells <- cell_groupings[cell_groupings$cell_group_name %in% cnv_regions$cell_group, ]

# 
malignant_cells_need <- malignant_cells$cell
malignant_cells_need
write.table(malignant_cells_need,"UCEC_EPI_cells_E.txt")



















