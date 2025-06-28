#BRCA####
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
setwd("G:\\pan_result\\gpu_result_tcga\\BRCA\\SCT")

list.files(path = "./GSM6592048_M1/")


img = Read10X_Image("./GSM6592048_M1/spatial/", 
                    image.name = "tissue_hires_image.png")
spe = Load10X_Spatial(data.dir = "./GSM6592048_M1/",
                      filename = "GSM6592049_M2_filtered_feature_bc_matrix.h5",
                      assay = "Spatial", 
                      slice = "BRCA",
                      image = img)# 
head(spe@meta.data,2)
# spe@images$slice1@scale.factors$lowres = spe@images$slice1@scale.factors$hires

#SpatialFeaturePlot FeaturePlot
SpatialFeaturePlot(spe, features = "nFeature_Spatial")

T0<-spe

mt.genes <- grep(pattern = "^MT-", x = rownames(T0), value = TRUE)

T0$percent.mito <- (Matrix::colSums(T0@assays[["Spatial"]][mt.genes, ])/Matrix::colSums(T0@assays$Spatial@counts))*100



T0 <- subset(T0,
             # features =genes_to_keep,  #
             subset = nFeature_Spatial > 300 #& percent.mito < 30 #spots
)
plot1 <- VlnPlot(T0, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(T0, features = "nCount_Spatial",pt.size.factor = 3) + 
  theme(legend.position = "right")
wrap_plots(plot1, plot2)

#NormalizeData
T0 <- SCTransform(T0, assay = "Spatial", verbose = FALSE)

T0 <- RunPCA(T0, assay = "SCT", verbose = FALSE)
T0 <- FindNeighbors(T0, reduction = "pca", dims = 1:30)
T0 <- FindClusters(T0, verbose = FALSE)
T0 <- RunUMAP(T0, reduction = "pca", dims = 1:30)

p1 <- DimPlot(T0, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(T0, label = TRUE, label.size = 3)
p1 + p2

T0 <- FindSpatiallyVariableFeatures(T0, assay = "SCT",
                                       features = VariableFeatures(T0)[1:1000],
                                       selection.method = "moransi")
#top.features <- head(SpatiallyVariableFeatures(T0, selection.method = "moransi"), 6)

top.features <- head(T0@commands[["FindSpatiallyVariableFeatures.SCT"]]@params[["features"]], 6)

#  Seurat 

# expression_matrix <- as.matrix(T0[["SCT"]]@data)

SpatialFeaturePlot(T0, features = top.features, 
                   pt.size.factor = 3,
                   ncol = 3, alpha = c(0.1, 1))

head(T0)
BiocManager::install("MCPcounter")
library(MCPcounter)


#HNSC####
#1#
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)

library(jsonlite)
library(png)
library(tidyverse)
library(ggpubr)
setwd("G:/pan_result/gpu_result_tcga/HNSC/SCT")

list.files(path = "./GSM6339633/")


img = Read10X_Image("./GSM6339633/Spatial/", 
                    image.name = "tissue_hires_image.png")
spe = Load10X_Spatial(data.dir = "./GSM6339633/",
                      filename = "GSM6339633_s3_filtered_feature_bc_matrix.h5",
                      assay = "Spatial", 
                      slice = "hnsc",
                      image = img)# 


head(spe@meta.data,2)
# spe@images$slice1@scale.factors$lowres = spe@images$slice1@scale.factors$hires

#SpatialFeaturePlot FeaturePlot
SpatialFeaturePlot(spe, features = "nFeature_Spatial", pt.size = 3)

T0<-spe

mt.genes <- grep(pattern = "^MT-", x = rownames(T0), value = TRUE)

T0$percent.mito <- (Matrix::colSums(T0@assays[["Spatial"]]$counts[mt.genes, ])/Matrix::colSums(T0@assays[["Spatial"]]$counts))*100
#spot5
genes_to_keep <- setdiff(names(which(Matrix::rowSums(T0@assays[["Spatial"]]$counts)>5)),mt.genes)

#3 subsetspot
T0 <- subset(T0,
             features =genes_to_keep,  #
             subset = nFeature_Spatial > 300 & percent.mito < 30 #spots
)
plot1 <- VlnPlot(T0, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(T0, features = "nCount_Spatial",pt.size.factor = 3) + 
  theme(legend.position = "right")
wrap_plots(plot1, plot2)
# install.packages("Rcpp")
# BiocManager::install('glmGamPoi')
# library(glmGamPoi)
#NormalizeData
T0 <- SCTransform(T0, assay = "Spatial", verbose = FALSE)
#3
T0 <- RunPCA(T0, assay = "SCT", verbose = FALSE)
T0 <- FindNeighbors(T0, reduction = "pca", dims = 1:30)
T0 <- FindClusters(T0, verbose = FALSE)
T0 <- RunUMAP(T0, reduction = "pca", dims = 1:30)

p1 <- DimPlot(T0, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(T0, label = TRUE, label.size = 3,pt.size=4)
p1 + p2

ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339633\\GSM6339633_sc_weizhushi.pdf",p1 + p2,width=10,height = 10)

#2 marker 
markergene<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")
Marker2<-markergene$marker
# Marker2 = c("EPCAM",
#             "PECAM1","PLVAP",
#             "COL3A1","COL1A1","COL1A2",
#             "CD79A","CD79B","CD19",
#             "CD3D","CD3E","CD8A","CD4",
#             "C1QA","C1QB","CD163","CD1C"
# )
# list
p3<-DotPlot(T0,features=Marker2) + coord_flip()
p3
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339633\\GSM6339633_markergene.pdf",p3,width=10,height = 10)

p4<-VlnPlot(T0,features = Marker2,pt.size = 0,ncol = 5)
p4
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339633\\GSM6339633_markergene_vio.pdf",p4,width=10,height = 10)


sign_celltype <-read.csv("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339633\\zhushi.csv")
sign_celltype<-sign_celltype$cell


new.cluster.ids <- sign_celltype
names(new.cluster.ids) <- levels(T0)
ST_object <- RenameIdents(T0, new.cluster.ids)

p1 <- DimPlot(ST_object, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(ST_object, label = TRUE, label.size = 3,pt.size.factor = 5)
p1 + p2
SpatialFeaturePlot(ST_object, features = c("EPCAM","PECAM1","COL1A1",
                                           "CD79A","CD19","CD3D","CD3E","C1QA","CD1C") )

ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339633\\GSM6339633_zhushi_umap.pdf",p1,width=10,height = 10)

ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339633\\GSM6339633_zhushi_sct.pdf",p2,width=10,height = 10)

# 
custom_colors <- c("Fibroblasts" = "#E15759", "Epithelium" = "#F28E2B", "B_cell" ="#59A14F","Macrophage"="#8DD3C7")

# 
p2 <- SpatialDimPlot(ST_object, label = TRUE, label.size = 3, pt.size.factor = 5, cols = custom_colors)
p2 

ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339633\\GSM6339633_zhushi_sct_new.pdf",p2,width=10,height = 10)






gene<-read.csv("G:\\pan_result\\gpu_result_tcga\\HNSC\\297featuresweightme.csv",row.names = 1)
# gene<-colnames(gene)


#cluster####
label<-read.table("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_cluster.txt")
gene1<-gene[rownames(label)[which(label$clusters=="1")],]
gene2<-gene[rownames(label)[which(label$clusters=="2")],]
data<-cbind(t(gene1),t(gene2))
T_p<-c()
FC_result<-c()
for(i in 1:nrow(data)){
  T_result<-wilcox.test(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]),as.numeric(data[i,1:nrow(gene1)]),paired = F,var.equal = T)
  p1<-T_result$p.value#p
  # t1<-T_result[[1]]#T
  T_p<-c(T_p,p1)
  #T_value<-c(T_value,t1)
  FC<-mean(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]))/mean(as.numeric(data[i,1:nrow(gene1)]))
  ##fold change########
  FC_result<-c(FC_result,FC)
}
#logFC<-log2(FC_result)
#T_p<-T_p[2:length(T_p)]#NAP
#T_value<-T_value[2:length(T_value)]#NAT
#pFDR
#FDR_p<-p.adjust(T_p,method="fdr",length(T_p))#p
DEG_all<-data.frame(rownames(data),T_p,FC_result,stringsAsFactors = FALSE)

colnames(DEG_all)<-c("Gene","pval","FC") 
####logFC>1logFC<-1adjp<0.05
write.csv(DEG_all,file = "DEG_all.csv")
DEG<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC>2),]#1295
write.table(DEG,"G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\HNSC_DEGUP.txt")
DEG_down<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC<0.5),]#1295

expression_matrix <- ST_object@assays[["Spatial"]]$counts
class(expression_matrix)

gene_expression2<-0
#
DEG<-intersect(DEG$Gene,expression_matrix@Dimnames[[1]])
for (i in 1:length(DEG)) {
  gene_expression <- expression_matrix[DEG[i], ]
  gene_expression2<-gene_expression2+gene_expression
}


ST_object <- AddMetaData(ST_object, metadata = gene_expression2, col.name = "upgene")
p3<-SpatialFeaturePlot(ST_object, features = "upgene" ,pt.size=5)
p3
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339633\\GSM6339633_upgene.pdf",p3,width=10,height = 10)



#0-1
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

# 
normalized_data <- normalize(gene_expression2)
print(normalized_data)
ST_object <- AddMetaData(ST_object, metadata = normalized_data, col.name = "downgene")
p3<-SpatialFeaturePlot(ST_object, features = "downgene" ,pt.size=5)
p3
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339633\\GSM6339633_upgene_new.pdf",p3,width=10,height = 10)
#2#

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)

library(jsonlite)
library(png)
library(tidyverse)
library(ggpubr)
setwd("G:/pan_result/gpu_result_tcga/HNSC/SCT")

list.files(path = "./GSM6339635/")


img = Read10X_Image("./GSM6339635/Spatial/", 
                    image.name = "tissue_hires_image.png")
spe = Load10X_Spatial(data.dir = "./GSM6339635/",
                      filename = "GSM6339635_s5_filtered_feature_bc_matrix.h5",
                      assay = "Spatial", 
                      slice = "hnsc",
                      image = img)# 


head(spe@meta.data,2)
# spe@images$slice1@scale.factors$lowres = spe@images$slice1@scale.factors$hires

#SpatialFeaturePlot FeaturePlot
SpatialFeaturePlot(spe, features = "nFeature_Spatial", pt.size = 3)

T0<-spe

mt.genes <- grep(pattern = "^MT-", x = rownames(T0), value = TRUE)

T0$percent.mito <- (Matrix::colSums(T0@assays[["Spatial"]]$counts[mt.genes, ])/Matrix::colSums(T0@assays[["Spatial"]]$counts))*100
#spot5
genes_to_keep <- setdiff(names(which(Matrix::rowSums(T0@assays[["Spatial"]]$counts)>5)),mt.genes)

#3 subsetspot
T0 <- subset(T0,
             features =genes_to_keep,  #
             subset = nFeature_Spatial > 300 & percent.mito < 30 #spots
)
plot1 <- VlnPlot(T0, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(T0, features = "nCount_Spatial",pt.size.factor = 3) + 
  theme(legend.position = "right")
wrap_plots(plot1, plot2)
# install.packages("Rcpp")
# BiocManager::install('glmGamPoi')
# library(glmGamPoi)
#NormalizeData
T0 <- SCTransform(T0, assay = "Spatial", verbose = FALSE)
#3
T0 <- RunPCA(T0, assay = "SCT", verbose = FALSE)
T0 <- FindNeighbors(T0, reduction = "pca", dims = 1:30)
T0 <- FindClusters(T0, verbose = FALSE)
T0 <- RunUMAP(T0, reduction = "pca", dims = 1:30)

p1 <- DimPlot(T0, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(T0, label = TRUE, label.size = 3,pt.size=4)
p1 + p2

ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339635\\GSM6339635_sc_weizhushi.pdf",p1 + p2,width=10,height = 10)

#2 marker 
markergene<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")
Marker2<-markergene$marker
# Marker2 = c("EPCAM",
#             "PECAM1","PLVAP",
#             "COL3A1","COL1A1","COL1A2",
#             "CD79A","CD79B","CD19",
#             "CD3D","CD3E","CD8A","CD4",
#             "C1QA","C1QB","CD163","CD1C"
# )
# list
p3<-DotPlot(T0,features=Marker2) + coord_flip()
p3
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339635\\GSM6339635_markergene.pdf",p3,width=10,height = 10)

p4<-VlnPlot(T0,features = Marker2,pt.size = 0,ncol = 5)
p4
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339635\\GSM6339635_markergene_vio.pdf",p4,width=10,height = 10)


sign_celltype <-read.csv("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339635\\zhushi.csv")
sign_celltype<-sign_celltype$cell


new.cluster.ids <- sign_celltype
names(new.cluster.ids) <- levels(T0)
ST_object <- RenameIdents(T0, new.cluster.ids)

p1 <- DimPlot(ST_object, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(ST_object, label = TRUE, label.size = 3,pt.size.factor = 5)
p1 + p2
SpatialFeaturePlot(ST_object, features = c("EPCAM","PECAM1","COL1A1",
                                           "CD79A","CD19","CD3D","CD3E","C1QA","CD1C") )

ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339635\\GSM6339635_zhushi_umap.pdf",p1,width=10,height = 10)

ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339635\\GSM6339635_zhushi_sct.pdf",p2,width=10,height = 10)


# 
custom_colors <- c("Fibroblasts" = "#E15759", "Epithelium" = "#F28E2B", "B_cell" ="#59A14F","Macrophage"="#8DD3C7","Plasmocyte"="#BEBADA")

# 
p2 <- SpatialDimPlot(ST_object, label = TRUE, label.size = 3, pt.size.factor = 5, cols = custom_colors)
p2 
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339635\\GSM6339635_zhushi_sct_new.pdf",p2,width=10,height = 10)




gene<-read.csv("G:\\pan_result\\gpu_result_tcga\\HNSC\\297featuresweightme.csv",row.names = 1)
# gene<-colnames(gene)


#cluster####
label<-read.table("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_cluster.txt")
gene1<-gene[rownames(label)[which(label$clusters=="1")],]
gene2<-gene[rownames(label)[which(label$clusters=="2")],]
data<-cbind(t(gene1),t(gene2))
T_p<-c()
FC_result<-c()
for(i in 1:nrow(data)){
  T_result<-wilcox.test(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]),as.numeric(data[i,1:nrow(gene1)]),paired = F,var.equal = T)
  p1<-T_result$p.value#p
  # t1<-T_result[[1]]#T
  T_p<-c(T_p,p1)
  #T_value<-c(T_value,t1)
  FC<-mean(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]))/mean(as.numeric(data[i,1:nrow(gene1)]))
  ##fold change########
  FC_result<-c(FC_result,FC)
}
#logFC<-log2(FC_result)
#T_p<-T_p[2:length(T_p)]#NAP
#T_value<-T_value[2:length(T_value)]#NAT
#pFDR
#FDR_p<-p.adjust(T_p,method="fdr",length(T_p))#p
DEG_all<-data.frame(rownames(data),T_p,FC_result,stringsAsFactors = FALSE)

colnames(DEG_all)<-c("Gene","pval","FC") 
####logFC>1logFC<-1adjp<0.05
write.csv(DEG_all,file = "DEG_all.csv")
DEG<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC>2),]#1295
write.table(DEG,"G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\HNSC_DEGUP.txt")
DEG_down<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC<0.5),]#1295

expression_matrix <- ST_object@assays[["Spatial"]]$counts
class(expression_matrix)

gene_expression2<-0
#
DEG<-intersect(DEG$Gene,expression_matrix@Dimnames[[1]])
for (i in 1:length(DEG)) {
  gene_expression <- expression_matrix[DEG[i], ]
  gene_expression2<-gene_expression2+gene_expression
}


ST_object <- AddMetaData(ST_object, metadata = gene_expression2, col.name = "upgene")
p3<-SpatialFeaturePlot(ST_object, features = "upgene" ,pt.size=5)
p3
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339635\\GSM6339635_upgene.pdf",p3,width=10,height = 10)

#0-1
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

# 
normalized_data <- normalize(gene_expression2)
print(normalized_data)
ST_object <- AddMetaData(ST_object, metadata = normalized_data, col.name = "downgene")
p3<-SpatialFeaturePlot(ST_object, features = "downgene" ,pt.size=5)
p3

ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339635\\GSM6339635_upgene_new.pdf",p3,width=10,height = 10)

#3#
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)

library(jsonlite)
library(png)
library(tidyverse)
library(ggpubr)
setwd("G:/pan_result/gpu_result_tcga/HNSC/SCT")

list.files(path = "./GSM6339632/")


img = Read10X_Image("./GSM6339632/Spatial/", 
                    image.name = "tissue_hires_image.png")
spe = Load10X_Spatial(data.dir = "./GSM6339632/",
                      filename = "GSM6339632_s2_filtered_feature_bc_matrix.h5",
                      assay = "Spatial", 
                      slice = "hnsc",
                      image = img)# 


head(spe@meta.data,2)
# spe@images$slice1@scale.factors$lowres = spe@images$slice1@scale.factors$hires

#SpatialFeaturePlot FeaturePlot
SpatialFeaturePlot(spe, features = "nFeature_Spatial", pt.size = 3)

T0<-spe

mt.genes <- grep(pattern = "^MT-", x = rownames(T0), value = TRUE)

T0$percent.mito <- (Matrix::colSums(T0@assays[["Spatial"]]$counts[mt.genes, ])/Matrix::colSums(T0@assays[["Spatial"]]$counts))*100
#spot5
genes_to_keep <- setdiff(names(which(Matrix::rowSums(T0@assays[["Spatial"]]$counts)>5)),mt.genes)

#3 subsetspot
T0 <- subset(T0,
             features =genes_to_keep,  #
             subset = nFeature_Spatial > 300 & percent.mito < 30 #spots
)
plot1 <- VlnPlot(T0, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(T0, features = "nCount_Spatial",pt.size.factor = 3) + 
  theme(legend.position = "right")
wrap_plots(plot1, plot2)
# install.packages("Rcpp")
# BiocManager::install('glmGamPoi')
# library(glmGamPoi)
#NormalizeData
T0 <- SCTransform(T0, assay = "Spatial", verbose = FALSE)
#3
T0 <- RunPCA(T0, assay = "SCT", verbose = FALSE)
T0 <- FindNeighbors(T0, reduction = "pca", dims = 1:30)
T0 <- FindClusters(T0, verbose = FALSE)
T0 <- RunUMAP(T0, reduction = "pca", dims = 1:30)

p1 <- DimPlot(T0, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(T0, label = TRUE, label.size = 3,pt.size=4)
p1 + p2

ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339632\\GSM6339632_sc_weizhushi.pdf",p1 + p2,width=10,height = 10)

#2 marker 
markergene<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")
Marker2<-markergene$marker
# Marker2 = c("EPCAM",
#             "PECAM1","PLVAP",
#             "COL3A1","COL1A1","COL1A2",
#             "CD79A","CD79B","CD19",
#             "CD3D","CD3E","CD8A","CD4",
#             "C1QA","C1QB","CD163","CD1C"
# )
# list
p3<-DotPlot(T0,features=Marker2) + coord_flip()
p3
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339632\\GSM6339632_markergene.pdf",p3,width=10,height = 10)

p4<-VlnPlot(T0,features = Marker2,pt.size = 0,ncol = 5)
p4
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339632\\GSM6339632_markergene_vio.pdf",p4,width=10,height = 10)


sign_celltype <-read.csv("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339632\\zhushi.csv")
sign_celltype<-sign_celltype$cell


new.cluster.ids <- sign_celltype
names(new.cluster.ids) <- levels(T0)
ST_object <- RenameIdents(T0, new.cluster.ids)

p1 <- DimPlot(ST_object, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(ST_object, label = TRUE, label.size = 5,pt.size.factor = 5)
p1 + p2
SpatialFeaturePlot(ST_object, features = c("EPCAM","PECAM1","COL1A1",
                                           "CD79A","CD19","CD3D","CD3E","C1QA","CD1C") )

ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339632\\GSM6339632_zhushi_umap.pdf",p1,width=10,height = 10)

ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339632\\GSM6339632_zhushi_sct.pdf",p2,width=10,height = 10)

# 
custom_colors <- c("Fibroblasts" = "#E15759", "Epithelium" = "#4E79A7", "B_cell" ="#59A14F","Macrophage"="#FFCC33","Plasmocyte"="#C8BF2C")

# 
p2 <- SpatialDimPlot(ST_object, label = TRUE, label.size = 3, pt.size.factor = 5, cols = custom_colors)
p2 

ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339632\\GSM6339632_zhushi_sct_new.pdf",p2,width=10,height = 10)





gene<-read.csv("G:\\pan_result\\gpu_result_tcga\\HNSC\\297featuresweightme.csv",row.names = 1)
# gene<-colnames(gene)


#cluster####
label<-read.table("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_cluster.txt")
gene1<-gene[rownames(label)[which(label$clusters=="1")],]
gene2<-gene[rownames(label)[which(label$clusters=="2")],]
data<-cbind(t(gene1),t(gene2))
T_p<-c()
FC_result<-c()
for(i in 1:nrow(data)){
  T_result<-wilcox.test(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]),as.numeric(data[i,1:nrow(gene1)]),paired = F,var.equal = T)
  p1<-T_result$p.value#p
  # t1<-T_result[[1]]#T
  T_p<-c(T_p,p1)
  #T_value<-c(T_value,t1)
  FC<-mean(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]))/mean(as.numeric(data[i,1:nrow(gene1)]))
  ##fold change########
  FC_result<-c(FC_result,FC)
}
#logFC<-log2(FC_result)
#T_p<-T_p[2:length(T_p)]#NAP
#T_value<-T_value[2:length(T_value)]#NAT
#pFDR
#FDR_p<-p.adjust(T_p,method="fdr",length(T_p))#p
DEG_all<-data.frame(rownames(data),T_p,FC_result,stringsAsFactors = FALSE)

colnames(DEG_all)<-c("Gene","pval","FC") 
####logFC>1logFC<-1adjp<0.05
write.csv(DEG_all,file = "DEG_all.csv")
DEG<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC>2),]#1295
write.table(DEG,"G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\HNSC_DEGUP.txt")
DEG_down<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC<0.5),]#1295

expression_matrix <- ST_object@assays[["Spatial"]]$counts
class(expression_matrix)

gene_expression2<-0
#
DEG<-intersect(DEG$Gene,expression_matrix@Dimnames[[1]])
for (i in 1:length(DEG)) {
  gene_expression <- expression_matrix[DEG[i], ]
  gene_expression2<-gene_expression2+gene_expression
}


ST_object <- AddMetaData(ST_object, metadata = gene_expression2, col.name = "upgene")
p3<-SpatialFeaturePlot(ST_object, features = "upgene" ,pt.size=5)
p3
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339632\\GSM6339632_upgene.pdf",p3,width=10,height = 10)

#0-1
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

# 
normalized_data <- normalize(gene_expression2)
print(normalized_data)
ST_object <- AddMetaData(ST_object, metadata = normalized_data, col.name = "upgene")
p3<-SpatialFeaturePlot(ST_object, features = "upgene" ,pt.size=5)
p3
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339632\\GSM6339632_upgene_new.pdf",p3,width=10,height = 10)

#4#
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)

library(jsonlite)
library(png)
library(tidyverse)
library(ggpubr)
setwd("G:/pan_result/gpu_result_tcga/HNSC/SCT")

list.files(path = "./GSM6339631/")


img = Read10X_Image("./GSM6339631/Spatial/", 
                    image.name = "tissue_hires_image.png")
spe = Load10X_Spatial(data.dir = "./GSM6339631/",
                      filename = "GSM6339631_s1_filtered_feature_bc_matrix.h5",
                      assay = "Spatial", 
                      slice = "hnsc",
                      image = img)# 


head(spe@meta.data,2)
# spe@images$slice1@scale.factors$lowres = spe@images$slice1@scale.factors$hires

#SpatialFeaturePlot FeaturePlot
SpatialFeaturePlot(spe, features = "nFeature_Spatial", pt.size = 3)

T0<-spe

mt.genes <- grep(pattern = "^MT-", x = rownames(T0), value = TRUE)

T0$percent.mito <- (Matrix::colSums(T0@assays[["Spatial"]]$counts[mt.genes, ])/Matrix::colSums(T0@assays[["Spatial"]]$counts))*100
#spot5
genes_to_keep <- setdiff(names(which(Matrix::rowSums(T0@assays[["Spatial"]]$counts)>5)),mt.genes)

#3 subsetspot
T0 <- subset(T0,
             features =genes_to_keep,  #
             subset = nFeature_Spatial > 300 & percent.mito < 30 #spots
)
plot1 <- VlnPlot(T0, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(T0, features = "nCount_Spatial",pt.size.factor = 3) + 
  theme(legend.position = "right")
wrap_plots(plot1, plot2)
# install.packages("Rcpp")
# BiocManager::install('glmGamPoi')
# library(glmGamPoi)
#NormalizeData
T0 <- SCTransform(T0, assay = "Spatial", verbose = FALSE)
#3
T0 <- RunPCA(T0, assay = "SCT", verbose = FALSE)
T0 <- FindNeighbors(T0, reduction = "pca", dims = 1:30)
T0 <- FindClusters(T0, verbose = FALSE)
T0 <- RunUMAP(T0, reduction = "pca", dims = 1:30)

p1 <- DimPlot(T0, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(T0, label = TRUE, label.size = 3,pt.size=4)
p1 + p2

ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339631\\GSM6339631_sc_weizhushi.pdf",p1 + p2,width=10,height = 10)

#2 marker 
markergene<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")
Marker2<-markergene$marker
# Marker2 = c("EPCAM",
#             "PECAM1","PLVAP",
#             "COL3A1","COL1A1","COL1A2",
#             "CD79A","CD79B","CD19",
#             "CD3D","CD3E","CD8A","CD4",
#             "C1QA","C1QB","CD163","CD1C"
# )
# list
p3<-DotPlot(T0,features=Marker2) + coord_flip()
p3
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339631\\GSM6339631_markergene.pdf",p3,width=10,height = 10)

p4<-VlnPlot(T0,features = Marker2,pt.size = 0,ncol = 5)
p4
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339631\\GSM6339631_markergene_vio.pdf",p4,width=10,height = 10)


sign_celltype <-read.csv("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339631\\zhushi.csv")
sign_celltype<-sign_celltype$cell


new.cluster.ids <- sign_celltype
names(new.cluster.ids) <- levels(T0)
ST_object <- RenameIdents(T0, new.cluster.ids)

p1 <- DimPlot(ST_object, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(ST_object, label = TRUE, label.size = 5,pt.size.factor = 5)
p1 + p2
# SpatialFeaturePlot(ST_object, features = c("EPCAM","PECAM1","COL1A1",
#                                            "CD79A","CD19","CD3D","CD3E","C1QA","CD1C") )

ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339631\\GSM6339631_zhushi_umap.pdf",p1,width=10,height = 10)

ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339631\\GSM6339631_zhushi_sct.pdf",p2,width=10,height = 10)

# 
custom_colors <- c("Fibroblasts" = "#E15759", "Epithelium" = "#4E79A7", "B_cell" ="#59A14F","Macrophage"="#FFCC33","Plasmocyte"="#C8BF2C")

# 
p2 <- SpatialDimPlot(ST_object, label = TRUE, label.size = 3, pt.size.factor = 5, cols = custom_colors)
p2 
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339631\\GSM6339631_zhushi_sct_new.pdf",p2,width=10,height = 10)



gene<-read.csv("G:\\pan_result\\gpu_result_tcga\\HNSC\\297featuresweightme.csv",row.names = 1)
# gene<-colnames(gene)


#cluster####
label<-read.table("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_cluster.txt")
gene1<-gene[rownames(label)[which(label$clusters=="1")],]
gene2<-gene[rownames(label)[which(label$clusters=="2")],]
data<-cbind(t(gene1),t(gene2))
T_p<-c()
FC_result<-c()
for(i in 1:nrow(data)){
  T_result<-wilcox.test(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]),as.numeric(data[i,1:nrow(gene1)]),paired = F,var.equal = T)
  p1<-T_result$p.value#p
  # t1<-T_result[[1]]#T
  T_p<-c(T_p,p1)
  #T_value<-c(T_value,t1)
  FC<-mean(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]))/mean(as.numeric(data[i,1:nrow(gene1)]))
  ##fold change########
  FC_result<-c(FC_result,FC)
}
#logFC<-log2(FC_result)
#T_p<-T_p[2:length(T_p)]#NAP
#T_value<-T_value[2:length(T_value)]#NAT
#pFDR
#FDR_p<-p.adjust(T_p,method="fdr",length(T_p))#p
DEG_all<-data.frame(rownames(data),T_p,FC_result,stringsAsFactors = FALSE)

colnames(DEG_all)<-c("Gene","pval","FC") 
####logFC>1logFC<-1adjp<0.05
write.csv(DEG_all,file = "DEG_all.csv")
DEG<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC>2),]#1295
write.table(DEG,"G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\HNSC_DEGUP.txt")
DEG_down<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC<0.5),]#1295

expression_matrix <- ST_object@assays[["Spatial"]]$counts
class(expression_matrix)

gene_expression2<-0
#
DEG<-intersect(DEG$Gene,expression_matrix@Dimnames[[1]])
for (i in 1:length(DEG)) {
  gene_expression <- expression_matrix[DEG[i], ]
  gene_expression2<-gene_expression2+gene_expression
}


ST_object <- AddMetaData(ST_object, metadata = gene_expression2, col.name = "upgene")
p3<-SpatialFeaturePlot(ST_object, features = "upgene" ,pt.size=6)
p3
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339631\\GSM6339631_upgene.pdf",p3,width=10,height = 10)



#0-1
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

# 
normalized_data <- normalize(gene_expression2)
print(normalized_data)
ST_object <- AddMetaData(ST_object, metadata = normalized_data, col.name = "downgene")
p3<-SpatialFeaturePlot(ST_object, features = "downgene" ,pt.size=5)
p3
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\SCT\\GSM6339631\\GSM6339631_upgene.pdf",p3,width=10,height = 10)

