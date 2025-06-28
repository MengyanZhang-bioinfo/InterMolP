#BRCA####
rm(list=ls())
library(Seurat)
library(cowplot)
library(dplyr)
# library(glmGamPoi)
library(data.table)
setwd("G:/pan_result/gpu_result_tcga/BRCA/scRNA_new")
#  R.utils 
library(R.utils)
data<-readRDS("E:\\\\sample.obj.rds")


##---harmony---#
library(harmony)
data=RunHarmony(data,group.by="orig.ident",assay.use="RNA",max.inter.harmony=20)



###--cluster---#
#rm(list=ls())
#setwd("E:/test_death/test1/01.analysis221126_v2_2/01.Seurat.v2/merge")
library(Seurat)
library(ggthemes)
library(dplyr)
library(ggplot2)
library(clustree)

pc.num=1:50
data=RunTSNE(data, reduction="harmony", dims=pc.num) %>% 
  RunUMAP(reduction="harmony", dims=pc.num) %>%
  FindNeighbors(reduction="harmony", dims=pc.num)

data= FindClusters(data,resolution=0.6)


#length(levels(Idents(data)))
#table(data@active.ident) 
dir.create("harmony_0.4_20")
setwd("harmony_0.4_20")

sample_col=c("#59A14F","#4E79A7","#F28E2B","#E15759","#4E80AB",
             "#FFCC33","#C8BF2C","#8c310a","#228B22","#815c94",
             "#3692a8" ,"#c95043" ,"#b5574d","#548faf","#e97e33",
             "#76a695",'#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', 
             '#80B1D3', '#FDB462')

#sample_col=c("#BC3C29FF","#E18727FF","#20854EFF", "#7876B1FF" ,"#6F99ADFF","#FB8072", "#EDC948","#B2DF8A",  "#BEBADA",  "#BEBADA","#4E79A7", "#BAB0AC","#59A14F","#4E79A7")
##color
p=DimPlot(data,reduction = "umap",label=T,raster = FALSE)
p
ggsave("UMAP_cluster_harmony.pdf",p,width=9,height=9)

p3=DimPlot(data,reduction = "tsne",label=T,raster = FALSE)
p3
ggsave("TSNE_cluster_harmony.pdf",p3,width=9,height=9)
saveRDS(data,"Seurat.rds")

#--------
markergene<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")

markergene<-markergene$marker
library(stringr)
library(ggplot2)
# genes_to_check = unique(intersect(rownames(data),markergene))
genes_to_check = markergene
dir.create("shoudong")
setwd("shoudong")

P13 <- DotPlot(data, features = genes_to_check,
               assay='RNA' )  + coord_flip()
P13
P14 <- VlnPlot(object = data, features =genes_to_check,log =T )
P14
P15 <- FeaturePlot(object = data, features=genes_to_check )
P15

Sample="GSE161529"

pdf(paste(Sample,"_marker_DotPlot.pdf",sep=""),width=10,height=10)
print(P13)
dev.off()

pdf(paste(Sample,"_marker_VlnPlot.pdf",sep=""),width=20,height=20)
print(P14)
dev.off()

pdf(paste(Sample,"_marker_FeaturePlot.pdf",sep=""),width=10,height=10)
print(P15)
dev.off()
ggsave(filename = "mk.pdf",P13/P14/P15,height = 15,width = 10)

# 
unique(data@meta.data[["seurat_clusters"]])
manual_celltype<-read.csv("marker_celltype.csv")

new.cluster.ids <- manual_celltype$celltype
names(new.cluster.ids) <- levels(data)
seu.obj <- RenameIdents(data, new.cluster.ids)

p16 <- DimPlot(seu.obj, 
               reduction = "umap", 
               label = TRUE, 
               pt.size = 0.5) + NoLegend()
p16
p17 <- DimPlot(seu.obj, 
               reduction = "tsne", 
               label = TRUE, 
               pt.size = 0.5) + NoLegend()
p17

pdf(paste(Sample,"_marker_umap.pdf",sep=""),width=20,height=20)
print(p16)
dev.off()

pdf(paste(Sample,"_marker_tsnet.pdf",sep=""),width=10,height=10)
print(p17)
dev.off()
save(seu.obj,file = paste0(Sample,".seu.obj.Rdata"))


#COAD####
rm(list=ls())
library(Seurat)
library(cowplot)
library(dplyr)
# library(glmGamPoi)
library(data.table)
setwd("G:/pan_result/gpu_result_tcga/COAD/scRNA_new/GSE188711_RAW")
#  R.utils 
library(R.utils)

files <- c(list.files("G:/pan_result/gpu_result_tcga/COAD/scRNA_new/GSE188711_RAW/"))

# 
prefixes <- unique(sub("_[^.]+\\..+$", "", files))

base_path<-"G:/pan_result/gpu_result_tcga/COAD/scRNA_new/GSE188711_RAW"

# # 
# for (prefix in prefixes) {
#   # 
#   target_path <- file.path(base_path, prefix)
#   if (!dir.exists(target_path)) {
#     dir.create(target_path)
#   }
#   
#   # 
#   current_files <- list.files(base_path, pattern = paste0(prefix, "_"), full.names = TRUE)
#   
#   # 
#   for (file in current_files) {
#     # 
#     gunzip(file, remove = FALSE, overwrite = TRUE)
#     
#     # 
#     unzipped_file <- sub("^(.*)\\.gz$", "\\1", file)
#     
#     # 
#     if (grepl("barcodes", unzipped_file)) {
#       new_name <- "barcodes.tsv"
#     } else if (grepl("features", unzipped_file)) {
#       new_name <- "features.tsv"
#     } else if (grepl("matrix", unzipped_file)) {
#       new_name <- "matrix.mtx"
#     }
#     
#     # 
#     new_file_path <- file.path(target_path, new_name)
#     
#     # 
#     file.rename(unzipped_file, new_file_path)
#   }
# }


prefix<-prefixes[1]
file<-current_files[1]
# 
for (prefix in prefixes) {
  # 
  target_path <- file.path(base_path, prefix)
  if (!dir.exists(target_path)) {
    dir.create(target_path)
  }

  # 
  current_files <- list.files(base_path, pattern = paste0(prefix, "_"), full.names = TRUE)

  # 
  for (file in current_files) {


    # 
    # unzipped_file <- sub("^(.*)\\.gz$", "\\1", file)
    unzipped_file<-file
    # 
    if (grepl("barcodes", unzipped_file)) {
      new_name <- "barcodes.tsv.gz"
    } else if (grepl("features", unzipped_file)) {
      new_name <- "features.tsv.gz"
    } else if (grepl("matrix", unzipped_file)) {
      new_name <- "matrix.mtx.gz"
    }

    # 
    new_file_path <- file.path(target_path, new_name)

    # 
    file.rename(unzipped_file, new_file_path)
  }
}




library(Seurat)
library(cowplot)
library(dplyr)
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(data.table)
setwd("G:/pan_result/gpu_result_tcga/COAD/scRNA_new")
tmp = list.dirs('GSE188711_RAW/')[-1]
tmp
sce =CreateSeuratObject(counts =  Read10X(tmp) ,
                        #project = gsub('.gz','',gsub('^GSM[0-9]*_','') ) ,
                        min.cells = 5,
                        min.features = 300 )
GetAssay(sce,assay = "RNA")#17524   cells  22722  features 2
data<-sce

########
rownames(data)[grepl('^mt-',rownames(data),ignore.case = T)]
rownames(data)[grepl('^Rp[sl]',rownames(data),ignore.case = T)]
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
data <- subset(data, subset = nFeature_RNA > 200 & nCount_RNA > 2500 & percent.mt < 20)

##---NormalizeData--#
data <- NormalizeData(data, normalization.method =  "LogNormalize", 
                      scale.factor = 10000)
data=FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)

#all.genes <- rownames(data)
data <- ScaleData(data)

##---PCA--#

data=RunPCA(data, features = VariableFeatures(data),dims = 1:30)

##---harmony---#
library(harmony)
data=RunHarmony(data,group.by="orig.ident",assay.use="RNA",max.inter.harmony=20)



###--cluster---#
#rm(list=ls())
#setwd("E:/test_death/test1/01.analysis221126_v2_2/01.Seurat.v2/merge")
library(Seurat)
library(ggthemes)
library(dplyr)
library(ggplot2)
library(clustree)

pc.num=1:50
data=RunTSNE(data, reduction="harmony", dims=pc.num) %>% 
  RunUMAP(reduction="harmony", dims=pc.num) %>%
  FindNeighbors(reduction="harmony", dims=pc.num)

data= FindClusters(data,resolution=0.6)


#length(levels(Idents(data)))
#table(data@active.ident) 
dir.create("harmony_0.4_20")
setwd("harmony_0.4_20")

sample_col=c("#59A14F","#4E79A7","#F28E2B","#E15759","#4E80AB",
             "#FFCC33","#C8BF2C","#8c310a","#228B22","#815c94",
             "#3692a8" ,"#c95043" ,"#b5574d","#548faf","#e97e33",
             "#76a695",'#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', 
             '#80B1D3', '#FDB462')

#sample_col=c("#BC3C29FF","#E18727FF","#20854EFF", "#7876B1FF" ,"#6F99ADFF","#FB8072", "#EDC948","#B2DF8A",  "#BEBADA",  "#BEBADA","#4E79A7", "#BAB0AC","#59A14F","#4E79A7")
##color
p=DimPlot(data,reduction = "umap",label=T,raster = FALSE)
p
ggsave("UMAP_cluster_harmony.pdf",p,width=9,height=9)

p3=DimPlot(data,reduction = "tsne",label=T,raster = FALSE)
p3
ggsave("TSNE_cluster_harmony.pdf",p3,width=9,height=9)
saveRDS(data,"Seurat.rds")
##---cluster marker---#
clus=data.frame(data@active.ident)###cluster 
colnames(clus)=c("Cluster")
dir.create("Markers")
write.csv(clus,'Markers/Cluster.csv',quote=F)


#######---------- ---------

markergene<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")
data<-readRDS("Seurat.rds")
markergene<-markergene$marker
library(stringr)
library(ggplot2)
# genes_to_check = unique(intersect(rownames(data),markergene))
genes_to_check = markergene
dir.create("shoudong")
setwd("shoudong")

P13 <- DotPlot(data, features = genes_to_check,
               assay='RNA' )  + coord_flip()
P13
P14 <- VlnPlot(object = data, features =genes_to_check,log =T )
P14
P15 <- FeaturePlot(object = data, features=genes_to_check )
P15
Sample="GSE188711"

pdf(paste(Sample,"_marker_DotPlot.pdf",sep=""),width=10,height=10)
print(P13)
dev.off()

pdf(paste(Sample,"_marker_VlnPlot.pdf",sep=""),width=20,height=20)
print(P14)
dev.off()

pdf(paste(Sample,"_marker_FeaturePlot.pdf",sep=""),width=10,height=10)
print(P15)
dev.off()
ggsave(filename = "mk.pdf",P13/P14/P15,height = 15,width = 10)

# 
unique(data@meta.data[["seurat_clusters"]])
manual_celltype<-read.csv("marker_celltype.csv")

new.cluster.ids <- manual_celltype$celltype
names(new.cluster.ids) <- levels(data)
seu.obj <- RenameIdents(data, new.cluster.ids)

p16 <- DimPlot(seu.obj, 
               reduction = "umap", 
               label = TRUE, 
               pt.size = 0.5) + NoLegend()
p16
p17 <- DimPlot(seu.obj, 
               reduction = "tsne", 
               label = TRUE, 
               pt.size = 0.5) + NoLegend()
p17

pdf(paste(Sample,"_marker_umap.pdf",sep=""),width=20,height=20)
print(p16)
dev.off()

pdf(paste(Sample,"_marker_tsnet.pdf",sep=""),width=10,height=10)
print(p17)
dev.off()
save(seu.obj,file = paste0(Sample,".seu.obj.Rdata"))

#ESCA####
rm(list=ls())
library(Seurat)
library(cowplot)
library(dplyr)
# library(glmGamPoi)
library(data.table)
setwd("G:/pan_result/gpu_result_tcga/ESCA/scRNA_new")
#  R.utils 
library(R.utils)
load("H:\\\\database\\SCDTD\\data\\chenlong\\Esophageal Cancer\\disease\\GSE196756_geneplotneed.Rdata")
data<-data.combined.new

##---harmony---#
library(harmony)
data=RunHarmony(data,group.by="orig.ident",assay.use="RNA",max.inter.harmony=20)



###--cluster---#
#rm(list=ls())
#setwd("E:/test_death/test1/01.analysis221126_v2_2/01.Seurat.v2/merge")
library(Seurat)
library(ggthemes)
library(dplyr)
library(ggplot2)
library(clustree)

pc.num=1:50
data=RunTSNE(data, reduction="harmony", dims=pc.num) %>% 
  RunUMAP(reduction="harmony", dims=pc.num) %>%
  FindNeighbors(reduction="harmony", dims=pc.num)

data= FindClusters(data,resolution=0.6)


#length(levels(Idents(data)))
#table(data@active.ident) 
dir.create("harmony_0.4_20")
setwd("harmony_0.4_20")

sample_col=c("#59A14F","#4E79A7","#F28E2B","#E15759","#4E80AB",
             "#FFCC33","#C8BF2C","#8c310a","#228B22","#815c94",
             "#3692a8" ,"#c95043" ,"#b5574d","#548faf","#e97e33",
             "#76a695",'#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', 
             '#80B1D3', '#FDB462')

#sample_col=c("#BC3C29FF","#E18727FF","#20854EFF", "#7876B1FF" ,"#6F99ADFF","#FB8072", "#EDC948","#B2DF8A",  "#BEBADA",  "#BEBADA","#4E79A7", "#BAB0AC","#59A14F","#4E79A7")
##color
p=DimPlot(data,reduction = "umap",label=T,raster = FALSE)
p
ggsave("UMAP_cluster_harmony.pdf",p,width=9,height=9)

p3=DimPlot(data,reduction = "tsne",label=T,raster = FALSE)
p3
ggsave("TSNE_cluster_harmony.pdf",p3,width=9,height=9)
saveRDS(data,"Seurat.rds")

#--------
markergene<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")

markergene<-markergene$marker
library(stringr)
library(ggplot2)
# genes_to_check = unique(intersect(rownames(data),markergene))
genes_to_check = markergene
dir.create("shoudong")
setwd("shoudong")

P13 <- DotPlot(data, features = genes_to_check,
               assay='RNA' )  + coord_flip()
P13
P14 <- VlnPlot(object = data, features =genes_to_check,log =T )
P14
P15 <- FeaturePlot(object = data, features=genes_to_check )
P15

Sample="GSE196756"

pdf(paste(Sample,"_marker_DotPlot.pdf",sep=""),width=10,height=10)
print(P13)
dev.off()

pdf(paste(Sample,"_marker_VlnPlot.pdf",sep=""),width=20,height=20)
print(P14)
dev.off()

pdf(paste(Sample,"_marker_FeaturePlot.pdf",sep=""),width=10,height=10)
print(P15)
dev.off()
ggsave(filename = "mk.pdf",P13/P14/P15,height = 15,width = 10)

# 
unique(data@meta.data[["seurat_clusters"]])
manual_celltype<-read.csv("marker_celltype.csv")

new.cluster.ids <- manual_celltype$celltype
names(new.cluster.ids) <- levels(data)
seu.obj <- RenameIdents(data, new.cluster.ids)

p16 <- DimPlot(seu.obj, 
               reduction = "umap", 
               label = TRUE, 
               pt.size = 0.5) + NoLegend()
p16
p17 <- DimPlot(seu.obj, 
               reduction = "tsne", 
               label = TRUE, 
               pt.size = 0.5) + NoLegend()
p17

pdf(paste(Sample,"_marker_umap.pdf",sep=""),width=20,height=20)
print(p16)
dev.off()

pdf(paste(Sample,"_marker_tsnet.pdf",sep=""),width=10,height=10)
print(p17)
dev.off()
save(seu.obj,file = paste0(Sample,".seu.obj.Rdata"))


#HNSC####
rm(list=ls())
library(Seurat)
library(cowplot)
library(dplyr)
# library(glmGamPoi)
library(data.table)
setwd("G:/pan_result/gpu_result_tcga/HNSC/scRNA_new")
# read.table()txt.gz
meta<- read.table(gzfile("GSE181919_UMI_counts.txt.gz"),  header = TRUE, sep = "\t")
seurat_data<- fread("GSE181919_Barcode_metadata.txt",data.table = F)

meta_n = subset(meta,Tissue == "Normal")



# CreateSeuratObject()Seurat
data <- CreateSeuratObject(counts = meta,
                           min.features = 200,
                           min.cells = 3, 
                           project = "GSE181919")
#rm(seurat_obj)
GetAssay(data,assay = "RNA")#  cells  features 

########
rownames(data)[grepl('^mt-',rownames(data),ignore.case = T)]
rownames(data)[grepl('^Rp[sl]',rownames(data),ignore.case = T)]
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
data <- subset(data, subset = nFeature_RNA > 200 & nCount_RNA > 2500 & percent.mt < 20)

##---NormalizeData--#
data <- NormalizeData(data, normalization.method =  "LogNormalize", 
                      scale.factor = 10000)
data=FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)

#all.genes <- rownames(data)
data <- ScaleData(data)

##---PCA--#

data=RunPCA(data, features = VariableFeatures(data),dims = 1:30)




###--cluster---#
#rm(list=ls())
#setwd("E:/test_death/test1/01.analysis221126_v2_2/01.Seurat.v2/merge")
library(Seurat)
library(ggthemes)
library(dplyr)
library(ggplot2)
library(clustree)

pc.num=1:50
data=RunTSNE(data, reduction="pca", dims=pc.num) %>% 
  RunUMAP(reduction="pca", dims=pc.num) %>%
  FindNeighbors(reduction="pca", dims=pc.num)

data= FindClusters(data,resolution=0.6)


#length(levels(Idents(data)))
#table(data@active.ident) 
dir.create("harmony_0.4_20")
setwd("harmony_0.4_20")

sample_col=c("#59A14F","#4E79A7","#F28E2B","#E15759","#4E80AB",
             "#FFCC33","#C8BF2C","#8c310a","#228B22","#815c94",
             "#3692a8" ,"#c95043" ,"#b5574d","#548faf","#e97e33",
             "#76a695",'#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', 
             '#80B1D3', '#FDB462')

#sample_col=c("#BC3C29FF","#E18727FF","#20854EFF", "#7876B1FF" ,"#6F99ADFF","#FB8072", "#EDC948","#B2DF8A",  "#BEBADA",  "#BEBADA","#4E79A7", "#BAB0AC","#59A14F","#4E79A7")
##color
p=DimPlot(data,reduction = "umap",label=T,raster = FALSE)
p
p1=DimPlot(data,reduction = "umap",group.by="orig.ident",split.by="orig.ident",ncol=3,cols=sample_col,raster = FALSE)+NoLegend()
p1
ggsave("UMAP_Samples_harmony.pdf",p1,width=9,height=9)
ggsave("UMAP_cluster_harmony.pdf",p,width=9,height=9)
p2=DimPlot(data,reduction = "tsne",group.by="orig.ident",split.by="orig.ident",ncol=3,cols=sample_col,raster = FALSE)+NoLegend()
p2

p3=DimPlot(data,reduction = "tsne",label=T,raster = FALSE)
p3
ggsave("TSNE_Samples_harmony.pdf",p2,width=9,height=9)
ggsave("TSNE_cluster_harmony.pdf",p3,width=9,height=9)
saveRDS(data,"Seurat.rds")
##---cluster marker---#
clus=data.frame(data@active.ident)###cluster 
colnames(clus)=c("Cluster")
dir.create("Markers")
write.csv(clus,'Markers/Cluster.csv',quote=F)


#######---------- ---------

markergene<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")

markergene<-markergene$marker
library(stringr)
library(ggplot2)
genes_to_check = unique(intersect(rownames(data),markergene))

dir.create("shoudong")
setwd("shoudong")

P13 <- DotPlot(data, features = genes_to_check,
               assay='RNA' )  + coord_flip()
P13
P14 <- VlnPlot(object = data, features =genes_to_check,log =T )
P14
P15 <- FeaturePlot(object = data, features=genes_to_check )
P15
Sample="GSE181919"

pdf(paste(Sample,"_marker_DotPlot.pdf",sep=""),width=10,height=10)
print(P13)
dev.off()

pdf(paste(Sample,"_marker_VlnPlot.pdf",sep=""),width=20,height=20)
print(P14)
dev.off()

pdf(paste(Sample,"_marker_FeaturePlot.pdf",sep=""),width=10,height=10)
print(P15)
dev.off()
ggsave(filename = "mk.pdf",P13/P14/P15,height = 15,width = 10)

# 
unique(data@meta.data[["seurat_clusters"]])
manual_celltype<-read.csv("marker_celltype.csv")

new.cluster.ids <- manual_celltype$celltype
names(new.cluster.ids) <- levels(data)
seu.obj <- RenameIdents(data, new.cluster.ids)

p16 <- DimPlot(seu.obj, 
              reduction = "umap", 
              label = TRUE, 
              pt.size = 0.5) + NoLegend()
p16
p17 <- DimPlot(seu.obj, 
              reduction = "tsne", 
              label = TRUE, 
              pt.size = 0.5) + NoLegend()
p17

pdf(paste(Sample,"_marker_umap.pdf",sep=""),width=20,height=20)
print(p16)
dev.off()

pdf(paste(Sample,"_marker_tsnet.pdf",sep=""),width=10,height=10)
print(p17)
dev.off()
save(seu.obj,file = paste0(Sample,".seu.obj.Rdata"))


#######----------singleR ---------
##---singleR HCL human ---#

Sample="GSE181919"
type="umap"
library(Seurat)
library(dplyr)
library(ggplot2)
library(SingleR)
library(tibble)
library(celldex)

options(future.globals.maxSize= 1300 * 1024^2)
### read 10X data
#---- five human ref databse -----
load("E:\\\\\\HumanPrimaryCellAtlas_hpca.se_human.RData")
hpca.se <- get(load("E:\\\\HumanPrimaryCellAtlas_hpca.se_human.RData"))

#hpca.se <- HumanPrimaryCellAtlasData()
REF=hpca.se
ctrl=data
ctrl@meta.data$cell.type=Idents(ctrl)
test=as.SingleCellExperiment(ctrl)
## Annot
Anno=SingleR(test = test,
             ref = REF,
             labels = REF$label.main,
             method = "cluster",
             cluster = test$cell.type
)
Anno$cluster=rownames(Anno)
### extract anno 
fin=Anno %>% dplyr::tbl_df() %>% dplyr::select(cluster,labels)
new.cluster.ids=fin$labels
names(new.cluster.ids)=fin$cluster
ctrl=RenameIdents(object = ctrl,new.cluster.ids)
dir.create("singleR")
setwd("singleR")

write.table(fin,file="cluster_after.xls",sep="\t",row.names = FALSE,quote=F)


dir.create("singleR")
setwd("singleR")
p1=DimPlot(ctrl, reduction = "umap", pt.size = 0.5,label=TRUE)
pdf(paste(Sample,"_UMAP_Anno.pdf",sep=""),width=10,height=6)
print(p1)
# + ggtitle(label = "UMAP")
dev.off()
png(paste(Sample,"_UMAP_Anno.png",sep=""))
# + ggtitle(label = "UMAP")
print(p1)
dev.off()

p2=DimPlot(ctrl, reduction = "tsne", pt.size = 0.5,label=TRUE)
pdf(paste(Sample, "_TSNE_Anno.pdf",sep=""),width=10,height=6)
print(p2)
dev.off()
png(paste(Sample,"_TSNE_Anno.png",sep=""))
print(p1)
# + ggtitle(label = "UMAP")
dev.off()

ctrl@meta.data$cell_type=Idents(ctrl)
#saveRDS(ctrl, file = paste(Sample,"_Seurat_Anno_singleR.rds",sep=""))
pdf(paste(Sample,"_SingleR_heatmap.pdf",sep=""),width=8,height=5)
plotScoreHeatmap(Anno)
dev.off()
tiff(paste(Sample,"_SingleR_heatmap.tiff",sep=""))
plotScoreHeatmap(Anno)
dev.off()
write.csv(Idents(ctrl),file = paste(Sample,"_cluster_Anno.csv",sep=""),quote=F)

data.combined.new = ctrl

save(data.combined.new,file = paste(Sample,"_geneplotneed.Rdata",sep=""))


#KIPC####
rm(list=ls())
library(Seurat)
library(cowplot)
library(dplyr)
# library(glmGamPoi)
library(data.table)
setwd("G:/pan_result/gpu_result_tcga/KIRC/scRNA_new")
#  R.utils 
library(R.utils)
load("H:\\\\database\\SCDTD\\data\\chenlong\\Kidney Chromophobe\\disease\\HP\\GSE159115_geneplotneed.Rdata")
data<-data.combined.new

##---harmony---#
library(harmony)
data=RunHarmony(data,group.by="orig.ident",assay.use="RNA",max.inter.harmony=20)



###--cluster---#
#rm(list=ls())
#setwd("E:/test_death/test1/01.analysis221126_v2_2/01.Seurat.v2/merge")
library(Seurat)
library(ggthemes)
library(dplyr)
library(ggplot2)
library(clustree)

pc.num=1:50
data=RunTSNE(data, reduction="harmony", dims=pc.num) %>% 
  RunUMAP(reduction="harmony", dims=pc.num) %>%
  FindNeighbors(reduction="harmony", dims=pc.num)

data= FindClusters(data,resolution=0.6)


#length(levels(Idents(data)))
#table(data@active.ident) 
dir.create("harmony_0.4_20")
setwd("harmony_0.4_20")

sample_col=c("#59A14F","#4E79A7","#F28E2B","#E15759","#4E80AB",
             "#FFCC33","#C8BF2C","#8c310a","#228B22","#815c94",
             "#3692a8" ,"#c95043" ,"#b5574d","#548faf","#e97e33",
             "#76a695",'#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', 
             '#80B1D3', '#FDB462')

#sample_col=c("#BC3C29FF","#E18727FF","#20854EFF", "#7876B1FF" ,"#6F99ADFF","#FB8072", "#EDC948","#B2DF8A",  "#BEBADA",  "#BEBADA","#4E79A7", "#BAB0AC","#59A14F","#4E79A7")
##color
p=DimPlot(data,reduction = "umap",label=T,raster = FALSE)
p
ggsave("UMAP_cluster_harmony.pdf",p,width=9,height=9)

p3=DimPlot(data,reduction = "tsne",label=T,raster = FALSE)
p3
ggsave("TSNE_cluster_harmony.pdf",p3,width=9,height=9)
saveRDS(data,"Seurat.rds")
data<-readRDS("G:\\pan_result\\gpu_result_tcga\\KIRC\\scRNA_new\\harmony_0.4_20\\Seurat.rds")
#--------
setwd("G:\\pan_result\\gpu_result_tcga\\KIRC\\scRNA_new\\harmony_0.4_20\\shoudong")
markergene<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")

markergene<-markergene$marker
library(stringr)
library(ggplot2)
# genes_to_check = unique(intersect(rownames(data),markergene))
genes_to_check = markergene
dir.create("shoudong")
setwd("shoudong")

P13 <- DotPlot(data, features = genes_to_check,
               assay='RNA' )  + coord_flip()
P13
P14 <- VlnPlot(object = data, features =genes_to_check,log =T )
P14
P15 <- FeaturePlot(object = data, features=genes_to_check )
P15

Sample="GSE159115"

pdf(paste(Sample,"_marker_DotPlot.pdf",sep=""),width=10,height=10)
print(P13)
dev.off()

pdf(paste(Sample,"_marker_VlnPlot.pdf",sep=""),width=20,height=20)
print(P14)
dev.off()

pdf(paste(Sample,"_marker_FeaturePlot.pdf",sep=""),width=10,height=10)
print(P15)
dev.off()
ggsave(filename = "mk.pdf",P13/P14/P15,height = 15,width = 10)

# 
unique(data@meta.data[["seurat_clusters"]])
manual_celltype<-read.csv("marker_celltype.csv")

new.cluster.ids <- manual_celltype$celltype
names(new.cluster.ids) <- levels(data)
seu.obj <- RenameIdents(data, new.cluster.ids)

p16 <- DimPlot(seu.obj, 
               reduction = "umap", 
               label = TRUE, 
               pt.size = 0.5) + NoLegend()
p16
p17 <- DimPlot(seu.obj, 
               reduction = "tsne", 
               label = TRUE, 
               pt.size = 0.5) + NoLegend()
p17

pdf(paste(Sample,"_marker_umap.pdf",sep=""),width=20,height=20)
print(p16)
dev.off()

pdf(paste(Sample,"_marker_tsnet.pdf",sep=""),width=10,height=10)
print(p17)
dev.off()
save(seu.obj,file = paste0(Sample,".seu.obj.Rdata"))

#KIRP####
rm(list=ls())
library(Seurat)
library(cowplot)
library(dplyr)
# library(glmGamPoi)
library(data.table)
setwd("G:/pan_result/gpu_result_tcga/KIRP/scRNA_new")
#  R.utils 
library(R.utils)
load("H:\\\\database\\SCDTD\\data\\chenlong\\Kidney Chromophobe\\disease\\HP\\GSE152938_geneplotneed.Rdata")
data<-data.combined.new

##---harmony---#
library(harmony)
data=RunHarmony(data,group.by="orig.ident",assay.use="RNA",max.inter.harmony=20)



###--cluster---#
#rm(list=ls())
#setwd("E:/test_death/test1/01.analysis221126_v2_2/01.Seurat.v2/merge")
library(Seurat)
library(ggthemes)
library(dplyr)
library(ggplot2)
library(clustree)

pc.num=1:50
data=RunTSNE(data, reduction="harmony", dims=pc.num) %>% 
  RunUMAP(reduction="harmony", dims=pc.num) %>%
  FindNeighbors(reduction="harmony", dims=pc.num)

data= FindClusters(data,resolution=0.6)


#length(levels(Idents(data)))
#table(data@active.ident) 
dir.create("harmony_0.4_20")
setwd("harmony_0.4_20")

sample_col=c("#59A14F","#4E79A7","#F28E2B","#E15759","#4E80AB",
             "#FFCC33","#C8BF2C","#8c310a","#228B22","#815c94",
             "#3692a8" ,"#c95043" ,"#b5574d","#548faf","#e97e33",
             "#76a695",'#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', 
             '#80B1D3', '#FDB462')

#sample_col=c("#BC3C29FF","#E18727FF","#20854EFF", "#7876B1FF" ,"#6F99ADFF","#FB8072", "#EDC948","#B2DF8A",  "#BEBADA",  "#BEBADA","#4E79A7", "#BAB0AC","#59A14F","#4E79A7")
##color
p=DimPlot(data,reduction = "umap",label=T,raster = FALSE)
p
ggsave("UMAP_cluster_harmony.pdf",p,width=9,height=9)

p3=DimPlot(data,reduction = "tsne",label=T,raster = FALSE)
p3
ggsave("TSNE_cluster_harmony.pdf",p3,width=9,height=9)
saveRDS(data,"Seurat.rds")

#--------
markergene<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")

markergene<-markergene$marker
library(stringr)
library(ggplot2)
# genes_to_check = unique(intersect(rownames(data),markergene))
genes_to_check = markergene
dir.create("shoudong")
setwd("shoudong")

P13 <- DotPlot(data, features = genes_to_check,
               assay='RNA' )  + coord_flip()
P13
P14 <- VlnPlot(object = data, features =genes_to_check,log =T )
P14
P15 <- FeaturePlot(object = data, features=genes_to_check )
P15

Sample="GSE152938"

pdf(paste(Sample,"_marker_DotPlot.pdf",sep=""),width=10,height=10)
print(P13)
dev.off()

pdf(paste(Sample,"_marker_VlnPlot.pdf",sep=""),width=20,height=20)
print(P14)
dev.off()

pdf(paste(Sample,"_marker_FeaturePlot.pdf",sep=""),width=10,height=10)
print(P15)
dev.off()
ggsave(filename = "mk.pdf",P13/P14/P15,height = 15,width = 10)

# 
unique(data@meta.data[["seurat_clusters"]])
manual_celltype<-read.csv("marker_celltype.csv")

new.cluster.ids <- manual_celltype$celltype
names(new.cluster.ids) <- levels(data)
seu.obj <- RenameIdents(data, new.cluster.ids)

p16 <- DimPlot(seu.obj, 
               reduction = "umap", 
               label = TRUE, 
               pt.size = 0.5) + NoLegend()
p16
p17 <- DimPlot(seu.obj, 
               reduction = "tsne", 
               label = TRUE, 
               pt.size = 0.5) + NoLegend()
p17

pdf(paste(Sample,"_marker_umap.pdf",sep=""),width=20,height=20)
print(p16)
dev.off()

pdf(paste(Sample,"_marker_tsnet.pdf",sep=""),width=10,height=10)
print(p17)
dev.off()
save(seu.obj,file = paste0(Sample,".seu.obj.Rdata"))


#LIHC####
rm(list=ls())
library(Seurat)
library(cowplot)
library(dplyr)
# library(glmGamPoi)
library(data.table)
setwd("G:/pan_result/gpu_result_tcga/LIHC/scRNA_new/GSE112271_RAW")
#  R.utils 
library(R.utils)

files <- c(list.files("G:/pan_result/gpu_result_tcga/LIHC/scRNA_new/GSE112271_RAW/"))

# 
prefixes <- unique(sub("_[^.]+\\..+$", "", files))

base_path<-"G:/pan_result/gpu_result_tcga/LIHC/scRNA_new/GSE112271_RAW"

# # 
# for (prefix in prefixes) {
#   # 
#   target_path <- file.path(base_path, prefix)
#   if (!dir.exists(target_path)) {
#     dir.create(target_path)
#   }
#   
#   # 
#   current_files <- list.files(base_path, pattern = paste0(prefix, "_"), full.names = TRUE)
#   
#   # 
#   for (file in current_files) {
#     # 
#     gunzip(file, remove = FALSE, overwrite = TRUE)
#     
#     # 
#     unzipped_file <- sub("^(.*)\\.gz$", "\\1", file)
#     
#     # 
#     if (grepl("barcodes", unzipped_file)) {
#       new_name <- "barcodes.tsv"
#     } else if (grepl("features", unzipped_file)) {
#       new_name <- "features.tsv"
#     } else if (grepl("matrix", unzipped_file)) {
#       new_name <- "matrix.mtx"
#     }
#     
#     # 
#     new_file_path <- file.path(target_path, new_name)
#     
#     # 
#     file.rename(unzipped_file, new_file_path)
#   }
# }


prefix<-prefixes[1]
file<-current_files[1]
# 
for (prefix in prefixes) {
  # 
  target_path <- file.path(base_path, prefix)
  if (!dir.exists(target_path)) {
    dir.create(target_path)
  }
  
  # 
  current_files <- list.files(base_path, pattern = paste0(prefix, "_"), full.names = TRUE)
  
  # 
  for (file in current_files) {
    
    
    # 
    # unzipped_file <- sub("^(.*)\\.gz$", "\\1", file)
    unzipped_file<-file
    # 
    if (grepl("barcodes", unzipped_file)) {
      new_name <- "barcodes.tsv.gz"
    } else if (grepl("genes", unzipped_file)) {
      new_name <- "features.tsv.gz"
    } else if (grepl("matrix", unzipped_file)) {
      new_name <- "matrix.mtx.gz"
    }
    
    # 
    new_file_path <- file.path(target_path, new_name)
    
    # 
    file.rename(unzipped_file, new_file_path)
  }
}




library(Seurat)
library(cowplot)
library(dplyr)
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(data.table)
setwd("G:/pan_result/gpu_result_tcga/LIHC/scRNA_new")
tmp = list.dirs('GSE112271_RAW/')[-1]
tmp
sce =CreateSeuratObject(counts =  Read10X(tmp) ,
                        #project = gsub('.gz','',gsub('^GSM[0-9]*_','') ) ,
                        min.cells = 5,
                        min.features = 300 )
GetAssay(sce,assay = "RNA")#17524   cells  22722  features 2
data<-sce

########
rownames(data)[grepl('^mt-',rownames(data),ignore.case = T)]
rownames(data)[grepl('^Rp[sl]',rownames(data),ignore.case = T)]
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
data <- subset(data, subset = nFeature_RNA > 200 & nCount_RNA > 2500 & percent.mt < 20)

##---NormalizeData--#
data <- NormalizeData(data, normalization.method =  "LogNormalize", 
                      scale.factor = 10000)
data=FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)

#all.genes <- rownames(data)
data <- ScaleData(data)

##---PCA--#

data=RunPCA(data, features = VariableFeatures(data),dims = 1:30)

##---harmony---#
library(harmony)
data=RunHarmony(data,group.by="orig.ident",assay.use="RNA",max.inter.harmony=20)



###--cluster---#
#rm(list=ls())
#setwd("E:/test_death/test1/01.analysis221126_v2_2/01.Seurat.v2/merge")
library(Seurat)
library(ggthemes)
library(dplyr)
library(ggplot2)
library(clustree)

pc.num=1:50
data=RunTSNE(data, reduction="harmony", dims=pc.num) %>% 
  RunUMAP(reduction="harmony", dims=pc.num) %>%
  FindNeighbors(reduction="harmony", dims=pc.num)

data= FindClusters(data,resolution=0.6)


#length(levels(Idents(data)))
#table(data@active.ident) 
dir.create("harmony_0.4_20")
setwd("harmony_0.4_20")

sample_col=c("#59A14F","#4E79A7","#F28E2B","#E15759","#4E80AB",
             "#FFCC33","#C8BF2C","#8c310a","#228B22","#815c94",
             "#3692a8" ,"#c95043" ,"#b5574d","#548faf","#e97e33",
             "#76a695",'#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', 
             '#80B1D3', '#FDB462')

#sample_col=c("#BC3C29FF","#E18727FF","#20854EFF", "#7876B1FF" ,"#6F99ADFF","#FB8072", "#EDC948","#B2DF8A",  "#BEBADA",  "#BEBADA","#4E79A7", "#BAB0AC","#59A14F","#4E79A7")
##color
p=DimPlot(data,reduction = "umap",label=T,raster = FALSE)
p
ggsave("UMAP_cluster_harmony.pdf",p,width=9,height=9)

p3=DimPlot(data,reduction = "tsne",label=T,raster = FALSE)
p3
ggsave("TSNE_cluster_harmony.pdf",p3,width=9,height=9)
saveRDS(data,"Seurat.rds")
##---cluster marker---#
clus=data.frame(data@active.ident)###cluster 
colnames(clus)=c("Cluster")
dir.create("Markers")
write.csv(clus,'Markers/Cluster.csv',quote=F)


#######---------- ---------

markergene<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")

markergene<-markergene$marker
library(stringr)
library(ggplot2)
# genes_to_check = unique(intersect(rownames(data),markergene))
genes_to_check = markergene
dir.create("shoudong")
setwd("shoudong")

P13 <- DotPlot(data, features = genes_to_check,
               assay='RNA' )  + coord_flip()
P13
P14 <- VlnPlot(object = data, features =genes_to_check,log =T )
P14
P15 <- FeaturePlot(object = data, features=genes_to_check )
P15
Sample="GSE112271"

pdf(paste(Sample,"_marker_DotPlot.pdf",sep=""),width=10,height=10)
print(P13)
dev.off()

pdf(paste(Sample,"_marker_VlnPlot.pdf",sep=""),width=20,height=20)
print(P14)
dev.off()

pdf(paste(Sample,"_marker_FeaturePlot.pdf",sep=""),width=10,height=10)
print(P15)
dev.off()
ggsave(filename = "mk.pdf",P13/P14/P15,height = 15,width = 10)

# 
unique(data@meta.data[["seurat_clusters"]])
manual_celltype<-read.csv("marker_celltype.csv")

new.cluster.ids <- manual_celltype$celltype
names(new.cluster.ids) <- levels(data)
seu.obj <- RenameIdents(data, new.cluster.ids)

p16 <- DimPlot(seu.obj, 
               reduction = "umap", 
               label = TRUE, 
               pt.size = 0.5) + NoLegend()
p16
p17 <- DimPlot(seu.obj, 
               reduction = "tsne", 
               label = TRUE, 
               pt.size = 0.5) + NoLegend()
p17

pdf(paste(Sample,"_marker_umap.pdf",sep=""),width=20,height=20)
print(p16)
dev.off()

pdf(paste(Sample,"_marker_tsnet.pdf",sep=""),width=10,height=10)
print(p17)
dev.off()
save(seu.obj,file = paste0(Sample,".seu.obj.Rdata"))


#LUNG####
RDS<-readRDS("G:\\pan_result\\gpu_result_tcga\\LUNG\\scRNA_new\\GSE131907_Lung_Cancer_normalized_log2TPM_matrix.rds")
colnames(RDS)[1:5]
rownames(RDS)[1:5]
RDS <- RDS[, grep("LUNG_T", colnames(RDS))]

setwd("G:/pan_result/gpu_result_tcga/LUNG/scRNA_new")
# write.table(RDS,"RDS.txt")
meta<- RDS
seurat_data<- fread("GSE131907_Lung_Cancer_cell_annotation.txt",data.table = F)
# CreateSeuratObject()Seurat
data <- CreateSeuratObject(counts = RDS,
                           min.features = 200,
                           min.cells = 3, 
                           project = "GSE131907")
data@assays[["RNA"]]@layers[["counts"]]@Dimnames[[1]]<-rownames(data)
data@assays[["RNA"]]@layers[["counts"]]@Dimnames[[2]]<-colnames(data)

GetAssay(data,assay = "RNA")#  cells  features 

########
rownames(data)[grepl('^mt-',rownames(data),ignore.case = T)]
rownames(data)[grepl('^Rp[sl]',rownames(data),ignore.case = T)]
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
data <- subset(data, subset = nFeature_RNA > 200 & nCount_RNA > 2500 & percent.mt < 20)
##---NormalizeData--#
data <- NormalizeData(data, normalization.method =  "LogNormalize", 
                      scale.factor = 10000)
data=FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)

#all.genes <- rownames(data)
data <- ScaleData(data)

##---PCA--#

data=RunPCA(data, features = VariableFeatures(data),dims = 1:30)




###--cluster---#
#rm(list=ls())
#setwd("E:/test_death/test1/01.analysis221126_v2_2/01.Seurat.v2/merge")
library(Seurat)
library(ggthemes)
library(dplyr)
library(ggplot2)
library(clustree)

pc.num=1:50
data=RunTSNE(data, reduction="pca", dims=pc.num) %>% 
  RunUMAP(reduction="pca", dims=pc.num) %>%
  FindNeighbors(reduction="pca", dims=pc.num)

data= FindClusters(data,resolution=0.6)


#length(levels(Idents(data)))
#table(data@active.ident) 
dir.create("harmony_0.4_20")
setwd("harmony_0.4_20")

sample_col=c("#59A14F","#4E79A7","#F28E2B","#E15759","#4E80AB",
             "#FFCC33","#C8BF2C","#8c310a","#228B22","#815c94",
             "#3692a8" ,"#c95043" ,"#b5574d","#548faf","#e97e33",
             "#76a695",'#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', 
             '#80B1D3', '#FDB462')

#sample_col=c("#BC3C29FF","#E18727FF","#20854EFF", "#7876B1FF" ,"#6F99ADFF","#FB8072", "#EDC948","#B2DF8A",  "#BEBADA",  "#BEBADA","#4E79A7", "#BAB0AC","#59A14F","#4E79A7")
##color
p=DimPlot(data,reduction = "umap",label=T,raster = FALSE)
p
p1=DimPlot(data,reduction = "umap",group.by="orig.ident",split.by="orig.ident",ncol=3,cols=sample_col,raster = FALSE)+NoLegend()
p1
ggsave("UMAP_Samples_harmony.pdf",p1,width=9,height=9)
ggsave("UMAP_cluster_harmony.pdf",p,width=9,height=9)
p2=DimPlot(data,reduction = "tsne",group.by="orig.ident",split.by="orig.ident",ncol=3,cols=sample_col,raster = FALSE)+NoLegend()
p2

p3=DimPlot(data,reduction = "tsne",label=T,raster = FALSE)
p3
ggsave("TSNE_Samples_harmony.pdf",p2,width=9,height=9)
ggsave("TSNE_cluster_harmony.pdf",p3,width=9,height=9)
saveRDS(data,"Seurat.rds")
##---cluster marker---#
clus=data.frame(data@active.ident)###cluster 
colnames(clus)=c("Cluster")
dir.create("Markers")
write.csv(clus,'Markers/Cluster.csv',quote=F)


#######---------- ---------

markergene<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")

markergene<-markergene$marker
library(stringr)
library(ggplot2)
genes_to_check = unique(intersect(rownames(data),markergene))

dir.create("shoudong")
setwd("shoudong")

P13 <- DotPlot(data, features = genes_to_check,
               assay='RNA' )  + coord_flip()
P13
P14 <- VlnPlot(object = data, features =genes_to_check,log =T )
P14
P15 <- FeaturePlot(object = data, features=genes_to_check )
P15
Sample="GSE131907"

pdf(paste(Sample,"_marker_DotPlot.pdf",sep=""),width=10,height=10)
print(P13)
dev.off()

pdf(paste(Sample,"_marker_VlnPlot.pdf",sep=""),width=20,height=20)
print(P14)
dev.off()

pdf(paste(Sample,"_marker_FeaturePlot.pdf",sep=""),width=10,height=10)
print(P15)
dev.off()
ggsave(filename = "mk.pdf",P13/P14/P15,height = 15,width = 10)

# 
unique(data@meta.data[["seurat_clusters"]])
manual_celltype<-read.csv("marker_celltype.csv")

new.cluster.ids <- manual_celltype$celltype
names(new.cluster.ids) <- levels(data)
seu.obj <- RenameIdents(data, new.cluster.ids)

p16 <- DimPlot(seu.obj, 
               reduction = "umap", 
               label = TRUE, 
               pt.size = 0.5) + NoLegend()
p16
p17 <- DimPlot(seu.obj, 
               reduction = "tsne", 
               label = TRUE, 
               pt.size = 0.5) + NoLegend()
p17

pdf(paste(Sample,"_marker_umap.pdf",sep=""),width=20,height=20)
print(p16)
dev.off()

pdf(paste(Sample,"_marker_tsnet.pdf",sep=""),width=10,height=10)
print(p17)
dev.off()
save(seu.obj,file = paste0(Sample,".seu.obj.Rdata"))

#PAAD####
rm(list=ls())
library(Seurat)
library(cowplot)
library(dplyr)
# library(glmGamPoi)
library(data.table)
setwd("G:/pan_result/gpu_result_tcga/PAAD/scRNA_new/GSE155698_RAW")
#  R.utils 
library(R.utils)

files <- c(list.files("G:/pan_result/gpu_result_tcga/PAAD/scRNA_new/GSE155698_RAW/"))

# 
prefixes <- unique(sub("_[^.]+\\..+$", "", files))
yasuo<-sub(".*_(PDAC_TISSUE_\\S+).*", "\\1", files)
# sub.tar.gz
yasuo <- sub("\\.tar\\.gz$", "", yasuo)
base_path<-"G:/pan_result/gpu_result_tcga/PAAD/scRNA_new/GSE155698_RAW"

# # 
# for (prefix in prefixes) {
#   # 
#   target_path <- file.path(base_path, prefix)
#   if (!dir.exists(target_path)) {
#     dir.create(target_path)
#   }
#   
#   # 
#   current_files <- list.files(base_path, pattern = paste0(prefix, "_"), full.names = TRUE)
#   
#   # 
#   for (file in current_files) {
#     # 
#     gunzip(file, remove = FALSE, overwrite = TRUE)
#     
#     # 
#     unzipped_file <- sub("^(.*)\\.gz$", "\\1", file)
#     
#     # 
#     if (grepl("barcodes", unzipped_file)) {
#       new_name <- "barcodes.tsv"
#     } else if (grepl("features", unzipped_file)) {
#       new_name <- "features.tsv"
#     } else if (grepl("matrix", unzipped_file)) {
#       new_name <- "matrix.mtx"
#     }
#     
#     # 
#     new_file_path <- file.path(target_path, new_name)
#     
#     # 
#     file.rename(unzipped_file, new_file_path)
#   }
# }

#  tar.gz 
tar_files <- list.files(base_path, pattern = "*.tar.gz$", full.names = TRUE)
# 
tar_file<-tar_files[1]
i=1
for (i in 1:length(tar_files)) {
  tar_file<-tar_files[i]
  prefix<-prefixes[i]
  target_path <- file.path(base_path, prefix)
  if (!dir.exists(target_path)) {
    dir.create(target_path)
  }
  
  
  # GSM4710689_PDAC_TISSUE_1
  base_name <- tools::file_path_sans_ext(basename(tar_file))
  
  # 
  new_folder <- file.path(base_path, base_name)
  if (!dir.exists(new_folder)) {
    dir.create(new_folder)
  }

  #  tar.gz 

  untar(tar_file, exdir = new_folder)
  
  #  filtered_feature_bc_matrix 
  filtered_dir <- file.path(new_folder, yasuo[i], "filtered_feature_bc_matrix")
  
  if (dir.exists(filtered_dir)) {
    #  filtered_feature_bc_matrix 
    compressed_files <- list.files(filtered_dir, pattern = "\\.gz$", full.names = TRUE)
    
    # 
    for (file in compressed_files) {
      #  file.rename() 
      file.rename(file, file.path(target_path, basename(file)))
      cat("Moved:", basename(file), "to", target_path, "\n")
    }
  } else {
    cat("No filtered_feature_bc_matrix found in", tar_file, "\n")
  }
  
  # 
  unlink(new_folder, recursive = TRUE)
}



prefix<-prefixes[1]
file<-current_files[1]
# 
for (prefix in prefixes) {
  # 
  target_path <- file.path(base_path, prefix)
  if (!dir.exists(target_path)) {
    dir.create(target_path)
  }
  
  # 
  current_files <- list.files(base_path, pattern = paste0(prefix, "_"), full.names = TRUE)
  
  # 
  for (file in current_files) {
    
    
    # 
    # unzipped_file <- sub("^(.*)\\.gz$", "\\1", file)
    unzipped_file<-file
    # 
    if (grepl("barcodes", unzipped_file)) {
      new_name <- "barcodes.tsv.gz"
    } else if (grepl("features", unzipped_file)) {
      new_name <- "features.tsv.gz"
    } else if (grepl("matrix", unzipped_file)) {
      new_name <- "matrix.mtx.gz"
    }
    
    # 
    new_file_path <- file.path(target_path, new_name)
    
    # 
    file.rename(unzipped_file, new_file_path)
  }
}




library(Seurat)
library(cowplot)
library(dplyr)
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(data.table)
setwd("G:/pan_result/gpu_result_tcga/PAAD/scRNA_new")
tmp = list.dirs('GSE155698_RAW/')[-1]
tmp
sce =CreateSeuratObject(counts =  Read10X(tmp) ,
                        #project = gsub('.gz','',gsub('^GSM[0-9]*_','') ) ,
                        min.cells = 5,
                        min.features = 300 )
GetAssay(sce,assay = "RNA")#17524   cells  22722  features 2
data<-sce

########
rownames(data)[grepl('^mt-',rownames(data),ignore.case = T)]
rownames(data)[grepl('^Rp[sl]',rownames(data),ignore.case = T)]
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
data <- subset(data, subset = nFeature_RNA > 200 & nCount_RNA > 2500 & percent.mt < 20)

##---NormalizeData--#
data <- NormalizeData(data, normalization.method =  "LogNormalize", 
                      scale.factor = 10000)
data=FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)

#all.genes <- rownames(data)
data <- ScaleData(data)

##---PCA--#

data=RunPCA(data, features = VariableFeatures(data),dims = 1:30)

##---harmony---#
library(harmony)
data=RunHarmony(data,group.by="orig.ident",assay.use="RNA",max.inter.harmony=20)



###--cluster---#
#rm(list=ls())
#setwd("E:/test_death/test1/01.analysis221126_v2_2/01.Seurat.v2/merge")
library(Seurat)
library(ggthemes)
library(dplyr)
library(ggplot2)
library(clustree)

pc.num=1:50
data=RunTSNE(data, reduction="harmony", dims=pc.num) %>% 
  RunUMAP(reduction="harmony", dims=pc.num) %>%
  FindNeighbors(reduction="harmony", dims=pc.num)

data= FindClusters(data,resolution=0.6)


#length(levels(Idents(data)))
#table(data@active.ident) 
dir.create("harmony_0.4_20")
setwd("harmony_0.4_20")

sample_col=c("#59A14F","#4E79A7","#F28E2B","#E15759","#4E80AB",
             "#FFCC33","#C8BF2C","#8c310a","#228B22","#815c94",
             "#3692a8" ,"#c95043" ,"#b5574d","#548faf","#e97e33",
             "#76a695",'#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', 
             '#80B1D3', '#FDB462')

#sample_col=c("#BC3C29FF","#E18727FF","#20854EFF", "#7876B1FF" ,"#6F99ADFF","#FB8072", "#EDC948","#B2DF8A",  "#BEBADA",  "#BEBADA","#4E79A7", "#BAB0AC","#59A14F","#4E79A7")
##color
p=DimPlot(data,reduction = "umap",label=T,raster = FALSE)
p
ggsave("UMAP_cluster_harmony.pdf",p,width=9,height=9)

p3=DimPlot(data,reduction = "tsne",label=T,raster = FALSE)
p3
ggsave("TSNE_cluster_harmony.pdf",p3,width=9,height=9)
saveRDS(data,"Seurat.rds")
##---cluster marker---#
clus=data.frame(data@active.ident)###cluster 
colnames(clus)=c("Cluster")
dir.create("Markers")
write.csv(clus,'Markers/Cluster.csv',quote=F)


#######---------- ---------

markergene<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")

markergene<-markergene$marker
library(stringr)
library(ggplot2)
# genes_to_check = unique(intersect(rownames(data),markergene))
genes_to_check = markergene
dir.create("shoudong")
setwd("shoudong")

P13 <- DotPlot(data, features = genes_to_check,
               assay='RNA' )  + coord_flip()
P13
P14 <- VlnPlot(object = data, features =genes_to_check,log =T )
P14
P15 <- FeaturePlot(object = data, features=genes_to_check )
P15
Sample="GSE155698"

pdf(paste(Sample,"_marker_DotPlot.pdf",sep=""),width=10,height=10)
print(P13)
dev.off()

pdf(paste(Sample,"_marker_VlnPlot.pdf",sep=""),width=20,height=20)
print(P14)
dev.off()

pdf(paste(Sample,"_marker_FeaturePlot.pdf",sep=""),width=10,height=10)
print(P15)
dev.off()
ggsave(filename = "mk.pdf",P13/P14/P15,height = 15,width = 10)

# 
unique(data@meta.data[["seurat_clusters"]])
manual_celltype<-read.csv("marker_celltype.csv")

new.cluster.ids <- manual_celltype$celltype
names(new.cluster.ids) <- levels(data)
seu.obj <- RenameIdents(data, new.cluster.ids)

p16 <- DimPlot(seu.obj, 
               reduction = "umap", 
               label = TRUE, 
               pt.size = 0.5) + NoLegend()
p16
p17 <- DimPlot(seu.obj, 
               reduction = "tsne", 
               label = TRUE, 
               pt.size = 0.5) + NoLegend()
p17

pdf(paste(Sample,"_marker_umap.pdf",sep=""),width=20,height=20)
print(p16)
dev.off()

pdf(paste(Sample,"_marker_tsnet.pdf",sep=""),width=10,height=10)
print(p17)
dev.off()
save(seu.obj,file = paste0(Sample,".seu.obj.Rdata"))

#THCA####
rm(list=ls())
library(Seurat)
library(cowplot)
library(dplyr)
# library(glmGamPoi)
library(data.table)
setwd("G:/pan_result/gpu_result_tcga/THCA/scRNA_new")
#  R.utils 
library(R.utils)
load("H:\\\\database\\SCDTD\\data\\chenlong\\Thyroid Cancer\\disease\\GSE191288_geneplotneed.Rdata")
data<-data.combined.new

##---harmony---#
library(harmony)
data=RunHarmony(data,group.by="orig.ident",assay.use="RNA",max.inter.harmony=20)



###--cluster---#
#rm(list=ls())
#setwd("E:/test_death/test1/01.analysis221126_v2_2/01.Seurat.v2/merge")
library(Seurat)
library(ggthemes)
library(dplyr)
library(ggplot2)
library(clustree)

pc.num=1:50
data=RunTSNE(data, reduction="harmony", dims=pc.num) %>% 
  RunUMAP(reduction="harmony", dims=pc.num) %>%
  FindNeighbors(reduction="harmony", dims=pc.num)

data= FindClusters(data,resolution=0.6)


#length(levels(Idents(data)))
#table(data@active.ident) 
dir.create("harmony_0.4_20")
setwd("harmony_0.4_20")

sample_col=c("#59A14F","#4E79A7","#F28E2B","#E15759","#4E80AB",
             "#FFCC33","#C8BF2C","#8c310a","#228B22","#815c94",
             "#3692a8" ,"#c95043" ,"#b5574d","#548faf","#e97e33",
             "#76a695",'#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', 
             '#80B1D3', '#FDB462')

#sample_col=c("#BC3C29FF","#E18727FF","#20854EFF", "#7876B1FF" ,"#6F99ADFF","#FB8072", "#EDC948","#B2DF8A",  "#BEBADA",  "#BEBADA","#4E79A7", "#BAB0AC","#59A14F","#4E79A7")
##color
p=DimPlot(data,reduction = "umap",label=T,raster = FALSE)
p
ggsave("UMAP_cluster_harmony.pdf",p,width=9,height=9)

p3=DimPlot(data,reduction = "tsne",label=T,raster = FALSE)
p3
ggsave("TSNE_cluster_harmony.pdf",p3,width=9,height=9)
saveRDS(data,"Seurat.rds")
#--------
setwd("G:\\pan_result\\gpu_result_tcga\\THCA\\scRNA_new\\harmony_0.4_20")
markergene<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")

markergene<-markergene$marker
library(stringr)
library(ggplot2)
# genes_to_check = unique(intersect(rownames(data),markergene))
genes_to_check = markergene
dir.create("shoudong")
setwd("shoudong")

P13 <- DotPlot(data, features = genes_to_check,
               assay='RNA' )  + coord_flip()
P13
P14 <- VlnPlot(object = data, features =genes_to_check,log =T )
P14
P15 <- FeaturePlot(object = data, features=genes_to_check )
P15

Sample="GSE191288"

pdf(paste(Sample,"_marker_DotPlot.pdf",sep=""),width=10,height=10)
print(P13)
dev.off()

pdf(paste(Sample,"_marker_VlnPlot.pdf",sep=""),width=20,height=20)
print(P14)
dev.off()

pdf(paste(Sample,"_marker_FeaturePlot.pdf",sep=""),width=10,height=10)
print(P15)
dev.off()
ggsave(filename = "mk.pdf",P13/P14/P15,height = 15,width = 10)

# 
unique(data@meta.data[["seurat_clusters"]])
manual_celltype<-read.csv("marker_celltype.csv")

new.cluster.ids <- manual_celltype$celltype
names(new.cluster.ids) <- levels(data)
seu.obj <- RenameIdents(data, new.cluster.ids)

p16 <- DimPlot(seu.obj, 
               reduction = "umap", 
               label = TRUE, 
               pt.size = 0.5) + NoLegend()
p16
p17 <- DimPlot(seu.obj, 
               reduction = "tsne", 
               label = TRUE, 
               pt.size = 0.5) + NoLegend()
p17

pdf(paste(Sample,"_marker_umap.pdf",sep=""),width=20,height=20)
print(p16)
dev.off()

pdf(paste(Sample,"_marker_tsnet.pdf",sep=""),width=10,height=10)
print(p17)
dev.off()
save(seu.obj,file = paste0(Sample,".seu.obj.Rdata"))


#UCEC####
rm(list=ls())
library(Seurat)
library(cowplot)
library(dplyr)
# library(glmGamPoi)
library(data.table)
setwd("G:/pan_result/gpu_result_tcga/UCEC/scRNA_new")
#  R.utils 
library(R.utils)
load("H:\\\\database\\SCDTD\\data\\chenlong\\Ovarian cancer\\Endometrioid Cancer\\GSE173682_geneplotneed.Rdata")
data<-data.combined.new

##---harmony---#
library(harmony)
data=RunHarmony(data,group.by="orig.ident",assay.use="RNA",max.inter.harmony=20)



###--cluster---#
#rm(list=ls())
#setwd("E:/test_death/test1/01.analysis221126_v2_2/01.Seurat.v2/merge")
library(Seurat)
library(ggthemes)
library(dplyr)
library(ggplot2)
library(clustree)

pc.num=1:50
data=RunTSNE(data, reduction="harmony", dims=pc.num) %>% 
  RunUMAP(reduction="harmony", dims=pc.num) %>%
  FindNeighbors(reduction="harmony", dims=pc.num)

data= FindClusters(data,resolution=0.6)


#length(levels(Idents(data)))
#table(data@active.ident) 
dir.create("harmony_0.4_20")
setwd("harmony_0.4_20")

sample_col=c("#59A14F","#4E79A7","#F28E2B","#E15759","#4E80AB",
             "#FFCC33","#C8BF2C","#8c310a","#228B22","#815c94",
             "#3692a8" ,"#c95043" ,"#b5574d","#548faf","#e97e33",
             "#76a695",'#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', 
             '#80B1D3', '#FDB462')

#sample_col=c("#BC3C29FF","#E18727FF","#20854EFF", "#7876B1FF" ,"#6F99ADFF","#FB8072", "#EDC948","#B2DF8A",  "#BEBADA",  "#BEBADA","#4E79A7", "#BAB0AC","#59A14F","#4E79A7")
##color
p=DimPlot(data,reduction = "umap",label=T,raster = FALSE)
p
ggsave("UMAP_cluster_harmony.pdf",p,width=9,height=9)

p3=DimPlot(data,reduction = "tsne",label=T,raster = FALSE)
p3
ggsave("TSNE_cluster_harmony.pdf",p3,width=9,height=9)
saveRDS(data,"Seurat.rds")

#--------
markergene<-read.csv("G:\\pan_result\\gpu_result_tcga\\.csv")

markergene<-markergene$marker
library(stringr)
library(ggplot2)
# genes_to_check = unique(intersect(rownames(data),markergene))
genes_to_check = markergene
dir.create("shoudong")
setwd("shoudong")

P13 <- DotPlot(data, features = genes_to_check,
               assay='RNA' )  + coord_flip()
P13
P14 <- VlnPlot(object = data, features =genes_to_check,log =T )
P14
P15 <- FeaturePlot(object = data, features=genes_to_check )
P15

Sample="GSE173682"

pdf(paste(Sample,"_marker_DotPlot.pdf",sep=""),width=10,height=10)
print(P13)
dev.off()

pdf(paste(Sample,"_marker_VlnPlot.pdf",sep=""),width=20,height=20)
print(P14)
dev.off()

pdf(paste(Sample,"_marker_FeaturePlot.pdf",sep=""),width=10,height=10)
print(P15)
dev.off()
ggsave(filename = "mk.pdf",P13/P14/P15,height = 15,width = 10)

# 
unique(data@meta.data[["seurat_clusters"]])
manual_celltype<-read.csv("marker_celltype.csv")

new.cluster.ids <- manual_celltype$celltype
names(new.cluster.ids) <- levels(data)
seu.obj <- RenameIdents(data, new.cluster.ids)

p16 <- DimPlot(seu.obj, 
               reduction = "umap", 
               label = TRUE, 
               pt.size = 0.5) + NoLegend()
p16
p17 <- DimPlot(seu.obj, 
               reduction = "tsne", 
               label = TRUE, 
               pt.size = 0.5) + NoLegend()
p17

pdf(paste(Sample,"_marker_umap.pdf",sep=""),width=20,height=20)
print(p16)
dev.off()

pdf(paste(Sample,"_marker_tsnet.pdf",sep=""),width=10,height=10)
print(p17)
dev.off()
save(seu.obj,file = paste0(Sample,".seu.obj.Rdata"))




















