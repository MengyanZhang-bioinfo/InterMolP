

 
#install.packages("dplyr")
library(clusterProfiler)  
library(stringr)  
library(AnnotationDbi)
library(org.Hs.eg.db)  
library(DOSE)
library(ggplot2)  
library(ggrepel)  
library(dplyr)
library(R.utils) 
R.utils::setOption("clusterProfiler.download.method","auto")
 
x<-gene
id_list <- bitr(x, 
                fromType="SYMBOL", 
                toType=c("ENTREZID","ENSEMBL",'SYMBOL'),
                OrgDb="org.Hs.eg.db")
 
id_list <- na.omit(id_list)
#id_list<-rownames(gini30)
 

go <- enrichGO(gene = id_list$ENTREZID, 
               OrgDb = org.Hs.eg.db, # 
               keyType = "ENTREZID", # 
               ont = "ALL", # ,BP()/CC()/MF()/ALL()
               pAdjustMethod = "BH", # P,fdr
               pvalueCutoff = 0.05,qvalueCutoff = 0.2, # p/q
               readable = T # IDsymbol
)
go.res <- data.frame(go) # GO
saveRDS(go.res,"G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\go.res.Rds")
###
gobar<-barplot(go, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
godot<-dotplot(go,showCategory = 6,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
gobar
godot
ggsave("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\BRCA_gobar_plot.pdf",gobar,width = 10,height = 10)
ggsave("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\BRCA_godot_plot.pdf",godot,width = 10,height = 10)

kegg <- enrichKEGG(gene = id_list$ENTREZID, 
                   organism = "hsa", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.2)
kegg$Description
kegg.res<-data.frame(kegg)
dotplot(kegg,orderBy = "x",showCategory =20) 
barplot(kegg,orderBy = "x",showCategory =20) 
 

#COAD




#install.packages("dplyr")
library(clusterProfiler)  
library(stringr)  
library(AnnotationDbi)
library(org.Hs.eg.db)  
library(DOSE)
library(ggplot2)  
library(ggrepel)  
library(dplyr)
library(R.utils) 
R.utils::setOption("clusterProfiler.download.method","auto")
 
x<-gene
id_list <- bitr(x, 
                fromType="SYMBOL", 
                toType=c("ENTREZID","ENSEMBL",'SYMBOL'),
                OrgDb="org.Hs.eg.db")
 
id_list <- na.omit(id_list)
#id_list<-rownames(gini30)
 

go <- enrichGO(gene = id_list$ENTREZID, 
               OrgDb = org.Hs.eg.db, # 
               keyType = "ENTREZID", # 
               ont = "ALL", # ,BP()/CC()/MF()/ALL()
               pAdjustMethod = "BH", # P,fdr
               pvalueCutoff = 0.05,qvalueCutoff = 0.2, # p/q
               readable = T # IDsymbol
)
go.res <- data.frame(go) # GO
saveRDS(go.res,"G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\go.res.Rds")
###
gobar<-barplot(go, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
godot<-dotplot(go,showCategory = 6,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
gobar
godot
ggsave("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\COAD_gobar_plot.pdf",gobar,width = 10,height = 10)
ggsave("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\COAD_godot_plot.pdf",godot,width = 10,height = 10)

kegg <- enrichKEGG(gene = id_list$ENTREZID, 
                   organism = "hsa", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.2)
kegg$Description
kegg.res<-data.frame(kegg)
dotplot(kegg,orderBy = "x",showCategory =20) 
barplot(kegg,orderBy = "x",showCategory =20) 
 
#ESCA####

#install.packages("dplyr")
library(clusterProfiler) 
library(stringr)  
library(AnnotationDbi)
library(org.Hs.eg.db)  
library(DOSE)
library(ggplot2)  
library(ggrepel)  
library(dplyr)
library(R.utils) 
R.utils::setOption("clusterProfiler.download.method","auto")
 
x<-gene
id_list <- bitr(x, 
                fromType="SYMBOL", 
                toType=c("ENTREZID","ENSEMBL",'SYMBOL'),
                OrgDb="org.Hs.eg.db")
 
id_list <- na.omit(id_list)
#id_list<-rownames(gini30)
 

go <- enrichGO(gene = id_list$ENTREZID, 
               OrgDb = org.Hs.eg.db, # 
               keyType = "ENTREZID", # 
               ont = "ALL", # ,BP()/CC()/MF()/ALL()
               pAdjustMethod = "BH", # P,fdr
               pvalueCutoff = 0.05,qvalueCutoff = 0.2, # p/q
               readable = T # IDsymbol
)
go.res <- data.frame(go) # GO
saveRDS(go.res,"G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\go.res.Rds")
###
gobar<-barplot(go, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
godot<-dotplot(go,showCategory = 6,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
gobar
godot
ggsave("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\ESCA_gobar_plot.pdf",gobar,width = 10,height = 10)
ggsave("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\ESCA_godot_plot.pdf",godot,width = 10,height = 10)

kegg <- enrichKEGG(gene = id_list$ENTREZID, 
                   organism = "hsa", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.2)
kegg$Description
kegg.res<-data.frame(kegg)
dotplot(kegg,orderBy = "x",showCategory =20) 
barplot(kegg,orderBy = "x",showCategory =20) 
 
#HNSC####


 
#install.packages("dplyr")
library(clusterProfiler)  
library(stringr) 
library(AnnotationDbi)
library(org.Hs.eg.db)  
library(DOSE)
library(ggplot2)  
library(ggrepel)  
library(dplyr)
library(R.utils) 
R.utils::setOption("clusterProfiler.download.method","auto")
 
x<-gene
id_list <- bitr(x, 
                fromType="SYMBOL", 
                toType=c("ENTREZID","ENSEMBL",'SYMBOL'),
                OrgDb="org.Hs.eg.db")
 
id_list <- na.omit(id_list)
#id_list<-rownames(gini30)
 

go <- enrichGO(gene = id_list$ENTREZID, 
               OrgDb = org.Hs.eg.db, # 
               keyType = "ENTREZID", # 
               ont = "ALL", # ,BP()/CC()/MF()/ALL()
               pAdjustMethod = "BH", # P,fdr
               pvalueCutoff = 0.05,qvalueCutoff = 0.2, # p/q
               readable = T # IDsymbol
)
go.res <- data.frame(go) # GO
saveRDS(go.res,"G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\go.res.Rds")
###
gobar<-barplot(go, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
godot<-dotplot(go,showCategory = 6,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
gobar
godot
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_gobar_plot.pdf",gobar,width = 10,height = 10)
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_godot_plot.pdf",godot,width = 10,height = 10)

kegg <- enrichKEGG(gene = id_list$ENTREZID, 
                   organism = "hsa", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.2)
kegg$Description
kegg.res<-data.frame(kegg)
dotplot(kegg,orderBy = "x",showCategory =20) 
barplot(kegg,orderBy = "x",showCategory =20) 
 
#KIRC####


 
#install.packages("dplyr")
library(clusterProfiler)  
library(stringr)  
library(AnnotationDbi)
library(org.Hs.eg.db) 
library(DOSE)
library(ggplot2)  
library(ggrepel)  
library(dplyr)
library(R.utils) 
R.utils::setOption("clusterProfiler.download.method","auto")
 
x<-gene
id_list <- bitr(x, 
                fromType="SYMBOL", 
                toType=c("ENTREZID","ENSEMBL",'SYMBOL'),
                OrgDb="org.Hs.eg.db")
 
id_list <- na.omit(id_list)
#id_list<-rownames(gini30)
 

go <- enrichGO(gene = id_list$ENTREZID, 
               OrgDb = org.Hs.eg.db, # 
               keyType = "ENTREZID", # 
               ont = "ALL", # ,BP()/CC()/MF()/ALL()
               pAdjustMethod = "BH", # P,fdr
               pvalueCutoff = 0.05,qvalueCutoff = 0.2, # p/q
               readable = T # IDsymbol
)
go.res <- data.frame(go) # GO
saveRDS(go.res,"G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\go.res.Rds")
###
gobar<-barplot(go, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
godot<-dotplot(go,showCategory = 6,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
gobar
godot
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\KIRC_gobar_plot.pdf",gobar,width = 10,height = 10)
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\KIRC_godot_plot.pdf",godot,width = 10,height = 10)

kegg <- enrichKEGG(gene = id_list$ENTREZID, 
                   organism = "hsa", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.2)
kegg$Description
kegg.res<-data.frame(kegg)
dotplot(kegg,orderBy = "x",showCategory =20) 
barplot(kegg,orderBy = "x",showCategory =20) 
 
#KIRP####


 
#install.packages("dplyr")
library(clusterProfiler)  
library(stringr)  
library(AnnotationDbi)
library(org.Hs.eg.db)  
library(DOSE)
library(ggplot2) 
library(ggrepel)  
library(dplyr)
library(R.utils) 
R.utils::setOption("clusterProfiler.download.method","auto")
 
x<-gene
id_list <- bitr(x, 
                fromType="SYMBOL", 
                toType=c("ENTREZID","ENSEMBL",'SYMBOL'),
                OrgDb="org.Hs.eg.db")
 
id_list <- na.omit(id_list)
#id_list<-rownames(gini30)
 

go <- enrichGO(gene = id_list$ENTREZID, 
               OrgDb = org.Hs.eg.db, # 
               keyType = "ENTREZID", # 
               ont = "ALL", # ,BP()/CC()/MF()/ALL()
               pAdjustMethod = "BH", # P,fdr
               pvalueCutoff = 0.05,qvalueCutoff = 0.2, # p/q
               readable = T # IDsymbol
)
go.res <- data.frame(go) # GO
saveRDS(go.res,"G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\go.res.Rds")
###
gobar<-barplot(go, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
godot<-dotplot(go,showCategory = 6,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
gobar
godot
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\KIRP_gobar_plot.pdf",gobar,width = 10,height = 10)
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\KIRP_godot_plot.pdf",godot,width = 10,height = 10)

kegg <- enrichKEGG(gene = id_list$ENTREZID, 
                   organism = "hsa", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.2)
kegg$Description
kegg.res<-data.frame(kegg)
dotplot(kegg,orderBy = "x",showCategory =20) 
barplot(kegg,orderBy = "x",showCategory =20) 
 

#LIHC####


 
#install.packages("dplyr")
library(clusterProfiler)  
library(stringr)  
library(AnnotationDbi)
library(org.Hs.eg.db)  
library(DOSE)
library(ggplot2)  
library(ggrepel) 
library(dplyr)
library(R.utils) 
R.utils::setOption("clusterProfiler.download.method","auto")
 
x<-gene
id_list <- bitr(x, 
                fromType="SYMBOL", 
                toType=c("ENTREZID","ENSEMBL",'SYMBOL'),
                OrgDb="org.Hs.eg.db")
 
id_list <- na.omit(id_list)
#id_list<-rownames(gini30)
 

go <- enrichGO(gene = id_list$ENTREZID, 
               OrgDb = org.Hs.eg.db, # 
               keyType = "ENTREZID", # 
               ont = "ALL", # ,BP()/CC()/MF()/ALL()
               pAdjustMethod = "BH", # P,fdr
               pvalueCutoff = 0.05,qvalueCutoff = 0.2, # p/q
               readable = T # IDsymbol
)
go.res <- data.frame(go) # GO
saveRDS(go.res,"G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\go.res.Rds")
###
gobar<-barplot(go, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
godot<-dotplot(go,showCategory = 6,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
gobar
godot
ggsave("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\LIHC_gobar_plot.pdf",gobar,width = 10,height = 10)
ggsave("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\LIHC_godot_plot.pdf",godot,width = 10,height = 10)

kegg <- enrichKEGG(gene = id_list$ENTREZID, 
                   organism = "hsa", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.2)
kegg$Description
kegg.res<-data.frame(kegg)
dotplot(kegg,orderBy = "x",showCategory =20) 
barplot(kegg,orderBy = "x",showCategory =20) 
 
#LUNG####


 
#install.packages("dplyr")
library(clusterProfiler)  
library(stringr)  
library(AnnotationDbi)
library(org.Hs.eg.db)  
library(DOSE)
library(ggplot2)  
library(ggrepel)  
library(dplyr)
library(R.utils)
R.utils::setOption("clusterProfiler.download.method","auto")
 
x<-gene
id_list <- bitr(x, 
                fromType="SYMBOL", 
                toType=c("ENTREZID","ENSEMBL",'SYMBOL'),
                OrgDb="org.Hs.eg.db")
 
id_list <- na.omit(id_list)
#id_list<-rownames(gini30)
 

go <- enrichGO(gene = id_list$ENTREZID, 
               OrgDb = org.Hs.eg.db, # 
               keyType = "ENTREZID", # 
               ont = "ALL", # ,BP()/CC()/MF()/ALL()
               pAdjustMethod = "BH", # P,fdr
               pvalueCutoff = 0.05,qvalueCutoff = 0.2, # p/q
               readable = T # IDsymbol
)
go.res <- data.frame(go) # GO
saveRDS(go.res,"G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\go.res.Rds")
###
gobar<-barplot(go, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
godot<-dotplot(go,showCategory = 6,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
gobar
godot
ggsave("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\NSCLC_gobar_plot.pdf",gobar,width = 10,height = 10)
ggsave("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\NSCLC_godot_plot.pdf",godot,width = 10,height = 10)

kegg <- enrichKEGG(gene = id_list$ENTREZID, 
                   organism = "hsa", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.2)
kegg$Description
kegg.res<-data.frame(kegg)
dotplot(kegg,orderBy = "x",showCategory =20) 
barplot(kegg,orderBy = "x",showCategory =20) 
 

#PAAD####


 
#install.packages("dplyr")
library(clusterProfiler)  
library(stringr)  
library(AnnotationDbi)
library(org.Hs.eg.db)  
library(DOSE)
library(ggplot2)  
library(ggrepel)  
library(dplyr)
library(R.utils) 
R.utils::setOption("clusterProfiler.download.method","auto")

x<-gene
id_list <- bitr(x, 
                fromType="SYMBOL", 
                toType=c("ENTREZID","ENSEMBL",'SYMBOL'),
                OrgDb="org.Hs.eg.db")
 
id_list <- na.omit(id_list)
#id_list<-rownames(gini30)
 

go <- enrichGO(gene = id_list$ENTREZID, 
               OrgDb = org.Hs.eg.db, # 
               keyType = "ENTREZID", # 
               ont = "ALL", # ,BP()/CC()/MF()/ALL()
               pAdjustMethod = "BH", # P,fdr
               pvalueCutoff = 0.05,qvalueCutoff = 0.2, # p/q
               readable = T # IDsymbol
)
go.res <- data.frame(go) # GO
saveRDS(go.res,"G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\go.res.Rds")
###
gobar<-barplot(go, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
godot<-dotplot(go,showCategory = 6,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
gobar
godot
ggsave("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\PAAD_gobar_plot.pdf",gobar,width = 10,height = 10)
ggsave("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\PAAD_godot_plot.pdf",godot,width = 10,height = 10)

kegg <- enrichKEGG(gene = id_list$ENTREZID, 
                   organism = "hsa", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.2)
kegg$Description
kegg.res<-data.frame(kegg)
dotplot(kegg,orderBy = "x",showCategory =20) 
barplot(kegg,orderBy = "x",showCategory =20) 
 
#THCA####


 
#install.packages("dplyr")
library(clusterProfiler)  
library(stringr)  
library(AnnotationDbi)
library(org.Hs.eg.db)  
library(DOSE)
library(ggplot2)  
library(ggrepel)  
library(dplyr)
library(R.utils) 
R.utils::setOption("clusterProfiler.download.method","auto")
 
x<-gene
id_list <- bitr(x, 
                fromType="SYMBOL", 
                toType=c("ENTREZID","ENSEMBL",'SYMBOL'),
                OrgDb="org.Hs.eg.db")
 
id_list <- na.omit(id_list)
#id_list<-rownames(gini30)


go <- enrichGO(gene = id_list$ENTREZID, 
               OrgDb = org.Hs.eg.db, # 
               keyType = "ENTREZID", # 
               ont = "ALL", # ,BP()/CC()/MF()/ALL()
               pAdjustMethod = "BH", # P,fdr
               pvalueCutoff = 0.05,qvalueCutoff = 0.2, # p/q
               readable = T # IDsymbol
)
go.res <- data.frame(go) # GO
saveRDS(go.res,"G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\go.res.Rds")
###
gobar<-barplot(go, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
godot<-dotplot(go,showCategory = 6,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
gobar
godot
ggsave("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\THCA_gobar_plot.pdf",gobar,width = 10,height = 10)
ggsave("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\THCA_godot_plot.pdf",godot,width = 10,height = 10)

kegg <- enrichKEGG(gene = id_list$ENTREZID, 
                   organism = "hsa", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.2)
kegg$Description
kegg.res<-data.frame(kegg)
dotplot(kegg,orderBy = "x",showCategory =20) 
barplot(kegg,orderBy = "x",showCategory =20) 
 

#UCEC####


 
#install.packages("dplyr")
library(clusterProfiler)  
library(stringr)  
library(AnnotationDbi)
library(org.Hs.eg.db)  
library(DOSE)
library(ggplot2)  
library(ggrepel)  
library(dplyr)
library(R.utils) 
R.utils::setOption("clusterProfiler.download.method","auto")
 
x<-gene
id_list <- bitr(x, 
                fromType="SYMBOL", 
                toType=c("ENTREZID","ENSEMBL",'SYMBOL'),
                OrgDb="org.Hs.eg.db")
 
id_list <- na.omit(id_list)
#id_list<-rownames(gini30)
 

go <- enrichGO(gene = id_list$ENTREZID,
               OrgDb = org.Hs.eg.db, # 
               keyType = "ENTREZID", # 
               ont = "ALL", # ,BP()/CC()/MF()/ALL()
               pAdjustMethod = "BH", # P,fdr
               pvalueCutoff = 0.05,qvalueCutoff = 0.2, # p/q
               readable = T # IDsymbol
)
go.res <- data.frame(go) # GO
saveRDS(go.res,"G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\go.res.Rds")
###
gobar<-barplot(go, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
godot<-dotplot(go,showCategory = 6,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
gobar
godot
ggsave("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\UCEC_gobar_plot.pdf",gobar,width = 10,height = 10)
ggsave("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\UCEC_godot_plot.pdf",godot,width = 10,height = 10)

kegg <- enrichKEGG(gene = id_list$ENTREZID, 
                   organism = "hsa", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.2)
kegg$Description
kegg.res<-data.frame(kegg)
dotplot(kegg,orderBy = "x",showCategory =20) 
barplot(kegg,orderBy = "x",showCategory =20) 
