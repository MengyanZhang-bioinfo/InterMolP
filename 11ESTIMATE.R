#BRCA####
#0.1estimate####
library(estimate)
# 1cpm
dat = log2(edgeR::cpm(expr)+1)
write.table(dat, file = "BRCA_estimate_input.txt", sep = '\t', quote = F)
exp.file = "BRCA_estimate_input.txt"
in.gct.file = "ESTIMATE_BRCA_input.gct"
#  GCT 
outputGCT(exp.file, in.gct.file)
rt <- read.table("ESTIMATE_BRCA_input.gct", 
                 skip = 2, # 2
                 header = TRUE, 
                 sep = "\t")

# filterCommonGenes()common gene
filterCommonGenes(input.f = exp.file, # 
                  output.f = "ESTIMATE_BRCA_gene.gct", # gct
                  id = "GeneSymbol" # GeneSymbolGeneSymbolEntrezID"common_genes"
) 

#2esitimate

out.score.file = "ESTIMATE_score.gct"
# 
estimateScore(in.gct.file, # 
              out.score.file, # 
              platform = "illumina") # TCGA,illumina

# ./estimated_purity_plots/
# plotPurity(scores = out.score.file)
## ESTIMATE
ESTIMATE_score = read.table(out.score.file,
                            skip = 2,
                            header = T,
                            row.names = 1)
ESTIMATE_score = as.data.frame(t(ESTIMATE_score[, 2:ncol(ESTIMATE_score)])) # 21
ESTIMATE_score[1:6, 1:3]
write.table(ESTIMATE_score, file = "BRCA_estimate_score.txt", sep = '\t', quote = F)








#1.4.1estimate####
ESTIMATE_score<-read.table("G:\\pan_result\\gpu_result_tcga\\BRCA\\BRCA_estimate_score.txt")
#3
## 
colnames(ESTIMATE_score)
ESTIMATE_score$sample<-gsub("\\.","-",substr(rownames(ESTIMATE_score),1,15))
ESTIMATE_score2<-merge(ESTIMATE_score,clusters,by="sample", all = TRUE)
ESTIMATE_score2<-na.omit(ESTIMATE_score2)
ESTIMATE_score2$group<-ESTIMATE_score2$clusters
ESTIMATE_score2$group<-as.character(ESTIMATE_score2$group)
library(ggpubr)
library(ggsci)
library(patchwork)
p1 <- ggplot(ESTIMATE_score2, aes(x = group, y = StromalScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values=c("#80D1C8","#FFD4A9"))+
  labs(x = "", y = 'Stromal Score') +
  stat_compare_means() +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75),
        legend.position = 'none' # 
  )

p2 <- ggplot(ESTIMATE_score2, aes(x = group, y = ImmuneScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values=c("#80D1C8","#FFD4A9"))+
  labs(x = "", y = 'ImmuneScore Score') +
  stat_compare_means() +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75),
        legend.position = 'none' # 
  )

p3 <- ggplot(ESTIMATE_score2, aes(x = group, y = ESTIMATEScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values=c("#80D1C8","#FFD4A9"))+
  labs(x = "", y = 'ESTIMATE Score') +
  stat_compare_means() +
  theme_bw(base_size = 16)

p1 + p2 + p3 +  plot_layout(ncol = 3)

ppp<-p1 + p2 + p3 +  plot_layout(ncol = 3)

ggsave("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\BRCA_cluster_estimate_plot.pdf",ppp,width = 10,height = 10)







#2 cluster
gene1<-exp[rownames(clusters)[which(clusters$clusters=="1")],]
gene2<-exp[rownames(clusters)[which(clusters$clusters=="2")],]
# gene1<-gene[label$X[which(label$label=="1")],]
# gene2<-gene[label$X[which(label$label=="2")],]
data<-cbind(t(gene1),t(gene2))
T_p<-c()
FC_result<-c()
for(i in 1:nrow(data)){
  T_result<-wilcox.test(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]),as.numeric(data[i,1:nrow(gene1)]),paired = F,var.equal = T)
  #####################################t############
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

DEG<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC>2|DEG_all$FC<0.5),]#1295


test2<-rbind(biaoda1,biaoda2)
test2<-as.data.frame(t(test2))
test2


library(pheatmap)
table(cluster_sample$class)
annotation_col=data.frame(type=factor(rep(c('cluster1','cluster2'),c(nrow(biaoda1),nrow(biaoda2)))))
rownames(annotation_col)=colnames(test2)
aaa<-pheatmap(test2,annotation_col = annotation_col,show_colnames = F,show_rownames = F,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\BRCA_exp_heatmap_plot.pdf",aaa,width = 10,height = 10)


test3<-log2(test2+1)#
data_s<-as.data.frame(t(apply(test2,1,scale)))
colnames(data_s)<-colnames(test2)
aaa<-pheatmap(data_s,annotation_col = annotation_col,show_colnames = F,show_rownames = F,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa


test_gene<-test2[DEG$Gene,]
aaa<-pheatmap(test_gene,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\BRCA_exp_heatmap_deg_noscale_plot.pdf",aaa,width = 10,height = 10)
aaa<-pheatmap(test_gene,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\BRCA_exp_heatmap_deg_scalerow_plot.pdf",aaa,width = 10,height = 10)
aaa<-pheatmap(test_gene,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(20))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\BRCA_exp_heatmap_deg_scalerow20_plot.pdf",aaa,width = 10,height = 10)


test_gene_100<-1000*test_gene
aaa<-pheatmap(test_gene_100,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
#
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\BRCA_exp_heatmap_deg100_noscale50_plot.pdf",aaa,width = 10,height = 10)

aaa<-pheatmap(test_gene_100,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(20))## 
#
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\BRCA_exp_heatmap_deg100_scalerow20_plot.pdf",aaa,width = 10,height = 10)

aaa<-pheatmap(test_gene_100,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
#
ggsave("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\BRCA_exp_heatmap_deg100_scalerow50_plot.pdf",aaa,width = 10,height = 10)




test_gene_abs<-abs(test_gene)
aaa<-pheatmap(test_gene_abs,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
#
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\BRCA_exp_heatmap_degabs_plot.pdf",aaa,width = 10,height = 10)

#COAD####
#0.1estimate####
library(estimate)
# 1cpm
dat = log2(edgeR::cpm(expr)+1)
write.table(dat, file = "COAD_estimate_input.txt", sep = '\t', quote = F)
exp.file = "COAD_estimate_input.txt"
in.gct.file = "ESTIMATE_COAD_input.gct"
#  GCT 
outputGCT(exp.file, in.gct.file)
rt <- read.table("ESTIMATE_COAD_input.gct", 
                 skip = 2, # 2
                 header = TRUE, 
                 sep = "\t")

# filterCommonGenes()common gene
filterCommonGenes(input.f = exp.file, # 
                  output.f = "ESTIMATE_COAD_gene.gct", # gct
                  id = "GeneSymbol" # GeneSymbolGeneSymbolEntrezID"common_genes"
) 

#2esitimate

out.score.file = "ESTIMATE_score.gct"
# 
estimateScore(in.gct.file, # 
              out.score.file, # 
              platform = "illumina") # TCGA,illumina

# ./estimated_purity_plots/
# plotPurity(scores = out.score.file)
## ESTIMATE
ESTIMATE_score = read.table(out.score.file,
                            skip = 2,
                            header = T,
                            row.names = 1)
ESTIMATE_score = as.data.frame(t(ESTIMATE_score[, 2:ncol(ESTIMATE_score)])) # 21
ESTIMATE_score[1:6, 1:3]
write.table(ESTIMATE_score, file = "COAD_estimate_score.txt", sep = '\t', quote = F)

getwd()






#1.4.1estimate####
ESTIMATE_score<-read.table("G:\\pan_result\\gpu_result_tcga\\COAD\\COAD_estimate_score.txt")
#3
## 
ESTIMATE_score$sample<-substr(rownames(ESTIMATE_score),1,15)
ESTIMATE_score$group<-as.character(clusters$clusters)

library(ggpubr)
library(ggsci)
library(patchwork)
p1 <- ggplot(ESTIMATE_score, aes(x = group, y = StromalScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8"))+
  labs(x = "", y = 'Stromal Score') +
  stat_compare_means() +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75),
        legend.position = 'none' # 
  )

p2 <- ggplot(ESTIMATE_score, aes(x = group, y = ImmuneScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8"))+
  labs(x = "", y = 'ImmuneScore Score') +
  stat_compare_means() +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75),
        legend.position = 'none' # 
  )

p3 <- ggplot(ESTIMATE_score, aes(x = group, y = ESTIMATEScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8"))+
  labs(x = "", y = 'ESTIMATE Score') +
  stat_compare_means() +
  theme_bw(base_size = 16)

p1 + p2 + p3 +  plot_layout(ncol = 3)

ppp<-p1 + p2 + p3 +  plot_layout(ncol = 3)

ggsave("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\COAD_cluster_estimate_plot.pdf",ppp,width = 10,height = 10)







#2 cluster
gene1<-exp[rownames(clusters)[which(clusters$clusters=="1")],]
gene2<-exp[rownames(clusters)[which(clusters$clusters=="2")],]
# gene1<-gene[label$X[which(label$label=="1")],]
# gene2<-gene[label$X[which(label$label=="2")],]
data<-cbind(t(gene1),t(gene2))
T_p<-c()
FC_result<-c()
for(i in 1:nrow(data)){
  T_result<-wilcox.test(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]),as.numeric(data[i,1:nrow(gene1)]),paired = F,var.equal = T)
  #####################################t############
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

DEG<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC>2|DEG_all$FC<0.5),]#1295


test2<-rbind(biaoda1,biaoda2)
test2<-as.data.frame(t(test2))
test2


library(pheatmap)
table(cluster_sample$class)
annotation_col=data.frame(type=factor(rep(c('cluster1','cluster2'),c(nrow(biaoda1),nrow(biaoda2)))))
rownames(annotation_col)=colnames(test2)
aaa<-pheatmap(test2,annotation_col = annotation_col,show_colnames = F,show_rownames = F,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\COAD_exp_heatmap_plot.pdf",aaa,width = 10,height = 10)


test3<-log2(test2+1)#
data_s<-as.data.frame(t(apply(test2,1,scale)))
colnames(data_s)<-colnames(test2)
aaa<-pheatmap(data_s,annotation_col = annotation_col,show_colnames = F,show_rownames = F,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa


test_gene<-test2[DEG$Gene,]
aaa<-pheatmap(test_gene,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\COAD_exp_heatmap_deg_noscale_plot.pdf",aaa,width = 10,height = 10)
aaa<-pheatmap(test_gene,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\COAD_exp_heatmap_deg_scalerow_plot.pdf",aaa,width = 10,height = 10)
aaa<-pheatmap(test_gene,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(20))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\COAD_exp_heatmap_deg_scalerow20_plot.pdf",aaa,width = 10,height = 10)

aaa<-pheatmap(test_gene,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(10))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\COAD_exp_heatmap_deg_scalerow20_plot.pdf",aaa,width = 10,height = 10)

test_gene_100<-1000*test_gene
aaa<-pheatmap(test_gene_100,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
#
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\COAD_exp_heatmap_deg100_noscale50_plot.pdf",aaa,width = 10,height = 10)

aaa<-pheatmap(test_gene_100,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(20))## 
#
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\COAD_exp_heatmap_deg100_scalerow20_plot.pdf",aaa,width = 10,height = 10)

aaa<-pheatmap(test_gene_100,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
#
ggsave("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\COAD_exp_heatmap_deg100_scalerow50_plot.pdf",aaa,width = 10,height = 10)




test_gene_abs<-abs(test_gene)
aaa<-pheatmap(test_gene_abs,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
#
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\COAD_exp_heatmap_degabs_plot.pdf",aaa,width = 10,height = 10)

#ESCA####
#0.1estimate####
library(estimate)
# 1cpm
dat = log2(edgeR::cpm(expr)+1)
write.table(dat, file = "ESCA_estimate_input.txt", sep = '\t', quote = F)
exp.file = "ESCA_estimate_input.txt"
in.gct.file = "ESTIMATE_ESCA_input.gct"
#  GCT 
outputGCT(exp.file, in.gct.file)
rt <- read.table("ESTIMATE_ESCA_input.gct", 
                 skip = 2, # 2
                 header = TRUE, 
                 sep = "\t")

# filterCommonGenes()common gene
filterCommonGenes(input.f = exp.file, # 
                  output.f = "ESTIMATE_ESCA_gene.gct", # gct
                  id = "GeneSymbol" # GeneSymbolGeneSymbolEntrezID"common_genes"
) 

#2esitimate

out.score.file = "ESTIMATE_score.gct"
# 
estimateScore(in.gct.file, # 
              out.score.file, # 
              platform = "illumina") # TCGA,illumina

# ./estimated_purity_plots/
# plotPurity(scores = out.score.file)
## ESTIMATE
ESTIMATE_score = read.table(out.score.file,
                            skip = 2,
                            header = T,
                            row.names = 1)
ESTIMATE_score = as.data.frame(t(ESTIMATE_score[, 2:ncol(ESTIMATE_score)])) # 21
ESTIMATE_score[1:6, 1:3]
write.table(ESTIMATE_score, file = "ESCA_estimate_score.txt", sep = '\t', quote = F)








#1.4.1estimate####
ESTIMATE_score<-read.table("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\ESCA_estimate_score.txt")
#3
## 
ESTIMATE_score$sample<-substr(rownames(ESTIMATE_score),1,15)
ESTIMATE_score$group<-as.character(clusters$clusters)

library(ggpubr)
library(ggsci)
library(patchwork)
p1 <- ggplot(ESTIMATE_score, aes(x = group, y = StromalScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8"))+
  labs(x = "", y = 'Stromal Score') +
  stat_compare_means() +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75),
        legend.position = 'none' # 
  )

p2 <- ggplot(ESTIMATE_score, aes(x = group, y = ImmuneScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8"))+
  labs(x = "", y = 'ImmuneScore Score') +
  stat_compare_means() +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75),
        legend.position = 'none' # 
  )

p3 <- ggplot(ESTIMATE_score, aes(x = group, y = ESTIMATEScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8"))+
  labs(x = "", y = 'ESTIMATE Score') +
  stat_compare_means() +
  theme_bw(base_size = 16)

p1 + p2 + p3 +  plot_layout(ncol = 3)

ppp<-p1 + p2 + p3 +  plot_layout(ncol = 3)

ggsave("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\ESCA_cluster_estimate_plot.pdf",ppp,width = 10,height = 10)







#2 cluster
gene1<-exp[rownames(clusters)[which(clusters$clusters=="1")],]
gene2<-exp[rownames(clusters)[which(clusters$clusters=="2")],]
# gene1<-gene[label$X[which(label$label=="1")],]
# gene2<-gene[label$X[which(label$label=="2")],]
data<-cbind(t(gene1),t(gene2))
T_p<-c()
FC_result<-c()
for(i in 1:nrow(data)){
  T_result<-wilcox.test(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]),as.numeric(data[i,1:nrow(gene1)]),paired = F,var.equal = T)
  #####################################t############
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

DEG<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC>2|DEG_all$FC<0.5),]#1295


test2<-rbind(biaoda1,biaoda2)
test2<-as.data.frame(t(test2))
test2


library(pheatmap)
table(cluster_sample$class)
annotation_col=data.frame(type=factor(rep(c('cluster1','cluster2'),c(nrow(biaoda1),nrow(biaoda2)))))
rownames(annotation_col)=colnames(test2)
aaa<-pheatmap(test2,annotation_col = annotation_col,show_colnames = F,show_rownames = F,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\ESCA_exp_heatmap_plot.pdf",aaa,width = 10,height = 10)

aaa<-pheatmap(test2,annotation_col = annotation_col,show_colnames = F,show_rownames = F,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\ESCA_exp_heatmap_plot.pdf",aaa,width = 10,height = 10)


test3<-log2(test2+1)#
data_s<-as.data.frame(t(apply(test2,1,scale)))
colnames(data_s)<-colnames(test2)
aaa<-pheatmap(data_s,annotation_col = annotation_col,show_colnames = F,show_rownames = F,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa


test_gene<-test2[DEG$Gene,]
aaa<-pheatmap(test_gene,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\ESCA_exp_heatmap_deg_noscale_plot.pdf",aaa,width = 10,height = 10)
aaa<-pheatmap(test_gene,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\ESCA_exp_heatmap_deg_scalerow_plot.pdf",aaa,width = 10,height = 10)
aaa<-pheatmap(test_gene,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(20))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\ESCA_exp_heatmap_deg_scalerow20_plot.pdf",aaa,width = 10,height = 10)


test_gene_100<-1000*test_gene
aaa<-pheatmap(test_gene_100,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
#
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\ESCA_exp_heatmap_deg100_noscale50_plot.pdf",aaa,width = 10,height = 10)

aaa<-pheatmap(test_gene_100,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(20))## 
#
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\ESCA_exp_heatmap_deg100_scalerow20_plot.pdf",aaa,width = 10,height = 10)

aaa<-pheatmap(test_gene_100,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
#
ggsave("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\ESCA_exp_heatmap_deg100_scalerow50_plot.pdf",aaa,width = 10,height = 10)




test_gene_abs<-abs(test_gene)
aaa<-pheatmap(test_gene_abs,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
#
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\ESCA_exp_heatmap_degabs_plot.pdf",aaa,width = 10,height = 10)


#HNSC####
#0.1estimate####
library(estimate)
# 1cpm
dat = log2(edgeR::cpm(expr)+1)
write.table(dat, file = "HNSC_estimate_input.txt", sep = '\t', quote = F)
exp.file = "HNSC_estimate_input.txt"
in.gct.file = "ESTIMATE_HNSC_input.gct"
#  GCT 
outputGCT(exp.file, in.gct.file)
rt <- read.table("ESTIMATE_HNSC_input.gct", 
                 skip = 2, # 2
                 header = TRUE, 
                 sep = "\t")

# filterCommonGenes()common gene
filterCommonGenes(input.f = exp.file, # 
                  output.f = "ESTIMATE_HNSC_gene.gct", # gct
                  id = "GeneSymbol" # GeneSymbolGeneSymbolEntrezID"common_genes"
) 

#2esitimate

out.score.file = "ESTIMATE_score.gct"
# 
estimateScore(in.gct.file, # 
              out.score.file, # 
              platform = "illumina") # TCGA,illumina

# ./estimated_purity_plots/
# plotPurity(scores = out.score.file)
## ESTIMATE
ESTIMATE_score = read.table(out.score.file,
                            skip = 2,
                            header = T,
                            row.names = 1)
ESTIMATE_score = as.data.frame(t(ESTIMATE_score[, 2:ncol(ESTIMATE_score)])) # 21
ESTIMATE_score[1:6, 1:3]
write.table(ESTIMATE_score, file = "HNSC_estimate_score.txt", sep = '\t', quote = F)








#1.4.1estimate####
ESTIMATE_score<-read.table("G:\\pan_result\\gpu_result_tcga\\HNSC\\HNSC_estimate_score.txt")
#3
## 
ESTIMATE_score$sample<-substr(rownames(ESTIMATE_score),1,15)
ESTIMATE_score$group<-as.character(clusters$clusters)

library(ggpubr)
library(ggsci)
library(patchwork)
p1 <- ggplot(ESTIMATE_score, aes(x = group, y = StromalScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8"))+
  labs(x = "", y = 'Stromal Score') +
  stat_compare_means() +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75),
        legend.position = 'none' # 
  )

p2 <- ggplot(ESTIMATE_score, aes(x = group, y = ImmuneScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8"))+
  labs(x = "", y = 'ImmuneScore Score') +
  stat_compare_means() +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75),
        legend.position = 'none' # 
  )

p3 <- ggplot(ESTIMATE_score, aes(x = group, y = ESTIMATEScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8"))+
  labs(x = "", y = 'ESTIMATE Score') +
  stat_compare_means() +
  theme_bw(base_size = 16)

p1 + p2 + p3 +  plot_layout(ncol = 3)

ppp<-p1 + p2 + p3 +  plot_layout(ncol = 3)

ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_cluster_estimate_plot.pdf",ppp,width = 10,height = 10)







#2 cluster
gene1<-exp[rownames(clusters)[which(clusters$clusters=="1")],]
gene2<-exp[rownames(clusters)[which(clusters$clusters=="2")],]
# gene1<-gene[label$X[which(label$label=="1")],]
# gene2<-gene[label$X[which(label$label=="2")],]
data<-cbind(t(gene1),t(gene2))
T_p<-c()
FC_result<-c()
for(i in 1:nrow(data)){
  T_result<-wilcox.test(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]),as.numeric(data[i,1:nrow(gene1)]),paired = F,var.equal = T)
  #####################################t############
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

DEG<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC>2|DEG_all$FC<0.5),]#1295


test2<-rbind(biaoda1,biaoda2)
test2<-as.data.frame(t(test2))
test2


library(pheatmap)
table(cluster_sample$class)
annotation_col=data.frame(type=factor(rep(c('cluster1','cluster2'),c(nrow(biaoda1),nrow(biaoda2)))))
rownames(annotation_col)=colnames(test2)
aaa<-pheatmap(test2,annotation_col = annotation_col,show_colnames = F,show_rownames = F,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_exp_heatmap_plot.pdf",aaa,width = 10,height = 10)

aaa<-pheatmap(test2,annotation_col = annotation_col,show_colnames = F,show_rownames = F,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_exp_heatmap_plot.pdf",aaa,width = 10,height = 10)


test3<-log2(test2+1)#
data_s<-as.data.frame(t(apply(test2,1,scale)))
colnames(data_s)<-colnames(test2)
aaa<-pheatmap(data_s,annotation_col = annotation_col,show_colnames = F,show_rownames = F,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa


test_gene<-test2[DEG$Gene,]
aaa<-pheatmap(test_gene,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_exp_heatmap_deg_noscale_plot.pdf",aaa,width = 10,height = 10)
aaa<-pheatmap(test_gene,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_exp_heatmap_deg_scalerow_plot.pdf",aaa,width = 10,height = 10)
aaa<-pheatmap(test_gene,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(20))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_exp_heatmap_deg_scalerow20_plot.pdf",aaa,width = 10,height = 10)


test_gene_100<-1000*test_gene
aaa<-pheatmap(test_gene_100,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
#
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_exp_heatmap_deg100_noscale50_plot.pdf",aaa,width = 10,height = 10)

aaa<-pheatmap(test_gene_100,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(20))## 
#
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_exp_heatmap_deg100_scalerow20_plot.pdf",aaa,width = 10,height = 10)

aaa<-pheatmap(test_gene_100,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
#
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_exp_heatmap_deg100_scalerow50_plot.pdf",aaa,width = 10,height = 10)




test_gene_abs<-abs(test_gene)
aaa<-pheatmap(test_gene_abs,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
#
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_exp_heatmap_degabs_plot.pdf",aaa,width = 10,height = 10)


#KIRC####
#0.1estimate####
library(estimate)
# 1cpm
dat = log2(edgeR::cpm(expr)+1)
write.table(dat, file = "KIRC_estimate_input.txt", sep = '\t', quote = F)
exp.file = "KIRC_estimate_input.txt"
in.gct.file = "ESTIMATE_KIRC_input.gct"
#  GCT 
outputGCT(exp.file, in.gct.file)
rt <- read.table("ESTIMATE_KIRC_input.gct", 
                 skip = 2, # 2
                 header = TRUE, 
                 sep = "\t")

# filterCommonGenes()common gene
filterCommonGenes(input.f = exp.file, # 
                  output.f = "ESTIMATE_KIRC_gene.gct", # gct
                  id = "GeneSymbol" # GeneSymbolGeneSymbolEntrezID"common_genes"
) 

#2esitimate

out.score.file = "ESTIMATE_score.gct"
# 
estimateScore(in.gct.file, # 
              out.score.file, # 
              platform = "illumina") # TCGA,illumina

# ./estimated_purity_plots/
# plotPurity(scores = out.score.file)
## ESTIMATE
ESTIMATE_score = read.table(out.score.file,
                            skip = 2,
                            header = T,
                            row.names = 1)
ESTIMATE_score = as.data.frame(t(ESTIMATE_score[, 2:ncol(ESTIMATE_score)])) # 21
ESTIMATE_score[1:6, 1:3]
write.table(ESTIMATE_score, file = "KIRC_estimate_score.txt", sep = '\t', quote = F)








#1.4.1estimate####
ESTIMATE_score<-read.table("G:\\pan_result\\gpu_result_tcga\\KIRC\\KIRC_estimate_score.txt")
#3
## 
ESTIMATE_score$sample<-substr(rownames(ESTIMATE_score),1,15)
ESTIMATE_score$group<-as.character(clusters$clusters)

library(ggpubr)
library(ggsci)
library(patchwork)
p1 <- ggplot(ESTIMATE_score, aes(x = group, y = StromalScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+
  labs(x = "", y = 'Stromal Score') +
  stat_compare_means() +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75),
        legend.position = 'none' # 
  )

p2 <- ggplot(ESTIMATE_score, aes(x = group, y = ImmuneScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+
  labs(x = "", y = 'ImmuneScore Score') +
  stat_compare_means() +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75),
        legend.position = 'none' # 
  )

p3 <- ggplot(ESTIMATE_score, aes(x = group, y = ESTIMATEScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+
  labs(x = "", y = 'ESTIMATE Score') +
  stat_compare_means() +
  theme_bw(base_size = 16)

p1 + p2 + p3 +  plot_layout(ncol = 3)

ppp<-p1 + p2 + p3 +  plot_layout(ncol = 3)

ggsave("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\KIRC_cluster_estimate_plot.pdf",ppp,width = 10,height = 10)







#2 cluster
gene1<-exp[rownames(clusters)[which(clusters$clusters=="1")],]
gene2<-exp[rownames(clusters)[which(clusters$clusters=="2")],]
gene3<-exp[rownames(clusters)[which(clusters$clusters=="3")],]
# gene1<-gene[label$X[which(label$label=="1")],]
# gene2<-gene[label$X[which(label$label=="2")],]
data<-cbind(t(gene1),t(gene2),t(gene3))
class<-rep(c(1,2,3),c(nrow(gene1),nrow(gene2),nrow(gene3)))
T_p<-c()
FC_result<-c()
for(i in 1:nrow(data)){
  T_result<-kruskal.test(as.numeric(data[i,])~class)
  #####################################t############
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

DEG<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC>2|DEG_all$FC<0.5),]#1295

# 
# biaoda1<-exp[which(cluster_sample[which(cluster_sample$class=="1"),1]%in%yangben),]
# biaoda2<-exp[which(cluster_sample[which(cluster_sample$class=="2"),1]%in%yangben),]
# biaoda3<-exp[which(cluster_sample[which(cluster_sample$class=="3"),1]%in%yangben),]
# # biaoda1<-biaoda1[,-ncol(biaoda1)]
# biaoda2<-biaoda2[,-ncol(biaoda2)]
#biaoda3<-biaoda3[,-ncol(biaoda3)]


test2<-rbind(biaoda1,biaoda2,biaoda3)
test2<-as.data.frame(t(test2))
test2


library(pheatmap)
table(cluster_sample$class)
annotation_col=data.frame(type=factor(rep(c('cluster1','cluster2'),c(nrow(biaoda1),nrow(biaoda2)))))
annotation_col=data.frame(type=factor(rep(c('cluster1','cluster2','cluster3'),c(nrow(biaoda1),nrow(biaoda2),nrow(biaoda3)))))

rownames(annotation_col)=colnames(test2)
aaa<-pheatmap(test2,annotation_col = annotation_col,show_colnames = F,show_rownames = F,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\KIRC_exp_heatmap_plot.pdf",aaa,width = 10,height = 10)


test3<-log2(test2+1)#
data_s<-as.data.frame(t(apply(test2,1,scale)))
colnames(data_s)<-colnames(test2)
aaa<-pheatmap(data_s,annotation_col = annotation_col,show_colnames = F,show_rownames = F,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa


test_gene<-test2[DEG$Gene,]
aaa<-pheatmap(test_gene,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\KIRC_exp_heatmap_deg_noscale_plot.pdf",aaa,width = 10,height = 10)
aaa<-pheatmap(test_gene,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\KIRC_exp_heatmap_deg_scalerow_plot.pdf",aaa,width = 10,height = 10)
aaa<-pheatmap(test_gene,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(20))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\KIRC_exp_heatmap_deg_scalerow20_plot.pdf",aaa,width = 10,height = 10)


test_gene_100<-1000*test_gene
aaa<-pheatmap(test_gene_100,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
#
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\KIRC_exp_heatmap_deg100_noscale50_plot.pdf",aaa,width = 10,height = 10)

aaa<-pheatmap(test_gene_100,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(20))## 
#
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\KIRC_exp_heatmap_deg100_scalerow20_plot.pdf",aaa,width = 10,height = 10)

aaa<-pheatmap(test_gene_100,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
#
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\KIRC_exp_heatmap_deg100_scalerow50_plot.pdf",aaa,width = 10,height = 10)




test_gene_abs<-abs(test_gene)
aaa<-pheatmap(test_gene_abs,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
#
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\KIRC_exp_heatmap_degabs_plot.pdf",aaa,width = 10,height = 10)



#KIRP####
#0.1estimate####
library(estimate)
# 1cpm
dat = log2(edgeR::cpm(expr)+1)
write.table(dat, file = "KIRP_estimate_input.txt", sep = '\t', quote = F)
exp.file = "KIRP_estimate_input.txt"
in.gct.file = "ESTIMATE_KIRP_input.gct"
#  GCT 
outputGCT(exp.file, in.gct.file)
rt <- read.table("ESTIMATE_KIRP_input.gct", 
                 skip = 2, # 2
                 header = TRUE, 
                 sep = "\t")

# filterCommonGenes()common gene
filterCommonGenes(input.f = exp.file, # 
                  output.f = "ESTIMATE_KIRP_gene.gct", # gct
                  id = "GeneSymbol" # GeneSymbolGeneSymbolEntrezID"common_genes"
) 

#2esitimate

out.score.file = "ESTIMATE_score.gct"
# 
estimateScore(in.gct.file, # 
              out.score.file, # 
              platform = "illumina") # TCGA,illumina

# ./estimated_purity_plots/
# plotPurity(scores = out.score.file)
## ESTIMATE
ESTIMATE_score = read.table(out.score.file,
                            skip = 2,
                            header = T,
                            row.names = 1)
ESTIMATE_score = as.data.frame(t(ESTIMATE_score[, 2:ncol(ESTIMATE_score)])) # 21
ESTIMATE_score[1:6, 1:3]
write.table(ESTIMATE_score, file = "KIRP_estimate_score.txt", sep = '\t', quote = F)








#1.4.1estimate####
ESTIMATE_score<-read.table("G:\\pan_result\\gpu_result_tcga\\KIRP\\KIRP_estimate_score.txt")
#3
## 
ESTIMATE_score$sample<-substr(rownames(ESTIMATE_score),1,15)
ESTIMATE_score$group<-as.character(clusters$clusters)

library(ggpubr)
library(ggsci)
library(patchwork)
p1 <- ggplot(ESTIMATE_score, aes(x = group, y = StromalScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+
  labs(x = "", y = 'Stromal Score') +
  stat_compare_means() +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75),
        legend.position = 'none' # 
  )

p2 <- ggplot(ESTIMATE_score, aes(x = group, y = ImmuneScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+
  labs(x = "", y = 'ImmuneScore Score') +
  stat_compare_means() +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75),
        legend.position = 'none' # 
  )

p3 <- ggplot(ESTIMATE_score, aes(x = group, y = ESTIMATEScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+
  labs(x = "", y = 'ESTIMATE Score') +
  stat_compare_means() +
  theme_bw(base_size = 16)

p1 + p2 + p3 +  plot_layout(ncol = 3)

ppp<-p1 + p2 + p3 +  plot_layout(ncol = 3)

ggsave("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\KIRP_cluster_estimate_plot.pdf",ppp,width = 10,height = 10)







#2 cluster
gene1<-exp[rownames(clusters)[which(clusters$clusters=="1")],]
gene2<-exp[rownames(clusters)[which(clusters$clusters=="2")],]
gene3<-exp[rownames(clusters)[which(clusters$clusters=="3")],]
# gene1<-gene[label$X[which(label$label=="1")],]
# gene2<-gene[label$X[which(label$label=="2")],]
data<-cbind(t(gene1),t(gene2),t(gene3))
class<-rep(c(1,2,3),c(nrow(gene1),nrow(gene2),nrow(gene3)))
T_p<-c()
FC_result<-c()
for(i in 1:nrow(data)){
  T_result<-kruskal.test(as.numeric(data[i,])~class)
  #####################################t############
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

DEG<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC>2|DEG_all$FC<0.5),]#1295

# 
# biaoda1<-exp[which(cluster_sample[which(cluster_sample$class=="1"),1]%in%yangben),]
# biaoda2<-exp[which(cluster_sample[which(cluster_sample$class=="2"),1]%in%yangben),]
# biaoda3<-exp[which(cluster_sample[which(cluster_sample$class=="3"),1]%in%yangben),]
# # biaoda1<-biaoda1[,-ncol(biaoda1)]
# biaoda2<-biaoda2[,-ncol(biaoda2)]
#biaoda3<-biaoda3[,-ncol(biaoda3)]


test2<-rbind(biaoda1,biaoda2,biaoda3)
test2<-as.data.frame(t(test2))
test2


library(pheatmap)
table(cluster_sample$class)
annotation_col=data.frame(type=factor(rep(c('cluster1','cluster2'),c(nrow(biaoda1),nrow(biaoda2)))))
annotation_col=data.frame(type=factor(rep(c('cluster1','cluster2','cluster3'),c(nrow(biaoda1),nrow(biaoda2),nrow(biaoda3)))))

rownames(annotation_col)=colnames(test2)
aaa<-pheatmap(test2,annotation_col = annotation_col,show_colnames = F,show_rownames = F,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\KIRP_exp_heatmap_plot.pdf",aaa,width = 10,height = 10)


test3<-log2(test2+1)#
data_s<-as.data.frame(t(apply(test2,1,scale)))
colnames(data_s)<-colnames(test2)
aaa<-pheatmap(data_s,annotation_col = annotation_col,show_colnames = F,show_rownames = F,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa


test_gene<-test2[DEG$Gene,]
aaa<-pheatmap(test_gene,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\KIRP_exp_heatmap_deg_noscale_plot.pdf",aaa,width = 10,height = 10)
aaa<-pheatmap(test_gene,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\KIRP_exp_heatmap_deg_scalerow_plot.pdf",aaa,width = 10,height = 10)
aaa<-pheatmap(test_gene,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(20))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\KIRP_exp_heatmap_deg_scalerow20_plot.pdf",aaa,width = 10,height = 10)


test_gene_100<-1000*test_gene
aaa<-pheatmap(test_gene_100,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
#
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\KIRP_exp_heatmap_deg100_noscale50_plot.pdf",aaa,width = 10,height = 10)

aaa<-pheatmap(test_gene_100,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(20))## 
#
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\KIRP_exp_heatmap_deg100_scalerow20_plot.pdf",aaa,width = 10,height = 10)

aaa<-pheatmap(test_gene_100,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
#
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\KIRP_exp_heatmap_deg100_scalerow50_plot.pdf",aaa,width = 10,height = 10)




test_gene_abs<-abs(test_gene)
aaa<-pheatmap(test_gene_abs,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
#
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\KIRP_exp_heatmap_degabs_plot.pdf",aaa,width = 10,height = 10)



#LIHC####
#0.1estimate####
library(estimate)
# 1cpm
dat = log2(edgeR::cpm(expr)+1)
write.table(dat, file = "LIHC_estimate_input.txt", sep = '\t', quote = F)
exp.file = "LIHC_estimate_input.txt"
in.gct.file = "ESTIMATE_LIHC_input.gct"
#  GCT 
outputGCT(exp.file, in.gct.file)
rt <- read.table("ESTIMATE_LIHC_input.gct", 
                 skip = 2, # 2
                 header = TRUE, 
                 sep = "\t")

# filterCommonGenes()common gene
filterCommonGenes(input.f = exp.file, # 
                  output.f = "ESTIMATE_LIHC_gene.gct", # gct
                  id = "GeneSymbol" # GeneSymbolGeneSymbolEntrezID"common_genes"
) 

#2esitimate

out.score.file = "ESTIMATE_score.gct"
# 
estimateScore(in.gct.file, # 
              out.score.file, # 
              platform = "illumina") # TCGA,illumina

# ./estimated_purity_plots/
# plotPurity(scores = out.score.file)
## ESTIMATE
ESTIMATE_score = read.table(out.score.file,
                            skip = 2,
                            header = T,
                            row.names = 1)
ESTIMATE_score = as.data.frame(t(ESTIMATE_score[, 2:ncol(ESTIMATE_score)])) # 21
ESTIMATE_score[1:6, 1:3]
write.table(ESTIMATE_score, file = "LIHC_estimate_score.txt", sep = '\t', quote = F)








#1.4.1estimate####
ESTIMATE_score<-read.table("G:\\pan_result\\gpu_result_tcga\\LIHC\\LIHC_estimate_score.txt")
#3
## 
ESTIMATE_score$sample<-substr(rownames(ESTIMATE_score),1,15)
ESTIMATE_score$group<-as.character(clusters$clusters)

library(ggpubr)
library(ggsci)
library(patchwork)
p1 <- ggplot(ESTIMATE_score, aes(x = group, y = StromalScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+
  labs(x = "", y = 'Stromal Score') +
  stat_compare_means() +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75),
        legend.position = 'none' # 
  )

p2 <- ggplot(ESTIMATE_score, aes(x = group, y = ImmuneScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+
  labs(x = "", y = 'ImmuneScore Score') +
  stat_compare_means() +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75),
        legend.position = 'none' # 
  )

p3 <- ggplot(ESTIMATE_score, aes(x = group, y = ESTIMATEScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+
  labs(x = "", y = 'ESTIMATE Score') +
  stat_compare_means() +
  theme_bw(base_size = 16)

p1 + p2 + p3 +  plot_layout(ncol = 3)

ppp<-p1 + p2 + p3 +  plot_layout(ncol = 3)

ggsave("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\LIHC_cluster_estimate_plot.pdf",ppp,width = 10,height = 10)







#2 cluster
gene1<-exp[rownames(clusters)[which(clusters$clusters=="1")],]
gene2<-exp[rownames(clusters)[which(clusters$clusters=="2")],]
gene3<-exp[rownames(clusters)[which(clusters$clusters=="3")],]
# gene1<-gene[label$X[which(label$label=="1")],]
# gene2<-gene[label$X[which(label$label=="2")],]
data<-cbind(t(gene1),t(gene2),t(gene3))
class<-rep(c(1,2,3),c(nrow(gene1),nrow(gene2),nrow(gene3)))
T_p<-c()
FC_result<-c()
for(i in 1:nrow(data)){
  T_result<-kruskal.test(as.numeric(data[i,])~class)
  #####################################t############
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

DEG<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC>2|DEG_all$FC<0.5),]#1295

# 
# biaoda1<-exp[which(cluster_sample[which(cluster_sample$class=="1"),1]%in%yangben),]
# biaoda2<-exp[which(cluster_sample[which(cluster_sample$class=="2"),1]%in%yangben),]
# biaoda3<-exp[which(cluster_sample[which(cluster_sample$class=="3"),1]%in%yangben),]
# # biaoda1<-biaoda1[,-ncol(biaoda1)]
# biaoda2<-biaoda2[,-ncol(biaoda2)]
#biaoda3<-biaoda3[,-ncol(biaoda3)]


test2<-rbind(biaoda1,biaoda2,biaoda3)
test2<-as.data.frame(t(test2))
test2


library(pheatmap)
table(cluster_sample$class)
annotation_col=data.frame(type=factor(rep(c('cluster1','cluster2'),c(nrow(biaoda1),nrow(biaoda2)))))
annotation_col=data.frame(type=factor(rep(c('cluster1','cluster2','cluster3'),c(nrow(biaoda1),nrow(biaoda2),nrow(biaoda3)))))

rownames(annotation_col)=colnames(test2)
aaa<-pheatmap(test2,annotation_col = annotation_col,show_colnames = F,show_rownames = F,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\LIHC_exp_heatmap_plot.pdf",aaa,width = 10,height = 10)


test3<-log2(test2+1)#
data_s<-as.data.frame(t(apply(test2,1,scale)))
colnames(data_s)<-colnames(test2)
aaa<-pheatmap(data_s,annotation_col = annotation_col,show_colnames = F,show_rownames = F,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa


test_gene<-test2[DEG$Gene,]
aaa<-pheatmap(test_gene,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\LIHC_exp_heatmap_deg_noscale_plot.pdf",aaa,width = 10,height = 10)
aaa<-pheatmap(test_gene,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\LIHC_exp_heatmap_deg_scalerow_plot.pdf",aaa,width = 10,height = 10)
aaa<-pheatmap(test_gene,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(20))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\LIHC_exp_heatmap_deg_scalerow20_plot.pdf",aaa,width = 10,height = 10)


test_gene_100<-1000*test_gene
aaa<-pheatmap(test_gene_100,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
#
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\LIHC_exp_heatmap_deg100_noscale50_plot.pdf",aaa,width = 10,height = 10)

aaa<-pheatmap(test_gene_100,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(20))## 
#
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\LIHC_exp_heatmap_deg100_scalerow20_plot.pdf",aaa,width = 10,height = 10)

aaa<-pheatmap(test_gene_100,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
#
ggsave("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\LIHC_exp_heatmap_deg100_scalerow50_plot.pdf",aaa,width = 10,height = 10)




test_gene_abs<-abs(test_gene)
aaa<-pheatmap(test_gene_abs,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
#
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\LIHC_exp_heatmap_degabs_plot.pdf",aaa,width = 10,height = 10)


#LUNG####
#0.1estimate####
library(estimate)
# 1cpm
dat = log2(edgeR::cpm(expr)+1)
write.table(dat, file = "LUNG_estimate_input.txt", sep = '\t', quote = F)
exp.file = "LUNG_estimate_input.txt"
in.gct.file = "ESTIMATE_LUNG_input.gct"
#  GCT 
outputGCT(exp.file, in.gct.file)
rt <- read.table("ESTIMATE_LUNG_input.gct", 
                 skip = 2, # 2
                 header = TRUE, 
                 sep = "\t")

# filterCommonGenes()common gene
filterCommonGenes(input.f = exp.file, # 
                  output.f = "ESTIMATE_LUNG_gene.gct", # gct
                  id = "GeneSymbol" # GeneSymbolGeneSymbolEntrezID"common_genes"
) 

#2esitimate

out.score.file = "ESTIMATE_score.gct"
# 
estimateScore(in.gct.file, # 
              out.score.file, # 
              platform = "illumina") # TCGA,illumina

# ./estimated_purity_plots/
# plotPurity(scores = out.score.file)
## ESTIMATE
ESTIMATE_score = read.table(out.score.file,
                            skip = 2,
                            header = T,
                            row.names = 1)
ESTIMATE_score = as.data.frame(t(ESTIMATE_score[, 2:ncol(ESTIMATE_score)])) # 21
ESTIMATE_score[1:6, 1:3]
write.table(ESTIMATE_score, file = "LUNG_estimate_score.txt", sep = '\t', quote = F)








#1.4.1estimate####
ESTIMATE_score<-read.table("G:\\pan_result\\gpu_result_tcga\\LUNG\\LUNG_estimate_score.txt")
#3
## 
ESTIMATE_score$sample<-substr(rownames(ESTIMATE_score),1,15)
ESTIMATE_score$group<-as.character(clusters$clusters)

library(ggpubr)
library(ggsci)
library(patchwork)
p1 <- ggplot(ESTIMATE_score, aes(x = group, y = StromalScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+
  labs(x = "", y = 'Stromal Score') +
  stat_compare_means() +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75),
        legend.position = 'none' # 
  )

p2 <- ggplot(ESTIMATE_score, aes(x = group, y = ImmuneScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+
  labs(x = "", y = 'ImmuneScore Score') +
  stat_compare_means() +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75),
        legend.position = 'none' # 
  )

p3 <- ggplot(ESTIMATE_score, aes(x = group, y = ESTIMATEScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+
  labs(x = "", y = 'ESTIMATE Score') +
  stat_compare_means() +
  theme_bw(base_size = 16)

p1 + p2 + p3 +  plot_layout(ncol = 3)

ppp<-p1 + p2 + p3 +  plot_layout(ncol = 3)

ggsave("G:\\pan_result\\gpu_result_tcga\\LUNG\\450k\\LUNG_cluster_estimate_plot.pdf",ppp,width = 10,height = 10)







#2 cluster
gene1<-exp[rownames(clusters)[which(clusters$clusters=="1")],]
gene2<-exp[rownames(clusters)[which(clusters$clusters=="2")],]
gene3<-exp[rownames(clusters)[which(clusters$clusters=="3")],]
# gene1<-gene[label$X[which(label$label=="1")],]
# gene2<-gene[label$X[which(label$label=="2")],]
data<-cbind(t(gene1),t(gene2),t(gene3))
class<-rep(c(1,2,3),c(nrow(gene1),nrow(gene2),nrow(gene3)))
T_p<-c()
FC_result<-c()
for(i in 1:nrow(data)){
  T_result<-kruskal.test(as.numeric(data[i,])~class)
  #####################################t############
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

DEG<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC>2|DEG_all$FC<0.5),]#1295

# 
# biaoda1<-exp[which(cluster_sample[which(cluster_sample$class=="1"),1]%in%yangben),]
# biaoda2<-exp[which(cluster_sample[which(cluster_sample$class=="2"),1]%in%yangben),]
# biaoda3<-exp[which(cluster_sample[which(cluster_sample$class=="3"),1]%in%yangben),]
# # biaoda1<-biaoda1[,-ncol(biaoda1)]
# biaoda2<-biaoda2[,-ncol(biaoda2)]
#biaoda3<-biaoda3[,-ncol(biaoda3)]


test2<-rbind(biaoda1,biaoda2,biaoda3)
test2<-as.data.frame(t(test2))
test2


library(pheatmap)
table(cluster_sample$class)
annotation_col=data.frame(type=factor(rep(c('cluster1','cluster2'),c(nrow(biaoda1),nrow(biaoda2)))))
annotation_col=data.frame(type=factor(rep(c('cluster1','cluster2','cluster3'),c(nrow(biaoda1),nrow(biaoda2),nrow(biaoda3)))))

rownames(annotation_col)=colnames(test2)
aaa<-pheatmap(test2,annotation_col = annotation_col,show_colnames = F,show_rownames = F,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\LUNG\\450k\\LUNG_exp_heatmap_plot.pdf",aaa,width = 10,height = 10)


test3<-log2(test2+1)#
data_s<-as.data.frame(t(apply(test2,1,scale)))
colnames(data_s)<-colnames(test2)
aaa<-pheatmap(data_s,annotation_col = annotation_col,show_colnames = F,show_rownames = F,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa


test_gene<-test2[DEG$Gene,]
aaa<-pheatmap(test_gene,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\LUNG\\450k\\LUNG_exp_heatmap_deg_noscale_plot.pdf",aaa,width = 10,height = 10)
aaa<-pheatmap(test_gene,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\LUNG\\450k\\LUNG_exp_heatmap_deg_scalerow_plot.pdf",aaa,width = 10,height = 10)
aaa<-pheatmap(test_gene,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(20))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\LUNG\\450k\\LUNG_exp_heatmap_deg_scalerow20_plot.pdf",aaa,width = 10,height = 10)


test_gene_100<-1000*test_gene
aaa<-pheatmap(test_gene_100,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
#
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\LUNG\\450k\\LUNG_exp_heatmap_deg100_noscale50_plot.pdf",aaa,width = 10,height = 10)

aaa<-pheatmap(test_gene_100,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(20))## 
#
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\LUNG\\450k\\LUNG_exp_heatmap_deg100_scalerow20_plot.pdf",aaa,width = 10,height = 10)

aaa<-pheatmap(test_gene_100,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
#
ggsave("G:\\pan_result\\gpu_result_tcga\\LUNG\\450k\\LUNG_exp_heatmap_deg100_scalerow50_plot.pdf",aaa,width = 10,height = 10)




test_gene_abs<-abs(test_gene)
aaa<-pheatmap(test_gene_abs,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
#
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\LUNG\\450k\\LUNG_exp_heatmap_degabs_plot.pdf",aaa,width = 10,height = 10)

#PAAD####
#0.1estimate####
library(estimate)
# 1cpm
dat = log2(edgeR::cpm(expr)+1)
write.table(dat, file = "PAAD_estimate_input.txt", sep = '\t', quote = F)
exp.file = "PAAD_estimate_input.txt"
in.gct.file = "ESTIMATE_PAAD_input.gct"
#  GCT 
outputGCT(exp.file, in.gct.file)
rt <- read.table("ESTIMATE_PAAD_input.gct", 
                 skip = 2, # 2
                 header = TRUE, 
                 sep = "\t")

# filterCommonGenes()common gene
filterCommonGenes(input.f = exp.file, # 
                  output.f = "ESTIMATE_PAAD_gene.gct", # gct
                  id = "GeneSymbol" # GeneSymbolGeneSymbolEntrezID"common_genes"
) 

#2esitimate

out.score.file = "ESTIMATE_score.gct"
# 
estimateScore(in.gct.file, # 
              out.score.file, # 
              platform = "illumina") # TCGA,illumina

# ./estimated_purity_plots/
# plotPurity(scores = out.score.file)
## ESTIMATE
ESTIMATE_score = read.table(out.score.file,
                            skip = 2,
                            header = T,
                            row.names = 1)
ESTIMATE_score = as.data.frame(t(ESTIMATE_score[, 2:ncol(ESTIMATE_score)])) # 21
ESTIMATE_score[1:6, 1:3]
write.table(ESTIMATE_score, file = "PAAD_estimate_score.txt", sep = '\t', quote = F)








#1.4.1estimate####
ESTIMATE_score<-read.table("G:\\pan_result\\gpu_result_tcga\\PAAD\\PAAD_estimate_score.txt")
#3
## 
ESTIMATE_score$sample<-substr(rownames(ESTIMATE_score),1,15)
ESTIMATE_score$group<-as.character(clusters$clusters)

library(ggpubr)
library(ggsci)
library(patchwork)
p1 <- ggplot(ESTIMATE_score, aes(x = group, y = StromalScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8"))+
  labs(x = "", y = 'Stromal Score') +
  stat_compare_means() +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75),
        legend.position = 'none' # 
  )

p2 <- ggplot(ESTIMATE_score, aes(x = group, y = ImmuneScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8"))+
  labs(x = "", y = 'ImmuneScore Score') +
  stat_compare_means() +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75),
        legend.position = 'none' # 
  )

p3 <- ggplot(ESTIMATE_score, aes(x = group, y = ESTIMATEScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8"))+
  labs(x = "", y = 'ESTIMATE Score') +
  stat_compare_means() +
  theme_bw(base_size = 16)

p1 + p2 + p3 +  plot_layout(ncol = 3)

ppp<-p1 + p2 + p3 +  plot_layout(ncol = 3)

ggsave("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\PAAD_cluster_estimate_plot.pdf",ppp,width = 10,height = 10)







#2 cluster
gene1<-exp[rownames(clusters)[which(clusters$clusters=="1")],]
gene2<-exp[rownames(clusters)[which(clusters$clusters=="2")],]
# gene1<-gene[label$X[which(label$label=="1")],]
# gene2<-gene[label$X[which(label$label=="2")],]
data<-cbind(t(gene1),t(gene2))
T_p<-c()
FC_result<-c()
for(i in 1:nrow(data)){
  T_result<-wilcox.test(as.numeric(data[i,(nrow(gene1)+1):ncol(data)]),as.numeric(data[i,1:nrow(gene1)]),paired = F,var.equal = T)
  #####################################t############
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

DEG<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC>2|DEG_all$FC<0.5),]#1295


test2<-rbind(biaoda1,biaoda2)
test2<-as.data.frame(t(test2))
test2


library(pheatmap)
table(cluster_sample$class)
annotation_col=data.frame(type=factor(rep(c('cluster1','cluster2'),c(nrow(biaoda1),nrow(biaoda2)))))
rownames(annotation_col)=colnames(test2)
aaa<-pheatmap(test2,annotation_col = annotation_col,show_colnames = F,show_rownames = F,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\PAAD_exp_heatmap_plot.pdf",aaa,width = 10,height = 10)


test3<-log2(test2+1)#
data_s<-as.data.frame(t(apply(test2,1,scale)))
colnames(data_s)<-colnames(test2)
aaa<-pheatmap(data_s,annotation_col = annotation_col,show_colnames = F,show_rownames = F,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa


test_gene<-test2[DEG$Gene,]
aaa<-pheatmap(test_gene,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\PAAD_exp_heatmap_deg_noscale_plot.pdf",aaa,width = 10,height = 10)
aaa<-pheatmap(test_gene,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\PAAD_exp_heatmap_deg_scalerow_plot.pdf",aaa,width = 10,height = 10)
aaa<-pheatmap(test_gene,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(20))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\PAAD_exp_heatmap_deg_scalerow20_plot.pdf",aaa,width = 10,height = 10)


test_gene_100<-1000*test_gene
aaa<-pheatmap(test_gene_100,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
#
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\PAAD_exp_heatmap_deg100_noscale50_plot.pdf",aaa,width = 10,height = 10)

aaa<-pheatmap(test_gene_100,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(20))## 
#
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\PAAD_exp_heatmap_deg100_scalerow20_plot.pdf",aaa,width = 10,height = 10)

aaa<-pheatmap(test_gene_100,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
#
ggsave("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\PAAD_exp_heatmap_deg100_scalerow50_plot.pdf",aaa,width = 10,height = 10)




test_gene_abs<-abs(test_gene)
aaa<-pheatmap(test_gene_abs,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
#
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\PAAD_exp_heatmap_degabs_plot.pdf",aaa,width = 10,height = 10)




#THCA####
#0.1estimate####
library(estimate)
# 1cpm
dat = log2(edgeR::cpm(expr)+1)
write.table(dat, file = "THCA_estimate_input.txt", sep = '\t', quote = F)
exp.file = "THCA_estimate_input.txt"
in.gct.file = "ESTIMATE_THCA_input.gct"
#  GCT 
outputGCT(exp.file, in.gct.file)
rt <- read.table("ESTIMATE_THCA_input.gct", 
                 skip = 2, # 2
                 header = TRUE, 
                 sep = "\t")

# filterCommonGenes()common gene
filterCommonGenes(input.f = exp.file, # 
                  output.f = "ESTIMATE_THCA_gene.gct", # gct
                  id = "GeneSymbol" # GeneSymbolGeneSymbolEntrezID"common_genes"
) 

#2esitimate

out.score.file = "ESTIMATE_score.gct"
# 
estimateScore(in.gct.file, # 
              out.score.file, # 
              platform = "illumina") # TCGA,illumina

# ./estimated_purity_plots/
# plotPurity(scores = out.score.file)
## ESTIMATE
ESTIMATE_score = read.table(out.score.file,
                            skip = 2,
                            header = T,
                            row.names = 1)
ESTIMATE_score = as.data.frame(t(ESTIMATE_score[, 2:ncol(ESTIMATE_score)])) # 21
ESTIMATE_score[1:6, 1:3]
write.table(ESTIMATE_score, file = "THCA_estimate_score.txt", sep = '\t', quote = F)








#1.4.1estimate####
ESTIMATE_score<-read.table("G:\\pan_result\\gpu_result_tcga\\THCA\\THCA_estimate_score.txt")
#3
## 
ESTIMATE_score$sample<-substr(rownames(ESTIMATE_score),1,15)
ESTIMATE_score$group<-as.character(clusters$clusters)

library(ggpubr)
library(ggsci)
library(patchwork)
p1 <- ggplot(ESTIMATE_score, aes(x = group, y = StromalScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+
  labs(x = "", y = 'Stromal Score') +
  stat_compare_means() +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75),
        legend.position = 'none' # 
  )

p2 <- ggplot(ESTIMATE_score, aes(x = group, y = ImmuneScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+
  labs(x = "", y = 'ImmuneScore Score') +
  stat_compare_means() +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75),
        legend.position = 'none' # 
  )

p3 <- ggplot(ESTIMATE_score, aes(x = group, y = ESTIMATEScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+
  labs(x = "", y = 'ESTIMATE Score') +
  stat_compare_means() +
  theme_bw(base_size = 16)

p1 + p2 + p3 +  plot_layout(ncol = 3)

ppp<-p1 + p2 + p3 +  plot_layout(ncol = 3)

ggsave("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\THCA_cluster_estimate_plot.pdf",ppp,width = 10,height = 10)







#2 cluster
gene1<-exp[rownames(clusters)[which(clusters$clusters=="1")],]
gene2<-exp[rownames(clusters)[which(clusters$clusters=="2")],]
gene3<-exp[rownames(clusters)[which(clusters$clusters=="3")],]
# gene1<-gene[label$X[which(label$label=="1")],]
# gene2<-gene[label$X[which(label$label=="2")],]
data<-cbind(t(gene1),t(gene2),t(gene3))
class<-rep(c(1,2,3),c(nrow(gene1),nrow(gene2),nrow(gene3)))
T_p<-c()
FC_result<-c()
for(i in 1:nrow(data)){
  T_result<-kruskal.test(as.numeric(data[i,])~class)
  #####################################t############
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

DEG<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC>2|DEG_all$FC<0.5),]#1295

# 
# biaoda1<-exp[which(cluster_sample[which(cluster_sample$class=="1"),1]%in%yangben),]
# biaoda2<-exp[which(cluster_sample[which(cluster_sample$class=="2"),1]%in%yangben),]
# biaoda3<-exp[which(cluster_sample[which(cluster_sample$class=="3"),1]%in%yangben),]
# # biaoda1<-biaoda1[,-ncol(biaoda1)]
# biaoda2<-biaoda2[,-ncol(biaoda2)]
#biaoda3<-biaoda3[,-ncol(biaoda3)]


test2<-rbind(biaoda1,biaoda2,biaoda3)
test2<-as.data.frame(t(test2))
test2


library(pheatmap)
table(cluster_sample$class)
annotation_col=data.frame(type=factor(rep(c('cluster1','cluster2'),c(nrow(biaoda1),nrow(biaoda2)))))
annotation_col=data.frame(type=factor(rep(c('cluster1','cluster2','cluster3'),c(nrow(biaoda1),nrow(biaoda2),nrow(biaoda3)))))

rownames(annotation_col)=colnames(test2)
aaa<-pheatmap(test2,annotation_col = annotation_col,show_colnames = F,show_rownames = F,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\THCA_exp_heatmap_plot.pdf",aaa,width = 10,height = 10)


test3<-log2(test2+1)#
data_s<-as.data.frame(t(apply(test2,1,scale)))
colnames(data_s)<-colnames(test2)
aaa<-pheatmap(data_s,annotation_col = annotation_col,show_colnames = F,show_rownames = F,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa


test_gene<-test2[DEG$Gene,]
aaa<-pheatmap(test_gene,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\THCA_exp_heatmap_deg_noscale_plot.pdf",aaa,width = 10,height = 10)
aaa<-pheatmap(test_gene,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\THCA_exp_heatmap_deg_scalerow_plot.pdf",aaa,width = 10,height = 10)
aaa<-pheatmap(test_gene,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(20))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\THCA_exp_heatmap_deg_scalerow20_plot.pdf",aaa,width = 10,height = 10)


test_gene_100<-1000*test_gene
aaa<-pheatmap(test_gene_100,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
#
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\THCA_exp_heatmap_deg100_noscale50_plot.pdf",aaa,width = 10,height = 10)

aaa<-pheatmap(test_gene_100,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(20))## 
#
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\THCA_exp_heatmap_deg100_scalerow20_plot.pdf",aaa,width = 10,height = 10)

aaa<-pheatmap(test_gene_100,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
#
ggsave("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\THCA_exp_heatmap_deg100_scalerow50_plot.pdf",aaa,width = 10,height = 10)




test_gene_abs<-abs(test_gene)
aaa<-pheatmap(test_gene_abs,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
#
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\THCA_exp_heatmap_degabs_plot.pdf",aaa,width = 10,height = 10)


#UCEC####
#0.1estimate####
library(estimate)
# 1cpm
dat = log2(edgeR::cpm(expr)+1)
write.table(dat, file = "UCEC_estimate_input.txt", sep = '\t', quote = F)
exp.file = "UCEC_estimate_input.txt"
in.gct.file = "ESTIMATE_UCEC_input.gct"
#  GCT 
outputGCT(exp.file, in.gct.file)
rt <- read.table("ESTIMATE_UCEC_input.gct", 
                 skip = 2, # 2
                 header = TRUE, 
                 sep = "\t")

# filterCommonGenes()common gene
filterCommonGenes(input.f = exp.file, # 
                  output.f = "ESTIMATE_UCEC_gene.gct", # gct
                  id = "GeneSymbol" # GeneSymbolGeneSymbolEntrezID"common_genes"
) 

#2esitimate

out.score.file = "ESTIMATE_score.gct"
# 
estimateScore(in.gct.file, # 
              out.score.file, # 
              platform = "illumina") # TCGA,illumina

# ./estimated_purity_plots/
# plotPurity(scores = out.score.file)
## ESTIMATE
ESTIMATE_score = read.table(out.score.file,
                            skip = 2,
                            header = T,
                            row.names = 1)
ESTIMATE_score = as.data.frame(t(ESTIMATE_score[, 2:ncol(ESTIMATE_score)])) # 21
ESTIMATE_score[1:6, 1:3]
write.table(ESTIMATE_score, file = "UCEC_estimate_score.txt", sep = '\t', quote = F)








#1.4.1estimate####
ESTIMATE_score<-read.table("G:\\pan_result\\gpu_result_tcga\\UCEC\\UCEC_estimate_score.txt")
#3
## 
ESTIMATE_score$sample<-substr(rownames(ESTIMATE_score),1,15)
ESTIMATE_score$group<-as.character(clusters$clusters)

library(ggpubr)
library(ggsci)
library(patchwork)
p1 <- ggplot(ESTIMATE_score, aes(x = group, y = StromalScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+
  labs(x = "", y = 'Stromal Score') +
  stat_compare_means() +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75),
        legend.position = 'none' # 
  )

p2 <- ggplot(ESTIMATE_score, aes(x = group, y = ImmuneScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+
  labs(x = "", y = 'ImmuneScore Score') +
  stat_compare_means() +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75),
        legend.position = 'none' # 
  )

p3 <- ggplot(ESTIMATE_score, aes(x = group, y = ESTIMATEScore, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values=c("#FFD4A9", "#80D1C8","pink"))+
  labs(x = "", y = 'ESTIMATE Score') +
  stat_compare_means() +
  theme_bw(base_size = 16)

p1 + p2 + p3 +  plot_layout(ncol = 3)

ppp<-p1 + p2 + p3 +  plot_layout(ncol = 3)

ggsave("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\UCEC_cluster_estimate_plot.pdf",ppp,width = 10,height = 10)







#2 cluster
gene1<-exp[rownames(clusters)[which(clusters$clusters=="1")],]
gene2<-exp[rownames(clusters)[which(clusters$clusters=="2")],]
gene3<-exp[rownames(clusters)[which(clusters$clusters=="3")],]
# gene1<-gene[label$X[which(label$label=="1")],]
# gene2<-gene[label$X[which(label$label=="2")],]
data<-cbind(t(gene1),t(gene2),t(gene3))
class<-rep(c(1,2,3),c(nrow(gene1),nrow(gene2),nrow(gene3)))
T_p<-c()
FC_result<-c()
for(i in 1:nrow(data)){
  T_result<-kruskal.test(as.numeric(data[i,])~class)
  #####################################t############
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

DEG<-DEG_all[which(DEG_all$pval < 0.05&DEG_all$FC>2|DEG_all$FC<0.5),]#1295

# 
# biaoda1<-exp[which(cluster_sample[which(cluster_sample$class=="1"),1]%in%yangben),]
# biaoda2<-exp[which(cluster_sample[which(cluster_sample$class=="2"),1]%in%yangben),]
# biaoda3<-exp[which(cluster_sample[which(cluster_sample$class=="3"),1]%in%yangben),]
# # biaoda1<-biaoda1[,-ncol(biaoda1)]
# biaoda2<-biaoda2[,-ncol(biaoda2)]
#biaoda3<-biaoda3[,-ncol(biaoda3)]


test2<-rbind(biaoda1,biaoda2,biaoda3)
test2<-as.data.frame(t(test2))
test2


library(pheatmap)
table(cluster_sample$class)
annotation_col=data.frame(type=factor(rep(c('cluster1','cluster2'),c(nrow(biaoda1),nrow(biaoda2)))))
annotation_col=data.frame(type=factor(rep(c('cluster1','cluster2','cluster3'),c(nrow(biaoda1),nrow(biaoda2),nrow(biaoda3)))))

rownames(annotation_col)=colnames(test2)
aaa<-pheatmap(test2,annotation_col = annotation_col,show_colnames = F,show_rownames = F,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\UCEC_exp_heatmap_plot.pdf",aaa,width = 10,height = 10)


test3<-log2(test2+1)#
data_s<-as.data.frame(t(apply(test2,1,scale)))
colnames(data_s)<-colnames(test2)
aaa<-pheatmap(data_s,annotation_col = annotation_col,show_colnames = F,show_rownames = F,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa


test_gene<-test2[DEG$Gene,]
aaa<-pheatmap(test_gene,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\UCEC_exp_heatmap_deg_noscale_plot.pdf",aaa,width = 10,height = 10)
aaa<-pheatmap(test_gene,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\UCEC_exp_heatmap_deg_scalerow_plot.pdf",aaa,width = 10,height = 10)
aaa<-pheatmap(test_gene,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(20))## 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\UCEC_exp_heatmap_deg_scalerow20_plot.pdf",aaa,width = 10,height = 10)


test_gene_100<-1000*test_gene
aaa<-pheatmap(test_gene_100,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
#
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\UCEC_exp_heatmap_deg100_noscale50_plot.pdf",aaa,width = 10,height = 10)

aaa<-pheatmap(test_gene_100,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(20))## 
#
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\UCEC_exp_heatmap_deg100_scalerow20_plot.pdf",aaa,width = 10,height = 10)

aaa<-pheatmap(test_gene_100,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
#
ggsave("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\UCEC_exp_heatmap_deg100_scalerow50_plot.pdf",aaa,width = 10,height = 10)




test_gene_abs<-abs(test_gene)
aaa<-pheatmap(test_gene_abs,annotation_col = annotation_col,show_colnames = F,show_rownames = T,
              cluster_col = FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))## 
#
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\UCEC_exp_heatmap_degabs_plot.pdf",aaa,width = 10,height = 10)







