
#BRCA####

library(ggplot2) #
library(dplyr) # 
#library(ggstatsplot) 
brca_tnm<-read.table("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\TCGA.BRCA.sampleMap_BRCA_clinicalMatrix",sep = "\t",header = T)
brca_tnm<-data.frame(Sample.ID=brca_tnm$sampleID,T=brca_tnm$pathologic_T,N=brca_tnm$pathologic_N,M=brca_tnm$pathologic_M)
sigmoid_450k<-data.frame(Sample.ID=rownames(exp),exp)
brca_tnm$Sample.ID<-gsub("-",".",brca_tnm$Sample.ID)
sigmoid_450k$Sample.ID<-substr(sigmoid_450k$Sample.ID,1,15)


brca_tnm_450k<-merge(brca_tnm,sigmoid_450k,by="Sample.ID")
brca_tnm_450k<-na.omit(brca_tnm_450k)
str(brca_tnm_450k)
table(brca_tnm_450k$T)
table(brca_tnm_450k$N)
table(brca_tnm_450k$M)

brca_tnm_450k$T<-gsub(".*T1.*", "T1", brca_tnm_450k$T)
brca_tnm_450k$T<-gsub(".*T2.*", "T2", brca_tnm_450k$T)
brca_tnm_450k$T<-gsub(".*T3.*", "T3", brca_tnm_450k$T)
brca_tnm_450k$T<-gsub(".*T4.*", "T4", brca_tnm_450k$T)
brca_tnm_450k$N<-gsub(".*N0.*", "N0", brca_tnm_450k$N)
brca_tnm_450k$N<-gsub(".*N1.*", "N1", brca_tnm_450k$N)
brca_tnm_450k$N<-gsub(".*N2.*", "N2", brca_tnm_450k$N)
brca_tnm_450k$N<-gsub(".*N3.*", "N3", brca_tnm_450k$N)
brca_tnm_450k<-brca_tnm_450k[-which(brca_tnm_450k$N==""),]
# brca_tnm_450k<-brca_tnm_450k[-which(brca_tnm_450k$T=="TX"),]
# brca_tnm_450k<-brca_tnm_450k[-which(brca_tnm_450k$N==""),]
brca_tnm_450k<-brca_tnm_450k[-which(brca_tnm_450k$N=="NX"),]


cluster_sample<-data.frame(Sample.ID=rownames(clusters),class=clusters$clusters)
cluster_sample$Sample.ID<-substr(cluster_sample$Sample.ID,1,15)
brca_tnm_450k
brca_tn_450k<-merge(brca_tnm_450k,cluster_sample,by="Sample.ID")
#T    
brca_t_450k<-data.frame(Tfenqi=brca_tn_450k$T,cluster=brca_tn_450k$class,brca_tn_450k[,5:(ncol(brca_tnm_450k)-1)])
brca_t_450k<-na.omit(brca_t_450k)
brca_t_450k$cluster<-as.character(brca_t_450k$cluster)
aaa<-brca_t_450k %>% 
  ggplot(aes(x = cluster, fill = Tfenqi)) + 
  geom_bar(position = position_fill()) + 
  scale_fill_brewer(palette = 'Set3') + #          
  theme_classic() + #        
  labs(y = 'Percent') + 
  coord_flip() #          
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\BRCA_cluster+T100_plot.pdf",aaa,width = 10,height = 10)



#N    
brca_n_450k<-data.frame(Nfenqi=brca_tn_450k$N,cluster=brca_tn_450k$class,brca_tn_450k[,5:(ncol(brca_tnm_450k)-1)])
brca_n_450k<-na.omit(brca_n_450k)
table(brca_n_450k$Nfenqi)
# brca_n_450k<-brca_n_450k[-which(brca_n_450k$Nfenqi=="NX"),]
brca_n_450k$cluster<-as.character(brca_n_450k$cluster)
aaa<-brca_n_450k %>% 
  ggplot(aes(x = cluster, fill = Nfenqi)) + #x          clarity，          color（J  H）
  geom_bar(position = position_fill()) + 
  scale_fill_brewer(palette = 'Set3') + #          
  theme_classic() + #        
  labs(y = 'Percent') + #    y      ‘Percent’
  coord_flip() #          
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\BRCA_cluster+N100_plot.pdf",aaa,width = 10,height = 10)




#stage#
library(tidytree)
colnames(label)[1]<-"Sample.ID"
label$Sample.ID<-substr(label$Sample.ID,1,15)
brca_stage_450<-left_join(cluster_sample,label,by='Sample.ID')
a<-brca_stage_450$label
a[which(a=="1")]<-"stage I"
a[which(a=="2"|a=="3")]<-"stage II&III"
brca_stage_450$stage<-a
colnames(brca_stage_450)[2]<-"cluster"
brca_stage_450$cluster<-as.character(brca_stage_450$cluster)
aaa<-brca_stage_450 %>% 
  ggplot(aes(x = cluster, fill = stage)) + 
  geom_bar(position = position_fill()) + 
  scale_fill_brewer(palette = 'Set3') + 
  theme_classic() + 
  labs(y = 'Percent') + 
  coord_flip() 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\BRCA_cluster+stage100_plot.pdf",aaa,width = 10,height = 10)


library(data.table)
LUADsurvival=fread(paste0("G:\\dnameth\\BRCA","\\","TCGA-","BRCA",".GDC_phenotype.tsv.gz"),header=TRUE)
LUADsurvival<-as.data.frame(LUADsurvival)
rownames(LUADsurvival)<-LUADsurvival[,1]
LUADsurvival<-LUADsurvival[,-1]
rownames(LUADsurvival)=gsub("-","\\.",rownames(LUADsurvival))
unique(LUADsurvival$clinical_stage)
unique(LUADsurvival$tumor_stage)
str(LUADsurvival$clinical_stage)

matches_stage_i <-grepl("\\bStage I(A|B)?\\b", LUADsurvival$clinical_stage)
matches_stage_ii <-grepl("\\bStage II(A|B|C)?\\b", LUADsurvival$clinical_stage)
matches_stage_iii <-grepl("\\bStage III(A|B|C)?\\b", LUADsurvival$clinical_stage)

matches_stage_i <-grepl("\\bstage i(a|b)?\\b", LUADsurvival$tumor_stage)
matches_stage_ii <-grepl("\\bstage ii(a|b|c)?\\b", LUADsurvival$tumor_stage)
matches_stage_iii <-grepl("\\bstage iii(a|b|c)?\\b", LUADsurvival$tumor_stage)
LUADsurvival$tumor_stage[matches_stage_i]<-"1"
LUADsurvival$tumor_stage[matches_stage_ii]<-"2"
LUADsurvival$tumor_stage[matches_stage_iii]<-"3"
LUADstage1=LUADsurvival[LUADsurvival$tumor_stage==1,]
LUADstage2=LUADsurvival[LUADsurvival$tumor_stage==2,]
LUADstage3=LUADsurvival[LUADsurvival$tumor_stage==3,]

labeltostage<-data.table(Sample.ID=c(rownames(LUADstage1),rownames(LUADstage2),rownames(LUADstage3)),stage=c(rep(c("Stage I","Stage II","Stage III"),c(nrow(LUADstage1),nrow(LUADstage2),nrow(LUADstage3)))))


# labeltostage<-read.csv("label2stage.csv")
library(tidytree)
colnames(label)[1]<-"Sample.ID"
colnames(labeltostage)[1]<-"Sample.ID"
labeltostage$Sample.ID<-substr(labeltostage$Sample.ID,1,15)

# brca_stage_450<-left_join(cluster_sample,label,by='Sample.ID')
brca_stage_450<-left_join(brca_stage_450,labeltostage,by='Sample.ID')
# brca_stage_450<-brca_stage_450[,-c(4,5)]
colnames(brca_stage_450)[5]<-"stage"
colnames(brca_stage_450)[2]<-"cluster"
# colnames(brca_stage_450)[3]<-"stage"
# brca_stage_450$stage<-paste("stage",brca_stage_450$stage)
brca_stage_450$cluster<-paste("Cluster",brca_stage_450$cluster)


# brca_stage_450$stage[1:162]<-"Stage I"
# brca_stage_450$stage[163:231]<-"Stage II"
# brca_stage_450$stage[232:277]<-"Stage III"
# brca_stage_450$stage[which(brca_stage_450$stage=="stage 2")]<-brca_stage_450$stage2[which(brca_stage_450$stage=="stage 2")]
# brca_stage_450$stage[which(brca_stage_450$stage=="stage 1")]<-"Stage I"
# brca_stage_450$stage[which(brca_stage_450$stage=="STAGE II")]<-"Stage II"
# brca_stage_450$stage[which(brca_stage_450$stage=="STAGE III")]<-"Stage III"
# devtools::install_github("davidsjoberg/ggsankey")
library(tidyverse)
library(ggsankey)
library(ggplot2)
library(cols4all)
head(brca_stage_450)
df <- brca_stage_450 %>%
  make_long(cluster, stage)

df$node <- factor(df$node,levels = c(rev(unique(brca_stage_450$stage)),
                                     rev(unique(brca_stage_450$cluster))))


c4a_gui()
mycol <- c4a('rainbow_wh_rd',53)
mycol <- sample(mycol,length(mycol)) 
mycol<-c("#FFD4A9", "#80D1C8","pink","#7FB974","#5072BF","#AE84BB")

p <- ggplot(df, aes(x = x, next_x = next_x,
                    node = node, next_node = next_node,
                    fill = node,
                    label = node)) +
  scale_fill_manual(values = mycol) 

p1 <- p +
  geom_sankey(flow.alpha = 0.5, 
              smooth = 8, 
              width = 0.12) + 
  geom_sankey_text(size = 3.2, 
                   color = "black") +
  theme_void() +
  theme(legend.position = 'none') 
p1
ggsave("G:\\pan_result\\gpu_result_tcga\\BRCA\\450k\\BRCA_cluster+stage_sankey_plot.pdf",p1,width = 10,height = 10)




#COAD####

 
library(ggplot2)  
library(dplyr) 
#library(ggstatsplot)  

brca_tnm<-read.table("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\TCGA.COAD.sampleMap_COAD_clinicalMatrix",sep = "\t",header = T)
brca_tnm<-data.frame(Sample.ID=brca_tnm$sampleID,T=brca_tnm$pathologic_T,N=brca_tnm$pathologic_N,M=brca_tnm$pathologic_M)
sigmoid_450k<-data.frame(Sample.ID=rownames(exp),exp)
brca_tnm$Sample.ID<-gsub("-",".",brca_tnm$Sample.ID)
sigmoid_450k$Sample.ID<-substr(sigmoid_450k$Sample.ID,1,15)


brca_tnm_450k<-merge(brca_tnm,sigmoid_450k,by="Sample.ID")
brca_tnm_450k<-na.omit(brca_tnm_450k)
str(brca_tnm_450k)
brca_tnm_450k$T<-gsub(".*T4.*", "T4", brca_tnm_450k$T)
brca_tnm_450k$N<-gsub(".*N1.*", "N1", brca_tnm_450k$N)
brca_tnm_450k$N<-gsub(".*N2.*", "N2", brca_tnm_450k$N)
brca_tnm_450k<-brca_tnm_450k[-which(brca_tnm_450k$T=="Tis"),]


cluster_sample<-data.frame(Sample.ID=rownames(clusters),class=clusters$clusters)
cluster_sample$Sample.ID<-substr(cluster_sample$Sample.ID,1,15)
brca_tnm_450k
brca_tn_450k<-merge(brca_tnm_450k,cluster_sample,by="Sample.ID")
#T 
brca_t_450k<-data.frame(Tfenqi=brca_tn_450k$T,cluster=brca_tn_450k$class,brca_tn_450k[,5:(ncol(brca_tnm_450k)-1)])
brca_t_450k<-na.omit(brca_t_450k)
brca_t_450k$cluster<-as.character(brca_t_450k$cluster)
aaa<-brca_t_450k %>% 
  ggplot(aes(x = cluster, fill = Tfenqi)) +  
  geom_bar(position = position_fill()) + 
  scale_fill_brewer(palette = 'Set3') +  
  theme_classic() +  
  labs(y = 'Percent') + 
  coord_flip()  
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\COAD_cluster+T100_plot.pdf",aaa,width = 10,height = 10)



#N 
brca_n_450k<-data.frame(Nfenqi=brca_tn_450k$N,cluster=brca_tn_450k$class,brca_tn_450k[,5:(ncol(brca_tnm_450k)-1)])
brca_n_450k<-na.omit(brca_n_450k)
table(brca_n_450k$Nfenqi)
# brca_n_450k<-brca_n_450k[-which(brca_n_450k$Nfenqi=="NX"),]
brca_n_450k$cluster<-as.character(brca_n_450k$cluster)
aaa<-brca_n_450k %>% 
  ggplot(aes(x = cluster, fill = Nfenqi)) +  
  geom_bar(position = position_fill()) + 
  scale_fill_brewer(palette = 'Set3') +  
  theme_classic() +  
  labs(y = 'Percent') +  
  coord_flip()  
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\COAD_cluster+N100_plot.pdf",aaa,width = 10,height = 10)

#
library(ggplot2)
library(dplyr)
library(stringr)
t_cluster1<-as.data.frame(table(brca_t_450k[which(brca_t_450k$cluster==1),1]))
t_cluster1$Var1<-paste0("cluster1_",t_cluster1$Var1)
t_cluster1$Proportion<-t_cluster1$Freq/sum(t_cluster1$Freq)
t_cluster2<-as.data.frame(table(brca_t_450k[which(brca_t_450k$cluster==2),1]))
t_cluster2$Var1<-paste0("cluster2_",t_cluster2$Var1)
t_cluster2$Proportion<-t_cluster2$Freq/sum(t_cluster2$Freq)
t_cluster3<-as.data.frame(table(brca_t_450k[which(brca_t_450k$cluster==3),1]))
t_cluster3$Var1<-paste0("cluster3_",t_cluster3$Var1)
t_cluster3$Proportion<-t_cluster3$Freq/sum(t_cluster3$Freq)




t_cluster<-rbind(t_cluster1,t_cluster2,t_cluster3)
str(t_cluster)
p <- ggplot(t_cluster) + 
  geom_col(
    aes(
      x = reorder(Var1, Proportion),
      y = Proportion,
      fill = Proportion
    ),
    position = "dodge2",
    show.legend = TRUE,
    alpha = .9
  )
p
p <- p + coord_polar()
ggsave("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\COAD_cluster+T_flower_plot.pdf",p,width = 10,height = 10)

# p+geom_hline(
#   aes(yintercept = y), 
#   data.frame(y = c(0:3) * 1000),
#   color = "lightgrey"
# )+
#   geom_segment(data = t_cluster,  #  
#     aes(
#       x = reorder(str_wrap(region, 5), sum_length),
#       y = 0,
#       xend = reorder(str_wrap(region, 5), sum_length),
#       yend = 3000
#     ),
#     linetype = "dashed",
#     color = "gray12"
#   ) +
#   scale_y_continuous(
#     limits = c(-1500, 3500),
#     expand = c(0, 0),
#     breaks = c(0, 1000, 2000, 3000)
#   ) + 
#   scale_fill_gradientn(
#     "Amount of Tracks",
#     colours = c( "#6C5B7B","#C06C84","#F67280","#F8B195")
#   ) +
#   guides(
#     fill = guide_colorsteps(
#       barwidth = 15, barheight = .5, title.position = "top", title.hjust = .5
#     )
#   ) + theme_minimal()+
#   theme(
#     axis.title = element_blank(),
#     axis.ticks = element_blank(),
#     axis.text.y = element_blank(),
#     axis.text.x = element_text(color = "gray12", size = 12),
#     legend.position = "bottom"
#   )

#stage#
library(tidytree)
colnames(label)[1]<-"Sample.ID"
label$Sample.ID<-substr(label$Sample.ID,1,15)
brca_stage_450<-left_join(cluster_sample,label,by='Sample.ID')
a<-brca_stage_450$label
a[which(a=="1")]<-"stage I"
a[which(a=="2"|a=="3")]<-"stage II&III"
brca_stage_450$stage<-a
colnames(brca_stage_450)[2]<-"cluster"
brca_stage_450$cluster<-as.character(brca_stage_450$cluster)
aaa<-brca_stage_450 %>% 
  ggplot(aes(x = cluster, fill = stage)) +  
  geom_bar(position = position_fill()) + 
  scale_fill_brewer(palette = 'Set3') +  
  theme_classic() +  
  labs(y = 'Percent') +  
  coord_flip() 
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\COAD_cluster+stage100_plot.pdf",aaa,width = 10,height = 10)

 
library(data.table)
LUADsurvival=fread(paste0("G:\\dnameth\\COAD","\\","TCGA-","COAD",".GDC_phenotype.tsv.gz"),header=TRUE)
LUADsurvival<-as.data.frame(LUADsurvival)
rownames(LUADsurvival)<-LUADsurvival[,1]
LUADsurvival<-LUADsurvival[,-1]
rownames(LUADsurvival)=gsub("-","\\.",rownames(LUADsurvival))
unique(LUADsurvival$clinical_stage)
unique(LUADsurvival$tumor_stage)
matches_stage_i <-grepl("\\bstage i(a|b)?\\b", LUADsurvival$tumor_stage)
matches_stage_ii <-grepl("\\bstage ii(a|b|c)?\\b", LUADsurvival$tumor_stage)
matches_stage_iii <-grepl("\\bstage iii(a|b|c)?\\b", LUADsurvival$tumor_stage)
LUADsurvival$tumor_stage[matches_stage_i]<-"1"
LUADsurvival$tumor_stage[matches_stage_ii]<-"2"
LUADsurvival$tumor_stage[matches_stage_iii]<-"3"
LUADstage1=LUADsurvival[LUADsurvival$tumor_stage==1,]
LUADstage2=LUADsurvival[LUADsurvival$tumor_stage==2,]
LUADstage3=LUADsurvival[LUADsurvival$tumor_stage==3,]

labeltostage<-data.table(Sample.ID=c(rownames(LUADstage1),rownames(LUADstage2),rownames(LUADstage3)),stage=c(rep(c("Stage I","Stage II","Stage III"),c(nrow(LUADstage1),nrow(LUADstage2),nrow(LUADstage3)))))


# labeltostage<-read.csv("label2 stage.csv")
library(tidytree)
colnames(label)[1]<-"Sample.ID"
colnames(labeltostage)[1]<-"Sample.ID"
labeltostage$Sample.ID<-substr(labeltostage$Sample.ID,1,15)

# brca_stage_450<-left_join(cluster_sample,label,by='Sample.ID')
brca_stage_450<-left_join(brca_stage_450,labeltostage,by='Sample.ID')
# brca_stage_450<-brca_stage_450[,-c(4,5)]
colnames(brca_stage_450)[5]<-"stage"
colnames(brca_stage_450)[2]<-"cluster"
# colnames(brca_stage_450)[3]<-"stage"
# brca_stage_450$stage<-paste("stage",brca_stage_450$stage)
brca_stage_450$cluster<-paste("Cluster",brca_stage_450$cluster)
# brca_stage_450$stage[which(brca_stage_450$stage=="stage 2")]<-brca_stage_450$stage2[which(brca_stage_450$stage=="stage 2")]
# brca_stage_450$stage[which(brca_stage_450$stage=="stage 1")]<-"Stage I"
# brca_stage_450$stage[which(brca_stage_450$stage=="STAGE II")]<-"Stage II"
# brca_stage_450$stage[which(brca_stage_450$stage=="STAGE III")]<-"Stage III"
# devtools::install_github("davidsjoberg/ggsankey")
library(tidyverse)
library(ggsankey)
library(ggplot2)
library(cols4all)
head(brca_stage_450)
df <- brca_stage_450 %>%
  make_long(cluster, stage)
 
df$node <- factor(df$node,levels = c(rev(unique(brca_stage_450$stage)),
                                     rev(unique(brca_stage_450$cluster))))
 
c4a_gui()
mycol <- c4a('rainbow_wh_rd',53)
mycol <- sample(mycol,length(mycol))  
mycol<-c("#FFD4A9", "#80D1C8","pink","#7FB974","#5072BF","#AE84BB")
 
p <- ggplot(df, aes(x = x, next_x = next_x,
                    node = node, next_node = next_node,
                    fill = node,
                    label = node)) +  
  scale_fill_manual(values = mycol)  

p1 <- p +
  geom_sankey(flow.alpha = 0.5, 
              smooth = 8,  
              width = 0.12) +  
  geom_sankey_text(size = 3.2,  
                   color = "black") +
  theme_void() +
  theme(legend.position = 'none')  
p1
ggsave("G:\\pan_result\\gpu_result_tcga\\COAD\\450k\\COAD_cluster+stage_sankey_plot.pdf",p1,width = 10,height = 10)
#ESCA####

 
library(ggplot2)  
library(dplyr)  
#library(ggstatsplot) 

brca_tnm<-read.table("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\TCGA.ESCA.sampleMap_ESCA_clinicalMatrix",sep = "\t",header = T)
brca_tnm<-data.frame(Sample.ID=brca_tnm$sampleID,T=brca_tnm$pathologic_T,N=brca_tnm$pathologic_N,M=brca_tnm$pathologic_M)
sigmoid_450k<-data.frame(Sample.ID=rownames(exp),exp)
brca_tnm$Sample.ID<-gsub("-",".",brca_tnm$Sample.ID)
sigmoid_450k$Sample.ID<-substr(sigmoid_450k$Sample.ID,1,15)


brca_tnm_450k<-merge(brca_tnm,sigmoid_450k,by="Sample.ID")
brca_tnm_450k<-na.omit(brca_tnm_450k)
str(brca_tnm_450k)
table(brca_tnm_450k$T)
table(brca_tnm_450k$N)

brca_tnm_450k$T<-gsub(".*T4.*", "T4", brca_tnm_450k$T)
# brca_tnm_450k$N<-gsub(".*N1.*", "N1", brca_tnm_450k$N)
# brca_tnm_450k$N<-gsub(".*N2.*", "N2", brca_tnm_450k$N)
brca_tnm_450k<-brca_tnm_450k[-which(brca_tnm_450k$N=="NX"),]


cluster_sample<-data.frame(Sample.ID=rownames(clusters),class=clusters$clusters)
cluster_sample$Sample.ID<-substr(cluster_sample$Sample.ID,1,15)
brca_tnm_450k
brca_tn_450k<-merge(brca_tnm_450k,cluster_sample,by="Sample.ID")
#T 
brca_t_450k<-data.frame(Tfenqi=brca_tn_450k$T,cluster=brca_tn_450k$class,brca_tn_450k[,5:(ncol(brca_tnm_450k)-1)])
brca_t_450k<-na.omit(brca_t_450k)
brca_t_450k$cluster<-as.character(brca_t_450k$cluster)
aaa<-brca_t_450k %>% 
  ggplot(aes(x = cluster, fill = Tfenqi)) +  
  geom_bar(position = position_fill()) + 
  scale_fill_brewer(palette = 'Set3') +  
  theme_classic() +  
  labs(y = 'Percent') +  
  coord_flip()  
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\ESCA_cluster+T100_plot.pdf",aaa,width = 10,height = 10)



#N 
brca_n_450k<-data.frame(Nfenqi=brca_tn_450k$N,cluster=brca_tn_450k$class,brca_tn_450k[,5:(ncol(brca_tnm_450k)-1)])
brca_n_450k<-na.omit(brca_n_450k)
table(brca_n_450k$Nfenqi)
# brca_n_450k<-brca_n_450k[-which(brca_n_450k$Nfenqi=="NX"),]
brca_n_450k$cluster<-as.character(brca_n_450k$cluster)
aaa<-brca_n_450k %>% 
  ggplot(aes(x = cluster, fill = Nfenqi)) +  
  geom_bar(position = position_fill()) + 
  scale_fill_brewer(palette = 'Set3') +  
  theme_classic() + 
  labs(y = 'Percent') +  
  coord_flip()  
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\ESCA_cluster+N100_plot.pdf",aaa,width = 10,height = 10)

# 
library(ggplot2)
library(dplyr)
library(stringr)
t_cluster1<-as.data.frame(table(brca_t_450k[which(brca_t_450k$cluster==1),1]))
t_cluster1$Var1<-paste0("cluster1_",t_cluster1$Var1)
t_cluster1$Proportion<-t_cluster1$Freq/sum(t_cluster1$Freq)
t_cluster2<-as.data.frame(table(brca_t_450k[which(brca_t_450k$cluster==2),1]))
t_cluster2$Var1<-paste0("cluster2_",t_cluster2$Var1)
t_cluster2$Proportion<-t_cluster2$Freq/sum(t_cluster2$Freq)
t_cluster3<-as.data.frame(table(brca_t_450k[which(brca_t_450k$cluster==3),1]))
t_cluster3$Var1<-paste0("cluster3_",t_cluster3$Var1)
t_cluster3$Proportion<-t_cluster3$Freq/sum(t_cluster3$Freq)




t_cluster<-rbind(t_cluster1,t_cluster2,t_cluster3)
str(t_cluster)
p <- ggplot(t_cluster) + 
  geom_col(
    aes(
      x = reorder(Var1, Proportion),
      y = Proportion,
      fill = Proportion
    ),
    position = "dodge2",
    show.legend = TRUE,
    alpha = .9
  )
p
p <- p + coord_polar()
ggsave("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\ESCA_cluster+T_flower_plot.pdf",p,width = 10,height = 10)

# p+geom_hline(
#   aes(yintercept = y), 
#   data.frame(y = c(0:3) * 1000),
#   color = "lightgrey"
# )+
#   geom_segment(data = t_cluster,  #  
#     aes(
#       x = reorder(str_wrap(region, 5), sum_length),
#       y = 0,
#       xend = reorder(str_wrap(region, 5), sum_length),
#       yend = 3000
#     ),
#     linetype = "dashed",
#     color = "gray12"
#   ) +
#   scale_y_continuous(
#     limits = c(-1500, 3500),
#     expand = c(0, 0),
#     breaks = c(0, 1000, 2000, 3000)
#   ) + 
#   scale_fill_gradientn(
#     "Amount of Tracks",
#     colours = c( "#6C5B7B","#C06C84","#F67280","#F8B195")
#   ) +
#   guides(
#     fill = guide_colorsteps(
#       barwidth = 15, barheight = .5, title.position = "top", title.hjust = .5
#     )
#   ) + theme_minimal()+
#   theme(
#     axis.title = element_blank(),
#     axis.ticks = element_blank(),
#     axis.text.y = element_blank(),
#     axis.text.x = element_text(color = "gray12", size = 12),
#     legend.position = "bottom"
#   )

#stage#
library(tidytree)
colnames(label)[1]<-"Sample.ID"
label$Sample.ID<-substr(label$Sample.ID,1,15)
brca_stage_450<-left_join(cluster_sample,label,by='Sample.ID')
a<-brca_stage_450$label
a[which(a=="1")]<-"stage I"
a[which(a=="2"|a=="3")]<-"stage II&III"
brca_stage_450$stage<-a
colnames(brca_stage_450)[2]<-"cluster"
brca_stage_450$cluster<-as.character(brca_stage_450$cluster)
aaa<-brca_stage_450 %>% 
  ggplot(aes(x = cluster, fill = stage)) +  
  geom_bar(position = position_fill()) + 
  scale_fill_brewer(palette = 'Set3') + 
  theme_classic() +  
  labs(y = 'Percent') +  
  coord_flip()  
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\ESCA_cluster+stage100_plot.pdf",aaa,width = 10,height = 10)

 
library(data.table)
LUADsurvival=fread(paste0("G:\\dnameth\\ESCA","\\","TCGA-","ESCA",".GDC_phenotype.tsv.gz"),header=TRUE)
LUADsurvival<-as.data.frame(LUADsurvival)
rownames(LUADsurvival)<-LUADsurvival[,1]
LUADsurvival<-LUADsurvival[,-1]
rownames(LUADsurvival)=gsub("-","\\.",rownames(LUADsurvival))
unique(LUADsurvival$clinical_stage)
unique(LUADsurvival$tumor_stage)
matches_stage_i <-grepl("\\bstage i(a|b)?\\b", LUADsurvival$tumor_stage)
matches_stage_ii <-grepl("\\bstage ii(a|b|c)?\\b", LUADsurvival$tumor_stage)
matches_stage_iii <-grepl("\\bstage iii(a|b|c)?\\b", LUADsurvival$tumor_stage)
LUADsurvival$tumor_stage[matches_stage_i]<-"1"
LUADsurvival$tumor_stage[matches_stage_ii]<-"2"
LUADsurvival$tumor_stage[matches_stage_iii]<-"3"
LUADstage1=LUADsurvival[LUADsurvival$tumor_stage==1,]
LUADstage2=LUADsurvival[LUADsurvival$tumor_stage==2,]
LUADstage3=LUADsurvival[LUADsurvival$tumor_stage==3,]

labeltostage<-data.table(Sample.ID=c(rownames(LUADstage1),rownames(LUADstage2),rownames(LUADstage3)),stage=c(rep(c("Stage I","Stage II","Stage III"),c(nrow(LUADstage1),nrow(LUADstage2),nrow(LUADstage3)))))


# labeltostage<-read.csv("label2 stage.csv")
library(tidytree)
colnames(label)[1]<-"Sample.ID"
colnames(labeltostage)[1]<-"Sample.ID"
labeltostage$Sample.ID<-substr(labeltostage$Sample.ID,1,15)

# brca_stage_450<-left_join(cluster_sample,label,by='Sample.ID')
brca_stage_450<-left_join(brca_stage_450,labeltostage,by='Sample.ID')
# brca_stage_450<-brca_stage_450[,-c(4,5)]
colnames(brca_stage_450)[5]<-"stage"
colnames(brca_stage_450)[2]<-"cluster"
# colnames(brca_stage_450)[3]<-"stage"
# brca_stage_450$stage<-paste("stage",brca_stage_450$stage)
brca_stage_450$cluster<-paste("Cluster",brca_stage_450$cluster)
# brca_stage_450$stage[which(brca_stage_450$stage=="stage 2")]<-brca_stage_450$stage2[which(brca_stage_450$stage=="stage 2")]
# brca_stage_450$stage[which(brca_stage_450$stage=="stage 1")]<-"Stage I"
# brca_stage_450$stage[which(brca_stage_450$stage=="STAGE II")]<-"Stage II"
# brca_stage_450$stage[which(brca_stage_450$stage=="STAGE III")]<-"Stage III"
# devtools::install_github("davidsjoberg/ggsankey")
library(tidyverse)
library(ggsankey)
library(ggplot2)
library(cols4all)
head(brca_stage_450)
df <- brca_stage_450 %>%
  make_long(cluster, stage)
 
df$node <- factor(df$node,levels = c(rev(unique(brca_stage_450$stage)),
                                     rev(unique(brca_stage_450$cluster))))
 
c4a_gui()
mycol <- c4a('rainbow_wh_rd',53)
mycol <- sample(mycol,length(mycol))  
mycol<-c("#FFD4A9", "#80D1C8","pink","#7FB974","#5072BF","#AE84BB")
 
p <- ggplot(df, aes(x = x, next_x = next_x,
                    node = node, next_node = next_node,
                    fill = node,
                    label = node)) + 
  scale_fill_manual(values = mycol)  

p1 <- p +
  geom_sankey(flow.alpha = 0.5,  
              smooth = 8, 
              width = 0.12) +  
  geom_sankey_text(size = 3.2,  
                   color = "black") +
  theme_void() +
  theme(legend.position = 'none')  
p1
ggsave("G:\\pan_result\\gpu_result_tcga\\ESCA\\450k\\ESCA_cluster+stage_sankey_plot.pdf",p1,width = 10,height = 10)

#HNSC####

 
library(ggplot2)  
library(dplyr)  
#library(ggstatsplot)  

brca_tnm<-read.table("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\TCGA.HNSC.sampleMap_HNSC_clinicalMatrix",sep = "\t",header = T)
brca_tnm<-data.frame(Sample.ID=brca_tnm$sampleID,T=brca_tnm$pathologic_T,N=brca_tnm$pathologic_N,M=brca_tnm$pathologic_M)
sigmoid_450k<-data.frame(Sample.ID=rownames(exp),exp)
brca_tnm$Sample.ID<-gsub("-",".",brca_tnm$Sample.ID)
sigmoid_450k$Sample.ID<-substr(sigmoid_450k$Sample.ID,1,15)


brca_tnm_450k<-merge(brca_tnm,sigmoid_450k,by="Sample.ID")
brca_tnm_450k<-na.omit(brca_tnm_450k)
str(brca_tnm_450k)
table(brca_tnm_450k$T)
table(brca_tnm_450k$N)
table(brca_tnm_450k$M)

brca_tnm_450k$T<-gsub(".*T1.*", "T1", brca_tnm_450k$T)
brca_tnm_450k$T<-gsub(".*T2.*", "T2", brca_tnm_450k$T)
brca_tnm_450k$T<-gsub(".*T3.*", "T3", brca_tnm_450k$T)
brca_tnm_450k$T<-gsub(".*T4.*", "T4", brca_tnm_450k$T)
brca_tnm_450k$N<-gsub(".*N1.*", "N1", brca_tnm_450k$N)
brca_tnm_450k$N<-gsub(".*N2.*", "N2", brca_tnm_450k$N)
brca_tnm_450k<-brca_tnm_450k[-which(brca_tnm_450k$N==""),]
brca_tnm_450k<-brca_tnm_450k[-which(brca_tnm_450k$T=="TX"),]
brca_tnm_450k<-brca_tnm_450k[-which(brca_tnm_450k$N==""),]
brca_tnm_450k<-brca_tnm_450k[-which(brca_tnm_450k$N=="NX"),]


cluster_sample<-data.frame(Sample.ID=rownames(clusters),class=clusters$clusters)
cluster_sample$Sample.ID<-substr(cluster_sample$Sample.ID,1,15)
brca_tnm_450k
brca_tn_450k<-merge(brca_tnm_450k,cluster_sample,by="Sample.ID")
#T 
brca_t_450k<-data.frame(Tfenqi=brca_tn_450k$T,cluster=brca_tn_450k$class,brca_tn_450k[,5:(ncol(brca_tnm_450k)-1)])
brca_t_450k<-na.omit(brca_t_450k)
brca_t_450k$cluster<-as.character(brca_t_450k$cluster)
aaa<-brca_t_450k %>% 
  ggplot(aes(x = cluster, fill = Tfenqi)) + 
  geom_bar(position = position_fill()) + 
  scale_fill_brewer(palette = 'Set3') +  
  theme_classic() +  
  labs(y = 'Percent') +  
  coord_flip()  
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_cluster+T100_plot.pdf",aaa,width = 10,height = 10)



#N
brca_n_450k<-data.frame(Nfenqi=brca_tn_450k$N,cluster=brca_tn_450k$class,brca_tn_450k[,5:(ncol(brca_tnm_450k)-1)])
brca_n_450k<-na.omit(brca_n_450k)
table(brca_n_450k$Nfenqi)
# brca_n_450k<-brca_n_450k[-which(brca_n_450k$Nfenqi=="NX"),]
brca_n_450k$cluster<-as.character(brca_n_450k$cluster)
aaa<-brca_n_450k %>% 
  ggplot(aes(x = cluster, fill = Nfenqi)) +  
  geom_bar(position = position_fill()) + 
  scale_fill_brewer(palette = 'Set3') +  
  theme_classic() +  
  labs(y = 'Percent') +  
  coord_flip()  
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_cluster+N100_plot.pdf",aaa,width = 10,height = 10)




#stage#
library(tidytree)
colnames(label)[1]<-"Sample.ID"
label$Sample.ID<-substr(label$Sample.ID,1,15)
brca_stage_450<-left_join(cluster_sample,label,by='Sample.ID')
a<-brca_stage_450$label
a[which(a=="1")]<-"stage I"
a[which(a=="2"|a=="3")]<-"stage II&III"
brca_stage_450$stage<-a
colnames(brca_stage_450)[2]<-"cluster"
brca_stage_450$cluster<-as.character(brca_stage_450$cluster)
aaa<-brca_stage_450 %>% 
  ggplot(aes(x = cluster, fill = stage)) +  
  geom_bar(position = position_fill()) + 
  scale_fill_brewer(palette = 'Set3') +  
  theme_classic() +  
  labs(y = 'Percent') +  
  coord_flip()  
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_cluster+stage100_plot.pdf",aaa,width = 10,height = 10)

 
library(data.table)
LUADsurvival=fread(paste0("G:\\dnameth\\HNSC","\\","TCGA-","HNSC",".GDC_phenotype.tsv.gz"),header=TRUE)
LUADsurvival<-as.data.frame(LUADsurvival)
rownames(LUADsurvival)<-LUADsurvival[,1]
LUADsurvival<-LUADsurvival[,-1]
rownames(LUADsurvival)=gsub("-","\\.",rownames(LUADsurvival))
unique(LUADsurvival$clinical_stage)
unique(LUADsurvival$tumor_stage)
str(LUADsurvival$clinical_stage)

matches_stage_i <-grepl("\\bStage I(A|B)?\\b", LUADsurvival$clinical_stage)
matches_stage_ii <-grepl("\\bStage II(A|B|C)?\\b", LUADsurvival$clinical_stage)
matches_stage_iii <-grepl("\\bStage III(A|B|C)?\\b", LUADsurvival$clinical_stage)

matches_stage_i <-grepl("\\bstage i(a|b)?\\b", LUADsurvival$tumor_stage)
matches_stage_ii <-grepl("\\bstage ii(a|b|c)?\\b", LUADsurvival$tumor_stage)
matches_stage_iii <-grepl("\\bstage iii(a|b|c)?\\b", LUADsurvival$tumor_stage)
LUADsurvival$tumor_stage[matches_stage_i]<-"1"
LUADsurvival$tumor_stage[matches_stage_ii]<-"2"
LUADsurvival$tumor_stage[matches_stage_iii]<-"3"
LUADstage1=LUADsurvival[LUADsurvival$tumor_stage==1,]
LUADstage2=LUADsurvival[LUADsurvival$tumor_stage==2,]
LUADstage3=LUADsurvival[LUADsurvival$tumor_stage==3,]

labeltostage<-data.table(Sample.ID=c(rownames(LUADstage1),rownames(LUADstage2),rownames(LUADstage3)),stage=c(rep(c("Stage I","Stage II","Stage III"),c(nrow(LUADstage1),nrow(LUADstage2),nrow(LUADstage3)))))


# labeltostage<-read.csv("label2 stage.csv")
library(tidytree)
colnames(label)[1]<-"Sample.ID"
colnames(labeltostage)[1]<-"Sample.ID"
labeltostage$Sample.ID<-substr(labeltostage$Sample.ID,1,15)

# brca_stage_450<-left_join(cluster_sample,label,by='Sample.ID')
brca_stage_450<-left_join(brca_stage_450,labeltostage,by='Sample.ID')
# brca_stage_450<-brca_stage_450[,-c(4,5)]
colnames(brca_stage_450)[5]<-"stage"
colnames(brca_stage_450)[2]<-"cluster"
# colnames(brca_stage_450)[3]<-"stage"
# brca_stage_450$stage<-paste("stage",brca_stage_450$stage)
brca_stage_450$cluster<-paste("Cluster",brca_stage_450$cluster)

# HNSC
# brca_stage_450$stage[1:162]<-"Stage I"
# brca_stage_450$stage[163:231]<-"Stage II"
# brca_stage_450$stage[232:277]<-"Stage III"
# brca_stage_450$stage[which(brca_stage_450$stage=="stage 2")]<-brca_stage_450$stage2[which(brca_stage_450$stage=="stage 2")]
# brca_stage_450$stage[which(brca_stage_450$stage=="stage 1")]<-"Stage I"
# brca_stage_450$stage[which(brca_stage_450$stage=="STAGE II")]<-"Stage II"
# brca_stage_450$stage[which(brca_stage_450$stage=="STAGE III")]<-"Stage III"
# devtools::install_github("davidsjoberg/ggsankey")
library(tidyverse)
library(ggsankey)
library(ggplot2)
library(cols4all)
head(brca_stage_450)
df <- brca_stage_450 %>%
  make_long(cluster, stage)
 
df$node <- factor(df$node,levels = c(rev(unique(brca_stage_450$stage)),
                                     rev(unique(brca_stage_450$cluster))))
 
c4a_gui()
mycol <- c4a('rainbow_wh_rd',53)
mycol <- sample(mycol,length(mycol))  
mycol<-c("#FFD4A9", "#80D1C8","pink","#7FB974","#5072BF","#AE84BB")
 
p <- ggplot(df, aes(x = x, next_x = next_x,
                    node = node, next_node = next_node,
                    fill = node,
                    label = node)) +  
  scale_fill_manual(values = mycol) 

p1 <- p +
  geom_sankey(flow.alpha = 0.5,  
              smooth = 8,  
              width = 0.12) + 
  geom_sankey_text(size = 3.2,  
                   color = "black") +
  theme_void() +
  theme(legend.position = 'none')  
p1
ggsave("G:\\pan_result\\gpu_result_tcga\\HNSC\\450k\\HNSC_cluster+stage_sankey_plot.pdf",p1,width = 10,height = 10)

#KIRC####

 
library(ggplot2)  
library(dplyr)  
#library(ggstatsplot)  
brca_tnm<-read.table("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\TCGA.KIRC.sampleMap_KIRC_clinicalMatrix",sep = "\t",header = T)
brca_tnm<-data.frame(Sample.ID=brca_tnm$sampleID,T=brca_tnm$pathologic_T,N=brca_tnm$pathologic_N,M=brca_tnm$pathologic_M)
sigmoid_450k<-data.frame(Sample.ID=rownames(exp),exp)
brca_tnm$Sample.ID<-gsub("-",".",brca_tnm$Sample.ID)
sigmoid_450k$Sample.ID<-substr(sigmoid_450k$Sample.ID,1,15)


brca_tnm_450k<-merge(brca_tnm,sigmoid_450k,by="Sample.ID")
brca_tnm_450k<-na.omit(brca_tnm_450k)
str(brca_tnm_450k)
table(brca_tnm_450k$T)
table(brca_tnm_450k$N)
table(brca_tnm_450k$M)

brca_tnm_450k$T<-gsub(".*T1.*", "T1", brca_tnm_450k$T)
brca_tnm_450k$T<-gsub(".*T2.*", "T2", brca_tnm_450k$T)
brca_tnm_450k$T<-gsub(".*T3.*", "T3", brca_tnm_450k$T)
# brca_tnm_450k$N<-gsub(".*N1.*", "N1", brca_tnm_450k$N)
# brca_tnm_450k$N<-gsub(".*N2.*", "N2", brca_tnm_450k$N)
brca_tnm_450k<-brca_tnm_450k[-which(brca_tnm_450k$N=="NX"),]




cluster_sample<-data.frame(Sample.ID=rownames(clusters),class=clusters$clusters)
cluster_sample$Sample.ID<-substr(cluster_sample$Sample.ID,1,15)
brca_tnm_450k
brca_tn_450k<-merge(brca_tnm_450k,cluster_sample,by="Sample.ID")
#T 
brca_t_450k<-data.frame(Tfenqi=brca_tn_450k$T,cluster=brca_tn_450k$class,brca_tn_450k[,5:(ncol(brca_tnm_450k)-1)])
brca_t_450k<-na.omit(brca_t_450k)
brca_t_450k$cluster<-as.character(brca_t_450k$cluster)
aaa<-brca_t_450k %>% 
  ggplot(aes(x = cluster, fill = Tfenqi)) +  
  geom_bar(position = position_fill()) + 
  scale_fill_brewer(palette = 'Set3') +  
  theme_classic() +  
  labs(y = 'Percent') +  
  coord_flip()  
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\KIRC_cluster+T100_plot.pdf",aaa,width = 10,height = 10)



#N 
brca_n_450k<-data.frame(Nfenqi=brca_tn_450k$N,cluster=brca_tn_450k$class,brca_tn_450k[,5:(ncol(brca_tnm_450k)-1)])
brca_n_450k<-na.omit(brca_n_450k)
table(brca_n_450k$Nfenqi)
# brca_n_450k<-brca_n_450k[-which(brca_n_450k$Nfenqi=="NX"),]
brca_n_450k$cluster<-as.character(brca_n_450k$cluster)
aaa<-brca_n_450k %>% 
  ggplot(aes(x = cluster, fill = Nfenqi)) +  
  geom_bar(position = position_fill()) + 
  scale_fill_brewer(palette = 'Set3') +  
  theme_classic() +  
  labs(y = 'Percent') +  
  coord_flip()  
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\KIRC_cluster+N100_plot.pdf",aaa,width = 10,height = 10)




#stage#
library(tidytree)
colnames(label)[1]<-"Sample.ID"
label$Sample.ID<-substr(label$Sample.ID,1,15)
brca_stage_450<-left_join(cluster_sample,label,by='Sample.ID')
a<-brca_stage_450$label
a[which(a=="1")]<-"stage I"
a[which(a=="2"|a=="3")]<-"stage II&III"
brca_stage_450$stage<-a
colnames(brca_stage_450)[2]<-"cluster"
brca_stage_450$cluster<-as.character(brca_stage_450$cluster)
aaa<-brca_stage_450 %>% 
  ggplot(aes(x = cluster, fill = stage)) +  
  geom_bar(position = position_fill()) + 
  scale_fill_brewer(palette = 'Set3') +  
  theme_classic() +  
  labs(y = 'Percent') +  
  coord_flip()  
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\KIRC_cluster+stage100_plot.pdf",aaa,width = 10,height = 10)

 
library(data.table)
LUADsurvival=fread(paste0("G:\\dnameth\\KIRC","\\","TCGA-","KIRC",".GDC_phenotype.tsv.gz"),header=TRUE)
LUADsurvival<-as.data.frame(LUADsurvival)
rownames(LUADsurvival)<-LUADsurvival[,1]
LUADsurvival<-LUADsurvival[,-1]
rownames(LUADsurvival)=gsub("-","\\.",rownames(LUADsurvival))
unique(LUADsurvival$clinical_stage)
unique(LUADsurvival$tumor_stage)
str(LUADsurvival$clinical_stage)

# matches_stage_i <-grepl("\\bStage I(A|B)?\\b", LUADsurvival$clinical_stage)
# matches_stage_ii <-grepl("\\bStage II(A|B|C)?\\b", LUADsurvival$clinical_stage)
# matches_stage_iii <-grepl("\\bStage III(A|B|C)?\\b", LUADsurvival$clinical_stage)

matches_stage_i <-grepl("\\bstage i(a|b)?\\b", LUADsurvival$tumor_stage)
matches_stage_ii <-grepl("\\bstage ii(a|b|c)?\\b", LUADsurvival$tumor_stage)
matches_stage_iii <-grepl("\\bstage iii(a|b|c)?\\b", LUADsurvival$tumor_stage)
LUADsurvival$tumor_stage[matches_stage_i]<-"1"
LUADsurvival$tumor_stage[matches_stage_ii]<-"2"
LUADsurvival$tumor_stage[matches_stage_iii]<-"3"
LUADstage1=LUADsurvival[LUADsurvival$tumor_stage==1,]
LUADstage2=LUADsurvival[LUADsurvival$tumor_stage==2,]
LUADstage3=LUADsurvival[LUADsurvival$tumor_stage==3,]

labeltostage<-data.table(Sample.ID=c(rownames(LUADstage1),rownames(LUADstage2),rownames(LUADstage3)),stage=c(rep(c("Stage I","Stage II","Stage III"),c(nrow(LUADstage1),nrow(LUADstage2),nrow(LUADstage3)))))


# labeltostage<-read.csv("label2 stage.csv")
library(tidytree)
colnames(label)[1]<-"Sample.ID"
colnames(labeltostage)[1]<-"Sample.ID"
labeltostage$Sample.ID<-substr(labeltostage$Sample.ID,1,15)

# brca_stage_450<-left_join(cluster_sample,label,by='Sample.ID')
brca_stage_450<-left_join(brca_stage_450,labeltostage,by='Sample.ID')
# brca_stage_450<-brca_stage_450[,-c(4,5)]
colnames(brca_stage_450)[5]<-"stage"
colnames(brca_stage_450)[2]<-"cluster"
# colnames(brca_stage_450)[3]<-"stage"
# brca_stage_450$stage<-paste("stage",brca_stage_450$stage)
brca_stage_450$cluster<-paste("Cluster",brca_stage_450$cluster)

# 
# brca_stage_450$stage[1:162]<-"Stage I"
# brca_stage_450$stage[163:231]<-"Stage II"
# brca_stage_450$stage[232:277]<-"Stage III"
# brca_stage_450$stage[which(brca_stage_450$stage=="stage 2")]<-brca_stage_450$stage2[which(brca_stage_450$stage=="stage 2")]
# brca_stage_450$stage[which(brca_stage_450$stage=="stage 1")]<-"Stage I"
# brca_stage_450$stage[which(brca_stage_450$stage=="STAGE II")]<-"Stage II"
# brca_stage_450$stage[which(brca_stage_450$stage=="STAGE III")]<-"Stage III"
# devtools::install_github("davidsjoberg/ggsankey")
library(tidyverse)
library(ggsankey)
library(ggplot2)
library(cols4all)
head(brca_stage_450)
df <- brca_stage_450 %>%
  make_long(cluster, stage)
 
df$node <- factor(df$node,levels = c(rev(unique(brca_stage_450$stage)),
                                     rev(unique(brca_stage_450$cluster))))
 
c4a_gui()
mycol <- c4a('rainbow_wh_rd',53)
mycol <- sample(mycol,length(mycol))  
mycol<-c("#FFD4A9", "#80D1C8","pink","#7FB974","#5072BF","#AE84BB")

p <- ggplot(df, aes(x = x, next_x = next_x,
                    node = node, next_node = next_node,
                    fill = node,
                    label = node)) +  
  scale_fill_manual(values = mycol)  

p1 <- p +
  geom_sankey(flow.alpha = 0.5,  
              smooth = 8,  
              width = 0.12) +  
  geom_sankey_text(size = 3.2, 
                   color = "black") +
  theme_void() +
  theme(legend.position = 'none')  
p1
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRC\\450k\\KIRC_cluster+stage_sankey_plot.pdf",p1,width = 10,height = 10)

#KIRP####

 
library(ggplot2)  
library(dplyr)  
#library(ggstatsplot)  

brca_tnm<-read.table("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\TCGA.KIRP.sampleMap_KIRP_clinicalMatrix",sep = "\t",header = T)
brca_tnm<-data.frame(Sample.ID=brca_tnm$sampleID,T=brca_tnm$pathologic_T,N=brca_tnm$pathologic_N,M=brca_tnm$pathologic_M)
sigmoid_450k<-data.frame(Sample.ID=rownames(exp),exp)
brca_tnm$Sample.ID<-gsub("-",".",brca_tnm$Sample.ID)
sigmoid_450k$Sample.ID<-substr(sigmoid_450k$Sample.ID,1,15)


brca_tnm_450k<-merge(brca_tnm,sigmoid_450k,by="Sample.ID")
brca_tnm_450k<-na.omit(brca_tnm_450k)
str(brca_tnm_450k)
table(brca_tnm_450k$T)
table(brca_tnm_450k$N)
table(brca_tnm_450k$M)
brca_tnm_450k$T<-gsub(".*T1.*", "T1", brca_tnm_450k$T)
brca_tnm_450k$T<-gsub(".*T2.*", "T2", brca_tnm_450k$T)

brca_tnm_450k$T<-gsub(".*T3.*", "T3", brca_tnm_450k$T)
brca_tnm_450k$N<-gsub(".*N1.*", "N1", brca_tnm_450k$N)
brca_tnm_450k$N<-gsub(".*N2.*", "N2", brca_tnm_450k$N)
brca_tnm_450k<-brca_tnm_450k[-which(brca_tnm_450k$N=="NX"),]
brca_tnm_450k<-brca_tnm_450k[-which(brca_tnm_450k$N=="[Discrepancy]"),]

cluster_sample<-data.frame(Sample.ID=rownames(clusters),class=clusters$clusters)
cluster_sample$Sample.ID<-substr(cluster_sample$Sample.ID,1,15)
brca_tnm_450k
brca_tn_450k<-merge(brca_tnm_450k,cluster_sample,by="Sample.ID")
#T 
brca_t_450k<-data.frame(Tfenqi=brca_tn_450k$T,cluster=brca_tn_450k$class,brca_tn_450k[,5:(ncol(brca_tnm_450k)-1)])
brca_t_450k<-na.omit(brca_t_450k)
brca_t_450k$cluster<-as.character(brca_t_450k$cluster)
aaa<-brca_t_450k %>% 
  ggplot(aes(x = cluster, fill = Tfenqi)) +  
  geom_bar(position = position_fill()) + 
  scale_fill_brewer(palette = 'Set3') +  
  theme_classic() +  
  labs(y = 'Percent') +  
  coord_flip()  
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\KIRP_cluster+T100_plot.pdf",aaa,width = 10,height = 10)



#N 
brca_n_450k<-data.frame(Nfenqi=brca_tn_450k$N,cluster=brca_tn_450k$class,brca_tn_450k[,5:(ncol(brca_tnm_450k)-1)])
brca_n_450k<-na.omit(brca_n_450k)
table(brca_n_450k$Nfenqi)
# brca_n_450k<-brca_n_450k[-which(brca_n_450k$Nfenqi=="NX"),]
brca_n_450k$cluster<-as.character(brca_n_450k$cluster)
aaa<-brca_n_450k %>% 
  ggplot(aes(x = cluster, fill = Nfenqi)) +  
  geom_bar(position = position_fill()) + 
  scale_fill_brewer(palette = 'Set3') +  
  theme_classic() +  
  labs(y = 'Percent') +  
  coord_flip()  
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\KIRP_cluster+N100_plot.pdf",aaa,width = 10,height = 10)

# p+geom_hline(
#   aes(yintercept = y), 
#   data.frame(y = c(0:3) * 1000),
#   color = "lightgrey"
# )+
#   geom_segment(data = t_cluster,  # 
#     aes(
#       x = reorder(str_wrap(region, 5), sum_length),
#       y = 0,
#       xend = reorder(str_wrap(region, 5), sum_length),
#       yend = 3000
#     ),
#     linetype = "dashed",
#     color = "gray12"
#   ) +
#   scale_y_continuous(
#     limits = c(-1500, 3500),
#     expand = c(0, 0),
#     breaks = c(0, 1000, 2000, 3000)
#   ) + 
#   scale_fill_gradientn(
#     "Amount of Tracks",
#     colours = c( "#6C5B7B","#C06C84","#F67280","#F8B195")
#   ) +
#   guides(
#     fill = guide_colorsteps(
#       barwidth = 15, barheight = .5, title.position = "top", title.hjust = .5
#     )
#   ) + theme_minimal()+
#   theme(
#     axis.title = element_blank(),
#     axis.ticks = element_blank(),
#     axis.text.y = element_blank(),
#     axis.text.x = element_text(color = "gray12", size = 12),
#     legend.position = "bottom"
#   )

#stage#
library(tidytree)
colnames(label)[1]<-"Sample.ID"
label$Sample.ID<-substr(label$Sample.ID,1,15)
brca_stage_450<-left_join(cluster_sample,label,by='Sample.ID')
a<-brca_stage_450$label
a[which(a=="1")]<-"stage I"
a[which(a=="2"|a=="3")]<-"stage II&III"
brca_stage_450$stage<-a
colnames(brca_stage_450)[2]<-"cluster"
brca_stage_450$cluster<-as.character(brca_stage_450$cluster)
aaa<-brca_stage_450 %>% 
  ggplot(aes(x = cluster, fill = stage)) +  
  geom_bar(position = position_fill()) + 
  scale_fill_brewer(palette = 'Set3') +  
  theme_classic() +  
  labs(y = 'Percent') +  
  coord_flip()  
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\KIRP_cluster+stage100_plot.pdf",aaa,width = 10,height = 10)

 
library(data.table)
LUADsurvival=fread(paste0("G:\\dnameth\\KIRP","\\","TCGA-","KIRP",".GDC_phenotype.tsv.gz"),header=TRUE)
LUADsurvival<-as.data.frame(LUADsurvival)
rownames(LUADsurvival)<-LUADsurvival[,1]
LUADsurvival<-LUADsurvival[,-1]
rownames(LUADsurvival)=gsub("-","\\.",rownames(LUADsurvival))
unique(LUADsurvival$clinical_stage)
unique(LUADsurvival$tumor_stage)
matches_stage_i <-grepl("\\bstage i(a|b)?\\b", LUADsurvival$tumor_stage)
matches_stage_ii <-grepl("\\bstage ii(a|b|c)?\\b", LUADsurvival$tumor_stage)
matches_stage_iii <-grepl("\\bstage iii(a|b|c)?\\b", LUADsurvival$tumor_stage)
LUADsurvival$tumor_stage[matches_stage_i]<-"1"
LUADsurvival$tumor_stage[matches_stage_ii]<-"2"
LUADsurvival$tumor_stage[matches_stage_iii]<-"3"
LUADstage1=LUADsurvival[LUADsurvival$tumor_stage==1,]
LUADstage2=LUADsurvival[LUADsurvival$tumor_stage==2,]
LUADstage3=LUADsurvival[LUADsurvival$tumor_stage==3,]

labeltostage<-data.table(Sample.ID=c(rownames(LUADstage1),rownames(LUADstage2),rownames(LUADstage3)),stage=c(rep(c("Stage I","Stage II","Stage III"),c(nrow(LUADstage1),nrow(LUADstage2),nrow(LUADstage3)))))


# labeltostage<-read.csv("label2 stage.csv")
library(tidytree)
colnames(label)[1]<-"Sample.ID"
colnames(labeltostage)[1]<-"Sample.ID"
labeltostage$Sample.ID<-substr(labeltostage$Sample.ID,1,15)

# brca_stage_450<-left_join(cluster_sample,label,by='Sample.ID')
brca_stage_450<-left_join(brca_stage_450,labeltostage,by='Sample.ID')
# brca_stage_450<-brca_stage_450[,-c(4,5)]
colnames(brca_stage_450)[5]<-"stage"
colnames(brca_stage_450)[2]<-"cluster"
# colnames(brca_stage_450)[3]<-"stage"
# brca_stage_450$stage<-paste("stage",brca_stage_450$stage)
brca_stage_450$cluster<-paste("Cluster",brca_stage_450$cluster)
# brca_stage_450$stage[which(brca_stage_450$stage=="stage 2")]<-brca_stage_450$stage2[which(brca_stage_450$stage=="stage 2")]
# brca_stage_450$stage[which(brca_stage_450$stage=="stage 1")]<-"Stage I"
# brca_stage_450$stage[which(brca_stage_450$stage=="STAGE II")]<-"Stage II"
# brca_stage_450$stage[which(brca_stage_450$stage=="STAGE III")]<-"Stage III"
# devtools::install_github("davidsjoberg/ggsankey")
library(tidyverse)
library(ggsankey)
library(ggplot2)
library(cols4all)
head(brca_stage_450)
df <- brca_stage_450 %>%
  make_long(cluster, stage)
 
df$node <- factor(df$node,levels = c(rev(unique(brca_stage_450$stage)),
                                     rev(unique(brca_stage_450$cluster))))
 
c4a_gui()
mycol <- c4a('rainbow_wh_rd',53)
mycol <- sample(mycol,length(mycol)) 
mycol<-c("#FFD4A9", "#80D1C8","pink","#7FB974","#5072BF","#AE84BB")
 
p <- ggplot(df, aes(x = x, next_x = next_x,
                    node = node, next_node = next_node,
                    fill = node,
                    label = node)) +  
  scale_fill_manual(values = mycol)  

p1 <- p +
  geom_sankey(flow.alpha = 0.5,  
              smooth = 8,  
              width = 0.12) +  
  geom_sankey_text(size = 3.2,  
                   color = "black") +
  theme_void() +
  theme(legend.position = 'none')  
p1
ggsave("G:\\pan_result\\gpu_result_tcga\\KIRP\\450k\\KIRP_cluster+stage_sankey_plot.pdf",p1,width = 10,height = 10)
#LIHC####


library(ggplot2)  
library(dplyr)  
#library(ggstatsplot)  

brca_tnm<-read.table("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\TCGA.LIHC.sampleMap_LIHC_clinicalMatrix",sep = "\t",header = T)
brca_tnm<-data.frame(Sample.ID=brca_tnm$sampleID,T=brca_tnm$pathologic_T,N=brca_tnm$pathologic_N,M=brca_tnm$pathologic_M)
sigmoid_450k<-data.frame(Sample.ID=rownames(exp),exp)
brca_tnm$Sample.ID<-gsub("-",".",brca_tnm$Sample.ID)
sigmoid_450k$Sample.ID<-substr(sigmoid_450k$Sample.ID,1,15)


brca_tnm_450k<-merge(brca_tnm,sigmoid_450k,by="Sample.ID")
brca_tnm_450k<-na.omit(brca_tnm_450k)
str(brca_tnm_450k)
table(brca_tnm_450k$T)
table(brca_tnm_450k$N)
table(brca_tnm_450k$M)

# brca_tnm_450k$T<-gsub(".*T1.*", "T1", brca_tnm_450k$T)
brca_tnm_450k$T<-gsub(".*T2.*", "T2", brca_tnm_450k$T)
brca_tnm_450k$T<-gsub(".*T3.*", "T3", brca_tnm_450k$T)
# brca_tnm_450k$T<-gsub(".*T4.*", "T4", brca_tnm_450k$T)
# brca_tnm_450k$N<-gsub(".*N1.*", "N1", brca_tnm_450k$N)
# brca_tnm_450k$N<-gsub(".*N2.*", "N2", brca_tnm_450k$N)
brca_tnm_450k<-brca_tnm_450k[-which(brca_tnm_450k$N==""),]
brca_tnm_450k<-brca_tnm_450k[-which(brca_tnm_450k$T=="[Discrepancy]"),]
# brca_tnm_450k<-brca_tnm_450k[-which(brca_tnm_450k$N==""),]
brca_tnm_450k<-brca_tnm_450k[-which(brca_tnm_450k$N=="NX"),]


cluster_sample<-data.frame(Sample.ID=rownames(clusters),class=clusters$clusters)
cluster_sample$Sample.ID<-substr(cluster_sample$Sample.ID,1,15)
brca_tnm_450k
brca_tn_450k<-merge(brca_tnm_450k,cluster_sample,by="Sample.ID")
#T 
brca_t_450k<-data.frame(Tfenqi=brca_tn_450k$T,cluster=brca_tn_450k$class,brca_tn_450k[,5:(ncol(brca_tnm_450k)-1)])
brca_t_450k<-na.omit(brca_t_450k)
brca_t_450k$cluster<-as.character(brca_t_450k$cluster)
aaa<-brca_t_450k %>% 
  ggplot(aes(x = cluster, fill = Tfenqi)) +  
  geom_bar(position = position_fill()) + 
  scale_fill_brewer(palette = 'Set3') +  
  theme_classic() +  
  labs(y = 'Percent') +  
  coord_flip()  
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\LIHC_cluster+T100_plot.pdf",aaa,width = 10,height = 10)



#N 
brca_n_450k<-data.frame(Nfenqi=brca_tn_450k$N,cluster=brca_tn_450k$class,brca_tn_450k[,5:(ncol(brca_tnm_450k)-1)])
brca_n_450k<-na.omit(brca_n_450k)
table(brca_n_450k$Nfenqi)
# brca_n_450k<-brca_n_450k[-which(brca_n_450k$Nfenqi=="NX"),]
brca_n_450k$cluster<-as.character(brca_n_450k$cluster)
aaa<-brca_n_450k %>% 
  ggplot(aes(x = cluster, fill = Nfenqi)) +  
  geom_bar(position = position_fill()) + 
  scale_fill_brewer(palette = 'Set3') +  
  theme_classic() +  
  labs(y = 'Percent') +  
  coord_flip()  
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\LIHC_cluster+N100_plot.pdf",aaa,width = 10,height = 10)




#stage#
library(tidytree)
colnames(label)[1]<-"Sample.ID"
label$Sample.ID<-substr(label$Sample.ID,1,15)
brca_stage_450<-left_join(cluster_sample,label,by='Sample.ID')
a<-brca_stage_450$label
a[which(a=="1")]<-"stage I"
a[which(a=="2"|a=="3")]<-"stage II&III"
brca_stage_450$stage<-a
colnames(brca_stage_450)[2]<-"cluster"
brca_stage_450$cluster<-as.character(brca_stage_450$cluster)
aaa<-brca_stage_450 %>% 
  ggplot(aes(x = cluster, fill = stage)) +  
  geom_bar(position = position_fill()) + 
  scale_fill_brewer(palette = 'Set3') +  
  theme_classic() +  
  labs(y = 'Percent') +  
  coord_flip()  
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\LIHC_cluster+stage100_plot.pdf",aaa,width = 10,height = 10)

 
library(data.table)
LUADsurvival=fread(paste0("G:\\dnameth\\LIHC","\\","TCGA-","LIHC",".GDC_phenotype.tsv.gz"),header=TRUE)
LUADsurvival<-as.data.frame(LUADsurvival)
rownames(LUADsurvival)<-LUADsurvival[,1]
LUADsurvival<-LUADsurvival[,-1]
rownames(LUADsurvival)=gsub("-","\\.",rownames(LUADsurvival))
unique(LUADsurvival$clinical_stage)
unique(LUADsurvival$tumor_stage)
str(LUADsurvival$clinical_stage)

matches_stage_i <-grepl("\\bStage I(A|B)?\\b", LUADsurvival$clinical_stage)
matches_stage_ii <-grepl("\\bStage II(A|B|C)?\\b", LUADsurvival$clinical_stage)
matches_stage_iii <-grepl("\\bStage III(A|B|C)?\\b", LUADsurvival$clinical_stage)

matches_stage_i <-grepl("\\bstage i(a|b)?\\b", LUADsurvival$tumor_stage)
matches_stage_ii <-grepl("\\bstage ii(a|b|c)?\\b", LUADsurvival$tumor_stage)
matches_stage_iii <-grepl("\\bstage iii(a|b|c)?\\b", LUADsurvival$tumor_stage)
LUADsurvival$tumor_stage[matches_stage_i]<-"1"
LUADsurvival$tumor_stage[matches_stage_ii]<-"2"
LUADsurvival$tumor_stage[matches_stage_iii]<-"3"
LUADstage1=LUADsurvival[LUADsurvival$tumor_stage==1,]
LUADstage2=LUADsurvival[LUADsurvival$tumor_stage==2,]
LUADstage3=LUADsurvival[LUADsurvival$tumor_stage==3,]

labeltostage<-data.table(Sample.ID=c(rownames(LUADstage1),rownames(LUADstage2),rownames(LUADstage3)),stage=c(rep(c("Stage I","Stage II","Stage III"),c(nrow(LUADstage1),nrow(LUADstage2),nrow(LUADstage3)))))


# labeltostage<-read.csv("label2 stage.csv")
library(tidytree)
colnames(label)[1]<-"Sample.ID"
colnames(labeltostage)[1]<-"Sample.ID"
labeltostage$Sample.ID<-substr(labeltostage$Sample.ID,1,15)

# brca_stage_450<-left_join(cluster_sample,label,by='Sample.ID')
brca_stage_450<-left_join(brca_stage_450,labeltostage,by='Sample.ID')
# brca_stage_450<-brca_stage_450[,-c(4,5)]
colnames(brca_stage_450)[5]<-"stage"
colnames(brca_stage_450)[2]<-"cluster"
# colnames(brca_stage_450)[3]<-"stage"
# brca_stage_450$stage<-paste("stage",brca_stage_450$stage)
brca_stage_450$cluster<-paste("Cluster",brca_stage_450$cluster)

# # LIHC
# brca_stage_450$stage[1:162]<-"Stage I"
# brca_stage_450$stage[163:231]<-"Stage II"
# brca_stage_450$stage[232:277]<-"Stage III"
# brca_stage_450$stage[which(brca_stage_450$stage=="stage 2")]<-brca_stage_450$stage2[which(brca_stage_450$stage=="stage 2")]
# brca_stage_450$stage[which(brca_stage_450$stage=="stage 1")]<-"Stage I"
# brca_stage_450$stage[which(brca_stage_450$stage=="STAGE II")]<-"Stage II"
# brca_stage_450$stage[which(brca_stage_450$stage=="STAGE III")]<-"Stage III"
# devtools::install_github("davidsjoberg/ggsankey")
library(tidyverse)
library(ggsankey)
library(ggplot2)
library(cols4all)
head(brca_stage_450)
df <- brca_stage_450 %>%
  make_long(cluster, stage)
 
df$node <- factor(df$node,levels = c(rev(unique(brca_stage_450$stage)),
                                     rev(unique(brca_stage_450$cluster))))

c4a_gui()
mycol <- c4a('rainbow_wh_rd',53)
mycol <- sample(mycol,length(mycol))  
mycol<-c("#FFD4A9", "#80D1C8","pink","#7FB974","#5072BF","#AE84BB")
 
p <- ggplot(df, aes(x = x, next_x = next_x,
                    node = node, next_node = next_node,
                    fill = node,
                    label = node)) +  
  scale_fill_manual(values = mycol)  

p1 <- p +
  geom_sankey(flow.alpha = 0.5,  
              smooth = 8,  
              width = 0.12) +  
  geom_sankey_text(size = 3.2,  
                   color = "black") +
  theme_void() +
  theme(legend.position = 'none')  
p1
ggsave("G:\\pan_result\\gpu_result_tcga\\LIHC\\450k\\LIHC_cluster+stage_sankey_plot.pdf",p1,width = 10,height = 10)


#LUNG####

 
library(ggplot2) 
library(dplyr)  
#library(ggstatsplot)  

brca_tnm<-read.table("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\TCGA.NSCLC.sampleMap_NSCLC_clinicalMatrix",sep = "\t",header = T)
brca_tnm<-data.frame(Sample.ID=brca_tnm$sampleID,T=brca_tnm$pathologic_T,N=brca_tnm$pathologic_N,M=brca_tnm$pathologic_M)
sigmoid_450k<-data.frame(Sample.ID=rownames(exp),exp)
brca_tnm$Sample.ID<-gsub("-",".",brca_tnm$Sample.ID)
sigmoid_450k$Sample.ID<-substr(sigmoid_450k$Sample.ID,1,15)


brca_tnm_450k<-merge(brca_tnm,sigmoid_450k,by="Sample.ID")
brca_tnm_450k<-na.omit(brca_tnm_450k)
str(brca_tnm_450k)
table(brca_tnm_450k$T)
table(brca_tnm_450k$N)
table(brca_tnm_450k$M)

brca_tnm_450k$T<-gsub(".*T1.*", "T1", brca_tnm_450k$T)
brca_tnm_450k$T<-gsub(".*T2.*", "T2", brca_tnm_450k$T)
brca_tnm_450k$T<-gsub(".*T3.*", "T3", brca_tnm_450k$T)
brca_tnm_450k$T<-gsub(".*T4.*", "T4", brca_tnm_450k$T)
brca_tnm_450k$N<-gsub(".*N1.*", "N1", brca_tnm_450k$N)
brca_tnm_450k$N<-gsub(".*N2.*", "N2", brca_tnm_450k$N)
brca_tnm_450k<-brca_tnm_450k[-which(brca_tnm_450k$N==""),]
brca_tnm_450k<-brca_tnm_450k[-which(brca_tnm_450k$T=="TX"),]
brca_tnm_450k<-brca_tnm_450k[-which(brca_tnm_450k$N==""),]
brca_tnm_450k<-brca_tnm_450k[-which(brca_tnm_450k$N=="NX"),]


cluster_sample<-data.frame(Sample.ID=rownames(clusters),class=clusters$clusters)
cluster_sample$Sample.ID<-substr(cluster_sample$Sample.ID,1,15)
brca_tnm_450k
brca_tn_450k<-merge(brca_tnm_450k,cluster_sample,by="Sample.ID")
#T 
brca_t_450k<-data.frame(Tfenqi=brca_tn_450k$T,cluster=brca_tn_450k$class,brca_tn_450k[,5:(ncol(brca_tnm_450k)-1)])
brca_t_450k<-na.omit(brca_t_450k)
brca_t_450k$cluster<-as.character(brca_t_450k$cluster)
aaa<-brca_t_450k %>% 
  ggplot(aes(x = cluster, fill = Tfenqi)) +  
  geom_bar(position = position_fill()) + 
  scale_fill_brewer(palette = 'Set3') +  
  theme_classic() +  
  labs(y = 'Percent') +  
  coord_flip()  
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\NSCLC_cluster+T100_plot.pdf",aaa,width = 10,height = 10)



#N 
brca_n_450k<-data.frame(Nfenqi=brca_tn_450k$N,cluster=brca_tn_450k$class,brca_tn_450k[,5:(ncol(brca_tnm_450k)-1)])
brca_n_450k<-na.omit(brca_n_450k)
table(brca_n_450k$Nfenqi)
# brca_n_450k<-brca_n_450k[-which(brca_n_450k$Nfenqi=="NX"),]
brca_n_450k$cluster<-as.character(brca_n_450k$cluster)
aaa<-brca_n_450k %>% 
  ggplot(aes(x = cluster, fill = Nfenqi)) +  
  geom_bar(position = position_fill()) + 
  scale_fill_brewer(palette = 'Set3') +  
  theme_classic() +  
  labs(y = 'Percent') +  
  coord_flip()  
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\NSCLC_cluster+N100_plot.pdf",aaa,width = 10,height = 10)




#stage#
library(tidytree)
colnames(label)[1]<-"Sample.ID"
label$Sample.ID<-substr(label$Sample.ID,1,15)
brca_stage_450<-left_join(cluster_sample,label,by='Sample.ID')
a<-brca_stage_450$label
a[which(a=="1")]<-"stage I"
a[which(a=="2"|a=="3")]<-"stage II&III"
brca_stage_450$stage<-a
colnames(brca_stage_450)[2]<-"cluster"
brca_stage_450$cluster<-as.character(brca_stage_450$cluster)
aaa<-brca_stage_450 %>% 
  ggplot(aes(x = cluster, fill = stage)) +  
  geom_bar(position = position_fill()) + 
  scale_fill_brewer(palette = 'Set3') +  
  theme_classic() +  
  labs(y = 'Percent') +  
  coord_flip()  
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\NSCLC_cluster+stage100_plot.pdf",aaa,width = 10,height = 10)

 
library(data.table)
LUADsurvival1=fread(paste0("G:\\dnameth\\LUSC","\\","TCGA-","LUSC",".GDC_phenotype.tsv.gz"),header=TRUE)
LUADsurvival2=fread(paste0("G:\\dnameth\\LUAD","\\","TCGA-","LUAD",".GDC_phenotype.tsv.gz"),header=TRUE)
LUADsurvival<-rbind(LUADsurvival1,LUADsurvival2)
LUADsurvival<-as.data.frame(LUADsurvival)
rownames(LUADsurvival)<-LUADsurvival[,1]
LUADsurvival<-LUADsurvival[,-1]
rownames(LUADsurvival)=gsub("-","\\.",rownames(LUADsurvival))
unique(LUADsurvival$clinical_stage)
unique(LUADsurvival$tumor_stage)
str(LUADsurvival$clinical_stage)

matches_stage_i <-grepl("\\bStage I(A|B)?\\b", LUADsurvival$clinical_stage)
matches_stage_ii <-grepl("\\bStage II(A|B|C)?\\b", LUADsurvival$clinical_stage)
matches_stage_iii <-grepl("\\bStage III(A|B|C)?\\b", LUADsurvival$clinical_stage)

matches_stage_i <-grepl("\\bstage i(a|b)?\\b", LUADsurvival$tumor_stage)
matches_stage_ii <-grepl("\\bstage ii(a|b|c)?\\b", LUADsurvival$tumor_stage)
matches_stage_iii <-grepl("\\bstage iii(a|b|c)?\\b", LUADsurvival$tumor_stage)
LUADsurvival$tumor_stage[matches_stage_i]<-"1"
LUADsurvival$tumor_stage[matches_stage_ii]<-"2"
LUADsurvival$tumor_stage[matches_stage_iii]<-"3"
LUADstage1=LUADsurvival[LUADsurvival$tumor_stage==1,]
LUADstage2=LUADsurvival[LUADsurvival$tumor_stage==2,]
LUADstage3=LUADsurvival[LUADsurvival$tumor_stage==3,]

labeltostage<-data.table(Sample.ID=c(rownames(LUADstage1),rownames(LUADstage2),rownames(LUADstage3)),stage=c(rep(c("Stage I","Stage II","Stage III"),c(nrow(LUADstage1),nrow(LUADstage2),nrow(LUADstage3)))))


# labeltostage<-read.csv("label2 stage.csv")
library(tidytree)
colnames(label)[1]<-"Sample.ID"
colnames(labeltostage)[1]<-"Sample.ID"
labeltostage$Sample.ID<-substr(labeltostage$Sample.ID,1,15)

# brca_stage_450<-left_join(cluster_sample,label,by='Sample.ID')
brca_stage_450<-left_join(brca_stage_450,labeltostage,by='Sample.ID')
# brca_stage_450<-brca_stage_450[,-c(4,5)]
colnames(brca_stage_450)[5]<-"stage"
colnames(brca_stage_450)[2]<-"cluster"
# colnames(brca_stage_450)[3]<-"stage"
# brca_stage_450$stage<-paste("stage",brca_stage_450$stage)
brca_stage_450$cluster<-paste("Cluster",brca_stage_450$cluster)

# NSCLC
# brca_stage_450$stage[1:162]<-"Stage I"
# brca_stage_450$stage[163:231]<-"Stage II"
# brca_stage_450$stage[232:277]<-"Stage III"
# brca_stage_450$stage[which(brca_stage_450$stage=="stage 2")]<-brca_stage_450$stage2[which(brca_stage_450$stage=="stage 2")]
# brca_stage_450$stage[which(brca_stage_450$stage=="stage 1")]<-"Stage I"
# brca_stage_450$stage[which(brca_stage_450$stage=="STAGE II")]<-"Stage II"
# brca_stage_450$stage[which(brca_stage_450$stage=="STAGE III")]<-"Stage III"
# devtools::install_github("davidsjoberg/ggsankey")
library(tidyverse)
library(ggsankey)
library(ggplot2)
library(cols4all)
head(brca_stage_450)
df <- brca_stage_450 %>%
  make_long(cluster, stage)

df$node <- factor(df$node,levels = c(rev(unique(brca_stage_450$stage)),
                                     rev(unique(brca_stage_450$cluster))))
 
c4a_gui()
mycol <- c4a('rainbow_wh_rd',53)
mycol <- sample(mycol,length(mycol))  
mycol<-c("#FFD4A9", "#80D1C8","pink","#7FB974","#5072BF","#AE84BB")
 
p <- ggplot(df, aes(x = x, next_x = next_x,
                    node = node, next_node = next_node,
                    fill = node,
                    label = node)) +  
  scale_fill_manual(values = mycol)  

p1 <- p +
  geom_sankey(flow.alpha = 0.5,  
              smooth = 8,  
              width = 0.12) +  
  geom_sankey_text(size = 3.2,  
                   color = "black") +
  theme_void() +
  theme(legend.position = 'none')  
p1
ggsave("G:\\pan_result\\gpu_result_tcga\\NSCLC\\450k\\NSCLC_cluster+stage_sankey_plot.pdf",p1,width = 10,height = 10)

#PAAD####

 
library(ggplot2)  
library(dplyr)  
#library(ggstatsplot)  

brca_tnm<-read.table("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\TCGA.PAAD.sampleMap_PAAD_clinicalMatrix",sep = "\t",header = T)
brca_tnm<-data.frame(Sample.ID=brca_tnm$sampleID,T=brca_tnm$pathologic_T,N=brca_tnm$pathologic_N,M=brca_tnm$pathologic_M)
sigmoid_450k<-data.frame(Sample.ID=rownames(exp),exp)
brca_tnm$Sample.ID<-gsub("-",".",brca_tnm$Sample.ID)
sigmoid_450k$Sample.ID<-substr(sigmoid_450k$Sample.ID,1,15)


brca_tnm_450k<-merge(brca_tnm,sigmoid_450k,by="Sample.ID")
brca_tnm_450k<-na.omit(brca_tnm_450k)
str(brca_tnm_450k)
table(brca_tnm_450k$T)
table(brca_tnm_450k$N)
table(brca_tnm_450k$M)

brca_tnm_450k$T<-gsub(".*T1.*", "T1", brca_tnm_450k$T)
brca_tnm_450k$T<-gsub(".*T2.*", "T2", brca_tnm_450k$T)
brca_tnm_450k$T<-gsub(".*T3.*", "T3", brca_tnm_450k$T)
brca_tnm_450k$T<-gsub(".*T4.*", "T4", brca_tnm_450k$T)
brca_tnm_450k$N<-gsub(".*N1.*", "N1", brca_tnm_450k$N)
brca_tnm_450k$N<-gsub(".*N2.*", "N2", brca_tnm_450k$N)
brca_tnm_450k<-brca_tnm_450k[-which(brca_tnm_450k$N==""),]
brca_tnm_450k<-brca_tnm_450k[-which(brca_tnm_450k$T=="TX"),]
brca_tnm_450k<-brca_tnm_450k[-which(brca_tnm_450k$N==""),]
brca_tnm_450k<-brca_tnm_450k[-which(brca_tnm_450k$N=="NX"),]


cluster_sample<-data.frame(Sample.ID=rownames(clusters),class=clusters$clusters)
cluster_sample$Sample.ID<-substr(cluster_sample$Sample.ID,1,15)
brca_tnm_450k
brca_tn_450k<-merge(brca_tnm_450k,cluster_sample,by="Sample.ID")
#T 
brca_t_450k<-data.frame(Tfenqi=brca_tn_450k$T,cluster=brca_tn_450k$class,brca_tn_450k[,5:(ncol(brca_tnm_450k)-1)])
brca_t_450k<-na.omit(brca_t_450k)
brca_t_450k$cluster<-as.character(brca_t_450k$cluster)
aaa<-brca_t_450k %>% 
  ggplot(aes(x = cluster, fill = Tfenqi)) +  
  geom_bar(position = position_fill()) + 
  scale_fill_brewer(palette = 'Set3') +  
  theme_classic() +  
  labs(y = 'Percent') +  
  coord_flip()  
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\PAAD_cluster+T100_plot.pdf",aaa,width = 10,height = 10)



#N 
brca_n_450k<-data.frame(Nfenqi=brca_tn_450k$N,cluster=brca_tn_450k$class,brca_tn_450k[,5:(ncol(brca_tnm_450k)-1)])
brca_n_450k<-na.omit(brca_n_450k)
table(brca_n_450k$Nfenqi)
# brca_n_450k<-brca_n_450k[-which(brca_n_450k$Nfenqi=="NX"),]
brca_n_450k$cluster<-as.character(brca_n_450k$cluster)
aaa<-brca_n_450k %>% 
  ggplot(aes(x = cluster, fill = Nfenqi)) +  
  geom_bar(position = position_fill()) + 
  scale_fill_brewer(palette = 'Set3') +  
  theme_classic() +  
  labs(y = 'Percent') +  
  coord_flip()  
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\PAAD_cluster+N100_plot.pdf",aaa,width = 10,height = 10)




#stage#
library(tidytree)
colnames(label)[1]<-"Sample.ID"
label$Sample.ID<-substr(label$Sample.ID,1,15)
brca_stage_450<-left_join(cluster_sample,label,by='Sample.ID')
a<-brca_stage_450$label
a[which(a=="1")]<-"stage I"
a[which(a=="2"|a=="3")]<-"stage II&III"
brca_stage_450$stage<-a
colnames(brca_stage_450)[2]<-"cluster"
brca_stage_450$cluster<-as.character(brca_stage_450$cluster)
aaa<-brca_stage_450 %>% 
  ggplot(aes(x = cluster, fill = stage)) +  
  geom_bar(position = position_fill()) + 
  scale_fill_brewer(palette = 'Set3') +  
  theme_classic() +  
  labs(y = 'Percent') +  
  coord_flip()  
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\PAAD_cluster+stage100_plot.pdf",aaa,width = 10,height = 10)

 
library(data.table)
LUADsurvival=fread(paste0("G:\\dnameth\\PAAD","\\","TCGA-","PAAD",".GDC_phenotype.tsv.gz"),header=TRUE)
LUADsurvival<-as.data.frame(LUADsurvival)
rownames(LUADsurvival)<-LUADsurvival[,1]
LUADsurvival<-LUADsurvival[,-1]
rownames(LUADsurvival)=gsub("-","\\.",rownames(LUADsurvival))
unique(LUADsurvival$clinical_stage)
unique(LUADsurvival$tumor_stage)
str(LUADsurvival$clinical_stage)

matches_stage_i <-grepl("\\bStage I(A|B)?\\b", LUADsurvival$clinical_stage)
matches_stage_ii <-grepl("\\bStage II(A|B|C)?\\b", LUADsurvival$clinical_stage)
matches_stage_iii <-grepl("\\bStage III(A|B|C)?\\b", LUADsurvival$clinical_stage)

matches_stage_i <-grepl("\\bstage i(a|b)?\\b", LUADsurvival$tumor_stage)
matches_stage_ii <-grepl("\\bstage ii(a|b|c)?\\b", LUADsurvival$tumor_stage)
matches_stage_iii <-grepl("\\bstage iii(a|b|c)?\\b", LUADsurvival$tumor_stage)
LUADsurvival$tumor_stage[matches_stage_i]<-"1"
LUADsurvival$tumor_stage[matches_stage_ii]<-"2"
LUADsurvival$tumor_stage[matches_stage_iii]<-"3"
LUADstage1=LUADsurvival[LUADsurvival$tumor_stage==1,]
LUADstage2=LUADsurvival[LUADsurvival$tumor_stage==2,]
LUADstage3=LUADsurvival[LUADsurvival$tumor_stage==3,]

labeltostage<-data.table(Sample.ID=c(rownames(LUADstage1),rownames(LUADstage2),rownames(LUADstage3)),stage=c(rep(c("Stage I","Stage II","Stage III"),c(nrow(LUADstage1),nrow(LUADstage2),nrow(LUADstage3)))))


# labeltostage<-read.csv("label2 stage.csv")
library(tidytree)
colnames(label)[1]<-"Sample.ID"
colnames(labeltostage)[1]<-"Sample.ID"
labeltostage$Sample.ID<-substr(labeltostage$Sample.ID,1,15)

# brca_stage_450<-left_join(cluster_sample,label,by='Sample.ID')
brca_stage_450<-left_join(brca_stage_450,labeltostage,by='Sample.ID')
# brca_stage_450<-brca_stage_450[,-c(4,5)]
colnames(brca_stage_450)[5]<-"stage"
colnames(brca_stage_450)[2]<-"cluster"
# colnames(brca_stage_450)[3]<-"stage"
# brca_stage_450$stage<-paste("stage",brca_stage_450$stage)
brca_stage_450$cluster<-paste("Cluster",brca_stage_450$cluster)

#PAAD
# brca_stage_450$stage[1:162]<-"Stage I"
# brca_stage_450$stage[163:231]<-"Stage II"
# brca_stage_450$stage[232:277]<-"Stage III"
# brca_stage_450$stage[which(brca_stage_450$stage=="stage 2")]<-brca_stage_450$stage2[which(brca_stage_450$stage=="stage 2")]
# brca_stage_450$stage[which(brca_stage_450$stage=="stage 1")]<-"Stage I"
# brca_stage_450$stage[which(brca_stage_450$stage=="STAGE II")]<-"Stage II"
# brca_stage_450$stage[which(brca_stage_450$stage=="STAGE III")]<-"Stage III"
# devtools::install_github("davidsjoberg/ggsankey")
library(tidyverse)
library(ggsankey)
library(ggplot2)
library(cols4all)
head(brca_stage_450)
df <- brca_stage_450 %>%
  make_long(cluster, stage)
 
df$node <- factor(df$node,levels = c(rev(unique(brca_stage_450$stage)),
                                     rev(unique(brca_stage_450$cluster))))
 
c4a_gui()
mycol <- c4a('rainbow_wh_rd',53)
mycol <- sample(mycol,length(mycol))  
mycol<-c("#FFD4A9", "#80D1C8","pink","#7FB974","#5072BF","#AE84BB")
 
p <- ggplot(df, aes(x = x, next_x = next_x,
                    node = node, next_node = next_node,
                    fill = node,
                    label = node)) +  
  scale_fill_manual(values = mycol)  

p1 <- p +
  geom_sankey(flow.alpha = 0.5,  
              smooth = 8,  
              width = 0.12) +  
  geom_sankey_text(size = 3.2,  
                   color = "black") +
  theme_void() +
  theme(legend.position = 'none')  
p1
ggsave("G:\\pan_result\\gpu_result_tcga\\PAAD\\450k\\PAAD_cluster+stage_sankey_plot.pdf",p1,width = 10,height = 10)


#THCA####

 
library(ggplot2)  
library(dplyr)  
#library(ggstatsplot)  

brca_tnm<-read.table("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\TCGA.THCA.sampleMap_THCA_clinicalMatrix",sep = "\t",header = T)
brca_tnm<-data.frame(Sample.ID=brca_tnm$sampleID,T=brca_tnm$pathologic_T,N=brca_tnm$pathologic_N,M=brca_tnm$pathologic_M)
sigmoid_450k<-data.frame(Sample.ID=rownames(exp),exp)
brca_tnm$Sample.ID<-gsub("-",".",brca_tnm$Sample.ID)
sigmoid_450k$Sample.ID<-substr(sigmoid_450k$Sample.ID,1,15)


brca_tnm_450k<-merge(brca_tnm,sigmoid_450k,by="Sample.ID")
brca_tnm_450k<-na.omit(brca_tnm_450k)
str(brca_tnm_450k)
table(brca_tnm_450k$T)
table(brca_tnm_450k$N)
table(brca_tnm_450k$M)

brca_tnm_450k$T<-gsub(".*T1.*", "T1", brca_tnm_450k$T)
brca_tnm_450k$T<-gsub(".*T2.*", "T2", brca_tnm_450k$T)
brca_tnm_450k$T<-gsub(".*T3.*", "T3", brca_tnm_450k$T)
brca_tnm_450k$T<-gsub(".*T4.*", "T4", brca_tnm_450k$T)
brca_tnm_450k$N<-gsub(".*N1.*", "N1", brca_tnm_450k$N)
brca_tnm_450k$N<-gsub(".*N2.*", "N2", brca_tnm_450k$N)
brca_tnm_450k<-brca_tnm_450k[-which(brca_tnm_450k$N==""),]
brca_tnm_450k<-brca_tnm_450k[-which(brca_tnm_450k$T=="TX"),]
brca_tnm_450k<-brca_tnm_450k[-which(brca_tnm_450k$N==""),]
brca_tnm_450k<-brca_tnm_450k[-which(brca_tnm_450k$N=="NX"),]


cluster_sample<-data.frame(Sample.ID=rownames(clusters),class=clusters$clusters)
cluster_sample$Sample.ID<-substr(cluster_sample$Sample.ID,1,15)
brca_tnm_450k
brca_tn_450k<-merge(brca_tnm_450k,cluster_sample,by="Sample.ID")
#T 
brca_t_450k<-data.frame(Tfenqi=brca_tn_450k$T,cluster=brca_tn_450k$class,brca_tn_450k[,5:(ncol(brca_tnm_450k)-1)])
brca_t_450k<-na.omit(brca_t_450k)
brca_t_450k$cluster<-as.character(brca_t_450k$cluster)
aaa<-brca_t_450k %>% 
  ggplot(aes(x = cluster, fill = Tfenqi)) +  
  geom_bar(position = position_fill()) + 
  scale_fill_brewer(palette = 'Set3') +  
  theme_classic() +  
  labs(y = 'Percent') +  
  coord_flip()  
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\THCA_cluster+T100_plot.pdf",aaa,width = 10,height = 10)



#N 
brca_n_450k<-data.frame(Nfenqi=brca_tn_450k$N,cluster=brca_tn_450k$class,brca_tn_450k[,5:(ncol(brca_tnm_450k)-1)])
brca_n_450k<-na.omit(brca_n_450k)
table(brca_n_450k$Nfenqi)
# brca_n_450k<-brca_n_450k[-which(brca_n_450k$Nfenqi=="NX"),]
brca_n_450k$cluster<-as.character(brca_n_450k$cluster)
aaa<-brca_n_450k %>% 
  ggplot(aes(x = cluster, fill = Nfenqi)) +  
  geom_bar(position = position_fill()) + 
  scale_fill_brewer(palette = 'Set3') +  
  theme_classic() +  
  labs(y = 'Percent') +  
  coord_flip()  
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\THCA_cluster+N100_plot.pdf",aaa,width = 10,height = 10)




#stage#
library(tidytree)
colnames(label)[1]<-"Sample.ID"
label$Sample.ID<-substr(label$Sample.ID,1,15)
brca_stage_450<-left_join(cluster_sample,label,by='Sample.ID')
a<-brca_stage_450$label
a[which(a=="1")]<-"stage I"
a[which(a=="2"|a=="3")]<-"stage II&III"
brca_stage_450$stage<-a
colnames(brca_stage_450)[2]<-"cluster"
brca_stage_450$cluster<-as.character(brca_stage_450$cluster)
aaa<-brca_stage_450 %>% 
  ggplot(aes(x = cluster, fill = stage)) +  
  geom_bar(position = position_fill()) + 
  scale_fill_brewer(palette = 'Set3') +  
  theme_classic() +  
  labs(y = 'Percent') +  
  coord_flip()  
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\THCA_cluster+stage100_plot.pdf",aaa,width = 10,height = 10)

 
library(data.table)
LUADsurvival=fread(paste0("G:\\dnameth\\THCA","\\","TCGA-","THCA",".GDC_phenotype.tsv.gz"),header=TRUE)
LUADsurvival<-as.data.frame(LUADsurvival)
rownames(LUADsurvival)<-LUADsurvival[,1]
LUADsurvival<-LUADsurvival[,-1]
rownames(LUADsurvival)=gsub("-","\\.",rownames(LUADsurvival))
unique(LUADsurvival$clinical_stage)
unique(LUADsurvival$tumor_stage)
str(LUADsurvival$clinical_stage)

matches_stage_i <-grepl("\\bStage I(A|B)?\\b", LUADsurvival$clinical_stage)
matches_stage_ii <-grepl("\\bStage II(A|B|C)?\\b", LUADsurvival$clinical_stage)
matches_stage_iii <-grepl("\\bStage III(A|B|C)?\\b", LUADsurvival$clinical_stage)

matches_stage_i <-grepl("\\bstage i(a|b)?\\b", LUADsurvival$tumor_stage)
matches_stage_ii <-grepl("\\bstage ii(a|b|c)?\\b", LUADsurvival$tumor_stage)
matches_stage_iii <-grepl("\\bstage iii(a|b|c)?\\b", LUADsurvival$tumor_stage)
LUADsurvival$tumor_stage[matches_stage_i]<-"1"
LUADsurvival$tumor_stage[matches_stage_ii]<-"2"
LUADsurvival$tumor_stage[matches_stage_iii]<-"3"
LUADstage1=LUADsurvival[LUADsurvival$tumor_stage==1,]
LUADstage2=LUADsurvival[LUADsurvival$tumor_stage==2,]
LUADstage3=LUADsurvival[LUADsurvival$tumor_stage==3,]

labeltostage<-data.table(Sample.ID=c(rownames(LUADstage1),rownames(LUADstage2),rownames(LUADstage3)),stage=c(rep(c("Stage I","Stage II","Stage III"),c(nrow(LUADstage1),nrow(LUADstage2),nrow(LUADstage3)))))


# labeltostage<-read.csv("label2stage.csv")
library(tidytree)
colnames(label)[1]<-"Sample.ID"
colnames(labeltostage)[1]<-"Sample.ID"
labeltostage$Sample.ID<-substr(labeltostage$Sample.ID,1,15)

# brca_stage_450<-left_join(cluster_sample,label,by='Sample.ID')
brca_stage_450<-left_join(brca_stage_450,labeltostage,by='Sample.ID')
# brca_stage_450<-brca_stage_450[,-c(4,5)]
colnames(brca_stage_450)[5]<-"stage"
colnames(brca_stage_450)[2]<-"cluster"
# colnames(brca_stage_450)[3]<-"stage"
# brca_stage_450$stage<-paste("stage",brca_stage_450$stage)
brca_stage_450$cluster<-paste("Cluster",brca_stage_450$cluster)

# THCA
# brca_stage_450$stage[1:162]<-"Stage I"
# brca_stage_450$stage[163:231]<-"Stage II"
# brca_stage_450$stage[232:277]<-"Stage III"
# brca_stage_450$stage[which(brca_stage_450$stage=="stage 2")]<-brca_stage_450$stage2[which(brca_stage_450$stage=="stage 2")]
# brca_stage_450$stage[which(brca_stage_450$stage=="stage 1")]<-"Stage I"
# brca_stage_450$stage[which(brca_stage_450$stage=="STAGE II")]<-"Stage II"
# brca_stage_450$stage[which(brca_stage_450$stage=="STAGE III")]<-"Stage III"
# devtools::install_github("davidsjoberg/ggsankey")
library(tidyverse)
library(ggsankey)
library(ggplot2)
library(cols4all)
head(brca_stage_450)
df <- brca_stage_450 %>%
  make_long(cluster, stage)
 
df$node <- factor(df$node,levels = c(rev(unique(brca_stage_450$stage)),
                                     rev(unique(brca_stage_450$cluster))))
 
c4a_gui()
mycol <- c4a('rainbow_wh_rd',53)
mycol <- sample(mycol,length(mycol))  
mycol<-c("#FFD4A9", "#80D1C8","pink","#7FB974","#5072BF","#AE84BB")
 
p <- ggplot(df, aes(x = x, next_x = next_x,
                    node = node, next_node = next_node,
                    fill = node,
                    label = node)) +  
  scale_fill_manual(values = mycol)  

p1 <- p +
  geom_sankey(flow.alpha = 0.5,  
              smooth = 8,  
              width = 0.12) +  
  geom_sankey_text(size = 3.2,  
                   color = "black") +
  theme_void() +
  theme(legend.position = 'none')  
p1
ggsave("G:\\pan_result\\gpu_result_tcga\\THCA\\450k\\THCA_cluster+stage_sankey_plot.pdf",p1,width = 10,height = 10)


#UCEC####

 
library(ggplot2)  
library(dplyr)  
#library(ggstatsplot)  

brca_tnm<-read.table("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\TCGA.UCEC.sampleMap_UCEC_clinicalMatrix",sep = "\t",header = T)
brca_tnm<-data.frame(Sample.ID=brca_tnm$sampleID,T=brca_tnm$pathologic_T,N=brca_tnm$pathologic_N,M=brca_tnm$pathologic_M)
sigmoid_450k<-data.frame(Sample.ID=rownames(exp),exp)
brca_tnm$Sample.ID<-gsub("-",".",brca_tnm$Sample.ID)
sigmoid_450k$Sample.ID<-substr(sigmoid_450k$Sample.ID,1,15)


brca_tnm_450k<-merge(brca_tnm,sigmoid_450k,by="Sample.ID")
brca_tnm_450k<-na.omit(brca_tnm_450k)
str(brca_tnm_450k)
table(brca_tnm_450k$T)
table(brca_tnm_450k$N)
table(brca_tnm_450k$M)

brca_tnm_450k$T<-gsub(".*T1.*", "T1", brca_tnm_450k$T)
brca_tnm_450k$T<-gsub(".*T2.*", "T2", brca_tnm_450k$T)
brca_tnm_450k$T<-gsub(".*T3.*", "T3", brca_tnm_450k$T)
brca_tnm_450k$T<-gsub(".*T4.*", "T4", brca_tnm_450k$T)
brca_tnm_450k$N<-gsub(".*N1.*", "N1", brca_tnm_450k$N)
brca_tnm_450k$N<-gsub(".*N2.*", "N2", brca_tnm_450k$N)
brca_tnm_450k<-brca_tnm_450k[-which(brca_tnm_450k$N==""),]
brca_tnm_450k<-brca_tnm_450k[-which(brca_tnm_450k$T=="TX"),]
brca_tnm_450k<-brca_tnm_450k[-which(brca_tnm_450k$N==""),]
brca_tnm_450k<-brca_tnm_450k[-which(brca_tnm_450k$N=="NX"),]


cluster_sample<-data.frame(Sample.ID=rownames(clusters),class=clusters$clusters)
cluster_sample$Sample.ID<-substr(cluster_sample$Sample.ID,1,15)
brca_tnm_450k
brca_tn_450k<-merge(brca_tnm_450k,cluster_sample,by="Sample.ID")
#T 
brca_t_450k<-data.frame(Tfenqi=brca_tn_450k$T,cluster=brca_tn_450k$class,brca_tn_450k[,5:(ncol(brca_tnm_450k)-1)])
brca_t_450k<-na.omit(brca_t_450k)
brca_t_450k$cluster<-as.character(brca_t_450k$cluster)
aaa<-brca_t_450k %>% 
  ggplot(aes(x = cluster, fill = Tfenqi)) +  
  geom_bar(position = position_fill()) + 
  scale_fill_brewer(palette = 'Set3') +  
  theme_classic() +  
  labs(y = 'Percent') +  
  coord_flip()  
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\UCEC_cluster+T100_plot.pdf",aaa,width = 10,height = 10)



#N 
brca_n_450k<-data.frame(Nfenqi=brca_tn_450k$N,cluster=brca_tn_450k$class,brca_tn_450k[,5:(ncol(brca_tnm_450k)-1)])
brca_n_450k<-na.omit(brca_n_450k)
table(brca_n_450k$Nfenqi)
# brca_n_450k<-brca_n_450k[-which(brca_n_450k$Nfenqi=="NX"),]
brca_n_450k$cluster<-as.character(brca_n_450k$cluster)
aaa<-brca_n_450k %>% 
  ggplot(aes(x = cluster, fill = Nfenqi)) +  
  geom_bar(position = position_fill()) + 
  scale_fill_brewer(palette = 'Set3') +  
  theme_classic() +  
  labs(y = 'Percent') +  
  coord_flip()  
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\UCEC_cluster+N100_plot.pdf",aaa,width = 10,height = 10)




#stage#
library(tidytree)
colnames(label)[1]<-"Sample.ID"
label$Sample.ID<-substr(label$Sample.ID,1,15)
brca_stage_450<-left_join(cluster_sample,label,by='Sample.ID')
a<-brca_stage_450$label
a[which(a=="1")]<-"stage I"
a[which(a=="2"|a=="3")]<-"stage II&III"
brca_stage_450$stage<-a
colnames(brca_stage_450)[2]<-"cluster"
brca_stage_450$cluster<-as.character(brca_stage_450$cluster)
aaa<-brca_stage_450 %>% 
  ggplot(aes(x = cluster, fill = stage)) +  
  geom_bar(position = position_fill()) + 
  scale_fill_brewer(palette = 'Set3') +  
  theme_classic() +  
  labs(y = 'Percent') +  
  coord_flip()  
aaa
ggsave("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\UCEC_cluster+stage100_plot.pdf",aaa,width = 10,height = 10)


library(data.table)
LUADsurvival=fread(paste0("G:\\dnameth\\UCEC","\\","TCGA-","UCEC",".GDC_phenotype.tsv.gz"),header=TRUE)
LUADsurvival<-as.data.frame(LUADsurvival)
rownames(LUADsurvival)<-LUADsurvival[,1]
LUADsurvival<-LUADsurvival[,-1]
rownames(LUADsurvival)=gsub("-","\\.",rownames(LUADsurvival))
unique(LUADsurvival$clinical_stage)
unique(LUADsurvival$tumor_stage)
str(LUADsurvival$clinical_stage)

matches_stage_i <-grepl("\\bStage I(A|B)?\\b", LUADsurvival$clinical_stage)
matches_stage_ii <-grepl("\\bStage II(A|B|C)?\\b", LUADsurvival$clinical_stage)
matches_stage_iii <-grepl("\\bStage III(A|B|C)?\\b", LUADsurvival$clinical_stage)

matches_stage_i <-grepl("\\bstage i(a|b)?\\b", LUADsurvival$tumor_stage)
matches_stage_ii <-grepl("\\bstage ii(a|b|c)?\\b", LUADsurvival$tumor_stage)
matches_stage_iii <-grepl("\\bstage iii(a|b|c)?\\b", LUADsurvival$tumor_stage)
LUADsurvival$tumor_stage[matches_stage_i]<-"1"
LUADsurvival$tumor_stage[matches_stage_ii]<-"2"
LUADsurvival$tumor_stage[matches_stage_iii]<-"3"
LUADstage1=LUADsurvival[LUADsurvival$tumor_stage==1,]
LUADstage2=LUADsurvival[LUADsurvival$tumor_stage==2,]
LUADstage3=LUADsurvival[LUADsurvival$tumor_stage==3,]

labeltostage<-data.table(Sample.ID=c(rownames(LUADstage1),rownames(LUADstage2),rownames(LUADstage3)),stage=c(rep(c("Stage I","Stage II","Stage III"),c(nrow(LUADstage1),nrow(LUADstage2),nrow(LUADstage3)))))


# labeltostage<-read.csv("label2 stage.csv")
library(tidytree)
colnames(label)[1]<-"Sample.ID"
colnames(labeltostage)[1]<-"Sample.ID"
labeltostage$Sample.ID<-substr(labeltostage$Sample.ID,1,15)

# brca_stage_450<-left_join(cluster_sample,label,by='Sample.ID')
brca_stage_450<-left_join(brca_stage_450,labeltostage,by='Sample.ID')
# brca_stage_450<-brca_stage_450[,-c(4,5)]
colnames(brca_stage_450)[5]<-"stage"
colnames(brca_stage_450)[2]<-"cluster"
# colnames(brca_stage_450)[3]<-"stage"
# brca_stage_450$stage<-paste("stage",brca_stage_450$stage)
brca_stage_450$cluster<-paste("Cluster",brca_stage_450$cluster)

# UCEC
# brca_stage_450$stage[1:162]<-"Stage I"
# brca_stage_450$stage[163:231]<-"Stage II"
# brca_stage_450$stage[232:277]<-"Stage III"
# brca_stage_450$stage[which(brca_stage_450$stage=="stage 2")]<-brca_stage_450$stage2[which(brca_stage_450$stage=="stage 2")]
# brca_stage_450$stage[which(brca_stage_450$stage=="stage 1")]<-"Stage I"
# brca_stage_450$stage[which(brca_stage_450$stage=="STAGE II")]<-"Stage II"
# brca_stage_450$stage[which(brca_stage_450$stage=="STAGE III")]<-"Stage III"
# devtools::install_github("davidsjoberg/ggsankey")
library(tidyverse)
library(ggsankey)
library(ggplot2)
library(cols4all)
head(brca_stage_450)
df <- brca_stage_450 %>%
  make_long(cluster, stage)
 
df$node <- factor(df$node,levels = c(rev(unique(brca_stage_450$stage)),
                                     rev(unique(brca_stage_450$cluster))))
 
c4a_gui()
mycol <- c4a('rainbow_wh_rd',53)
mycol <- sample(mycol,length(mycol))  
mycol<-c("#FFD4A9", "#80D1C8","pink","#7FB974","#5072BF","#AE84BB")
 
p <- ggplot(df, aes(x = x, next_x = next_x,
                    node = node, next_node = next_node,
                    fill = node,
                    label = node)) +  
  scale_fill_manual(values = mycol)  

p1 <- p +
  geom_sankey(flow.alpha = 0.5,  
              smooth = 8,  
              width = 0.12) +  
  geom_sankey_text(size = 3.2,  
                   color = "black") +
  theme_void() +
  theme(legend.position = 'none') 
p1
ggsave("G:\\pan_result\\gpu_result_tcga\\UCEC\\450k\\UCEC_cluster+stage_sankey_plot.pdf",p1,width = 10,height = 10)


