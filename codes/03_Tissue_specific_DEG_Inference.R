library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(GOSemSim)
library(DOSE)

immune.combined <- subset(immune.combined, subset= stim %in% c("BM","mPB","PB","thymus"))
sub1 <- subset(immune.combined, subset= stim %in% c("BM","mPB","PB"))
##EryMk lineage
sub1 <- subset(sub1, subset= celltype %in% c("HSC","MPP#1","MPP#2","LMPP#1","LMPP#2","BMEP","EryP","MkP","Ma/Eo/BaP"))
##B lineage
sub1 <- subset(sub1, subset= celltype %in% c("HSC","MPP#1","MPP#2","LMPP#1","LMPP#2","CLP","pro-B","pre-B"))
##T lineage
sub1 <- subset(sub1, subset= celltype %in% c("HSC","MPP#1","MPP#2","LMPP#1","LMPP#2","ETP-like","Thy#1","Thy#2","Thy#3"))
##My lineage
sub1 <- subset(sub1, subset= celltype %in% c("HSC","MPP#1","MPP#2","LMPP#1","LMPP#2","GMP","CDP","pre-pDC"))
sub1 <- subset(immune.combined, subset= ((celltype=="HSC"|celltype=="MPP #1"|celltype=="MPP #2"|celltype=="LMPP #1"|celltype=="GMP"|celltype=="pre-cDC")&(stim=="BM"|stim=="PB"|stim=="mPB"))
               |celltype=="pre-pDC")
##BM mPB PB all lineages
sub1 <- subset(sub1, subset= celltype %in% c("HSC","MPP#1","MPP#2","LMPP#1","LMPP#2","BMEP","EryP","MkP","Ma/Eo/BaP","GMP","CDP","pre-pDC","CLP","pro-B","pre-B"))
sub1 <- subset(sub1_all, subset= stim %in% c("BM","mPB"))
sub1 <- subset(sub1, subset= stim %in% c("BM","mPB","PB"))
sub1$tissue_celltype <- paste(sub1$stim, sub1$celltype, sep="_")
SaveH5Seurat(sub1,filename="/home/yushiya/data/cd34/data/fig/immune_BMPB.h5seurat", overwrite = TRUE)
Convert("/home/yushiya/data/cd34/data/fig/immune_BMPB.h5seurat", dest = "h5ad", overwrite = TRUE)


#plot of cluster percentage changing
library(reshape2)
library("ggalluvial")
#total results with x=stim fill=celltype
clus_percent <- table(immune.combined@meta.data[["stim"]],immune.combined@meta.data[["celltype"]])
clus_percent <- apply(clus_percent,1,function(x) prop.table(x))
clus_percent <- t(clus_percent[order(as.numeric(rownames(clus_percent))),])
Clus=colnames(clus_percent)
clus_percent=data.frame(t(clus_percent), Clus)
clus_percent=melt(clus_percent, id='Clus')
names(clus_percent)[2]='Stim'
clus_percent$Clus <- factor(clus_percent$Clus,levels = c('HSC','MPP#1',"MPP#2",'LMPP#1',"LMPP#2",'BMEP','EryP','MkP',"Ma/Eo/BaP",'GMP','CDP','pre-pDC','CLP','pro-B','pre-B',"ETP-like","Thy#1","Thy#2","Thy#3"))
p1=ggplot(clus_percent,
          aes(x=Stim, y=value*100, fill=Clus, stratum = Clus, alluvium = Clus)) +
  geom_bar(stat='identity', width=0.45) +
  geom_alluvium() +
  geom_stratum(width=0.45, size=0.1) +
  scale_fill_manual(values = color_all) +
  labs(x='Samples', y='Relative Abundance (%)')+
  #scale_y_continuous(expand=c(0, 0))+
  theme(axis.text.x=element_text(angle=45, hjust=1),
        panel.border = element_rect(color = "black", fill=NA),
        text = element_text(size = 12),
        panel.background = element_blank(),
        panel.grid = element_line(colour = "grey90"))# + theme_bw()
ggsave("/home/yushiya/data/cd34/data/fig/fig1-clusP.pdf", plot=p1, width=5, height=4.5)
#lineage results with x=celltype fill=stim
clus_percent <- table(sub1@meta.data[["celltype"]],sub1@meta.data[["stim"]])
clus_percent <- apply(clus_percent,2,function(x) prop.table(x))
clus_percent <- apply(clus_percent,1,function(x) prop.table(x))
clus_percent <- t(clus_percent[order(as.numeric(rownames(clus_percent))),])
clus_percent <- na.omit(clus_percent)
Tissue=colnames(clus_percent)
clus_percent=melt(clus_percent, id='Tissue')
names(clus_percent)[2]='Tissue'
names(clus_percent)[1]='Celltype'
clus_percent$Celltype <- factor(clus_percent$Celltype,levels = c('HSC','MPP#1',"MPP#2",'LMPP#1',"LMPP#2",'BMEP','EryP','MkP',"Ma/Eo/BaP",'GMP','CDP','pre-pDC','CLP','pro-B','pre-B'))
clus_percent$Celltype <- factor(clus_percent$Celltype,levels = c('HSC','MPP#1',"MPP#2",'LMPP#1',"LMPP#2",'BMEP','EryP','MkP',"Ma/Eo/BaP"))
clus_percent$Celltype <- factor(clus_percent$Celltype,levels = c('HSC','MPP#1',"MPP#2",'LMPP#1',"LMPP#2",'GMP','CDP','pre-pDC'))
clus_percent$Celltype <- factor(clus_percent$Celltype,levels = c('HSC','MPP#1',"MPP#2",'LMPP#1','LMPP#2','ETP-like',"Thy#1","Thy#2","Thy#3"))
clus_percent$Celltype <- factor(clus_percent$Celltype,levels = c('HSC','MPP#1',"MPP#2",'LMPP#1',"LMPP#2",'CLP','pro-B','pre-B'))
p1=ggplot(clus_percent,
          aes(x=Celltype, y=value*100, fill=Tissue)) +
  geom_bar(stat='identity', width=0.45) +
  #geom_alluvium() +
  #geom_stratum(width=0.45, size=0.1) +
  scale_fill_manual(values = color_stim) +
  labs(x='Celltypes', y='Relative Abundance (%)')+
  scale_y_continuous(expand=c(0, 0))+
  theme(axis.text.x=element_text(hjust=1, angle=45, vjust=1),
        text = element_text(size = 12),
        panel.background = element_blank(),
        panel.grid = element_line(colour = "grey90"),
        axis.line = element_line(colour = "grey90",size = rel(1)))# + theme_bw()
ggsave("/home/yushiya/data/cd34/data/fig/fig2-My-clusP.pdf", plot=p1, width=4, height=2.5)
ggsave("/home/yushiya/data/cd34/data/fig/fig2-clusP-BMmPBPB-B.pdf", plot=p1, width=5, height=2.5)

#plot of cell cycle Phase percentage changing
#total results with x=Phase, fill=celltype
phase_per <- table(immune.combined@meta.data[["celltype"]],immune.combined@meta.data[["Phase"]])
phase_per <- apply(phase_per,2,function(x) prop.table(x))
Clus=colnames(phase_per)
phase_per=melt(phase_per, id='Clus')
names(phase_per)[1]='Celltype'
names(phase_per)[2]='Phase'
phase_per$Phase <- factor(phase_per$Phase, levels=c("G1","S","G2M"))
phase_per$Celltype <- factor(phase_per$Celltype,levels = c('HSC','MPP #1',"MPP #2",'BMEP','EryP','MkP','LMPP #1','LMPP #2','GMP','pre-cDC','pre-pDC',"ETP","mThy","cThy #1","cThy #2",'pre-pro B','pro-B','pre-B'))
p1=ggplot(phase_per,
          aes(x=Phase, y=value*100, fill=Celltype, stratum = Celltype, alluvium = Celltype)) +
  geom_bar(stat='identity', width=0.45) +
  geom_alluvium() +
  geom_stratum(width=0.45, size=0.1) +
  scale_fill_manual(values = color_all) +
  labs(x='Phase', y='Relative Abundance (%)')+
  scale_y_continuous(expand=c(0, 0))+
  theme(axis.text.x=element_text(hjust=0.5),
        text = element_text(size = 18),
        panel.background = element_blank(),
        panel.grid = element_line(colour = "grey90"))# + theme_bw()
ggsave("/home/yushiya/data/cd34/manu_fig/v4/fig1-cycleP.pdf", plot=p1, width=6, height=6)
#lineage results with x=Celltype, fill=Phase
phase_per <- table(sub1@meta.data[["Phase"]],sub1@meta.data[["celltype"]])
phase_per <- apply(phase_per,2,function(x) prop.table(x))
Clus=colnames(phase_per)
phase_per=melt(phase_per, id='Clus')
phase_per <- na.omit(phase_per)
names(phase_per)[2]='Celltype'
names(phase_per)[1]='Phase'
phase_per$Phase <- factor(phase_per$Phase, levels=c("G1","S","G2M"))
phase_per$Celltype <- factor(phase_per$Celltype,levels = c('HSC','MPP #1',"MPP #2",'LMPP #1','LMPP #2'))
phase_per$Celltype <- factor(phase_per$Celltype,levels = c('HSC','MPP #1',"MPP #2",'MEP','EryP','MkP'))
phase_per$Celltype <- factor(phase_per$Celltype,levels = c('HSC','MPP #1',"MPP #2",'LMPP #2','CLP','pro-B',"pre-B"))
phase_per$Celltype <- factor(phase_per$Celltype,levels = c('HSC','MPP #1',"MPP #2",'LMPP #1','ETP','mThy',"cThy #1","cThy #2"))
phase_per$Celltype <- factor(phase_per$Celltype,levels = c('HSC','MPP #1',"MPP #2",'LMPP #1','GMP','CDP',"pre-pDC"))
p1=ggplot(phase_per,
          aes(x=Celltype, y=value*100, fill=Phase)) +
  geom_bar(stat='identity', width=0.45) +
  #geom_alluvium() +
  #geom_stratum(width=0.45, size=0.1) +
  scale_fill_manual(values = color_cycle) +
  labs(x='Celltypes', y='Relative Abundance (%)')+
  scale_y_continuous(expand=c(0, 0))+
  theme(axis.text.x=element_text(hjust=0.5, angle=45, vjust=0.5),
        text = element_text(size = 12),
        panel.background = element_blank(),
        panel.grid = element_line(colour = "grey90"),
        axis.line = element_line(colour = "grey90",size = rel(1)))# + theme_bw()
ggsave("/home/yushiya/cd34/fig/fig-My-cycleP.pdf", plot=p1, width=4, height=2.5)

#cell density via pseudotime
#density plot
sub1 <- subset(sub1, subset= stim=="BM"|stim=="PB"|stim=="mPB")
type_dis<-as.data.frame(colnames(sub1))
type_dis$Tissue <- sub1$stim
type_dis$Pseudotime <- sub1$pseudotime
names(type_dis)[1]='ID'
p1=ggplot(type_dis, aes(x = Pseudotime, fill = Tissue))+ geom_density(alpha=0.7, size=0.3)+
  theme_classic()+
  scale_fill_manual(values=color_stim %in% c("BM","PB"))+
  scale_x_continuous(limits = c(0,30))+
  theme(axis.text.x=element_text(hjust=0.5),
        text = element_text(size = 16),
        panel.background = element_blank(),
        panel.grid = element_line(colour = "black"),
        axis.line = element_line(colour = "black",size = rel(1),arrow = arrow(angle = 30,length = unit(0.1,"inches"))))
ggsave(file="/home/yushiya/fig/fig4-density-T.pdf", plot=p1, width=6, height=3)
p1=ggplot(type_dis, aes(x = Pseudotime, color=Tissue))+ geom_density(size=1)+
  theme_classic()+
  scale_color_manual(values=color_stim)+
  scale_x_continuous(limits = c(0,30))+
  theme(axis.text.x=element_text(hjust=0.5),
        text = element_text(size = 16),
        panel.background = element_blank(),
        panel.grid = element_line(colour = "black"),
        axis.line = element_line(colour = "black",size = rel(1),arrow = arrow(angle = 30,length = unit(0.1,"inches"))))
ggsave(file="/home/yushiya/data/cd34/data/fig/fig4-T-density.pdf", plot=p1, width=7, height=3)
ggsave(file="/home/yushiya/data/cd34/data/fig/fig2-density-BMmPBPB.pdf", plot=p1, width=8, height=3)
ggsave(file="/home/yushiya/data/cd34/data/fig/fig2-density-BMmPBPB-Mye.pdf", plot=p1, width=8, height=3)

#violin plot
p <- VlnPlot(sub1, features = "pseudotime", pt.size = 0, split.by = "stim",cols=color_stim,combine = FALSE,y.max=30)
p <- VlnPlot(sub1, features = "pseudotime", pt.size = 0, cols=color_all,combine = FALSE,y.max=30)
p1=p[[1]] + coord_flip() 
ggsave(file="/home/yushiya/data/cd34/manu_fig/v4/fig1-Vlndensity-my.pdf", plot=p1, width=6, height=6)

#calculate DEGs and GOs among stims in each celltype
library(pheatmap)
library(clusterProfiler, lib.loc='/home/panyuwen/miniconda3/envs/clusterprofiler/lib/R/library')
#library(org.Hs.eg.db, lib.loc='/home/panyuwen/miniconda3/envs/clusterprofiler/lib/R/library')
library(org.Hs.eg.db)
library(enrichplot, lib.loc='/home/panyuwen/miniconda3/envs/clusterprofiler/lib/R/library')
library(GOSemSim, lib.loc='/home/panyuwen/miniconda3/envs/clusterprofiler/lib/R/library')
library(DOSE, lib.loc='/home/panyuwen/miniconda3/envs/clusterprofiler/lib/R/library')
sub1 <- subset(immune.combined, subset= stim=="BM"|stim=="mPB"|stim=="PB")
Idents(sub1) <- sub1$stim
cp_stim.markers <- FindAllMarkers(sub1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.405)
write.csv(cp_stim.markers, file=paste("/home/yushiya/data/cd34/manu_fig/v4/fig2_DEG_BMmPBPB-1.csv",sep=""))
cp_stim.markers <- read.csv("/home/yushiya/data/cd34/manu_fig/v4/fig2_DEG_PBmPB-1.csv")
cp_stim.markers <- subset(cp_stim.markers, subset= p_val_adj<0.01)
sub1 <- subset(immune.combined, subset= stim=="BM"|stim=="mPB"|stim=="PB")
select="pre-pDC"
sub_select <- subset(sub1, subset= celltype==select)
Idents(sub_select) <- sub_select$stim
cp_stim.markers <- FindAllMarkers(sub_select, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1)
write.csv(cp_stim.markers, file=paste("/home/yushiya/data/cd34/manu_fig/v4/fig2_",select,".csv",sep=""))
top_genes <- cp_stim.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top_genes <- top_genes[,"gene"] %>% unique()
colnames(top_genes) <- "gene"
sub11<- subset(sub_select, subset= stim=="BM")
sub12<- subset(sub_select, subset= stim=="mPB")
sub13<- subset(sub_select, subset= stim=="PB")
sub14<- subset(sub_select, subset= stim=="thymus")
mean_exp <- apply(sub11@assays[["RNA"]]@scale.data[top_genes$gene,],1,mean)
top_genes$BM <- mean_exp
mean_exp <- apply(sub12@assays[["RNA"]]@scale.data[top_genes$gene,],1,mean)
top_genes$mPB <- mean_exp
mean_exp <- apply(sub13@assays[["RNA"]]@scale.data[top_genes$gene,],1,mean)
top_genes$PB <- mean_exp
mean_exp <- apply(sub14@assays[["RNA"]]@scale.data[top_genes$gene,],1,mean)
top_genes$thymus <- mean_exp
top_genes<-data.frame(top_genes)
rownames(top_genes) <- top_genes$gene
top_genes <- top_genes[,-1]
p1=pheatmap(top_genes, cluster_cols = F, #cluster_rows = F,  
            #clustering_method = "median",
            fontsize_row = 10, fontsize_col= 12, angle_col = 0,
            col=colorRampPalette(c("#000066","#339999" ,"yellow"))(50))
ggsave(file=paste("/home/yushiya/data/cd34/manu_fig/v4/fig2-DEG_",select,".pdf",sep=""), plot=p1, width=3.5, height=4)
BMID <- row.names(subset(cp_stim.markers, subset= cluster=="BM"))
BM_gene <- bitr(BMID, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
mPBID <- row.names(subset(cp_stim.markers, subset= cluster=="mPB"))
mPB_gene <- bitr(mPBID, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
PBID <- row.names(subset(cp_stim.markers, subset= cluster=="PB"))
PB_gene <- bitr(PBID, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
thymusID <- row.names(subset(cp_stim.markers, subset= cluster=="thymus"))
thymus_gene <- bitr(thymusID, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
#cp = list(BM.gene=BM_gene$ENTREZID, mPB_gene=mPB_gene$ENTREZID,PB_gene=PB_gene$ENTREZID)
cp = list(BM.gene=BM_gene$ENTREZID, mPB_gene=mPB_gene$ENTREZID,PB_gene=PB_gene$ENTREZID, thymus.gene=thymus_gene$ENTREZID)
go.p <- compareCluster(cp,
                       fun = "enrichGO",
                       OrgDb = "org.Hs.eg.db",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.01
)
go.p1 <- simplify(go.p,cutoff=0.4,by="p.adjust",select_fun=min)  #去除冗余
p1=ggplot(go.p1, aes(Cluster, Description), showCategory=6) +
  geom_point(aes(color=p.adjust, size=GeneRatio))+
  theme_classic()+
  theme(axis.text.x =element_text(angle = 45,vjust = 0.5,hjust = 0.5),
        text = element_text(size = 15),
        panel.background = element_blank(),
        panel.grid = element_line(colour = "grey90"))#+theme_bw()
ggsave(file=paste("/home/yushiya/data/cd34/manu_fig/v4/fig2-GO_",select,".pdf",sep=""), plot=p1, width=8, height=6)
p1=cnetplot(go.p) 
ggsave(file=paste("/home/yushiya/data/cd34/manu_fig/v4/fig2-cnet_",select,".pdf",sep=""), plot=p1, width=8, height=7)
