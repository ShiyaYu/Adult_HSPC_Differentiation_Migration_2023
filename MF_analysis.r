library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(devtools)
library(slingshot)
library(RColorBrewer)

MF <- readRDS("/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/MF.rds")
immune.combined <- readRDS("/home/yushiya/cd34/fig/HSPC.rds")
#map and annotate MF data to reference
immune.combined <- RunUMAP(immune.combined, dims = 1:40, reduction = "pca", return.model = TRUE)
all.anchors <- FindTransferAnchors(reference = immune.combined, query = MF,
                                        dims = 1:30, reference.reduction = "pca")
MF <- MapQuery(anchorset = all.anchors, reference = immune.combined, query = MF,
                           refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap")
MF$predicted.celltype <- factor(MF$predicted.celltype, levels=c("HSC","MPP #1","MPP #2","MEP","EryP","MkP","LMPP #1","LMPP #2","GMP","CDP","pre-pDC","ETP","Thy #1","Thy #2","CLP","pro-B"))
color_MF <- c("#550000", "#AA3939", "#CC6F66","#FFBD54","#B97802","#915900","#E79492","#FFF176","#CECA68", "#B2A13F","#DDA0DD","#98ADC4","#7F689D","#A5D6A7", "#47AF50")
color_stim <- c("#DB5C25","#F3B747","#649541","#9966CC")
color_stim <- c("#DB5C25","#649541","#9966CC")
p1=DimPlot(MF, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
        cols=color_MF,repel = TRUE,raster = T)
ggsave("/home/yushiya/cd34/fig/fig6-MF-umap.pdf", plot=p1, width=6, height=4)
MF$celltype <- MF$predicted.celltype
DefaultAssay(MF) <- "RNA"
p1=FeaturePlot(MF, reduction = "ref.umap", features = c("AVP","CRHBP","MPO","IRF8","GATA2","PPBP"), max.cutoff = 3, cols = c("#C5E9F5", "#981C12"), raster = T)
ggsave("/home/yushiya/cd34/fig/fig6-MF-Fea.pdf", plot=p1, width=6, height=7)

#merge MF and healthy samples
immune.combined$UMAP_1 <- immune.combined@reductions[["umap"]]@cell.embeddings[,1]
immune.combined$UMAP_2 <- immune.combined@reductions[["umap"]]@cell.embeddings[,2]
MF$UMAP_1 <- MF@reductions[["umap"]]@cell.embeddings[,1]
MF$UMAP_2 <- MF@reductions[["umap"]]@cell.embeddings[,2]
total <- merge(immune.combined,MF)
total$stim <- factor(total$stim, levels=c("BM","mPB","PB","MF"))
total$celltype <- factor(total$celltype, levels = c('HSC','MPP #1',"MPP #2",'MEP','EryP','MkP','LMPP #1','GMP','CDP',"pre-pDC","ETP","CLP","pro-B"))
Idents(total) <- total$celltype
total <- ScaleData(total, verbose = FALSE)
total <- RunPCA(total,npcs = 30, verbose = FALSE)
total <- FindNeighbors(total, reduction = "pca", dims = 1:30)
total <- FindClusters(total,resolution = 0.8)
total <- RunUMAP(total, reduction = "pca", dims = 1:30)
total@reductions[["umap"]]@cell.embeddings[,1] <- total$UMAP_1
total@reductions[["umap"]]@cell.embeddings[,2] <- total$UMAP_2
DimPlot(total, reduction = "umap", label = TRUE)
total <- subset(total, subset= stim=="BM"|stim=="PB"|stim=="mPB"|stim=="MF")
total <- subset(total, subset= stim=="BM"|stim=="PB"|stim=="MF")
total$stim <- factor(total$stim, levels=c("BM","PB","MF"))
total <- subset(total, subset= celltype=="HSC"|celltype=="MPP #1"|celltype=="MPP #2"|celltype=="MEP"|celltype=="EryP"|celltype=="MkP"|celltype=="LMPP #1"|celltype=="GMP"|celltype=="CDP"|celltype=="pre-pDC"|celltype=="ETP"|celltype=="CLP"|celltype=="pro-B")
saveRDS(total, file = "/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/MF_healthy.rds")
total <- readRDS(file = "/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/MF_healthy.rds")

#calculate Migratory score
library(ggpubr)
positive_regu_selected <- c("MPP1","MSN","CD99","RAC2","MAPK1","APP","PLCB1","ZNF580","TGFB1","CALR",
                            "TPGS1","C1QBP","CRK","PYCARD","MAPK3","SPN","ADAM10","RIN3","NCKAP1L","DNM1L",
                            "CD9","WNK1","VEGFB","DNAJC4","SPI1","MDK","RHOG","CD81","ADAM8",
                            "CAMK1D","ANXA1","DOCK8","LYN","TGS1","HOXA7","GPSM3","RIPOR2","CCL28","PTGER4",
                            "RHOH","CD47","SELENOK","INPP5D","ITGA4","ANXA4","RTN4","ADAM17","GCSAML","MIA3")
MF<-AddModuleScore(object = MF, features = list(Migratory_Score=positive_regu_selected), name="Migratory_Score")
VlnPlot(MF, features = c("Migratory_Score1"), group.by="predicted.celltype", pt.size = 0,combine = FALSE)
miGene <- as.data.frame(colnames(total))
miGene$stim <- total$stim
miGene$celltype <- total$celltype
miGene$celltype <- factor(miGene$celltype,levels = c('HSC','MPP #1',"MPP #2",'MEP','EryP','MkP','LMPP #1','GMP','CDP',"pre-pDC","ETP","CLP","pro-B"))
miGene$stim <- factor(miGene$stim, levels=c("BM","mPB","PB","MF"))
names(miGene)[1]='ID'
miGene$gene='Migratory_Score1'
miGene$Expression <- total$Migratory_Score1
p1=ggboxplot(miGene, x="celltype", y="Expression", color="stim", bxp.errorbar = T,
          palette = color_stim, add="boxplot", short.panel.labs = TRUE, ncol=4, add.params=list(size=0.2))+
  theme(axis.text.x =element_text(angle = 45,hjust = 1),
        plot.title = element_text(hjust = 0))
ggsave("/home/yushiya/cd34/fig/fig6-Migr-all.pdf", plot=p1, width=11, height=4)
#Expression of migration related genes
FeaturePlot(total,reduction = "umap", features = c("CXCR4"), max.cutoff = 2, cols = c("#C9E9FF","#981C12"), split.by = "stim", ncol=2) +
  theme(legend.position = "right")
VlnPlot(total, features = "CXCR4", pt.size = 0,combine = FALSE, split.by="stim", cols=color_stim)

#calculate cell type distribution
library(reshape2)
library("ggalluvial")
clus_percent <- table(total$celltype,total$stim)
clus_percent <- apply(clus_percent,1,function(x) prop.table(x))
clus_percent <- t(clus_percent[order(as.numeric(rownames(clus_percent))),])
clus_percent <- na.omit(clus_percent)
Stim=colnames(clus_percent)
clus_percent=melt(clus_percent, id='Stim')
names(clus_percent)[2]='Stim'
names(clus_percent)[1]='Celltype'
clus_percent$Celltype <- factor(clus_percent$Celltype,levels = c('HSC','MPP #1',"MPP #2",'MEP','EryP','MkP','LMPP #1','GMP','CDP',"pre-pDC","ETP","CLP","pro-B"))
clus_percent$Stim <- factor(clus_percent$Stim, levels=c("BM","mPB","PB","MF"))
p1=ggplot(clus_percent,
       aes(x=Celltype, y=value*100, fill=Stim)) +
  geom_bar(stat='identity', width=0.45) +
  #geom_alluvium() +
  #geom_stratum(width=0.45, size=0.1) +
  scale_fill_manual(values = color_stim) +
  labs(x='Celltypes', y='Relative Abundance (%)')+
  scale_y_continuous(expand=c(0, 0))+
  theme(axis.text.x=element_text(hjust=0.5, angle=45, vjust=0.5),
        text = element_text(size = 12),
        panel.background = element_blank(),
        panel.grid = element_line(colour = "grey90"),
        axis.line = element_line(colour = "grey90",size = rel(1)))
ggsave("/home/yushiya/cd34/fig/fig6-MF-clusP-all.pdf", plot=p1, width=6, height=3)

#calculate lineage scores
celltype_markers <- read.csv("/home/yushiya/cd34/fig/immune_rpca_celltype_markers.csv",stringsAsFactors = F)
LT_HSC_genes <- read_xlsx("/home/yushiya/cd34/fig/LT_HSC_genes.xlsx",col_names = T)
HSC_gene <- subset(celltype_markers, subset= cluster=="HSC")
HSC_gene <- HSC_gene[HSC_gene$gene %in% LT_HSC_genes$`Gene Symbol`,]
HSC_gene <- HSC_gene %>% top_n(n = 50, wt = avg_log2FC)
total<-AddModuleScore(object = total, features = list(HSC_Score=HSC_gene$gene), name="HSC_Score")
FeaturePlot(total,reduction = "umap", features = c("HSC_Score1"), max.cutoff = 2, cols = c("#C9E9FF","#981C12"), split.by = "stim") +
  theme(legend.position = "right")
p1=VlnPlot(total, features = "HSC_Score1", pt.size = 0,combine = FALSE, split.by="stim", cols=color_stim)
ggsave("/home/yushiya/cd34/fig/fig6-MF-HSC-score-all.pdf", plot=p1[[1]], width=7, height=3)
EryP_gene <- subset(celltype_markers, subset= cluster=="EryP")
EryP_gene <- EryP_gene %>% top_n(n = 50, wt = avg_log2FC)
total<-AddModuleScore(object = total, features = list(EryP_Score=EryP_gene$gene), name="EryP_Score")
p1=VlnPlot(total, features = "EryP_Score1", pt.size = 0,combine = FALSE, split.by="stim", cols=color_stim)
ggsave("/home/yushiya/cd34/fig/fig6-MF-Ery-score-all.pdf", plot=p1[[1]], width=7, height=3)
MkP_gene <- subset(celltype_markers, subset= cluster=="MkP")
MkP_gene <- MkP_gene %>% top_n(n = 50, wt = avg_log2FC)
total<-AddModuleScore(object = total, features = list(MkP_Score=MkP_gene$gene), name="MkP_Score")
p1=VlnPlot(total, features = "MkP_Score1", pt.size = 0,combine = FALSE, split.by="stim", cols=color_stim)
ggsave("/home/yushiya/cd34/fig/fig6-MF-Mk-score-all.pdf", plot=p1[[1]], width=7, height=3)
My_gene <- subset(celltype_markers, subset= cluster=="pre-cDC" | cluster=="pre-pDC")
My_gene <- My_gene %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
total<-AddModuleScore(object = total, features = list(My_Score=My_gene$gene), name="My_Score")
p1=VlnPlot(total, features = "My_Score1", pt.size = 0,combine = FALSE, split.by="stim", cols=color_stim)
ggsave("/home/yushiya/cd34/fig/fig6-MF-My-score-all.pdf", plot=p1[[1]], width=7, height=3)
B_gene <- subset(celltype_markers, subset= cluster=="pro-B"|cluster=="pre-B")
B_gene <- B_gene %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
total<-AddModuleScore(object = total, features = list(B_Score=B_gene$gene), name="B_Score")
p1=VlnPlot(total, features = "B_Score1", pt.size = 0,combine = FALSE, split.by="stim", cols=color_stim)
ggsave("/home/yushiya/cd34/fig/fig6-MF-B-score-all.pdf", plot=p1[[1]], width=7, height=3)
B_gene <- subset(celltype_markers, subset= cluster=="pro-B"|cluster=="pre-B")
B_gene <- B_gene %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
total<-AddModuleScore(object = total, features = list(B_Score=B_gene$gene), name="B_Score")
p1=VlnPlot(total, features = "B_Score1", pt.size = 0,combine = FALSE, split.by="stim", cols=color_stim)
ggsave("/home/yushiya/cd34/fig/fig6-MF-B-score-all.pdf", plot=p1[[1]], width=7, height=3)
T_gene <- subset(celltype_markers, subset= cluster=="Thy #1" | cluster=="Thy #2"| cluster=="Thy #3")
T_gene <- T_gene %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
total<-AddModuleScore(object = total, features = list(T_Score=T_gene$gene), name="T_Score")
p1=VlnPlot(total, features = "T_Score1", pt.size = 0,combine = FALSE, split.by="stim", cols=color_stim)
ggsave("/home/yushiya/cd34/fig/fig6-MF-T-score-all.pdf", plot=p1[[1]], width=7, height=3)


#calculate DEGs
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(GOSemSim)
library(DOSE)
celltype.markers <- FindAllMarkers(MF, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1)
write.csv(celltype.markers, "/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/MF_celltype_markers.csv")
Idents(total) <- total$stim
stim.markers <- FindAllMarkers(total, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1)
write.csv(stim.markers, "/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/MF_all_stim_markers.csv")
type <- c('HSC','MPP #1',"MPP #2",'MEP','EryP','MkP','LMPP #1','GMP','CDP',"pre-pDC","ETP","CLP","pro-B")
for (i in 1:length(type)) {
  sub1 <- subset(total, subset= celltype==type[i])
  Idents(sub1) <- sub1$stim
  cp_stim.markers <- FindAllMarkers(sub1, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1)
  write.csv(cp_stim.markers, paste("/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/cpstim_",type[i],".csv",sep=""))
}
stim.markers <- read.csv("/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/cpstim_MkP.csv",stringsAsFactors = F)
top_genes <- celltype.markers %>% group_by(cluster) %>% top_n(n = 40, wt = avg_log2FC)
top_genes <- subset(top_genes, cluster=="EryP")
sub1 <- subset(total, subset= celltype=="EryP")
DefaultAssay(sub1) <- "RNA"
sub1 <- ScaleData(sub1, verbose = FALSE)
DoHeatmap(sub1, features = top_genes$gene, group.by="stim",
          assay="RNA",
          group.colors=color_stim) + 
  scale_fill_gradientn(colors=c("blue","white","firebrick3"))
#GO analysis
BMID <- subset(stim.markers, subset= cluster=="BM")
BM_gene <- bitr(BMID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
mPBID <- subset(stim.markers, subset= cluster=="mPB")
mPB_gene <- bitr(mPBID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
PBID <- subset(stim.markers, subset= cluster=="PB")
PB_gene <- bitr(PBID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
MFID <- subset(stim.markers, subset= cluster=="MF")
MF_gene <- bitr(MFID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
cp = list(BM.gene=BM_gene$ENTREZID, mPB.gene=mPB_gene$ENTREZID, PB.gene=PB_gene$ENTREZID)
cp = list(BM.gene=BM_gene$ENTREZID, mPB.gene=mPB_gene$ENTREZID, PB.gene=PB_gene$ENTREZID, MF.gene=MF_gene$ENTREZID)
go.p <- compareCluster(cp,
                       fun = "enrichGO",
                       OrgDb = "org.Hs.eg.db",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.01
)
go.p1 <- simplify(go.p,cutoff=0.4,by="p.adjust",select_fun=min)  #去除冗余
ggplot(go.p1, aes(Cluster, Description), showCategory=6) +
  geom_point(aes(color=p.adjust, size=GeneRatio))+
  theme_classic()+
  theme(axis.text.x =element_text(angle = 45,vjust = 0.5,hjust = 0.5),
        text = element_text(size = 15),
        panel.background = element_blank(),
        panel.grid = element_line(colour = "grey90"))
#GSEA analysis
library(ggupset)
library(forcats)
library(ggstance)
stim.markers <- read.csv("/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/cpstim_MEP.csv",stringsAsFactors = F)
MF_genelist <- subset(stim.markers, subset= cluster=="MF")
gene <- MF_genelist$gene
gene <- bitr(gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
gene <- dplyr::distinct(gene,SYMBOL,.keep_all = T)
colnames(MF_genelist)[8] <- "SYMBOL"
MF_genelist <- MF_genelist %>% inner_join(gene, by="SYMBOL")
geneList= MF_genelist$avg_log2FC 
names(geneList)= MF_genelist$ENTREZID
geneList=sort(geneList,decreasing = T)
#gene_high <- geneList[abs(geneList)>0.5]
#gmt <- "https://wikipathways-data.wmcloud.org/current/gmt/wikipathways-20230410-gmt-Homo_sapiens.gmt"
#wp <- read.gmt.wp(gmt)
gse <- gseWP(geneList, organism = "Homo sapiens")
p1=upsetplot(gse) 
ggsave("/home/yushiya/cd34/fig/fig6-MF-gse-MEP-upset.pdf", plot=p1, width=9, height=4.5)
#p1=ridgeplot(gse)
#ggsave("/home/yushiya/cd34/fig/fig6-MF-gse-HSC-ridge.pdf", plot=p1, width=8, height=3.5)
p1=gseaplot2(gse, geneSetID = c(2,3,4,6), pvalue_table = TRUE, color=c("blue"))
ggsave("/home/yushiya/cd34/fig/fig6-MF-gse-MEP.pdf", plot=p1, width=11, height=6)
p1=ggplot(gse, aes(NES, fct_reorder(Description, NES),fill=qvalues), showCategory=10) + 
  geom_barh(stat='identity') + 
  scale_fill_continuous(low='#DB2D2D', high='#F1EC17') + 
  theme_minimal() + ylab(NULL)
ggsave("/home/yushiya/cd34/fig/fig6-MF-gse-MEP-bar.pdf", plot=p1, width=6.5, height=5)
saveRDS(gse, "/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/MF_gse.rds")

#gsea in all celltypes
gse@result$cluster="HSC"
gsea_all <- gse@result
type <- c('HSC','MPP #1',"MPP #2",'MEP','EryP','MkP','LMPP #1','GMP','CDP',"pre-pDC","ETP","CLP","pro-B")
for (i in 1:length(type)) {
  stim.markers <- read.csv(paste("/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/cpstim_",type[i],".csv",sep=""),stringsAsFactors = F)
  MF_genelist <- subset(stim.markers, subset= cluster=="MF")
  gene <- MF_genelist$gene
  gene <- bitr(gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
  gene <- dplyr::distinct(gene,SYMBOL,.keep_all = T)
  colnames(MF_genelist)[8] <- "SYMBOL"
  MF_genelist <- MF_genelist %>% inner_join(gene, by="SYMBOL")
  geneList= MF_genelist$avg_log2FC 
  names(geneList)= MF_genelist$ENTREZID
  geneList=sort(geneList,decreasing = T)
  gse <- gseWP(geneList, organism = "Homo sapiens")
  gse@result$cluster=type[i]
  gsea_all <- rbind(gsea_all,gse@result)
}
write.csv(gsea_all, "/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/gsea_all.csv")
gsea_all <- read.csv("/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/gsea_all.csv")
top <- gsea_all %>% group_by(cluster) %>% top_n(n = -5, wt = pvalue)
top_gsea <- top$Description %>% unique()
mtx_gsea <- data.frame(top_gsea)
colnames(mtx_gsea)[1] <- "Description"
for (i in 1:length(type)) {
  temp <- subset(gsea_all, subset= cluster==type[i])
  temp <- temp[,c("Description","enrichmentScore")]
  mtx_gsea <- mtx_gsea %>% left_join(temp, by="Description")
  mtx_gsea$enrichmentScore[is.na(mtx_gsea$enrichmentScore)] <- 0
  colnames(mtx_gsea)[i+1] <- type[i]
}
rownames(mtx_gsea) <- top_gsea
mtx_gsea <- mtx_gsea[,-1]
p1=pheatmap(mtx_gsea, cluster_cols = F, #cluster_rows = F,  
         clustering_method = "average",
         #scale = "row",
         #border=F,
         #col=colorRampPalette(c("navy","firebrick3"))(50),
         fontsize_row = 12, fontsize_col= 12, angle_col = 45)
ggsave("/home/yushiya/cd34/fig/fig6-MF-gse-all-heatmap.pdf", plot=p1, width=10, height=4.5)


#DEGs among MF and BM, mPB, PB
Idents(total) <- total$stim
stim.markers <- FindAllMarkers(total, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1)
write.csv(stim.markers, "/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/MF_all_stim_markers.csv")
stim.markers <- read.csv("/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/MF_all_stim_markers.csv", stringsAsFactors = F)
HSC_markers <- read.csv("/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/cpstim_HSC.csv", stringsAsFactors = F)
MPP1_markers <- read.csv("/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/cpstim_MPP #1.csv", stringsAsFactors = F)
MPP2_markers <- read.csv("/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/cpstim_MPP #2.csv", stringsAsFactors = F)
MEP_markers <- read.csv("/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/cpstim_MEP.csv", stringsAsFactors = F)
EryP_markers <- read.csv("/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/cpstim_EryP.csv", stringsAsFactors = F)
MkP_markers <- read.csv("/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/cpstim_MkP.csv", stringsAsFactors = F)
stim.markers <- rbind(HSC_markers,MPP1_markers,MPP2_markers,MEP_markers,EryP_markers,MkP_markers)
#find chemokines in DEGs
chemokines <- read.table("/home/yushiya/cd34/cytokines.txt")
chemo_genes <- stim.markers[stim.markers$gene %in% chemokines$V1,]
chemo_BM <- subset(chemo_genes, subset= cluster=="BM")
chemo_BM <- chemo_BM[order(chemo_BM$p_val_adj),]
chemo_BM <- chemo_BM[1:3,]
chemo_MF <- subset(chemo_genes, subset= cluster=="MF")
chemo_MF <- chemo_MF[order(chemo_MF$p_val_adj),]
chemo_MF <- chemo_MF[1:21,]
chemo <- rbind(chemo_BM,chemo_MF)
gene <- chemo$gene %>% unique()
#heatmap of DEGs in HSPC and Ery/Mk lineage
total <- subset(total, subset= stim %in% c("mPB","MF"))
total$stim_celltype <- paste(total$stim,total$celltype,sep="_")
sub1 <- subset(total, subset= celltype=="HSC" | celltype=="MPP #1" |celltype=="MPP #2" |celltype=="MEP" |celltype=="EryP" |celltype=="MkP")
cellInfo <- data.frame(stim_celltype=sub1$stim_celltype)
mtx <- data.frame(sub1@assays[["RNA"]]@data[gene,]) 
colnames(mtx) <- rownames(cellInfo)
top_exp <- sapply(split(rownames(cellInfo), cellInfo$stim_celltype),
                    function(cells) rowMeans(mtx[gene,cells]))
lev<-c("BM_HSC","BM_MPP #1","BM_MPP #2","BM_MEP","BM_EryP","BM_MkP",
       "mPB_HSC","mPB_MPP #1","mPB_MPP #2","mPB_MEP","mPB_EryP","mPB_MkP",
       "PB_HSC","PB_MPP #1","PB_MPP #2","PB_MEP","PB_EryP","PB_MkP",
       "MF_HSC","MF_MPP #1","MF_MPP #2","MF_MEP","MF_EryP","MF_MkP")
lev<-c("BM_HSC","mPB_HSC","PB_HSC","MF_HSC","BM_MPP #1","mPB_MPP #1","PB_MPP #1","MF_MPP #1",
       "BM_MPP #2","mPB_MPP #2","PB_MPP #2","MF_MPP #2","BM_MEP","mPB_MEP","PB_MEP","MF_MEP",
       "BM_EryP","mPB_EryP","PB_EryP","MF_EryP","BM_MkP","mPB_MkP","PB_MkP","MF_MkP")
lev<-c("BM_HSC","mPB_HSC","PB_HSC","MF_HSC","BM_MPP #1","mPB_MPP #1","PB_MPP #1","MF_MPP #1",
       "BM_MPP #2","mPB_MPP #2","PB_MPP #2","MF_MPP #2","BM_MEP","mPB_MEP","PB_MEP","MF_MEP",
       "BM_EryP","mPB_EryP","PB_EryP","MF_EryP","BM_MkP","mPB_MkP","PB_MkP","MF_MkP")
lev<-c("mPB_HSC","MF_HSC","mPB_MPP #1","MF_MPP #1",
       "mPB_MPP #2","MF_MPP #2","mPB_MEP","MF_MEP",
       "mPB_EryP","MF_EryP","mPB_MkP","MF_MkP")
top_exp <- top_exp[,lev] 
annotation_col <- data.frame(Stim = factor(rep(c("BM", "mPB","PB","MF"), c(6))))
annotation_col <- data.frame(Stim = factor(rep(c("mPB","MF"), c(6))))
annotation_col <- data.frame(Stim = factor(rep(c("BM", "PB","MF"), c(6))))
rownames(annotation_col) <- lev
annotation_colors =list(Stim=c(BM="#DB5C25",mPB="#F3B747",PB="#649541",MF="#9966CC"))
annotation_colors =list(Stim=c(mPB="#F3B747",MF="#9966CC"))
p1=pheatmap(top_exp, cluster_cols = F, #cluster_rows = F,  
            clustering_method = "average",
            annotation_col = annotation_col,
            annotation_colors = annotation_colors,
            scale = "row",
            #border=F,
            col=colorRampPalette(c("navy","white" ,"firebrick3"))(50),
            fontsize_row = 12, fontsize_col= 12, angle_col = 45)
ggsave("/home/yushiya/cd34/fig/fig6-MF-chemo-all-1.pdf", plot=p1, width=8, height=4)


#cellchat with BM stromal cells
library(CellChat)
library(svglite)
sub_sc <- readRDS(file = "/home/yushiya/cd34/data/stromal_sub.rds")
sub_tec <- readRDS(file = "/home/yushiya/cd34/data/TEC_sub.rds")
MF <- subset(total, subset= stim=="MF")
Idents(MF) <- MF$predicted.celltype
ETP <- subset(MF, subset= predicted.celltype=="ETP")
sum <- merge(MF,sub_sc)
sum <- merge(ETP,c(sub_sc,sub_tec))
cellchat <- createCellChat(sum)
cellchat@idents <- factor(cellchat@idents, levels=c("BM_SC","BM_Endo","HSC","MPP #1","MPP #2","LMPP #1","MEP","EryP","MkP","GMP","CDP","pre-pDC","ETP","CLP","pro-B"))
cellchat@idents <- factor(cellchat@idents, levels=c("ETP","BM_SC","BM_Endo","Thymic_Epi","Thymic_Mesen","Thymic_Endo"))
CellChatDB <- CellChatDB.human 
CellChatDB.use <- CellChatDB 
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
#识别细胞组中过度表达的配体或受体
cellchat <- identifyOverExpressedGenes(cellchat)
#识别过度表达的配体受体相互作用
cellchat <- identifyOverExpressedInteractions(cellchat)
#将基因表达数据投射到PPI网络上
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# 如果在某些细胞群中只有少数细胞，则过滤掉细胞间的通信
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
dev.off()
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat2[c(1:2),] <- mat[c(1:2),]
netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = "Interaction weights/strength")
ggsave("/home/yushiya/cd34/fig/fig6-cellchat-all-circle.pdf", plot=p1, width=6, height=6)
# 显示从某些细胞组到其他细胞组的所有显著的相互作用（L-R 对）
p1=netVisual_bubble(cellchat, sources.use = c(1:2), targets.use = c(3:14), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = c(2:6), targets.use = c(1), remove.isolate = FALSE)
ggsave("/home/yushiya/cd34/fig/fig6-cellchat-LRpairs.pdf", plot=p1, width=8, height=7)
pathways.show <- c("CXCL") 
netVisual_aggregate(cellchat, signaling = pathways.show, sources.use = c(1,2), targets.use = c(3:14))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds", sources.use = c(1:2), targets.use = c(3:14))
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # “netP”:推断出的信号通路的细胞间通信网络
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 18, height = 8, font.size = 14)
netAnalysis_contribution(cellchat, signaling = pathways.show)
ggsave("/home/yushiya/cd34/fig/fig6-cellchat-CXCL-circle.pdf", plot=p1, width=6, height=6)
ggsave("/home/yushiya/cd34/fig/fig6-cellchat-CXCL-heatmap.pdf", plot=p1, width=5, height=4)
ggsave("/home/yushiya/cd34/fig/fig6-cellchat-CXCL-sig.pdf", plot=p1, width=10, height=5)
ggsave("/home/yushiya/cd34/fig/fig6-cellchat-CXCL-ctb.pdf", plot=p1, width=3.5, height=4)
cellchat <- readRDS(file = "/home/yushiya/cd34/fig/cellchat_MF.rds")
saveRDS(cellchat, file = "/home/yushiya/cd34/fig/cellchat_MF.rds")





