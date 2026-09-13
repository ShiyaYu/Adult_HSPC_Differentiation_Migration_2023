library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)

#add stromal cell and endothelial cell
sc.data <- Read10X(data.dir="/home/yushiya/cd34/data/stromal_GSE190965")
sc <- CreateSeuratObject(counts = sc.data, project = "SC")
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
sc <- subset(sc, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000)
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sc)
sc <- ScaleData(sc, features = all.genes)
sc <- RunPCA(sc, features = VariableFeatures(object = sc))
VizDimLoadings(sc, dims = 1:5, reduction = "pca")
ElbowPlot(sc)
DimPlot(sc, reduction = "pca")
sc <- RunUMAP(sc, reduction = "pca", dims = 1:15)
sc <- FindNeighbors(sc, reduction = "pca", dims = 1:15)
sc <- FindClusters(sc, resolution = 0.8)
color <- c(brewer.pal(12,"Paired"),brewer.pal(12,"Set3"))
DimPlot(sc, reduction = "umap", label = T, repel=T, cols=color)
ggsave("/home/yushiya/cd34/fig/fig-sc-umap.pdf", plot=p1, width=6, height=4)
p1=FeaturePlot(sc, features = c("CXCL12","VCAN","LEPR","PECAM1","ICAM2"), max.cutoff = 3, cols = c("#89B7D3", "red"),ncol=3)
ggsave("/home/yushiya/cd34/fig/fig-sc-FeaP.pdf", plot=p1, width=8, height=4)
#saveRDS(sc, file = "/home/yushiya/cd34/data/stromal.rds")
#select cluster 0,3,8 as stromal cell, cluster 6 as endothelial cell
sub_sc <- subset(sc, subset= seurat_clusters=="0" | seurat_clusters=="3" | seurat_clusters=="8" | seurat_clusters=="6")
sc_gene <- c("NNMT","IFITM3","DCN","CXCL12","CFD","APOE","PTGDS","LEPR","TF","CHL1","VCAN","IGFBP5","PECAM1","ICAM2","LST1")
sub_sc$seurat_clusters <- factor(sub_sc$seurat_clusters, levels=c("0","3","8","6"))
DoHeatmap(sub_sc, features = sc_gene, group.by="seurat_clusters",
          group.colors=c("#A6CEE3","#33A02C","#CAB2D6","#FDBF6F")) + 
  scale_fill_gradientn(colors=c("blue","white","firebrick3"))
ggsave("/home/yushiya/cd34/fig/fig-sc-doHeat.pdf", plot=p1, width=12, height=4)
sub_sc <- RenameIdents(sub_sc,'0'="stromal cell",'3'="stromal cell",'8'="stromal cell",'6'="endothelial cell")
sub_sc <- RenameIdents(sub_sc,"stromal cell"="BM_SC","endothelial cell"="BM_Endo")
sub_sc$celltype <- Idents(sub_sc)
p1=DimPlot(sub_sc, reduction = "umap", label = F, cols=c("#1F78B4","#FDBF6F"))
ggsave("/home/yushiya/cd34/fig/fig-sc-sub-umap.pdf", plot=p1, width=6, height=4)
p1=FeaturePlot(sub_sc, features = c("CXCL12","VCAN","LEPR","PECAM1","ICAM2"), max.cutoff = 3, cols = c("#89B7D3", "red"),ncol=3)
ggsave("/home/yushiya/cd34/fig/fig-sc_sub-FeaP.pdf", plot=p1, width=8, height=4)
p1=VlnPlot(sub_sc, features = c("APP"), cols=c("#1F78B4","#FDBF6F"), pt.size = 0)
ggsave("/home/yushiya//fig/fig3-cellchat-APP-1.pdf", plot=p1, width=3, height=2.3)
#saveRDS(sub_sc, file = "/home/yushiya/cd34/data/stromal_sub.rds")
sc <- readRDS("/home/yushiya/cd34/data/stromal.rds")
sub_sc <- readRDS(file = "/home/yushiya/data/cd34/data/stromal_sub.rds")

#cellchat
#BM stromal cells
library(RColorBrewer)
allBM <- readRDS("/home/yushiya/data/cd34/data/2024_BM_CODEX/GSE253355_Normal_Bone_Marrow_Atlas_Seurat_SB_v2.rds")
DimPlot(allBM)
sub_sc <- subset(allBM, subset= cluster_anno_l2 %in% c("Fibro-MSC","APOD+ MSC","Osteoblast","Osteo-MSC","THY1+ MSC","Adipo-MSC","SEC","AEC"))
Idents(sub_sc) <- factor(Idents(sub_sc), levels = c("Fibro-MSC","APOD+ MSC","Osteoblast","Osteo-MSC","THY1+ MSC","Adipo-MSC","SEC","AEC"))
p1=VlnPlot(sub_sc, features = "TGM2", pt.size = 0, cols = brewer.pal(8, "Paired"))
ggsave("/home/yushiya/data/cd34/data/fig/fig3-cellchat-TGM2.pdf", plot=p1, width=5, height=3)
p1=VlnPlot(sum, features = "ADGRG1", pt.size = 0, cols = color_stim, split.by = "stim")
ggsave("/home/yushiya/data/cd34/data/fig/fig3-cellchat-ADGRG1.pdf", plot=p1, width=7, height=2.5)
#Idents()
sub1 <- subset(immune.combined, subset= celltype %in% c("HSC","MPP#1","MPP#2","LMPP#1","LMPP#2","BMEP","EryP","MkP","Ma/Eo/BaP","GMP","CDP","pre-pDC","CLP","pro-B","pre-B"))
sub1 <- subset(sub1, subset= stim=="PB" | stim=="BM" | stim=="mPB" | stim=="SP")
#devtools::install_github("sqjin/CellChat")
library(CellChat)
library(svglite)
#Idents(sub1) <- sub1$stim_celltype
sum <- merge(sub_sc,sub1)
#sum <- merge(sum,sub_tec)
sum <- readRDS("/home/yushiya/data/cd34/data/fig/immune_stromal-1.rds")
Idents(sum) <- sum$celltype
sum <- subset(sum, subset= (tissue=="PB") | (celltype %in% c("Fibro-MSC","APOD+ MSC","Osteoblast","Osteo-MSC","THY1+ MSC","Adipo-MSC","SEC","AEC")) )
cellchat <- createCellChat(sum)
cellchat@idents <- factor(cellchat@idents, levels=c("Fibro-MSC","APOD+ MSC","Osteoblast","Osteo-MSC","THY1+ MSC","Adipo-MSC","SEC","AEC","HSC","MPP#1","MPP#2","LMPP#1","LMPP#2","BMEP","EryP","MkP","Ma/Eo/BaP","GMP","CDP","pre-pDC","CLP","pro-B","pre-B"))
cellchat@idents <- factor(cellchat@idents, levels=c("BM_SC","BM_Endo","HSC","MPP#1","MPP#2","LMPP#1","LMPP#2","BMEP","EryP","MkP","Ma/Eo/BaP","GMP","CDP","pre-pDC","CLP","pro-B","pre-B"))
cellchat@idents <- factor(cellchat@idents, levels=c("Thymic_Epi","Thymic_Mesen","Thymic_Endo","HSC","MPP #1","MPP #2","LMPP #1","LMPP #2","ETP","Thy #1","Thy #2","Thy #3"))
cellchat@idents <- factor(cellchat@idents, levels=c("BM_SC","BM_Endo","Thymic_Epi","Thymic_Endo","BM_ETP","mPB_ETP","PB_ETP","thymus_ETP"))
cellchat@idents <- factor(cellchat@idents, levels=c("BM_SC","BM_Endo","BM_HSC","BM_MPP #1","BM_MPP #2","BM_LMPP #1","BM_LMPP #2",
                                                    "mPB_HSC","mPB_MPP #1","mPB_MPP #2","mPB_LMPP #1","mPB_LMPP #2",
                                                    "PB_HSC","PB_MPP #1","PB_MPP #2","PB_LMPP #1","PB_LMPP #2"))
CellChatDB <- CellChatDB.human 
CellChatDB.use <- CellChatDB 
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
#future::plan("multicore", workers = 10) 
#识别细胞组中过度表达的配体或受体
cellchat <- identifyOverExpressedGenes(cellchat)
#识别过度表达的配体受体相互作�?
cellchat <- identifyOverExpressedInteractions(cellchat)
#计算通信概率并推断细胞通信网络
#cellchat <- projectData(cellchat, PPI.human)
#cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- computeCommunProb(cellchat, type = "triMean")
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
mat2[c(1:8),c(9:23)] <- mat[c(1:8),c(9:23)]
mat2[c(1:2),c(3:17)] <- mat[c(1:2),c(3:17)]
mat2[c(1:2),c(6:9)] <- mat[c(1:2),c(6:9)]
mat_new <- mat[c(1:8),c(9:13,18:20)]
mat_new <- mat[c(1:8),c(9:13,21:23)]
sum <- as.matrix(rowSums(mat_new))
colnames(sum) <- "HSPC"
m <- diag(c(0, 0), nrow = 9, ncol = 9)
rownames(m) <- colnames(m) <- c(
  "Fibro-MSC", "APOD+ MSC", "Osteoblast",
  "Osteo-MSC", "THY1+ MSC", "Adipo-MSC", "SEC", "AEC", "HSPC"
)
m[,9] <- c(sum,0)
groupSize <- groupSize[1:9]
netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = "Interaction weights/strength")
netVisual_circle(m, vertex.weight = groupSize, weight.scale = T, title.name = "Sum of Interaction weights")
ggsave("/home/yushiya/fig/fig3-cellchat-SC-BM-circle.pdf", plot=p1, width=6, height=6)
# 显示从某些细胞组到其他细胞组的所有显著的相互作用（L-R 对）
p1=netVisual_bubble(cellchat_mPB, sources.use = c(1:8), targets.use = c(9:23), remove.isolate = FALSE)
p1=netVisual_bubble(cellchat, sources.use = c(1:2), targets.use = c(3:17), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = c(1:8), targets.use = c(9:23), remove.isolate = FALSE)
p1=netVisual_bubble(cellchat_BM, sources.use = c(1:8), targets.use = c(9:23),signaling = c("APP","CD34","ADGRE","NOTCH","CSF"), remove.isolate = TRUE)
netVisual_bubble(cellchat_BM, sources.use = c(1:8), targets.use = c(9:23),signaling = c("CD34","ADGRE"), remove.isolate = TRUE)
LRparirs <- data.frame(pairs_3)
LRparirs <- data.frame(c("CD55_ADGRE5","ICAM1_ITGAL","ICAM1_SPN","PECAM1_CD38","SELE_GLG1","TGM2_ADGRG1"))
LRparirs <- data.frame(c("CXCL12_CXCR4","MDK_NCL","APP_CD74","APP_SORL1"))
colnames(LRparirs) <- "interaction_name"
p1=netVisual_bubble(cellchat_mPB, sources.use = c(1:8), targets.use = c(9:23),pairLR.use = LRparirs, remove.isolate = TRUE)
netVisual_bubble(cellchat_PB, sources.use = c(1:8), targets.use = c(9:23),pairLR.use = LRparirs, remove.isolate = TRUE, max.quantile=0.9)
p1=netVisual_bubble(cellchat_mPB, sources.use = c(7,8), targets.use = c(9:23), pairLR.use = LRparirs,remove.isolate = FALSE)
p1=netVisual_bubble(cellchat_BM, sources.use = c(5), targets.use = c(9:23), pairLR.use = LRparirs,remove.isolate = FALSE)
p1=netVisual_bubble(cellchat, sources.use = c(1:2), targets.use = c(3:17),signaling = c("APP","CXCL","MK"), remove.isolate = TRUE)
p1=netVisual_bubble(cellchat, sources.use = c(1), targets.use = c(3:17),signaling = c("CXCL","FN1","APP"), remove.isolate = TRUE)
ggsave("/home/yushiya/data/cd34/data/fig/fig3-cellchat-SC-BM-LRpairs.pdf", plot=p1, width=8, height=12)
ggsave("/home/yushiya/data/cd34/data/fig/fig3-cellchat-LRpairs-all-nor-mPB-selected.pdf", plot=p1, width=16, height=5)
ggsave("/home/yushiya/data/cd34/data/fig/fig3-cellchat-LRpairs-all-nor-mPB-selected-1.pdf", plot=p1, width=7, height=2.5)
ggsave("/home/yushiya/data/cd34/data/fig/fig3-cellchat-SC-BM-LRpairs-selected.pdf", plot=p1, width=5, height=2.5)
ggsave("/home/yushiya/cd34/fig/fig3-cellchat-LRpairs-mPB.pdf", plot=p1, width=8, height=10)
p1=plotGeneExpression(cellchat, signaling = c("CXCL","APP","SELPLG","ANGPTL"),color.use = c("#550000","#AA3939","#CC6F66","#E79492","#FBD2CE","#1F78B4","#FDBF6F"))
p1=plotGeneExpression(cellchat, signaling = c("PTN","CD99","MIF","CXCL","NOTCH","CCL"),color.use = c("#1F78B4","#FDBF6F","#66C2A5","#FC8D62","#8DA0CB","#DB5C25","#F3B747","#649541","#4C82C5"))
ggsave("/home/yushiya/cd34/fig/fig-cellchat-geneExp-3.pdf", plot=p1, width=5, height=6)
saveRDS(cellchat, file = "/home/yushiya/data/cd34/data/fig/cellchat_SC_all_nor_BM.rds")
#"CXCL","APP","SELPLG","ANGPTL","ADGRE5"
pathways.show <- c("COLLAGEN") 
mat <- cellchat_BM@netP[["prob"]][,,pathways.show]
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat2[c(1:8),c(9:23)] <- mat[c(1:8),c(9:23)]
netVisual_circle(mat2, weight.scale = T, edge.weight.max = max(mat), title.name = "Interaction weights/strength")
netVisual_aggregate(cellchat, signaling = pathways.show, sources.use = c(1,8), targets.use = c(9:23))
netVisual_aggregate(cellchat_BM, signaling = pathways.show, sources.use = c(1:8), targets.use = c(9:23))
netVisual_chord_cell(cellchat_BM, signaling = pathways.show, title.name = paste0(pathways.show, " signaling network"))
netVisual_heatmap(cellchat_BM, signaling = pathways.show, color.heatmap = "Reds", sources.use = c(1,8), targets.use = c(9:23))
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # “netP�?:推断出的信号通路的细胞间通信网络
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, color.heatmap = "BuGn", width = 18, height = 8, font.size.title = 18, font.size = 14)
netAnalysis_contribution(cellchat_BM, signaling = pathways.show)
ggsave("/home/yushiya/data/cd34/data/fig/fig3-cellchat-SC-BM-CXCL-sig.pdf", plot=p1, width=6, height=6)
ggsave("/home/yushiya/data/cd34/data/fig/fig3-cellchat-COLLAGEN-BM-heatmap.pdf", plot=p1, width=9, height=5)
ggsave("/home/yushiya/data/cd34/data/fig/fig3-cellchat-SELPLG-mPB-sig.pdf", plot=p1, width=10, height=5)
ggsave("/home/yushiya/data/cd34/data/fig/fig3-cellchat-CXCL-ctb.pdf", plot=p1, width=10, height=5)
#cellchat <- readRDS(file = "/home/yushiya/data/cd34/data/fig/cellchat_all.rds")
cellchat_BM <- readRDS(file = "/home/yushiya/data/cd34/data/fig/cellchat_all_nor_BM.rds")
cellchat_mPB <- readRDS(file = "/home/yushiya/data/cd34/data/fig/cellchat_all_nor_mPB.rds")
cellchat_PB <- readRDS(file = "/home/yushiya/data/cd34/data/fig/cellchat_all_nor_PB.rds")
cellchat <- readRDS(file = "/home/yushiya/data/cd34/data/fig/cellchat_SC_all_nor_BM.rds")
#saveRDS(cellchat, file = "/home/yushiya/data/cd34/data/fig/cellchat_all-1.rds")
saveRDS(cellchat, file = "/home/yushiya/data/cd34/data/fig/cellchat_all_nor.rds")

#filter CCC signals
cellchat_BM <- readRDS(file = "/home/yushiya/data/cd34/data/fig/cellchat_all_nor_BM.rds")
p1=netVisual_bubble(cellchat_BM, sources.use = c(1:8), targets.use = c(9:23), remove.isolate = FALSE)
ccc_mtx <- data.frame(p1[["data"]][["source.target"]])
colnames(ccc_mtx)[1] <- "source_target"
ccc_mtx$interaction_name <- p1[["data"]][["interaction_name_2"]]
ccc_mtx$prob <- p1[["data"]][["prob"]]
ccc_mtx$pval <- p1[["data"]][["pval"]]
#ccc_mtx <- ccc_mtx[(log(ccc_mtx$pval)<0.01),]
var_df <- ccc_mtx %>%
  group_by(interaction_name) %>%
  summarise(prob_var = var(prob), .groups = "drop")
write.csv(var_df,file='/home/yushiya/data/cd34/data/fig/cellchat_BM_var.csv')
var_cellphoneDB <- apply(re_draw , 1, function(x) var(x, na.rm = TRUE))
var_cellphoneDB <- data.frame(var_cellphoneDB)
var_cellphoneDB$interaction <- rownames(re_draw)
write.csv(var_cellphoneDB,file='/home/yushiya/data/cd34/data/fig/cellphoneDB_BM_var.csv')

#CCC interaction statistics
library(reshape2)
library("ggalluvial")
library(RColorBrewer)
stas <- data.frame(c(3,4,6,3),c(9,1,8,7))
colnames(stas) <- c("non-migration","migration")
rownames(stas) <- c("Fibro","Osteo","Mesenchymal","EC")
stas <- t(stas)
Type=colnames(stas)
stas=melt(stas, id='Type')
names(stas)[1]='Type'
names(stas)[2]='Location'
stas$Location <- factor(stas$Location,levels = c("Fibro","Osteo","Mesenchymal","EC"))
ggplot(stas,
       aes(x=Location, y=value, fill=Type)) +
  geom_bar(stat='identity', width=0.45) +
  #geom_alluvium() +
  #geom_stratum(width=0.45, size=0.1) +
  scale_fill_manual(values = brewer.pal(9,"Set1")[2:3]) +
  labs(x='Location', y='CCC numbers')+
  scale_y_continuous(expand=c(0, 0))+
  theme(axis.text.x=element_text(hjust=1, angle=45, vjust=1),
        text = element_text(size = 18),
        panel.background = element_blank(),
        #panel.grid = element_line(colour = "grey90"),
        axis.line = element_line(colour = "black",size = rel(1)))# + theme_bw()


#spatial cellchat results
library(CellChat)
spatial_cellchat <- readRDS("/home/zhangyutao/jinlab/3.20240913_AgingNeutro/5.ST/1.CCC/singleCell_cellChat/sp_compare_Cd34/cellChat_spatial_y_cd34.rds")
netVisual_bubble(spatial_cellchat, sources.use = c(3), remove.isolate = FALSE)
netVisual_bubble(spatial_cellchat, sources.use = c(3), targets.use = c(12), remove.isolate = FALSE)
pathways.show <- c("ICAM")
netVisual_aggregate(spatial_cellchat, signaling = pathways.show, sources.use = c(3), targets.use = c(1:13))
netVisual_heatmap(spatial_cellchat, signaling = pathways.show, color.heatmap = "Reds", sources.use = c(3), targets.use = c(12))
spatial_cellchat <- netAnalysis_computeCentrality(spatial_cellchat, slot.name = "netP") # “netP�?:推断出的信号通路的细胞间通信网络
netAnalysis_signalingRole_network(spatial_cellchat, signaling = pathways.show, width = 18, height = 8, font.size.title = 16, font.size = 14)



#Nichnet
library(nichenetr) # Please update to v2.0.4
library(SeuratObject)
library(tidyverse)
sum <- readRDS("/home/yushiya/data/cd34/data/fig/immune_stromal.rds")
Idents(sum) <- sum$celltype
#sum <- subset(sum, subset= (tissue=="BM") | (celltype %in% c("Fibro-MSC","APOD+ MSC","Osteoblast","Osteo-MSC","THY1+ MSC","Adipo-MSC","SEC","AEC")) )
Idents(sum) <- factor(Idents(sum), levels=c("Fibro-MSC", "APOD+ MSC", "Osteoblast", "Osteo-MSC", "THY1+ MSC", "Adipo-MSC", "SEC", "AEC", "HSC","MPP#1","MPP#2","LMPP#1","LMPP#2","BMEP","EryP","MkP","Ma/Eo/BaP","GMP","CDP","pre-pDC","CLP","pro-B","pre-B"))
organism <- "human"
if(organism == "human"){
  lr_network <- readRDS("/home/yushiya/data/reference/Nichnet_database/lr_network_human_21122021.rds")
  ligand_target_matrix <- readRDS("/home/yushiya/data/reference/Nichnet_database/ligand_target_matrix_nsga2r_final.rds")
  weighted_networks <- readRDS("/home/yushiya/data/reference/Nichnet_database/weighted_networks_nsga2r_final.rds")
} else if(organism == "mouse"){
  lr_network <- readRDS("/home/yushiya/data/reference/Nichnet_database/lr_network_mouse_21122021.rds")
  ligand_target_matrix <- readRDS("/home/yushiya/data/reference/Nichnet_database/ligand_target_matrix_nsga2r_final_mouse.rds")
  weighted_networks <- readRDS("/home/yushiya/data/reference/Nichnet_database/weighted_networks_nsga2r_final_mouse.rds")
}
lr_network <- lr_network %>% distinct(from, to)
head(lr_network)
#Define a set of potential ligands for reciever and sender celltypes
receiver = "HSC"
expressed_genes_receiver <- get_expressed_genes(receiver, sum, pct = 0.05)
all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()
sender_celltypes <- c("Fibro-MSC", "APOD+ MSC", "Osteoblast", "Osteo-MSC", "THY1+ MSC", "Adipo-MSC", "SEC", "AEC")
# Use lapply to get the expressed genes of every sender cell type separately here
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, sum, 0.05)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()
potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 
#Define the gene set of interest
condition_oi <-  "PB"
condition_reference <- "BM"
seurat_obj_receiver <- subset(sum, idents = receiver)
DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
                                  ident.1 = condition_oi, ident.2 = condition_reference,
                                  group.by = "tissue",
                                  min.pct = 0.05) %>% rownames_to_column("gene")
geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
#Define the background genes
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
#Perform NicheNet ligand activity analysis
ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)
ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
#select top 30 genes
best_upstream_ligands <- ligand_activities %>% top_n(5, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)
best_upstream_ligands <- c("IL17F","NCAM1","OCLN","GSTP1","CD96")
best_upstream_ligands <- c("CXCL12","MDK","SELE","PODXL2","PECAM1","APP","ICAM1","CD55")
best_upstream_ligands <- c("CXCL12","MDK","APP")
vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)
make_heatmap_ggplot(vis_ligand_aupr,
                    "Prioritized ligands", "Ligand activity", 
                    legend_title = "AUPR", color = "darkorange") + 
  theme(axis.text.x.top = element_blank())
#Infer target genes and receptors of top-ranked ligands
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()
active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33) 
order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))
vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])
make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                    color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 
vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 
make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                    y_name = "Ligands", x_name = "Receptors",  
                    color = "mediumvioletred", legend_title = "Prior interaction potential")+ 
  theme(axis.text.y.top = element_blank())
#Infer signaling paths beween ligand(s) and target(s) of interest
ligands_oi <- "ICAM1" # this can be a list of multiple ligands if required
targets_oi <- c("SPN")
active_signaling_network <- get_ligand_signaling_path(ligands_all = ligands_oi,
                                                      targets_all = targets_oi,
                                                      weighted_networks = weighted_networks,
                                                      ligand_tf_matrix = ligand_target_matrix,
                                                      top_n_regulators = 4,
                                                      minmax_scaling = TRUE) 
graph_min_max <- diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network,
                                                   ligands_all = ligands_oi, targets_all = targets_oi,
                                                   sig_color = "indianred", gr_color = "steelblue")
save(graph_min_max, file = "/home/yushiya/data/cd34/data/fig/nichnet_data/graph_ICAM1.RData")
graph_svg <- DiagrammeRsvg::export_svg(DiagrammeR::render_graph(graph_min_max, layout = "tree", output = "graph"))
p1=cowplot::ggdraw() + cowplot::draw_image(charToRaw(graph_svg))
ggsave("/home/yushiya/data/cd34/data/fig/nichnet_data/fig_ICAM1.pdf", plot=p1, width=5, height=5)
data_source_network <- infer_supporting_datasources(signaling_graph_list = active_signaling_network,
                                                    lr_network = lr_network, sig_network = sig_network, gr_network = gr_network)
CCC_pairs_BM <- read.csv("/home/yushiya/data/cd34/data/fig/CCC_pairs_BM.txt", header = F)
CCC_pairs_PB <- read.csv("/home/yushiya/data/cd34/data/fig/CCC_pairs_PB.txt", header = F)
pairs <- union(CCC_pairs_BM$V1,CCC_pairs_PB$V1)
best_upstream_ligands <- str_split(pairs, "_", simplify = TRUE)[, 1] %>% unique()
#Visualizing results in all HSPC celltypes
type <- c("HSC","MPP#1","MPP#2","LMPP#1","LMPP#2","BMEP","EryP","MkP","Ma/Eo/BaP","GMP","CDP","pre-pDC","CLP","pro-B","pre-B")
for (i in 2:length(type)) {
  receiver = type[i]
  expressed_genes_receiver <- get_expressed_genes(receiver, sum, pct = 0.05)
  all_receptors <- unique(lr_network$to)  
  expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
  potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()
  sender_celltypes <- c("Fibro-MSC", "APOD+ MSC", "Osteoblast", "Osteo-MSC", "THY1+ MSC", "Adipo-MSC", "SEC", "AEC")
  list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, sum, 0.05)
  expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()
  potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 
  condition_oi <-  "PB"
  condition_reference <- "BM"
  seurat_obj_receiver <- subset(sum, idents = receiver)
  DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
                                    ident.1 = condition_oi, ident.2 = condition_reference,
                                    group.by = "tissue",
                                    min.pct = 0.05) %>% rownames_to_column("gene")
  geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
  geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
  background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
  ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                                 background_expressed_genes = background_expressed_genes,
                                                 ligand_target_matrix = ligand_target_matrix,
                                                 potential_ligands = potential_ligands)
  ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
  #best_upstream_ligands <- c("CXCL12","MDK","SELE","PODXL2","PECAM1","APP","ICAM1","CD55")
  vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
    column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)
  make_heatmap_ggplot(vis_ligand_aupr,
                      "Prioritized ligands", "Ligand activity", 
                      legend_title = "AUPR", color = "darkorange") + 
    theme(axis.text.x.top = element_blank())
  #Infer target genes and receptors of top-ranked ligands
  best_upstream_ligands1 <- intersect(best_upstream_ligands,colnames(ligand_target_matrix))
  active_ligand_target_links_df <- best_upstream_ligands1 %>%
    lapply(get_weighted_ligand_target_links,
           geneset = geneset_oi,
           ligand_target_matrix = ligand_target_matrix,
           n = 100) %>%
    bind_rows() %>% drop_na()
  active_ligand_target_links <- prepare_ligand_target_visualization(
    ligand_target_df = active_ligand_target_links_df,
    ligand_target_matrix = ligand_target_matrix,
    cutoff = 0.33) 
  order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
  order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))
  vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])
  vis_ligand_target <- data.frame(vis_ligand_target)
  if (type[i]=="Ma/Eo/BaP") {
    write.table(vis_ligand_target, paste("/home/yushiya/data/cd34/data/fig/nichnet_data/ligand_target_Ma.txt",sep=''),col.names = T, row.names = T)
  }else {
    write.table(vis_ligand_target, paste("/home/yushiya/data/cd34/data/fig/nichnet_data/ligand_target_",type[i],".txt",sep=''),col.names = T, row.names = T)
  }
  vis_ligand_target$ligand <- rownames(vis_ligand_target)
  vis_ligand_target <- melt(vis_ligand_target, id="ligand")
  rownames(vis_ligand_target) <- paste(vis_ligand_target$ligand, vis_ligand_target$variable, sep="_")
  vis_ligand_target$pairs <- paste(vis_ligand_target$ligand, vis_ligand_target$variable, sep="_")
  vis_ligand_target <- subset(vis_ligand_target, subset= (value > 0))
  vis_ligand_target <- vis_ligand_target[,c(4,3)]
  ligand_target_add <- full_join(ligand_target_add, vis_ligand_target, by="pairs")
  #make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
  #                    color = "purple", legend_title = "Regulatory potential") +
  #  scale_fill_gradient2(low = "whitesmoke",  high = "purple")
  ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
    best_upstream_ligands, expressed_receptors,
    lr_network, weighted_networks$lr_sig) 
  vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
    ligand_receptor_links_df,
    best_upstream_ligands,
    order_hclust = "both")
  vis_ligand_receptor_network <- data.frame(vis_ligand_receptor_network)
  if (type[i]=="Ma/Eo/BaP") {
    write.table(vis_ligand_receptor_network, paste("/home/yushiya/data/cd34/data/fig/nichnet_data/ligand_receptor_network_Ma.txt",sep=''),col.names = T, row.names = T)
  }else {
    write.table(vis_ligand_receptor_network, paste("/home/yushiya/data/cd34/data/fig/nichnet_data/ligand_receptor_network_",type[i],".txt",sep=''),col.names = T, row.names = T)
  }
  vis_ligand_receptor_network$receptor <- rownames(vis_ligand_receptor_network)
  vis_ligand_receptor_network <- melt(vis_ligand_receptor_network, id="receptor")
  rownames(vis_ligand_receptor_network) <- paste(vis_ligand_receptor_network$variable, vis_ligand_receptor_network$receptor, sep="_")
  vis_ligand_receptor_network$pairs <- paste(vis_ligand_receptor_network$variable, vis_ligand_receptor_network$receptor, sep="_")
  vis_ligand_receptor_network <- subset(vis_ligand_receptor_network, subset= (value > 0))
  vis_ligand_receptor_network <- vis_ligand_receptor_network[,c(4,3)]
  ligand_receptor_network_add <- full_join(ligand_receptor_network_add, vis_ligand_receptor_network, by="pairs")
  #make_heatmap_ggplot(t(vis_ligand_receptor_network), 
  #                    y_name = "Ligands", x_name = "Receptors",  
  #                    color = "mediumvioletred", legend_title = "Prior interaction potential")+ 
  #  theme(axis.text.y.top = element_blank())
}
ligand_target_add[is.na(ligand_target_add)] <- 0
rownames(ligand_target_add) <- ligand_target_add$pairs
ligand_target_add <- ligand_target_add[,c(2:16)]
colnames(ligand_target_add) <- type
write.table(ligand_target_add, paste("/home/yushiya/data/cd34/data/fig/nichnet_data/ligand_target_add.txt",sep=''),col.names = T, row.names = T)
ligand_target_add <- read.table("/home/yushiya/data/cd34/data/fig/nichnet_data/ligand_target_add.txt")

ligand_receptor_network_add[is.na(ligand_receptor_network_add)] <- 0
rownames(ligand_receptor_network_add) <- ligand_receptor_network_add$pairs
ligand_receptor_network_add <- ligand_receptor_network_add[,c(2:16)]
colnames(ligand_receptor_network_add) <- type
write.table(ligand_receptor_network_add, paste("/home/yushiya/data/cd34/data/fig/nichnet_data/ligand_receptor_network_add.txt",sep=''),col.names = T, row.names = T)
ligand_receptor_network_add$pairs <- rownames(ligand_receptor_network_add)
ligand_receptor_network_add=melt(ligand_receptor_network_add, id='pairs')
ligand_receptor_network_add <- subset(ligand_receptor_network_add, value != 0)
ggplot(data = ligand_receptor_network_add, mapping = aes_string(x = "variable",y = "pairs")) + 
  geom_point(mapping = aes_string(color = "value", size=1)) + 
  labs(x = "Features", y = "pairs") + 
  theme_linedraw() +
  #scale_color_manual(values = "RdYlBu") +
  scale_color_distiller(palette = "RdYlBu") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
pairs_3 <- intersect(CCC_pairs_BM$V1,ligand_receptor_network_add$pairs)
write.table(pairs_3, "/home/yushiya/data/cd34/data/fig/CCC_pairs_3.txt", row.names = F, col.names = F, quote = F)
pairs_3 <- read.table("/home/yushiya/data/cd34/data/fig/CCC_pairs_3.txt")

#multinichenet 
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(nichenetr)
library(multinichenetr)
sum <- readRDS("/home/yushiya/data/cd34/data/fig/immune_stromal.rds")
Idents(sum) <- sum$celltype
sum <- subset(sum, subset= (tissue %in% c("BM","PB")) | (celltype %in% c("Fibro-MSC","APOD+ MSC","Osteoblast","Osteo-MSC","THY1+ MSC","Adipo-MSC","SEC","AEC")) )
organism = "human"
options(timeout = 120)
if(organism == "human"){
  lr_network_all = 
    readRDS(
      "/home/yushiya/data/reference/Nichnet_database/lr_network_human_allInfo_30112033.rds"
    ) %>% 
    mutate(
      ligand = convert_alias_to_symbols(ligand, organism = organism), 
      receptor = convert_alias_to_symbols(receptor, organism = organism))
  lr_network_all = lr_network_all  %>% 
    mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 
  lr_network = lr_network_all %>% 
    distinct(ligand, receptor)
  ligand_target_matrix = readRDS(
    "/home/yushiya/data/reference/Nichnet_database/ligand_target_matrix_nsga2r_final.rds"
  )
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]
} else if(organism == "mouse"){
  lr_network_all = readRDS(
    "/home/yushiya/data/reference/Nichnet_database/lr_network_mouse_allInfo_30112033.rds"
  ) %>% 
    mutate(
      ligand = convert_alias_to_symbols(ligand, organism = organism), 
      receptor = convert_alias_to_symbols(receptor, organism = organism))
  lr_network_all = lr_network_all  %>% 
    mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 
  lr_network = lr_network_all %>% 
    distinct(ligand, receptor)
  ligand_target_matrix = readRDS(
    "/home/yushiya/data/reference/Nichnet_database/ligand_target_matrix_nsga2r_final_mouse.rds"
  )
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]
}
#sce = readRDS(url(
#  "https://zenodo.org/record/8010790/files/sce_subset_misc.rds"
#))
sum@meta.data[(sum$celltype %in% c("Fibro-MSC","APOD+ MSC","Osteoblast","Osteo-MSC","THY1+ MSC","Adipo-MSC","SEC","AEC")),][,"tissue"] <- sample(c("BM", "PB"), size = length(sum@meta.data[(sum$celltype %in% c("Fibro-MSC","APOD+ MSC","Osteoblast","Osteo-MSC","THY1+ MSC","Adipo-MSC","SEC","AEC")),][,"tissue"]), replace = TRUE)
sum@meta.data[(sum$celltype %in% c("Fibro-MSC","APOD+ MSC","Osteoblast","Osteo-MSC","THY1+ MSC","Adipo-MSC","SEC","AEC")),][,"ori_stim"] <- sample(c("BM_SC1", "BM_SC2","BM_SC3"), size = length(sum@meta.data[(sum$celltype %in% c("Fibro-MSC","APOD+ MSC","Osteoblast","Osteo-MSC","THY1+ MSC","Adipo-MSC","SEC","AEC")),][,"tissue"]), replace = TRUE)
sum$celltype <- make.names(sum$celltype)
sce = Seurat::as.SingleCellExperiment(sum, assay = "RNA")
sce = alias_to_symbol_SCE(sce, "human") %>% makenames_SCE()
#Define metadata
sample_id = "ori_stim"
group_id = "tissue"
celltype_id = "celltype"
#batches = "ori_stim"
batches = NA
covariates = NA
#Define the contrasts of interest
contrasts_oi = c("'BM-PB','PB-BM'")  
contrast_tbl = tibble(contrast = c("BM-PB","PB-BM"), group = c("BM","PB"))
#Define the sender and receiver cell types of interest
#senders_oi <- c("Fibro-MSC","APOD+ MSC","Osteoblast","Osteo-MSC","THY1+ MSC","Adipo-MSC","SEC","AEC")
#receivers_oi <- c("HSC","MPP#1","MPP#2","LMPP#1","LMPP#2","BMEP","EryP","MkP","Ma/Eo/BaP","GMP","CDP","pre-pDC","CLP","pro-B","pre-B")
senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% 
            c(senders_oi, receivers_oi)
          ]
#core analysis
min_cells = 10
abundance_info = get_abundance_info(
  sce = sce, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, receivers_oi = receivers_oi, 
  batches = batches
)
abundance_info$abund_plot_sample
abundance_df_summarized = abundance_info$abundance_data %>% 
  mutate(keep = as.logical(keep)) %>% 
  group_by(group_id, celltype_id) %>% 
  summarise(samples_present = sum((keep)))
celltypes_absent_one_condition = abundance_df_summarized %>% 
  filter(samples_present == 0) %>% pull(celltype_id) %>% unique() 
# find truly condition-specific cell types by searching for cell types 
# truely absent in at least one condition
celltypes_present_one_condition = abundance_df_summarized %>% 
  filter(samples_present >= 2) %>% pull(celltype_id) %>% unique() 
# require presence in at least 2 samples of one group so 
# it is really present in at least one condition
condition_specific_celltypes = intersect(
  celltypes_absent_one_condition, 
  celltypes_present_one_condition)
total_nr_conditions = SummarizedExperiment::colData(sce)[,group_id] %>% 
  unique() %>% length() 
absent_celltypes = abundance_df_summarized %>% 
  filter(samples_present < 2) %>% 
  group_by(celltype_id) %>% 
  count() %>% 
  filter(n == total_nr_conditions) %>% 
  pull(celltype_id)
print("condition-specific celltypes:")
## [1] "condition-specific celltypes:"
print(condition_specific_celltypes)
## character(0)
print("absent celltypes:")
## [1] "absent celltypes:"
print(absent_celltypes)
## character(0)
analyse_condition_specific_celltypes = FALSE
if(analyse_condition_specific_celltypes == TRUE){
  senders_oi = senders_oi %>% setdiff(absent_celltypes)
  receivers_oi = receivers_oi %>% setdiff(absent_celltypes)
} else {
  senders_oi = senders_oi %>% 
    setdiff(union(absent_celltypes, condition_specific_celltypes))
  receivers_oi = receivers_oi %>% 
    setdiff(union(absent_celltypes, condition_specific_celltypes))
}
sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% 
            c(senders_oi, receivers_oi)
          ]
min_sample_prop = 0.50
fraction_cutoff = 0.05
frq_list = get_frac_exprs(
  sce = sce, 
  sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id, 
  batches = batches, 
  min_cells = min_cells, 
  fraction_cutoff = fraction_cutoff, min_sample_prop = min_sample_prop)
genes_oi = frq_list$expressed_df %>% 
  filter(expressed == TRUE) %>% pull(gene) %>% unique() 
sce = sce[genes_oi, ]
abundance_expression_info = process_abundance_expression_info(
  sce = sce, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, receivers_oi = receivers_oi, 
  lr_network = lr_network, 
  #batches = batches, 
  frq_list = frq_list, 
  abundance_info = abundance_info)
abundance_expression_info$celltype_info$pb_df %>% head()
abundance_expression_info$celltype_info$pb_df_group %>% head()
DE_info = get_DE_info(
  sce = sce, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  batches = batches, covariates = covariates, 
  contrasts_oi = contrasts_oi, 
  min_cells = min_cells, 
  expressed_df = frq_list$expressed_df)
DE_info$celltype_de$de_output_tidy %>% head()
DE_info$hist_pvals
empirical_pval = FALSE
if(empirical_pval == TRUE){
  DE_info_emp = get_empirical_pvals(DE_info$celltype_de$de_output_tidy)
  celltype_de = DE_info_emp$de_output_tidy_emp %>% select(-p_val, -p_adj) %>% 
    rename(p_val = p_emp, p_adj = p_adj_emp)
} else {
  celltype_de = DE_info$celltype_de$de_output_tidy
} 
sender_receiver_de = combine_sender_receiver_de(
  sender_de = celltype_de,
  receiver_de = celltype_de,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network
)
sender_receiver_de %>% head()
p_val_threshold = 0.05
p_val_adj = FALSE 
geneset_assessment = contrast_tbl$contrast %>% 
  lapply(
    process_geneset_data, 
    celltype_de, logFC_threshold, p_val_adj, p_val_threshold
  ) %>% 
  bind_rows() 
geneset_assessment_adjustedPval = contrast_tbl$contrast %>% 
  lapply(
    process_geneset_data, 
    celltype_de, logFC_threshold, p_val_adj = TRUE, p_val_threshold
  ) %>% 
  bind_rows() 
top_n_target = 250
verbose = TRUE
cores_system = 8
n.cores = min(cores_system, celltype_de$cluster_id %>% unique() %>% length()) 
ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(
  get_ligand_activities_targets_DEgenes(
    receiver_de = celltype_de,
    receivers_oi = intersect(receivers_oi, celltype_de$cluster_id %>% unique()),
    ligand_target_matrix = ligand_target_matrix,
    logFC_threshold = logFC_threshold,
    p_val_threshold = p_val_threshold,
    p_val_adj = p_val_adj,
    top_n_target = top_n_target,
    verbose = verbose, 
    n.cores = n.cores
  )
))
ligand_activities_targets_DEgenes$ligand_activities %>% head(20)
ligand_activity_down = FALSE
sender_receiver_tbl = sender_receiver_de %>% distinct(sender, receiver)
metadata_combined = SummarizedExperiment::colData(sce) %>% tibble::as_tibble()
if(!is.na(batches)){
  grouping_tbl = metadata_combined[,c(sample_id, group_id, batches)] %>% 
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group",batches)
} else {
  grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% 
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group")
}
prioritization_tables = suppressMessages(generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  scenario = "regular", # all prioritization criteria will be weighted equally
  fraction_cutoff = fraction_cutoff, 
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender,
  ligand_activity_down = ligand_activity_down
))
prioritization_tables$group_prioritization_tbl %>% head(20)
lr_target_prior_cor = lr_target_prior_cor_inference(
  receivers_oi = prioritization_tables$group_prioritization_tbl$receiver %>% unique(), 
  abundance_expression_info = abundance_expression_info, 
  celltype_de = celltype_de, 
  grouping_tbl = grouping_tbl, 
  prioritization_tables = prioritization_tables, 
  ligand_target_matrix = ligand_target_matrix, 
  logFC_threshold = logFC_threshold, 
  p_val_threshold = p_val_threshold, 
  p_val_adj = p_val_adj
)
path = "/home/yushiya/data/cd34/data/fig/nichnet_data/"
multinichenet_output = list(
  celltype_info = abundance_expression_info$celltype_info,
  celltype_de = celltype_de,
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de =  sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  prioritization_tables = prioritization_tables,
  grouping_tbl = grouping_tbl,
  lr_target_prior_cor = lr_target_prior_cor
) 
multinichenet_output = make_lite_output(multinichenet_output)
save = TRUE
if(save == TRUE){
  saveRDS(multinichenet_output, paste0(path, "multinichenet_test_output_1.rds"))
}
#draw results
prioritized_tbl_oi_all = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  top_n = 50, 
  rank_per_group = FALSE
)
prioritized_tbl_oi = 
  multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_oi_all$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>% 
  left_join(prioritized_tbl_oi_all)
prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0
senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique()) %>% sort()
colors_sender = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
colors_receiver = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
circos_list = make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)



#cellphoneDB prepare data for python
library(SeuratDisk)
sum <- readRDS("/home/yushiya/data/cd34/data/fig/immune_stromal.rds")
Idents(sum) <- sum$celltype
markers <- FindAllMarkers(sum, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(markers,file='/home/yushiya/data/cd34/data/fig/cellphoneDB_data/immune_stromal_DEG.csv')
markers <- read.csv("/home/yushiya/data/cd34/data/fig/cellphoneDB_data/immune_stromal_DEG.csv")
markers_meta <- markers[,c(7,8)]
colnames(markers_meta)[1] <- "cell_type"
write.table(markers_meta,file="/home/yushiya/data/cd34/data/fig/cellphoneDB_data/immune_stromal_DEG_meta.txt", row.names = F, sep = "\t")
sum@assays[["RNA"]]@scale.data <- matrix()
SaveH5Seurat(sum,filename="/home/yushiya/data/cd34/data/fig/cellphoneDB_data/immune_stromal.h5seurat", overwrite = TRUE)
Convert("/home/yushiya/data/cd34/data/fig/cellphoneDB_data/immune_stromal.h5seurat", dest = "h5ad", overwrite = TRUE)
#subset BM,PB
sum <- readRDS("/home/yushiya/data/cd34/data/fig/immune_stromal.rds")
Idents(sum) <- sum$celltype
sum <- subset(sum, subset= (tissue %in% c("BM")) | (celltype %in% c("Fibro-MSC","APOD+ MSC","Osteoblast","Osteo-MSC","THY1+ MSC","Adipo-MSC","SEC","AEC")))
sum@assays[["RNA"]]@scale.data <- matrix()
SaveH5Seurat(sum,filename="/home/yushiya/data/cd34/data/fig/cellphoneDB_data/immune_stromal_BM.h5seurat", overwrite = TRUE)
Convert("/home/yushiya/data/cd34/data/fig/cellphoneDB_data/immune_stromal_BM.h5seurat", dest = "h5ad", overwrite = TRUE)
write.table(sum$celltype,file="/home/yushiya/data/cd34/data/fig/cellphoneDB_data/immune_stromal_BM_meta.txt", sep = "\t", col.names = F, quote = F)
sum <- subset(sum, subset= (tissue %in% c("PB")) | (celltype %in% c("Fibro-MSC","APOD+ MSC","Osteoblast","Osteo-MSC","THY1+ MSC","Adipo-MSC","SEC","AEC")))
sum@assays[["RNA"]]@scale.data <- matrix()
SaveH5Seurat(sum,filename="/home/yushiya/data/cd34/data/fig/cellphoneDB_data/immune_stromal_PB.h5seurat", overwrite = TRUE)
Convert("/home/yushiya/data/cd34/data/fig/cellphoneDB_data/immune_stromal_PB.h5seurat", dest = "h5ad", overwrite = TRUE)
write.table(sum$celltype,file="/home/yushiya/data/cd34/data/fig/cellphoneDB_data/immune_stromal_PB_meta.txt", sep = "\t", col.names = F, quote = F)

#read results
library(reshape2)
re_method1 <- read.csv("/home/yushiya/data/cd34/data/fig/cellphoneDB_data/method1_BM/simple_analysis_interaction_scores_04_09_2025_141541.txt", sep = "\t")
re_method1 <- read.csv("/home/yushiya/data/cd34/data/fig/cellphoneDB_data/method2_BM/statistical_analysis_interaction_scores_04_09_2025_160442.txt", sep = "\t")
re_method1 <- read.csv("/home/yushiya/data/cd34/data/fig/cellphoneDB_data/method1_PB/simple_analysis_interaction_scores_04_09_2025_163255.txt", sep = "\t")
re_method1 <- read.csv("/home/yushiya/data/cd34/data/fig/cellphoneDB_data/method2_PB/statistical_analysis_interaction_scores_04_09_2025_171843.txt", sep = "\t")
re_method1 <- read.csv("/home/yushiya/data/cd34/data/fig/cellphoneDB_data/method3/degs_analysis_interaction_scores_03_16_2025_172734.txt", sep = "\t")
re_method_sub1 <- re_method1[,grepl("^Adipo",colnames(re_method1))]
re_method_sub1 <- re_method_sub1[,c(4:7,9:16,21:23)]
re_method_sub1 <- re_method_sub1[,c(6,9,10,7,8,1,4,12,11,5,2,14,3,15,13)]
re_method_sub2 <- re_method1[,grepl("^AEC",colnames(re_method1))]
re_method_sub2 <- re_method_sub2[,c(4:7,9:16,21:23)]
re_method_sub2 <- re_method_sub2[,c(6,9,10,7,8,1,4,12,11,5,2,14,3,15,13)]
re_method_sub3 <- re_method1[,grepl("^SEC",colnames(re_method1))]
re_method_sub3 <- re_method_sub3[,c(4:7,9:16,21:23)]
re_method_sub3 <- re_method_sub3[,c(6,9,10,7,8,1,4,12,11,5,2,14,3,15,13)]
re_method_sub4 <- re_method1[,grepl("^APOD",colnames(re_method1))]
re_method_sub4 <- re_method_sub4[,c(4:7,9:16,21:23)]
re_method_sub4 <- re_method_sub4[,c(6,9,10,7,8,1,4,12,11,5,2,14,3,15,13)]
re_method_sub5 <- re_method1[,grepl("^Fibro",colnames(re_method1))]
re_method_sub5 <- re_method_sub5[,c(4:7,9:16,21:23)]
re_method_sub5 <- re_method_sub5[,c(6,9,10,7,8,1,4,12,11,5,2,14,3,15,13)]
re_method_sub6 <- re_method1[,grepl("^THY",colnames(re_method1))]
re_method_sub6 <- re_method_sub6[,c(4:7,9:16,21:23)]
re_method_sub6 <- re_method_sub6[,c(6,9,10,7,8,1,4,12,11,5,2,14,3,15,13)]
re_method_sub7 <- re_method1[,grepl("^Osteo.MSC",colnames(re_method1))]
re_method_sub7 <- re_method_sub7[,c(4:7,9:16,21:23)]
re_method_sub7 <- re_method_sub7[,c(6,9,10,7,8,1,4,12,11,5,2,14,3,15,13)]
re_method_sub8 <- re_method1[,grepl("^Osteoblast",colnames(re_method1))]
re_method_sub8 <- re_method_sub8[,c(4:7,9:16,21:23)]
re_method_sub8 <- re_method_sub8[,c(6,9,10,7,8,1,4,12,11,5,2,14,3,15,13)]
re_method_sub <- cbind(re_method1[,c(1:13)],re_method_sub1,re_method_sub2,re_method_sub3,re_method_sub4,re_method_sub5,re_method_sub6,re_method_sub7,re_method_sub8)
re_method_sub <- re_method_sub[(rowMeans(re_method_sub[,c(14:133)]) > 0),]
#re_method_sub <- re_method_sub[(apply(re_method_sub[,c(14:133)],1,sd) > 10),]
write.table(re_method_sub, "/home/yushiya/data/cd34/data/fig/cellphoneDB_data/method1_PB/simple_analysis_filtered.txt", sep="\t")
write.table(re_method_sub, "/home/yushiya/data/cd34/data/fig/cellphoneDB_data/method2_PB/statistical_analysis_filtered.txt", sep="\t")
write.table(re_method_sub, "/home/yushiya/data/cd34/data/fig/cellphoneDB_data/method3/deg_analysis_filtered.txt", sep="\t")
re_draw <- re_method_sub[,c(14:133)]
rownames(re_draw) <- re_method_sub$interacting_pair
re_draw$pairs <- rownames(re_draw)
selected_pairs <- c("CD55_ADGRE5","ICAM1_ITGAL","ICAM1_SPN","PECAM1_CD38","SELE_GLG1","TGM2_ADGRG1")
selected_pairs <- c("CXCL12_CXCR4","APP_CD74","APP_SORL1")
selected_pairs <- c("CXCL12_CXCR4","APP_CD74","APP_SORL1","MDK_NCL")
re_draw <- re_draw[(re_draw$pairs %in% selected_pairs),]
re_draw <- re_draw[,c(16:45)]
re_draw <- re_draw[,c(1:15)]
re_draw <- re_draw[,c(76:90)]
re_draw$pairs <- selected_pairs
re_draw=melt(re_draw, id='pairs')
re_draw_filtered <- subset(re_draw, value != 0)
p1=ggplot(data = re_draw_filtered, mapping = aes_string(x = "variable",y = "pairs")) + 
  geom_point(mapping = aes_string(color = "value", size=1)) + 
  labs(x = "Features", y = "pairs") + 
  theme_linedraw() +
  #scale_color_manual(values = "RdYlBu") +
  scale_color_distiller(palette = "RdYlBu") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
ggsave("/home/yushiya/data/cd34/data/fig/fig3-cellphoneDB-allpairs-PB.pdf", plot=p1, width=25, height=40)
ggsave("/home/yushiya/data/cd34/data/fig/fig3-cellphoneDB-allpairs-PB_selected.pdf", plot=p1, width=25, height=6)
ggsave("/home/yushiya/data/cd34/data/fig/fig3-cellphoneDB-allpairs-PB_selected-1.pdf", plot=p1, width=8, height=3)
ggsave("/home/yushiya/data/cd34/data/fig/fig3-cellphoneDB-allpairs-BM_selected-3.pdf", plot=p1, width=5, height=2.5)
VlnPlot(sum, features = "CD38", pt.size = 0)
VlnPlot(sum, features = "PECAM1", pt.size = 0)
VlnPlot(sum, features = "SELE", pt.size = 0)
VlnPlot(sum, features = "GLG1", pt.size = 0)
VlnPlot(sum, features = "CD44", pt.size = 0)
VlnPlot(sum, features = "SELL", pt.size = 0)
VlnPlot(sum, features = "PODXL2", pt.size = 0) ICAM1_SPN
VlnPlot(sum, features = "ICAM1", pt.size = 0)
VlnPlot(sum, features = "SPN", pt.size = 0)
VlnPlot(sum, features = "CD55", pt.size = 0)
VlnPlot(sum, features = "ADGRE5", pt.size = 0)
VlnPlot(sum, features = "APP", pt.size = 0)
VlnPlot(sum, features = "SORL1", pt.size = 0)
re_draw_selected <- subset(re_draw, pairs %in% c("SELE_GLG1","SELE_CD44","PECAM1_CD38"," ICAM1_SPN","PODXL_SELL","CXCL12_CXCR4","ICAM1_SPN","CD55_ADGRE5","APP_SORL1","APP_CD74"))
p1=ggplot(data = re_draw_selected, mapping = aes_string(x = "variable",y = "pairs")) + 
  geom_point(mapping = aes_string(color = "value", size=1)) + 
  labs(x = "Features", y = "pairs") + 
  theme_linedraw() +
  #scale_color_manual(values = "RdYlBu") +
  scale_color_distiller(palette = "RdYlBu") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
ggsave("/home/yushiya/data/cd34/data/fig/fig3-cellphoneDB-allpairs-selected.pdf", plot=p1, width=25, height=5)

#Find overlap signal of three functions
p1=netVisual_bubble(cellchat_BM, sources.use = c(1:8), targets.use = c(9:23), remove.isolate = FALSE)
cellchat_pairs <- paste(p1[["data"]]$ligand,p1[["data"]]$receptor,sep="_")
cellchat_pairs <- unique(cellchat_pairs)
cellphoneDB_pairs <- re_method_sub$interacting_pair
pairs <- intersect(cellchat_pairs,cellphoneDB_pairs)
write.table(pairs, "/home/yushiya/data/cd34/data/fig/CCC_pairs_BM.txt", row.names = F, col.names = F, quote = F)




#heatmap for comparison of cellchat
library(pheatmap)
cellchat_BM <- readRDS(file = "/home/yushiya/cd34/fig/cellchat_BM.rds")
cellchat_mPB <- readRDS(file = "/home/yushiya/cd34/fig/cellchat_mPB.rds")
cellchat_PB <- readRDS(file = "/home/yushiya/cd34/fig/cellchat_PB.rds")
#"CXCL","SELPLG","ADGRE5"
pathways.show <- c("ADGRE5") 
mtx <- rbind(cellchat_BM@netP[["prob"]][1,,pathways.show],cellchat_mPB@netP[["prob"]][1,,pathways.show],cellchat_PB@netP[["prob"]][1,,pathways.show])
mtx <- rbind(cellchat_BM@netP[["prob"]][2,,pathways.show],cellchat_mPB@netP[["prob"]][2,,pathways.show],cellchat_PB@netP[["prob"]][2,,pathways.show])
rownames(mtx) <- c("BM","mPB","PB")
annotation_row <- data.frame(Stim = factor(rep(c("BM", "mPB","PB"), c(1,1,1))))
rownames(annotation_row) <- c("BM","mPB","PB")
annotation_colors =list(Stim=c(BM="#DB5C25",mPB="#F3B747",PB="#649541"))
p1=pheatmap(mtx, cluster_cols = F, cluster_rows = F,
            annotation_row = annotation_row,
            annotation_colors = annotation_colors,
            col=colorRampPalette(c("white", "#F35137","#530000"))(50),
            fontsize_row = 12, fontsize_col= 12, angle_col = 45)
ggsave("/home/yushiya/cd34/fig/fig3-cellchat-comp-ADGRE5.pdf", plot=p1, width=8, height=2.2)
