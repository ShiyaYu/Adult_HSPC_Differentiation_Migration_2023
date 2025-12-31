library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(devtools)
library(slingshot)
library(RColorBrewer)

MF1.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/Sample02_cellranger")
MF1 <- CreateSeuratObject(counts = MF1.data, project = "MF")
MF1$stim <- 'MF1'
MF1[["percent.mt"]] <- PercentageFeatureSet(MF1, pattern = "^MT-")
VlnPlot(MF1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MF1 <- subset(MF1, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(MF1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MF1 <- NormalizeData(MF1, normalization.method = "LogNormalize", scale.factor = 10000)
MF1 <- FindVariableFeatures(MF1, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(MF1)
MF1 <- ScaleData(MF1, features = all.genes)

MF2.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/Sample03_cellranger")
MF2 <- CreateSeuratObject(counts = MF2.data, project = "MF")
MF2$stim <- 'MF2'
MF2[["percent.mt"]] <- PercentageFeatureSet(MF2, pattern = "^MT-")
VlnPlot(MF2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MF2 <- subset(MF2, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(MF2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MF2 <- NormalizeData(MF2, normalization.method = "LogNormalize", scale.factor = 10000)
MF2 <- FindVariableFeatures(MF2, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(MF2)
MF2 <- ScaleData(MF2, features = all.genes)

MF3.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/Sample04_cellranger")
MF3 <- CreateSeuratObject(counts = MF3.data, project = "MF")
MF3$stim <- 'MF3'
MF3[["percent.mt"]] <- PercentageFeatureSet(MF3, pattern = "^MT-")
VlnPlot(MF3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MF3 <- subset(MF3, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(MF3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MF3 <- NormalizeData(MF3, normalization.method = "LogNormalize", scale.factor = 10000)
MF3 <- FindVariableFeatures(MF3, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(MF3)
MF3 <- ScaleData(MF3, features = all.genes)

MF4.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/Sample05_cellranger")
MF4 <- CreateSeuratObject(counts = MF4.data, project = "MF")
MF4$stim <- 'MF4'
MF4[["percent.mt"]] <- PercentageFeatureSet(MF4, pattern = "^MT-")
VlnPlot(MF4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MF4 <- subset(MF4, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(MF4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MF4 <- NormalizeData(MF4, normalization.method = "LogNormalize", scale.factor = 10000)
MF4 <- FindVariableFeatures(MF4, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(MF4)
MF4 <- ScaleData(MF4, features = all.genes)

MF5.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/Sample07_cellranger")
MF5 <- CreateSeuratObject(counts = MF5.data, project = "MF")
MF5$stim <- 'MF5'
MF5[["percent.mt"]] <- PercentageFeatureSet(MF5, pattern = "^MT-")
VlnPlot(MF5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MF5 <- subset(MF5, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(MF5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MF5 <- NormalizeData(MF5, normalization.method = "LogNormalize", scale.factor = 10000)
MF5 <- FindVariableFeatures(MF5, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(MF5)
MF5 <- ScaleData(MF5, features = all.genes)

MF6.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/Sample08_cellranger")
MF6 <- CreateSeuratObject(counts = MF6.data, project = "MF")
MF6$stim <- 'MF6'
MF6[["percent.mt"]] <- PercentageFeatureSet(MF6, pattern = "^MT-")
VlnPlot(MF6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MF6 <- subset(MF6, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(MF6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MF6 <- NormalizeData(MF6, normalization.method = "LogNormalize", scale.factor = 10000)
MF6 <- FindVariableFeatures(MF6, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(MF6)
MF6 <- ScaleData(MF6, features = all.genes)

MF7.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/Sample10_cellranger")
MF7 <- CreateSeuratObject(counts = MF7.data, project = "MF")
MF7$stim <- 'MF7'
MF7[["percent.mt"]] <- PercentageFeatureSet(MF7, pattern = "^MT-")
VlnPlot(MF7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MF7 <- subset(MF7, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(MF7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MF7 <- NormalizeData(MF7, normalization.method = "LogNormalize", scale.factor = 10000)
MF7 <- FindVariableFeatures(MF7, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(MF7)
MF7 <- ScaleData(MF7, features = all.genes)

MF8.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/Sample11_cellranger")
MF8 <- CreateSeuratObject(counts = MF8.data, project = "MF")
MF8$stim <- 'MF8'
MF8[["percent.mt"]] <- PercentageFeatureSet(MF8, pattern = "^MT-")
VlnPlot(MF8, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MF8 <- subset(MF8, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(MF8, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MF8 <- NormalizeData(MF8, normalization.method = "LogNormalize", scale.factor = 10000)
MF8 <- FindVariableFeatures(MF8, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(MF8)
MF8 <- ScaleData(MF8, features = all.genes)

MF9.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/Sample12_cellranger")
MF9 <- CreateSeuratObject(counts = MF9.data, project = "MF")
MF9$stim <- 'MF9'
MF9[["percent.mt"]] <- PercentageFeatureSet(MF9, pattern = "^MT-")
VlnPlot(MF9, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MF9 <- subset(MF9, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(MF9, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MF9 <- NormalizeData(MF9, normalization.method = "LogNormalize", scale.factor = 10000)
MF9 <- FindVariableFeatures(MF9, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(MF9)
MF9 <- ScaleData(MF9, features = all.genes)

MF10.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/Sample14_cellranger")
MF10 <- CreateSeuratObject(counts = MF10.data, project = "MF")
MF10$stim <- 'MF10'
MF10[["percent.mt"]] <- PercentageFeatureSet(MF10, pattern = "^MT-")
VlnPlot(MF10, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MF10 <- subset(MF10, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(MF10, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MF10 <- NormalizeData(MF10, normalization.method = "LogNormalize", scale.factor = 10000)
MF10 <- FindVariableFeatures(MF10, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(MF10)
MF10 <- ScaleData(MF10, features = all.genes)

MF11.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/Sample15_cellranger")
MF11 <- CreateSeuratObject(counts = MF11.data, project = "MF")
MF11$stim <- 'MF11'
MF11[["percent.mt"]] <- PercentageFeatureSet(MF11, pattern = "^MT-")
VlnPlot(MF11, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MF11 <- subset(MF11, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(MF11, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MF11 <- NormalizeData(MF11, normalization.method = "LogNormalize", scale.factor = 10000)
MF11 <- FindVariableFeatures(MF11, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(MF11)
MF11 <- ScaleData(MF11, features = all.genes)

MF12.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/Sample16_cellranger")
MF12 <- CreateSeuratObject(counts = MF12.data, project = "MF")
MF12$stim <- 'MF12'
MF12[["percent.mt"]] <- PercentageFeatureSet(MF12, pattern = "^MT-")
VlnPlot(MF12, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MF12 <- subset(MF12, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(MF12, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MF12 <- NormalizeData(MF12, normalization.method = "LogNormalize", scale.factor = 10000)
MF12 <- FindVariableFeatures(MF12, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(MF12)
MF12 <- ScaleData(MF12, features = all.genes)

MF13.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/Sample19_cellranger")
MF13 <- CreateSeuratObject(counts = MF13.data, project = "MF")
MF13$stim <- 'MF13'
MF13[["percent.mt"]] <- PercentageFeatureSet(MF13, pattern = "^MT-")
VlnPlot(MF13, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MF13 <- subset(MF13, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(MF13, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MF13 <- NormalizeData(MF13, normalization.method = "LogNormalize", scale.factor = 10000)
MF13 <- FindVariableFeatures(MF13, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(MF13)
MF13 <- ScaleData(MF13, features = all.genes)

MF14.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/Sample20_cellranger")
MF14 <- CreateSeuratObject(counts = MF14.data, project = "MF")
MF14$stim <- 'MF14'
MF14[["percent.mt"]] <- PercentageFeatureSet(MF14, pattern = "^MT-")
VlnPlot(MF14, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MF14 <- subset(MF14, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(MF14, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MF14 <- NormalizeData(MF14, normalization.method = "LogNormalize", scale.factor = 10000)
MF14 <- FindVariableFeatures(MF14, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(MF14)
MF14 <- ScaleData(MF14, features = all.genes)

MF15.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/Sample21_cellranger")
MF15 <- CreateSeuratObject(counts = MF15.data, project = "MF")
MF15$stim <- 'MF15'
MF15[["percent.mt"]] <- PercentageFeatureSet(MF15, pattern = "^MT-")
VlnPlot(MF15, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MF15 <- subset(MF15, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(MF15, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MF15 <- NormalizeData(MF15, normalization.method = "LogNormalize", scale.factor = 10000)
MF15 <- FindVariableFeatures(MF15, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(MF15)
MF15 <- ScaleData(MF15, features = all.genes)

MF.anchors <- FindIntegrationAnchors(object.list = c(MF1,MF2,MF3,MF4,MF5,MF6,MF7,MF8,MF9,MF10,MF11,MF12,MF13,MF14,MF15), dims = 1:20)
MF <- IntegrateData(anchorset = MF.anchors, dims = 1:20)
DefaultAssay(MF) <- "integrated"
MF@meta.data$ori_stim <- MF@meta.data$stim
MF@meta.data$stim <- "MF"
# Run the standard workflow for visualization and clustering
MF <- ScaleData(MF, verbose = FALSE)
MF <- RunPCA(MF, npcs = 30, verbose = FALSE)
ElbowPlot(MF)
#Clustering
MF <- FindNeighbors(MF, reduction = "pca", dims = 1:30)
MF <- FindClusters(MF, resolution = 0.8)
# umap
MF <- RunUMAP(MF, reduction = "pca", dims = 1:30)
DimPlot(MF, reduction = "umap", group.by = "ori_stim")
DimPlot(MF, reduction = "umap", label = TRUE)
saveRDS(MF, file = "/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/MF.rds")
saveRDS(MF, file = "/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/MF_mapped.rds")
saveRDS(MF, file = "/home/yushiya/data3/cd34/data/fig/MF_mapped.rds")

MF <- readRDS("/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/MF.rds")
MF <- readRDS("/home/yushiya/data/cd34/data/fig/MF_mapped.rds")
immune.combined <- readRDS("/home/yushiya/data/cd34/data/immune.combined-33-celltype.rds")
#immune.combined <- readRDS("/home/yushiya/cd34/fig/immune_combined-3I-2-celltype.rds")
#map and annotate MF data to reference
immune.combined1 <- RunUMAP(immune.combined, dims = 1:40, reduction = "pca", return.model = TRUE)
immune.combined1@reductions[["umap"]]@cell.embeddings <- immune.combined@reductions[["umap"]]@cell.embeddings
immune.combined1@reductions[["umap"]]@misc[["model"]][["embedding"]] <- immune.combined@reductions[["umap"]]@cell.embeddings
immune.combined1 <- FindVariableFeatures(immune.combined1, selection.method = "vst", nfeatures = 2000)
DimPlot(immune.combined1, reduction = "umap", label = T)
all.anchors <- FindTransferAnchors(reference = immune.combined1, query = MF,
                                   dims = 1:30, reference.reduction = "pca")
MF <- MapQuery(anchorset = all.anchors, reference = immune.combined1, query = MF,refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap")
MF$predicted.celltype <- factor(MF$predicted.celltype, levels=c('HSC','MPP#1',"MPP#2",'LMPP#1',"LMPP#2",'BMEP','EryP','MkP',"Ma/Eo/BaP",'GMP','CDP','pre-pDC','CLP','pro-B','pre-B',"ETP-like","Thy#1","Thy#2","Thy#3"))
MF <- subset(MF, subset=predicted.celltype %in% c('HSC','MPP#1',"MPP#2",'LMPP#1',"LMPP#2",'BMEP','EryP','MkP',"Ma/Eo/BaP",'GMP','CDP','pre-pDC','CLP','pro-B','pre-B',"ETP-like"))
color_MF <- c("#550000", "#AA3939", "#CC6F66","#E79492","#FBD2CE","#F2B662","#F38F44","#B97802","#915900","#49AFB8","#3465A0","#2C2A50","#A5D6A7", "#47AF50", "#2E7D32","#8290BC")
color_stim <- c("#DB5C25","#F3B747","#649541","#BC80BD","#4F3763")
color_stim <- c("#DB5C25","#BC80BD","#F3B747","#4F3763")
p1=DimPlot(MF, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
           cols=color_MF, repel = TRUE,raster = T, raster.dpi = c(300,300))
ggsave("/home/yushiya/data/cd34/data/fig/fig6-MF-umap.pdf", plot=p1, width=6, height=4)
MF$celltype <- MF$predicted.celltype
DefaultAssay(MF) <- "RNA"
p1=FeaturePlot(MF, reduction = "ref.umap", features = c("AVP","CRHBP","MPO","IRF8","GATA2","PPBP"), max.cutoff = 3, cols = c("#C5E9F5", "#981C12"), raster = T)
ggsave("/home/yushiya/data3/cd34/data/fig/fig6-MF-Fea.pdf", plot=p1, width=6, height=7)
saveRDS(MF,"/home/yushiya/data/cd34/data/fig/MF_mapped.rds")


#merge MF and healthy samples
immune.combined$UMAP_1 <- immune.combined@reductions[["umap"]]@cell.embeddings[,1]
immune.combined$UMAP_2 <- immune.combined@reductions[["umap"]]@cell.embeddings[,2]
MF$UMAP_1 <- MF@reductions[["ref.umap"]]@cell.embeddings[,1]
MF$UMAP_2 <- MF@reductions[["ref.umap"]]@cell.embeddings[,2]
total <- merge(immune.combined,MF)
total <- subset(total, subset= stim %in% c("BM","mPB","PB","SP","MF"))
total$stim <- factor(total$stim, levels=c("BM","mPB","PB","SP","MF"))
total <- subset(total, subset= celltype %in% c('HSC','MPP#1',"MPP#2",'LMPP#1',"LMPP#2",'BMEP','EryP','MkP',"Ma/Eo/BaP",'GMP','CDP','pre-pDC','CLP','pro-B','pre-B',"ETP-like"))
total$celltype <- factor(total$celltype, levels = c('HSC','MPP#1',"MPP#2",'LMPP#1',"LMPP#2",'BMEP','EryP','MkP',"Ma/Eo/BaP",'GMP','CDP','pre-pDC','CLP','pro-B','pre-B',"ETP-like"))
Idents(total) <- total$celltype
total <- NormalizeData(total, normalization.method = "LogNormalize", scale.factor = 10000)
total <- FindVariableFeatures(total, selection.method = "vst", nfeatures = 2000)
total <- ScaleData(total, verbose = FALSE)
total <- RunPCA(total, npcs = 50, verbose = FALSE)
total <- RunUMAP(total, reduction = "pca", dims = 1:30)
total <- FindNeighbors(total, reduction = "pca", dims = 1:30)
total <- FindClusters(total, resolution = 0.5)
total@reductions[["umap"]]@cell.embeddings[,1] <- total$UMAP_1
total@reductions[["umap"]]@cell.embeddings[,2] <- total$UMAP_2
DimPlot(total, reduction = "umap", label = TRUE, group.by = "celltype")
total <- subset(total, subset= stim=="BM"|stim=="mPB"|stim=="SP"|stim=="MF")
#total <- subset(total, subset= stim=="BM"|stim=="PB"|stim=="MF")
total$stim <- factor(total$stim, levels=c("BM","SP","mPB","MF"))
saveRDS(total, file = "/home/yushiya/data/cd34/data/fig/MF_merged.rds")
total <- readRDS(file = "/home/yushiya/data/cd34/data/fig/MF_merged.rds")

#calculate Migratory score
library(ggpubr)
positive_regu_selected <- c("MPP1","MSN","CD99","RAC2","MAPK1","APP","PLCB1","ZNF580","TGFB1","CALR",
                            "TPGS1","C1QBP","CRK","PYCARD","MAPK3","SPN","ADAM10","RIN3","NCKAP1L","DNM1L",
                            "CD9","WNK1","VEGFB","DNAJC4","SPI1","MDK","RHOG","CD81","ADAM8",
                            "CAMK1D","ANXA1","DOCK8","LYN","TGS1","HOXA7","GPSM3","RIPOR2","CCL28","PTGER4",
                            "RHOH","CD47","SELENOK","INPP5D","ITGA4","ANXA4","RTN4","ADAM17","GCSAML","MIA3")
total<-AddModuleScore(object = total, features = list(Migratory_Score=positive_regu_selected), name="Migratory_Score")
VlnPlot(total, features = c("Migratory_Score1"), group.by="celltype", pt.size = 0,combine = FALSE)
miGene <- as.data.frame(colnames(total))
miGene$Tissue <- total$stim
miGene$celltype <- total$celltype
miGene$celltype <- factor(miGene$celltype,levels = c('HSC','MPP#1',"MPP#2",'LMPP#1',"LMPP#2",'BMEP','EryP','MkP',"Ma/Eo/BaP",'GMP','CDP','pre-pDC','CLP','pro-B','pre-B',"ETP-like"))
#miGene$Tissue <- factor(miGene$Tissue, levels=c("BM","mPB","PB","SP","MF"))
miGene$Tissue <- factor(miGene$Tissue, levels=c("BM","SP","mPB","MF"))
names(miGene)[1]='ID'
miGene$gene='Migratory_Score1'
miGene$Migr_Score <- total$Migratory_Score1
p1=ggboxplot(miGene, x="celltype", y="Migr_Score", color="Tissue", bxp.errorbar = T,
          palette = color_stim, add="boxplot", short.panel.labs = TRUE, ncol=4, add.params=list(size=0.2))+
  stat_compare_means(label = "p.signif",method = "t.test",ref.group = ".all.")+
  theme(axis.text.x =element_text(angle = 45,hjust = 1),
        plot.title = element_text(hjust = 0))
ggsave("/home/yushiya/data/cd34/data/fig/fig6-Migr-all-1.pdf", plot=p1, width=11, height=4)
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
clus_percent=melt(clus_percent, id='Tissue')
names(clus_percent)[2]='Tissue'
names(clus_percent)[1]='Celltype'
clus_percent$Celltype <- factor(clus_percent$Celltype,levels = c('HSC','MPP#1',"MPP#2",'LMPP#1',"LMPP#2",'BMEP','EryP','MkP',"Ma/Eo/BaP",'GMP','CDP','pre-pDC','CLP','pro-B','pre-B',"ETP-like"))
clus_percent$Stim <- factor(clus_percent$Stim, levels=c("BM","mPB","PB","SP","MF"))
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
        axis.line = element_line(colour = "grey90",size = rel(1)))
ggsave("/home/yushiya/data/cd34/data/fig/fig6-MF-clusP-all.pdf", plot=p1, width=6, height=3)

#calculate lineage scores
library(readxl)
celltype_markers <- read.csv("/home/yushiya/data/cd34/data/fig/DEG_celltype.csv",stringsAsFactors = F)
LT_HSC_genes <- read.csv("//home/yushiya/data/cd34/data/fig/LT_HSC_genes.csv")
HSC_gene <- subset(celltype_markers, subset= cluster=="HSC")
HSC_gene <- HSC_gene[HSC_gene$gene %in% LT_HSC_genes$Gene.Symbol,]
HSC_gene <- HSC_gene %>% top_n(n = 50, wt = avg_log2FC)
total<-AddModuleScore(object = total, features = list(HSC_Score=HSC_gene$gene), name="HSC_Score")
p1=VlnPlot(total, features = c("HSC_Score1"), cols=color_stim, pt.size = 0, split.by = "stim", group.by = "celltype") 
VlnPlot(MF, features = c("HSC_Score1"), pt.size = 0, group.by = "celltype")
ggsave(file="/home/yushiya/data/cd34/data/fig/fig6-MF-Score-HSC-1.pdf", plot=p1, width=8, height=3)
EryP_gene <- subset(celltype_markers, subset= cluster=="EryP")
EryP_gene <- EryP_gene %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
total<-AddModuleScore(object = total, features = list(EryP_Score=EryP_gene$gene), name="EryP_Score")
p1=VlnPlot(total, features = c("EryP_Score1"), cols=color_stim, pt.size = 0, split.by = "stim", group.by = "celltype")
ggsave(file="/home/yushiya/data/cd34/data/fig/fig6-MF-Score-EryP-1.pdf", plot=p1, width=8, height=3)
MkP_gene <- subset(celltype_markers, subset= cluster=="MkP")
MkP_gene <- MkP_gene %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
total<-AddModuleScore(object = total, features = list(MkP_Score=MkP_gene$gene), name="MkP_Score")
p1=VlnPlot(total, features = c("MkP_Score1"), cols=color_stim, pt.size = 0, split.by = "stim", group.by = "celltype")
ggsave(file="/home/yushiya/data/cd34/data/fig/fig6-MF-Score-Mk-1.pdf", plot=p1, width=8, height=3)
My_gene <- subset(celltype_markers, subset= cluster=="CDP" | cluster=="pre-pDC")
My_gene <- My_gene %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
total<-AddModuleScore(object = total, features = list(My_Score=My_gene$gene), name="My_Score")
p1=VlnPlot(total, features = c("My_Score1"), cols=color_stim, pt.size = 0, split.by = "stim", group.by = "celltype") 
ggsave(file="/home/yushiya/data/cd34/data/fig/fig6-MF-Score-My-1.pdf", plot=p1, width=8, height=3)
B_gene <- subset(celltype_markers, subset= cluster=="pro-B"|cluster=="pre-B")
B_gene <- B_gene %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
total<-AddModuleScore(object = total, features = list(B_Score=B_gene$gene), name="B_Score")
p1=VlnPlot(total, features = c("B_Score1"), cols=color_stim, pt.size = 0, split.by = "stim", group.by = "celltype")
ggsave(file="/home/yushiya/data/cd34/data/fig/fig6-MF-Score-B-1.pdf", plot=p1, width=8, height=3)
T_gene <- subset(celltype_markers, subset= cluster=="Thy#1" | cluster=="Thy#2"| cluster=="Thy#3")
T_gene <- T_gene %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
total<-AddModuleScore(object = total, features = list(T_Score=T_gene$gene), name="T_Score")
p1=VlnPlot(total, features = c("T_Score1"), cols=color_stim, pt.size = 0, split.by = "stim", group.by = "celltype")
ggsave(file="/home/yushiya/data/cd34/data/fig/fig6-MF-Score-T-1.pdf", plot=p1, width=8, height=3)


#calculate DEGs
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(GOSemSim)
library(DOSE)
celltype.markers <- FindAllMarkers(MF, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(celltype.markers, "/home/yushiya/data3/cd34/data/fig/MF_celltype_markers.csv")
Idents(total) <- total$stim
stim.markers <- FindAllMarkers(total, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(stim.markers, "/home/yushiya/data3/cd34/data/fig/MF_all_stim_markers.csv")
type <- c('HSC','MPP#1',"MPP#2",'LMPP#1',"LMPP#2",'BMEP','EryP','MkP','GMP','CDP','pre-pDC','CLP')
type <- c("Ma/Eo/Ba cell")
type <- c("CLP","GMP-like","CDP","pre-cDC","pre-pDC","ETP")
for (i in 1:length(type)) {
  sub1 <- subset(total, subset= celltype==type[i])
  Idents(sub1) <- sub1$stim
  cp_stim.markers <- FindAllMarkers(sub1, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
  write.csv(cp_stim.markers, paste("/home/yushiya/data/cd34/cd34/data/fig/MF_cpstim_",type[i],".csv",sep=""))
}
stim.markers <- read.csv("/home/yushiya/data3/cd34/data/fig/cpstim_MkP.csv",stringsAsFactors = F)
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
stim.markers <- read.csv("/home/yushiya/data/cd34/data/fig/MF_cpstim_HSC.csv",stringsAsFactors = F)
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
#gse@result$cluster="HSC"
#gsea_all <- gse@result
gsea_all <- ""
type <- c("HSC","MPP#1","MPP#2","LMPP#1","LMPP#2","BMEP","MkP","CLP","GMP","CDP","pre-pDC")
type <- c("HSC")
type <- c("Ma/Eo/Ba cell")
for (i in 1:length(type)) {
  stim.markers <- read.csv(paste("/home/yushiya/data/cd34/data/fig/MF_cpstim_",type[i],".csv",sep=""),stringsAsFactors = F)
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
gsea_all <- gsea_all[c(2:length(gsea_all$ID)),]
write.csv(gsea_all, "/home/yushiya/data/cd34/data/fig/gsea_all.csv")
gsea_all <- read.csv("/home/yushiya/data/cd34/data/fig/gsea_all.csv")
top <- gsea_all %>% group_by(cluster) %>% top_n(n = -5, wt = pvalue)
top_gsea <- top$Description %>% unique()
mtx_gsea <- data.frame(top_gsea)
colnames(mtx_gsea)[1] <- "Description"
type <- c("HSC","MPP#1","MPP#2","LMPP#1","LMPP#2","BMEP","MkP","CLP","GMP","CDP","pre-pDC")
for (i in 1:length(type)) {
  temp <- subset(gsea_all, subset= cluster==type[i])
  temp <- temp[,c("Description","enrichmentScore")]
  temp$enrichmentScore <- as.numeric(temp$enrichmentScore)
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
ggsave("/home/yushiya/data/cd34/data/fig/fig6-MF-gse-all-heatmap.pdf", plot=p1, width=9, height=5)


#DEGs among MF and BM, mPB, PB
#Idents(total) <- total$stim
#stim.markers <- FindAllMarkers(total, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1)
#write.csv(stim.markers, "/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/MF_all_stim_markers.csv")
#stim.markers <- read.csv("/home/yushiya/data/cd34/data/2020_Mead_PBMC/MF/MF_all_stim_markers.csv", stringsAsFactors = F)
HSC_markers <- read.csv("/home/yushiya/data/cd34/data/fig/MF_cpstim_HSC.csv", stringsAsFactors = F)
MPP1_markers <- read.csv("/home/yushiya/data/cd34/data/fig/MF_cpstim_MPP#1.csv", stringsAsFactors = F)
MPP2_markers <- read.csv("/home/yushiya/data/cd34/data/fig/MF_cpstim_MPP#2.csv", stringsAsFactors = F)
BMEP_markers <- read.csv("/home/yushiya/data/cd34/data/fig/MF_cpstim_BMEP.csv", stringsAsFactors = F)
EryP_markers <- read.csv("/home/yushiya/data/cd34/data/fig/MF_cpstim_EryP.csv", stringsAsFactors = F)
MkP_markers <- read.csv("/home/yushiya/data/cd34/data/fig/MF_cpstim_MkP.csv", stringsAsFactors = F)
Ma_markers <- read.csv("/home/yushiya/data/cd34/data/fig/MF_cpstim_Ma.csv", stringsAsFactors = F)
stim.markers <- rbind(HSC_markers,MPP1_markers,MPP2_markers,BMEP_markers,EryP_markers,MkP_markers,Ma_markers)
#find chemokines in DEGs
chemokines <- read.table("/home/yushiya/cd34/cytokines.txt")
chemo_genes <- stim.markers[stim.markers$gene %in% chemokines$V1,]
#chemo_BM <- subset(chemo_genes, subset= cluster=="BM")
#chemo_BM <- chemo_BM[order(chemo_BM$p_val_adj),]
#chemo_BM <- chemo_BM[1:5,]
#chemo_MF <- subset(chemo_genes, subset= cluster=="MF")
#chemo_MF <- chemo_MF[order(chemo_MF$p_val_adj),]
#chemo_MF <- chemo_MF[1:21,]
#chemo <- rbind(chemo_BM,chemo_MF)
chemo <- subset(chemo_genes, subset= cluster %in% c("mPB","SP","MF"))
chemo <- chemo %>% group_by(cluster) %>% top_n(n = -5, wt = p_val_adj)
gene <- chemo$gene %>% unique()
gene <- c("IL18","TNFSF10","CXCR4","FLT3","IFNGR1","LAP3","TGFB1","CMTM6","CCL5","PF4","PPBP","CXCL8","CXCL2","CXCL3","MIF","IL1B")
#heatmap of DEGs in HSPC and Ery/Mk lineage
sub1 <- subset(total, subset= stim %in% c("BM","mPB","SP","MF"))
sub1$stim_celltype <- paste(sub1$stim,sub1$celltype,sep="_")
sub1 <- subset(sub1, subset= celltype=="HSC" |celltype=="MPP#1" |celltype=="MPP#2" |celltype=="BMEP" |celltype=="EryP" |celltype=="MkP" |celltype=="Ma/Eo/BaP")
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
lev<-c("BM_HSC","SP_HSC","mPB_HSC","MF_HSC",
       "BM_MPP#1","SP_MPP#1","mPB_MPP#1","MF_MPP#1","BM_MPP#2","SP_MPP#2","mPB_MPP#2","MF_MPP#2",
       "BM_BMEP","SP_BMEP","mPB_BMEP","MF_BMEP",
       "BM_EryP","SP_EryP","mPB_EryP","MF_EryP","BM_MkP","SP_MkP","mPB_MkP","MF_MkP",
       "BM_Ma/Eo/BaP","SP_Ma/Eo/BaP","mPB_Ma/Eo/BaP","MF_Ma/Eo/BaP")
top_exp <- top_exp[,lev] 
annotation_col <- data.frame(Tissue = factor(rep(c("BM", "mPB","PB","MF"), c(6))))
annotation_col <- data.frame(Tissue = factor(rep(c("BM","SP","mPB","MF"), c(7))))
annotation_col <- data.frame(Tissue = factor(rep(c("BM", "PB","MF"), c(6))))
rownames(annotation_col) <- lev
annotation_colors =list(Tissue=c(BM="#DB5C25",mPB="#F3B747",PB="#649541",MF="#4F3763"))
annotation_colors =list(Tissue=c(BM="#DB5C25",SP="#BC80BD",mPB="#F3B747",MF="#4F3763"))
p1=pheatmap(top_exp, cluster_cols = F, #cluster_rows = F,  
            clustering_method = "average",
            annotation_col = annotation_col,
            annotation_colors = annotation_colors,
            scale = "row",
            #border=F,
            col=colorRampPalette(c("navy","white" ,"firebrick3"))(50),
            fontsize_row = 12, fontsize_col= 12, angle_col = 45)
ggsave("/home/yushiya/data/cd34/data/fig/fig6-MF-chemo-1.pdf", plot=p1, width=9, height=4)


#cellchat with BM stromal cells
library(CellChat)
library(svglite)
sub_sc <- readRDS(file = "/home/yushiya/data/cd34/data/stromal_sub.rds")
sub_tec <- readRDS(file = "/home/yushiya/data/cd34/data/TEC_sub.rds")
#MF <- subset(total, subset= stim=="MF")
Idents(MF) <- MF$predicted.celltype
MF <- subset(MF, subset=predicted.celltype %in% c('HSC','MPP#1',"MPP#2",'LMPP#1',"LMPP#2",'BMEP','EryP','MkP',"Ma/Eo/BaP",'GMP','CDP','pre-pDC'))
#mPB,SP
sub1 <- subset(immune.combined, subset= tissue %in% 'mPB')
sub1 <- subset(immune.combined, subset= tissue %in% 'SP')
sub1 <- subset(MF, subset=predicted.celltype %in% c('HSC','MPP#1',"MPP#2",'LMPP#1',"LMPP#2",'BMEP','EryP','MkP',"Ma/Eo/BaP",'GMP','CDP','pre-pDC','CLP','pro-B','pre-B',))
#ETP <- subset(MF, subset= predicted.celltype=="ETP")
sum <- merge(MF,sub_sc)
sum <- merge(sub1,sub_sc)
#sum <- merge(ETP,c(sub_sc,sub_tec))
cellchat <- createCellChat(sum)
cellchat@idents <- factor(cellchat@idents, levels=c("BM_SC","BM_Endo",'HSC','MPP#1',"MPP#2",'LMPP#1',"LMPP#2",'BMEP','EryP','MkP',"Ma/Eo/BaP",'GMP','CDP','pre-pDC',"CLP","pro-B","pre-B","ETP-like"))
#cellchat@idents <- factor(cellchat@idents, levels=c("Fibro-MSC","APOD+ MSC","Osteoblast","Osteo-MSC","THY1+ MSC","Adipo-MSC","SEC","AEC",'HSC','MPP#1',"MPP#2",'LMPP#1',"LMPP#2",'BMEP','EryP','MkP',"Ma/Eo/BaP",'GMP','CDP','pre-pDC'))
#cellchat@idents <- factor(cellchat@idents, levels=c("BM_SC","BM_Endo",'HSC','MPP#1',"MPP#2",'LMPP#1',"LMPP#2",'BMEP','EryP','MkP',"Ma/Eo/BaP",'GMP','CDP','pre-pDC'))
#cellchat@idents <- factor(cellchat@idents, levels=c("ETP","BM_SC","BM_Endo","Thymic_Epi","Thymic_Mesen","Thymic_Endo"))
CellChatDB <- CellChatDB.human 
CellChatDB.use <- CellChatDB 
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
#识别细胞组中过度表达的配体或受体
cellchat <- identifyOverExpressedGenes(cellchat)
#识别过度表达的配体受体相互作用
cellchat <- identifyOverExpressedInteractions(cellchat)
#将基因表达数据投射到PPI网络上
#cellchat <- projectData(cellchat, PPI.human)
#cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- computeCommunProb(cellchat, type = "triMean")
# 如果在某些细胞群中只有少数细胞，则过滤掉细胞间的通信
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
#par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
mat <- cellchat@net$weight
#par(mfrow = c(3,4), xpd=TRUE)
#dev.off()
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat2[c(1:2),c(3:14)] <- mat[c(1:2),c(3:14)]
netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = "Interaction weights/strength")
ggsave("/home/yushiya/data/cd34/data/fig/fig6-cellchat-all-circle.pdf", plot=p1, width=6, height=6)
# 显示从某些细胞组到其他细胞组的所有显著的相互作用（L-R 对）
p1=netVisual_bubble(cellchat, sources.use = c(1:8), targets.use = c(9:20), remove.isolate = FALSE)
ggsave("/home/yushiya/data/cd34/data/fig/fig6-cellchat-LRpairs-1.pdf", plot=p1, width=12, height=14)
p1=netVisual_bubble(cellchat, sources.use = c(5), targets.use = c(9:20), signaling = c("CXCL","MK","APP"), remove.isolate = FALSE)
ggsave("/home/yushiya/data/cd34/data/fig/fig6-cellchat_sig_select-1.pdf", plot=p1, width=4.5, height=3)
pathways.show <- c("MK") 
netVisual_aggregate(cellchat, signaling = pathways.show, sources.use = c(1,2), targets.use = c(3:14))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds", sources.use = c(1:2), targets.use = c(3:14))
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # “netP”:推断出的信号通路的细胞间通信网络
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 18, height = 8, font.size = 14, font.size.title = 18)
netAnalysis_contribution(cellchat, signaling = pathways.show)
ggsave("/home/yushiya/fig/fig6-cellchat-CXCL-circle.pdf", plot=p1, width=6, height=6)
ggsave("/home/yushiya/fig/fig6-cellchat-CXCL-heatmap.pdf", plot=p1, width=6, height=10)
ggsave("/home/yushiya/fig/fig6-cellchat-CXCL-role.pdf", plot=p1, width=9, height=5)
ggsave("/home/yushiya/fig/fig6-cellchat-CXCL-ctb.pdf", plot=p1, width=3.5, height=4)
cellchat <- readRDS(file = "/home/yushiya/data/cd34/data/fig/cellchat_MF.rds")
saveRDS(cellchat, file = "/home/yushiya/data/cd34/data/fig/cellchat_MF-1.rds")

#comparison of cellchat from MF and mPB
cellchat_MF <- readRDS(file = "/home/yushiya/data/cd34/data/fig/cellchat_MF.rds")
cellchat_mPB <- readRDS(file = "/home/yushiya/data/cd34/data/fig/cellchat_MF_mPB.rds")


#meta data of MF
Idents(MF) <- MF$ori_stim
MF <- RenameIdents(MF, "MF1"=12,"MF2"=7,"MF3"=14,"MF4"=21,"MF5"=8,"MF6"=13,"MF7"=20,
                   "MF8"=11,"MF9"=7,"MF10"=13,"MF11"=6,"MF12"=6,"MF13"=14.7,"MF14"=22,"MF15"=22.4)
MF$SP_size <- Idents(MF)
MF$SP_size <- factor(MF$SP_size, levels=c(6,7,8,11,12,13,14,14.7,20,21,22,22.4))
Idents(MF) <- MF$ori_stim
MF <- RenameIdents(MF, "MF1"=3,"MF2"=3,"MF3"=3,"MF4"="2-3","MF5"=3,"MF6"=3,"MF7"=3,
                   "MF8"=3,"MF9"=3,"MF10"=3,"MF11"="2-3","MF12"=2,"MF13"=2,"MF14"=3,"MF15"=3)
MF$Fibrosis_grade <- Idents(MF)
MF$Fibrosis_grade <- factor(MF$Fibrosis_grade, levels=c("2","2-3","3"))
#scores
miGene <- as.data.frame(colnames(MF))
miGene$Tissue <- MF$stim
miGene$celltype <- MF$predicted.celltype
miGene$Fibrosis_grade <- MF$Fibrosis_grade
miGene$SP_size <- MF$SP_size
miGene$celltype <- factor(miGene$celltype,levels = c('HSC','MPP#1',"MPP#2",'LMPP#1',"LMPP#2",'BMEP','EryP','MkP',"Ma/Eo/BaP",'GMP','CDP','pre-pDC','CLP','pro-B','pre-B',"ETP-like"))
#miGene$Tissue <- factor(miGene$Tissue, levels=c("BM","mPB","PB","SP","MF"))
miGene$Tissue <- factor(miGene$Tissue, levels=c("BM","SP","mPB","MF"))
miGene$Fibrosis_grade <- factor(miGene$Fibrosis_grade, levels=c("2","2-3","3"))
miGene$SP_size <- factor(miGene$SP_size, levels=c(6,7,8,11,12,13,14,14.7,20,21,22,22.4))
names(miGene)[1]='ID'
miGene$Migr_Score <- MF$Migratory_Score1
p1=ggboxplot(miGene, x="celltype", y="Migr_Score", color="Fibrosis_grade", bxp.errorbar = T,
             palette = c("#BCBDDC","#6A51A3","#3F007D"), add="boxplot", short.panel.labs = TRUE, ncol=4, add.params=list(size=0.2))+
  stat_compare_means(label = "p.signif",method = "t.test",ref.group = ".all.")+
  theme(axis.text.x =element_text(angle = 45,hjust = 1),
        plot.title = element_text(hjust = 0))
ggsave("/home/yushiya/data/cd34/data/fig/fig6-Migr-fibro-1.pdf", plot=p1, width=11, height=4)
p1=ggboxplot(miGene, x="celltype", y="Migr_Score", color="celltype", bxp.errorbar = T,
             palette = color_MF, add="boxplot", short.panel.labs = TRUE, ncol=4, add.params=list(size=0.2))+
  stat_compare_means(label = "p.signif",method = "t.test",ref.group = ".all.")+
  theme(axis.text.x =element_text(angle = 45,hjust = 1),
        plot.title = element_text(hjust = 0))
ggsave("/home/yushiya/data/cd34/data/fig/fig6-Migr-celltype-1.pdf", plot=p1, width=11, height=5)
p1=ggboxplot(miGene, x="SP_size", y="Migr_Score", color="SP_size", bxp.errorbar = T,
             add="boxplot", short.panel.labs = TRUE, ncol=4, add.params=list(size=0.2))+
  stat_compare_means(label = "p.signif",method = "t.test",ref.group = ".all.")+
  theme(axis.text.x =element_text(angle = 45,hjust = 1),
        plot.title = element_text(hjust = 0))
ggsave("/home/yushiya/data/cd34/data/fig/fig6-Migr-SPsize-1.pdf", plot=p1, width=11, height=4.5)



