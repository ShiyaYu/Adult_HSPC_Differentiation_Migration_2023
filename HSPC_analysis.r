library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(devtools)
library(slingshot)
library(RColorBrewer)
#library(DoubletFinder)
#library(circlize)
library(SeuratDisk)
library(hdf5r)
library(loomR)

thymus1.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2020_Parekh_thymus/thymus1/thymus1-c7/outs/filtered_feature_bc_matrix")
thymus1 <- CreateSeuratObject(counts = thymus1.data, project = "CD34")
thymus1$stim <- 'thymus1'
thymus1[["percent.mt"]] <- PercentageFeatureSet(thymus1, pattern = "^MT-")
VlnPlot(thymus1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
thymus1 <- subset(thymus1, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(thymus1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
thymus1 <- NormalizeData(thymus1, normalization.method = "LogNormalize", scale.factor = 10000)
thymus1 <- FindVariableFeatures(thymus1, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(thymus1)
thymus1 <- ScaleData(thymus1, features = all.genes)
thymus1 <- RunPCA(thymus1, features = VariableFeatures(object = thymus1))
VizDimLoadings(thymus1, dims = 1:5, reduction = "pca")
ElbowPlot(thymus1)
DimPlot(thymus1, reduction = "pca")
thymus1 <- RunUMAP(thymus1, reduction = "pca", dims = 1:13)
thymus1 <- FindNeighbors(thymus1, reduction = "pca", dims = 1:13)
thymus1 <- FindClusters(thymus1, resolution = 0.8)
DimPlot(thymus1, reduction = "umap", label = T)
VlnPlot(thymus1, features = "CD34")
FeaturePlot(thymus1, features = c("CD34"), max.cutoff = 3, cols = c("grey", "red"))
table(Idents(subset(thymus1, subset= CD34>0)))/table(Idents(thymus1))
markers.to.plot <- c("CD34","BAALC","MEIS1","MEF2C","IRF8","CEBPA","MPO","DTX1","TYROBP","CD7","CD2","CD1A","CD44","BCL11B","PTCRA","RAG1","CD8A","CD8B","PTPRC")
plot3=DotPlot(thymus1, features = markers.to.plot, cols = c("grey", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave("Dotplot.pdf", plot=plot3, width=12, height=6)
saveRDS(thymus1, file = "/home/yushiya/data/cd34/data/2020_Parekh_thymus/thymus1-c7.rds")
thymus1 <- readRDS(file = "/home/yushiya/data/cd34/data/2020_Parekh_thymus/thymus1-c7.rds")
thymus1.loom <- as.loom(thymus1, filename = "/home/yushiya/data/cd34/data/2020_Parekh_thymus/thymus1-c7.loom", verbose = FALSE)
thymus1.loom$close_all()

thymus2.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2020_Parekh_thymus/thymus2/thymus2-c7/outs/filtered_feature_bc_matrix")
thymus2 <- CreateSeuratObject(counts = thymus2.data, project = "CD34")
thymus2$stim <- 'thymus2'
thymus2[["percent.mt"]] <- PercentageFeatureSet(thymus2, pattern = "^MT-")
VlnPlot(thymus2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
thymus2 <- subset(thymus2, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(thymus2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
thymus2 <- NormalizeData(thymus2, normalization.method = "LogNormalize", scale.factor = 10000)
thymus2 <- FindVariableFeatures(thymus2, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(thymus2)
thymus2 <- ScaleData(thymus2, features = all.genes)
thymus2 <- RunPCA(thymus2, features = VariableFeatures(object = thymus2))
VizDimLoadings(thymus2, dims = 1:5, reduction = "pca")
ElbowPlot(thymus2)
DimPlot(thymus2, reduction = "pca")
thymus2 <- RunUMAP(thymus2, reduction = "pca", dims = 1:13)
thymus2 <- FindNeighbors(thymus2, reduction = "pca", dims = 1:13)
thymus2 <- FindClusters(thymus2, resolution = 0.8)
DimPlot(thymus2, reduction = "umap", label=T)
VlnPlot(thymus2, features = "CD34")
FeaturePlot(thymus2, features = c("CD34"), max.cutoff = 3, cols = c("grey", "red"))
table(Idents(subset(thymus2, subset= CD34>0)))/table(Idents(thymus2))
markers.to.plot <- c("CD34","BAALC","MEIS1","MEF2C","IRF8","CEBPA","MPO","DTX1","TYROBP","CD7","CD2","CD1A","CD44","BCL11B","PTCRA","RAG1","CD8A","CD8B","PTPRC")
plot3=DotPlot(thymus2, features = markers.to.plot, cols = c("grey", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave("Dotplot.pdf", plot=plot3, width=12, height=6)
saveRDS(thymus2, file = "/home/yushiya/data/cd34/data/2020_Parekh_thymus/thymus2-c7.rds")
thymus2 <- readRDS(file = "/home/yushiya/data/cd34/data/2020_Parekh_thymus/thymus2-c7.rds")
thymus2.loom <- as.loom(thymus2, filename = "/home/yushiya/data/cd34/data/2020_Parekh_thymus/thymus2-c7.loom", verbose = FALSE)
thymus2.loom$close_all()

thymus3.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2020_Parekh_thymus/thymus3/thymus3-c7/outs/filtered_feature_bc_matrix")
thymus3 <- CreateSeuratObject(counts = thymus3.data, project = "CD34")
thymus3$stim <- 'thymus3'
thymus3[["percent.mt"]] <- PercentageFeatureSet(thymus3, pattern = "^MT-")
VlnPlot(thymus3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
thymus3 <- subset(thymus3, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(thymus3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
thymus3 <- NormalizeData(thymus3, normalization.method = "LogNormalize", scale.factor = 10000)
thymus3 <- FindVariableFeatures(thymus3, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(thymus3)
thymus3 <- ScaleData(thymus3, features = all.genes)
thymus3 <- RunPCA(thymus3, features = VariableFeatures(object = thymus3))
VizDimLoadings(thymus3, dims = 1:5, reduction = "pca")
ElbowPlot(thymus3)
DimPlot(thymus3, reduction = "pca", label=T)
thymus3 <- RunUMAP(thymus3, reduction = "pca", dims = 1:13)
thymus3 <- FindNeighbors(thymus3, reduction = "pca", dims = 1:13)
thymus3<- FindClusters(thymus3, resolution = 0.8)
DimPlot(thymus3, reduction = "umap",label=T)
FeaturePlot(thymus3, features = c("CD34"), max.cutoff = 3, cols = c("grey", "red"))
VlnPlot(thymus3, features = "CD34")
table(Idents(subset(thymus3, subset= CD34>0)))/table(Idents(thymus3))
markers.to.plot <- c("CD34","BAALC","MEIS1","MEF2C","IRF8","CEBPA","MPO","DTX1","TYROBP","CD7","CD2","CD1A","CD44","BCL11B","PTCRA","RAG1","CD8A","CD8B","PTPRC")
plot3=DotPlot(thymus3, features = markers.to.plot, cols = c("grey", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave("Dotplot.pdf", plot=plot3, width=12, height=6)
saveRDS(thymus3, file = "/home/yushiya/data/cd34/data/2020_Parekh_thymus/thymus3-c7.rds")
thymus3 <- readRDS(file = "/home/yushiya/data/cd34/data/2020_Parekh_thymus/thymus3-c7.rds")
thymus3.loom <- as.loom(thymus3, filename = "/home/yushiya/data/cd34/data/2020_Parekh_thymus/thymus3-c7.loom", verbose = FALSE)
thymus3.loom$close_all()

thymus.anchors <- FindIntegrationAnchors(object.list = c(thymus1,thymus2,thymus3), dims = 1:20)
thymus.combined <- IntegrateData(anchorset = thymus.anchors, dims = 1:20)
DefaultAssay(thymus.combined) <- "integrated"
thymus.combined@meta.data$ori_stim <- thymus.combined@meta.data$stim
thymus.combined@meta.data$stim <- "thymus"
# Run the standard workflow for visualization and clustering
thymus.combined <- ScaleData(thymus.combined, verbose = FALSE)
thymus.combined <- RunPCA(thymus.combined, npcs = 30, verbose = FALSE)
ElbowPlot(thymus.combined)
#Clustering
thymus.combined <- FindNeighbors(thymus.combined, reduction = "pca", dims = 1:30)
thymus.combined <- FindClusters(thymus.combined, resolution = 0.8)
# umap
thymus.combined <- RunUMAP(thymus.combined, reduction = "pca", dims = 1:30)
thymus.combined <- CellCycleScoring(thymus.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
RidgePlot(thymus.combined, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
# Visualization
DimPlot(thymus.combined, reduction = "umap", group.by = "ori_stim")
DimPlot(thymus.combined, reduction = "umap", split.by = "ori_stim")
DimPlot(thymus.combined, reduction = "umap", label = TRUE)
DimPlot(thymus.combined, reduction = "umap", group.by = "Phase")
#set cell identity
Idents(thymus.combined) <- 'seurat_clusters'
DimPlot(thymus.combined, cells.highlight = colnames(subset(thymus.combined,subset = seurat_clusters=="0")), cols.highlight = "blue")
table(thymus.combined@meta.data[["stim"]],thymus.combined@meta.data[["seurat_clusters"]])
#Find cell markers
DefaultAssay(thymus.combined) <- "RNA"
thymus.combined.markers <- FindAllMarkers(thymus.combined, assay="RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top20 <- thymus.combined.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
plot1=DoHeatmap(thymus.combined, features = top20$gene) + NoLegend()
ggsave("top10_markers.pdf", plot=plot1, width=15, height=25)
write.csv(thymus.combined.markers, "top20_diff_genes.csv", row.names = F)
FeaturePlot(thymus.combined, features = c("CD34"), max.cutoff = 3, cols = c("grey", "red"))
markers.to.plot <- c("CD34","CD7","CD2","CD1A","CD44","MEF2C","MME","BCL11A","CD3","CD4","CD8A","CD8B","PTPRC","CD3G","TYROBP","CD19","CD27","GZMB","HOXA9","FXYD2","SH3TC1","CCR9","VPREB1","MS4A1","IL3RA")
markers.to.plot <- c("CD34","BAALC","MEIS1","MEF2C","IRF8","CEBPA","MPO","DTX1","TYROBP","CD7","CD2","CD1A","CD44","BCL11B","PTCRA","RAG1","CD8A","CD8B","PTPRC")
plot3=DotPlot(thymus.combined, features = markers.to.plot, cols = c("grey", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave("Dotplot.pdf", plot=plot3, width=13, height=6)
VlnPlot(thymus.combined, features = "CD34")
table(Idents(subset(thymus.combined, subset= CD34>0)))/table(Idents(thymus.combined))
VlnPlot(thymus.combined, features = "CD34", pt.size = 0, combine = FALSE)
thymus.combined <- FindVariableFeatures(thymus.combined, selection.method = "vst", nfeatures = 2000)
thymus.combined <- subset(thymus.combined, subset= seurat_clusters=='0' | seurat_clusters=='1' | seurat_clusters=='2' | seurat_clusters=='3' | seurat_clusters=='4' | seurat_clusters=='5' | seurat_clusters=='7')
saveRDS(thymus.combined, file = "/home/yushiya/data/cd34/data/2020_Parekh_thymus/thymus_combined-c7.rds")
thymus.combined <- readRDS(file = "/home/yushiya/data/cd34/data/2020_Parekh_thymus/thymus_combined-c7.rds")

PB_LD1.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2020_bioRxiv_BMPBSP/LD1/LD1-PB/outs/filtered_feature_bc_matrix")
PB_LD1 <- CreateSeuratObject(counts = PB_LD1.data, project = "CD34")
PB_LD1$stim <- 'PB_LD1'
PB_LD1[["percent.mt"]] <- PercentageFeatureSet(PB_LD1, pattern = "^MT-")
VlnPlot(PB_LD1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
PB_LD1 <- subset(PB_LD1, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(PB_LD1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
PB_LD1 <- NormalizeData(PB_LD1, normalization.method = "LogNormalize", scale.factor = 10000)
PB_LD1 <- FindVariableFeatures(PB_LD1, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(PB_LD1)
PB_LD1 <- ScaleData(PB_LD1, features = all.genes)
PB_LD1 <- RunPCA(PB_LD1, features = VariableFeatures(object = PB_LD1))
VizDimLoadings(PB_LD1, dims = 1:5, reduction = "pca")
ElbowPlot(PB_LD1)
DimPlot(PB_LD1, reduction = "pca")
PB_LD1 <- RunUMAP(PB_LD1, reduction = "pca", dims = 1:17)
PB_LD1 <- FindNeighbors(PB_LD1, reduction = "pca", dims = 1:17)
PB_LD1 <- FindClusters(PB_LD1, resolution = 0.8)
DimPlot(PB_LD1, reduction = "umap",label=T)
FeaturePlot(PB_LD1, features = c("IGHM","GATA2"), max.cutoff = 3, cols = c("grey", "red"))
VlnPlot(PB_LD1, features = "CD34")
table(Idents(subset(PB_LD1, subset= CD34>0)))/table(Idents(PB_OD1))
markers.to.plot <- c("AVP","CRHBP","MEG3","EMCN","THY1","CD34","CD38","PTPRC","FLT3","KIT","GATA1","GATA2","ITGA2B","KLF1","TFRC","HBD","CD36","APOE","CA1","AHSP","GYPA","HEMGN","SELP","GP1BA","GP9","PF4","TPSAB1","HDC","MS4A2","EPX","CSF3R","CSF1R","CD33","MPO",
                     "ELANE","AZU1","PRTN3","LYZ","IRF8","IL3RA","SPIB","CD68","CD14","MME","IL7R","VPREB1","VPREB3","CD79A","CD79B","IGHM","IGHD","FCER2")
plot3=DotPlot(PB_LD1, features = markers.to.plot, cols = c("grey", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave("Dotplot.pdf", plot=plot3, width=13, height=6)
saveRDS(PB_LD1, file = "/home/yushiya/data/cd34/data/2020_bioRxiv_BMPBSP/PB_LD1.rds")
PB_LD1 <- readRDS(file = "/home/yushiya/data/cd34/data/2020_bioRxiv_BMPBSP/PB_LD1.rds")
PB_LD1.loom <- as.loom(PB_LD1, filename = "/home/yushiya/data/cd34/data/2020_bioRxiv_BMPBSP/PB_LD1.loom", verbose = FALSE)
PB_LD1.loom$close_all()

PB_LD2.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2020_bioRxiv_BMPBSP/LD2/LD2-PB/outs/filtered_feature_bc_matrix")
PB_LD2 <- CreateSeuratObject(counts = PB_LD2.data, project = "CD34")
PB_LD2$stim <- 'PB_LD2'
PB_LD2[["percent.mt"]] <- PercentageFeatureSet(PB_LD2, pattern = "^MT-")
VlnPlot(PB_LD2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
PB_LD2 <- subset(PB_LD2, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(PB_LD2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
PB_LD2 <- NormalizeData(PB_LD2, normalization.method = "LogNormalize", scale.factor = 10000)
PB_LD2 <- FindVariableFeatures(PB_LD2, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(PB_LD2)
PB_LD2 <- ScaleData(PB_LD2, features = all.genes)
PB_LD2 <- RunPCA(PB_LD2, features = VariableFeatures(object = PB_LD2))
VizDimLoadings(PB_LD2, dims = 1:5, reduction = "pca")
ElbowPlot(PB_LD2)
DimPlot(PB_LD2, reduction = "pca")
PB_LD2 <- RunUMAP(PB_LD2, reduction = "pca", dims = 1:17)
PB_LD2 <- FindNeighbors(PB_LD2, reduction = "pca", dims = 1:17)
PB_LD2 <- FindClusters(PB_LD2, resolution = 0.8)
DimPlot(PB_LD2, reduction = "umap",label=T)
FeaturePlot(PB_LD2, features = c("IGHM","GATA2"), max.cutoff = 3, cols = c("grey", "red"))
VlnPlot(PB_LD2, features = "CD34")
table(Idents(subset(PB_LD2, subset= CD34>0)))/table(Idents(PB_LD2))
markers.to.plot <- c("AVP","CRHBP","MEG3","EMCN","THY1","CD34","CD38","PTPRC","FLT3","KIT","GATA1","GATA2","ITGA2B","KLF1","TFRC","HBD","CD36","APOE","CA1","AHSP","GYPA","HEMGN","SELP","GP1BA","GP9","PF4","TPSAB1","HDC","MS4A2","EPX","CSF3R","CSF1R","CD33","MPO",
                     "ELANE","AZU1","PRTN3","LYZ","IRF8","IL3RA","SPIB","CD68","CD14","MME","IL7R","VPREB1","VPREB3","CD79A","CD79B","IGHM","IGHD","FCER2")
plot3=DotPlot(PB_LD2, features = markers.to.plot, cols = c("grey", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave("Dotplot.pdf", plot=plot3, width=13, height=6)
saveRDS(PB_LD2, file = "/home/yushiya/data/cd34/data/2020_bioRxiv_BMPBSP/PB_LD2.rds")
PB_LD2 <- readRDS(file = "/home/yushiya/data/cd34/data/2020_bioRxiv_BMPBSP/PB_LD2.rds")
PB_LD2.loom <- as.loom(PB_LD2, filename = "/home/yushiya/data/cd34/data/2020_bioRxiv_BMPBSP/PB_LD2.loom", verbose = FALSE)
PB_LD2.loom$close_all()

PB_LD3.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2020_bioRxiv_BMPBSP/LD3/LD3-PB/outs/filtered_feature_bc_matrix")
PB_LD3 <- CreateSeuratObject(counts = PB_LD3.data, project = "CD34")
PB_LD3$stim <- 'PB_LD3'
PB_LD3[["percent.mt"]] <- PercentageFeatureSet(PB_LD3, pattern = "^MT-")
VlnPlot(PB_LD3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
PB_LD3 <- subset(PB_LD3, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(PB_LD3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
PB_LD3 <- NormalizeData(PB_LD3, normalization.method = "LogNormalize", scale.factor = 10000)
PB_LD3 <- FindVariableFeatures(PB_LD3, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(PB_LD3)
PB_LD3 <- ScaleData(PB_LD3, features = all.genes)
PB_LD3 <- RunPCA(PB_LD3, features = VariableFeatures(object = PB_LD3))
VizDimLoadings(PB_LD3, dims = 1:5, reduction = "pca")
ElbowPlot(PB_LD3)
DimPlot(PB_LD3, reduction = "pca")
PB_LD3 <- RunUMAP(PB_LD3, reduction = "pca", dims = 1:17)
PB_LD3 <- FindNeighbors(PB_LD3, reduction = "pca", dims = 1:17)
PB_LD3 <- FindClusters(PB_LD3, resolution = 0.8)
DimPlot(PB_LD3, reduction = "umap",label=T)
FeaturePlot(PB_LD3, features = c("IGHM","GATA2"), max.cutoff = 3, cols = c("grey", "red"))
VlnPlot(PB_LD3, features = "CD34")
table(Idents(subset(PB_LD3, subset= CD34>0)))/table(Idents(PB_LD2))
markers.to.plot <- c("AVP","CRHBP","MEG3","EMCN","THY1","CD34","CD38","PTPRC","FLT3","KIT","GATA1","GATA2","ITGA2B","KLF1","TFRC","HBD","CD36","APOE","CA1","AHSP","GYPA","HEMGN","SELP","GP1BA","GP9","PF4","TPSAB1","HDC","MS4A2","EPX","CSF3R","CSF1R","CD33","MPO",
                     "ELANE","AZU1","PRTN3","LYZ","IRF8","IL3RA","SPIB","CD68","CD14","MME","IL7R","VPREB1","VPREB3","CD79A","CD79B","IGHM","IGHD","FCER2")
plot3=DotPlot(PB_LD3, features = markers.to.plot, cols = c("grey", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave("Dotplot.pdf", plot=plot3, width=13, height=6)
saveRDS(PB_LD3, file = "/home/yushiya/data/cd34/data/2020_bioRxiv_BMPBSP/PB_LD3.rds")
PB_LD3 <- readRDS(file = "/home/yushiya/data/cd34/data/2020_bioRxiv_BMPBSP/PB_LD3.rds")
PB_LD3.loom <- as.loom(PB_LD3, filename = "/home/yushiya/data/cd34/data/2020_bioRxiv_BMPBSP/PB_LD3.loom", verbose = FALSE)
PB_LD3.loom$close_all()

PB_LD4.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2020_bioRxiv_BMPBSP/LD4/LD4-PB/outs/filtered_feature_bc_matrix")
PB_LD4 <- CreateSeuratObject(counts = PB_LD4.data, project = "CD34")
PB_LD4$stim <- 'PB_LD4'
PB_LD4[["percent.mt"]] <- PercentageFeatureSet(PB_LD4, pattern = "^MT-")
VlnPlot(PB_LD4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
PB_LD4 <- subset(PB_LD4, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(PB_LD4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
PB_LD4 <- NormalizeData(PB_LD4, normalization.method = "LogNormalize", scale.factor = 10000)
PB_LD4 <- FindVariableFeatures(PB_LD4, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(PB_LD4)
PB_LD4 <- ScaleData(PB_LD4, features = all.genes)
PB_LD4 <- RunPCA(PB_LD4, features = VariableFeatures(object = PB_LD4))
VizDimLoadings(PB_LD4, dims = 1:5, reduction = "pca")
ElbowPlot(PB_LD4)
DimPlot(PB_LD4, reduction = "pca")
PB_LD4 <- RunUMAP(PB_LD4, reduction = "pca", dims = 1:17)
PB_LD4 <- FindNeighbors(PB_LD4, reduction = "pca", dims = 1:17)
PB_LD4 <- FindClusters(PB_LD4, resolution = 0.8)
DimPlot(PB_LD4, reduction = "umap",label=T)
FeaturePlot(PB_LD4, features = c("IGHM","GATA2"), max.cutoff = 3, cols = c("grey", "red"))
VlnPlot(PB_LD4, features = "CD34")
table(Idents(subset(PB_LD4, subset= CD34>0)))/table(Idents(PB_LD2))
markers.to.plot <- c("AVP","CRHBP","MEG3","EMCN","THY1","CD34","CD38","PTPRC","FLT3","KIT","GATA1","GATA2","ITGA2B","KLF1","TFRC","HBD","CD36","APOE","CA1","AHSP","GYPA","HEMGN","SELP","GP1BA","GP9","PF4","TPSAB1","HDC","MS4A2","EPX","CSF3R","CSF1R","CD33","MPO",
                     "ELANE","AZU1","PRTN3","LYZ","IRF8","IL3RA","SPIB","CD68","CD14","MME","IL7R","VPREB1","VPREB3","CD79A","CD79B","IGHM","IGHD","FCER2")
plot3=DotPlot(PB_LD4, features = markers.to.plot, cols = c("grey", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave("Dotplot.pdf", plot=plot3, width=13, height=6)
saveRDS(PB_LD4, file = "/home/yushiya/data/cd34/data/2020_bioRxiv_BMPBSP/PB_LD4.rds")
PB_LD4 <- readRDS(file = "/home/yushiya/data/cd34/data/2020_bioRxiv_BMPBSP/PB_LD4.rds")
PB_LD4.loom <- as.loom(PB_LD4, filename = "/home/yushiya/data/cd34/data/2020_bioRxiv_BMPBSP/PB_LD4.loom", verbose = FALSE)
PB_LD4.loom$close_all()

PB.anchors <- FindIntegrationAnchors(object.list = c(PB_LD1,PB_LD2,PB_LD3,PB_LD4), dims = 1:30)
PB.combined <- IntegrateData(anchorset = PB.anchors, dims = 1:30)
DefaultAssay(PB.combined) <- "integrated"
PB.combined@meta.data$ori_stim <- PB.combined@meta.data$stim
PB.combined@meta.data$stim <- "PB"
# Run the standard workflow for visualization and clustering
PB.combined <- ScaleData(PB.combined, verbose = FALSE)
PB.combined <- RunPCA(PB.combined, npcs = 30, verbose = FALSE)
ElbowPlot(PB.combined)
#Clustering
PB.combined <- FindNeighbors(PB.combined, reduction = "pca", dims = 1:25)
PB.combined <- FindClusters(PB.combined, resolution = 0.8)
# umap
PB.combined <- RunUMAP(PB.combined, reduction = "pca", dims = 1:25)
# Visualization
DimPlot(PB.combined, reduction = "umap", group.by = "stim")
DimPlot(PB.combined, reduction = "umap", split.by = "ori_stim", ncol=2)
DimPlot(PB.combined, reduction = "umap", label = TRUE)
DimPlot(PB.combined, reduction = "umap", group.by = "cd34_clus", label=TRUE)
table(PB.combined@meta.data[["stim"]],PB.combined@meta.data[["seurat_clusters"]])
DefaultAssay(PB.combined) <- "RNA"
PB.combined.markers <- FindAllMarkers(PB.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top20 <- PB.combined.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
plot1=DoHeatmap(PB.combined, features = top10$gene) + NoLegend()
ggsave("top10_markers.pdf", plot=plot1, width=11, height=9)
write.csv(PB.combined.markers, "top20_diff_genes.csv", row.names = F)
FeaturePlot(PB.combined, features = c("CD34"), max.cutoff = 3, cols = c("grey", "red"))
FeaturePlot(PB.combined, features = c("PLEK","HBB","MPO","SPIB","CD79A","DNTT"), max.cutoff = 3, cols = c("grey", "blue"))
VlnPlot(PB.combined, features = "CD34")
table(Idents(subset(PB.combined, subset= CD34>0)))/table(Idents(PB.combined))
markers.to.plot <- c("AVP","CRHBP","MEG3","EMCN","THY1","CD34","CD38","PTPRC","FLT3","KIT","GATA1","GATA2","ITGA2B","KLF1","TFRC","HBD","CD36","APOE","CA1","AHSP","GYPA","HEMGN","SELP","GP1BA","GP9","PF4","TPSAB1","HDC","MS4A2","EPX","CSF3R","CSF1R","CD33","MPO",
                     "ELANE","AZU1","PRTN3","LYZ","IRF8","IL3RA","SPIB","CD68","CD14","MME","IL7R","VPREB1","VPREB3","RAG2","EBF1","PAX5","CD24","CD19","MS4A1","CD79A","CD79B","IGHM","IGHD","FCER2","IGHA1","IGHG1","SDC1","TNFRSF17","CD7","CD8A","CD8B","CD3D",
                     "CD3E","CD3G","CD247","KLRB1","NCR3","KLRD1","GZMH")
markers.to.plot <- c("AVP","CRHBP","MEG3","EMCN","THY1","CD34","CD38","PTPRC","FLT3","KIT","GATA1","GATA2","ITGA2B","KLF1","TFRC","HBD","CD36","APOE","CA1","AHSP","GYPA","HEMGN","SELP","GP1BA","GP9","PF4","TPSAB1","HDC","MS4A2","EPX","CSF3R","CSF1R","CD33","MPO",
                     "ELANE","AZU1","PRTN3","LYZ","IRF8","IL3RA","SPIB","CD68","CD14","MME","IL7R","VPREB1","VPREB3","CD79A","CD79B","IGHM","IGHD","FCER2")
plot3=DotPlot(PB.combined, features = markers.to.plot, cols = c("grey", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave("Dotplot.pdf", plot=plot3, width=16, height=6)
PB.combined <- subset(PB.combined, subset= seurat_clusters=="0" | seurat_clusters=="1" | seurat_clusters=="2" | seurat_clusters=="3" | seurat_clusters=="4" | seurat_clusters=="5" | seurat_clusters=="6" | seurat_clusters=="7" | seurat_clusters=="8"|seurat_clusters=="9"|seurat_clusters=="10"|seurat_clusters=="11"|seurat_clusters=="12"|seurat_clusters=="13"|seurat_clusters=="14"|seurat_clusters=="15"|seurat_clusters=="16" | seurat_clusters=="17" | seurat_clusters=="18"|seurat_clusters=="19"|seurat_clusters=="20"|seurat_clusters=="21"|seurat_clusters=="23"|seurat_clusters=="24")
saveRDS(PB.combined, file = "/home/yushiya/data/cd34/data/2020_bioRxiv_BMPBSP/PB.combined_1234.rds")
PB.combined <- readRDS(file = "/home/yushiya/data/cd34/data/2020_bioRxiv_BMPBSP/PB.combined_1234.rds")

BM1.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2019_Blood_BM/BM1/BM1-c7/outs/filtered_feature_bc_matrix")
BM1 <- CreateSeuratObject(counts = BM1.data, project = "CD34")
BM1$stim <- 'BM1'
BM1[["percent.mt"]] <- PercentageFeatureSet(BM1, pattern = "^MT-")
VlnPlot(BM1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
BM1 <- subset(BM1, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(BM1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
BM1 <- NormalizeData(BM1, normalization.method = "LogNormalize", scale.factor = 10000)
BM1 <- FindVariableFeatures(BM1, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(BM1)
BM1 <- ScaleData(BM1, features = all.genes)
BM1 <- RunPCA(BM1, features = VariableFeatures(object = BM1))
VizDimLoadings(BM1, dims = 1:5, reduction = "pca")
ElbowPlot(BM1)
DimPlot(BM1, reduction = "pca")
BM1 <- RunUMAP(BM1, reduction = "pca", dims = 1:17)
BM1 <- FindNeighbors(BM1, reduction = "pca", dims = 1:17)
BM1 <- FindClusters(BM1, resolution = 0.8)
DimPlot(BM1, reduction = "umap",label=T)
FeaturePlot(BM1, features = c("CD34"), max.cutoff = 3, cols = c("grey", "red"))
VlnPlot(BM1, features = "CD34")
table(Idents(subset(BM1, subset= CD34>0)))/table(Idents(BM1))
markers.to.plot <- c("AVP","CRHBP","MEG3","EMCN","THY1","CD34","CD38","PTPRC","FLT3","KIT","GATA1","GATA2","ITGA2B","KLF1","TFRC","HBD","CD36","APOE","CA1","AHSP","GYPA","HEMGN","SELP","GP1BA","GP9","PF4","TPSAB1","HDC","MS4A2","EPX","CSF3R","CSF1R","CD33","MPO",
                     "ELANE","AZU1","PRTN3","LYZ","IRF8","IL3RA","SPIB","CD68","CD14","MME","IL7R","VPREB1","VPREB3","RAG2","EBF1","PAX5","CD24","CD19","MS4A1","CD79A","CD79B","IGHM","IGHD")
plot3=DotPlot(BM1, features = markers.to.plot, cols = c("grey", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave("Dotplot.pdf", plot=plot3, width=13, height=6)
#BM1 <- subset(BM1, subset= seurat_clusters=="0" | seurat_clusters=="1" | seurat_clusters=="2" | seurat_clusters=="3" | seurat_clusters=="4" | seurat_clusters=="5" | seurat_clusters=="6" | seurat_clusters=="7" | seurat_clusters=="8"|seurat_clusters=="9"|seurat_clusters=="11"|seurat_clusters=="13"|seurat_clusters=="15")
saveRDS(BM1, file = "/home/yushiya/data/cd34/data/2019_Blood_BM/BM1_c7.rds")
BM1 <- readRDS(file = "/home/yushiya/data/cd34/data/2019_Blood_BM/BM1_c7.rds")
BM1.loom <- as.loom(BM1, filename = "/home/yushiya/data/cd34/data/2019_Blood_BM/BM1_c7.loom", verbose = FALSE)
BM1.loom$close_all()

BM2.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2019_Blood_BM/BM2/BM2-c7/outs/filtered_feature_bc_matrix")
BM2 <- CreateSeuratObject(counts = BM2.data, project = "CD34")
BM2$stim <- 'BM2'
BM2[["percent.mt"]] <- PercentageFeatureSet(BM2, pattern = "^MT-")
VlnPlot(BM2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
BM2 <- subset(BM2, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(BM2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
BM2 <- NormalizeData(BM2, normalization.method = "LogNormalize", scale.factor = 10000)
BM2 <- FindVariableFeatures(BM2, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(BM2)
BM2 <- ScaleData(BM2, features = all.genes)
BM2 <- RunPCA(BM2, features = VariableFeatures(object = BM2))
VizDimLoadings(BM2, dims = 1:5, reduction = "pca")
ElbowPlot(BM2)
DimPlot(BM2, reduction = "pca")
BM2 <- RunUMAP(BM2, reduction = "pca", dims = 1:17)
BM2 <- FindNeighbors(BM2, reduction = "pca", dims = 1:17)
BM2 <- FindClusters(BM2, resolution = 0.8)
DimPlot(BM2, reduction = "umap", label=T)
FeaturePlot(BM2, features = c("CD34"), max.cutoff = 3, cols = c("grey", "red"))
VlnPlot(BM2, features = "CD34")
table(Idents(subset(BM2, subset= CD34>0)))/table(Idents(BM2))
markers.to.plot <- c("AVP","CRHBP","MEG3","EMCN","THY1","CD34","CD38","PTPRC","FLT3","KIT","GATA1","GATA2","ITGA2B","KLF1","TFRC","HBD","CD36","APOE","CA1","AHSP","GYPA","HEMGN","SELP","GP1BA","GP9","PF4","TPSAB1","HDC","MS4A2","EPX","CSF3R","CSF1R","CD33","MPO",
                     "ELANE","AZU1","PRTN3","LYZ","IRF8","IL3RA","SPIB","CD68","CD14","MME","IL7R","VPREB1","VPREB3","RAG2","EBF1","PAX5","CD24","CD19","MS4A1","CD79A","CD79B","IGHM","IGHD")
plot3=DotPlot(BM2, features = markers.to.plot, cols = c("grey", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave("Dotplot.pdf", plot=plot3, width=13, height=6)
#BM2 <- subset(BM2, subset= seurat_clusters=="0" | seurat_clusters=="1" | seurat_clusters=="2" | seurat_clusters=="3" | seurat_clusters=="4" | seurat_clusters=="5" | seurat_clusters=="6" | seurat_clusters=="7" | seurat_clusters=="8"|seurat_clusters=="9"|seurat_clusters=="12"|seurat_clusters=="14")
saveRDS(BM2, file = "/home/yushiya/data/cd34/data/2019_Blood_BM/BM2_c7.rds")
BM2 <- readRDS(file = "/home/yushiya/data/cd34/data/2019_Blood_BM/BM2_c7.rds")
BM2.loom <- as.loom(BM2, filename = "/home/yushiya/data/cd34/data/2019_Blood_BM/BM2_c7.loom", verbose = FALSE)
BM2.loom$close_all()

BM3.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2019_Blood_BM/BM3/BM3-c7/outs/filtered_feature_bc_matrix")
BM3 <- CreateSeuratObject(counts = BM3.data, project = "CD34")
BM3$stim <- 'BM3'
BM3[["percent.mt"]] <- PercentageFeatureSet(BM3, pattern = "^MT-")
VlnPlot(BM3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
BM3 <- subset(BM3, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(BM3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
BM3 <- NormalizeData(BM3, normalization.method = "LogNormalize", scale.factor = 10000)
BM3 <- FindVariableFeatures(BM3, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(BM3)
BM3 <- ScaleData(BM3, features = all.genes)
BM3 <- RunPCA(BM3, features = VariableFeatures(object = BM3))
VizDimLoadings(BM3, dims = 1:5, reduction = "pca")
ElbowPlot(BM3)
DimPlot(BM3, reduction = "pca")
BM3 <- RunUMAP(BM3, reduction = "pca", dims = 1:17)
BM3 <- FindNeighbors(BM3, reduction = "pca", dims = 1:17)
BM3 <- FindClusters(BM3, resolution = 0.9)
DimPlot(BM3, reduction = "umap", label=T)
FeaturePlot(BM3, features = c("CD34"), max.cutoff = 3, cols = c("grey", "red"))
VlnPlot(BM3, features = "CD34")
table(Idents(subset(BM3, subset= CD34>0)))/table(Idents(BM3))
markers.to.plot <- c("AVP","CRHBP","MEG3","EMCN","THY1","CD34","CD38","PTPRC","FLT3","KIT","GATA1","GATA2","ITGA2B","KLF1","TFRC","HBD","CD36","APOE","CA1","AHSP","GYPA","HEMGN","SELP","GP1BA","GP9","PF4","TPSAB1","HDC","MS4A2","EPX","CSF3R","CSF1R","CD33","MPO",
                     "ELANE","AZU1","PRTN3","LYZ","IRF8","IL3RA","SPIB","CD68","CD14","MME","IL7R","VPREB1","VPREB3","RAG2","EBF1","PAX5","CD24","CD19","MS4A1","CD79A","CD79B","IGHM","IGHD")
plot3=DotPlot(BM3, features = markers.to.plot, cols = c("grey", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave("Dotplot.pdf", plot=plot3, width=13, height=6)
#BM3 <- subset(BM3, subset= seurat_clusters=="0" | seurat_clusters=="1" | seurat_clusters=="2" | seurat_clusters=="3" | seurat_clusters=="4" | seurat_clusters=="5" | seurat_clusters=="6" | seurat_clusters=="7" | seurat_clusters=="8"|seurat_clusters=="9")
saveRDS(BM3, file = "/home/yushiya/data/cd34/data/2019_Blood_BM/BM3_c7.rds")
BM3 <- readRDS(file = "/home/yushiya/data/cd34/data/2019_Blood_BM/BM3_c7.rds")
BM3.loom <- as.loom(BM3, filename = "/home/yushiya/data/cd34/data/2019_Blood_BM/BM3_c7.loom", verbose = FALSE)
BM3.loom$close_all()

BM4.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2019_Blood_BM/BM4/BM4-c7/outs/filtered_feature_bc_matrix")
BM4 <- CreateSeuratObject(counts = BM4.data, project = "CD34")
BM4$stim <- 'BM4'
BM4[["percent.mt"]] <- PercentageFeatureSet(BM4, pattern = "^MT-")
VlnPlot(BM4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
BM4 <- subset(BM4, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(BM4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
BM4 <- NormalizeData(BM4, normalization.method = "LogNormalize", scale.factor = 10000)
BM4 <- FindVariableFeatures(BM4, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(BM4)
BM4 <- ScaleData(BM4, features = all.genes)
BM4 <- RunPCA(BM4, features = VariableFeatures(object = BM4))
VizDimLoadings(BM4, dims = 1:5, reduction = "pca")
ElbowPlot(BM4)
DimPlot(BM4, reduction = "pca")
BM4 <- RunUMAP(BM4, reduction = "pca", dims = 1:17)
BM4 <- FindNeighbors(BM4, reduction = "pca", dims = 1:17)
BM4 <- FindClusters(BM4, resolution = 0.9)
DimPlot(BM4, reduction = "umap", label=T)
FeaturePlot(BM4, features = c("CD34"), max.cutoff = 3, cols = c("grey", "red"))
VlnPlot(BM4, features = "CD34")
table(Idents(subset(BM4, subset= CD34>0)))/table(Idents(BM4))
markers.to.plot <- c("AVP","CRHBP","MEG3","EMCN","THY1","CD34","CD38","PTPRC","FLT3","KIT","GATA1","GATA2","ITGA2B","KLF1","TFRC","HBD","CD36","APOE","CA1","AHSP","GYPA","HEMGN","SELP","GP1BA","GP9","PF4","TPSAB1","HDC","MS4A2","EPX","CSF3R","CSF1R","CD33","MPO",
                     "ELANE","AZU1","PRTN3","LYZ","IRF8","IL3RA","SPIB","CD68","CD14","MME","IL7R","VPREB1","VPREB3","RAG2","EBF1","PAX5","CD24","CD19","MS4A1","CD79A","CD79B","IGHM","IGHD")
plot3=DotPlot(BM4, features = markers.to.plot, cols = c("grey", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave("Dotplot.pdf", plot=plot3, width=13, height=6)
#BM4 <- subset(BM4, subset= seurat_clusters=="0" | seurat_clusters=="1" | seurat_clusters=="2" | seurat_clusters=="3" | seurat_clusters=="4" | seurat_clusters=="5" | seurat_clusters=="6" | seurat_clusters=="7" | seurat_clusters=="8"|seurat_clusters=="9")
saveRDS(BM4, file = "/home/yushiya/data/cd34/data/2019_Blood_BM/BM4_c7.rds")
BM4 <- readRDS(file = "/home/yushiya/data/cd34/data/2019_Blood_BM/BM4_c7.rds")
BM4.loom <- as.loom(BM4, filename = "/home/yushiya/data/cd34/data/2019_Blood_BM/BM4_c7.loom", verbose = FALSE)
BM4.loom$close_all()

bmmc.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2020_Qin/bmmc-c7/outs/filtered_feature_bc_matrix")
bmmc <- CreateSeuratObject(counts = bmmc.data, project = "CD34")
bmmc$stim <- 'bmmc'
bmmc[["percent.mt"]] <- PercentageFeatureSet(bmmc, pattern = "^MT-")
VlnPlot(bmmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
bmmc[["percent.CD34"]] <- PercentageFeatureSet(bmmc, pattern = "CD34")
VlnPlot(bmmc, features = c("percent.CD34"))
bmmc <- subset(bmmc, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(bmmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
bmmc <- NormalizeData(bmmc, normalization.method = "LogNormalize", scale.factor = 10000)
bmmc <- FindVariableFeatures(bmmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(bmmc)
bmmc <- ScaleData(bmmc, features = all.genes)
bmmc <- RunPCA(bmmc, features = VariableFeatures(object = bmmc))
VizDimLoadings(bmmc, dims = 1:5, reduction = "pca")
ElbowPlot(bmmc)
DimPlot(bmmc, reduction = "pca")
bmmc <- RunUMAP(bmmc, reduction = "pca", dims = 1:30)
bmmc <- RunTSNE(bmmc, reduction = "pca", dims = 1:30)
bmmc <- FindNeighbors(bmmc, reduction = "pca", dims = 1:30)
bmmc <- FindClusters(bmmc, resolution = 0.9)
DimPlot(bmmc, reduction = "umap", label=T)
DimPlot(bmmc, reduction = "tsne", label=T)
DimPlot(bmmc,reduction = "umap", group.by = "cd34_clus", label=TRUE)
DimPlot(bmmc,reduction = "tsne", group.by = "cd34_clus", label=TRUE)
FeaturePlot(bmmc, features = c("CD34"), max.cutoff = 3, cols = c("grey", "red"))
markers.to.plot <- c("AVP","CRHBP","MEG3","EMCN","THY1","CD34","CD38","PTPRC","FLT3","KIT","GATA1","GATA2","ITGA2B","KLF1","TFRC","HBD","CD36","APOE","CA1","AHSP","GYPA","HEMGN","SELP","GP1BA","GP9","PF4","TPSAB1","HDC","MS4A2","EPX","CSF3R","CSF1R","CD33","MPO",
                     "ELANE","AZU1","PRTN3","LYZ","IRF8","IL3RA","SPIB","CD68","CD14","MME","IL7R","VPREB1","VPREB3","RAG2","EBF1","PAX5","CD24","CD19","MS4A1","CD79A","CD79B","IGHM","IGHD")
plot3=DotPlot(bmmc, features = markers.to.plot, cols = c("grey", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave("Dotplot.pdf", plot=plot3, width=13, height=6)
saveRDS(bmmc, file = "/home/yushiya/data/cd34/data/2019_Blood_BM/bmmc_c7.rds")
bmmc <- readRDS(file = "/home/yushiya/data/cd34/data/2019_Blood_BM/bmmc_c7.rds")
bmmc.loom <- as.loom(bmmc, filename = "/home/yushiya/data/cd34/data/2019_Blood_BM/bmmc_c7.loom", verbose = FALSE)
bmmc.loom$close_all()

BM.anchors <- FindIntegrationAnchors(object.list = c(BM1,BM2,BM3,BM4,bmmc), dims = 1:20)
BM.combined <- IntegrateData(anchorset = BM.anchors, dims = 1:20)
DefaultAssay(BM.combined) <- "integrated"
BM.combined$ori_stim <- BM.combined@meta.data$stim
BM.combined@meta.data$stim <- "BM"
# Run the standard workflow for visualization and clustering
BM.combined <- ScaleData(BM.combined, verbose = FALSE)
BM.combined <- RunPCA(BM.combined, npcs = 30, verbose = FALSE)
ElbowPlot(BM.combined)
#Clustering
BM.combined <- FindNeighbors(BM.combined, reduction = "pca", dims = 1:28)
BM.combined <- FindClusters(BM.combined, resolution = 0.9)
# umap
BM.combined <- RunUMAP(BM.combined, reduction = "pca", dims = 1:28)
#tsne
BM.combined <- RunTSNE(BM.combined, reduction = "pca", dims = 1:28)
# Visualization
DimPlot(BM.combined, reduction = "umap", group.by = "ori_stim")
DimPlot(BM.combined, reduction = "umap", split.by = "ori_stim", ncol=3)
DimPlot(BM.combined, reduction = "umap", label = TRUE)
DimPlot(BM.combined, reduction = "tsne", label = TRUE)
DimPlot(BM.combined, reduction = "tsne", group.by = "ori_stim")
DimPlot(BM.combined, reduction = "umap", group.by = "cd34_clus", label=TRUE)
DimPlot(BM.combined, reduction = "tsne", group.by = "cd34_clus", label=TRUE)
table(BM.combined@meta.data[["ori_stim"]],BM.combined@meta.data[["seurat_clusters"]])
DefaultAssay(BM.combined) <- "RNA"
BM.combined.markers <- FindAllMarkers(BM.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top20 <- BM.combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
plot1=DoHeatmap(BM.combined, features = top20$gene, assay = "integrated") + NoLegend()
ggsave("top10_markers.pdf", plot=plot1, width=16, height=14)
write.csv(BM.combined.markers, "top20_diff_genes.csv", row.names = F)
FeaturePlot(BM.combined, features = c("CD34"), max.cutoff = 3, cols = c("grey", "red"))
FeaturePlot(BM.combined, features = c("PLEK","HBB","MPO","SPIB","CD79A","DNTT"), max.cutoff = 3, cols = c("grey", "blue"))
VlnPlot(BM.combined, features = "CD34")
table(Idents(subset(BM.combined, subset= CD34>0)))/table(Idents(BM.combined))
markers.to.plot <- c("AVP","CRHBP","MEG3","EMCN","THY1","CD34","CD38","PTPRC","FLT3","KIT","GATA1","GATA2","ITGA2B","KLF1","TFRC","HBD","CD36","APOE","CA1","AHSP","GYPA","HEMGN","SELP","GP1BA","GP9","PF4","TPSAB1","HDC","MS4A2","EPX","CSF3R","CSF1R","CD33","MPO",
                     "ELANE","AZU1","PRTN3","LYZ","IRF8","IL3RA","SPIB","CD68","CD14","MME","IL7R","VPREB1","VPREB3","RAG2","EBF1","PAX5","CD24","CD19","MS4A1","CD79A","CD79B","IGHM","IGHD","FCER2","IGHA1","IGHG1","SDC1","TNFRSF17","CD7","CD8A","CD8B","CD3D",
                     "CD3E","CD3G","CD247","KLRB1","NCR3","KLRD1","GZMH")
markers.to.plot <- c("AVP","CRHBP","MEG3","EMCN","THY1","CD34","CD38","PTPRC","FLT3","KIT","GATA1","GATA2","ITGA2B","KLF1","TFRC","HBD","CD36","APOE","CA1","AHSP","GYPA","HEMGN","SELP","GP1BA","GP9","PF4","TPSAB1","HDC","MS4A2","EPX","CSF3R","CSF1R","CD33","MPO",
                     "ELANE","AZU1","PRTN3","LYZ","IRF8","IL3RA","SPIB","CD68","CD14","MME","IL7R","VPREB1","VPREB3","RAG2","EBF1","PAX5","CD24","CD19","MS4A1","CD79A","CD79B","IGHM","IGHD")
plot3=DotPlot(BM.combined, features = markers.to.plot, cols = c("grey", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave("Dotplot.pdf", plot=plot3, width=13, height=6)
BM.combined <- subset(BM.combined, subset= seurat_clusters=='0' | seurat_clusters=='1' | seurat_clusters=='2' | seurat_clusters=='3' | seurat_clusters=='4' | seurat_clusters=='5' | seurat_clusters=='6' | seurat_clusters=='7' | seurat_clusters=='8' | seurat_clusters=='9' | seurat_clusters=='10'| seurat_clusters=='11' |seurat_clusters=='12' | seurat_clusters=='13' | seurat_clusters=='14'| seurat_clusters=='15'| seurat_clusters=='16'| seurat_clusters=='17'| seurat_clusters=='18' | seurat_clusters=='19' | seurat_clusters=='20'| seurat_clusters=='21' | seurat_clusters=='22'| seurat_clusters=='23' | seurat_clusters=='24')
saveRDS(BM.combined, file = "/home/yushiya/data/cd34/data/2019_Blood_BM/BM_combined-c7_ori.rds")
BM.combined <- readRDS(file = "/home/yushiya/data/cd34/data/2019_Blood_BM/BM_combined-c7_ori.rds")

mPB1.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2020_Mead_PBMC/sample1/mPB1/outs/filtered_feature_bc_matrix")
mPB1 <- CreateSeuratObject(counts = mPB1.data, project = "CD34")
mPB1$stim <- 'mPB1'
mPB1[["percent.mt"]] <- PercentageFeatureSet(mPB1, pattern = "^MT-")
VlnPlot(mPB1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
mPB1 <- subset(mPB1, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(mPB1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
mPB1 <- NormalizeData(mPB1, normalization.method = "LogNormalize", scale.factor = 10000)
mPB1 <- FindVariableFeatures(mPB1, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(mPB1)
mPB1 <- ScaleData(mPB1, features = all.genes)
mPB1 <- RunPCA(mPB1, features = VariableFeatures(object = mPB1))
VizDimLoadings(mPB1, dims = 1:5, reduction = "pca")
ElbowPlot(mPB1)
DimPlot(mPB1, reduction = "pca")
mPB1 <- RunUMAP(mPB1, reduction = "pca", dims = 1:26)
mPB1 <- FindNeighbors(mPB1, reduction = "pca", dims = 1:26)
mPB1 <- FindClusters(mPB1, resolution = 0.8)
DimPlot(mPB1, reduction = "umap",label=T)
FeaturePlot(mPB1, features = c("IGHM","GATA2"), max.cutoff = 3, cols = c("grey", "red"))
VlnPlot(mPB1, features = "CD34")
table(Idents(subset(mPB1, subset= CD34>0)))/table(Idents(mPB2))
markers.to.plot <- c("AVP","CRHBP","MEG3","EMCN","THY1","CD34","CD38","PTPRC","FLT3","KIT","GATA1","GATA2","ITGA2B","KLF1","TFRC","HBD","CD36","APOE","CA1","AHSP","GYPA","HEMGN","SELP","GP1BA","GP9","PF4","TPSAB1","HDC","MS4A2","EPX","CSF3R","CSF1R","CD33","MPO",
                     "ELANE","AZU1","PRTN3","LYZ","IRF8","IL3RA","SPIB","CD68","CD14","MME","IL7R","VPREB1","VPREB3","RAG2","EBF1","PAX5","CD24","CD19","MS4A1","CD79A","CD79B","IGHM","IGHD")
plot3=DotPlot(mPB1, features = markers.to.plot, cols = c("grey", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave("Dotplot.pdf", plot=plot3, width=13, height=6)
saveRDS(mPB1, file = "/home/yushiya/data/cd34/data/2020_Mead_PBMC/mPB1.rds")
mPB1 <- readRDS(file = "/home/yushiya/data/cd34/data/2020_Mead_PBMC/mPB1.rds")
mPB1.loom <- as.loom(mPB1, filename = "/home/yushiya/data/cd34/data/2020_Mead_PBMC/mPB1.loom", verbose = FALSE)
mPB1.loom$close_all()

mPB2.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2020_Mead_PBMC/sample6/mPB2/outs/filtered_feature_bc_matrix")
mPB2 <- CreateSeuratObject(counts = mPB2.data, project = "CD34")
mPB2$stim <- 'mPB2'
mPB2[["percent.mt"]] <- PercentageFeatureSet(mPB2, pattern = "^MT-")
VlnPlot(mPB2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
mPB2 <- subset(mPB2, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(mPB2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
mPB2 <- NormalizeData(mPB2, normalization.method = "LogNormalize", scale.factor = 10000)
mPB2 <- FindVariableFeatures(mPB2, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(mPB2)
mPB2 <- ScaleData(mPB2, features = all.genes)
mPB2 <- RunPCA(mPB2, features = VariableFeatures(object = mPB2))
VizDimLoadings(mPB2, dims = 1:5, reduction = "pca")
ElbowPlot(mPB2)
DimPlot(mPB2, reduction = "pca")
mPB2 <- RunUMAP(mPB2, reduction = "pca", dims = 1:26)
mPB2 <- FindNeighbors(mPB2, reduction = "pca", dims = 1:26)
mPB2 <- FindClusters(mPB2, resolution = 0.8)
DimPlot(mPB2, reduction = "umap",label=T)
FeaturePlot(mPB2, features = c("IGHM","GATA2"), max.cutoff = 3, cols = c("grey", "red"))
VlnPlot(mPB2, features = "CD34")
table(Idents(subset(mPB2, subset= CD34>0)))/table(Idents(mPB2))
markers.to.plot <- c("AVP","CRHBP","MEG3","EMCN","THY1","CD34","CD38","PTPRC","FLT3","KIT","GATA1","GATA2","ITGA2B","KLF1","TFRC","HBD","CD36","APOE","CA1","AHSP","GYPA","HEMGN","SELP","GP1BA","GP9","PF4","TPSAB1","HDC","MS4A2","EPX","CSF3R","CSF1R","CD33","MPO",
                     "ELANE","AZU1","PRTN3","LYZ","IRF8","IL3RA","SPIB","CD68","CD14","MME","IL7R","VPREB1","VPREB3","RAG2","EBF1","PAX5","CD24","CD19","MS4A1","CD79A","CD79B","IGHM","IGHD")
plot3=DotPlot(mPB2, features = markers.to.plot, cols = c("grey", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave("Dotplot.pdf", plot=plot3, width=13, height=6)
saveRDS(mPB2, file = "/home/yushiya/data/cd34/data/2020_Mead_PBMC/mPB2.rds")
mPB2 <- readRDS(file = "/home/yushiya/data/cd34/data/2020_Mead_PBMC/mPB2.rds")
mPB2.loom <- as.loom(mPB2, filename = "/home/yushiya/data/cd34/data/2020_Mead_PBMC/mPB2.loom", verbose = FALSE)
mPB2.loom$close_all()

mPB3.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2020_Mead_PBMC/sample9/mPB3/outs/filtered_feature_bc_matrix")
mPB3 <- CreateSeuratObject(counts = mPB3.data, project = "CD34")
mPB3$stim <- 'mPB3'
mPB3[["percent.mt"]] <- PercentageFeatureSet(mPB3, pattern = "^MT-")
VlnPlot(mPB3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
mPB3 <- subset(mPB3, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(mPB3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
mPB3 <- NormalizeData(mPB3, normalization.method = "LogNormalize", scale.factor = 10000)
mPB3 <- FindVariableFeatures(mPB3, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(mPB3)
mPB3 <- ScaleData(mPB3, features = all.genes)
mPB3 <- RunPCA(mPB3, features = VariableFeatures(object = mPB3))
VizDimLoadings(mPB3, dims = 1:5, reduction = "pca")
ElbowPlot(mPB3)
DimPlot(mPB3, reduction = "pca")
mPB3 <- RunUMAP(mPB3, reduction = "pca", dims = 1:26)
mPB3 <- FindNeighbors(mPB3, reduction = "pca", dims = 1:26)
mPB3 <- FindClusters(mPB3, resolution = 0.8)
DimPlot(mPB3, reduction = "umap",label=T)
FeaturePlot(mPB3, features = c("IGHM","GATA2"), max.cutoff = 3, cols = c("grey", "red"))
VlnPlot(mPB3, features = "CD34")
table(Idents(subset(mPB3, subset= CD34>0)))/table(Idents(mPB3))
markers.to.plot <- c("AVP","CRHBP","MEG3","EMCN","THY1","CD34","CD38","PTPRC","FLT3","KIT","GATA1","GATA2","ITGA2B","KLF1","TFRC","HBD","CD36","APOE","CA1","AHSP","GYPA","HEMGN","SELP","GP1BA","GP9","PF4","TPSAB1","HDC","MS4A2","EPX","CSF3R","CSF1R","CD33","MPO",
                     "ELANE","AZU1","PRTN3","LYZ","IRF8","IL3RA","SPIB","CD68","CD14","MME","IL7R","VPREB1","VPREB3","RAG2","EBF1","PAX5","CD24","CD19","MS4A1","CD79A","CD79B","IGHM","IGHD")
plot3=DotPlot(mPB3, features = markers.to.plot, cols = c("grey", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave("Dotplot.pdf", plot=plot3, width=13, height=6)
saveRDS(mPB3, file = "/home/yushiya/data/cd34/data/2020_Mead_PBMC/mPB3.rds")
mPB3 <- readRDS(file = "/home/yushiya/data/cd34/data/2020_Mead_PBMC/mPB3.rds")
mPB3.loom <- as.loom(mPB3, filename = "/home/yushiya/data/cd34/data/2020_Mead_PBMC/mPB3.loom", verbose = FALSE)
mPB3.loom$close_all()

mPB4.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2020_Mead_PBMC/sample13/mPB4/outs/filtered_feature_bc_matrix")
mPB4 <- CreateSeuratObject(counts = mPB4.data, project = "CD34")
mPB4$stim <- 'mPB4'
mPB4[["percent.mt"]] <- PercentageFeatureSet(mPB4, pattern = "^MT-")
VlnPlot(mPB4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
mPB4 <- subset(mPB4, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(mPB4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
mPB4 <- NormalizeData(mPB4, normalization.method = "LogNormalize", scale.factor = 10000)
mPB4 <- FindVariableFeatures(mPB4, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(mPB4)
mPB4 <- ScaleData(mPB4, features = all.genes)
mPB4 <- RunPCA(mPB4, features = VariableFeatures(object = mPB4))
VizDimLoadings(mPB4, dims = 1:5, reduction = "pca")
ElbowPlot(mPB4)
DimPlot(mPB4, reduction = "pca")
mPB4 <- RunUMAP(mPB4, reduction = "pca", dims = 1:26)
mPB4 <- FindNeighbors(mPB4, reduction = "pca", dims = 1:26)
mPB4 <- FindClusters(mPB4, resolution = 0.8)
DimPlot(mPB4, reduction = "umap",label=T)
FeaturePlot(mPB4, features = c("IGHM","GATA2"), max.cutoff = 3, cols = c("grey", "red"))
VlnPlot(mPB4, features = "CD34")
table(Idents(subset(mPB4, subset= CD34>0)))/table(Idents(BM_OD4))
markers.to.plot <- c("AVP","CRHBP","MEG3","EMCN","THY1","CD34","CD38","PTPRC","FLT3","KIT","GATA1","GATA2","ITGA2B","KLF1","TFRC","HBD","CD36","APOE","CA1","AHSP","GYPA","HEMGN","SELP","GP1BA","GP9","PF4","TPSAB1","HDC","MS4A2","EPX","CSF3R","CSF1R","CD33","MPO",
                     "ELANE","AZU1","PRTN3","LYZ","IRF8","IL3RA","SPIB","CD68","CD14","MME","IL7R","VPREB1","VPREB3","RAG2","EBF1","PAX5","CD24","CD19","MS4A1","CD79A","CD79B","IGHM","IGHD")
plot3=DotPlot(mPB4, features = markers.to.plot, cols = c("grey", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave("Dotplot.pdf", plot=plot3, width=13, height=6)
saveRDS(mPB4, file = "/home/yushiya/data/cd34/data/2020_Mead_PBMC/mPB4.rds")
mPB4 <- readRDS(file = "/home/yushiya/data/cd34/data/2020_Mead_PBMC/mPB4.rds")
mPB4.loom <- as.loom(mPB4, filename = "/home/yushiya/data/cd34/data/2020_Mead_PBMC/mPB4.loom", verbose = FALSE)
mPB4.loom$close_all()

mPB.anchors <- FindIntegrationAnchors(object.list = c(mPB1,mPB2,mPB3,mPB4), dims = 1:20)
mPB.combined <- IntegrateData(anchorset = mPB.anchors, dims = 1:20)
DefaultAssay(mPB.combined) <- "integrated"
mPB.combined$ori_stim <- mPB.combined@meta.data$stim
mPB.combined@meta.data$stim <- "mPB"
# Run the standard workflow for visualization and clustering
mPB.combined <- ScaleData(mPB.combined, verbose = FALSE)
mPB.combined <- RunPCA(mPB.combined, npcs = 30, verbose = FALSE)
ElbowPlot(mPB.combined)
#Clustering
mPB.combined <- FindNeighbors(mPB.combined, reduction = "pca", dims = 1:28)
mPB.combined <- FindClusters(mPB.combined, resolution = 0.9)
# umap
mPB.combined <- RunUMAP(mPB.combined, reduction = "pca", dims = 1:28)
#tsne
mPB.combined <- RunTSNE(mPB.combined, reduction = "pca", dims = 1:28)
# Visualization
DimPlot(mPB.combined, reduction = "umap", group.by = "ori_stim")
DimPlot(mPB.combined, reduction = "umap", split.by = "ori_stim", ncol=2)
DimPlot(mPB.combined, reduction = "umap", label = TRUE)
DimPlot(mPB.combined, reduction = "tsne", label = TRUE)
DimPlot(mPB.combined, reduction = "tsne", group.by = "ori_stim")
DimPlot(mPB.combined, reduction = "umap", group.by = "cd34_clus", label=TRUE)
DimPlot(mPB.combined, reduction = "tsne", group.by = "cd34_clus", label=TRUE)
table(mPB.combined@meta.data[["ori_stim"]],mPB.combined@meta.data[["seurat_clusters"]])
DefaultAssay(mPB.combined) <- "RNA"
mPB.combined.markers <- FindAllMarkers(mPB.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top20 <- mPB.combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
plot1=DoHeatmap(mPB.combined, features = top20$gene, assay = "integrated") + NoLegend()
ggsave("top10_markers.pdf", plot=plot1, width=16, height=14)
write.csv(mPB.combined.markers, "top20_diff_genes.csv", row.names = F)
FeaturePlot(mPB.combined, features = c("CD34"), max.cutoff = 3, cols = c("grey", "red"))
FeaturePlot(mPB.combined, features = c("PLEK","HBB","MPO","SPIB","CD79A","DNTT"), max.cutoff = 3, cols = c("grey", "blue"))
VlnPlot(mPB.combined, features = "CD34")
table(Idents(subset(mPB.combined, subset= CD34>0)))/table(Idents(mPB.combined))
markers.to.plot <- c("AVP","CRHBP","MEG3","EMCN","THY1","CD34","CD38","PTPRC","FLT3","KIT","GATA1","GATA2","ITGA2B","KLF1","TFRC","HBD","CD36","APOE","CA1","AHSP","GYPA","HEMGN","SELP","GP1BA","GP9","PF4","TPSAB1","HDC","MS4A2","EPX","CSF3R","CSF1R","CD33","MPO",
                     "ELANE","AZU1","PRTN3","LYZ","IRF8","IL3RA","SPIB","CD68","CD14","MME","IL7R","VPREB1","VPREB3","RAG2","EBF1","PAX5","CD24","CD19","MS4A1","CD79A","CD79B","IGHM","IGHD","FCER2","IGHA1","IGHG1","SDC1","TNFRSF17","CD7","CD8A","CD8B","CD3D",
                     "CD3E","CD3G","CD247","KLRB1","NCR3","KLRD1","GZMH")
markers.to.plot <- c("AVP","CRHBP","MEG3","EMCN","THY1","CD34","CD38","PTPRC","FLT3","KIT","GATA1","GATA2","ITGA2B","KLF1","TFRC","HBD","CD36","APOE","CA1","AHSP","GYPA","HEMGN","SELP","GP1BA","GP9","PF4","TPSAB1","HDC","MS4A2","EPX","CSF3R","CSF1R","CD33","MPO",
                     "ELANE","AZU1","PRTN3","LYZ","IRF8","IL3RA","SPIB","CD68","CD14","MME","IL7R","VPREB1","VPREB3","RAG2","EBF1","PAX5","CD24","CD19","MS4A1","CD79A","CD79B","IGHM","IGHD")
plot3=DotPlot(mPB.combined, features = markers.to.plot, cols = c("grey", "red"), dot.scale = 8) + 
  RotatedAxis()
ggsave("Dotplot.pdf", plot=plot3, width=13, height=6)
mPB.combined <- subset(mPB.combined, subset= seurat_clusters=='0' | seurat_clusters=='1' | seurat_clusters=='2' | seurat_clusters=='3' | seurat_clusters=='4' | seurat_clusters=='5' | seurat_clusters=='6' | seurat_clusters=='7' | seurat_clusters=='8' | seurat_clusters=='9' | seurat_clusters=='10'| seurat_clusters=='11' |seurat_clusters=='12' | seurat_clusters=='13' | seurat_clusters=='14'| seurat_clusters=='15'| seurat_clusters=='16'| seurat_clusters=='17'| seurat_clusters=='18' | seurat_clusters=='19' | seurat_clusters=='20'| seurat_clusters=='21' | seurat_clusters=='22'| seurat_clusters=='23' | seurat_clusters=='24')
saveRDS(mPB.combined, file = "/home/yushiya/data/cd34/data/2020_Mead_PBMC/mPB.combined-c7.rds")
mPB.combined <- readRDS(file = "/home/yushiya/data/cd34/data/2020_Mead_PBMC/mPB.combined-c7.rds")

#Integrate by Seurat CCA
immune.anchors <- FindIntegrationAnchors(object.list = c(thymus.combined,PB.combined,BM.combined,mPB.combined), assay = c("RNA","RNA","integrated","integrated"), k.anchor = 10, k.score = 50, k.filter = 100, anchor.features = 2000, max.features = 200, dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)
DefaultAssay(immune.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 50, verbose = FALSE)
ElbowPlot(immune.combined)
#Clustering
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:40)
immune.combined <- FindClusters(immune.combined, reduction = "pca", resolution = 0.95)
#umap
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:40)
immune.combined <- RunTSNE(immune.combined, reduction = "pca", dims = 1:30)
#cell cycle 
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
immune.combined <- CellCycleScoring(immune.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
immune.combined$Phase <- factor(immune.combined$Phase, levels=c("G1","S","G2M"))
#colours
#HSPCs "#550000", "#AA3939", "#CC6F66", "#E79492","#FBD2CE"
#EryMk "#FFBD54","#B97802","#915900"
#B "#A5D6A7", "#47AF50", "#2E7D32",
#T "#DDA0DD","#98ADC4","#7F689D","#54689A",
#My "#FFF176","#CECA68", "#B2A13F",
color_all <- c("#550000", "#AA3939", "#CC6F66","#FFBD54","#B97802","#915900","#E79492","#FBD2CE","#FFF176","#CECA68", "#B2A13F","#DDA0DD","#98ADC4","#7F689D","#54689A","#A5D6A7", "#47AF50", "#2E7D32")
color_BMPB_lineages <- c("#550000", "#AA3939", "#CC6F66","#E79492","#FBD2CE","#FFBD54","#B97802","#915900","#FFF176","#CECA68", "#B2A13F","#A5D6A7", "#47AF50", "#2E7D32")
color_stim <- c("#DB5C25","#F3B747","#649541","#4C82C5")
color_cycle <- c("#EA63A2","#FDD685","#52C6EC")
#visualization
DimPlot(immune.combined, reduction = "umap", label = T, repel=T)
DimPlot(immune.combined, reduction = "umap", label = T, cols=color_all, repel=T)
ggsave("/home/yushiya/cd34/fig/fig1-integrated.pdf", plot=p1, width=6, height=4)
DimPlot(immune.combined, reduction = "umap", group.by="stim", cols=color_stim)
ggsave("/home/yushiya/cd34/fig/fig1-groupBYstim.pdf", plot=p1, width=6, height=4)
DimPlot(immune.combined, reduction = "umap", split.by="stim", ncol=2)
ggsave("/home/yushiya/cd34/fig/fig1-splitBYstim.pdf", plot=p1, width=16, height=4)
p1=DimPlot(immune.combined, reduction = "umap", group.by="Phase", cols=color_cycle)
ggsave("/home/yushiya/cd34/fig/fig1-groupBYphase.pdf", plot=p1, width=6, height=4)
p1=DimPlot(sub1, reduction = "umap", split.by="stim", ncol=1, cols=color_all)
ggsave("/home/yushiya/cd34/fig/fig2-BMmPBPB.pdf", plot=p1, width=6, height=10)
p1=DimPlot(immune.combined, reduction = "umap", split.by="stim", ncol=1, cols=color_all)
ggsave("/home/yushiya/cd34/fig/fig6-allUMAP.pdf", plot=p1, width=6, height=13)
DimPlot(immune.combined, reduction = "umap", split.by="ori_stim", ncol=2)
DimPlot(immune.combined, reduction = "tsne", label = T)
DimPlot(immune.combined, reduction = "tsne", group.by="stim")
DefaultAssay(immune.combined) <- "RNA"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
table(immune.combined@meta.data[["stim"]],immune.combined@meta.data[["seurat_clusters"]])
VlnPlot(immune.combined, features = "CD34", combine = FALSE)
p1=FeaturePlot(immune.combined, features = c("CD34"), max.cutoff = 3, cols = c("#FFE9AE", "#981C12"))
ggsave("/home/yushiya/data/cd34/manu_fig/fig1-CD34.pdf", plot=p1, width=6, height=4)
p1=FeaturePlot(immune.combined, features = c("AVP","MPO","IRF8","CD7","CD3D","CD1A","CD79A","GATA2","PF4","PPBP","CCR7"),pt.size=0.001,max.cutoff = 3, cols = c("#C5E9F5", "#981C12"), order=T)
ggsave("/home/yushiya/cd34/fig/fig1-Feamarkers-1.pdf", plot=p1, width=10, height=5.5)
immune.combined.markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(immune.combined.markers, "/home/yushiya/data/cd34/manu_fig/v4/celltype_markers.csv")
Idents(immune.combined) <- immune.combined$seurat_clusters
immune.combined.markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(immune.combined.markers, "/home/yushiya/data/cd34/manu_fig/v4/seurat-cluster_markers.csv")
markers.to.plot <- c("CD34","AVP","CRHBP","MEG3","EMCN","KIT","GATA1","GATA2","TFRC","HBD","CD36","PPBP","PF4","CSF3R","CSF1R","MPO",
                     "ELANE","LYZ","IRF8","IL3RA","SPIB","CD68","CCR7","NKG7","IL7R","CD19","VPREB1","EBF1","PAX5","MS4A1","CD79A",
                     "CD7","CD3D","CD3E","RAG1","RAG2","CD1A")
plot3=DotPlot(immune.combined, features = markers.to.plot, cols = c(c("grey","red"),c("grey","red"),c("grey", "red")), dot.scale = 8) + 
  RotatedAxis()
ggsave("Dotplot.pdf", plot=plot3, width=16, height=7)
ggsave("/home/yushiya/cd34/fig/fig1-markers.pdf", plot=plot3, width=12, height=6.5)
immune.combined <- RenameIdents(immune.combined,'0'='HSC','3'='HSC','6'='MPP #1','7'='MPP #1','1'='MPP #2','2'='MPP #2','13'='MPP #2','5'='BMEP','21'='BMEP','12'='EryP','24'='MkP',
                                '4'='LMPP #1','15'='LMPP #2','18'='GMP','23'='pre-cDC','22'='pre-pDC','8'='ETP','10'='mThy','11'='cThy #1','17'='cThy #2','19'='cThy #2',
                                '14'='pre-pro B','9'='pro-B','16'='pre-B','20'='pre-B')
immune.combined$celltype <- factor(immune.combined$celltype, levels=c("HSC","MPP #1","MPP #2","MEP","EryP","MkP","LMPP #1","LMPP #2","GMP","CDP","pre-pDC","ETP","Thy #1","Thy #2","Thy #3","CLP","pro-B","pre-B"))
immune.combined$celltype <- factor(immune.combined$celltype, levels=c("HSC","MPP #1","MPP #2","LMPP #1","LMPP #2","MEP","EryP","MkP","GMP","CDP","pre-pDC","ETP","Thy #1","Thy #2","Thy #3","CLP","pro-B","pre-B"))
saveRDS(immune.combined, file = "/home/yushiya/cd34/fig/HSPC.rds")

#monocle3
dyn.load("/home/wangjl/local/gdal-2.4.4/libgdal.so.20.5.4")
library(monocle3)
##创建CDS对象并预处理数据
data <- GetAssayData(immune.combined, assay = 'RNA', slot = 'counts')
cell_metadata <- immune.combined@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
save(cds,file="/home/yushiya/cd34/fig/immune.combined.cds-2.RData")
load(file="/home/yushiya/cd34/fig/immune.combined.cds-2.RData")
#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)
#umap降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('cds.umap')
##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(immune.combined, reduction = "umap")
#int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
## Monocle3聚类分区
cds <- cluster_cells(cds)
cds@clusters@listData[["UMAP"]][["clusters"]] <- immune.combined$celltype
## 识别轨迹
cds <- learn_graph(cds)
p<-plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
              label_branch_points = FALSE)
# cds <- order_cells(cds) 存在bug，使用辅助线选择root细胞
p + geom_vline(xintercept = seq(-7,-6,0.25)) + geom_hline(yintercept = seq(-2,-1,0.25))
embed <- data.frame(Embeddings(immune.combined, reduction = "umap"))
embed <- subset(embed, UMAP_1 > -7 & UMAP_1 < -6.75 & UMAP_2 > -1.75 & UMAP_2 < -1.5)
root.cell <- rownames(embed)
cds <- order_cells(cds, root_cells = root.cell)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
              label_leaves = FALSE,  label_branch_points = FALSE)
plot_cells(cds, color_cells_by = "celltype", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = FALSE)
immune.combined$pseudotime <- pseudotime(cds)
ggsave("/home/yushiya/cd34/fig/fig1-pseu.pdf", plot=p1, width=6, height=4)
#find genes through pseudotime
sub_cds <- cds[,colData(cds)$stim %in% c("BM","mPB","PB")]
sub_cds <- sub_cds[,colData(sub_cds)$celltype %in% c("HSC","MPP #1","MPP #2","LMPP #1","LMPP #2","MEP","EryP","MkP","CLP","pro-B","pre-B","GMP","CDP","pre-pDC")]
plot_cells(sub_cds, color_cells_by="celltype")
#pr_graph_test_res <- graph_test(sub_cds, neighbor_graph="knn", cores=40)
#pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))
#write.csv(pr_graph_test_res, "/home/yushiya/cd34/fig/peu_gene.csv")
cds_pr_test_res <- graph_test(sub_cds, neighbor_graph="principal_graph", cores=80)  #DEGs along trajectory
write.csv(cds_pr_test_res, "/home/yushiya/cd34/fig/peu_tra_gene.csv")
cds_pr_test_res <- read.csv("/home/yushiya/cd34/fig/peu_tra_gene.csv", row.names = 1)
pr_deg_ids <- row.names(subset(cds_pr_test_res, q_value < 0.01 & morans_I > 0.2))
#heatmap along peudotime
plot_matrix <- exprs(sub_cds)[match(pr_deg_ids,#挑选我们前面确定需要展示的基因
                                   rownames(rowData(sub_cds))),
                             order(pseudotime(sub_cds))]#按照拟时对基因排序
plot_matrix <- t(apply(plot_matrix,1,function(x){smooth.spline(x,df=3)$y}))
plot_matrix <- t(apply(plot_matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(plot_matrix) <- pr_deg_ids
dim(plot_matrix)
source('/home/yushiya/software/monocle3_heatmap.R')
plot_matrix_combin <- cutColumn_Means(plot_matrix,cut = 50)
##排序设置
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}
p1 <- pheatmap::pheatmap(plot_matrix_combin, 
                   useRaster = T,
                   cluster_cols=FALSE, 
                   cluster_rows=T, 
                   show_rownames=F, 
                   show_colnames=F, 
                   clustering_method = "ward.D2",
                   cutree_rows=7,
                   filename=NA,
                   border_color = NA,
                   fontsize_row = 8,
                   color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
                   clustering_callback = callback)
#行注释
annotation_row <- data.frame(Cluster=factor(cutree(p1$tree_row, 7)))
row.names(annotation_row) <- rownames(plot_matrix_combin)
rowcolor <- brewer.pal(7,"Paired")
names(rowcolor) <- c("1","2","3","4","5","7") #类型颜色
#注释颜色设置
ann_colors <- list(Cluster=rowcolor) #颜色设置
p2 <- pheatmap::pheatmap(plot_matrix_combin, 
                   cluster_cols=FALSE, 
                   cluster_rows=T, 
                   show_rownames=F, 
                   show_colnames=F, 
                   clustering_method = "ward.D2",
                   filename=NA,
                   cutree_rows=7,
                   border_color = NA,
                   fontsize_row = 8,
                   color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
                   #annotation_colors=ann_colors,
                   annotation_row = annotation_row,
                   clustering_callback = callback,
                   annotation_names_col = F,
                   annotation_names_row = F,
                   main="Pseudotime")
###首先提取热图中各个module的基因
module_gene <- as.data.frame(cutree(p2$tree_row, k=7))
colnames(module_gene) <- "Module"
module_gene$gene <- rownames(module_gene)
#我们这里展示GO结果
Module_GO=data.frame()
for (i in unique(module_gene$Module)) {
  data=filter(module_gene,module_gene$Module==i)
  go <- enrichGO(gene= data$gene,
                 OrgDb= org.Hs.eg.db,
                 keyType= 'SYMBOL',
                 ont= "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff= 0.05,
                 qvalueCutoff= 0.05,
                 readable= TRUE)
  go_res=go@result
  
  if (dim(go_res)[1] != 0) {
    go_res$cluster=i
    Module_GO=rbind(Module_GO,go_res)
  }
}
Module_GO <- Module_GO[which(Module_GO$qvalue <= 0.05),]
write.csv(Module_GO, file = '/home/yushiya/cd34/fig/Module_GO_pheatmap.csv')
top_Module_GO <- Module_GO %>% group_by(cluster) %>% top_n(n = -10, wt = p.adjust)
Module_GO <- Module_GO[,c("ID","Description","qvalue","cluster")]
#展示需要的基因
gene <- c("S100A9","FCER1G","CXCL5","CTSG","CCR7","AZU1","PRTN3","CXCR4","GPR183","GAS6","LGMN","SMPD3","CXCR3")
p3 <- pheatmap::pheatmap(plot_matrix_combin, 
                         cluster_cols=FALSE, 
                         cluster_rows=T, 
                         show_rownames=T, 
                         show_colnames=F, 
                         clustering_method = "ward.D2",
                         filename=NA,
                         border_color = NA,
                         fontsize_row = 8,
                         color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
                         #annotation_colors=ann_colors,
                         annotation_row = annotation_row,
                         clustering_callback = callback,
                         annotation_names_col = F,
                         annotation_names_row = F,
                         main="Pseudotime")
source('/home/yushiya/software/add.flag.R')
library(grid)
p4=add.flag(p3,kept.labels = gene,repel.degree = 0.2)
ggsave("/home/yushiya/cd34/fig/fig2-tra-heatmap.pdf", plot=p4, width=6, height=7)
#find modules
gene_module_df <- find_gene_modules(sub_cds[pr_deg_ids,], resolution=c(10^seq(-8,-3)))
gene_module_df <- find_gene_modules(sub_cds[pr_deg_ids,], resolution=c(10^seq(-9,-4)))
table(gene_module_df$module)
write.csv(gene_module_df, "/home/yushiya/cd34/fig/peu_tra_module-2.csv")
gene_module_df <- read.csv("/home/yushiya/cd34/fig/peu_tra_module-2.csv", row.names = 1)
sub_cds1 <- sub_cds[,colData(sub_cds)$celltype %in% c("HSC","MPP #1","MPP #2","LMPP #1","LMPP #2","MEP","EryP","MkP","CLP","pro-B","pre-B","GMP","CDP","pre-pDC")]
new_type <- paste(colData(sub_cds)$stim,colData(sub_cds)$celltype,sep="_")
cell_group_df <- tibble::tibble(cell=row.names(colData(sub_cds)), 
                                cell_group=new_type)
agg_mat <- aggregate_gene_expression(sub_cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
#colnames(agg_mat) <- factor(colnames(agg_mat), levels=c("HSC","MPP #1","MPP #2","LMPP #1","LMPP #2","MEP","EryP","MkP","CLP","pro-B","pre-B","GMP","CDP","pre-pDC"))
p<-pheatmap::pheatmap(agg_mat,cluster_rows = T, cluster_cols = T,
                   scale="column", clustering_method="ward.D2")
ggsave("/home/yushiya/cd34/fig/fig2-tra-module-heatmap-1.pdf", plot=p, width=8, height=4)
#GO of modules
ID <- subset(gene_module_df, subset= module=="1")
ego1 <- enrichGO(gene         = ID$id,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
barplot(ego1, showCategory=10, drop=T)
dotplot(ego1, showCategory=10)
#GO terms in all modules
module_gene1 <- gene_module_df[,c(1,2)]
colnames(module_gene1) <- c("gene","Module")
rownames(module_gene1) <- module_gene1$gene
Module_GO1=data.frame()
for (i in unique(module_gene1$Module)) {
  data=filter(module_gene1,module_gene1$Module==i)
  go <- enrichGO(gene= data$gene,
                 OrgDb= org.Hs.eg.db,
                 keyType= 'SYMBOL',
                 ont= "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff= 0.05,
                 qvalueCutoff= 0.05,
                 readable= TRUE)
  go_res=go@result
  if (dim(go_res)[1] != 0) {
    go_res$cluster=i
    Module_GO1=rbind(Module_GO1,go_res)
  }
}
Module_GO1 <- Module_GO1[which(Module_GO1$qvalue <= 0.05),]
write.csv(Module_GO1, file = '/home/yushiya/cd34/fig/Module_GO_monocle3.csv')
top_Module_GO1 <- Module_GO1 %>% group_by(cluster) %>% top_n(n = -10, wt = p.adjust)
selected_GO <- top_Module_GO1[c(1,2,3,4,11,12,21,22,31,33,35,41,42,29),]
mtx_go <- data.frame(selected_GO$Description)
colnames(mtx_go)[1] <- "Description"
for (i in 1:5) {
  temp <- subset(selected_GO, subset= cluster==i)
  temp <- temp[,c("Description","GeneRatio")]
  split_num <- data.frame(t(data.frame(strsplit(temp$GeneRatio,"/"))))
  temp$GeneRatio <- as.numeric(split_num$X1)/as.numeric(split_num$X2)
  #temp$GeneRatio <- strsplit()
  mtx_go <- mtx_go %>% left_join(temp, by="Description")
  mtx_go$GeneRatio[is.na(mtx_go$GeneRatio)] <- 0
  colnames(mtx_go)[i+1] <- paste("Module",i)
}
rownames(mtx_go) <- selected_GO$Description
mtx_go <- mtx_go[,-1]
p1 <- pheatmap(mtx_go, cluster_cols = F, #cluster_rows = F,  
            clustering_method = "average",
            scale = "column",
            #border=F,
            #col=colorRampPalette(c("navy","firebrick3"))(50),
            fontsize_row = 12, fontsize_col= 12, angle_col = 45)
ggsave("/home/yushiya/cd34/fig/fig2-tra-module-GO.pdf", plot=p1, width=10, height=4.5)


##EryMk lineage
sub1 <- subset(immune.combined, subset= celltype=="HSC"|celltype=="MPP #1"|celltype=="MPP #2"|celltype=="BMEP"|celltype=="EryP"|celltype=="MkP")
##B lineage
sub1 <- subset(immune.combined, subset= celltype=="HSC"|celltype=="MPP #1"|celltype=="MPP #2"|celltype=="LMPP #2"|celltype=="CLP"|celltype=="pro-B"|celltype=="pre-B")
##T lineage
sub1 <- subset(immune.combined, subset= celltype=="HSC"|celltype=="MPP #1"|celltype=="MPP #2"|celltype=="LMPP #1"|celltype=="LMPP #2"|celltype=="ETP"|celltype=="Thy #1"|celltype=="Thy #2"|celltype=="Thy #3")
sub1 <- subset(immune.combined, subset= ((celltype=="HSC"|celltype=="MPP #1"|celltype=="MPP #2"|celltype=="LMPP #1")&(stim=="BM"|stim=="PB"|stim=="mPB"))
               |celltype=="ETP"
               |((celltype=="mThy"|celltype=="cThy #1"|celltype=="cThy #2")&stim=="thymus"))
##My lineage
sub1 <- subset(immune.combined, subset= celltype=="HSC"|celltype=="MPP #1"|celltype=="MPP #2"|celltype=="LMPP #1"|celltype=="GMP"|celltype=="pre-cDC"|celltype=="pre-pDC")
sub1 <- subset(immune.combined, subset= ((celltype=="HSC"|celltype=="MPP #1"|celltype=="MPP #2"|celltype=="LMPP #1"|celltype=="GMP"|celltype=="pre-cDC")&(stim=="BM"|stim=="PB"|stim=="mPB"))
               |celltype=="pre-pDC")
##BM mPB PB all lineages
sub1 <- subset(immune.combined, subset= celltype=="HSC"|celltype=="MPP #1"|celltype=="MPP #2"|celltype=="MEP"|celltype=="EryP"|celltype=="MkP"|celltype=="LMPP #1"|celltype=="LMPP #2"|celltype=="GMP"|celltype=="CDP"|celltype=="pre-pDC"|celltype=="CLP"|celltype=="pro-B"|celltype=="pre-B")
sub1$celltype <- factor(sub1$celltype, levels= c("HSC","MPP #1","MPP #2","LMPP #1","LMPP #2","MEP","EryP","MkP","GMP","CDP","pre-pDC","CLP","pro-B","pre-B"))


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
clus_percent$Clus <- factor(clus_percent$Clus,levels = c('HSC','MPP #1',"MPP #2",'BMEP','EryP','MkP','LMPP #1','LMPP #2','GMP','pre-cDC','pre-pDC',"ETP","mThy","cThy #1","cThy #2",'pre-pro B','pro-B','pre-B'))
p1=ggplot(clus_percent,
          aes(x=Stim, y=value*100, fill=Clus, stratum = Clus, alluvium = Clus)) +
  geom_bar(stat='identity', width=0.45) +
  geom_alluvium() +
  geom_stratum(width=0.45, size=0.1) +
  scale_fill_manual(values = color_all) +
  labs(x='Samples', y='Relative Abundance (%)')+
  scale_y_continuous(expand=c(0, 0))+
  theme(axis.text.x=element_text(hjust=0.5),
        text = element_text(size = 18),
        panel.background = element_blank(),
        panel.grid = element_line(colour = "grey90"))# + theme_bw()
ggsave("/home/yushiya/data/cd34/manu_fig/v4/fig1-clusP.pdf", plot=p1, width=6, height=6)
#lineage results with x=celltype fill=stim
clus_percent <- table(sub1@meta.data[["celltype"]],sub1@meta.data[["stim"]])
clus_percent <- apply(clus_percent,1,function(x) prop.table(x))
clus_percent <- t(clus_percent[order(as.numeric(rownames(clus_percent))),])
clus_percent <- na.omit(clus_percent)
Stim=colnames(clus_percent)
clus_percent=melt(clus_percent, id='Stim')
names(clus_percent)[2]='Stim'
names(clus_percent)[1]='Celltype'
clus_percent$Celltype <- factor(clus_percent$Celltype,levels = c('HSC','MPP #1',"MPP #2",'LMPP #1','LMPP #2'))
clus_percent$Celltype <- factor(clus_percent$Celltype,levels = c('HSC','MPP #1',"MPP #2",'MEP','EryP','MkP'))
clus_percent$Celltype <- factor(clus_percent$Celltype,levels = c('HSC','MPP #1',"MPP #2",'LMPP #2','CLP','pro-B',"pre-B"))
clus_percent$Celltype <- factor(clus_percent$Celltype,levels = c('HSC','MPP #1',"MPP #2",'LMPP #1','ETP','mThy',"cThy #1","cThy #2"))
clus_percent$Celltype <- factor(clus_percent$Celltype,levels = c('HSC','MPP #1',"MPP #2",'LMPP #1','GMP','CDP',"pre-pDC"))
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
        axis.line = element_line(colour = "grey90",size = rel(1)))# + theme_bw()
ggsave("/home/yushiya/cd34/fig/fig-My-clusP.pdf", plot=p1, width=4, height=2.5)

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
type_dis$stim <- sub1$stim
type_dis$pseudotime <- sub1$pseudotime
names(type_dis)[1]='ID'
p1=ggplot(type_dis, aes(x = pseudotime, fill = stim))+ geom_density(alpha=0.7, size=0.3)+
  theme_classic()+
  scale_fill_manual(values=color_stim)+
  scale_x_continuous(limits = c(0,30))+
  theme(axis.text.x=element_text(hjust=0.5),
        text = element_text(size = 16),
        panel.background = element_blank(),
        panel.grid = element_line(colour = "black"),
        axis.line = element_line(colour = "black",size = rel(1),arrow = arrow(angle = 30,length = unit(0.1,"inches"))))
ggsave(file="/home/yushiya/cd34/fig/fig1-density-B.pdf", plot=p1, width=6, height=3)
p1=ggplot(type_dis, aes(x = pseudotime, color=stim))+ geom_density(size=1)+
  theme_classic()+
  scale_color_manual(values=color_stim)+
  scale_x_continuous(limits = c(0,30))+
  theme(axis.text.x=element_text(hjust=0.5),
        text = element_text(size = 16),
        panel.background = element_blank(),
        panel.grid = element_line(colour = "black"),
        axis.line = element_line(colour = "black",size = rel(1),arrow = arrow(angle = 30,length = unit(0.1,"inches"))))
ggsave(file="/home/yushiya/cd34/fig/fig-HSPC-density.pdf", plot=p1, width=6, height=3)

#violin plot
p <- VlnPlot(sub1, features = "pseudotime", pt.size = 0, split.by = "stim",cols=color_stim,combine = FALSE,y.max=30)
p1=p[[1]] + coord_flip() 
ggsave(file="/home/yushiya/data/cd34/manu_fig/v4/fig1-Vlndensity-my.pdf", plot=p1, width=6, height=6)

#calculate DEGs and GOs among stims in each celltype
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(GOSemSim)
library(DOSE)
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

#Find markers
type <- c("MPP #1","MPP #2","LMPP #1","LMPP #2","BMEP","EryP","MkP","CLP","pro-B","GMP")
type <- c("ETP","pre-cDC","pre-pDC")
for (i in 1:length(type)) {
  #sub1 <- sub1 <- subset(immune.combined, subset= celltype==type[i] & (stim=="BM"|stim=="PB"|stim=="mPB"))
  sub1 <- sub1 <- subset(immune.combined, subset= celltype==type[i])
  Idents(sub1) <- sub1$stim
  cp_stim.markers <- FindAllMarkers(sub1, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1)
  write.csv(cp_stim.markers, paste("/home/yushiya/cd34/fig/cpstim_",type[i],".csv",sep=""))
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
  saveRDS(go.p,paste("/home/yushiya/cd34/fig/GO_",type[i],".rds",sep=""))
}
for (i in 1:length(type)) {
  go.p <- readRDS(paste("/home/yushiya/cd34/fig/GO_",type[i],".rds",sep=""))
  go.p1 <- simplify(go.p,cutoff=0.4,by="p.adjust",select_fun=min)  #去除冗余
  p1=ggplot(go.p1, aes(Cluster, Description), showCategory=6) +
    geom_point(aes(color=p.adjust, size=GeneRatio))+
    theme_classic()+
    theme(axis.text.x =element_text(angle = 45,vjust = 0.5,hjust = 0.5),
          text = element_text(size = 15),
          panel.background = element_blank(),
          panel.grid = element_line(colour = "grey90"))
  ggsave(file=paste("/home/yushiya/cd34/fig/GO_",type[i],".pdf",sep=""), plot=p1, width=8, height=7)
}

#stemness and entropy by SCENT CCAT
library(SCENT)
data(net13Jun12)
#data(dataChu)
sub1 <- subset(immune.combined, subset= stim=="PB")
phenoImmune.v <- as.character(sub1$celltype)
#immune.m <- data.frame(sub1@assays$RNA@counts)
mtx <- immune.combined[["RNA"]]@data
immune.m <- immune.combined[["RNA"]]@data  #log normalized data
#immune.m <- log2(immune.m+1)
geneID <- rownames(immune.m)
Entrez_ID <- bitr(geneID, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
temp <- Entrez_ID$ENTREZID[match(geneID,Entrez_ID$SYMBOL)]
tempID <- which(is.na(temp))
immune.m <- immune.m[-tempID,]
temp <- na.omit(temp)
rownames(immune.m) <- temp
immune.m <- data.matrix(immune.m)
ccat.v <- CompCCAT(exp = immune.m, ppiA = net13Jun12.m)
immune.combined$ccat <- ccat.v
p1=VlnPlot(immune.combined, features = c("ccat"), cols=color_all, pt.size = 0) +
  geom_boxplot(width=.2,col="black",fill="white")+ NoLegend()
ggsave(file="/home/yushiya/cd34/fig/fig2-Score-ccat.pdf", plot=p1, width=8, height=3.5)

#Migratory score
library(ggpubr)
positive_regu_selected <- c("MPP1","MSN","CD99","RAC2","MAPK1","APP","PLCB1","ZNF580","TGFB1","CALR",
                            "TPGS1","C1QBP","CRK","PYCARD","MAPK3","SPN","ADAM10","RIN3","NCKAP1L","DNM1L",
                            "CD9","WNK1","VEGFB","DNAJC4","SPI1","MDK","RHOG","CD81","ADAM8",
                            "CAMK1D","ANXA1","DOCK8","LYN","TGS1","HOXA7","GPSM3","RIPOR2","CCL28","PTGER4",
                            "RHOH","CD47","SELENOK","INPP5D","ITGA4","ANXA4","RTN4","ADAM17","GCSAML","MIA3")
immune.combined<-AddModuleScore(object = immune.combined, features = list(Migratory_Score=positive_regu_selected), name="Migratory_Score")
sub1 <- subset(immune.combined, subset= celltype=="HSC"|celltype=="MPP #1"|celltype=="MPP #2"|celltype=="LMPP #1"|celltype=="LMPP #2"|celltype=="BMEP"|celltype=="EryP"|celltype=="MkP"|celltype=="pre-pro B"|celltype=="pro-B"|celltype=="pre-B"|celltype=="GMP"|celltype=="pre-cDC")
type <- c("MPP #1","MPP #2","LMPP #1","LMPP #2","BMEP","EryP","MkP","pre-pro B","pro-B","pre-B","GMP","pre-cDC")
clus="HSC"
for (i in 1:length(type)) {
  clus=type[i]
  sub1 <- subset(immune.combined, subset= celltype==clus)
  sub1 <- subset(sub1, subset= stim=="BM"|stim=="PB"|stim=="mPB")
  miGene <- as.data.frame(colnames(sub1))
  miGene$stim <- sub1$stim
  miGene$celltype <- sub1$celltype
  miGene$stim <- factor(miGene$stim, levels=c("BM","mPB","PB","thymus"))
  names(miGene)[1]='ID'
  miGene$gene='Migratory_Score1'
  miGene$Expression <- sub1$Migratory_Score1
  #plot with different facet
  p1=ggboxplot(miGene, x="stim", y="Expression", color = "stim",
               palette = "nejm", short.panel.labs = T, ncol=7)+
    #stat_compare_means(method = "anova", label.y=0.5, vjust=1)+
    #stat_compare_means(label = "p.signif",method = "t.test", ref.group = ".all.")+
    stat_compare_means(comparisons = list(c("BM","mPB"),c("BM","PB"),c("mPB","PB")),
                       #stat_compare_means(comparisons = list(c("BM","mPB"),c("BM","PB"),c("mPB","PB"),c("BM","thymus"),c("mPB","thymus"),c("PB","thymus")),
                       label = "p.signif",method = "t.test")+
    ggtitle(clus)+
    theme(axis.text.x =element_text(angle = 45,vjust = 0.5,hjust = 0.5),
          plot.title = element_text(hjust = 0.5))
  ggsave(file=paste("/home/yushiya/data/cd34/manu_fig/v4/fig3-migr_",clus,".pdf",sep=""), plot=p1, width=3.5, height=4.5)
}
sub1 <- subset(immune.combined, subset= celltype==clus)
sub1 <- subset(sub1, subset= stim=="BM"|stim=="PB"|stim=="mPB")
miGene <- as.data.frame(colnames(sub1))
miGene$stim <- sub1$stim
miGene$celltype <- sub1$celltype
miGene$stim <- factor(miGene$stim, levels=c("BM","mPB","PB","thymus"))
names(miGene)[1]='ID'
miGene$gene='Migratory_Score1'
miGene$Expression <- sub1$Migratory_Score1
#plot with different facet
ggboxplot(miGene, x="stim", y="Expression", color = "stim",
          palette = "nejm", short.panel.labs = T, ncol=7)+
  #stat_compare_means(method = "anova", label.y=0.5, vjust=1)+
  #stat_compare_means(label = "p.signif",method = "t.test", ref.group = ".all.")+
  stat_compare_means(comparisons = list(c("BM","mPB"),c("BM","PB"),c("mPB","PB")),
                     #stat_compare_means(comparisons = list(c("BM","mPB"),c("BM","PB"),c("mPB","PB"),c("BM","thymus"),c("mPB","thymus"),c("PB","thymus")),
                     label = "p.signif",method = "t.test")+
  #ggtitle(clus)+
  theme(axis.text.x =element_text(angle = 45,vjust = 0.5,hjust = 0.5),
        plot.title = element_text(hjust = 0.5))
ggsave(file=paste("/home/yushiya/data/cd34/manu_fig/v4/fig3-migr_",clus,".pdf",sep=""), plot=p1, width=3.5, height=4.5)
#boxplot with lineages in the same facet
p1=ggboxplot(miGene, x="celltype", y="Expression", color="stim", bxp.errorbar = T,
             add="boxplot", short.panel.labs = TRUE, ncol=4, add.params=list(size=0.2))+
  #stat_compare_means(comparisons = list(c("HSC","MPP #1"),c("MPP #1","MPP #2"),c("MPP #2","BMEP"),c("BMEP","EryP"),c("BMEP","MkP")),
  #stat_compare_means(comparisons = list(c("HSC","MPP #1"),c("MPP #1","MPP #2"),c("MPP #2","LMPP #2"),c("LMPP #2","pre-pro B"),c("LMPP #2","pro-B"),c("pre-pro B","pro-B"),c("pro-B","pre-B")), 
  #stat_compare_means(comparisons = list(c("HSC","MPP #1"),c("MPP #1","MPP #2"),c("MPP #2","LMPP #1"),c("LMPP #1","ETP"),c("ETP","mThy"),c("mThy","cThy #1"),c("mThy","cThy #2"),c("cThy #1","cThy #2")),  
  #stat_compare_means(comparisons = list(c("HSC","MPP #1"),c("MPP #1","MPP #2"),c("MPP #2","LMPP #1"),c("LMPP #1","GMP"),c("GMP","pre-cDC"),c("GMP","pre-pDC"),c("pre-cDC","pre-pDC")),  
  #                   label = "p.signif",method = "t.test")+
  scale_color_manual(values = color_stim) +
  theme(axis.text.x =element_text(angle = 45,hjust = 1),
        plot.title = element_text(hjust = 0))
ggsave(file="/home/yushiya/cd34/fig/fig3-migr-BMPB-lineages.pdf", plot=p1, width=9, height=3)
ggsave(file="/home/yushiya/cd34/fig/fig5-migr-ETP.pdf", plot=p1, width=3, height=3)
#boxplot with lineages in the same facet not split by stim
p1=ggboxplot(miGene, x="celltype", y="Expression", color="celltype", bxp.errorbar = T,
             add="boxplot", short.panel.labs = TRUE, ncol=4, add.params=list(size=0.2))+
  #stat_compare_means(comparisons = list(c("HSC","MPP #1"),c("MPP #1","MPP #2"),c("MPP #2","BMEP"),c("BMEP","EryP"),c("BMEP","MkP")),
  #stat_compare_means(comparisons = list(c("HSC","MPP #1"),c("MPP #1","MPP #2"),c("MPP #2","LMPP #2"),c("LMPP #2","pre-pro B"),c("LMPP #2","pro-B"),c("pre-pro B","pro-B"),c("pro-B","pre-B")), 
  #stat_compare_means(comparisons = list(c("HSC","MPP #1"),c("MPP #1","MPP #2"),c("MPP #2","LMPP #1"),c("LMPP #1","ETP"),c("ETP","mThy"),c("mThy","cThy #1"),c("mThy","cThy #2"),c("cThy #1","cThy #2")),  
  #stat_compare_means(comparisons = list(c("HSC","MPP #1"),c("MPP #1","MPP #2"),c("MPP #2","LMPP #1"),c("LMPP #1","GMP"),c("GMP","pre-cDC"),c("GMP","pre-pDC"),c("pre-cDC","pre-pDC")),  
  #                   label = "p.signif",method = "t.test")+
  scale_color_manual(values = color_BMPB_lineages) +
  theme(axis.text.x =element_text(angle = 45,hjust = 1),
        plot.title = element_text(hjust = 0))
ggsave(file="/home/yushiya/cd34/fig/fig3-migr-BMPB-lineages-all.pdf", plot=p1, width=9, height=4)


#cell proliferation score
#GO:0072089 stem cell proliferation QuickGO-annotations-1680317212768-20230401.tsv
#GO:0008283 cell population proliferation QuickGO-annotations-1680319286203-20230401.tsv
proliferation <- data.table::fread("/home/yushiya/QuickGO-annotations-1680317212768-20230401.tsv")
proliferation <- data.table::fread("/home/yushiya/QuickGO-annotations-1680319286203-20230401.tsv")
proliferation <- proliferation$SYMBOL
proliferation <- immune.combined[["RNA"]]@data[rownames(immune.combined[["RNA"]]@data) %in% proliferation,]
proliferation <- rownames(proliferation)
immune.combined<-AddModuleScore(object = immune.combined, features = list(Stem_Proliferation_Score=proliferation), name="Stem_Proliferation_Score")
FeaturePlot(immune.combined, features = c("Proliferation_Score1"), max.cutoff = 4, cols = c("#000066","yellow"),split.by="stim")
VlnPlot(immune.combined, features = c("Stem_Proliferation_Score1"), cols=color_stim,split.by="stim", pt.size = 0,combine = FALSE)

#HSC multilineage scores
sub1 <- subset(immune.combined, subset= celltype=="HSC")
sub1 <- subset(sub1, subset= stim=="BM"|stim=="PB"|stim=="mPB")
sub1$stim <- factor(sub1$stim, levels=c("BM","mPB","PB"))
Idents(sub1) <- sub1$stim
p1=VlnPlot(sub1, features = c("HSC_Score1"), cols=color_stim, pt.size = 0) +
  geom_boxplot(width=.2,col="black",fill="white")+ NoLegend()
ggsave(file="/home/yushiya/cd34/fig/fig2-HSC_Score-HSC.pdf", plot=p1, width=4, height=3)
p1=VlnPlot(sub1, features = c("EryP_Score1"), cols=color_stim, pt.size = 0) +
  geom_boxplot(width=.2,col="black",fill="white")+ NoLegend()
ggsave(file="/home/yushiya/cd34/fig/fig2-HSC_Score-EryP.pdf", plot=p1, width=4, height=3)
p1=VlnPlot(sub1, features = c("MkP_Score1"), cols=color_stim, pt.size = 0) +
  geom_boxplot(width=.2,col="black",fill="white")+ NoLegend()
ggsave(file="/home/yushiya/cd34/fig/fig2-HSC_Score-MkP.pdf", plot=p1, width=4, height=3)
p1=VlnPlot(sub1, features = c("My_Score1"), cols=color_stim, pt.size = 0) +
  geom_boxplot(width=.2,col="black",fill="white")+ NoLegend()
ggsave(file="/home/yushiya/cd34/fig/fig2-HSC_Score-My.pdf", plot=p1, width=4, height=3)
p1=VlnPlot(sub1, features = c("B_Score1"), cols=color_stim, pt.size = 0) +
  geom_boxplot(width=.2,col="black",fill="white")+ NoLegend()
ggsave(file="/home/yushiya/cd34/fig/fig2-HSC_Score-B.pdf", plot=p1, width=4, height=3)
p1=VlnPlot(sub1, features = c("T_Score1"), cols=color_stim, pt.size = 0) +
  geom_boxplot(width=.2,col="black",fill="white")+ NoLegend()
ggsave(file="/home/yushiya/cd34/fig/fig2-HSC_Score-T.pdf", plot=p1, width=4, height=3)

#HSC score, GMP score
library(readxl)   #read xlsx files
celltype_markers <- read.csv("/home/yushiya/cd34/fig/immune_rpca_celltype_markers.csv",stringsAsFactors = F)
LT_HSC_genes <- read_xlsx("/home/yushiya/cd34/fig/LT_HSC_genes.xlsx",col_names = T)
HSC_gene <- subset(celltype_markers, subset= cluster=="HSC")
HSC_gene <- HSC_gene[HSC_gene$gene %in% LT_HSC_genes$`Gene Symbol`,]
HSC_gene <- HSC_gene %>% top_n(n = 50, wt = avg_log2FC)
immune.combined<-AddModuleScore(object = immune.combined, features = list(HSC_Score=HSC_gene$gene), name="HSC_Score")
p1=VlnPlot(immune.combined, features = c("HSC_Score1"), cols=color_all, pt.size = 0) +
  geom_boxplot(width=.2,col="black",fill="white")+ NoLegend()
ggsave(file="/home/yushiya/cd34/fig/fig2-Score-HSC-1.pdf", plot=p1, width=8, height=3.5)
EryP_gene <- subset(celltype_markers, subset= cluster=="EryP")
EryP_gene <- EryP_gene %>% top_n(n = 50, wt = avg_log2FC)
immune.combined<-AddModuleScore(object = immune.combined, features = list(EryP_Score=EryP_gene$gene), name="EryP_Score")
p1=VlnPlot(immune.combined, features = c("EryP_Score1"), cols=color_all, pt.size = 0) +
  geom_boxplot(width=.2,col="black",fill="white")+ NoLegend()
ggsave(file="/home/yushiya/cd34/fig/fig2-Score-EryP.pdf", plot=p1, width=8, height=3.5)
MkP_gene <- subset(celltype_markers, subset= cluster=="MkP")
MkP_gene <- MkP_gene %>% top_n(n = 50, wt = avg_log2FC)
immune.combined<-AddModuleScore(object = immune.combined, features = list(MkP_Score=MkP_gene$gene), name="MkP_Score")
p1=VlnPlot(immune.combined, features = c("MkP_Score1"), cols=color_all, pt.size = 0) +
  geom_boxplot(width=.2,col="black",fill="white")+ NoLegend()
ggsave(file="/home/yushiya/cd34/fig/fig2-Score-MkP.pdf", plot=p1, width=8, height=3.5)
My_gene <- subset(celltype_markers, subset= cluster=="CDP" | cluster=="pre-pDC")
My_gene <- My_gene %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
immune.combined<-AddModuleScore(object = immune.combined, features = list(My_Score=My_gene$gene), name="My_Score")
p1=VlnPlot(immune.combined, features = c("My_Score1"), cols=color_all, pt.size = 0) +
  geom_boxplot(width=.2,col="black",fill="white")+ NoLegend()
ggsave(file="/home/yushiya/cd34/fig/fig2-Score-My.pdf", plot=p1, width=8, height=3.5)
B_gene <- subset(celltype_markers, subset= cluster=="pro-B"|cluster=="pre-B")
B_gene <- B_gene %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
immune.combined<-AddModuleScore(object = immune.combined, features = list(B_Score=B_gene$gene), name="B_Score")
p1=VlnPlot(immune.combined, features = c("B_Score1"), cols=color_all, pt.size = 0) +
  geom_boxplot(width=.2,col="black",fill="white")+ NoLegend()
ggsave(file="/home/yushiya/cd34/fig/fig2-Score-B.pdf", plot=p1, width=8, height=3.5)
T_gene <- subset(celltype_markers, subset= cluster=="Thy #1" | cluster=="Thy #2"| cluster=="Thy #3")
T_gene <- T_gene %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
immune.combined<-AddModuleScore(object = immune.combined, features = list(T_Score=T_gene$gene), name="T_Score")
p1=VlnPlot(immune.combined, features = c("T_Score1"), cols=color_all, pt.size = 0) +
  geom_boxplot(width=.2,col="black",fill="white")+ NoLegend()
ggsave(file="/home/yushiya/cd34/fig/fig2-Score-T.pdf", plot=p1, width=8, height=3.5)

#compare HSPCs: HSC, MPP #1, MPP #2, LMPP #1, LMPP #2
sub1 <- subset(immune.combined, subset= (celltype=="HSC"|celltype=="MPP #1"|celltype=="MPP #2"|celltype=="LMPP #1"|celltype=="LMPP #2")&(stim=="BM"|stim=="PB"|stim=="mPB"))
color_sub <- c("#550000", "#AA3939", "#CC6F66","grey","grey","grey","#E79492", "#FBD2CE","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey")
p1=DimPlot(immune.combined, reduction = "umap", label = F, cols=color_sub)
ggsave("/home/yushiya/cd34/fig/fig-HSPC-umap.pdf", plot=p1, width=6, height=4)
color_sub1 <- c("#550000", "#AA3939", "#CC6F66","#E79492", "#FBD2CE")
p1=DimPlot(sub1, reduction = "umap", label = F, cols=color_sub1, split.by = "stim", ncol=1)
ggsave("/home/yushiya/cd34/fig/fig-HSPC-split.pdf", plot=p1, width=6, height=10)
score <- data.frame(sub1$stim)
colnames(score) <- "stim"
score$celltype <- sub1$celltype
score$HSC_Score <- sub1$HSC_Score1
score$Ery_Score <- sub1$EryP_Score1
score$Mk_Score <- sub1$MkP_Score1
score$My_Score <- sub1$My_Score1
score$B_Score <- sub1$B_Score1
score$T_Score <- sub1$T_Score1
p1=ggplot(score, aes(x=celltype,y=T_Score, fill=stim))+
  geom_violin(position = position_dodge(0.9), width=1) +
  geom_boxplot(width=.2,position = position_dodge(0.9), outlier.shape = NA) +
  theme_classic()+
  scale_fill_manual(values=color_stim)+
  theme(text = element_text(size = 18),
        axis.text.x=element_text(angle=45,vjust=0.5))
ggsave(file="/home/yushiya/cd34/fig/fig-HSPC-T_Score.pdf", plot=p1, width=7, height=3)
cpstim_HSC <- read.csv("/home/yushiya/cd34/fig/cpstim_HSC.csv",stringsAsFactors = F)
cpstim_MPP1 <- read.csv("/home/yushiya/cd34/fig/cpstim_MPP #1.csv",stringsAsFactors = F)
cpstim_MPP2 <- read.csv("/home/yushiya/cd34/fig/cpstim_MPP #2.csv",stringsAsFactors = F)
cpstim_LMPP1 <- read.csv("/home/yushiya/cd34/fig/cpstim_LMPP #1.csv",stringsAsFactors = F)
cpstim_LMPP2 <- read.csv("/home/yushiya/cd34/fig/cpstim_LMPP #2.csv",stringsAsFactors = F)
common_gene <- cpstim_HSC[(cpstim_HSC$gene %in% cpstim_MPP1$gene),]
common_gene <- common_gene[common_gene$gene %in% cpstim_MPP2$gene,]
common_gene <- common_gene[common_gene$gene %in% cpstim_LMPP1$gene,]
common_gene <- common_gene[common_gene$gene %in% cpstim_LMPP2$gene,]
top20 <- common_gene %>% group_by(cluster) %>% top_n(n=13, wt=avg_log2FC)
top20 <- top20[,"gene"] %>% unique()
top20 <- top20[-"RPL27A",]
sub1$stim_celltype <- paste(sub1$stim,sub1$celltype,sep="_")
cellInfo <- data.frame(stim_celltype=sub1$stim_celltype)
mtx <- data.frame(sub1@assays[["RNA"]]@data) 
colnames(mtx) <- rownames(cellInfo)
top10_exp <- sapply(split(rownames(cellInfo), cellInfo$stim_celltype),
                    function(cells) rowMeans(mtx[top20$gene,cells]))
lev<-c("BM_HSC","BM_MPP #1","BM_MPP #2","BM_LMPP #1","BM_LMPP #2",
       "mPB_HSC","mPB_MPP #1","mPB_MPP #2","mPB_LMPP #1","mPB_LMPP #2",
       "PB_HSC","PB_MPP #1","PB_MPP #2","PB_LMPP #1","PB_LMPP #2")
top10_exp <- top10_exp[,lev] 
annotation_col <- data.frame(Stim = factor(rep(c("BM", "mPB","PB"), c(5,5,5))))
rownames(annotation_col) <- lev
annotation_colors =list(Stim=c(BM="#DB5C25",mPB="#F3B747",PB="#649541"))
p1=pheatmap(top10_exp, cluster_cols = F, #cluster_rows = F,  
            clustering_method = "average",
            annotation_col = annotation_col,
            annotation_colors = annotation_colors,
            scale = "row",
            #border=F,
            col=colorRampPalette(c("navy","white" ,"firebrick3"))(50),
            fontsize_row = 12, fontsize_col= 12, angle_col = 45)
ggsave("/home/yushiya/cd34/fig/fig-HSPC-DEG.pdf", plot=p1, width=7, height=7)
HSC_go <- readRDS("/home/yushiya/cd34/fig/GO_HSC.rds")
MPP1_go <- readRDS("/home/yushiya/cd34/fig/GO_MPP #1.rds")
BMID <- subset(common_gene, subset= cluster=="BM")
BM_gene <- bitr(BMID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
mPBID <- subset(common_gene, subset= cluster=="mPB")
mPB_gene <- bitr(mPBID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
PBID <- subset(common_gene, subset= cluster=="PB")
PB_gene <- bitr(PBID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
cp = list(BM.gene=BM_gene$ENTREZID, mPB.gene=mPB_gene$ENTREZID, PB.gene=PB_gene$ENTREZID)
go.p <- compareCluster(cp,
                       fun = "enrichGO",
                       OrgDb = "org.Hs.eg.db",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.01
)
go.p1 <- simplify(go.p,cutoff=0.8,by="p.adjust",select_fun=min)
p1=ggplot(go.p1, aes(Cluster, Description), showCategory=5) +
  geom_point(aes(color=p.adjust, size=GeneRatio))+
  theme_classic()+
  theme(axis.text.x =element_text(angle = 45,vjust = 0.5,hjust = 0.5),
        text = element_text(size = 15),
        panel.background = element_blank(),
        panel.grid = element_line(colour = "grey90"))#+theme_bw()
ggsave("/home/yushiya/cd34/fig/fig-HSPC-GO.pdf", plot=p1, width=7, height=7)

#compare B lineage: CLP, pro-B, pre-B
sub1 <- subset(immune.combined, subset= (celltype=="CLP"|celltype=="pro-B"|celltype=="pre-B")&(stim=="BM"|stim=="PB"|stim=="mPB"))
sub1 <- subset(sub1, subset= celltype %in% c("CLP","pro-B"))
stim.markers <- FindAllMarkers(sub1, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1)
top20 <- stim.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top20_combined <- rbind(subset(top20, subset = cluster=="CLP")[,"gene"],subset(top20, subset = cluster=="pro-B")[,"gene"]) %>% unique()
sub11<- subset(sub1, subset= celltype=="CLP")
sub12<- subset(sub1, subset= celltype=="pro-B")
mean_exp <- apply(sub11@assays[["RNA"]]@scale.data[top20_combined$gene,],1,mean)
top20_combined$CLP <- mean_exp
mean_exp <- apply(sub12@assays[["RNA"]]@scale.data[top20_combined$gene,],1,mean)
top20_combined$proB <- mean_exp
top20_combined<-data.frame(top20_combined)
rownames(top20_combined) <- top20_combined$gene
top20_combined <- top20_combined[,-1]
p1=pheatmap(top20_combined, cluster_cols = F, #cluster_rows = F,  
            #clustering_method = "median",
            fontsize_row = 7, fontsize_col= 10, angle_col = 45,
            col=colorRampPalette(c("#000066","#339999" ,"yellow"))(50))
ggsave("/home/yushiya/cd34/fig/fig-B-DEG.pdf", plot=p1, width=3, height=4)
BMID <- subset(stim.markers, subset= cluster=="CLP")
BM_gene <- bitr(BMID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
mPBID <- subset(stim.markers, subset= cluster=="pro-B")
mPB_gene <- bitr(mPBID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
cp = list(CLP=BM_gene$ENTREZID, proB=mPB_gene$ENTREZID)
go.p <- compareCluster(cp,
                       fun = "enrichGO",
                       OrgDb = "org.Hs.eg.db",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.01
)
go.p1 <- simplify(go.p,cutoff=0.8,by="p.adjust",select_fun=min)
p1=ggplot(go.p1, aes(Cluster, Description), showCategory=5) +
  geom_point(aes(color=p.adjust, size=GeneRatio))+
  theme_classic()+
  theme(axis.text.x =element_text(angle = 45,vjust = 0.5,hjust = 0.5),
        text = element_text(size = 15),
        panel.background = element_blank(),
        panel.grid = element_line(colour = "grey90"))#+theme_bw()
ggsave("/home/yushiya/cd34/fig/fig-B-GO.pdf", plot=p1, width=7, height=5)

#compare EryMk lineage: MEP, EryP, MkP
sub1 <- subset(immune.combined, subset= (celltype=="MEP"|celltype=="EryP"|celltype=="MkP")&(stim=="BM"|stim=="PB"|stim=="mPB"))
color_sub <- c("grey","grey","grey","#FFBD54","#B97802","#915900","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey")
p1=DimPlot(immune.combined, reduction = "umap", label = F, cols=color_sub)
ggsave("/home/yushiya/cd34/fig/fig-EryMk-umap.pdf", plot=p1, width=6, height=4)
cpstim_MEP <- read.csv("/home/yushiya/cd34/fig/cpstim_BMEP.csv",stringsAsFactors = F)
cpstim_Ery <- read.csv("/home/yushiya/cd34/fig/cpstim_EryP.csv",stringsAsFactors = F)
cpstim_Mk <- read.csv("/home/yushiya/cd34/fig/cpstim_MkP.csv",stringsAsFactors = F)
common_gene <- cpstim_MEP[(cpstim_MEP$gene %in% cpstim_Ery$gene),]
common_gene <- common_gene[common_gene$gene %in% cpstim_Mk$gene,]
top20 <- common_gene %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
top20 <- top20[,"gene"] %>% unique()
sub1$stim_celltype <- paste(sub1$stim,sub1$celltype,sep="_")
cellInfo <- data.frame(stim_celltype=sub1$stim_celltype)
mtx <- data.frame(sub1@assays[["RNA"]]@data[top20$gene,]) 
colnames(mtx) <- rownames(cellInfo)
top10_exp <- sapply(split(rownames(cellInfo), cellInfo$stim_celltype),
                    function(cells) rowMeans(mtx[top20$gene,cells]))
lev<-c("BM_MEP","BM_EryP","BM_MkP",
       "mPB_MEP","mPB_EryP","mPB_MkP",
       "PB_MEP","PB_EryP","PB_MkP")
top10_exp <- top10_exp[,lev] 
annotation_col <- data.frame(Stim = factor(rep(c("BM", "mPB","PB"), c(3,3,3))))
rownames(annotation_col) <- lev
annotation_colors =list(Stim=c(BM="#DB5C25",mPB="#F3B747",PB="#649541"))
p1=pheatmap(top10_exp, cluster_cols = F, #cluster_rows = F,  
            clustering_method = "average",
            annotation_col = annotation_col,
            annotation_colors = annotation_colors,
            scale = "row",
            #border=F,
            col=colorRampPalette(c("navy","white" ,"firebrick3"))(50),
            fontsize_row = 12, fontsize_col= 12, angle_col = 45)
ggsave("/home/yushiya/cd34/fig/fig-EryMk-DEG.pdf", plot=p1, width=7, height=7)
BMID <- subset(common_gene, subset= cluster=="BM")
BM_gene <- bitr(BMID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
mPBID <- subset(common_gene, subset= cluster=="mPB")
mPB_gene <- bitr(mPBID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
PBID <- subset(common_gene, subset= cluster=="PB")
PB_gene <- bitr(PBID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
cp = list(BM.gene=BM_gene$ENTREZID, mPB.gene=mPB_gene$ENTREZID, PB.gene=PB_gene$ENTREZID)
go.p <- compareCluster(cp,
                       fun = "enrichGO",
                       OrgDb = "org.Hs.eg.db",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.01
)
go.p1 <- simplify(go.p,cutoff=0.4,by="p.adjust",select_fun=min)
p1=ggplot(go.p1, aes(Cluster, Description), showCategory=5) +
  geom_point(aes(color=p.adjust, size=GeneRatio))+
  theme_classic()+
  theme(axis.text.x =element_text(angle = 45,vjust = 0.5,hjust = 0.5),
        text = element_text(size = 15),
        panel.background = element_blank(),
        panel.grid = element_line(colour = "grey90"))#+theme_bw()
ggsave("/home/yushiya/cd34/fig/fig-EryMk-GO.pdf", plot=p1, width=7, height=7)

#compare My lineage: GMP,CDP,pre-pDC
sub1 <- subset(immune.combined, subset= celltype=="GMP"|celltype=="CDP"|celltype=="pre-pDC")
color_sub <- c("grey","grey","grey","grey","grey","grey","grey","grey","#FFF176","#CECA68","#B2A13F","grey","grey","grey","grey","grey","grey","grey")
p1=DimPlot(immune.combined, reduction = "umap", label = F, cols=color_sub)
ggsave("/home/yushiya/cd34/fig/fig-My-umap.pdf", plot=p1, width=6, height=4)
cp_stim.markers <- read.csv("/home/yushiya/cd34/fig/cpstim_pre-pDC.csv",stringsAsFactors = F)
top_genes <- cp_stim.markers %>% group_by(cluster) %>% top_n(n = 6, wt = avg_log2FC)
top_genes <- top_genes[,"gene"] %>% unique()
colnames(top_genes) <- "gene"
sub11<- subset(sub1, subset= stim=="BM")
sub12<- subset(sub1, subset= stim=="mPB")
sub13<- subset(sub1, subset= stim=="PB")
sub14<- subset(sub1, subset= stim=="thymus")
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
            scale = "row",
            fontsize_row = 10, fontsize_col= 12, angle_col = 0,
            col=colorRampPalette(c("#000066","#339999" ,"yellow"))(50))
ggsave("/home/yushiya/cd34/fig/fig-My-DEG-prepDC.pdf", plot=p1, width=4, height=5)
BMID <- subset(cp_stim.markers, subset= cluster=="BM")
BM_gene <- bitr(BMID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
mPBID <- subset(cp_stim.markers, subset= cluster=="mPB")
mPB_gene <- bitr(mPBID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
PBID <- subset(cp_stim.markers, subset= cluster=="PB")
PB_gene <- bitr(PBID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
thymusID <- subset(cp_stim.markers, subset= cluster=="thymus")
thymus_gene <- bitr(thymusID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
cp = list(BM.gene=BM_gene$ENTREZID, mPB.gene=mPB_gene$ENTREZID, PB.gene=PB_gene$ENTREZID, thymus.gene=thymus_gene$ENTREZID)
go.p <- compareCluster(cp,
                       fun = "enrichGO",
                       OrgDb = "org.Hs.eg.db",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.01
)
go.p1 <- simplify(go.p,cutoff=0.8,by="p.adjust",select_fun=min)
p1=ggplot(go.p1, aes(Cluster, Description), showCategory=5) +
  geom_point(aes(color=p.adjust, size=GeneRatio))+
  theme_classic()+
  theme(axis.text.x =element_text(angle = 45,vjust = 0.5,hjust = 0.5),
        text = element_text(size = 15),
        panel.background = element_blank(),
        panel.grid = element_line(colour = "grey90"))#+theme_bw()
ggsave("/home/yushiya/cd34/fig/fig-My-GO-prepDC.pdf", plot=p1, width=9, height=7)
sub1 <- subset(sub1, subset= celltype=="pre-pDC")
sub1$stim <- factor(sub1$stim, levels=c("BM","mPB","PB","thymus"))
Idents(sub1) <- sub1$stim
p1=VlnPlot(sub1, features = c("My_Score1"), cols=color_stim, pt.size = 0) +
  geom_boxplot(width=.2,col="black",fill="white")+ NoLegend()
ggsave(file="/home/yushiya/cd34/fig/fig-My-Score-My.pdf", plot=p1, width=4, height=3)
p1=VlnPlot(sub1, features = c("B_Score1"), cols=color_stim, pt.size = 0) +
  geom_boxplot(width=.2,col="black",fill="white")+ NoLegend()
ggsave(file="/home/yushiya/cd34/fig/fig-My-Score-B.pdf", plot=p1, width=4, height=3)
p1=VlnPlot(sub1, features = c("T_Score1"), cols=color_stim, pt.size = 0) +
  geom_boxplot(width=.2,col="black",fill="white")+ NoLegend()
ggsave(file="/home/yushiya/cd34/fig/fig-My-Score-T.pdf", plot=p1, width=4, height=3)



#chemokines
sub1 <- subset(immune.combined, subset= stim=="BM"|stim=="PB"|stim=="mPB")
Idents(sub1) <- sub1$stim
cp_stim.markers <- FindAllMarkers(sub1, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1)
write.csv(cp_stim.markers, file=paste("/home/yushiya/data/cd34/manu_fig/v4/fig2-DEG-BMPB.csv",sep=""))
p1=FeaturePlot(sub1, features = c("CXCR4"), max.cutoff = 3, cols = c("#C9E9FF","#981C12"), split.by="stim",order=T, pt.size = 0.001) + 
  theme(legend.position = "right")
ggsave(file="/home/yushiya/cd34/fig/fig2-CXCR4-fea.pdf", plot=p1, width=10, height=3)
p1=FeaturePlot(sub1, features = c("RPL31"), max.cutoff = 3, cols = c("#C9E9FF","#981C12"), split.by="stim") + 
  theme(legend.position = "right")
ggsave(file="/home/yushiya/cd34/fig/fig2-RPL31-fea.pdf", plot=p1, width=10, height=3)
p1=VlnPlot(sub1, features = "CXCR4", pt.size = 0,combine = FALSE, split.by="stim", cols=color_stim)
ggsave(file="/home/yushiya/data/cd34/manu_fig/v4/fig3-CXCXR4-vln.pdf", plot=p1[[1]], width=8, height=2.5)
p1=FeaturePlot(immune.combined, features = c("CCR7"), max.cutoff = 3, cols = c("#C9E9FF","#981C12"), split.by="stim", order=T) + 
  theme(legend.position = "right")
ggsave(file="/home/yushiya/cd34/fig/fig3-CCR7-fea-1.pdf", plot=p1, width=12, height=3)
p1=FeaturePlot(immune.combined, features = c("CCR9"), max.cutoff = 3, cols = c("#C9E9FF","#981C12"), split.by="stim", order=T) + 
  theme(legend.position = "right")
ggsave(file="/home/yushiya/cd34/fig/fig3-CCR9-fea-1.pdf", plot=p1, width=12, height=3)


#compare pre-pro B and pro-B
sub1 <- subset(immune.combined, subset= celltype=="pre-pro B"|celltype=="pro-B")
cp_stim.markers <- FindAllMarkers(sub1, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1)
write.csv(cp_stim.markers, file=paste("/home/yushiya/data/cd34/manu_fig/v4/fig2-compare-Bpros.csv",sep=""))
VlnPlot(sub1, features = "CXCR4", pt.size = 0,combine = FALSE, cols=color_stim)
miGene <- as.data.frame(colnames(sub1))
miGene$stim <- sub1$stim
miGene$celltype <- sub1$celltype
names(miGene)[1]='ID'
miGene$gene="CXCR4"
miGene$Expression <- sub1@assays[["integrated"]]@scale.data["CXCR4",]
p1=ggboxplot(miGene, x="celltype", y="Expression", fill="celltype", bxp.errorbar = T,
             palette = "jco", add="boxplot", short.panel.labs = TRUE, ncol=4, add.params=list(size=0.2))+
  stat_compare_means(comparisons = list(c("pre-pro B","pro-B")),
                     label = "p.signif",method = "t.test")+
  theme(axis.text.x =element_text(angle = 45,hjust = 1),
        plot.title = element_text(hjust = 0))
ggsave(file="/home/yushiya/data/cd34/manu_fig/v4/fig1-compareB-CXCR4.pdf", plot=p1, width=4, height=4)

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
#saveRDS(sub_sc, file = "/home/yushiya/cd34/data/stromal_sub.rds")
sc <- readRDS("/home/yushiya/cd34/data/stromal.rds")
sub_sc <- readRDS(file = "/home/yushiya/cd34/data/stromal_sub.rds")

#cellchat
sub1 <- subset(immune.combined, subset= celltype=="HSC"|celltype=="MPP #1"|celltype=="MPP #2"|celltype=="MEP"|celltype=="EryP"|celltype=="MkP"|celltype=="LMPP #1"|celltype=="LMPP #2"|celltype=="GMP"|celltype=="CDP"|celltype=="pre-pDC"|celltype=="CLP"|celltype=="pro-B"|celltype=="pre-B")
sub1 <- subset(sub1, subset= stim=="PB" | stim=="BM" | stim=="mPB")
#devtools::install_github("sqjin/CellChat")
library(CellChat)
library(svglite)
Idents(sub1) <- sub1$stim_celltype
sum <- merge(sub_sc,sub1)
sum <- merge(sum,sub_tec)
cellchat <- createCellChat(sum)
cellchat@idents <- factor(cellchat@idents, levels=c("BM_SC","BM_Endo","HSC","MPP #1","MPP #2","MEP","EryP","MkP","LMPP #1","LMPP #2","GMP","CDP","pre-pDC","ETP","Thy #1","Thy #2","Thy #3","CLP","pro-B","pre-B"))
cellchat@idents <- factor(cellchat@idents, levels=c("BM_SC","BM_Endo","HSC","MPP #1","MPP #2","MEP","EryP","MkP","LMPP #1","LMPP #2","GMP","CDP","pre-pDC","CLP","pro-B","pre-B"))
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
mat2[c(3:5),c(6:9)] <- mat[c(3:5),c(6:9)]
mat2[c(1:2),c(6:9)] <- mat[c(1:2),c(6:9)]
mat2[c(1:2),c(3:16)] <- mat[c(1:2),c(3:16)]
netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = "Interaction weights/strength")
ggsave("/home/yushiya/cd34/fig/fig5-cellchat-all-circle.pdf", plot=p1, width=6, height=6)
# 显示从某些细胞组到其他细胞组的所有显著的相互作用（L-R 对）
netVisual_bubble(cellchat, sources.use = c(1:2), targets.use = c(3:17), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = c(1:2), targets.use = c(6:9), remove.isolate = FALSE)
p1=netVisual_bubble(cellchat, sources.use = c(1), targets.use = c(6:9),signaling = c("SELPLG","PTN","PTPRM","PECAM1","MIF","CD99","CXCL","ANGPTL"), remove.isolate = TRUE, sort.by.target = T)
p1=netVisual_bubble(cellchat, sources.use = c(4), targets.use = c(6:9), remove.isolate = FALSE,thresh=0.01)
p1=netVisual_bubble(cellchat, sources.use = c(4), targets.use = c(6,8,9),signaling = c("PTN","CD99"), remove.isolate = TRUE)
p1=netVisual_bubble(cellchat, sources.use = c(4), targets.use = c(6,8,9),signaling = c("PTN","CD99","MIF","CXCL","MIF","NOTCH"), remove.isolate = TRUE)
p1=netVisual_bubble(cellchat, sources.use = c(1:2), targets.use = c(3:17),signaling = c("CXCL","SELPLG","PECAM1","ADGRE5"), remove.isolate = TRUE)
p1=netVisual_bubble(cellchat, sources.use = c(1), targets.use = c(3:17),signaling = c("CXCL","FN1","APP"), remove.isolate = TRUE)
ggsave("/home/yushiya/cd34/fig/fig3-cellchat-LRpairs-BM_SC-1.pdf", plot=p1, width=5, height=2.5)
ggsave("/home/yushiya/cd34/fig/fig5-cellchat-LRpairs-Thymic_Mesen.pdf", plot=p1, width=4.5, height=3.5)
ggsave("/home/yushiya/cd34/fig/fig3-cellchat-LRpairs-mPB.pdf", plot=p1, width=8, height=10)
p1=plotGeneExpression(cellchat, signaling = c("CXCL","APP","SELPLG","ANGPTL"),color.use = c("#550000","#AA3939","#CC6F66","#E79492","#FBD2CE","#1F78B4","#FDBF6F"))
p1=plotGeneExpression(cellchat, signaling = c("PTN","CD99","MIF","CXCL","NOTCH","CCL"),color.use = c("#1F78B4","#FDBF6F","#66C2A5","#FC8D62","#8DA0CB","#DB5C25","#F3B747","#649541","#4C82C5"))
ggsave("/home/yushiya/cd34/fig/fig-cellchat-geneExp-3.pdf", plot=p1, width=5, height=6)
saveRDS(cellchat, file = "/home/yushiya/cd34/fig/cellchat_all_BM.rds")
#"CXCL","APP","SELPLG","ANGPTL","ADGRE5"
pathways.show <- c("CXCL") 
netVisual_aggregate(cellchat, signaling = pathways.show, sources.use = c(1,2), targets.use = c(3:17))
netVisual_aggregate(cellchat, signaling = pathways.show, sources.use = c(1:5), targets.use = c(6:9))
netVisual_chord_cell(cellchat, signaling = pathways.show, title.name = paste0(pathways.show, " signaling network"))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds", sources.use = c(1,2), targets.use = c(3:16))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds", sources.use = c(1:5), targets.use = c(6:9))
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # “netP”:推断出的信号通路的细胞间通信网络
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 18, height = 8, font.size = 14)
netAnalysis_contribution(cellchat, signaling = pathways.show)
ggsave("/home/yushiya/cd34/fig/fig3-cellchat-CXCL-PB-circle.pdf", plot=p1, width=6, height=6)
ggsave("/home/yushiya/cd34/fig/fig3-cellchat-CXCL-mPB-heatmap.pdf", plot=p1, width=6, height=5)
ggsave("/home/yushiya/cd34/fig/fig3-cellchat-SELPLG-mPB-sig.pdf", plot=p1, width=10, height=5)
ggsave("/home/yushiya/cd34/fig/fig3-cellchat-SELPLG-mPB-ctb.pdf", plot=p1, width=3, height=4)
cellchat <- readRDS(file = "/home/yushiya/cd34/fig/cellchat_mPB.rds")
#saveRDS(cellchat, file = "/home/yushiya/cd34/fig/cellchat_ETP-1.rds")

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


#Analyze ETP
saveRDS(sub1, "/home/yushiya/cd34/fig/ETP.rds")
sub1 <- readRDS("/home/yushiya/cd34/fig/ETP.rds")
sub1 <- subset(immune.combined, subset= celltype=="ETP")
sub1$stim <- factor(sub1$stim,levels = c('BM','mPB','PB','thymus'))
Idents(sub1) <- sub1$stim
sub1 <- RenameIdents(sub1, "BM"="TSP_BM", "mPB"="TSP_mPB", "PB"="TSP_PB","thymus"="ETP")
sub1$celltype1 <- Idents(sub1)
sub1$stim_celltype <- paste(sub1$stim,sub1$celltype,sep="_")
DefaultAssay(sub1) <- "integrated"
sub1 <- RunUMAP(sub1, reduction = "pca", dims = 1:25)
sub1 <- FindNeighbors(sub1, reduction = "pca", dims = 1:25)
sub1 <- FindClusters(sub1, resolution = 0.5)
DimPlot(sub1, reduction = "umap", label = TRUE, cols=brewer.pal(6,"Paired"))
p1=DimPlot(sub1, reduction = "umap", group.by="stim", cols=color_stim)
ggsave("/home/yushiya/cd34/fig/fig4-ETP-umap.pdf", plot=p1, width=6, height=4)
DefaultAssay(sub1) <- "RNA"
sub1 <- ScaleData(sub1, verbose = FALSE)
table(sub1@meta.data[["stim"]],sub1@meta.data[["seurat_clusters"]])
markers1 <- FindAllMarkers(sub1, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1)
sub1 <- RenameIdents(sub1,'1'="c1",'4'="c2",'3'="c3",'0'="c4",'2'="c5",'5'="c6")
FeaturePlot(sub1, features = c("pseudotime"), min.cutoff = 20,max.cutoff = 25,cols = c("blue", "red"),keep.scale=NULL)
FeaturePlot(sub1, features = c("CD44","CD7","CD2","IL7R"), cols = c("#C9E9FF", "red"), order=T,keep.scale=NULL)
ggsave("/home/yushiya/cd34/fig/fig4-ETP-pseu.pdf", plot=p1, width=6, height=4)
VlnPlot(sub1, features = "pseudotime", pt.size = 0, cols=brewer.pal(6,"Paired"))
VlnPlot(sub1, features = "pseudotime", pt.size = 0, cols=color_stim)
ETP_gene <- c("AVP","CRHBP","NKAIN2","AREG",     #HSPC
              "MPO","	S100A10","IRF8","AZU1",
              "GATA2","PF4",
              "VPREB1","EBF1","MS4A1","LEF1","DNTT",
              "CD7","CD3D","CD3E","TCF7")
DoHeatmap(sub1, features = ETP_gene, 
          group.colors=color_stim) + 
  scale_fill_gradientn(colors=c("blue","white","firebrick3"))
ggsave("/home/yushiya/cd34/fig/fig4-ETP-heatmap.pdf", plot=p1, width=6, height=4)
VlnPlot(sub1, features = "My_Score1", pt.size = 0, cols=color_stim)
#monocle 2
library(monocle)
library(igraph)
library(ggsci)
library(viridis)
library(pheatmap)
library(grid)
expr_matrix <- as(as.matrix(sub1@assays$RNA@counts), 'sparseMatrix')
p_data <- sub1@meta.data
p_data$celltype <- sub1$stim
f_data <- data.frame(gene_short_name = row.names(sub1),row.names = row.names(sub1))
pd <- new('AnnotatedDataFrame', data = p_data)
fd <- new('AnnotatedDataFrame', data = f_data)
cds <- newCellDataSet(expr_matrix, phenoData = pd,featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
#cluster 高变基因
deg.cluster <- read.csv("/home/yushiya/cd34/fig/cpstim_ETP.csv",stringsAsFactors = F)
express_genes <- subset(deg.cluster,p_val_adj<0.01)
express_genes <- express_genes %>% group_by(cluster) %>% top_n(n = -500, wt = p_val)
express_genes <- express_genes$gene
cds <- setOrderingFilter(cds,express_genes)
#“dpFeature”选择高变基因
sub1 <- FindVariableFeatures(object = sub1)
expressed_genes<- VariableFeatures(sub1)
cds <- detectGenes(cds, min_expr = 0.1) 
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 10)) 
diff <-differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr="~stim",cores=30)
deg <- subset(diff, qval < 0.01)
deg <-deg[order(deg$qval,decreasing=F),]
ordergene <-row.names(deg)[order(deg$qval)][1:2000]
#ordergene <- rownames(deg)
cds <- setOrderingFilter(cds,ordergene) 
plot_ordering_genes(cds)
#order cells
cds <- reduceDimension(cds,max_components = 2,method = 'DDRTree')
source("/home/yushiya/software/order_cells.R")
cds <- orderCells(cds)
cds <- orderCells(cds, root_state = 3)
p1=plot_cell_trajectory(cds,color_by="Pseudotime",size=1,show_backbone=TRUE)
ggsave("/home/yushiya/cd34/fig/fig4-ETP-pseu_mo2_pseu.pdf", plot=p1, width=5, height=4)
p1=plot_cell_trajectory(cds,color_by="stim",size=1,show_backbone=TRUE)+ scale_color_manual(values = color_stim)
ggsave("/home/yushiya/cd34/fig/fig4-ETP-pseu_mo2.pdf", plot=p1, width=5, height=4)
plot_cell_trajectory(cds,color_by="State",size=1,show_backbone=TRUE)
plot_cell_trajectory(cds, color_by = "stim") + facet_wrap("~stim", nrow = 1)
plot_cell_trajectory(cds, color_by =c("CD44","CD7","CD2","IL7R")) +scale_color_gsea()
colnames(pData(cds))
pData(cds)$CD44 = log2(exprs(cds)['CD44',]+1)
p1=plot_cell_trajectory(cds, color_by ="CD44") + scale_color_gsea()
pData(cds)$CD7 =log2(exprs(cds)['CD7',]+1)
p2=plot_cell_trajectory(cds, color_by ="CD7") +scale_color_gsea()
pData(cds)$CD2 =log2(exprs(cds)['CD2',]+1)
p3=plot_cell_trajectory(cds, color_by ="CD2") +scale_color_gsea()
pData(cds)$IL7R =log2(exprs(cds)['IL7R',]+1)
p4=plot_cell_trajectory(cds, color_by ="IL7R") +scale_color_gsea()
library(patchwork)
p=p1+p2+p3+p4
ggsave("/home/yushiya/cd34/fig/fig4-ETP-pseu_mo2_markers.pdf", plot=p, width=6, height=6)
s.genes <-c("IL7R","CD84")
plot_genes_in_pseudotime(cds[s.genes,], color_by = "stim")
sub1$pseu_mo2 <- cds$Pseudotime
save(cds,file="/home/yushiya/cd34/fig/ETP.cds.monocle2.RData")
#DEG
load("/home/yushiya/cd34/fig/ETP.cds.monocle2.RData")
etp_monocle <- cds
peu_gene <- differentialGeneTest(etp_monocle,fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = 20)
write.csv(peu_gene,file='/home/yushiya/cd34/fig/ETP_peu_gene.csv')
peu_gene <- read.csv(file='/home/yushiya/cd34/fig/ETP_peu_gene.csv', row.names = 1)
peu_gene <- peu_gene[which(peu_gene$qval<0.01 & peu_gene$num_cells_expressed>100),]
peu_gene %>% arrange(qval)  -> peu_gene#按照qval排个序
peu_gene <- peu_gene[1:500,] #这里我们取前100个基因演示
source('/home/yushiya/software/monocle2_heatmap.R')
p11 <- plot_pseudotime_heatmap(etp_monocle[peu_gene$gene_short_name,],
                        num_clusters = 3,
                        cores = 20, 
                        show_rownames = T,return_heatmap =T,
                        hmcols = viridis(256),
                        use_gene_short_name = T)
###首先提取热图中各个module的基因
module_gene <- as.data.frame(cutree(p11$tree_row, k=3))
colnames(module_gene) <- "Module"
module_gene$gene <- rownames(module_gene)
Module_GO=data.frame()
for (i in unique(module_gene$Module)) {
  data=filter(module_gene,module_gene$Module==i)
  go <- enrichGO(gene= data$gene,
                 OrgDb= org.Hs.eg.db,
                 keyType= 'SYMBOL',
                 ont= "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff= 0.05,
                 qvalueCutoff= 0.05,
                 readable= TRUE)
  go_res=go@result
  
  if (dim(go_res)[1] != 0) {
    go_res$cluster=i
    Module_GO=rbind(Module_GO,go_res)
  }
}
Module_GO <- Module_GO[which(Module_GO$qvalue <= 0.05),]
write.csv(Module_GO, file = '/home/yushiya/cd34/fig/ETP_Module_GO.csv')
top_Module_GO <- Module_GO %>% group_by(cluster) %>% top_n(n = -10, wt = p.adjust)
gene <- c("IL7R", "CAMK4","DOCK2","CDK6","CD3D","CD3E","CD3G","BCL11B","LCK","CD2","SOX4","TCF7",
          "CD74","HLA-A","HLA-DRA","HLA-DRB1","B2M",
          "RPL11","RPS21","RPS14","RPL7")
source('/home/yushiya/software/add.flag.R')
p <- add.flag(p11,kept.labels = gene,repel.degree = 0.2)
ggsave("/home/yushiya/cd34/fig/fig4-ETP-pseu_heatmap.pdf", plot=p, width=5, height=6)




#monocle 3
dyn.load("/home/wangjl/local/gdal-2.4.4/libgdal.so.20.5.4")
library(monocle3)
data <- GetAssayData(sub1, assay = 'RNA', slot = 'counts')
cell_metadata <- sub1@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
save(cds,file="/home/yushiya/cd34/fig/ETP.cds.RData")
load(file="/home/yushiya/cd34/fig/ETP.cds.RData")
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, preprocess_method = "PCA")
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('cds.umap')
##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(sub1, reduction = "umap")
cds@int_colData$reducedDims$UMAP <- int.embed
## Monocle3聚类分区
cds <- cluster_cells(cds)
cds@clusters@listData[["UMAP"]][["clusters"]] <- sub1$stim
## 识别轨迹
cds <- learn_graph(cds)
p<-plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
              label_branch_points = FALSE)
cds <- order_cells(cds)
p1=plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = FALSE)
ggsave("/home/yushiya/cd34/fig/fig4-ETP-pseu_new.pdf", plot=p1, width=6, height=4)
p1=plot_cells(cds, genes=c("CD44","CD7","CD2","IL7R"),
              show_trajectory_graph=T,
              label_cell_groups=FALSE,
              label_leaves=FALSE, label_branch_points=F)
ggsave("/home/yushiya/cd34/fig/fig4-ETP-pseu_new_markers.pdf", plot=p1, width=6, height=4)
sub1$pseudotime_new <- pseudotime(cds)
#Expression via pseudotime
pseu_exp<-as.data.frame(colnames(sub1))
pseu_exp$stim <- sub1$stim
pseu_exp$celltype <- sub1$celltype1
pseu_exp$pseudotime <- sub1$pseu_mo2
names(pseu_exp)[1]='ID'
pseu_exp$Migratory_Score <- sub1$Migratory_Score1
pseu_exp$HSC_Score <- sub1$HSC_Score1
pseu_exp$B_Score <- sub1$B_Score1
pseu_exp$T_Score <- sub1$T_Score1
pseu_exp$My_Score <- sub1$My_Score1
p4=ggplot(pseu_exp, aes(x=pseudotime,y=T_Score))+
  geom_point(size=0.3,aes(colour=celltype))+
  stat_smooth(mapping = aes(colour=celltype))+
  #facet_wrap(~celltype)+#theme_bw()
  #scale_x_continuous(limits = c(20,25))+
  scale_color_manual(values = color_stim) +
  theme_classic()+
  theme(axis.text.x=element_text(hjust=0.5),
        text = element_text(size = 18),
        panel.background = element_blank(),
        panel.grid = element_line(colour = "grey90"),
        axis.line = element_line(colour = "grey60",size = rel(1),arrow = arrow(angle = 30,length = unit(0.1,"inches"))))
library(patchwork)
p=p1+p2+p3+p4
ggsave("/home/yushiya/cd34/fig/fig4-ETP-scores.pdf", plot=p, width=10, height=5)
#Expression in vlnplot splited by stim
VlnPlot(sub1, features = "T_Score1", pt.size = 0, split.by = "stim",cols=color_stim)
#DEGs of TFs and surface markers of ETP
library(pheatmap)
ETP_markers <- read.csv("/home/yushiya/cd34/fig/cpstim_ETP.csv",stringsAsFactors = F)
cell_surface <- read.csv("/home/yushiya/data/cd34/cell_surface_protein.csv")
cs.ETP_markers <- ETP_markers[ETP_markers$gene %in% cell_surface$ENTREZ.gene.symbol,]
TF <- read.table("/home/yushiya/data/cd34/TF.txt", sep="\t", header=T)
tf.ETP_markers <- ETP_markers[ETP_markers$gene %in% TF$Symbol,]
selected_ETP.markers <- ETP_markers[(ETP_markers$gene %in% TF$Symbol | ETP_markers$gene %in% cell_surface$ENTREZ.gene.symbol),]
BM_ETP_gene <- subset(selected_ETP.markers, subset=cluster=="BM")
PB_ETP_gene <- subset(selected_ETP.markers, subset=cluster=="PB")
mPB_ETP_gene <- subset(selected_ETP.markers, subset=cluster=="mPB")
thymus_ETP_gene <- subset(selected_ETP.markers, subset=cluster=="thymus")
DEG_ETP <- rbind(BM_ETP_gene[1:5,],mPB_ETP_gene[1:5,],PB_ETP_gene[1:5,],thymus_ETP_gene[1:15,],selected_ETP.markers[selected_ETP.markers$gene %in% c("CD44","CD2"),])
DEG_ETP <- data.frame(DEG_ETP$gene %>% unique())
colnames(DEG_ETP)[1] <- "gene"
sub11<- subset(sub1, subset= stim=="BM")
sub12<- subset(sub1, subset= stim=="mPB")
sub13<- subset(sub1, subset= stim=="PB")
sub14<- subset(sub1, subset= stim=="thymus")
mean_exp <- apply(sub11@assays[["RNA"]]@scale.data[DEG_ETP$gene,],1,mean)
DEG_ETP$BM <- mean_exp
mean_exp <- apply(sub12@assays[["RNA"]]@scale.data[DEG_ETP$gene,],1,mean)
DEG_ETP$mPB <- mean_exp
mean_exp <- apply(sub13@assays[["RNA"]]@scale.data[DEG_ETP$gene,],1,mean)
DEG_ETP$PB <- mean_exp
mean_exp <- apply(sub14@assays[["RNA"]]@scale.data[DEG_ETP$gene,],1,mean)
DEG_ETP$thymus <- mean_exp
rownames(DEG_ETP) <- DEG_ETP$gene
DEG_ETP <- DEG_ETP[,-1]
p1=pheatmap(DEG_ETP, cluster_cols = F, #cluster_rows = F,  
            #clustering_method = "median",
            scale = "row",
            fontsize_row = 10, fontsize_col= 12, angle_col = 0,
            col=colorRampPalette(c("#000066","#339999" ,"yellow"))(50))
ggsave("/home/yushiya/cd34/fig/fig4-ETP-DEG-selected.pdf", plot=p1, width=3.5, height=5)
#TSP and ETP markers
Idents(sub1) <- sub1$stim
sub1 <- RenameIdents(sub1, "BM"="TSP", "PB"="TSP","mPB"="TSP", "thymus"="ETP")
sub1$newtype <- Idents(sub1)
TSP_ETP.markers <- FindAllMarkers(sub1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.405)
TSPID <- subset(TSP_ETP.markers, subset= cluster=="TSP")
TSP_gene <- bitr(TSPID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
ETPID <- subset(TSP_ETP.markers, subset= cluster=="ETP")
ETP_gene <- bitr(ETPID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
cp = list(TSP.gene=TSP_gene$ENTREZID, ETP.gene=ETP_gene$ENTREZID)
#GO analysis
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(GOSemSim)
library(DOSE)
ETP_markers <- read.csv("/home/yushiya/cd34/fig/cpstim_ETP.csv",stringsAsFactors = F)
BMID <- subset(ETP_markers, subset= cluster=="BM")
BM_gene <- bitr(BMID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
PBID <- subset(ETP_markers, subset= cluster=="PB")
PB_gene <- bitr(PBID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
mPBID <- subset(ETP_markers, subset= cluster=="mPB")
mPB_gene <- bitr(mPBID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
thymusID <- subset(ETP_markers, subset= cluster=="thymus")
thymus_gene <- bitr(thymusID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
cp = list(BM.gene=BM_gene$ENTREZID, thymus.gene=thymus_gene$ENTREZID)
cp = list(mPB.gene=mPB_gene$ENTREZID, thymus.gene=thymus_gene$ENTREZID)
cp = list(PB.gene=PB_gene$ENTREZID, thymus.gene=thymus_gene$ENTREZID)
TSP_gene <- rbind(BM_gene,mPB_gene,PB_gene)
TSP_gene <- TSP_gene$ENTREZID %>% unique()
cp = list(TSP.gene=TSP_gene, ETP.gene=thymus_gene$ENTREZID)
cp = list(BM.gene=BM_gene$ENTREZID,mPB.gene=mPB_gene$ENTREZID, PB.gene=PB_gene$ENTREZID)
go.p <- compareCluster(cp,
                       fun = "enrichGO",
                       OrgDb = "org.Hs.eg.db",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.01
)
go.p1 <- simplify(go.p,cutoff=0.4,by="p.adjust",select_fun=min)  #去除冗余
p1=ggplot(go.p1, aes(Cluster, Description), showCategory=5) +
  geom_point(aes(color=p.adjust, size=GeneRatio))+
  theme_classic()+
  theme(axis.text.x =element_text(angle = 45,vjust = 0.5,hjust = 0.5),
        text = element_text(size = 15),
        panel.background = element_blank(),
        panel.grid = element_line(colour = "grey90"))#+theme_bw()
ggsave("/home/yushiya/cd34/fig/fig4-GO-TSPs.pdf", plot=p1, width=7, height=4)
A <- data.frame(go.p1@compareClusterResult[["Description"]])
colnames(A)[1] <- "Description"
A$GeneRatio <- go.p1@compareClusterResult[["GeneRatio"]]
split_num <- data.frame(t(data.frame(strsplit(A$GeneRatio,"/"))))
A$GeneRatio_num <- as.numeric(split_num$X1)/as.numeric(split_num$X2)
A$p.adjust <- as.numeric(go.p1@compareClusterResult[["p.adjust"]])
A$p.adjust <- formatC(A$p.adjust, format = "e", digits = 2) 
A$group <- go.p1@compareClusterResult[["Cluster"]]
A$GeneRatio_num[21:36] <- -A$GeneRatio_num[21:36]   #TSP ETP
A_draw <- A[c(1:4,21:26),]  #TSP ETP
A$GeneRatio_num[16:31] <- -A$GeneRatio_num[16:31]   #PB thymus
A_draw <- A[c(1:5,17:21),]  #PB thymus
A$GeneRatio_num[20:35] <- -A$GeneRatio_num[20:35]   #BM thymus
A_draw <- A[c(1:5,20:25),]  #BM thymus
A$GeneRatio_num[20:35] <- -A$GeneRatio_num[20:35]   #mPB thymus
A_draw <- A[c(1:4,20:25),]  #mPB thymus
p1=ggplot(A_draw,aes(reorder(Description, GeneRatio_num),GeneRatio_num,fill=group))+
  geom_col()+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color="black",size=10),
        axis.line.x = element_line(color='black'),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = 'none')+
  coord_flip()+
  geom_segment(aes(y=0, yend=0,x=0,xend=11))+
  geom_text(data = A_draw[which(A_draw$GeneRatio_num>0),],aes(x=Description, y=-0.01, label=Description),
            hjust=1, size=4)+
  geom_text(data = A_draw[which(A_draw$GeneRatio_num<0),],aes(x=Description, y=0.01, label=Description),
            hjust=0, size=4)+
  geom_text(data = A_draw[which(A_draw$GeneRatio_num>0),],aes(label=p.adjust),
            hjust=-0.1, size=3, color='red')+
  geom_text(data = A_draw[which(A_draw$GeneRatio_num<0),],aes(label=p.adjust),
            hjust=1.1, size=3, color="red")+
  scale_fill_manual(values = c("#98ADC4",
                               "#DDA0DD"))+
  scale_x_discrete(expand = expansion(mult = c(0,0)))+
  ylim(-0.1, 0.1)+
  labs(x='', y='GeneRatio')
#color TSP, ETP, "#DDA0DD","#98ADC4","#7F689D","#54689A"
ggsave("/home/yushiya/cd34/fig/fig4-GO-TSP-ETP.pdf", plot=p1, width=8.5, height=5.5)


#thymus stromal cells
Convert("/home/yushiya/cd34/data/TEC/GSE147520_all_cells.h5ad", dest = "h5seurat", overwrite = F)
tec <- LoadH5Seurat("/home/yushiya/cd34/data/TEC/GSE147520_all_cells.h5seurat",meta.data = T)
VlnPlot(tec, features = c("n_genes", "n_counts", "percent_mito"), ncol = 3)
Idents(tec) <- tec$cell_types
color <- c(brewer.pal(12,"Set3"))
p1=DimPlot(tec, reduction = "umap", label = T, repel=T, cols=color, raster=T)
ggsave("/home/yushiya/cd34/fig/fig5-tec-umap.pdf", plot=p1, width=6, height=4)
sub_tec <- subset(tec, subset= samples=="Postnatal (6 days)"|samples=="Postnatal (10 months)"|samples=="Adult (25 yo)")
sub_tec <- subset(tec, subset= cell_types=="Epithelium-1" | cell_types=="Epithelium-2" | cell_types=="Epithelium-3"|cell_types=="Mesenchyme" |cell_types=="Endo-1 (venous)"|cell_types=="Endo-2 (arterial)"|cell_types=="Endo-3 (venous)"|cell_types=="Endo-4 (lymph.)" )
sub_tec <- RenameIdents(sub_tec,"Epithelium-1"="Thymic_Epi","Epithelium-2"="Thymic_Epi","Epithelium-3"="Thymic_Epi","Mesenchyme"="Thymic_Mesen","Endo-1 (venous)"="Thymic_Endo","Endo-2 (arterial)"="Thymic_Endo","Endo-3 (venous)"="Thymic_Endo","Endo-4 (lymph.)"="Thymic_Endo")
p1=DimPlot(sub_tec, reduction = "umap", label = T, repel=T, cols=c("#66C2A5","#FC8D62","#8DA0CB"), raster=T)
ggsave("/home/yushiya/cd34/fig/fig5-tec-sub_umap.pdf", plot=p1, width=6, height=4)
p1=FeaturePlot(sub_tec, features = c("EPCAM","KRT8","PDGFRA","LUM","PECAM1","ACKR1"), cols = c("#89B7D3", "red"),ncol=2, raster=T)
ggsave("/home/yushiya/cd34/fig/fig5-tec-FeaP.pdf", plot=p1, width=6, height=8)
tec.markers <- FindAllMarkers(sub_tec, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1)
write.csv(tec.markers, file=paste("/home/yushiya/data/cd34/manu_fig/v4/fig5_DEG_tec.csv",sep=""))
tec.markers <- read.csv("/home/yushiya/data/cd34/manu_fig/v4/fig5_DEG_tec.csv")
FeaturePlot(sub_tec, features = c("CCL19","CCL21","CCL2","CCL8"), cols = c("#89B7D3", "red"),ncol=2, raster=T, label=T)
p1=DotPlot(sub_tec, features = c("ccl25","CCL19","CCL21","CCL2","CCL13","CCL11","CCL8","CCL14","CCL23"), cols = c("blue","red"), dot.scale = 8) + 
  RotatedAxis()
ggsave("/home/yushiya/cd34/fig/fig5-tec-sub_dot.pdf", plot=p1, width=7, height=3)
sub_tec$new_celltype <- Idents(sub_tec)
sub_tec <- readRDS("/home/yushiya/cd34/data/TEC_sub.rds")
sub_tec <- subset(sub_tec, subset= new_celltype=="Thymic_Epi" | new_celltype=="Thymic_Endo")
#saveRDS(sub_tec, file = "/home/yushiya/cd34/data/TEC_sub.rds")
cytokines <- read.csv("/home/yushiya/data/cd34/cytokines.txt",header=F)
tec_cyto <- tec.markers[tec.markers$gene %in% cytokines$V1,]
all_thy <- merge(sub1,sub_tec)
Idents(all_thy) <- factor(Idents(all_thy),levels=c("TSP_BM","TSP_mPB","TSP_PB","ETP","Thymic_Epi","Thymic_Mesen","Thymic_Endo"))
p1=DotPlot(all_thy, features = c("CCL19","CCL21","CCL25","CXCL14","CCL2","CCL14","TNFSF10","TNFRSF4"), cols = c("blue","red"), dot.scale = 8) + 
  RotatedAxis()
ggsave("/home/yushiya/cd34/fig/fig5-tec-sub_dot1.pdf", plot=p1, width=7, height=3)


sc.data <- Read10X(data.dir="/home/yushiya/data/blood_atlas/2022_Science_ImmuneMap/A29/BMA/A29_BMA/outs/multi/count/raw_feature_bc_matrix")
sc <- CreateSeuratObject(counts = sc.data, project = "SC")
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
sc <- subset(sc, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
sc <- LoadH5Seurat("/home/yushiya/data/blood_atlas/2022_Science_ImmuneMap/global.h5seurat",meta.data = T) 

