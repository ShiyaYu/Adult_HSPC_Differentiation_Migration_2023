library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)

#Integrate by Seurat CCA
thymus.combined <- readRDS(file = "/home/yushiya/data/cd34/data/thymus_combined-c7-add-1.rds")
SP.combined <- readRDS(file = "/home/yushiya/data/cd34/data/SP.combined.rds")
PB.combined <- readRDS(file = "/home/yushiya/data/cd34/data/PB.combined.rds")
BM.combined <- readRDS(file = "/home/yushiya/data/cd34/data/BM_combined.rds")
mPB.combined <- readRDS(file = "/home/yushiya/data/cd34/data/mPB.combined.rds")
immune.combined <- readRDS("/home/yushiya/data/cd34/data/immune.combined-celltype.rds")
#,k.anchor = 10, k.score = 100, k.filter=50
immune.anchors <- FindIntegrationAnchors(object.list = c(BM.combined,mPB.combined,PB.combined,SP.combined,thymus.combined), assay = c("integrated","RNA","RNA","RNA","integrated"), k.filter = 50, dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)
DefaultAssay(immune.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 50, verbose = FALSE)
ElbowPlot(immune.combined)
#Clustering
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:26)
immune.combined <- FindClusters(immune.combined, reduction = "pca", resolution = 0.95)
immune.combined <- FindClusters(immune.combined, reduction = "pca", resolution = 2.5)
#umap
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:26)
immune.combined <- RunTSNE(immune.combined, reduction = "pca", dims = 1:30)
#cell cycle 
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
immune.combined <- CellCycleScoring(immune.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
immune.combined$Phase <- factor(immune.combined$Phase, levels=c("G1","S","G2M"))
#colours
#HSPCs "#550000", "#AA3939", "#CC6F66", "#E79492","#FBD2CE"
#EryMk "#F2B662","#F38F44","#B97802","#915900"
#B "#A5D6A7", "#47AF50", "#2E7D32",
#T "#8290BC","#816CAA","#7B2C7B","#472349",
#My "#49AFB8","#3790B3","#3465A0","#2C2A50"
#color_all <- c("#550000", "#AA3939", "#CC6F66","#E79492","#FBD2CE","#F2B662","#F38F44","#B97802","#915900","#49AFB8","#3790B3","#3465A0","#2C2A50","#A5D6A7", "#47AF50", "#2E7D32","#8290BC","#816CAA","#7B2C7B","#472349")
color_all <- c("#550000", "#AA3939", "#CC6F66","#E79492","#FBD2CE","#F2B662","#F38F44","#B97802","#915900","#49AFB8","#3465A0","#2C2A50","#A5D6A7", "#47AF50", "#2E7D32","#8290BC","#816CAA","#7B2C7B","#472349")
color_BMPB_lineages <- c("#550000", "#AA3939", "#CC6F66","#E79492","#FBD2CE","#F2B662","#F38F44","#B97802","#915900","#49AFB8","#3465A0","#2C2A50","#A5D6A7", "#47AF50", "#2E7D32")
color_all_mye <- c("#550000", "#AA3939", "#CC6F66","#E79492","#FBD2CE","#49AFB8","#3465A0","#2C2A50")
color_all_B <- c("#550000", "#AA3939", "#CC6F66","#E79492","#FBD2CE","#A5D6A7", "#47AF50", "#2E7D32")
color_T_lin <- c("#550000", "#AA3939", "#CC6F66","#E79492","#FBD2CE","grey","grey","grey","grey","grey","grey","grey","grey", "grey", "grey","#8290BC","#816CAA","#7B2C7B","#472349")
color_stim <- c("#DB5C25","#F3B747","#649541","#AF86BA","#4C82C5")
color_stim <- c("#DB5C25","#F3B747","#649541","#4C82C5")
color_stim <- c("#DB5C25","#649541")
color_cycle <- c("#EA63A2","#FDD685","#52C6EC")
#visualization
DimPlot(immune.combined, reduction = "umap", label = T, repel=T)
DimPlot(immune.combined, reduction = "umap", split.by="stim", ncol=3)
DimPlot(immune.combined, reduction = "umap", split.by="ori_stim", ncol=8)
p1=DimPlot(immune.combined, reduction = "umap", label = T, cols=color_all, repel=T,raster.dpi = c(300,300))
ggsave("/home/yushiya/data/cd34/data/fig/fig1-integrated.pdf", plot=p1, width=6, height=4)
p1=DimPlot(immune.combined, reduction = "umap", split.by="stim", cols=color_stim, ncol=2, raster.dpi = c(200,200))
ggsave("/home/yushiya/data/cd34/data/fig/fig1-colorBYstim.pdf", plot=p1, width=5, height=6)
p1=DimPlot(immune.combined, reduction = "umap", split.by="stim", ncol=3)
ggsave("/home/yushiya/data/cd34/data/fig/fig1-splitBYstim.pdf", plot=p1, width=4, height=10)
p1=DimPlot(immune.combined, reduction = "umap", group.by="Phase", cols=color_cycle)
ggsave("/home/yushiya/data/cd34/data/fig/fig1-groupBYphase.pdf", plot=p1, width=6, height=4)
p1=DimPlot(immune.combined, reduction = "umap", cols=color_T_lin, label=T, raster.dpi = c(300,300))
ggsave("/home/yushiya/data/cd34/data/fig/fig4-T_umap.pdf", plot=p1, width=6, height=4)
p1=DimPlot(immune.combined, reduction = "umap", split.by="stim", ncol=1, cols=color_all)
ggsave("/home/yushiya/data/cd34/data/fig/fig6-allUMAP.pdf", plot=p1, width=6, height=13)
DimPlot(immune.combined, reduction = "umap", split.by="ori_stim", ncol=2)
DimPlot(immune.combined, reduction = "tsne", label = T)
DimPlot(immune.combined, reduction = "tsne", group.by="stim")
DefaultAssay(immune.combined) <- "RNA"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
table(immune.combined@meta.data[["stim"]],immune.combined@meta.data[["seurat_clusters"]])
VlnPlot(immune.combined, features = "CD34", combine = FALSE)
p1=FeaturePlot(immune.combined, features = c("CD34"), max.cutoff = 3, cols = c("#FFE9AE", "#981C12"))
ggsave("/home/yushiya/data/cd34/manu_fig/fig1-CD34.pdf", plot=p1, width=6, height=4)
p1=FeaturePlot(immune.combined, features = c("AVP","MPO","IRF8","CCR7","CD7","CD1A","CD79A","GATA2","PPBP"), max.cutoff = 3, cols = c("#C5E9F5", "#981C12"))
p1=FeaturePlot(immune.combined, features = c("AVP","GATA2","CD36","PPBP","HDC","MPO","LYZ","IRF8","CCR7","VPREB1","CD44","CD7","CD3D","RAG2","CD1A"),
               max.cutoff = 3, cols = c("#C5E9F5", "#981C12", "#981C12"), order=T, ncol=5,
               raster=TRUE, raster.dpi = c(300,300))
ggsave("/home/yushiya/data/cd34/data/fig/fig1-Feamarkers-c1.pdf", plot=p1, width=15, height=7)
immune.combined.markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(immune.combined.markers, "/home/yushiya/data/cd34/manu_fig/v4/celltype_markers.csv")
Idents(immune.combined) <- immune.combined$seurat_clusters
immune.combined.markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(immune.combined.markers, "/home/yushiya/data/cd34/manu_fig/v4/seurat-cluster_markers.csv")
markers.to.plot <- c("EMCN","MEG3","AVP","CRHBP","GATA2","GATA1","HBD","CD36","PPBP","PF4","TPSAB1","HDC","MPO",
                     "ELANE","LYZ","IRF8","CD68","CCR7","NKG7","IL7R","CD19","VPREB1","PAX5",
                     "CD7","CD3D","RAG2","CD1A")
plot3=DotPlot(immune.combined, features = markers.to.plot, cols = c(c("grey","red"),c("grey","red"),c("grey", "red")), dot.scale = 8) + 
  RotatedAxis()
plot3=DotPlot(immune.combined, features = markers.to.plot, cols = c(c("#150788","#F4E83A")), dot.scale = 6) + 
  RotatedAxis()
plot3=DotPlot(immune.combined, features = markers.to.plot, dot.scale = 6) +  
  #scale_color_gradientn(colors = c("#F0F921FF","#FDC926FF","#FA9E3BFF","#ED7953FF","#ED7953FF","#BD3786FF","#9C179EFF","#7301A8FF","#47039FFF","#0D0887FF")) +
  scale_color_gradientn(colors = viridis::plasma(10)) +
  theme(panel.border = element_rect(color = "black", fill=NA),
        panel.background = element_blank()) + 
  RotatedAxis()
ggsave("Dotplot.pdf", plot=plot3, width=16, height=7)
ggsave("/home/yushiya/data/cd34/data/fig/fig1-markers_c2.pdf", plot=plot3, width=10, height=5)
saveRDS(immune.combined, file = "/home/yushiya/data/cd34/data/immune.combined-celltype.rds")
#save h5ad
immune.combined$tissue_celltype <- paste(immune.combined$stim, immune.combined$celltype, sep="_")
immune.combined$celltype <- as.character(immune.combined$celltype)
SaveH5Seurat(immune.combined,filename="/home/yushiya/data/cd34/data/fig/immune_combined.h5seurat", overwrite = TRUE)
Convert("/home/yushiya/data/cd34/data/fig/immune_combined.h5seurat", dest = "h5ad", overwrite = TRUE)
