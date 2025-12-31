library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(devtools)
library(slingshot)
library(RColorBrewer)

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

#monocle3
library(Matrix)
library(ggplot2)
library(monocle3)
##创建CDS对象并预处理数据
data <- GetAssayData(immune.combined, assay = 'RNA', slot = 'counts')
cell_metadata <- immune.combined@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
save(cds,file="/home/yushiya/data/cd34/data/fig/immune.combined.cds.RData")
load(file="/home/yushiya/data/cd34/data/fig/immune.combined.cds.RData")
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
cds <- learn_graph(cds,learn_graph_control=list(geodesic_distance_ratio=0.5))
p<-plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
              label_branch_points = FALSE)
pdf("/home/yushiya/data/cd34/data/fig/fig1-monocle.pdf", width = 6, height = 4, dpi=200)  
plot()
dev.off()
# cds <- order_cells(cds) 存在bug，使用辅助线选择root细胞
p + geom_vline(xintercept = seq(5,7,0.25)) + geom_hline(yintercept = seq(4.5,5.5,0.25))
embed <- data.frame(Embeddings(immune.combined, reduction = "umap"))
embed <- subset(embed, UMAP_1 > 5.5 & UMAP_1 < 6.5 & UMAP_2 > 4.2 & UMAP_2 < 5.5)
root.cell <- rownames(embed)
cds <- order_cells(cds, root_cells = root.cell)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = FALSE)
plot_cells(cds, color_cells_by = "celltype", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = FALSE)
immune.combined$pseudotime <- pseudotime(cds)
pseudotime <- read.csv("/home/yushiya/data/cd34/data/fig/pseudotime.csv")
immune.combined$pseudotime <- pseudotime$x
ggsave("/home/yushiya/data/cd34/data/cd34/fig/fig1-pseu.pdf", plot=p1, width=6, height=4)
#plot HSC genes along pseudotime
genes <- c("CD34", "AVP", "MEG3", "ID1")
lineage_cds <- cds[rowData(cds)$gene_short_name %in% genes,]
plot_genes_in_pseudotime(lineage_cds,color_cells_by="pseudotime",min_expr=2)
plot_cells(cds, genes=c("SPP1","AVP","MEG3","ID1"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE)
#find genes through pseudotime
sub_cds <- cds[,colData(cds)$stim %in% c("BM","mPB","PB")]
sub_cds <- sub_cds[,colData(sub_cds)$celltype %in% c("HSC","MPP#1","MPP#2","LMPP#1","LMPP#2","BMEP","EryP","MkP","Ma/Eo/BaP","GMP","CDP","pre-pDC","CLP","pro-B","pre-B")]
plot_cells(sub_cds, color_cells_by="celltype")
#pr_graph_test_res <- graph_test(sub_cds, neighbor_graph="knn", cores=40)
#pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))
#write.csv(pr_graph_test_res, "/home/yushiya/cd34/fig/peu_gene.csv")
cds_pr_test_res <- graph_test(sub_cds, neighbor_graph="principal_graph", cores=80)  #DEGs along trajectory
write.csv(cds_pr_test_res, "/home/yushiya/data/cd34/data/fig/peu_tra_gene_all.csv")
cds_pr_test_res <- read.csv("/home/yushiya/data/cd34/data/fig/peu_tra_gene_all.csv", row.names = 1)
pr_deg_ids <- row.names(subset(cds_pr_test_res, q_value < 0.01 & morans_I > 0.2))
cell_group_df <- tibble::tibble(cell=row.names(colData(sub_cds)), 
                                cell_group=colData(sub_cds)$celltype)
agg_mat <- aggregate_gene_expression(sub_cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
write.table(t(agg_mat), "/home/yushiya/data/cd34/data/fig/fig2-peu_tra_module-aggmat-BM-1.csv", sep = "\t")
agg_mat <- read.table("/home/yushiya/data/cd34/data/fig/fig2-peu_tra_module-aggmat-BM-1.csv", sep="\t")
#colnames(agg_mat) <- factor(colnames(agg_mat), levels=c("HSC","MPP #1","MPP #2","LMPP #1","LMPP #2","MEP","EryP","MkP","CLP","pro-B","pre-B","GMP","CDP","pre-pDC"))
p<-pheatmap::pheatmap(agg_mat,cluster_rows = T, cluster_cols = T,
                      scale="column", clustering_method="ward.D2")
ggsave("/home/yushiya/data/cd34/data/fig/fig2-tra-module-heatmap-BM.pdf", plot=p, width=6, height=4)
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
write.csv(Module_GO1, file = '/home/yushiya/data/cd34/data/fig/fig2-Module_GO_BM-1.csv')

#correlation of module genes in BM, mPB, PB
library(pheatmap)
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
gene_module_BM <- read.csv("/home/yushiya/data/cd34/data/fig/fig2-peu_tra_module-BM-1.csv", row.names = 1)
gene_module_mPB <- read.csv("/home/yushiya/data/cd34/data/fig/fig2-peu_tra_module-mPB.csv", row.names = 1)
gene_module_PB <- read.csv("/home/yushiya/data/cd34/data/fig/fig2-peu_tra_module-PB.csv", row.names = 1)
go_module_BM <- read.csv("/home/yushiya/data/cd34/data/fig/fig2-Module_GO_BM-1.csv", row.names = 1)
go_module_mPB <- read.csv("/home/yushiya/data/cd34/data/fig/fig2-Module_GO_mPB.csv", row.names = 1)
go_module_PB <- read.csv("/home/yushiya/data/cd34/data/fig/fig2-Module_GO_PB.csv", row.names = 1)
gene_module_BM$module <- stringr::str_c("Module ", gene_module_BM$module)
gene_module_mPB$module <- stringr::str_c("Module ", gene_module_mPB$module)
gene_module_PB$module <- stringr::str_c("Module ", gene_module_PB$module)
# 提取每个module的id
get_ids_by_module <- function(df) {
  df %>%
    group_by(module) %>%
    summarise(ids = list(id)) %>%
    mutate(module = as.character(module))
}
grouped_ids1 <- get_ids_by_module(gene_module_BM)
grouped_ids2 <- get_ids_by_module(gene_module_mPB)
# 计算Jaccard相似性
calculate_jaccard <- function(ids1, ids2) {
  intersection <- length(intersect(unlist(ids1), unlist(ids2)))
  union <- length(union(unlist(ids1), unlist(ids2)))
  intersection / union
}
# 生成所有module对并计算相似性
modules1 <- grouped_ids1$module
modules2 <- grouped_ids2$module
# 创建模块对的组合
module_combinations <- expand.grid(BM = modules1, mPB = modules2)
# 计算每个模块对的相似性
module_combinations$jaccard_similarity <- apply(module_combinations, 1, function(row) {
  ids1 <- grouped_ids1 %>% filter(module == row[1]) %>% pull(ids)
  ids2 <- grouped_ids2 %>% filter(module == row[2]) %>% pull(ids)
  if (length(ids1) == 0 || length(ids2) == 0) {
    return(NA)
  }
  calculate_jaccard(ids1, ids2)
})
# 查看结果
print(module_combinations)
# 将相似性结果转换为矩阵形式
similarity_matrix <- module_combinations %>%
  na.omit() %>%
  pivot_wider(names_from = mPB, values_from = jaccard_similarity) %>%
  select(-BM) %>%
  as.matrix()
rownames(similarity_matrix) <- unique(module_combinations$BM)
# 绘制热图
p<-pheatmap(
  similarity_matrix,
  color = viridis::viridis(100),
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  main = "Module Similarity Heatmap"
)
ggsave("/home/yushiya/data/cd34/data/fig/fig2-tra-module-similarity-BMmPB.pdf", plot=p, width=5, height=4)
#3d figures
library(plotly)
# 提取每个module的id
get_ids_by_module <- function(df) {
  df %>%
    group_by(module) %>%
    summarise(ids = list(id)) %>%
    mutate(module = as.character(module))
}
grouped_ids1 <- get_ids_by_module(gene_module_BM)
grouped_ids2 <- get_ids_by_module(gene_module_mPB)
grouped_ids3 <- get_ids_by_module(gene_module_PB)
# 计算Jaccard相似性
calculate_jaccard <- function(ids1, ids2) {
  intersection <- length(intersect(unlist(ids1), unlist(ids2)))
  union <- length(union(unlist(ids1), unlist(ids2)))
  intersection / union
}
# 生成所有module对并计算相似性
modules1 <- grouped_ids1$module
modules2 <- grouped_ids2$module
modules3 <- grouped_ids3$module
# 创建模块对的组合
module_combinations <- expand.grid(BM = modules1, mPB = modules2, PB = modules3)
# 计算每个模块对的相似性
module_combinations$jaccard_similarity_BM_mPB <- apply(module_combinations, 1, function(row) {
  ids1 <- grouped_ids1 %>% filter(module == row[1]) %>% pull(ids)
  ids2 <- grouped_ids2 %>% filter(module == row[2]) %>% pull(ids)
  if (length(ids1) == 0 || length(ids2) == 0) {
    return(NA)
  }
  calculate_jaccard(ids1, ids2)
})
module_combinations$jaccard_similarity_BM_PB <- apply(module_combinations, 1, function(row) {
  ids1 <- grouped_ids1 %>% filter(module == row[1]) %>% pull(ids)
  ids3 <- grouped_ids3 %>% filter(module == row[3]) %>% pull(ids)
  if (length(ids1) == 0 || length(ids3) == 0) {
    return(NA)
  }
  calculate_jaccard(ids1, ids3)
})
module_combinations$jaccard_similarity_mPB_PB <- apply(module_combinations, 1, function(row) {
  ids2 <- grouped_ids2 %>% filter(module == row[2]) %>% pull(ids)
  ids3 <- grouped_ids3 %>% filter(module == row[3]) %>% pull(ids)
  if (length(ids2) == 0 || length(ids3) == 0) {
    return(NA)
  }
  calculate_jaccard(ids2, ids3)
})
# 计算三个数组的交集相似性
module_combinations$jaccard_similarity_all <- apply(module_combinations, 1, function(row) {
  ids1 <- grouped_ids1 %>% filter(module == row[1]) %>% pull(ids)
  ids2 <- grouped_ids2 %>% filter(module == row[2]) %>% pull(ids)
  ids3 <- grouped_ids3 %>% filter(module == row[3]) %>% pull(ids)
  if (length(ids1) == 0 || length(ids2) == 0 || length(ids3) == 0) {
    return(NA)
  }
  # 计算三者的交集和并集
  intersection <- length(intersect(intersect(unlist(ids1), unlist(ids2)), unlist(ids3)))
  union <- length(union(union(unlist(ids1), unlist(ids2)), unlist(ids3)))
  intersection / union
})
# 查看结果
print(module_combinations)
#module_combinations <- module_combinations %>%
#  filter(!is.na(jaccard_similarity_all) & jaccard_similarity_all != 0)
#module_combinations$BM <- as.numeric(as.character(module_combinations$BM))
#module_combinations$mPB <- as.numeric(as.character(module_combinations$mPB))
#module_combinations$PB <- as.numeric(as.character(module_combinations$PB))
module_combinations$jaccard_similarity_all <- as.numeric(module_combinations$jaccard_similarity_all)
plot_ly(module_combinations, x = ~jaccard_similarity_BM_PB, y = ~jaccard_similarity_BM_mPB, z = ~jaccard_similarity_mPB_PB, 
        color = ~jaccard_similarity_all, 
        colors = c("#0D0887FF", "#47039FFF", "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF", "#ED7953FF", "#FA9E3BFF")) %>%
  add_markers(size = 5) %>%  # 设置点的大小
  layout(
    scene = list(
      xaxis = list(title = 'BM PB Similarity', range = c(0, 0.7), 
                   titlefont = list(size = 17),  # 设置x轴标题的字体大小
                   tickfont = list(size = 12)),    # 设置x轴刻度标签的字体大小),
      yaxis = list(title = 'BM mPB Similarity', range = c(0, 0.6), 
                   titlefont = list(size = 17),  
                   tickfont = list(size = 12)),
      zaxis = list(title = 'mPB PB Similarity', range = c(0, 0.6), 
                   titlefont = list(size = 17),  
                   tickfont = list(size = 12)),
      coloraxis = list(colorbar = list(title = 'Jaccard Similarity'))
    )
  )

ggsave("/home/yushiya/data/cd34/data/fig/fig2-tra-module-similarity-3d.pdf", plot=p1, width=5, height=4)

#calculate migration score by GO modules
#BM CSF3R/SELL/IL1B/CD74/AIF1/ANXA1/MPP1
#mPB CSF3R/SELL/CD74/AIF1/ANXA1/MDK 	
#PB CSF3R/SELL/IL1B/ANXA1/MDK/CTSG/FOXJ1/AZU1/PRTN3/ELANE/MPP1
migr_genes <- c("CSF3R","SELL","IL1B","CD74","AIF1","ANXA1","MPP1",
                "MDK","CTSG","FOXJ1","AZU1","PRTN3","ELANE")
migr_genes <- c("CSF3R","SELL","IL1B","CD74","AIF1","ANXA1","MPP1")
sub_cds_BM <- cds[,colData(cds)$stim %in% c("BM")]
sub_cds_BM <- sub_cds_BM[,colData(sub_cds_BM)$celltype %in% c("HSC","MPP#1","MPP#2","LMPP#1","LMPP#2","BMEP","EryP","MkP","Ma/Eo/BaP","GMP","CDP","pre-pDC","CLP","pro-B","pre-B")]
genes_module=gene_module_BM %>% filter(module %in% c("Module 3"))
genes_module <- genes_module[,c("Module 1")]
p1=plot_cells(sub_cds_BM,
              genes=migr_genes,
              label_cell_groups=FALSE,
              show_trajectory_graph=FALSE,
              rasterize = TRUE)
p1=plot_cells(sub_cds_BM,
              genes=gene_module_BM %>% filter(module %in% c(3,5)),
              label_cell_groups=FALSE,
              show_trajectory_graph=FALSE,
              rasterize = TRUE)
ggsave("/home/yushiya/data/cd34/data/fig/fig2-pseu-umap-BM-module3.pdf", plot=p1, width=6, height=4)
sub_cds_BM_lin <- sub_cds_BM[rowData(sub_cds_BM)$gene_short_name %in% migr_genes,]
p1=plot_genes_in_pseudotime(sub_cds_BM_lin,
                            #color_cells_by="celltype",
                            min_expr=0.5)
ggsave("/home/yushiya/data/cd34/data/fig/fig2-pseu-BM-migrgenes.pdf", plot=p1, width=5, height=12)
#migratory gene in module
migr_genes <- c("CSF3R","SELL","IL1B","CD74","AIF1","ANXA1","MPP1","MDK")
sub1 <- AddModuleScore(object = sub1, features = list(Migr_score=migr_genes), name="Migr_score")
VlnPlot(sub1, features = c("Migr_score1"), cols = color_stim, split.by="tissue", pt.size = 0)
#chemotaxis gene in module
chemo_gene <- c("CSF3R","CSF1","SELL","IL1B","PF4V1","CXCL5","CD74","AIF1","ANXA1","MDK","CCL5","CXADR","MPP1")
chemo_gene <- c("PF4V1","PF4","PPBP","CMTM5","THBS1","CCL5","ITGB3")
sub1 <- AddModuleScore(object = sub1, features = list(chemo_score=chemo_gene), name="chemo_score")
VlnPlot(sub1, features = c("chemo_score1"), cols = color_stim, split.by="tissue", pt.size = 0)
#adhesion gene in module
adhesion_gene <- c("IL1B","IGFBP2","ANXA1","MDK","ANGPT1","ZBTB16","TESPA1","CTSG","SKAP1")
sub1 <- AddModuleScore(object = sub1, features = list(adhesion_score=adhesion_gene), name="adhesion_score")
VlnPlot(sub1, features = c("adhesion_score1"), cols = color_stim, split.by="tissue", pt.size = 0)
#differentiation genes in DEG
#draw boxplot
miGene <- as.data.frame(colnames(sub1))
miGene$tissue <- sub1$tissue
miGene$celltype <- sub1$celltype
miGene$tissue <- factor(miGene$tissue, levels=c("BM","mPB","PB"))
names(miGene)[1]='ID'
miGene$Migratory_Score <- sub1$Migr_score1
miGene$Adhesion_Score <- sub1$adhesion_score1
miGene$Differentiation_Score <- sub1$Diff_Score1
miGene$Chemotaxis_Score <- sub1$chemo_score1 
miGene <- subset(miGene, subset= celltype %in% c("HSC","MPP#1","MPP#2","LMPP#1","LMPP#2","GMP","CDP","pre-pDC"))
miGene <- subset(miGene, subset= celltype %in% c("HSC","MPP#1","MPP#2","LMPP#1","LMPP#2","BMEP","EryP","MkP","Ma/Eo/BaP"))
miGene <- subset(miGene, subset= celltype %in% c("HSC","MPP#1","MPP#2","LMPP#1","LMPP#2","CLP","pro-B","pre-B"))
#boxplot with lineages in the same facet
p1=ggboxplot(miGene, x="celltype", y="Differentiation_Score", color="tissue", bxp.errorbar = T,
             add="boxplot", short.panel.labs = TRUE, ncol=4, add.params=list(size=0.2))+
  stat_compare_means(label = "p.signif",method = "t.test",ref.group = ".all.")+
  scale_color_manual(values = color_stim) +
  theme(axis.text.x =element_text(angle = 45,hjust = 1),
        plot.title = element_text(hjust = 0))
ggsave(file="/home/yushiya/data/cd34/data/fig/fig3-migr-BMmPBPB.pdf", plot=p1, width=7, height=3)
ggsave(file="/home/yushiya/data/cd34/data/fig/fig3-adhe-BMmPBPB.pdf", plot=p1, width=7, height=3)
ggsave(file="/home/yushiya/data/cd34/data/fig/fig3-diff-BMmPBPB.pdf", plot=p1, width=7, height=3)
ggsave(file="/home/yushiya/data/cd34/data/fig/fig3-diff-BMmPBPB-B.pdf", plot=p1, width=5, height=3)
#boxplot with lineages in the same facet not split by stim
p1=ggboxplot(miGene, x="celltype", y="Differentiation_Score", color="celltype", bxp.errorbar = T,
             add="boxplot", short.panel.labs = TRUE, ncol=4, add.params=list(size=0.2))+
  stat_compare_means(label = "p.signif",method = "t.test",ref.group = ".all.")+
  scale_color_manual(values = color_all_B) +
  theme(axis.text.x =element_text(angle = 45,hjust = 1),
        plot.title = element_text(hjust = 0))
ggsave(file="/home/yushiya/data/cd34/data/fig/fig3-migr-BMmPBPB-lineages-all.pdf", plot=p1, width=9, height=4)
ggsave(file="/home/yushiya/data/cd34/data/fig/fig3-adhe-BMmPBPB-lineages-all.pdf", plot=p1, width=9, height=4)
ggsave(file="/home/yushiya/data/cd34/data/fig/fig3-diff-BMmPBPB-lineages-all.pdf", plot=p1, width=9, height=4)
ggsave(file="/home/yushiya/data/cd34/data/fig/fig3-diff-BMmPBPB-lineages-all-B.pdf", plot=p1, width=5, height=4)
#plot through pseudotime
sub1_ery <- subset(sub1, subset = celltype %in% c("HSC","MPP#1","MPP#2","LMPP#1","LMPP#2","BMEP","EryP","MkP","Ma/Eo/BaP"))
sub1_mye <- subset(sub1, subset = celltype %in% c("HSC","MPP#1","MPP#2","LMPP#1","LMPP#2","GMP","CDP","pre-pDC"))
sub1_B <- subset(sub1, subset = celltype %in% c("HSC","MPP#1","MPP#2","LMPP#1","LMPP#2","CLP","pro-B","pre-B"))
mtx <- as.data.frame(sub1_B$celltype)
colnames(mtx)[1] <- "celltype"
mtx$tissue <- sub1_B$tissue
mtx$pseudotime <- sub1_B$pseudotime
mtx$CSF3R <- sub1_B@assays[["RNA"]]@data["CSF3R",]
mtx$CD74 <- sub1_B@assays[["RNA"]]@data["CD74",]
mtx$SELL <- sub1_B@assays[["RNA"]]@data["SELL",]
mtx$ANXA1 <- sub1_B@assays[["RNA"]]@data["ANXA1",]
mtx$MPP1 <- sub1_B@assays[["RNA"]]@data["MPP1",]
mtx$IL1B <- sub1_B@assays[["RNA"]]@data["IL1B",]
mtx$AIF1 <- sub1_B@assays[["RNA"]]@data["AIF1",]
mtx$CXCR4 <- sub1@assays[["RNA"]]@data["CXCR4",]
mtx$MDK <- sub1@assays[["RNA"]]@data["MDK",]
mtx$migr_exp <- sub1$m1_exp1
col <- c(color_BMPB_lineages,color_stim)
col <- c("#550000", "#AA3939", "#CC6F66", "#E79492", "#FBD2CE","#F2B662", "#F38F44", "#B97802","#DB5C25", "#F3B747", "#649541")
p1=ggplot(mtx, aes(x=pseudotime,y=IL1B))+
  geom_point(size=0.3,aes(colour=celltype),alpha = 0)+
  stat_smooth(mapping = aes(colour=tissue))+
  #facet_wrap(~celltype)+#theme_bw()
  #scale_y_continuous(limits = c(-0.5,0.5))+
  scale_color_manual(values = col) +
  #scale_color_manual(values = col) +
  theme_classic()+
  theme(axis.text.x=element_text(hjust=0.5),
        text = element_text(size = 18),
        #panel.border = element_rect(color = "black", fill=NA),
        panel.background = element_blank())
#panel.grid = element_line(colour = "grey60"),
#axis.line = element_line(colour = "grey60",size = rel(1),arrow = arrow(angle = 30,length = unit(0.1,"inches"))))
ggsave(file=paste("/home/yushiya/data/cd34/data/fig/fig2-pseu-B-IL1B.pdf.pdf",sep=""), plot=p1, width=8, height=4)



#Venn plot among gene modules
library(VennDiagram)
module_genes_BM <- gene_module_BM %>% filter(module == "Module 3") %>% pull(id)
module_genes_mPB <- gene_module_mPB %>% filter(module == "Module 6") %>% pull(id)
module_genes_PB <- gene_module_PB %>% filter(module == "Module 2") %>% pull(id)
# 计算每个集合的比例
#total <- length(union(union(module_genes_BM, module_genes_mPB),module_genes_PB))
#prob_BM <- length(module_genes_BM) / total
#prob_mPB <- length(module_genes_mPB) / total
#prob_PB <- length(module_genes_PB) / total
p1 <- draw.triple.venn(
  area1 = length(module_genes_BM),
  area2 = length(module_genes_mPB),
  area3 = length(module_genes_PB),
  n12 = length(intersect(module_genes_BM, module_genes_mPB)),
  n13 = length(intersect(module_genes_BM, module_genes_PB)),
  n23 = length(intersect(module_genes_mPB, module_genes_PB)),
  n123 = length(intersect(intersect(module_genes_BM, module_genes_mPB), module_genes_PB)),
  category = c("BM_Module3", "mPB_Module6", "PB_Module2"),
  col = c("#DB5C25","#F3B747","#649541"),
  fill = c("#DB5C25","#F3B747","#649541"),
  alpha = c(0.85, 0.85, 0.85),
  #label.col = c("black", "black", "black"),
  cex = 1.5,
  cat.cex = 1.2,
  margin = 0.05,
  lty = "blank"  #不显示圆圈颜色
)
ggsave("/home/yushiya/data/cd34/data/fig/fig2-tra-module-overlap-7.pdf", plot=p1, width=4.5, height=4)

#draw module exp through pseudotime
sub1_BM <- subset(sub1, subset= tissue=="BM")
sub1_BM <- AddModuleScore(object = sub1_BM, features = list(BM_Module=module_genes_BM), name="BM_Module")
sub1_mPB <- subset(sub1, subset= tissue=="mPB")
sub1_mPB <- AddModuleScore(object = sub1_mPB, features = list(mPB_Module=module_genes_mPB), name="mPB_Module")
sub1_PB <- subset(sub1, subset= tissue=="PB")
sub1_PB <- AddModuleScore(object = sub1_PB, features = list(PB_Module=module_genes_PB), name="PB_Module")
mtx <- as.data.frame(sub1_PB$celltype)
colnames(mtx)[1] <- "celltype"
mtx$pseudotime <- sub1_PB$pseudotime
mtx$tissue <- sub1_PB$tissue
mtx$PB_Module_2 <- sub1_PB$PB_Module1
#col <- c(color_BMPB_lineages,color_stim)
#col <- c("#550000", "#AA3939", "#CC6F66", "#E79492", "#FBD2CE","#F2B662", "#F38F44", "#B97802","#915900","#DB5C25", "#F3B747", "#649541")
p1=ggplot(mtx, aes(x=pseudotime,y=PB_Module_2))+
  geom_point(size=0.3,aes(colour=celltype),alpha = 0)+
  stat_smooth(mapping = aes(colour=tissue))+
  #facet_wrap(~celltype)+#theme_bw()
  #scale_y_continuous(limits = c(-0.5,0.5))+
  #scale_color_manual(values = col) +
  #scale_color_manual(values = col) +
  theme_classic()+
  theme(axis.text.x=element_text(hjust=0.5),
        text = element_text(size = 18),
        #panel.border = element_rect(color = "black", fill=NA),
        panel.background = element_blank())


#compare module differences among tissues
agg_mat_BM <- read.table("/home/yushiya/data/cd34/data/fig/fig2-peu_tra_module-aggmat-BM-1.csv", sep="\t")
agg_mat_mPB <- read.table("/home/yushiya/data/cd34/data/fig/fig2-peu_tra_module-aggmat-mPB.csv", sep="\t")
agg_mat_PB <- read.table("/home/yushiya/data/cd34/data/fig/fig2-peu_tra_module-aggmat-PB.csv", sep="\t")
cp_aggmat <- data.frame(agg_mat_BM$Module.3)
cp_aggmat$v2 <- agg_mat_mPB$Module.6
cp_aggmat$v3 <- agg_mat_PB$Module.2
colnames(cp_aggmat) <- c("BM_Module3", "mPB_Module6", "PB_Module2")
#colnames(cp_aggmat) <- c("BM_Module5", "mPB_Module2", "mPB_Module7", "PB_Module5")
rownames(cp_aggmat) <- rownames(agg_mat_BM)
cp_aggmat <- cp_aggmat[c(6,11,12,7,8,5,2,14),]  #mye
cp_aggmat <- cp_aggmat[c(6,11,12,7,8,1,4,10,9),]  #ery
cp_aggmat <- cp_aggmat[c(6,11,12,7,8,3,15,13),]  #B
cp_aggmat <- cp_aggmat[c(6,11,12,7,8,1,4,10,9,5,2,14,3,15,13),]
p<-pheatmap::pheatmap(t(cp_aggmat),cluster_rows = F, cluster_cols = T,
                      #scale="column", 
                      clustering_method="ward.D2")
ggsave("/home/yushiya/data/cd34/data/fig/fig2-tra-module-cpmodule-B.pdf", plot=p, width=6, height=2.3)

#draw module expression among tissues
library(dplyr)
library(tidyr)
library(reshape2)
library("ggalluvial")
agg_mat <- read.table("/home/yushiya/data/cd34/data/fig/fig2-peu_tra_module-aggmat-1.csv",sep=",")
rownames(agg_mat) <- agg_mat$V1
agg_mat <- data.frame(t(agg_mat[,2:4]))
#mat <- data.frame(agg_mat$Module.1)
#colnames(mat)[1] <- "Exp."
#mat$tissue_celltype <- agg_mat$V1
#mat <- separate(agg_mat, col = V1, into = c("tissue", "celltype"), sep = "_")
mat=melt(agg_mat, id='V1')
colnames(mat)[1] <- "tissue"
mat$value <- as.numeric(mat$value)
mat <- mat %>% mutate(value = round(value, 2))
mat1 <- subset(mat, subset= variable=="Module.6")
#mat1 <- mat1[c(1,3),]
p1=ggplot(data = mat1, aes(x = tissue, y = value, fill = tissue)) + 
  geom_bar(stat = 'identity') +
  labs(title = "Module 6") +
  scale_fill_manual(values=color_stim) +
  geom_text(aes(label = value), vjust=0) +
  theme_classic() +
  theme(axis.title = element_text(color = "black", size = 12),
        axis.text = element_text(size = 10, color = "black"),
        plot.title = element_text(size = 14, hjust = .5, color = "black"))
ggsave("/home/yushiya/data/cd34/data/fig/fig2-tra-module6-mPB-exp.pdf", plot=p1, width=6, height=2)


#draw module go
color_all <- c("#550000", "#AA3939", "#CC6F66","#E79492","#FBD2CE","#F2B662","#F38F44","#B97802","#915900","#49AFB8","#3465A0","#2C2A50","#A5D6A7", "#47AF50", "#2E7D32","#8290BC","#816CAA","#7B2C7B","#472349")
Module_GO_BM <- read.csv('/home/yushiya/data/cd34/data/fig/fig2-Module_GO_BM-1.csv')
Module_GO_mPB <- read.csv('/home/yushiya/data/cd34/data/fig/fig2-Module_GO_mPB.csv')
Module_GO_PB <- read.csv('/home/yushiya/data/cd34/data/fig/fig2-Module_GO_PB.csv')
mytheme <- theme(axis.text.x = element_text(hjust = 0.5,size = 12), 
                 ## 删去y轴label
                 axis.text.y = element_blank(),
                 ## 删去y轴刻度线
                 axis.ticks.y = element_blank(), 
                 ## 删去x\y轴标题
                 #axis.title.x = element_blank(), 
                 #axis.title.y = element_blank(), 
                 axis.title = element_text(size = 20),
                 axis.title.x = element_text(size = 16),  # 修改横轴标题字体大小
                 axis.title.y = element_text(size = 16),  # 修改纵轴标题字体大小
                 plot.title = element_text(hjust = 0.5,size =  18),
                 legend.position = "none")
sub_GO <- subset(Module_GO_PB, subset= cluster==2)
sub_GO$log_pvalue <- -log(as.numeric(sub_GO$pvalue))
sub_GO_top <- sub_GO[c(2,3,5,6,7),]   #BM module 1
sub_GO_top <- sub_GO[c(1,2,3,4,7),]   #BM module 2
sub_GO_top <- sub_GO[c(4,7,8,9,123),]   #BM module 3
sub_GO_top <- sub_GO[c(1,2,4,5,60),]   #BM module 4
sub_GO_top <- sub_GO[c(1,2,3,5,9),]   #BM module 5
sub_GO_top <- sub_GO[c(1,2,3,4,5),]   #BM module 6
sub_GO_top <- sub_GO[c(1,4,8,10,12),]   #BM module 7
sub_GO_top <- sub_GO[c(1,2,3,6,9),]   #BM module 8
sub_GO_top <- sub_GO[c(1,2,3,5,6),]   #mPB module 1
sub_GO_top <- sub_GO[c(1,3,6,7,60),]   #mPB module 2
sub_GO_top <- sub_GO[c(2,5,6,8,68),]   #mPB module 3
sub_GO_top <- sub_GO[c(1,2,4,5,7),]   #mPB module 4
sub_GO_top <- sub_GO[c(1,2,4,6,7),]   #mPB module 5
sub_GO_top <- sub_GO[c(1,4,11,20,61),]   #mPB module 6
sub_GO_top <- sub_GO[c(1,4,9,17,24),]   #mPB module 7
sub_GO_top <- sub_GO[c(1,3,6,7,87),]   #PB module 1
sub_GO_top <- sub_GO[c(1,3,6,9,38),]   #PB module 2
sub_GO_top <- sub_GO[c(1,2,5,6,7),]   #PB module 3
sub_GO_top <- sub_GO[c(1,2,6,9,11),]   #PB module 4
sub_GO_top <- sub_GO[c(1,2,5,10,13),]   #PB module 5
sub_GO_top <- sub_GO[c(1,3,4,7,8),]   #PB module 6
sub_GO_top$Description <- factor(sub_GO_top$Description,levels = rev(sub_GO_top$Description))
p1=ggplot(data = sub_GO_top, aes(x = Description, y = log_pvalue)) +
  geom_bar(stat = "identity", width = 0.8, fill = "#649541",alpha = 0.8) + #绘制条形图
  geom_text(aes(y = 0, #控制文本标签起始位置
                label = Description),
            size = 6,hjust = 0) + #hjust = 0左对齐
  coord_flip() + theme_bw() + mytheme +
  labs(x = "GO Terms", y = "-log p value", title = "PB Module 2 GO Term Enrichment") 
ggsave("/home/yushiya/data/cd34/data/fig/fig2-tra-module-GO-PB-2.pdf", plot=p1, width=5, height=4.5)

#box plot of p value in leukocyte migration among tissues
#BM 3.227947e-03 CSF3R/SELL/IL1B/CD74/AIF1/ANXA1/MPP1
#mPB 3.882165e-05 CSF3R/SELL/CD74/AIF1/ANXA1/MDK
#PB 3.587438e-05 CSF3R/SELL/IL1B/ANXA1/MDK/CTSG/FOXJ1/AZU1/PRTN3/ELANE/MPP1
#cell-cell adhesion 1.434874e-04,6.372566e-08,2.485559e-04
#cytokine production 3.085333e-05,7.066359e-08,6.345075e-05
#chemotaxis 1.139707e-03,7.576118e-04,3.420411e-05
dat <- data.frame(c(3.227947e-03,3.882165e-05,3.587438e-05))
dat <- data.frame(c(1.434874e-04,6.372566e-08,2.485559e-04))
dat <- data.frame(c(3.085333e-05,7.066359e-08,6.345075e-05))
dat <- data.frame(c(1.139707e-03,7.576118e-04,3.420411e-05))
colnames(dat)[1] <- "pvalue"
dat$log_pvalue <- -log(dat$pvalue)
dat$bar <- c("BM_Module3","mPB_Module6","PB_Module2")
p1=ggplot(dat,aes(bar,log_pvalue))+
  geom_col(aes(fill=bar)) + 
  scale_fill_manual(values = color_stim) + 
  theme_classic() +
  theme(axis.text.x=element_text(angle = 45,hjust = 1),
        text = element_text(size = 18, color = "black"),
        #panel.border = element_rect(color = "black", fill=NA),
        panel.background = element_blank())
ggsave("/home/yushiya/data/cd34/data/fig/fig2-tra-module-chemo.pdf", plot=p1, width=5, height=3.5)
#add module genes
migr_genes <- c("CSF3R","SELL","IL1B","CD74","AIF1","ANXA1","MPP1")
sub1 <- AddModuleScore(object = sub1, features = list(module_exp=migr_genes), name="m1_exp")


#chord plot
library(dplyr)
library(stringr)
library(GOplot)
cds_pr_test_res <- read.csv("/home/yushiya/data/cd34/data/fig/peu_tra_gene_BM.csv", row.names = 1)
sub_GO <- subset(Module_GO_BM, subset= cluster==4)
sub_GO$log_pvalue <- -log(as.numeric(sub_GO$pvalue))
sub_GO_top <- sub_GO[c(4,7,8,9,123),]   #BM module 3
sub_GO_top <- sub_GO[c(1,2,3,5,9),]   #BM module 5
sub_GO_top <- sub_GO[c(1,2,4,5,60),]   #BM module 4
sub_GO_top <- sub_GO[c(1,4,11,14,61),]   #mPB module 6
sub_GO_top <- sub_GO[c(1,3,6,7,9),]   #PB module 2
go <- sub_GO_top[,c(2,3,11,6,9)]
colnames(go)<-c( 'ID', 'term','category','adj_pval','genes')
go <- go %>% mutate(genes = str_replace_all(genes, "/", ","))
#go <- go %>% mutate(genes = str_c('"', genes, '"'))
fc <- data.frame(cds_pr_test_res$morans_I)
fc$ID <- rownames(cds_pr_test_res)
colnames(fc)<-c("logFC",'ID')
circ <- circle_dat(go,fc)
genes_draw <- circ$genes %>% unique()
process_draw <- circ$term %>% unique()
chord <-chord_dat(data = circ, genes = genes_draw,process = process_draw)
p1=GOChord(chord, space = 0.02, gene.order = 'logFC', 
           gene.space = 0.3, gene.size = 4, border.size=0.1,
           nlfc = 0,
           ribbon.col = brewer.pal(5,"Paired"))
ggsave("/home/yushiya/data/cd34/data/fig/fig2-tra-module-GO-chord-BM-4.pdf", plot=p1, width=5, height=5.8)


#plot module genes
library(ComplexHeatmap)
library(grid)
#migration_gene <- c('CXCR4','CSF3R','CSF1','SELL','ADD2','IL1B','PF4V1','CXCL5','CD74','AIF1','ANXA1','MDK','CCL5','ITGA2B','ITGB3','CXADR','MPP1')
#ribosome_gene <- c('RPL11','RPS8','RPL5','RPS7','RPL14','RPL24','RPL35A','RPS3A','RPS23','RPS14','NPM1')
#prolif_gene <- c('LEF1','IL7R','MZB1','LST1','CARD11','RAG2','FLT3','CORO1A','TYROBP')
#adhesion_gene <- c('S100A10','SPTA1','CD36','ANGPT1')
migration_gene <- c('CSF3R','CSF1','SELL','ADD2','IL1B','PF4V1','CXCL5','CD74','AIF1','ANXA1','MDK','CCL5','ITGA2B','ITGB3','CXADR','MPP1')
ribosome_gene <- c('RPL11','RPL5','RPL24','RPS3A','RPS23','RPS14','NPM1')
prolif_gene <- c('LEF1','IL7R','MZB1','LST1','CARD11','RAG2','FLT3','CORO1A','TYROBP')
adhesion_gene <- c('S100A10','SPTA1','CD36','ANGPT1')
gene <- c(migration_gene,ribosome_gene,prolif_gene,adhesion_gene)
gene_module_df <- read.csv("/home/yushiya/data/cd34/data/fig/fig2-peu_tra_module-all.csv", row.names = 1)
#gene_module_1 <- subset(gene_module_df, subset= module==1)
type="HSC"
deg <- read.csv(paste("/home/yushiya/data/cd34/data/fig/cpstim_BMPB_",type,".csv",sep=""))
deg <- subset(deg, subset= p_val < 0.01)
module_genes_BM <- gene_module_BM %>% filter(module == "Module 3") %>% pull(id)
module_genes_mPB <- gene_module_mPB %>% filter(module == "Module 6") %>% pull(id)
module_genes_PB <- gene_module_PB %>% filter(module == "Module 2") %>% pull(id)
module_genes <- union(union(gene_module_BM$id,gene_module_mPB$id),gene_module_PB$id)
gene <- intersect(module_genes, deg$gene)
#sub1 <- subset(sub1_all, subset= tissue %in% c("BM","mPB","PB"))
sub11 <- subset(sub1, subset= celltype==type)
cellInfo <- data.frame(tissue=sub11$tissue)
mtx <- data.frame(sub11@assays[["RNA"]]@data[gene,]) 
colnames(mtx) <- rownames(cellInfo)
top_exp <- sapply(split(rownames(cellInfo), cellInfo$tissue),
                  function(cells) rowMeans(mtx[gene,cells]))
df_scaled <- t(scale(t(top_exp)))
anno_genes <- c("ANXA1","MDK","SELL","CXCR4","CXCL5","CCL5",
                #"CRHBP","AVP","AREG",
                "RPS18","RPL15","RPL9",
                'SPTA1','CD36','ANGPT1',
                'IL7R','RAG2','FLT3','CORO1A','TYROBP')
anno_genes <- c("ANXA1","MDK","SELL","CXCR4","CXCL5","CCL5","MPP1",
                "CRHBP",
                "RPL15","RPL35","RPS25",
                'IL7R','RAG2','FLT3','CORO1A','TYROBP',
                "JUN","JUND")
anno_genes <- c("CSF3R","SELL","IL1B","CD74","AIF1","ANXA1","MPP1")  #BM migration
anno_genes <- intersect(intersect(module_genes_BM, module_genes_mPB), module_genes_PB)
p11 <-pheatmap::pheatmap(df_scaled, 
                         cluster_cols=FALSE, 
                         cluster_rows=T, 
                         show_rownames=T, 
                         show_colnames=F, 
                         clustering_method = "ward.D2",
                         filename=NA,
                         border_color = NA,
                         fontsize_row = 8,
                         color=colorRampPalette(c("#4979B6","white","#D93429"))(100),
                         #annotation_colors=ann_colors,
                         #annotation_row = annotation_row,
                         #clustering_callback = callback,
                         annotation_names_col = F,
                         annotation_names_row = F)
source('/home/yushiya/code/add.flag.R')
p <- add.flag(p11,kept.labels = anno_genes,repel.degree = 0.2)
ggsave(paste("/home/yushiya/data/cd34/data/fig/fig2-DEG-mark_",type,".pdf",sep=""), plot=p, width=3.5, height=4.5)


#Module score
gene_module_df <- read.csv("/home/yushiya/data/cd34/data/fig/fig2-peu_tra_module-all.csv", row.names = 1)
m1 <- subset(gene_module_df, subset= module=="1")
pdf("/home/yushiya/data/cd34/data/fig/fig2-tra-module-umap.pdf", width=8, height=6)
plot_cells(sub_cds1,
           genes=m1$id,
           #genes=gene_module_df %>% filter(module %in% c(1:6)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
dev.off()
sub1 <- AddModuleScore(object = sub1, features = list(m1_exp=m1$id), name="m1_exp")
p1=VlnPlot(sub1, features = c("m1_exp1"), cols = color_stim, split.by="tissue", pt.size = 0)
ggsave("/home/yushiya/data/cd34/data/fig/fig2-tra-module-m1exp.pdf", plot=p1, width=5, height=2.5)
#migratory gene in module
migr_gene <- c("CXCR4","CSF3R","CSF1","SELL","ADD2","IL1B","PF4V1","PF4","PPBP","CXCL5","CD74","AIF1","ANXA1","MDK","CCL5","ITGA2B","ITGB3","CXADR","MPP1")
migr_gene <- c("CSF3R","SELL","IL1B","CD74","AIF1","ANXA1")
sub1 <- AddModuleScore(object = sub1, features = list(Migr_score=migr_gene), name="Migr_score")
VlnPlot(sub1, features = c("Migr_score1"), cols = color_stim, split.by="tissue", pt.size = 0)
#chemotaxis gene in module
chemo_gene <- c("CSF3R","CSF1","SELL","IL1B","PF4V1","CXCL5","CD74","AIF1","ANXA1","MDK","CCL5","CXADR","MPP1")
sub1 <- AddModuleScore(object = sub1, features = list(chemo_score=chemo_gene), name="chemo_score")
VlnPlot(sub1, features = c("chemo_score1"), cols = color_stim, split.by="tissue", pt.size = 0)
#plot through pseudotime
sub1 <- subset(sub1, subset = celltype %in% c("HSC","MPP#1","MPP#2","LMPP#1","LMPP#2","BMEP","EryP","MkP","Ma/Eo/BaP"))
mtx <- as.data.frame(sub1$celltype)
colnames(mtx)[1] <- "celltype"
mtx$tissue <- sub1$tissue
mtx$pseudotime <- sub1$pseudotime
mtx$CD74 <- sub1@assays[["RNA"]]@data["CD74",]
mtx$SELL <- sub1@assays[["RNA"]]@data["SELL",]
mtx$ANXA1 <- sub1@assays[["RNA"]]@data["ANXA1",]
mtx$CXCR4 <- sub1@assays[["RNA"]]@data["CXCR4",]
mtx$MDK <- sub1@assays[["RNA"]]@data["MDK",]
mtx$MPP1 <- sub1@assays[["RNA"]]@data["MPP1",]
mtx$IL1B <- sub1@assays[["RNA"]]@data["IL1B",]
mtx$AIF1 <- sub1@assays[["RNA"]]@data["AIF1",]
mtx$Migr_score <- sub1$Migr_score1
col <- c(color_BMPB_lineages,color_stim)
col <- c("#550000", "#AA3939", "#CC6F66", "#E79492", "#FBD2CE","#F2B662", "#F38F44", "#B97802","#915900","#DB5C25", "#F3B747", "#649541")
p1=ggplot(mtx, aes(x=pseudotime,y=Migr_score))+
  geom_point(size=0.3,aes(colour=celltype),alpha = 0)+
  stat_smooth(mapping = aes(colour=tissue))+
  #facet_wrap(~celltype)+#theme_bw()
  #scale_y_continuous(limits = c(-0.5,0.5))+
  scale_color_manual(values = col) +
  #scale_color_manual(values = col) +
  theme_classic()+
  theme(axis.text.x=element_text(hjust=0.5),
        text = element_text(size = 18),
        #panel.border = element_rect(color = "black", fill=NA),
        panel.background = element_blank())
#panel.grid = element_line(colour = "grey60"),
#axis.line = element_line(colour = "grey60",size = rel(1),arrow = arrow(angle = 30,length = unit(0.1,"inches"))))
ggsave(file=paste("/home/yushiya/data/cd34/data/fig/fig2-pseu-SELL.pdf.pdf",sep=""), plot=p1, width=8, height=4)


#prepare for python SCENIC and read results
library(Seurat)     
library(SCopeLoomR) 
library(AUCell)     
library(SCENIC)      
library(dplyr)      
library(KernSmooth)    
library(RColorBrewer)  
library(plotly) 
sub1 <- subset(immune.combined, subset= celltype %in% c("HSC","MPP#1","MPP#2","LMPP#1","LMPP#2","BMEP","EryP","MkP","Ma/Eo/BaP","GMP","CDP","pre-pDC","CLP","pro-B","pre-B"))
sub1 <- subset(sub1, subset= stim %in% c("BM","mPB","PB"))
#提取基因表达矩阵
#normalized_data_matrix <- GetAssayData(sub1, assay = "RNA", slot = "data")
raw_count_matrix <- GetAssayData(sub1, assay = "RNA", slot = "counts")
#过滤，保留至少在 1% 细胞中表达的基因，仅保留总体表达量高于 3% 细胞总数的基因。
#gene_expression_sums <- Matrix::colSums(raw_count_matrix)
#threshold_sum <- ncol(raw_count_matrix) * 0.03
#genes_over_threshold <- gene_expression_sums > threshold_sum
genes_expressed_in_cells <- Matrix::rowSums(raw_count_matrix > 0) > (ncol(raw_count_matrix) * 0.01)
filtered_genes <- genes_expressed_in_cells
exprMat_filtered <- raw_count_matrix[filtered_genes, ]
#生成loom
loom <- SCopeLoomR::build_loom(
  file.name = "/home/yushiya/data/cd34/data/fig/SCENIC_immume_BMmPBPB.loom",
  dgem = exprMat_filtered,
  default.embedding = NULL
)
loom$close()
#https://www.jianshu.com/p/b78687101562
#python code
#sh /home/yushiya/code/pyscenic_from_loom.sh -i SCENIC_immume_BMmPBPB.loom -n 10
#Read pyscenic results
library(optparse)
op_list <- list(
  make_option(c("-l", "--input_loom"), type = "character", default = NULL, action = "store", help = "The input of aucell loom file",metavar="rds"),
  make_option(c("-m", "--input_meta"), type = "character", default = NULL, action = "store", help = "The metadata of Seurat object",metavar="idents"),
  make_option(c("-c", "--celltype"), type = "character", default = NULL, action = "store", help = "The colname of metadata to calculate RSS",metavar="lab
el")
)
parser <- OptionParser(option_list = op_list)
opt = parse_args(parser)
#read aucell results
loom <- open_loom("/home/yushiya/code/aucell.loom")
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
close_loom(loom)
#read metadata
meta <- sub1@meta.data
cellinfo <- meta[,c("celltype","nFeature_RNA","nCount_RNA")]
colnames(cellinfo)=c('celltype', 'nGene' ,'nUMI')
cellTypes <-  as.data.frame(subset(cellinfo,select = 'celltype'))
selectedResolution <- "celltype"
#calculate regulon特异性分数(Regulon Specificity Score, RSS) in celltypes 
sub_regulonAUC <- regulonAUC
rss <- calcRSS(AUC=getAUC(sub_regulonAUC),
               cellAnnotation=cellTypes[colnames(sub_regulonAUC),
                                        selectedResolution])
rss=na.omit(rss)
try({
  rssPlot <- plotRSS(rss)
  save(regulonAUC,rssPlot,regulons,file='regulon_RSS.Rdata')
})
saveRDS(rss,paste0("/home/yushiya/data/cd34/data/fig/SCENIC_immume_BMmPBPB_rss.rds"))
rss <- readRDS("/home/yushiya/data/cd34/data/fig/SCENIC_immume_BMmPBPB_rss.rds")
plotRSS_oneSet(rss, setName = "MPP#2")
# 计算每个细胞组中各调控子(regulon)的平均活性，并将这些平均活性值存储在一个矩阵中
# cellsPerGroup这里得到是不同细胞群中的样本列表
# function(x)rowMeans(getAUC(sub_regulonAUC)[,x])可以计算每个细胞群的regulon平均AUC值
cellTypes <- data.frame(row.names = colnames(sub1), 
                        celltype = sub1$celltype)
cellsPerGroup <- split(rownames(cellTypes),cellTypes[,"celltype"])
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(x) 
                                    rowMeans(getAUC(sub_regulonAUC)[,x]))
#range(regulonActivity_byGroup)
# 对结果进行归一化
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T))
regulonActivity_byGroup_Scaled <- regulonActivity_byGroup_Scaled[,c(1:15)]
#展示转录因子平均活性(全部)
library(ComplexHeatmap)
library(circlize)
Heatmap(
  regulonActivity_byGroup_Scaled,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),
                                            rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = TRUE,
  row_names_gp                 = gpar(fontsize = 12),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE
)
#基于平均活性将得到的scaled Data进行不同组别的差异分析
library(dplyr) 
rss=regulonActivity_byGroup_Scaled
head(rss)
df = do.call(rbind,
             lapply(1:ncol(rss), function(i){
               dat= data.frame(
                 path  = rownames(rss), # 当前regulon的名称
                 cluster = colnames(rss)[i], # 当前cluster的名称
                 sd.1 = rss[,i], # 当前cluster中每个调控因子的值
                 sd.2 = apply(rss[,-i], 1, median)  #除了当前cluster之外的所有cluster 中该调控因子的中位值
               )
             }))
df$fc = df$sd.1 - df$sd.2
top5 <- df %>% 
  group_by(cluster) %>% 
  top_n(15, fc)
op5_unique <- top5[!duplicated(top5$path), ]
rowcn = data.frame(path = top5_unique$cluster) 
rowcn <- data.frame(rowcn[c(1:24),])
n = rss[top5_unique$path,] 
n <- n[c(1:24),]
breaksList = seq(-1.5, 1.5, by = 0.1)
colors <- colorRampPalette(c("#336699", "white", "tomato"))(length(breaksList))
pdf("TFs_output.pdf", width = 6, height = 10)
pheatmap(n,
         #annotation_row = rowcn,
         color = colors,
         cluster_rows = T,
         cluster_cols = FALSE,
         show_rownames = T,
         #gaps_col = cumsum(table(annCol$Type)),  # 使用排序后的列分割点
         #gaps_row = cumsum(table(annRow$Methods)), # 行分割
         fontsize_row = 10,
         fontsize_col = 12,
         annotation_names_row = FALSE)
dev.off()
#展示单个转录因子
library(SummarizedExperiment)
regulonsToPlot = c("CEBPB(+)","MAFF(+)","MECOM(+)","NFKB1(+)","NFKB2(+)","PPARD(+)")
regulonsToPlot %in% row.names(sub_regulonAUC)
sub1@meta.data = cbind(sub1@meta.data ,
                       t(assay(sub_regulonAUC[regulonsToPlot,])))

# Vis
DotPlot(sub1, features = unique(regulonsToPlot)) + RotatedAxis()
RidgePlot(sub1, features = regulonsToPlot , ncol = 3, cols=color_all) 
VlnPlot(sub1, features = regulonsToPlot, split.by = "tissue", cols=color_stim, pt.size = 0)
FeaturePlot(sub1,features = regulonsToPlot)
#through pseudotime
mtx <- as.data.frame(sub1$celltype)
colnames(mtx)[1] <- "celltype"
mtx$tissue <- sub1$tissue
mtx$pseudotime <- sub1$pseudotime
mtx$CEBPB <- sub1@meta.data[,"CEBPB(+)"]
mtx$MAFF <- sub1@meta.data[,"MAFF(+)"]
mtx$MECOM <- sub1@meta.data[,"MECOM(+)"]
mtx$NFKB1 <- sub1@meta.data[,"NFKB1(+)"]
mtx$NFKB2 <- sub1@meta.data[,"NFKB2(+)"]
mtx$PPARD <- sub1@meta.data[,"PPARD(+)"]
col <- c(color_BMPB_lineages,color_stim)
#col <- c("#550000", "#AA3939", "#CC6F66", "#E79492", "#FBD2CE","#F2B662", "#F38F44", "#B97802","#DB5C25", "#F3B747", "#649541")
p1=ggplot(mtx, aes(x=pseudotime,y=PPARD))+
  geom_point(size=0.3,aes(colour=celltype),alpha = 0)+
  stat_smooth(mapping = aes(colour=tissue))+
  #facet_wrap(~celltype)+#theme_bw()
  #scale_y_continuous(limits = c(-0.5,0.5))+
  scale_color_manual(values = col) +
  #scale_color_manual(values = col) +
  theme_classic()+
  theme(axis.text.x=element_text(hjust=0.5),
        text = element_text(size = 18),
        #panel.border = element_rect(color = "black", fill=NA),
        panel.background = element_blank())
#panel.grid = element_line(colour = "grey60"),
#axis.line = element_line(colour = "grey60",size = rel(1),arrow = arrow(angle = 30,length = unit(0.1,"inches"))))
ggsave(file=paste("/home/yushiya/data/cd34/data/fig/fig2-pseu-B-IL1B.pdf.pdf",sep=""), plot=p1, width=8, height=4)



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

#Find markers
type <- c("HSC","MPP#1","MPP#2","LMPP#1","LMPP#2","BMEP","EryP","MkP","GMP","CDP","pre-pDC","CLP","pro-B","pre-B")
type <- c("Ma/Eo/BaP")
type <- c("pre-pDC","CLP","pro-B","pre-B")
type <- c("HSC")
type <- c("ETP","pre-pDC")
for (i in 1:length(type)) {
  #sub1 <- sub1 <- subset(immune.combined, subset= celltype==type[i] & (stim=="BM"|stim=="PB"|stim=="mPB"))
  sub11 <- subset(sub1, subset= celltype==type[i])
  Idents(sub11) <- sub11$stim
  cp_stim.markers <- FindAllMarkers(sub11, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
  write.csv(cp_stim.markers, paste("/home/yushiya/data/cd34/data/fig/cpstim_mPBPB_",type[i],".csv",sep=""))
  #write.csv(cp_stim.markers, paste("/home/yushiya/data/cd34/data/fig/cpstim_BMmPB_Ma.csv",sep=""))
  #BMID <- subset(cp_stim.markers, subset= cluster=="BM")
  #BM_gene <- bitr(BMID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
  mPBID <- subset(cp_stim.markers, subset= cluster=="mPB")
  mPB_gene <- bitr(mPBID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
  PBID <- subset(cp_stim.markers, subset= cluster=="PB")
  PB_gene <- bitr(PBID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
  #SPID <- subset(cp_stim.markers, subset= cluster=="SP")
  #SP_gene <- bitr(SPID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
  #thymusID <- row.names(subset(cp_stim.markers, subset= cluster=="thymus"))
  #thymus_gene <- bitr(thymusID, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
  #cp = list(BM.gene=BM_gene$ENTREZID, mPB_gene=mPB_gene$ENTREZID,PB_gene=PB_gene$ENTREZID)
  #cp = list(BM.gene=BM_gene$ENTREZID, mPB_gene=mPB_gene$ENTREZID,PB_gene=PB_gene$ENTREZID, SP.gene=SP_gene$ENTREZID)
  cp = list(mPB_gene=mPB_gene$ENTREZID,PB_gene=PB_gene$ENTREZID)
  go.p <- compareCluster(cp,
                         fun = "enrichGO",
                         OrgDb = "org.Hs.eg.db",
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.01
  )
  saveRDS(go.p,paste("/home/yushiya/data/cd34/data/fig/GO_mPBPB_",type[i],".rds",sep=""))
  #saveRDS(go.p,paste("/home/yushiya/data/cd34/data/fig/GO_BMmPB_Ma.rds",sep=""))
}
for (i in 1:length(type)) {
  go.p <- readRDS(paste("/home/yushiya/data/cd34/data/fig/GO_BMPB_",type[i],".rds",sep=""))
  go.p <- readRDS(paste("/home/yushiya/data/cd34/data/fig/GO_BMmPB_",type[i],".rds",sep=""))
  #go.p <- readRDS(paste("/home/yushiya/data/cd34/data/fig/GO_BMPB_Ma.rds",sep=""))
  go.p1 <- simplify(go.p,cutoff=0.6,by="p.adjust",select_fun=min)  #去除冗余
  #go.p1@compareClusterResult <- subset(go.p1@compareClusterResult, subset= Cluster %in% c("BM.gene","mPB_gene","PB_gene"))
  go.p@compareClusterResult <- go.p@compareClusterResult[c(1,3,4,9,13,16,47,48,49,50,61,366),]  #BM mPB HSC
  go.p@compareClusterResult <- go.p@compareClusterResult[c(1,2,3,4,5,6,399,403,406,410,420,427),]  #BM PB BMEP
  go.p@compareClusterResult <- go.p@compareClusterResult[c(1,2,3,5,6,7,233,234,235,237,239,240),]  #BM mPB BMEP
  go.p@compareClusterResult <- go.p@compareClusterResult[c(2,3,4,5,6,8,53,55,56,59,65,68),]  #BM PB GMP
  go.p@compareClusterResult <- go.p@compareClusterResult[c(2,3,4,6,8,9,70,71,74,77,78,409),]  #BM mPB GMP
  p1=ggplot(go.p, aes(Cluster, Description), showCategory=6) +
    geom_point(aes(color=p.adjust, size=GeneRatio))+
    theme_classic()+
    theme(axis.text.x =element_text(angle = 45,vjust = 0.5,hjust = 0.5),
          text = element_text(size = 15),
          panel.background = element_blank(),
          panel.grid = element_line(colour = "grey90"))
  p1=ggplot(go.p1, aes(Cluster, Description), showCategory=6) +
    geom_point(aes(color=p.adjust, size=GeneRatio))+
    theme_classic()+
    theme(axis.text.x =element_text(angle = 45,vjust = 0.5,hjust = 0.5),
          text = element_text(size = 15),
          panel.background = element_blank(),
          panel.grid = element_line(colour = "grey90"))
  ggsave(file=paste("/home/yushiya/data/cd34/data/fig/fig2-GO_BMmPB_",type[i],".pdf",sep=""), plot=p1, width=7, height=5)
  ggsave(file=paste("/home/yushiya/data/cd34/data/fig/fig2-GO_BMPB_",type[i],".pdf",sep=""), plot=p1, width=7, height=5)
  #ggsave(file=paste("/home/yushiya/data/cd34/data/fig/fig2-GO_BMPB_Ma.pdf",sep=""), plot=p1, width=8, height=5)
}
#draw DEG heatmap and GO plots
sub1 <- subset(immune.combined, subset= tissue %in% c("BM","mPB","PB"))
type <- c("HSC","MPP#1","MPP#2","LMPP#1","LMPP#2","BMEP","EryP","MkP","GMP","CDP","CLP","pro-B","GMP")
type <- c("CLP")
type <- c("Ma/Eo/BaP")
sub11 <- subset(sub1, subset= celltype==type)
cp_stim.markers <- read.csv(paste("/home/yushiya/data/cd34/data/fig/cpstim_BMPB_",type,".csv",sep=""))
#DEG_mmt <- read.csv(paste("/home/yushiya/data/cd34/data/fig/DEG_mmt_",type,".csv",sep=""))
#DEG_mmt <- subset(DEG_mmt, subset= de_pval < 0.01)
#top_genes <- DEG_mmt %>% top_n(n = -20, wt = de_pval)
top_genes <- cp_stim.markers %>% group_by(cluster) %>% top_n(n = -12, wt = p_val)
top_genes <- cp_stim.markers %>% group_by(cluster) %>% top_n(n = 12, wt = avg_log2FC)
top_genes <- top_genes[,"gene"] %>% unique()
colnames(top_genes) <- "gene"
new_gene <- c("CSF3R","SELL","IL1B","AIF1","ANXA1")
top_genes <- rbind(top_genes,
                   data.frame(gene = new_gene))
cellInfo <- data.frame(tissue=sub11$tissue)
mtx <- data.frame(sub11@assays[["RNA"]]@data[top_genes$gene,]) 
colnames(mtx) <- rownames(cellInfo)
top_exp <- sapply(split(rownames(cellInfo), cellInfo$tissue),
                  function(cells) rowMeans(mtx[top_genes$gene,cells]))
p1=pheatmap(top_exp, cluster_cols = F, #cluster_rows = F,  
            #clustering_method = "median",
            scale = "row",
            fontsize_row = 10, fontsize_col= 12, angle_col = 0,
            col=viridis::viridis(100))
ggsave(file=paste("/home/yushiya/data/cd34/data/fig/fig2-DEG_BMmPBPB_",type,".pdf",sep=""), plot=p1, width=4.5, height=5)
ggsave(file=paste("/home/yushiya/data/cd34/data/fig/fig2-DEG_mmt_",type,".pdf",sep=""), plot=p1, width=3.5, height=5)
#proliferation id 6194/3726/7128/972/2355/3115/3133/1026/64115/58515/3113/1191/7940/604/10892/57162/3615/5621/4860/4208/1051
#differentiation id 6659/60/57492/1385/6688/6777/56998/5925/29117/861/5530/5590/6599/4192/57178/3965/6603/8289/196528/301/11221/7040/3123/9760/5591/2309/387/8546
#adhesion id 3122/972/3127/8013/1326/3115/3133/5329/64115/9308/3119/58515/3108/3113/960/84807/288/604/3117/10892/7704/2113/399/4860/10725/50848
prolif_gene <- c("6194","3726",'7128',"972","2355","3115","3133","1026","64115","58515","3113","1191","7940","604","10892","57162","3615","5621","4860","4208","1051")
prolif_gene <- bitr(prolif_gene, fromType="ENTREZID", toType=c("SYMBOL"), OrgDb="org.Hs.eg.db")
sub1 <- AddModuleScore(object = sub1, features = list(Prolif_Score=prolif_gene$SYMBOL), name="Prolif_Score")
diff_gene <- c("6659","60","57492","1385","6688","6777","56998","5925","29117","861","5530","5590","6599","4192","57178","3965","6603","8289","196528","301","11221","7040","3123","9760","5591","2309","387","8546")
diff_gene <- bitr(diff_gene, fromType="ENTREZID", toType=c("SYMBOL"), OrgDb="org.Hs.eg.db")
sub1 <- AddModuleScore(object = sub1, features = list(Diff_Score=diff_gene$SYMBOL), name="Diff_Score")
VlnPlot(sub1, features = "Diff_Score1", split.by = "tissue", pt.size = 0, cols = color_stim)
adhe_gene <- c("3122","972","3127","8013","1326","3115","3133","5329","64115","9308","3119","58515","3108","3113","960","84807","288","604","3117","10892","7704","2113","399","4860","10725","50848")
adhe_gene <- bitr(adhe_gene, fromType="ENTREZID", toType=c("SYMBOL"), OrgDb="org.Hs.eg.db")
#adhe_gene <- c("CSF3R","IL1B","CD74","ANXA1","MPP1")
sub1 <- AddModuleScore(object = sub1, features = list(Adhe_Score=adhe_gene$SYMBOL), name="Adhesion_Score")
VlnPlot(sub1, features = "Adhesion_Score1", split.by = "tissue", pt.size = 0, cols = color_stim)
#heatmap of genes
genes <- c(migr_gene,diff_gene$SYMBOL, adhe_gene$SYMBOL)
genes <- unique(genes)
genes <- genes[!grepl("^HLA", genes)]
sub11 <- subset(sub1, subset= celltype=="HSC")
cellInfo <- data.frame(tissue=sub11$tissue)
mtx <- data.frame(sub11@assays[["RNA"]]@data[genes,]) 
colnames(mtx) <- rownames(cellInfo)
top_exp <- sapply(split(rownames(cellInfo), cellInfo$tissue),
                  function(cells) rowMeans(mtx[genes,cells]))
pheatmap(top_exp, cluster_cols = F, #cluster_rows = F,  
         clustering_method = "median",
         scale = "row",
         fontsize_row = 10, fontsize_col= 12, angle_col = 0,
         col=viridis::viridis(100))
#relationships between adhesion score and differention score
df <- data.frame(sub1$tissue)
colnames(df) <- "Tissue"
df$Celltype <- sub1$celltype
df$Differentiation_Score <- sub1$Diff_Score1
df$Adhesion_Score <- sub1$adhesion_score1
df$Migration_Score <- sub1$Migr_score1
med_df <- df %>%
  group_by(Tissue, Celltype) %>%
  summarise(
    Diff_med = median(Differentiation_Score,  na.rm = TRUE),
    Migr_med = median(Migration_Score, na.rm = TRUE),
    Adh_med = median(Adhesion_Score, na.rm = TRUE),
    .groups = "drop"
  )
p1=ggplot(med_df, aes(x = Diff_med, y = Migr_med, color = Tissue)) +
  geom_point(size = 3) +
  scale_color_manual(values = color_stim) +
  labs(x = "Diff_Score median",
       y = "Migr_Score median",
       title = "Tissue-CellType median scores") +
  theme_bw()
ggsave(file=paste("/home/yushiya/data/cd34/data/fig/fig2-relation_diff_migr1.pdf",sep=""), plot=p1, width=5.5, height=3)
#difference among lineages
Idents(sub1) <- sub1$celltype
sub1 <- RenameIdents(sub1, "HSC"="HSCs","MPP#1"="HSCs","MPP#2"="HSCs","LMPP#1"="HSCs","LMPP#2"="HSCs",
                     "BMEP"="Ery/Mk","EryP"="Ery/Mk","MkP"="Ery/Mk","Ma/Eo/BaP"="Ery/Mk",
                     "GMP"="Mye","CDP"="Mye","pre-pDC"="Mye",
                     "CLP"="B","pre-B"="B","pro-B"="B")
VlnPlot(sub1, features = "ANXA1", pt.size = 0,combine = FALSE, split.by = "tissue", cols=color_stim)
VlnPlot(sub1, features = "SELL", pt.size = 0,combine = FALSE, split.by = "tissue", cols=color_stim)
VlnPlot(sub1, features = "Migr_Score1", pt.size = 0,combine = FALSE, split.by = "tissue", cols=color_stim)


#GO results of momento
type <- c("HSC","MPP#1","MPP#2","LMPP#1","LMPP#2","BMEP","EryP","MkP","CLP","pro-B","pre-B","GMP","CDP","pre-pDC")
type <- c("Ma/Eo/BaP")
for (i in 1:length(type)) {
  DEG_mmt <- read.csv(paste("/home/yushiya/data/cd34/data/fig/DEG_mmt_",type[i],".csv",sep=""))
  #DEG_mmt <- read.csv(paste("/home/yushiya/data/cd34/data/fig/DEG_mmt_Ma.csv",sep=""))
  DEG_mmt <- subset(DEG_mmt, subset= de_pval < 0.001)
  BMID <- subset(DEG_mmt, subset= de_coef > 0)
  PBID <- subset(DEG_mmt, subset= de_coef < 0)
  cp = list(BM_gene=BMID$gene, PB_gene=PBID$gene)
  go.p <- compareCluster(cp,
                         fun = "enrichGO",
                         OrgDb = "org.Hs.eg.db",
                         ont = "BP",
                         keyType = "SYMBOL",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.01
  )
  saveRDS(go.p,paste("/home/yushiya/data/cd34/data/fig/GO_mmt_",type[i],".rds",sep=""))
  #saveRDS(go.p,paste("/home/yushiya/data/cd34/data/fig/GO_mmt_Ma.rds",sep=""))
}
#View DEG and GO results of momento
type <- c("CDP")
DEG_mmt <- read.csv(paste("/home/yushiya/data/cd34/data/fig/DEG_mmt_",type,".csv",sep=""))
DEG_mmt <- read.csv(paste("/home/yushiya/data/cd34/data/fig/DEG_mmt_allBMPB.csv",sep=""))
DEG_mmt <- subset(DEG_mmt, subset= de_pval < 0.001)
go <- readRDS(paste("/home/yushiya/data/cd34/data/fig/GO_mmt_",type,".rds",sep=""))
go <- readRDS(paste("/home/yushiya/data/cd34/data/fig/GO_",type,".rds",sep=""))
go.p1 <- simplify(go,cutoff=0.5,by="p.adjust",select_fun=min)
p1=ggplot(go, aes(Cluster, Description)) +
  geom_point(aes(color=p.adjust, size=GeneRatio))+
  theme_classic()+
  theme(axis.text.x =element_text(angle = 45,vjust = 0.5,hjust = 0.5),
        text = element_text(size = 15),
        panel.background = element_blank(),
        panel.grid = element_line(colour = "grey90"))#+theme_bw()
ggsave(file=paste("/home/yushiya/data/cd34/data/fig/fig2-GO_mmt_",type,".pdf",sep=""), plot=p1, width=7, height=5)
View(go@compareClusterResult)
#DEG plots of all cells seperated by GO terms
#ribosome biogenesis,MRTO4/RPL11/ZNF593/EBNA1BP2/ERI3/RPS8/ZNHIT6/RPL5/RPS27/RRP15/RPS7/NOL10/WDR43/XPO1/NIFK/WDR12/XRCC5/RPL14/DHX30/NSUN3/RPL24/RPL35A/RPS3A/TMA16/FRG1/DROSHA/MTREX/NSA2/RPS23/WDR36/NHP2/RPL7L1/RPF2/RPS12/TFB1M/GTF2H5/PRKDC/RBIS/RCL1/RPL7A/MRPS2/RPP38/RPS24/EXOSC1/PDCD11/RPS13/METTL15/WDR74/DDX10/RPL6/RPLP0/POP5/GTF3A/EXOSC8/CUL4A/CINP/NOP10/RPS27L/EFL1/TRAF7/RSL1D1/RPS15A/EXOSC6/METTL16/NUP88/C1QBP/RPL26/UTP6/AATF/RPL38/RBFA/WDR18/RPS15/RPS28/RPS19/NOP53/PIH1D1/RPS5/NOP56/XRN2/RPS21/DDX17/GNL3L/PIN4/RPL10/DKC1
#regulation of hemopoiesis,TNFSF4/ID2/CTNNB1/CD74/CD83/HSPA1A/HSPA1B/HLA-DRA/KLF10/VSIR/CTR9/ZBTB16/UBASH3B/ETS1/TSC22D1/PNP/NFKBIA/FOS/BATF/SOCS1/EVI2B/MALT1/JUNB/NFKBID/ZFP36/CEBPB
#leukocyte cell-cell adhesion,LAPTM5/F11R/TNFSF4/PELI1/SELENOK/CBLB/CD74/CD83/HLA-E/HLA-DRA/HLA-DRB5/HLA-DQB1/HLA-DMA/EZR/MAP3K8/VSIR/ZBTB16/ETS1/PNP/SOCS1/NFAT5/MALT1/NFKBID/PRNP/CEBPB
#regulation of leukocyte differentiation,TNFSF4/ID2/CTNNB1/CD74/CD83/HLA-DRA/KLF10/VSIR/ZBTB16/UBASH3B/PNP/FOS/BATF/SOCS1/EVI2B/MALT1/JUNB/NFKBID/CEBPB
#regulation of cytokine production,CD84/NLRP3/EIF2AK3/RFTN1/IL1RAP/RBM47/HMGB2/PDE4D/IFNGR1/IRF8/SLC7A5/CYBA/BRCA1/MALT1/AZU1/ELANE/TYROBP/CEBPB/IL17RA
#leukocyte migration,TNFRSF14/CDC42/CSF3R/S100A9/FCER1G/LYST/ADAM17/SOS1/IL1B/ITGA4/MYD88/OXSR1/FER/CSF1R/CD74/STK10/DUSP1/NEDD9/AIF1/GPSM3/TREM1/FYN/RAC1/MYO1G/RABGEF1/LYN/PLEC/DOCK8/B4GALT1/ANXA1/NINJ1/CAMK1D/ITGB1/ADAM8/SPI1/JAML/WNK1/EPS8/SELPLG/P2RX4/GAS6/LGALS3/RIN3/ADAM10/LGALS9/CD300A/CNN2/ICAM1/HCK/ITGB2/MSN
gene <- c("MRTO4","RPL11","ZNF593","EBNA1BP2","ERI3","RPS8","ZNHIT6","RPL5","RPS27","RRP15","RPS7","NOL10",'WDR43',"XPO1","NIFK","WDR12","XRCC5",'RPL14',
          "TNFSF4","CTNNB1","CD83","HSPA1A","HSPA1B","KLF10","VSIR",'CTR9',"ZBTB16","UBASH3B","ETS1","TSC22D1","PNP",'NFKBIA',"FOS","BATF",'SOCS1',"EVI2B",'MALT1',"JUNB","NFKBID",'ZFP36',"CEBPB",
          "LAPTM5","F11R","TNFSF4","PELI1","SELENOK","CBLB","CD74","CD83","EZR",'MAP3K8',"VSIR","ZBTB16",'ETS1',"SOCS1","NFAT5","MALT1","PRNP","CEBPB",
          "CD84","NLRP3","EIF2AK3","RFTN1","IL1RAP",'RBM47',"HMGB2","PDE4D",'IFNGR1','IRF8',"SLC7A5","CYBA","BRCA1","MALT1","AZU1","ELANE","TYROBP","IL17RA",
          "CD74","STK10","DUSP1","NEDD9","AIF1","GPSM3","TREM1",'FYN','RAC1',"MYO1G","RABGEF1","LYN",'PLEC',"DOCK8","B4GALT1","ANXA1","NINJ1","CAMK1D","ITGB1","ADAM8","SPI1","JAML","EPS8","SELPLG","WNK1","P2RX4",'GAS6','LGALS3',"RIN3")
gene <- c(#"MRTO4","RPL11","ZNF593","EBNA1BP2","ERI3","RPS8","ZNHIT6","RPL5","RPS27","RRP15","RPS7","NOL10",'WDR43',"XPO1","NIFK","WDR12","XRCC5",'RPL14',
  "CD83","HSPA1A","HSPA1B","KLF10",'NFKBIA',"FOS","EVI2B","JUNB",'ZFP36',
  "TNFSF4","PELI1","SELENOK","SOCS1","NFAT5",
  "IL1RAP","HMGB2","BRCA1",'RBM47','IRF8',"SLC7A5","CYBA",
  "CD74","STK10","DUSP1","AIF1","GPSM3",'PLEC',"DOCK8","NINJ1","ADAM8","WNK1","P2RX4","ANXA1","CAMK1D","EPS8","SELPLG",'RAC1',"ITGB1")
sub1$stim_celltype <- paste(sub1$stim,sub1$celltype,sep="_")
cellInfo <- data.frame(stim_celltype=sub1$stim_celltype)
mtx <- data.frame(sub1@assays[["RNA"]]@data[gene,]) 
colnames(mtx) <- rownames(cellInfo)
top_exp <- sapply(split(rownames(cellInfo), cellInfo$stim_celltype),
                  function(cells) rowMeans(mtx[gene,cells]))
lev<-c("BM_HSC","BM_MPP#1","BM_MPP#2","BM_LMPP#1","BM_LMPP#2","BM_BMEP","BM_EryP","BM_MkP","BM_Ma/Eo/BaP","BM_CLP","BM_pro-B","BM_pre-B","BM_GMP","BM_CDP","BM_pre-pDC",
       "PB_HSC","PB_MPP#1","PB_MPP#2","PB_LMPP#1","PB_LMPP#2","PB_BMEP","PB_EryP","PB_MkP","PB_Ma/Eo/BaP","PB_CLP","PB_pro-B","PB_pre-B","PB_GMP","PB_CDP","PB_pre-pDC")
lev<-c("BM_HSC","PB_HSC","BM_MPP#1","PB_MPP#1","BM_MPP#2","PB_MPP#2","BM_LMPP#1","PB_LMPP#1","BM_LMPP#2","PB_LMPP#2",
       "BM_BMEP","PB_BMEP","BM_EryP","PB_EryP","BM_MkP","PB_MkP","BM_Ma/Eo/BaP","PB_Ma/Eo/BaP",
       "PB_pre-B","BM_GMP","PB_GMP","BM_CDP","PB_CDP","BM_pre-pDC","PB_pre-pDC","BM_CLP","PB_CLP","BM_pro-B","PB_pro-B","BM_pre-B")
top_exp <- top_exp[,lev] 
annotation_col <- data.frame(Tissue = factor(rep(c("BM","PB"), c(15))))
annotation_col <- data.frame(Tissue = factor(rep(c("BM","PB"), c(15,15))))
rownames(annotation_col) <- lev
annotation_colors =list(Tissue=c(BM="#DB5C25",PB="#649541"))
p1=pheatmap(top_exp, cluster_cols = F, cluster_rows = F,  
            clustering_method = "average",
            annotation_col = annotation_col,
            annotation_colors = annotation_colors,
            scale = "row",
            #border=F,
            col=colorRampPalette(c("navy","white" ,"firebrick3"))(50),
            fontsize_row = 6, fontsize_col= 12, angle_col = 45)
ggsave(file="/home/yushiya/data/cd34/data/fig/fig2-allHeatmap-1.pdf", plot=p1, width=9, height=4)


#GO plots
type <- c("HSC","MPP#1","MPP#2","LMPP#1","LMPP#2","BMEP","EryP","MkP","GMP","CDP","CLP","pro-B","GMP")
type <- c("GMP")
go <- readRDS(paste("/home/yushiya/data/cd34/data/fig/GO_",type,".rds",sep=""))
View(go@compareClusterResult)
go_mtx <- go@compareClusterResult
go_mtx$celltype <- type
go_mtx <- go_mtx[c(3,17,23,25),]  #"MPP#1"
go_mtx <- go_mtx[c(6,9),]  #"MPP#2"
go_mtx <- go_mtx[c(5,7,8,12),]  #"LMPP#1"
go_mtx <- go_mtx[c(29,37,46),]  #"LMPP#2"
go_mtx <- go_mtx[c(69,95,105,116,118),]  #"BMEP"
go_mtx <- go_mtx[c(137,145,149,155,156,201,202),]  #"EryP"
go_mtx <- go_mtx[c(8,29,112,113,115,125,131,141),]  #"GMP"
go_combined <- data.frame()
go_combined <- rbind(go_combined,go_mtx)
migrID <- c(301,4267,57118,23236,11151,928,1843,4282,972,351,4192,566,5657,1991,5478,5479,811,1130,3146,3146,1511,708,6280,30817,6279,5734,2920,3553)
migr_gene <- bitr(migrID, fromType="ENTREZID", toType=c("SYMBOL"), OrgDb="org.Hs.eg.db")
immune.combined<-AddModuleScore(object = immune.combined, features = list(Migr_Score=migr_gene$SYMBOL), name="Migr_Score")
adhID <- c(284,6659,3123,8631,1021,972,3122,960,3117,11151,91663,6281,4192,3119,6280,6279,960,3113,3659,301,9314,3553)
6280/6659/3123/6279/960/3119/3113/3659/301/9314/3553/3117
adh_gene <- bitr(adhID, fromType="ENTREZID", toType=c("SYMBOL"), OrgDb="org.Hs.eg.db")
immune.combined<-AddModuleScore(object = immune.combined, features = list(Adh_Score=adh_gene$SYMBOL), name="Adhesion_Score")
VlnPlot(sub1, features = "Adhesion_Score1", split.by = "tissue", pt.size = 0)
VlnPlot(immune.combined, features = "Migr_Score1", split.by = "tissue", pt.size = 0)
VlnPlot(immune.combined, features = "Migratory_Score1", split.by = "tissue", pt.size = 0)


#stemness by cytotrace
library(CytoTRACE2)
expression_data <- immune.combined[["RNA"]]@counts
cytotrace2_result <- cytotrace2(expression_data, species = "human")
write.csv(cytotrace2_result, file=paste("/home/yushiya/data/cd34/data/fig/cytotrace2_result.csv",sep=""))
cytotrace2_result <- read.csv("/home/yushiya/data/cd34/data/fig/cytotrace2_result.csv")
immune.combined$CytoTRACE2_Score <- cytotrace2_result$CytoTRACE2_Score
immune.combined$CytoTRACE2_Potency <- cytotrace2_result$CytoTRACE2_Potency
VlnPlot(immune.combined, features = "CytoTRACE2_Score", split.by = "tissue", pt.size = 0)
DimPlot(immune.combined, group.by = "CytoTRACE2_Potency")
p1=VlnPlot(immune.combined, features = c("CytoTRACE2_Score"), cols = color_all, pt.size = 0) +
  geom_boxplot(width=.2,col="black",fill="white")+ NoLegend()
ggsave(file="/home/yushiya/data/cd34/data/fig/fig2-cytotrace.pdf", plot=p1, width=8, height=3.5)
#boxplot in different facet
mtx <- as.data.frame(colnames(immune.combined))
mtx$tissue <- immune.combined$tissue
mtx$celltype <- immune.combined$celltype
mtx$CytoTRACE2_Score <- immune.combined$CytoTRACE2_Score
mtx$ccat_Score <- immune.combined$ccat
mtx$HSC_Score <- immune.combined$HSC_Score1
mtx <- subset(mtx, subset= tissue %in% c("BM","mPB","PB"))
mtx <- subset(mtx, subset= celltype %in% c("HSC","MPP#1","MPP#2","LMPP#1","LMPP#2","BMEP","EryP","MkP","Ma/Eo/BaP","GMP","CDP","pre-pDC","CLP","pro-B","pre-B","ETP-like"))
mtx$tissue <- factor(mtx$tissue, levels=c("BM","mPB","PB"))
names(mtx)[1]='ID'
#plot with different facet
p1=ggboxplot(mtx, x="tissue", y="HSC_Score", color = "tissue",
             short.panel.labs = T, ncol=7)+
  facet_wrap(~celltype, ncol=5) +
  scale_color_manual(values = color_stim) +
  #ggtitle(clus)+
  theme(axis.text.x =element_text(angle = 45,vjust = 1,hjust = 1),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(color = "black", fill=NA))
ggsave(file=paste("/home/yushiya/data/cd34/data/fig/fig2-cytotrace2.pdf.pdf",sep=""), plot=p1, width=6, height=7)
ggsave(file=paste("/home/yushiya/data/cd34/data/fig/fig2-HSC-score-2.pdf.pdf",sep=""), plot=p1, width=5.5, height=7)

allBM <- readRDS("/home/yushiya/data/cd34/data/2024_BM_CODEX/GSE253355_Normal_Bone_Marrow_Atlas_Seurat_SB_v2.rds")
DimPlot(allBM)
expression_data1 <- allBM[["RNA"]]@counts
cytotrace2_result1 <- cytotrace2(expression_data1, species = "human")
allBM$CytoTRACE2_Score <- cytotrace2_result1$CytoTRACE2_Score
VlnPlot(allBM, features = "CytoTRACE2_Score", pt.size = 0)
write.csv(cytotrace2_result1, file=paste("/home/yushiya/data/cd34/data/2024_BM_CODEX/cytotrace2_result.csv",sep=""))

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
ccat.v <- read.csv("/home/yushiya/data/cd34/data/fig/ccat_result.csv")
immune.combined$ccat <- ccat.v$x
p1=VlnPlot(immune.combined, features = c("ccat"), pt.size = 0) +
  geom_boxplot(width=.2,col="black",fill="white")+ NoLegend()
ggsave(file="/home/yushiya/data/cd34/data/fig/fig2-ccat.pdf", plot=p1, width=8, height=3.5)
VlnPlot(immune.combined, features = "ccat", split.by = "tissue", pt.size = 0)


#Migratory score
library(ggpubr)
positive_regu_selected <- c("MPP1","MSN","CD99","RAC2","MAPK1","APP","PLCB1","ZNF580","TGFB1","CALR",
                            "TPGS1","C1QBP","CRK","PYCARD","MAPK3","SPN","ADAM10","RIN3","NCKAP1L","DNM1L",
                            "CD9","WNK1","VEGFB","DNAJC4","SPI1","MDK","RHOG","CD81","ADAM8",
                            "CAMK1D","ANXA1","DOCK8","LYN","TGS1","HOXA7","GPSM3","RIPOR2","CCL28","PTGER4",
                            "RHOH","CD47","SELENOK","INPP5D","ITGA4","ANXA4","RTN4","ADAM17","GCSAML","MIA3")
immune.combined<-AddModuleScore(object = immune.combined, features = list(Migratory_Score=positive_regu_selected), name="Migratory_Score")
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
miGene$tissue <- sub1$tissue
miGene$celltype <- sub1$celltype
#miGene$stim <- factor(miGene$stim, levels=c("BM","mPB","PB","thymus"))
miGene$stim <- factor(miGene$tissue, levels=c("BM","mPB","PB"))
names(miGene)[1]='ID'
miGene$gene='Migratory_Score1'
miGene$Migratory_Score <- sub1$Migratory_Score1
miGene$Adhesion_Score <- sub1$Adhesion_Score1
miGene$Migratory_Score <- sub1$Migr_score1
miGene$Chemotaxis_Score <- sub1$chemo_score1
#miGene$gene='MDK'
#miGene$Expression <- sub1@assays[["RNA"]]@data["MDK",]
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
p1=ggboxplot(miGene, x="celltype", y="Migratory_Score", color="tissue", bxp.errorbar = T,
             add="boxplot", short.panel.labs = TRUE, ncol=4, add.params=list(size=0.2))+
  stat_compare_means(label = "p.signif",method = "t.test",ref.group = ".all.")+
  scale_color_manual(values = color_stim) +
  theme(axis.text.x =element_text(angle = 45,hjust = 1),
        plot.title = element_text(hjust = 0))
ggsave(file="/home/yushiya/data/cd34/data/fig/fig3-chemo-BMmPBPB.pdf", plot=p1, width=9, height=3.5)
ggsave(file="/home/yushiya/data/cd34/data/fig/fig3-migr-BMPB.pdf", plot=p1, width=7, height=3)
p1=ggboxplot(miGene, x="celltype", y="Migratory_Score", color="celltype", bxp.errorbar = T,
             add="boxplot", short.panel.labs = TRUE, ncol=4, add.params=list(size=0.2))+
  scale_color_manual(values = color_T) +
  stat_compare_means(label = "p.signif",method = "t.test",ref.group = "TSP_BM") +
  theme(axis.text.x =element_text(angle = 45,hjust = 1),
        plot.title = element_text(hjust = 0))
ggsave(file="/home/yushiya/data/cd34/data/fig/fig5-migr-ETP3.pdf", plot=p1, width=4, height=4)
#boxplot with lineages in the same facet not split by stim
p1=ggboxplot(miGene, x="celltype", y="Migratory_Score", color="celltype", bxp.errorbar = T,
             add="boxplot", short.panel.labs = TRUE, ncol=4, add.params=list(size=0.2))+
  #stat_compare_means(comparisons = list(c("HSC","MPP #1"),c("MPP #1","MPP #2"),c("MPP #2","BMEP"),c("BMEP","EryP"),c("BMEP","MkP")),
  #stat_compare_means(comparisons = list(c("HSC","MPP #1"),c("MPP #1","MPP #2"),c("MPP #2","LMPP #2"),c("LMPP #2","pre-pro B"),c("LMPP #2","pro-B"),c("pre-pro B","pro-B"),c("pro-B","pre-B")), 
  #stat_compare_means(comparisons = list(c("HSC","MPP #1"),c("MPP #1","MPP #2"),c("MPP #2","LMPP #1"),c("LMPP #1","ETP"),c("ETP","mThy"),c("mThy","cThy #1"),c("mThy","cThy #2"),c("cThy #1","cThy #2")),  
  #stat_compare_means(comparisons = list(c("HSC","MPP #1"),c("MPP #1","MPP #2"),c("MPP #2","LMPP #1"),c("LMPP #1","GMP"),c("GMP","pre-cDC"),c("GMP","pre-pDC"),c("pre-cDC","pre-pDC")),  
  #                   label = "p.signif",method = "t.test")+
  scale_color_manual(values = color_all) +
  theme(axis.text.x =element_text(angle = 45,hjust = 1),
        plot.title = element_text(hjust = 0))
ggsave(file="/home/yushiya/data/cd34/data/cd34/fig/fig3-migr-BMPB-lineages-all.pdf", plot=p1, width=9, height=4)


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
#library(readxl)   #read xlsx files
celltype_markers <- read.csv("/home/yushiya/data/cd34/data/fig/DEG_celltype.csv",stringsAsFactors = F)
LT_HSC_genes <- read.csv("/home/yushiya/data/cd34/data/fig/LT_HSC_genes.csv")
HSC_gene <- subset(celltype_markers, subset= cluster=="HSC")
HSC_gene <- HSC_gene[HSC_gene$gene %in% LT_HSC_genes$Gene.Symbol,]
HSC_gene <- HSC_gene %>% top_n(n = 50, wt = avg_log2FC)
immune.combined<-AddModuleScore(object = immune.combined, features = list(HSC_Score=HSC_gene$gene), name="HSC_Score")
p1=VlnPlot(immune.combined, features = c("HSC_Score1"), cols=color_all, pt.size = 0) +
  geom_boxplot(width=.2,col="black",fill="white")+ NoLegend()
ggsave(file="/home/yushiya/data/cd34/data/fig/fig2-Score-HSC.pdf", plot=p1, width=8, height=3.5)
EryP_gene <- subset(celltype_markers, subset= cluster=="EryP")
EryP_gene <- EryP_gene %>% top_n(n = 50, wt = avg_log2FC)
immune.combined<-AddModuleScore(object = immune.combined, features = list(EryP_Score=EryP_gene$gene), name="EryP_Score")
p1=VlnPlot(immune.combined, features = c("EryP_Score1"), cols=color_all, pt.size = 0) +
  geom_boxplot(width=.2,col="black",fill="white")+ NoLegend()
ggsave(file="/home/yushiya/data/cd34/data/fig/fig2-Score-EryP.pdf", plot=p1, width=8, height=3.5)
MkP_gene <- subset(celltype_markers, subset= cluster=="MkP")
MkP_gene <- MkP_gene %>% top_n(n = 50, wt = avg_log2FC)
immune.combined<-AddModuleScore(object = immune.combined, features = list(MkP_Score=MkP_gene$gene), name="MkP_Score")
p1=VlnPlot(immune.combined, features = c("MkP_Score1"), cols=color_all, pt.size = 0) +
  geom_boxplot(width=.2,col="black",fill="white")+ NoLegend()
ggsave(file="/home/yushiya/data/cd34/data/fig/fig2-Score-MkP.pdf", plot=p1, width=8, height=3.5)
My_gene <- subset(celltype_markers, subset= cluster=="CDP" | cluster=="pre-pDC")
My_gene <- My_gene %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
immune.combined<-AddModuleScore(object = immune.combined, features = list(My_Score=My_gene$gene), name="My_Score")
p1=VlnPlot(immune.combined, features = c("My_Score1"), cols=color_all, pt.size = 0) +
  geom_boxplot(width=.2,col="black",fill="white")+ NoLegend()
ggsave(file="/home/yushiya/data/cd34/data/fig/fig2-Score-My.pdf", plot=p1, width=8, height=3.5)
B_gene <- subset(celltype_markers, subset= cluster=="pro-B"|cluster=="pre-B")
B_gene <- B_gene %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
immune.combined<-AddModuleScore(object = immune.combined, features = list(B_Score=B_gene$gene), name="B_Score")
p1=VlnPlot(immune.combined, features = c("B_Score1"), cols=color_all, pt.size = 0) +
  geom_boxplot(width=.2,col="black",fill="white")+ NoLegend()
ggsave(file="/home/yushiya/data/cd34/data/fig/fig2-Score-B.pdf", plot=p1, width=8, height=3.5)
T_gene1 <- subset(celltype_markers, subset= cluster=="Thy#1"| cluster=="Thy#3")
T_gene1 <- T_gene1 %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
T_gene2 <- subset(celltype_markers, subset= cluster=="Thy#2")
T_gene2 <- T_gene2 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
T_gene <- rbind(T_gene1,T_gene2) %>% unique()
immune.combined<-AddModuleScore(object = immune.combined, features = list(T_Score=T_gene$gene), name="T_Score")
p1=VlnPlot(immune.combined, features = c("T_Score1"), cols=color_all, pt.size = 0) +
  geom_boxplot(width=.2,col="black",fill="white")+ NoLegend()
ggsave(file="/home/yushiya/data/cd34/data/fig/fig2-Score-T.pdf", plot=p1, width=8, height=3.5)
#plot scores in 3D figure
library(scatterplot3d)
mtx <- as.data.frame(immune.combined$celltype)
colnames(mtx)[1] <- "celltype"
mtx$tissue <- immune.combined$tissue
mtx$HSC_Score <- immune.combined$HSC_Score1
mtx$EryP_Score <- immune.combined$EryP_Score1
mtx$MkP_Score <- immune.combined$MkP_Score1
mtx$B_Score <- immune.combined$B_Score1
mtx$T_Score <- immune.combined$T_Score1
mtx$My_Score <- immune.combined$My_Score1
mtx$CytoTRACE2_Score <- immune.combined$CytoTRACE2_Score
mtx$ccat <- immune.combined$ccat
mtx1 <- subset(mtx, subset= celltype=="HSC")
scatterplot3d(mtx1$My_gene,mtx1$EryP_gene,mtx1$B_gene,
              xlim = c(-0.5, 0.5),ylim = c(-0.5, 0.5),zlim = c(-0.5, 0.5),
              pch=16,  highlight.3d=TRUE)
mtx <- subset(mtx, subset= tissue %in% c("BM","mPB","PB"))
mtx <- subset(mtx, subset= celltype %in% c("HSC","MPP#1","MPP#2","LMPP#1","LMPP#2","BMEP","EryP","MkP","Ma/Eo/BaP","GMP","CDP","pre-pDC","CLP","pro-B","pre-B","ETP-like"))
mtx$tissue <- factor(mtx$tissue, levels=c("BM","mPB","PB"))
#plot with different facet
p1=ggboxplot(mtx, x="tissue", y="My_Score", color = "tissue",
             short.panel.labs = T, ncol=7)+
  facet_wrap(~celltype, ncol=5) +
  scale_color_manual(values = color_stim) +
  #ggtitle(clus)+
  theme(axis.text.x =element_text(angle = 45,vjust = 0.5,hjust = 0.5),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(color = "black", fill=NA))
ggsave(file=paste("/home/yushiya/data/cd34/data/fig/fig2-Score-My1.pdf.pdf",sep=""), plot=p1, width=6, height=7)
#plot through pseudotime
mtx$pseudotime <- immune.combined$pseudotime
col <- c(color_all,color_stim)
p1=ggplot(mtx, aes(x=pseudotime,y=HSC_Score))+
  geom_point(size=0.3,aes(colour=celltype),alpha = 0)+
  stat_smooth(mapping = aes(colour=tissue))+
  #facet_wrap(~celltype)+#theme_bw()
  #scale_y_continuous(limits = c(-0.5,0.5))+
  #scale_color_manual(values = color_stim) +
  scale_color_manual(values = col) +
  theme_classic()+
  theme(axis.text.x=element_text(hjust=0.5),
        text = element_text(size = 18),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill=NA))
#panel.grid = element_line(colour = "grey60"),
#axis.line = element_line(colour = "grey60",size = rel(1),arrow = arrow(angle = 30,length = unit(0.1,"inches"))))
ggsave(file=paste("/home/yushiya/data/cd34/data/fig/fig2-pseuScore-HSC-1.pdf.pdf",sep=""), plot=p1, width=7, height=6, dpi=100)

#plot gene expression through pseudotime
mtx <- as.data.frame(immune.combined$pseudotime)
colnames(mtx)[1] <- "pseudotime"
#mtx$GATA2 <- sub1@assays[["RNA"]]@data["GATA2",]
#mtx$VPREB1 <- sub1@assays[["RNA"]]@data["VPREB1",]
#mtx$MPO <- sub1@assays[["RNA"]]@data["MPO",]
#mtx$CD3D <- sub1@assays[["RNA"]]@data["CD3D",]
mtx$tissue <- immune.combined$tissue
mtx$HSC_Score <- immune.combined$HSC_Score1
mtx$EryP_Score <- immune.combined$EryP_Score1
mtx$MkP_Score <- immune.combined$MkP_Score1
mtx$B_Score <- immune.combined$B_Score1
mtx$T_Score <- immune.combined$T_Score1
mtx$Mye_Score <- immune.combined$My_Score1
p1=ggplot(mtx, aes(x=pseudotime,y=Mye_Score))+
  geom_point(size=0.3,aes(colour=tissue),alpha = 0)+
  stat_smooth(mapping = aes(colour=tissue))+
  #facet_wrap(~celltype)+#theme_bw()
  #scale_y_continuous(limits = c(0,1.5))+
  scale_color_manual(values = color_stim) +
  #scale_color_manual(values = col) +
  theme_classic()+
  theme(axis.text.x=element_text(hjust=0.5),
        text = element_text(size = 18),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill=NA))
#color_stim <- c("BM"="#DB5C25","mPB"="#F3B747","PB"="#649541","SP"="#AF86BA","thymus"="#4C82C5")
ggsave(file=paste("/home/yushiya/data/cd34/data/fig/fig1-pseuScore-HSC_Score.pdf",sep=""), plot=p1, width=7, height=6, dpi=100)
ggsave(file=paste("/home/yushiya/data/cd34/data/fig/fig1-pseuScore-Ery_Score.pdf",sep=""), plot=p1, width=7, height=6, dpi=100)
ggsave(file=paste("/home/yushiya/data/cd34/data/fig/fig1-pseuScore-Mk_Score.pdf",sep=""), plot=p1, width=7, height=6, dpi=100)
ggsave(file=paste("/home/yushiya/data/cd34/data/fig/fig1-pseuScore-Mye_Score.pdf",sep=""), plot=p1, width=7, height=6, dpi=100)
ggsave(file=paste("/home/yushiya/data/cd34/data/fig/fig1-pseuScore-B_Score.pdf",sep=""), plot=p1, width=7, height=6, dpi=100)
ggsave(file=paste("/home/yushiya/data/cd34/data/fig/fig1-pseuScore-T_Score.pdf",sep=""), plot=p1, width=7, height=6, dpi=100)


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
p1=FeaturePlot(sub1, features = c("CXCR4"), max.cutoff = 3, cols = c("#C9E9FF","#981C12"), split.by="stim",order=T, raster = T, raster.dpi = c(200,200)) + 
  theme(legend.position = "right")
ggsave(file="/home/yushiya/data/cd34/data/fig/fig2-CXCR4-fea.pdf", plot=p1, width=10, height=3)
p1=FeaturePlot(sub1, features = c("ANXA1"), max.cutoff = 3, cols = c("#C9E9FF","#C9E9FF","#981C12"), split.by="stim", raster = T, raster.dpi = c(200,200)) + 
  theme(legend.position = "right")
ggsave(file="/home/yushiya/data/cd34/data/fig/fig2-ANXA1-fea.pdf", plot=p1, width=10, height=3)
p1=FeaturePlot(sub1, features = c("CRHBP"), max.cutoff = 3, cols = c("#C9E9FF","#C9E9FF","#981C12"), split.by="stim",order=T, raster = T, raster.dpi = c(200,200)) + 
  theme(legend.position = "right")
ggsave(file="/home/yushiya/data/cd34/data/fig/fig2-CRHBP-fea.pdf", plot=p1, width=10, height=3)
p1=FeaturePlot(sub1, features = c("SELL"), max.cutoff = 3, cols = c("#C9E9FF","#C9E9FF","#981C12"), split.by="stim",order=T, raster = T, raster.dpi = c(200,200)) + 
  theme(legend.position = "right")
ggsave(file="/home/yushiya/data/cd34/data/fig/fig2-SELL-fea.pdf", plot=p1, width=10, height=3)
p1=FeaturePlot(sub1, features = c("RPL31"), max.cutoff = 3, cols = c("#C9E9FF","#C9E9FF","#C9E9FF","#C9E9FF","#C9E9FF","#C9E9FF","#C9E9FF","#C9E9FF","#981C12"), split.by="stim", raster = T, raster.dpi = c(200,200)) + 
  theme(legend.position = "right")
ggsave(file="/home/yushiya/data/cd34/data/fig/fig2-RPL-fea.pdf", plot=p1, width=10, height=3)
VlnPlot(sub1, features = "CXCR4", pt.size = 0,combine = FALSE, split.by="stim", cols=color_stim)
ggsave(file="/home/yushiya/data/cd34/data/fig/fig3-cellchat-CXCXR4.pdf", plot=p1[[1]], width=8, height=2.5)
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
#识别过度表达的配体受体相互作用
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
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # “netP”:推断出的信号通路的细胞间通信网络
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
spatial_cellchat <- netAnalysis_computeCentrality(spatial_cellchat, slot.name = "netP") # “netP”:推断出的信号通路的细胞间通信网络
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


#Analyze ETP
sub1 <- readRDS("/home/yushiya/data/cd34/data/fig/ETP_mono_2.rds")
sub1 <- subset(immune.combined, subset= celltype=="ETP-like")
sub1 <- subset(sub1, subset= tissue %in% c('BM','mPB','PB','thymus'))
sub1$stim <- factor(sub1$stim,levels = c('BM','mPB','PB','thymus'))
Idents(sub1) <- sub1$stim
sub1 <- RenameIdents(sub1, "BM"="TSP_BM", "mPB"="TSP_mPB", "PB"="TSP_PB","thymus"="ETP")
color_T <- c("#DB5C25","#F3B747","#649541","#4C82C5","#A6D96A","#8DA0CB","#E78AC3","#66C2A5","#FFD92F","#E5C494","#FC8D62")
color_T1 <- c("#4C82C5","#A6D96A")
DefaultAssay(sub1) <- "integrated"
#calculate cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
sub1 <- CellCycleScoring(sub1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
sub1$Phase <- factor(sub1$Phase, levels=c("G1","S","G2M"))
#delete cell cycle influence
sub1 <- FindVariableFeatures(sub1, selection.method = "vst")
#sub1$CC.Difference <- sub1$S.Score - sub1$G2M.Score
#sub1 <- ScaleData(sub1, vars.to.regress = "CC.Difference", features = rownames(sub1))
sub1 <- ScaleData(sub1, 
                  vars.to.regress = c("S.Score", "G2M.Score"), 
                  features = rownames(sub1))
sub1 <- RunPCA(sub1, npcs = 50, verbose = FALSE)
sub1 <- RunPCA(sub1, features = c(s.genes, g2m.genes))
DimPlot(sub1,group.by = "Phase",reduction = "pca",
        cols = pal_npg("nrc", alpha = 0.7)(3))
#Run umap
sub1 <- RunUMAP(sub1, reduction = "pca", dims = 1:16)
sub1 <- FindNeighbors(sub1, reduction = "pca", dims = 1:16)
sub1 <- FindClusters(sub1, resolution = 0.43)
DimPlot(sub1, reduction = "umap", label = TRUE)
DimPlot(sub1, reduction = "umap", label = TRUE, group.by = "tissue")
DimPlot(sub1, reduction = "umap", label = TRUE, group.by = "Phase")
DimPlot(sub1, reduction = "umap", label = TRUE, split.by = "integrated_snn_res.0.5", ncol=3)
p1=DimPlot(sub1, reduction = "umap", label = TRUE, cols=color_T)
ggsave("/home/yushiya/data/cd34/data/fig/fig4-ETP-umap-2.pdf", plot=p1, width=6, height=4)
sub1$stim_celltype <- paste(sub1$stim,sub1$seurat_clusters,sep="_")
Idents(sub1) <- sub1$stim_celltype
sub1 <- RenameIdents(sub1, "BM_5"="TSP_BM", "mPB_5"="TSP_mPB", "PB_5"="TSP_PB","thymus_5"="TSP_thymus",
                     "BM_2"="ETP1", "mPB_2"="ETP1", "PB_2"="ETP1","thymus_2"="ETP1",
                     "BM_3"="ETP2", "mPB_3"="ETP2", "PB_3"="ETP2","thymus_3"="ETP2",
                     "BM_1"="ETP3", "mPB_1"="ETP3", "PB_1"="ETP3","thymus_1"="ETP3",
                     "BM_0"="ETP4", "mPB_0"="ETP4", "PB_0"="ETP4","thymus_0"="ETP4",
                     "BM_4"="ETP5", "mPB_4"="ETP5", "PB_4"="ETP5","thymus_4"="ETP5",
                     "BM_6"="ETP6", "mPB_6"="ETP6", "PB_6"="ETP6","thymus_6"="ETP6")
sub1$T_type <- Idents(sub1)
sub1 <- RenameIdents(sub1, "ETP1"="ETP_IGLL1","ETP3"="ETP_TRBC1","ETP2"="ETP_RPS21")
sub1$T_type1 <- Idents(sub1)
sub1$T_type1 <- factor(sub1$T_type1, levels=c("TSP_BM","TSP_mPB","TSP_PB","TSP_thymus","ETP_IGLL1","ETP_TRBC1","ETP_RPS21"))
DefaultAssay(sub1) <- "integrated"
DefaultAssay(sub1) <- "RNA"
sub1 <- ScaleData(sub1, verbose = FALSE)
Idents(sub1) <- sub1$integrated_snn_res.0.5
table(sub1@meta.data[["stim"]],sub1@meta.data[["seurat_clusters"]])
#markers <- FindAllMarkers(sub1, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
markers <- FindAllMarkers(sub1, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(markers,file='/home/yushiya/data/cd34/data/fig/ETP_markers.csv')
sub1 <- RenameIdents(sub1,'1'="c1",'4'="c2",'3'="c3",'0'="c4",'2'="c5",'5'="c6")
FeaturePlot(sub1, features = c("pseudotime"), min.cutoff = 20,max.cutoff = 25,cols = c("blue", "red"),keep.scale=NULL)
FeaturePlot(sub1, features = c("CD44","CD7","CD2","IL7R"), cols = c("#C9E9FF", "red"), order=T,keep.scale=NULL)
ggsave("/home/yushiya/cd34/fig/fig4-ETP-pseu.pdf", plot=p1, width=6, height=4)
VlnPlot(sub1, features = "pseudotime", pt.size = 0, cols=brewer.pal(6,"Paired"))
VlnPlot(sub1, features = "pseudotime", pt.size = 0, cols=color_stim)
VlnPlot(sub1, features = "CD44", pt.size = 0, cols=color_T)
ETP_gene <- c("AVP","CRHBP","NKAIN2","AREG",     #HSPC
              "MPO","	S100A10","IRF8","AZU1",
              "GATA2","PF4",
              "VPREB1","EBF1","MS4A1","LEF1","DNTT",
              "CD7","CD3D","CD3E","TCF7")
DoHeatmap(sub1, features = ETP_gene, 
          group.colors=color_T) + 
  scale_fill_gradientn(colors=c("blue","white","firebrick3"))
ggsave("/home/yushiya/cd34/fig/fig4-ETP-heatmap.pdf", plot=p1, width=6, height=4)
VlnPlot(sub1, features = "My_Score1", pt.size = 0, cols=color_T)
markers.to.plot <- c("HOPX","CRHBP","MEIS1","IGHM","CD74","TSC22D1","JCHAIN","IL7R","CD44",
                     "CD7","CD3D","BCL11B","TRBC2")
plot3=DotPlot(sub1, features = markers.to.plot, dot.scale = 6) +  
  #scale_color_gradientn(colors = c("#F0F921FF","#FDC926FF","#FA9E3BFF","#ED7953FF","#ED7953FF","#BD3786FF","#9C179EFF","#7301A8FF","#47039FFF","#0D0887FF")) +
  scale_color_gradientn(colors = viridis::plasma(10)) +
  theme(panel.border = element_rect(color = "black", fill=NA),
        panel.background = element_blank()) + 
  RotatedAxis()
ggsave("/home/yushiya/data/cd34/data/fig/fig4-ETP-dot-markers.pdf", plot=plot3, width=7, height=3.5)
#saveRDS(sub1, "/home/yushiya/data/cd34/data/fig/ETP_regress.rds")
saveRDS(sub1, "/home/yushiya/data/cd34/data/fig/ETP_mono_1.rds")
SaveH5Seurat(sub1,filename="/home/yushiya/data/cd34/data/fig/immune_ETP.h5seurat", overwrite = TRUE)
Convert("/home/yushiya/data/cd34/data/fig/immune_ETP.h5seurat", dest = "h5ad", overwrite = TRUE)

#sub ETP only select clus 5,2
sub_etp <- subset(sub1, subset= seurat_clusters %in% c("5","2"))
sub_etp$tissue_clus <- paste(sub_etp$tissue, sub_etp$seurat_clusters, sep="_")
Idents(sub_etp) <- sub_etp$tissue_clus
sub_etp <- RenameIdents(sub_etp, "BM_5"="TSP_BM", "mPB_5"="TSP_mPB", "PB_5"="TSP_PB", "thymus_5"="TSP_thymus",
                        "BM_2"="ETP1", "mPB_2"="ETP1", "PB_2"="ETP1", "thymus_2"="ETP1")
sub_etp$new_type <- Idents(sub_etp)
markers <- FindAllMarkers(sub_etp, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
saveRDS(sub_etp, "/home/yushiya/data/cd34/data/fig/ETP_sub.rds")

#rename clusters
immune.combined$celltype_T <- immune.combined$celltype
meta <- immune.combined@meta.data
meta_T <- sub1@meta.data
etp.idx <- which(meta$celltype_T == "ETP-like")
meta$celltype_T <- as.character(meta$celltype_T)
new.type <- meta_T[rownames(meta)[etp.idx], "T_type", drop = TRUE]
new.type[is.na(new.type)] <- "ETP6"
meta$celltype_T[etp.idx] <- new.type
meta$celltype_T <- factor(meta$celltype_T)
immune.combined@meta.data <- meta
#draw umap
DimPlot(immune.combined, reduction = "umap", label = TRUE)
saveRDS(sub1, "/home/yushiya/data/cd34/data/immune.combined-33-celltype-1")
#rename
meta <- immune.combined@meta.data %>% 
  mutate(celltype_T = recode(celltype_T,
                             "ETP2"  = "Thy#1",
                             "ETP4" = "Thy#1",
                             "ETP6" = "Thy#1",
                             "ETP1" = "ETP-like",
                             "ETP3" = "ETP-like",
                             "ETP5" = "ETP-like"))
immune.combined@meta.data <- meta
sub1 <- subset(sub1, subset= T_type %in% c("TSP_BM","TSP_mPB","TSP_PB","TSP_thymus","ETP1","ETP3","ETP5"))
sub1 <- RenameIdents(sub1, "ETP3"="ETP2", "ETP5"="ETP3")
Idents(sub1) <- factor(Idents(sub1), levels=c("TSP_BM","TSP_mPB","TSP_PB","TSP_thymus","ETP1","ETP2","ETP3"))
saveRDS(sub1, "/home/yushiya/data/cd34/data/fig/ETP_mono_2.rds")

#compare TSP_thymus and ETP1
sub_sub <- subset(sub1, subset= T_type %in% c("TSP_thymus","ETP1"))
sub_sub_markers <- FindAllMarkers(sub_sub, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
library(clusterProfiler, lib.loc='/home/panyuwen/miniconda3/envs/clusterprofiler/lib/R/library')
library(org.Hs.eg.db)
library(enrichplot, lib.loc='/home/panyuwen/miniconda3/envs/clusterprofiler/lib/R/library')
library(GOSemSim, lib.loc='/home/panyuwen/miniconda3/envs/clusterprofiler/lib/R/library')
library(DOSE, lib.loc='/home/panyuwen/miniconda3/envs/clusterprofiler/lib/R/library')
TSP_thymusID <- subset(sub_sub_markers, subset= cluster=="TSP_thymus")
TSP_thymus_gene <- bitr(TSP_thymusID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
ETP1ID <- subset(sub_sub_markers, subset= cluster=="ETP1")
ETP1_gene <- bitr(ETP1ID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
cp = list(TSP_thymus_gene=TSP_thymus_gene$ENTREZID,ETP1.gene=ETP1_gene$ENTREZID)
go.p <- compareCluster(cp,
                       fun = "enrichGO",
                       OrgDb = "org.Hs.eg.db",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.01
)
go.p1 <- simplify(go.p,cutoff=0.6,by="p.adjust",select_fun=min)  #去除冗余
p1=ggplot(go.p1, aes(Cluster, Description), showCategory=5) +
  geom_point(aes(color=p.adjust, size=GeneRatio))+
  theme_classic()+
  theme(axis.text.x =element_text(angle = 45,vjust =1,hjust = 1),
        text = element_text(size = 15),
        panel.background = element_blank(),
        panel.grid = element_line(colour = "grey90"))#+theme_bw()
#DEG heatmap
top_genes <- sub_sub_markers %>% group_by(cluster) %>% top_n(n = -11, wt = p_val)
top_genes <- data.frame(top_genes[!grepl("MT-",top_genes$gene),])
top_genes <- top_genes$gene %>% unique()
cellInfo <- data.frame(celltype=sub_sub$T_type)
mtx <- data.frame(sub_sub@assays[["RNA"]]@data) 
colnames(mtx) <- rownames(cellInfo)
top_exp <- sapply(split(rownames(cellInfo), cellInfo$celltype),
                  function(cells) rowMeans(mtx[top_genes,cells]))
top_exp <- top_exp[,c(4,5,6)]
p1=pheatmap(top_exp, cluster_cols = F, #cluster_rows = F,  
            #clustering_method = "median",
            scale = "row",
            fontsize_row = 10, fontsize_col= 12, angle_col = 45,
            col=colorRampPalette(c("#000066","#339999" ,"yellow"))(50))
ggsave("/home/yushiya/data/cd34/data/fig/fig4-DEG-TSPthymus-ETP1.pdf", plot=p1, width=3.5, height=5)
#boxplot of lineage scores
VlnPlot(sub_sub, features = c("My_Score1", "B_Score1", "T_Score1"), ncol = 3, pt.size = 0)
p1=VlnPlot(sub_sub, features = c("CXCR4","CD74","CD44"), ncol = 3, pt.size = 0, cols = color_T1)
ggsave(file=paste("/home/yushiya/data/cd34/data/fig/fig4-TSP_ETP-exp.pdf.pdf",sep=""), plot=p1, width=7, height=3)
mtx <- as.data.frame(sub_sub$T_type)
colnames(mtx)[1] <- "T_type"
mtx$HSC_Score <- sub_sub$HSC_Score1
mtx$B_Score <- sub_sub$B_Score1
mtx$T_Score <- sub_sub$T_Score1
mtx$My_Score <- sub_sub$My_Score1
mtx$pseudotime <- sub_sub$pseudotime
mtx$T_type <- factor(mtx$T_type, levels=c("TSP_thymus","ETP1"))
#plot with different facet
p1=ggboxplot(mtx, x="T_type", y="pseudotime", color = "T_type",
             short.panel.labs = T, ncol=7)+
  scale_color_manual(values = color_T1) +
  #stat_compare_means(label = "p.signif",method = "t.test", comparisons = c("TSP_thymus","ETP1"))+
  stat_compare_means(label = "p.signif",method = "t.test",label.x.npc="center")+
  theme(axis.text.x =element_text(angle = 45,vjust = 1,hjust = 1),
        plot.title = element_text(hjust = 0.5))
ggsave(file=paste("/home/yushiya/data/cd34/data/fig/fig4-TSP_ETP-Score-pseu.pdf.pdf",sep=""), plot=p1, width=3, height=4)



#DEGs of TSP subsets and ETP subsets
sub_tsp <- subset(sub1, subset= T_type %in% c("TSP_BM","TSP_mPB","TSP_PB","TSP_thymus","ETP1"))
TSP_sub_markers <- FindAllMarkers(sub_tsp, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(markers,file='/home/yushiya/data/cd34/data/fig/TSP_sub_markers.csv')
TSP_sub_markers <- read.csv('/home/yushiya/data/cd34/data/fig/TSP_sub_markers.csv')
sub_etp <- subset(sub1, subset= T_type %in% c("ETP1","ETP2","ETP3","ETP4","ETP5","ETP6"))
ETP_sub_markers <- FindAllMarkers(sub_etp, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(ETP_sub_markers,file='/home/yushiya/data/cd34/data/fig/ETP_sub_markers.csv')
ETP_sub_markers <- read.csv('/home/yushiya/data/cd34/data/fig/ETP_sub_markers.csv')
#GO terms of cluster
#top_gene <- markers %>% group_by(cluster) %>% top_n(n = -200, wt = p_val)
#top_gene_clus <- subset(markers, subset= cluster=="TSP_BM")
top_gene_clus <- subset(TSP_sub_markers, subset= cluster=="TSP_thymus")
ego <- enrichGO(gene         = top_gene_clus$gene,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
ego1 <- simplify(ego,cutoff=0.6,by="p.adjust",select_fun=min)
dotplot(ego, showCategory = 6) + 
  labs(title = "ETP1 GO Term Enrichment") 
ggsave(file=paste("/home/yushiya/data/cd34/data/fig/fig4-GO-TSP_BM.pdf",sep=""), plot=p1, width=4, height=3)

#GO barplot
GO_TSP <- readRDS("/home/yushiya/data/cd34/data/fig/GO-TSPs-new.rds")
GO_ETP <- readRDS("/home/yushiya/data/cd34/data/fig/GO-ETPs-new.rds")
mytheme <- theme(axis.text.x = element_text(hjust = 0.5,size = 12), 
                 ## 删去y轴label
                 axis.text.y = element_blank(),
                 ## 删去y轴刻度线
                 axis.ticks.y = element_blank(), 
                 ## 删去x\y轴标题
                 #axis.title.x = element_blank(), 
                 #axis.title.y = element_blank(), 
                 axis.title = element_text(size = 20),
                 axis.title.x = element_text(size = 16),  # 修改横轴标题字体大小
                 axis.title.y = element_text(size = 16),  # 修改纵轴标题字体大小
                 plot.title = element_text(hjust = 0.5,size =  18),
                 legend.position = "none")
#sub_GO <- GO_TSP@compareClusterResult
sub_GO <- GO_ETP@compareClusterResult
#sub_GO <- subset(sub_GO, subset= Cluster=="TSP_BM.gene")
sub_GO <- ego@result
#sub_GO <- subset(sub_GO, subset= Cluster=="ETP1.gene")
sub_GO$log_pvalue <- -log(as.numeric(sub_GO$pvalue))
sub_GO_top <- sub_GO[c(1:5),]   #TSP_BM
sub_GO_top <- sub_GO[c(1,3,5,6,78),]   #ETP1
sub_GO_top <- sub_GO[c(1,4,25,26,32),]   #TSP_thymus
sub_GO_top$Description <- factor(sub_GO_top$Description,levels = rev(sub_GO_top$Description))
p1=ggplot(data = sub_GO_top, aes(x = Description, y = log_pvalue)) +
  geom_bar(stat = "identity", width = 0.8, fill = "#4C82C5",alpha = 0.8) + #绘制条形图
  geom_text(aes(y = 0, #控制文本标签起始位置
                label = Description),
            size = 6,hjust = 0) + #hjust = 0左对齐
  coord_flip() + theme_classic() + mytheme +
  scale_y_continuous(limits = c(0, 30)) + # 去掉留白
  labs(x = "GO Terms", y = "-log p value") 
ggsave("/home/yushiya/data/cd34/data/fig/fig4-GO-TSP_thymus.pdf", plot=p1, width=6, height=4.5)

#chemotaxis genes expression
chemo_genes <- c("PDE4B","PPIB","CALR","PPIA","LEF1","HMGB1","CXCR4","ADGRE2","CD74","ELANE","HMGB2","RNASE2","LSP1",
                 "CD44","LGALS1","PYCARD","IL1B","LYN","TSPAN32","ITGA4","ZMIZ1")
sub_etp <- subset(sub1, subset= T_type %in% c("TSP_BM","TSP_mPB","TSP_PB","TSP_thymus","ETP1"))
cellInfo <- data.frame(Celltype=sub_etp$T_type)
mtx <- data.frame(sub_etp@assays[["RNA"]]@data[chemo_genes,]) 
colnames(mtx) <- rownames(cellInfo)
chemo_exp <- sapply(split(rownames(cellInfo), cellInfo$Celltype),
                    function(cells) rowMeans(mtx[chemo_genes,cells]))
lev<-c("TSP_BM","TSP_mPB","TSP_PB","TSP_thymus","ETP1")
chemo_exp <- chemo_exp[,lev] 
#annotation_col <- data.frame(Stim = factor(rep(c("BM", "mPB","PB"), c(3,3,3))))
#annotation_colors =list(Stim=c(BM="#DB5C25",mPB="#F3B747",PB="#649541"))
p1=pheatmap(chemo_exp, cluster_cols = F, #cluster_rows = F,  
         clustering_method = "average",
         #annotation_col = annotation_col,
         #annotation_colors = annotation_colors,
         scale = "row",
         #border=F,
         col=colorRampPalette(c("navy","white" ,"firebrick3"))(50),
         fontsize_row = 12, fontsize_col= 12, angle_col = 45)
ggsave("/home/yushiya/data/cd34/data/fig/fig4-DEG-migr-TSPs.pdf", plot=p1, width=5, height=4.5)


#gene expression through pseudotime
mtx <- as.data.frame(sub_etp$new_type)
colnames(mtx)[1] <- "celltype"
mtx$pseudotime <- sub_etp$pseudotime
mtx$CALR <- sub_etp@assays[["RNA"]]@data["CALR",]
mtx$CXCR4 <- sub_etp@assays[["RNA"]]@data["CXCR4",]
mtx$CD74 <- sub_etp@assays[["RNA"]]@data["CD74",]
mtx$CD3E <- sub_etp@assays[["RNA"]]@data["CD3E",]
mtx$CD2 <- sub_etp@assays[["RNA"]]@data["CD2",]
mtx$CD7 <- sub_etp@assays[["RNA"]]@data["CD7",]
mtx$CD44 <- sub_etp@assays[["RNA"]]@data["CD44",]
#col <- c(color_BMPB_lineages,color_stim)
p1=ggplot(mtx, aes(x=pseudotime,y=CD44))+
  geom_point(size=0.3,aes(colour=celltype),alpha = 0)+
  stat_smooth(mapping = aes(colour=celltype), se = FALSE, linewidth=1.5)+
  #facet_wrap(~celltype)+#theme_bw()
  #scale_y_continuous(limits = c(-0.5,0.5))+
  scale_color_manual(values = color_T) +
  theme_classic()+
  theme(axis.text.x=element_text(hjust=0.5),
        text = element_text(size = 18),
        #panel.border = element_rect(color = "black", fill=NA),
        panel.background = element_blank())
#panel.grid = element_line(colour = "grey60"),
#axis.line = element_line(colour = "grey60",size = rel(1),arrow = arrow(angle = 30,length = unit(0.1,"inches"))))
ggsave(file=paste("/home/yushiya/data/cd34/data/fig/fig4-ETP-pseu-CD44.pdf",sep=""), plot=p1, width=5, height=4)


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
##过滤低质量细胞，过滤低于1%细胞中检出的基因 
cds <- detectGenes(cds, min_expr = 0.1)
#differentialGeneTest 找高变基因
expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 10))
diff <- differentialGeneTest(cds[expressed_genes,],
                             fullModelFormulaStr="~T_type",
                             cores=10) 
dergene <- subset(diff, qval < 0.01)
dergene <- dergene[order(dergene$pval,decreasing=F),]
ordergene <-row.names(dergene[1:3500,])
cds <- setOrderingFilter(cds, ordergene)
#cluster 高变基因
deg.cluster <- read.csv("/home/yushiya/data/cd34/data/fig/ETP_markers.csv",stringsAsFactors = F)
express_genes <- subset(deg.cluster,p_val_adj<0.01)
express_genes <- express_genes %>% group_by(cluster) %>% top_n(n = -120, wt = p_val)
express_genes <- express_genes$gene
cds <- setOrderingFilter(cds,express_genes)
#“dpFeature”选择高变基因
sub1 <- FindVariableFeatures(object = sub1)
expressed_genes<- VariableFeatures(sub1)
cds <- detectGenes(cds, min_expr = 0.1) 
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 10)) 
diff <-differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr="~T_type",cores=30)
deg <- subset(diff, qval < 0.01)
deg <-deg[order(deg$qval,decreasing=F),]
ordergene <-row.names(deg)[order(deg$qval)][1:1500]
#ordergene <- rownames(deg)
cds <- setOrderingFilter(cds,ordergene) 
plot_ordering_genes(cds)
#order cells
cds <- reduceDimension(cds,max_components = 2,method = 'DDRTree')
source("/home/yushiya/code/order_cells.R")
cds <- orderCells(cds)
cds <- orderCells(cds, root_state = 5)
p1=plot_cell_trajectory(cds,color_by="Pseudotime",size=1,show_backbone=TRUE)
ggsave("/home/yushiya/data/cd34/data/fig/fig4-ETP-pseu_mo2_pseu-1.pdf", plot=p1, width=5, height=4)
p1=plot_cell_trajectory(cds,color_by="T_type1",size=1,show_backbone=TRUE)+ scale_color_manual(values = color_T)
plot_cell_trajectory(cds,color_by="integrated_snn_res.0.6",size=1,show_backbone=TRUE)
plot_cell_trajectory(cds,color_by="T_type",size=1,show_backbone=TRUE)
ggsave("/home/yushiya/data/cd34/data/fig/fig4-ETP-pseu_mo2-1.pdf", plot=p1, width=5, height=4)
plot_cell_trajectory(cds,color_by="State",size=1,show_backbone=TRUE)
plot_cell_trajectory(cds, color_by = "integrated_snn_res.0.6") + facet_wrap("~integrated_snn_res.0.6", nrow = 2)
p1=plot_cell_trajectory(cds, color_by = "T_type1") + facet_wrap("~T_type1", nrow = 2)+ scale_color_manual(values = color_T)
ggsave("/home/yushiya/data/cd34/data/fig/fig4-ETP-pseu_mo2-split-1.pdf", plot=p1, width=5, height=4)
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
ggsave("/home/yushiya/data/cd34/data/fig/fig4-ETP-pseu_mo2_markers-1.pdf", plot=p, width=6, height=6)
s.genes <-c("IL7R","CD84")
plot_genes_in_pseudotime(cds[s.genes,], color_by = "stim")
sub1$pseu_mo2 <- cds$Pseudotime
save(cds,file="/home/yushiya/data/cd34/data/fig/ETP.cds.12-1.RData")
#DEG
load("/home/yushiya/data/cd34/data/fig/ETP.cds.3.RData")
etp_monocle <- cds
peu_gene <- differentialGeneTest(etp_monocle,fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = 20)
write.csv(peu_gene,file='/home/yushiya/data/cd34/data/fig/ETP_peu_gene.csv')
peu_gene <- read.csv(file='/home/yushiya/data/cd34/data/fig/ETP_peu_gene.csv', row.names = 1)
peu_gene <- peu_gene[which(peu_gene$qval<0.01 & peu_gene$num_cells_expressed>100),]
peu_gene %>% arrange(qval)  -> peu_gene#按照qval排个序
peu_gene <- peu_gene[1:1000,] #这里我们取前100个基因演示
source('/home/yushiya/code/monocle2_heatmap.R')
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
write.csv(Module_GO, file = '/home/yushiya/data/cd34/data/fig/fig5-ETP_Module_GO.csv')
top_Module_GO <- Module_GO %>% group_by(cluster) %>% top_n(n = -10, wt = p.adjust)
gene <- c("IL7R", "CAMK4","DOCK2","CDK6","CD3D","CD3E","CD3G","BCL11B","LCK","CD2","SOX4","TCF7",
          "CD74","HLA-A","HLA-DRA","HLA-DRB1","B2M",
          "RPL11","RPS21","RPS14","RPL7")
gene <- c("FOXP1","PDE4D","HLA-DPB1","IGHM","HLA-DRB1","ELF1","VAV3","BCL2","PDE4B","CD81","CD79B","PTPRC","AREG",
          "RPS8","RPL7A","RPS13","RPS28")
source('/home/yushiya/code/add.flag.R')
p <- add.flag(p11,kept.labels = gene,repel.degree = 0.2)
ggsave("/home/yushiya/data/cd34/data/fig/fig4-ETP-pseu_heatmap.pdf", plot=p, width=5, height=6)

#ETP subsets through pseudotime
#expr.df <- data.frame(sub1$T_type1)
expr.df <- data.frame(sub_etp$new_type)
colnames(expr.df)[1] <- "T_type"
expr.df$pseudotime <- sub_etp$pseudotime
p1=ggplot(expr.df, aes(y= T_type, x=pseudotime, color=T_type))+
  geom_jitter(size=1, alpha=1)+
  scale_color_manual(values = color_T) +
  theme_bw()+
  labs(title='Cell types', y='', x='Pseudotime')&NoLegend()
ggsave("/home/yushiya/data/cd34/data/fig/fig4-ETP-pseu_celltype.pdf", plot=p1, width=4, height=4)
expr.df$exp <- sub1@assays[["RNA"]]@data["CD44",]
p1=ggplot(expr.df, aes(x=pseudotime,y=exp))+
  geom_point(size=0.3,aes(colour=T_type),alpha = 0)+
  stat_smooth(mapping = aes(colour=T_type))+
  #facet_wrap(~celltype)+#theme_bw()
  #scale_y_continuous(limits = c(-0.5,0.5))+
  #scale_color_manual(values = color_stim) +
  scale_color_manual(values = color_T) +
  theme_classic()+
  theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size = 12),
        panel.border = element_blank(), axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))+
  #scale_color_manual(values = colors_protein)+
  coord_cartesian(ylim = c(0,2))+
  labs(title='', y='smoothed mean', x='Pseudotime')
ggsave("/home/yushiya/data/cd34/data/fig/fig4-ETP-pseu_CD44.pdf", plot=p1, width=5, height=3)


#monocle 3
#dyn.load("/home/wangjl/local/gdal-2.4.4/libgdal.so.20.5.4")
#library(monocle3)
library(Matrix)
library(ggplot2)
library(monocle3,lib.loc='/home/panyuwen/miniconda3/envs/monocle3/lib/R/library')
library(leidenbase, lib.loc='/home/panyuwen/miniconda3/envs/monocle3/lib/R/library') 
library(spdep, lib.loc='/home/panyuwen/miniconda3/envs/monocle3/lib/R/library') 
library(pbmcapply, lib.loc='/home/panyuwen/miniconda3/envs/monocle3/lib/R/library') 
sub1 <- readRDS("/home/yushiya/data/cd34/data/fig/ETP_mono.rds")
data <- GetAssayData(sub1, assay = 'RNA', slot = 'counts')
cell_metadata <- sub1@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
save(cds,file="/home/yushiya/data/cd34/data/fig/ETP_mono3_clus.cds.RData")
load(file="/home/yushiya/data/cd34/data/fig/ETP_mono3_clus.cds.RData")
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, preprocess_method = "PCA")
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('cds.umap')
##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(sub1, reduction = "umap")
cds@int_colData$reducedDims$UMAP <- int.embed
## Monocle3聚类分区
cds <- cluster_cells(cds)
#cds@clusters@listData[["UMAP"]][["clusters"]] <- sub1$integrated_snn_res.0.5
## 识别轨迹
cds <- learn_graph(cds)
p<-plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
              label_branch_points = FALSE)
ggsave("/home/yushiya/data/cd34/data/fig/fig4-ETP-pseu_branch_clus.pdf", plot=p, width=6, height=4)
p1 <- p + geom_vline(xintercept = seq(-8,-7,0.25)) + geom_hline(yintercept = seq(1.5,2.5,0.25))
embed <- data.frame(Embeddings(sub1, reduction = "umap"))
embed <- subset(embed, UMAP_1 > -7.5 & UMAP_1 < -7 & UMAP_2 > 1.5 & UMAP_2 < 2)
root.cell <- rownames(embed)
cds <- order_cells(cds, root_cells = root.cell)
#cds <- order_cells(cds)
p1=plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
              label_leaves = FALSE,  label_branch_points = FALSE)
ggsave("/home/yushiya/data/cd34/data/fig/fig4-ETP-pseu_branch_clus-1.pdf", plot=p1, width=6.5, height=4)
p1=plot_cells(cds, genes=c("CD44","CD7","CD2","IL7R"),
              show_trajectory_graph=T,
              label_cell_groups=FALSE,
              label_leaves=FALSE, label_branch_points=F)
ggsave("/home/yushiya/data/cd34/data/fig/fig4-ETP-pseu_branch_clus_markers.pdf", plot=p1, width=6, height=4)
sub1$pseudotime <- pseudotime(cds)
#find gene modules
cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=80)  #DEGs along trajectory
write.csv(cds_pr_test_res, "/home/yushiya/data/cd34/data/fig/ETP_peu_tra_gene-1.csv")
cds_pr_test_res <- read.csv("/home/yushiya/data/cd34/data/fig/ETP_peu_tra_gene-1.csv", row.names = 1)
pr_deg_ids <- row.names(subset(cds_pr_test_res, q_value < 0.01 & morans_I > 0.2))
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=c(10^seq(-8,-1)))
table(gene_module_df$module)
write.csv(gene_module_df, "/home/yushiya/data/cd34/data/fig/ETP-peu_tra_module-1.csv")
gene_module_df <- read.csv("/home/yushiya/data/cd34/data/fig/ETP-peu_tra_module-1.csv", row.names = 1)
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$integrated_snn_res.0.5)
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=Idents(sub1))
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
#write.table(t(agg_mat), "/home/yushiya/data/cd34/data/fig/ETP-peu_tra_module-aggmat.csv", sep = "\t")
#agg_mat <- read.table("/home/yushiya/data/cd34/data/fig/ETP-peu_tra_module-aggmat.csv", sep="\t")
p<-pheatmap::pheatmap(agg_mat,cluster_rows = T, cluster_cols = T,
                      scale="column", clustering_method="ward.D2")
ggsave("/home/yushiya/data/cd34/data/fig/ETP-tra-module-heatmap-1.pdf", plot=p, width=6, height=4)
#GO terms in all modules
module_gene <- gene_module_df[,c(1,2)]
colnames(module_gene) <- c("gene","Module")
rownames(module_gene) <- module_gene$gene
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
write.csv(Module_GO, file = '/home/yushiya/data/cd34/data/fig/ETP-Module_GO-1.csv')
#draw GO modules
mytheme <- theme(axis.text.x = element_text(hjust = 0.5,size = 12), 
                 ## 删去y轴label
                 axis.text.y = element_blank(),
                 ## 删去y轴刻度线
                 axis.ticks.y = element_blank(), 
                 ## 删去x\y轴标题
                 #axis.title.x = element_blank(), 
                 #axis.title.y = element_blank(), 
                 axis.title = element_text(size = 20),
                 axis.title.x = element_text(size = 16),  # 修改横轴标题字体大小
                 axis.title.y = element_text(size = 16),  # 修改纵轴标题字体大小
                 plot.title = element_text(hjust = 0.5,size =  18),
                 legend.position = "none")
sub_GO <- subset(Module_GO, subset= cluster==4)
sub_GO$log_pvalue <- -log(as.numeric(sub_GO$pvalue))
sub_GO_top <- sub_GO[c(1:5),]   #BM module 
sub_GO_top$Description <- factor(sub_GO_top$Description,levels = rev(sub_GO_top$Description))
p1=ggplot(data = sub_GO_top, aes(x = Description, y = log_pvalue)) +
  geom_bar(stat = "identity", width = 0.8, fill = "firebrick3",alpha = 0.8) + #绘制条形图
  geom_text(aes(y = 0, #控制文本标签起始位置
                label = Description),
            size = 6,hjust = 0) + #hjust = 0左对齐
  coord_flip() + theme_bw() + mytheme +
  labs(x = "GO Terms", y = "-log p value", title = "Module 4 GO Term Enrichment") 
ggsave("/home/yushiya/data/cd34/data/fig/fig4-ETP-tra-module-GO-4.pdf", plot=p1, width=5, height=4.5)



#Expression via pseudotime
pseu_exp<-as.data.frame(colnames(sub1))
pseu_exp$stim <- sub1$stim
pseu_exp$celltype <- sub1$T_type1
pseu_exp$pseudotime <- sub1$pseu_mo2
#pseu_exp$celltype <- sub1$celltype
#pseu_exp$pseudotime <- sub1$pseudotime
names(pseu_exp)[1]='ID'
pseu_exp$Migratory_Score <- sub1$Migratory_Score1
pseu_exp$HSC_Score <- sub1$HSC_Score1
pseu_exp$B_Score <- sub1$B_Score1
pseu_exp$T_Score <- sub1$T_Score1
pseu_exp$My_Score <- sub1$My_Score1
pseu_exp$Expression <- sub1@assays[["RNA"]]@data["NCL",]
p4=ggplot(pseu_exp, aes(x=pseudotime,y=T_Score))+
  geom_point(size=0.7,aes(colour=celltype))+
  #stat_smooth(mapping = aes(colour=celltype),method = "loess")+
  #facet_wrap(~celltype)+#theme_bw()
  #scale_x_continuous(limits = c(20,25))+
  scale_color_manual(values = color_T) +
  theme_classic()+
  theme(axis.text.x=element_text(hjust=0.5),
        text = element_text(size = 18),
        panel.background = element_blank(),
        panel.grid = element_line(colour = "grey90"),
        axis.line = element_line(colour = "grey60",size = rel(1),arrow = arrow(angle = 30,length = unit(0.1,"inches"))))
ggplot(pseu_exp, aes(x=pseudotime,y=Expression))+
  geom_point(size=0.3,aes(colour=stim))+
  stat_smooth(mapping = aes(colour=stim))+
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
ggsave("/home/yushiya/data/cd34/data/fig/fig4-ETP-scores-1.pdf", plot=p, width=10, height=5)
#Expression in vlnplot splited by stim
VlnPlot(sub1, features = "T_Score1", pt.size = 0, split.by = "stim",cols=color_stim)
#DEGs of TFs and surface markers of ETP
library(pheatmap)
ETP_markers <- read.csv("/home/yushiya/data/cd34/data/fig/cpstim_ETP.csv",stringsAsFactors = F)
cell_surface <- read.csv("/home/yushiya/data/cd34/cell_surface_protein.csv")
cs.ETP_markers <- ETP_markers[ETP_markers$gene %in% cell_surface$ENTREZ.gene.symbol,]
TF <- read.table("/home/yushiya/data/cd34/TF.txt", sep="\t", header=T)
tf.ETP_markers <- ETP_markers[ETP_markers$gene %in% TF$Symbol,]
selected_ETP.markers <- ETP_markers[(ETP_markers$gene %in% TF$Symbol | ETP_markers$gene %in% cell_surface$ENTREZ.gene.symbol),]
#TSP subsets
top_genes <- TSP_sub_markers %>% group_by(cluster) %>% top_n(n = -5, wt = p_val)
top_genes <- TSP_sub_markers %>% group_by(cluster) %>% top_n(n = 8, wt = avg_log2FC)
top_genes <- data.frame(top_genes[!grepl("MT-",top_genes$gene),])
top_genes <- top_genes$gene %>% unique()
cellInfo <- data.frame(celltype=sub_etp$T_type)
mtx <- data.frame(sub_etp@assays[["RNA"]]@data) 
colnames(mtx) <- rownames(cellInfo)
top_exp <- sapply(split(rownames(cellInfo), cellInfo$celltype),
                  function(cells) rowMeans(mtx[top_genes,cells]))
top_exp <- top_exp[,c(1:4)]
p1=pheatmap(top_exp, cluster_cols = F, #cluster_rows = F,  
            #clustering_method = "median",
            scale = "row",
            fontsize_row = 10, fontsize_col= 12, angle_col = 45,
            col=colorRampPalette(c("#000066","#339999" ,"yellow"))(50))
ggsave("/home/yushiya/data/cd34/data/fig/fig4-TSP-DEG-2.pdf", plot=p1, width=3.5, height=5)
#ETP subsets
sub1_etp <- subset(sub1, subset= T_type %in% c("ETP1","ETP2","ETP3","ETP4","ETP5","ETP6"))
ETP_sub_markers <- FindAllMarkers(sub1_etp, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
top_genes <- ETP_sub_markers %>% group_by(cluster) %>% top_n(n = -5, wt = p_val)
top_genes <- ETP_sub_markers %>% group_by(cluster) %>% top_n(n = 8, wt = avg_log2FC)
top_genes <- data.frame(top_genes[!grepl("MT-",top_genes$gene),])
top_genes <- top_genes$gene %>% unique()
cellInfo <- data.frame(celltype=sub_etp$T_type)
mtx <- data.frame(sub_etp@assays[["RNA"]]@data) 
colnames(mtx) <- rownames(cellInfo)
top_exp <- sapply(split(rownames(cellInfo), cellInfo$celltype),
                  function(cells) rowMeans(mtx[top_genes,cells]))
top_exp <- top_exp[,c(5:10)]
p1=pheatmap(top_exp, cluster_cols = F, #cluster_rows = F,  
            #clustering_method = "median",
            scale = "row",
            fontsize_row = 10, fontsize_col= 12, angle_col = 45,
            col=colorRampPalette(c("#000066","#339999" ,"yellow"))(50))
ggsave("/home/yushiya/data/cd34/data/fig/fig4-ETP-DEG-2.pdf", plot=p1, width=3.5, height=5)
#TSP and ETP markers
Idents(sub1) <- sub1$stim
sub1 <- RenameIdents(sub1, "TSP_BM"="TSP", "TSP_PB"="TSP","TSP_mPB"="TSP", "TSP_thymus"="TSP","ETP1"="ETP","ETP2"="ETP","ETP3"="ETP")
sub1$newtype <- Idents(sub1)
TSP_ETP.markers <- FindAllMarkers(sub1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.405)
TSPID <- subset(TSP_ETP.markers, subset= cluster=="TSP")
TSP_gene <- bitr(TSPID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
ETPID <- subset(TSP_ETP.markers, subset= cluster=="ETP")
ETP_gene <- bitr(ETPID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
cp = list(TSP.gene=TSP_gene$ENTREZID, ETP.gene=ETP_gene$ENTREZID)
#GO analysis
library(clusterProfiler, lib.loc='/home/panyuwen/miniconda3/envs/clusterprofiler/lib/R/library')
library(org.Hs.eg.db)
library(enrichplot, lib.loc='/home/panyuwen/miniconda3/envs/clusterprofiler/lib/R/library')
library(GOSemSim, lib.loc='/home/panyuwen/miniconda3/envs/clusterprofiler/lib/R/library')
library(DOSE, lib.loc='/home/panyuwen/miniconda3/envs/clusterprofiler/lib/R/library')
markers <- read.csv("/home/yushiya/data/cd34/data/fig/ETP_markers.csv",stringsAsFactors = F)
#TSP subsets
TSP_markers <- subset(markers, subset= cluster %in% c("TSP_BM","TSP_mPB","TSP_PB","TSP_thymus"))
BMID <- subset(top_gene, subset= cluster=="TSP_BM")
BM_gene <- bitr(BMID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
PBID <- subset(markers, subset= cluster=="TSP_PB")
PB_gene <- bitr(PBID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
mPBID <- subset(markers, subset= cluster=="TSP_mPB")
mPB_gene <- bitr(mPBID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
thymusID <- subset(markers, subset= cluster=="TSP_thymus")
thymus_gene <- bitr(thymusID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
cp = list(TSP_BM.gene=BM_gene$ENTREZID,TSP_mPB.gene=mPB_gene$ENTREZID,TSP_PB.gene=PB_gene$ENTREZID,TSP_thymus.gene=thymus_gene$ENTREZID)
#ETP subsets
ETP_markers <- subset(markers, subset= cluster %in% c("ETP1","ETP2","ETP3"))
ETP1ID <- subset(ETP_sub_markers, subset= cluster=="ETP1")
ETP1_gene <- bitr(ETP1ID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
ETP2ID <- subset(ETP_sub_markers, subset= cluster=="ETP2")
ETP2_gene <- bitr(ETP2ID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
ETP3ID <- subset(ETP_sub_markers, subset= cluster=="ETP3")
ETP3_gene <- bitr(ETP3ID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
ETP4ID <- subset(ETP_sub_markers, subset= cluster=="ETP4")
ETP4_gene <- bitr(ETP4ID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
ETP5ID <- subset(ETP_sub_markers, subset= cluster=="ETP5")
ETP5_gene <- bitr(ETP5ID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
ETP6ID <- subset(ETP_sub_markers, subset= cluster=="ETP6")
ETP6_gene <- bitr(ETP6ID$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
cp = list(ETP1.gene=ETP1_gene$ENTREZID,ETP2.gene=ETP2_gene$ENTREZID,ETP3.gene=ETP3_gene$ENTREZID,ETP4.gene=ETP4_gene$ENTREZID,ETP5.gene=ETP5_gene$ENTREZID,ETP6.gene=ETP6_gene$ENTREZID)
cp = list(ETP1.gene=ETP1_gene$ENTREZID,ETP2.gene=ETP3_gene$ENTREZID,ETP3.gene=ETP5_gene$ENTREZID)
go.p <- compareCluster(cp,
                       fun = "enrichGO",
                       OrgDb = "org.Hs.eg.db",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.01
)
go.p1 <- simplify(go.p,cutoff=0.6,by="p.adjust",select_fun=min)  #去除冗余
p1=ggplot(go.p1, aes(Cluster, Description), showCategory=5) +
  geom_point(aes(color=p.adjust, size=GeneRatio))+
  theme_classic()+
  theme(axis.text.x =element_text(angle = 45,vjust =1,hjust = 1),
        text = element_text(size = 15),
        panel.background = element_blank(),
        panel.grid = element_line(colour = "grey90"))#+theme_bw()
saveRDS(go.p, "/home/yushiya/data/cd34/data/fig/GO-TSPs-new.rds")
saveRDS(go.p, "/home/yushiya/data/cd34/data/fig/GO-ETPs-new-1.rds")
ggsave("/home/yushiya/data/cd34/data/fig/fig4-GO-TSPs-new.pdf", plot=p1, width=9, height=4.5)
ggsave("/home/yushiya/data/cd34/data/fig/fig4-GO-ETPs-new-1.pdf", plot=p1, width=8.5, height=4.5)
A1 <- go.p@compareClusterResult[c(1,8,9,10,11,122,123,128,129),]
A <- data.frame(A1[["Description"]])
colnames(A)[1] <- "Description"
A$p.adjust <- as.numeric(A1[["p.adjust"]])
A$p.adjust <- -log(A$p.adjust)
A$group <- A1[["Cluster"]]
A$p.adjust[6:9] <- -A$p.adjust[6:9]   #TSP ETP
A_draw <- A[c(1:5,6:9),]  #TSP ETP
A$GeneRatio_num[16:31] <- -A$GeneRatio_num[16:31]   #PB thymus
A_draw <- A[c(1:5,17:21),]  #PB thymus
A$GeneRatio_num[20:35] <- -A$GeneRatio_num[20:35]   #BM thymus
A_draw <- A[c(1:5,20:25),]  #BM thymus
A$GeneRatio_num[20:35] <- -A$GeneRatio_num[20:35]   #mPB thymus
A_draw <- A[c(1:4,20:25),]  #mPB thymus
p1=ggplot(A_draw,aes(reorder(Description, p.adjust),p.adjust,fill=group))+
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
  geom_segment(aes(y=0, yend=0,x=0,xend=9))+
  geom_text(data = A_draw[which(A_draw$p.adjust>0),],aes(x=Description, y=-0.01, label=Description),
            hjust=1, size=4)+
  geom_text(data = A_draw[which(A_draw$p.adjust<0),],aes(x=Description, y=0.01, label=Description),
            hjust=0, size=4)+
  geom_text(data = A_draw[which(A_draw$p.adjust>0),],aes(label=p.adjust),
            hjust=-0.1, size=3, color='red')+
  geom_text(data = A_draw[which(A_draw$p.adjust<0),],aes(label=p.adjust),
            hjust=1.1, size=3, color="red")+
  scale_fill_manual(values = c("#98ADC4",
                               "#DDA0DD"))+
  scale_x_discrete(expand = expansion(mult = c(0,0)))+
  ylim(-25, 25)+
  labs(x='', y='-log(p_value)')
#color TSP, ETP, "#DDA0DD","#98ADC4","#7F689D","#54689A"
ggsave("/home/yushiya/data/cd34/data/fig/fig4-GO-TSP-ETP.pdf", plot=p1, width=8.5, height=5.5)


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


#add thymus data
t1.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2020_thymus_Lavaert/Thy1/Thy1/outs/filtered_feature_bc_matrix")
t1 <- CreateSeuratObject(counts = t1.data, project = "thymocyte1")
t1$ori_stim <- 'Thy1'
t1$stim <- 'thymus'
t1[["percent.mt"]] <- PercentageFeatureSet(t1, pattern = "^MT-")
VlnPlot(t1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
t1 <- subset(t1, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(t1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
t1 <- NormalizeData(t1, normalization.method = "LogNormalize", scale.factor = 10000)
t1 <- FindVariableFeatures(t1, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(t1)
t1 <- ScaleData(t1, features = all.genes)

t2.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2020_thymus_Lavaert/Thy2/Thy2/outs/filtered_feature_bc_matrix")
t2 <- CreateSeuratObject(counts = t2.data, project = "thymocyte1")
t2$ori_stim <- 'Thy2'
t2$stim <- 'thymus'
t2[["percent.mt"]] <- PercentageFeatureSet(t2, pattern = "^MT-")
VlnPlot(t2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
t2 <- subset(t2, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(t2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
t2 <- NormalizeData(t2, normalization.method = "LogNormalize", scale.factor = 10000)
t2 <- FindVariableFeatures(t2, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(t2)
t2 <- ScaleData(t2, features = all.genes)

t3.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2020_thymus_Lavaert/Thy3/Thy3/outs/filtered_feature_bc_matrix")
t3 <- CreateSeuratObject(counts = t3.data, project = "thymocyte1")
t3$ori_stim <- 'Thy3'
t3$stim <- 'thymus'
t3[["percent.mt"]] <- PercentageFeatureSet(t3, pattern = "^MT-")
VlnPlot(t3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
t3 <- subset(t3, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(t3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
t3 <- NormalizeData(t3, normalization.method = "LogNormalize", scale.factor = 10000)
t3 <- FindVariableFeatures(t3, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(t3)
t3 <- ScaleData(t3, features = all.genes)

t4.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2020_thymus_Lavaert/Thy4/Thy4/outs/filtered_feature_bc_matrix")
t4 <- CreateSeuratObject(counts = t4.data, project = "thymocyte1")
t4$ori_stim <- 'Thy4'
t4$stim <- 'thymus'
t4[["percent.mt"]] <- PercentageFeatureSet(t4, pattern = "^MT-")
VlnPlot(t4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
t4 <- subset(t4, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(t4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
t4 <- NormalizeData(t4, normalization.method = "LogNormalize", scale.factor = 10000)
t4 <- FindVariableFeatures(t4, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(t4)
t4 <- ScaleData(t4, features = all.genes)

t5.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2020_thymus_Lavaert/Thy5/Thy5/outs/filtered_feature_bc_matrix")
t5 <- CreateSeuratObject(counts = t5.data, project = "thymocyte1")
t5$ori_stim <- 'Thy5'
t5$stim <- 'thymus'
t5[["percent.mt"]] <- PercentageFeatureSet(t5, pattern = "^MT-")
VlnPlot(t5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
t5 <- subset(t5, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(t5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
t5 <- NormalizeData(t5, normalization.method = "LogNormalize", scale.factor = 10000)
t5 <- FindVariableFeatures(t5, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(t5)
t5 <- ScaleData(t5, features = all.genes)

t6.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2020_thymus_Lavaert/Thy6/Thy6/outs/filtered_feature_bc_matrix")
t6 <- CreateSeuratObject(counts = t6.data, project = "thymocyte1")
t6$ori_stim <- 'Thy6'
t6$stim <- 'thymus'
t6[["percent.mt"]] <- PercentageFeatureSet(t6, pattern = "^MT-")
VlnPlot(t6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
t6 <- subset(t6, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(t6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
t6 <- NormalizeData(t6, normalization.method = "LogNormalize", scale.factor = 10000)
t6 <- FindVariableFeatures(t6, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(t6)
t6 <- ScaleData(t6, features = all.genes)

t7.data <- Read10X(data.dir="/home/yushiya/data/cd34/data/2020_thymus_Lavaert/Thy7/Thy7/outs/filtered_feature_bc_matrix")
t7 <- CreateSeuratObject(counts = t7.data, project = "thymocyte1")
t7$ori_stim <- 'Thy7'
t7$stim <- 'thymus'
t7[["percent.mt"]] <- PercentageFeatureSet(t7, pattern = "^MT-")
VlnPlot(t7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
t7 <- subset(t7, subset = nCount_RNA > 1000 & nFeature_RNA > 200 & percent.mt < 10)
VlnPlot(t7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
t7 <- NormalizeData(t7, normalization.method = "LogNormalize", scale.factor = 10000)
t7 <- FindVariableFeatures(t7, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(t7)
t7 <- ScaleData(t7, features = all.genes)

thy.anchors <- FindIntegrationAnchors(object.list = c(t3,t4,t5,t6,t7), dims = 1:20)
thy <- IntegrateData(anchorset = thy.anchors, dims = 1:20)
DefaultAssay(thy) <- "integrated"
# Run the standard workflow for visualization and clustering
thy <- ScaleData(thy, verbose = FALSE)
thy <- RunPCA(thy, npcs = 30, verbose = FALSE)
ElbowPlot(thy)
thy <- FindVariableFeatures(thy, selection.method = "vst", nfeatures = 2000)
#Clustering
thy <- FindNeighbors(thy, reduction = "pca", dims = 1:30)
thy <- FindClusters(thy, resolution = 0.8)
# umap
thy <- RunUMAP(thy, reduction = "pca", dims = 1:30)
DimPlot(thy, reduction = "umap", split.by = "ori_stim")
DimPlot(thy, reduction = "umap", label = TRUE)
saveRDS(thy, file = "/home/yushiya/data/cd34/data/2020_thymus_Lavaert/thy_combined.rds")
thy <- readRDS("/home/yushiya/data/cd34/data/2020_thymus_Lavaert/thy_combined.rds")

immune.combined <- RunUMAP(immune.combined, dims = 1:40, reduction = "pca", return.model = TRUE)
immune.combined <- FindVariableFeatures(immune.combined, selection.method = "vst", nfeatures = 2000)
thy.anchors <- FindTransferAnchors(reference = immune.combined, query = thy,
                                   dims = 1:30, reference.reduction = "pca")
thy <- MapQuery(anchorset = thy.anchors, reference = immune.combined, query = thy,
                refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap")
#thy$predicted.celltype <- factor(thy$predicted.celltype, levels=c("HSC","MPP #1","MPP #2","MEP","EryP","MkP","LMPP #1","LMPP #2","GMP","CDP","pre-pDC","ETP","Thy #1","Thy #2","CLP","pro-B"))
DimPlot(thy, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
        repel = TRUE,raster = T)
thy$celltype <- thy$predicted.celltype
thymocyte.combined <- readRDS("/home/yushiya/data/cd34/data/2022_bioRxiv_thymus_7w-3yr/thymocyte_combined.rds")
anchors <- FindTransferAnchors(reference = immune.combined, query = thymocyte.combined,
                               dims = 1:30, reference.reduction = "pca")
thymocyte.combined <- MapQuery(anchorset = anchors, reference = immune.combined, query = thymocyte.combined,
                               refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap")
DimPlot(sub1, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
        repel = TRUE,raster = T)
sub1 <- subset(thymocyte.combined, subset= predicted.celltype=="ETP")
thymocyte.combined$celltype <- thymocyte.combined$predicted.celltype
saveRDS(thymocyte.combined,"/home/yushiya/data/cd34/data/2022_bioRxiv_thymus_7w-3yr/thymocyte_combined.rds")
thy_add <- merge(immune.combined,c(thy,thymocyte.combined))
thy_add <- subset(thy_add, subset= celltype=="ETP")
DimPlot(thy_add, reduction = "ref.umap", label = TRUE)

etp_list <- SplitObject(thy_add, split.by = "orig.ident")
etp.anchors <- FindIntegrationAnchors(object.list = etp_list, dims = 1:20)
etp <- IntegrateData(anchorset = etp.anchors, dims = 1:20)
DefaultAssay(etp) <- "integrated"
etp <- ScaleData(etp, verbose = FALSE)
etp <- RunPCA(etp, npcs = 30, verbose = FALSE)
ElbowPlot(etp)
etp <- FindVariableFeatures(etp, selection.method = "vst", nfeatures = 2000)
#Clustering
etp <- FindNeighbors(etp, reduction = "pca", dims = 1:18)
etp <- FindClusters(etp, resolution = 0.5)
etp <- RunUMAP(etp, reduction = "pca", dims = 1:18)
DimPlot(etp, reduction = "umap", split.by = "orig.ident")
DimPlot(etp, reduction = "umap", split.by = "stim")
DimPlot(etp, reduction = "umap", label = TRUE)
DefaultAssay(etp) <- "RNA"
etp <- ScaleData(etp, verbose = FALSE)
saveRDS(etp,"/home/yushiya/cd34/fig/Thy_ETP.rds")

markers <- FindAllMarkers(etp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, file = '/home/yushiya/cd34/fig/Thy_markers.csv')
top10 <- markers %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)
top10_gene <- top10$gene %>% unique()
cellInfo <- data.frame(celltype=etp$integrated_snn_res.0.5)
mtx <- data.frame(etp@assays[["RNA"]]@data[top10_gene,]) 
colnames(mtx) <- rownames(cellInfo)
top10_exp <- sapply(split(rownames(cellInfo), cellInfo$celltype),
                    function(cells) rowMeans(mtx[top10_gene,cells]))
top10_exp <- na.omit(top10_exp)
p1 <- pheatmap(top10_exp, #cluster_cols = F, #cluster_rows = F,  
               #clustering_method = "median",
               scale = "row",
               fontsize_row = 10, fontsize_col= 12, angle_col = 45,
               col=colorRampPalette(c("#000066","#339999" ,"yellow"))(50))
ggsave("/home/yushiya/cd34/fig/new/thy_heatmap.pdf", plot=p1, width=7, height=8)
Thy_GO=data.frame()
for (i in unique(markers$cluster)) {
  data=subset(markers, subset= cluster==i)
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
    Thy_GO=rbind(Thy_GO,go_res)
  }
}
Thy_GO <- Thy_GO[which(Thy_GO$qvalue <= 0.05),]
write.csv(Thy_GO, file = '/home/yushiya/cd34/fig/Thy_GO.csv')
top_GO <- Thy_GO %>% group_by(cluster) %>% top_n(n = -10, wt = p.adjust)
selected_GO <- top_GO[c(1,15,17,18,23,26,35,36,47,56,65,71,72,76,90,97,107,109,117,118,125,132,129,135,136,138,139,141,144,147),]
selected_GO <- top_GO[c(1,4,11,21,23,29,30,31,32,42,52,61,64,68,74,87,95,97,109,110,112,113,115,116,120,27,104,124),]
mtx_go <- data.frame(selected_GO$Description)
colnames(mtx_go)[1] <- "Description"
for (i in 0:11) {
  temp <- subset(selected_GO, subset= cluster==i)
  temp <- temp[,c("Description","GeneRatio")]
  split_num <- data.frame(t(data.frame(strsplit(temp$GeneRatio,"/"))))
  temp$GeneRatio <- as.numeric(split_num$X1)/as.numeric(split_num$X2)
  #temp$GeneRatio <- strsplit()
  mtx_go <- mtx_go %>% left_join(temp, by="Description")
  mtx_go$GeneRatio[is.na(mtx_go$GeneRatio)] <- 0
  colnames(mtx_go)[i+2] <- paste("Cluster",i)
}
mtx_go <- mtx_go[-c(9,13,14,10,11,16,17,8,12,21,23,25,19,20,27,28,6,24),]
rownames(mtx_go) <- mtx_go$Description
mtx_go <- mtx_go[,-1]
pheatmap(mtx_go, #cluster_cols = F, #cluster_rows = F,  
         clustering_method = "average",
         scale = "column",
         #border=F,
         #col=colorRampPalette(c("navy","firebrick3"))(50),
         fontsize_row = 12, fontsize_col= 12, angle_col = 45)

#Monocle find lineage
dyn.load("/home/wangjl/local/gdal-2.4.4/libgdal.so.20.5.4")
library(monocle3)
data <- GetAssayData(etp, assay = 'RNA', slot = 'counts')
cell_metadata <- etp@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
save(cds,file="/home/yushiya/cd34/fig/new/Thy_ETP.cds.RData")
load(file="/home/yushiya/data/blood_atlas/2022_Science_ImmuneMap/hema.combined.cds.RData")
#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, preprocess_method = "PCA")
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(etp, reduction = "umap")
cds@int_colData$reducedDims$UMAP <- int.embed
cds <- cluster_cells(cds)
cds@clusters@listData[["UMAP"]][["clusters"]] <- etp$integrated_snn_res.0.5
## 识别轨迹
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = FALSE, rasterize = T)
#DEGs along trajectory
cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=80)  
write.csv(cds_pr_test_res, "/home/yushiya/cd34/fig/new/ETP_peu_tra_gene.csv")
cds_pr_test_res <- read.csv("/home/yushiya/cd34/fig/new/ETP_peu_tra_gene.csv", row.names = 1)
pr_deg_ids <- row.names(subset(cds_pr_test_res, q_value < 0.01 & morans_I > 0.2))
#find modules
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=c(10^seq(-7,-2)))
table(gene_module_df$module)
write.csv(gene_module_df, "/home/yushiya/cd34/fig/new/ETP_peu_tra_module.csv")
gene_module_df <- read.csv("/home/yushiya/cd34/fig/new/ETP_peu_tra_module.csv", row.names = 1)
#new_type <- paste(colData(EM_cds)$stim,colData(EM_cds)$celltype,sep="_")
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=etp$integrated_snn_res.0.5)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,#cluster_rows = T, cluster_cols = T,
                   scale="column", clustering_method="ward.D2")
#GO terms in all modules
module_gene <- gene_module_df[,c(1,2)]
colnames(module_gene) <- c("gene","Module")
rownames(module_gene) <- module_gene$gene
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
write.csv(Module_GO, file = '/home/yushiya/cd34/fig/new/ETP_Module_GO_monocle3.csv')
top_Module_GO <- Module_GO %>% group_by(cluster) %>% top_n(n = -10, wt = p.adjust)
selected_Module_GO <- top_Module_GO[c(1,2,7,13,16,21,24,31,44),]
mtx_go <- data.frame(selected_Module_GO$Description)
colnames(mtx_go)[1] <- "Description"
for (i in 1:5) {
  temp <- subset(selected_Module_GO, subset= cluster==i)
  temp <- temp[,c("Description","GeneRatio")]
  split_num <- data.frame(t(data.frame(strsplit(temp$GeneRatio,"/"))))
  temp$GeneRatio <- as.numeric(split_num$X1)/as.numeric(split_num$X2)
  #temp$GeneRatio <- strsplit()
  mtx_go <- mtx_go %>% left_join(temp, by="Description")
  mtx_go$GeneRatio[is.na(mtx_go$GeneRatio)] <- 0
  colnames(mtx_go)[i+1] <- paste("Module",i)
}
mtx_go <- mtx_go[-c(4),]
rownames(mtx_go) <- mtx_go$Description
mtx_go <- mtx_go[,-1]
pheatmap(mtx_go, #cluster_cols = F, #cluster_rows = F,  
         clustering_method = "average",
         scale = "column",
         #border=F,
         #col=colorRampPalette(c("navy","firebrick3"))(50),
         fontsize_row = 12, fontsize_col= 12, angle_col = 45)

#calculate cell type distribution
library(reshape2)
library("ggalluvial")
clus_percent <- table(etp$integrated_snn_res.0.5,etp$stim)
clus_percent <- replace(clus_percent, clus_percent < 20, 0) 
clus_percent <- apply(clus_percent,2,function(x) prop.table(x))
clus_percent <- apply(clus_percent,1,function(x) prop.table(x))
clus_percent <- t(clus_percent[order(as.numeric(rownames(clus_percent))),])
clus_percent <- na.omit(clus_percent)
Stim=colnames(clus_percent)
clus_percent=melt(clus_percent, id='Stim')
names(clus_percent)[2]='Stim'
names(clus_percent)[1]='Cluster'
clus_percent$Stim <- factor(clus_percent$Stim, levels=c("BM","mPB","PB","thymus"))
clus_percent$Cluster <- factor(clus_percent$Cluster, levels=c("0","1","2","3","4","5","6","7","8","9","10","11"))
ggplot(clus_percent,
       aes(x=Cluster, y=value*100, fill=Stim)) +
  geom_bar(stat='identity', width=0.45) +
  #geom_alluvium() +
  #geom_stratum(width=0.45, size=0.1) +
  scale_fill_manual(values = color_stim) +
  labs(x='Cluster', y='Relative Abundance (%)')+
  scale_y_continuous(expand=c(0, 0))+
  theme(axis.text.x=element_text(hjust=0.5, angle=45, vjust=0.5),
        text = element_text(size = 12),
        panel.background = element_blank(),
        panel.grid = element_line(colour = "grey90"),
        axis.line = element_line(colour = "grey90",size = rel(1)))

sub_etp <- subset(etp, subset= integrated_snn_res.0.5=="11")
sub_etp <- subset(etp, subset= stim=="BM" | stim=="thymus")
sub_etp <- subset(etp, subset= integrated_snn_res.0.5=="2")
Idents(sub_etp) <- sub_etp$stim
markers_11 <- FindAllMarkers(sub_etp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, file = '/home/yushiya/cd34/fig/new/markers_2.csv')
markers_11 <- read.csv('/home/yushiya/cd34/fig/new/markers_11.csv')
top10 <- markers_11 %>% group_by(cluster) %>% top_n(n=20, wt=avg_log2FC)
top10_gene <- top10$gene %>% unique()
cellInfo <- data.frame(celltype=sub_etp$stim)
mtx <- data.frame(sub_etp@assays[["RNA"]]@data[top10_gene,]) 
colnames(mtx) <- rownames(cellInfo)
top10_exp <- sapply(split(rownames(cellInfo), cellInfo$celltype),
                    function(cells) rowMeans(mtx[top10_gene,cells]))
top10_exp <- na.omit(top10_exp)
pheatmap(top10_exp, cluster_cols = F, cluster_rows = F,  
         #clustering_method = "median",
         scale = "row",
         fontsize_row = 10, fontsize_col= 12, angle_col = 45,
         col=colorRampPalette(c("#000066","#339999" ,"yellow"))(50))
BMID <- row.names(subset(markers_11, subset= cluster=="BM"))
BM_gene <- bitr(BMID, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
mPBID <- row.names(subset(markers_11, subset= cluster=="mPB"))
mPB_gene <- bitr(mPBID, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
PBID <- row.names(subset(markers_11, subset= cluster=="PB"))
PB_gene <- bitr(PBID, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
thymusID <- row.names(subset(markers_11, subset= cluster=="thymus"))
thymus_gene <- bitr(thymusID, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
cp = list(BM.gene=BM_gene$ENTREZID, thymus.gene=thymus_gene$ENTREZID)
cp = list(BM.gene=BM_gene$ENTREZID, mPB_gene=mPB_gene$ENTREZID,PB_gene=PB_gene$ENTREZID, thymus.gene=thymus_gene$ENTREZID)
go.p <- compareCluster(cp,
                       fun = "enrichGO",
                       OrgDb = "org.Hs.eg.db",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.01
)
saveRDS(go.p, "/home/yushiya/cd34/fig/new/GO_2.rds")
go.p <- readRDS("/home/yushiya/cd34/fig/new/GO_11.rds")
go.p1 <- simplify(go.p,cutoff=0.4,by="p.adjust",select_fun=min)  #去除冗余
go.p1 <- go.p
go.p1@compareClusterResult <- go.p1@compareClusterResult[c(1,2,4,26,32,59,73,93,106,121,207,216,264,356,122,342),]
go.p1@compareClusterResult <- go.p1@compareClusterResult[c(1,13,30,133,137,188,199,204,224,228,231,342,327,349,350,331,356,367,427,443,446,535,525,302,548),]
ggplot(go.p1, aes(Cluster, Description), showCategory=15) +
  geom_point(aes(color=p.adjust, size=GeneRatio))+
  theme_classic()+
  theme(axis.text.x =element_text(angle = 45,vjust = 0.5,hjust = 0.5),
        text = element_text(size = 15),
        panel.background = element_blank(),
        panel.grid = element_line(colour = "grey90"))#+theme_bw()

#cellchat
library(CellChat)
library(svglite)
sub_tec <- readRDS("/home/yushiya/data/cd34/data/TEC_sub.rds")
sub1 <- readRDS("/home/yushiya/data/cd34/data/fig/ETP.rds")
sub1 <- subset(sub1, subset= T_type %in% c("TSP_BM","TSP_PB","TSP_thymus","ETP1"))
#Idents(sub_etp) <- sub_etp$stim
sum <- merge(sub_etp,sub_tec)
#sum <- merge(sub1,sub_tec)
Idents(sum) <- factor(Idents(sum),levels=c("TSP_BM","TSP_mPB","TSP_PB","TSP_thymus","ETP_IGLL1","ETP_TRBC1","ETP_RPS21","Thymic_Epi","Thymic_Mesen","Thymic_Endo"))
Idents(sum) <- factor(Idents(sum),levels=c("TSP_BM","TSP_mPB","TSP_PB","TSP_thymus","ETP1","Thymic_Epi","Thymic_Mesen","Thymic_Endo"))
p1=DotPlot(sum, features = c("CCL19","CCL21","CCL25","CXCL12","CXCL14","CCL2","CCL14","TNFSF10","TNFRSF4"), cols = c("blue","red"), dot.scale = 8) + 
  RotatedAxis()
DotPlot(sum, features = c("CCL19","CCL21","CCL25","CXCL12","CXCL14","CCL2","CCL14","TNFSF10","TNFRSF4"), cols = c("#3E4A89FF","#FDE725FF"), dot.scale = 8) + 
  RotatedAxis()
ggsave("/home/yushiya/data/cd34/data/fig/fig5-ETP-dot3.pdf", plot=p1, width=7, height=3.5)
cellchat <- createCellChat(sum)
#cellchat@idents <- factor(cellchat@idents, levels=c("Thymic_Epi","Thymic_Mesen","Thymic_Endo","TSP_BM","TSP_mPB","TSP_PB","TSP_thymus","ETP_IGLL1","ETP_TRBC1","ETP_RPS21"))
cellchat@idents <- factor(cellchat@idents, levels=c("Thymic_Epi","Thymic_Mesen","Thymic_Endo","TSP_BM","TSP_mPB","TSP_PB","TSP_thymus","ETP1"))
CellChatDB <- CellChatDB.human 
CellChatDB.use <- CellChatDB 
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)  
#识别细胞组中过度表达的配体或受体
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
mat <- cellchat@net$weight
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat2[c(1:3),c(4:10)] <- mat[c(1:3),c(4:10)]
mat2[c(1:3),c(4:8)] <- mat[c(1:3),c(4:8)]
netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = "Interaction weights/strength")
# 显示从某些细胞组到其他细胞组的所有显著的相互作用（L-R 对）
p1=netVisual_bubble(cellchat, sources.use = c(1:3), targets.use = c(4:8), remove.isolate = FALSE)
p1=netVisual_bubble(cellchat, sources.use = c(2), targets.use = c(4:8), signaling = c("CD99","MIF","CXCL","NOTCH","CCL","COLLAGEN"), remove.isolate = FALSE)
plotGeneExpression(cellchat, signaling = c("CXCL","MIF","PTN","NOTCH","CD99"))
color_ETP <- c("#66C2A5","#FC8D62","#8DA0CB","#DB5C25","#F3B747","#649541","#4C82C5","#A6D96A","#8DA0CB","#E78AC3")
p1=plotGeneExpression(cellchat, signaling = c("PTN","CD99","MIF","CXCL","NOTCH","CCL","COLLAGEN","ADGRE"),color.use=color_ETP)
#ggsave("/home/yushiya/data/cd34/fig_23/fig3-cellchat_select.pdf", plot=p1, width=5, height=2.5)
ggsave("/home/yushiya/data/cd34/data/fig/fig5-cellchat-ETP-LRpairs3.pdf", plot=p1, width=6, height=9)
ggsave("/home/yushiya/data/cd34/data/fig/fig5-cellchat-ETP-LRpairs-selected3.pdf", plot=p1, width=4, height=5)
ggsave("/home/yushiya/data/cd34/data/fig/fig5-cellchat-ETP-LRpairs-exp3.pdf", plot=p1, width=4, height=8)
saveRDS(cellchat, file = "/home/yushiya/data/cd34/data/fig/cellchat_ETP3.rds")
#"CXCL","APP","SELPLG","ANGPTL","ADGRE5"
pathways.show <- c("CXCL") 
netVisual_aggregate(cellchat, signaling = pathways.show, sources.use = c(1,2), targets.use = c(3:17))
netVisual_aggregate(cellchat, signaling = pathways.show, sources.use = c(1:5), targets.use = c(6:9))
netVisual_chord_cell(cellchat, signaling = pathways.show, title.name = paste0(pathways.show, " signaling network"))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds", sources.use = c(1,2), targets.use = c(3:16))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds", sources.use = c(1:5), targets.use = c(6:9))
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # “netP”:推断出的信号通路的细胞间通信网络
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 18, height = 8, font.size = 12, font.size.title = 18)
netAnalysis_contribution(cellchat, signaling = pathways.show)
ggsave("/home/yushiya/cd34/fig/fig5-cellchat-CXCL-PB-circle.pdf", plot=p1, width=6, height=6)
ggsave("/home/yushiya/cd34/fig/fig5-ETP-cellchat-CXCL-heatmap-3.pdf", plot=p1, width=6, height=10)
cellchat <- readRDS(file = "/home/yushiya/data/cd34/data/fig/cellchat_ETP3.rds")

#thymus visium data
library(SeuratDisk)
library(reticulate)
use_condaenv(condaenv = "/home/yushiya/miniconda3/envs/cell2location_cuda", required = TRUE)
ad <- import("anndata")
adata_vis_162 <- ad$read_h5ad('/home/yushiya/data/cd34/data/2024_thymus_Spatial/adata_CMA_162.h5ad')
colnames(adata_vis_162$X) <- adata_vis_162$var


