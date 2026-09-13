library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(Matrix)
library(ggplot2)
library(monocle3)

#Initialize CDS Object
immune.combined <- readRDS("./data/immune.combined-celltype.rds")
expr_data <- GetAssayData(immune.combined, assay = 'RNA', slot = 'counts')
gene_anno <- data.frame(gene_short_name = rownames(expr_data), row.names = rownames(expr_data))
cds <- new_cell_data_set(
  expr_data,
  cell_metadata = immune.combined@meta.data,
  gene_metadata = gene_anno
)
#Preprocess & Align Embeddings
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, preprocess_method = "PCA")
#Transfer UMAP from Seurat
cds@int_colData$reducedDims$UMAP <- Embeddings(immune.combined, reduction = "umap")
cds <- cluster_cells(cds)
cds@clusters@listData[["UMAP"]][["clusters"]] <- immune.combined$celltype
#Learn Graph & Order Cells
cds <- learn_graph(cds, learn_graph_control=list(geodesic_distance_ratio=0.5))
# Define Root Cell mathematically (e.g., based on HSC cluster UMAP coordinates)
embed <- data.frame(Embeddings(immune.combined, reduction = "umap"))
root.cell <- rownames(subset(embed, UMAP_1 > 5.5 & UMAP_1 < 6.5 & UMAP_2 > 4.2 & UMAP_2 < 5.5))
cds <- order_cells(cds, root_cells = root.cell)
# Append pseudotime back to Seurat object
immune.combined$pseudotime <- pseudotime(cds)
# Trajectory Differential Expression (Graph Test)
sub_cds <- cds[, colData(cds)$stim %in% c("BM","mPB","PB")]
cds_pr_test_res <- graph_test(sub_cds, neighbor_graph="principal_graph", cores=8)
write.csv(cds_pr_test_res, "./results/Pseudotime_DEGs.csv")

#Gene Modules
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
# 计算Jaccard相似�?
calculate_jaccard <- function(ids1, ids2) {
  intersection <- length(intersect(unlist(ids1), unlist(ids2)))
  union <- length(union(unlist(ids1), unlist(ids2)))
  intersection / union
}
# 生成所有module对并计算相似�?
modules1 <- grouped_ids1$module
modules2 <- grouped_ids2$module
# 创建模块对的组合
module_combinations <- expand.grid(BM = modules1, mPB = modules2)
# 计算每个模块对的相似�?
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
# 计算Jaccard相似�?
calculate_jaccard <- function(ids1, ids2) {
  intersection <- length(intersect(unlist(ids1), unlist(ids2)))
  union <- length(union(unlist(ids1), unlist(ids2)))
  intersection / union
}
# 生成所有module对并计算相似�?
modules1 <- grouped_ids1$module
modules2 <- grouped_ids2$module
modules3 <- grouped_ids3$module
# 创建模块对的组合
module_combinations <- expand.grid(BM = modules1, mPB = modules2, PB = modules3)
# 计算每个模块对的相似�?
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
# 计算三个数组的交集相似�?
module_combinations$jaccard_similarity_all <- apply(module_combinations, 1, function(row) {
  ids1 <- grouped_ids1 %>% filter(module == row[1]) %>% pull(ids)
  ids2 <- grouped_ids2 %>% filter(module == row[2]) %>% pull(ids)
  ids3 <- grouped_ids3 %>% filter(module == row[3]) %>% pull(ids)
  if (length(ids1) == 0 || length(ids2) == 0 || length(ids3) == 0) {
    return(NA)
  }
  # 计算三者的交集和并�?
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
# 计算每个集合的比�?
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
  lty = "blank"  #不显示圆圈颜�?
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
                 ## 删去x\y轴标�?
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
  geom_bar(stat = "identity", width = 0.8, fill = "#649541",alpha = 0.8) + #绘制条形�?
  geom_text(aes(y = 0, #控制文本标签起始位置
                label = Description),
            size = 6,hjust = 0) + #hjust = 0左对�?
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
