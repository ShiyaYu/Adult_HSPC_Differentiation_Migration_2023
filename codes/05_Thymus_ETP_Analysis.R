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
                 ## 删去x\y轴标�?
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
  geom_bar(stat = "identity", width = 0.8, fill = "#4C82C5",alpha = 0.8) + #绘制条形�?
  geom_text(aes(y = 0, #控制文本标签起始位置
                label = Description),
            size = 6,hjust = 0) + #hjust = 0左对�?
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
#differentialGeneTest 找高变基�?
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
peu_gene %>% arrange(qval)  -> peu_gene#按照qval排个�?
peu_gene <- peu_gene[1:1000,] #这里我们取前100个基因演�?
source('/home/yushiya/code/monocle2_heatmap.R')
p11 <- plot_pseudotime_heatmap(etp_monocle[peu_gene$gene_short_name,],
                               num_clusters = 3,
                               cores = 20, 
                               show_rownames = T,return_heatmap =T,
                               hmcols = viridis(256),
                               use_gene_short_name = T)
###首先提取热图中各个module的基�?
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
                 ## 删去x\y轴标�?
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
  geom_bar(stat = "identity", width = 0.8, fill = "firebrick3",alpha = 0.8) + #绘制条形�?
  geom_text(aes(y = 0, #控制文本标签起始位置
                label = Description),
            size = 6,hjust = 0) + #hjust = 0左对�?
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
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # “netP�?:推断出的信号通路的细胞间通信网络
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 18, height = 8, font.size = 12, font.size.title = 18)
netAnalysis_contribution(cellchat, signaling = pathways.show)
ggsave("/home/yushiya/cd34/fig/fig5-cellchat-CXCL-PB-circle.pdf", plot=p1, width=6, height=6)
ggsave("/home/yushiya/cd34/fig/fig5-ETP-cellchat-CXCL-heatmap-3.pdf", plot=p1, width=6, height=10)
cellchat <- readRDS(file = "/home/yushiya/data/cd34/data/fig/cellchat_ETP3.rds")
