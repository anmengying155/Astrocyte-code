if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install("singscore")
BiocManager::install("GSVA")
BiocManager::install("GSEABase")
BiocManager::install("limma")
 
install.packages("devtools")
library(devtools)
devtools::install_github('dviraran/SingleR')
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
setwd("D:/File/document/单细胞测序代码/step6")             #设置工作目录
 
#读取文件，并对重复基因取均值
rt=read.table("geneMatrix.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)



#将矩阵转换为Seurat对象，并对数据进行过滤
pbmc <- CreateSeuratObject(counts = data,project = "seurat", min.cells = 3, min.features = 50, names.delim = "_",)

#使用PercentageFeatureSet函数计算线粒体基因的百分比
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^mt-")
pdf(file="04.featureViolin.pdf",width=10,height=6)           #保存基因特征小提琴图
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
pbmc <- subset(x = pbmc, subset = nFeature_RNA > 50 & percent.mt < 5)    #对数据进行过滤
 
#测序深度的相关性绘图
pdf(file="04.featureCor.pdf",width=10,height=6)              #保存基因特征相关性图
plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5)
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",,pt.size=1.5)
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#对数据进行标准化
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#提取那些在细胞间变异系数较大的基因
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 1500)
#输出特征方差图
top10 <- head(x = VariableFeatures(object = pbmc), 10)
pdf(file="04.featureVar.pdf",width=10,height=6)              #保存基因特征方差图
plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
##PCA分析
pbmc=ScaleData(pbmc)                     #PCA降维之前的标准预处理步骤
pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc))     #PCA分析
 
#绘制每个PCA成分的相关基因
pdf(file="05.pcaGene.pdf",width=10,height=8)
VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca",nfeatures = 20)
dev.off()
 
#主成分分析图形
pdf(file="05.PCA.pdf",width=6.5,height=6)
DimPlot(object = pbmc, reduction = "pca")
dev.off()
 
#主成分分析热图
pdf(file="05.pcaHeatmap.pdf",width=10,height=8)
DimHeatmap(object = pbmc, dims = 1:4, cells = 500, balanced = TRUE,nfeatures = 30,ncol=2)#1：4是热图分成几个，ncol是分为几列。
dev.off()
 
#每个PC的p值分布和均匀分布
pbmc <- JackStraw(object = pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20)
pdf(file="05.pcaJackStraw.pdf",width=8,height=6)
JackStrawPlot(object = pbmc, dims = 1:20)
dev.off()
#
#
##TSNE聚类分析
pcSelect=20
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)                #计算邻接距离
pbmc <- FindClusters(object = pbmc, resolution = 0.5)                  #对细胞分组,优化标准模块化
pbmc <- RunTSNE(object = pbmc, dims = 1:pcSelect)                      #TSNE聚类
pdf(file="06.TSNE.pdf",width=6.5,height=6)
DimPlot(object = pbmc, reduction = "tsne", label = TRUE, pt.size = 2)    #TSNE可视化
dev.off()
write.table(pbmc$seurat_clusters,file="06.tsneCluster.txt",quote=F,sep="\t",col.names=F)
 
##寻找差异表达的特征
logFCfilter=0.5
adjPvalFilter=0.05
pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter)
sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_logFC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.markers,file="06.markers.xls",sep="\t",row.names=F,quote=F)
head(pbmc.markers)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#绘制marker在各个cluster的热图
pdf(file="06.tsneHeatmap.pdf",width=12,height=9)
DoHeatmap(object = pbmc, features = top10$gene) + NoLegend()
dev.off()
 
#绘制marker的小提琴图
pdf(file="06.markerViolin.pdf",width=10,height=6)
VlnPlot(object = pbmc, features = c("IGLL5", "MBOAT1"))
dev.off()
 
#绘制marker在各个cluster的散点图
pdf(file="06.markerScatter.pdf",width=10,height=6)
FeaturePlot(object = pbmc, features = c("IGLL5", "MBOAT1"),cols = c("green", "red"))
dev.off()
 
#绘制marker在各个cluster的气泡图
pdf(file="06.markerBubble.pdf",width=12,height=6)
cluster10Marker=c("MBOAT1", "NFIB", "TRPS1", "SOX4", "CNN3", "PIM2", "MZB1", "MS4A1", "ELK2AP", "IGLL5")
DotPlot(object = pbmc, features = cluster10Marker)
dev.off()
 
#
#
#这一章的代码从这里开始
 
 
 ###################################07.注释细胞类型###################################
library(SingleR)
library(celldex)
reference <- celldex::MouseRNAseqData()
# 获取表达矩阵
counts <- GetAssayData(pbmc, slot = "counts")
# 获取Seurat中的聚类结果作为标签
clusters <- pbmc@meta.data$seurat_clusters
# 加载小鼠参考数据集（假设你使用celldex包）
mouse_reference <- celldex::MouseRNAseqData()
# 使用SingleR进行注释
singler_result <- SingleR(test = counts, 
                          ref = mouse_reference,  # 使用celldex包中的小鼠参考数据集
                          labels = mouse_reference$label.main)  # 这里使用Seurat的聚类信息作为标签

# 获取SingleR注释结果
clusterAnn = singler_result$labels
# 将注释结果写入文件
write.table(clusterAnn, file="07.clusterAnn.txt", quote=FALSE, sep="\t", col.names=FALSE)
write.table(singler_result, file="07.cellAnn.txt", quote=FALSE, sep="\t", col.names=FALSE)


#准备monocle分析需要的文件
monocle.matrix=as.matrix(pbmc@assays$RNA@data)
monocle.matrix=cbind(id=row.names(monocle.matrix),monocle.matrix)
write.table(monocle.matrix,file="07.monocleMatrix.txt",quote=F,sep="\t",row.names=F)
monocle.sample=as.matrix(pbmc@meta.data)
monocle.sample=cbind(id=row.names(monocle.sample),monocle.sample)
write.table(monocle.sample,file="07.monocleSample.txt",quote=F,sep="\t",row.names=F)
monocle.geneAnn=data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
monocle.geneAnn=cbind(id=row.names(monocle.geneAnn),monocle.geneAnn)
write.table(monocle.geneAnn,file="07.monocleGene.txt",quote=F,sep="\t",row.names=F)
write.table(singler$other,file="07.monocleClusterAnn.txt",quote=F,sep="\t",col.names=F)
write.table(sig.markers,file="07.monocleMarkers.txt",sep="\t",row.names=F,quote=F)
 

# 1. 提取原始count矩阵
monocle.matrix = as.matrix(GetAssayData(pbmc, slot = "counts"))
monocle.matrix = cbind(id = rownames(monocle.matrix), monocle.matrix)
write.table(monocle.matrix, file="07.monocleMatrix.txt", quote=F, sep="\t", row.names=F)
# 2. 提取meta信
monocle.sample = as.matrix(pbmc@meta.data)
monocle.sample = cbind(id = rownames(monocle.sample), monocle.sample)
write.table(monocle.sample, file="07.monocleSample.txt", quote=F, sep="\t", row.names=F)
# 3. 基因注释信息
monocle.geneAnn = data.frame(gene_short_name = rownames(GetAssayData(, slot = "counts")),
                             row.names = rownames(GetAssayData(pbmc, slot = "counts")))
monocle.geneAnn = cbind(id = rownames(monocle.geneAnn), monocle.geneAnn)
write.table(monocle.geneAnn, file="07.monocleGene.txt", quote=F, sep="\t", row.names=F)
# 4. 其他辅助文件（如果需要）
write.table(singler$other, file="07.monocleClusterAnn.txt", quote=F, sep="\t", col.names=F)
write.table(sig.markers, file="07.monocleMarkers.txt", sep="\t", row.names=F, quote=F)



 ###################################07.注释细胞类型###################################
library(SingleR)
library(celldex)
reference <- celldex::MouseRNAseqData()
# 获取表达矩阵
counts <- GetAssayData(seurat_obj, slot = "counts")
# 获取Seurat中的聚类结果作为标签
clusters <- seurat_obj@meta.data$seurat_clusters
# 加载小鼠参考数据集（假设你使用celldex包）
mouse_reference <- celldex::MouseRNAseqData()
# 使用SingleR进行注释
singler_result <- SingleR(test = counts, 
                          ref = mouse_reference,  # 使用celldex包中的小鼠参考数据集
                          labels = mouse_reference$label.main)  # 这里使用Seurat的聚类信息作为标签

# 获取SingleR注释结果
clusterAnn = singler_result$labels
seurat_obj$singler_labels <- singler_result$labels
# 将注释结果写入文件
write.table(clusterAnn, file="07.clusterAnn.txt", quote=FALSE, sep="\t", col.names=FALSE)
write.table(singler_result, file="07.cellAnn.txt", quote=FALSE, sep="\t", col.names=FALSE)

# -------------------------------
# UMAP 可视化 MR 注释
# -------------------------------
pdf(file="Labels.pdf",width=10,height=6)
DimPlot(seurat_obj, reduction = "umap", group.by = "singler_labels", label = FALSE, pt.size = 2) + 
  ggtitle("Cell Types - MR Labels") + coord_fixed() +
  theme(text = element_text(size = 12))  # 调整字体大小
dev.off()
##寻找差异表达的特征
logFCfilter=0.5
adjPvalFilter=0.05
seurat_obj.markers <- FindAllMarkers(object = seurat_obj,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter)
head(seurat_obj.markers$cluster)
unique(seurat_obj.markers$cluster)
# 查看seurat_obj.markers的前几行
head(seurat_obj.markers)
# 看一眼列名
colnames(seurat_obj.markers)
# 统计一下符合条件的基因数目
sum(abs(as.numeric(as.vector(seurat_obj.markers$avg_log2FC))) > logFCfilter)
sum(as.numeric(as.vector(seurat_obj.markers$p_val_adj)) < adjPvalFilter)

sig.markers=seurat_obj.markers[(abs(as.numeric(as.vector(seurat_obj.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(seurat_obj.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.markers,file="06.markers.xls",sep="\t",row.names=F,quote=F)
head(seurat_obj.markers)
top10 <- seurat_obj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#绘制marker在各个cluster的热图
pdf(file="06.tsneHeatmap.pdf",width=12,height=9)
DoHeatmap(object =seurat_obj, features = top10$gene) + NoLegend()
dev.off()
#绘制marker的小提琴图
pdf(file="06.markerViolin.pdf",width=10,height=6)
VlnPlot(object = seurat_obj, features = c("IGLL5", "MBOAT1"))
dev.off()
#绘制marker在各个cluster的散点图
pdf(file="06.markerScatter.pdf",width=10,height=6)
FeaturePlot(object = seurat_obj, features = c("IGLL5", "MBOAT1"),cols = c("green", "red"))
dev.off()
#绘制marker在各个cluster的气泡图
pdf(file="06.markerBubble.pdf",width=12,height=6)
cluster10Marker=c("MBOAT1", "NFIB", "TRPS1", "SOX4", "CNN3", "PIM2", "MZB1", "MS4A1", "ELK2AP", "IGLL5")
DotPlot(object = seurat_obj, features = cluster10Marker)
dev.off()
 


# 1. 提取原始count矩阵
monocle.matrix = as.matrix(GetAssayData(seurat_obj, slot = "counts"))
monocle.matrix = cbind(id = rownames(monocle.matrix), monocle.matrix)
write.table(monocle.matrix, file="07.monocleMatrix.txt", quote=F, sep="\t", row.names=F)
# 2. 提取meta信
monocle.sample = as.matrix(seurat_obj@meta.data)
monocle.sample = cbind(id = rownames(monocle.sample), monocle.sample)
write.table(monocle.sample, file="07.monocleSample.txt", quote=F, sep="\t", row.names=F)
# 3. 基因注释信息
monocle.geneAnn = data.frame(gene_short_name = rownames(GetAssayData(seurat_obj, slot = "counts")),
                             row.names = rownames(GetAssayData(seurat_obj, slot = "counts")))
monocle.geneAnn = cbind(id = rownames(monocle.geneAnn), monocle.geneAnn)
write.table(monocle.geneAnn, file="07.monocleGene.txt", quote=F, sep="\t", row.names=F)
# 4. 其他辅助文件（如果需要）
write.table(clusterAnn, file="07.monocleClusterAnn.txt", quote=F, sep="\t", col.names=F)
write.table(sig.markers, file="07.monocleMarkers.txt", sep="\t", row.names=F, quote=F)
