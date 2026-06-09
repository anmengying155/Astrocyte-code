#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("singscore")
#BiocManager::install("GSVA")
#BiocManager::install("GSEABase")
#BiocManager::install("limma")
#BiocManager的包安装命令可能会改变，需要去BiocManager官网查询最新安装方法，短时间内不会改变的。
#install.packages("devtools")
#library(devtools)
#devtools::install_github('dviraran/SingleR')
#这个包作用为对细胞的注释，对富集分析作用很关键，在这一篇章可以先不安装，在以后的go富集分析才需要。
正文代码开始：`#读取数据
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)

setwd("数据所在目录")             #设置工作目录

#读取文件，并对重复基因取均值
rt=read.table("geneMatrix.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)#如果基因存在多行中，为了取得结果就会取均值

#将矩阵转换为Seurat对象，并对数据进行过滤
pbmc <- CreateSeuratObject(counts = data,project = "seurat", min.cells = 3, min.features = 50, names.delim = "_",)#min.cells为基因存在样本最小数，需要根据实际情况选择，min.features = 50基因最小存在细胞数
#使用PercentageFeatureSet函数计算线粒体基因的百分比
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
pdf(file="04.featureViolin.pdf",width=10,height=6)           #保存对于基因特征的小提琴图
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
pbmc <- subset(x = pbmc, subset = nFeature_RNA > 50 & percent.mt < 7)    #对数据进行过滤，线粒体比例大于7%进行删掉，因为线粒体基因比例过高的细胞，会干扰细胞分群

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
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)#top10是对基因
CombinePlots(plots = list(plot1, plot2))
dev.off()`
