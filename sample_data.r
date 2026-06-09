library(Seurat)
## 
# 加载Cell Ranger输出的过滤后的基因表达矩阵
data_dir <- "D:/单细胞测序/Seq2_3_Comparison/summary (2)/summary/1_Cellranger_result/A0_EV/filtered_feature_bc_matrix"
seurat_obj <- Read10X(data.dir = data_dir)
# 创建Seurat对象
seurat_obj<- CreateSeuratObject(counts = seurat_obj)
seurat_obj
# 提取基因表达矩阵
expression_matrix <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
##dim(expression_matrix)
##str(expression_matrix)
head(expression_matrix)
sample_expression_matrix<- rowSums(expression_matrix)/ncol(expression_matrix)
# 查看新的样本表达矩阵
head(sample_expression_matrix)
sample_expression_matrix<-t(t(sample_expression_matrix))
colnames(sample_expression_matrix)<-c("A0_EV")
colnames(sample_expression_matrix)
A0_EV=sample_expression_matrix
head(A0_EV)
##
# 加载Cell Ranger输出的过滤后的基因表达矩阵
data_dir <- "D:/单细胞测序/Seq2_3_Comparison/summary (2)/summary/1_Cellranger_result/A1_EV/filtered_feature_bc_matrix"
seurat_obj <- Read10X(data.dir = data_dir)
# 创建Seurat对象
seurat_obj<- CreateSeuratObject(counts = seurat_obj)
seurat_obj
# 提取基因表达矩阵
expression_matrix <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
str(expression_matrix)
head(expression_matrix)
sample_expression_matrix<- rowSums(expression_matrix)/ncol(expression_matrix)
# 查看新的样本表达矩阵
head(sample_expression_matrix)
sample_expression_matrix<-t(t(sample_expression_matrix))
colnames(sample_expression_matrix)<-c("A1_EV")
colnames(sample_expression_matrix)
A_EV=sample_expression_matrix
head(A1_EV)
##
# 加载Cell Ranger输出的过滤后的基因表达矩阵
data_dir <- "D:/单细胞测序/Seq2_3_Comparison/summary (2)/summary/1_Cellranger_result/A2_EV/filtered_feature_bc_matrix"
seurat_obj <- Read10X(data.dir = data_dir)
# 创建Seurat对象
seurat_obj<- CreateSeuratObject(counts = seurat_obj)
seurat_obj
# 提取基因表达矩阵
expression_matrix <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
str(expression_matrix)
head(expression_matrix)
sample_expression_matrix<- rowSums(expression_matrix)/ncol(expression_matrix)
# 查看新的样本表达矩阵
head(sample_expression_matrix)
sample_expression_matrix<-t(t(sample_expression_matrix))
colnames(sample_expression_matrix)<-c("A2_EV")
colnames(sample_expression_matrix)
A2_EV=sample_expression_matrix
head(A2_EV)
unique(A2_EV)
##
# 加载Cell Ranger输出的过滤后的基因表达矩阵
data_dir <- "D:/单细胞测序/Seq2_3_Comparison/summary (2)/summary/1_Cellranger_result/DBCO/filtered_feature_bc_matrix"
seurat_obj <- Read10X(data.dir = data_dir)
# 创建Seurat对象
seurat_obj<- CreateSeuratObject(counts = seurat_obj)
seurat_obj
# 提取基因表达矩阵
expression_matrix <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
str(expression_matrix)
head(expression_matrix)
sample_expression_matrix<- rowSums(expression_matrix)/ncol(expression_matrix)
# 查看新的样本表达矩阵
head(sample_expression_matrix)
sample_expression_matrix<-t(t(sample_expression_matrix))
colnames(sample_expression_matrix)<-c("DBCO")
colnames(sample_expression_matrix)
DBCO=sample_expression_matrix
head(DBCO)
unique(DBCO)
##
# 加载Cell Ranger输出的过滤后的基因表达矩阵
data_dir <- "D:/单细胞测序/Seq2_3_Comparison/summary (2)/summary/1_Cellranger_result/Stroke/filtered_feature_bc_matrix"
seurat_obj <- Read10X(data.dir = data_dir)
# 创建Seurat对象
seurat_obj<- CreateSeuratObject(counts = seurat_obj)
seurat_obj
# 提取基因表达矩阵
expression_matrix <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
str(expression_matrix)
head(expression_matrix)
sample_expression_matrix<- rowSums(expression_matrix)/ncol(expression_matrix)
# 查看新的样本表达矩阵
head(sample_expression_matrix)
sample_expression_matrix<-t(t(sample_expression_matrix))
colnames(sample_expression_matrix)<-c("Stroke")
colnames(sample_expression_matrix)
Stroke=sample_expression_matrix
head(Stroke)
unique(Stroke)
##
# 加载Cell Ranger输出的过滤后的基因表达矩阵
data_dir <- "D:/单细胞测序/Seq2_3_Comparison/summary (2)/summary/1_Cellranger_result/T_EV/filtered_feature_bc_matrix"
seurat_obj <- Read10X(data.dir = data_dir)
# 创建Seurat对象
seurat_obj<- CreateSeuratObject(counts = seurat_obj)
seurat_obj
# 提取基因表达矩阵
expression_matrix <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
str(expression_matrix)
head(expression_matrix)
sample_expression_matrix<-rowSums(expression_matrix)/ncol(expression_matrix)
# 查看新的样本表达矩阵
head(sample_expression_matrix)
sample_expression_matrix<-t(t(sample_expression_matrix))
colnames(sample_expression_matrix)<-c("T_EV")
colnames(sample_expression_matrix)
T_EV=sample_expression_matrix
head(T_EV)
unique(T_EV)

##
# 将它们合并成一个新的数据框（保留列名）
merged_matrix <- cbind(A0_EV, A1_EV, A2_EV,DBCO,Stroke,T_EV)

# 检查合并后的矩阵
print(merged_matrix)
