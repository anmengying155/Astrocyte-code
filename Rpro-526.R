{
  # 合并所有Seurat对象
  combined_seurat <- merge(x = seurat_list[[1]], 
                           y = seurat_list[-1], 
                           add.cell.ids = names(seurat_list))
  
  
  
  combined_seurat<-NormalizeData(combined_seurat,
                                 normalization.method="LogNormalize",
                                 scale.factor=10000)
  
  #鉴定2000个高变基因(数量可人为设置，一般是2K)
  combined_seurat<-FindVariableFeatures(combined_seurat,
                                        selection.method="vst",
                                        nfeatures=2000)
  #找出10个高变的基因
  top10<-head(VariableFeatures(combined_seurat),10)
  
  #对高可变基因进行可视化
  plot1<-VariableFeaturePlot(combined_seurat)
  plot2<-LabelPoints(plot=plot1,points=top10,
                     repel=TRUE);plot2
  ggsave("Variablegenes.png",plot=plot2,width=8,height=5,dpi=300)
  
  combined_seurat<-ScaleData(combined_seurat,features= VariableFeatures(object=combined_seurat))
  
  combined_seurat<-RunPCA(combined_seurat,features=VariableFeatures(object=combined_seurat),)
  
  
  combined_seurat <- RunHarmony(combined_seurat, group.by = 'orig.ident')
  
  DimPlot(combined_seurat,reduction="harmony",group.by="orig.ident")
  
  ElbowPlot(combined_seurat,ndims=50)
  ggsave("ElbowPlot.png",width=8,height=6,dpi=300)
  
  
  
  combined_seurat<-FindNeighbors(combined_seurat,dims=1:10,
                                 reduction="harmony",#单样本需要pca
                                 verbose=FALSE)
  combined_seurat <- FindClusters(
    combined_seurat, 
    resolution = 0.2
  )
  
  head(Idents(combined_seurat), 5)
  
  library(raster)
  library(grid)
  #假设电脑的显卡非常高级的话，可以不用PCA降维，直接UMAP
  #Umap方式
  combined_seurat <- RunUMAP(combined_seurat, reduction = "harmony", dims = 1:20)
  DimPlot(combined_seurat,label=T,raster=FALSE)
  ggsave("UMAP_plot.pdf", plot = UMAP_1, width = 20, height = 20) 
  table(Idents(combined_seurat))
  # 删除群体9到13
  combined_seurat1 <- subset(combined_seurat, idents = setdiff(levels(Idents(combined_seurat)), c("9", "10", "11", "12", "13")))
  DimPlot(combined_seurat1,label=T,raster=FALSE)
  
  saveRDS(combined_seurat, file = "/home/zbh1/combined_seurat.rds")
  saveRDS(my_data, file = "my_data.rds")
}#数据处理
{
  library(scMetabolism)
  library(ggplot2)
  library(rsvd)
  countexp.Seurat<-sc.metabolism.Seurat(obj = reaname, method = "AUCell", 
                                        imputation = F, ncores = 2, 
                                        metabolism.type = "KEGG")
  metabolism <- reaname@assays$METABOLISM$score
  Metabolismfig<-DotPlot.metabolism(
    obj = countexp.Seurat, 
    pathway =metabolism_pathways, 
    phenotype = "celltype", 
    norm = "y"
  )
  
  pdf(" ",family = 'sans', width = 12, height = 6)
  print( Metabolismfig)
  dev.off()
  
}#scmetabolism
{
  library(Startrac)
  library(ggplot2)
  library(tictoc)
  library(ggpubr)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(circlize)
  library(tidyverse)
  library(sscVis)
  R_oe<-calTissueDist(dat,
                      byPatient=F,
                      colname.cluster="celltype",
                      colname.patient="sample",
                      colname.tissue="RS",
                      method="chisq",
                      min.rowSum=0)
  R_oe
  
  library(ComplexHeatmap)
  library(grid)
  
  # 假设 col_fun 根据 "+/-" 定义颜色
  # 定义颜色映射函数
  # 定义固定颜色映射函数
  col_fun <- function(x) {
    ifelse(x > 1, "#FF9999", "#99CCFF")  # 大于1为红色，小于等于1为蓝色
  }
  
  
  # 自定义新的 Heatmap
  cell_fun = function(j, i, x, y, width, height, fill) {
    label <- ifelse(R_oe[i, j] > 1, "+", "+/-")  # 定义符号
    label_color <- ifelse(R_oe[i, j] > 1, "#FF9999", "#99CCFF")  # 定义符号颜色
    grid.text(label, x, y, gp = gpar(fontsize = 8, col = label_color))  # 绘制
  }
  
  # 绘制热图
  Heatmap(
    as.matrix(R_oe),,
    width = unit(4, "cm"),
    show_heatmap_legend = FALSE,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    row_names_side = "right",
    show_column_names = TRUE,
    show_row_names = TRUE,
    col = col_fun,  # 使用颜色映射函数
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 10),
    heatmap_legend_param = list(
      title = "Category",
      at = c(0.5, 1.5),  # 图例中的刻度
      labels = c("+/- (R_o/e ≤ 1)", "+ (R_o/e > 1)"),  # 图例的描述
      legend_gp = gpar(fill = c("blue", "red"))  # 图例颜色
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      # 根据 R_oe 值定义符号
      label <- ifelse(R_oe[i, j] > 1, "+", "+/-")
      # 绘制符号
      grid.text(label, x, y, gp = gpar(fontsize = 8))
    }
  )
  
  }#Ro/e
{  
  library(CellChat)
  library(ggalluvial)
  library(Seurat)
  { 
    reaname<-reanamenew
    celltype_counts <- table(reaname$celltype)
    table(Idents(reaname))
    # 筛选出细胞数大于等于1000的celltype
    valid_celltypes <- names(celltype_counts[celltype_counts >= 600])
    
    # 筛选Seurat对象，保留这些celltype
    combined_low <- subset(combined_low, cells = WhichCells(combined_low, idents = valid_celltypes))
    combined_high <- subset(combined_high, cells = WhichCells(combined_high, idents = valid_celltypes))
    
    # 对combined_high进行数据标准化
    combined_high <- NormalizeData(combined_high,
                                   normalization.method="LogNormalize",
                                   scale.factor=10000)
    
    # 鉴定2000个高变基因(数量可人为设置，一般是2K)
    combined_high <- FindVariableFeatures(combined_high,
                                          selection.method="vst",
                                          nfeatures=2000)
    
    # 对高变基因进行数据缩放
    combined_high <- ScaleData(combined_high, features=rownames(combined_high))
    
    # 对combined_low进行数据标准化
    combined_low <- NormalizeData(combined_low,
                                  normalization.method="LogNormalize",
                                  scale.factor=10000)
    
    # 鉴定2000个高变基因(数量可人为设置，一般是2K)
    combined_low <- FindVariableFeatures(combined_low,
                                         selection.method="vst",
                                         nfeatures=2000)
    
    # 对高变基因进行数据缩放
    combined_low <- ScaleData(combined_low, features=rownames(combined_low))
    cellchat <- createCellChat(object = combined_high,
                               meta = combined_high@meta.data,
                               group.by = "celltype")
    cellchat
    CellChatDB <- CellChatDB.human  
    showDatabaseCategory(CellChatDB)
    CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling", "ECM-Receptor" ,"Cell-Cell Contact"))
    cellchat@DB <- CellChatDB.use
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- projectData(cellchat, PPI.human)
    cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE) 
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    chat_high <- cellchat
    
    
    
    cellchat <- createCellChat(object = combined_low,
                               meta = combined_low@meta.data,
                               group.by = "celltype")
    cellchat
    CellChatDB <- CellChatDB.human  
    showDatabaseCategory(CellChatDB)
    CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling", "ECM-Receptor" ,"Cell-Cell Contact"))
    cellchat@DB <- CellChatDB.use
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- projectData(cellchat, PPI.human)
    cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE) 
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    chat_low <- cellchat
    # 保存 CellChat 对象为 RDS 文件
    saveRDS(chat_high, file = "chat_high.rds")
    # 保存 CellChat 对象为 RDS 文件
    saveRDS(chat_low, file = "chat_low.rds")
    
  }#cellchatsteps
  
  
  }#cellchat
{
  # 提取"T cell"的细胞
  t_cell_subset1<-t_cell_subset
  
  t_cell_subset <- subset(tn_cell_subset, idents = "Tcell")
  t_cell_subset <- subset(t_cell_subset,subset = nFeature_RNA > 800 & nFeature_RNA < 6000 & percent.mt < 40)
  
  t_cell_subset <- NormalizeData(t_cell_subset,
                                 normalization.method="LogNormalize",
                                 scale.factor=10000)
  
  # 鉴定2000个高变基因（数量可人为设置，一般是2K）
  t_cell_subset <- FindVariableFeatures(t_cell_subset,
                                        selection.method="vst",
                                        nfeatures=2000)
  t_cell_subset <- ScaleData(t_cell_subset, features=VariableFeatures(object=t_cell_subset))
  
  t_cell_subset <- RunPCA(t_cell_subset, features=VariableFeatures(object=t_cell_subset), verbose=FALSE)
  
  ElbowPlot(t_cell_subset, ndims=50)
  t_cell_subset <- FindNeighbors(t_cell_subset, dims=1:30, reduction="harmony", verbose=FALSE)
  
  t_cell_subset <- FindClusters(t_cell_subset, resolution=1, verbose=FALSE)
  t_cell_subset <- RunUMAP(t_cell_subset, reduction="harmony", dims=1:30)
  UMAP_1 <- DimPlot(t_cell_subset, label=TRUE); 
  print(UMAP_1)
  table(Idents(t_cell_subset))
  
  FeaturePlot(t_cell_subset, features = c("nFeature_RNA", 'nCount_RNA',"percent.mt"))
  DimPlot(t_cell_subset, reduction = "umap", label = T, pt.size = 0.3)
  markers_4clusters <- FindAllMarkers(t_cell_subset, only.pos = TRUE)
  FeaturePlot(t_cell_subset,features = 'PLCG2',label = F)
  
  
  {# 1. 构建 marker 列表
    CD4_markers_list <- list(
      Texaust = c('TIGIT', 'CTLA4',' PDCD1', 'HAVCR2', 'LAG3'),
      CD4    = c("CD4"),
      CD8 = c('CD8A','CD8B'),
      Treg   = c("TNFRSF4", "BATF", "TNFRSF18", "FOXP3", "IL2RA", "IKZF2"),
      naive  = c("CCR7", "SELL", "CD5",'GPR183'),
      Tfh    = c("CXCR5", "BCL6", "ICA1", "TOX", "TOX2", "IL6ST"),
      Tcm    = c('HSPA6', 'MT1E', 'MT1F'),
      Tem    = c('RPS19')
    )
    
    # 2. 展平 gene 向量 + 分组信息
    marker_genes <- unlist(CD4_markers_list)
    marker_groups <- rep(names(CD4_markers_list), times = sapply(CD4_markers_list, length))
    
    # 创建 data.frame 便于匹配分组
    marker_df <- data.frame(
      gene = marker_genes,
      group = marker_groups
    )
    dp <- DotPlot(t_cell_subset, features = marker_genes)
    
    # 5. 添加分组信息（将 group 加入数据框）
    dp$data$group <- marker_df$group[match(dp$data$features.plot, marker_df$gene)]
    
    # 6. 重新绘图，并按功能 marker 分组分面展示
    library(ggplot2)
    p <- ggplot(dp$data, aes(x = id, y = features.plot)) +
      geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
      scale_color_gradient(low = "lightgrey", high = "blue") +
      facet_wrap(~ factor(group, levels = names(CD4_markers_list)), scales = "free_y") +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_text(size = 10, face = "bold")
      ) +
      labs(x = "Cluster", y = "Gene", color = "Avg. Expression", size = "Pct. Exp")
    
    print(p)
  }
  
  
  genes_to_check <- c('CD16a','CD4','HNK1','CXCR6','CCR7','SELL','CD69','CXCL13','BAG3','PLCG2','GPR183','STMN1',
                      'CD8A','CD8B','CX3CR1','PDCD1','GNLY','KLRD1','GZMB','GZMK','IFIT3','PLCG2','SLC4A10','FOXP3','NKG7','CD19','CD3D','FCGR3A','CX3CR1','CXCR4','IFNG','PLCG2','CD160','STMN1','CTLA4','LAG3','TRGC1',' TRGC2')
  genes_to_check=str_to_upper(genes_to_check)
  myeloidscRNA<-SetIdent(t_cell_subset,value="RNA_snn_res.1")
  current_idents <- as.numeric(levels(Idents(myeloidscRNA)))
  # 按数字排序
  sorted_idents <- sort(current_idents)
  # 设置排序后的Idents
  myeloidscRNA <- SetIdent(myeloidscRNA, 
                           cells = Cells(myeloidscRNA), 
                           value = factor(Idents(myeloidscRNA), levels = sorted_idents))
  library(stringr)
  table(Idents(t_cell_subset))
  p=DotPlot(myeloidscRNA,features=unique(genes_to_check),
            assay='RNA')+coord_flip()
  p
  
  
  
  new.cluster.ids <- c('CD4T_CCR7','CD4T_CCR7','CD8T_GZMK','CD8T_GZMK','Treg','CD8T_GNLY','CD4T_CCR7','CD4T_CD69','Treg','CD8T_PDCD1','CD8T_GZMK','CD4T_PLCG2','CD8T_STMN1','CD8T_PLCG2')
  length(new.cluster.ids)
  t_cell_subset1<-t_cell_subset
  names(new.cluster.ids) <- levels(t_cell_subset1)
  t_cell_subset1 <- RenameIdents(t_cell_subset1, new.cluster.ids)
  # 将 active.ident 添加到 meta.data 中作为新的列
  t_cell_subset1@meta.data$celltype <- t_cell_subset1@active.ident
  
  table(Idents(t_cell_subset))
  table(Idents(t_cell_subset1))
  
  celltype_colors <- c(
    "CD8T_GZMK" = "#F4A460",     # 金黄色
    "CD4T_CD69" = "#FFA07A",    # 浅珊瑚
    "CD8T_GNLY" = "#FFE4B5",    # 浅小麦色
    "CD8T_PDCD1" = "#FFA500",  # 橙色
    "CD8T_PLCG2" = "mistyrose",     # 浅桃色
    
    "CD8T_STMN1" = "#ADD8E6",     # 浅蓝色
    "CD4T_CCR7" = "#F08080",   # 蒂芙尼蓝
    "CD4T_PLCG2" = "#B22222",     # 钢蓝色
    "CD4T_CXCL13" = "#5F9EA0",   # 灰绿色
    
    "CD8T_CD69" = "#FFDEAD",          # 纳瓦荷白
    "Treg" = "#D8BFD8"        # 爱丽丝蓝
  )
  
  
  DimPlot(t_cell_subset1, reduction = "umap", label = F, pt.size = 0.5)+scale_color_manual(values = celltype_colors)
  genes_to_check <- c('CD4','CCR7','SELL','CD69','CXCL13','BAG3','PLCG2','GPR183','STMN1',
                      'CD8A','CD8B','CX3CR1','PDCD1','GNLY','KLRD1','GZMK','IFIT3','PLCG2','SLC4A10','FOXP3')
  genes_to_check=str_to_upper(genes_to_check)
  myeloidscRNA<-SetIdent(t_cell_subset1,value="celltype")
  current_idents <- as.character(levels(Idents(myeloidscRNA)))
  # 按数字排序
  sorted_idents <- sort(current_idents)
  # 设置排序后的Idents
  myeloidscRNA <- SetIdent(myeloidscRNA, 
                           cells = Cells(myeloidscRNA), 
                           value = factor(Idents(myeloidscRNA), levels = sorted_idents))
  library(stringr)
  p=DotPlot(myeloidscRNA,features=unique(genes_to_check),
            assay='RNA')+coord_flip()
  p
  FeaturePlot(t_cell_subset1,features = 'CD4',label = F)
  FeaturePlot(t_cell_subset1,features = 'CD8A',label = F)
  FeaturePlot(t_cell_subset,features = "RPS19",label = F)
  DimPlot(t_cell_subset,reduction = "umap", label = T, pt.size = 0.5)
  # 比较 cluster 3 和 cluster 7
  markers_3_vs_7 <- FindMarkers(t_cell_subset, ident.1 = 0, ident.2 = 1, only.pos = TRUE)
  
  
}#Tcell
{celltype_colors <- c(
  "CD8T_GZMK" = "#F4A460",     # 金黄色
  "CD4T_CD69" = "#FFA07A",    # 浅珊瑚
  "CD8T_GNLY" = "#FFE4B5",    # 浅小麦色
  "CD8T_PDCD1" = "#FFA500",  # 橙色
  "CD8T_PLCG2" = "mistyrose",     # 浅桃色
  
  "CD8T_STMN1" = "#ADD8E6",     # 浅蓝色
  "CD4T_CCR7" = "#F08080",   # 蒂芙尼蓝
  "CD4T_PLCG2" = "#B22222",     # 钢蓝色
  "CD4T_CXCL13" = "#5F9EA0",   # 灰绿色
  
  "CD8T_CD69" = "#FFDEAD",          # 纳瓦荷白
  "Treg" = "#D8BFD8"        # 爱丽丝蓝
)
  
  table(Idents(t_cell_subset1))
  p<-DimPlot(t_cell_subset1, reduction = "umap", label = F, pt.size = 0.5)+
    scale_color_manual(values = celltype_colors)
  p
  pdf("/home/zbh1/ZBH-SC/UMAPTcell.pdf",family = 'sans', width = 5.4, height = 3.61)
  print(p)
  dev.off()
  genes_to_check <- c('CD4','CD8A','CD8B','FOXP3','PLCG2','GZMK','GNLY','STMN1',
                      'PDCD1','CD69','CCR7')
  p<-DotPlot(object = t_cell_subset1, features = genes_to_check)
  
  cell_order <- c("CD4T_CCR7", "CD4T_CD69", "CD4T_PLCG2", 
                  "CD8T_PDCD1", "CD8T_STMN1", "CD8T_GNLY", "CD8T_GZMK", "CD8T_PLCG2",'Treg')
  # 重新排序
  p$data$id <- factor(p$data$id, levels = cell_order)
  
  p1 <- ggplot(p$data, aes(x = features.plot, y = id)) + 
    geom_point(aes(size = pct.exp, color = avg.exp.scaled))+
    scale_size(range = c(4, 6)) + 
    theme_classic() + 
    scale_color_gradient2(low = "#BDBADB",  mid = "white", high = "#FBB463") +  # 暖色调渐变
    theme(
      axis.text.x = element_text(angle = 90, face = 1, size = 10, 
                                 family = "sans", hjust = 1, vjust = 0.5, color = "black"),
      axis.text.y = element_text(size = 10, face = 1, 
                                 family = "sans", color = "black"),
      legend.text = element_text(size = 8, face = 1, family = "sans"),
      legend.title = element_text(size = 10, face = 1, family = "sans"),
      legend.position = 'top',
      strip.placement = "outside",
      strip.text.x = element_text(size = 9, family = "sans"),
      axis.title = element_blank()
    ) +
    guides(colour = guide_colourbar(title.vjust = 0.9, title.hjust = 0)) +
    labs(size = "Percent Expressed", color = "Average Expression") +
    ggplot2::ggtitle("DotPlot of Marker Genes Across Cell Types")
  
  p1
  pdf("/home/zbh1/ZBH-SC/DotplotTcell.pdf",family = 'sans', width = 7.07, height = 4.09)
  print(p1)
  dev.off()
}#T fig