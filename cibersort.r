
# 定义一个函数，检查并安装包
install_if_missing <- function(packages) {
  for (package in packages) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package)
    }
  }
}

# 安装所有必要的包
install_if_missing(c("devtools", "CIBERSORT", "ggplot2", "dplyr", "ggthemes", "pheatmap", "tibble", "tidyr", "ggpubr", "ggsci"))

# 加载这些包
library(devtools)
devtools::install_github("Moonerss/CIBERSORT")
library(CIBERSORT)
library(ggplot2)
library(dplyr)
library(ggthemes)
library(pheatmap)
library(tibble)
library(tidyr)
library(ggpubr)
library(ggsci)



remove(list = ls()) #一键清空
#读取LM22文件（免疫细胞特征基因文件）
# 强制读取数据

sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
mixture_file <- system.file("extdata", "exampleForLUAD.txt", package = "CIBERSORT")
data(LM22)
data(mixed_expr)
data <- read.table(file="D:/单细胞测序/Seq2_3_Comparison/summary (2)/summary/1_Cellranger_result/A0_EV/A0_EV.expression.txt/expression_addgeneid.txt",header=F, sep="\t",check.names=F,quote="")
data_filter = data
# 将数据框转换为矩阵
data_matrix <- as.matrix(data)
head(mixed_expr[1:3],n=3)
head(data_matrix[1:3],n=3)
head(data_filter[1:3],n=3)
str(data_filter)
dim(data)
dim(data_filter)
dataresults <- cibersort(sig_matrix = LM22, mixture_file = data_filter,perm = 1000,QN = F)
# perm置换次数=1000
# QN如果是芯片设置为T，如果是测序就设置为F
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggsci)
res <- data.frame(dataresults[,1:22])%>%
  mutate(group = c(rep('shControl',4),rep('shPHF8',8)))%>%
  rownames_to_column("sample")%>%
  pivot_longer(cols = colnames(.)[2:23],
               names_to = "cell.type",
               values_to = 'value')
head(res,n=6)
## # A tibble: 6 × 4
##   sample      group     cell.type                   value
##   <chr>       <chr>     <chr>                       <dbl>
## 1 shControl.1 shControl B.cells.naive              0     
## 2 shControl.1 shControl B.cells.memory             0.0404
## 3 shControl.1 shControl Plasma.cells               0.0163
## 4 shControl.1 shControl T.cells.CD8                0.129 
## 5 shControl.1 shControl T.cells.CD4.naive          0     
## 6 shControl.1 shControl T.cells.CD4.memory.resting 0

ggplot(res,aes(cell.type,value,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 14,face = "italic",colour = 'black'),
        axis.text.y = element_text(face = "italic",size = 14,colour = 'black'))+
  scale_fill_nejm()+
  stat_compare_means(aes(group = group),label = "p.format",size=3,method = "kruskal.test")
