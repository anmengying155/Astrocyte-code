import pandas as pd

# 读取 Excel 文件
df = pd.read_excel('D:\File\document\单细胞测序代码\AUcell\参考功能基因集.xlsx')  # 替换为你的 Excel 文件路径

# 创建一个空的字符串列表，用于存储 .gmt 格式的每一行
gmt_lines = []

# 遍历每一行，格式化为 .gmt 格式
for _, row in df.iterrows():
    gene_set_name = row['Functional gene sets']  # 假设列名是 'Gene Set Name'
    genes = row['Genes']  # 假设列名是 'Genes'
    
    # 将基因集名称和基因按制表符分隔
    gmt_line = f"{gene_set_name}\t{gene_set_name}\t" + "\t".join(genes.split(','))
    gmt_lines.append(gmt_line)

# 将结果保存为 .gmt 文件
with open('referencegeneset.gmt', 'w') as f:
    f.write("\n".join(gmt_lines))

print("转换完成，.gmt 文件已保存为 'referencegeneset.gmt'")
