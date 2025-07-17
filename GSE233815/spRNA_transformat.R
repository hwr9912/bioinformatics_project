rm(list = ls())
library(Seurat)
# 设置工作路径
setwd("/Volumes/MacPassport/project/bioinfo/GSE233815")
getwd()

# 读取 .rds 文件
spatial_mcao <- readRDS("data/spatial/seurat_1stSpatial.rds")
spatial_sham <- readRDS("data/spatial/seurat_2ndSpatial.rds")
spatial_all <- list(spatial_mcao[[1]], spatial_mcao[[2]], spatial_mcao[[3]],
                    spatial_sham[[1]], spatial_sham[[2]])

seurat2hdf <- function(seurat_obj, name="data/"){
  # 1. 提取表达数据
  expr_data <- as.matrix(GetAssayData(seurat_obj, slot = "counts")) # 原始counts
  
  # 2. 提取细胞元数据 (cell metadata)
  cell_meta <- seurat_obj@meta.data
  
  # 3. 提取基因元数据 (feature metadata)
  # 在 Seurat v5 中，特征元数据通常直接在 Assay 的 features 槽位中
  # 否则，可能在 seurat_obj@assays$RNA@meta.features
  gene_meta <- GetAssay(seurat_obj, assay = "Spatial")@meta.features
  
  # 4. 提取降维嵌入 (例如 PCA 或 UMAP 坐标)
  # 注意：确保这些降维结果已经计算过了
  # pca_coords <- Embeddings(object = seurat_obj, reduction = "pca")
  # umap_coords <- Embeddings(object = seurat_obj, reduction = "umap")
  
  # 5. 提取聚类信息 (如果已经计算了聚类)
  # 例如，默认的 Seurat 聚类结果通常在 meta.data 中，列名为 'seurat_clusters'
  clusters <- seurat_obj@meta.data$seurat_clusters
  
  # --- 保存这些数据到中间文件 ---
  # 可以保存为 CSV 文件，Python 很容易读取
  write.csv(expr_data, file = paste0(name, "expression_data.csv"), row.names = TRUE)
  write.csv(cell_meta, file = paste0(name, "cell_metadata.csv"), row.names = TRUE)
  write.csv(gene_meta, file = paste0(name, "gene_metadata.csv"), row.names = TRUE) # 如果有的话
  
  # 对于降维和聚类，可以单独保存
  # write.csv(pca_coords, file = paste0("pca_coords.csv"), row.names = TRUE)
  # write.csv(umap_coords, file = paste0("umap_coords.csv"), row.names = TRUE)
  # write.csv(data.frame(clusters = clusters), file = paste0("clusters.csv"), row.names = TRUE)
}

fname <- c("d1","d3","d7","sham","d7b")
for (i in seq_along(fname)) {
  seurat2hdf(seurat_obj = spatial_all[[i]], name = paste0("data/spatial/GSE233813_csv/", fname[[i]], "_"))
}
