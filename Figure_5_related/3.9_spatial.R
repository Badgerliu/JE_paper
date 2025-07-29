setwd("../../../zhanglab/lh/3.9_spatial_GEO")
setwd("../3.9_spatial_GEO")
setwd("../")
test <- Read10X_Image("B2/spatial/",image.name = "image.png")

sptaialST <- test@image

count <- Read10X("B2")
seurat <- CreateSeuratObject(count)
seurat@images$spatial <- test

SpatialDimPlot(seurat,images = "spatial",pt.size.factor=1,image.alpha=0.4)
saveRDS(seurat,"processed_data/B2_seurat.Rds")

seurat <- SCTransform(seurat,assay = "RNA", verbose = FALSE)
seurat@assays$RNA@key <- "rna_"
FetchData(seurat, vars=c("GZMH","GZMB"))
SpatialFeaturePlot(seurat, features = c("KRT10", "KRT5"),pt.size.factor=10)
seurat <- RunPCA(seurat, assay = "SCT", verbose = FALSE)
seurat <- FindNeighbors(seurat, reduction = "pca", dims = 1:30)
seurat <- FindClusters(seurat, verbose = FALSE,resolution = 1.6)
seurat <- RunUMAP(seurat, reduction = "pca", dims = 1:30)
SpatialDimPlot(seurat, label = TRUE, label.size = 3,pt.size.factor=10)


seurat <- FindSpatiallyVariableFeatures(seurat, assay = "SCT", features = VariableFeatures(seurat)[1:1000],
                                       selection.method = "moransi")
top.features <- head(SpatiallyVariableFeatures(seurat, selection.method = "moransi"), 50)

de_markers2 <- FindAllMarkers(seurat )
top.markers <- de_markers2%>% group_by(cluster) %>% top_n(n = 20,wt = avg_log2FC)

SpatialFeaturePlot(seurat, features =head(top.markers$gene),pt.size.factor=10)


markerGeneName <- markerGene$Gene

load("../../database/gene-related/geneinfo_2022.rda")
geneInfo <- geneinfo_2022%>%
  dplyr::select(symbol, symbol_mouse)
mapper = function(df, value_col, name_col) setNames(df[value_col]%>%unlist, df[name_col]%>%unlist)
humansymbol2mousesymbol = mapper(geneinfo_2022,"symbol","symbol_mouse")

convertFun <- function(symbols){
  converted_symbols = symbols %>% as.character %>%humansymbol2mousesymbol[.]
  return(converted_symbols)
}

humanGeneMouse <- convertFun(markerGeneName)

SpatialFeaturePlot(seurat, features =humanGeneMouse,pt.size.factor=10,ncol = 5,image.alpha = 0)
ggsave("result/markerGene.pdf",width = 25,height = 20)


varGene <- VariableFeatures(gingival_v2_d0)
varGeneHuman <- convertFun(varGene)
intersectGene <- intersect(varGeneHuman,top.markers$gene)
SpatialFeaturePlot(seurat, features =intersectGene[1:40],pt.size.factor=10,ncol = 5,image.alpha = 0)
ggsave("result/spatialGene_1.pdf",width = 25,height = 20)
SpatialFeaturePlot(seurat, features =intersectGene[41:79],pt.size.factor=10,ncol = 5,image.alpha = 0)
ggsave("result/spatialGene_2.pdf",width = 25,height = 20)
# 
# 
# 
# 
# SpatialFeaturePlot(brain, features = top.features, ncol = 3, alpha = c(0.1, 1))
# 
# #DimPlot(seurat, reduction = "umap", label = TRUE)
# seurat@images$spatial@key <- "rna_"
# 
# RenameCells(seurat,new.names = paste0(Cells(x = seurat), 
#                                       "_", "reference"))




DefaultAssay(gingival_v2_d0) <- "SCT"

anchors <- FindTransferAnchors(reference = seurat, query = gingival_v2_d0, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass, prediction.assay = TRUE,
                                  weight.reduction = cortex[["pca"]], dims = 1:30)