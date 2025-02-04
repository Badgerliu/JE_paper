# 2024.4.12-----------------------------
# read h5ad file -------
library(ComplexHeatmap)
library(clusterProfiler)
st <- readRDS("important_processed_data/st_seurat_v3.Rds")
scMerge <- zellkonverter::readH5AD("important_processed_data/4.7_scMerge.h5ad")
#debug(as.Seurat)
#undebug(as.Seurat)
scMergeSeurat <- as.Seurat(scMerge,counts = "X",data = "logcounts")

scMergeSeurat <- NormalizeData(scMergeSeurat, normalization.method = "LogNormalize", scale.factor = 10000)


#split(colnames(scMergeSeurat))
#scMergeSeurat@images
scCol = strsplit(colnames(scMergeSeurat),"_")%>%lapply(`[` ,1) %>%unlist
coordMerge <- st@images$spatial@coordinates[scCol,]
newImage <- st@images
newImage$spatial@coordinates <- coordMerge
rownames(newImage$spatial@coordinates ) <- colnames(scMergeSeurat)

scMergeSeurat@images <- newImage
scMergeSeurat@images$spatial@assay <- "originalexp"
scMergeSeurat@images$spatial@key <- "originalexp_"
#debug(SpatialDimPlot)
#undebug(SpatialDimPlot)

#Images(scMergeSeurat)
#debug(Images)
coord_merge <- unique(coordst[,c("tissue", "row", "col", "imagerow", "imagecol", "spot"),])
rownames(coord_merge) <- coord_merge$spot
coord_merge <-coord_merge[as.character(scMergeSeurat$point),]
#GetTissueCoordinates
coord_merge
#coord_merge <- coordst%>%column_to_rownames("spot")
newImage$spatial@coordinates[c("tissue", "row","col","imagerow","imagecol")] <- coord_merge[c("tissue", "row", "col", "imagerow", "imagecol")]
scMergeSeurat@images <- newImage
Images(scMergeSeurat,"originalexp")

#scMergeSeurat@assays$originalexp@data
scMergeSeurat@images$spatial@coordinates[c("row","col")] <- scMergeSeurat@reductions$spatial_bk@cell.embeddings
SpatialDimPlot(scMergeSeurat)
coordst1 <- scMergeSeurat@images$spatial@coordinates
newImage <- st@images
imgCoor <- newImage$spatial@coordinates
lm1 <- lm(imagecol~row, data=imgCoor)
imgColPred <- predict(lm1,coordst1)
coordst1$imagecol <- imgColPred
lm2 <- lm(imagerow~col, data=imgCoor)
imgrowPred <- predict(lm2,coordst1)
coordst1$imagerow <- imgrowPred

scMergeSeurat@images$spatial@coordinates <- coordst1
scMergeSeurat@images$spatial@image <- st@images$spatial@image 
st <- readRDS("important_processed_data/st_seurat_v3.Rds")
st <- rotateSeuratImage(st,slide="spatial",rotation="L90")
SpatialFeaturePlot(st,c("Krt5"),pt.size.factor=5,image.alpha=0.4,ncol = 4)&scale_fill_viridis_c(option ="D")
scMergeSeurat@images$spatial@image <- st@images$spatial@image 

scMergeSeurat <- rotateSeuratImage(scMergeSeurat,slide="spatial",rotation="Hf")
SpatialFeaturePlot(scMergeSeurat,c("Krt5"),pt.size.factor=5,image.alpha=0.4,ncol = 4)&scale_fill_viridis_c(option ="D")

scMergeSeurat <- rotateSeuratImage(scMergeSeurat,slide="spatial",rotation="R90")

saveRDS(scMergeSeurat,"important_processed_data/6.18_merge_3day_6sample.Rds")

SpatialFeaturePlot(scMergeSeurat,c("pattern0", "pattern1", "pattern2", "pattern3", "pattern4", 
                       "pattern5", "pattern6", "pattern7"),pt.size.factor=5,image.alpha=0.4,ncol = 4)&scale_fill_viridis_c(option ="D")


scMergeSeurat@images$spatial@image <- gingivalST@images$spatial@image
scMergeSeurat <- rotateSeuratImage(scMergeSeurat,slide="spatial",rotation="R90")
#scMergeSeurat2@images$spatial@scale.factors$lowres <- 0.58
SpatialFeaturePlot(scMergeSeurat,c("Krt5"),pt.size.factor=5,image.alpha=0.4,ncol = 4)&scale_fill_viridis_c(option ="D")

SpatialFeaturePlot(scMergeSeurat,c("pattern0", "pattern1", "pattern2", "pattern3", "pattern4", 
                                   "pattern5", "pattern6", "pattern7"),pt.size.factor=5,image.alpha=0.4,ncol = 4)&scale_fill_viridis_c(option ="D")

ggsave("results/FigureD/pattern_all.pdf",width = 10,height = 8)
ggsave("results/FigureD/pattern_all.png",width = 10,height = 8)

SpatialFeaturePlot(scMergeSeurat,c("Dsg1a","Clip4","Klf3","Dlx3"),pt.size.factor=5,image.alpha=0.4,ncol = 4)&scale_fill_viridis_c(option ="D")

ggsave("results/FigureD/coor_feature.pdf",width = 10,height = 4)
ggsave("results/FigureD/coor_feature.png",width = 10,height = 4)
SpatialFeaturePlot(scMergeSeurat,c("Fam25c"),pt.size.factor=5,image.alpha=0.4,ncol = 4)&scale_fill_viridis_c(option ="D")


SpatialFeaturePlot(scMergeSeurat,c("pattern1"),pt.size.factor=5,image.alpha=0.4,ncol = 1)&scale_fill_viridis_c(option ="D")
ggsave("results/FigureD/pattern1.pdf",width = 4,height = 4)
ggsave("results/FigureD/pattern1.png",width = 4,height = 4)

patternCorrDf <- read.csv("results/patternDE/corr_mean_pattern1.csv",row.names = 1)
patternCorrDf <- sort(patternCorrDf,decreasing = TRUE)
ComplexHeatmap::Heatmap(patternCorrDf,,show_row_names = F,cluster_rows = F)

genesSelect <- c("Dsg1a","Clip4","Klf3","Dlx3")
index = which(rownames(patternCorrDf) %in% genesSelect)
anno = anno_mark(at = index, labels = genesSelect, which = "column",side="bottom")
#Heatmap(patternCorrDf, show_row_names = F,cluster_rows = F) + rowAnnotation(mark = anno)
col_fun = colorRamp2(range(patternCorrDf), hcl_palette = "Blues 3", reverse = TRUE)
pdf("results/FigureD/CorrHM.pdf",width = 10,height = 1.3)
Heatmap(t(patternCorrDf), show_column_names  = F,show_row_names = F, col = col_fun,
        cluster_columns = F,bottom_annotation = columnAnnotation(mark = anno)) 
dev.off()



#===Figure D dotplot----------------------------
filterGene <- read.csv("4.6_spatial_pattern/mannual_pattern_filter.csv")
filterGeneList <- split(filterGene$Gene,filterGene$Class)
goCompare <- compareCluster(geneCluster = filterGeneList, fun = enrichGO,OrgDb = "org.Mm.eg.db",ont="ALL",keyType="SYMBOL")
dotplot(goCompare)
write.csv(goCompare@compareClusterResult,"results/patternDE/GO_result.csv")

cnetplot(goCompare)
ggsave("results/patternDE/GO_cnetplot.pdf",width = 8,height = 8)
genes <- goCompare@compareClusterResult

GOselect <- read.table("processed_data/GO_select")%>%unlist

clusterRes <- goCompare@compareClusterResult
resultFull <- clusterRes%>%filter(ID%in%GOselect)
goCompare@compareClusterResult <- resultFull
dotplot(goCompare,showCategory=26)
ggsave("results/FigureD/GO_full.pdf",width = 7,height = 8)

GOselectFil <- read.table("processed_data/GO_select_main")%>%unlist
resultMain <- clusterRes%>%filter(ID%in%GOselectFil)
goCompare@compareClusterResult <- resultMain
dotplot(goCompare,showCategory=15)
ggsave("results/FigureD/GO_main.pdf",width = 6,height = 5)


goplot(goCompare)

# D0GO <- enrichGO(filterGeneList$D0,OrgDb = "org.Mm.eg.db",ont="ALL",keyType="SYMBOL")
# D3GO <- enrichGO(filterGeneList$D3,OrgDb = "org.Mm.eg.db",ont="ALL",keyType="SYMBOL")
# D5GO <- enrichGO(filterGeneList$D5,OrgDb = "org.Mm.eg.db",ont="ALL",keyType="SYMBOL")
# goplot(D0GO)

#=== Featureplot
scMergeSeurat$day <- NA
scMergeSeurat$day[scMergeSeurat$orig.ident%in%c("D0","D0_K5_GFP")] ="D0"
scMergeSeurat$day[scMergeSeurat$orig.ident%in%c("D3","D3_K5_GFP")] ="D3"
scMergeSeurat$day[scMergeSeurat$orig.ident%in%c("D5","D5_K5_GFP")] ="D5"
d0 <- scMergeSeurat[,scMergeSeurat$day=="D0"]
d3 <- scMergeSeurat[,scMergeSeurat$day=="D3"]
d5 <- scMergeSeurat[,scMergeSeurat$day=="D5"]
p1 <- SpatialFeaturePlot(d0,c("S100a9"),pt.size.factor=5,image.alpha=0.4)&scale_fill_viridis_c(option ="D")
p2 <-  SpatialFeaturePlot(d3,c("S100a9"),pt.size.factor=5,image.alpha=0.4)&scale_fill_viridis_c(option ="D")
p3 <-  SpatialFeaturePlot(d5,c("S100a9"),pt.size.factor=5,image.alpha=0.4)&scale_fill_viridis_c(option ="D")
p = p1|p2|p3
p
ggsave("results/FigureD/S100a9.pdf",width = 9,height = 3)
ggsave("results/FigureD/S100a9.png",width = 9,height = 3)
