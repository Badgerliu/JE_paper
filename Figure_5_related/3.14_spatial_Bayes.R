library(BayesSpace)
library(Seurat)
seurat_ca[]
sce = as.SingleCellExperiment(seurat)
coordinates <- read.csv("../3.9_spatial_GEO/B2/spatial/tissue_positions_list.csv",header = F)
rownames(coordinates) <- NULL
coordinates <- coordinates%>%column_to_rownames("V1")
#rownames(coordinates) <- coordinates$V1
colnames(coordinates) <- colnames(test@coordinates)
coordinates <- coordinates[colnames(sce),]

colData(sce) = cbind(colData(sce), coordinates)
sce=sce[,coordinates$tissue!=0]

sce = spatialPreprocess(sce, platform = "Visium", skip.PCA = T, log.normalize = F)
sce = spatialCluster(sce, nrep = 1000, burn.in = 100, q = 10) #quickly cluster via BayesSpace
clusterPlot(sce) #plot via BayesSpace

sce_orig = as.SingleCellExperiment(seurat)

colData(sce_orig) = cbind(colData(sce_orig), coordinates)
sce_orig=sce_orig[,coordinates$tissue!=0]
sce_orig = spatialPreprocess(sce_orig, platform = "Visium", skip.PCA = T, log.normalize = F)
sce_orig = spatialCluster(sce_orig, nrep = 1000, burn.in = 100, q = 10) #quickly cluster via BayesSpace

sce <- spatialEnhance(sce_orig, q=10, platform="Visium", d=7,
                                    model="t", gamma=3,
                                    jitter_prior=0.3, jitter_scale=3.5,
                                    nrep=1000, burn.in=100,
                                    save.chain=TRUE)

clusterPlot(sce)

susetIdent <- c("1","2","9","5","6","10")



markers <- c("KRT5", "KRT10", "HMGB2")
marker2 <- rownames(sce)[11:20]
sce <- enhanceFeatures(sce, sce_orig,model = "lm",
                                     feature_names=rownames(sce),
                                     nrounds=100)
sce <- enhanceFeatures(sce, sce_orig,model = "lm",
                       feature_names=rownames(sce),assay.type = "counts",
                       nrounds=100)
sce_orig@assays@data$counts[1:10,1:10]

count_pred <-  sce@assays@data$counts
logcount_pred <- sce@assays@data$logcounts



sce_sub <- sce[,sce@colData$spatial.cluster%in%susetIdent]
saveRDS(sce_sub,"processed_data/spatial_bayes_subset.Rds")

count_sub <- sce_sub@assays@data$counts



featurePlot(sce, c("HMGB2","KRT10","KRT5"))
featurePlot(sce, c("KRT10"))
featurePlot(sce, c("KRT14"))
#test <- Read10X_Image("B2/spatial/",image.name = "image.png")

clusterPlot(sce_sub)
ggsave("result/bayes_spatil.pdf")
