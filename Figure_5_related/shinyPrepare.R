
scMergeSeurat <- FindVariableFeatures(scMergeSeurat,nfeatures = 4000)
selectGene <- c("Krt5","Trp63")

scMergeSeuratDiet <- scMergeSeurat[c(VariableFeatures(scMergeSeurat),selectGene),]
# Sample 1000 unique values
#unique_values <- coord_merge %>% distinct()

index <- coord_merge%>%
  rownames_to_column("rowname")%>%
  group_by(spot) %>%
  sample_n(1) %>%
  ungroup()
subsetIndex <-  rownames(coord_merge)%in%index$rowname 

scMergeSeuratDiet <- scMergeSeuratDiet[,subsetIndex]
library(Seurat)
SpatialFeaturePlot(scMergeSeuratDiet,"Krt5",pt.size.factor=5,image.alpha=0.4)&scale_fill_viridis_c(option ="D")
saveRDS(scMergeSeuratDiet,"shinyapp/data/spatialSeurat.Rds")


write.csv(markers,"important_processed_data/4.15_marker_spatial.csv")

# label stability------------
label_stability <- read.table("4.5_cytospace/processed_data/4.5_cytospace_results/Seurat_cellfracs.txt",row.names = 1)
label_stability_val <- apply(label_stability, 1, max, na.rm=TRUE)
label_stability2 <- read.table("4.5_cytospace/processed_data/v2_cytoscape/Seurat_cellfracs.txt",row.names = 1)
label_stability_val2 <- apply(label_stability2, 1, max, na.rm=TRUE)

scMergeSeuratDiet$stability1 <- label_stability_val[scMergeSeuratDiet$point]
scMergeSeuratDiet$stability2 <- label_stability_val2[scMergeSeuratDiet$point]
p1 <- SpatialFeaturePlot(scMergeSeuratDiet,"stability1",pt.size.factor=5,image.alpha=0.4)&scale_fill_viridis_c(option ="D")
p2 <- SpatialFeaturePlot(scMergeSeuratDiet,"stability2",pt.size.factor=5,image.alpha=0.4)&scale_fill_viridis_c(option ="D")
p1+p2
