# cytpscape file prepare----------------------------------
setwd("../3.16_cytoscpace/")
library(tidyverse)

seurat2 <- readRDS("../3.9_spatial_GEO/processed_data/B2_seurat.Rds")
seurat <- readRDS("../3.9_spatial_GEO/processed_data/B2_seurat.Rds")
gingival_v2_d0Mat <- as.matrix(gingival_v2_d0@assays$RNA@counts)
gingival_v2_d0Mat <- as.data.frame(gingival_v2_d0Mat)
gingival_v2_d0Mat <- gingival_v2_d0Mat%>%rownames_to_column("GENES")
data.table::fwrite(gingival_v2_d0Mat,"input/scRNA_gingivald0v2.txt",sep = "\t",row.names = F)
unique(gingival_v2_d0@meta.data$celltype)
ctFile <- as.data.frame(gingival_v2_d0@meta.data$celltype)
rownames(ctFile) <- colnames(gingival_v2_d0)
ctFile <- ctFile%>%rownames_to_column("Cell IDs")
colnames(ctFile) <- c("Cell IDs","CellType")
write.table(ctFile,"input/celltype.txt",row.names = F,sep = "\t")
seurat@images$spatial@assay <- "RNA"
seurat@images$spatial@key <- "rna_"
coordinates <- coordinates[colnames(seurat),]

#colData(sce) = cbind(colData(sce), coordinates)
cells=colnames(seurat)[coordinates$tissue!=0]
stSub <- subset(seurat,cells=rownames(seurat@images$spatial@coordinates))

stFile <- as.matrix(stSub@assays$RNA@counts)
mouseGene <- convertFunH2M(rownames(stFile))
stFile <- stFile[!is.na(mouseGene),]
mouseGene <- na.omit(mouseGene)
rownames(stFile) <- mouseGene
stFile <- as.data.frame(stFile)
stFile <- stFile%>%rownames_to_column("GENES")
data.table::fwrite(stFile,"input/stMat_2.txt",sep = "\t")
dim(coordinates)
#coordFile <- coordinates[c(2,3)]
coordFile <- seurat@images$spatial@coordinates[c(2,3)]
coordFile <- coordFile%>%rownames_to_column("SpotID")
#coordinates <- coordinates%>%rownames_to_column("SpotID")

stFileByes <- as.data.frame(count_sub)
mouseGene <- convertFunH2M(rownames(stFileByes))
stFileByes <- stFileByes[!is.na(mouseGene),]
mouseGene <- na.omit(mouseGene)
#remove duplicate
subsetGene <- names(table(mouseGene))[(table(mouseGene)<2)]
logicSubset=mouseGene%in%subsetGene
stFileByes <- stFileByes[logicSubset,]
mouseGene <- mouseGene[logicSubset]
rownames(stFileByes) <- mouseGene

stFileByes <- stFileByes%>%rownames_to_column("GENES")
data.table::fwrite(stFileByes,"input/stMat_bayes.txt",sep = "\t")

sceCol <- sce_sub@colData
coord_bayes <- sceCol[c(5,6)]%>%as.data.frame()
coord_bayes <- coord_bayes%>%rownames_to_column("SpotID")
data.table::fwrite(coord_bayes,"input/stCoord_bayes.txt",sep = "\t")
coord_bayes <- coord_bayes[c(5,6)]


data.table::fwrite(coordFile,"input/stCoord_2.txt",sep = "\t")
data.table::fwrite(coordinates,"processed_data/fullCoord.txt",sep = "\t")
library(SeuratDisk)
SaveH5Seurat(seurat,"processed_data/spatial_GSE.h5Seurat",overwrite = T)
Convert("processed_data/spatial_GSE.h5Seurat",dest="h5ad",,overwrite = T)

DefaultAssay(seurat) <- "SCT"
mouseGene <- convertFunH2M(rownames(seurat@assays$SCT))
seurat_mouse<- seurat[!is.na(mouseGene),]
mouseGene <- na.omit(mouseGene)
mouseGene <- as.character(mouseGene)
RenameGenesSeurat <- function(obj, newnames,assay) { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  RNA <- obj@assays[[assay]]
  
  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
    #if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}
seurat_mouse <- RenameGenesSeurat(seurat_mouse,mouseGene,"SCT")
mouseGeneScale <- convertFunH2M(rownames(seurat@assays$SCT@scale.data))
mouseGeneScale <- na.omit(mouseGeneScale)
mouseGeneScale <- as.character(mouseGeneScale)
seurat_mouse_bk <- seurat_mouse
seurat_mouse <- seurat_mouse_bk

#seurat_mouse <- DietSeurat(seurat_mouse)
seurat_mouse@assays$SCT@counts@Dimnames[[1]] <- mouseGene
seurat_mouse@assays$SCT@data@Dimnames[[1]] <- mouseGene
rownames(seurat_mouse@assays$SCT@scale.data) <- mouseGeneScale
rownames(seurat_mouse@assays$RNA@scale.data) <- mouseGeneScale
rownames(seurat_mouse@assays$SCT@scale.data) <- mouseGeneScale
DefaultAssay(seurat_mouse) <- "RNA"
seurat_mouse_sce <- as.SingleCellExperiment(seurat_mouse)
zellkonverter::writeH5AD(seurat_mouse_sce,"processed_data/spatial_GSE2.h5ad")
SaveH5Seurat(seurat_mouse,"processed_data/spatial_GSE.h5Seurat",overwrite = T)
Convert("processed_data/spatial_GSE.h5Seurat",dest="h5ad",,overwrite = T)
