gingival_sc_1 <- readRDS("../../liuhuan/disk2/JE_epi_communication/inputdata/JE_combined_epi_harmony_MAGIC_renamed.rds")
gingival_sc_2 <- readRDS("../../liuhuan/disk2/JE_epi_communication/inputdata/JE_combined_epi_celltype_with_epi_20240305.rds")


FeaturePlot(gingival_sc_1,"Odam",split.by = "orig.ident")

FeaturePlot(gingival_sc_2,"Odam",split.by = "orig.ident")

unique(gingival_sc_1$orig.ident)
unique(gingival_sc_2$orig.ident)


gingival_count <- as.matrix(gingival_sc_1@assays$RNA@counts)
gingival_count <- as.data.frame(gingival_count)
gingival_count <- gingival_count%>%rownames_to_column("GENES")
data.table::fwrite(gingival_count,"4.5_cytospace/input/gingival_v1.txt",sep = "\t",row.names = F)


ctFile <- as.data.frame(gingival_sc_1@meta.data$celltype)
rownames(ctFile) <- colnames(gingival_sc_1)
ctFile <- ctFile%>%rownames_to_column("Cell IDs")
colnames(ctFile) <- c("Cell IDs","CellType")
write.table(ctFile,"4.5_cytospace/input/celltype.txt",row.names = F,sep = "\t")


#== run results processed by cytospace----------------------------------
location <- read.csv("../4.5_cytospace/processed_data/4.5_cytospace_results/assigned_locations.csv")


RNA <- gingival_sc_1@assays$RNA@counts
RNA_st <- RNA[,location$OriginalCID]
colnames(RNA_st) <- location$UniqueCID
gingivalST@assays$RNA <- CreateAssayObject(RNA_st)
gingivalST <- CreateSeuratObject(RNA_st)

metaST <- gingival_sc_1@meta.data[location$OriginalCID,]
rownames(metaST) <- location$UniqueCID
gingivalST@meta.data <- metaST

coordst1 <- location[c(5,6)]

locationID <- coordst[,c("spot","idents")]
spotID <- locationID$idents
names(spotID) <- locationID$spot

location$idents=spotID[location$SpotID]

gingivalST$idents <- location$idents

SpatialDimPlot(gingivalST,group.by = "idents",pt.size.factor = 5)
ggsave("results/dimplot_st.pdf")
SpatialDimPlot(gingivalST,group.by = "orig.ident",pt.size.factor = 5)
mycolor<-colorRampPalette(brewer.pal(8,'Spectral'))(length(unique(metaBar$idents)))


metaBar <- gingivalST@meta.data[c("orig.ident","idents")]
ggplot(metaBar, aes(x=orig.ident, fill=idents)) + 
  geom_bar(position = "fill")+theme_classic()+scale_fill_manual(values = mycolor)+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1,size = 10,face = "bold"),title = NULL)
ggsave("results/barplot_cluster.pdf")
ggplot(metaBar, aes(x=idents, fill=orig.ident)) + 
  geom_bar(position = "fill")+theme_classic()+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1,size = 10,face = "bold"),title = NULL)
ggsave("results/barplot_ident.pdf")
saveRDS(gingivalST,"../important_processed_data/4.5_v3_gingival_1st.Rds")
gingivalST <- readRDS("../important_processed_data/4.5_v3_gingival_1st.Rds")

STv3 <- readRDS("important_processed_data/st_seurat_v3.Rds")

img <- STv3@images$spatial
imgCoor <- img@coordinates
lm1 <- lm(imagecol~row, data=imgCoor)
imgColPred <- predict(lm1,coordst1)
coordst1$imagecol <- imgColPred
lm2 <- lm(imagerow~col, data=imgCoor)
imgrowPred <- predict(lm2,coordst1)
coordst1$imagerow <- imgrowPred
coordst1$tissue=1
coordst1 <- coordst1[colnames(imgCoor)]
rownames(coordst1) <- location$UniqueCID
img@coordinates <- coordst1
gingivalST@images$spatial <- img
Idents(gingivalST) <- gingivalST$celltype
SpatialDimPlot(gingivalST,alpha = 1,pt.size.factor = 5)
ggsave("results/cluster_dimplot.pdf")

SpatialDimPlot(gingivalST,group.by = "orig.ident",alpha = 1,pt.size.factor = 5)
ggsave("results/cluster_dimplot_orig.ident.pdf")


#==  prepare-----------------
st <- data.table::fread("input/stMat_bayes.txt")
geneInterset <- intersect(st$GENES,gingival_count$GENES)

gingival_sc_2 <- gingival_sc_2[,gingival_sc_2$orig.ident%in%c("D0_K5_GFP", "D3_K5_GFP", "D5_K5_GFP")]
gingival_sc_2 <- gingival_sc_2[geneInterset,]

gingival_count <- as.matrix(gingival_sc_2@assays$RNA@counts)
gingival_count <- as.data.frame(gingival_count)
gingival_count <- gingival_count%>%rownames_to_column("GENES")
data.table::fwrite(gingival_count,"input/gingival_v2.txt",sep = "\t",row.names = F)


ctFile <- as.data.frame(gingival_sc_2@meta.data$celltype)
rownames(ctFile) <- colnames(gingival_sc_2)
ctFile <- ctFile%>%rownames_to_column("Cell IDs")
colnames(ctFile) <- c("Cell IDs","CellType")
write.table(ctFile,"input/celltype_v2.txt",row.names = F,sep = "\t")

#== data----------------------

location <- read.csv("processed_data/v2_cytoscape/assigned_locations.csv")


RNA <- gingival_sc_2@assays$RNA
RNA_st <- RNA@counts[,location$OriginalCID]
colnames(RNA_st) <- location$UniqueCID
#gingivalST2@assays$RNA <-  CreateAssayObject(RNA_st)
gingivalST2 <- CreateSeuratObject(RNA_st)

metaST <- gingival_sc_2@meta.data[location$OriginalCID,]
rownames(metaST) <- location$UniqueCID
gingivalST@meta.data <- metaST

coordst1 <- location[c(5,6)]

locationID <- coordst[,c("spot","idents")]
spotID <- locationID$idents
names(spotID) <- locationID$spot

location$idents=spotID[location$SpotID]

gingivalST2$idents <- location$idents

img <- STv3@images$spatial
imgCoor <- img@coordinates
lm1 <- lm(imagecol~row, data=imgCoor)
imgColPred <- predict(lm1,coordst1)
coordst1$imagecol <- imgColPred
lm2 <- lm(imagerow~col, data=imgCoor)
imgrowPred <- predict(lm2,coordst1)
coordst1$imagerow <- imgrowPred
coordst1$tissue=1
coordst1 <- coordst1[colnames(imgCoor)]
rownames(coordst1) <- location$UniqueCID
img@coordinates <- coordst1
gingivalST2@images$spatial <- img

metaST <- gingival_sc_2@meta.data[location$OriginalCID,]
rownames(metaST) <- location$UniqueCID
gingivalST2@meta.data <- metaST


SpatialDimPlot(gingivalST2,group.by = "idents",pt.size.factor = 5)
ggsave("results/dimplot_st.pdf")
SpatialDimPlot(gingivalST2,group.by = "celltype",pt.size.factor = 5)
SpatialDimPlot(gingivalST2,group.by = "orig.ident",pt.size.factor = 5)
ggsave("results/v2_dimplot_st_ident.pdf")


metaBar <- gingivalST2@meta.data[c("orig.ident","idents")]
ggplot(metaBar, aes(x=orig.ident, fill=idents)) + 
  geom_bar(position = "fill")+theme_classic()+scale_fill_manual(values = mycolor)+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1,size = 10,face = "bold"),title = NULL)
ggsave("results/v2barplot_cluster.pdf")
ggplot(metaBar, aes(x=idents, fill=orig.ident)) + 
  geom_bar(position = "fill")+theme_classic()+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1,size = 10,face = "bold"),title = NULL)
ggsave("results/v2_barplot_ident.pdf")
saveRDS(gingivalST2,"../important_processed_data/4.5_v3_gingival_2st.Rds")


mycolor<-colorRampPalette(brewer.pal(8,'Spectral'))(length(unique(metaBar$idents)))