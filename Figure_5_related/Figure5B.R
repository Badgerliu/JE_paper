## Fig B plot-----------------------------
library(Seurat)
library(viridis)
#setwd("../")
source("../utils/Utils.R")
st <- readRDS("important_processed_data/st_seurat_v3.Rds")

st@images$spatial@scale.factors$lowres <- 0.595

st <- rotateSeuratImage(st,slide="spatial",rotation="R90")
SpatialFeaturePlot(st,"Krt5",pt.size.factor=5,image.alpha=0.4)

SpatialDimPlot(st,group.by = "celltype",pt.size.factor=5,image.alpha=0.4)
SpatialDimPlot(st,group.by = "bayesIdent",pt.size.factor=5,image.alpha=0.4)

# mycol_9 <- colorRampPalette(brewer.pal(8,'Set1'))(9)
# 
# SpatialDimPlot(st,group.by = "bayesIdent",pt.size.factor=5,image.alpha=0.4) + scale_fill_manual(values=mycol_9)
# SpatialDimPlot(st,group.by = "bayesIdent",pt.size.factor=5,image.alpha=0.4) + scale_fill_manual(values=mycol_9)
# SpatialDimPlot(st,group.by = "bayesIdent",pt.size.factor=5,image.alpha=0.4) + scale_fill_manual(values=mycol_9)
mycol_4 <- colorRampPalette(brewer.pal(4,'Set2'))(4)
SpatialDimPlot(st,group.by = "celltype",pt.size.factor=5,image.alpha=0.4)+ scale_fill_manual(values=mycol_4)
ggsave("results/FigureB/4.11_sc_celltype.pdf",width = 6,height = 6)
ggsave("results/FigureB/4.11_sc_celltype.png",width = 6,height = 6)


mycol_set2 <- colorRampPalette(brewer.pal(8,'Set2'))(9)
SpatialDimPlot(st,group.by = "bayesIdent",pt.size.factor=5,image.alpha=0.4) + scale_fill_manual(values=mycol_set2)
ggsave("results/FigureB/4.11_SpatialDimplot_set2.pdf",width = 6,height = 6)
# mycol_paired <- colorRampPalette(brewer.pal(8,'Paired'))(9)
# SpatialDimPlot(st,group.by = "bayesIdent",pt.size.factor=5,image.alpha=0.4) + scale_fill_manual(values=mycol_paired)

ggsave("results/FigureB/4.11_SpatialDimplot.pdf",width = 6,height = 6)
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "RdBu")))

SpatialFeaturePlot(st,"Krt5",pt.size.factor=5,image.alpha=0.4) +scale_fill_gradientn( colours=SpatialColors(n=100))
SpatialFeaturePlot(st,"Krt5",pt.size.factor=5,image.alpha=0.4)+scale_fill_viridis_c(option ="A")
SpatialFeaturePlot(st,"Krt5",pt.size.factor=5,image.alpha=0.4)+scale_fill_viridis_c(option ="D")

#== supp celltype heatmap-----------------------------
hm <- table(st$bayesIdent,st$celltype)%>%as.matrix()%>%as.data.frame()

library(reshape2)
library(circlize)
wide_df <- dcast(hm, Var1 ~ Var2, value.var = "Freq")%>%column_to_rownames("Var1")
wide_df <- t(scale(t(wide_df)))
pdf("results/FigureB/supp_heatmap.pdf",width = 6,height = 4)
ComplexHeatmap::Heatmap(t(wide_df),col = colorRamp2(c(-2, 0, 2), hcl_palette="Viridis"))
dev.off()

#== Figure C-----------------------
SpatialFeaturePlot(st,"Odam",pt.size.factor=5,image.alpha=0.4)+scale_fill_viridis_c(option ="D")
ggsave("results/FigureC/Odam.pdf",width = 6,height = 6)

SpatialFeaturePlot(st,"Dsg1a",pt.size.factor=5,image.alpha=0.4)+scale_fill_viridis_c(option ="D")
ggsave("results/FigureC/Dsg1a.pdf",width = 6,height = 6)
SpatialFeaturePlot(st,"Slc7a11",pt.size.factor=5,image.alpha=0.4)+scale_fill_viridis_c(option ="D")
ggsave("results/FigureC/Slc7a11.pdf",width = 6,height = 6)
ggsave("results/FigureC/Slc7a11.png",width = 6,height = 6)

SpatialFeaturePlot(st,"Lama5",pt.size.factor=5,image.alpha=0.4)+scale_fill_viridis_c(option ="D")
ggsave("results/FigureC/Lama5.pdf",width = 6,height = 6)
ggsave("results/FigureC/Lama5.png",width = 6,height = 6)

SpatialFeaturePlot(st,"Hmgb2",pt.size.factor=5,image.alpha=0.4)+scale_fill_viridis_c(option ="D")
ggsave("results/FigureC/Hmgb2.pdf",width = 6,height = 6)
ggsave("results/FigureC/Hmgb2.png",width = 6,height = 6)
SpatialFeaturePlot(st,"Hmgb2",pt.size.factor=5,image.alpha=0.4)+scale_fill_viridis_c(option ="D")
ggsave("results/FigureC/Hmgb2.pdf",width = 6,height = 6)
ggsave("results/FigureC/Hmgb2.png",width = 6,height = 6)

SpatialFeaturePlot(st,"Ggt1",pt.size.factor=5,image.alpha=0.4)+scale_fill_viridis_c(option ="D")
ggsave("results/FigureC/Ggt1.pdf",width = 6,height = 6)
ggsave("results/FigureC/Ggt1.png",width = 6,height = 6)

Idents(st) <- st$bayesIdent
markers <- FindAllMarkers(st)
specificity_base = 0.001
markers$specificity = ((markers$pct.1 + specificity_base) / (markers$pct.2 + specificity_base)) * markers$avg_log2FC
cluster_levels = c("IBL",  "JE_progenitor", "EBL", "sulcular_epi", "basal layer","gingiva progenitor", "spinous+granular","cornified", 
                   "stratum corneum")

comparisonGene <- markers%>%
  group_by(cluster)%>%
  arrange(desc(specificity))%>%
  top_n(2)%>%
  arrange(fct_relevel(cluster, cluster_levels))%>%
  ungroup()%>%
  select(gene)%>%
  unlist()
levels(st) <- cluster_levels

p <- DotPlot(st,assay = "RNA",features = comparisonGene)+theme(
  axis.text.x = element_text(angle = 70,size = 10,face = "bold",hjust  = 1)
)+ scale_color_viridis_c()
ggsave("results/FigureC/Supp_sigGene.pdf",width = 10,height = 6)

saveRDS(st,"important_processed_data/st_seurat_v3.Rds")
