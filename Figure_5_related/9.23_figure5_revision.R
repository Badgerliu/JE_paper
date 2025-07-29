library(RColorBrewer)
library(viridis)
st <- readRDS("important_processed_data/st_seurat_v3.Rds")
mycol_set2 <- colorRampPalette(brewer.pal(8,'Set2'))(9)
SpatialDimPlot(st,group.by = "bayesIdent",pt.size.factor=5,image.alpha=0.4,stroke = NA) + scale_fill_manual(values=mycol_set2)
ggsave("results/9.23_figure_revision/bayes_ident.pdf",width = 6,height = 6)
SpatialDimPlot(st,group.by = "bayesIdent",pt.size.factor=5,image.alpha=0.4) + scale_fill_manual(values=mycol_set2)
ggsave("results/9.23_figure_revision/bayes_ident.png",width = 6,height = 6)


scMerge <- readRDS('important_processed_data/6.18_merge_3day_6sample.Rds')
scMerge$day <- NA
scMerge$day[scMerge$orig.ident%in%c("D0","D0_K5_GFP")] ="D0"
scMerge$day[scMerge$orig.ident%in%c("D3","D3_K5_GFP")] ="D3"
scMerge$day[scMerge$orig.ident%in%c("D5","D5_K5_GFP")] ="D5"

d0 <- scMerge[,scMerge$day=="D0"]
d3 <- scMerge[,scMerge$day=="D3"]
d5 <- scMerge[,scMerge$day=="D5"]
p1 <- SpatialFeaturePlot(d0,c("S100a9"),pt.size.factor=5,image.alpha=0.4,stroke = NA)&scale_fill_viridis_c(option ="D")
p1
ggsave(plot = p1,"results/9.23_figure_revision/S100a9_d0.pdf",width = 6,height = 6)
p2 <-  SpatialFeaturePlot(d3,c("S100a9"),pt.size.factor=5,image.alpha=0.4, stroke = NA)&scale_fill_viridis_c(option ="D")
p2
ggsave(plot = p2,"results/9.23_figure_revision/S100a9_d3.pdf",width = 6,height = 6)
p3 <-  SpatialFeaturePlot(d5,c("S100a9"),pt.size.factor=5,image.alpha=0.4, stroke = NA)&scale_fill_viridis_c(option ="D")
p3
ggsave(plot = p3,"results/9.23_figure_revision/S100a9_d5.pdf",width = 6,height = 6)



go_meta <- read.csv("processed_data/6.18_aucell_3_stage_meta.csv",row.names = 1)
scMerge@meta.data[colnames(go_meta)] <- go_meta
SpatialFeaturePlot(scMerge,c("GOBP_WOUND_HEALING"),pt.size.factor=5,image.alpha=0.4, stroke = NA)&scale_fill_viridis_c(option ="D")
ggsave("results/9.23_figure_revision/go_wound_whole.pdf",width = 6,height = 6)


d0 <- scMerge[,scMerge$day=="D0"]
d3 <- scMerge[,scMerge$day=="D3"]
d5 <- scMerge[,scMerge$day=="D5"]

p1 <- SpatialFeaturePlot(d0,c("GOBP_WOUND_HEALING"),pt.size.factor=5,image.alpha=0.4,stroke = NA)&scale_fill_viridis_c(option ="D")
p2 <-  SpatialFeaturePlot(d3,c("GOBP_WOUND_HEALING"),pt.size.factor=5,image.alpha=0.4,stroke = NA)&scale_fill_viridis_c(option ="D")
p3 <-  SpatialFeaturePlot(d5,c("GOBP_WOUND_HEALING"),pt.size.factor=5,image.alpha=0.4,stroke = NA)&scale_fill_viridis_c(option ="D")
p = p1|p2|p3
p
ggsave("results/9.23_figure_revision/go_wound_split.pdf",width = 9,height = 4)

