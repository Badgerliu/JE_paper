go_meta <- read.csv("processed_data/6.18_aucell_3_stage_meta.csv",row.names = 1)
scMergeSeurat@meta.data[colnames(go_meta)] <- go_meta
SpatialFeaturePlot(scMergeSeurat,c("GOBP_WOUND_HEALING"),pt.size.factor=5,image.alpha=0.4)&scale_fill_viridis_c(option ="D")
ggsave("results/6.18_go_project/6.18_go_wound_whole.pdf",width = 5,height = 4)

scMergeSeurat$day <- NA
scMergeSeurat$day[scMergeSeurat$orig.ident%in%c("D0","D0_K5_GFP")] ="D0"
scMergeSeurat$day[scMergeSeurat$orig.ident%in%c("D3","D3_K5_GFP")] ="D3"
scMergeSeurat$day[scMergeSeurat$orig.ident%in%c("D5","D5_K5_GFP")] ="D5"

d0 <- scMergeSeurat[,scMergeSeurat$day=="D0"]
d3 <- scMergeSeurat[,scMergeSeurat$day=="D3"]
d5 <- scMergeSeurat[,scMergeSeurat$day=="D5"]

p1 <- SpatialFeaturePlot(d0,c("GOBP_WOUND_HEALING"),pt.size.factor=5,image.alpha=0.4)&scale_fill_viridis_c(option ="D")
p2 <-  SpatialFeaturePlot(d3,c("GOBP_WOUND_HEALING"),pt.size.factor=5,image.alpha=0.4)&scale_fill_viridis_c(option ="D")
p3 <-  SpatialFeaturePlot(d5,c("GOBP_WOUND_HEALING"),pt.size.factor=5,image.alpha=0.4)&scale_fill_viridis_c(option ="D")
p = p1|p2|p3
p
ggsave("results/6.18_go_project/6.18_go_wound_split.pdf",width = 9,height = 4)


SpatialFeaturePlot(scMergeSeurat,c("GOBP_WOUND_HEALING"),pt.size.factor=5,image.alpha=0.4)&scale_fill_viridis_c(option ="H")

p1 <- SpatialFeaturePlot(d0,c("GOBP_WOUND_HEALING"),pt.size.factor=5,image.alpha=0.4)&scale_fill_viridis_c(option ="H")
p2 <-  SpatialFeaturePlot(d3,c("GOBP_WOUND_HEALING"),pt.size.factor=5,image.alpha=0.4)&scale_fill_viridis_c(option ="H")
p3 <-  SpatialFeaturePlot(d5,c("GOBP_WOUND_HEALING"),pt.size.factor=5,image.alpha=0.4)&scale_fill_viridis_c(option ="H")
p = p1|p2|p3
p
ggsave("results/6.18_go_project/6.18_go_wound_split_colorH.pdf",width = 9,height = 4)


SpatialFeaturePlot(scMergeSeurat,c("pattern2"),pt.size.factor=5,image.alpha=0.4)&scale_fill_viridis_c(option ="D")
ggsave("results/6.18_plot_mod/6.18_pattern2.pdf",width = 5,height = 4)
