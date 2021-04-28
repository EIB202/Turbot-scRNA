library(Seurat)
library(ggplot2)

#Read data
###############################################
HK_0d<-readRDS(file = "dataset/HK_0d.rds")
HK_0.5d<-readRDS(file = "dataset/HK_0.5d.rds")
HK_3d<-readRDS(file = "dataset/HK_3d.rds")
HK_7d<-readRDS(file = "dataset/HK_7d.rds")
HK_14d<-readRDS(file = "dataset/HK_14d.rds")
HK_21d<-readRDS(file = "dataset/HK_21d.rds")

#Annotation
###############################################
HK_0d$stim <- "0"
HK_0.5d$stim <- "0.5"
HK_3d$stim <- "3"
HK_7d$stim <- "7"
HK_14d$stim <- "14"
HK_21d$stim <- "21"

#Integrate data
###############################################
hk_combined.anchors <- FindIntegrationAnchors(object.list = list(HK_0d, HK_0.5d, HK_3d, HK_7d, HK_14d,HK_21d), dims = 1:30)
hk.combined <- IntegrateData(anchorset = hk_combined.anchors, dims = 1:30)
DefaultAssay(hk.combined) <- "integrated"

#Dimensional reduction and cluster cells
###############################################
hk.combined <- ScaleData(hk.combined, verbose = FALSE)
hk.combined <- RunPCA(hk.combined, npcs = 30, verbose = FALSE)

hk.combined <- FindNeighbors(hk.combined, dims = 1:20)
hk.combined <- FindClusters(hk.combined, resolution = 0.5)
hk.combined <- RunUMAP(hk.combined, dims = 1:20)
DimPlot(hk.combined, reduction = "umap",label = T,label.size = 5,repel = T)+NoLegend()+
  theme(axis.line = element_blank())+
  theme(axis.ticks = element_line(size=3, colour = "black"))+
  theme(axis.text= element_text(size = 20, face = "bold"))+
  theme(axis.title = element_blank())

#Dotplot of conserved marker genes
###############################################
markers.to.plot <- c("cd3e","cd3gd","trac","cd79a", "cd79b","ighm-m","ight-m","cd74a", "cd80/86", "mhcII-1", "mhcII-2", "cd209a","fcer1a1","mpeg1.2",'tnfb',"mpx","mmp9","cxcr1","rorc")

DotPlot(hk.combined,assay = "RNA",features = markers.to.plot, cols = c("gray", "springgreen4"), dot.scale = 8,
        col.min = -1, col.max = 2, scale.max = 100) + 
  RotatedAxis()+theme(axis.line = element_line(size=1.2, colour = "black"))+
  labs(y="Clusters")+theme(axis.title.y = element_text(size = 20, face = "bold"))+
  theme(axis.ticks = element_line(size=1, colour = "black"))+
  theme(axis.text= element_text(size = 15, face = "bold"))+
  theme(axis.text.x = element_text(size = 15, face = "bold.italic"))+
  theme(axis.title.x = element_blank())

#Find markers
###############################################
logFCfilter=0.5
adjPvalFilter=0.05
hk.markers <- FindAllMarkers(object = hk.combined.filter,
                             only.pos = FALSE,
                             assay = "RNA",
                             min.pct = 0.25,
                             logfc.threshold = logFCfilter)
sig.markers=hk.markers[(abs(as.numeric(as.vector(hk.markers$avg_logFC)))>logFCfilter & as.numeric(as.vector(hk.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.markers,file="hk.markers.xls",sep="\t",row.names=F,quote=F)

#Rename clusters
###############################################
new.cluster.ids <- c("B cells","Neutrophils","T cells","Neutrophils","T cells",
                     "Non-immune cells","Neutrophils","DCs","T cells","Neutrophils",
                     "Neutrophils","B cells","B cells","B cells","Doublelets","B cells",
                     "B cells","Non-immune cells","Doublelets","MΦ","Doublelets","Non-immune cells",
                     "Non-immune cells","T cells","Doublelets","DCs","DCs")
names(new.cluster.ids) <- levels(hk.combined)
hk.combined <- RenameIdents(hk.combined, new.cluster.ids)

hk.combined.filter<-subset(hk.combined,idents =c('T cells','B cells','Neutrophils','DCs',"MΦ"))

DimPlot(hk.combined.filter, reduction = "umap",label = TRUE,label.size = 5,repel = T)+NoLegend()+
  theme(axis.line = element_line(size=1.2, colour = "black"))+
  theme(axis.ticks = element_line(size=1, colour = "black"))+
  theme(axis.text= element_text(size = 15, face = "bold"))+
  theme(axis.title = element_text(size = 15, face = "bold"))

saveRDS(hk.combined.filter, "hk.combined.filter.rds")


#Heatmap of Top10 genes
###############################################
top10 <- hk.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
color <- colorRampPalette(c("white","springgreen4"))(100)
DoHeatmap(object = hk.combined.filter, features = top10$gene, size =5,assay = "RNA",slot = "data")+
  scale_fill_gradientn(colors=color,guide="colorbar")+
  theme(axis.text= element_text(size = 10,colour = "black",face = "bold.italic"))+
  theme(legend.text = element_text(colour="black", size = 12),legend.title =  element_text(colour="black", size = 12, face = "bold"))+
  theme(legend.position="right")

#Pathway socre(datasets of spleen, gill and posterior intestine should be integrated additionally)
###############################################
combined<-merge(x=hk.combined.filter,y=c(spleen.combined.filter,gill.combined.filter,PI.combined.filter))
combined$stim <- factor(combined$stim, levels = c("0",'0.5','3','7','14','21'))

ISG<-read.table("Response to type I interferon.txt",header = T,check.names=F,sep = "\t")
ISG<-list(ISG$Gene)

combined <- AddModuleScore( object = combined, features = ISG, ctrl = 100, name = 'ISG',assay = "RNA")

VlnPlot(combined,features = "ISG1",split.by = "stim",pt.size = 0,assay = "RNA",y.max = 0.6)+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.line = element_line(size=1.2, colour = "black"))+
  theme(axis.ticks = element_line(size=1, colour = "black"))+
  theme(axis.text= element_text(size = 15, face = "bold",angle = 0))+
  theme(axis.text.x = element_blank())+
  theme(axis.title.x = element_blank())+
  stat_summary(fun.data = mean_sdl,geom = "pointrange",color="black",size=1,position = position_dodge(0.905))


