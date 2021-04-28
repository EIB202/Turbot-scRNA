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

#Subset T cells
###############################################
T_cells_hk_mock<-subset(HK_0d, idents =c('3','7'))
T_cells_hk_mock$stim <- "0"
T_cells_hk_0.5d<-subset(HK_0.5d, idents =c('1','3'))
T_cells_hk_0.5d$stim <- "0.5"
T_cells_hk_3d<-subset(HK_3d, idents =c('5','6'))
T_cells_hk_3d$stim <- "3"
T_cells_hk_7d<-subset(HK_7d, idents =c('3','4','6'))
T_cells_hk_7d$stim <- "7"
T_cells_hk_14d<-subset(HK_14d, idents =c('2','3',"15"))
T_cells_hk_14d$stim <- "14"
T_cells_hk_21d<-subset(HK_21d, idents =c('0','3','13'))
T_cells_hk_21d$stim <- "21"

#Integrate T cells of head kidney
###############################################
T_cells_hk.anchors <- FindIntegrationAnchors(object.list = list(T_cells_hk_mock, T_cells_hk_0.5d, T_cells_hk_3d, T_cells_hk_7d, T_cells_hk_14d,T_cells_hk_21d), dims = 1:30)
T_cells_hk.combined <- IntegrateData(anchorset = T_cells_hk.anchors, dims = 1:30)
DefaultAssay(T_cells_hk.combined) <- "integrated"

#Dimensional reduction and cluster cells
###############################################
T_cells_hk.combined <- ScaleData(T_cells_hk.combined, verbose = FALSE)
T_cells_hk.combined <- RunPCA(T_cells_hk.combined, npcs = 30, verbose = FALSE)
T_cells_hk.combined <- FindNeighbors(T_cells_hk.combined, dims = 1:20)
T_cells_hk.combined <- FindClusters(T_cells_hk.combined, resolution = 0.5)
T_cells_hk.combined <- RunUMAP(T_cells_hk.combined, dims = 1:20)
DimPlot(T_cells_hk.combined, reduction = "umap",label = TRUE,label.size = 5)+NoLegend()+
  theme(axis.line = element_line(size=1.2, colour = "black"))+
  theme(axis.ticks = element_line(size=1, colour = "black"))+
  theme(axis.text= element_text(size = 15, face = "bold"))+
  theme(axis.title = element_text(size = 15, face = "bold"))

#Find markers
###############################################
logFCfilter=0.3
adjPvalFilter=0.05
T_cells_hk.markers <- FindAllMarkers(object = T_cells_hk.combined,
                                 assay = "RNA",
                                 only.pos = FALSE,
                                 min.pct = 0.25,
                                 logfc.threshold = logFCfilter)
sig.markers=T_cells_hk.markers[(abs(as.numeric(as.vector(T_cells_hk.markers$avg_logFC)))>logFCfilter & as.numeric(as.vector(T_cells_hk.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.markers,file="T_cells_hk_combined.markers.xls",sep="\t",row.names=F,quote=F)

#Rename clusters
###############################################
new.cluster.ids <- c("CD4+ naive T","CD4+Foxp3- Tregs","CD4-CD8- T","CD8+ CTLs","NK cells","CD8+ naive T","Doublelets1",
                     "Doublelets2","Doublelets3","ILCs3-1-like","CD4+Foxp3+ Tregs","Doublelets4","Doublelets5")
names(new.cluster.ids) <- levels(T_cells_hk.combined)
T_cells_hk.combined <- RenameIdents(T_cells_hk.combined, new.cluster.ids)

T_cells_hk.combined.filter<-subset(T_cells_hk.combined,idents =c('CD4+ naive T','CD4+Foxp3- Tregs','CD4-CD8- T',
                                                                 'CD8+ CTLs','NK cells','CD8+ naive T','ILCs3-1-like',"CD4+Foxp3+ Tregs"))

DimPlot(T_cells_hk.combined.filter, reduction = "umap",label = T,label.size = 5,repel = T)+NoLegend()+
  theme(axis.line = element_blank())+
  theme(axis.ticks = element_line(size=3, colour = "black"))+
  theme(axis.text= element_text(size = 15, face = "bold"))+
  theme(axis.title = element_blank())

markers.to.plot <- c("cd3gd","cd3e","trac","trbc1","trbc2","cd4-1","cd4-2","cd8a","cd8b","nkl.3",
                     "cd226","gzma","gzmb","faslg","foxp3","ctla4","ikzf4","il10","rorc","tnfa")

DotPlot(T_cells_hk.combined.filter,assay = "RNA",features = markers.to.plot,cols = c("gray", "springgreen4"), dot.scale = 8,
        col.min = -2.5, col.max = 2.5, scale.by = "size", scale.max = 100) + 
  RotatedAxis()+theme(axis.line = element_line(size=1.2, colour = "black"))+
  labs(y="Clusters")+theme(axis.title.y = element_text(size = 20, face = "bold"))+
  theme(axis.ticks = element_line(size=1, colour = "black"))+
  theme(axis.text= element_text(size = 15, face = "bold"))+
  theme(axis.text.x = element_text(size = 15, face = "bold.italic"))+
  theme(axis.title.x = element_blank())+
  theme(legend.position=("right"))

saveRDS(T_cells_hk.combined.filter,"T_cells_hk.combined.filter.rds")

#Average expression analysis(example of CD8+ CTLs)
###############################################
CTLs<-subset(T_cells_hk.combined.filter,idents =c("CD8+ CTLs"))
Idents(CTLs) <- "stim"
avg.CTLs <- as.data.frame(log1p(AverageExpression(CTLs, verbose = FALSE)$RNA))
write.table(avg.CTLs, file = "avg.CTLs.txt",sep = "\t",quote = F)





