library(Seurat)
library(limma)
library(dplyr)
library(magrittr)
library(ggplot2)

#Read data
###############################################
matrix.hk <- read.table("dataset/HK_CTRL.txt",sep="\t",header=T,check.names=F)
matrix.hk=as.matrix(matrix.hk)
rownames(matrix.hk)=matrix.hk[,1]
exp.hk=matrix.hk [,2:ncol(matrix.hk)]
dimnames=list(rownames(exp.hk),colnames(exp.hk))
data.hk=matrix(as.numeric(as.matrix(exp.hk)),nrow=nrow(exp.hk),dimnames=dimnames)
data.hk=avereps(data.hk)

#Create Seurat project and QC
###############################################
hk <- CreateSeuratObject(counts = data.hk, project = "Head Kidney", min.cells = 3, min.features = 200)
hk[["percent.mt"]] <- PercentageFeatureSet(hk, pattern = "^MT-")

VlnPlot(hk, features = c("nFeature_RNA"))+NoLegend()+
  theme(plot.title = element_blank())+
  labs(y="Number of genes")+theme(axis.title.y = element_text(size = 15, face = "bold"))+
  theme(axis.line = element_line(size=1.2, colour = "black"))+
  theme(axis.ticks = element_line(size=1, colour = "black"))+
  theme(axis.text.y = element_text(size = 15, face = "bold"))+
  theme(axis.text.x = element_blank())+
  theme(axis.title.x = element_blank())

VlnPlot(hk, features = c("nCount_RNA"))+NoLegend()+
  theme(plot.title = element_blank())+
  labs(y="Number of counts")+theme(axis.title.y = element_text(size = 15, face = "bold"))+
  theme(axis.line = element_line(size=1.2, colour = "black"))+
  theme(axis.ticks = element_line(size=1, colour = "black"))+
  theme(axis.text.y = element_text(size = 15, face = "bold"))+
  theme(axis.text.x = element_blank())+
  theme(axis.title.x = element_blank())

VlnPlot(hk, features = c("percent.mt"))+NoLegend()+
  theme(plot.title = element_blank())+
  labs(y="Percent of MT")+theme(axis.title.y = element_text(size = 15, face = "bold"))+
  theme(axis.line = element_line(size=1.2, colour = "black"))+
  theme(axis.ticks = element_line(size=1, colour = "black"))+
  theme(axis.text.y = element_text(size = 15, face = "bold"))+
  theme(axis.text.x = element_blank())+
  theme(axis.title.x = element_blank())

hk <- subset(hk, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
hk <- NormalizeData(hk, normalization.method = "LogNormalize", scale.factor = 10000)
hk <- FindVariableFeatures(hk, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(object = hk)

#Dimensional reduction and cluster cells to extract T/B cells for step 3
###############################################
all.genes <- rownames(hk)
hk <- ScaleData(hk, features = all.genes)
hk <- RunPCA(hk, features = VariableFeatures(object = hk))

hk <- FindNeighbors(hk, dims = 1:20)
hk <- FindClusters(hk, resolution = 0.5)
hk <- RunUMAP(hk, dims = 1:20)
DimPlot(hk, reduction = "umap",label = T)+NoLegend()+
  theme(axis.line = element_line(size=1.2, colour = "black"))+
  theme(axis.ticks = element_line(size=1, colour = "black"))+
  theme(axis.text= element_text(size = 15, face = "bold"))+
  theme(axis.title = element_text(size = 15, face = "bold"))

markers.to.plot <- c("cd3e","cd3gd","trac","cd79a", "cd79b","ighm-m","ight-m","cd74a", "cd80/86", "mhcII-1", "mhcII-2", "cd209a","fcer1a1","mpeg1.2","tnfa","mpx","mmp9","cxcr1","rorc")
DotPlot(hk,features = markers.to.plot, cols = c("gray", "springgreen4"), dot.scale = 8) + 
  RotatedAxis()+theme(axis.line = element_line(size=1.2, colour = "black"))+
  labs(y="Clusters")+theme(axis.title.y = element_text(size = 20, face = "bold"))+
  theme(axis.ticks = element_line(size=1, colour = "black"))+
  theme(axis.text= element_text(size = 15, face = "bold"))+
  theme(axis.text.x = element_text(size = 15, face = "bold.italic"))+
  theme(axis.title.x = element_blank())

saveRDS(hk,"HK_CTRL.rds")
