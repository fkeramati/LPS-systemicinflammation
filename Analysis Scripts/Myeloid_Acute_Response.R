suppressMessages(library(Seurat))
suppressMessages(library(scater))
suppressMessages(library(scDblFinder))
suppressMessages(library(HGNChelper))
suppressMessages(library(SeuratWrappers))
suppressMessages(library(harmony))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(ggthemes))
suppressMessages(library(miloR))
suppressMessages(library(ggrastr))

# You have to generate the following Seurat object by runing the first part of the analysis script, presented in this depository.
endo = readRDS("All_Donor_TimePoint_Tissue_logNorm_MNN.rds")
endo = subset(endo, subset = TimePoint == "D7", invert = T)

# Myeloid Lineage #
myelo_tmp = subset(endo, idents = c("HSC", "GMP", "Pro Mono", "Monocyte", "Non Cla Mono", "cDC", 
                                    "Inflam Mono", "Inflam cDC"))
count.tmp = GetAssayData(myelo_tmp, slot = "counts")
myelo = CreateSeuratObject(counts = count.tmp, project = "myelo", min.cells = 0, min.features = 0)
myelo@meta.data$Tissue = unlist(strsplit(names(myelo$orig.ident), split = "_"))[c(T,F,F,F)]
myelo@meta.data$Donor = unlist(strsplit(names(myelo$orig.ident), split = "_"))[c(F,T,F,F)]
myelo@meta.data$TimePoint = unlist(strsplit(names(myelo$orig.ident), split = "_"))[c(F,F,T,F)]
myelo$CellType = endo$CellType
myelo$orig.ident = endo$orig.ident
myelo$percent.mt = endo$percent.mt
myelo$percent.rb = endo$percent.rb
rm(myelo_tmp, count.tmp)

# Normalize and Scale
myelo = NormalizeData(myelo) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData()
s.genes = cc.genes.updated.2019$s.genes
g2m.genes = cc.genes.updated.2019$g2m.genes
myelo = CellCycleScoring(myelo, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)

# PCA, UMAP visualization and clustering
myelo = RunPCA(myelo, ndims.print = 1:30, verbose = T) %>% RunUMAP(dims = 1:30) %>% 
  FindNeighbors(dims = 1:30, k.param = 20) %>% FindClusters(resolution = 1)

DimPlot(myelo, label = TRUE, pt.size = 0.1)
DimPlot(myelo, label = F, pt.size = 0.1, group.by = "Tissue", shuffle = T)
DimPlot(myelo, label = F, pt.size = 0.1, group.by = "Donor", shuffle = T)
DimPlot(myelo, label = F, pt.size = 0.1, group.by = "TimePoint", shuffle = T)
DimPlot(myelo, label = T, group.by = "CellType")

mark = FindAllMarkers(myelo, only.pos = T, min.pct = 0.3, min.diff.pct = 0.1)

# Clusters 13 is not a myeloid cluster (remove it)
myelo = subset(myelo, idents = 13, invert = T)

count.tmp = GetAssayData(myelo, slot = "counts")
myelo2 = CreateSeuratObject(counts = count.tmp, project = "myelo", min.cells = 0, min.features = 0)
myelo2@meta.data$Tissue = unlist(strsplit(names(myelo2$orig.ident), split = "_"))[c(T,F,F,F)]
myelo2@meta.data$Donor = unlist(strsplit(names(myelo2$orig.ident), split = "_"))[c(F,T,F,F)]
myelo2@meta.data$TimePoint = unlist(strsplit(names(myelo2$orig.ident), split = "_"))[c(F,F,T,F)]
myelo2$CellType = myelo$CellType
myelo2$orig.ident = myelo$orig.ident
myelo2$percent.mt = myelo$percent.mt
myelo2$percent.rb = myelo$percent.rb
rm(count.tmp)

# Normalize and Scale
myelo2 = NormalizeData(myelo2) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData()
s.genes = cc.genes.updated.2019$s.genes
g2m.genes = cc.genes.updated.2019$g2m.genes
myelo2 = CellCycleScoring(myelo2, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)

# PCA, UMAP visualization and clustering
myelo2 = RunPCA(myelo2, ndims.print = 1:30, verbose = T) %>% RunUMAP(dims = 1:30) %>% 
  FindNeighbors(dims = 1:30, k.param = 20) %>% FindClusters(resolution = 1)

DimPlot(myelo2, label = TRUE, pt.size = 0.1)
DimPlot(myelo2, label = F, pt.size = 0.1, group.by = "Tissue", shuffle = T)
DimPlot(myelo2, label = F, pt.size = 0.1, group.by = "Donor", shuffle = T)
DimPlot(myelo2, label = F, pt.size = 0.1, group.by = "TimePoint", shuffle = T)
myelo2$TisDnr = paste(myelo2$Tissue, myelo2$Donor)
saveRDS(myelo2, "Myeloid_D0_4h_NoBatchCorrect.rds")

########################
### Batch Correction ###
########################
# There is a donor and sequencing batch effect, correct for it using MNN
myelo_mnn = RunFastMNN(object.list = SplitObject(myelo2, split.by = "TisDnr"), features = 2000, verbose = T)
myelo_mnn = RunUMAP(myelo_mnn, dims = 1:30, reduction = "mnn", min.dist = 0.3)
myelo_mnn = FindNeighbors(myelo_mnn, dims = 1:30, reduction = "mnn", k.param = 20)
myelo_mnn = FindClusters(myelo_mnn, resolution = 2)
myelo_mnn@assays$RNA@scale.data = myelo2@assays$RNA@scale.data

myelo_mnn = FindSubCluster(myelo_mnn, cluster = 6, graph.name = "RNA_snn", resolution = 0.2)
DimPlot(myelo_mnn, label = T, group.by = "sub.cluster")
Idents(myelo_mnn) = myelo_mnn$sub.cluster

myelo_mnn = FindSubCluster(myelo_mnn, cluster = 7, graph.name = "RNA_snn", resolution = 0.3)
DimPlot(myelo_mnn, label = T, group.by = "sub.cluster")
Idents(myelo_mnn) = myelo_mnn$sub.cluster

myelo_mnn = FindSubCluster(myelo_mnn, cluster = 10, graph.name = "RNA_snn", resolution = 0.2)
DimPlot(myelo_mnn, label = T, group.by = "sub.cluster")
Idents(myelo_mnn) = myelo_mnn$sub.cluster

myelo_mnn = FindSubCluster(myelo_mnn, cluster = 13, graph.name = "RNA_snn", resolution = 0.3)
DimPlot(myelo_mnn, label = T, group.by = "sub.cluster")
Idents(myelo_mnn) = myelo_mnn$sub.cluster

myelo_mnn = FindSubCluster(myelo_mnn, cluster = "13_0", graph.name = "RNA_snn", resolution = 0.3)
DimPlot(myelo_mnn, label = T, group.by = "sub.cluster")
Idents(myelo_mnn) = myelo_mnn$sub.cluster

myelo_mnn = FindSubCluster(myelo_mnn, cluster = 14, graph.name = "RNA_snn", resolution = 0.2)
DimPlot(myelo_mnn, label = T, group.by = "sub.cluster")
Idents(myelo_mnn) = myelo_mnn$sub.cluster

myelo_mnn = FindSubCluster(myelo_mnn, cluster = 17, graph.name = "RNA_snn", resolution = 0.3)
DimPlot(myelo_mnn, label = T, group.by = "sub.cluster")
Idents(myelo_mnn) = myelo_mnn$sub.cluster

myelo_mnn = FindSubCluster(myelo_mnn, cluster = 19, graph.name = "RNA_snn", resolution = 0.2)
DimPlot(myelo_mnn, label = T, group.by = "sub.cluster")
Idents(myelo_mnn) = myelo_mnn$sub.cluster

myelo_mnn = RenameIdents(myelo_mnn, "0" = "Inflam Pro Mono", "1" = "Inflam Pro Mono", "2" = "Pro Mono", "3" = "Pro Mono", "4" = "Non Cla Mono",
                         "5" = "Monocyte", "6_0" = "Inflam Pro Mono", "6_1" = "Inflam Mono", "7_0" = "Monocyte", "7_1" = "Inter Mono", "8" = "Monocyte", "9" = "Monocyte",
                         "10_0" = "Pro Mono", "10_1" = "Cycl Pro Mono", "11" = "Pro Mono", "12" = "Inflam Mono", "13_0_0" = "GMP", 
                         "13_0_1" = "Inflam GMP", "13_1" = "GMP", "13_2" = "Cycl Pro Mono", "13_3" = "Cycl Pro Mono", "14_0" = "LT-HSC", "14_1" = "Inflam LT-HSC",
                         "15" = "Inflam Pro Mono", "16" = "ST-HSC", "17_0" = "Monocyte", "17_1" = "cDC", "18" = "cDC", 
                         "19_0" = "Inflam Cycl Pro Mono", "19_1" = "Inflam GMP", "20" = "Inflam cDC", "21" = "Prog cDC")
df = data.frame(Cell = colnames(myelo_mnn), Ident = as.character(Idents(myelo_mnn)), TimePoint = myelo_mnn$TimePoint)
df[which(df$Ident == "ST-HSC" & df$TimePoint == "4h"),"Ident"] = "Inflam ST-HSC"
Idents(myelo_mnn) = df$Ident
Idents(myelo_mnn) = factor(Idents(myelo_mnn), levels = c("LT-HSC","ST-HSC","GMP","Cycl Pro Mono","Pro Mono","Monocyte",
                                                         "Inter Mono","Non Cla Mono","Prog cDC","cDC","Inflam LT-HSC",
                                                         "Inflam ST-HSC","Inflam GMP","Inflam Cycl Pro Mono",
                                                         "Inflam Pro Mono","Inflam Mono","Inflam cDC"))
DimPlot(myelo_mnn, label = T, label.size = 3) + NoLegend()
saveRDS(myelo_mnn, "Myeloid_D0_4h_MNN.rds")