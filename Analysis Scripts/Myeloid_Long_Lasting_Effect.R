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
  FindNeighbors(dims = 1:30, k.param = 20) %>% FindClusters(resolution = 1.5)

DimPlot(myelo, label = TRUE, pt.size = 0.1)
DimPlot(myelo, label = F, pt.size = 0.1, group.by = "Tissue", shuffle = T)
DimPlot(myelo, label = F, pt.size = 0.1, group.by = "Donor", shuffle = T)
DimPlot(myelo, label = F, pt.size = 0.1, group.by = "TimePoint", shuffle = T)
DimPlot(myelo, label = T, group.by = "CellType")

mark = FindAllMarkers(myelo, only.pos = T, min.pct = 0.3, min.diff.pct = 0.1)

# Clusters 20, 22 are not a myeloid clusters (remove it)
myelo = subset(myelo, idents = c(20,22), invert = T)

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
saveRDS(myelo2, "Myeloid_All_TimePoint_NoBatchCorrect.rds")

# Subsample Donor Sbj3 PBMC D0 and D7 as they are exceedingly high in number (distorting UMAP)
tmp1 = subset(myelo2, subset = orig.ident == "PBMC_Sbj3_D0")
tmp2 = subset(myelo2, subset = orig.ident == "PBMC_Sbj3_D7")

tmp3 = tmp1[, sample(colnames(tmp1), size = 1500, replace=F)]
tmp4 = tmp2[, sample(colnames(tmp2), size = 1500, replace=F)]

cells1 = setdiff(colnames(tmp1), colnames(tmp3))
cells2 = setdiff(colnames(tmp2), colnames(tmp4))

myelo2 = subset(myelo2, cells = c(cells1, cells2), invert = T)

# There is a donor and sequencing batch effect, correct for it using MNN
myelo_mnn = RunFastMNN(object.list = SplitObject(myelo2, split.by = "TisDnr"), features = 1000, verbose = T)
myelo_mnn = RunUMAP(myelo_mnn, dims = 1:30, reduction = "mnn", min.dist = 0.4)
myelo_mnn = FindNeighbors(myelo_mnn, dims = 1:30, reduction = "mnn", k.param = 20)
myelo_mnn = FindClusters(myelo_mnn, resolution = 1.5)
myelo_mnn@assays$RNA@scale.data = myelo2@assays$RNA@scale.data
saveRDS(myelo_mnn, "Myeloid_All_TimePoint_MNN.rds")
