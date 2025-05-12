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

# Lymphoid B Lineage #
lympho_tmp = subset(endo, idents = c("Pro B", "Pre B", "B Naive", "B Mem", "Plasma", "pDC"))
count.tmp = GetAssayData(lympho_tmp, slot = "counts")
lympho = CreateSeuratObject(counts = count.tmp, project = "lympho", min.cells = 0, min.features = 0)
lympho@meta.data$Tissue = unlist(strsplit(names(lympho$orig.ident), split = "_"))[c(T,F,F,F)]
lympho@meta.data$Donor = unlist(strsplit(names(lympho$orig.ident), split = "_"))[c(F,T,F,F)]
lympho@meta.data$TimePoint = unlist(strsplit(names(lympho$orig.ident), split = "_"))[c(F,F,T,F)]
lympho$CellType = endo$CellType
lympho$orig.ident = endo$orig.ident
lympho$percent.mt = endo$percent.mt
lympho$percent.rb = endo$percent.rb
rm(lympho_tmp, count.tmp)

# Normalize and Scale
lympho = NormalizeData(lympho) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData()
s.genes = cc.genes.updated.2019$s.genes
g2m.genes = cc.genes.updated.2019$g2m.genes
lympho = CellCycleScoring(lympho, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)

# PCA, UMAP visualization and clustering
lympho = RunPCA(lympho, ndims.print = 1:30, verbose = T) %>% RunUMAP(dims = 1:30) %>% 
  FindNeighbors(dims = 1:30, k.param = 20) %>% FindClusters(resolution = 1)

DimPlot(lympho, label = TRUE, pt.size = 0.1)
DimPlot(lympho, label = F, pt.size = 0.1, group.by = "Tissue", shuffle = T)
DimPlot(lympho, label = F, pt.size = 0.1, group.by = "Donor", shuffle = T)
DimPlot(lympho, label = F, pt.size = 0.1, group.by = "TimePoint", shuffle = T)
DimPlot(lympho, label = T, group.by = "CellType")

saveRDS(lympho, "LymphoidB_D0_4h_NoBatchCorrect.rds")

########################
### Batch Correction ###
########################
# There is a donor and sequencing batch effect, correct for it using MNN
lympho$TisDnr = paste(lympho$Tissue, lympho$Donor)
lympho_mnn = RunFastMNN(object.list = SplitObject(lympho, split.by = "TisDnr"), features = 3000, verbose = T)
lympho_mnn = RunUMAP(lympho_mnn, dims = 1:10, reduction = "mnn", min.dist = 0.3)
lympho_mnn = FindNeighbors(lympho_mnn, dims = 1:10, reduction = "mnn", k.param = 20)
lympho_mnn = FindClusters(lympho_mnn, resolution = 1.5)
lympho_mnn@assays$RNA@scale.data = lympho@assays$RNA@scale.data
saveRDS(lympho_mnn, "LymphoidB_D0_4h_MNN.rds")