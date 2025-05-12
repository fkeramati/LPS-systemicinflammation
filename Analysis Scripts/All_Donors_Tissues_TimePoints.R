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

#####################
##### IMPORTANT #####
# In order to run the following code you need to download the single-cell RNA sequencing results from GEO depository and put them in appropriate order.
# I have put files in the following order: 
#BM
#---Sbj1
#------D0
#------4h
#------D7
#---Sbj2
#------D0
#------4h
#------D7
#---Sbj3
#------D0
#------4h
#------D7
#PBMC
#---Sbj1
#------D0
#------4h
#------D7
#---Sbj3
#------D0
#------4h
#------D7
#---Sbj4
#------D0
#------4h
#------D7
# If you change the ordering or organization of the files, you should change the code accordingly.
####################

########################
### Helper functions ###
########################
# Helper function to Normalize & Scale & Reduce Dimentionality of data & Cluster
normVis = function(obj){
  # Normalize and Scale
  obj = NormalizeData(obj, verbose = F) %>% FindVariableFeatures(verbose = F) %>% 
    ScaleData(vars.to.regress = c("nCount_RNA"), verbose = F)
  # PCA, UMAP dimension reduction and clustering
  set.seed(1234)
  obj = RunPCA(obj, verbose = F) %>% suppressWarnings(RunUMAP(dims = 1:15, verbose = F)) %>% 
    FindNeighbors(dims = 1:15, k.param = 30, verbose = F) %>% FindClusters(resolution = 1, verbose = F)
  return(obj)
}

# Helper function to remove low quality cells, 
# automatically assign cell types to clusters, and
# remove potential doublets
removeDoublets = function (sample, path, cluster = 4){
  print(paste("Processing:", sample))
  # Read in data and generate Seurat object
  obj.data = Read10X(data.dir = path)
  obj = CreateSeuratObject(counts = obj.data, project = sample, min.cells = 3, min.features = 50)
  obj = RenameCells(obj, add.cell.id = sample)
  # Remove really low/high gene expressing and very low quality cells
  obj[["percent.mt"]] = PercentageFeatureSet(obj, pattern = "MT-")
  obj = subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 30 & nCount_RNA < 50000)
  # Identify Doublets
  set.seed(1234)
  sce = suppressMessages(scDblFinder(GetAssayData(obj, slot="counts"), clusters = cluster, verbose = F))
  obj$scDblFinder.class = sce$scDblFinder.class
  # Use sctype to automatically assign cell types to clusters
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  # Remove mito and rb genes
  mt.genes = rownames(obj)[grep("^MT-", rownames(obj))]
  rb.genes = rownames(obj)[grep("^RP[LS]", rownames(obj))]
  genes = c(mt.genes, rb.genes)
  genes.keep = setdiff(rownames(obj), genes)
  obj = subset(obj, features = genes.keep)
  rm(genes, genes.keep, mt.genes, rb.genes)
  # Remove non-coding genes
  # For this step you need to have access to a gtf file containing the annotation of genes.
  # You can find this file in the github depository of our paper.
  genes = read.table("genes.gtf", sep = "\t")
  genes = genes[genes$V3 == "gene", ]
  df_genes = data.frame(gene_id = unlist(strsplit(genes$V9, ";| "))[c(F,T,F,F,F,F,F,F,F,F,F,F,F,F)],
                        gene_name = unlist(strsplit(genes$V9, ";| "))[c(F,F,F,F,F,F,F,T,F,F,F,F,F,F)],
                        gene_biotype = unlist(strsplit(genes$V9, ";| "))[c(F,F,F,F,F,F,F,F,F,F,F,F,F,T)])
  df_genes = df_genes[df_genes$gene_biotype == "protein_coding",]
  obj = subset(obj, features = df_genes$gene_name)
  rm(genes, df_genes)
  obj = normVis(obj)
  tissue = "Immune system" 
  # prepare gene sets
  # For this step you need simplified excel file that contains marker genes for major immune cell types.
  # You can find this file in the github depository of our paper.
  gs_list = gene_sets_prepare("Simplified_ScTypeDB.xlsx", tissue)
  # get cell-type by cell matrix
  es.max = sctype_score(scRNAseqData = obj[["RNA"]]@scale.data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
  # merge by cluster
  cL_resutls = do.call("rbind", lapply(unique(obj@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(obj@meta.data[obj@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(obj@meta.data$seurat_clusters==cl)), 10)}))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
  sctype_scores = unique(sctype_scores)
  sctype_scores = sctype_scores[order(sctype_scores$cluster),]
  Tissue = unlist(strsplit(colnames(obj)[1], "_"))[1]
  tnk_cls = which(sctype_scores$type %in% c("CD4 T", "CD8 T", "NK")) - 1
  mye_cls = which(sctype_scores$type %in% c("Monocytes","HSC","GMP","DC")) - 1
  bpl_cls = which(sctype_scores$type %in% c("B", "Plasma", "PDC")) - 1
  tnk = subset(obj, idents = tnk_cls)
  tmp.data = GetAssayData(tnk, slot = "counts")
  tnk = CreateSeuratObject(counts = tmp.data)
  tnk = normVis(tnk)
  tnk$scDblFinder.class = obj$scDblFinder.class
  tbl = as.data.frame.matrix(table(tnk$seurat_clusters, tnk$scDblFinder.class))
  tbl = tbl[tbl$singlet < tbl$doublet,]
  if (nrow(tbl) > 0){
    tnk_cells = union(colnames(subset(tnk, idents = rownames(tbl))), colnames(subset(tnk, subset = scDblFinder.class == "doublet")))
  } else {
    tnk_cells = colnames(subset(tnk, subset = scDblFinder.class == "doublet"))
  }
  mye = subset(obj, idents = mye_cls)
  tmp.data = GetAssayData(mye, slot = "counts")
  mye = CreateSeuratObject(counts = tmp.data)
  mye = normVis(mye)
  mye$scDblFinder.class = obj$scDblFinder.class
  tbl = as.data.frame.matrix(table(mye$seurat_clusters, mye$scDblFinder.class))
  tbl = tbl[tbl$singlet < tbl$doublet,]
  if (nrow(tbl) > 0){
    mye_cells = union(colnames(subset(mye, idents = rownames(tbl))), colnames(subset(mye, subset = scDblFinder.class == "doublet")))
  } else {
    mye_cells = colnames(subset(mye, subset = scDblFinder.class == "doublet"))
  }
  bpl = subset(obj, idents = bpl_cls)
  tmp.data = GetAssayData(bpl, slot = "counts")
  bpl = CreateSeuratObject(counts = tmp.data)
  bpl = normVis(bpl)
  bpl$scDblFinder.class = obj$scDblFinder.class
  tbl = as.data.frame.matrix(table(bpl$seurat_clusters, bpl$scDblFinder.class))
  tbl = tbl[tbl$singlet < tbl$doublet,]
  if (nrow(tbl) > 0){
    bpl_cells = union(colnames(subset(bpl, idents = rownames(tbl))), colnames(subset(bpl, subset = scDblFinder.class == "doublet")))
  } else {
    bpl_cells = colnames(subset(bpl, subset = scDblFinder.class == "doublet"))
  }
  ery_cells = character()
  if (Tissue == "BM"){
    ery_cls = which(sctype_scores$type %in% c("Ery","Platelets")) - 1
    ery = subset(obj, idents = ery_cls)
    tmp.data = GetAssayData(ery, slot = "counts")
    ery = CreateSeuratObject(counts = tmp.data)
    ery = normVis(ery)
    ery$scDblFinder.class = obj$scDblFinder.class
    tbl = as.data.frame.matrix(table(ery$seurat_clusters, ery$scDblFinder.class))
    tbl = tbl[tbl$singlet < tbl$doublet,]
    if (nrow(tbl) > 0){
      ery_cells = union(colnames(subset(ery, idents = rownames(tbl))), colnames(subset(ery, subset = scDblFinder.class == "doublet")))
    } else {
      ery_cells = colnames(subset(ery, subset = scDblFinder.class == "doublet"))
    }
  }
  cells = union(union(union(tnk_cells,mye_cells),bpl_cells),ery_cells)
  return(cells)
}

# Do the doublet removal for each sample
# This is the first step of doublet removal and we remove doublets identified by scDblFinder.
df = data.frame(Sample = c("BM_Sbj1_D0","BM_Sbj1_4h","BM_Sbj1_D7","BM_Sbj2_D0","BM_Sbj2_4h","BM_Sbj2_D7",
                           "BM_Sbj3_D0","BM_Sbj3_4h","BM_Sbj3_D7","PBMC_Sbj4_D0","PBMC_Sbj4_4h","PBMC_Sbj4_D7",
                           "PBMC_Sbj1_D0","PBMC_Sbj1_4h","PBMC_Sbj1_D7","PBMC_Sbj3_D0","PBMC_Sbj3_4h","PBMC_Sbj3_D7"),
                Path = c("./BM/Sbj1/D0/filtered_feature_bc_matrix/",
                         "./BM/Sbj1/4h/filtered_feature_bc_matrix/",
                         "./BM/Sbj1/D7/filtered_feature_bc_matrix/",
                         "./BM/Sbj2/D0/filtered_feature_bc_matrix/",
                         "./BM/Sbj2/4h/filtered_feature_bc_matrix/",
                         "./BM/Sbj2/D7/filtered_feature_bc_matrix/",
                         "./BM/Sbj3/D0/filtered_feature_bc_matrix/",
                         "./BM/Sbj3/4h/filtered_feature_bc_matrix/",
                         "./BM/Sbj3/D7/filtered_feature_bc_matrix/",
                         "./PBMC/Sbj1/D0/filtered_feature_bc_matrix/",
                         "./PBMC/Sbj1/4h/filtered_feature_bc_matrix/",
                         "./PBMC/Sbj1/D7/filtered_feature_bc_matrix/",
                         "./PBMC/Sbj3/D0/filtered_feature_bc_matrix/",
                         "./PBMC/Sbj3/4h/filtered_feature_bc_matrix/",
                         "./PBMC/Sbj3/D7/filtered_feature_bc_matrix/",
                         "./PBMC/Sbj4/D0/filtered_feature_bc_matrix/",
                         "./PBMC/Sbj4/4h/filtered_feature_bc_matrix/",
                         "./PBMC/Sbj4/D7/filtered_feature_bc_matrix/"))
cells = list()
for (i in 1:nrow(df)){
  cells[[df$Sample[i]]] = removeDoublets(sample = df$Sample[i], path = df$Path[i])
}
saveRDS(cells, "Cells_All.rds")
rm(i, removeDoublets)

# To refine more and remove the doublets more precisely, at this step we only check Day 0 (baseline) samples
# and visualize and cluster them by cell type and identify potential doublets within each lineage and further
# remove them.
###########################
df = data.frame(Sample = c("BM_Sbj1_D0","BM_Sbj2_D0","BM_Sbj3_D0","PBMC_Sbj4_D0","PBMC_Sbj1_D0","PBMC_Sbj3_D0"),
                Path = c("./BM/Sbj1/D0/filtered_feature_bc_matrix/",
                         "./BM/Sbj2/D0/filtered_feature_bc_matrix/",
                         "./BM/Sbj3/D0/filtered_feature_bc_matrix/",
                         "./PBMC/Sbj4/D0/filtered_feature_bc_matrix/",
                         "./PBMC/Sbj1/D0/filtered_feature_bc_matrix/",
                         "./PBMC/Sbj3/D0/filtered_feature_bc_matrix/"))

obj = list()
for (i in 1:nrow(df)){
  obj.data = Read10X(data.dir = df$Path[i])
  obj[[i]] = CreateSeuratObject(counts = obj.data, project = df$Sample[i], min.cells = 3, min.features = 50)
  obj[[i]] = RenameCells(obj[[i]], add.cell.id = df$Sample[i])
  if (i == 6){
    endo = merge(obj[[1]], c(obj[[2]],obj[[3]],obj[[4]],obj[[5]],obj[[6]]))
    rm(obj,obj.data)
  }
}
rm(i,df)

# Keep only identified high quality singlets
cells = readRDS("Cells_All.rds")
endo = subset(endo, cells = unlist(cells), invert = T)

# Remove low quality cells
endo[["percent.mt"]] = PercentageFeatureSet(endo, pattern = "MT-")
endo[["percent.rb"]] = PercentageFeatureSet(endo, pattern = "^RP[LS]")
endo = subset(endo, subset = nFeature_RNA > 350 & nFeature_RNA < 7500 & percent.mt < 15 & nCount_RNA < 50000)

# Remove mito and rb genes
mt.genes = rownames(endo)[grep("^MT-", rownames(endo))]
rb.genes = rownames(endo)[grep("^RP[LS]", rownames(endo))]
genes = c(mt.genes, rb.genes)
genes.keep = setdiff(rownames(endo), genes)
endo = subset(endo, features = genes.keep)
rm(genes, genes.keep, mt.genes, rb.genes)

# Remove non-coding genes
genes = read.table("genes.gtf", sep = "\t")
genes = genes[genes$V3 == "gene", ]
df_genes = data.frame(gene_id = unlist(strsplit(genes$V9, ";| "))[c(F,T,F,F,F,F,F,F,F,F,F,F,F,F)],
                      gene_name = unlist(strsplit(genes$V9, ";| "))[c(F,F,F,F,F,F,F,T,F,F,F,F,F,F)],
                      gene_biotype = unlist(strsplit(genes$V9, ";| "))[c(F,F,F,F,F,F,F,F,F,F,F,F,F,T)])
df_genes = df_genes[df_genes$gene_biotype == "protein_coding",]
endo = subset(endo, features = df_genes$gene_name)
rm(genes, df_genes)

endo = normVis(endo)

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
tissue = "Immune system" 
# prepare gene sets
gs_list = gene_sets_prepare("Simplified_ScTypeDB.xlsx", tissue)
# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = endo[["RNA"]]@scale.data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(endo@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(endo@meta.data[endo@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(endo@meta.data$seurat_clusters==cl)), 10)}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
sctype_scores = unique(sctype_scores)
sctype_scores = sctype_scores[order(sctype_scores$cluster),]

tnk_cls = which(sctype_scores$type %in% c("CD4 T", "CD8 T", "NK")) - 1
mye_cls = which(sctype_scores$type %in% c("Monocytes","HSC","GMP","DC")) - 1
bpl_cls = which(sctype_scores$type %in% c("B", "Plasma", "PDC")) - 1
ery_cls = which(sctype_scores$type %in% c("Ery","Platelets")) - 1
tnk = subset(endo, idents = tnk_cls)
tmp.data = GetAssayData(tnk, slot = "counts")
tnk = CreateSeuratObject(counts = tmp.data)
tnk = normVis(tnk)
mye = subset(endo, idents = mye_cls)
tmp.data = GetAssayData(mye, slot = "counts")
mye = CreateSeuratObject(counts = tmp.data)
mye = normVis(mye)
bpl = subset(endo, idents = bpl_cls)
tmp.data = GetAssayData(bpl, slot = "counts")
bpl = CreateSeuratObject(counts = tmp.data)
bpl = normVis(bpl)
ery = subset(endo, idents = ery_cls)
tmp.data = GetAssayData(ery, slot = "counts")
ery = CreateSeuratObject(counts = tmp.data)
ery = normVis(ery)

# Clusters 12 (granulocytes), 13 (doublet with B), 14 (doublet with T) in mye object should be removed 
# Clusters 10 (doublet with NK), 12 (doublet with B), 11 (doublet with Monos) in ery object should be removed
# Platelets also attach to T cells, so I will remove cells expressing CD3E
# Cluster 14 (doublet with T) in bpl object should be removed
cells1 = colnames(subset(mye, idents = c(12,13,14)))
cells2 = colnames(subset(ery, idents = c(10,11,12)))
cells2 = c(cells2, colnames(subset(ery, subset = CD3E > 1, slot = "counts")))
cells3 = colnames(subset(bpl, idents = c(14)))
cells_all = c(cells1, cells2, cells3)
cells[["Other"]] = cells_all
saveRDS(cells, "Cells_All.rds")

# To refine more and remove the doublets more precisely, at this step we only check 4 hours (inflammatory) samples
# and visualize and cluster them by cell type and identify potential doublets within each lineage and further
# remove them.
############################
df = data.frame(Sample = c("BM_Sbj1_4h","BM_Sbj2_4h","BM_Sbj3_4h","PBMC_Sbj4_4h","PBMC_Sbj1_4h","PBMC_Sbj3_4h"),
                Path = c("./BM/Sbj1/4h/filtered_feature_bc_matrix/",
                         "./BM/Sbj2/4h/filtered_feature_bc_matrix/",
                         "./BM/Sbj3/4h/filtered_feature_bc_matrix/",
                         "./PBMC/Sbj4/4h/filtered_feature_bc_matrix/",
                         "./PBMC/Sbj1/4h/filtered_feature_bc_matrix/",
                         "./PBMC/Sbj3/4h/filtered_feature_bc_matrix/"))

obj = list()
for (i in 1:nrow(df)){
  obj.data = Read10X(data.dir = df$Path[i])
  obj[[i]] = CreateSeuratObject(counts = obj.data, project = df$Sample[i], min.cells = 3, min.features = 50)
  obj[[i]] = RenameCells(obj[[i]], add.cell.id = df$Sample[i])
  if (i == 6){
    endo = merge(obj[[1]], c(obj[[2]],obj[[3]],obj[[4]],obj[[5]],obj[[6]]))
    rm(obj,obj.data)
  }
}
rm(i,df)

# Keep only identified high quality singlets
cells = readRDS("Cells_All.rds")
endo = subset(endo, cells = unlist(cells), invert = T)

# Remove low quality cells
endo[["percent.mt"]] = PercentageFeatureSet(endo, pattern = "MT-")
endo[["percent.rb"]] = PercentageFeatureSet(endo, pattern = "^RP[LS]")
endo = subset(endo, subset = nFeature_RNA > 350 & nFeature_RNA < 7500 & percent.mt < 15 & nCount_RNA < 50000)

# Remove mito and rb genes
mt.genes = rownames(endo)[grep("^MT-", rownames(endo))]
rb.genes = rownames(endo)[grep("^RP[LS]", rownames(endo))]
genes = c(mt.genes, rb.genes)
genes.keep = setdiff(rownames(endo), genes)
endo = subset(endo, features = genes.keep)
rm(genes, genes.keep, mt.genes, rb.genes)

# Remove non-coding genes
genes = read.table("genes.gtf", sep = "\t")
genes = genes[genes$V3 == "gene", ]
df_genes = data.frame(gene_id = unlist(strsplit(genes$V9, ";| "))[c(F,T,F,F,F,F,F,F,F,F,F,F,F,F)],
                      gene_name = unlist(strsplit(genes$V9, ";| "))[c(F,F,F,F,F,F,F,T,F,F,F,F,F,F)],
                      gene_biotype = unlist(strsplit(genes$V9, ";| "))[c(F,F,F,F,F,F,F,F,F,F,F,F,F,T)])
df_genes = df_genes[df_genes$gene_biotype == "protein_coding",]
endo = subset(endo, features = df_genes$gene_name)
rm(genes, df_genes)

endo = normVis(endo)

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
tissue = "Immune system" 
# prepare gene sets
gs_list = gene_sets_prepare("Simplified_ScTypeDB.xlsx", tissue)
# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = endo[["RNA"]]@scale.data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(endo@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(endo@meta.data[endo@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(endo@meta.data$seurat_clusters==cl)), 10)}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
sctype_scores = unique(sctype_scores)
sctype_scores = sctype_scores[order(sctype_scores$cluster),]

tnk_cls = which(sctype_scores$type %in% c("CD4 T", "CD8 T", "NK")) - 1
mye_cls = which(sctype_scores$type %in% c("Monocytes","HSC","GMP","DC")) - 1
bpl_cls = which(sctype_scores$type %in% c("B", "Plasma", "PDC")) - 1
ery_cls = which(sctype_scores$type %in% c("Ery","Platelets")) - 1
tnk = subset(endo, idents = tnk_cls)
tmp.data = GetAssayData(tnk, slot = "counts")
tnk = CreateSeuratObject(counts = tmp.data)
tnk = normVis(tnk)
mye = subset(endo, idents = mye_cls)
tmp.data = GetAssayData(mye, slot = "counts")
mye = CreateSeuratObject(counts = tmp.data)
mye = normVis(mye)
bpl = subset(endo, idents = bpl_cls)
tmp.data = GetAssayData(bpl, slot = "counts")
bpl = CreateSeuratObject(counts = tmp.data)
bpl = normVis(bpl)
ery = subset(endo, idents = ery_cls)
tmp.data = GetAssayData(ery, slot = "counts")
ery = CreateSeuratObject(counts = tmp.data)
ery = normVis(ery)

# Clusters 10 (granulocytes) in mye object should be removed 
# Clusters 11 (doublet with B), 6 (doublet with T) in ery object should be removed
# Cluster 15 (doublet with B) in tnk object should be removed
cells1 = colnames(subset(mye, idents = c(10)))
cells2 = colnames(subset(ery, idents = c(6, 11)))
cells3 = colnames(subset(tnk, idents = c(15)))
cells_all = c(cells1, cells2)
cells[["Other"]] = c(cells[["Other"]], cells_all)
saveRDS(cells, "Cells_All.rds")



# To refine more and remove the doublets more precisely, at this step we only check Day 7 (long-term) samples
# and visualize and cluster them by cell type and identify potential doublets within each lineage and further
# remove them.
###########################
###########################
df = data.frame(Sample = c("BM_Sbj1_D7","BM_Sbj2_D7","BM_Sbj3_D7","PBMC_Sbj4_D7","PBMC_Sbj1_D7","PBMC_Sbj3_D7"),
                Path = c("./BM/Sbj1/D7/filtered_feature_bc_matrix/",
                         "./BM/Sbj2/D7/filtered_feature_bc_matrix/",
                         "./BM/Sbj3/D7/filtered_feature_bc_matrix/",
                         "./PBMC/Sbj4/D7/filtered_feature_bc_matrix/",
                         "./PBMC/Sbj1/D7/filtered_feature_bc_matrix/",
                         "./PBMC/Sbj3/D7/filtered_feature_bc_matrix/"))

obj = list()
for (i in 1:nrow(df)){
  obj.data = Read10X(data.dir = df$Path[i])
  obj[[i]] = CreateSeuratObject(counts = obj.data, project = df$Sample[i], min.cells = 3, min.features = 50)
  obj[[i]] = RenameCells(obj[[i]], add.cell.id = df$Sample[i])
  if (i == 6){
    endo = merge(obj[[1]], c(obj[[2]],obj[[3]],obj[[4]],obj[[5]],obj[[6]]))
    rm(obj,obj.data)
  }
}
rm(i,df)

# Keep only identified high quality singlets
cells = readRDS("Cells_All.rds")
endo = subset(endo, cells = unlist(cells), invert = T)

# Remove low quality cells
endo[["percent.mt"]] = PercentageFeatureSet(endo, pattern = "MT-")
endo[["percent.rb"]] = PercentageFeatureSet(endo, pattern = "^RP[LS]")
endo = subset(endo, subset = nFeature_RNA > 350 & nFeature_RNA < 7500 & percent.mt < 15 & nCount_RNA < 50000)

# Remove mito and rb genes
mt.genes = rownames(endo)[grep("^MT-", rownames(endo))]
rb.genes = rownames(endo)[grep("^RP[LS]", rownames(endo))]
genes = c(mt.genes, rb.genes)
genes.keep = setdiff(rownames(endo), genes)
endo = subset(endo, features = genes.keep)
rm(genes, genes.keep, mt.genes, rb.genes)

# Remove non-coding genes
genes = read.table("genes.gtf", sep = "\t")
genes = genes[genes$V3 == "gene", ]
df_genes = data.frame(gene_id = unlist(strsplit(genes$V9, ";| "))[c(F,T,F,F,F,F,F,F,F,F,F,F,F,F)],
                      gene_name = unlist(strsplit(genes$V9, ";| "))[c(F,F,F,F,F,F,F,T,F,F,F,F,F,F)],
                      gene_biotype = unlist(strsplit(genes$V9, ";| "))[c(F,F,F,F,F,F,F,F,F,F,F,F,F,T)])
df_genes = df_genes[df_genes$gene_biotype == "protein_coding",]
endo = subset(endo, features = df_genes$gene_name)
rm(genes, df_genes)

endo = normVis(endo)

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
tissue = "Immune system" 
# prepare gene sets
gs_list = gene_sets_prepare("Simplified_ScTypeDB.xlsx", tissue)
# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = endo[["RNA"]]@scale.data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(endo@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(endo@meta.data[endo@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(endo@meta.data$seurat_clusters==cl)), 10)}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
sctype_scores = unique(sctype_scores)
sctype_scores = sctype_scores[order(sctype_scores$cluster),]

tnk_cls = which(sctype_scores$type %in% c("CD4 T", "CD8 T", "NK")) - 1
mye_cls = which(sctype_scores$type %in% c("Monocytes","HSC","GMP","DC")) - 1
bpl_cls = which(sctype_scores$type %in% c("B", "Plasma", "PDC")) - 1
ery_cls = which(sctype_scores$type %in% c("Ery","Platelets")) - 1
tnk = subset(endo, idents = tnk_cls)
tmp.data = GetAssayData(tnk, slot = "counts")
tnk = CreateSeuratObject(counts = tmp.data)
tnk = normVis(tnk)
mye = subset(endo, idents = mye_cls)
tmp.data = GetAssayData(mye, slot = "counts")
mye = CreateSeuratObject(counts = tmp.data)
mye = normVis(mye)
bpl = subset(endo, idents = bpl_cls)
tmp.data = GetAssayData(bpl, slot = "counts")
bpl = CreateSeuratObject(counts = tmp.data)
bpl = normVis(bpl)
ery = subset(endo, idents = ery_cls)
tmp.data = GetAssayData(ery, slot = "counts")
ery = CreateSeuratObject(counts = tmp.data)
ery = normVis(ery)

# Clusters 6+8 (doublet with T), 12 (doublet with Monos), 11 (doublet with B) in ery object should be removed
# Clusters 12 (doublet with B), 11 (granulocytes) in mye object should be removed
cells1 = colnames(subset(ery, idents = c(6,8,11,12)))
cells2 = colnames(subset(mye, idents = c(11,12)))
cells_all = c(cells1, cells2)
cells[["Other"]] = c(cells[["Other"]], cells_all)
saveRDS(cells, "Cells_All.rds")



# Now that we have checked all samples both one by one and samples from the same time points, we will
# generate a merged object out of all sample (from all time points and tissues). Do the final check of
# the quality and save this object for downstream analysis.
############################
############################
rm(list = ls())
df = data.frame(Sample = c("BM_Sbj1_D0","BM_Sbj1_4h","BM_Sbj1_D7","BM_Sbj2_D0","BM_Sbj2_4h","BM_Sbj2_D7",
                           "BM_Sbj3_D0","BM_Sbj3_4h","BM_Sbj3_D7","PBMC_Sbj4_D0","PBMC_Sbj4_4h","PBMC_Sbj4_D7",
                           "PBMC_Sbj1_D0","PBMC_Sbj1_4h","PBMC_Sbj1_D7","PBMC_Sbj3_D0","PBMC_Sbj3_4h","PBMC_Sbj3_D7"),
                Path = c("./BM/Sbj1/D0/filtered_feature_bc_matrix/",
                         "./BM/Sbj1/4h/filtered_feature_bc_matrix/",
                         "./BM/Sbj1/D7/filtered_feature_bc_matrix/",
                         "./BM/Sbj2/D0/filtered_feature_bc_matrix/",
                         "./BM/Sbj2/4h/filtered_feature_bc_matrix/",
                         "./BM/Sbj2/D7/filtered_feature_bc_matrix/",
                         "./BM/Sbj3/D0/filtered_feature_bc_matrix/",
                         "./BM/Sbj3/4h/filtered_feature_bc_matrix/",
                         "./BM/Sbj3/D7/filtered_feature_bc_matrix/",
                         "./PBMC/Sbj4/D0/filtered_feature_bc_matrix/",
                         "./PBMC/Sbj4/4h/filtered_feature_bc_matrix/",
                         "./PBMC/Sbj4/D7/filtered_feature_bc_matrix/",
                         "./PBMC/Sbj1/D0/filtered_feature_bc_matrix/",
                         "./PBMC/Sbj1/4h/filtered_feature_bc_matrix/",
                         "./PBMC/Sbj1/D7/filtered_feature_bc_matrix/",
                         "./PBMC/Sbj3/D0/filtered_feature_bc_matrix/",
                         "./PBMC/Sbj3/4h/filtered_feature_bc_matrix/",
                         "./PBMC/Sbj3/D7/filtered_feature_bc_matrix/"))
obj = list()
for (i in 1:nrow(df)){
  obj.data = Read10X(data.dir = df$Path[i])
  obj[[i]] = CreateSeuratObject(counts = obj.data, project = df$Sample[i], min.cells = 3, min.features = 50)
  obj[[i]] = RenameCells(obj[[i]], add.cell.id = df$Sample[i])
  if (i == 18){
    endo = merge(obj[[1]], c(obj[[2]],obj[[3]],obj[[4]],obj[[5]],obj[[6]],obj[[7]],obj[[8]],obj[[9]],
                             obj[[10]],obj[[11]],obj[[12]],obj[[13]],obj[[14]],obj[[15]],obj[[16]],obj[[17]],obj[[18]]))
    rm(obj,obj.data)
  }
}

# Keep only identified high quality singlets
cells = readRDS("Cells_All.rds")
endo = subset(endo, cells = unlist(cells), invert = T)

# Remove low quality cells
endo[["percent.mt"]] = PercentageFeatureSet(endo, pattern = "MT-")
endo[["percent.rb"]] = PercentageFeatureSet(endo, pattern = "^RP[LS]")
endo = subset(endo, subset = nFeature_RNA > 350 & nFeature_RNA < 7500 & percent.mt < 20 & nCount_RNA < 50000)

# Remove non-coding genes
genes = read.table("genes.gtf", sep = "\t")
genes = genes[genes$V3 == "gene", ]
df_genes = data.frame(gene_id = unlist(strsplit(genes$V9, ";| "))[c(F,T,F,F,F,F,F,F,F,F,F,F,F,F)],
                      gene_name = unlist(strsplit(genes$V9, ";| "))[c(F,F,F,F,F,F,F,T,F,F,F,F,F,F)],
                      gene_biotype = unlist(strsplit(genes$V9, ";| "))[c(F,F,F,F,F,F,F,F,F,F,F,F,F,T)])
df_genes = df_genes[df_genes$gene_biotype == "protein_coding",]
endo = subset(endo, features = df_genes$gene_name)
rm(genes, df_genes)

# Normalize and Scale
endo = NormalizeData(endo, normalization.method = "LogNormalize", scale.factor = 10000)
endo = FindVariableFeatures(endo, selection.method = "vst", nfeatures = 5000)
endo = ScaleData(endo)

# PCA, UMAP visualization and clustering
set.seed(1234)
endo = RunPCA(endo, ndims.print = 1:30, verbose = F)
ElbowPlot(endo, ndims = 50)
endo = RunUMAP(endo, dims = 1:15)
endo = FindNeighbors(endo, dims = 1:15, k.param = 30)
endo = FindClusters(endo, resolution = 2.5)
endo$Tissue = unlist(strsplit(endo$orig.ident, split = "_"))[c(T,F,F)]
endo$Donor = unlist(strsplit(endo$orig.ident, split = "_"))[c(F,T,F)]
endo$TimePoint = unlist(strsplit(endo$orig.ident, split = "_"))[c(F,F,T)]
DimPlot(endo, label = TRUE, pt.size = 0.01)
DimPlot(endo, label = F, pt.size = 0.01, group.by = "Donor", shuffle = T)
DimPlot(endo, label = F, pt.size = 0.01, group.by = "TimePoint", shuffle = T)
DimPlot(endo, label = F, pt.size = 0.01, group.by = "Tissue", shuffle = T)

# Remove mito and rb genes
mt.genes = rownames(endo)[grep("^MT-", rownames(endo))]
rb.genes = rownames(endo)[grep("^RP[LS]", rownames(endo))]
genes = c(mt.genes, rb.genes)
genes.keep = setdiff(rownames(endo), genes)
endo = subset(endo, features = genes.keep)
rm(genes, genes.keep, mt.genes, rb.genes)

# There are some clusters that are specifically enriched in high MT% cells
FeaturePlot(endo, order = T, min.cutoff = "10", max.cutoff = "q99", raster = F, label = T, "percent.mt")
endo = FindSubCluster(endo, cluster = 15, graph.name = "RNA_snn", resolution = 0.1)
Idents(endo) = endo$sub.cluster
# Remove those clusters
endo = subset(endo, idents = c("15_1",28,50,38,51,30,47,25,22,23,18,45), invert = T)
endo = subset(endo, subset = percent.mt < 10)

# Normalize and Scale
endo = NormalizeData(endo, normalization.method = "LogNormalize", scale.factor = 10000)
endo = FindVariableFeatures(endo, selection.method = "vst", nfeatures = 5000)
endo = ScaleData(endo)

# PCA, UMAP visualization and clustering
set.seed(1234)
endo = RunPCA(endo, ndims.print = 1:30, verbose = F)
ElbowPlot(endo, ndims = 50)
endo = RunUMAP(endo, dims = 1:15)
endo = FindNeighbors(endo, dims = 1:15, k.param = 30)
endo = FindClusters(endo, resolution = 2.5)
endo$Tissue = unlist(strsplit(endo$orig.ident, split = "_"))[c(T,F,F)]
endo$Donor = unlist(strsplit(endo$orig.ident, split = "_"))[c(F,T,F)]
endo$TimePoint = unlist(strsplit(endo$orig.ident, split = "_"))[c(F,F,T)]
DimPlot(endo, label = TRUE, pt.size = 0.01)
DimPlot(endo, label = F, pt.size = 0.01, group.by = "Donor", shuffle = T)
DimPlot(endo, label = F, pt.size = 0.01, group.by = "TimePoint", shuffle = T)
DimPlot(endo, label = F, pt.size = 0.01, group.by = "Tissue", shuffle = T)

# save the raw (not batch corrected file for downstream analysis)
saveRDS(endo, "All_Donor_TimePoint_Tissue_logNorm_NoBatchCorrect.rds")
cells = colnames(endo)
# save the list of all qualified (high quality and singlets) cells.
saveRDS(cells, "Qualified_Cells_All.rds")



########################
### Batch correction ###
########################
# Correct for Sequencing Batch Differences
endo_mnn = RunFastMNN(object.list = SplitObject(endo, split.by = "Tissue"), features = 2000, verbose = T)
endo_mnn = RunUMAP(endo_mnn, dims = 1:30, reduction = "mnn", min.dist = 0.4)
endo_mnn = FindNeighbors(endo_mnn, dims = 1:30, reduction = "mnn", k.param = 20)
endo_mnn = FindClusters(endo_mnn, resolution = 1.5)
DimPlot(endo_mnn, label = T)
myTheme = theme_bw() + theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
                             panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))

DimPlot(endo_mnn, label = F, group.by = "TimePoint", shuffle = T, raster = F, pt.size = 0.001) + scale_color_manual(values=c("#e6ab02", "#33bbee", "#ee3377")) +
  myTheme

# Assign cell types
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system" 
# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)
# get cell-type by cell matrix
VariableFeatures(endo_mnn) = VariableFeatures(endo)
endo_mnn = ScaleData(endo_mnn)
es.max = sctype_score(scRNAseqData = endo_mnn[["RNA"]]@scale.data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(endo_mnn@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(endo_mnn@meta.data[endo_mnn@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(endo_mnn@meta.data$seurat_clusters==cl)), 10)}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
sctype_scores = unique(sctype_scores)
sctype_scores = sctype_scores[order(sctype_scores$cluster),]

endo_mnn = RenameIdents(endo_mnn, "0" = "CD4 Naive", "1" = "CD8 Naive", "2" = "CD4 Naive", "3" = "B Naive", "4" = "Monocyte",
                        "5" = "CD4 Mem", "6" = "CD4 Mem", "7" = "NK", "8" = "Monocyte", "9" = "Inflam Mono", "10" = "B Mem",
                        "11" = "Inflam CD4 Naive", "12" = "CD8 Mem", "13" = "B Naive", "14" = "CD8 Mem", "15" = "CD4 Mem", "16" = "Inflam CD8 Naive",
                        "17" = "CD4 Mem", "18" = "Pro Mono", "19" = "CD8 Mem", "20" = "Non Cla Mono", "21" = "Pre B", "22" = "HSC",
                        "23" = "Late Ery", "24" = "GMP", "25" = "Late Ery", "26" = "Early Ery", "27" = "Inflam NK", "28" = "Pre B",
                        "29" = "pDC", "30" = "cDC", "31" = "Plasma", "32" = "Platelet", "33" = "Pro B", "34" = "Pro B",
                        "35" = "Inflam Mono", "36" = "Late Ery", "37" = "Inflam cDC", "38" = "Late Ery", "39" = "Plasma", "40" = "CD4 Mem")
DimPlot(endo_mnn, label = T)
endo_mnn$CellType = Idents(endo_mnn)
# Save the batch correction Seurat object of all samples (from all time ponts and tissues) for downstream analysis.
saveRDS(endo_mnn, "All_Donor_TimePoint_Tissue_logNorm_MNN.rds")
