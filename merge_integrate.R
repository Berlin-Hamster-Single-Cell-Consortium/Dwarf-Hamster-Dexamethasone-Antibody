library(ggplot2)
library(stringr)
library(Seurat)
library(HDF5Array)
library(DelayedArray)
library(DelayedMatrixStats)
require(gridExtra)

data <- Read10X(data.dir = c("Hamster1/DexAb_Sample1/outs/per_sample_outs/U1/count/sample_feature_bc_matrix/", "Hamster2/DexAb_Sample2/outs/per_sample_outs/U1/count/sample_feature_bc_matrix/", "Hamster3/DexAb_Sample3/outs/per_sample_outs/U1/count/sample_feature_bc_matrix/"), gene.column = 2)
Untr_1 = CreateSeuratObject(counts = data$`Gene Expression`, project="Untr_1", min.cells=5, min.features=1000)
Untr_1@meta.data$orig.ident <- "Untr_1"

data <- Read10X(data.dir = c("Hamster1/DexAb_Sample1/outs/per_sample_outs/U2/count/sample_feature_bc_matrix/", "Hamster2/DexAb_Sample2/outs/per_sample_outs/U2/count/sample_feature_bc_matrix/", "Hamster3/DexAb_Sample3/outs/per_sample_outs/U2/count/sample_feature_bc_matrix/"), gene.column = 2)
Untr_2 = CreateSeuratObject(counts = data$`Gene Expression`, project="Untr_2", min.cells=5, min.features=1000)
Untr_2@meta.data$orig.ident <- "Untr_2"

data <- Read10X(data.dir = c("Hamster1/DexAb_Sample1/outs/per_sample_outs/U3/count/sample_feature_bc_matrix/", "Hamster2/DexAb_Sample2/outs/per_sample_outs/U3/count/sample_feature_bc_matrix/", "Hamster3/DexAb_Sample3/outs/per_sample_outs/U3/count/sample_feature_bc_matrix/"), gene.column = 2)
Untr_3 = CreateSeuratObject(counts = data$`Gene Expression`, project="Untr_3", min.cells=5, min.features=1000)
Untr_3@meta.data$orig.ident <- "Untr_3"

data <- Read10X(data.dir = c("Hamster1/DexAb_Sample1/outs/per_sample_outs/D1/count/sample_feature_bc_matrix/", "Hamster2/DexAb_Sample2/outs/per_sample_outs/D1/count/sample_feature_bc_matrix/", "Hamster3/DexAb_Sample3/outs/per_sample_outs/D1/count/sample_feature_bc_matrix/"), gene.column = 2)
dex_1 = CreateSeuratObject(counts = data$`Gene Expression`, project="dex_1", min.cells=5, min.features=1000)
dex_1@meta.data$orig.ident <- "dex_1"

data <- Read10X(data.dir = c("Hamster1/DexAb_Sample1/outs/per_sample_outs/D2/count/sample_feature_bc_matrix/", "Hamster2/DexAb_Sample2/outs/per_sample_outs/D2/count/sample_feature_bc_matrix/", "Hamster3/DexAb_Sample3/outs/per_sample_outs/D2/count/sample_feature_bc_matrix/"), gene.column = 2)
dex_2 = CreateSeuratObject(counts = data$`Gene Expression`, project="dex_2", min.cells=5, min.features=1000)
dex_2@meta.data$orig.ident <- "dex_2"

data <- Read10X(data.dir = c("Hamster1/DexAb_Sample1/outs/per_sample_outs/D3/count/sample_feature_bc_matrix/", "Hamster2/DexAb_Sample2/outs/per_sample_outs/D3/count/sample_feature_bc_matrix/", "Hamster3/DexAb_Sample3/outs/per_sample_outs/D3/count/sample_feature_bc_matrix/"), gene.column = 2)
dex_3 = CreateSeuratObject(counts = data$`Gene Expression`, project="dex_3", min.cells=5, min.features=1000)
dex_3@meta.data$orig.ident <- "dex_3"

data <- Read10X(data.dir = c("Hamster1/DexAb_Sample1/outs/per_sample_outs/AB1/count/sample_feature_bc_matrix/", "Hamster2/DexAb_Sample2/outs/per_sample_outs/AB1/count/sample_feature_bc_matrix/", "Hamster3/DexAb_Sample3/outs/per_sample_outs/AB1/count/sample_feature_bc_matrix/"), gene.column = 2)
ab_1 = CreateSeuratObject(counts = data$`Gene Expression`, project="ab_1", min.cells=5, min.features=1000)
ab_1@meta.data$orig.ident <- "ab_1"

data <- Read10X(data.dir = c("Hamster1/DexAb_Sample1/outs/per_sample_outs/AB2/count/sample_feature_bc_matrix/", "Hamster2/DexAb_Sample2/outs/per_sample_outs/AB2/count/sample_feature_bc_matrix/", "Hamster3/DexAb_Sample3/outs/per_sample_outs/AB2/count/sample_feature_bc_matrix/"), gene.column = 2)
ab_2 = CreateSeuratObject(counts = data$`Gene Expression`, project="ab_2", min.cells=5, min.features=1000)
ab_2@meta.data$orig.ident <- "ab_2"

data <- Read10X(data.dir = c("Hamster1/DexAb_Sample1/outs/per_sample_outs/AB3/count/sample_feature_bc_matrix/", "Hamster2/DexAb_Sample2/outs/per_sample_outs/AB3/count/sample_feature_bc_matrix/", "Hamster3/DexAb_Sample3/outs/per_sample_outs/AB3/count/sample_feature_bc_matrix/"), gene.column = 2)
ab_3 = CreateSeuratObject(counts = data$`Gene Expression`, project="ab_3", min.cells=5, min.features=1000)
ab_3@meta.data$orig.ident <- "ab_3"

data <- Read10X(data.dir = c("Hamster1/DexAb_Sample1/outs/per_sample_outs/ABD1/count/sample_feature_bc_matrix/", "Hamster2/DexAb_Sample2/outs/per_sample_outs/ABD1/count/sample_feature_bc_matrix/", "Hamster3/DexAb_Sample3/outs/per_sample_outs/ABD1/count/sample_feature_bc_matrix/"), gene.column = 2)
dexab_1 = CreateSeuratObject(counts = data$`Gene Expression`, project="dexab_1", min.cells=5, min.features=1000)
dexab_1@meta.data$orig.ident <- "dexab_1"

data <- Read10X(data.dir = c("Hamster1/DexAb_Sample1/outs/per_sample_outs/ABD2/count/sample_feature_bc_matrix/", "Hamster2/DexAb_Sample2/outs/per_sample_outs/ABD2/count/sample_feature_bc_matrix/", "Hamster3/DexAb_Sample3/outs/per_sample_outs/ABD2/count/sample_feature_bc_matrix/"), gene.column = 2)
dexab_2 = CreateSeuratObject(counts = data$`Gene Expression`, project="dexab_2", min.cells=5, min.features=1000)
dexab_2@meta.data$orig.ident <- "dexab_2"

data <- Read10X(data.dir = c("Hamster1/DexAb_Sample1/outs/per_sample_outs/ABD3/count/sample_feature_bc_matrix/", "Hamster2/DexAb_Sample2/outs/per_sample_outs/ABD3/count/sample_feature_bc_matrix/", "Hamster3/DexAb_Sample3/outs/per_sample_outs/ABD3/count/sample_feature_bc_matrix/"), gene.column = 2)
dexab_3 = CreateSeuratObject(counts = data$`Gene Expression`, project="dexab_3", min.cells=5, min.features=1000)
dexab_3@meta.data$orig.ident <- "dexab_3"



hamster_all <- merge(Untr_1, y = c(Untr_2, Untr_3, ab_1, ab_2, ab_3, dex_1, dex_2, dex_3, dexab_1, dexab_2, dexab_3), add.cell.ids = c("Untr_1", "Untr_2", "Untr_3", "ab_1", "ab_2", "ab_3", "dex_1", "dex_2", "dex_3", "dexab_1", "dexab_2", "dexab_3"), project = "DexAb")
saveRDS(hamster_all, "./DexAb_combined.rds")

DefaultAssay(hamster_all) <- "RNA"


## integrate by hamster to remove batch effects

hamster.list <- SplitObject(hamster_all, split.by = "orig.ident")
for (i in 1:length(hamster.list)){
  hamster.list[[i]] <- SCTransform(hamster.list[[i]], verbose =T)
}

hamster.features <- SelectIntegrationFeatures(object.list = hamster.list, nfeatures = 3000)
hamster.list <- PrepSCTIntegration(object.list = hamster.list, anchor.features = hamster.features, 
                                   verbose = T)

hamster.anchors <- FindIntegrationAnchors(object.list = hamster.list, normalization.method = "SCT", 
                                          anchor.features = hamster.features, verbose = T)
hamster.integrated <- IntegrateData(anchorset = hamster.anchors, normalization.method = "SCT", 
                                    verbose = T)

## run dimensional reductions
#   PCA
hamster.integrated<- RunPCA(hamster.integrated, verbose = FALSE)
#   UMAP
hamster.integrated<- RunUMAP(hamster.integrated, dims = 1:30, verbose = FALSE)

hamster.integrated <- FindNeighbors(hamster.integrated, dims = 1:30)
hamster.integrated <- FindClusters(hamster.integrated, resolution = 0.5)


saveRDS(hamster.integrated, "./DexAb_combined_integrated.rds")

