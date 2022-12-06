# Data processing

library(parallel)
library(future)
library(Seurat)
library(Matrix)
options(future.globals.maxSize=256*1024^3)
plan(multiprocess)
library(ggplot2)
library(cowplot)
library(magrittr)
library(dplyr)
library(purrr)

# Analysis parameters
mc <- 32        # number of given cores
pcaDims <- 50   # number of PCA components to compute

source("../functions.R")

# Use Seurat built in list of cell cycle related genes
ccGenes <- cc.genes.updated.2019
# Change gene names to mm genes
ccGenes <- lapply(ccGenes, function(x)  {paste0(substring(x, 1, 1), tolower(substring(x, 2, nchar(x))))})

# Read raw data (gene/cell count matrix from cellranger, filtered: use only detected cellular barcodes)
samples <- dir(path = ".", 
               pattern = "filtered_feature_bc_matrix$", 
               full.names = T, 
               recursive = T, 
               include.dirs = T)
samplesRawData <- readRawData(samples = samples)
names(samplesRawData) <- substr(samples, 3, nchar(samples) - 27)

# Split features into RNA, HTOs and ADT counts
rnaCounts <- lapply(samplesRawData, function(x) x[grepl("ENSMUSG", rownames(x)), ])
htosCounts <- lapply(samplesRawData, function(x) x[grepl("HTO", rownames(x)), ])
protCounts <- lapply(samplesRawData, function(x) x[!grepl("ENSMUSG|HTO", rownames(x)), ])

# Get annotation df from samples genes identifiers
annot <- getAnnot(samples)

# Divide cell cycle genes list into markers of G2/M phase and markers of S phase
sGenes <- annot[annot$Gene.name %in% ccGenes$s.genes, "Gene.stable.ID"]
sGenes <- sGenes[!is.na(sGenes)]
g2mGenes <- annot[annot$Gene.name %in% ccGenes$g2m.genes, "Gene.stable.ID"]
g2mGenes <- g2mGenes[!is.na(g2mGenes)]

# Set up Seurat objects
seuObjectsHashtag <- lapply(seq_along(rnaCounts), function(i) {
  CreateSeuratObject(counts = rnaCounts[[i]], 
                     project = names(rnaCounts)[i])
})
seuObjects <- seuObjectsHashtag
names(seuObjects) <- names(rnaCounts)

# Analyze percentage of mitochondrial genes in cells and no. of genes in cells
mitoFeatures <- annot[grep(pattern = "^mt-", annot$Gene.name), ][, 1]
percentMito <- mclapply(seq_along(seuObjects), function(i) {
  mitoReads <- Matrix::colSums(x = GetAssayData(object = seuObjects[[i]], 
                               slot = 'counts')[rownames(seuObjects[[i]]) %in% mitoFeatures, ])
  totalReads <- Matrix::colSums(x = GetAssayData(object = seuObjects[[i]],
                                                 slot = 'counts'))
  pctMito <- mitoReads / totalReads
  pctMito
})

# Create metrics used in QC
# percent.mito
# log10GenesPerUMI
seuObjects <- mclapply(seq_along(seuObjects), function(i) {
  seuObjects[[i]][['percent.mito']] <- percentMito[[i]]
  seuObjects[[i]]
}, mc.cores = mc)
names(seuObjects) <- names(samplesRawData)

seuObjects <- mclapply(seq_along(seuObjects), function(i) {
  seuObjects[[i]][['log10GenesPerUMI']] <- log10(seuObjects[[i]]$nFeature_RNA) / log10(seuObjects[[i]]$nCount_RNA)
  seuObjects[[i]]
}, mc.cores = mc)
names(seuObjects) <- names(samplesRawData)

# Add HTO data as a new assay independent from RNA
for(i in 1:length(seuObjects)) {
  seuObjects[[i]][["HTO"]] <- CreateAssayObject(counts = htosCounts[[i]])
}

# Normalize HTO data using centered log-ratio (CLR) transformation and demultiplex samples
for(i in 1:length(seuObjects)) {
  seuObjects[[i]] <- NormalizeData(seuObjects[[i]], assay = "HTO", normalization.method = "CLR")
  seuObjects[[i]] <- HTODemux(seuObjects[[i]], assay = "HTO", positive.quantile = 0.99)  
}

# Remove negative and doublet cells from the object
for(i in 1:length(seuObjects)) {
  Idents(seuObjects[[i]]) <- "HTO_classification.global"
}

seuObjectsAllCells <- seuObjects
seuObjects <- lapply(seuObjects, function(x) subset(x, idents = "Singlet"))
names(seuObjects) <- names(seuObjects)
names(seuObjects) <- gsub("full_nec", "singlets", names(seuObjects))

# Adjust cutoff threshold for replicate 3 (lower quality of the single-cell library) 
rep3singletsFiltered <- sapply(paste0("HTO", c(1:6)), function(hto) {
  hto.cells <- Cells(seuObjects$mm_rna_prot_gam_sex_rep3_singlets)[seuObjects$mm_rna_prot_gam_sex_rep3_singlets$HTO_maxID == hto]
  hto.data <- seuObjects$mm_rna_prot_gam_sex_rep3_singlets@assays$HTO@data
  hto.cells[(apply(hto.data[rownames(hto.data) != hto, hto.cells], 2, max) < 0.4)]
}) %>% unlist 

# Filter replicate 3
seuObjects$mm_rna_prot_gam_sex_rep3_singlets <- subset(seuObjects$mm_rna_prot_gam_sex_rep3_singlets, cells=rep3singletsFiltered)
names(seuObjects)[3] <- gsub("singlets", "singlets_filtered", names(seuObjects)[3])

# Permissive cell filtering settings
seuObjects <- lapply(seq_along(seuObjects), function(i) {
  seuObjects[[i]] <- subset(x = seuObjects[[i]], subset = (nFeature_RNA > 300) &
                              (percent.mito < 0.1) &
                              (log10GenesPerUMI > 0.8))
  seuObjects[[i]]
})
names(seuObjects) <- names(samplesRawData)

# Gene filtering - keep genes present in 5 or more cells
seuObjects <- lapply(seq_along(seuObjects), function(i) {
  counts <- GetAssayData(object = seuObjects[[i]], slot = "counts")

  nonzero <- counts > 0
  
  # Sums all TRUE values and returns TRUE if more than 5 TRUE values per gene
  keepGenes <- Matrix::rowSums(nonzero) >= 5
  
  # Only keeping those genes expressed in more than 5 cells
  filteredCounts <- counts[keepGenes, ]
  
  # Reassign to Seurat object
  seuObjects[[i]] <- CreateSeuratObject(filteredCounts, meta.data = seuObjects[[i]]@meta.data)
  seuObjects[[i]]
})
names(seuObjects) <- names(samplesRawData)


# Standard integration workflow
seuObjects <- lapply(seq_along(seuObjects), function(i) {
  # Normalize data
  seuObjects[[i]] <- NormalizeData(seuObjects[[i]])
  seuObjects[[i]]
})
names(seuObjects) <- names(samplesRawData)

# Cell cycle scoring for QC purposes
seuObjects <- lapply(seq_along(seuObjects), function(i) {
  # Normalize data
  seuObjects[[i]] <- CellCycleScoring(object = seuObjects[[i]], 
                                     s.features = sGenes, 
                                     g2m.features = g2mGenes, 
                                     set.ident = TRUE)
  seuObjects[[i]]@meta.data$CC.Difference <- (seuObjects[[i]]@meta.data$S.Score - seuObjects[[i]]@meta.data$G2M.Score)
  seuObjects[[i]]
})
names(seuObjects) <- names(samplesRawData)

# Find variable features
seuObjects <- lapply(seq_along(seuObjects), function(i) {
  seuObjects[[i]] <- FindVariableFeatures(seuObjects[[i]], selection.method = "vst", nfeatures = 2000)
  seuObjects[[i]]
})
names(seuObjects) <- names(samplesRawData)

# Add protein (ADT) data to a new assay
for(i in 1:length(seuObjects)) {
  seuObjects[[i]][["ADT"]] <- CreateAssayObject(counts = protCounts[[i]][, colnames(seuObjects[[i]])])
}  

# Select features and anchors used for integration procedure
features <- SelectIntegrationFeatures(object.list = seuObjects)
anchors <- FindIntegrationAnchors(object.list = seuObjects, anchor.features = features)

allFeaturesStandard <- lapply(seuObjects, row.names) %>% Reduce(intersect, .) 

# Integrate data
samplesIntegratedStandard <- IntegrateData(anchorset = anchors,
                                           features.to.integrate = allFeaturesStandard,
                                           verbose = T)

# Set integrated data as default for downstream analysis
DefaultAssay(samplesIntegratedStandard) <- "integrated"

# Scale data
samplesIntegratedStandard <- ScaleData(samplesIntegratedStandard, verbose = FALSE)

# Run PCA 
samplesIntegratedStandard <- RunPCA(object = samplesIntegratedStandard, 
                                       npcs = pcaDims, 
                                       verbose = FALSE)

# Run UMAP on 30 PCs
samplesIntegratedStandard <- RunUMAP(object = samplesIntegratedStandard,
                                        reduction = "pca",
                                        dims=1:30)

# Cluster data (resolution chosen based on results of intermediary analysis)
samplesIntegratedStandard <- FindNeighbors(object = samplesIntegratedStandard,
                                              dims = 1:30)
samplesIntegratedStandard <- FindClusters(object = samplesIntegratedStandard,
                                             resolution = 1.1)

# Assign metadata based on demultiplexing results
samplesIntegratedStandard$day <- "D0"
samplesIntegratedStandard$day[samplesIntegratedStandard$HTO_classification %in% c("HTO3", "HTO4")] <- "D14"
samplesIntegratedStandard$day[samplesIntegratedStandard$HTO_classification %in% c("HTO5", "HTO6")] <- "D21"

samplesIntegratedStandard$sex <- "female"
samplesIntegratedStandard$sex[samplesIntegratedStandard$HTO_classification %in% c("HTO2", "HTO4", "HTO6")] <- "male"

samplesIntegratedStandard$condition <- "ctrl"
samplesIntegratedStandard$condition[samplesIntegratedStandard$HTO_classification %in% c("HTO3", "HTO4", "HTO5", "HTO6")] <- "tumor"

samplesIntegratedStandard$replicate <- "R1"
samplesIntegratedStandard$replicate[samplesIntegratedStandard$orig.ident == "mm_rna_prot_gam_sex_rep2_full_nec"] <- "R2"
samplesIntegratedStandard$replicate[samplesIntegratedStandard$orig.ident == "mm_rna_prot_gam_sex_rep3_full_nec"] <- "R3"

# Preprocess ADT data
seuObjectsADT <- seuObjects 

seuObjectsADT <- lapply(seuObjectsADT, function(x) {
  DefaultAssay(x) <- "ADT"
  x[["RNA"]] <- NULL
  x <- NormalizeData(x, assay = "ADT", normalization.method = "CLR")
  x <- FindVariableFeatures(x, selection.method = "vst", verbose=T)
  x <- ScaleData(x)
  x <- RunPCA(x, npcs = 5, approx = F)
  x
})
names(seuObjectsADT) <- names(seuObjects)

# Perform independent integration of ADT data
repAnchorsADT <- FindIntegrationAnchors(object.list = seuObjectsADT,
                                          assay = rep("ADT", length(seuObjectsADT)),
                                          anchor.features = rownames(seuObjectsADT[[1]]@assays$ADT),
                                          dims = c(1:5),
                                          verbose = T)

allFeaturesStandardADT <- lapply(seuObjectsADT, row.names) %>% Reduce(intersect, .) 

seuObjectIntegratedADT <- IntegrateData(anchorset = repAnchorsADT,
                                           features.to.integrate = allFeaturesStandardADT,
                                           new.assay.name = "integrated_ADT",
                                           dims = c(1:5),
                                           verbose = T)

# Check if cells in the same order in two objects
(seuObjectIntegratedADT@assays$integrated_ADT@data %>% colnames) %>% identical(colnames(samplesIntegratedStandard@assays$integrated@data))

# Copy integrated_ADT data to samplesIntegratedStandard object
samplesIntegratedStandard[["integrated_ADT"]] <- CreateAssayObject(data = seuObjectIntegratedADT@assays$integrated_ADT@data)
samplesIntegratedStandard <- ScaleData(samplesIntegratedStandard, assay = "integrated_ADT", verbose = FALSE)

# Rename protein names in ADT assays to avoid names conflict with integrated_ADT assay
rownames(samplesIntegratedStandard@assays$ADT@counts) <- paste0(rownames(samplesIntegratedStandard@assays$ADT), "-raw")
rownames(samplesIntegratedStandard@assays$ADT@data) <- paste0(rownames(samplesIntegratedStandard@assays$ADT), "-raw")

# Find markers
samplesIntegratedStandardMarkers_1.1 <- FindAllMarkers(object = samplesIntegratedStandard, 
                                                       assay = "RNA",
                                                       only.pos = TRUE, 
                                                       min.pct = 0.50, 
                                                       logfc.threshold = 0.1)



