## Running DecontX on each sample, and integration again
### Seurat preprocessing and integration
### Seurat preprocessing and integration
suppressPackageStartupMessages({
    ## basic use
    library(Seurat)
    library(stringr)
    library(dplyr)
    library(patchwork)
    library(ggplot2)
    library(future)
    library(tidydr)
    library(harmony)
    library(scCustomize)
    ## those used for decontamination
    library(decontX)
    library(SingleCellExperiment)
    ## used for ploting the clustering structer
    library(clustree)
    })

## Provide enough memory for preprocessing  
options(future.globals.maxSize = 500*1024^3)

## Set directory
dir = "/oak/stanford/projects/kibr/Reorganizing/Projects/Andi/vascular_atlas/"
setwd(dir)

## Loading previously integrated files
integrated <- readRDS("./Results/Revision_2/02.integrated_object_revision.rds")
ls()
print(integrated)

## Removing those already been predicted as doublets
integrated = subset(integrated, subset = scDblFinder.class == "singlet")
print(integrated)

###################################################################################
###############-- Integration on the doublets removed object --####################
###################################################################################
## Generating UMAP before decontamination
DefaultAssay(integrated) <- "RNA"
## As the data already split by individualID, do not need to split again
# integrated[["RNA"]] <- split(split[["RNA"]], f = split$individualID)

integrated <- SCTransform(integrated,assay = "RNA",new.assay.name = "SCT",vars.to.regress = c("percent.mt"),
        layer = "counts",verbose = F,method = "glmGamPoi") %>% 
        RunPCA(ndims=30,verbose = F) %>% 
        FindNeighbors(dims = 1:30,verbose = F) %>% 
        RunUMAP(dims = 1:30, reduction.name="umap.rna",reduction.key = "rnaUMAP_",verbose = F)

integrated <- IntegrateLayers(
                    object = integrated, method = HarmonyIntegration,
                    orig.reduction = "pca", new.reduction = "harmony",
                    verbose = FALSE
                    )

integrated <- FindNeighbors(integrated, reduction = "harmony", dims = 1:30)
# integrated <- FindClusters(obj, resolution = c(0.5,1,1.5,2))
integrated <- RunUMAP(integrated, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

##################################################################################
###############-- decontX decontamination before annotation --####################
##################################################################################
## Extract the RNA assay counts layer for decontX decontamination
DefaultAssay(integrated) = "RNA"
### Need to JoinLayer again for matrix extraction
integrated = JoinLayers(integrated)
counts <- LayerData(object = integrated, layer = "counts")

## Generating new singlecellexperiment object
sce <- SingleCellExperiment(list(counts = counts),
                            colData = integrated@meta.data)

## running decontX on each sample 
sce <- decontX(sce,batch = sce$sampleID)

print(range(assay(sce, "counts")))
print(range(assay(sce,"decontXcounts")))

##############################################################################################
###############-- Running integration and clustering on denoised matrix --####################
##############################################################################################
## Create new assay on original seurat object
integrated[["decontX"]] <- CreateAssayObject(data = assay(sce, "decontXcounts"))
DefaultAssay(integrated) = "decontX"
integrated[["decontX"]]$counts = integrated[["decontX"]]$data
print(integrated)

## Get the metadata back to the seurat object
integrated$decontX_contamination = sce$decontX_contamination

## Now, need to split the data
DefaultAssay(integrated) <- "decontX"
integrated[["decontX"]] <- split(integrated[["decontX"]], f = integrated$individualID)

## Running integration and clustering
integrated <- SCTransform(integrated,assay = "decontX",new.assay.name = "SCT",vars.to.regress = c("percent.mt"),layer = "counts",verbose = F,method = "glmGamPoi") %>% 
        RunPCA(ndims=30,verbose = F) %>% 
        FindNeighbors(dims = 1:30,verbose = F) %>% 
        RunUMAP(dims = 1:30, reduction.name="umap.decontX",reduction.key = "rnaUMAP_",verbose = F)

integrated <- IntegrateLayers(
  object = integrated, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

integrated <- FindNeighbors(integrated, reduction = "harmony", dims = 1:30)
integrated <- FindClusters(integrated, resolution = c(0.1,0.3,0.5,0.7,1,1.5,2))
integrated <- RunUMAP(integrated, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony.denoised")

##############################################################################################
###############-- Plotting the UMAP of original and denosied object --########################
##############################################################################################
umap_theme <- theme_dr()+theme(panel.grid.major = element_blank(), 
                                            panel.grid.minor = element_blank(),
                                            panel.background = element_blank(), 
                                            axis.line = element_line(colour = "black"))
## Original UMAP
p1 <- DimPlot(integrated, label = T, repel = TRUE, pt.size = 1,
            reduction = "umap.harmony",group.by = "individualID")+umap_theme+ggtitle("DonorID")
p2 <- FeaturePlot(integrated, label = T, repel = TRUE, 
              reduction = "umap.harmony",features = "decontX_contamination")+umap_theme+ggtitle("Contamination_level")
p3 <- FeaturePlot(integrated, repel = TRUE, 
              reduction = "umap.harmony",features = "nCount_RNA")+umap_theme+ggtitle("nCount_RNA")
p4 <- FeaturePlot(integrated, repel = TRUE, 
              reduction = "umap.harmony",features = "nFeature_RNA")+umap_theme+ggtitle("nFeature_RNA")                                          

ggsave(filename = "./Results/Revision_2/Figures/0_integrated_UMAP_before_decontX.pdf",
    patchwork::wrap_plots(p1,p2,p3,p4, ncol = 2),
      scale = 1, width = 13, height = 10)

## New UMAP
p1 <- DimPlot(integrated, label = T, repel = TRUE, pt.size = 1,
            reduction = "umap.harmony.denoised",group.by = "individualID")+umap_theme+ggtitle("DonorID")
p2 <- FeaturePlot(integrated, label = T, repel = TRUE, 
              reduction = "umap.harmony.denoised",features = "decontX_contamination")+umap_theme+ggtitle("Contamination_level")
p3 <- FeaturePlot(integrated, repel = TRUE, 
              reduction = "umap.harmony.denoised",features = "nCount_RNA")+umap_theme+ggtitle("nCount_RNA")
p4 <- FeaturePlot(integrated, repel = TRUE, 
              reduction = "umap.harmony.denoised",features = "nFeature_RNA")+umap_theme+ggtitle("nFeature_RNA")                                          

ggsave(filename = "./Results/Revision_2/Figures/0_integrated_UMAP_after_decontX.pdf",
    patchwork::wrap_plots(p1,p2,p3,p4, ncol = 2),
      scale = 1, width = 13, height = 10)

##################################################################################
###############-- Saving Object for downstream analysis --########################
##################################################################################
saveRDS(integrated,file = "./Results/Revision_2/03_1.denoised_object_dbl_rm.rds")


#### Saving as h5 file for scanpy pipeline
# DefaultAssay(integrated) = "decontX"
# integrated = JoinLayers(integrated)
# integrated[["SCT"]] = NULL
# integrated[["RNA"]]$counts = round(LayerData(integrated,layer = "counts",assay = "decontX"))
# ### remove unwanted assays
# # integrated[["decontX"]] = NULL
# # integrated[["SCT"]] = NULL

# ## saving as h5ad
# DefaultAssay(integrated) <- "RNA"
# integrated = JoinLayers(integrated)
# integrated <- NormalizeData(integrated) %>% 
#                 FindVariableFeatures(nfeatures = 1000,verbose =F) %>%
#                 ScaleData(verbose =F) %>%
#                 RunPCA(ndims=30,verbose =F) %>%
#                 RunUMAP(dims = 1:30,reduction.name="umap.rna",reduction.key = "rnaUMAP_",verbose =F)

# integrated[["RNA"]] <- as(object = integrated[["RNA"]], Class = "Assay")
# integrated[["RNA"]] # check the version of the seurat assay

# integrated[["RNA"]]$scale.data = NULL
# integrated[["RNA"]]$data = NULL

# sceasy::convertFormat(integrated, from="seurat", to="anndata",
#                        outFile='./Results/Revision_2/03.integrated_object_decontX.h5ad')