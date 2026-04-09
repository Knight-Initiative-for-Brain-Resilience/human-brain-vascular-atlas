### Seurat preprocessing and integration
suppressPackageStartupMessages({
    library(sceasy)
    library(Seurat)
    library(stringr)
    library(dplyr)
    library(patchwork)
    library(ggplot2)
    library(future)
    library(tidydr)
    library(harmony)
    library(scCustomize)
    library(scDblFinder)
    })

# Provide enough memory for preprocessing  
options(future.globals.maxSize = 540*1024^3)

## Set the working directory
# dir = "/oak/stanford/projects/kibr/Reorganizing/Projects/Andi/vascular_atlas/"
# setwd(dir)

## load in the metadata for cellranger files
meta_df <- read.csv("./Results/preprocessing_files_1203.csv")
meta_df$sample_path <- paste(meta_df$samplepath,"/outs/",sep = "")
## Correct the path
meta_df$sample_mtx_path = gsub("processed","cr_data",meta_df$sample_mtx_path)
meta_df$sample_path = gsub("processed","cr_data",meta_df$sample_mtx_path)

## Provide the preprocessing as a function
object_preprocessing <- function(data.dir){
    t1 = Sys.time()
    print(t1)

    ## Read in the count matrix
    obj.data <-Read10X(data.dir = data.dir)

    ## Create seurat object
    obj <- CreateSeuratObject(counts = obj.data, min.cells = 3, min.features = 200)
    print(paste("Raw dimension:",dim(obj)))

    ## Initial filteration based on number of features in each cell
    obj <- subset(
        x = obj,
        subset = nFeature_RNA < 10000 &
        nFeature_RNA > 200)
    print(paste("Dimensions after initial filteration:",dim(obj)))

    ##############################################################################
    # calculate k means n = 2 based on mitocondril reads -------------------------
    ##############################################################################
    DefaultAssay(obj) <- "RNA"

    # calculating MT reads
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

    # then keep the group with lower mean valus 
    mt_kmeans <- kmeans(obj$percent.mt,centers = 2)
    print(mt_kmeans$centers)
    print(mt_kmeans$size)

    obj$mt_keep <- mt_kmeans$cluster 
    if(mt_kmeans$centers[1] > mt_kmeans$centers[2]){
        obj <- subset(x = obj,subset = mt_keep == 2)
        }else{
            obj <- subset(x = obj,subset = mt_keep == 1)
            }
    cat(paste("MT_range after filteration: ", range(obj$percent.mt)))

    ## Doublet detection using scDblFinder
    ## More conservative than DoubletFinder
    sce = scDblFinder(LayerData(obj),clusters=FALSE)
    obj$scDblFinder.class = sce$scDblFinder.class
    obj$scDblFinder.score = sce$scDblFinder.score

    # get time
    t2 = Sys.time()
    print(t2-t1)

    #return object
    return(obj)
  }

################################################################################
######################## create the object(s) list     #########################
################################################################################
## Creating an empty object list
object.list <- list()
i = 1
# read in data and preprocessing the data based on the parameters
for (i in 1:length(meta_df$sampleID)) {
  print(paste("Sample",meta_df[i,]$sampleID,"on processing. /"))
  # generating two paths
  data.dir <- meta_df[i,]$sample_path
  # preproceesing each object
  object.list[[i]] <- object_preprocessing(
    data.dir = data.dir
  )
  print(paste("Sample",meta_df[i,]$sampleID,"fisnished! /"))
}
print(object.list)

# ################################################################################
# ############ Doublet removal based on DoubletFinder results ####################
# ################################################################################
for (i in 1:length(object.list)){
    #cat(i)
    meta <- object.list[[i]]@meta.data
    meta$doublet <- ifelse(meta[,7] == "Doublet",T,F)
    # deleting the two columns contain redundant information
    meta <- meta[,-c(6:7)]
    object.list[[i]]@meta.data <- meta
    
    object.list[[i]] <- subset(object.list[[i]],subset = doublet == F)
}

saveRDS(object.list,file = "./Results/Revision_2/01.merged_objects_list.rds")
# object.list <- readRDS("./Results/Revision_2/01.merged_objects_list.rds")

################################################################################
### Merging all individual object to a merged object for downstream analysis ###
################################################################################
## removing data and scale.data because of merging time. 
for (i in 1:length(object.list)){
  DefaultAssay(object.list[[i]]) <- "RNA"
  # object.list[[i]][["RNA"]]$data <- NULL
  # object.list[[i]][["RNA"]]$scale.data <- NULL
}

# merging and adding sample id
object.merge <- merge(object.list[[1]], y = c(object.list[[2]], object.list[[3]], object.list[[4]], object.list[[5]], object.list[[6]], object.list[[7]],
                                             object.list[[8]], object.list[[9]], object.list[[10]], object.list[[11]], object.list[[12]], object.list[[13]],object.list[[14]], object.list[[15]],
                                             object.list[[16]], object.list[[17]], object.list[[18]], object.list[[19]], object.list[[20]], object.list[[21]],
                                             object.list[[22]], object.list[[23]], object.list[[24]], object.list[[25]], object.list[[26]], object.list[[27]],
                                             object.list[[28]], object.list[[29]], object.list[[30]], object.list[[31]], object.list[[32]], object.list[[33]],
                                             object.list[[34]], object.list[[35]], object.list[[36]], object.list[[37]], object.list[[38]], object.list[[39]],
                                             object.list[[40]], object.list[[41]], object.list[[42]], object.list[[43]], object.list[[44]], object.list[[45]],
                                             object.list[[46]], object.list[[47]], object.list[[48]], object.list[[49]], object.list[[50]], object.list[[51]],
                                             object.list[[52]], object.list[[53]], object.list[[54]], object.list[[55]], object.list[[56]], object.list[[57]],
                                             object.list[[58]], object.list[[59]], object.list[[60]], object.list[[61]], object.list[[62]], object.list[[63]],
                                             object.list[[64]], object.list[[65]], object.list[[66]], object.list[[67]], object.list[[68]], object.list[[69]],
                                             object.list[[70]], object.list[[71]], object.list[[72]], object.list[[73]], object.list[[74]], object.list[[75]],
                                             object.list[[76]], object.list[[77]], object.list[[78]], object.list[[79]], object.list[[80]], object.list[[81]],
                                             object.list[[82]], object.list[[83]], object.list[[84]], object.list[[85]], object.list[[86]], object.list[[87]],
                                             object.list[[88]], object.list[[89]], object.list[[90]], object.list[[91]], object.list[[92]], object.list[[93]],
                                             object.list[[94]], object.list[[95]], object.list[[96]], object.list[[97]], object.list[[98]], object.list[[99]],
                                             object.list[[100]], object.list[[101]], object.list[[102]], object.list[[103]], object.list[[104]], object.list[[105]],
                                             object.list[[106]], object.list[[107]], object.list[[108]], object.list[[109]], object.list[[110]], object.list[[111]],
                                             object.list[[112]], object.list[[113]], object.list[[114]], object.list[[115]], object.list[[116]], object.list[[117]],
                                             object.list[[118]], object.list[[119]], object.list[[120]], object.list[[121]], object.list[[122]], object.list[[123]],
                                             object.list[[124]], object.list[[125]], object.list[[126]], object.list[[127]], object.list[[128]], object.list[[129]],
                                             object.list[[130]], object.list[[131]], object.list[[132]], object.list[[133]], object.list[[134]], object.list[[135]],
                                             object.list[[136]], object.list[[137]], object.list[[138]], object.list[[139]], object.list[[140]], object.list[[141]],
                                             object.list[[142]], object.list[[143]], object.list[[144]], object.list[[145]], object.list[[146]]),
                      add.cell.ids = meta_df$sampleID)

######################################################################
#################### Organizing all the metadata #####################  
######################################################################
meta <- object.merge@meta.data

## adding donor information
temp <- str_split_fixed(rownames(meta),"_",5)
meta$individualID <- paste(temp[,1],temp[,2], sep = "_")
meta$sampleID <- paste(temp[,1],temp[,2],temp[,3],temp[,4],sep = "_")
ID <- match(meta$sampleID, meta_df$sampleID)
## Region information
meta$regionID <- meta_df[ID,]$regionID
meta$region_name <- meta_df[ID,]$region_name
meta$region_abb <- meta_df[ID,]$region_short
meta$sorting_date <- meta_df[ID,]$sorting_date
meta$sequencing_run <- meta_df[ID,]$Sequencing_run
## Additional region information
df <- read.csv("./data/region.csv")
df$brain_region <- paste(df$regionID, df$region_abb,sep = "_")
meta$brain_region <- paste(meta$regionID, meta$region_abb,sep = "_")
id = match(meta$brain_region,df$brain_region)
meta$region_layer <- df[id,]$Region_layer_4
table(meta$region_layer)
## Age at death
meta$ageatdeath <- NA
meta[meta$individualID == "UW_7029",]$ageatdeath <- 30
meta[meta$individualID == "LB_4008",]$ageatdeath <- 35
meta[meta$individualID == "LB_9770",]$ageatdeath <- 40
meta[meta$individualID == "Stanford_12052023",]$ageatdeath <- 49
## Sex
meta$sex <- ifelse(meta$individualID %in% c("UW_7029","LB_9770"), "M","F")
## assign back the meta data
object.merge@meta.data <- meta

###########################
### Saving object again ###
###########################
saveRDS(object.merge, file = "./Results/Revision_2/01.merged_objects.rds")

# object.merge <- readRDS("./Results/Revision/01.merged_objects.rds")

######################################################################
################ Prepare for Harmony integration #####################  
######################################################################
## Merging the layers and split the layers based on batch effect
object.merge <- JoinLayers(object.merge)
object.merge[["RNA"]] <- split(object.merge[["RNA"]], f = object.merge$individualID)

DefaultAssay(object.merge) <- "RNA"
object.merge <- SCTransform(object.merge,assay = "RNA",new.assay.name = "SCT",vars.to.regress = c("percent.mt"),layer = "counts",verbose = T,method = "glmGamPoi") %>% 
        RunPCA(ndims=30) %>% 
        FindNeighbors(dims = 1:30) %>% 
        RunUMAP(dims = 1:30, reduction.name="umap.rna",reduction.key = "rnaUMAP_")

######################################################################
###########-- FIRST UMAP of unintegrated overall raw data--###########
######################################################################
umap_theme <- theme_dr()+theme(panel.grid.major = element_blank(), 
                                            panel.grid.minor = element_blank(),
                                            panel.background = element_blank(), 
                                            axis.line = element_line(colour = "black"))

p1 <- DimPlot(object.merge, label = T, repel = TRUE, 
              reduction = "umap.rna",group.by = "individualID")+umap_theme+ggtitle("IndividualID")
p2 <- DimPlot(object.merge, label = T, repel = TRUE, 
              reduction = "umap.rna",group.by = "scDblFinder.class")+umap_theme+ggtitle("Doublets Prediction")
p3 <- FeaturePlot(object.merge, repel = TRUE, 
              reduction = "umap.rna",features = "nCount_RNA")+umap_theme+ggtitle("nCount_RNA")
p4 <- FeaturePlot(object.merge, repel = TRUE, 
              reduction = "umap.rna",features = "nFeature_RNA")+umap_theme+ggtitle("nFeature_RNA")

ggsave(filename = "./Results/Revision_2/Figures/0_merged.object_UMAP_individual_revision.pdf",
    patchwork::wrap_plots(p1,p2,p3,p4, ncol = 2),
      scale = 1, width = 10, height = 10)

##################################################################################
############################## Harmony integration ###############################
##################################################################################
integrated <- IntegrateLayers(
  object = object.merge, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

integrated <- FindNeighbors(integrated, reduction = "harmony", dims = 1:30)
integrated <- FindClusters(integrated, resolution = c(0.5,1,1.5,2))
integrated <- RunUMAP(integrated, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

## Plotting
p1 <- DimPlot(integrated, label = T, repel = TRUE, 
              reduction = "umap.harmony",group.by = "individualID")+umap_theme+ggtitle("IndividualID")
p2 <- DimPlot(integrated, label = T, repel = TRUE, 
              reduction = "umap.harmony",group.by = "scDblFinder.class")+umap_theme+ggtitle("Doublets Prediction")
p3 <- FeaturePlot(integrated, repel = TRUE, 
              reduction = "umap.harmony",features = "nCount_RNA")+umap_theme+ggtitle("nCount_RNA")
p4 <- FeaturePlot(integrated, repel = TRUE, 
              reduction = "umap.harmony",features = "nFeature_RNA")+umap_theme+ggtitle("nFeature_RNA")

ggsave(filename = "./Results/Revision_2/Figures/0_merged.object_UMAP_individual_layered_harmony.pdf",
    patchwork::wrap_plots(p1,p2,p3,p4, ncol = 2),
      scale = 1, width = 10, height = 10)

##################################################################################
################################# Saving Object ##################################
##################################################################################
print(integrated)
saveRDS(integrated,file = "./Results/Revision_2/02.integrated_object_revision.rds")
