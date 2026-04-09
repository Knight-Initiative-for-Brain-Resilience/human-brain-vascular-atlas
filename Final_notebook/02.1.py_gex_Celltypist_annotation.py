### Doing cell type annotation using celltypist
# --- imports (as you had) ---
import os, gc, psutil, time
import scanpy as sc
import anndata as ad
import celltypist
from celltypist import models
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from rich import print

print(psutil.virtual_memory())
gc.collect()
PATH="/oak/stanford/projects/kibr/Reorganizing/Projects/Andi/vascular_atlas"
os.chdir(PATH)

# --- load query data ---
# adata = sc.read_h5ad(PATH + "/Results/GEX/merged_anndata_normalized_rm_doublets_400.h5ad")
# adata = sc.read_h5ad(PATH + "/Results/Revision_2/02.integrated_object_revision.h5ad")
adata = sc.read_h5ad(PATH + "/Results/Revision_2/03.integrated_object_decontX.h5ad")
print(adata)

# Saving count data
adata.raw = adata.copy()
adata.layers["counts"] = adata.X.copy()
print(adata)

# Normalizing to median total counts
sc.pp.normalize_per_cell(adata, counts_per_cell_after=10**4)  # normalize to 10,000 counts per cell
sc.pp.log1p(adata)

# In case need to do the normalization again and need to convert the sparse matrix into array
adata.X = adata.X.toarray()

#### Actuall cell type annotation command using Celltypist
########################################################################
############## Loading model from celltypist_AY ########################
########################################################################
model = models.Model.load(model = './data/celltypist_model/celltypist_model_AY.pkl')
print("Model cell types:", model.cell_types)

## Prediction using AY data
t_start = time.time()
prediction = celltypist.annotate(adata, model = model, majority_voting=True)

## Transform the prediction to adata to get the full output
prediction_adata = prediction.to_adata()
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds")

## copy the results to the original AnnData object
adata.obs["celltypist_cell_label_AY"] = prediction_adata.obs.loc[adata.obs.index, "majority_voting"]
adata.obs["celltypist_conf_score_AY"] = prediction_adata.obs.loc[adata.obs.index, "conf_score"]

print(pd.crosstab(adata.obs["celltypist_cell_label_AY"],adata.obs["sampleID"]))

########################################################################
############ Loading model from celltypist_Garcia ######################
########################################################################
model = models.Model.load(model = './data/celltypist_model/celltypist_model_Garcia.pkl')
print("Model cell types:", model.cell_types)

## Prediction using AY data
t_start = time.time()
prediction = celltypist.annotate(adata, model = model, majority_voting=True)

## Transform the prediction to adata to get the full output
prediction_adata = prediction.to_adata()
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds")

## copy the results to the original AnnData object
adata.obs["celltypist_cell_label_Garcia"] = prediction_adata.obs.loc[adata.obs.index, "majority_voting"]
adata.obs["celltypist_conf_score_Garcia"] = prediction_adata.obs.loc[adata.obs.index, "conf_score"]

print(pd.crosstab(adata.obs["celltypist_cell_label_Garcia"],adata.obs["sampleID"]))

########################################################################
# ############ Loading model from celltypist_Walchli #####################
# ########################################################################
# model = models.Model.load(model = './data/celltypist_model/celltypist_model_Walchli.pkl')
# print("Model cell types:", model.cell_types)

# ## Prediction using AY data
# t_start = time.time()
# prediction = celltypist.annotate(adata, model = model, majority_voting=True)

# ## Transform the prediction to adata to get the full output
# prediction_adata = prediction.to_adata()
# t_end = time.time()
# print(f"Time elapsed: {t_end - t_start} seconds")

# ## copy the results to the original AnnData object
# adata.obs["celltypist_cell_label_Walchli"] = prediction_adata.obs.loc[adata.obs.index, "majority_voting"]
# adata.obs["celltypist_conf_score_Walchli"] = prediction_adata.obs.loc[adata.obs.index, "conf_score"]

# print(pd.crosstab(adata.obs["celltypist_cell_label_Walchli"],adata.obs["sampleID"]))

## converting the matrix into sparse matrix
adata.X = csr_matrix(adata.X)

## show the memory usage after cleaning
print(adata)
print(psutil.virtual_memory())

meta = adata.obs
meta_path = PATH+"/Results/Revision_2/celltypist_annotation_decontX.csv"
meta.to_csv(meta_path)

## Saving the annotated data
adata = adata.raw.to_adata()
results_file = PATH+"/Results/Revision_2/03.integrated_object_decontX.h5ad"

adata.write_h5ad(
    results_file,
    #compression=hdf5plugin.FILTERS["zstd"]
)

