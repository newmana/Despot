import MENDER
import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.metrics import *

# Using leiden to do dimension reduction and clustering
def Spatial_Cluster_MENDER(adata, n_top_genes=2000, n_comps=50):
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor='seurat', inplace=True)
    ncp = n_comps if adata.shape[1] > n_comps else adata.shape[1] - 1
    use_hvg = True if adata.shape[1] > n_comps else False
    sc.pp.pca(adata, n_comps=ncp, use_highly_variable=use_hvg, svd_solver='arpack')
    sc.pp.neighbors(adata)
    # use umap and leiden for clustering
    sc.tl.umap(adata)
    sc.tl.leiden(adata, key_added='clusters')
    # input parameters of MENDER
    scale = 6

    # main body of MENDER
    msm = MENDER.MENDER_single(
        adata,
        # determine which cell state to use
        # we use the cell state got by Leiden
        ct_obs='ct'
    )
    msm.estimate_radius()
    # set the MENDER parameters

    msm.set_MENDER_para(
        # default of n_scales is 6
        n_scales=scale,

        # for single cell data, nn_mode is set to 'radius'
        nn_mode='radius',

        # default of n_scales is 15 um (see the manuscript for why).
        # MENDER also provide a function 'estimate_radius' for estimating the radius
        # this para in 3D might be smaller than 2D
        nn_para=7,

    )
    # construct the context representation
    msm.run_representation(

        # the number of processings
    )

    # set the spatial clustering parameter
    # positive values for the expected number of domains
    # negative values for the clustering resolution
    msm.run_clustering_normal(-0.5)
    return msm.adata_MENDER