import pkg_resources
from utils.check import Check_Requirements
import subprocess
from utils.common import *


def MENDER_install():
    package = "MENDER"
    print(f"Dependencies will be installed when Using {package} for the first time.")
    # download SPROD handle python dependencies
    if "somender" not in {pkg.key for pkg in pkg_resources.working_set}:
        # handle python dependencies
        py_req = Check_Requirements({"louvain", "scikit-image", "numpy", "igraph", "scikit-learn", "umap-learn"})
        python = sys.executable
        py_ins = subprocess.check_call([python, '-m', 'pip', 'install', "SOMENDER"], stdout=subprocess.DEVNULL)
        if py_ins+py_req == 0:
            print(f"{package} installation succeeded.")
            return 0
        else:
            print(f"{package} installation failed.")
            exit(-1)


def Spatial_Cluster_MENDER(adata, n_comps=50):
    import MENDER
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
        ct_obs='clusters'
    )
    msm.estimate_radius()
    # set the MENDER parameters

    msm.set_MENDER_para(n_scales=scale, nn_mode='radius',nn_para=10)
    # construct the context representation
    msm.run_representation()

    # set the spatial clustering parameter
    # positive values for the expected number of domains
    # negative values for the clustering resolution
    msm.run_clustering_normal(-0.5)
    return msm.adata_MENDER