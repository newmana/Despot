from utils.common import *
import pkg_resources
from utils.check import Check_Requirements
import subprocess

def Harmony_install():
    print("Dependencies will be installed when Using Harmony for the first time.")
    # download SPROD handle python dependencies
    if "harmonypy" not in {pkg.key for pkg in pkg_resources.working_set}:
        # handle python dependencies
        py_req = Check_Requirements({"scanpy"})
        python = sys.executable
        py_ins = subprocess.check_call([python, '-m', 'pip', 'install', "harmonypy"], stdout=subprocess.DEVNULL)
        if py_ins+py_req == 0:
            print("Harmony installation succeeded.")
            return 0
        else:
            print("Harmony installation failed.")
            exit(-1)

            
# Using leiden to do dimension reduction and clustering
def Embed_Harmony(adata: ad.AnnData, key: str="batch", n_comps:int=50, max_iter:int=20) -> ad.AnnData:
    if key in adata.obs.columns:
      sc.pp.normalize_total(adata, inplace=True)
      sc.pp.log1p(adata)
      sc.pp.highly_variable_genes(adata, flavor='seurat', inplace=True)
      ncp = n_comps if adata.shape[1] > n_comps else adata.shape[1] - 1
      use_hvg = True if adata.shape[1] > n_comps else False
      sc.pp.pca(adata, n_comps=ncp, use_highly_variable=use_hvg, svd_solver='arpack')
      sc.external.pp.harmony_integrate(adata, key=key,max_iter_harmony=max_iter)
      sc.pp.neighbors(adata, use_rep="X_pca_harmony")
      # use umap and leiden for clustering
      sc.tl.umap(adata)
      sc.tl.tsne(adata)
    else:
      print("No batches in smdFiles, skip embedding.")
    return adata