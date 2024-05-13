from utils.common import *


# Using leiden to do dimension reduction and clustering
def Spatial_Cluster_Analysis(adata, n_top_genes=2000, n_comps=50):
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
    return adata


# obsm = pd.DataFrame()
# for idx, model in enumerate(models):
#     losses = []
#     with torch.no_grad():
#         for i, x in enumerate(data_loader0):
#             x = x.to(device).view(-1, count_size)
#             x_reconst, mu, log_var = model(x)
#             reconst_loss = F.binary_cross_entropy(x_reconst, x, reduction='mean')
#             kl_div = - 0.5 * torch.sum(1 + log_var - mu.pow(2) - log_var.exp())
#
#             loss = np.float((reconst_loss + kl_div).cpu())
#             losses.append(loss)
#     obsm[uni_types[idx]] = losses