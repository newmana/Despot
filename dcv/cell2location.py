from utils.common import *
from utils.io import *
import matplotlib as mpl
import pkg_resources
from utils.check import Check_Requirements
import subprocess
import cell2location

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text

results_folder = './temp'

# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'


def cell2location_install():
    # download cell2location and handle python dependencies
    if "cell2location" not in {pkg.key for pkg in pkg_resources.working_set}:
        print("Dependencies will be installed when Using cell2location for the first time.")

        # handle python dependencies
        py_req = Check_Requirements({"pyro-ppl", "scvi-tools", "black", "flake8", "pytest", "pytest-cov",
                                     "isort", "pre-commit", "jinja2", "leidenalg"})

        python = sys.executable
        py_ins = subprocess.check_call([python, '-m', 'pip', 'install', "cell2location"], stdout=subprocess.DEVNULL)

        if py_req + py_ins == 0:
            print("cell2location installed successfully.")
            return 0
        else:
            print("cell2location installation failed.")
            exit(-1)
    else:
        print("cell2location has been installed correctly.")
        return 0

def Cell2Location_pp_sc(smdFile, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12):
    adata_sc = Load_smd_sc_to_AnnData(smdFile)
    adata_sc.obs_names_make_unique()
    adata_sc.var_names_make_unique()
    # Use ENSEMBL as gene IDs to make sure IDs are unique and correctly matched
    adata_sc.var['SYMBOL'] = adata_sc.var.index
    plt.ion()
    selected = cell2location.utils.filtering.filter_genes(adata_sc,
                            cell_count_cutoff=cell_count_cutoff,
                            cell_percentage_cutoff2=cell_percentage_cutoff2,
                            nonz_mean_cutoff=nonz_mean_cutoff)
    plt.ioff()
    # filter the object
    adata_sc = adata_sc[:, selected].copy()
    return adata_sc


def Cell2Location_rg_sc(adata_sc, max_epoches=250, batch_size=2500, train_size=1, lr=0.002,
                        num_samples=1000, use_gpu=True):
    from cell2location.models import RegressionModel

    # prepare anndata for the regression model
    cell2location.models.RegressionModel.setup_anndata(adata=adata_sc,
                            # cell type, covariate used for constructing signatures
                            labels_key='annotation')

    # create and train the regression model
    mod = RegressionModel(adata_sc)
    # mod.view_anndata_setup()

    # Use all data for training (validation not implemented yet, train_size=1)
    mod.train(max_epochs=max_epoches, batch_size=batch_size, train_size=train_size, lr=lr)

    # plot ELBO loss history during training, removing first 20 epochs from the plot
    # mod.plot_history(20)

    # In this section, we export the estimated cell abundance (summary of the posterior distribution).
    adata_sc = mod.export_posterior(
        adata_sc, sample_kwargs={'num_samples': num_samples, 'batch_size': batch_size, 'use_gpu': use_gpu}
    )

    # export estimated expression in each cluster
    if 'means_per_cluster_mu_fg' in adata_sc.varm.keys():
        inf_aver = adata_sc.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                        for i in adata_sc.uns['mod']['factor_names']]].copy()
    else:
        inf_aver = adata_sc.var[[f'means_per_cluster_mu_fg_{i}'
                                        for i in adata_sc.uns['mod']['factor_names']]].copy()
    inf_aver.columns = adata_sc.uns['mod']['factor_names']
    return inf_aver

def Cell2Location_pp_sp(smdFile:str, h5data:str) -> ad.AnnData:
    info = smdInfo(smdFile)
    adata_sp = Load_smd_to_AnnData(smdFile, h5data=h5data)
    if h5data == 'SPCS_mat':
        adata_sp.X = np.array(np.exp(adata_sp.X) - 1)
    adata_sp.var_names_make_unique()
    adata_sp.obs['sample'] = info.name

    # rename genes to ENSEMBL
    adata_sp.var['SYMBOL'] = adata_sp.var_names

    # find mitochondria-encoded (MT) genes
    adata_sp.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_sp.var['SYMBOL']]

    # remove MT genes for spatial mapping (keeping their counts in the object)
    adata_sp.obsm['MT'] = adata_sp[:, adata_sp.var['MT_gene'].values].X.toarray()
    adata_sp = adata_sp[:, ~adata_sp.var['MT_gene'].values]
    return adata_sp


def Cell2Location_rg_sp(adata_sp, inf_aver, N_cells_per_location=30, detection_alpha=20,
                        max_epoches=10000, batch_size=None, train_size=1, lr=0.002,
                        num_samples=1000, use_gpu=True):
    # do spatial mapping
    # find shared genes and subset both anndata and reference signatures
    intersect = np.intersect1d(adata_sp.var['gene_name'], inf_aver.index)
    adata_sp.var_names = adata_sp.var['gene_name']
    adata_sp.var_names_make_unique()
    adata_sp = adata_sp[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()
    # prepare anndata for cell2location model
    cell2location.models.Cell2location.setup_anndata(adata=adata_sp, batch_key="sample")

    # create and train the model
    mod = cell2location.models.Cell2location(
        adata_sp, cell_state_df=inf_aver,
        # the expected average cell abundance: tissue-dependent
        # hyper-prior which can be estimated from paired histology:
        N_cells_per_location=N_cells_per_location,
        # hyperparameter controlling normalisation of
        # within-experiment variation in RNA detection:
        detection_alpha=detection_alpha
    )
    # mod.view_anndata_setup()
    mod.train(max_epochs=max_epoches,
              # train using full data (batch_size=None)
              batch_size=batch_size,
              # use all data points in training because
              # we need to estimate cell abundance at all locations
              train_size=train_size,
              use_gpu=use_gpu)

    # plot ELBO loss history during training, removing first 100 epochs from the plot
    # mod.plot_history(1000)

    # In this section, we export the estimated cell abundance (summary of the posterior distribution).
    adata_sp = mod.export_posterior(
        adata_sp, sample_kwargs={'num_samples': num_samples, 'batch_size': mod.adata.n_obs, 'use_gpu': use_gpu}
    )

    # add 5% quantile, representing confident cell abundance, 'at least this amount is present',
    # to adata.obs with nice names for plotting
    adata_sp.obs[adata_sp.uns['mod']['factor_names']] = adata_sp.obsm['q05_cell_abundance_w_sf']
    return adata_sp


def Cell2Location_run(smdFile, h5data, sc_max_epoches=250, sc_batch_size=2500, sc_train_size=1, sc_lr=0.002, sc_num_samples=1000,
                      N_cells_per_location=30, detection_alpha=20,
                      sp_max_epoches=10000, sp_batch_size=2500, sp_train_size=1, sp_lr=0.002, sp_num_samples=1000, use_gpu=True):
    adata_sc = Cell2Location_pp_sc(smdFile)
    adata_sp = Cell2Location_pp_sp(smdFile, h5data=h5data)
    inf_aver = Cell2Location_rg_sc(adata_sc, sc_max_epoches, sc_batch_size, sc_train_size,sc_lr,sc_num_samples, use_gpu)
    adata_sp = Cell2Location_rg_sp(adata_sp, inf_aver,
                                   N_cells_per_location=N_cells_per_location,
                                   detection_alpha=detection_alpha,
                                   max_epoches=sp_max_epoches,
                                   batch_size=sp_batch_size,
                                   train_size=sp_train_size,
                                   lr=sp_lr,
                                   num_samples=sp_num_samples,
                                   use_gpu=use_gpu)
    weight = adata_sp.obsm['q05_cell_abundance_w_sf']
    weight.columns = adata_sp.uns['mod']['factor_names']
    adata_sp.obsm['Cell2Location'] = weight
    return adata_sp




