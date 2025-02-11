import subprocess
import time
import pkg_resources
from typing import Union, List
from utils.io import *
from utils.preprocess import *
from utils.geo import Despot_Find_bestGroup


def Pip_cluVisual(smdFile, h5data, imgPath):
    adata = Load_smd_to_AnnData(smdFile, h5data)
    figa, axsa = plt.subplots(2, 3, figsize=(9, 6), constrained_layout=True)
    sc.pl.spatial(adata, img_key="lowres", color="Seurat", size=1.5, ax=axsa[0, 0], show=False)
    sc.pl.spatial(adata, img_key="lowres", color="Giotto", size=1.5, ax=axsa[0, 1], show=False)
    sc.pl.spatial(adata, img_key="lowres", color="ground_truth", size=1.5, ax=axsa[0, 2], show=False)
    sc.pl.spatial(adata, img_key="lowres", color="BayesSpace", size=1.5, ax=axsa[1, 0], show=False)
    sc.pl.spatial(adata, img_key="lowres", color="Squidpy", size=1.5, ax=axsa[1, 1], show=False)
    sc.pl.spatial(adata, img_key="lowres", color="stlearn", size=1.5, ax=axsa[1, 2], show=False)
    figa.savefig(imgPath)


def Pip_estimate(smdFile, h5data='matrix', method='SpatialDE', force=False):
    method_list = ['SpatialDE', 'SpaGCN', 'SPARK', 'Seurat', 'BayesSpace', 'Giotto', 'Squidpy']
    with h5.File(smdFile, 'r') as f:
        if method not in method_list:
            print("Clustering method:{0} has not been supported yet, default `leiden`.".format(method))
            method = 'SpatialDE'
        if (h5data + '/features/is_HVG/' + method in f) and (force == False):
            print("Estimating with " + method + " has done, skip it.")
            return
    if method == 'SpatialDE':
        from est.Estim_SpatialDE import HDG_Analysis_by_SpatialDE
        adata = Load_smd_to_AnnData(smdFile, h5data)
        adata.uns['HDGs'] = HDG_Analysis_by_SpatialDE(adata)
        Save_smd_from_SpatialDE(smdFile, adata, h5data)
    elif method == 'SpaGCN':
        from est.Estim_SpaGCN import HDG_Analysis_by_SpaGCN
        adata = Load_smd_to_AnnData(smdFile, h5data)
        adata = HDG_Analysis_by_SpaGCN(adata)
        Save_smd_from_SpaGCN_est(smdFile, adata, h5data)
    elif method == 'SPARK':
        # Using R scripts
        os.system("Rscript est/Estim_SPARK.R")
    elif method == 'trendsceek':
        # Using R scripts
        os.system("Rscript est/Estim_trendsceek.R")
    elif method == 'Giotto':
        # Using R scripts
        os.system("Rscript est/Estim_Giotto.R")


def Pip_benchmark(smdFile, h5data):
    adata = Load_smd_to_AnnData(smdFile, h5data)
    clu = pd.DataFrame(adata.obs[['Seurat', 'SpaGCN', 'stlearn', 'BayesSpace', 'leiden', 'Giotto']], dtype='category')
    from sklearn import metrics
    ground_truth = adata.obs['ground_truth'].astype('category')
    benchmarks = pd.DataFrame(columns=clu.columns)
    NMI = [metrics.normalized_mutual_info_score(clu[method], ground_truth) for method in clu.columns]
    ARI = [metrics.adjusted_rand_score(clu[method], ground_truth) for method in clu.columns]
    benchmarks.loc['NMI'] = NMI
    benchmarks.loc['ARI'] = ARI
    benchmarks.to_csv("temps/151673_" + h5data + ".csv", index=True, header=True)


def Create_h5datas(cfg: Union[dict, List[dict]]) -> List[str]:
    # whether you need Decontamination
    if isinstance(cfg, dict):
        cfgs = [cfg]
    all_h5datas = []
    for i,c in enumerate(cfgs):
        h5datas = []
        Decont = c['Decontamination']
        smdFile = c['smdFile']
        if type(Decont) != list:
            Decont = [Decont]
        with h5.File(smdFile, 'r') as f:
            for dec in Decont:
                if dec == "SpotClean":
                    h5data = 'SpotClean_mat'
                elif dec == "SPCS":
                    h5data = "SPCS_mat"
                elif dec == "SPROD":
                    h5data = "SPROD_mat"
                elif dec == "none":
                    h5data = "matrix"
                else:
                    h5data = dec
                if h5data in f:
                    h5datas.append(h5data)
                else:
                    print("{0} can't be found in smdFile.".format(h5data))
            if i==0:
                all_h5datas.extend(h5datas)
            else:
                all_h5datas = list(set(all_h5datas).intersection(h5datas))
    return all_h5datas


def Despot_Estimate(smdFile, cfg, force=False):
    # whether need Decontamination
    h5datas = Create_h5datas(cfg)

    methods = cfg['Estimation']
    for h5data in h5datas:
        for method in methods:
            Pip_estimate(smdFile, h5data, method=method, force=force)


def Despot_Benchmark(smdFile, cfg, mode: str = "cluster", force: bool = False):
    mode_list = ['cluster', 'estimate', 'deconvolution']
    h5datas = []
    # whether need Decontamination
    Decont = cfg['Decontamination']
    for dec in Decont:
        if dec == "SpotClean":
            h5data = 'SpotClean_mat'
            h5datas.append(h5data)
        elif dec == "SPCS":
            h5data = "SPCS_mat"
            h5datas.append(h5data)
        elif dec == "none":
            h5data = "matrix"
            h5datas.append(h5data)
        else:
            print("no matched Decontamination method, default `none`")

    if len(h5datas) == 0:
        h5datas.append("matrix")

    bm_done = True
    with h5.File(smdFile, 'r') as f:
        if mode not in mode_list:
            print("Benchmark mode:{0} has not been supported yet".format(mode))
            return
        for h5data in h5datas:
            if (h5data + '/benchmark/' + mode not in f) or (force is True):
                bm_done = False
                break
    if bm_done is True:
        print("Benchmark with " + mode + " has done, skip it.")
        return
    if mode == "cluster":
        print("Computing all cluster benchmarks in smdFile.")
        os.system("Rscript clu/Benchmark_cluster.R")
    elif mode == "estimate":
        print("Computing all estimate benchmarks in smdFile.")
        os.system("Rscript est/Benchmark_estimate.R")


# we provide 2 methods for decontamination
def Pip_decont(smdFile, cfg, method='none', force=False, ):
    # checking whether need to do decontamination
    method_list = ['none', 'SpotClean', 'SPCS', 'SPROD']
    with h5.File(smdFile, 'r') as f:
        if method not in method_list:
            print("Decontamination method:{0} has not been supported yet, default `none`.".format(method))
            return
        if (method + "_mat" in f) and (force == False):
            print(method + " has done, skip it.")
            return

    # do decontamination
    if method == 'SpotClean':
        print("Decontamination method: SpotClean.")
        if cfg['platform'] == 'ST':
            print("SpotClean only support 10X series. But your platform is ST")
        else:
            os.system("Rscript dct/Decont_SpotClean.R")
    elif method == 'SPCS':
        print("Decontamination method: SPCS.")
        os.system("Rscript dct/Decont_SPCS.R")
    elif method == 'SPROD':
        print("Decontamination method: SPROD.")
        from dct.Decont_SPROD import sprod_install, sprod_pp, sprod_run, sprod_save
        sprod_install(cfg['pythonPath'])
        sprod_pp(smdFile, tempdir='temps')
        sprod_run(sprod_dir='temps', out_dir='temp_result', pythonPath=cfg['pythonPath'])
        sprod_save(smdFile, out_dir='temp_result')
    elif method == 'none':
        print("Decontamination method: none.")
        return


# we provide 7 inline methods for clustering and supports Squidpy and Giotto
def Pip_cluster(smdFile, cfg, h5data, method='leiden', force=False, tif=None):
    method_list = ['leiden', 'SpaGCN', 'stlearn', 'Seurat', 'BayesSpace', 'Giotto', 'Squidpy', 'SEDR', 'MENDER', "BASS"]
    platform = cfg['platform']
    hires = cfg['load_hires']
    with h5.File(smdFile, 'r') as f:
        if method not in method_list:
            print("Clustering method:{0} has not been supported yet, default `leiden`.".format(method))
            method = 'leiden'
        if (h5data + '/idents/' + method in f) and (force == False):
            print("Clustering with " + method + " has done, skip it.")
            return
        if method == "BASS":
            if (f"{h5data}/idents/BASS_c" in f) and (force == False) and (f"{h5data}/idents/BASS_z" in f):
                print(f"Clustering with {method} has done, skip it.")
                return
    if method == 'leiden':
        from clu.Cluster_leiden import Spatial_Cluster_Analysis
        adata = Load_smd_to_AnnData(smdFile, h5data)
        adata = Remove_mito_genes(adata)
        adata = Filter_genes(adata)
        adata = Spatial_Cluster_Analysis(adata)
        Save_smd_from_leiden(smdFile, adata, h5data)
        # adata = Differential_Expression_Analysis(adata)
    elif method == 'MENDER':
        from clu.Cluster_MENDER import MENDER_install, Spatial_Cluster_MENDER
        MENDER_install()
        adata = Load_smd_to_AnnData(smdFile, h5data)
        adata = Remove_mito_genes(adata)
        adata = Filter_genes(adata)
        adata = Spatial_Cluster_MENDER(adata)
        Save_smd_from_leiden(smdFile, adata, h5data)
        Save_smd_from_MENDER(smdFile, adata, h5data)
    elif method == 'SpaGCN':
        from clu.Cluster_SpaGCN import spagcn_install
        spagcn_install()
        from clu.Cluster_SpaGCN import Spatial_Cluster_Analysis_SpaGCN
        adata = Load_smd_to_AnnData(smdFile, h5data)
        adata = Spatial_Cluster_Analysis_SpaGCN(adata)
        Save_smd_from_SpaGCN_clu(smdFile, adata, h5data)
        # adata = Differential_Expression_Analysis(adata)
    elif method == "SEDR":
        from clu.Cluster_SEDR import sedr_install
        sedr_install()
        from clu.Cluster_SEDR import sedr_run
        adata = Load_smd_to_AnnData(smdFile, h5data)
        adata = sedr_run(adata, n_clusters=20)
        Save_smd_from_SEDR(smdFile, adata, h5data)
    elif method == 'BayesSpace':
        # Using R scripts
        os.system("Rscript clu/Cluster_BayesSpace.R")
    elif method == 'BASS':
        # Using R scripts
        os.system("Rscript clu/Cluster_BASS.R")
    elif method == "Seurat":
        # Using R scripts
        os.system("Rscript clu/Cluster_Seurat.R")
    elif method == "stlearn":
        if platform in ["10X", "10X_Visium", "Visium", "visium", "10X Visium"]:
            from clu.Cluster_stlearn import Spatial_Cluster_Analysis_stSME
            stdata = Load_smd_to_stData(smdFile, count=h5data, hires=hires)
            stdata = Preprocess_stLearn(stdata)
            stdata = Spatial_Cluster_Analysis_stSME(stdata)
            Save_smd_from_stlearn(smdFile, stdata, h5data)
    elif method == "Giotto":
        os.system("Rscript clu/Cluster_Giotto.R")
    elif method == "Squidpy":
        from clu.Cluster_leiden import Spatial_Cluster_Analysis
        from clu.Cluster_Squidpy import Spatial_Cluster_Analysis_squidpy
        adata = Load_smd_to_AnnData(smdFile, h5data)
        if tif is None:
            print("Squidpy needs input image for cluster.Default leiden.")
            adata = Remove_mito_genes(adata)
            adata = Filter_genes(adata)
            adata = Spatial_Cluster_Analysis(adata)
            Save_smd_from_leiden(smdFile, adata, h5data)
        else:
            adata = Spatial_Cluster_Analysis_squidpy(adata, tif)
            Save_smd_from_Squidpy_clu(smdFile, adata, h5data)


def Pip_deconv(smdFile, cfg, h5data='matrix', method="stereoScope", force=False, standard_size=10):
    method_list = ['CARD', 'Cell2Location', 'SPOTlight', 'SPOTlight_vae', 'SPOTlight_es', 'spacexr', 'spacexr_es',
                   'StereoScope', 'Seurat', 'Seurat_es', 'Giotto', 'StereoScope_es', 'StereoScope_na']
    pythonPath = cfg['pythonPath']
    with h5.File(smdFile, 'r') as f:
        if method not in method_list:
            print("Deconvolution method:{0} has not been supported yet, default `stereoScope`.".format(method))
            method = 'stereoScope'
        if (h5data + '/deconv/' + method in f) and (force == False):
            print("Deconvolution with " + method + " has done, skip it.")
            return
    if method == 'spacexr':
        do_vae = True
        with h5.File(smdFile, 'r') as f:
            if ('sc-ref' in f) and (force == False):
                do_vae = False
        if do_vae:
            from dcv.spacexr import spacexr_pp_VAE
            spacexr_pp_VAE(smdFile, standard_size=standard_size)
        os.system("Rscript dcv/Deconv_spacexr.R")
    elif method == 'spacexr_es':
        from dcv.spacexr import spacexr_pp_EasySample
        spacexr_pp_EasySample(smdFile)
        os.system("Rscript dcv/Deconv_spacexr_es.R")
    elif method == 'CARD':
        os.system("Rscript dcv/Deconv_CARD.R")
    elif method == 'SPOTlight':
        # Using R scripts
        os.system("Rscript dcv/Deconv_SPOTlight.R")
    elif method == 'SPOTlight_es':
        do_es = True
        with h5.File(smdFile, 'r') as f:
            if ('sc-ref-es' in f) and (force == False):
                do_es = False
        if do_es:
            from dcv.spacexr import spacexr_pp_EasySample
            spacexr_pp_EasySample(smdFile, standard_size=standard_size)
        # Using EasyEnsample for SPOTlight reference
        os.system("Rscript dcv/Deconv_SPOTlight_es.R")
    elif method == 'SPOTlight_vae':
        do_vae = True
        with h5.File(smdFile, 'r') as f:
            if ('sc-ref' in f) and (force == False):
                do_vae = False
        if do_vae:
            from dcv.spacexr import spacexr_pp_VAE
            spacexr_pp_VAE(smdFile, standard_size=standard_size)
        # Using VAE for SPOTLight reference
        os.system("Rscript dcv/Deconv_SPOTlight_vae.R")
    elif method == 'StereoScope':
        from dcv.stereoScope import stereoscope_install
        stereoscope_install()
        from dcv.stereoScope import StereoScope_pp_VAE, StereoScope_run, StereoScope_pp_EasySample
        StereoScope_pp_VAE(smdFile, tempdir="temps", h5data=h5data,  standard_size=standard_size)
        StereoScope_run(stereo_dir="temps/stsc_temp", pythonPath=pythonPath, out_dir="temps/stsc_res")
        Wfile = os.listdir("temps/stsc_res/spt_data")[0]
        Wfile = "temps/stsc_res/spt_data/" + Wfile
        Save_smd_from_StereoScope(smdFile, Wfile, h5data, name='StereoScope')
        os.remove(Wfile)
    elif method == 'StereoScope_es':
        from dcv.stereoScope import stereoscope_install
        stereoscope_install()
        # stereoScope_es represents using easy sample in single-cell data
        from dcv.stereoScope import StereoScope_pp_VAE, StereoScope_run, StereoScope_pp_EasySample
        StereoScope_pp_EasySample(smdFile, h5data=h5data, standard_size=standard_size)
        StereoScope_run(stereo_dir="temps/stsc_temp", pythonPath=pythonPath, out_dir="temps/stsc_res")
        Wfile = os.listdir("temps/stsc_res/spt_data")[0]
        Wfile = "temps/stsc_res/spt_data/" + Wfile
        Save_smd_from_StereoScope(smdFile, Wfile, h5data, name='StereoScope_es')
        os.remove(Wfile)
    elif method == 'StereoScope_na':
        from dcv.stereoScope import stereoscope_install
        stereoscope_install()
        # stereoScope_na represents no preprocessing in single-cell data
        from dcv.stereoScope import StereoScope_pp_na, StereoScope_run
        StereoScope_pp_na(smdFile, h5data=h5data)
        StereoScope_run(stereo_dir="temps/stsc_temp", pythonPath=pythonPath, out_dir="temps/stsc_res")
        Wfile = os.listdir("temps/stsc_res/spt_data")[0]
        Wfile = "temps/stsc_res/spt_data/" + Wfile
        Save_smd_from_StereoScope(smdFile, Wfile, h5data, name='StereoScope_na')
        os.remove(Wfile)
    elif method == 'Cell2Location':
        from dcv.cell2location import cell2location_install
        cell2location_install()
        from dcv.cell2location import Cell2Location_run
        try:
            adata = Cell2Location_run(smdFile, h5data)
            Save_smd_from_Cell2Location(smdFile, adata, h5data)
        except:
            print("The distribution of data may not suitable for Cell2Location.")
    elif method == 'Seurat':
        state = subprocess.call(["Rscript", "dcv/Deconv_Seurat.R"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if state == 0:
            print("Deconvolution using Seurat finished.")
    elif method == 'Giotto':
        os.system("Rscript dcv/Deconv_Giotto.R")


def Despot_Decont(smdFile, cfg, force=False):
    methods = cfg['Decontamination']
    for method in methods:
        Pip_decont(smdFile, method=method, force=force, cfg=cfg)


def Despot_Cluster(smdFile, cfg, force=False):
    # if you need Decontamination
    h5datas = Create_h5datas(cfg)

    # decode dataPath
    dataPath = cfg['dataPath']
    methods = cfg['Clustering']
    tif = dataPath + '/' + cfg['imgPath'] + '/' + cfg['fullres_path']
    for h5data in h5datas:
        for method in methods:
            Pip_cluster(smdFile, cfg, h5data, method=method, tif=tif, force=force)


def Despot_Deconv(smdFile, cfg=None, force=False):
    sptinfo = smdInfo(smdFile)
    if cfg is None:
        cfg = sptinfo.configs
    do_scRNA_seq = False
    with h5.File(smdFile, 'r') as f:
        if 'scRNA_seq' not in f or force is True:
            print("scRNA_seq data not in smdFile, Generating it...")
            do_scRNA_seq = True
        else:
            print("Found scRNA_seq data in smdFile.")
    if do_scRNA_seq:
        print("Using scRNA-seq provided by users...")
        state = subprocess.check_call(["Rscript", "sc/Generate_scRNA-seq.R"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if state == 0:
            print("Generating scRNA-seq data finished.")
    # whether need Decontamination
    h5datas = Create_h5datas(cfg)

    methods = cfg['Deconvolution']
    f = open("log.txt", 'a')
    for method in methods:
        for h5data in h5datas:
            if h5data != "matrix" and method == "Cell2Location":
                print("Denoising by SPROD doesn't follow the GammaPoisson suppose in "
                      "Cell2Location, skip it.")
                continue
            print("Using {0} to do deconvolution in {1}.".format(method, h5data))
            start = time.time()
            Pip_deconv(smdFile, cfg, h5data, method=method, force=force)
            end = time.time()
            print("method {0} using: {1}min.".format(method, (end - start) / 60), file=f)
    f.close()


# do deconvolution for certain subset
def Despot_SubsetDeconv(smdFile, h5data="matrix", clu_mtd=None, domain=None) -> None:
    adata = Load_smd_to_AnnData(smdFile, h5data)
    adata_sub = adata[adata.obs[clu_mtd] == domain]
    idx = list(adata_sub.obs_names)
    subset_name = clu_mtd + "_" + str(domain)
    Save_smd_from_Subset(smdFile, h5data=h5data, subset_name=clu_mtd + "_" + str(domain), refer_idx=idx, force=False)
    Save_smd_from_configs(smdFile, {"Decontamination": subset_name})
    Despot_Deconv(smdFile)


def Despot_Ensemble(smdFile, beta=1, greedy=1, spmat_subset=None, clu_subset=None, dcv_subset=None) -> None:
    Despot_Find_bestGroup(smdFile, beta, greedy, spmat_subset, clu_subset, dcv_subset)


def Despot_Embedding(smdFile: Union[str, List[str]], cfg: Union[dict, List[dict]]) -> None:
    from ebd.Embed_Harmony import Embed_Harmony
    h5datas = Create_h5datas(cfg)
    for h5data in h5datas:
        adata = Load_smd_to_AnnData(smdFile, h5data)
        # embedding by harmony and hpca
        adata = Embed_Harmony(adata)
        Save_smd_from_Harmony(smdFiles=smdFile, adata=adata, h5data=h5data, name="Harmony")


def Despot_Run(cfg):
    name = cfg['name']
    platform = cfg['platform']
    cfg_json = json.dumps(cfg)
    smdFile = cfg['smdFile']
    # set running configs
    Save_json_configs(cfg, "params.json")
    # need hires?
    print(f"=========Despot Info============")
    print("smdFile:{0}".format(smdFile))
    print("dataset name:{0}".format(name))
    print("platform:{0}".format(platform))
    print("Using hires img: {0}".format(cfg['load_hires']))
    print(f"=========Despot Start===========")
    SMD_init(smdFile=smdFile)
    Save_smd_from_configs(smdFile, items=cfg)
    Despot_Decont(smdFile, cfg)
    Despot_Cluster(smdFile, cfg)
    Despot_Deconv(smdFile, cfg)
    Despot_Ensemble(smdFile)
    print(f"=========Despot Finish==========")
    print("smdFile:{0}".format(smdFile))
    print("dataset name:{0}".format(name))
    print("platform:{0}".format(platform))
    print("Using hires img: {0}".format(cfg['load_hires']))
    info = smdInfo(smdFile)
    return info