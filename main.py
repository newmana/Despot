from utils.check import *
print("Checking requirements...")
Check_Environments()
Check_Requirements({"anndata", "h5py","matplotlib", "numpy", "pandas", "scanpy", "scipy", "torch", "rtree", "geopandas", "esda"})
print("Importing requirements...")
import json
import shutil
from utils.io import *
from utils.api import Create_h5datas, Despot_Decont, Despot_Cluster, Despot_Deconv, Despot_Ensemble
from utils.geo import Despot_Find_bestGroup, Despot_group_correlation,\
    Show_self_correlation, Show_bash_best_group, Gen_venn
print("All requirements imported...")


def SMD_init(smdFile: str, force: bool = False):
    # whether the file exists
    if os.path.isfile(smdFile) and force is False:
        # check the corrections of sptFile
        if SMD_check_corrections(smdFile):
            return

    print("Initializing smdFile...")
    cmd = "Rscript utils/Init.R"
    os.system(cmd)
    print("Done.")


def Replicate_run(cfg:dict, iterations=1):
    smdFile = cfg['smdFile']
    cfg['smdFile'] = f"{smdFile.split('.')[0]}-{iterations}.h5smd"
    return cfg

if __name__ == "__main__":
    cfg_list = os.listdir('configs')
    for cfg_name in cfg_list:
        if cfg_name == 'osmFISH.json':
            cfg_path = 'configs/' + cfg_name
            shutil.copy(cfg_path, dst="params.json")
            cfg = Load_json_configs("params.json")

            # set python path
            pythonPath = sys.executable
            cfg["pythonPath"] = pythonPath

            # set R lib path
            cfg["R_library_Path"] = 'sptranr/R'

            # set working dictionary
            cfg["working_dir"] = os.path.dirname(os.path.abspath(__file__))

            # set venv
            cfg["venv"] = sys.prefix.split('/')[-1]

            for i in range(1):
                smdFile0 = cfg['smdFile']
                name = cfg['name']
                platform = cfg['platform']
                cfg['smdFile'] = f"{smdFile0.split('.')[0]}-{i}.h5smd"
                cfg_json = json.dumps(cfg)
                smdFile = cfg['smdFile']
                # set running configs
                Save_json_configs(cfg, "params.json")
                Save_smd_from_configs(smdFile, items=cfg)
                # need hires?
                hires = cfg['load_hires']
                print(f"=========Despot {i} Info============")
                print("smdFile:{0}".format(smdFile))
                print("dataset name:{0}".format(name))
                print("platform:{0}".format(platform))
                print("Using hires img: {0}".format(hires))
                print(f"=========Despot {i} Start===========")
                SMD_init(smdFile=smdFile)
                Despot_Decont(smdFile, cfg)
                Despot_Cluster(smdFile, cfg)
                Despot_Deconv(smdFile, cfg)
                Despot_Ensemble(smdFile)
                print(f"=========Despot {i} Finish==========")