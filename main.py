from utils.check import *
print("Checking requirements...")
Check_Environments()
Check_Requirements({"anndata", "h5py","matplotlib", "numpy","collection", "pandas", "scanpy", "scipy", "torch", "rtree", "geopandas", "esda"})
print("Importing requirements...")
import json
import shutil
from utils.io import *
from utils.api import Despot_Decont, Despot_Cluster, Despot_Deconv, Despot_Ensemble
print("All requirements imported...")


def Replicate_run(cfg:dict, iterations=1):
    smdFile = cfg['smdFile']
    cfg['smdFile'] = f"{smdFile.split('.')[0]}-{iterations}.h5smd"
    return cfg

if __name__ == "__main__":
    cfg_list = os.listdir('configs')
    for cfg_name in cfg_list:
        if cfg_name == 'V1_mouse_brain.json':
            cfg_path = f'configs/{cfg_name}'
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
            smdFile0 = cfg['smdFile']
            for i in range(5):
                name = cfg['name']
                platform = cfg['platform']
                cfg['smdFile'] = f"{smdFile0.split('.')[0]}-{i}.h5smd"
                cfg_json = json.dumps(cfg)
                smdFile = cfg['smdFile']
                # set running configs
                Save_json_configs(cfg, "params.json")
                # need hires?
                hires = cfg['load_hires']
                print(f"=========Despot {i} Info============")
                print("smdFile:{0}".format(smdFile))
                print("dataset name:{0}".format(name))
                print("platform:{0}".format(platform))
                print("Using hires img: {0}".format(hires))
                print(f"=========Despot {i} Start===========")
                SMD_init(smdFile=smdFile)
                Save_smd_from_configs(smdFile, items=cfg)
                Despot_Decont(smdFile, cfg)
                Despot_Cluster(smdFile, cfg)
                Despot_Deconv(smdFile, cfg)
                Despot_Ensemble(smdFile)
                print(f"=========Despot {i} Finish==========")