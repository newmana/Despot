def getMoranHP(platform:str) -> float:
  dict = {
    "10X Visium": 0.5,
    "10X_Visium": 0.5,
    "ST": 0.5,
    "MERFISH": 0.05,
    "Slide-seq": 0.1,
    "Stereo-seq":0.1,
    "osmFISH": 0.05,
  }
  return dict[platform]