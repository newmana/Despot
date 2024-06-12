# Despot

Cell-type-specific domains are the anatomical domains in spatially resolved transcriptome (SRT) tissues where particular cell types are enriched coincidentally. Precisely detecting these domains is essential for the biological understanding of tissue architectures and cell co-localizations. However, using existing computational methods is challenging to detect specific domains with low-proportion cell types, which are partly overlapped with or even inside other cell-type-specific domains. Here, we propose an approach called Despot that synthesizes segmentation and deconvolution as an ensemble to generate cell-type patterns, detect low-proportion cell-type specific domains, and display these domains intuitively. Experimental evaluation showed that Despot enabled us to discover the co-localizations between cancer-associated fibroblasts and immune-related cells that indicate potential tumor microenvironment (TME) domains in given slices, which was obscured by previous computational methods. We further elucidated the identified domains and found that Srgn may be a critical TME marker in SRT slices. By deciphering T-cell-specific domains in breast cancer tissues, Despot also revealed that the proportions of exhausted T cells were significantly more enormous in invasive than in ductal carcinoma. 

![overview](overview.png)

## Installation

Despot is implemented by Python 3.9.7 and R 4.1.3 (tested on Ubuntu 20.04, MacOS Catalina). 

**For Python 3.9:**

First we create a new virtual environment, and activate it. Then we install the basic python packages in `requirements.txt`

```shell
conda create --name Despot python=3.9
conda activate Despot
pip install -r requirements.txt
```

Time to install the requirements depends on the operating system.

**For R 4.1.3:**

Run this code in R Console for Despot basic requirements

```R
install.packages(c("BiocManager", "Matrix","stringr"))
BiocManager::install(c("rhdf5", "png", "rjson","data.table","SingleCellExperiment", "optparse", "distances", "gtools", "umap", "dplyr", "mvtnorm", "qvalue", "ComplexHeatmap", "tidyverse"))
```

## Usage

### Execution in Shell/Bash

We run the `main.py` to execute the configures in `/configs` using python:

```shell
python main.py
```

The required softwares will be installed automatically. Users also can install them manually using the URLs in `source.md `.  The smdFiles generated are stored in `h5smds`.

**Config Settings**

In default, we run the configures in `/configs`. Users can set their own configs with the guidance below:

```json
{
  "dataPath": "path/to/your/SRT data/input",
  "imgPath": "path/to/your/image/input under dataPath",
  "platform": "10X_Visium",  /*10X_Visium, ST, Slide-seq*/	
  "dataSpices": "Human",	/*Human, Mice*/
  "name": "your_sample_name",
  "filter_matrix": true,  /*Whether the matrix is filtered*/
  "smdFile": "path/to/output/smdFiles",
  "ground_truth": "metadata.csv", /*SRT ground truth file*/
  "ground_name": "ground_truth", /*SRT ground truth name*/
  "scDataPath": "path/to/your/single-cell data/input",
  "scType": "Any", /*Any, tsv, mtx, txt, h5ad*/
  "scBarcodes": "count_matrix_barcodes.tsv",
  "scFeatures": "count_matrix_genes.tsv",
  "scMatrix": "count_matrix_sparse.mtx",
  "scGroundTruth": "metadata.csv",
  "scGroundName": "celltype_minor",
  "load_hires": true,
  "load_fullres": false,
  "fullres_path": "path/to/your/fullres_datapath",
  "Decontamination": ["none","SPCS","SpotClean", "SPROD"...],
  "Clustering": ["BayesSpace", "SpaGCN", "leiden", "stlearn", "SEDR"...],
  "Deconvolution": ["CARD","Cell2Location","SPOTlight", "spacexr","StereoScope"...],
}
```

Despot will load these configs first, check for the environments, and run for outputs in smdFiles.

### Execution in Python console

`TODO`

## Description

### Cell-type Specific Domains

Conventional spatial domains are expected to have high intra-cluster similarity within each domain and low inter-cluster similarity between different domains. Cell-type-specific domains are further expected to be spatial heterogenous domains where the proportions of given cell types are significantly higher than other domains. In addition, cell-type-specific domains is suitable for SRT slices in both spot-level and single-cell resolution. Compared with conventional domains, cell-type-specific domains mainly focus on the location of potential cell types and relax restrictions on regional segmentation. Here is an example:



<img src="domain_example.png" alt="cts" style="zoom:33%;" />

Figure 1A, conventional spatial domains separate the tissue into $D_1$ , $D_2$, $D_3$ and $D_4$. Each pair of them has no overlaps. Figure 1B, the specific domains of cell-type $c_1$, $c_2$, $c_3$, $c_4$ and $c_5$ are identified as $D_1^*$, $D_2^*$, $D_3^*$, $D_4^*$ and $D_5^*$. In particular, the $D_1^*$ specific domain is totally overlapped with D2, and they are inside the  $D_3^*$ specific domain. The $D_4^*$ specific domain is partly overlapped $D_5^*$, and their overlapped domain is marked by yellow.



### Despot Algorithm

Despot aims to precisely detect cell-type specific domains in sequencing-based spatial transcriptomics data via ensemble learning. We developed a data representation called *smd* as the core data structure of Despot. In addition to spatial modalities analysis, Despot provides a variational inference module to derive single-cell patterns. Despot integrates several segmentation and deconvolution methods as an ensemble to detect cell-type-specific domains in the modality integration module. After detecting cell-type-specific domains, Despot examines the upregulated genes of these domains, deciphers the cell-type co-localizations, and finally visualizes them using 3D landscapes.



### Data Respresentation

Despot designs a spatial multi-omics data representation format called *smd* based on HDF5 binary data. A *smd* file has interfaces with Python and R infrastructures, containing *spContexts* and a series of *scMatrixs* and *spMatrixs*, whose core matrices are column compressed. 

![smd](smd.png)

### 3D Landscape

Despot designs 3D landscapes to display the detected cell-type-specific domains. The generated landscapes curve the outlines of cell-type-specific domains in tissues and display the overlapped domains in layered structures. Despot uses the Alpha-shape algorithm to curve the contour lines in tissues and uses dots with colors to represent the layers of cell-type-specific domains. To map the layered structures onto contour lines, Despot randomly selects dots of 10 percent for each layer and creates waterfall lines to link the dots and outlined regions together. 3DLandscapes can be displayed using the function `Show_3DLandscape`, which receives an output SMD file from Despot. Here is an example:

<img src="3DLandscape.png" alt="3DLandscape" style="zoom: 33%;" />

**Figure 3**. A demo of a 3D Landscape, corresponding to Figure S1C. Key elements like cell-type specific domains, slice images, waterfall lines, contour lines, and co-localization statuses are marked out by red arrays. Slice Image is at the bottom of the 3D Landscape, contour lines are drawn in the slice image. Cell-type-specific domains are plotted in hierarchical layers, which are linked with contours in the slice image by the waterfall lines. The 3D landscape intuitively displays the cell co-localization status. In particular, domain brown and domain blue are totally overlapped, domain blue is inside the domain yellow, domain purple, and domain green are partly overlapped, and domain pink is adjacent to domain purple. 

## APIs

`TODO`

