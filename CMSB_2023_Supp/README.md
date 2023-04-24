# Bolteau et al. CMSB 2023 Supplementary

This folder contains supplementary information for the paper "Inferring Boolean networks from single-cell human embryo datasets: proof of concept with trophectoderm maturation" submitted for the CMSB 2023 conference.


## `figs` folder
In this folder, we provide the figures present in the paper.


## `medium_late_TE_discrimination` folder

This folder contains the results obtained for the discrimination of the medium and late
trophectoderm stages. We divide the folder into three subfolders,  each corresponding to each step of the framework. 

### PKN reconstruction step
- List of the transcription factor used as input for pyBRAvo
- The reconstructed PKN from pyBRAvo: Cytoscape visualization and PDF
- The reduced PKN: Cytoscape visualization and PDF

### Experimental design reconstruction
- Experimental design setup containing input, intermediate and readout genes
- MIDAS file for medium and late TE used by Caspo for the BNs inference step

### BNs inference
- Inferred BNs for medium and late TE

## `scRNAseq_data` folder

Here, we provide the scRNAseq data as an expression matrix. The first matrix (`scRNAseq_data.zip`) is the initial matrix composed of the four datasets from Meistermann et al. The second one (`Petropoulos2016_scRNAseq_data.zip`) is only composed of the Petropoulos2016 dataset. The latter was used for the paper.