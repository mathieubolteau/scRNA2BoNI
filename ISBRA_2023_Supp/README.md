# Supplementary Material: Inferring Boolean Networks from Single-Cell Human Embryo Datasets

This folder contains supplementary information for the paper "Inferring Boolean networks from single-cell human embryo datasets" of ISBRA 2023 conference.



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

### Raw results
- Method raw results archive

## `scRNAseq_data` folder

Here, we provide the scRNAseq data as an expression matrix. The first matrix (`scRNAseq_data.zip`) is the initial matrix composed of the four datasets from Meistermann et al. The second one (`Petropoulos2016_scRNAseq_data.zip`) is only composed of the Petropoulos2016 dataset. The latter was used for the paper.

## `supplementary_material` folder
We provide supplementary methods and results complementary to the main paper.

## Contact

Mathieu Bolteau 
* [@mathieubolteau](https://github.com/mathieubolteau)
* [mathieu.bolteau@ls2n.fr](mailto:mathieu.bolteau@ls2n.fr)
* Nantes Université, Ecole Centrale Nantes, CNRS, LS2N, UMR 6004, F-44000, Nantes