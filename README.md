# DESpace2_manuscript
This repository contains the code used for the analyses presented in the manuscript “DESpace2: detection of differential spatial patterns in spatial omics data”.

The *./Analyses* folder contains the following subdirectories:

- `01_preprocessing`: quality control and filtering;
- `real_data_analyses_ARTISTA`: analyses of the ARTISTA real dataset;
- `sensitivity_analyses_ARTISTA`: sensitivity analyses of the ARTISTA simulation dataset;
- `simulation_ARTISTA` and `simulation_LIBD`: simulation studies based on the ARTISTA and LIBD datasets, respectively. 

Within each simulation subfolder, scripts are organized to first generate the simulated datasets, then apply clustering methods (Banksy or BayesSpace), and finally run all methods. `highly_abundant` or `lowly_abundant` refer to scenarios in which randomly selected spatial domain(s) show higher or lower abundance relative to the rest of the tissue.

The *./Figures* folder contains the code used to generate all figures included in the manuscript and Supplementary Information. Input data are available on [Zenodo](https://zenodo.org/records/15744023).
