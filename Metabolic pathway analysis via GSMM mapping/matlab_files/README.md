# File provenance

The scripts ```count_reac.m```, ```count_unique.m```, ```create_super_rxn_set.m```, ```weighted_union.m``` were obtained from the code generated in Swiss IBD Cohort Investigators *et al.* , 2019. DOI: 10.1038/s41591-018-0308-z. 
They can be found at https://figshare.com/s/a34f96698ca6fcd36ac2.

The file ```microbiome.m``` was heavily adapted from the file ```macpherson_analysis.m``` of the same source. In order to run this script, it is assumed that the ```data``` folder at the root of this repo contains a folder named AGORA2_models organized as follows:
```
data/
|_...
|_ AGORA2_models/
|  |_AGORA2_TableS1.xlsx
|  |_mat/
```
where ```mat``` contains the MATLAB files for the AGORA2 models of interest. AGORA2_TableS1.xlsx corresponds to Table S1 of Heinken et al., 2023 (DOI: 10.1038/s41587-022-01628-0) which has been extracted as a sheet out of the Excel file containing the supplementary material, and for which the first three rows were removed in order to obtain a dataframe-like structure.

The files in the folder 'subsystem_annotation' were written by Joerg Stelling to obtain map each reaction to a metabolic subsystem. We used the level 3 of the KEGG mapping (see ```python_files/parse_subsystem_assignation.ipynb```) and mapping was obtained on June 29th, 2023. The script ```fAnnotateReactions.m``` also needs to be provided with the path to the MATLAB AGORA2 models, currently set to ```../data/AGORA2_models/mat/```. The results of the subsystem annotation script are currently saved in ```../data/subsystem_assignation/``` provided such a folder exists.