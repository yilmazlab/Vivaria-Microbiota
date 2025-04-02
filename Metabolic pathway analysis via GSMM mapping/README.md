# Environments

- The python environment used is specified in python_files/requirements.txt. This environment can be installed by running
``` pip install -r requirements.txt ```
- The renv for processing the initial data can be found in renv_process_data.zip
- The renv used when fitting the models can be found in the fit_model_scripts folder\\

Both those renv need the installation of the ```renv``` R package to function. To install a renv, copy the ```activate.R``` file into a folder named renv. In the parent folder, copy the ```renv.lock```. Provided the ```renv``` package isinstalled, running R in the parent folder should automatically bootstrap the renv. Running ```renv::restore()``` will then install the necessary packages with the correct versions.

- MATLAB scripts were run with MATLAB R2022b

# Data

- Raw data from the analysis can be obtained here:
- The AGORA2 models and Supplementary Table S1 of the Heinken et al., 2023 (DOI: 10.1038/s41587-022-01628-0) should be available in the ```data``` folder in order to run steps 2 and 3 of the analysis. See ```matlab_files/README.md``` for detailed instructions on how to obtain and set up those files.

# Running the code

### Data folders

The scripts assume that the data is present in a ```data``` folder within the repo. All references to file in the current code should be relative to this folder and should thus be portable between systems, with the exception of the results obtained when fitting the models.

```
data\
|
|_AGORA2_models\
|_ original_files\
|    |   contains the unprocessed data files.
|    |_class_analysis_16_10\
|    |   contains the additional unprocessed data files from the supervised cluster  in Figure 3
|_processed_files\
|    |   contains the files derived from the data files used to generate the inputs to model fitting.
|    |_Classes_October_16\
|    |   contains the derived files from the supervised cluster analysis from Figure 3
|_results_files\
|    |   contains the files obtained after fitting the models
|_subsystem_assignation\
|    |  contains the output files from the subsystem assignation scripts
|_tables\
|    |   generated tables will be saved here
|_figures\
|   |   generated figures will be saved here
```

### Main scripts and pipeline

1. Generate the taxonomy and otu files ```taxonomy.csv, otumat.mat, otumat.rds```.
- ```R_files/process_original_rds.RMD```
- *Note: We could not separately intialize a renv for this file but package specifications can be found in the renv_process_data.zip folder. The renv from fit_model_script does not contain packages for the R_files files due to compatibility problems when both renv specifications were present in the project.*

2. Generate the normalized reaction abundance
- ```matlab_files/microbiome.m```
- *Note: See ```matlab_files/README.md``` for instructions on how to run this script*

3. Exploratory analysis of the initial data  (in the folder ```python_files/exploratory_analysis```):
- ```process_metadata.ipynb```: match the given metadata with the samples in the dataset and curate the chip annotation. The matching is necessary as the metadata file was obtained after the samples and contains irrelevant rows. The chip annotation is curated by replacing each chip number Chip117vX by Chip117.
- ```mapped_species_AGORA2.ipynb```: evaluate the proportion of the total reads which can be mapped to a model. Generates the data for Extended Figure 6.d. Computes the same mapping as the one done in the MATLAB script.
- ```studying_class_imbalance.ipynb```: explore the repartition of the data between groups 


4. Generate the subsystem mapping
- ```matlab_files/subsystem_annotation/fAnnotateReactions.m```: generate the mapping
- *Note: Subsystem mapping was obtained on June 29th, 2023. Mappings may change with time as the script pulls information from online databases.*
- ```python_files/parse_subsystem_assignation.ipynb```: process the mapping for use in downstream analysis.

5. Prepare the input files for model fitting
- ```python_files/create_input_dataset.ipynb```

6. Fit the models
- Scripts in ```R_files/fit_reaction models```  
- Instructions on how to create the renv and run this script are given in ```fit_model_scripts/README.md```
- *Note: each model takes on average 2.30 mn to fit and thre 7604 reaction-specific models.  This script was thus run in parallel on a cluster. You may need to adapt it to your own infrastructure*

7. Process multiple cluster outputs (if run in parallel)
- ```python_files/results_files/process_cluster_output.ipynb```

8. Analyze results 
- ```python_files/results_files/result_Human_SPF_Wild.ipynb``` generate the data for Figure 3.d and Extended Figures 6f, 6g and 6h: 
- ```python_files/results_files/vivarium_ranef.ipynb``` investigates the significance of the random effect modelling the origin of the samples.

9. Generation of Supplementary Figure 6e, Supplementary Table 8, and additional tables and statistics
- ```python_files/supplementary_analyses/supplementary_analyses.ipynb```


