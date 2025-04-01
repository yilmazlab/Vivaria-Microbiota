
# On any computer

Please make sure you have R version 4.2.1 installed. You can then follow steps 2 and 3 as detailed for the cluster to initialize the renv. You will not need to load any modules and may instead renv globally in your R installation if you wish to do so. Be aware that fitting all reaction-specific models on a single machine is very strongly discouraged as each model takes around 2 minutes and 30 seconds to fit and they are 7604 such models to fit.

# On the cluster

Once you have reached step 6. of the root README file.

Make sure you have a Python and a R module loaded. R 4.2.1 was used to produce the results described in the manuscript.

1. ```scp``` the necessary files to your folder of choice on the scratch. You will need to ```scp```:
    - fit_model_scripts/renv.lock
    - fit_model_scripts/.Rprofile
    - fit_model_scripts/renv/activate.R
    - fit_model_scripts/fit_reaction_models.R
    - fit_model_scripts/create_sh_files.py
    - data/processed_files/abund_and_meta.csv
    - data/rxn_names.csv


2. initiate the renv in the scratch
    1. make sure you are loading the correct R module
    2. create a renv folder and move the script activate.R to this folder
    3. start R by typing ```R```: this should bootstrap the renv. you should see something like this:
        ```# Bootstrapping renv 0.17.3
            Downloading renv 0.17.3 ... OK (downloaded source)
            Installing renv 0.17.3 ... Done!
            Successfully installed and loaded renv 0.17.3.
            Project 'your_project' loaded. [renv 0.17.3]
            One or more packages recorded in the lockfile are not installed.
            Use `renv::status()` for more details.
        ```
    4. In the R console, type:\
         ```renv::restore()```  
         This should install the necessary packages with the correct version. This  can take up to 30-40 minutes if those packages are not already cached.
         Then quit the R console using ```q()```.
    5. Next you need to install the renv package in a local library so that Rscript is able to load the package and the renv. This is done via lines 8-9 of the ```fit_reaction_models.R```. For those commands to work, create a folder named ```Rlib_421``` (you can choose another name but make sure to update l. 8 accordingly). Go back inside the R console and install the ```renv``` package in the library directory you just created by running:\
        ```install.packages("renv", lib = "./Rlib_421/")```
    5. Quit the R console using ```q()```, the environment installation should now be complete.
    

3. Test the environment with a fast and small run of the script:\
    ```Rscript --vanilla fit_reaction_models.R -i abund_and_metadata.csv -c rxn_names.csv -f results_files -b 2 -e 5```
    Remove the generated results file before running the scripts on the entire dataset.

4. In practice one should be able to run a signle bash file and do without the ```create_sh.py``` script. However, in case your cluster does not support it, you can use this Python script along with the ```matser_submit.sh``` bash script to generate bash scripts sending single jobs to the cluster (be sure to adapt the submission command in ```master_submit.sh``` to your infrastructure). Note that reaction indices are encoded in ```create_sh.py```.
    1. create a ```submit2cluster``` folder and move the create_sh.py script in it
    2. move to the folder you just created and invoke the python script:\
        ```python create_sh.py```
    3. move back to the root script and call the master submit script

Fitting all the reaction-specific models parallelized in this manner takes about 2 hours and 30 minutes.

5. Copy your results file back to your local computer and run the analysis notebooks.

