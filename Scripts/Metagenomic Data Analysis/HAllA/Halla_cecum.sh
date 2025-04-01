#!/bin/sh
#!/bin/bash
#SBATCH --job-name="Halla2"
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=96G
#SBATCH --partition=bdw

##SBATCH --output=/path/to/outfile
##SBATCH --error=/path/to/errfile

# For array jobs
# Array job containing 100 tasks, run max 10 tasks at the same time
##SBATCH --array=1-100%10

#### Your shell commands below this line ####

## QIIME pipeline shell script ##

export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

source activate halla
python -m venv halla-env
source halla-env/bin/activate
pip install matplotlib fontTools
cd "Library/CloudStorage/Dropbox/Ongoing Analysis/Wild-Type Mice vs Wild Mice Analysis/2022/Halla_cecum_metabolites_microbiota"

halla -x cecum_microbiota.txt -y cecum_metabolites.txt -o Microbiota_vs_Metabolomics -m spearman --fdr_alpha 0.05 --fnr_thresh 0.2 --cbar_label 'Pairwise Spearman' 
hallagnostic -i Microbiota_vs_Metabolomics/ -n 100
hallagram \
    -i Microbiota_vs_Metabolomics/ \
    --cbar_label 'Pairwise Spearman' \
    --x_dataset_label 'Metabolites' \
    --y_dataset_label 'Microbiota' \
    -n 50 \
    --block_border_width 1


hallagram \
    -i Cecum_metabolites_microbiota/ \
    --cbar_label 'Pairwise Spearman' \
    --x_dataset_label Metabolites \
    --y_dataset_label 'Microbial species' \
    --block_border_width 2