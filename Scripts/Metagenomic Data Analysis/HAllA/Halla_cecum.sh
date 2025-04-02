#!/bin/sh
#!/bin/bash

# SLURM job configuration
#SBATCH --job-name="Halla2"          # Job name
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --time=24:00:00              # Maximum runtime
#SBATCH --mem-per-cpu=96G            # Memory per CPU
#SBATCH --partition=bdw              # Partition to submit the job

## Optional output and error file paths
##SBATCH --output=/path/to/outfile    # Path to save standard output
##SBATCH --error=/path/to/errfile     # Path to save standard error

# For array jobs
# Array job containing 100 tasks, running a maximum of 10 tasks at the same time
##SBATCH --array=1-100%10

#### Shell commands start below ####

# Set locale to avoid potential encoding issues
export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

# Create and activate a Python virtual environment
python -m venv halla-env
source halla-env/bin/activate

# Install required Python packages
pip install matplotlib fontTools

# Navigate to the working directory containing input data
cd "Library/CloudStorage/Dropbox/Ongoing Analysis/Wild-Type Mice vs Wild Mice Analysis/2022/Halla_cecum_metabolites_microbiota"

# Run HAllA analysis with specified parameters
halla -x cecum_microbiota.txt -y cecum_metabolites.txt -o Microbiota_vs_Metabolomics \
      -m spearman --fdr_alpha 0.05 --fnr_thresh 0.2 --cbar_label 'Pairwise Spearman'

# Run HAllA diagnostics with 100 permutations
hallagnostic -i Microbiota_vs_Metabolomics/ -n 100

# Generate a hallagram visualization for the first dataset
hallagram \
    -i Microbiota_vs_Metabolomics/ \
    --cbar_label 'Pairwise Spearman' \
    --x_dataset_label 'Metabolites' \
    --y_dataset_label 'Microbiota' \
    -n 50 \
    --block_border_width 1

# Generate a hallagram visualization for the second dataset
hallagram \
    -i Cecum_metabolites_microbiota/ \
    --cbar_label 'Pairwise Spearman' \
    --x_dataset_label Metabolites \
    --y_dataset_label 'Microbial species' \
    --block_border_width 2