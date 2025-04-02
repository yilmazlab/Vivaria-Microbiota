#!/bin/sh
#!/bin/bash

# SLURM job configuration
#SBATCH --job-name="Halla2"          # Job name for SLURM scheduler
#SBATCH --nodes=1                    # Number of nodes to allocate
#SBATCH --time=24:00:00              # Maximum runtime for the job
#SBATCH --mem-per-cpu=96G            # Memory allocated per CPU
#SBATCH --partition=bdw              # Partition to submit the job

## Optional output and error file paths
##SBATCH --output=/path/to/outfile    # Uncomment to specify standard output file
##SBATCH --error=/path/to/errfile     # Uncomment to specify standard error file

# For array jobs
# Array job containing 100 tasks, running a maximum of 10 tasks at the same time
##SBATCH --array=1-100%10            # Uncomment to enable array jobs

#### Shell commands start below ####

# Set locale to avoid potential encoding issues
export LC_CTYPE=en_US.UTF-8          # Set character encoding to UTF-8
export LC_ALL=en_US.UTF-8            # Set locale to UTF-8

# Create and activate a Python virtual environment
python -m venv halla-env             # Create a virtual environment named 'halla-env'
source halla-env/bin/activate        # Activate the virtual environment

# Install required Python packages
pip install matplotlib fontTools     # Install necessary Python libraries

# Navigate to the working directory containing input data
cd "Library/CloudStorage/Dropbox/Ongoing Analysis/Wild-Type Mice vs Wild Mice Analysis/2022/Halla_cecum_metabolites_microbiota"

# Run HAllA analysis with specified parameters
halla -x cecum_microbiota.txt -y cecum_metabolites.txt -o Microbiota_vs_Metabolomics \
      -m spearman --fdr_alpha 0.05 --fnr_thresh 0.2 --cbar_label 'Pairwise Spearman'
# -x: Input file for microbiota data
# -y: Input file for metabolites data
# -o: Output directory for results
# -m: Correlation method (Spearman)
# --fdr_alpha: False discovery rate threshold
# --fnr_thresh: False negative rate threshold
# --cbar_label: Label for the color bar in visualizations

# Run HAllA diagnostics with 100 permutations
hallagnostic -i Microbiota_vs_Metabolomics/ -n 100
# -i: Input directory containing HAllA results
# -n: Number of permutations for diagnostics

# Generate a hallagram visualization for the first dataset
hallagram \
    -i Microbiota_vs_Metabolomics/ \
    --cbar_label 'Pairwise Spearman' \
    --x_dataset_label 'Metabolites' \
    --y_dataset_label 'Microbiota' \
    -n 50 \
    --block_border_width 1
# -i: Input directory containing HAllA results
# --cbar_label: Label for the color bar
# --x_dataset_label: Label for the x-axis dataset
# --y_dataset_label: Label for the y-axis dataset
# -n: Number of top associations to display
# --block_border_width: Width of block borders in the visualization

# Generate a hallagram visualization for the second dataset
hallagram \
    -i Cecum_metabolites_microbiota/ \
    --cbar_label 'Pairwise Spearman' \
    --x_dataset_label Metabolites \
    --y_dataset_label 'Microbial species' \
    --block_border_width 2
# Similar parameters as above, but for a different dataset