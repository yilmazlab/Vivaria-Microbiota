#!/bin/sh
#!/bin/bash
#SBATCH --mail-user=bahtiyar.yilmaz@dbmr.unibe.ch  # Email for job notifications
#SBATCH --mail-type=end  # Notify when the job ends
#SBATCH --job-name="GutEco1"  # Job name
#SBATCH --nodes=1  # Number of nodes
#SBATCH --time=6:00:00  # Maximum runtime
#SBATCH --mem-per-cpu=32G  # Memory per CPU
#SBATCH --partition=epyc2  # Partition to run the job

##SBATCH --output=~/StomaMetagenomics/SlurmFiles  # Uncomment to specify output file
##SBATCH --error=~/StomaMetagenomics/SlurmFiles  # Uncomment to specify error file

## QIIME pipeline shell script ##
export LC_CTYPE=en_US.UTF-8  # Set locale for character encoding
export LC_ALL=en_US.UTF-8  # Set locale for all settings

# Navigate to the working directory
cd /storage/homefs/terziev/Wild_Mice_Analysis/QIIME2_WildMice/Combined_data

# Merge feature tables from multiple datasets
singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime feature-table merge  \
--i-tables table_2013v1.qza \
--i-tables table_2013v2.qza \
--i-tables table_2015.qza \
--i-tables table_data.qza \
--i-tables table_GutEco1.qza \
--i-tables table_GutEcorep.qza \
--i-tables table_PatientsYoungSPF.qza \
--i-tables table_Chip68_126.qza \
--o-merged-table table_Combined_data.qza

# Merge representative sequences from multiple datasets
singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime feature-table merge-seqs   \
--i-data rep-seqs_2013v1.qza \
--i-data rep-seqs_2013v2.qza \
--i-data rep-seqs_2015.qza \
--i-data rep-seqs_data.qza \
--i-data rep-seqs_GutEco1.qza \
--i-data rep-seqs_GutEcorep.qza \
--i-data rep-seqs_PatientsYoungSPF.qza \
--i-data rep-seqs_Chip68_126.qza \
--o-merged-data rep-seqs_Combined_data.qza

# Summarize the merged feature table
singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime feature-table summarize \
--i-table table_Combined_data.qza \
--o-visualization table_Combined_data.qzv \
--m-sample-metadata-file master_mappingfile.tsv

# Tabulate the merged representative sequences
singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime feature-table tabulate-seqs \
--i-data rep-seqs_Combined_data.qza \
--o-visualization rep-seqs_Combined_data.qzv

# Generate a phylogenetic tree from the merged representative sequences
singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences rep-seqs_Combined_data.qza  \
--o-alignment aligned-rep-seqs_Combined_data.qza \
--o-masked-alignment masked-aligned-rep-seqs_Combined_data.qza  \
--o-tree unrooted-tree_Combined_data.qza  \
--o-rooted-tree rooted-tree_Combined_data.qza

# Classify taxonomy using a pre-trained classifier
singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime feature-classifier classify-sklearn \
--i-classifier silva-138-99-nb-classifier.qza  \
--i-reads rep-seqs_Combined_data.qza \
--o-classification taxonomy_Combined_data.qza

# Tabulate the taxonomy classification results
singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simgqiime metadata tabulate \
--m-input-file taxonomy_Combined_data.qza \
--o-visualization taxonomy_Combined_data.qzv

# Generate a taxa bar plot visualization
singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simgqiime taxa barplot \
          --i-table table_Combined_data.qza \
          --i-taxonomy taxonomy_Combined_data.qza \
          --m-metadata-file After_Combined_data.tsv \
          --o-visualization taxa-bar-plots_Combined_data.qzv

# Tabulate taxonomy classification results again (duplicate command, consider removing)
singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime metadata tabulate \
--m-input-file taxonomy_Combined_data.qza \
--o-visualization taxonomy_Combined_data.qzv

# Generate another taxa bar plot visualization (duplicate command, consider removing)
singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime taxa barplot \
          --i-table table_Combined_data.qza \
          --i-taxonomy taxonomy_Combined_data.qza \
          --m-metadata-file master_mappingfile.tsv \
          --o-visualization taxa-bar-plots_Combined_data.qzv
