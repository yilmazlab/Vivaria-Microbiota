#!/bin/sh
#!/bin/bash
#SBATCH --mail-user=bahtiyar.yilmaz@dbmr.unibe.ch
#SBATCH --mail-type=end
#SBATCH --job-name="GutEco1"
#SBATCH --nodes=1
#SBATCH --time=6:00:00
#SBATCH --mem-per-cpu=32G
#SBATCH --partition=epyc2

##SBATCH --output=~/StomaMetagenomics/SlurmFiles
##SBATCH --error=~/StomaMetagenomics/SlurmFiles

####For array jobs
###Array job containing 10 tasks, run max 20 tasks at the same time
####SBATCH --array=1-43%20


## QIIME pipeline shell script ##

# singularity help /software/singularity/containers/qiime2-2018.8-1.debian9.simg
export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

cd /storage/homefs/terziev/Wild_Mice_Analysis/QIIME2_WildMice/Combined_2022

singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime metadata tabulate \
--m-input-file taxonomy_Combined_2022.qza \
--o-visualization taxonomy_Combined_2022.qzv

singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime taxa barplot \
          --i-table table_Combined_2022.qza \
          --i-taxonomy taxonomy_Combined_2022.qza \
          --m-metadata-file After_Chip126_2022.tsv \
          --o-visualization taxa-bar-plots_Combined_2022.qzv
