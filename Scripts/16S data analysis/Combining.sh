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

singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime feature-table merge  \
--i-tables table_2013v1.qza \
--i-tables table_2013v2.qza \
--i-tables table_2015.qza \
--i-tables table_2022.qza \
--i-tables table_GutEco1.qza \
--i-tables table_GutEcorep.qza \
--i-tables table_PatientsYoungSPF.qza \
--i-tables table_Chip68_126.qza \
--o-merged-table table_Combined_2022.qza

singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime feature-table merge-seqs   \
--i-data rep-seqs_2013v1.qza \
--i-data rep-seqs_2013v2.qza \
--i-data rep-seqs_2015.qza \
--i-data rep-seqs_2022.qza \
--i-data rep-seqs_GutEco1.qza \
--i-data rep-seqs_GutEcorep.qza \
--i-data rep-seqs_PatientsYoungSPF.qza \
--i-data rep-seqs_Chip68_126.qza \
--o-merged-data rep-seqs_Combined_2022.qza

singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime feature-table summarize \
--i-table table_Combined_2022.qza \
--o-visualization table_Combined_2022.qzv \
--m-sample-metadata-file After_Chip126_2022.tsv

singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime feature-table tabulate-seqs \
--i-data rep-seqs_Combined_2022.qza \
--o-visualization rep-seqs_Combined_2022.qzv

singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences rep-seqs_Combined_2022.qza  \
--o-alignment aligned-rep-seqs_Combined_2022.qza \
--o-masked-alignment masked-aligned-rep-seqs_Combined_2022.qza  \
--o-tree unrooted-tree_Combined_2022.qza  \
--o-rooted-tree rooted-tree_Combined_2022.qza

singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime feature-classifier classify-sklearn \
--i-classifier gg-13-8-99-nb-classifier.qza  \
--i-reads rep-seqs_Combined_2022.qza \
--o-classification taxonomy_Combined_2022.qza

singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simgqiime metadata tabulate \
--m-input-file taxonomy_Combined_2022.qza \
--o-visualization taxonomy_Combined_2022.qzv

singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simgqiime taxa barplot \
          --i-table table_Combined_2022.qza \
          --i-taxonomy taxonomy_Combined_2022.qza \
          --m-metadata-file After_Combined_2022.tsv \
          --o-visualization taxa-bar-plots_Combined_2022.qzv
