#!/bin/sh
#!/bin/bash
#SBATCH --job-name="Chip147"
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=64G
#SBATCH --partition=all

##SBATCH --output=/path/to/outfile
##SBATCH --error=/path/to/errfile

# For array jobs
# Array job containing 100 tasks, run max 10 tasks at the same time
##SBATCH --array=1-100%10


#### Your shell commands below this line ####

## QIIME pipeline shell script ##

# singularity help /software/singularity/containers/qiime2-2018.8-1.debian9.simg

cd ~/Wild_Mice_Analysis/QIIME2_WildMice/Chip147/

singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime tools import \
--type MultiplexedSingleEndBarcodeInSequence \
--input-path Chip147_Diets.fastq.gz \
--output-path multiplexed-seqs_Chip147.qza

singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime cutadapt demux-single \
--i-seqs multiplexed-seqs_Chip147.qza \
--m-barcodes-file Chip147.tsv \
--m-barcodes-column BarcodesSequence \
--p-error-rate 0 \
--o-per-sample-sequences demultiplexed-seqs_Chip147.qza \
--o-untrimmed-sequences untrimmed_Chip147.qza \
--verbose

singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime cutadapt trim-single \
--i-demultiplexed-sequences demultiplexed-seqs_Chip147.qza \
--p-front ATTAGATACCCYGGTAGTCC \
--p-error-rate 0 \
--o-trimmed-sequences trimmed-seqs_Chip147.qza \
--verbose

singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg  qiime demux summarize \
--i-data trimmed-seqs_Chip147.qza \
--o-visualization trimmed-seqs_Chip147.qzv

singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime dada2 denoise-single \
--i-demultiplexed-seqs demultiplexed-seqs_Chip147.qza  \
--p-trunc-len 120  \
--o-table table_Chip147.qza \
--o-representative-sequences rep-seqs_Chip147.qza  \
--o-denoising-stats denoising-stats_Chip147.qza

singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime feature-table summarize \
--i-table table_Chip147.qza \
--o-visualization table_Chip147.qzv \
--m-sample-metadata-file Chip147.tsv

singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime feature-table tabulate-seqs \
--i-data rep-seqs_Chip147.qza \
--o-visualization rep-seqs_Chip147.qzv

singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences rep-seqs_Chip147.qza  \
--o-alignment aligned-rep-seqs_Chip147.qza \
--o-masked-alignment masked-aligned-rep-seqs_Chip147.qza  \
--o-tree unrooted-tree_Chip147.qza  \
--o-rooted-tree rooted-tree_Chip147.qza

  singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree_Chip147.qza \
  --i-table table_Chip147.qza \
  --p-sampling-depth 5000 \
  --m-metadata-file Chip147.tsv \
  --output-dir core-metrics-results_Chip147

  singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results_Chip147/faith_pd_vector.qza \
  --m-metadata-file Chip147.tsv \
  --o-visualization core-metrics-results_Chip147/faith-pd-group-significance_Chip147.qzv

singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results_Chip147/evenness_vector.qza \
  --m-metadata-file Chip147.tsv \
  --o-visualization core-metrics-results_Chip147/evenness-group-significance_Chip147.qzv

singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime feature-classifier classify-sklearn \
--i-classifier gg-13-8-99-nb-classifier.qza  \
--i-reads rep-seqs_Chip147.qza \
 --o-classification taxonomy_Chip147.qza

singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg  qiime metadata tabulate \
          --m-input-file taxonomy_Chip147.qza \
          --o-visualization taxonomy_Chip147.qzv

singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime taxa barplot \
          --i-table table_Chip147.qza \
          --i-taxonomy taxonomy_Chip147.qza \
          --m-metadata-file Chip147.tsv \
          --o-visualization taxa-bar-plots_Chip147.qzv
