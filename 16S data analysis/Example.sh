#!/bin/sh
#!/bin/bash
#SBATCH --job-name="Example"  # Job name
#SBATCH --nodes=1  # Number of nodes
#SBATCH --time=48:00:00  # Maximum runtime (48 hours)
#SBATCH --mem-per-cpu=64G  # Memory per CPU
#SBATCH --partition=all  # Partition to submit the job to

##SBATCH --output=/path/to/outfile  # Uncomment to specify output file
##SBATCH --error=/path/to/errfile  # Uncomment to specify error file

# For array jobs
# Array job containing 100 tasks, run max 10 tasks at the same time
##SBATCH --array=1-100%10  # Uncomment to enable array jobs

#### Your shell commands below this line ####

## QIIME pipeline shell script ##

# singularity help /software/singularity/containers/qiime2-2018.8-1.debian9.simg
# Navigate to the working directory
cd ~/Wild_Mice_Analysis/QIIME2_WildMice/Example/

# Step 1: Import data into QIIME2
singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime tools import \
--type MultiplexedSingleEndBarcodeInSequence \
--input-path Example.fastq.gz \
--output-path multiplexed-seqs_Example.qza

# Step 2: Demultiplex sequences using barcodes
singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime cutadapt demux-single \
--i-seqs multiplexed-seqs_Example.qza \
--m-barcodes-file Example.tsv \
--m-barcodes-column BarcodesSequence \
--p-error-rate 0 \
--o-per-sample-sequences demultiplexed-seqs_Example.qza \
--o-untrimmed-sequences untrimmed_Example.qza \
--verbose

# Step 3: Trim primers from sequences
singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime cutadapt trim-single \
--i-demultiplexed-sequences demultiplexed-seqs_Example.qza \
--p-front ATTAGATACCCYGGTAGTCC \
--p-error-rate 0.1 \
--o-trimmed-sequences trimmed-seqs_Example.qza \
--verbose

# Step 4: Summarize demultiplexed data
singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg  qiime demux summarize \
--i-data trimmed-seqs_Example.qza \
--o-visualization trimmed-seqs_Example.qzv

# Step 5: Perform DADA2 denoising
singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime dada2 denoise-single \
--i-demultiplexed-seqs demultiplexed-seqs_Example.qza  \
--p-trunc-len 120  \
--o-table table_Example.qza \
--o-representative-sequences rep-seqs_Example.qza  \
--o-denoising-stats denoising-stats_Example.qza

# Step 6: Summarize feature table
singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime feature-table summarize \
--i-table table_Example.qza \
--o-visualization table_Example.qzv \
--m-sample-metadata-file Example.tsv

# Step 7: Tabulate representative sequences
singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime feature-table tabulate-seqs \
--i-data rep-seqs_Example.qza \
--o-visualization rep-seqs_Example.qzv

# Step 8: Generate phylogenetic tree
singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences rep-seqs_Example.qza  \
--o-alignment aligned-rep-seqs_Example.qza \
--o-masked-alignment masked-aligned-rep-seqs_Example.qza  \
--o-tree unrooted-tree_Example.qza  \
--o-rooted-tree rooted-tree_Example.qza

# Step 9: Perform core diversity metrics analysis
singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime diversity core-metrics-phylogenetic \
--i-phylogeny rooted-tree_Example.qza \
--i-table table_Example.qza \
--p-sampling-depth 5000 \
--m-metadata-file Example.tsv \
--output-dir core-metrics-results_Example

# Step 10: Analyze alpha diversity (Faith's Phylogenetic Diversity)
singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime diversity alpha-group-significance \
--i-alpha-diversity core-metrics-results_Example/faith_pd_vector.qza \
--m-metadata-file Example.tsv \
--o-visualization core-metrics-results_Example/faith-pd-group-significance_Example.qzv

# Step 11: Analyze alpha diversity (Evenness)
singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime diversity alpha-group-significance \
--i-alpha-diversity core-metrics-results_Example/evenness_vector.qza \
--m-metadata-file Example.tsv \
--o-visualization core-metrics-results_Example/evenness-group-significance_Example.qzv

# Step 12: Assign taxonomy using a pre-trained classifier
singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime feature-classifier classify-sklearn \
--i-classifier silva-138-99-nb-classifier.qza  \
--i-reads rep-seqs_Example.qza \
--o-classification taxonomy_Example.qza

# Step 13: Tabulate taxonomy
singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg  qiime metadata tabulate \
--m-input-file taxonomy_Example.qza \
--o-visualization taxonomy_Example.qzv

# Step 14: Generate taxa bar plots
singularity exec /software/singularity/containers/qiime2-2018.8-1.debian9.simg qiime taxa barplot \
--i-table table_Example.qza \
--i-taxonomy taxonomy_Example.qza \
--m-metadata-file Example.tsv \
--o-visualization taxa-bar-plots_Example.qzv
