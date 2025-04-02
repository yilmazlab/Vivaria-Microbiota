#!/bin/sh
#!/bin/bash
#SBATCH --mail-user=bahtiyar.yilmaz@dbmr.unibe.ch  # Email address for job notifications
#SBATCH --mail-type=end                            # Notify when the job ends
#SBATCH --job-name="inStrainW"                     # Job name
#SBATCH --nodes=1                                  # Number of nodes
#SBATCH --time=24:00:00                            # Maximum runtime
#SBATCH --mem-per-cpu=11G                          # Memory per CPU
#SBATCH --cpus-per-task=9                          # Number of CPUs per task
#SBATCH --partition=bdw                            # Partition to submit the job

##SBATCH --output=/path/to/outfile                 # Uncomment to specify output file
##SBATCH --error=/path/to/errfile                  # Uncomment to specify error file

# For array jobs
# Array job containing 10 tasks, run max 20 tasks at the same time
#SBATCH --array=1-4%4

#### Your shell commands below this line ####
#### Tutorial: https://github.com/bxlab/metaWRAP/blob/master/Usage_tutorial.md

# Set locale settings
export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

# Load required modules
module load vital-it
module add UHTS/Aligner/bowtie2/2.3.4.1
module add SequenceAnalysis/GenePrediction/prodigal/2.6.3
module add UHTS/Analysis/Mash/2.0
module add UHTS/Analysis/mummer/4.0.0beta1
module add UHTS/Analysis/fastANI/1.1
module add UHTS/Analysis/samtools/1.10

# Change to the working directory
cd /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics

# Process samples for different regions
# Each block runs bowtie2 for alignment and inStrain for profiling

# Process Crick samples
bowtie2 -p 10 -x /storage/homefs/terziev/UHGGv1/Reps_bt2/UHGG_reps.fasta.bt2 \
  -1 /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/RAW_READS/Crick${SLURM_ARRAY_TASK_ID}_R1.fastq.gz \
  -2 /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/RAW_READS/Crick${SLURM_ARRAY_TASK_ID}_R2.fastq.gz \
  > /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/samfiles/Crick${SLURM_ARRAY_TASK_ID}.sam
inStrain profile /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/samfiles/Crick${SLURM_ARRAY_TASK_ID}.sam \
  /storage/homefs/terziev/UHGGv1/UHGG_reps.fasta \
  -o /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/Instrain/Crick${SLURM_ARRAY_TASK_ID}.IS \
  -p 24 -g /storage/homefs/terziev/UHGGv1/UHGG_reps.genes.fna \
  -s /storage/homefs/terziev/UHGGv1/UHGG_reps.stb --database_mode

# Process Greece samples
bowtie2 -p 10 -x /storage/homefs/terziev/UHGGv1/Reps_bt2/UHGG_reps.fasta.bt2 \
  -1 /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/RAW_READS/Greece${SLURM_ARRAY_TASK_ID}_R1.fastq.gz \
  -2 /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/RAW_READS/Greece${SLURM_ARRAY_TASK_ID}_R2.fastq.gz \
  > /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/samfiles/Greece${SLURM_ARRAY_TASK_ID}.sam
inStrain profile /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/samfiles/Greece${SLURM_ARRAY_TASK_ID}.sam \
  /storage/homefs/terziev/UHGGv1/UHGG_reps.fasta \
  -o /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/Instrain/Greece${SLURM_ARRAY_TASK_ID}.IS \
  -p 24 -g /storage/homefs/terziev/UHGGv1/UHGG_reps.genes.fna \
  -s /storage/homefs/terziev/UHGGv1/UHGG_reps.stb --database_mode

# Process Poland samples
bowtie2 -p 10 -x /storage/homefs/terziev/UHGGv1/Reps_bt2/UHGG_reps.fasta.bt2 \
  -1 /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/RAW_READS/Poland${SLURM_ARRAY_TASK_ID}_R1.fastq.gz \
  -2 /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/RAW_READS/Poland${SLURM_ARRAY_TASK_ID}_R2.fastq.gz \
  > /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/samfiles/Poland${SLURM_ARRAY_TASK_ID}.sam
inStrain profile /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/samfiles/Poland${SLURM_ARRAY_TASK_ID}.sam \
  /storage/homefs/terziev/UHGGv1/UHGG_reps.fasta \
  -o /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/Instrain/Poland${SLURM_ARRAY_TASK_ID}.IS \
  -p 24 -g /storage/homefs/terziev/UHGGv1/UHGG_reps.genes.fna \
  -s /storage/homefs/terziev/UHGGv1/UHGG_reps.stb --database_mode

# Process RF samples
bowtie2 -p 10 -x /storage/homefs/terziev/UHGGv1/Reps_bt2/UHGG_reps.fasta.bt2 \
  -1 /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/RAW_READS/RF${SLURM_ARRAY_TASK_ID}_R1.fastq.gz \
  -2 /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/RAW_READS/RF${SLURM_ARRAY_TASK_ID}_R2.fastq.gz \
  > /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/samfiles/RF${SLURM_ARRAY_TASK_ID}.sam
inStrain profile /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/samfiles/RF${SLURM_ARRAY_TASK_ID}.sam \
  /storage/homefs/terziev/UHGGv1/UHGG_reps.fasta \
  -o /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/Instrain/RF${SLURM_ARRAY_TASK_ID}.IS \
  -p 24 -g /storage/homefs/terziev/UHGGv1/UHGG_reps.genes.fna \
  -s /storage/homefs/terziev/UHGGv1/UHGG_reps.stb --database_mode

# Process Sweden samples
bowtie2 -p 10 -x /storage/homefs/terziev/UHGGv1/Reps_bt2/UHGG_reps.fasta.bt2 \
  -1 /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/RAW_READS/Sweden${SLURM_ARRAY_TASK_ID}_R1.fastq.gz \
  -2 /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/RAW_READS/Sweden${SLURM_ARRAY_TASK_ID}_R2.fastq.gz \
  > /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/samfiles/Sweden${SLURM_ARRAY_TASK_ID}.sam
inStrain profile /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/samfiles/Sweden${SLURM_ARRAY_TASK_ID}.sam \
  /storage/homefs/terziev/UHGGv1/UHGG_reps.fasta \
  -o /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/Instrain/Sweden${SLURM_ARRAY_TASK_ID}.IS \
  -p 24 -g /storage/homefs/terziev/UHGGv1/UHGG_reps.genes.fna \
  -s /storage/homefs/terziev/UHGGv1/UHGG_reps.stb --database_mode

# Process Wild mice samples
bowtie2 -p 10 -x /storage/homefs/terziev/UHGGv1/Reps_bt2/UHGG_reps.fasta.bt2 \
  -1 /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/RAW_READS/W${SLURM_ARRAY_TASK_ID}_R1.fastq.gz \
  -2 /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/RAW_READS/W${SLURM_ARRAY_TASK_ID}_R2.fastq.gz \
  > /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/samfiles/W${SLURM_ARRAY_TASK_ID}.sam
inStrain profile /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/samfiles/W${SLURM_ARRAY_TASK_ID}.sam \
  /storage/homefs/terziev/UHGGv1/UHGG_reps.fasta \
  -o /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/Instrain/W${SLURM_ARRAY_TASK_ID}.IS \
  -p 24 -g /storage/homefs/terziev/UHGGv1/UHGG_reps.genes.fna \
  -s /storage/homefs/terziev/UHGGv1/UHGG_reps.stb --database_mode

# Process SK samples
bowtie2 -p 10 -x /storage/homefs/terziev/UHGGv1/Reps_bt2/UHGG_reps.fasta.bt2 \
  -1 /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/RAW_READS/SK${SLURM_ARRAY_TASK_ID}_R1.fastq.gz \
  -2 /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/RAW_READS/SK${SLURM_ARRAY_TASK_ID}_R2.fastq.gz \
  > /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/samfiles/SK${SLURM_ARRAY_TASK_ID}.sam
inStrain profile /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/samfiles/SK${SLURM_ARRAY_TASK_ID}.sam \
  /storage/homefs/terziev/UHGGv1/UHGG_reps.fasta \
  -o /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/Instrain/SK${SLURM_ARRAY_TASK_ID}.IS \
  -p 24 -g /storage/homefs/terziev/UHGGv1/UHGG_reps.genes.fna \
  -s /storage/homefs/terziev/UHGGv1/UHGG_reps.stb --database_mode

# Process Human samples
bowtie2 -p 10 -x /storage/homefs/terziev/UHGGv1/Reps_bt2/UHGG_reps.fasta.bt2 \
  -1 /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/RAW_READS/Human/Human${SLURM_ARRAY_TASK_ID}_R1.fastq.gz \
  -2 /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/RAW_READS/Human/Human${SLURM_ARRAY_TASK_ID}_R2.fastq.gz \
  > /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/samfiles/Human${SLURM_ARRAY_TASK_ID}.sam
inStrain profile /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/samfiles/Human${SLURM_ARRAY_TASK_ID}.sam \
  /storage/homefs/terziev/UHGGv1/UHGG_reps.fasta \
  -o /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/Instrain/Human${SLURM_ARRAY_TASK_ID}.IS \
  -p 24 -g /storage/homefs/terziev/UHGGv1/UHGG_reps.genes.fna \
  -s /storage/homefs/terziev/UHGGv1/UHGG_reps.stb --database_mode
