#!/bin/sh
#!/bin/bash
#SBATCH --job-name="Hu6"
#SBATCH --mail-user=bahtiyar.yilmaz@dbmr.unibe.ch
#SBATCH --mail-type=end
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=16
#SBATCH --partition=epyc2

##SBATCH --output=/path/to/outfile
##SBATCH --error=/path/to/errfile

# For array jobs
# Array job containing 100 tasks, run max 10 tasks at the same time
##SBATCH --array=5-6%2

# humann_databases --download chocophlan full  humann3_chocophlan_database
# humann_databases --download uniref uniref90_diamond  humann3_UniRef90
# conda install -c bioconda diamond=0.9.36
# metaphlan --install --bowtie2db bowtie2db

export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8
source activate mpa
cd /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/RAW_READS

humann --input /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/Concatenate/Human6.fastq.gz --output /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/Humann3/Human6 --nucleotide-database /storage/homefs/terziev/humann3_chocophlan_database/chocophlan --bowtie2 /storage/homefs/terziev/bowtie2db/bowtie2-build --threads 20
