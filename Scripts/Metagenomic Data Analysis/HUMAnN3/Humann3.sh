#!/bin/sh
#!/bin/bash
#SBATCH --job-name="DownloadDatabase"
#SBATCH --nodes=1
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=24
#SBATCH --partition=epyc2

##SBATCH --output=/path/to/outfile
##SBATCH --error=/path/to/errfile

#For array jobs
#Array job containing 10 tasks, run max 10tasks at the same time
#SBATCH --array=1-37%20

# humann_databases --download chocophlan full  humann3_chocophlan_database
# humann_databases --download uniref uniref90_diamond  humann3_UniRef90
# conda install -c bioconda diamond=0.9.36
# metaphlan --install --bowtie2db bowtie2db
source activate mpa
cd /storage/workspaces/dbmr_mg/nbmisi/Bahti_stoma

cat /storage/homefs/terziev/StomaMetagenomics/CLEAN_READS/SF1_1.fastq /storage/homefs/terziev/StomaMetagenomics/CLEAN_READS/SF1_2.fastq > /storage/homefs/terziev/StomaMetagenomics/CLEAN_READS/SF1_12.fastq
humann --input /storage/homefs/terziev/StomaMetagenomics/CLEAN_READS/SF1_12.fastq  --output /storage/research/dbmr_yilmaz_gut_evo/Stoma/RAW_READS/Humann3/SF1 --nucleotide-database /storage/homefs/terziev/humann3_chocophlan_database/chocophlan --bowtie2  /storage/homefs/terziev/bowtie2db/bowtie2-build \

for i in {1..43};
 do cat /storage/homefs/terziev/StomaMetagenomics/CLEAN_READS/SF${i}_1.fastq /storage/homefs/terziev/StomaMetagenomics/CLEAN_READS/SF${i}_2.fastq > /storage/workspaces/dbmr_mg/nbmisi/Bahti_stoma/SF${i}_12.fastq
done

for i in {1..37}
    do humann -i SF${i}_12.fastq -o Output/SF${i} --threads 20 --nucleotide-database /storage/homefs/terziev/humann3_chocophlan_database/chocophlan --bowtie2 storage/homefs/terziev/bowtie2db/bowtie2-build
done


cd /storage/research/dbmr_yilmaz_gut_evo/Stoma/Humann3\ Output
# Joining the tables
humann_join_tables --input  genefamilies --output genefamilies_merged.tsv
humann_join_tables --input  pathcoverage --output pathcoverage.tsv
humann_join_tables --input  pathabundance --output pathabundance.tsv

# Normalizing RPKs to relative abundance
humann_renorm_table --input genefamilies_merged.tsv --output genefamilies_merged-cpm.tsv --units relab --update-snames
humann_renorm_table --input pathcoverage.tsv --output pathcoverage_merged-cpm.tsv --units relab --update-snames
humann_renorm_table --input pathabundance.tsv --output pathabundance_merged-cpm.tsv --units relab --update-snames

humann_regroup_table --input genefamilies_merged-cpm.tsv  --output genefamilies_merged-rxn-cpm.tsv --groups uniref90_rxn
humann_regroup_table --input pathcoverage_merged-cpm.tsv  --output pathcoverage_merged-rxn-cpm.tsv --groups uniref90_rxn
humann_regroup_table --input pathabundance_merged-cpm.tsv  --output pathabundance_merged-rxn-cpm.tsv --groups uniref90_rxn

# Attaching names to features
humann_rename_table --input genefamilies_merged-rxn-cpm.tsv --output genefamilies_merged-rxn-cpm-named.tsv --names metacyc-rxn
humann_rename_table --input pathcoverage_merged-rxn-cpm.tsv --output pathcoverage_merged-rxn-cpm-named.tsv --names metacyc-rxn
humann_rename_table --input pathabundance_merged-rxn-cpm.tsv --output pathabundance_merged-rxn-cpm-named.tsv --names metacyc-rxn

## 
