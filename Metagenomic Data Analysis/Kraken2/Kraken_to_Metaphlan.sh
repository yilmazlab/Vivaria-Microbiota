#!/bin/sh
#!/bin/bash
#SBATCH --mail-user=bahtiyar.yilmaz@dbmr.unibe.ch  # Email for job notifications
#SBATCH --mail-type=end                            # Notify when the job ends
#SBATCH --job-name="Kraken2"                       # Job name
#SBATCH --nodes=1                                  # Number of nodes
#SBATCH --time=24:00:00                            # Maximum runtime
#SBATCH --mem-per-cpu=8G                           # Memory per CPU
#SBATCH --cpus-per-task=24                         # Number of CPUs per task
#SBATCH --partition=epyc2                          # Partition to run the job

##SBATCH --output=~/StomaMetagenomics/SlurmFiles   # Uncomment to specify output file
##SBATCH --error=~/StomaMetagenomics/SlurmFiles    # Uncomment to specify error file

export LC_CTYPE=en_US.UTF-8                        # Set locale for character encoding
export LC_ALL=en_US.UTF-8                          # Set locale for all settings
source activate kraken2                            # Activate the Kraken2 environment

# Build a database for plants
kraken2-build --download-taxonomy --db kraken2_plant_db  # Download taxonomy for the plant database
kraken2-build --download-library plant --db kraken2_plant_db --threads 24  # Download plant library
kraken2-build --build --db kraken2_plant_db --no-masking  # Build the plant database

# Build a custom database for other organisms
# Note: Plant DB is kept separate due to errors when combining
kraken2-build --build --db CustomDB                # Initialize the custom database
kraken2-build --build --db CustomDB/library        # Build the library structure
kraken2-build --download-library bacteria --db CustomDB  # Add bacteria library
kraken2-build --download-library plasmid --db CustomDB   # Add plasmid library
kraken2-build --download-library archaea --db CustomDB   # Add archaea library
kraken2-build --download-library viral --db CustomDB     # Add viral library
kraken2-build --download-library human --db CustomDB     # Add human library
kraken2-build --download-library fungi --db CustomDB     # Add fungi library
kraken2-build --download-library protozoa --db CustomDB  # Add protozoa library
kraken2-build --build --max-db-size 96000000000 --db CustomDB  # Build the custom database with a max size
kraken2-inspect --db /storage/research/dbmr_yilmaz_gut_evo/CustomDB/ | head 10  # Inspect the database

# Uncomment the following line to use a prebuilt Kraken2 library
# wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20210127.tar.gz

source activate kraken2                            # Reactivate Kraken2 environment

mkdir kraken_analysis                              # Create a directory for Kraken2 analysis results

# Run Kraken2 analysis for W samples
for i in {1..122};
      do kraken2 --use-names --db /home/ubelix/dbmr/terziev/CustomDB --threads 10 --minimum-base-quality 0 --confidence 0.2 --report-zero-counts --report kraken_analysis/W${i}.kreport --paired --gzip-compressed W${i}_R1.fastq.gz  W${i}_R2.fastq.gz  > kraken_analysis/W${i}.kraken
done

# Run Kraken2 analysis for SK samples
for i in {1..18};
      do kraken2 --use-names --db /home/ubelix/dbmr/terziev/CustomDB --threads 10 --minimum-base-quality 0  --confidence 0.2 --report-zero-counts  --report kraken_analysis/SK${i}.kreport --paired --gzip-compressed SK${i}_R1.fastq.gz  SK${i}_R2.fastq.gz  > kraken_analysis/SK${i}.kraken
done

# Run Kraken2 analysis for Bahti_H samples
for i in {1..50};
     do kraken2 --use-names --db /home/ubelix/dbmr/terziev/CustomDB --threads 10 --minimum-base-quality 0 --confidence 0.2 --report-zero-counts --report kraken_analysis/Bahti_H${i}.kreport --paired --gzip-compressed Bahti_H${i}_R1.fastq.gz  Bahti_H${i}_R2.fastq.gz  > kraken_analysis/Bahti_H${i}.kraken
  done

# Run Kraken2 analysis for Human samples
for i in {1..16};
      do kraken2 --use-names --db /storage/research/dbmr_yilmaz_gut_evo/CustomKrakenDB --threads 10 --minimum-base-quality 0 --confidence 0.2 --report-zero-counts --report kraken_analysis/report/Human${i}.kreport --paired --gzip-compressed /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/RAW_READS/Human/Human${i}_R1.fastq.gz /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/RAW_READS/Human/Human${i}_R2.fastq.gz  > kraken_analysis/Human${i}.kraken
done

cd /storage/research/dbmr_yilmaz_gut_evo/Wild_Mice_Analysis/Metagenomics/  # Change to analysis directory
# mkdir bracken_analysis                           # Uncomment to create a directory for Bracken analysis

# Run Bracken analysis for various sample groups
# Process Bahti_H samples
for i in {1..50}; do
   bracken -d /storage/research/dbmr_yilmaz_gut_evo/CustomKrakenDB \
           -i kraken_analysis/report/Bahti_H${i}.kreport \
           -l S -o bracken_analysis/Bahti_H${i}.bracken.txt
done

# Process W samples
for i in {1..122}; do
   bracken -d /storage/research/dbmr_yilmaz_gut_evo/CustomKrakenDB \
           -i kraken_analysis/report/W${i}.kreport \
           -l S -o bracken_analysis/W${i}.bracken.txt
done

# Process additional W samples
for i in {1308..1351}; do
   bracken -d /storage/research/dbmr_yilmaz_gut_evo/CustomKrakenDB \
           -i kraken_analysis/report/W${i}.kreport \
           -l S -o bracken_analysis/W${i}.bracken.txt
done

# Process SK samples
for i in {1..18}; do
   bracken -d /storage/research/dbmr_yilmaz_gut_evo/CustomKrakenDB \
           -i kraken_analysis/report/SK${i}.kreport \
           -l S -o bracken_analysis/SK${i}.bracken.txt
done

# Process smaller sample groups (RF, Poland, Greece, Crick, Sweden)
for i in {1..4}; do
   bracken -d /storage/research/dbmr_yilmaz_gut_evo/CustomKrakenDB \
           -i kraken_analysis/report/RF${i}.kreport \
           -l S -o bracken_analysis/RF${i}.bracken.txt
done

for i in {1..4}; do
   bracken -d /storage/research/dbmr_yilmaz_gut_evo/CustomKrakenDB \
           -i kraken_analysis/report/Poland${i}.kreport \
           -l S -o bracken_analysis/Poland${i}.bracken.txt
done

for i in {1..4}; do
   bracken -d /storage/research/dbmr_yilmaz_gut_evo/CustomKrakenDB \
           -i kraken_analysis/report/Greece${i}.kreport \
           -l S -o bracken_analysis/Greece${i}.bracken.txt
done

for i in {1..4}; do
   bracken -d /storage/research/dbmr_yilmaz_gut_evo/CustomKrakenDB \
           -i kraken_analysis/report/Crick${i}.kreport \
           -l S -o bracken_analysis/Crick${i}.bracken.txt
done

for i in {1..4}; do
   bracken -d /storage/research/dbmr_yilmaz_gut_evo/CustomKrakenDB \
           -i kraken_analysis/report/Sweden${i}.kreport \
           -l S -o bracken_analysis/Sweden${i}.bracken.txt
done

# Convert Kraken reports to Metaphlan format
# Process Bahti_H samples
for i in {1..50}; do
   python kreport2mpa.py -r /Users/bahti/Library/CloudStorage/Dropbox/Ongoing\ Analysis/Wild-Type\ Mice\ vs\ Wild\ Mice\ Analysis/2022/Metagenomic/report/Bahti_H${i}.kreport \
                         --percentages --display-header \
                         -o /Users/bahti/Library/CloudStorage/Dropbox/Ongoing\ Analysis/Wild-Type\ Mice\ vs\ Wild\ Mice\ Analysis/2022/Metagenomic/report/mpa/Bahti_H${i}.mpa.txt
done

# Process W samples
for i in {1..122}; do
   python kreport2mpa.py -r /Users/bahti/Library/CloudStorage/Dropbox/Ongoing\ Analysis/Wild-Type\ Mice\ vs\ Wild\ Mice\ Analysis/2022/Metagenomic/report/W${i}.kreport \
                         --percentages --display-header \
                         -o /Users/bahti/Library/CloudStorage/Dropbox/Ongoing\ Analysis/Wild-Type\ Mice\ vs\ Wild\ Mice\ Analysis/2022/Metagenomic/report/mpa/W${i}.mpa.txt
done

# Process additional W samples
for i in {1308..1351}; do
   python kreport2mpa.py -r /Users/bahti/Library/CloudStorage/Dropbox/Ongoing\ Analysis/Wild-Type\ Mice\ vs\ Wild\ Mice\ Analysis/2022/Metagenomic/report/W${i}.kreport \
                         --percentages --display-header \
                         -o /Users/bahti/Library/CloudStorage/Dropbox/Ongoing\ Analysis/Wild-Type\ Mice\ vs\ Wild\ Mice\ Analysis/2022/Metagenomic/report/mpa/W${i}.mpa.txt
done

# Process SK samples
for i in {1..18}; do
   python kreport2mpa.py -r /Users/bahti/Library/CloudStorage/Dropbox/Ongoing\ Analysis/Wild-Type\ Mice\ vs\ Wild\ Mice\ Analysis/2022/Metagenomic/report/SK${i}.kreport \
                         --percentages --display-header \
                         -o /Users/bahti/Library/CloudStorage/Dropbox/Ongoing\ Analysis/Wild-Type\ Mice\ vs\ Wild\ Mice\ Analysis/2022/Metagenomic/report/mpa/SK${i}.mpa.txt
done

# Process smaller sample groups (RF, Poland, Greece, Crick, Sweden)
for i in {1..4}; do
   python kreport2mpa.py -r /Users/bahti/Library/CloudStorage/Dropbox/Ongoing\ Analysis/Wild-Type\ Mice\ vs\ Wild\ Mice\ Analysis/2022/Metagenomic/report/RF${i}.kreport \
                         --percentages --display-header \
                         -o /Users/bahti/Library/CloudStorage/Dropbox/Ongoing\ Analysis/Wild-Type\ Mice\ vs\ Wild\ Mice\ Analysis/2022/Metagenomic/report/mpa/RF${i}.mpa.txt
done

for i in {1..4}; do
   python kreport2mpa.py -r /Users/bahti/Library/CloudStorage/Dropbox/Ongoing\ Analysis/Wild-Type\ Mice\ vs\ Wild\ Mice\ Analysis/2022/Metagenomic/report/Poland${i}.kreport \
                         --percentages --display-header \
                         -o /Users/bahti/Library/CloudStorage/Dropbox/Ongoing\ Analysis/Wild-Type\ Mice\ vs\ Wild\ Mice\ Analysis/2022/Metagenomic/report/mpa/Poland${i}.mpa.txt
done

for i in {1..4}; do
   python kreport2mpa.py -r /Users/bahti/Library/CloudStorage/Dropbox/Ongoing\ Analysis/Wild-Type\ Mice\ vs\ Wild\ Mice\ Analysis/2022/Metagenomic/report/Greece${i}.kreport \
                         --percentages --display-header \
                         -o /Users/bahti/Library/CloudStorage/Dropbox/Ongoing\ Analysis/Wild-Type\ Mice\ vs\ Wild\ Mice\ Analysis/2022/Metagenomic/report/mpa/Greece${i}.mpa.txt
done

for i in {1..4}; do
   python kreport2mpa.py -r /Users/bahti/Library/CloudStorage/Dropbox/Ongoing\ Analysis/Wild-Type\ Mice\ vs\ Wild\ Mice\ Analysis/2022/Metagenomic/report/Crick${i}.kreport \
                         --percentages --display-header \
                         -o /Users/bahti/Library/CloudStorage/Dropbox/Ongoing\ Analysis/Wild-Type\ Mice\ vs\ Wild\ Mice\ Analysis/2022/Metagenomic/report/mpa/Crick${i}.mpa.txt
done

for i in {1..4}; do
   python kreport2mpa.py -r /Users/bahti/Library/CloudStorage/Dropbox/Ongoing\ Analysis/Wild-Type\ Mice\ vs\ Wild\ Mice\ Analysis/2022/Metagenomic/report/Sweden${i}.kreport \
                         --percentages --display-header \
                         -o /Users/bahti/Library/CloudStorage/Dropbox/Ongoing\ Analysis/Wild-Type\ Mice\ vs\ Wild\ Mice\ Analysis/2022/Metagenomic/report/mpa/Sweden${i}.mpa.txt
done

# Process Human samples
for i in {1..16}; do
   python kreport2mpa.py -r /Users/bahti/Library/CloudStorage/Dropbox/Ongoing\ Analysis/Wild-Type\ Mice\ vs\ Wild\ Mice\ Analysis/2022/Metagenomic/report/Human${i}.kreport \
                         --percentages --display-header \
                         -o /Users/bahti/Library/CloudStorage/Dropbox/Ongoing\ Analysis/Wild-Type\ Mice\ vs\ Wild\ Mice\ Analysis/2022/Metagenomic/report/mpa/Human${i}.mpa.txt
done

# Combine all Metaphlan files into a single file
python combine_mpa.py -i Bahti_H1.mpa.txt Bahti_H2.mpa.txt Bahti_H3.mpa.txt Bahti_H4.mpa.txt Bahti_H5.mpa.txt Bahti_H6.mpa.txt Bahti_H7.mpa.txt Bahti_H8.mpa.txt Bahti_H9.mpa.txt Bahti_H10.mpa.txt Bahti_H11.mpa.txt Bahti_H12.mpa.txt Bahti_H13.mpa.txt Bahti_H14.mpa.txt Bahti_H15.mpa.txt Bahti_H16.mpa.txt Bahti_H17.mpa.txt Bahti_H18.mpa.txt Bahti_H19.mpa.txt Bahti_H20.mpa.txt Bahti_H21.mpa.txt Bahti_H22.mpa.txt Bahti_H23.mpa.txt Bahti_H24.mpa.txt Bahti_H25.mpa.txt Bahti_H26.mpa.txt Bahti_H27.mpa.txt Bahti_H28.mpa.txt Bahti_H29.mpa.txt Bahti_H30.mpa.txt Bahti_H31.mpa.txt Bahti_H32.mpa.txt Bahti_H33.mpa.txt Bahti_H34.mpa.txt Bahti_H35.mpa.txt Bahti_H36.mpa.txt Bahti_H37.mpa.txt Bahti_H38.mpa.txt Bahti_H39.mpa.txt Bahti_H40.mpa.txt Bahti_H41.mpa.txt Bahti_H42.mpa.txt Bahti_H43.mpa.txt Bahti_H44.mpa.txt Bahti_H45.mpa.txt Bahti_H46.mpa.txt Bahti_H47.mpa.txt Bahti_H48.mpa.txt Bahti_H49.mpa.txt Bahti_H50.mpa.txt Crick1.mpa.txt Crick2.mpa.txt Crick3.mpa.txt Crick4.mpa.txt Greece1.mpa.txt Greece2.mpa.txt Greece3.mpa.txt Greece4.mpa.txt Human1.mpa.txt Human2.mpa.txt Human3.mpa.txt Human4.mpa.txt Human5.mpa.txt Human6.mpa.txt Human7.mpa.txt Human8.mpa.txt Human9.mpa.txt Human10.mpa.txt Human11.mpa.txt Human12.mpa.txt Human13.mpa.txt Human14.mpa.txt Human15.mpa.txt Human16.mpa.txt Poland1.mpa.txt Poland2.mpa.txt Poland3.mpa.txt Poland4.mpa.txt RF1.mpa.txt RF2.mpa.txt RF3.mpa.txt RF4.mpa.txt SK1.mpa.txt SK2.mpa.txt SK3.mpa.txt SK4.mpa.txt SK5.mpa.txt SK6.mpa.txt SK7.mpa.txt SK8.mpa.txt SK9.mpa.txt SK10.mpa.txt SK11.mpa.txt SK12.mpa.txt SK13.mpa.txt SK14.mpa.txt SK15.mpa.txt SK16.mpa.txt SK17.mpa.txt SK18.mpa.txt Sweden1.mpa.txt Sweden2.mpa.txt Sweden3.mpa.txt Sweden4.mpa.txt W1.mpa.txt W2.mpa.txt W3.mpa.txt W4.mpa.txt W5.mpa.txt W6.mpa.txt W7.mpa.txt W8.mpa.txt W9.mpa.txt W10.mpa.txt W11.mpa.txt W12.mpa.txt W13.mpa.txt W14.mpa.txt W15.mpa.txt W16.mpa.txt W17.mpa.txt W18.mpa.txt W19.mpa.txt W20.mpa.txt W21.mpa.txt W22.mpa.txt W23.mpa.txt W24.mpa.txt W25.mpa.txt W26.mpa.txt W27.mpa.txt W28.mpa.txt W29.mpa.txt W30.mpa.txt W31.mpa.txt W32.mpa.txt W33.mpa.txt W34.mpa.txt W35.mpa.txt W36.mpa.txt W37.mpa.txt W38.mpa.txt W39.mpa.txt W40.mpa.txt W41.mpa.txt W42.mpa.txt W43.mpa.txt W44.mpa.txt W45.mpa.txt W46.mpa.txt W47.mpa.txt W48.mpa.txt W49.mpa.txt W50.mpa.txt W51.mpa.txt W52.mpa.txt W53.mpa.txt W54.mpa.txt W55.mpa.txt W56.mpa.txt W57.mpa.txt W58.mpa.txt W59.mpa.txt W60.mpa.txt W61.mpa.txt W62.mpa.txt W63.mpa.txt W64.mpa.txt W65.mpa.txt W66.mpa.txt W67.mpa.txt W68.mpa.txt W69.mpa.txt W70.mpa.txt W71.mpa.txt W72.mpa.txt W73.mpa.txt W74.mpa.txt W75.mpa.txt W76.mpa.txt W77.mpa.txt W78.mpa.txt W79.mpa.txt W80.mpa.txt W81.mpa.txt W82.mpa.txt W83.mpa.txt W84.mpa.txt W85.mpa.txt W86.mpa.txt W87.mpa.txt W88.mpa.txt W89.mpa.txt W90.mpa.txt W91.mpa.txt W92.mpa.txt W93.mpa.txt W94.mpa.txt W95.mpa.txt W96.mpa.txt W97.mpa.txt W98.mpa.txt W99.mpa.txt W100.mpa.txt W101.mpa.txt W102.mpa.txt W103.mpa.txt W104.mpa.txt W105.mpa.txt W106.mpa.txt W107.mpa.txt W108.mpa.txt W109.mpa.txt W110.mpa.txt W111.mpa.txt W112.mpa.txt W113.mpa.txt W114.mpa.txt W115.mpa.txt W116.mpa.txt W117.mpa.txt W118.mpa.txt W119.mpa.txt W120.mpa.txt W121.mpa.txt W122.mpa.txt W1308.mpa.txt W1309.mpa.txt W1310.mpa.txt W1311.mpa.txt W1312.mpa.txt W1313.mpa.txt W1314.mpa.txt W1315.mpa.txt W1316.mpa.txt W1317.mpa.txt W1318.mpa.txt W1319.mpa.txt W1320.mpa.txt W1321.mpa.txt W1322.mpa.txt W1323.mpa.txt W1324.mpa.txt W1325.mpa.txt W1326.mpa.txt W1327.mpa.txt W1328.mpa.txt W1329.mpa.txt W1330.mpa.txt W1331.mpa.txt W1332.mpa.txt W1333.mpa.txt W1334.mpa.txt W1335.mpa.txt W1336.mpa.txt W1337.mpa.txt W1338.mpa.txt W1339.mpa.txt W1340.mpa.txt W1341.mpa.txt W1342.mpa.txt W1343.mpa.txt W1344.mpa.txt W1345.mpa.txt W1346.mpa.txt W1347.mpa.txt W1348.mpa.txt W1349.mpa.txt W1350.mpa.txt W1351.mpa.txt -o wildmicestudy.combined.mpa.txt
