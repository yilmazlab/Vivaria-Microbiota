# Vivaria Microbiota

This repository contains the scripts and data used in our manuscript published in ... (202..). The study focuses on analyzing microbiota data using state-of-the-art bioinformatics pipelines. 

## Manuscript
Link to manuscript: [Insert Link Here]

## Overview of the Pipelines
The following pipelines and tools were utilized in this study:

1. **QIIME2** ([Documentation](https://qiime2.org)) - For microbial community analysis.
2. **Phyloseq** ([Documentation](https://joey711.github.io/phyloseq/)) - For microbiome data analysis and visualization in R.
3. **MaAsLin2** ([Documentation](https://huttenhower.sph.harvard.edu/maaslin/)) - For multivariable association analysis.
4. **Kraken2** ([Documentation](https://ccb.jhu.edu/software/kraken2/)) - For taxonomic classification of metagenomic sequences.
5. **Braken** ([Documentation](https://ccb.jhu.edu/software/bracken/)) - For abundance estimation from Kraken2 outputs.
6. **InStrain** ([Documentation](https://instrain.readthedocs.io/en/latest/)) - For strain-level analysis of metagenomic data.
7. **HAllA** ([Documentation](https://huttenhower.sph.harvard.edu/halla/)) - For discovering multi-omic associations.
8. **metaWRAP** ([Documentation](https://github.com/bxlab/metaWRAP)) - For metagenomic data analysis and binning.
9. **MIMOSA2** ([Documentation](https://borenstein-lab.github.io/MIMOSA2shiny/)) - For metabolic modeling and analysis.

## Data and Code Availability
- **Metagenomic Sequencing Data**: All data generated in this study has been deposited in ... with BioProject Accession (Accession: [Insert Accession Here]).
- **16S rRNA Amplicon Sequencing Dataset**: Processed files and sample details can be downloaded from [Figshare](https://figshare.com/s/6b9e2c30c49bef0790ef) (DOI: [10.6084/m9.figshare.24523324](https://doi.org/10.6084/m9.figshare.24523324)).

- **An R Shiny app** was created to facilitate quick and efficient exploration of the meta-omics dataset from this study. The app is accessible at the following URLs:
  - OMICS data: [https://yilmazlab.shinyapps.io/vivaria_app/](https://yilmazlab.shinyapps.io/vivaria_app/)
  - Vivaria variants: [https://yilmazlab.shinyapps.io/vivaria_app2/](https://yilmazlab.shinyapps.io/vivaria_app2/)

> **Note**: No new code was developed in this study. We followed the developers' recommendations for each tool with minor optimizations.

## How to apply similar the analysis to your studies
1. Clone this repository:
   ```bash
   git clone https://github.com/your-repo/Vivaria-Microbiota.git
   cd Vivaria-Microbiota
   ```
2. Follow the documentation for each pipeline listed above to set up the required tools and dependencies.
3. Use the provided data and scripts to replicate the representative analysis.

## Citation
If you use this repository or data in your research, please cite our manuscript:
[Insert Citation Here]

## Contact
For questions or issues, please contact [Prof.Dr.Bahtiyar Yilmaz : bahtiyar.yilmaz@unibe.ch].

## License and third-party code
This repository is released under the MIT License. The MIT License applies to materials authored in this repository, including wrapper scripts, configuration files, and example workflows.
Please note that this repository makes use of third-party tools and published pipelines (e.g., QIIME2, MaasLin2, Kraken2, Humann3 etc.). These tools retain their original licenses and are not re-licensed under MIT by inclusion here. Users must comply with the respective licenses of these external software packages.



