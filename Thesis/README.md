# Thesis: Bridging Genomics and Preparedness: Quality Control Metrics and Analysis for Emerging and Circulating Avian Influenza in 2024-2025
This thesis has been approved and published May 2025 for Christie Rose's conclusion of her M.S. Degree in Bioinformatics and Molecular Biochemistry. The thesis can be accessed at the [GWScholarSpace](https://scholarspace.library.gwu.edu/concern/gw_etds/6q182k87k?locale=en) website, but not available for open access until Spring 2026. If you would like to access the document please emial me at crwood230@gmail.com.

## Abstract
"The ongoing spread of highly pathogenic avian influenza (HPAI) H5N1 across a growing range of mammalian hosts emphasizes the urgent need for high-quality genomic data. This study presents a comprehensive quality control (QC) analysis of over 3,000 H5N1 genome assemblies collected from diverse host organisms across the United States, with a particular focus on cattle, poultry, and domestic animals. Leveraging a cloud-based computational pipeline from HIVE, sequence data were evaluated using standardized QC metrics, including GC content and Phred quality scores, to assess data integrity and identify potential sources of sequencing error or contamination. Regional differences in sequencing quality, such as elevated GC content and variable Phred scores in samples from California cattle, highlight inconsistencies in upstream sample handling and laboratory protocols. In addition to QC assessments, Sankey diagrams were used to visualize potential clonal diversity and segment-specific alignment across nearest neighbor pipeline-determined reference genomes. The results of this work provide a method to evaluate H5N1 sequences and address current gaps in influenza genomic data. These findings contribute to enhanced biosurveillance capacity and potential diagnostic uses and offer a scalable framework for evaluating future zoonotic threats."

## Thesis Code
The code within this folder was utilized to create the datatables for my thesis paper. The data tables used the same schema as those found on data.argosdb.org. It created three tables:
- H5N1_assemblyQC
- H5N1_ngsQC
- H5N1_biosampleMeta

## Data Collection
All of the data used in this dissertation was from the USDA BioProject PRJNA1102327. This includes all surveillance reported by the USDA and its affiliated labs. The organisms that are included, but not limited to, were Chickens, Duck, Dairy Cattle, Domestic Cats, Turkey, Mountain Lion, Canada Goose, Snow Goose, Racoon, Fox, Alpaca, Bobcat, Raven, and more. The data retreived, in addition to their QC metrics, can be found in the tables names APHIS_H5N1_biosampleMeta-march29, APHIS_H5N1_ngsQC-march29, and APHIS_H5N1_assemblyQC-march29 under the datasets folder in the repo.

## Observable Figures
The code that was used to create the figures can be found on Observable. 
The most current renderings of the image can be [found here](https://observablehq.com/d/26eb67bd1cc77b6a).
