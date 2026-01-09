# Christie Woodside for HIVE-Lab
All public code and projects created by me for the HIVE Lab at George Washington University SMHS. 

## FDA-ARGOS Project
FDA-ARGOS database updates may help researchers rapidly validate diagnostic tests and use qualified genetic sequences to support future product development

As of September 2021, Embleema and George Washington University have been conducting bioinformatic research and system development, focusing on expanding the FDA-ARGOS database. This project expands datasets publicly available in FDA-ARGOS, improves quality control by developing quality matrix tools and scoring approaches that will allow the mining of public sequence databases, and identifies high-quality sequences for upload to the FDA-ARGOS database as regulatory-grade sequences. Building on expansions during the COVID-19 pandemic, this project aims to further improve the utility of the FDA-ARGOS database as a key tool for medical countermeasure development and validation.

The FDA-ARGOS Project GitHub Repository can be found [here](https://github.com/FDA-ARGOS/data.argosdb).

### Location
- Code files can be found in the [lib](https://github.com/FDA-ARGOS/data.argosdb/tree/main/lib) section of the FDA-ARGOS GitHub. 
- Most files created and edited by me can be found in the FDA-ARGOS [HIVE3](https://github.com/FDA-ARGOS/data.argosdb/tree/main/lib/HIVE3) folder.

## CensuScope Reference Databases

CensuScope is a census-based metagenomic analysis tool on the GW HIVE Platform used for rapid detection and classification of organisms in next-generation sequencing (NGS) metagenomic datasets. It performs census-based read generation and short-read alignment using BLAST and Bowtie, producing standardized taxonomic reports at the species level or higher resolution. CensuScope is routinely used across multiple HIVE Lab projects for contamination detection, microbiome analysis, and exploratory metagenomic profiling.

Methodology details:
https://pmc.ncbi.nlm.nih.gov/articles/PMC4218995/

CensuScope requires curated **reference database objects** for read mapping. This repository contains the database preparation and curation code used to support CensuScope analyses across multiple HIVE Lab research efforts.

This includes **slimNT**, a streamlined and curated alternative to the full NCBI NT database developed by the HIVE Lab (https://github.com/GW-HIVE/slimNT
). slimNT reduces computational overhead while preserving representative genomic diversity, enabling scalable and reproducible metagenomic analysis on the HIVE Platform. slimNT and related reference databases are shared across projects to ensure consistency while supporting project-specific performance requirements.

This repository also contains workflows to construct curated **reference genome databases** derived from NCBI Reference Sequence collections. Because NCBI reference genome databases are not directly downloadable as static files, custom scripts were developed to programmatically retrieve genomes flagged as reference genomes, assemble them into consolidated database objects, and prepare them for use within CensuScope.

Together, these reference databases provide a reusable, reproducible foundation for CensuScope-based metagenomic analyses across ongoing and future HIVE Lab projects.

## Master's Thesis
My thesis titled _Bridging Genomics and Preparedness: Quality Control Metrics and Analysis for Emerging and Circulating Avian Influenza in 2024-2025_ was completed and accepted in May 2025 for my Master's Program in Bioinformatics and Molecular Biochemistry.
All the data and code can be found in the folder titled [Thesis](https://github.com/cwoodside1278/Christie--HIVE-Lab/tree/main/Thesis).
