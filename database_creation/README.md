# CensuScope Reference Databases
## Overview

CensuScope is a census-based metagenomic analysis tool available on the GW HIVE Platform. It is designed and optimized for the rapid detection and classification of components within next-generation sequencing (NGS) metagenomic datasets. CensuScope generates census-based read files and performs short-read alignment using BLAST and Bowtie, producing standardized reports of detected organisms at the species level or higher taxonomic resolution. Common applications include contamination detection, microbiome profiling, and exploratory metagenomic analysis.

More details on the CensuScope methodology can be found here:
[https://pmc.ncbi.nlm.nih.gov/articles/PMC4218995/](https://pmc.ncbi.nlm.nih.gov/articles/PMC4218995/)

CensuScope requires a database object that serves as the reference for read mapping. This repository supports database preparation and curation workflows used by CensuScope on the HIVE Platform.

# Reference Genome Database Curation (This Repository)

The code in this repository supports the manual construction of a curated reference genome database compatible with CensuScope.

NCBI BLAST reference sequence databases such as refseq_genomes and refseq_reference_genomes are not directly downloadable as standalone objects. Instead, they are typically accessed dynamically through NCBI tools using specific flags. To support reproducible and offline workflows on the HIVE Platform, this repository provides scripts that:

Access the NCBI API

1. Identify genomes annotated with the “reference genome” flag
2. Download the corresponding genome assemblies
3. Manually concatenate them into a single, usable database file for CensuScope

The resulting database mirrors the intent of NCBI’s reference genome collections while enabling local control, transparency, and reproducibility within the HIVE environment.

# slimNT Database

Our lab has developed slimNT, a curated and streamlined alternative to the full NCBI NT database, which is available at:
[https://github.com/GW-HIVE/slimNT](https://github.com/GW-HIVE/slimNT)

Metagenomic analysis against the full NT database (currently containing hundreds of thousands of genomes) can be both time-consuming and computationally expensive, particularly during the indexing step. slimNT addresses this challenge by strategically subsetting NT to create a compact, user-defined database tailored to specific analytical needs.

The slimNT database is derived from Representative Proteomes (RPs) and Reference Proteome Groups (RPGs) provided by the Protein Information Resource (PIR). Users may select reference proteomes and viral reference proteomes based on configurable cutoff thresholds:

* Higher cutoff values generate larger, more comprehensive databases that improve detection sensitivity at the cost of increased computational requirements.

* Lower cutoff values produce smaller databases that enable faster querying but may preferentially detect more distantly related sequences.

This flexibility allows users to balance speed and accuracy based on project goals and available resources.

To ensure portability and reproducibility, slimNT is built using a Nextflow-based pipeline that automates database aggregation and construction. Detailed installation and usage instructions are provided in the slimNT repository.


# Intended Use

The databases generated through these workflows are intended for use with CensuScope to support:

* Rapid assessment of metagenomic sample composition

* Detection of sample contamination

* Microbiome and environmental sequencing analysis

* Exploratory pathogen and taxonomic profiling

Please see the associated repositories for implementation details and pipeline code.
