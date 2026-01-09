#!/bin/bash
# Updated December 5, 2025 by Christie Rose
# Script concatenates all of the genome assembly (fasta) files that were collected and unzipped in the preious steps. They are all concatenated into 
# one file that will be the database file. The compression step of this file needs to be run using run_compression.sh. 

source "$(dirname "$0")/config.sh"

logstepstart "Starting Step 5: Concatenating the Results"

cd "$OUTDIR" || exit 1

log "Creating missing_fna.txt (if empty_list.txt is present)..."
# Create missing_fna.txt
touch genomes/empty_list2.txt
cat genomes/empty_list.txt genomes/empty_list2.txt > missing_fna.txt

log "Concatenating .fna files into refseq_database_${VERSION}.fa..."
# Concatenate
cat genomes/*.fna > refseq_database_${VERSION}.fa    #added the test_ so I can compare outputs


logstepend "Step 5 completed successfully. Please run 'run_compression.sh' to complete final step."
