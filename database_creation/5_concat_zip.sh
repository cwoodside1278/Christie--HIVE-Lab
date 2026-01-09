#!/bin/bash

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
