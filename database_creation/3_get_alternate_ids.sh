#!/bin/bash

source "$(dirname "$0")/config.sh"

LOGDIR=/data/christie/refdata/logs
mkdir -p "$LOGDIR"



logstepstart "Starting Step 3: Getting Alternate IDs"

cd "$OUTDIR/genomes" || exit 1

log "Downloading asssembly summary..."
# Download and process assembly summary
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt
awk -F '\t' '{print $1, $18}' assembly_summary_refseq.txt > ref_gen_map.txt
rm assembly_summary_refseq.txt     

log "Processing empty list..."
# Process empty list
awk 'NR==FNR {empty[$0]; next} {for (i=1; i<=NF; i++) if ($i in empty) {for (j=1; j<=NF; j++) if (j != i) printf "%s", $j; print ""}}' 'empty_list.txt' ref_gen_map.txt > alt_ids.txt

logstepend "Step 3 completed successfully"
