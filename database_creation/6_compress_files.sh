#!/bin/bash
# Updated and added June 6, 2025
# Created so that we can run the compression step separately from the pipeline

source "$(dirname "$0")/config.sh"

logstepstart "Starting Step 6: Zipping the Results"

cd "$OUTDIR" || exit 1

log "Finding refseq_database_${VERSION}.fa..."

log "Compressing final database for refseq_database_${VERSION}.fa..."
gzip refseq_database_${VERSION}.fa   #now has versioning capabilities

logstepend "Step 6- Compression of RefSeq Database (HIVE version) was completed successfully"


