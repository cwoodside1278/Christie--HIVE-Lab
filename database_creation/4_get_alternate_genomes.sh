#!/bin/bash
# Updated December 5, 2025 by Christie Rose
# Script process the alternate genomes that were found and unzips them to be added to the database. This will only run if alt_ids.txt exists
# and is not empty.

source "$(dirname "$0")/config.sh"

logstepstart "Starting Step 4: Getting Alternate Genomes"

cd "$OUTDIR/genomes" || exit 1

# Only proceed if alt_ids.txt exists and is not empty
if [ -s alt_ids.txt ]; then
  log "Processing alternate genomes from alt_ids.txt..."
  total_alts=$(wc -l < alt_ids.txt)
  current=0
  failed_downloads=0

  # Set up download logging
  exec 3>&1 4>&2
  trap 'exec 2>&4 1>&3' 0 1 2 3
  exec 1>"../../logs/4_download_progress.log" 2>&1

  while IFS= read -r i; do
    ((current++))
    i=$(echo "$i" | tr -d '[:space:]')

    log "Processing alternate genome $current/$total_alts: $i"

    # Check if zip file already exists
    if [ -f "${i}.zip" ]; then
      log "File already exists for ${i}, skipping download"
      continue
    fi

    url="https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/$i/download?include_annotation_type=CDS_FASTA&filename=$i.zip"

    if ! curl -f -OJX GET "$url" -H "Accept: application/zip" 2>&1; then
      log "Failed to download: $i"
      echo "$i" >> "../../logs/4_failed_downloads.txt"
      ((failed_downloads++))
      continue
    fi
  done < "alt_ids.txt"

  # Restore stdout/stderr before changing logs
  exec 1>&3 2>&4

  log "Download phase complete. Total: $total_alts, Failed: $failed_downloads"

  # Only proceed with extraction if we have zip files
  if ls *.zip 1> /dev/null 2>&1; then
    log "Extracting alternate genome files..."

    # Set up extraction logging
    exec 3>&1 4>&2
    trap 'exec 2>&4 1>&3' 0 1 2 3
    exec 1>"../../logs/4_extract_progress.log" 2>&1

    for zip_file in *.zip; do
      genome_id="${zip_file%.zip}"
      
      # Skip extraction if the FNA file already exists and isn't empty
      if [ -f "${genome_id}.fna" ] && [ -s "${genome_id}.fna" ]; then
        log "FNA file already exists for ${genome_id}, skipping extraction"
        continue
      fi

      log "Extracting: $zip_file"
      if ! unzip -p "$zip_file" "*.fna" > "${zip_file%.zip}.fna"; then
        log "Error extracting file: $zip_file"
        echo "${zip_file%.zip}" >> "../../logs/4_extraction_failed.txt"
        continue
      fi
    done

    # Restore stdout/stderr
    exec 1>&3 2>&4
  else
    # Restore stdout/stderr
    exec 1>&3 2>&4
    log "No zip files found to extract"
  fi

  log "Processing empty files..."
  # Process empty files
  find . -name "*.fna" -size 0 > empty_list2.txt
  sed -i 's/^\.\///;s/\.fna$//' empty_list2.txt
  find . -name "*.fna" -size 0 -delete

else 
  log "No alternate genomes to process (alt_ids.txt is empty or missing)"
fi

logstepend "Step 4 completed successfully"
