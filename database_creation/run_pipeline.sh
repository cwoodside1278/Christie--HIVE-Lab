#!/usr/bin/env bash
# Created December 5, 2025 by Christie Rose

# Run this file via command line to run the database creation pipeline. There is version control, so please add either text or number after the 
# --version flag. A hard coded filepath location needs to be updated when you use these scripts, see below. 


set -euo pipefail
mkdir -p logs


# Parse command line args
VERSION=""

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --version)
      VERSION="$2"
      shift 2
      ;;
    *)
      echo "Unknown argument: $1" >&2
      echo "Usage: sbatch run_pipeline.sh --version <version>" >&2
      exit 1
      ;;
  esac
done

if [[ -z "$VERSION" ]]; then
  echo "ERROR: --version is required." >&2
  echo "Usage: sbatch run_pipeline.sh --version <version>" >&2
  exit 1
fi


BACKUP_DIR=""

while [[ "$#" -gt 0 ]]; do
  case $1 in
    --backup-dir) BACKUP_DIR="$2"; shift ;;
    #--version)    VERSION="$2"; shift ;;
    *) echo "Unknown arg: $1"; exit 1 ;;
  esac
  shift
done


#if [[ -z "$VERSION" ]]; then
  #echo "ERROR: you must supply --version <string> (e.g. 1.0)"; exit 1
#fi

export BACKUP_DIR
export VERSION

#new
TS="$(date +%Y%m%d_%H%M%S)"
exec > >(tee -a "logs/job_${VERSION:-unknown}_${TS}.out") 2> >(tee -a "logs/job_${VERSION:-unknown}_${TS}.err" >&2)


set -e

cd /data/USER/refdata || exit 1  # MANUALLY INPUT YOUR DIRECTORY PATH


mkdir -p logs

chmod +x run_pipeline.sh
chmod +x ./pipeline/*.sh
#second option:
#chmod +x run_pipeline.sh || true
#chmod +x ./pipeline/*.sh || true

echo "Reference Genome Database Pipeline started at $(date)"
echo "Using backup directory: ${BACKUP_DIR:-None}"
echo "Creating file version: $VERSION"

./pipeline/2_get_genomes.sh
./pipeline/3_get_alternate_ids.sh
./pipeline/4_get_alternate_genomes.sh
./pipeline/5_concat_zip.sh

echo "Reference Genome Pipeline completed at $(date)"
