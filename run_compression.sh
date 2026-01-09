#!/usr/bin/env bash
set -euo pipefail
mkdir -p logs



# Parse command line args
VERSION=""
BACKUP_DIR=""


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


#BACKUP_DIR=""
while [[ "$#" -gt 0 ]]; do
  case $1 in
    --backup-dir) BACKUP_DIR="$2"; shift ;;
    *) echo "Unknown arg: $1"; exit 1 ;;
  esac
  shift
done

export BACKUP_DIR
export VERSION

#new
#new
TS="$(date +%Y%m%d_%H%M%S)"
exec > >(tee -a "logs/compress_${VERSION:-unknown}_${TS}.out") 2> >(tee -a "logs/compress_${VERSION:-unknown}_${TS}.err" >&2)



set -e
cd /data/christie/refdata || exit 1

mkdir -p logs

chmod +x run_compression.sh      #updated
chmod +x ./pipeline/*.sh

echo "Compression of file started at $(date)"
echo "Using backup directory: ${BACKUP_DIR:-None}"
echo "Compressing the file for version: $VERSION"
#updated so it only runs the compression script

./pipeline/6_compress_files.sh

echo "Compression done, pipeline completed at $(date)"

