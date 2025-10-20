#!/usr/bin/env bash
set -euo pipefail

# ./find_missing_biosamples.sh \
#   /Users/christiewoodside/Desktop/ngs_gca.tsv \
#   /Users/christiewoodside/Desktop/ngsQC_ARGOS_updated.tsv \
#   biosample

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: $0 biosample_gca.tsv biosampleMeta_ARGOS_updated.tsv [biosample_column_name]" >&2
  echo "Example: $0 biosample_gca.tsv biosampleMeta_ARGOS_updated.tsv biosample_id" >&2
  exit 1
fi

file_a="$1"   # biosample_gca.tsv
file_b="$2"   # biosampleMeta_ARGOS_updated.tsv
col="${3:-}"  # optional explicit column name (e.g., biosample_id)

awk -F'\t' -v OFS='\t' -v col="$col" \
    -v out_rows="missing_in_ARGOS_rows.tsv" \
    -v out_ids="missing_in_ARGOS_ids.txt" '
function trim(s){ gsub(/^[ \t\r\n]+|[ \t\r\n]+$/, "", s); return s }
function detect_idx(header,     i,n,low,arr,cands)
{
  split(header, arr, FS)
  # Candidate column names (case-insensitive). Add/remove as needed.
  cands["biosample"]=1
  cands["biosample_id"]=1
  cands["biosampleid"]=1
  cands["biosample accession"]=1
  cands["biosample_accession"]=1
  cands["sample_accession"]=1

  # If user specified a column, prefer that exact name.
  if (col != "") {
    for (i=1; i<=length(arr); i++) {
      n = trim(arr[i]); low = tolower(n)
      if (low == tolower(col)) return i
    }
    return 0
  }

  # Auto-detect from candidates
  for (i=1; i<=length(arr); i++) {
    n = trim(arr[i]); low = tolower(n)
    if (low in cands) return i
  }
  return 0
}

# First file read is file_b (ARGOS updated)
FNR==1 {
  if (NR==1) {
    b_idx = detect_idx($0)
    if (!b_idx) {
      print "ERROR: Could not find BioSample column in " FILENAME > "/dev/stderr"
      exit 1
    }
  } else {
    a_idx = detect_idx($0)
    if (!a_idx) {
      print "ERROR: Could not find BioSample column in " FILENAME > "/dev/stderr"
      exit 1
    }
    a_header = $0
    print a_header > out_rows
  }
  next
}

# Build set of BioSample IDs seen in ARGOS (file_b)
NR==FNR {
  key = trim($(b_idx))
  if (key != "") seen[key] = 1
  next
}

# For rows in biosample_gca (file_a): if BioSample not in seen, output row + collect ID
{
  key = trim($(a_idx))
  if (!(key in seen)) {
    print $0 >> out_rows
    missing[key] = 1
  }
}

END {
  # Write unique IDs
  for (k in missing) print k > out_ids
}
' "$file_b" "$file_a"

echo "Wrote:"
echo "  - missing_in_ARGOS_rows.tsv"
echo "  - missing_in_ARGOS_ids.txt"
