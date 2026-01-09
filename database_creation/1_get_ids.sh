#!/bin/bash
# Updated December 5, 2025 by Christie Rose

### No longer applicable, this script is slimNT specific. It was left in here to allow for minimal changes in the pipeline to accomodate the time crunch
# of updating the code to accomodate reference sequence database creation. Ignore this script. 

source "$(dirname "$0")/config.sh"

logstepstart "Starting Step 1: Getting IDs"

mkdir -p "$OUTDIR/genomes"
cd "$OUTDIR" || exit 1

log "Downloading mapping file..."
# Download and process mapping file
wget -qO - 'https://rest.uniprot.org/proteomes/stream?compressed=true&fields=upid%2Corganism%2Corganism_id%2Cgenome_assembly&format=tsv&query=%28*%29' | gzip -d > mapping.txt


log "Creating whitelist of eukaryotes to include..."
WHITELIST_FILE="../eukaryotes_whitelist.txt"

if [ ! -f "$WHITELIST_FILE" ]; then
  log "Did not find whitelist file: creating default eukaryotes whitelist at $WHITELIST_FILE"
  cat > "$WHITELIST_FILE" << EOF
Caenorhabditis elegans
Drosophila melanogaster
Canis lupus
Felis catus
Gallus gallus
Mus musculus
Rattus norvegicus
Macaca mulatta
Danio rerio
Chlorocebus sabaeus
Arabidopsis thaliana
Saccharomyces cerevisiae
Fusarium sporotrichioides
Komagataella phaffii (strain GS115 / ATCC 20864) (Yeast) (Pichia pastoris)
Xenopus laevis
Pan troglodytes
Bos taurus
Sus scrofa
Candida albicans (strain SC5314 / ATCC MYA-2876) (Yeast)
E coli K12
Escherichia coli (strain K12)
Shigella flexneri (301 / Serotype 2a)
Cricetulus griseus
Aspergillus fumigatus (strain ATCC MYA-4609 / CBS 101355 / FGSC A1100 / Af293) (Neosartorya fumigata)
Pseudomonas aeruginosa (strain ATCC 15692 / DSM 22644 / CIP 104116 / JCM 14847 / LMG 12228 / 1C / PRS 101 / PAO1)
Pseudomonas syringae pv. tomato (strain ATCC BAA-871 / DC3000)
Alkalihalophilus pseudofirmus (strain ATCC BAA-2126 / JCM 17055 / OF4) (Bacillus pseudofirmus)
Bacillus cereus (strain ATCC 14579 / DSM 31 / CCUG 7414 / JCM 2152 / NBRC 15305 / NCIMB 9373 / NCTC 2599 / NRRL B-3711)
Bacillus subtilis (strain 168)
Lactiplantibacillus plantarum (strain ATCC BAA-793 / NCIMB 8826 / WCFS1)
Geobacillus kaustophilus (strain HTA426)
Mycoplasma genitalium (strain ATCC 33530 / DSM 19775 / NCTC 10195 / G37)
Mycoplasma mycoides subsp. mycoides SC (strain CCUG 32753 / NCTC 10114 / PG1)
EOF
  log "Default whitelist created. You can edit $WHITELIST_FILE to customize organisms."
  else
    log "Using existing eukaryotes whitelist: $WHITELIST_FILE"
fi

log "Combining proteome files..."
# Combine proteome files and extract IDs
cat <(wget -O - "${RPG75}") <(wget -O - "${VIRUS95}") > full_proteomes.txt

log "Filtering organisms..."
grep -v "Euk/" full_proteomes.txt > non_eukaryotes.txt

# For the eukaryotes, only keep those in the whitelist
grep "Euk/" full_proteomes.txt > all_eukaryotes.txt

# Process each whitelisted organism
> whitelisted_eukaryotes.txt
while IFS= read -r org; do
  [[ -z "$org" || "$org" =~ ^[[:space:]]*# ]] && continue
  grep "$org" all_eukaryotes.txt >> whitelisted_eukaryotes.txt || true
done < "$WHITELIST_FILE"

# Combine the results, the line below is old code
cat non_eukaryotes.txt whitelisted_eukaryotes.txt > filtered_proteomes.txt

# Extract IDs from filtered proteomes 
grep '^>' filtered_proteomes.txt | awk '{ sub(/^>/, ""); print $1 }' > ids.txt  #commented out sept 8


log "Creating mapped.db..."
# Create mapped.db
awk 'NR==FNR{a[$1];next} $1 in a {print $NF > "mapped.db"}' ids.txt mapping.txt


logstepend "Step 1 completed successfully"
