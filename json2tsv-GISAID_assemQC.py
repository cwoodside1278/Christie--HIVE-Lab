#Christie Woodside
#Oct ober 2024

'''This code takes in a folder of assemblyQC jsons and spits out the json information as a tsv. 
Since our GISAID assemblies are uploaded there is no recorded assembly_genome_id. This mean the 
other assembly codes will not work. This allos for a folder of JSONs to be formatted as is with no extra info'''

#Command line argument options
'''python3 json2tsv-GISAID_assemQC.py /Users/christiewoodside/desktop/argos/end_sept_push/assembly/ /Users/christiewoodside/desktop/argos/end_sept_push/assembly/justflu_assembly.tsv 
TSV file saved as /Users/christiewoodside/desktop/argos/end_sept_push/assembly/justflu_assembly.tsv'''

'''python3 json2tsv-GISAID_assemQC.py --json_folder /Users/christiewoodside/desktop/argos/oct_H5N1/all/ 
--output_tsv /Users/christiewoodside/desktop/argos/oct_H5N1/H5N1_texas_assembly_v1.tsv'''


import os
import pandas as pd
import json
import argparse
# import sys
# from Bio import Entrez
# import re
# import time

# sleeptime = 0.12
# sleeptime_withtoken = 0.11
# sleeptime_notoken = 0.34
# sep = '\t'

schema_keys = [
    "organism_name", 
    "infraspecific_name",
    "assembled_genome_acc",
    "genome_assembly_id",
    "representative_genome_acc",
    "representative_genome_org",
    "representative_genome_uniprot_acc",
    "lineage",
    "taxonomy_id",
    "bco_id",
    "schema_version",
    "analysis_platform",
    "analysis_platform_object_id",
    "assembly_file_source",    
    "genomic_section", 
    "num_chromosomes",
    "num_genes",    
    "assembly_gc_content",    
    "length",
    "size_gaps",
    "size_contigs",
    "contig_percentile",
    "contig_momentum",
    "coverage_contigs",
    "cnt_contigs",
    "cnt_gaps",
    "gap_percentile",
    "genome_coverage",
    "n50",
    "n75", 
    "n90",
    "n95",
    "l50",     
    "l75",
    "l90",
    "l95",
    "phred_average",
    "count_major_mutations",
    "count_major_indels",
    "mutation_momentum",
    "indels_momentum",
    "major_mutation_momentum",
    "major_indels_momentum",
    "alignment_anisotropy",
    "overhang_momentum",
    "aligned_momentum",
    "entropic_momentum",
    "reads_unaligned",
    "reads_aligned",
    "percent_reads_unaligned",
    "percent_reads_aligned",
    "assembly_level",
    "rpkm"
]

def listify(d, key_order):
    l = []
    for key in key_order:
        
        if key == "schema_version":
            value = 'v1.6'

        elif key == "bco_id":
            value = 'ARGOS_000086'

        elif key == "genome_assembly_id":  # Change here
            # Extract the EPI ID part after the '|' character from assembled_genome_acc
            value = d.get("assembled_genome_acc", '')  # Get value from the original key
            if '|' in value:
                value = value.split('|')[-1]  # Get the part after '|'
        elif key == "assembled_genome_acc":
            value = ''
        
        elif key == "assembly_file_source":
            value = 'GISAID'
        
        else:
            value = d.get(key, '')
        
        if isinstance(value, float):
            value = f"{value:.4f}"

        l.append(value)
    return l

def process_json_file(json_file):
    # Read JSON data from file
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    # Access the relevant part of the JSON
    if "refseq" in data:
        data = data["refseq"]
    
    # Convert each JSON object to a list based on schema keys
    rows = [listify(item, schema_keys) for item in data]
    
    return rows

def json_to_tsv(json_folder, output_tsv):
    # List all JSON files in the folder
    json_files = [os.path.join(json_folder, file) for file in os.listdir(json_folder) if file.endswith('.json')]
    
    # Initialize an empty list to collect rows from all files
    all_rows = []
    
    # Process each JSON file
    for json_file in json_files:
        all_rows.extend(process_json_file(json_file))
    
    # Convert the list of rows to a DataFrame
    df = pd.DataFrame(all_rows, columns=schema_keys)

    # Count occurrences of genome_assembly_id
    counts = df['genome_assembly_id'].value_counts().to_dict()
    # Map the counts back to the num_chromosomes column
    df['num_chromosomes'] = df['genome_assembly_id'].map(counts)
    
    # Save DataFrame to TSV file
    df.to_csv(output_tsv, sep='\t', index=False)
    
    print(f"\nTSV file saved as {output_tsv}\n")

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--json_folder", help="Folder containing JSON files")
    parser.add_argument("--output_tsv", help="Output TSV file")
    # parser.add_argument('--email', help='Email address associated with the NCBI account', type=str)
    # parser.add_argument('--api',
    #                     help='API key associated with the NCBI account (optional)',
    #                     type=str) 
    # if len(sys.argv) <= 1:
    #     sys.argv.append("--help")

    options = parser.parse_args()
    # Entrez.email = options.email
    # Entrez.api_key = options.api

    # # Set the API key if provided
    # if options.api:
    #     Entrez.api_key = options.api

    
    # Convert JSON to TSV
    json_to_tsv(options.json_folder, options.output_tsv)
