import pandas as pd
import json
import argparse

"""Only takes in ONE JSON FILE. Not a folder"""
'''python3 JSON2tsv-just_assem_out.py /Users/christiewoodside/desktop/argos/Campylobactercoli_test/all/coli_originalsample-qcAll.json /Users/christiewoodside/desktop/argos/Campylobactercoli_test/all/coli_original_assembly.tsv'''

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
        value = d.get(key, '')  # Default to empty string if the key doesn't exist
        
        # Optional: Format numeric values to 4 decimal places
        if isinstance(value, float):
            value = f"{value:.4f}"

        l.append(value)
    return l

def json_to_tsv(json_file, output_tsv):
    # Read JSON data from file
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    # Access the relevant part of the JSON
    if "refseq" in data:
        data = data["refseq"]
    
    # Convert each JSON object to a list based on schema keys
    rows = [listify(item, schema_keys) for item in data]
    
    # Convert the list of rows to a DataFrame
    df = pd.DataFrame(rows, columns=schema_keys)
    
    # Save DataFrame to TSV file
    df.to_csv(output_tsv, sep='\t', index=False)
    
    print(f"TSV file saved as {output_tsv}")

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument("json_file", help="Input JSON file")
    parser.add_argument("output_tsv", help="Output TSV file")
    args = parser.parse_args()
    
    # Convert JSON to TSV
    json_to_tsv(args.json_file, args.output_tsv)
