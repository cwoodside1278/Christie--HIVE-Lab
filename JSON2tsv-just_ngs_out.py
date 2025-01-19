import pandas as pd
import json
import argparse

'''Spits out the contents of the NGS JSON file as is. Does not pull from EUtils or any of the IDs. 
Only need to use one JSON not a file'''

'''python3 JSON2tsv-just_ngs_out.py /Users/christiewoodside/desktop/argos/dec_3/o48075-qcNGS.json /Users/christiewoodside/desktop/argos/dec_3/ngs/sudanEbolavirusTEST_ngs.tsv /Users/christiewoodside/Desktop/ARGOS/code/HIVE3/columns_ngs_justNGS_code.json
'''

#Don't even use these
schema_key = [
    "organism_name",
    "infraspecific_name",
    "lineage",
    "representative_genome_org",
    "representative_genome_acc",
    "representative_genome_uniprot_acc",
    "genome_assembly_id",
    "taxonomy_id",
    "bco_id",
    "schema_version",
    "analysis_platform",
    "analysis_platform_object_id",
    "strain",
    "bioproject",
    "biosample",
    "sra_run_id",
    "ngs_read_file_name",
    "ngs_read_file_source",
    "ngs_gc_content",
    "avg_phred_score",
    "min_read_length",
    "num_reads",
    "num_reads_unique",
    "avg_read_length",
    "max_read_length",
    "max_duplicate_read",
    "strategy",
    "instrument",
    "id_method",
    "codon_table",
    "percent_coding",
    "percent_not_coding",
    "complexity_percent",
    "non_complex_percent",
    "stdev_quality",
    "avg_quality_a",
    "avg_quality_t",
    "avg_quality_g",
    "avg_quality_c",
    "count_a",
    "count_c",
    "count_g",
    "count_t",
    "count_n",
    "percent_a",
    "percent_c",
    "percent_g",
    "percent_t",
    "count_all_wn",
    "count_all"
]

def listify(d, key_order, header_map):
    l = []
    #print(d, '\n')
    for key in key_order:
        #print(key, "-----")
        # Use the header_map to translate key if it exists
        mapped_key = header_map.get(key,key)  # If no mapping, use the original key
        #print(mapped_key, '\n')

        # First, try to get the value from the top-level dictionary
        value = d.get(mapped_key, '')
        #print(value, '\n')

        # If the key belongs to the nested 'bases' dictionary, look for it inside 'bases'
        if value == '' and mapped_key in ['count_a', 'count_c', 'count_g', 'count_t', 'count_n',
                                          'avg_quality_a', 'avg_quality_c', 'avg_quality_g', 'avg_quality_t', 
                                          'percent_a', 'percent_c', 'percent_g', 'percent_t', 'percent_n']:
            # Access the 'bases' dictionary if it exists
            bases_data = d.get('bases', {})
            #print(bases_data, '\n')
            value = bases_data.get(mapped_key, '')
            #print("value:  ", value)

            # if mapped_key == "analysis_platform":
            #     value = 'HIVE3'

        # Optional: Format numeric values to 4 decimal places
        if isinstance(value, float):
            value = f"{value:.4f}"

        l.append(value)
    return l


# def listify(d, key_order, header_map):
#     l = []
#     for key in key_order:
#         # Use the header_map to translate key if it exists
#         mapped_key = header_map.get(key, key)
        
#         # Get the value from the JSON object using the mapped key
#         value = d.get(mapped_key, '')  # Default to empty string if the key doesn't exist
        
#         # Optional: Format numeric values to 4 decimal places
#         if isinstance(value, float):
#             value = f"{value:.4f}"
        
#         l.append(value)
#     return l


def json_to_tsv(json_file, output_tsv, columns_data):
    # Read JSON data from file
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    # Access the relevant part of the JSON (assuming `top_level` points to the correct section)
    if columns_data["top_level"] in data:
        data = data[columns_data["top_level"]]
    
    # Extract the column names and the header map
    schema_keys = columns_data["columns"]
    header_map = columns_data.get("header_map", {})
    #print(header_map)
    
    # Convert each JSON object to a list based on schema keys
    rows = [listify(item, schema_keys, header_map) for item in data]
    #print(rows, '\n')
    
    # Convert the list of rows to a DataFrame
    df = pd.DataFrame(rows, columns=schema_keys)
    
    # Save DataFrame to TSV file
    df.to_csv(output_tsv, sep='\t', index=False)
    
    print(f"TSV file saved as {output_tsv}")


# def json_to_tsv(json_file, output_tsv, columns_data):
#     # Read JSON data from file
#     with open(json_file, 'r') as f:
#         data = json.load(f)
    
#     # Access the relevant part of the JSON (assuming `top_level` points to the correct section)
#     if columns_data["top_level"] in data:
#         data = data[columns_data["top_level"]]
    
#     # Extract the column names and the header map
#     schema_keys = columns_data["columns"]
#     header_map = columns_data.get("header_map", {})
    
#     # Convert each JSON object to a list based on schema keys
#     rows = [listify(item, schema_keys, header_map) for item in data]
    
#     # Convert the list of rows to a DataFrame
#     df = pd.DataFrame(rows, columns=schema_keys)
    
#     # Save DataFrame to TSV file
#     df.to_csv(output_tsv, sep='\t', index=False)
    
#     print(f"TSV file saved as {output_tsv}")


# if __name__ == "__main__":
#     # Set up argument parser
#     parser = argparse.ArgumentParser()
#     parser.add_argument("json_file", help="Input JSON file")
#     parser.add_argument("output_tsv", help="Output TSV file")
#     args = parser.parse_args()
    
#     # Convert JSON to TSV
#     json_to_tsv(args.json_file, args.output_tsv)

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument("json_file", help="Input JSON file")
    parser.add_argument("output_tsv", help="Output TSV file")
    parser.add_argument("columns_json", help="Path to the JSON file that contains column mappings")
    args = parser.parse_args()
    
    # Load the columns data (columns and header mapping) from JSON
    with open(args.columns_json, 'r') as f:
        columns_data = json.load(f)
    
    # Convert JSON to TSV
    json_to_tsv(args.json_file, args.output_tsv, columns_data)
