#Newly Updated json2tsv-assemQC.py file as of Sept 23, 2024
# Updated for schema v1.6 August 06, 2024
# Christie Woodside. Acknowledgements: Hadley
# Github handle: cwoodside1278
#
'''This is for GISAID/external downloader submission'''

"""json to tsv code and formatted schema v1.6 assemblyQC HIVE3 (qcAll.json) and GISAID assemblies.

****Required the bs_text file used in biosampleQC, also tested this one all biosamples from one org. Not sure how great it will work with multiple orgs so double check always.

This code will take the JSON file QC output for assemblyQC from schema v1.6 stored in a local folder 
and reformat it into a combined tsv. The output tsv will be pasted and aligned with the 
columns/headers for assemblyQC_HIVE (BCO ID ARGOS_000012) found in the data.argosdb.org dataset."""

"""only takes in the folder of jsons not one individual json as of June 24, 2024"""



import json, glob, os
import csv
import argparse
import sys
from Bio import Entrez
import re
import time
import xmltodict

__version__ = "1.1.0"
__status__ = "Development"

sleeptime_withtoken = 0.11
sleeptime_notoken = 0.34
argos_schema_version = 'v1.6'
sep = '\t'

# Note that this sleeptime is if authentication (a token) is provided. Use 0.34 if not authenticated.
# Suggesting 0.11 instead of 0.1 just to be safe! Getting banned by NCBI is a huge pain.
#This code was taken from biosample_datagrabber_v2.py

'''Example entry: 
python3 json2tsv-GISAID_assemQC.py -t /Users/christiewoodside/desktop/argos/end_sept_push/assembly/test_assembly.tsv --schema /Users/christiewoodside/desktop/argos/end_sept_push/assembly 
--email christie.woodside@email.gwu.edu --bs /Users/christiewoodside/desktop/argos/end_sept_push/biosample/flu_BS.txt  '''

def usr_args():
    """
    Functional arguments for process
    """
    parser = argparse.ArgumentParser(prog="argosdb", usage="%(prog)s [options]")
    # schema version. ex: v1.6
    parser.add_argument("-v", "--version", action="version", version="%(prog)s " + __version__)
    # The name of the output tsv file (and include .tsv extension)
    parser.add_argument("-t", "--tsv", help="TSV file to create.")
    # input assemblyQC (qcAll) json file directory
    parser.add_argument("-s", "--schema", required=True, help="Root JSON schema to parse. (Input JSON file)")
    #So NCBI doesn't ban us
    parser.add_argument('--email', help='Email address associated with the NCBI account', type=str)
    #the biosample file
    parser.add_argument('--bs', help='Text file containing biosample IDs (one per line)', type=str, required=True)

    if len(sys.argv) <= 1:
        sys.argv.append("--help")

    options = parser.parse_args()
    Entrez.email = options.email
    return options


def flatten_json(y):
    '''Flattening the JSON input'''
    out = {}

    def flatten(x, name=""):
        if isinstance(x, dict):
            for a in x:
                flatten(x[a], name + a + "_")
        elif isinstance(x, list):
            for i, a in enumerate(x):
                flatten(a, name + str(i) + "_")
        else:
            out[name[:-1]] = x

    flatten(y)
    return out
#__________________________________________________________________________________________________________________________________

def read_biosample_ids(file_path):
    """Read biosample IDs from the text file"""
    with open(file_path, 'r') as f:
        biosample_ids = [line.strip() for line in f]
        #print('      read_biosample:  ', biosample_ids)
    return biosample_ids

def match_biosample_to_json(biosample_id, json_files):
    """Match the biosample ID to the correct JSON file based on its filename."""
    for json_file in json_files:
        if biosample_id in json_file:
            return json_file
    return None  # No matching file found

def getGISAID(G_id):
    '''Get the GISAID id from the file itself'''
    #print('getGISAID id input:   ', G_id)
    if '|' in G_id:
        result = G_id.split('|', 1)[1].strip()
        return result
    

def bsDataGet(bs_term, sleeptime):
    '''uses biosample file to get information Outputs the assembly id that is needed'''
    search = Entrez.esearch(db='biosample', term=bs_term, retmode='xml')
    time.sleep(sleeptime)
    record = Entrez.read(search)
    time.sleep(sleeptime)
    bs_id = record['IdList'][0]
    info = Entrez.esummary(db='biosample', id=bs_id)
    time.sleep(sleeptime)

    record = Entrez.read(info)
    time.sleep(sleeptime)
    r = record['DocumentSummarySet']['DocumentSummary'][0]
    sd = r['SampleData']
    sd_json = xmltodict.parse(sd)['BioSample']

    attr = sd_json['Attributes']['Attribute']
    attr_set = {}
    for att in attr:
        attr_set[att['@attribute_name']] = att['#text']

        if att['@attribute_name'] == 'isolate':
            isolate_value = att['#text']
            attr_set['infraspecific_name'] = isolate_value
            #print("Isolate: ", isolate_value)

    attr_set['organism_name'] = sd_json['Description']['Organism']['OrganismName']
    #print("org name:     ", attr_set['organism_name'])
    attr_set['taxonomy_id'] = sd_json['Description']['Organism']['@taxonomy_id']

    lineage = getLin(bs_term, sleeptime) #some samples may not have a genome and therefore no lineage. Look it up on NCBI Taxonomy
    attr_set['lineage'] = lineage
    # attr_set['lineage'] = sd_json['Description']['Organism']
    # print(sd_json)
    #print(sd_json['Attributes']['Attribute'])

    return attr_set


def getLin(l_term, sleeptime):
    search = Entrez.esearch(db='nucleotide', term=l_term, retmode='xml', idtype="acc")
    time.sleep(sleeptime)
    record = Entrez.read(search)
    time.sleep(sleeptime)
    
    if record["IdList"]:
        assembly_record_id = record['IdList'][0]
        info = Entrez.efetch(db='nucleotide', id=assembly_record_id, rettype="gb", retmode="xml")
        time.sleep(sleeptime)
        record = Entrez.read(info)
        time.sleep(sleeptime)
        record_dict = record[0]
        return record_dict.get("GBSeq_taxonomy", "")
    return ""

def make_tsv(options, biosample_ids):
    """
    This function writes the data to a tsv file.
    """
    columns_data = json.load(open("./columns_assembly.json", "r"))
    json_files = glob.glob(os.path.join(options.schema, '*.json'))

    def get_data_from_flat_item(flat_item, key, attr_set):
        """Retrieve the data based on the key."""
        if key in attr_set:
            return attr_set[key]
        if key == "genome_assembly_id":
            return getGISAID(flat_item.get("assembled_genome_acc", ""))
        if key == "assembled_genome_acc":
            return ""
        elif key == "assembly_file_source":
            return 'GISAID'
        elif key == "schema_version":
            return 'v1.6'      #manually adding this in here
        elif key == "bco_id":
            return 'ARGOS_000012'   #manually adding this in
        else:
            return flat_item.get(columns_data["header_map"].get(key, key), "")

    
    id_counts = {}
    for schema in glob.glob(os.path.join(options.schema, '*.json')):
        with open(schema, "r") as jsonfile:
            data = json.load(jsonfile)
            for item in data[columns_data["top_level"]]:
                flat_item = flatten_json(item)
                analysis_platform_object_id = flat_item.get("analysis_platform_object_id", "")
                if analysis_platform_object_id:
                    if analysis_platform_object_id in id_counts:
                        id_counts[analysis_platform_object_id] += 1
                    else:
                        id_counts[analysis_platform_object_id] = 1

    with open(options.tsv, "w", newline="") as tsvfile:
        writer = csv.writer(tsvfile, delimiter="\t")
        # Write the header row
        writer.writerow(columns_data["columns"])

        for biosample_id in biosample_ids:
            # Get the biosample data
            print(f"Processing biosample: {biosample_id}")
            attr_set = bsDataGet(biosample_id, sleeptime_withtoken)
            
            # Find the correct JSON file for the biosample
            matched_json = match_biosample_to_json(biosample_id, json_files)
            if matched_json:
                print(f"Matching JSON file found: {matched_json}")
                with open(matched_json, "r") as jsonfile:
                    data = json.load(jsonfile)
                    for item in data[columns_data["top_level"]]:
                        flat_item = flatten_json(item)

                        # Generate the row, pulling data from flat_item and attr_set
                        row = [get_data_from_flat_item(flat_item, key, attr_set) for key in columns_data["columns"]]
                        writer.writerow(row)
            else:
                print(f"No matching JSON file found for biosample: {biosample_id}")


def main():
    """
    Main function
    """

    options = usr_args()
    biosample_ids = read_biosample_ids(options.bs)

    for bs_term in biosample_ids:
        bsDataGet(bs_term, sleeptime_withtoken)
        # if bs_id:
        #     print(f"BS ID for {bs_term}: {bs_id}")

    Entrez.email = options.email
    make_tsv(options, biosample_ids)


# ______________________________________________________________________________#
if __name__ == "__main__":
    options = usr_args()
    Entrez.email = options.email
    main()