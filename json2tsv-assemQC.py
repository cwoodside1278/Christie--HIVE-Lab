# Updated for schema v1.6 August 20, 2024
# Christie Woodside. Acknowledgements: Hadley
# Github handle: cwoodside1278
#
"""json to tsv code and formatted schema v1.6 assemblyQC HIVE3. (qcAll.json)
This code will take the JSON file QC output for assemblyQC from schema v1.6 stored in a local folder and reformat it into a combined tsv. The output tsv will be pasted and aligned with the 
columns/headers for assemblyQC_HIVE (BCI ID ARGOS_000012) found in the data.argosdb.org dataset. This code was used for the datapush to ARGOSdb"""
"""only takes in the folder of jsons not one individual json as of June 24, 2024"""

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json, glob, os
import csv
import argparse
import sys
from Bio import Entrez
import re
import time

__version__ = "1.1.0"
__status__ = "Development"

sleeptime_withtoken = 0.11 # seconds
sleeptime_notoken = 0.34 # seconds
argos_schema_version = 'v1.6'
sep = '\t'

# Note that this sleeptime is if authentication (a token) is provided. Use 0.34 if not authenticated.
# Suggesting 0.11 instead of 0.1 just to be safe! Getting banned by NCBI is a huge pain.
#This code was taken from biosample_datagrabber_v2.py

'''Example entry: 
python3 json2tsv-assemQC.py -t test_assembly.tsv --schema /Users/christiewoodside/desktop/argos/data_push_june/all --email christie.woodside@gwu.edu'''

def usr_args():
    """
    functional arguments for process
    """
    parser = argparse.ArgumentParser()
    # set usages options
    parser = argparse.ArgumentParser(prog="argosdb", usage="%(prog)s [options]")
    # schema version. ex: v1.6
    parser.add_argument("-v", "--version", action="version", version="%(prog)s " + __version__)
    # The name of the output tsv file (and include .tsv extension)
    parser.add_argument("-t", "--tsv", help="tsv file to create.")
    # input assemblyQC (qcAll) json file directory
    parser.add_argument(
        "-s",
        "--schema",
        required=True,
        # type = argparse.FileType('r'),
        help="Root json schema to parse. (Input json file)",
    )
    parser.add_argument('--email',
                        help='email address associated with the NCBI account',
                        type=str)
    

    # Print usage message if no args are supplied.
    if len(sys.argv) <= 1:
        sys.argv.append("--help")

    options = parser.parse_args()
    Entrez.email = options.email
    return options


def flatten_json(y):
    '''Flattening the json input'''
    out = {}

    def flatten(x, name=""):
        if isinstance(x, dict):
            for a in x:
                flatten(x[a], name + a + "_")
        elif isinstance(x, list):
            i = 0
            for a in x:
                flatten(a, name + str(i) + "_")
                i += 1
        else:
            out[name[:-1]] = x

    flatten(y)
    return out

#from biosample_datagrabber_v2.py on GitHub but used to grab assembly genome accession
def bsDataGet(as_term, sleeptime): #removed bco_id value
    #Get the genome assembly id, not biosample id (reusing code) -----------------------------------------------------------------------
    search = Entrez.esearch(db = 'nucleotide', term = as_term, retmode='xml', idtype="acc")
    time.sleep(sleeptime)
    record = Entrez.read(search)
    time.sleep(sleeptime)
    #print(f"LINE 101: {record}\n")
    assembly_id = ""
    if record["IdList"]:
        #print(f"LINE 93: {record}\n")
        assembly_record_id = record['IdList'][0]
        info = Entrez.efetch(db = 'nucleotide', id = assembly_record_id, rettype="gb", retmode="xml")
        time.sleep(sleeptime)
        record = Entrez.read(info)
        time.sleep(sleeptime)
        record_dict = record[0]
        
        #This is to get the genome assembly id
        for xref in record_dict["GBSeq_xrefs"]:
            if xref.get("GBXref_dbname") == "Assembly":
                assembly_id = xref["GBXref_id"]
                return assembly_id
            
#This si to get the organism name
def getOrg(o_term, sleeptime): #removed bco_id value
    search = Entrez.esearch(db = 'nucleotide', term = o_term, retmode='xml', idtype="acc")
    time.sleep(sleeptime)
    record = Entrez.read(search)
    time.sleep(sleeptime)
    
    org_id = ""
    if record["IdList"]:
        assembly_record_id = record['IdList'][0] #it is because I copy and pasted
        info = Entrez.efetch(db = 'nucleotide', id = assembly_record_id, rettype="gb", retmode="xml")
        time.sleep(sleeptime)
        record = Entrez.read(info)
        time.sleep(sleeptime)
        record_dict = record[0]
        org_id= record_dict["GBSeq_organism"]
        return org_id
    
#to get the taxonomy id
def getTax(t_term, sleeptime):
    search = Entrez.esearch(db = 'nucleotide', term = t_term, retmode='xml', idtype="acc")
    time.sleep(sleeptime)
    record = Entrez.read(search)
    time.sleep(sleeptime)
    
    tax_id = ""
    if record["IdList"]:
        assembly_record_id = record['IdList'][0] #it is because I copy and pasted
        info = Entrez.efetch(db = 'nucleotide', id = assembly_record_id, rettype="gb", retmode="xml")
        time.sleep(sleeptime)
        record = Entrez.read(info)
        time.sleep(sleeptime)
        record_dict = record[0]
        
        for feature in record_dict["GBSeq_feature-table"]:
            for qualifier in feature['GBFeature_quals']:
                if qualifier['GBQualifier_name'] == 'db_xref':
                    q= qualifier['GBQualifier_value']
                    tax_id = q.replace("taxon:","")
                    #print(f"db_xref: {tax_id}")
                    return tax_id
                
#to get the lineage
def getLin(l_term, sleeptime): #removed bco_id value
    #Get the genome assembly id -----------------------------------------------------------------------
    search = Entrez.esearch(db = 'nucleotide', term = l_term, retmode='xml', idtype="acc")
    time.sleep(sleeptime)
    record = Entrez.read(search)
    time.sleep(sleeptime)
    lin_id = ""
    if record["IdList"]:
        assembly_record_id = record['IdList'][0] #copy and pasted
        info = Entrez.efetch(db = 'nucleotide', id = assembly_record_id, rettype="gb", retmode="xml")
        time.sleep(sleeptime)
        record = Entrez.read(info)
        time.sleep(sleeptime)
        record_dict = record[0]
        lin_id = record_dict["GBSeq_taxonomy"]
        return lin_id

def getGene(g_term, sleeptime):
    search = Entrez.esearch(db='nucleotide', term=g_term, retmode='xml', idtype="acc")
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
        
        # Get the comment field from the record
        comment = record_dict.get("GBSeq_comment", "")

        # Extract genome-annotation-data section directly
        sections = {}
        # Extract genome annotation data
        annotation_data = re.search(r'##Genome-Annotation-Data-START##(.*?)##Genome-Annotation-Data-END##', comment, re.DOTALL)
        if annotation_data:
            annotation_data_str = annotation_data.group(1).strip()
            # Parse the annotation data directly
            parsed_data = {}
            lines = annotation_data_str.split(';')
            for line in lines:
                if '::' in line:
                    key, value = line.split('::', 1)
                    key = key.strip()
                    value = value.strip()
                    parsed_data[key] = value
            # Extract the desired value
            genes_total = parsed_data.get('Genes (total)', 'Not Found')
            return genes_total
        
def getBP(bp_term, sleeptime): #removed bco_id value
    #Get the BioProject associated with these samples -----------------------------------------------------------------------
    search = Entrez.esearch(db = 'nucleotide', term = bp_term, retmode='xml', idtype="acc")
    time.sleep(sleeptime)
    record = Entrez.read(search)
    time.sleep(sleeptime)
    bp_id = ""

    if record["IdList"]:
        assembly_record_id = record['IdList'][0]
        info = Entrez.efetch(db = 'nucleotide', id = assembly_record_id, rettype="gb", retmode="xml")
        time.sleep(sleeptime)
        record = Entrez.read(info)
        time.sleep(sleeptime)
        record_dict = record[0]
        
        #This is to get the genome assembly id
        for xref in record_dict["GBSeq_xrefs"]:
            if xref.get("GBXref_dbname") == "BioProject":
                bp_id = xref["GBXref_id"]
                return bp_id
        

def make_tsv(options):
    """
    This function writes the data to a tsv file
    """
    columns_data = json.load(open("./columns_assembly.json", "r")) #must have columns_assembly.json file stored locally on your computer to work. It is specific to each GC json output
    with open(options.tsv, "w", newline="") as tsvfile:
        writer = csv.writer(tsvfile, delimiter="\t")

        # Write the header row
        writer.writerow(columns_data["columns"])
         #take in a folder of jsons to be combined together as a tsv
        json_files = glob.glob(os.path.join(options.schema, '*.json')) 
    
        for schema in json_files:
            with open(schema, "r") as jsonfile:
                data = json.load(jsonfile)

                for item in data[columns_data["top_level"]]:
                    flat_item = flatten_json(item)
                    row = []
                    for key in columns_data["columns"]:
                        if key == "genome_assembly_id":
                            # Fetch the assembly_id using bsDataGet
                            a_id = flat_item.get("assembled_genome_acc", "")
                            assembly_id = bsDataGet(a_id, sleeptime_withtoken)
                            row.append(assembly_id if assembly_id else "-")

                        elif key == "organism_name":
                            o_id = flat_item.get("assembled_genome_acc", "")#it needs the assembled genome accession to access the json correctly and grab the org name
                            org_id = getOrg(o_id, sleeptime_withtoken)
                            row.append(org_id if org_id else "-")

                        elif key == "taxonomy_id":
                            t_id = flat_item.get("assembled_genome_acc", "")#it needs the assembled genome accession to access the json correctly and grab the org name
                            t_id = getTax(t_id, sleeptime_withtoken)
                            row.append(t_id if t_id else "-")

                        elif key == "lineage":
                            l_id = flat_item.get("assembled_genome_acc", "")#it needs the assembled genome accession to access the json correctly and grab the org name
                            l_id = getLin(l_id, sleeptime_withtoken)
                            row.append(l_id if l_id else "-")

                        elif key == "num_genes":
                            l_id = flat_item.get("assembled_genome_acc", "")#it needs the assembled genome accession to access the json correctly and grab the org name
                            l_id = getGene(l_id, sleeptime_withtoken)
                            row.append(l_id if l_id else "-")

                        elif key == "schema_version":
                            row.append('v1.6')
                        else:
                            if key in columns_data["header_map"]:
                                key = columns_data["header_map"][key]
                            row.append(flat_item.get(key, ""))
                    writer.writerow(row)


def main():
    """
    Main function
    """

    options = usr_args()
    Entrez.email = options.email
    make_tsv(options)


# ______________________________________________________________________________#
if __name__ == "__main__":
    options = usr_args()
    Entrez.email = options.email
    main()