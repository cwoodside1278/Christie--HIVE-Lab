#Christie Woodside Jan 19, 2025
# Github handle: cwoodside1278
#
'''For APHIS H5N1 orgs specifically'''
"""json to tsv code and formatted schema v1.6 assemblyQC HIVE3. (qcAll.json)
This code will take the JSON file QC output for assemblyQC from schema v1.6 stored in a local folder and reformat it into a combined tsv. The output tsv will be pasted and aligned with the 
columns/headers for assemblyQC_HIVE (BCO ID ARGOS_000012) found in the data.argosdb.org dataset. This code was used for the datapush to ARGOSdb"""
"""only takes in the folder of jsons not one individual json as of June 24, 2024"""


import json, glob, os
import csv
import argparse
import sys
from Bio import Entrez
import re
import time

__version__ = "1.1.0"
__status__ = "Development"

sleeptime = 0.12
sleeptime_withtoken = 0.12
sleeptime_notoken = 0.34
argos_schema_version = 'v1.6'
sep = '\t'

# Note that this sleeptime is if authentication (a token) is provided. Use 0.34 if not authenticated.
# Suggesting 0.11 instead of 0.1 just to be safe! Getting banned by NCBI is a huge pain.
#This code was taken from biosample_datagrabber_v2.py

'''python3 APHIS_assemQC.py -t /Users/christiewoodside/desktop/thesis/H5N1/jan19/all/H5N1_BATCH1_assembly.tsv --schema /Users/christiewoodside/desktop/thesis/H5N1/jan19/all --email christie.woodside@email.gwu.edu --api bfbde99c962d228023e8d62a078bdb12d108'''

'''python3 APHIS_assemQC.py --email christie.woodside@email.gwu.edu --api bfbde99c962d228023e8d62a078bdb12d108 -t /Users/christiewoodside/desktop/thesis/H5N1/jan19/all/H5N1_BATCH1_assembly.tsv --schema /Users/christiewoodside/desktop/thesis/H5N1/jan19/all'''

columns_data = {
  "top_level": "refseq",
  "columns": [
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
    "percent_reads_unaligned",
    "percent_reads_aligned",
    "reads_aligned",
    "assembly_level",
    "rpkm"
  ],
  "header_map": {
    "bco_id": "id"
  }
}


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

    parser.add_argument('--api',
                        help='API key associated with the NCBI account (optional)',
                        type=str) 
    #--api bfbde99c962d228023e8d62a078bdb12d108
    #CW as of Oct 1 bfbde99c962d228023e8d62a078bdb12d108

    if len(sys.argv) <= 1:
        sys.argv.append("--help")

    options = parser.parse_args()
    Entrez.email = options.email
    #Entrez.api_key = options.api

    # Set the API key if provided
    if options.api:
        Entrez.api_key = options.api

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
def bsDataGet(as_term, sleeptime):
    '''Outputs the assembly id that is needed'''
    #print('as_term', as_term)
    search = Entrez.esearch(db='nucleotide', term=as_term, retmode='xml', idtype="acc")
    #print("heyyyyyyy")
    time.sleep(sleeptime)
    record = Entrez.read(search)
    time.sleep(sleeptime)
    
    if record["IdList"]:
        assembly_record_id = record['IdList'][0]
        #print(assembly_record_id)
        info = Entrez.efetch(db='nucleotide', id=assembly_record_id, rettype="gb", retmode="xml")
        time.sleep(sleeptime)
        record = Entrez.read(info)
        time.sleep(sleeptime)
        record_dict = record[0]
        #print(record_dict)
        
        for xref in record_dict["GBSeq_xrefs"]:
            if xref.get("GBXref_dbname") == "Assembly":
                return xref["GBXref_id"]
    return ""

def getLevel(l_term, sleeptime):
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
        
        # Extracting the GBSeq_definition value
        definition = record_dict.get("GBSeq_definition", "")
        if ',' not in definition:
            if 'chromosome' in definition:
                return 'chromosome'
            else:
                return ''
            
        if ',' in definition:
            # Split the definition by the comma and strip any leading/trailing spaces
            result = definition.split(',', 1)[1].strip()
            #print('assem level result: ',result)
            if result == 'whole genome shotgun sequence':
                return "contig whole genome shotgun sequence"
            elif result == 'complete sequence':
                return "complete sequence"
            elif result == 'complete genome':
                return "complete genome"
            elif result == 'contig':
                return 'contig'
            elif "plasmid" in definition:
                return "chromosome"
            elif 'chromosome' in definition:
                return "chromosome"
            # if result == '':
            #     return ''
            return result
        return ""
    return "No ID found"

def getGenomicSection(l_term, sleeptime):
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

        # Extracting the GBSeq_definition value
        definition = record_dict.get("GBSeq_definition", "")
        if ',' not in definition:
            if 'chromosome' in definition:
                return 'chromosome'
            else:
                return ''
            
        if ',' in definition:
            result = definition.split(',', 1)[0].strip()
            print('genomic result: ',result)
            # # if 'plasmid' in definition:
            # #     return 'plasmid'
            # if 'unitig_0_quiver_pilon' in definition:
            #     return 'contig'
            
            if 'Influenza A virus' in definition:
                new = result.split('))', 1)[1].strip()
                #print("\n\n",new)
                return new
            return result
    return "No ID found"

def getOrg(o_term, sleeptime):
    search = Entrez.esearch(db='nucleotide', term=o_term, retmode='xml', idtype="acc")
    #print(search, "\n\n")
    time.sleep(sleeptime)
    record = Entrez.read(search)
    #print(record)
    time.sleep(sleeptime)
    
    if record["IdList"]:
        assembly_record_id = record['IdList'][0]
        info = Entrez.efetch(db='nucleotide', id=assembly_record_id, rettype="gb", retmode="xml")
        time.sleep(sleeptime)
        record = Entrez.read(info)
        time.sleep(sleeptime)
        record_dict = record[0]
        return record_dict.get("GBSeq_organism", "")
    return ""

def getTax(t_term, sleeptime):
    search = Entrez.esearch(db='nucleotide', term=t_term, retmode='xml', idtype="acc")
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
        
        for feature in record_dict["GBSeq_feature-table"]:
            for qualifier in feature['GBFeature_quals']:
                if qualifier['GBQualifier_name'] == 'db_xref':
                    if "taxon:" not in qualifier['GBQualifier_value']:
                        return ''
                    return qualifier['GBQualifier_value'].replace("taxon:", "")
    return ""

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

#Doesn't work:
def getInfra(l_term, sleeptime):
    '''This gets that specific name for the influenza A like
    (A/cattle/Texas/24-009290-006/2024(H5N1))'''
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

        # Extracting the GBSeq_definition value
        definition = record_dict.get("GBSeq_definition", "")
        if ',' not in definition:
            return ''
            
        if ',' in definition:
            result = definition.split(',', 1)[0].strip()
            #print('genomic result: ',result)
            
            if 'Influenza A virus' in definition:
                new = result.split('virus', 1)[1].strip()
                #print("\n\n",new)
                infra = new.split('segment', 1)[0].strip()
                #print("\n\n",infra)
                return infra
            return result
    return "No ID found"

def getGene(g_term, sleeptime):
    """
    Fetches the total number of genes from the genome annotation data.
    """
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
        
        comment = record_dict.get("GBSeq_comment", "")
        if not comment:
            print(f"No comment field found for accession {assembly_record_id}")
            return "None Found"

        annotation_data = re.search(r'##Genome-Annotation-Data-START##(.*?)##Genome-Annotation-Data-END##', comment, re.DOTALL)
        
        if annotation_data:
            annotation_data_str = annotation_data.group(1).strip()
            parsed_data = {}
            for line in annotation_data_str.split(';'):
                if '::' in line:
                    key, value = line.split('::', 1)
                    parsed_data[key.strip()] = value.strip()
            return parsed_data.get('Genes (total)', 'Not Found')
        else:
            print(f"No annotation data found for accession {assembly_record_id}")
            return "None Found"
    else:
        print(f"No records found for term {g_term}")
        return "None Found"

def make_tsv(options):
    """
    This function writes the data to a tsv file.
    """
    #columns_data = json.load(open("./columns_assembly.json", "r"))

    def get_data_from_flat_item(flat_item, key):
        """Retrieve the data based on the key."""
        if key == "genome_assembly_id":
            return assembly_value
        elif key == "organism_name":
            return getOrg(flat_item.get("assembled_genome_acc", ""), sleeptime_withtoken) or ""
        elif key == "assembly_file_source":
            return 'NCBI'
        elif key == "taxonomy_id":
            return getTax(flat_item.get("assembled_genome_acc", ""), sleeptime_withtoken) or ""
        elif key == "lineage":
            return getLin(flat_item.get("assembled_genome_acc", ""), sleeptime_withtoken) or ""
        elif key == "infraspecific_name":
            return getInfra(flat_item.get("assembled_genome_acc", ""), sleeptime_withtoken) or ""
        elif key == "num_genes":
            return getGene(flat_item.get("assembled_genome_acc", ""), sleeptime_withtoken) or ""
        elif key == "schema_version":
            return 'v1.6'      #manually adding this in here
        elif key == "bco_id":
            return 'ARGOS_000086'   #manually adding this in
        elif key == "assembly_level":
            return getLevel(flat_item.get("assembled_genome_acc", ""), sleeptime_withtoken) or ""
        elif key == "num_chromosomes":
            analysis_platform_object_id = flat_item.get("analysis_platform_object_id", "")
            return id_counts.get(analysis_platform_object_id, 0)
        elif key == 'genomic_section':
            return getGenomicSection(flat_item.get("assembled_genome_acc", ""), sleeptime_withtoken) or ""
        else:
            return flat_item.get(columns_data["header_map"].get(key, key), "")

    # First, count the occurrences of each analysis_platform_object_id
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
        # Process each JSON file
        for schema in glob.glob(os.path.join(options.schema, '*.json')):
            with open(schema, "r") as jsonfile:
                data = json.load(jsonfile)
                assembly_value = data['assembly']
                #assembly_value = data.get("assembly", "")
                #print("here:                        ",assembly_value)
                top_level_data = data.get(columns_data["top_level"], [])
                for item in top_level_data:
                    # #print(f"Processing item: {item}")
                    # flat_item = flatten_json(item)
                    # row = []

                    # # Extract `assembly` value
                    # assembly_value = data.get("assembly", "")
                    
                # for item in data[columns_data["top_level"]]:
                    flat_item = flatten_json(item)
                    row = [get_data_from_flat_item(flat_item, key) for key in columns_data["columns"]]
                    # Extract `assembly` value
                    assembly_value = data.get("assembly", "")
                    writer.writerow(row)

def main():
    """
    Main function
    """

    options = usr_args()
    Entrez.email = options.email
    Entrez.api_key = options.api
    make_tsv(options)


# ______________________________________________________________________________#
if __name__ == "__main__":
    options = usr_args()
    Entrez.email = options.email
    Entrez.api_key = options.api
    main()