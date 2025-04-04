#!/bin/bash

#input: bash assemblyQC_APHIS_datatable.sh /Users/christiewoodside/Desktop/Thesis/H5N1/mar19/all/ /Users/christiewoodside/Desktop/Thesis/H5N1/mar19/all/test_newcode.tsv

# Created March 20, 2025
#This code is used to create the assemblyQC_APHIS* tables for my thesis.. This was because there were server issues with NCBI and the API is just faster


# Check if jq is installed
if ! command -v jq &> /dev/null; then
    echo "jq is required but not installed. Install jq and try again."
    exit 1
fi

# Input directory and output file (output file is optional, default will be used if not provided)
input_dir="$1"
output_file="${2:-'combined_data.tsv'}"  # Use provided output file, or default to 'ncbi_data.tsv'

# Check if the input directory is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <input_directory> [<output_file>]"
    exit 1
fi

# Check if the input directory exists
if [ ! -d "$input_dir" ]; then
    echo "Input directory $input_dir does not exist. Exiting."
    exit 1
fi

# Ensure output file is empty before writing
> "$output_file"

#---------------------------------------------------------------------------------------------------------------


# Write headers to the TSV file manually (this is specifically for assemblyQC table)
echo -e "organism_name\tinfraspecific_name\tassembled_genome_acc\tgenome_assembly_id\trepresentative_genome_acc\trepresentative_genome_org\trepresentative_genome_uniprot_acc\tlineage\ttaxonomy_id\tbco_id\tschema_version\tanalysis_platform\tanalysis_platform_object_id\tassembly_file_source\tgenomic_section\tnum_chromosomes\tnum_genes\tassembly_gc_content\tlength\tsize_gaps\tsize_contigs\tcontig_percentile\tcontig_momentum\tcoverage_contigs\tcnt_contigs\tcoverage_gaps\tcnt_gaps\tgap_percentile\tgenome_coverage\tn50\tn75\tn90\tn95\tl50\tl75\tl90\tl95\tphred_average\tcount_major_mutations\tcount_major_indels\tmutation_momentum\tindels_momentum\tmajor_mutation_momentum\tmajor_indels_momentum\talignment_anisotropy\toverhang_momentum\taligned_momentum\tentropic_momentum\treads_unaligned\tpercent_reads_unaligned\tpercent_reads_aligned\treads_aligned\tassembly_level\trpkm" > "$output_file"

# Loop through all JSON files in the provided folder (nested directory structure) matching *-qcNGS.json
for json_file in "$input_dir"*-qcAll.json; do
    echo "Grabbing values from $json_file ..."
    echo ""

    #Will be the same genome assembly ID for all rows so putting it here
    GAID=$(jq -r '.assembly // "NA"' "$json_file")
    #echo "$GAID"

    jq -c '.refseq[]' "$json_file" | while read -r entry; do
        # Extract assembled_genome_acc which is the nucleotide ID
        nucleotide=$(echo "$entry" | jq -r '.assembled_genome_acc // "NA"')
        #echo "$nucleotide"

        # Extract values specific to this assembled_genome_acc
        analysis_platform_object_id=$(echo "$entry" | jq -r '.analysis_platform_object_id // "NA"')
        analysis_platform=$(echo "$entry" | jq -r '.analysis_platform // "NA"')
        bco_id=$(echo "$entry" | jq -r '.id // "NA"')
        length=$(echo "$entry" | jq -r '.length // "NA"')
        size_gaps=$(echo "$entry" | jq -r '.size_gaps // "NA"')
        size_contigs=$(echo "$entry" | jq -r '.size_contigs // "NA"')
        genome_coverage=$(echo "$entry" | jq -r '.genome_coverage // "NA"')
        coverage_contigs=$(echo "$entry" | jq -r '.coverage_contigs // "NA"')
        coverage_gaps=$(echo "$entry" | jq -r '.coverage_gaps // "NA"')
        cnt_contigs=$(echo "$entry" | jq -r '.cnt_contigs // "NA"')
        cnt_gaps=$(echo "$entry" | jq -r '.cnt_gaps // "NA"')
        contig_percentile=$(echo "$entry" | jq -r '.contig_percentile // "NA"')
        gap_percentile=$(echo "$entry" | jq -r '.gap_percentile // "NA"')
        contig_momentum=$(echo "$entry" | jq -r '.contig_momentum // "NA"')

        n50=$(echo "$entry" | jq -r '.n50 // "NA"')
        l50=$(echo "$entry" | jq -r '.l50 // "NA"')
        n75=$(echo "$entry" | jq -r '.n75 // "NA"')
        l75=$(echo "$entry" | jq -r '.l75 // "NA"')
        n90=$(echo "$entry" | jq -r '.n90 // "NA"')
        l90=$(echo "$entry" | jq -r '.l90 // "NA"')
        n95=$(echo "$entry" | jq -r '.n95 // "NA"')
        l95=$(echo "$entry" | jq -r '.l95 // "NA"')

        assembly_gc_content=$(echo "$entry" | jq -r '.assembly_gc_content // "NA"')
        phred_average=$(echo "$entry" | jq -r '.phred_average // "NA"')
        count_major_mutations=$(echo "$entry" | jq -r '.count_major_mutations // "NA"')
        count_major_indels=$(echo "$entry" | jq -r '.count_major_indels // "NA"')

        mutation_momentum=$(echo "$entry" | jq -r '.mutation_momentum // "NA"')
        indels_momentum=$(echo "$entry" | jq -r '.indels_momentum // "NA"')
        major_mutation_momentum=$(echo "$entry" | jq -r '.major_mutation_momentum // "NA"')
        major_indels_momentum=$(echo "$entry" | jq -r '.major_indels_momentum // "NA"')

        alignment_anisotropy=$(echo "$entry" | jq -r '.alignment_anisotropy // "NA"')
        overhang_momentum=$(echo "$entry" | jq -r '.overhang_momentum // "NA"')
        aligned_momentum=$(echo "$entry" | jq -r '.aligned_momentum // "NA"')
        entropic_momentum=$(echo "$entry" | jq -r '.entropic_momentum // "NA"')

        reads_unaligned=$(echo "$entry" | jq -r '.reads_unaligned // "NA"')
        reads_aligned=$(echo "$entry" | jq -r '.reads_aligned // "NA"')
        percent_reads_aligned=$(echo "$entry" | jq -r '.percent_reads_aligned // "NA"')
        percent_reads_unaligned=$(echo "$entry" | jq -r '.percent_reads_unaligned // "NA"')
        rpkm=$(echo "$entry" | jq -r '.rpkm // "NA"')


        # Grabbing Metadata from the NCBI APIS below---------------------------------------------------------------------------------------------


        # Query the assembly database using eutils API and $GAID which is the assemblyID
        SEARCH_RESULT=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&term=$GAID&retmode=json&api_key=bfbde99c962d228023e8d62a078bdb12d108")
        ASSEM_ID=$(echo "$SEARCH_RESULT" | jq -r '.esearchresult.idlist[0] // empty')
        #echo ""
        #echo "ASSEMBLY API Response: $SEARCH_RESULT"

        if [[ -n "$ASSEM_ID" ]]; then

            # Query the ASSEMBLY metadata to get more information to fill out the table
            ASSEM_METADATA=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&id=$ASSEM_ID&retmode=xml&api_key=bfbde99c962d228023e8d62a078bdb12d108")
            
            # Decode HTML entities in the ASSEM_METADATA
            DECODED_ASSEM_METADATA=$(echo "$ASSEM_METADATA" | sed 's/&lt;/</g; s/&gt;/>/g; s/&amp;/&/g')
            echo "$DECODED_ASSEM_METADATA" >> assembly_metadata.txt

            # Values for the Table
            TAXID=$(echo "$DECODED_ASSEM_METADATA" | xmllint --xpath "//*[local-name()='Taxid']/text()" - 2>/dev/null)
            #ORGANISM=$(echo "$DECODED_ASSEM_METADATA" | xmllint --xpath "//*[local-name()='Organism']/text()" - 2>/dev/null) #in my other assemblyQC file this is 'SpeciesName'
            ORGANISM=$(echo "$DECODED_ASSEM_METADATA" | xmllint --xpath "//*[local-name()='Organism']/text()" - 2>/dev/null | sed 's/ (.*)//')
            ASSEMBLY_STATUS=$(echo "$DECODED_ASSEM_METADATA" | xmllint --xpath "//*[local-name()='AssemblyStatus']/text()" - 2>/dev/null)
            

            # Extract the content of the Meta tag (which includes CDATA)
            META_CONTENT=$(echo "$DECODED_ASSEM_METADATA" | xmllint --xpath "//*[local-name()='Meta']/text()" - 2>/dev/null)
            # Remove CDATA wrapper by replacing the CDATA tags with their content
            CLEAN_META_CONTENT=$(echo "$META_CONTENT" | sed 's/<!\[CDATA\[\(.*\)\]\]>/\1/')
            # Now parse the cleaned META content
            # Ensure the content is valid XML by wrapping it in a proper root element
            CLEAN_META_XML="<root>$CLEAN_META_CONTENT</root>"
            # Extract chromosome count from the cleaned and wrapped META content
            CHROMOSOME_COUNT=$(echo "$CLEAN_META_XML" | xmllint --xpath 'string(//Stat[@category="chromosome_count" and @sequence_tag="all"])' -)
            #echo "chromosome counts: $CHROMOSOME_COUNT"
            

            # Print extracted values
            # echo "Taxid: $TAXID"
            # echo "Organism: $ORGANISM"
            # echo "Assembly Status: $ASSEMBLY_STATUS"
        fi





        # Query the nucleotide database using eutils API and $nucleotide which is the nucleotide assecion iD
        SEARCH_RESULT=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term=$nucleotide&retmode=json&api_key=bfbde99c962d228023e8d62a078bdb12d108")
        NUC_ID=$(echo "$SEARCH_RESULT" | jq -r '.esearchresult.idlist[0] // empty')
        #echo "NUC API Response: $SEARCH_RESULT"

        if [[ -n "$NUC_ID" ]]; then
            # Query the NUCLEOTIDE metadata to get more information to fill out the table
            NUC_METADATA=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=$NUC_ID&retmode=xml&api_key=bfbde99c962d228023e8d62a078bdb12d108")
            #echo "$NUC_METADATA"

            # Decode HTML entities in the NUC_METADATA -> some tags/f;ags were causing parsing errors which is why these are so long
            #DECODED_NUC_METADATA=$(echo "$NUC_METADATA" | sed 's/&lt;/</g; s/&gt;/>/g; s/&amp;/&/g' | sed 's/<GBFeature_location><[0-9]*\.\.[0-9]*<\/GBFeature_location>//g')
            DECODED_NUC_METADATA=$(echo "$NUC_METADATA" | sed 's/&lt;/</g; s/&gt;/>/g; s/&amp;/&/g' | sed -e 's/<GBFeature_location><[0-9]*\.\.[0-9]*<\/GBFeature_location>//g' -e 's/<GBFeature_location>complement(<[^<]*<\/GBFeature_location>//g')

            #gives error outputs
            echo "$DECODED_NUC_METADATA" | xmllint --noout -
            #echo "$DECODED_NUC_METADATA" >> nucleotide_metadata.txt


            # Extract values to add to the table--------
            #LINEAGE=$(echo "$DECODED_NUC_METADATA" | xmllint --xpath "string(//GBSeq_taxonomy)" - 2>/dev/null)
                #This one is more direct/specific than the above one
            LINEAGE=$(echo "$DECODED_NUC_METADATA" | xmllint --xpath "string(//GBSeq/GBSeq_taxonomy)" - 2>/dev/null)
            GENES_TOTAL=$(echo "$DECODED_NUC_METADATA" | grep -o '; Genes (total) :: [0-9,]\+' | sed 's/; Genes (total) :: //')
            GENE_SEC=$(echo "$DECODED_NUC_METADATA" | xmllint --xpath 'string(//GBSeq_definition)' - 2>/dev/null)

            #echo " "
            #echo "Raw GBSeq Definition: '$GENE_SEC'"
            #echo "Raw Lineage: '$LINEAGE'"


            #How I am getting the infraspecific name for the table
            STRAIN_INFRA=""
            if [[ -z "$GENE_SEC" ]]; then
                STRAIN_INFRA="  "
                echo "No GBSEq_Definition (for Genome_Section) found"
            
            elif [[ $GENE_SEC != *","* ]]; then
                echo " "
                echo "         THERE IS STRAIN NAME OH NO -------"
                STRAIN_INFRA=" "

            #if there is a , in the GBSeq_definition name (which there usually is for my flu  strains)
            else
                RESULT=$(echo $GENE_SEC | cut -d ',' -f1 | xargs)  # Extract first part before comma
                #echo "Infraspecific name before , result: $RESULT"
                echo "$RESULT"
                #echo ""

                if [[ "$GENE_SEC" == *"Influenza A virus"* ]]; then
                        NEW=$(echo "$RESULT" | sed -E 's/.*virus//;s/segment.*//;s/^[[:space:]]*//;s/[[:space:]]*$//')
                        STRAIN_INFRA="$NEW"
                else
                    STRAIN_INFRA="$RESULT"
                fi

            fi
            # echo "Infraspecific/Strain Name: $STRAIN_INFRA"
            # echo ""



            #This is how I am getting the genomic section still
            genomic_section=""

            if [[ -z "$GENE_SEC" ]]; then
                genomic_section=""
            elif [[ "$GENE_SEC" != *","* ]]; then  # If no comma in definition
                if [[ "$GENE_SEC" == *"chromosome"* ]]; then
                    genomic_section="chromosome"
                else
                    genomic_section=""
                fi
            else
                # Extract the first part before the comma
                result=$(echo "$GENE_SEC" | cut -d ',' -f1 | xargs)

                # Check for specific substrings
                if [[ "$GENE_SEC" == *"unitig_0_quiver_pilon"* ]]; then
                    genomic_section="contig"
                elif [[ "$GENE_SEC" == *"Influenza A virus"* ]]; then
                    new=$(echo "$result" | sed -E 's/.*\)\)//;s/^[[:space:]]*//;s/[[:space:]]*$//')  # Extract after "))"
                    genomic_section="$new"
                else
                    genomic_section="$result"
                fi
            fi

            #echo "The genomic section for this strain: $genomic_section"



        fi

         # Append extracted values to the output TSV file
        echo -e "$ORGANISM\t$STRAIN_INFRA\t$nucleotide\t$GAID\t\t\t\t$LINEAGE\t$TAXID\t\tv1.6\t$analysis_platform\t$analysis_platform_object_id\tNCBI\t${genomic_section:-'NA'}\t$CHROMOSOME_COUNT\t$GENES_TOTAL\t$assembly_gc_content\t$length\t$size_gaps\t$size_contigs\t$contig_percentile\t$contig_momentum\t$coverage_contigs\t$cnt_contigs\t$coverage_gaps\t$cnt_gaps\t$gap_percentile\t$genome_coverage\t$n50\t$n75\t$n90\t$n95\t$l50\t$l75\t$l90\t$l95\t$phred_average\t$count_major_mutations\t$count_major_indels\t$mutation_momentum\t$indels_momentum\t$major_mutation_momentum\t$major_indels_momentum\t$alignment_anisotropy\t$overhang_momentum\t$aligned_momentum\t$entropic_momentum\t$reads_unaligned\t$percent_reads_unaligned\t$percent_reads_aligned\t$reads_aligned\t$ASSEMBLY_STATUS\t$rpkm" >> "$output_file"

        # Print values to check if they are extracted correctly
        # echo "Extracted assembled_genome_acc: $nucleotide"
    done
done

# Indicate script completion
echo "Processing complete. Output written to $output_file"
