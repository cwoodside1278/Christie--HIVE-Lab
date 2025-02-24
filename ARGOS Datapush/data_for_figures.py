import pandas as pd
import argparse


'''python data_for_figures.py -tsv1 /Users/christiewoodside/Desktop/ngsQC_HIVE3-Jan31_fixedinstrument.tsv -tsv2 /Users/christiewoodside/Desktop/assemblyQC_HIVE3.tsv -out /Users/christiewoodside/Desktop/test.tsv'''
def compute_gc_averages(tsv1, tsv2, output_file):
    # Read the TSV files
    df_ngs = pd.read_csv(tsv1, sep='\t')
    df_assembly = pd.read_csv(tsv2, sep='\t')
    
    # Ensure necessary columns exist
    required_ngs_cols = {'genome_assembly_id', 'organism_name', 'instrument', 'ngs_read_file_name', 'ngs_gc_content', 'avg_phred_score'}
    required_assembly_cols = {'genome_assembly_id', 'assembly_gc_content', 'lineage'}
    
    if not required_ngs_cols.issubset(df_ngs.columns):
        raise ValueError(f"The NGS TSV file must contain columns: {required_ngs_cols}")
    if not required_assembly_cols.issubset(df_assembly.columns):
        raise ValueError(f"The Assembly TSV file must contain columns: {required_assembly_cols}")
    
    # Filter only Illumina rows
    df_ngs = df_ngs[df_ngs['instrument'].str.contains('Illumina', case=False, na=False)]
    
    # Extract SRR ID from ngs_read_file_name (assuming format SRRxxxxx_1.fastq, _2.fastq, etc.)
    # df_ngs['SRR_id'] = df_ngs['ngs_read_file_name'].str.extract(r'(SRR\d+)')
    df_ngs['SRR_id'] = df_ngs['ngs_read_file_name'].str.extract(r'([ES]RR\d+)')

    
    # Convert percentages to numeric values
    df_ngs['ngs_gc_content'] = df_ngs['ngs_gc_content'].str.rstrip('%').astype(float)
    df_assembly['assembly_gc_content'] = df_assembly['assembly_gc_content'].str.rstrip('%').astype(float)

    # Extract the first word from the lineage column
    #df_assembly['Family'] = df_assembly['lineage'].str.split().str[0]
    # Extract the first word from the lineage column, removing any trailing semicolon
    df_assembly['Family'] = df_assembly['lineage'].str.split().str[0].str.rstrip(';')

    #print(df_ngs[df_ngs['genome_assembly_id'] == 'GCA_000865725.1'][['genome_assembly_id', 'ngs_read_file_name', 'SRR_id', 'ngs_gc_content', 'avg_phred_score']])

    # Compute averages per genome_assembly_id and SRR_id
    grouped_ngs = df_ngs.groupby(['genome_assembly_id', 'SRR_id']).agg(
        average_ngs_gc_content=('ngs_gc_content', 'mean'),
        phred_average=('avg_phred_score', 'mean')
    ).reset_index()
    #print(grouped_ngs[grouped_ngs['genome_assembly_id'] == 'GCA_000865725.1'])
    
    # Compute final NGS averages per genome_assembly_id
    final_ngs = grouped_ngs.groupby('genome_assembly_id', as_index=False).agg(
        average_ngs_gc_content=('average_ngs_gc_content', 'mean'),
        phred_average=('phred_average', 'mean')
    )
    
    # Compute the average assembly GC content per genome_assembly_id
    #assembly_avg = df_assembly.groupby('genome_assembly_id', as_index=False)['assembly_gc_content'].mean()
    assembly_avg = df_assembly.groupby('genome_assembly_id', as_index=False)[['assembly_gc_content', 'Family']].first()

    #print(df_ngs[df_ngs['organism_name'].isna()])

    merged_df = pd.merge(final_ngs, assembly_avg, on='genome_assembly_id', how='inner') #f you only want rows that had Illumina reads, you could change the merge type to inner join
    #merged_df = pd.merge(final_ngs, assembly_avg, on='genome_assembly_id', how='outer')
    #merged_df = pd.merge(df_ngs[['genome_assembly_id', 'organism_name']].drop_duplicates(), merged_df, on='genome_assembly_id', how='left')
    merged_df = pd.merge(df_ngs[['genome_assembly_id', 'organism_name']].drop_duplicates(), merged_df, on='genome_assembly_id', how='outer')

    # Rounding to one decimal place
    merged_df['average_ngs_gc_content'] = merged_df['average_ngs_gc_content'].round(1)
    merged_df['phred_average'] = merged_df['phred_average'].round(1)
    merged_df['assembly_gc_content'] = merged_df['assembly_gc_content'].round(1)


    # Reorder columns
    merged_df = merged_df[['organism_name', 'assembly_gc_content', 'phred_average', 'average_ngs_gc_content', 'Family', 'genome_assembly_id']]
    
    # Save output
    merged_df.to_csv(output_file, sep='\t', index=False)
    print(f"\n Output saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute average GC content per genome assembly ID.")
    parser.add_argument("-tsv1", required=True, help="Path to the NGS TSV")
    parser.add_argument("-tsv2", required=True, help="Path to the Assembly TSV file")
    parser.add_argument("-out", required=True, help="Path to save the output TSV file.")
    
    args = parser.parse_args()
    compute_gc_averages(args.tsv1, args.tsv2, args.out)
