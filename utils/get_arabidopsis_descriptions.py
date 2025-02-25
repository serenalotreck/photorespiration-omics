"""
Get metadata associated with TAIR ID's for a list of genes.

Author: Serena G. Lotreck
"""
from os.path import abspath
import pandas as pd


def get_arabidopsis_descriptions(gene_list, gene2GO, metadata_path):
    """
    Get arabidopsis gene metadata.

    parameters:
        gene_list, list of str: genes to check
        gene2GO, df: has columns object_name, GO_term, and GO_ID
        metadata_path, str: path to the file containing Arabidopsis metadata
            (Araport11_functional_descriptions_20241231.txt)

    returns:
        gene_df, pandas df: genes with their associated data
    """
    # Read in the file
    metadata = pd.read_csv(abspath(metadata_path), sep='\t', header=0, encoding='Windows-1252')

    # Make a new column that can match the candidate genes
    metadata['name_base'] = metadata['name'].str.split('.').str[0]

    # Add GO terms to metadata
    go_df = gene2GO[['object_name', 'GO_term', 'GO_ID']].groupby('object_name').agg(list)
    metadata = pd.merge(metadata, go_df, how='outer', left_on='name_base', right_on='object_name')

    # Make gene_list into a Series
    gene_list = pd.Series(gene_list, name='name_base') # To match with metadata df

    # Subset by gene list; left join to preserve IDs for novel genes/chimeras
    gene_df = pd.merge(gene_list, metadata, how='left', on='name_base').set_index(['name_base', 'name'])

    return gene_df