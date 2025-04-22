"""
Get metadata associated with TAIR ID's for a list of genes.

Author: Serena G. Lotreck
"""
from os.path import abspath
import pandas as pd


def get_arabidopsis_descriptions(gene_list, metadata_paths, gene2GO=None):
    """
    Get arabidopsis gene metadata.

    parameters:
        gene_list, list of str: genes to check
        metadata_paths, dict: keys are annotation names, values are paths to
            a file containing Arabidopsis metadata. Assumes all dfs have the
            same column names
        gene2GO, df: has columns object_name, GO_term, and GO_ID, provide if GO
            terms should be included in output

    returns:
        gene_df, pandas df: genes with their associated data
    """
    # Process metadata
    metadata_dfs = {}
    for ann, metadata_path in metadata_paths.items():
        # Read in the file
        metadata = pd.read_csv(abspath(metadata_path), sep='\t', header=0, encoding='Windows-1252')
        metadata_dfs[ann] = metadata
    # If more than one, merge
    if len(metadata_dfs) == 2:
        all_metadata = pd.merge(
        list(metadata_dfs.values())[0], list(metadata_dfs.values())[1],
        left_on='name', right_on='name',
        suffixes=(f'_{list(metadata_dfs.keys())[0]}', f'_{list(metadata_dfs.keys())[1]}'))
    else:
        all_metadata = metadata_dfs.values()[0]
    # Make a new column that can match the candidate genes
    all_metadata['gene_id'] = all_metadata['name'].str.split('.').str[0]

    # Add GO terms to metadata
    if gene2GO is not None:
        go_df = gene2GO[['object_name', 'GO_term', 'GO_ID']].groupby('object_name').agg(list)
        all_metadata = pd.merge(all_metadata, go_df, how='outer', left_on='gene_id', right_on='object_name')

    # Make gene_list into a Series
    gene_list = pd.Series(gene_list, name='gene_id') # To match with metadata df

    # Subset by gene list; left join to preserve IDs for novel genes/chimeras
    gene_df = pd.merge(gene_list, all_metadata, how='left', on='gene_id').set_index(['gene_id', 'name'])

    return gene_df