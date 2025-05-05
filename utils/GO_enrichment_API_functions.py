"""
Function sfor making and processing calls to the PANTHER API GO tools.

Author: Serena G. Lotreck
"""
import requests
import pandas as pd
import numpy as np


def getPANTHER(query, checkPageInfo=True):
    """
    Get GO terms from a mmultiple-page search result.
    """
    # Check number of pages
    r = requests.get(query, headers={ "Accept" : "application/json"})
    if not r.ok:
      r.raise_for_status()
      sys.exit()
    responseBody = r.json()
    overall_jsons = responseBody['results']

    return overall_jsons

def processGOenrichments(enrichments, data, conditions_semantic):
    """
    Make a summary dataframe with GO terms from each aspect that are enriched
    in each group, and print a summary of the numbers.
    """
    processed_go = {
        'group': [],
        'aspect': [],
        'term': [],
        'GOid': [],
        'p_value_fdr': [],
        'associated_gene_IDs': [],
        'plus_minus': []
    }
    for group, res in enrichments.items():
        for aspect in res:
            for termDict in res[aspect]['result']:
                if termDict['fdr'] < 0.05:
                    processed_go['group'].append(group)
                    processed_go['aspect'].append(aspect)
                    processed_go['term'].append(termDict['term']['label'])
                    try:
                        processed_go['GOid'].append(termDict['term']['id'])
                    except KeyError:
                        processed_go['GOid'].append(np.nan)
                    processed_go['p_value_fdr'].append(termDict['fdr'])
                    try:
                        processed_go['associated_gene_IDs'].append(termDict['input_list']['mapped_ids'])
                        processed_go['plus_minus'].append('+')
                    except KeyError:
                        processed_go['associated_gene_IDs'].append(np.nan)
                        processed_go['plus_minus'].append('-')
                        

    go_results = pd.DataFrame(processed_go)
    go_results = go_results.set_index(['group', 'aspect']).sort_index()

    for group in data:
        print(f'\nFor comparison group {conditions_semantic[group]}, there are...')
        for aspect in ['biological_process', 'molecular_function', 'cellular_component']:
            try:
                print(f'{len(go_results.loc[(group, aspect)])} {aspect} GO terms')
            except KeyError:
                print(f'0 {aspect} GO terms')
        print(f'... enriched, for a total of {len(data[group])} genes in the comparison.\n')

    return go_results


