"""
Convert the paired read files in a directory to a samplesheet for input to
the nf-core/rnaseq pipeline.

Author: Serena G. Lotreck
"""
import argparse
from os import listdir, extsep
from os.path import abspath
import json
import pandas as pd
import warnings


def main(data_dir, semantic_names, skip_missing, out_loc):

    print('\nReading in data...')

    # Read in the semantic_names if provided
    if semantic_names != '':
        with open(semantic_names) as f:
            semantic_dict = json.load(f)

    # Set up dataframe structure
    samples = {'sample': [],
            'fastq_1': [],
            'fastq_2': [],
            'strandedness':[]}

    # Get the pairs of filenames
    print('\nPairing filenames...')
    all_fastq = []
    for f in listdir(data_dir):
        try:
            if f.split(extsep, 1)[1] == 'fastq.gz':
                all_fastq.append(f)
        except IndexError:
            continue
    all_fastq.sort()
    first_even = True
    for i, elt in enumerate(all_fastq):
        if (first_even and i % 2 == 0) or (not first_even and i % 2 != 0):
            # Check for the paired files
            if elt.split('_')[:3] != all_fastq[i+1].split('_')[:3]:
                if not skip_missing:
                    assert elt.split('_')[:3] != all_fastq[i+1].split('_')[:3], (
                        f'File {elt} is missing its paired counterpart')
                else:
                    warning = (f'File {elt} is missing its paired counterpart, '
                        'this sample will be skipped')
                    warnings.warn(warning)
                first_even = not first_even
                continue
            # Semantic name conversion if requested
            if semantic_names != '':
                samp_name = semantic_dict[elt.split('_')[0]]
            else:
                samp_name = elt.split('_')[0]
            samples['sample'].append(samp_name)
            samples['fastq_1'].append(elt)
            samples['fastq_2'].append(all_fastq[i+1])
            samples['strandedness'].append('auto')

    # Make dataframe
    print('\nMaking dataframe...')
    samplesheet_df = pd.DataFrame(samples)
    print(f'\nSnapshot of samplesheet:\n{samplesheet_df.head()}')

    # Save
    print(f'\nSaving as {out_loc}/samplesheet.csv...')
    samplesheet_df.to_csv(f'{out_loc}/samplesheet.csv', index=False)
    print('\Done!')

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Make samplesheet')

    parser.add_argument('data_dir', type=abspath,
            help='Directory containing RNAseq files to use')
    parser.add_argument('-semantic_names', type=str, default='',
            help='Optionally, provide a mapping between the sample '
            'names as present in data_dir and a semantic naming to '
            'use in the `sample` column of the samplesheet. '
            'Provide a json where keys are the string that appears before the '
            'first hyphen in the filename, and values are the semantic names')
    parser.add_argument('--skip_missing', action='store_true',
            help='Whether or not to skip files that are missing pairs')
    parser.add_argument('out_loc', type=abspath,
            help='Where to write output samplesheet')

    args = parser.parse_args()

    if args.semantic_names != '':
        args.semantic_names = abspath(args.semantic_names)

    main(args.data_dir, args.semantic_names, args.skip_missing, args.out_loc)
