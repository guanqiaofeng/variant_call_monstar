#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
  Copyright (C) 2021,  Ontario Institute for Cancer Research

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  Authors:
    Guanqiao Feng
"""

import xml.etree.ElementTree as ET
import pandas as pd
import argparse

def create_vcf_header(root):
    # Extract attributes from variant-report element
    curation_version_id = root.get('curation-version-id')
    disease = root.get('disease')
    flowcell_analysis = root.get('flowcell-analysis')
    gender = root.get('gender')
    human_genome_assembly = root.get('human-genome-assembly')
    specimen_id = root.get('specimen-id')
    study = root.get('study')
    test_request = root.get('test-request')
    test_type = root.get('test-type')

    # Extract attributes from sample element (assuming there is only one sample)
    sample = root.find('.//sample')
    name = sample.get('name')
    bait_set = sample.get('bait-set')

    # Construct the headers
    headers = [
        '##fileformat=VCFv4.3',
        f'##curation-version-id=\"{curation_version_id}\"',
        f'##disease="{disease}"',
        f'##flowcell-analysis="{flowcell_analysis}"',
        f'##gender="{gender}"',
        f'##human-genome-assembly="{human_genome_assembly}"',
        f'##specimen-id="{specimen_id}"',
        f'##study="{study}"',
        f'##test-request="{test_request}"',
        f'##test-type="{test_type}"',
        f'##name="{name}"',
        f'##bait-set="{bait_set}"',
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">',
        '##INFO=<ID=SOMATIC-OR-GERMLINE,Type=String,Description="somatic or germline">',
        '##INFO=<ID=ZYGOSITY,Type=String,Description="sample zygosity">'
    ]

    return headers

def extract_data_from_short_variant(short_variant):
    return {
        # Extract useful information from the xml file
        'AF': short_variant.get('allele-fraction'),
        'ALT': short_variant.get('alternate-sequence'),
        '#CHROM': short_variant.get('chromosome'),
        'DP': short_variant.get('depth'),
        'POS': short_variant.get('genomic-start'),
        'REF': short_variant.get('reference-sequence'),
        'SOMATIC-or-GERMLINE': short_variant.get('germline-status'),
        'ZYGOSITY': short_variant.get('tumor-zygosity')
    }

def chr_pos_sort(df):
    # Define a custom sort order for chromosome numbers
    chrom_order = {'chr' + str(i): i for i in range(1, 23)}
    chrom_order['chrX'] = 23
    chrom_order['chrY'] = 24

    # Apply the custom sort order
    df['#CHROM'] = df['#CHROM'].map(chrom_order)
    df['POS'] = df['POS'].astype(int)

    # Sort the DataFrame first by '#CHROM' and then by 'POS'
    sorted_df = df.sort_values(by=['#CHROM', 'POS'])

    # Convert the chromosome numbers back to their original format
    sorted_df['#CHROM'] = 'chr' + sorted_df['#CHROM'].astype(str)

    return sorted_df

# Process dataframe by adding "ID", "QUAL", "FILTER" with default value "." and "INFO" with "DP=#;AF=#"
# where DP means depth and AF means allele frequency
def process_dataframe(df):

    df['ID'] = '.'
    df['QUAL'] = '.'
    df['FILTER'] = '.'
    df['INFO'] = 'DP=' + df['DP'].astype(str) + ';AF=' + df['AF'].astype(str) + ';SOMATIC-OR-GERMLINE=' + df['SOMATIC-or-GERMLINE'] + ';ZYGOSITY=' + df['ZYGOSITY']
    df.drop(['DP', 'AF'], axis=1, inplace=True)
    desired_order = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    df = df[desired_order]
    sorted_df = chr_pos_sort(df.copy())

    return sorted_df

def main():
    # Create the parser
    parser = argparse.ArgumentParser(description="Process some files.")

    # Add input file argument
    parser.add_argument('-i', '--input', type=str, required=True, help="Input file name")

    # Add output file argument
    parser.add_argument('-o', '--output', type=str, required=True, help="Output file name")

    # Parse the arguments
    args = parser.parse_args()

    # Use the parsed input and output file names
    input_file_name = args.input
    output_file_name = args.output

    tree = ET.parse(input_file_name)
    root = tree.getroot()

    # Create VCF headers
    vcf_headers = create_vcf_header(root)

    # Extract information from xml
    data = []

    for short_variant in root.findall('.//short-variant'):
        short_variant_data = extract_data_from_short_variant(short_variant)
        data.append(short_variant_data)
    df = pd.DataFrame(data)

    # Process the DataFrame
    df_processed = process_dataframe(df)

    # Write headers and data to VCF file
    with open(output_file_name, 'w') as f:
        for header_line in vcf_headers:
            f.write(f"{header_line}\n")
        df_processed.to_csv(f, sep='\t', index=False)

if __name__ == "__main__":
    main()
