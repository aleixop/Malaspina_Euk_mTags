#!/usr/bin/env python3

__author__ = "Aleix Obiol"
__email__ = "obiol@icm.csic.es"

import sys
import argparse
import os.path

def get_args():
    """
    Parse and return arguments passed in.
    """
    parser = argparse.ArgumentParser(
        description = "Make OTU table out of a usearch filtered map by script `mitags_consensus_tax.py`")
    # Add arguments
    parser.add_argument("input",
                        help = "input filtered map")
    parser.add_argument("-i",
                        "--sample_identifier",
                        help = "label or character preceding sample name in each read, should be at the end of the field/column. Example: for a read like  `READ123;sample=Sample1`, you should write `-i 'sample='`",
                        required = True)
    parser.add_argument("-o",
                        "--output_file",
                        help = "path to output file where otu table should be written. It defaults to `otu_table.txt`",
                        required = False,
                        default = 'otu_table.txt')

    # Array for all arguments passed to script
    args = parser.parse_args()
    # Assign args to variables
    input = args.input
    sample_identifier = args.sample_identifier
    output_file = args.output_file
    # Return all variable values
    return input, sample_identifier, output_file

def get_samples(input, sample_identifier):
    """
    Input: usearch filtered map by script `mitags_consensus_tax.py`.
    Output: list of all sample names
    """
    samples = []
    with open(input, 'r') as uc_file:
        for line in uc_file:
            line = line.strip().split()
            sample = line[0].split(sample_identifier)[-1]
            if sample[-1] == ';': # patch to remove a ';' in my data
                sample = sample[:-1]
            if sample not in samples:
                samples.append(sample)
    return sorted(samples)

def make_otu_table(input, sample_identifier, samples):
    """
    Input: usearch filtered map by script `mitags_consensus_tax.py`.
    Output: dictionary with otu as keys, each containing another dictionary with samples as keys and number of reads as values.
    """
    dict_otu_table = {}
    with open(input, 'r') as uc_file:
        for line in uc_file:
            line = line.strip().split()
            otu = line[-1]
            sample = line[0].split(sample_identifier)[-1]
            if sample[-1] == ';':
                sample = sample[:-1] # patch to remove a ';' in my data
            if otu in dict_otu_table:
                dict_otu_table[otu][sample] += 1
            else:
                dict_otu_table[otu] = {}
                for sam in samples:
                    if sam == sample:
                        dict_otu_table[otu][sam] = 1
                    else:
                        dict_otu_table[otu][sam] = 0
    return dict_otu_table

def main():
    input, sample_identifier, output_file = get_args()
    if os.path.isfile(output_file): # check if file already exists and stop script if it does
        sys.exit('Stopping... The output file "'+ output_file + '" already exists. Remove it or specify another output file and try again.')
    samples = get_samples(input, sample_identifier)
    dict_otu_table = make_otu_table(input, sample_identifier, samples)
    with open(output_file, 'a') as file:
        print('OTUId', '\t'.join(samples), sep = '\t', file = file)
        for otu in dict_otu_table:
            print(otu, '\t'.join(str(x) for x in dict_otu_table[otu].values()), sep = '\t', file = file)

if __name__ == '__main__':
    main()
