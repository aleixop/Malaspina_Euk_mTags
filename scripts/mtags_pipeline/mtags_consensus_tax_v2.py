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
        description = "Make consensus taxonomy out of a usearch tophits map")
    # Add arguments
    parser.add_argument("input",
                        help = "input file in usearch's UC format.")
    parser.add_argument("-t",
                        "--tax_separator",
                        help = "character separating taxonomic levels.",
                        required = True)
    parser.add_argument("-s",
                        "--tax_sense",
                        choices = ['asc', 'desc'],
                        help = "sense of taxonomic levels in your database. 'asc' for lower to higher levels (e.g. ID_Diatomea_Stramenopiles_SAR_Eukaryota), 'desc' for higher to lower levels (e.g. Eukaryota_SAR_Stramenopiles_Diatomea_ID).",
                        required = True)
    parser.add_argument("-p",
                        "--pair_separator",
                        help = "pair (forward & reverse) character separator. Use this argument to remove redundancies from your dataset (i.e. reads that are represented for both forward and reverse pairs).",
                        required = False,
                        default = None)
    parser.add_argument("-o",
                        "--output_file",
                        help = "path to output file where filtered map should be written. It defaults to `filtered_map.uc`",
                        required = False,
                        default = 'filtered_map.uc')
    # Array for all arguments passed to script
    args = parser.parse_args()
    # Assign args to variables
    input = args.input
    tax_separator = args.tax_separator
    tax_sense = args.tax_sense
    pair_separator = args.pair_separator
    outfile = args.output_file
    # Return all variable values
    return input, tax_separator, tax_sense, pair_separator, outfile

def read_uc(input, tax_separator, tax_sense):
    """
    Read uc file and return a dict with read id as key and tax hits as value.
    Input: usearch's tophits uc file.
    Output: dict with read id as key and taxonomy hits as values.
    """
    with open(input, 'r') as uc_file:
        dict_taxs = {}
        for line in uc_file:
            line = line.strip().split()
            is_hit = line[0] == 'H' # check if line is for a hit or for a no hit ('H' vs 'N', respectively)
            if is_hit:
                read_id = line[8] # take read id, located in 9th column of the file
                if tax_sense == 'asc':
                    taxonomy = line[9].split(tax_separator) # take taxonomy column and split it
                elif tax_sense == 'desc':
                    taxonomy = line[9].split(tax_separator)[::-1] # take taxonomy and reverse order
                try:
                    dict_taxs[read_id]['hits'] += 1 # sum hits for each sequence
                    dict_taxs[read_id]['taxonomy'].append(taxonomy) # add taxonomy to taxonomy dict
                except KeyError: # fires when a read_id is read for the first time
                    percentage_identity = line[3] # take percentage_identity to the database
                    cigar_alignment = line[7]
                    dict_taxs[read_id] = {'hits': 1,
                                          'taxonomy': [taxonomy],
                                          'perc_id': percentage_identity,
                                          'alignment': cigar_alignment}
    return dict_taxs

def make_consensus_taxonomy(dict_taxs):
    """
    Read dict_taxs and make consensus taxonomy from all hits.
    Input: dict_taxs.
    Output: dict_consensus.
    """
    dict_consensus = {}
    for read_id in dict_taxs:
        hits = dict_taxs[read_id]['hits']
        taxonomy = dict_taxs[read_id]['taxonomy']
        first_taxonomy = taxonomy[0] # take first taxonomy list of the list, will be used to fill tax fields
        length_taxonomy = len(first_taxonomy)
        if hits == 1: # only one taxonomy, no problem here
            dict_consensus[read_id] = first_taxonomy
        elif hits > 1: # more than one hit, a consensus has to be reached
            dict_consensus[read_id] = [] # open empty list, it will be filled after
            check_list = [] # list to store the output of tax comparisons
            for taxs in zip(*taxonomy): # create an iterable that gives all entries for each position of the taxonomy lists
                check = len(set(taxs)) == 1 # if the iterable is length 1, all taxs agree. If it's larger, there's not an agreement
                check_list.append(check) # append comparison results to previosuly opened list
            try:
                first_true = check_list.index(True) # look for first position where there is a True
                if first_true != 0:
                    dict_consensus[read_id].extend(first_true * ['NA']) # add as many NA's as Falses are before first True
                    dict_consensus[read_id].extend(first_taxonomy[first_true:]) # from first True on, add the taxonomy
                else: # if first true is 0, taxs should be exactly the same
                    dict_consensus[read_id] = first_taxonomy
            except ValueError: # it fires when there's not a single True, so there's no consensus. The read is totally ambiguous
                dict_consensus[read_id].extend(length_taxonomy*['NA'])
    total_reads = len(dict_consensus)
    return dict_consensus, length_taxonomy, total_reads

def remove_redundancies(dict_consensus, pair_separator):
    """
    Read dict_consensus and merge reads represented twice by forward & reverse pairs into one.
    Input: dict_consensus.
    Output: dict_no_redundancies.
    """
    dict_pairs = {}
    dict_no_redundancies = {}
    redundant_reads = 0
    for read_id in sorted(dict_consensus.keys()):
        taxonomy = dict_consensus[read_id]
        read_no_pair = read_id.split(pair_separator)[0] # remove pair info
        rest_of_read = read_id[(read_id.index(pair_separator)+1):] # take everything rightwards to first pair separator
        try: # if passes this, means there are 2 pairs representing the same sequence
            dict_pairs[read_no_pair]['taxonomy'].append(taxonomy) # add taxonomy of second pair
            redundant_reads += 1
            new_taxonomy = []
            previous_read_id = dict_pairs[read_no_pair]['original_read']
            del dict_no_redundancies[previous_read_id]
            for taxs in zip(*dict_pairs[read_no_pair]['taxonomy']):
                taxs = list(taxs)
                tax1 = taxs[0]
                tax2 = taxs[1]
                if tax1 == tax2:
                    new_taxonomy += [tax1]
                elif tax1 != tax2:
                    if 'NA' in taxs:
                        taxs.remove('NA')
                        new_taxonomy += taxs
                    else:
                        new_taxonomy += ['new_NA']
            if 'new_NA' in new_taxonomy:
                new_taxonomy = new_taxonomy[::-1] # reverse the order of taxonomy to check for 'new_NA' in the higher tax rank
                na_position = new_taxonomy.index('new_NA')
                new_taxonomy = new_taxonomy[:na_position] + (len(new_taxonomy)-na_position)*['NA']
                new_taxonomy = new_taxonomy[::-1]
                new_read_id = read_no_pair + pair_separator + 'X' + rest_of_read[1:] # remove number of pair in the original read
                dict_no_redundancies[new_read_id] = new_taxonomy
            else:
                dict_no_redundancies[read_id] = new_taxonomy
        except KeyError:
            dict_pairs[read_no_pair] = {'original_read' : read_id,
                                        'taxonomy' : [taxonomy]}
            dict_no_redundancies[read_id] = taxonomy
    final_reads = len(dict_no_redundancies)
    return dict_no_redundancies, redundant_reads, final_reads

def main():
    # Take arguments passed in
    input, tax_separator, tax_sense, pair_separator, outfile = get_args()
    if os.path.isfile(outfile): # check if file already exists and stop script if it does
        sys.exit('Stopping... The output file "'+ outfile + '" already exists. Remove it or specify another output file and try again.')
    print("# Reading " + '"' + input + '"' + "...")
    print("# The specified tax separator is " + '"' + tax_separator + '"')
    if tax_sense == 'asc':
        tax_sense_long = 'ascending'
    elif tax_sense == 'desc':
        tax_sense_long = 'descending'
    print("# The specified tax sense is " + '"' + tax_sense_long + '"')
    dict_taxs = read_uc(input, tax_separator, tax_sense)
    print("# Making consensus taxonomy...")
    dict_consensus, length_taxonomy, total_reads = make_consensus_taxonomy(dict_taxs)
    dict_tax_levels = {}
    for i in range(0,length_taxonomy+1):
        dict_tax_levels[i] = 0
    if pair_separator:
        print("# Removing redundancies...")
        dict_no_redundancies, redundant_reads, final_reads = remove_redundancies(dict_consensus, pair_separator)
        for read_id in dict_no_redundancies:
            if tax_sense == 'asc':
                taxonomy = dict_no_redundancies[read_id]
            elif tax_sense == 'desc':
                taxonomy = dict_no_redundancies[read_id][::-1] # reverse order back to original
            NA_count = taxonomy.count('NA')
            dict_tax_levels[NA_count] += 1
            try:
                perc_id = dict_taxs[read_id]['perc_id']
                cigar_alignment = dict_taxs[read_id]['alignment']
            except KeyError:
                perc_id = '-'
                cigar_alignment = '-'
            with open(outfile, 'a') as file:
                print(read_id, perc_id, cigar_alignment, tax_separator.join(taxonomy), sep = '\t', file = file)
    else:
        for read_id in dict_consensus:
            if tax_sense == 'asc':
                taxonomy = dict_consensus[read_id]
            elif tax_sense == 'desc':
                taxonomy = dict_consensus[read_id][::-1] # reverse order back to original
            NA_count = taxonomy.count('NA')
            dict_tax_levels[NA_count] += 1
            perc_id = dict_taxs[read_id]['perc_id']
            cigar_alignment = dict_taxs[read_id]['alignment']
            with open(outfile, 'a') as file:
                print(read_id, perc_id, cigar_alignment, tax_separator.join(taxonomy), sep = '\t', file = file)
    print("# A filtered map was written in " + '"' + outfile + '"')
    print("# Here's a summary of your data:")
    print("## Original reads: " + str(total_reads))
    if pair_separator:
        print("## Reads represented by forward & reverse pairs: " + str(redundant_reads))
        print("## Reads in final filtered map: " + str(final_reads))
    else:
        print("## Reads in final filtered map: " + str(total_reads))
    print("### Number of reads for each taxonomic level definition:")
    for NA_count in sorted(dict_tax_levels.keys()):
        print("#### Level " + str(NA_count) + ': ' + str(dict_tax_levels[NA_count]))
    print("# All done!")

if __name__ == '__main__':
    main()
