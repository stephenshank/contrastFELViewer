import json
from csv import DictReader

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np


class Alignment:
    def __init__(self, path):
        records = SeqIO.parse(path, 'fasta')
        nucleotides = []
        self.headers = []
        for record in records:
            nucleotides.append(list(str(record.seq)))
            self.headers.append(record.id)
        self.sequence_data = np.array(nucleotides, dtype='<U1')
        self.number_of_sequences = self.sequence_data.shape[0]
        self.number_of_sites = self.sequence_data.shape[1]

    def remove_all_gap_columns(self):
        number_of_gaps = np.sum(self.sequence_data == '-', axis=0)
        all_gaps = number_of_gaps == self.number_of_sequences
        indices = np.arange(self.number_of_sites)[~all_gaps]
        print('removing ', np.arange(self.number_of_sites)[all_gaps])
        self.sequence_data = self.sequence_data[:, ~all_gaps]
        return indices

    def write(self, path):
        with open(path, 'w') as output_file:
            for i, header in enumerate(self.headers):
                output_file.write('>' + header + '\n')
                output_file.write(''.join(self.sequence_data[i, :]))
                output_file.write('\n')

    def sequence(self, header):
        index = self.headers.index(header)
        return self.sequence_data[index, :]

    def translate(self):
        amino_acids = []
        for i in range(self.number_of_sequences):
            nucleotides = np.array(self.sequence_data[i, :], dtype='<U1')
            for i in range(0, len(nucleotides), 3):
                codon = str(''.join(nucleotides[i:i+3]))
                if codon.count('-') < 3:
                    codon = codon.replace('-', 'N')
                    nucleotides[i] = codon[0]
                    nucleotides[i+1] = codon[1]
                    nucleotides[i+2] = codon[2]
            translated = Seq(''.join(nucleotides)).translate(gap='-')
            amino_acids.append(translated)
        self.sequence_data = np.array(amino_acids, dtype='<U1')


def remove_all_gap_columns(input, output_fasta, output_json):
    alignment = Alignment(input)
    indices = alignment.remove_all_gap_columns()
    count_gaps = np.sum(alignment.sequence_data == '-', axis=1)
    no_gap = np.argmin(count_gaps)
    assert count_gaps[no_gap] == 0
    header = alignment.headers[no_gap]
    alignment.write(output_fasta)
    
    result = {"indices": indices.tolist(), "no_gaps": header}
    with open(output_json, 'w') as output_file:
        json.dump(result, output_file, indent=4)


def get_plot_data(hyphy_indices, added_alignment, reference, hyphy_data_path, output):
    with open(hyphy_indices) as hyphy_indices_file:
        hyphy_indices = json.load(hyphy_indices_file)
    no_gap_header = hyphy_indices["no_gaps"]
    hyphy_indices = np.array(hyphy_indices["indices"], dtype=np.int)
    added_alignment = Alignment(added_alignment)
    reference = SeqIO.to_dict(SeqIO.parse(reference, 'fasta'))
    with open(hyphy_data_path) as hyphy_data_file:
        hyphy_data = json.load(hyphy_data_file)
    output_dict = {}
    no_gap_seq = added_alignment.sequence(no_gap_header)
    no_gap_indices = np.arange(added_alignment.number_of_sites)[no_gap_seq != '-']

    headers = [row[0] for row in hyphy_data["MLE"]["headers"]]
    if 'Cells' in hyphy_data_path:
        categories = [
            "alpha",
            "beta (LI)",
            "beta (JI)",
            "beta (background)",
            "beta (D)"
        ]
    else:
        categories = [
            "alpha",
            "beta (T_cells)",
            "beta (Monocytes)",
            "beta (background)",
            "beta (Plasma)"
        ]

    for category in categories:
        index = headers.index(category)
        data = np.array([row[index] for row in hyphy_data["MLE"]["content"]["0"]])
        percentile = np.percentile(data, 95)
        data[data > percentile] = percentile
        hyphy_output = np.zeros(added_alignment.number_of_sites)
        hyphy_output[no_gap_indices] = data[hyphy_indices]
        output_dict[category] = hyphy_output.tolist()

    hxb2 = added_alignment.sequence('HXB2_GP120')
    hxb2_map = np.arange(added_alignment.number_of_sites)[hxb2 != '-']
    hxb2_output = []
    with open('data/input/gp120_annotations.csv') as hxb2_file:
        reader = DictReader(hxb2_file)
        for row in reader:
            hxb2_index = int(row['fastaIndex'])
            if hxb2_index < len(hxb2_map):
                hxb2_output.append({
                    'annotation': row['protein feature'],
                    'index': int(hxb2_map[hxb2_index])+1
                })
    output_dict['hxb2'] = hxb2_output
    pairs = [('LI', 'JI'), ('LI', 'D'), ('JI', 'D')]
    pdb = added_alignment.sequence('3JWO')
    pdb_ungapped = pdb[pdb != '-']
    pdb_map = np.arange(added_alignment.number_of_sites)[pdb != '-']
    invert = lambda index, indicial_map: int(np.argmax(indicial_map == index))
    all_categories = ['overall'] + pairs
    all_selected_sites = []
    for category in all_categories:
        all_sites = enumerate(hyphy_data['MLE']['content']['0'])
        if isinstance(category, tuple):
            item1 = category[0]
            item2 = category[1]
            header = 'P-value for %s vs %s' % (item1, item2)
            if not header in headers:
                header = 'P-value for %s vs %s' % (item2, item1)
            output_key = '%s vs %s' % (item1, item2)
        else:
            header = 'P-value (overall)'
            output_key = 'Overall'
        index = headers.index(header)
        selected_sites = [
            {'original': i, 'pvalue': val[index], 'category': output_key.replace('_', ' ')} for i, val in all_sites if val[index] < .1
        ]

        for selected_site in selected_sites:
            print('original index', selected_site['original'])
            ungapped_index = invert(selected_site['original'], hyphy_indices)
            print('ungapped index', ungapped_index)
            added_index = no_gap_indices[ungapped_index]
            print('added index', added_index)
            pdb_index = invert(added_index, pdb_map)
            print('pdb index', pdb_index)
            pdb_character = pdb_ungapped[pdb_index]
            pdb_added_character = pdb[added_index]
            if pdb_added_character != '-':
                print('added:', added_index, pdb_added_character, ', pdb:', pdb_index, pdb_character)
                assert pdb_character == pdb_added_character
                selected_site['pdb'] = int(pdb_index+1)
                selected_site['pdbchar'] = pdb_character
            selected_site['index'] = int(added_index+1)
        all_selected_sites += selected_sites
        
        output_dict['siteAnnotations'] = sorted(
            hxb2_output + all_selected_sites,
            key=lambda x: x['index']
        )
    with open(output, 'w') as output_file:
        json.dump(output_dict, output_file, indent=4)


def translate(input, output):
    alignment = Alignment(input)
    alignment.translate()
    alignment.write(output)


def bundle_json(fasta_input, newick_input, json_input, json_output):
    with open(fasta_input) as file:
        fasta = file.read()

    with open(newick_input) as file:
        newick = file.read()

    with open(json_input) as file:
        hyphy = json.load(file)

    with open('data/input/3jwo.pdb') as file:
        structure = file.read()

    with open(json_output, 'w') as file:
        output = {
            "fasta": fasta,
            "newick": newick,
            "hyphy": hyphy,
            "structure": structure
        }
        json.dump(output, file, indent=4)

