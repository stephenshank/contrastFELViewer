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
            translated = Seq(''.join(self.sequence_data[i, :])).translate(gap='-')
            amino_acids.append(translated)
        self.sequence_data = np.array(amino_acids, dtype='<U1')


def remove_all_gap_columns(input, output_fasta, output_json):
    alignment = Alignment(input)
    indices = alignment.remove_all_gap_columns()
    count_gaps = np.sum(alignment.sequence_data == '-', axis=0)
    no_gap = np.argmin(count_gaps)
    assert count_gaps[no_gap] == 0
    header = alignment.headers[no_gap]
    alignment.write(output_fasta)
    
    result = {"indices": indices.tolist(), "no_gaps": header}
    with open(output_json, 'w') as output_file:
        json.dump(result, output_file)


def get_plot_data(hyphy_indices, added_alignment, reference, hyphy_data, output):
    with open(hyphy_indices) as hyphy_indices_file:
        hyphy_indices = json.load(hyphy_indices_file)
    no_gap_header = hyphy_indices["no_gaps"]
    hyphy_indices = np.array(hyphy_indices["indices"], dtype=np.int)
    added_alignment = Alignment(added_alignment)
    reference = SeqIO.to_dict(SeqIO.parse(reference, 'fasta'))
    with open(hyphy_data) as hyphy_data_file:
        hyphy_data = json.load(hyphy_data_file)
    output_dict = {}
    no_gap_seq = added_alignment.sequence(no_gap_header)
    no_gap_indices = np.arange(added_alignment.number_of_sites)[no_gap_seq != '-']

    headers = [row[0] for row in hyphy_data["MLE"]["headers"]]
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
        hyphy_output = np.zeros(added_alignment.number_of_sites)
        hyphy_output[no_gap_indices] = data[hyphy_indices]
        output_dict[category] = hyphy_output.tolist()

    hxb2 = added_alignment.sequence('HXB2_GP120')
    hxb2_map = np.arange(added_alignment.number_of_sites)[hxb2 != '-']
    hxb2_output = []
    with open('data/gp120_annotations.csv') as hxb2_file:
        reader = DictReader(hxb2_file)
        for row in reader:
            hxb2_index = int(row['fastaIndex'])
            if hxb2_index < len(hxb2_map):
                hxb2_output.append({
                    'annotation': row['protein feature'],
                    'index': int(hxb2_map[hxb2_index])+1
                })
    output_dict['hxb2'] = hxb2_output

    with open(output, 'w') as output_file:
        json.dump(output_dict, output_file)
    

def translate(input, output):
    alignment = Alignment(input)
    alignment.translate()
    alignment.write(output)


def bundle_json(input, output, patient_id):
    with open(input) as file:
        fasta = file.read()

    newick_path = 'data/%s_cFEL/%s.new' % (patient_id, patient_id)
    with open(newick_path) as file:
        newick = file.read()

    fel_json_path = 'data/%s_cFEL/%s_mappedIndices.json' % (patient_id, patient_id)
    with open(fel_json_path) as file:
        hyphy = json.load(file)

    with open('data/3jwo.pdb') as file:
        structure = file.read()

    with open(output, 'w') as file:
        output = {
            "fasta": fasta,
            "newick": newick,
            "hyphy": hyphy,
            "structure": structure
        }
        json.dump(output, file)

