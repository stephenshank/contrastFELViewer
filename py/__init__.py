import json

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
    alignment.write(output_fasta)
    
    result = {"indices": indices.tolist(), "no_gaps": no_gap}
    with open(output_json, 'w') as output_file:
        json.dump(indices.tolist(), output_file)


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

    fel_json_path = 'data/%s_cFEL/%s.fna.FEL.json' % (patient_id, patient_id)
    with open(fel_json_path) as file:
        hyphy = file.read()

    indices_path = "data/%s_cFEL/%s_noGaps.json" % (patient_id, patient_id)
    with open(indices_path) as file:
        indices = np.array(json.load(file))

    hxb2_annotations = pd.read_csv('data/gp120_annotations.csv')

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

