import json

import pandas as pd
from Bio import SeqIO
import numpy as np


def remove_all_gap_columns(input, output_fasta, output_json):
    records = SeqIO.parse(input, 'fasta')
    nucleotides = []
    headers = []
    for record in records:
        nucleotides.append(list(str(record.seq)))
        headers.append(record.id)
    sequence_data = np.array(nucleotides, dtype='<U1')
    number_of_sequences = sequence_data.shape[0]
    number_of_sites = sequence_data.shape[1]
    all_gaps = np.sum(sequence_data == '-', axis=0) == number_of_sequences

    with open(output_fasta, 'w') as output_file:
        for i, header in enumerate(headers):
            output_file.write('>' + header + '\n')
            output_file.write(''.join(sequence_data[i, ~all_gaps]))
            output_file.write('\n')
    indices = np.arange(number_of_sites)[~all_gaps]

    count_gaps = np.sum(sequence_data[:, ~all_gaps] == '-', axis=1)
    no_gap = np.argmin(count_gaps)
    assert count_gaps[no_gap] == 0
    with open(output_json, 'w') as output_file:
        result = {"indices": indices.tolist(), "no_gaps": no_gap}
        json.dump(indices.tolist(), output_file)


def translate(input, output):
    records = list(SeqIO.parse(input, 'fasta'))
    for record in records:
        record.seq = record.seq.translate(gap='-')
    SeqIO.write(records, output, 'fasta')


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

