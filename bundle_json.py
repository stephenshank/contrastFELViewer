import json


with open('../data/P01_cFEL/P01_AA.fasta') as file:
    fasta = file.read()

with open('../data/P01_cFEL/P01.new') as file:
    newick = file.read()

with open('../data/P01_cFEL/P01.fna.FEL.json') as file:
    hyphy = file.read()

with open('../data/3jwo.pdb') as file:
    structure = file.read()

with open('public/P01.json', 'w') as file:
    output = {
        "fasta": fasta,
        "newick": newick,
        "hyphy": hyphy,
        "structure": structure
    }
    json.dump(output, file)
