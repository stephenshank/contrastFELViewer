from py import *


THAI_PATIENT_IDS = ['P01', 'P02', 'P13']
GUT_PATIENT_IDS = ['P19', 'P21', 'P29']
GUT_MODELS = ["Cells", "SLI"]

rule unpack_thai:
  input:
    tar="data/thai.tar.gz"
  output:
    expand("data/input/{patient_id}.fasta", patient_id=THAI_PATIENT_IDS),
    expand("data/input/{patient_id}.new", patient_id=THAI_PATIENT_IDS),
    expand("data/input/{patient_id}.fna.FEL.json", patient_id=THAI_PATIENT_IDS),
    "data/input/3jwo.pdb",
    "data/input/gp120_annotations.csv",
    "data/input/ref_struct.fasta"
  shell:
    "tar -xvzf data/thai.tar.gz -C data/input"

rule unpack_gut:
  input:
    tar="data/gut.tar.gz"
  output:
    expand("data/input/{patient_id}_Env_{model}_noBL.fna", patient_id=GUT_PATIENT_IDS, model=GUT_MODELS),
    expand("data/input/{patient_id}_Env_{model}_noBL.fna.FEL.json", patient_id=GUT_PATIENT_IDS, model=GUT_MODELS)
  shell:
    "tar -xvzf data/gut.tar.gz -C data/input"

rule translate:
  input:
    fasta="data/input/{patient_id}_Env_{model}_noBL.fasta"
  output:
    fasta="data/{patient_id}_Env_{model}_noBL/AA.fasta"
  run:
    translate(input.fasta, output.fasta)

rule remove_all_gap_columns:
  input:
    fasta=rules.translate.output.fasta
  output:
    fasta="data/{patient_id}_Env_{model}_noBL/noGaps.fasta",
    json="data/{patient_id}_Env_{model}_noBL/noGaps.json"
  run:
    remove_all_gap_columns(input.fasta, output.fasta, output.json)

rule added_alignment:
  input:
    reference="data/input/ref_struct.fasta",
    patient=rules.remove_all_gap_columns.output.fasta
  output:
    fasta="data/{patient_id}_Env_{model}_noBL/added.fasta"
  shell:
    "mafft --add {input.reference} {input.patient} > {output.fasta}"

rule get_plot_data:
  input:
    json=rules.remove_all_gap_columns.output.json,
    fasta=rules.added_alignment.output.fasta,
    reference="data/input/ref_struct.fasta",
    hyphy="data/input/{patient_id}_Env_{model}_noBL.fna.FEL.json"
  output:
    json="data/{patient_id}_Env_{model}_noBL/mappedIndices.json"
  run:
    get_plot_data(input.json, input.fasta, input.reference, input.hyphy, output.json)

rule json:
  input:
    fasta=rules.added_alignment.output.fasta,
    json=rules.get_plot_data.output.json,
    newick="data/input/{patient_id}_Env_{model}_noBL.new"
  output:
    json="data/{patient_id}_Env_{model}_noBL/dashboard.json"
  run:
    bundle_json(input.fasta, input.newick, input.json, output.json)

rule gut:
  input:
    expand(rules.json.output.json, patient_id=THAI_PATIENT_IDS, model=GUT_MODELS)

