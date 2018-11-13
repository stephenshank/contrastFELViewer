from py import *


rule translate:
  input:
    fasta="data/input/{patient_id}.fasta"
  output:
    fasta="data/{patient_id}/AA.fasta"
  run:
    translate(input.fasta, output.fasta)

rule remove_all_gap_columns:
  input:
    fasta=rules.translate.output.fasta
  output:
    fasta="data/{patient_id}/noGaps.fasta",
    json="data/{patient_id}/noGaps.json"
  run:
    remove_all_gap_columns(input.fasta, output.fasta, output.json)

rule added_alignment:
  input:
    reference="data/input/ref_struct.fasta",
    patient=rules.remove_all_gap_columns.output.fasta
  output:
    fasta="data/{patient_id}/added.fasta"
  shell:
    "mafft --add {input.reference} {input.patient} > {output.fasta}"

rule get_plot_data:
  input:
    json=rules.remove_all_gap_columns.output.json,
    fasta=rules.added_alignment.output.fasta,
    reference="data/input/ref_struct.fasta",
    hyphy="data/input/{patient_id}.fna.FEL.json"
  output:
    json="data/{patient_id}/mappedIndices.json"
  run:
    get_plot_data(input.json, input.fasta, input.reference, input.hyphy, output.json)

rule json:
  input:
    fasta=rules.added_alignment.output.fasta,
    json=rules.get_plot_data.output.json
  output:
    json="data/{patient_id}/dashboard.json"
  run:
    bundle_json(input.fasta, output.json, wildcards.patient_id)

rule all:
  input:
    expand(rules.json.output.json, patient_id=['P01', 'P02', 'P13'])

