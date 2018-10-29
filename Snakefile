from py import *


rule remove_all_gap_columns:
  input:
    fasta="data/{patient_id}_cFEL/{patient_id}.fasta"
  output:
    fasta="data/{patient_id}_cFEL/{patient_id}_noGaps.fasta",
    json="data/{patient_id}_cFEL/{patient_id}_noGaps.json"
  run:
    remove_all_gap_columns(input.fasta, output.fasta, output.json)

rule translate:
  input:
    fasta=rules.remove_all_gap_columns.output.fasta
  output:
    fasta="data/{patient_id}_cFEL/{patient_id}_AA.fasta"
  run:
    translate(input.fasta, output.fasta)

rule alignment:
  input:
    reference="data/ref_struct.fasta",
    patient=rules.translate.output.fasta
  output:
    fasta="data/{patient_id}_cFEL/{patient_id}_displayed.fasta"
  shell:
    "mafft --add {input.reference} {input.patient} > {output.fasta}"

rule json:
  input:
    fasta=rules.alignment.output.fasta,
    json=rules.remove_all_gap_columns.output.json
  output:
    json="data/{patient_id}_cFEL/{patient_id}.json"
  run:
    bundle_json(input.fasta, output.json, wildcards.patient_id)

rule all:
  input:
    expand(rules.json.output.json, patient_id=['P01', 'P02', 'P13'])

