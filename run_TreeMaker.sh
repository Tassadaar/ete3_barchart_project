#!/bin/bash

source activate treeProject

python TreeMaker.py \
  -t input_files/1st_scenario.tree \
  -n input_files/toy.fasta \
  -f fasta -m relative \
  -s fymink,garp

conda deactivate