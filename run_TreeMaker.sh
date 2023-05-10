#!/bin/bash

source activate treeProject

python TreeMaker.py -t input_files/nematode_CH.tree -n input_files/nematode.fasta -f fasta

conda deactivate