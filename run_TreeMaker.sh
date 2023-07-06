#!/bin/bash

source activate treeProject

python TreeMaker.py \
  -t Martijn_et_al_2019/alphaproteobacteria_untreated.aln.treefile \
  -n Martijn_et_al_2019/alphaproteobacteria_untreated.aln \
  -f fasta -m relative \
  -s fymink,garp \
  -g Dechloromonas_aromatica_RCB,Pseudomonas_aeruginosa_PA7

conda deactivate