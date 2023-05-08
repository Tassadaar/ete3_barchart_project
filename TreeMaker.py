"""
This program makes a simple unrooted tree from a newwick string, and show said tree in GUI.
"""

import argparse
from ete3 import PhyloTree
from Bio import Phylo
from FaceMaker import FaceMaker

# specify tags
#parser = argparse.ArgumentParser(description='Tree making')
#parser.add_argument('-n', '--filename', required=True)
#parser.add_argument('-f', '--format', required=True)
#args = parser.parse_args()

# variable initiations
newwick_tree = "input_files/nematode_CH.tree"
file_location = "input_files/nematode.fasta" #args.filename
file_format = "fasta" #args.format
face_dict = FaceMaker(file_location, file_format).make_face()

# tree "growing"
tree = PhyloTree(newwick_tree)

# print(unrooted_tree)
for node in tree.traverse():
    if node.name != "":
        face = face_dict[node.name]
        node.add_face(face, column=0, position="branch-right")

Phylo.draw(tree)
