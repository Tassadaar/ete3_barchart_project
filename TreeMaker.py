"""
This program makes a simple unrooted tree from a newwick string, and show said tree in GUI.
"""

import argparse
from ete3 import Tree
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
tree = Tree(newwick_tree)

# add faces to the tree
for node in tree.traverse():
    if node.is_leaf():
        face = face_dict[node.name]
        node.add_face(face, column=1, position="branch-right")

# render to "file name"
tree.render("test.png", units="px", h=1920, w=1080)
