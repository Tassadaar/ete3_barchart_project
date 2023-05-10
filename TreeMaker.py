"""
This program makes a simple unrooted tree from a newwick string, and render said tree to PNG image.
"""

import argparse
from ete3 import Tree, faces, TreeStyle
from FaceMaker import FaceMaker

# specify tags
parser = argparse.ArgumentParser(description='Tree making')

parser.add_argument('-t', '--tree', required=True)
parser.add_argument('-n', '--filename', required=True)
parser.add_argument('-f', '--format', required=True)

args = parser.parse_args()

newwick_tree = args.tree
file_location = args.filename
file_format = args.format

# making faces
faceMaker = FaceMaker(file_location, file_format)
face_dict = faceMaker.make_barchartface()

# tree "growing"
tree = Tree(newwick_tree)

# tree styling
tree_style = TreeStyle()
tree_style.show_scale = False  # do not show scale


# layout function
def layout_fn(node):
    if node.is_leaf():
        face = face_dict[node.name]
        faces.add_face_to_node(face=face, node=node, column=1, position="aligned")


# render to "file name"
tree.render("test.png", units="px", h=2000, w=2500, dpi=70,  tree_style=tree_style, layout=layout_fn)
