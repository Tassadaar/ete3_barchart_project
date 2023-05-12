"""
This program makes a simple unrooted tree from a newwick string, and render said tree to PNG image.
"""

import argparse
from ete3 import Tree, faces, TreeStyle
from FrequencyCalculator import FrequencyCalculator
from FaceMaker import FaceMaker

# specify tags
parser = argparse.ArgumentParser(description="Tree making")

parser.add_argument("-t", "--tree", required=True)
parser.add_argument("-n", "--filename", required=True)
parser.add_argument("-f", "--format", required=True)
parser.add_argument("-m", "--mode", type=str, default="normal")

args = parser.parse_args()

newwick_tree = args.tree
file_location = args.filename
file_format = args.format
mode = args.mode

# get taxa list
taxa_list = FrequencyCalculator(file_location, file_format).get_taxa_list()

# making faces
faceMaker = FaceMaker()
face_dict = {}
fymink_dict = {}
garp_dict = {}
others_dict = {}

for taxon in taxa_list:
    face_dict[taxon.name] = faceMaker.make_barchartface(taxon.freq_dict)
    fymink_dict[taxon.name] = faceMaker.make_barchartface(taxon.fymink_freq_dict)
    garp_dict[taxon.name] = faceMaker.make_barchartface(taxon.garp_freq_dict)
    others_dict[taxon.name] = faceMaker.make_barchartface(taxon.others_freq_dict)


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


try:
    if mode == "normal":
        # render to "file name"
        tree.render("test.png", units="px", h=2000, w=2500, dpi=70,  tree_style=tree_style, layout=layout_fn)
    elif mode == "special":
        print("Get pranked, sucker!")
    else:
        raise ValueError("Invalid tag value for mode, make sure to check the list of valid tags and check spelling!")
except ValueError as e:
    print(f"Error: {e}")
