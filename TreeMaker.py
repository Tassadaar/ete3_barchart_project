"""
This program makes a simple unrooted tree from a newwick string, and render said tree to PNG image.
"""

import argparse
from ete3 import Tree, faces, TreeStyle, BarChartFace
from FrequencyCalculator import FrequencyCalculator

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
taxa_list = FrequencyCalculator(file_location, file_format).get_taxa_list(mode)

# tree "growing"
tree = Tree(newwick_tree)

# tree styling
tree_style = TreeStyle()
tree_style.show_scale = False  # do not show scale


# layout function
def layout_fn(node):
    if node.is_leaf():
        taxon = taxa_list[node.name]
        face = BarChartFace(
            values=list(taxon.freq_dict.values()),
            labels=list(taxon.freq_dict.keys()),
            colors= ["blue" for key in taxon.freq_dict.keys()],
            max_value = 0.2,
            width=100,
            height=50
        )
        faces.add_face_to_node(face=face, node=node, column=1, position="aligned")


def layout_sp(node):
    if node.is_leaf(): # && node.name = "HONEYBEE":
        taxon = taxa_list[node.name]
        fymink_face = BarChartFace(
            values=list(taxon.fymink_freq_dict.values()),
            labels=list(taxon.fymink_freq_dict.keys()),
            colors= ["blue" for key in taxon.fymink_freq_dict.keys()],
            max_value = 0.2,
            width=100,
            height=50
        )
        print(taxon.freq_dict.values())
        #garp_face = garp_dict[node.name]
        #others_face = others_dict[node.name]
        #print(fymink_face.values)
        faces.add_face_to_node(face=fymink_face, node=node, column=1, position="aligned")
        #faces.add_face_to_node(face=garp_face, node=node, column=2, position="aligned")
        #faces.add_face_to_node(face=others_face, node=node, column=3, position="aligned")


try:
    if mode == "normal":
        # render to "file name"
        tree.render("test.png", units="px", h=2000, w=2500, dpi=70,  tree_style=tree_style, layout=layout_fn)
    elif mode == "special":
        tree.render("test.png", units="px", h=2000, w=2500, dpi=70,  tree_style=tree_style, layout=layout_sp)
    else:
        raise ValueError("Invalid tag value for mode, make sure to check the list of valid tags and check spelling!")
except ValueError as e:
    print(f"Error: {e}")
