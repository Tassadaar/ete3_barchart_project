"""
This program makes a simple unrooted tree from a newwick string, and render said tree to PNG image.
"""

import argparse
import sys
from Bio import SeqIO
from TaxonHolder import Taxon
from ete3 import Tree, faces, TreeStyle, BarChartFace

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

try:
    if mode == "normal" or mode == "special":
        pass
    else:
        raise ValueError("Invalid tag value for mode, make sure to check the list of valid tags and check spelling!")
except ValueError as e:
    print(f"Error: {e}")
    sys.exit()

taxa_dict = {}  # dictionary to store taxa

# read in fasta and parse, then update
for seq_record in SeqIO.parse(file_location, file_format):
    new_taxon = Taxon(seq_record.id, seq_record.seq)
    new_taxon.calculate_all_amino_acid_frequencies()
    if mode == "special":
        new_taxon.calculate_fymink_garp_frequencies()
    taxa_dict[seq_record.id] = new_taxon  # get taxa dict

# tree "growing"
tree = Tree(newwick_tree)

# tree styling
tree_style = TreeStyle()
tree_style.show_scale = False  # do not show scale


def make_face(freq_dict):
    face = BarChartFace(
        values=list(freq_dict.values()),
        labels=list(freq_dict.keys()),
        colors=["blue" for key in freq_dict.keys()],
        max_value=0.2,
        width=100,
        height=50
    )

    return face


# layout function
def make_layout():
    if mode == "normal":
        def layout(node):
            if node.is_leaf():
                taxon = taxa_dict[node.name]
                face = make_face(taxon.freq_dict)
                faces.add_face_to_node(face=face, node=node, column=1, position="aligned")

    elif mode == "special":
        def layout(node):
            if node.is_leaf():
                taxon = taxa_dict[node.name]

                i = 1
                for freq_dict in [taxon.fymink_freq_dict, taxon.garp_freq_dict, taxon.other_freq_dict]:
                    face = make_face(freq_dict)
                    faces.add_face_to_node(face=face, node=node, column=i, position="aligned")
                    i += 1
    else:
        raise ValueError("Invalid tag value for mode, make sure to check the list of valid tags and check spelling!")

    return layout


# render tree
tree.render("test.png", units="px", h=2000, w=2500, dpi=70,  tree_style=tree_style, layout=make_layout())
