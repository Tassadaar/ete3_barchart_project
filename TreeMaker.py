"""
This program makes a simple unrooted tree from a newwick string, make BarChartFace for each leaf node in layout and
render said tree to a PNG image.

mandatory arguments: -t or --tree, newwick tree
                     -n or --filename, path to alignment file
                     -f or --format, format of said alignment file (only tested for fasta)
optional arguments: -m or --mode, type of BarChartFace to load (all or fymink/garp subsets)
"""

import argparse
import sys
from Bio import SeqIO
from TaxonHolder import Taxon
from ete3 import Tree, faces, TreeStyle, BarChartFace

# specify tags
parser = argparse.ArgumentParser(description="Tree making")

parser.add_argument("-t", "--tree", required=True)
parser.add_argument("-n", "--file", required=True)
parser.add_argument("-f", "--format", required=True)
parser.add_argument("-m", "--mode", type=str, default="normal")

args = parser.parse_args()

mode = args.mode

modes = ["normal", "special", "inverted", "fgInverted"]

try:
    if mode in modes:
        pass
    else:
        raise ValueError("Invalid tag value for mode, make sure to check the list of valid tags and check spelling!")
except ValueError as e:
    print(f"Error: {e}")
    sys.exit()

taxa_dict = {}  # dictionary to store taxa
all_seq = ""  # string to hold all the sequences in the alignment

# read in fasta and parse, then update
for seq_record in SeqIO.parse(args.file, args.format):
    new_taxon = Taxon(seq_record.id, seq_record.seq)
    all_seq += seq_record.seq
    new_taxon.set_all_amino_acid_frequencies()
    if mode == "special":
        new_taxon.set_fymink_garp_frequencies()
    taxa_dict[seq_record.id] = new_taxon  # get taxa dict

if mode == "inverted" or mode == "fgInverted":
    all_seq = all_seq.replace("-", "")
    all_amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    unsorted_avg_freq_dict = {aa: all_seq.count(aa) / len(all_seq) for aa in all_amino_acids}
    sorted_keys = sorted(unsorted_avg_freq_dict.keys())
    sorted_avg_freq_dict = {key: unsorted_avg_freq_dict[key] for key in sorted_keys}

    if mode == "inverted":
        for name, taxon in taxa_dict.items():
            taxon.set_frequency_deviations(sorted_avg_freq_dict)

    if mode == "fgInverted":
        for name, taxon in taxa_dict.items():
            taxon.set_fymink_garp_frequencies()
            taxon.set_fymink_garp_frequency_deviations(sorted_avg_freq_dict)

# tree "growing"
tree = Tree(args.tree)

# tree styling
tree_style = TreeStyle()
tree_style.show_scale = False  # do not show scale


# layout function
def layout(node):
    if node.is_leaf():
        taxon = taxa_dict[node.name]
        dict_list = []
        max_value = 0.2
        if mode == "normal":
            dict_list.append(taxon.freq_dict)
        elif mode == "special":
            dict_list.extend([taxon.fymink_freq_dict, taxon.garp_freq_dict, taxon.other_freq_dict])
        elif mode == "inverted":
            dict_list.append(taxon.all_freq_deviation_dict)
            max_value = 0.05
        elif mode == "fgInverted":
            dict_list.extend([taxon.fymink_freq_deviation_dict,
                              taxon.garp_freq_deviation_dict,
                              taxon.other_freq_deviation_dict])
            max_value = 0.05

        i = 1
        for freq_dict in dict_list:
            face = BarChartFace(
                values=[abs(x) for x in freq_dict.values()],
                labels=list(freq_dict.keys()),
                colors=["blue" if f > 0 else "red" for f in freq_dict.values()],
                width=100,
                height=50
            )
            face.max_value = max_value
            faces.add_face_to_node(face=face, node=node, column=i, position="aligned")
            i += 1


# render tree
tree.render("test.png", units="px", h=2000, w=2500, dpi=70,  tree_style=tree_style, layout=layout)
