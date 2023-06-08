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
from Taxon import Taxon
from ete3 import Tree, faces, TreeStyle, BarChartFace

# specify tags
parser = argparse.ArgumentParser(description="Tree making")

parser.add_argument("-t", "--tree", required=True)
parser.add_argument("-n", "--file", required=True)
parser.add_argument("-f", "--format", required=True)
parser.add_argument("-o", "--output", type=str, default="barchart_tree")
parser.add_argument("-s", "--subset", type=str, default="none")
parser.add_argument("-m", "--frequency_type", type=str, default="absolute")

args = parser.parse_args()

subset = args.subset
frequency_type = args.frequency_type

subsets = [word.upper() for word in subset.split(",")]  # a list assumed to have two items
frequency_types = ["absolute", "relative"]

all_amino_acids = "ACDEFGHIKLMNPQRSTVWY"

try:
    if len(subsets) != 2:
        raise ValueError("Subsets contain less or more than 2 groupings, make sure to check format!")
    for subset in subsets:
        if not set(subset).issubset(set(all_amino_acids)):
            raise ValueError("Subsets contain invalid amino acid(s), make sure to check spelling!")
    if frequency_type not in frequency_types:
        raise ValueError("Invalid tag for frequency type, make sure to check the list of valid tags and spelling!")
except ValueError as e:
    print(f"Error: {e}")
    sys.exit()

taxa_dict = {}  # dictionary to store taxa
all_seq = ""  # string to hold all the sequences in the alignment

# read in fasta and parse, then update
for seq_record in SeqIO.parse(args.file, args.format):
    new_taxon = Taxon(seq_record.id, seq_record.seq)
    all_seq += seq_record.seq
    new_taxon.set_aa_abs_freq()
    taxa_dict[seq_record.id] = new_taxon  # get taxa dict

# calculate relative frequencies if specified
if frequency_type == "relative":
    # calculate average frequency
    all_seq = all_seq.replace("-", "")
    avg_freq_dict = {aa: all_seq.count(aa) / len(all_seq) for aa in all_amino_acids}

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

        if subset == "none":

            if frequency_type == "absolute":
                dict_list.append(taxon.get_aa_abs_freq())
            elif frequency_type == "relative":
                dict_list.append(taxon.get_all_relative_freq(avg_freq_dict))
                max_value = 0.05

        else:
            taxon.set_subset_abs_freq(subsets)

            if frequency_type == "absolute":
                dict_list = taxon.get_subset_abs_freq()
            elif frequency_type == "relative":
                dict_list = taxon.get_subset_relative_freq(subsets, avg_freq_dict)
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
tree.render(args.output + ".png", units="px", h=2000, w=2500, dpi=70,  tree_style=tree_style, layout=layout)
