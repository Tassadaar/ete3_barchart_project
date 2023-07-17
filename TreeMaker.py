"""
This program makes a simple unrooted tree from a newwick string, \
        make BarChartFace for each leaf node in layout and render said tree to a PNG image.

mandatory arguments: -t or --tree, newick tree
                     -n or --filename, path to alignment file
                     -f or --format, format of said alignment file (only tested for fasta)
optional arguments: -m or --mode, type of BarChartFace to load (all or fymink/garp subsets)
"""

import argparse
import sys
from Bio import SeqIO
from Taxon import Taxon
from ete3 import Tree, faces, TreeStyle, BarChartFace, TextFace


# argument parsing validation functions

# check if subset is in the right format and if it contains valid amino acids
def validate_subsets(subset):
    # ensure subsets are processed in upper case
    subsets = [aa_group.upper() for aa_group in subset.split(",")]

    # raise error if number of subsets is not exactly 2
    if len(subsets) != 2:
        raise argparse.ArgumentTypeError("There must be exactly two \
                                          subsets of amino acids")
    # raise error if any subset includes a non-aminoacid char
    if any(char not in "ACDEFGHIKLMNPQRSTVWY" for aa_group in subsets for char in aa_group):
        raise argparse.ArgumentTypeError("At least one of the subsets \
                                          contains invalid amino acid(s)")
    # return uppercase subsets if they passed the checks
    return subsets


# check if frequency_type is valid
def validate_frequency(freq_type):
    # check if freq_type is an allowed one
    if freq_type not in ["absolute", "relative"]:
        raise argparse.ArgumentTypeError(
            "Invalid frequency type, only 'absolute' and 'relative' are allowed"
        )
    # return freq_type if it passed the check
    return freq_type


# check if the provided outgroup is valid
def validate_outgroup(outgroup_reps, leaves):
    outgroup_reps = outgroup_reps.split(",")
    # if any given representative is not found,
    if any(rep not in leaves for rep in outgroup_reps):
        # raise an error
        raise argparse.ArgumentTypeError("Given outgroup not found in tree, \
                                          make sure to check spelling!")
    # if no errors are found,
    # return the representatives as a list
    return outgroup_reps


# specify options, disable for debugging
parser = argparse.ArgumentParser(description="Tree making")

parser.add_argument("-t", "--tree", required=True)
parser.add_argument("-n", "--file", required=True)
parser.add_argument("-f", "--format", required=True)
parser.add_argument("-o", "--output", type=str, default="tree")
parser.add_argument("-s", "--subsets", type=validate_subsets, nargs='?')
parser.add_argument("-m", "--frequency_type", type=validate_frequency, default="absolute")
parser.add_argument("-g", "--outgroup_reps", type=str, nargs='?')
parser.add_argument("-c", "--show_chi2_score", type=bool, default=False)

args = parser.parse_args()


# encapsulated barchart face function
def get_barchart_face(freq_dict, max_value):
    face = BarChartFace(
        values=[abs(x) for x in freq_dict.values()],
        labels=[" " for x in freq_dict.keys()],
        label_fsize=9,  # this value dictates scaling if bar widths are uniform
        colors=["blue" if f > 0 else "red" for f in freq_dict.values()],
        width=40,  # when below a certain threshold, all the bar widths are scaled to be uniform
        height=50,
        max_value=max_value,
    )
    return face


# layout function
def layout(node):
    # exit the function without returning anything
    # if the node is not a leaf node
    if not node.is_leaf():
        return

    taxon = taxa_dict[node.name]
    dict_list = []

    # calculate frequency values to be displayed
    # if no subsets are given (default)
    if args.subsets is None:

        dict_list.append(taxon.display_freqs)
        if args.frequency_type == "absolute":
            # dict_list.append(taxon.freq_dict)
            max_value = 0.20
        elif args.frequency_type == "relative":
            # dict_list.append(taxon.get_all_relative_freq(avg_freq_dict))
            max_value = 0.05

    # if subsets are specified
    else:
        taxon.set_subset_abs_freq(subsets)

        if args.frequency_type == "absolute":
            dict_list = taxon.get_subset_abs_freq()
            max_value = 0.20
        elif args.frequency_type == "relative":
            dict_list = taxon.get_subset_relative_freq(subsets, avg_freq_dict)
            max_value = 0.05

    # determine barchartfaces to be displayed
    i = 1
    for freq_dict in dict_list:
        face = get_barchart_face(freq_dict, max_value)

        if node.name == tree.get_leaf_names()[-1]:
            face.labels = list(freq_dict.keys())

        # ensure a healthy width of gap between the tree and the faces
        if i == 1:
            face.margin_left = 50
        face.margin_right = 10

        if subsets[0] != "NONE" and i != len(dict_list):
            face.scale_fsize = 1  # this ensures that only one set of scale is shown for all columns
        faces.add_face_to_node(face=face, node=node, column=i, position="aligned")
        i += 1

    # determine textfaces to be displayed
    if args.show_chi2_score is True:
        text_face = TextFace(
            taxon.calculate_chi_square(avg_freq_dict)
        )
        text_face.margin_left = 50
        faces.add_face_to_node(face=text_face, node=node, column=i, position="aligned")


# if user specified outgroup taxa in the flags then root accordingly
def root(tree, outgroup_reps):
    # check if the outgroup is only one taxon
    if len(outgroup_reps) == 1:
        tree.set_outgroup(outgroup_reps[0])  # just make that one taxon the outgroup

        return tree

    # if not, make it monophyletic
    common_ancestor = tree.get_common_ancestor(*outgroup_reps)

    if not common_ancestor.is_root():
        tree.set_outgroup(common_ancestor)

        return tree

    # this is a workaround for how ete3 handles unrooted trees,
    # as you cannot reroot to the current "root"
    # the user needs to provide a proper ingroup based on biological information
    print("\nCommon ancestor is root, taking a detour")
    ingroup = input("\nEnter an ingroup taxon: ")

    # reject non-leaf inputs
    while ingroup not in tree.get_leaf_names():
        print("You entered a taxon that is not a part of the tree, check spelling!")
        ingroup = input("\nEnter an ingroup taxon: ")

    # reject outgroup inputs
    while ingroup in outgroup_reps:
        print("You entered an outgroup taxon, check spelling!")
        ingroup = input("\nEnter an ingroup taxon: ")

    tree.set_outgroup(ingroup)
    # get the desired monophyletic outgroup
    common_ancestor = tree.get_common_ancestor(*outgroup_reps)
    tree.set_outgroup(common_ancestor)

    return tree


def main(args):
    # load newick tree
    tree = Tree(args.tree)  # tree "growing"
    leaves = tree.get_leaf_names()

    # validate given outgroup representatives if they are given
    if args.outgroup_reps is not None:
        # check if they are given as the script expects
        try:
            outgroup_reps = validate_outgroup(args.outgroup_reps, leaves)
        except argparse.ArgumentTypeError as e:
            print(e)
            sys.exit()
        # then root the tree
        tree = root(tree, outgroup_reps)

    # set tree style
    tree.ladderize()
    tree_style = TreeStyle()
    tree_style.show_scale = False  # do not show scale

    # load fasta and store Taxon objects
    taxa_dict = {}  # dictionary to store taxa
    all_seq = ""  # string to hold all the sequences in the alignment
    for seq_record in SeqIO.parse(args.file, args.format):
        new_taxon = Taxon(seq_record.id, seq_record.seq)
        all_seq += new_taxon.seq
        taxa_dict[seq_record.id] = new_taxon

    # calculate mean frequencies of entire alignment
    # if relative frequencies or show_chi2_score is invoked
    if args.frequency_type == "relative" or args.show_chi2_score is True:
        avg_freq_dict = {aa: all_seq.count(aa) / len(all_seq) for aa in all_amino_acids}

    # determine what frequencies we want to display
    for tax in taxa_dict.values():

        if args.subsets is None:
            if args.frequency_type == 'absolute':
                tax.set_display_freqs('absolute')
            elif args.frequency_type == 'relative':
                tax.set_display_freqs('relative', avg_freq_dict)

        elif args.subsets:
            if args.frequency_type == "absolute":
                tax.set_subset_abs_freq(subsets)
            elif args.frequency_type == "relative":
                tax.set_subset_rel_freq(subsets, avg_freq_dict)

    # render tree
    tree.render(
        file_name=args.output + ".png",
        units="px", h=200 * len(leaves),
        tree_style=tree_style,
        layout=layout
    )


# Guard against undesired invocation upon import
if __name__ == "__main__":
    all_amino_acids = "ACDEFGHIKLMNPQRSTVWY"

    tree = None
    subsets = None
    args.frequency_type = None
    args.show_chi2_score = None
    taxa_dict = None
    avg_freq_dict = None

    main(args)
