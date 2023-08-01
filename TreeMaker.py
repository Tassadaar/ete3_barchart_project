"""
This program makes a simple unrooted tree from a newwick string, make BarChartFace for each leaf node in layout and
render said tree to a PNG image.

mandatory arguments: -t or --tree, newick tree
                     -n or --filename, path to alignment file
                     -f or --format, format of said alignment file (only tested for fasta)
optional arguments: -m or --mode, type of BarChartFace to load (all or fymink/garp subsets)
"""

import argparse
import sys
import os
from Bio import SeqIO
from Taxon import Taxon
from ete3 import Tree, faces, TreeStyle, BarChartFace, TextFace


# check if subset is in the right format and if it contains valid amino acids
def validate_subsets(subset):
    subsets = [aa_group.upper() for aa_group in subset.split(",")]

    if len(subsets) != 2:
        raise argparse.ArgumentTypeError(f"Subsets does not have exactly 2 groupings: {subset}")

    if any(char not in "ACDEFGHIKLMNPQRSTVWY" for aa_group in subsets for char in aa_group):
        raise argparse.ArgumentTypeError(f"At least one subset contains invalid amino acid(s)")

    return subsets


# check if frequency_type is valid
def validate_frequency(freq_type):
    if freq_type not in ["absolute", "relative"]:
        raise argparse.ArgumentTypeError(
            "Invalid frequency_type, only \"absolute\" or \"relative\" are allowed!"
        )
    return freq_type


# check if the provided outgroup is valid
def validate_outgroup(outgroup_reps, leaves):
    outgroup_reps = outgroup_reps.split(",")

    if any(rep not in leaves for rep in outgroup_reps):
        raise argparse.ArgumentTypeError("Invalid outgroup, make sure to check spelling!")

    return outgroup_reps


"""
the above functions are for argument validation
"""


# layout function
def layout(node):
    # only consider leaf nodes
    if not node.is_leaf():
        return

    # retrieve Taxon() object from memory
    taxon = taxa_dict[node.name]

    # taxon.display_freqs holds a list of dictionaries,
    # where each dictionary is 'aa' : frequency
    # one dict per subset group
    i = 1
    for freq_dict in taxon.display_freqs:

        # make a face for each subset
        face = get_barchart_face(freq_dict, taxon.display_max_value)

        # by default, all faces have empty labels
        # here we make sure the bottom taxon in the tree
        # do have annotated amino acid labels
        root = node.get_tree_root()
        if node.name == root.get_leaf_names()[-1]:
            face.labels = list( freq_dict.keys() )

        # ensure a healthy margin between the tree and the first face
        if i == 1:
            face.margin_left = 50
        # a smaller margin for between faces
        face.margin_right = 10

        # suppress y-axis labels for faces that are not the last face
        if len(freq_dict) > 1 and i != len(taxon.display_freqs):
            face.scale_fsize = 1

        # display face
        faces.add_face_to_node(face=face, node=node, column=i, position="aligned")
        i += 1

    # display per taxon chi2 score if desired by the user
    if args.show_chi2_score:
        face = TextFace(taxon.chi_square_score)
        face.margin_left = 50
        faces.add_face_to_node(face=face, node=node, column=i, position="aligned")


def get_barchart_face(freq_dict, max_value):
    face = BarChartFace(
        values = [ abs(x) for x in freq_dict.values() ],
        labels = [ ' ' for i in range(len(freq_dict)) ],
        label_fsize = 9,  # label font size
        colors = [ "blue" if f > 0 else "red" for f in freq_dict.values() ],
        width = 40,  # bar widths
        height = 20,
        max_value = max_value,
    )
    # face.inner_background.color = 'lightgrey'

    return face


# if user specified outgroup taxa in the flags then root accordingly
def root(tree, outgroup_reps):
    # check if the outgroup is only one taxon
    if len(outgroup_reps) == 1:
        # just make that one taxon the outgroup
        tree.set_outgroup(outgroup_reps[0])
        return tree

    # if requested outgroup has multiple taxa,
    # get their common ancestor node
    common_ancestor = tree.get_common_ancestor(*outgroup_reps)

    # if common ancestor is not the root node,
    # setting the outgroup with common ancestor is sufficient
    if not common_ancestor.is_root():
        tree.set_outgroup(common_ancestor)
        return tree

    # if it is the root node, the script needs to know from the user
    # at least one in-group taxon to know what is the outgroup
    # this is necessary as setting the outgroup with the root node doesn't work
    print("Outgroup common ancestor is root node, additional information is required")
    ingroup_taxon = input("\nEnter an ingroup taxon: ")

    # user input error-handling:
    # reject non-leaf inputs
    while ingroup_taxon not in tree.get_leaf_names():
        print("You entered a taxon that is not a part of the tree, check spelling!")
        ingroup_taxon = input("\nEnter an ingroup taxon: ")

    # reject outgroup inputs
    while ingroup_taxon in outgroup_reps:
        print("You entered an outgroup taxon, check spelling!")
        ingroup_taxon = input("\nEnter an ingroup taxon: ")

    # setting the ougroup with the ingroup-taxon ensures that
    # the real outgroup common ancestor is not the root node
    tree.set_outgroup(ingroup_taxon)
    common_ancestor = tree.get_common_ancestor(*outgroup_reps)
    tree.set_outgroup(common_ancestor)

    return tree


def main(args):
    # needs to be global because layout function needs it
    global taxa_dict
 
    # load tree and get leaf names
    tree = Tree(args.tree, format=1)
    leaves = tree.get_leaf_names()

    # if outgroups are given by the user,
    if args.outgroup_reps is not None:
        # validate outgroup_reps user input
        try:
            outgroup_reps = validate_outgroup(args.outgroup_reps, leaves)
        except argparse.ArgumentTypeError as e:
            print(e)
            sys.exit()
        # and root the tree
        tree = root(tree, outgroup_reps)

    # collect all Taxon() objects in memory
    taxa_dict = {}
    # collect all alignment sequences into a single string in memory
    all_seq = ""

    # parse fasta and generate Taxon() objects into memory
    for seq_record in SeqIO.parse(args.file, args.format):
        new_taxon = Taxon(seq_record.id, seq_record.seq)
        all_seq += new_taxon.seq
        taxa_dict[seq_record.id] = new_taxon

    # calculate mean frequencies if relative frequencies
    # or show_chi2_score is desired by the user
    if args.frequency_type == "relative" or args.show_chi2_score is True:
        # calculate mean frequencies
        avg_freq_dict = {aa: all_seq.count(aa) / len(all_seq) for aa in "ACDEFGHIKLMNPQRSTVWY"}

    # determine which frequencies to display
    for taxon in taxa_dict.values():

        # 4 possibilities:
        if args.subsets is None:
            ## 1 - no subsets, absolute frequencies
            ## the display frequencies are absolute by default
            ## when Taxon() objects are initialized
            ## no further action required

            ## 2 - no subsets, relative frequencies
            ## the display frequencies need to be reset
            ## to relative frequencies
            if args.frequency_type == "relative":
                taxon.set_all_relative_freq(avg_freq_dict)
                taxon.display_max_value = 0.05

        else:
            ## 3 - subsets, absolute frequencies
            ## the display frequencies need to be reset
            ## to a list of absolute frequencies dictionaries
            taxon.set_subset_abs_freq(args.subsets)

            ## 4 - subsets, relative frequencies
            ## the display frequencies need to be reset
            ## to a list of relative frequency dictionaries
            if args.frequency_type == "relative":
                taxon.set_subset_relative_freq(args.subsets, avg_freq_dict)
                taxon.display_max_value = 0.05

        if args.show_chi2_score is True:
            taxon.calculate_chi_square(avg_freq_dict)

    # tree styling
    tree.ladderize(direction=1)
    tree_style = TreeStyle()
    tree_style.show_scale = False  # do not show scale

    # render tree
    tree.render(
        file_name = args.output + '.png',
        units = "px",
        h = 200 * len(leaves),
        tree_style = tree_style,
        layout = layout
    )


# guard against undesired invocation upon import
if __name__ == "__main__":

    # specify options, disable for debugging
    parser = argparse.ArgumentParser(description="Tree making")

    parser.add_argument("-t", "--tree", required=True)
    parser.add_argument("-n", "--file", required=True)
    parser.add_argument("-f", "--format", required=True)
    parser.add_argument("-o", "--output", type=str, default="tree")
    parser.add_argument("-s", "--subsets", type=validate_subsets, nargs="?")
    parser.add_argument("-m", "--frequency_type", type=validate_frequency, default="absolute")
    parser.add_argument("-g", "--outgroup_reps", type=str, nargs="?")
    parser.add_argument("-c", "--show_chi2_score", action='store_true')

    # not sure why this is here
    # taxa_dict = None

    args = parser.parse_args()
    main(args)
