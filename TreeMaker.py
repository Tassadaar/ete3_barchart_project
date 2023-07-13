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
from Bio import SeqIO
from Taxon import Taxon
from ete3 import Tree, faces, TreeStyle, BarChartFace, TextFace


def main(args):
    global tree, subsets, frequency_type, chi2_score, taxa_dict, avg_freq_dict

    tree = Tree(args.tree)  # tree "growing"
    leaves = tree.get_leaf_names()

    try:
        outgroup_reps = validate_outgroup(args.outgroup_reps, leaves)
    except argparse.ArgumentTypeError as e:
        print(e)
        sys.exit()

    frequency_type = args.frequency_type
    subsets = args.subset  # a list assumed to have two items
    chi2_score = args.show_chi2_score

    taxa_dict = {}  # dictionary to store taxa
    all_seq = ""  # string to hold all the sequences in the alignment

    # read in fasta and parse, then update
    for seq_record in SeqIO.parse(args.file, args.format):
        new_taxon = Taxon(seq_record.id, seq_record.seq)
        all_seq += seq_record.seq
        taxa_dict[seq_record.id] = new_taxon  # get taxa dict

    # calculate relative frequencies if specified
    if frequency_type == "relative" or chi2_score is True:
        # calculate average frequency
        all_seq = all_seq.replace("-", "")
        avg_freq_dict = {aa: all_seq.count(aa) / len(all_seq) for aa in all_amino_acids}

    # tree rooting
    if outgroup_reps[0] != "NONE":  # check if rooting is required
        tree = root(tree, outgroup_reps)

    # tree styling
    tree.ladderize()
    tree_style = TreeStyle()
    tree_style.show_scale = False  # do not show scale

    # render tree
    tree.render(
        file_name=args.output + ".png",
        units="px", h=200 * len(leaves),
        tree_style=tree_style,
        layout=layout
    )


# layout function
def layout(node):

    if not node.is_leaf():
        return

    taxon = taxa_dict[node.name]
    dict_list = []
    max_value = 0.2

    if subsets[0] == "NONE":

        if frequency_type == "absolute":
            dict_list.append(taxon.freq_dict)
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
            labels=[" " for x in freq_dict.keys()],
            label_fsize=9,  # this value dictates scaling if bar widths are uniform
            colors=["blue" if f > 0 else "red" for f in freq_dict.values()],
            width=40,  # when below a certain threshold, all the bar widths are scaled to be uniform
            height=50,
            max_value=max_value,
        )

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

    if chi2_score is True:
        text_face = TextFace(
            taxon.calculate_chi_square(
                avg_freq_dict
            )
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

    # this is a workaround for how ete3 handles unrooted trees, as you cannot reroot to the current "root"
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
    common_ancestor = tree.get_common_ancestor(*outgroup_reps)  # get the desired monophyletic outgroup
    tree.set_outgroup(common_ancestor)

    return tree


"""
the following functions are for argument validation
"""


# check if subset is in the right format and if it contains valid amino acids
def validate_subset(subset):
    allowed_characters = set(all_amino_acids)
    subsets = subset.upper().split(",")

    if len(subsets) != 2:
        raise argparse.ArgumentTypeError(f"Subset contains less or more than 2 groupings: {subset}")
    elif not set(subsets[0]).issubset(allowed_characters) or not set(subsets[1]).issubset(allowed_characters):
        raise argparse.ArgumentTypeError(f"Subsets contain invalid amino acid(s)")

    return subsets


# check if frequency_type is valid
def validate_frequency(freq_type):
    allowed_types = ["absolute", "relative"]

    if freq_type not in allowed_types:
        raise argparse.ArgumentTypeError(
            "Invalid tag for frequency type, make sure to check the list of valid tags and spelling!"
        )

    return freq_type


# check if the provided outgroup is valid
def validate_outgroup(outgroup_reps, leaves):
    outgroup_reps = outgroup_reps.split(",")

    if outgroup_reps[0] != "NONE":

        if outgroup_reps[0] not in leaves:
            raise argparse.ArgumentTypeError("Invalid outgroup, make sure to check spelling!")

        elif len(outgroup_reps) > 1:

            if outgroup_reps[1] not in leaves:
                raise argparse.ArgumentTypeError("Invalid outgroup, make sure to check spelling!")

    return outgroup_reps


# Guard against undesired invocation upon import
if __name__ == "__main__":
    all_amino_acids = "ACDEFGHIKLMNPQRSTVWY"

    tree = None
    subsets = None
    frequency_type = None
    chi2_score = None
    taxa_dict = None
    avg_freq_dict = None

    # specify options, disable for debugging
    parser = argparse.ArgumentParser(description="Tree making")

    parser.add_argument("-t", "--tree", required=True)
    parser.add_argument("-n", "--file", required=True)
    parser.add_argument("-f", "--format", required=True)
    parser.add_argument("-o", "--output", type=str, default="tree")
    parser.add_argument("-s", "--subset", type=validate_subset)
    parser.add_argument("-m", "--frequency_type", type=validate_frequency, default="absolute")
    parser.add_argument("-g", "--outgroup_reps", type=str, default="none")
    parser.add_argument("-c", "--show_chi2_score", type=bool, default=False)

    # emulating commandline arguments for debugging, disable for normal execution
    # sys.argv = [
    #             "TreeMaker.py",
    #             "-t", "Martijn_et_al_2019/alphaproteobacteria_untreated.aln.treefile",
    #             "-n", "Martijn_et_al_2019/alphaproteobacteria_untreated.aln",
    #             "-f", "fasta",
    #             "-s", "fymink,garp",
    #             "-m", "relative",
    #             "-g", "Dechloromonas_aromatica_RCB,Pseudomonas_aeruginosa_PA7"
    # ]

    arguments = parser.parse_args()
    main(arguments)
