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


def main():
    # specify tags
    parser = argparse.ArgumentParser(description="Tree making")

    parser.add_argument("-t", "--tree", required=True)
    parser.add_argument("-n", "--file", required=True)
    parser.add_argument("-f", "--format", required=True)
    parser.add_argument("-o", "--output", type=str, default="tree")
    parser.add_argument("-s", "--subset", type=str, default="none")
    parser.add_argument("-m", "--frequency_type", type=str, default="absolute")
    parser.add_argument("-g", "--outgroup_reps", type=str, default="none")

    args = parser.parse_args()

    tree = Tree(args.tree)  # tree "growing"
    outgroup_reps = [word.upper() for word in args.outgroup_reps.split(",")]
    subset = args.subset
    frequency_type = args.frequency_type

    subsets = [word.upper() for word in subset.split(",")]  # a list assumed to have two items
    frequency_types = ["absolute", "relative"]
    all_amino_acids = "ACDEFGHIKLMNPQRSTVWY"

    # exception handling for flags
    try:
        # check if subset is in the right format and if it contains valid amino acids
        if len(subsets) != 2:
            raise ValueError("Subsets contain less or more than 2 groupings, make sure to check format!")

        for subset in subsets:

            if not set(subset).issubset(set(all_amino_acids)):
                raise ValueError("Subsets contain invalid amino acid(s), make sure to check spelling!")

        # check if frequency_type is valid
        if frequency_type not in frequency_types:
            raise ValueError("Invalid tag for frequency type, make sure to check the list of valid tags and spelling!")

        # check if the provided outgroup is valid
        if outgroup_reps[0] != "NONE":

            if outgroup_reps[0] not in tree.get_leaf_names():
                raise ValueError("Invalid outgroup, make sure to check spelling!")

            elif len(outgroup_reps) > 1:

                if outgroup_reps[1] not in tree.get_leaf_names():
                    raise ValueError("Invalid outgroup, make sure to check spelling!")

    except ValueError as e:
        print(f"Error: {e}")  # print error message
        sys.exit()  # terminate the program

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

    # tree rooting
    if outgroup_reps[0] != "NONE":  # check if rooting is required
        tree = root(tree, outgroup_reps)

    # experiments
    tree = root(tree, ["NEMATODE", "LOCUST"])
    tree.get_tree_root().unroot()

    # tree styling
    tree.ladderize()
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


def root(tree, outgroup_reps):
    if len(outgroup_reps) > 1:  # check if the outgroup is only one taxon
        # if not, make it monophyletic
        common_ancestor = tree.get_common_ancestor(outgroup_reps[0], outgroup_reps[1])
        children = common_ancestor.get_children()

        # this is a workaround for how ete3 handles unrooted trees, as you cannot reroot to the current "root"
        if common_ancestor.is_root():
            print("Common ancestor is root, taking a detour")
            fringe_case = True

            # find the first ingroup taxon for rooting
            for child in children:
                leaves = child.get_leaf_names()

                if all(rep not in leaves for rep in outgroup_reps):
                    tree.set_outgroup(leaves[0])
                    common_ancestor = tree.get_common_ancestor(outgroup_reps[0], outgroup_reps[1])
                    fringe_case = False
                    break

            # if all immediate children of the root contains an outgroup taxon
            while fringe_case:
                ingroup = input("Enter a leaf of the ingroup: ")

                if ingroup in common_ancestor.get_leaf_names() and ingroup not in outgroup_reps:
                    tree.set_outgroup(ingroup)
                    common_ancestor = tree.get_common_ancestor(outgroup_reps[0], outgroup_reps[1])
                    break
                else:
                    print("Invalid ingroup, check spelling!")

        tree.set_outgroup(common_ancestor)
    else:
        tree.set_outgroup(outgroup_reps[0])  # just make that one taxon the outgroup

    return tree

if __name__ == "__main__":
    main()
