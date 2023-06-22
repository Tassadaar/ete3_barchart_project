"""
This is a variant of the root() function in TreeMaker which asks for user prompt when it encounters a scenario 3.
A better solution was realized, but the implementation of the user prompt is worth archiving.
"""


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

                if type(ingroup) != str:
                    print("Invalid ingroup, please input a string!")
                elif ingroup in common_ancestor.get_leaf_names() and ingroup not in outgroup_reps:
                    tree.set_outgroup(ingroup)
                    common_ancestor = tree.get_common_ancestor(outgroup_reps[0], outgroup_reps[1])
                    break
                else:
                    print("Invalid ingroup, check spelling!")

        tree.set_outgroup(common_ancestor)
    else:
        tree.set_outgroup(outgroup_reps[0])  # just make that one taxon the outgroup

    return tree
