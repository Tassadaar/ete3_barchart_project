"""
This program makes a simple unrooted tree from a newwick string, and show said tree in GUI.
"""

from ete3 import Tree
from ete3 import BarChartFace
from FrequencyCalculator import FrequencyCalculator

# variable initiations
#newwick_tree = "input_files/nematode_CH.tree"

# tree "growing"
#unrooted_tree = Tree(newwick_tree)
# print(unrooted_tree)
#unrooted_tree.show()

print(FrequencyCalculator("input_files/nematode.fasta", "fasta").calculate_frequency())
