from unittest import TestCase
from ete3 import Tree
from TreeMaker import root


class Test(TestCase):
    test_tree = "(((HONEYBEE:1,KILLERBEE:1)1:1,NEMATODE:1)1:0.5,((FRUITFLY:1,LOCUST:1)1:1,(BRINESHR:1," \
                "((CHICKEN:1,SEAURCHIN:1)1:1,ALLOMYCES:1)1:1)1:1)1:0.5);"

    result_tree = "((((CHICKEN, SEAURCHIN), ALLOMYCES), BRINESHR), (((KILLERBEE, HONEYBEE), NEMATODE), " \
                  "(FRUITFLY, LOCUST)));"

    # Scenario 1: the tree is unrooted
    def test_root_scenario_1(self):
        input_tree = Tree(self.test_tree)
        input_tree.unroot()
        result_tree = Tree(self.result_tree)
        root(input_tree, ["BRINESHR", "CHICKEN"])
        rf_distance = input_tree.robinson_foulds(result_tree)[0]
        self.assertEqual(0, rf_distance, f"The distance should have been 0, but it is {rf_distance}")

    # Scenario 2: the tree is rooted but only one child contains outgroups
    def test_root_scenario_2(self):
        input_tree = Tree(self.test_tree)
        result_tree = Tree(self.result_tree)
        root(input_tree, ["BRINESHR", "CHICKEN"])
        rf_distance = input_tree.robinson_foulds(result_tree)[0]
        self.assertEqual(0, rf_distance, f"The distance should have been 0, but it is {rf_distance}")

    # Scenario 3: the tree is rooted but both children contain outgroups
    def test_root_scenario_3(self):
        input_tree = Tree(self.test_tree)
        result_tree = Tree(self.result_tree)
        root(input_tree, ["NEMATODE", "LOCUST"])
        rf_distance = input_tree.robinson_foulds(result_tree)[0]
        self.assertEqual(0, rf_distance, f"The distance should have been 0, but it is {rf_distance}")




