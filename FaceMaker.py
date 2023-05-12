"""
This program makes BarChartFaces and store them in a list

Returns: dictionary(faces)
"""

from ete3 import BarChartFace
from FrequencyCalculator import FrequencyCalculator


class FaceMaker:

    @staticmethod
    def make_barchartface(freq_dict):

        face = BarChartFace(list(freq_dict.values()))
        face.width = 100
        face.height = 50
        face.colors = ["blue" for key in list(freq_dict.keys())]  # ensure it doesn't run out of colors
        face.labels = list(freq_dict.keys())
        face.min_value = None
        face.max_value = 0.2

        return face


