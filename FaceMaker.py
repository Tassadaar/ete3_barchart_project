"""
This program makes BarChartFaces and store them in a list

Returns: dictionary(faces)
"""

from ete3 import BarChartFace
from FrequencyCalculator import FrequencyCalculator


class FaceMaker:
    file_location = None
    file_format = None

    def __init__(self, file_location, file_format):
        self.file_location = file_location
        self.file_format = file_format

    def make_barchartface(self):
        taxa_list = FrequencyCalculator(self.file_location, self.file_format).calculate_frequency()
        face_dict = {}

        for taxon in taxa_list:
            face = BarChartFace(list(taxon.freq_dict.values()))
            face.width = 100
            face.height = 50
            face.colors = ["blue" for key in list(taxon.freq_dict.keys())]  # ensure it doesn't run out of colors
            face.labels = list(taxon.freq_dict.keys())
            face.min_value = None
            face.max_value = 0.2

            face_dict[taxon.name] = face

        return face_dict
