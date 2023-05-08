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

    def make_face(self):
        taxa_list = FrequencyCalculator(self.file_location, self.file_format).calculate_frequency()
        face_dict = {}

        for taxon in taxa_list:
            new_face = BarChartFace(taxon.freq_list, 100, 50, (0, 0, 255), taxon.aa_list, None, 1)
            face_dict[taxon.name] = new_face

        return face_dict
