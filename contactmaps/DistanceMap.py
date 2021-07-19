# ======================================================================
#  AUTHOR:  Ben Rafferty, Purdue University
#  Copyright (c) 2010  Purdue Research Foundation
#
#  See the file "LICENSE.txt" for information on usage and
#  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
# ======================================================================

from GenericMap import GenericMap
import numpy
from utils import residue_array, calc_residue_dist

class DistanceMap(GenericMap):
    """
    Represents a distance map for a given protein model.  Initializer requires
    a Bio.PDB.Model object.  Optionally, a title may be given.  The boolean
    arguments "contour" and "color" will enable or disable contour lines, and
    select between color or grayscale output.  The interpolation argument
    enables or disables interpolation.  With interpolation disabled, the
    output will have a discrete color value for each entry in the matrix.  If
    enabled, the output will appear smoothed or blended.
    """
    def __init__(self, model, title="", contour=True, color=True, \
                 interpolate=True):
        if color: cmap = "jet"
        else: cmap = "gray"
        if interpolate: interpolation="bilinear"
        else: interpolation="nearest"
        GenericMap.__init__(self, model, title=title, \
                            xlabel="Amino acid index", \
                            ylabel="Amino acid index", \
                            cmap=cmap, contour=contour, \
                            interpolation=interpolation, \
                            colorbarlabel=r"Distance (\u00c5)")

    def calc_matrix(self, model):
        """
        Used internally to compute the matrix data when the object is
        initialized.
        """
        residues = residue_array(model)
        matrix = numpy.zeros((len(residues), len(residues)), numpy.float)
        for row, r1 in enumerate(residues):
            for col, r2 in enumerate(residues[:row]):
                val = calc_residue_dist(r1, r2)
                matrix[row, col] = val
                matrix[col, row] = val
        return matrix

    def get_maximum(self):
        """
        Return the maximum distance currently present in the matrix data.
        """
        return numpy.max(self.matrix)

    def saturate(self, val):
        """
        Modify the distance matrix such that all elements greater than the
        given value are levelled off.  Can be used to force a constant scale
        across multiple distance maps.
        """
        (n, m) = numpy.shape(self.matrix)
        for i in range(n):
            for j in range(m):
                if self.matrix[i][j] > val:
                    self.matrix[i][j] = val

    def set_color(b):
        """
        Use this method to select between color and grayscale outputs. Input
        argument is assumed to be a boolean.
        """
        if b: self.cmap = "jet"
        else: self.cmap = "gray"

    def set_interpolate(b):
        """
        Use this method to turn on/off image interpolation.  Input argument is
        assumed to be a boolean."
        """
        if b: self.interpolation = "bilinear"
        else: self.interpolation = "nearest"

    def set_contour(b):
        """
        Use this method to turn on/off contour lines.  Input argument is
        assumed to be a boolean."
        """

