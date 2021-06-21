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

class ContactMap(GenericMap):
    """
    Represents a contact map for a given protein model.  Initializer requires
    a Bio.PDB.Model object.  Optionally, a title may be given.  The default
    contact threshold of 10 Angstrom may be overridden.
    """
    def __init__(self, model, title="", threshold=8.0):
        GenericMap.__init__(self, model, title=title, \
                            xlabel="Amino acid index", \
                            ylabel="Amino acid index", \
                            cmap="binary", colorbar=False, threshold=threshold)

    def calc_matrix(self, model, threshold):
        """
        Used internally to compute the matrix data when the object is
        initialized.
        """
        residues = residue_array(model)
        matrix = numpy.ones((len(residues), len(residues)), numpy.bool)
        for row, r1 in enumerate(residues):
            for col, r2 in enumerate(residues[:row]):
                val = calc_residue_dist(r1, r2) < threshold
                matrix[row, col] = val
                matrix[col, row] = val
        return matrix
