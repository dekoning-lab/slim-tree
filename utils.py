# ======================================================================
#  AUTHOR:  Ben Rafferty, Purdue University
#  Copyright (c) 2010  Purdue Research Foundation
#
#  See the file "LICENSE.txt" for information on usage and
#  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
# ======================================================================

"""
Collection of utility functions used in the contactmaps module.
"""

import Bio.PDB
import numpy
import os.path
import tempfile
import urllib.request, urllib.parse, urllib.error
import http.client
import signal

def calc_residue_dist(residue_one, residue_two) :
    """
    Returns the C-alpha distance between two residues.
    """
    diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
    return numpy.sqrt(numpy.vdot(diff_vector, diff_vector))

def count_residues(model):
    """
    Returns the number of residues contained by the given model by determining
    the length of the generator returned by Model.get_residues().
    """
    return sum(1 for e in model.get_residues())

def fetch_pdb(pdbid, filename):
    """
    Attemps to fetch from rcsb.org a pdb file representing the protein with the
    given ID.  If successful, saves the pdb to a temp file and returns the
    filename.
    """
    pdb_url = "https://www.rcsb.org/pdb/files/%s.pdb" % pdbid.lower()
    urllib.request.urlretrieve(pdb_url, filename=filename)

def get_structure(pdbfile):
    """
    For convenience - returns a Bio.PDB.Structure object constructed from the
    given pdb file.
    """
    root, ext = os.path.splitext(pdbfile)
    if ext.lower() != ".pdb":
        raise ValueError("Given file is not a .pdb file.")
    return Bio.PDB.PDBParser().get_structure(root, pdbfile)

def remove_hetero(model):
    """
    Removes all heterogeneous residues from a model.  Use this if your
    calculations should only be considered with the actual protein residues.
    """
    for chain in model:
        hetero_list = []
        for residue in chain:
            id = residue.get_id()
            if id[0] != " ":
                hetero_list.append(id)
        for id in hetero_list:
            chain.detach_child(id)

def residue_array(model):
    """
    Creates a numpy array containing all residues in the given model.
    Arrays can be convenient as opposed to the generator returned by 
    Model.get_residues()
    """
    a = numpy.zeros(count_residues(model), numpy.object)
    for i, res in enumerate(model.get_residues()): a[i] = res
    return a

def check_url(url):
    """
    Check if the given url is accessible.  Returns nothing if accessible,
    raise exception if not.  Use signals to implement timeout because
    timeout argument of HTTPConnection does not exist in Python 2.5.2
    """
    def handler(signum, frame):
        raise Exception("%s: Timeout" % url)
    signal.signal(signal.SIGALRM, handler)
    signal.alarm(10)
    conn = http.client.HTTPConnection(url)
    conn.request("HEAD", "/")
    res = conn.getresponse()
    signal.alarm(0)
    if res.status != 200:
        raise Exception("%s: Status %d returned." % (url, res.status))