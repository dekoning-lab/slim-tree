#! /usr/bin/env python

import optparse
import sys
import numpy as np

parser = optparse.OptionParser(usage=
"%prog input.pdb [-o output] [-t threshold (A)] [-m model_num] [-c chain_id]")
parser.add_option("-o", dest="outfile", default="contactmap.png",
                  help="Output file name (default: contactmap.png)")
parser.add_option("-t", dest="threshold", type="float", default=7.0,
                  help="Contact threshold in angstroms (default: 7.0)")
parser.add_option("-m", dest="model_num", default=0, type="int",
                  help="Model number from pdb file to be used (default: 0)")
parser.add_option("-c", dest="chain_id", default=" ", type="string",
                  help="Chain ID from pdb file to be used (default: None)")
options, args = parser.parse_args()

if len(args) != 1:
    parser.print_usage()
    sys.exit(1)

from ContactMap import ContactMap
from utils import get_structure 

structure = get_structure(args[0])
model = structure[options.model_num]
contactMap = ContactMap(model, threshold=options.threshold, chain_id = options.chain_id)
loc = np.where(contactMap.get_data())
loc = list(zip(loc[0], loc[1]))

def getkey(item):
    return item[0]
loc = sorted(loc, key = getkey)
np.savetxt(options.outfile+".csv", loc, fmt = "%d", delimiter= ',', newline = ",")


contactMap.print_figure(options.outfile)



