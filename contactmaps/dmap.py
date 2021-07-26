#! /usr/bin/env python

import optparse
import sys

parser = optparse.OptionParser(usage=
"%prog input.pdb [-o output] [-m model_num]")
parser.add_option("-o", dest="outfile", default="distancemap.png",
                  help="Output file name (default: distancemap.png)")
parser.add_option("-m", dest="model_num", default=0, type="int",
                  help="Model number from pdb file to be used (default: 0)")
options, args = parser.parse_args()

if len(args) != 1:
    parser.print_usage()
    sys.exit(1)

from DistanceMap import DistanceMap
from utils import get_structure 

structure = get_structure(args[0])
model = structure[options.model_num]

map = DistanceMap(model)
with open (options.outfile+".tsv", "w") as f:
    f.write(map.tabsep_data())
f.close()

map.print_figure(options.outfile + ".png")






