# -*- coding: iso-8859-15 -*-


import sys, time, os
import pymol

pdbId = sys.argv[1]
chain = sys.argv[2]
units = sys.argv[3].split(";")
bound = sys.argv[4]


print pdbId
print chain

