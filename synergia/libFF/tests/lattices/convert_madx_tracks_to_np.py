#!/usr/bin/env python

import sys
import os

import numpy as np

madxoutput = sys.argv[1]
print "reading file: ", madxoutput

basename = os.path.splitext(madxoutput)[0]
print "writing numpy data to ", basename+".npy"

f = open(madxoutput)
lines = f.readlines()
nlines = len(lines)
f.close()

# format after header and first observation is
# number turn x px y py t pt s e
trackdata = np.loadtxt(madxoutput, skiprows=nlines-16, usecols=(2,3,4,5,6,7))
np.save(basename+".npy", trackdata)
