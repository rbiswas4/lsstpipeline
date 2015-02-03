#!/usr/bin/env python

"""
simulate a set of SN from a simlib file

"""
import sncosmo
import simulate_lsst as sl

# Read a simlib file and create a set of observation tables 
#restrict to 10 SN from each field
lcs = sl.simulate_simlib(simlibfile='cache/snana_test.simlib', snmodelsource='salt2-extended', outfile='LC/lc', restrict=None)
# simulate SN






