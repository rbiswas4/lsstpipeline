#!/usr/bin/env python

"""
simulate a set of SN from a simlib file

"""
import simulate_lsst as sl

# Read a simlib file and create a set of observation tables 

simlibfname = 'cache/RBtest_opsim.simlib'
meta, obstables = sncosmo.read_snana_simlib(simlibfilename)
def prefixbandname(prefix, obstable):

    _bb = np.array(obstable['FLT'])
    _lsst = np.array([prefix]*len(_bb), dtype='S6')
    return map ( ''.join, zip(_lsst,_bb))
# Get a distribution of redshifts from the rate equation


# simulate SN






