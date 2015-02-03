#!/usr/bin/env python

"""
build the z distribution of light curves in directory

"""
import sncosmo
import glob
def buildzdist(dirname):

    lclist = glob.glob(dirname + '/*')

    z = [] 
    for lcfile in lclist:
        lc = sncosmo.read_lc(lcfile)
        zval  = lc.meta['z']
        z.append(zval)
    return z

z = buildzdist('LC')
import matplotlib.pyplot as plt
plt.hist(z, bins = 12, histtype='step')
plt.show()

