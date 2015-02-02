#! /usr/bin/python

"""
This is an initial driver program for the LSST SN Simulation pipeline. It is
a direct transliteration from the ipython notebook.
"""

import sncosmo
import simulate_lsst as sl
import numpy as np
import pylab as pl
from numpy.random import normal, uniform
from astropy.table import Table

print "sncosmo version : ", sncosmo.__version__

"""
The basis of the simulation is a "SIMLIB" file which is a list of observations 
in a particular field of the sky. This simlib file is therefore created from 
either the OPSIMs output (which is a database of observations made over the 
course of the survey) or from OpsimObs (which is a time ordered set of 
observations) made across the survey for Deep Drilling fields for a 
particular strategy.
"""

simlibfilename = '/home/kushner/src/lsstpipeline/cache/RBtest_opsim.simlib'
#simlibfilename = 'cache/RBtest_opsim.simlib'

# set up bandpasses for LSST
sl.getlsstbandpassobjs(plot=False)

def prefixbandname(prefix, obstable):
        
    _bb = np.array(obstable['FLT'])
    _lsst = np.array([prefix]*len(_bb), dtype='S6')
    return map ( ''.join, zip(_lsst,_bb))
#Read the simlib file, get the libids which index different fields
meta, obstables = sncosmo.read_snana_simlib(simlibfilename)
print "META\n"
print meta
print "\n\n OBSTABLE LIBIDS \n"
print obstables.keys()

#Choose a field: for example the one indexed by 519
minmjd = obstables[519]['MJD'].min()
maxmjd = obstables[519]['MJD'].max()
obstable = obstables[519]
print "Number of observations in obstable: ", len(obstable)
print "min and max MJD in the set: ", minmjd, maxmjd

#Set SNCosmo model and set its t0 aritificially to values I know are in the above table
model = sncosmo.Model(source='salt2-extended')
params = []

# For a longer/complete answer, we will not set the z by hand, we will 
# draw from a distribution, but we are just doing one SN. So this is a
# good way

numSN  = 2
zlow = 0.
zhigh = 1.2
zvals =  uniform(zlow, zhigh, size=numSN)

for z in zvals: 
    mabs = normal(-19.3, 0.3)
    print 'mabs', mabs
    model.set(z=z)
    model.set_source_peakabsmag(mabs, 'bessellb', 'ab')
    
    # RB: really should not be min, max but done like in catalogs
    p = {'z':z, 't0': uniform(minmjd, maxmjd), 'x1': normal(0., 1.), 'c': normal(0., 0.1)}
    x0 = model.get('x0')
    p['x0']=x0
    p['t0']= 50000.
    params.append(p)
    
params[1]['t0'] = 50100.
print params[0]['t0']
print params[1]['t0']

print model # This is the second model

# Manipulate Obstable to form we can use for SNCosmo realize_lcs
# Should make a function using this for future
def manipulateObsTable(obstable):
    """
    Manipulate obstable read in from a SNANA style SIMLIB file into required
    columns for SNCosmo to work. This is an inplace modification leading to 
    no returns, but the input table itself is modified


    Parameters
    ----------
    obstable: `~astropy.Table`, mandatory
    Table of observations corresponding to a single field having the columns
    'MJD', 'ZPTAVG', 'CCD_GAIN' 'SKYSKIG', 'FLT':
    """
    col = Table.Column(obstable['SEARCH'].size*['ab'], name='zpsys')
    obstable.add_column(col)
    obstable['MJD'].name =  'time'
    obstable['ZPTAVG'].name =  'zp'
    obstable['CCD_GAIN'].name =  'gain'
    obstable['SKYSIG'].name =  'skynoise'
    col = Table.Column(prefixbandname("LSST_", obstable), name='band')
    obstable.add_column(col)       

    return None

manipulateObsTable(obstable)
# Now display the table of observations
print str(obstable)

## Simulate the SN, plot and profile

relevantdata = sncosmo.realize_lcs(obstable, model, params, trim_observations=True)

from astropy.table import Table
print type(relevantdata)

x = Table(relevantdata[0])

print str(x)

print str(len(relevantdata[0]))

# Write out to a file as this is how we will do things in a larger set by looping through, and then read in the file

model.set(**params[0])

# fig_relevant = sncosmo.plot_lc(relevantdata[0], model=model)

#print "Close Window to continute"
#pl.show()
sncosmo.write_lc(Table(relevantdata[0]), fname='lc.dat', format='ascii')

# sncosmo.write_lc(Table(relevantdata[0]), fname='lc.dat.json', format='json')

# sncosmo.write_lc(Table(relevantdata[0]), fname='lc.dat.fits', format='snana') # fits

lc = sncosmo.read_lc('lc.dat', format='ascii')

fmodel = sncosmo.Model(source='salt2-extended')
for z in zvals:
    fmodel.set(z=z)
    res, fitmodel = sncosmo.fit_lc(relevantdata[0], fmodel, ['t0', 'x0', 'x1', 'c']) 

print res

print params[0]


