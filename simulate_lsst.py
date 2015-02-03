import sncosmo
import numpy as np
import os
from numpy.random import normal, uniform
from astropy.table import Table, vstack
from astropy.units import Unit
from sncosmo import Model
import matplotlib.pyplot as plt

def getlsstbandpassobjs(plot=False):
    """
    General utility to return a list of the baseline LSST bandpasses
    and register them as SNCosmo bandpasses accessible through strings 
    like 'LSSTu'.
    
    Parameters
    ----------
    plot : bool, optional, defaults to False
        plot filter functions obtained
    Returns
    -------
    None
    
    Examples
    --------
    >>> getlsstbandpassobjs() 

    """
    bandPassList = ['u', 'g', 'r', 'i', 'z', 'y']
    banddir = os.path.join(os.getenv('THROUGHPUTS_DIR'), 'baseline')
    lsstbands = []
    lsstbp = {}

    for band in bandPassList:
        # setup sncosmo bandpasses
        bandfname = banddir + "/total_" + band + '.dat'
        loadsncosmo = True
        if loadsncosmo:
            # register the LSST bands to the SNCosmo registry
            # Not needed for LSST, but useful to compare independent codes
            # Usually the next two lines can be merged,
            # but there is an astropy bug currently which affects only OSX.
            numpyband = np.loadtxt(bandfname)
            sncosmoband = sncosmo.Bandpass(wave=numpyband[:, 0],
                                           trans=numpyband[:, 1],
                                           wave_unit=Unit('nm'),
                                           name='LSST_' + band)
            sncosmo.registry.register(sncosmoband, force=True)
    if plot:
        filterfigs, filterax = plt.subplots()
        for band in bandPassList:
            b = sncosmo.get_bandpass('LSST_' + band)
            filterax.plot(b.wave, b.trans, '-k', lw=2.0)
    return None

def _prefixbands(prefix, obstable):
    _bb = np.array(obstable['FLT'])
    _lsst = np.array([prefix]*len(_bb), dtype='S6')
    return map ( ''.join, zip(_lsst,_bb))


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
    col = Table.Column(_prefixbandname("LSST_", obstable), name='band')
    obstable.add_column(col)       

    return None


def simulate_simlib(simlibfile, snmodelsource, outfile):
    """
    Simulate SN based on the simlibfile using SNCosmo SALT models


    Parameters
    ----------
    simlibfile :

    snmodelsource:

    outfile:

    Returns
    -------

    """
    from astropy.units import Unit
    from astropy.coordinates import SkyCoord

    # read the simlibfile into obstables
    meta, obstables = sncosmo.read_snana_simlib(simlibfilename)

    # set the SNCosmo model source
    dustmaproot = os.getenv('SIMS_DUSTMAPS_DIR')
    map_dir = os.path.join(dustmaproot, 'DustMaps')
    dust = sncosmo.CCM89Dust()
    model = Model(source="salt2-extended",
                  effects=[dust, dust],
                  effect_frames=['rest', 'obs'],
                  effect_names=['host', 'mw'])

    # Different fields in SIMLIB are indexed by libids
    libids = obstables.keys()
    lcs = []  
    for libid in libids:

        # Get the obstable corresponding to each field
        obstable = obstables[libid]
        manipulateObsTable(ObsTable)
        # Need Area from PixSize
        ra  =  obstable.meta['RA']
        dec  =  obstable.meta['DECL']
        skycoords = SkyCoord(ra, dec, unit='deg')
        t_mwebv = sncosmo.get_ebv_from_map(skycoords, mapdir=map_dir,
                                           interpolate=False)
        model.set(mwebv=t_mwebv) 
        params = []
        #col = Table.Column(obstable['SEARCH'].size*['ab'], name='zpsys')
        #obstable['FLT'].name =  'band'
        #obstable['MJD'].name =  'time'
        #obstable['ZPTAVG'].name =  'zp'
        #obstable['CCD_GAIN'].name =  'gain'
        #obstable['SKYSIG'].name =  'skynoise'
        #obstable.add_column(col)
        redshifts = list(sncosmo.zdist(0., 1.2, ratefunc=cosmoRate, time = rangemjd, area=1.))
        print 'num SN generated ', len(redshifts)
	for z in redshifts:
	    mabs = normal(-19.3, 0.3)
	    model.set(z=z)
	    model.set_source_peakabsmag(mabs, 'bessellb', 'ab')
	    x0 = model.get('x0')
            # RB: really should not be min, max but done like in catalogs
	    p = {'z':z, 't0':uniform(minmjd, maxmjd), 'x0':x0, 'x1': normal(0., 1.), 'c': normal(0., 0.1)}
	    params.append(p)
        print 'realizing SN'
        lcslib  =  sncosmo.realize_lcs(obstable, model, params)
        lcs.append(lcslib)
    alllcsintables = vstack(lcs) 
    print alllcsintables[0]
    print alllcsintables[MJD].size
    sncosmo.write_lc(alllcsintables, fname='simulatedlc.dat', format='ascii')


def cosmoRate(z, alpha=2.6e-5, beta=1.5, H0=70):
    """
    Returns the rate of SN in units of number of SN per Mpc^3 comoving volume
    per year for the redshift(s) z according to power law relations 
    alpha ( 1 + z) ** beta usually used in the literature.


    Parameters
    ---------
    z : float or array of floats, mandatory
        redshifts at which 
    alpha : float, optional, defaults to 2.65e-5
    beta : float, optional, defaults to 1.5 
    H0 : float optional, defaults to 70 
        Hubble constant in units of Km/s/Mpc
       
    Returns
    -------

    Examples
    --------
    """
    z = np.asarray (z)
    return alpha * (70.0 / H0 )* (1.0 + z)**beta
              


getlsstbandpassobjs()
simlibfilename = 'cache/RBtest_opsim.simlib'
#simulate_simlib(simlibfilename)
