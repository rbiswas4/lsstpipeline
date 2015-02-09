#!/usr/bin/env python

"""
build the z distribution of light curves in directory

"""
import numpy as np
import sncosmo
import glob
from astropy.cosmology import Planck13 as cosmo
from astropy.utils import lazyproperty

def buildzdist(dirname):

    lclist = glob.glob(dirname + '/*')

    z = [] 
    for lcfile in lclist:
        lc = sncosmo.read_lc(lcfile)
        zval  = lc.meta['z']
        z.append(zval)
    return z

class lc(object):

    def __init__(self, lc) :
        """
        Parameters
        ----------
        lc : `~astropy.Table`

        """
        self.lc = lc
        self.data = np.asarray(lc)

    @lazyproperty 
    def bandlcs(self):
        """
        group a light curve by bands and obtain the snr for each of these
        """
        lc = self.lc
        bandgroups = lc.group_by('band')
        return bandgroups

    @lazyproperty
    def usedbands(self):
        """
        array of band names that have measurements 
        """
        lc = self.lc
        return np.unique(np.asarray(lc['band']))

    def bandmaxsnr(self):
        """

        """
        bandgroups = self.bandlcs
        return map(lambda xx: max(xx['flux'] / xx['fluxerr']), bandgroups.groups)







def getbandlcs(lc):
    """
    group a light curve by bands and obtain the snr for each of these
    """
    bandgroups = lc.group_by('band')

    return bandgroups


if __name__ == '__main__':
    z = buildzdist('LC')
    import seaborn as sns
    sns.set()
    sns.set_style('darkgrid')
    import matplotlib.pyplot as plt
    plt.hist(z, bins = 12, histtype='step', lw=2.0)
    plt.show()

