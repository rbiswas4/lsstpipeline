Prerequisites:
==============

- gfortran, f2py, standard python packages
 - anaconda python: Things should work with other python 2.7.5 versions, but we 
   often work with anaconda python. So, this is well tested on anaconda python.
 - throughputs for LSST. You need the path to the throughputs directory:: 
  
    wget http://dev.lsstcorp.org/cgit/LSST/sims/throughputs.git/snapshot/throughputs-1.2.tar.gz
    tar -xzvf throughputs-1.2.tar.gz 
    cd throughputs-1.2
 - anaconda has astropy, but you may need to update it::

    conda update astropy
    pip install emcee triangle iminuit

 - SNcosmo (install following instructions, working with developer version: https://github.com/sncosmo/sncosmo)::

    git clone https://github.com/sncosmo/sncosmo.git 
    cd sncosmo
    python setup.py install --user
 
