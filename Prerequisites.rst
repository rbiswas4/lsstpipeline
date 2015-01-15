Prerequisites:
==============

- gfortran, f2py
- Python packages: 
  We mostly use anaconda python which conveniently comes packaged with several 
  packages that we will need like astropy. A complete list of requirements are 
  in SNCosmo. 

- anaconda has astropy, but you may need to update it. Also it is useful to have the optional packages which are not required for all SNCosmo functions, but will be useful for us::

   conda update astropy
   pip install emcee triangle iminuit

- throughputs for LSST. You need the path to the throughputs directory:: 
  
   wget http://dev.lsstcorp.org/cgit/LSST/sims/throughputs.git/snapshot/throughputs-1.2.tar.gz
   tar -xzvf throughputs-1.2.tar.gz 
   cd throughputs-1.2


- SNcosmo (install following instructions, working with developer version: https://github.com/sncosmo/sncosmo)::

   git clone https://github.com/sncosmo/sncosmo.git
   cd sncosmo
   python setup.py install --user
