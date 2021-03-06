Getting simulations from OPSIMOBS:
=================================
This README explains how to set up the pipeline and run it for a particular 
strategy. The strategy is defined by a control file as used in David Cinabro's 
Deep Field Arbitrary Cadence Generator. 

Prerequisites:
==============
 - gfortran, f2py, standard python packages
 - anaconda python: Things should work with other python 2.7.5 versions, but we 
   often work with anaconda python. So, this is well tested on anaconda python.
 - pyslalib: clone from https://github.com/scottransom/pyslalib and install 
 - throughputs (just wget http://dev.lsstcorp.org/cgit/LSST/sims/throughputs.git/snapshot/throughputs-1.2.tar.gz and     unatar it. You will require the path to this directory).  
 - opsimobs (forked from Lynne Jones' https://github.com/rhiannonlynne/OpsimObs , slightly modified at present, clone from https://github.com/rbiswas4/OpsimObs) 
 - anaconda has astropy, but you may need to update it::

    conda update astropy
    pip install emcee triangle iminuit
 - SNcosmo (install following instructions, working with developer version: https://github.com/sncosmo/sncosmo)::

    git clone https://github.com/sncosmo/sncosmo.git 
    cd sncosmo
    python setup.py install --user
 - SNANA (install following instructions http://das.sdss2.org/ge/sample/sdsssn/SNANA-PUBLIC/)
 
Setting up:
==========
 - Edit setup.sh to point to the correct directories 
 - Edit Makefiles if you are using other compilers. With gfortran, you should not have to edit this. Used to 
	the source files in the fortran directory which are essentially the files from David Cinabro's 
	Deep Field Arbitrary Cadence Generator: 
	(https://github.com/DarkEnergyScienceCollaboration/DeepFieldArbitraryCadenceGenerator)  
 - source setup.sh (assuming you are on bash, this step sets up the environment variables required by opsimobs, 
	compiles the fortran code into a library of functions that can be called from python

Running:
=======

 - python lsstpipeline --inputfile inputfilename
