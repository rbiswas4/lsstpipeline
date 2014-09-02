#!/usr/bin/env python

#class SimLSST(strategyfile, simname):
import strategy2opsimin
import os
import opsimObs.opsimObs as op
import opsim2simlib




#use environment variable setup by setup.sh to obtain the directory of opsims
opsimdir = os.environ.get("OPSIMOBS")


programdir = os.path.dirname(os.path.realpath(__file__))
mydir = os.getcwd() + "/"


#Get INPUT required for OPSIMOBS from strategy
def strategy2simlib(strategyfile , dataprefix, cache, simlibfile=None):


	#if simlibfile is None:
	#	simlibfile = dataprefix + "_opsim.simlib"
	opsiminfile  = dataprefix + "_" + "opsiminput"
	strategy2opsimin.lsstdfacg(act="I", dfacgf=strategyfile, outfile=cache + opsiminfile)

	#TODO: Even better change opsimObs, so that it need not be run from the 
	#directory
	os.chdir(opsimdir)
	op.getopsimouts( cache + opsiminfile,  "opsimSky.conf", 
		outfile=cache + dataprefix+ "_opsimout.dat")
	os.chdir(programdir)
	#Get simlib from output
	opsim2simlib.lsstsimlib(opsimf = cache + dataprefix +"_opsimout.dat", simlib = cache + simlibfile)
