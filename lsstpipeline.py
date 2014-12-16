#!/usr/bin/env python

#import strategy2opsimin
import os
import simlibio as sio
import OpsimObs.opsimObs as op
#import opsim2simlib
import utils 
import argparse


#use environment variable setup by setup.sh to obtain the directory of opsims
opsimdir = os.environ.get("OPSIMOBS")
mydir = os.getcwd() + "/"
programdir = os.path.dirname(os.path.realpath(__file__))

#Initialize a logstring
logstring = ""
logstring += utils.get_userinfo()


#parse command line arguments
parser = argparse.ArgumentParser(prog="lsstpipeline.py", 
                                description= """Runs the lsst pipeline based on several input files, but takes the jobs file as an argument """)
parser.add_argument("--inputfile", type=str, default=None, 
                    required=True, help="Absolute path to the input file for the pipeline")
args = parser.parse_args()
inputfile = args.inputfile

#Start 
#record on log file
logstring += __file__ + " on inputfile " + inputfile +" started. \n"
#Read input file and find out what to do
jobdict = utils.processfile(inputfile)

print jobdict
#Fill in some values from input file 
#cache
cache = jobdict["cache"] [0] + "/"
print cache
    # create the cache if it does not exist:
if not os.path.exists(cache):
    os.makedirs(cache)
    logstring += "cache created \n"

#Set the data file prefix
dataprefix = jobdict["runprefix"][0]
logstring += "found the data prefix " + dataprefix + "\n"

#get the simlib
simlibfilename = dataprefix +"_opsim.simlib"
strategyfilename = jobdict["createsimlib"][1]
if utils.str2bool(jobdict["createsimlib"][0]):
    logstring += "creating simlib file " + simlibfilename + " from stratgy file " + strategyfilename + "\n"
    sio.strategy2simlib(strategyfile=strategyfilename, cache=cache, dataprefix=dataprefix, simlibfile=simlibfilename) 
	
else:
    simlibfilename = strategyfilename
    logstring += "simlib file provided as "+ simlibfilename + "\n"
 
	
if jobdict['simulate_using'][0] == 'snana':
    utils.modifysnanafiles(cache, sndict)	
    os.system('sim_SNMix.pl Mix_sim.input')
elif jobdict['simulate_using'][0] == 'sncosmo':

print logstring
exit()
