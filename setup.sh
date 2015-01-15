#!/usr/bin/env bash

mydir=`pwd`
export OPSIMOBS=${HOME}/src/OpsimObs
export THROUGHPUTS_DIR=${HOME}/data/throughputs-1.2
#cd ${throughputsdir} 
#source ups/throughputs.sh
#cd ${mydir}
make strategy2opsimout
make opsim2simlib 

