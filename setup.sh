#!/usr/bin/env bash

mydir=`pwd`
export OPSIMOBS=${HOME}/src/OpsimObs
export throughputsdir=${THROUGHPUTSDIR}
#cd ${throughputsdir} 
#source ups/throughputs.sh
#cd ${mydir}
make strategy2opsimout
make opsim2simlib 

