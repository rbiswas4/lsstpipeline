#!/usr/bin/env bash

mydir=`pwd`
export OPSIMOBS=/home/rbiswas/src/LSST/opsimobs
export throughputsdir=/home/soft/throughputs-1.2
cd ${throughputsdir} 
source ups/throughputs.sh
cd ${mydir}
make strategy2opsimout
make opsim2simlib 

