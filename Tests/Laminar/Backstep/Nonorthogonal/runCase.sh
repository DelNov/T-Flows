#!/bin/bash

# Specify the number of processors here (for parallel mode): 
NUMBEROFCPUS=4

# Cleaning old files 
rm -f out* *.faces *-monit* readme
rm -f *.cfn *.dim *.pvtu *.vtu
rm -rf Sub*
rm -f Process
rm -f Generate
rm -f Divide
rm -f Convert

# Creating soft links from Binaries/ 
ln -i -s ../../../../Binaries/* . 

# Generating mesh
./Generate < generate.scr

# running the solver (in either serial or parallel)
echo "#===========================================================#"
echo "#  Please specify the run mode i.e. 'serial' or 'parallel'  #"
echo "#===========================================================#"

# Read the input from user
read RUN

if [ $RUN = "serial" ]
then
    echo "#-----------------------------------#"
    echo "# Running T-Flows in serial mode... #"
    echo "#-----------------------------------#"
    ./Process >| out_results_serialMode.dat &
    tail -100f out_results_serialMode.dat
elif [ $RUN = "parallel" ]
then
    # Decomposing the domain
    ./Divide < divide.scr
    echo "#-------------------------------------#"
    echo "# Running T-Flows in parallel mode... #"
    echo "#-------------------------------------#"
    mpirun -np $NUMBEROFCPUS ./Process >| out_results_parallelMode.dat & 
    tail -100f out_results_parallelMode.dat
else 
    echo "ERROR! Please specify one of the available run modes (i.e., 'serial' or 'parallel')"
    exit 1
fi

