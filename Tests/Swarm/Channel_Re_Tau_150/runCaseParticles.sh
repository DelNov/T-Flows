#!/bin/bash

# Specify the number of processors here (for parallel mode): 
NUMBEROFCPUS=4

# Creating a soft link for the fluid flow control file 
ln -i -s control_b_particles

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
    ./Process >| out_particleResults_serialMode.dat &
    tail -100f out_particleResults_serialMode.dat
elif [ $RUN = "parallel" ]
then
    # Decomposing the domain
    ./Divide < divide.scr
    echo "#-------------------------------------#"
    echo "# Running T-Flows in parallel mode... #"
    echo "#-------------------------------------#"
    mpirun -np $NUMBEROFCPUS ./Process >| out_particleResults_parallelMode.dat & 
    tail -100f out_particleResults_parallelMode.dat
else 
    echo "ERROR! Please specify one of the available run modes (i.e., 'serial' or 'parallel')"
    exit 1
fi
