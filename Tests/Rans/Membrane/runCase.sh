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
ln -i -s ../../../Binaries/* .

# Creating a soft link for the main control file
ln -i -s control.0 control 

# Extracting the meshes of all three domains
echo "#=========================#"
echo "#  Extracting the meshes  #"
echo "#=========================#"
gunzip *.gz

# Converting the mesh in T-Flows  
echo "#=============================================#"
echo "#  Converting the mesh into the local format  #"
echo "#=============================================#"
./Convert < convert.1.scr
./Convert < convert.2.scr
./Convert < convert.3.scr

# running the solver (in either serial or parallel)
echo "#===========================================================#"
echo "#  Please specify the run mode i.e. 'serial' or 'parallel'  #"
echo "#===========================================================#"

# Read the input from user
#read RUN

#if [ $RUN = "serial" ]
if [ $1 = serial ]
then
    echo "#-----------------------------------#"
    echo "# Running T-Flows in serial mode... #"
    echo "#-----------------------------------#"
    ./Process >| out_results_serialMode.dat &
    tail -100f out_results_serialMode.dat
#elif [ $RUN = "parallel" ]
elif [ $1 = parallel ]
then
    # Decomposing the domain
    ./Divide < divide.1.scr
    ./Divide < divide.2.scr
    ./Divide < divide.3.scr
    echo "#-------------------------------------#"
    echo "# Running T-Flows in parallel mode... #"
    echo "#-------------------------------------#"
    mpirun -np $NUMBEROFCPUS ./Process >| out_results_parallelMode.dat & 
    tail -100f out_results_parallelMode.dat
else 
    echo "ERROR! Please specify one of the available run modes (i.e., 'serial' or 'parallel')"
    exit 1
fi

