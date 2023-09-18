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

# Moving a possible readme file - lest overwriting something by mistake
mv readme README

# Creating soft links from Binaries/ 
ln -i -s ../../../../Binaries/* .

# Creating a soft link for the main control file
ln -i -s control.0 control 

# Constructing the mesh of both upper and lower domains (using gmsh)
echo "#=============================================#"
echo "#  Gmsh: Saving the mesh of the upper domain  #"
echo "#=============================================#"
gmsh  upper_dom.geo -3 -o upper_dom.msh
echo "#=============================================#"
echo "#  Gmsh: Saving the mesh of the lower domain  #"
echo "#=============================================#"
gmsh  lower_dom.geo -3 -o lower_dom.msh
echo "#  Mesh is successfully converted!  #"

# Converting the mesh in T-Flows  
echo "#=============================================#"
echo "#  Converting the mesh into the local format  #"
echo "#=============================================#"
./Convert < convert_upper_dom.scr
./Convert < convert_lower_dom.scr

# running the solver (in either serial or parallel)
echo "#===========================================================#"
echo "#  Please specify the run mode i.e. 'serial' or 'parallel'  #"
echo "#===========================================================#"

# Read the input from user
#read RUN

#if [ $RUN = "serial" ]
if [ $1 = serial ]
then
    echo "#-------------------------------------#"
    echo "#  Running T-Flows in serial mode...  #"
    echo "#-------------------------------------#"
    ./Process >| out_results_serialMode.dat &
    tail -100f out_results_serialMode.dat
#elif [ $RUN = "parallel" ]
elif [ $1 = parallel ]
then
    # Decomposing the domain
    ./Divide < divide_upper_dom.scr
    ./Divide < divide_lower_dom.scr
    echo "#---------------------------------------#"
    echo "#  Running T-Flows in parallel mode...  #"
    echo "#---------------------------------------#"
    mpirun -np $NUMBEROFCPUS ./Process >| out_results_parallelMode.dat & 
    tail -100f out_results_parallelMode.dat
else 
    echo "ERROR! Please specify one of the available run modes (i.e., 'serial' or 'parallel')"
    exit 1
fi

