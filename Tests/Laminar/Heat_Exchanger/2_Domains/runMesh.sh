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

# Constructing the mesh of both upper and lower domains (using gmsh)
echo "#===========================================#"
echo "#  Gmsh: Saving the mesh of the upper domain  #"
echo "#===========================================#"
gmsh  upper_dom.geo -3 -o upper_dom.msh
echo "#===========================================#"
echo "#  Gmsh: Saving the mesh of the lower domain  #"
echo "#===========================================#"
gmsh  lower_dom.geo -3 -o lower_dom.msh
echo "#  Mesh is successfully converted!  #"

# Converting the mesh in T-Flows  
echo "#=============================================#"
echo "#  Converting the mesh into the local format  #"
echo "#=============================================#"
./Convert < convert_upper_dom.scr
./Convert < convert_lower_dom.scr

# Visualize the mesh in Paraview
paraview upper_dom.vtu lower_dom.vtu & 
