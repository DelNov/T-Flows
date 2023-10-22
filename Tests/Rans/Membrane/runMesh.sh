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

# Visualize the mesh in Paraview 
paraview upper.vtu membrane.vtu lower.vtu &  
