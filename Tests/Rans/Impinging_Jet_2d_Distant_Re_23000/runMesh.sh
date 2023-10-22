#!/bin/bash

# Specify the number of processors here (for parallel mode): 
NUMBEROFCPUS=4

# Cleaning any old files 
rm -f out* *.faces *-monit* readme
rm -f *.cfn *.dim *.pvtu *.vtu
rm -rf Sub*
rm -f Process
rm -f Generate
rm -f Divide
rm -f Convert

# Creating soft links from Binaries/ 
ln -i -s ../../../Binaries/* . 

# Extracting the mesh
gunzip jet.neu.gz 

# Converting the mesh
./Convert < convert.scr 

# Visualize the mesh in Paraview
paraview jet.vtu & 
