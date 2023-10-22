#!/bin/bash

# Cleaning any old files 
rm -f out* *.faces *.monit readme
rm -f *.cfn *.dim *.pvtu *.vtu
rm -rf Sub*
rm -f Process
rm -f Generate
rm -f Divide
rm -f Convert

# Creating soft links from Binaries/ 
ln -i -s ../../../Binaries/* . 

# Generating mesh
./Generate < generate.scr

# Visualize the mesh in Paraview
paraview chan.vtu & 
