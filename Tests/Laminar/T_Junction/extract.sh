#!/bin/bash

OUT=$1

cat $OUT | grep 'particle:   1' | awk '{print $7 " " $8}' > particle_01_mpi.dat
cat $OUT | grep 'particle:   2' | awk '{print $7 " " $8}' > particle_02_mpi.dat
cat $OUT | grep 'particle:   3' | awk '{print $7 " " $8}' > particle_03_mpi.dat
cat $OUT | grep 'particle:   4' | awk '{print $7 " " $8}' > particle_04_mpi.dat
cat $OUT | grep 'particle:   5' | awk '{print $7 " " $8}' > particle_05_mpi.dat
cat $OUT | grep 'particle:   6' | awk '{print $7 " " $8}' > particle_06_mpi.dat
cat $OUT | grep 'particle:   7' | awk '{print $7 " " $8}' > particle_07_mpi.dat
cat $OUT | grep 'particle:   8' | awk '{print $7 " " $8}' > particle_08_mpi.dat
cat $OUT | grep 'particle:   9' | awk '{print $7 " " $8}' > particle_09_mpi.dat
cat $OUT | grep 'particle:  10' | awk '{print $7 " " $8}' > particle_10_mpi.dat
cat $OUT | grep 'particle:  11' | awk '{print $7 " " $8}' > particle_11_mpi.dat
cat $OUT | grep 'particle:  12' | awk '{print $7 " " $8}' > particle_12_mpi.dat
cat $OUT | grep 'particle:  13' | awk '{print $7 " " $8}' > particle_13_mpi.dat
cat $OUT | grep 'particle:  14' | awk '{print $7 " " $8}' > particle_14_mpi.dat
cat $OUT | grep 'particle:  15' | awk '{print $7 " " $8}' > particle_15_mpi.dat
cat $OUT | grep 'particle:  16' | awk '{print $7 " " $8}' > particle_16_mpi.dat
