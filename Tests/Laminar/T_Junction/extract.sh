#!/bin/bash
OUT="out_08_with_particles"
cat $OUT | grep ' 1 position'| awk '{print $3 " " $5}' > particle_01_test.dat
cat $OUT | grep ' 2 position'| awk '{print $3 " " $5}' > particle_02_test.dat
cat $OUT | grep ' 3 position'| awk '{print $3 " " $5}' > particle_03_test.dat
cat $OUT | grep ' 4 position'| awk '{print $3 " " $5}' > particle_04_test.dat
cat $OUT | grep ' 5 position'| awk '{print $3 " " $5}' > particle_05_test.dat
cat $OUT | grep ' 6 position'| awk '{print $3 " " $5}' > particle_06_test.dat
cat $OUT | grep ' 7 position'| awk '{print $3 " " $5}' > particle_07_test.dat
cat $OUT | grep ' 8 position'| awk '{print $3 " " $5}' > particle_08_test.dat
cat $OUT | grep ' 9 position'| awk '{print $3 " " $5}' > particle_09_test.dat
cat $OUT | grep '10 position'| awk '{print $3 " " $5}' > particle_10_test.dat
cat $OUT | grep '11 position'| awk '{print $3 " " $5}' > particle_11_test.dat
cat $OUT | grep '12 position'| awk '{print $3 " " $5}' > particle_12_test.dat
cat $OUT | grep '13 position'| awk '{print $3 " " $5}' > particle_13_test.dat
cat $OUT | grep '14 position'| awk '{print $3 " " $5}' > particle_14_test.dat
cat $OUT | grep '15 position'| awk '{print $3 " " $5}' > particle_15_test.dat
cat $OUT | grep '16 position'| awk '{print $3 " " $5}' > particle_16_test.dat
