bubbly flow simulation with mass transfer in a subchannel

How to perform simulation
1. gmsh -3 fuel_bundle.geo > gmsh.out
2. Convert < convert.scr > convert.out
3. Divide 
4. InitialBubbles/main.f90 -> ellipsoid_parameters.ini
5. mpiexec -np 16 ./Process < /dev/null > out &
      output: volume.out, massTransfer.out
6. Post-process
      post.sh

Note:
In control file, coefficients for Lee model are defined.
c_lee = 10000.0 for condensation
        20000.0 for vaporization

