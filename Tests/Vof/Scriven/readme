Scriven problem: vapor bubble growth in a superheated liquid
   conditions of simulation:
      water at 1 bar
   initial conditions:
      degree of superheat = 1.25 K
      initial bubble radius = 50 microns, defined in ellipsoid_parameters.ini
      initial temperature profile = analytical solution defined in dT0.75_R50microns.txt
   computational domain and cell:
      domain size [-0.00015:0.00015] in the all directions = cube of 300 microns
      cell size = 300 microns / 48 = 6.25 microns
   analytical solution:
      Scriven-dT1.25-R50microns.nb (Mathematica)
      beta = 4.06022
      radius = 2 * beta * (alpha * time_analytical)

1. Generate < generate.scr > generate.out
2. Divide
3. Compile, check User_Mod
4. Process
   % mpiexec -np 4 Process </dev/null > out &
      input : dT1.25_R50microns.txt, ellipsoid_parameters.ini, control
      output: bench-data.dat
5. Post-process
   % gnuplot radius.gnu (warning! beta is hard coded)
      output radius.png (radius calculated from bubble volume)
   % tec360 iso-tpr.lay


