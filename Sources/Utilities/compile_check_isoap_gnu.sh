rm -f *.o
rm -f *.mod
gfortran -c -cpp -fdefault-real-8 ../Process/Polyhedron_Mod.f90 
gfortran -c -cpp -fdefault-real-8 ../Process/Iso_Polygons_Mod.f90 
gfortran -c -cpp -fdefault-real-8 ../Process/Isoap_Mod.f90 
gfortran -o Check_Isoap *.o -cpp -fdefault-real-8 Check_Isoap.f90
ln -i -s Isoap_Data/phi-cube-07  phi
