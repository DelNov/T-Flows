// Gmsh project created on Sat May  1 11:02:16 2021
SetFactory("OpenCASCADE");

Cylinder(1) = {0, 0, 0, 0, 0, 2, 1.0, 2*Pi};

Physical Surface("bottom") = {3};
Physical Surface("top") = {2};
Physical Surface("side") = {1};

Physical Volume("volume") = {1};

Characteristic Length {1, 2} = 0.1;
