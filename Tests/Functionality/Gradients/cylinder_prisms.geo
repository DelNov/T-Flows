SetFactory("OpenCASCADE");

A =  2.0;
N = 20;
gridsize = A / N;

Circle(1) = {0, 0, 0, A/2, 0, 2*Pi};
Curve Loop(1) = {1};
Plane Surface(1) = {1};

surfaceVector[] = Extrude {0, 0, A} {
 Surface{1};
 Layers{N};

 Recombine;
};

Characteristic Length {2, 1} = gridsize;

Physical Surface("bottom") = {1};
Physical Surface("top") = {2};
Physical Surface("side") = {3};
Physical Volume("volume") = {1};
