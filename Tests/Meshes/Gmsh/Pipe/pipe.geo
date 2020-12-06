SetFactory("OpenCASCADE");

Circle(1) = {0, 0, 0, 0.5, 0, 2*Pi};
Curve Loop(1) = {1};
Plane Surface(1) = {1};
Transfinite Curve {1} = 60 Using Progression 1;

Extrude {0, 0, 8} {
  Surface{1}; Layers{120};
  Recombine;
}

Physical Surface("bottom") = {1};
Physical Surface("top") = {3};
Physical Surface("side_wall") = {2};
Physical Volume("vessel") = {1};
