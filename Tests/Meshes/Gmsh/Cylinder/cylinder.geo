SetFactory("OpenCASCADE");
Circle(1) = {0, 0, 0, 0.5, 0, 2*Pi};

Curve Loop(1) = {1};
Plane Surface(1) = {1};

Extrude {0, 0, 4} {
  Surface{1};
}

Transfinite Curve {3, 1} =  47 Using Progression 1;
Transfinite Curve {2}    =  60 Using Progression 1;
// Transfinite Curve {3, 1} = 22 Using Progression 1;
// Transfinite Curve {2}    =  7 Using Progression 1;

Physical Surface("bottom")    = {1};
Physical Surface("side_wall") = {2};
Physical Surface("top")       = {3};

Physical Volume("vessel") = {1};
