SetFactory("OpenCASCADE");

R_VESSEL = 2;
R_PIPES  = 0.3;

Cylinder(1) = {   0, 0, 0,    0,  0,  1,  R_VESSEL, 2*Pi};
Cylinder(2) = {-4.5, 0, 0.5,  3,  0,  0,  R_PIPES,  2*Pi};
Cylinder(3) = { 4.5, 0, 0.5, -3,  0,  0,  R_PIPES,  2*Pi};

// Unite them all
BooleanUnion{ Volume{1}; Delete; }{ Volume{3}; Volume{2}; Delete; }

//---------------------------------------
// Define boundary and volume conditions
//---------------------------------------
Physical Surface("bottom") = {3};
Physical Surface("top") = {4};
Physical Surface("walls") = {1, 5, 2};
Physical Surface("inlet") = {7};
Physical Surface("outlet") = {6};
Physical Volume("domain") = {1};

//--------------------------
// Characteristic mesh size
//--------------------------
Characteristic Length {7, 5, 1, 6, 2, 3, 4} = 0.04;

