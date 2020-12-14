SetFactory("OpenCASCADE");

R_VESSEL = 2;
R_PIPES  = 0.3;

//-----------------
// Define geometry
//-----------------
Cylinder(1) = {   0, 0, 0,    0,  0,  1,  R_VESSEL, 2*Pi};
Cylinder(2) = {-4.5, 0, 0.5,  3,  0,  0,  R_PIPES,  2*Pi};
Cylinder(3) = { 4.5, 0, 0.5, -3,  0,  0,  R_PIPES,  2*Pi};

// Unite them all
BooleanUnion{ Volume{1}; Delete; }{ Volume{3}; Volume{2}; Delete; }

//---------------------------------------
// Define boundary and volume conditions
//---------------------------------------
Physical Surface("inlet") = {7};
Physical Surface("outlet") = {6};
Physical Surface("bottom") = {3};
Physical Surface("top") = {4};
Physical Surface("vessel_walls") = {1};
Physical Surface("pipe_walls") = {2, 5};
Physical Volume("vessel") = {1};

//-----------------------------------------------------
// Characteristic mesh size (it is defined per points)
//-----------------------------------------------------
Characteristic Length {7, 5, 3, 1, 6} = 0.03;
Characteristic Length {2, 4} = 0.06;
