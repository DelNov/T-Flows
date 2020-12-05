//------------------------------------------
// Flow is in Y direction, periodicity in X
//------------------------------------------

// Constants
SURF_BOTTOM = 101;
SURF_TOP    = 102;
SURF_WEST   = 103;
SURF_EAST   = 104;
SURF_SOUTH  = 105;
SURF_NORTH  = 106;

//--------
// Points
//--------
 Point(1) = {0.0, 0.0, -0.20};
 Point(2) = {1.0, 0.0, -0.20};
 Point(3) = {0.0, 1.0, -0.20};
 Point(4) = {1.0, 1.0, -0.20};
 Point(5) = {0.0, 0.0, -0.00};
 Point(6) = {1.0, 0.0, -0.00};
 Point(7) = {0.0, 1.0, -0.00};
 Point(8) = {1.0, 1.0, -0.00};

//-------
// Lines
//-------
 Line(1) = {1, 2};
 Line(2) = {2, 4};
 Line(3) = {4, 3};
 Line(4) = {3, 1};
 Line(5) = {5, 6};
 Line(6) = {6, 8};
 Line(7) = {8, 7};
 Line(8) = {7, 5};
 Line(9) = {1, 5};
 Line(10) = {2, 6};
 Line(11) = {3, 7};
 Line(12) = {4, 8};

//----------
// Surfaces
//----------
Line Loop(13) = {1, 2, 3, 4};
Plane Surface(SURF_BOTTOM) = {13};
Line Loop(15) = {5, 6, 7, 8};
Plane Surface(SURF_TOP) = {15};
Line Loop(17) = {9, 5, -10, -1};
Plane Surface(SURF_WEST) = {17};
Line Loop(19) = {11, -7, -12, 3};
Plane Surface(SURF_EAST) = {19};
Line Loop(21) = {9, -8, -11, 4};
Plane Surface(SURF_SOUTH) = {21};
Line Loop(23) = {10, 6, -12, -2};
Plane Surface(SURF_NORTH) = {23};

//--------
// Volume
//--------
Surface Loop(25) = {SURF_BOTTOM, SURF_TOP, SURF_WEST, SURF_EAST, SURF_SOUTH, SURF_NORTH};
Volume(26) = {25};

//---------------------
// Boundary conditions
//---------------------
Physical Surface("bottom_wall") = {SURF_BOTTOM};
Physical Surface("top_wall")    = {SURF_TOP};
Physical Surface("inlet")       = {SURF_WEST};
Physical Surface("outlet")      = {SURF_EAST};
Physical Surface("periodic_x")  = {SURF_SOUTH, SURF_NORTH};

Physical Volume("lower_dom") = {26};

//------
// Mesh
//------
Transfinite Line {1, 3, 5, 7}    = 41;
Transfinite Line {9, 10, 11, 12} = 41;
Transfinite Line {2, 4, 8, 6}    = 41;
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";
