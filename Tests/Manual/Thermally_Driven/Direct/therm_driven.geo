A    =  1.0;  // length and height of the cavity
B    =  0.1;  // width of the cavity
NA   = 60;    // resolution in length and height
NB   =  3;    // resolution in width (periodic direction)
BUMP = 0.1;

// Define points
Point(1) = {0, 0, 0};  Point(2) = {A, 0, 0};
Point(3) = {A, 0, A};  Point(4) = {0, 0, A};

// Define lines
Line(1) = {1, 2};  Line(2) = {2, 3};
Line(3) = {3, 4};  Line(4) = {4, 1};

// Define front surface
Curve Loop(1) = {1, 2, 3, 4};  Plane Surface(1) = {1};

// Set lines to transfinite
// (+1 is because GMSH expects number of points here)
Transfinite Curve {1, 2, 3, 4} = NA+1  Using Bump BUMP;

// Set surface to be transfinite and quadrilateral
Transfinite Surface {1} = {1, 2, 3, 4};  Recombine Surface {1};

// Expand in spanwise direction
Extrude {0, B, 0} {Surface{1}; Layers {NB}; Recombine;}

// Define boundary conditions
Physical Surface("upper_wall", 27) = {21};
Physical Surface("left_wall", 28) = {17};
Physical Surface("right_wall", 29) = {25};
Physical Surface("lower_wall", 30) = {13};
Physical Surface("periodic_y", 31) = {1, 26};

// Define volume condition
Physical Volume("fluid", 32) = {1};
