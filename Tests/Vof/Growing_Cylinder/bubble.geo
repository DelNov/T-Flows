A  = 1.0;  // length and height of the cavity
B  = 0.2;  // width of the cavity
NA = 40;    // resolution in length and height
NB = 20;    // resolution in width (periodic direction)

// Define points
Point(1) = {-0.5*A, 0, -0.5*A};  Point(2) = {A*0.5, 0, -0.5*A};
Point(3) = {A*0.5, 0, A*0.5};  Point(4) = {-0.5*A, 0, 0.5*A};

// Define lines
Line(1) = {1, 2};  Line(2) = {2, 3};
Line(3) = {3, 4};  Line(4) = {4, 1};

// Define front surface
Curve Loop(1) = {1, 2, 3, 4};  Plane Surface(1) = {1};

// Set lines to transfinite
// (+1 is because GMSH expects number of points here)
Transfinite Curve {1, 2, 3, 4} = NA+1  Using Bump  1.0;

// Set surface to be transfinite and quadrilateral (comment for tetrahedrals)
Transfinite Surface {1} = {1, 2, 3, 4};  Recombine Surface {1};

// Expand in spanwise direction
Extrude {0, B, 0} {Surface{1}; Layers {NB}; Recombine;}

// Define boundary conditions
Physical Surface("upper_wall", 27) = {21};
Physical Surface("side_walls", 28) = {25, 17};
Physical Surface("lower_wall", 29) = {13};
Physical Surface("periodic_y", 30) = {1, 26};

// Define volume condition
Physical Volume("fluid", 31) = {1};
