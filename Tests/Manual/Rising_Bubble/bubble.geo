A  =   1.0;  // length and width of the pool
B  =   2.0;  // height of the pool
NA =  64;    // resolution in length and height
NB = 128;    // resolution in width (periodic direction)

// Define points
Point(1) = {0, 0, 0};  Point(2) = {A, 0, 0};
Point(3) = {A, A, 0};  Point(4) = {0, A, 0};

// Define lines
Line(1) = {1, 2};  Line(2) = {2, 3};
Line(3) = {3, 4};  Line(4) = {4, 1};

// Define front surface
Curve Loop(1) = {1, 2, 3, 4};  Plane Surface(1) = {1};

// Set lines to transfinite
// (+1 is because GMSH expects number of points here)
Transfinite Curve {1, 2, 3, 4} = NA+1  Using Bump  1.0;

// Set surface to be transfinite and quadrilateral
Transfinite Surface {1} = {1, 2, 3, 4};  Recombine Surface {1};

// Expand in spanwise direction
Extrude {0, 0, B} {Surface{1}; Layers {NB}; Recombine;}

// Define boundary conditions
Physical Surface("pool_walls", 28) = {21, 25, 17, 13, 1, 26};

// Define volume condition
Physical Volume("fluid", 31) = {1};
