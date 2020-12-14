// Specify mesh size
DELTA = 0.06;

// Points at the front
Point(1) = {-0.7071,  0,     0,  DELTA};
Point(2) = { 0.0,     0.5,   0,  DELTA};
Point(3) = { 0.0,     0.25,  0,  DELTA};
Point(4) = { 0.75,    0.25,  0,  DELTA};
Point(5) = { 0.5,     0.0,   0,  DELTA};

// Lines at the front
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 1};

// Surface at the front
Curve Loop(1) = {1, 2, 3, 4, 5};  Plane Surface(1) = {1};

// Create volume
Extrude {0, 0, 0.5} {
  Surface{1}; 
}

// Boundary conditions
Physical Surface("head") = {15};
Physical Surface("tail") = {27};
Physical Surface("symmetry") = {31};
Physical Surface("neck") = {19};
Physical Surface("body") = {23};
Physical Surface("periodic") = {32, 1};

// Volume condition
Physical Volume("arrow") = {1};
