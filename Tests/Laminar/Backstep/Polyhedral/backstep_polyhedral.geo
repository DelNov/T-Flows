
// Points in front (x-z) plane
Point(1) = { 0,  0,  0,  0.1};
Point(2) = { 2,  0,  0,  0.1};
Point(3) = { 2,  0, -1,  0.1};
Point(4) = { 3,  0, -1,  0.1};
Point(5) = {13,  0, -1,  0.1};
Point(6) = {13,  0,  1,  0.1};
Point(7) = { 3,  0,  1,  0.1};
Point(8) = { 0,  0,  1,  0.1};

// Lines in front (x-z) plane
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};

// Front surface in x-z plane
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Plane Surface(1) = {1};

// Characterstic mesh sizes at points
Characteristic Length {1, 2, 3, 4, 7, 8} = 0.05;
Characteristic Length {5, 6} = 0.15;

// Extrude in z-direction
Extrude {0, 1, 0} {
  Surface{1}; Layers{5}; Recombine;
}

// Boundary (and volume) conditions
Physical Surface("in") = {49};
Physical Surface("out") = {37};
Physical Surface("low_wall") = {21, 29, 33};
Physical Surface("back_wall") = {25};
Physical Surface("top_wall") = {45, 41};
Physical Surface("periodic_y") = {1, 50};
Physical Volume("back") = {1};

