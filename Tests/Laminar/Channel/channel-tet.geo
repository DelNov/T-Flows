L = 12;
W =  1;
H =  1;

DEL = 0.1;

// Points
Point(1) = {0, 0, 0, DEL};
Point(2) = {L, 0, 0, DEL};
Point(3) = {L, W, 0, DEL};
Point(4) = {0, W, 0, DEL};

// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Surfaces
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Volume
Extrude {0, 0, H} {
  Surface{1}; Layers{H/DEL+1}; Recombine;
}

// Volume and boundary conditions
Physical Volume("fluid") = {1};
Physical Surface("periodic_x") = {25, 17};
Physical Surface("periodic_z") = {26, 1};
Physical Surface("walls") = {21, 13};
