
Point(1) = {0,   0,   0, 0.1};
Point(2) = {1,   0,   0, 0.1};
Point(3) = {1,   0.5, 0, 0.1};
Point(4) = {0.5, 0.5, 0, 0.1};
Point(5) = {0.5, 1.0, 0, 0.1};
Point(6) = {0,   1.0, 0, 0.1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

Curve Loop(1) = {1, 2, 3, 4, 5, 6};  Plane Surface(1) = {1};

Extrude {0, 0, 0.5} {
  Surface{1};
}

Physical Surface("x_direction") = {21, 29, 37};
Physical Surface("y_direction") = {17, 25, 33};
Physical Surface("z_direction") = {38, 1};
Physical Volume("fluid") = {1};
