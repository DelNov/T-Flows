N_plane = 61;
N_depth = 15;

Point(1) = {0, 0, 0};
Point(2) = {0, 0, 1};
Point(3) = {1, 0, 1};
Point(4) = {1, 0, 0};
//+
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Transfinite Curve {1, 2, 3, 4} = N_plane Using Bump 0.1;

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Transfinite Surface {1} = {1, 2, 3, 4};
Recombine Surface "*";

Extrude {0, 0.25, 0} {
  Surface{1};  Layers {N_depth}; Recombine;
}

Physical Surface("LEFT") = {13};
Physical Surface("RIGHT", 27) = {21};
Physical Surface("SIDE") = {1, 26};
Physical Surface("BOTTOM") = {25};
Physical Surface("TOP") = {17};
Physical Volume("AIR", 31) = {1};
