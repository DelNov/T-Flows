N_PLANE = 61;
N_DEPTH =  5;
BUMP    =  0.1;

Point(1) = {0, 0, 0};
Point(2) = {0, 0, 1};
Point(3) = {1, 0, 1};
Point(4) = {1, 0, 0};
//+
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Transfinite Curve {1, 2, 3, 4} = N_PLANE Using Bump BUMP;

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Transfinite Surface {1} = {1, 2, 3, 4};
Recombine Surface "*";

Extrude {0, 0.25, 0} {
  Surface{1};  Layers {N_DEPTH}; Recombine;
}

Physical Surface("LEFT_WALL") = {13};
Physical Surface("RIGHT_WALL", 27) = {21};
Physical Surface("SIDE_WALLS") = {1, 26};
Physical Surface("BOTTOM_WALL") = {25};
Physical Surface("TOP_WALL") = {17};
Physical Volume("FLUID", 31) = {1};
