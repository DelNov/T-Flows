N_HEIGHT = 61;
N_THICK  = 16;
N_DEPTH =  5;
BUMP    =  0.1;
PROG    =  0.8;

Point(1) = { 1.0, 0, 0};
Point(2) = { 1.0, 0, 1};
Point(3) = { 1.1, 0, 1};
Point(4) = { 1.1, 0, 0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Transfinite Curve {1, 3} = N_HEIGHT Using Bump BUMP;
Transfinite Curve {4}    = N_THICK  Using Progression PROG;
Transfinite Curve {2}    = N_THICK  Using Progression 1.0/PROG;
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Transfinite Surface {1} = {1, 4, 3, 2};
Recombine Surface "*";

Extrude {0, 0.25, 0} {
  Surface{1};  Layers {N_DEPTH}; Recombine;
}

Physical Surface("LEFT_WALL") = {13};
Physical Surface("RIGHT_WALL", 27) = {21};
Physical Surface("SIDE_WALLS") = {1, 26};
Physical Surface("BOTTOM_WALL") = {25};
Physical Surface("TOP_WALL") = {17};
Physical Volume("COLD_SOLID", 31) = {1};
