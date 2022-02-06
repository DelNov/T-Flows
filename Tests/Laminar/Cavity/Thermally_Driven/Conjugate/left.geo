N_HEIGHT = 61;
N_THICK  = 16;
N_DEPTH =  5;
BUMP    =  0.1;
PROG    =  0.8;

Point(1) = {-0.1, 0, 0};
Point(2) = {-0.1, 0, 1};
Point(3) = { 0.0, 0, 1};
Point(4) = { 0.0, 0, 0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Transfinite Curve {1, 3} = N_HEIGHT Using Bump BUMP;
Transfinite Curve {2}    = N_THICK  Using Progression PROG;
Transfinite Curve {4}    = N_THICK  Using Progression 1.0/PROG;
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Transfinite Surface {1} = {1, 4, 3, 2};
Recombine Surface "*";

Extrude {0, 0.25, 0} {
  Surface{1};  Layers {N_DEPTH}; Recombine;
}

Physical Surface("LEFT") = {13};
Physical Surface("RIGHT", 27) = {21};
Physical Surface("SIDE") = {1, 26};
Physical Surface("BOTTOM") = {25};
Physical Surface("TOP") = {17};
Physical Volume("SOLID", 31) = {1};
