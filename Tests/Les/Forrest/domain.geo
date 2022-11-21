h = 20.0;

L = 9.6 * h;
W = 4.8 * h;
H = 3.0 * h;

NL = 97;
NW = 48;
NH = 30;

//------------------------------------------------------------------------------
// Points
//------------------------------------------------------------------------------
Point(1) = {0, 0, 0};
Point(2) = {L, 0, 0};
Point(3) = {L, W, 0};
Point(4) = {0, W, 0};

//------------------------------------------------------------------------------
// Lines
//------------------------------------------------------------------------------
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Transfinite Curve {1, 3} = NL Using Progression 1;
Transfinite Curve {4, 2} = NW Using Progression 1;

//------------------------------------------------------------------------------
// Surface
//------------------------------------------------------------------------------
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Transfinite Surface {1};
Recombine Surface{1};

//------------------------------------------------------------------------------
// Volume
//------------------------------------------------------------------------------
Extrude {0, 0, H} {
  Surface{1}; Layers{NH}; Recombine;
}

Physical Surface("BOTTOM", 27) = {1};
Physical Surface("TOP", 28) = {26};
Physical Surface("PERIODIC_X", 29) = {25, 17};
Physical Surface("PERIODIC_Y", 30) = {13, 21};
Physical Volume("DOMAIN", 31) = {1};
