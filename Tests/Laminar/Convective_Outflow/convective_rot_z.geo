//------------------------------------------------------------------------------
//
//  Grid to reproduce Bottaro's case on convective outflow conditions
//
//------------------------------------------------------------------------------

// Resolutions in streamwise, wall-normal and periodic direction
NX = 200;
NZ =  20;
NY =   3;

// Angle of rotation
ANGLE = 0.0;  // 45.0 * Pi / 180.0;

// The geometry is so simple that I hardly care to comment
Point(1) = { 0, 0, 0, 1.0};
Point(2) = {10, 0, 0, 1.0};
Point(3) = {10, 0, 1, 1.0};
Point(4) = { 0, 0, 1, 1.0};

// Abort;
Rotate {{0, 0, 1}, {0, 0, 0}, ANGLE} {
  Point{2}; Point{3}; Point{1}; Point{4}; 
}

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};  Plane Surface(1) = {1};

Transfinite Curve {1, 3} = NX+1 Using Progression 1;
Transfinite Curve {2, 4} = NZ+1 Using Progression 1;

Transfinite Surface {1} = {1, 2, 3, 4};
Recombine Surface {1};

Extrude {-5 * Sin(ANGLE), 5 * Cos(ANGLE), 0} {
  Surface{1};
  Layers{NY};
  Recombine;
}

Physical Surface("inlet", 27) = {25};
Physical Surface("outlet", 28) = {17};
Physical Surface("bottom_wall", 29) = {13};
Physical Surface("top_wall", 30) = {21};
Physical Surface("periodic", 31) = {1, 26};
Physical Volume("fluid", 32) = {1};

