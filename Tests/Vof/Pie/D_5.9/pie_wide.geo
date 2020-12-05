//--------------------------
//
// Set dimension parameters
//
//--------------------------
RB = 0.00295;       // radius of the bubble

NP = 6;             // number of parts

RI = RB *  2.5;     // inner diameter
RO = RB * 18.0;     // vessel diameter
LY = RB * 36.0;     // height of the domain

D_MIN = RB / 15.0;  // coarse, medium and fine: 10, 15 or 20
D_MAX = D_MIN * 12.0;

Printf("Bubble radius %g", RB);
Printf("Domain height %g", LY);

NY = LY / D_MIN;
NT = 2.0 * RI * Pi / D_MIN / NP;
NR = RI / D_MIN;
ANGLE = 2.0 * Pi / NP;

Printf("Number of axial layers %g", NY);
Printf("Number of radial layers %g", NR);
Printf("Number of tangential layers %g", NT);

//------------------------------------------
// Define the surface inside and extrude it
//------------------------------------------

// Points
Point(1) = {0,                 0,  0,               D_MIN};
Point(2) = {RI * Cos(ANGLE),   0,  RI * Sin(ANGLE), D_MIN};
Point(3) = {RI,                0,  0,               D_MIN};

// Lines
Line(1) = {1, 2};
Line(2) = {1, 3};
Circle(3) = {2, 1, 3};

Transfinite Curve {1, 2} = NR Using Progression 1;
Transfinite Curve {3}    = NT Using Progression 1;

// Surface
Curve Loop(1) = {1, 3, -2};
Plane Surface(1) = {1};
Recombine Surface{1};
Mesh.Algorithm = 8;

// Extrude
Extrude {0, LY, 0} {
  Surface {1};
  Layers{NY};         // layers measured in cells, not nodes
  Recombine;
}

//---------------------------------------------------
// Define the surface on the side and extrude it too
//---------------------------------------------------

// Points
Point(11) = {RO,  0,   0,  D_MAX};
Point(12) = {RO,  LY,  0,  D_MAX};

// Lines
Line(15) = {3, 11};
Line(16) = {10, 12};
Line(17) = {12, 11};

// Surface
Curve Loop(2) = {14, 16, 17, -15};
Plane Surface(21) = {2};
Recombine Surface{21};
Mesh.Algorithm = 8;

// Extrude
Extrude {{0, 1, 0}, {0, 0, 0}, -2.0 * Pi / NP} {
  Surface{21};
  Layers{NT-1};
  Recombine;
}

//---------------------------------------
// Define volume and boundary conditions
//---------------------------------------
Physical Volume("Fluid") = {1, 2};

Physical Surface("Walls") = {1, 20, 34, 42};
Physical Surface("Cylinder") = {38};
Physical Surface("Cuts") = {11, 43, 21, 19};
