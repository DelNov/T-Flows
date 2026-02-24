SetFactory("OpenCASCADE");

P   =  1.2;    // pitch
R   =  0.4;    // radius
H   =  0.4;    // height
X_0 =  0.0;
Y_0 =  0.0;
M   = 21;      // number of layers in height
N   = 21;      // number of nodes between rods

PNTS_CENT_OUT_RODS = 100;
PNTS_FRST_OUT_RODS = 200;
PNTS_LAST_OUT_RODS = 300;

ARC_CENT      = 1000;
ARCS_OUT      = 2000;
LINS_OUT_RODS = 3000;
N_FULL        = N * 2 * R * Pi / (P - 2*R);

//------------
//   Points
//------------

// Points on outer rods
angle_0 = Pi / 6;
For n In{ 1 : 6 }
  angle = angle_0 + (n-1) * Pi / 3;

  Point(PNTS_CENT_OUT_RODS + n) = {X_0 + P * Cos(angle),
                                   Y_0 + P * Sin(angle),
                                   0.0};

  Point(PNTS_LAST_OUT_RODS + n) = {X_0+P*Cos(angle) - R*Cos(angle-Pi/3),
                                   Y_0+P*Sin(angle) - R*Sin(angle-Pi/3),
                                   0.0};

  Point(PNTS_FRST_OUT_RODS + n) = {X_0+P*Cos(angle) - R*Cos(angle+Pi/3),
                                   Y_0+P*Sin(angle) - R*Sin(angle+Pi/3),
                                   0.0};
EndFor

//----------
//   Arcs
//----------

// Arcs on the central rod
Circle(ARC_CENT) = {X_0, Y_0, 0, R, 0, 2*Pi};

// Arcs on outer rods
For n In{ 1 : 6 }
  Circle(ARCS_OUT + n) = {PNTS_FRST_OUT_RODS + n,
                          PNTS_CENT_OUT_RODS + n,
                          PNTS_LAST_OUT_RODS + n};
EndFor

For i In{ 1 : 6 }
  j = i + 1;  If(j > 6) j = 1; EndIf
  Line(LINS_OUT_RODS + i) = {PNTS_LAST_OUT_RODS + i,
                             PNTS_FRST_OUT_RODS + j};
EndFor

// It would be nice to make this inside a loop
Curve Loop(1) = {2001, 3001,
                 2002, 3002,
                 2003, 3003,
                 2004, 3004,
                 2005, 3005,
                 2006, 3006};
Curve Loop(2) = {ARC_CENT};
Plane Surface(1) = {1, 2};

// Set transfinite curves in between rods
For n In{ 1 : 6 }
  Transfinite Curve {LINS_OUT_RODS + n} = N          Using Progression 1;
  Transfinite Curve {ARCS_OUT      + n} = N_FULL / 3 Using Progression 1;
EndFor

Transfinite Curve {ARC_CENT} = N_FULL Using Progression 1;

//------------
//   Volume
//------------
Extrude {0, 0, H} {
  Surface{1}; Layers{M};
  Recombine;
}

//-----------------------------------
//   Boundary and volume condition
//-----------------------------------
Physical Surface("top_and_bottom") = {1, 15};
Physical Surface("rods") = {2, 4, 6, 8, 10, 12, 14};
Physical Surface("dir_x") = {13, 7};
Physical Surface("dir_x_p60") = {3, 9};
Physical Surface("dir_x_m60") = {5, 11};
Physical Volume("all") = {1};
