//==============================================================================
//
//
// Extrusion mode
//
//
//------------------------------------------------------------------------------
HOT  = 1;
COLD = 2;

EXTRUDE = COLD;

//==============================================================================
//
//
// Constants
//
//
//------------------------------------------------------------------------------

R  =  0.5;
W  = 12.0;
WI =  3.0;
P  =  6;    // polygon faces
NP =  8;    // number of points per polygon
NI =  5;    // number of internal points in O-grid
ND = 41;    // 21;
NC = 41;    // 31;
NB = 41;    // 41;
NA = 41;    // 51;
NZ = 21;
NI =  8;

PNT_CENTRAL = 0;

P_IN_COLD  =  1000;  L_IN_COLD  =  1000;   S_IN_COLD  = 1000;
P_OUT_COLD =  2000;  L_OUT_COLD =  2000;   S_OUT_COLD = 2000;
P_IN_HOT   =  3000;  L_IN_HOT   =  3000;   S_IN_HOT   = 3000;
P_OUT_HOT  =  4000;  L_OUT_HOT  =  4000;   S_OUT_HOT  = 4000;

P_BND      =  5000;  L_BND =  5000;
P_SPIR     =  6000;
L_SPIR_BL  =  7000;
L_SPIR_A   =  8000;
L_SPIR_B   =  9000;
L_SPIR_C   = 10000;
L_SPIR_D   = 11000;
P_COLD     = 12000;  L_COLD     = 12000;  S_COLD = 12000;
P_HOT      = 13000;
L_TURN     = 14000;
L_HOT_BL   = 15000;  S_HOT_BL = 15000;

//==============================================================================
//
//
// Definition of macros
//
//
//------------------------------------------------------------------------------

//==============================================================================
   Macro NextFromP
//------------------------------------------------------------------------------
// Parameters: m, returns n
//------------------------------------------------------------------------------

  n = m + 1;
  If(n > P)
    n = 1;
  EndIf

  Return

//==============================================================================
   Macro EightPointsAroundCenter
//------------------------------------------------------------------------------
// Parameters: pnt_start, lin_start, xc, yc, angle_0, r
//------------------------------------------------------------------------------

  // Central point
  Point(pnt_start) = {xc, yc, 0.0};

  // Create bounding points
  For m In{ 1 : P }
    angle = angle_0 + (m-1) * 2 * Pi / P;
    Point(pnt_start+m) = {xc+r*Cos(angle), yc+r*Sin(angle), 0.0};
  EndFor

  // Create circular arcs
  For m In{ 1 : P }
    Call NextFromP;  // estimate next point
    Circle(lin_start+m) = {pnt_start+m, pnt_start, pnt_start+n};
    Transfinite Curve{lin_start+m} = NP Using Progression 1;
  EndFor

  // Create internal points
  For m In{ 1 : P }
    angle = angle_0 + (m-1) * 2 * Pi / P;
    Point(pnt_start + m + P) = {xc+0.7*r*Cos(angle), yc+0.7*r*Sin(angle), 0.0};
  EndFor

  // Create connections in between internal points
  For m In{ 1 : P }
    Call NextFromP;  // estimate next point
    Line(lin_start+m+P) = {pnt_start+m+P, pnt_start+n+P};
    Transfinite Curve{lin_start+m+P} = NP Using Progression 1;
  EndFor

  // Create connections between internal and external points
  For m In{ 1 : P }
    Line(lin_start+m+2*P) = {pnt_start+m, pnt_start+m+P};
    Transfinite Curve{lin_start+m+2*P} = NI Using Progression 1;
  EndFor

  // Create all the planar surfaces
  For m In{ 1 : P }
    n = m + 2*P + 1;
    If(m == P)
      n = m + P + 1;
    EndIf
    Curve Loop(sur_start+m) = {  lin_start+m,
                                 lin_start+n,
                               -(lin_start+m+P),
                               -(lin_start+m+2*P)};
    Plane Surface(sur_start+m) = {sur_start+m};
    Transfinite Surface{sur_start+m};
    Recombine Surface{sur_start+m};
  EndFor

  Curve Loop(sur_start+P+1) = {lin_start+P+1:lin_start+P+P};
  Plane Surface(sur_start+P+1) = {sur_start+P+1};
  Recombine Surface{sur_start+P+1};

  Return

//==============================================================================
//
//
// Planar geometry
//
//
//------------------------------------------------------------------------------

Point(PNT_CENTRAL)   = {0.0, 0.0, 0.0};
Point(PNT_CENTRAL+1) = {1.0, 0.0, 0.0};

// Call macro EightPointsAroundCenter to create inlet cold
pnt_start = P_IN_COLD;
lin_start = L_IN_COLD;
sur_start = S_IN_COLD;
xc        = -4.5;
yc        =  0.0;
angle_0   =  0.0;
r         =  R;
Call EightPointsAroundCenter;

// Call macro EightPointsAroundCenter to create outlet cold
pnt_start = P_OUT_COLD;
lin_start = L_OUT_COLD;
sur_start = S_OUT_COLD;
xc        = -0.5;
yc        =  0.0;
angle_0   =  0.0;
r         =  R;
Call EightPointsAroundCenter;

// Call macro EightPointsAroundCenter to create outlet cold
pnt_start = P_IN_HOT;
lin_start = L_IN_HOT;
sur_start = S_IN_HOT;
xc        = +1.0;
yc        =  0.0;
angle_0   =  0.0;
r         =  R;
Call EightPointsAroundCenter;

// Call macro EightPointsAroundCenter to create outlet cold
pnt_start = P_OUT_HOT;
lin_start = L_OUT_HOT;
sur_start = S_OUT_HOT;
xc        = -3.5;
yc        = -2.5;
angle_0   =  0.0;
r         =  R;
Call EightPointsAroundCenter;

//--------------------------
//  Define domain boundary
//--------------------------
m = P_BND;
Point(m) = { 0.0,  0.0,  0.0};  m++;
Point(m) = {-5.5,  0.0,  0.0};  m++;
Point(m) = { 0.0,  5.5,  0.0};  m++;
Point(m) = { 5.5,  0.0,  0.0};  m++;
Point(m) = { 0.0, -5.5,  0.0};  m++;

m = L_BND;
p = P_BND;
Circle(m) = {p+1, p, p+2}; m++;
Circle(m) = {p+2, p, p+3}; m++;
Circle(m) = {p+3, p, p+4}; m++;
Circle(m) = {p+4, p, p+1}; m++;

//---------------------------------

m = P_SPIR;
Point(m) = { 0.0,  5.0,   0.0};  m++;
Point(m) = { 0.0,  4.75,  0.0};  m++;
Point(m) = { 0.0,  4.25,  0.0};  m++;
Point(m) = { 0.0,  4.0,   0.0};  m++;

Point(m) = { 5.0,   0.0,  0.0};  m++;
Point(m) = { 4.75,  0.0,  0.0};  m++;
Point(m) = { 4.25,  0.0,  0.0};  m++;
Point(m) = { 4.0,   0.0,  0.0};  m++;

Point(m) = { 1.0, -4.0,   0.0};  m++;
Point(m) = { 1.0, -3.75,  0.0};  m++;
Point(m) = { 1.0, -3.25,  0.0};  m++;
Point(m) = { 1.0, -3.0,   0.0};  m++;

Point(m) = {-3.0,   0.0,  0.0};  m++;
Point(m) = {-2.75,  0.0,  0.0};  m++;
Point(m) = {-2.25,  0.0,  0.0};  m++;
Point(m) = {-2.0,   0.0,  0.0};  m++;

Point(m) = { 0.0,  3.0,   0.0};  m++;
Point(m) = { 0.0,  2.75,  0.0};  m++;
Point(m) = { 0.0,  2.25,  0.0};  m++;
Point(m) = { 0.0,  2.0,   0.0};  m++;

Point(m) = { 3.0,   0.0,  0.0};  m++;
Point(m) = { 2.75,  0.0,  0.0};  m++;
Point(m) = { 2.25,  0.0,  0.0};  m++;
Point(m) = { 2.0,   0.0,  0.0};  m++;

Point(m) = { 1.0, -2.0,   0.0};  m++;
Point(m) = { 1.0, -1.75,  0.0};  m++;
Point(m) = { 1.0, -1.25,  0.0};  m++;
Point(m) = { 1.0, -1.0,   0.0};  m++;

m = L_SPIR_A;
ms = m;
Circle(m)  = {P_IN_COLD+4, PNT_CENTRAL, P_SPIR};    m++;
Circle(m)  = {P_IN_COLD+1, PNT_CENTRAL, P_SPIR+3};  m++;
Circle(m)  = {P_SPIR,    PNT_CENTRAL, P_SPIR+4};  m++;
Circle(m)  = {P_SPIR+3,  PNT_CENTRAL, P_SPIR+7};  m++;

Ellipse(m) = {P_IN_COLD+2, PNT_CENTRAL, P_SPIR+2, P_SPIR+2};  m++;
Ellipse(m) = {P_IN_COLD+3, PNT_CENTRAL, P_SPIR+1, P_SPIR+1};  m++;
Ellipse(m) = {P_SPIR+2,  PNT_CENTRAL, P_SPIR+6, P_SPIR+6};  m++;
Ellipse(m) = {P_SPIR+1,  PNT_CENTRAL, P_SPIR+5, P_SPIR+5};  m++;
mc = m;
Transfinite Curve{ms:mc} = NA Using Progression 1;

m = L_SPIR_B;
ms = m;
Circle(m)  = {P_SPIR+ 4,  PNT_CENTRAL+1, P_SPIR+ 8};  m++;
Circle(m)  = {P_SPIR+ 7,  PNT_CENTRAL+1, P_SPIR+11};  m++;
Circle(m)  = {P_SPIR+ 8,  PNT_CENTRAL+1, P_SPIR+12};  m++;
Circle(m)  = {P_SPIR+11,  PNT_CENTRAL+1, P_SPIR+15};  m++;

Ellipse(m) = {P_SPIR+ 5, PNT_CENTRAL+1, P_SPIR+ 9, P_SPIR+ 9}; m++;
Ellipse(m) = {P_SPIR+ 6, PNT_CENTRAL+1, P_SPIR+10, P_SPIR+10}; m++;
Ellipse(m) = {P_SPIR+ 9, PNT_CENTRAL+1, P_SPIR+13, P_SPIR+13}; m++;
Ellipse(m) = {P_SPIR+10, PNT_CENTRAL+1, P_SPIR+14, P_SPIR+14}; m++;
mc = m;
Transfinite Curve{ms:mc} = NB Using Progression 1.0;

m = L_SPIR_C;
ms = m;
Circle(m) = {P_SPIR+12, PNT_CENTRAL, P_SPIR+16};  m++;
Circle(m) = {P_SPIR+15, PNT_CENTRAL, P_SPIR+19};  m++;
Circle(m) = {P_SPIR+16, PNT_CENTRAL, P_SPIR+20};  m++;
Circle(m) = {P_SPIR+19, PNT_CENTRAL, P_SPIR+23};  m++;

Ellipse(m) = {P_SPIR+13, PNT_CENTRAL, P_SPIR+17, P_SPIR+17}; m++;
Ellipse(m) = {P_SPIR+14, PNT_CENTRAL, P_SPIR+18, P_SPIR+18}; m++;
Ellipse(m) = {P_SPIR+17, PNT_CENTRAL, P_SPIR+21, P_SPIR+21}; m++;
Ellipse(m) = {P_SPIR+18, PNT_CENTRAL, P_SPIR+22, P_SPIR+22}; m++;
mc = m;
Transfinite Curve{ms:mc} = NC Using Progression 1.0;

m = L_SPIR_D;
ms = m;
Circle(m)  = {P_SPIR+20,  PNT_CENTRAL+1, P_SPIR+24};   m++;
Circle(m)  = {P_SPIR+23,  PNT_CENTRAL+1, P_SPIR+27};   m++;
Circle(m)  = {P_SPIR+24,  PNT_CENTRAL+1, P_OUT_COLD+4};  m++;
Circle(m)  = {P_SPIR+27,  PNT_CENTRAL+1, P_OUT_COLD+1};  m++;

Ellipse(m) = {P_SPIR+21, PNT_CENTRAL+1, P_SPIR+25, P_SPIR+25}; m++;
Ellipse(m) = {P_SPIR+22, PNT_CENTRAL+1, P_SPIR+26, P_SPIR+26}; m++;
Ellipse(m) = {P_SPIR+25, 1, P_OUT_COLD+5, P_OUT_COLD+5};  m++;
Ellipse(m) = {P_SPIR+26, 1, P_OUT_COLD+6, P_OUT_COLD+6};  m++;
mc = m;
Transfinite Curve{ms:mc} = ND Using Progression 1.0;

// Lines perependicular to cold spiral
m = L_SPIR_BL;

// These are in boundary layers
ms = m;
Line(m) = {P_SPIR+ 0, P_SPIR+ 1};  m++;
Line(m) = {P_SPIR+ 3, P_SPIR+ 2};  m++;
Line(m) = {P_SPIR+ 4, P_SPIR+ 5};  m++;
Line(m) = {P_SPIR+ 7, P_SPIR+ 6};  m++;
Line(m) = {P_SPIR+ 8, P_SPIR+ 9};  m++;
Line(m) = {P_SPIR+11, P_SPIR+10};  m++;
Line(m) = {P_SPIR+12, P_SPIR+13};  m++;
Line(m) = {P_SPIR+15, P_SPIR+14};  m++;
Line(m) = {P_SPIR+16, P_SPIR+17};  m++;
Line(m) = {P_SPIR+19, P_SPIR+18};  m++;
Line(m) = {P_SPIR+20, P_SPIR+21};  m++;
Line(m) = {P_SPIR+23, P_SPIR+22};  m++;
Line(m) = {P_SPIR+24, P_SPIR+25};  m++;
Line(m) = {P_SPIR+27, P_SPIR+26};  m++;
mc = m;
Transfinite Curve{ms:mc} = NP Using Progression 1.2;

ms = m;
Line(m) = {P_SPIR+ 2, P_SPIR+ 1};  m++;
Line(m) = {P_SPIR+ 6, P_SPIR+ 5};  m++;
Line(m) = {P_SPIR+10, P_SPIR+ 9};  m++;
Line(m) = {P_SPIR+14, P_SPIR+13};  m++;
Line(m) = {P_SPIR+18, P_SPIR+17};  m++;
Line(m) = {P_SPIR+22, P_SPIR+21};  m++;
Line(m) = {P_SPIR+26, P_SPIR+25};  m++;
mc = m;
Transfinite Curve{ms:mc} = NP Using Progression 1.0;

m = S_COLD;
ms = m;
Curve Loop(m) = {1001,  L_SPIR_A+4,  -7001,  -L_SPIR_A-1};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {1002,  L_SPIR_A+5,  -7014,  -L_SPIR_A-4};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {1003,  L_SPIR_A+0,   7000,  -L_SPIR_A-5};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {7001,  L_SPIR_A+6,  -7003,  -L_SPIR_A-3};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {7014,  L_SPIR_A+7,  -7015,  -L_SPIR_A-6};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {7000,  L_SPIR_A+7,  -7002,  -L_SPIR_A-2};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {7003,  L_SPIR_B+5,  -7005,  -L_SPIR_B-1};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {7015,  L_SPIR_B+4,  -7016,  -L_SPIR_B-5};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {7002,  L_SPIR_B+4,  -7004,  -L_SPIR_B-0};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {7005,  L_SPIR_B+7,  -7007,  -L_SPIR_B-3};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {7016,  L_SPIR_B+6,  -7017,  -L_SPIR_B-7};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {7004,  L_SPIR_B+6,  -7006,  -L_SPIR_B-2};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {7006,  L_SPIR_C+04, -7008, -L_SPIR_C-0};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {7017,  L_SPIR_C+04, -7018, -L_SPIR_C-5};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {7007,  L_SPIR_C+05, -7009, -L_SPIR_C-1};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {7008,  L_SPIR_C+06, -7010, -L_SPIR_C-2};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {7018,  L_SPIR_C+06, -7019, -L_SPIR_C-7};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {7009,  L_SPIR_C+07, -7011, -L_SPIR_C-3};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {7010,  L_SPIR_D+4, -7012, -L_SPIR_D-0};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {7019,  L_SPIR_D+4, -7020, -L_SPIR_D-5};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {7011,  L_SPIR_D+5, -7013, -L_SPIR_D-1};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {2006, -L_SPIR_D-3,  7013,  L_SPIR_D+7};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {2005, -L_SPIR_D-7,  7020,  L_SPIR_D+6};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {2004, -L_SPIR_D-6, -7012,  L_SPIR_D+2};  Plane Surface(m) = {m};  m++;
mc = m;

For m In{ ms : mc }
  Transfinite Surface{m};
  Recombine Surface{m};
EndFor

m = P_HOT;

Point(m) = {-5.25,  0,     0,  1.0};  m++;
Point(m) = { 0,     5.25,  0,  1.0};  m++;
Point(m) = { 5.25,  0,     0,  1.0};  m++;
Point(m) = { 1.00, -4.25,  0,  1.0};  m++;
Point(m) = {-3.25,  0,     0,  1.0};  m++;
Point(m) = { 0,     3.25,  0,  1.0};  m++;
Point(m) = { 3.25,  0,     0,  1.0};  m++;
Point(m) = { 1.0,  -2.25,  0,  1.0};  m++;
Point(m) = {-1.25,  0,     0,  1.0};  m++;

Point(m) = {-0.5-0.75/2, 0.75*Sqrt(3)/2.0,  0.0};  m++;
Point(m) = {-0.5+0.75/2, 0.75*Sqrt(3)/2.0,  0.0};  m++;

Point(m) = { 0.25,  0,     0,  1.0};  m++;
Point(m) = { 1,    -0.75,  0,  1.0};  m++;
Point(m) = { 1.75,  0,     0,  1.0};  m++;
Point(m) = { 0,     1.75,  0,  1.0};  m++;
Point(m) = {-1.75,  0,     0,  1.0};  m++;
Point(m) = { 1,    -2.75,  0,  1.0};  m++;
Point(m) = { 3.75,  0,     0,  1.0};  m++;
Point(m) = { 0,     3.75,  0,  1.0};  m++;
Point(m) = {-3.75,  0,     0,  1.0};  m++;

Point(m) = {-4.5-0.75/2, -0.75*Sqrt(3)/2.0,  0.0};  m++;
Point(m) = {-4.5+0.75/2, -0.75*Sqrt(3)/2.0,  0.0};  m++;

m = L_SPIR_A + 8;
Circle(m) = {P_HOT,    PNT_CENTRAL+0, P_HOT+01};  m++;
Circle(m) = {P_HOT+01, PNT_CENTRAL+0, P_HOT+02};  m++;
Transfinite Curve{m-2:m-1} = NA Using Progression 1;

m = L_SPIR_B + 8;
Circle(m) = {P_HOT+02, PNT_CENTRAL+1, P_HOT+03};  m++;
Circle(m) = {P_HOT+03, PNT_CENTRAL+1, P_HOT+04};  m++;
Transfinite Curve{m-2:m-1} = NB Using Progression 1;

m = L_SPIR_C + 8;
Circle(m) = {P_HOT+04, PNT_CENTRAL+0, P_HOT+05};  m++;
Circle(m) = {P_HOT+05, PNT_CENTRAL+0, P_HOT+06};  m++;
Transfinite Curve{m-2:m-1} = NC Using Progression 1;

m = L_SPIR_D + 8;
Circle(m) = {P_HOT+06, PNT_CENTRAL+1, P_HOT+07};  m++;
Circle(m) = {P_HOT+07, PNT_CENTRAL+1, P_HOT+08};  m++;
Transfinite Curve{m-2:m-1} = ND Using Progression 1;

m = L_TURN;
Circle(m) = {P_HOT+08, P_OUT_COLD, P_HOT+09};  m++;
Circle(m) = {P_HOT+09, P_OUT_COLD, P_HOT+10};  m++;
Circle(m) = {P_HOT+10, P_OUT_COLD, P_HOT+11};  m++;
Transfinite Curve{m-3:m-1} = NP Using Progression 1;

m = L_SPIR_D + 8 + 2;
Circle(m) = {P_HOT+11, PNT_CENTRAL+1, P_HOT+12};  m++;
Circle(m) = {P_HOT+12, PNT_CENTRAL+1, P_HOT+13};  m++;
Transfinite Curve{m-2:m-1} = ND Using Progression 1;

m = L_SPIR_C + 8 + 2;
Circle(m) = {P_HOT+13, PNT_CENTRAL+0, P_HOT+14};  m++;
Circle(m) = {P_HOT+14, PNT_CENTRAL+0, P_HOT+15};  m++;
Transfinite Curve{m-2:m-1} = NC Using Progression 1;

m = L_SPIR_B + 8 + 2;
Circle(m) = {P_HOT+15, PNT_CENTRAL+1, P_HOT+16};  m++;
Circle(m) = {P_HOT+16, PNT_CENTRAL+1, P_HOT+17};  m++;
Transfinite Curve{m-2:m-1} = NB Using Progression 1;

m = L_SPIR_A + 8 + 2;
Circle(m) = {P_HOT+17, PNT_CENTRAL+0, P_HOT+18};  m++;
Circle(m) = {P_HOT+18, PNT_CENTRAL+0, P_HOT+19};  m++;
Transfinite Curve{m-2:m-1} = NA Using Progression 1;

m = L_TURN + 3;
Circle(m) = {P_HOT+19, P_IN_COLD, P_HOT+21};  m++;
Circle(m) = {P_HOT+21, P_IN_COLD, P_HOT+20};  m++;
Circle(m) = {P_HOT+20, P_IN_COLD, P_HOT+00};  m++;
Transfinite Curve{m-3:m-1} = NP Using Progression 1;

m = L_HOT_BL;
ms = m;
Line(m) = {P_HOT,    P_IN_COLD+04};  m++;
Line(m) = {P_HOT+01, P_SPIR+00};  m++;
Line(m) = {P_HOT+02, P_SPIR+04};  m++;
Line(m) = {P_HOT+03, P_SPIR+08};  m++;
Line(m) = {P_HOT+04, P_SPIR+12};  m++;
Line(m) = {P_HOT+05, P_SPIR+16};  m++;
Line(m) = {P_HOT+06, P_SPIR+20};  m++;
Line(m) = {P_HOT+07, P_SPIR+24};  m++;
Line(m) = {P_HOT+08, P_OUT_COLD+4};  m++;
Line(m) = {P_HOT+09, P_OUT_COLD+3};  m++;
Line(m) = {P_HOT+10, P_OUT_COLD+2};  m++;
Line(m) = {P_HOT+11, P_OUT_COLD+1};  m++;
Line(m) = {P_HOT+12, P_SPIR+27};  m++;
Line(m) = {P_HOT+13, P_SPIR+23};  m++;
Line(m) = {P_HOT+14, P_SPIR+19};  m++;
Line(m) = {P_HOT+15, P_SPIR+15};  m++;
Line(m) = {P_HOT+16, P_SPIR+11};  m++;
Line(m) = {P_HOT+17, P_SPIR+07};  m++;
Line(m) = {P_HOT+18, P_SPIR+03};  m++;
Line(m) = {P_HOT+19, P_IN_COLD+01};  m++;
Line(m) = {P_HOT+21, P_IN_COLD+06};  m++;
Line(m) = {P_HOT+20, P_IN_COLD+05};  m++;
mc = m;
Transfinite Curve{ms:mc} = NP Using Progression 1.0/1.2;

// The following five are in inner, not boundary layers
ms = m;
Line(m) = {P_HOT+19, P_HOT+04};  m++;
Line(m) = {P_HOT+18, P_HOT+05};  m++;
Line(m) = {P_HOT+17, P_HOT+06};  m++;
Line(m) = {P_HOT+16, P_HOT+07};  m++;
Line(m) = {P_HOT+15, P_HOT+08};  m++;
mc = m;
Transfinite Curve{ms:mc} = NP Using Progression 1;

// These three are in outer layers
ms = m;
Line(m) = {P_BND+1, P_HOT};    m++;
Line(m) = {P_BND+2, P_HOT+1};  m++;
Line(m) = {P_BND+3, P_HOT+2};  m++;
mc = m;
Transfinite Curve{ms:mc} = 3 Using Progression 1;

m = S_HOT_BL;
ms = m;
Curve Loop(m) = {L_HOT_BL,     L_SPIR_A,    -L_HOT_BL- 1, -L_SPIR_A-8};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {L_HOT_BL+ 1,  L_SPIR_A+ 2, -L_HOT_BL- 2, -L_SPIR_A-9};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {L_HOT_BL+ 2,  L_SPIR_B,    -L_HOT_BL- 3, -L_SPIR_B-8};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {L_HOT_BL+ 3,  L_SPIR_B+ 2, -L_HOT_BL- 4, -L_SPIR_B-9};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {L_HOT_BL+ 4,  L_SPIR_C,    -L_HOT_BL- 5, -L_SPIR_C-8};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {L_HOT_BL+ 5,  L_SPIR_C+ 2, -L_HOT_BL- 6, -L_SPIR_C-9};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {L_HOT_BL+ 6,  L_SPIR_D,    -L_HOT_BL- 7, -L_SPIR_D-8};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {L_HOT_BL+ 7,  L_SPIR_D+2,  -L_HOT_BL- 8, -L_SPIR_D-9};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {L_HOT_BL+ 8, -L_OUT_COLD-3,-L_HOT_BL- 9, -L_TURN};      Plane Surface(m) = {m};  m++;
Curve Loop(m) = {L_HOT_BL+ 9, -L_OUT_COLD-2,-L_HOT_BL-10, -L_TURN-01};   Plane Surface(m) = {m};  m++;
Curve Loop(m) = {L_HOT_BL+10, -L_OUT_COLD-1,-L_HOT_BL-11, -L_TURN-02};   Plane Surface(m) = {m};  m++;
Curve Loop(m) = {L_HOT_BL+11, -L_SPIR_D-3,  -L_HOT_BL-12, -L_SPIR_D-10}; Plane Surface(m) = {m};  m++;
Curve Loop(m) = {L_HOT_BL+12, -L_SPIR_D-1,  -L_HOT_BL-13, -L_SPIR_D-11}; Plane Surface(m) = {m};  m++;
Curve Loop(m) = {L_HOT_BL+13, -L_SPIR_C-3,  -L_HOT_BL-14, -L_SPIR_C-10}; Plane Surface(m) = {m};  m++;
Curve Loop(m) = {L_HOT_BL+14, -L_SPIR_C-1,  -L_HOT_BL-15, -L_SPIR_C-11}; Plane Surface(m) = {m};  m++;
Curve Loop(m) = {L_HOT_BL+15, -L_SPIR_B-3,  -L_HOT_BL-16, -L_SPIR_B-10}; Plane Surface(m) = {m};  m++;
Curve Loop(m) = {L_HOT_BL+16, -L_SPIR_B-1,  -L_HOT_BL-17, -L_SPIR_B-11}; Plane Surface(m) = {m};  m++;
Curve Loop(m) = {L_HOT_BL+17, -L_SPIR_A-3,  -L_HOT_BL-18, -L_SPIR_A-10}; Plane Surface(m) = {m};  m++;
Curve Loop(m) = {L_HOT_BL+18, -L_SPIR_A-1,  -L_HOT_BL-19, -L_SPIR_A-11}; Plane Surface(m) = {m};  m++;
Curve Loop(m) = {L_HOT_BL+19, -L_IN_COLD-6, -L_HOT_BL-20, -L_TURN-3};    Plane Surface(m) = {m};  m++;
Curve Loop(m) = {L_HOT_BL+20, -L_IN_COLD-5, -L_HOT_BL-21, -L_TURN-4};    Plane Surface(m) = {m};  m++;
Curve Loop(m) = {L_HOT_BL+21, -L_IN_COLD-4, -L_HOT_BL,    -L_TURN-5};    Plane Surface(m) = {m};  m++;

// The following four are in inner, not boundary layers
Curve Loop(m) = {L_HOT_BL+22, L_SPIR_C+8, -L_HOT_BL-23, L_SPIR_A+11};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {L_HOT_BL+23, L_SPIR_C+9, -L_HOT_BL-24, L_SPIR_A+10};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {L_HOT_BL+24, L_SPIR_D+8, -L_HOT_BL-25, L_SPIR_B+11};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {L_HOT_BL+25, L_SPIR_D+9, -L_HOT_BL-26, L_SPIR_B+10};  Plane Surface(m) = {m};  m++;

// These two are in outer layers
Curve Loop(m) = {L_HOT_BL+27, 8008, -L_HOT_BL-28, -5000};  Plane Surface(m) = {m};  m++;
Curve Loop(m) = {L_HOT_BL+28, 8009, -L_HOT_BL-29, -5001};  Plane Surface(m) = {m};  m++;

mc = m;

For m In{ ms : mc }
  Transfinite Surface{m};
  Recombine Surface{m};
EndFor

Transfinite Curve {P_BND  :P_BND+1} = NA    Using Progression 1;
Transfinite Curve {P_BND+2:P_BND+3} = NA+10 Using Progression 1;

// Inner area
Curve Loop(1) = {L_SPIR_C+11,
                 L_HOT_BL+26,
                 L_TURN,
                 L_TURN+1,
                 L_TURN+2,
                 L_SPIR_D+10,
                 L_SPIR_D+11,
                 L_SPIR_C+10};
Curve Loop(2) = {P_IN_HOT+1:P_IN_HOT+6};
Plane Surface(1) = {1, 2};
Recombine Surface{1};

Curve Loop(3) = {L_HOT_BL+29,
                 L_SPIR_B+8,
                 L_SPIR_B+9,
                -L_HOT_BL-22,
                 L_TURN+03,
                 L_TURN+04,
                 L_TURN+05,
                -L_HOT_BL-27,
                -5003,
                -5002};
Curve Loop(4) = {P_OUT_HOT+1:P_OUT_HOT+6};
Plane Surface(2) = {3, 4};
Recombine Surface{2};

//==============================================================================
//
//
// Extrude to 3D
//
//
//------------------------------------------------------------------------------

//==============================================================================
//
// Extrude cold
//
//------------------------------------------------------------------------------
If(EXTRUDE == COLD)
  // Extrusion of the hot stream
  Extrude {0, 0, W} {
    Surface{S_IN_COLD+1 :S_IN_COLD +7};
    Surface{S_OUT_COLD+1:S_OUT_COLD+7};
    Surface{S_COLD:S_COLD+23};
    Layers{21}; Recombine;
  }
  Extrude {0, 0, -WI} {
    Surface{S_IN_COLD+1:S_IN_COLD+7};
    Layers{NI}; Recombine;
  }
  Extrude {0, 0, WI} {
    Surface{15325};
    Surface{15215};
    Surface{15237};
    Surface{15259};
    Surface{15281};
    Surface{15303};
    Surface{15357};
    Layers{NI}; Recombine;
  }
  Recursive Delete {
    Surface{1:2};
    Surface{S_IN_HOT+1:S_IN_HOT+7};
    Surface{S_OUT_HOT+1:S_OUT_HOT+7};
    Surface{S_HOT_BL:S_HOT_BL+27};
  }

  Physical Surface("membrane") = {15414, 15488, 15554, 15620, 15642, 15708,
                                  15774, 15884, 15246, 15224, 15202, 15832,
                                  15818, 15752, 15686, 15576, 15510, 15444,
                                  15378, 15148, 15126, 15104};

  Physical Surface("walls") = {15193, 15051, 15073, 15095, 15117, 15139, 15161,
                               15423, 15401, 15379, 15489, 15467, 15445, 15555,
                               15533, 15511, 15577, 15599, 15621, 15643, 15665,
                               15687, 15709, 15731, 15753, 15775, 15797, 15819,
                               15841, 15863, 15885, 16058, 16102, 16080, 16146,
                               16124, 16168, 15960, 15894, 15982, 15916, 15938,
                               16004, 12000, 12001, 12002, 12003, 12004, 12005,
                               12006, 12007, 12008, 12011, 12010, 12009, 12014,
                               12013, 12012, 12017, 12016, 12015, 12020, 12019,
                               12018, 12023, 12022, 12021, 2004, 2003, 2002,
                               2001, 2006, 2005, 2007};

  Physical Surface("cold_in") = {15995, 15973, 15951, 15929, 15907, 16017, 16049};

  Physical Surface("cold_out") = {16115, 16137, 16159, 16181, 16071, 16093, 16213};

  Physical Volume("cold") = {1:52};
EndIf

//==============================================================================
//
// Extrude hot
//
//------------------------------------------------------------------------------
If(EXTRUDE == HOT)
  Extrude {0, 0, W} {
    Surface{1:2};
    Surface{S_IN_HOT+1:S_IN_HOT+7};
    Surface{S_OUT_HOT+1:S_OUT_HOT+7};
    Surface{S_HOT_BL:S_HOT_BL+27};
    Layers{21}; Recombine;
  }
  Extrude {0, 0, -WI} {
    Surface{S_IN_HOT+1:S_IN_HOT+7};
    Layers{NI}; Recombine;
  }
  Extrude {0, 0, WI} {
    Surface{15369};
    Surface{15391};
    Surface{15413};
    Surface{15435};
    Surface{15457};
    Surface{15479};
    Surface{15511};
    Layers{NI}; Recombine;
  }
  Recursive Delete {
    Surface{S_IN_COLD+1 :S_IN_COLD +7};
    Surface{S_OUT_COLD+1:S_OUT_COLD+7};
    Surface{S_COLD:S_COLD+23};
  }

  Physical Surface("membrane") = {15524, 15546, 15568, 15590, 15612, 15634,
                                  15656, 15678, 15700, 15722, 15744, 15766,
                                  15788, 15810, 15832, 15854, 15876, 15898,
                                  15920, 15942, 15964, 15986};

  Physical Surface("cylinder") = {15158, 16126, 16104, 15154};

  Physical Surface("walls") = {16105, 15533, 16127, 15555, 15929, 16017, 16039,
                               15907, 15643, 15621, 15183, 15577, 15665, 16061,
                               15885, 16083, 15687, 15863, 15841, 15101, 15819,
                               15797, 15775, 15227, 15249, 15271, 15347, 15315,
                               15293, 15205, 15753, 15731, 15709, 15973, 15995,
                               15951, 15599, 16344, 16300, 16322, 16388, 16366,
                               16410, 2, 15026, 15000, 15021, 15020, 15019,
                               15018, 15022, 15004, 15002, 15001, 15027, 15005,
                               15023, 15017, 15006, 15024, 15016, 15015, 15025,
                               15007, 15014, 1, 15012, 15013, 15011, 15008,
                               15009, 15010, 15003, 4005, 4004, 4006, 4003,
                               4002, 4001, 4007, 16224, 16202, 16180, 16158,
                               16136, 16246};

  Physical Surface("hot_in") = {16237, 16215, 16193, 16171, 16149, 16259, 16291};

  Physical Surface("hot_out") = {16335, 16357, 16379, 16401, 16423, 16313, 16455};

  Physical Volume("hot") = {1:58};
EndIf

