//==============================================================================
//
// MATRIX OF CUBES
//
//==============================================================================

H = 0.015;    // 15 mm
h = 3.4*H;    // 51 mm

// Ns are in cells
N_CUBE = 40;  B_CUBE = 0.2;
N_HIGH = 48;  B_HIGH = 0.12;
N_REST = 30;  P_REST = 1.1;

POINT_DELTA = 100;
POINT_FLOOR =   0;
POINT_CUBE  = POINT_FLOOR + POINT_DELTA;
POINT_TOP   = POINT_CUBE  + POINT_DELTA;

// Beginnings of horizontal lines
LINE_DELTA   = 100;
LINE_H_FLOOR =   0;
LINE_H_CUBE  = LINE_H_FLOOR + LINE_DELTA;
LINE_H_TOP   = LINE_H_CUBE  + LINE_DELTA;
LINE_V_LOW   = LINE_H_TOP   + LINE_DELTA;
LINE_V_HIGH  = LINE_V_LOW   + LINE_DELTA;

// Beginnings of "z" surfaces
SURF_DELTA   = 100;
SURF_Z_FLOOR =   0;
SURF_Z_CUBE  = SURF_Z_FLOOR + SURF_DELTA;
SURF_Z_TOP   = SURF_Z_CUBE  + SURF_DELTA;
SURF_Y_LOW   = SURF_Z_TOP   + SURF_DELTA;
SURF_Y_HIGH  = SURF_Y_LOW   + SURF_DELTA;
SURF_X_LOW   = SURF_Y_HIGH  + SURF_DELTA;
SURF_X_HIGH  = SURF_X_LOW   + SURF_DELTA;

TINY = 1.0e-9;
HUGE = 1.0e+9;

//==============================================================================
//
// POINTS
//
//==============================================================================

For i_f In{ 1 : 3 }
  If(i_f == 1)  SP = POINT_FLOOR;  Z = 0;  EndIf
  If(i_f == 2)  SP = POINT_CUBE;   Z = H;  EndIf
  If(i_f == 3)  SP = POINT_TOP;    Z = h;  EndIf
  Point(SP+ 1) = {0,       0,         Z};
  Point(SP+ 2) = {2*H-H/2, 0,         Z};
  Point(SP+ 3) = {2*H+H/2, 0,         Z};
  Point(SP+ 4) = {4*H,     0,         Z};
  Point(SP+ 5) = {0,       2*H-H/2,   Z};
  Point(SP+ 6) = {2*H-H/2, 2*H-H/2,   Z};
  Point(SP+ 7) = {2*H+H/2, 2*H-H/2,   Z};
  Point(SP+ 8) = {4*H,     2*H-H/2,   Z};
  Point(SP+ 9) = {0,       2*H+H/2,   Z};
  Point(SP+10) = {2*H-H/2, 2*H+H/2,   Z};
  Point(SP+11) = {2*H+H/2, 2*H+H/2,   Z};
  Point(SP+12) = {4*H,     2*H+H/2,   Z};
  Point(SP+13) = {0,       4*H,       Z};
  Point(SP+14) = {2*H-H/2, 4*H,       Z};
  Point(SP+15) = {2*H+H/2, 4*H,       Z};
  Point(SP+16) = {4*H,     4*H,       Z};
EndFor

//==============================================================================
//
// LINES
//
//==============================================================================

//------------------
// Horizontal lines
//------------------
For i_f In{ 1 : 3 }
  If(i_f == 1)  SL = LINE_H_FLOOR;  SP = POINT_FLOOR;  EndIf
  If(i_f == 2)  SL = LINE_H_CUBE;   SP = POINT_CUBE;   EndIf
  If(i_f == 3)  SL = LINE_H_TOP;    SP = POINT_TOP;    EndIf
  // Longer lines
  Line(SL+ 1) = {SP +2, SP +1};
  Line(SL+ 2) = {SP +3, SP +4};
  Line(SL+ 3) = {SP +6, SP +5};
  Line(SL+ 4) = {SP +7, SP +8};
  Line(SL+ 5) = {SP+10, SP +9};
  Line(SL+ 6) = {SP+11, SP+12};
  Line(SL+ 7) = {SP+14, SP+13};
  Line(SL+ 8) = {SP+15, SP+16};
  Line(SL+ 9) = {SP+ 5, SP+ 1};
  Line(SL+10) = {SP+ 6, SP+ 2};
  Line(SL+11) = {SP+ 7, SP+ 3};
  Line(SL+12) = {SP+ 8, SP+ 4};
  Line(SL+13) = {SP+ 9, SP+13};
  Line(SL+14) = {SP+10, SP+14};
  Line(SL+15) = {SP+11, SP+15};
  Line(SL+16) = {SP+12, SP+16};
  // Shorter lines
  Line(SL+17) = {SP+ 2, SP+ 3};
  Line(SL+18) = {SP+ 6, SP+ 7};
  Line(SL+19) = {SP+10, SP+11};
  Line(SL+20) = {SP+14, SP+15};
  Line(SL+21) = {SP+ 5, SP+ 9};
  Line(SL+22) = {SP+ 6, SP+10};
  Line(SL+23) = {SP+ 7, SP+11};
  Line(SL+24) = {SP+ 8, SP+12};
EndFor

// Vertical lines
For i_f In{ 1 : 2 }
  If(i_f == 1)  SP = POINT_FLOOR;  SL = LINE_V_LOW;   EndIf
  If(i_f == 2)  SP = POINT_CUBE;   SL = LINE_V_HIGH;  EndIf
  For j In{ 1 : 4 }
    For i In{ 1 : 4 }
      FP = (j-1) * 4 + i;
      Printf("First point for defining vertical line %g", FP);
      SL = SL + 1;
      Line(SL) = {SP + FP, SP + FP + POINT_DELTA};
    EndFor
  EndFor
EndFor

// Define longer lines as transfinite
For i_f In{ 1 : 3 }
  If(i_f == 1)  SL = LINE_H_FLOOR;  EndIf
  If(i_f == 2)  SL = LINE_H_CUBE;   EndIf
  If(i_f == 3)  SL = LINE_H_TOP;    EndIf
  For l In{ 1 : 16 }
    Transfinite Curve {SL+l} = N_REST+1 Using Progression P_REST;
  EndFor
EndFor

// Define shorter lines as transfinite
For i_f In{ 1 : 3 }
  If(i_f == 1)  SL = LINE_H_FLOOR;  EndIf
  If(i_f == 2)  SL = LINE_H_CUBE;   EndIf
  If(i_f == 3)  SL = LINE_H_TOP;    EndIf
  For l In{ 17 : 24 }
    Transfinite Curve {SL+l} = N_CUBE+1 Using Bump B_CUBE;
  EndFor
EndFor

Transfinite Curve {LINE_V_LOW :LINE_V_LOW +16} = N_CUBE+1 Using Bump B_CUBE;
Transfinite Curve {LINE_V_HIGH:LINE_V_HIGH+16} = N_HIGH+1 Using Bump B_HIGH;

//==============================================================================
//
// SURFACES
//
//==============================================================================

//----------------------------------------------
// Define horizontal surfaces
// (These are a bit too messy for an algorithm)
//----------------------------------------------
For i_f In{ 1 : 3 }
  If(i_f == 1)  SL = LINE_H_FLOOR;  SS = SURF_Z_FLOOR;  EndIf
  If(i_f == 2)  SL = LINE_H_CUBE;   SS = SURF_Z_CUBE;   EndIf
  If(i_f == 3)  SL = LINE_H_TOP;    SS = SURF_Z_TOP;    EndIf
  Curve Loop(SS+1) = {SL+ 1, -(SL+ 3),   SL+10,  -(SL+ 9)}; Plane Surface(SS+1) = {SS+1};
  Curve Loop(SS+2) = {SL+17, -(SL+11), -(SL+18),   SL+10};  Plane Surface(SS+2) = {SS+2};
  Curve Loop(SS+3) = {SL+ 2, -(SL+12), -(SL+ 4),   SL+11};  Plane Surface(SS+3) = {SS+3};
  Curve Loop(SS+4) = {SL+ 3,   SL+21,  -(SL+ 5), -(SL+22)}; Plane Surface(SS+4) = {SS+4};
  Curve Loop(SS+5) = {SL+18,   SL+23,  -(SL+19), -(SL+22)}; Plane Surface(SS+5) = {SS+5};
  Curve Loop(SS+6) = {SL+ 4,   SL+24,  -(SL+ 6), -(SL+23)}; Plane Surface(SS+6) = {SS+6};
  Curve Loop(SS+7) = {SL+ 5,   SL+13,  -(SL+ 7), -(SL+14)}; Plane Surface(SS+7) = {SS+7};
  Curve Loop(SS+8) = {SL+19,   SL+15,  -(SL+20), -(SL+14)}; Plane Surface(SS+8) = {SS+8};
  Curve Loop(SS+9) = {SL+ 6,   SL+16,  -(SL+ 8), -(SL+15)}; Plane Surface(SS+9) = {SS+9};
EndFor

// Define horizontal surfaces to be transfinite
Printf("==== SETTING Z SURFACES AS TRANSFINITE ====");
For i_f In{ 1 : 3 }
  If(i_f == 1)  SP = POINT_FLOOR;  SS = SURF_Z_FLOOR;  EndIf
  If(i_f == 2)  SP = POINT_CUBE;   SS = SURF_Z_CUBE;   EndIf
  If(i_f == 3)  SP = POINT_TOP;    SS = SURF_Z_TOP;    EndIf
  For j In{ 1 : 4 }
    For i In{ 1 : 4 }
      FP = (j-1) * 4 + i;
      If(i < 4) If(j < 4)
        Printf("First point for Z surfaces %g", FP);
        SS = SS + 1;
        Transfinite Surface {SS} = {SP+FP, SP+FP+1, SP+FP+5, SP+FP+4};
      EndIf EndIf
    EndFor
  EndFor
EndFor

//----------------------------------
// Define surfaces in "y" direction
//----------------------------------
For i_f In{ 1 : 2 }
  If(i_f == 1)
    SL1 = LINE_H_FLOOR;  SL2 = LINE_H_CUBE;  SL3 = LINE_V_LOW;   SS = SURF_Y_LOW;
  EndIf
  If(i_f == 2)
    SL1 = LINE_H_CUBE;   SL2 = LINE_H_TOP;   SL3 = LINE_V_HIGH;  SS = SURF_Y_HIGH;
  EndIf
  Curve Loop(SS+ 1) = {SL1+01, SL3+01, -(SL2+01), -(SL3+02)};  Plane Surface(SS+ 1) = {SS+ 1};
  Curve Loop(SS+ 2) = {SL1+17, SL3+03, -(SL2+17), -(SL3+02)};  Plane Surface(SS+ 2) = {SS+ 2};
  Curve Loop(SS+ 3) = {SL1+02, SL3+04, -(SL2+02), -(SL3+03)};  Plane Surface(SS+ 3) = {SS+ 3};
  Curve Loop(SS+ 4) = {SL1+03, SL3+05, -(SL2+03), -(SL3+06)};  Plane Surface(SS+ 4) = {SS+ 4};
  Curve Loop(SS+ 5) = {SL1+18, SL3+07, -(SL2+18), -(SL3+06)};  Plane Surface(SS+ 5) = {SS+ 5};
  Curve Loop(SS+ 6) = {SL1+04, SL3+08, -(SL2+04), -(SL3+07)};  Plane Surface(SS+ 6) = {SS+ 6};
  Curve Loop(SS+ 7) = {SL1+05, SL3+09, -(SL2+05), -(SL3+10)};  Plane Surface(SS+ 7) = {SS+ 7};
  Curve Loop(SS+ 8) = {SL1+19, SL3+11, -(SL2+19), -(SL3+10)};  Plane Surface(SS+ 8) = {SS+ 8};
  Curve Loop(SS+ 9) = {SL1+06, SL3+12, -(SL2+06), -(SL3+11)};  Plane Surface(SS+ 9) = {SS+ 9};
  Curve Loop(SS+10) = {SL1+07, SL3+13, -(SL2+07), -(SL3+14)};  Plane Surface(SS+10) = {SS+10};
  Curve Loop(SS+11) = {SL1+20, SL3+15, -(SL2+20), -(SL3+14)};  Plane Surface(SS+11) = {SS+11};
  Curve Loop(SS+12) = {SL1+08, SL3+16, -(SL2+08), -(SL3+15)};  Plane Surface(SS+12) = {SS+12};
EndFor

// Define y surfaces to be transfinite
Printf("==== SETTING Y SURFACES AS TRANSFINITE ====");
For i_f In{ 1 : 2 }
  If(i_f == 1)
    SS = SURF_Y_LOW;   SP1 = POINT_FLOOR;  SP2 = POINT_CUBE;
  EndIf
  If(i_f == 2)
    SS = SURF_Y_HIGH;  SP1 = POINT_CUBE;   SP2 = POINT_TOP;
  EndIf
  For j In{ 1 : 4 }
    For i In{ 1 : 4 }
      FP = (j-1) * 4 + i;
      If(i < 4)
        SS = SS + 1;
        Printf("First point for Y surfaces %g is %g", SS, FP);
        Transfinite Surface {SS} = {SP1+FP, SP1+FP+1, SP2+FP+1, SP2+FP};
      EndIf
    EndFor
  EndFor
EndFor
//+

//----------------------------------
// Define surfaces in "x" direction
//----------------------------------
For i_f In{ 1 : 2 }
  If(i_f == 1)
    SS = SURF_X_LOW;   SL1 = LINE_H_FLOOR;  SL2 = LINE_H_CUBE;  SL3 = LINE_V_LOW;
  EndIf
  If(i_f == 2)
    SS = SURF_X_HIGH;  SL1 = LINE_H_CUBE;    SL2 = LINE_H_TOP;   SL3 = LINE_V_HIGH;
  EndIf
  Curve Loop(SS+ 1) = {SL1+ 9, SL3+ 1, -(SL2+ 9), -(SL3+ 5)};  Plane Surface(SS+ 1) = {SS+ 1};
  Curve Loop(SS+ 2) = {SL1+10, SL3+ 2, -(SL2+10), -(SL3+ 6)};  Plane Surface(SS+ 2) = {SS+ 2};
  Curve Loop(SS+ 3) = {SL1+11, SL3+ 3, -(SL2+11), -(SL3+ 7)};  Plane Surface(SS+ 3) = {SS+ 3};
  Curve Loop(SS+ 4) = {SL1+12, SL3+ 4, -(SL2+12), -(SL3+ 8)};  Plane Surface(SS+ 4) = {SS+ 4};
  Curve Loop(SS+ 5) = {SL1+21, SL3+ 9, -(SL2+21), -(SL3+ 5)};  Plane Surface(SS+ 5) = {SS+ 5};
  Curve Loop(SS+ 6) = {SL1+22, SL3+10, -(SL2+22), -(SL3+ 6)};  Plane Surface(SS+ 6) = {SS+ 6};
  Curve Loop(SS+ 7) = {SL1+23, SL3+11, -(SL2+23), -(SL3+ 7)};  Plane Surface(SS+ 7) = {SS+ 7};
  Curve Loop(SS+ 8) = {SL1+24, SL3+12, -(SL2+24), -(SL3+ 8)};  Plane Surface(SS+ 8) = {SS+ 8};
  Curve Loop(SS+ 9) = {SL1+13, SL3+13, -(SL2+13), -(SL3+ 9)};  Plane Surface(SS+ 9) = {SS+ 9};
  Curve Loop(SS+10) = {SL1+14, SL3+14, -(SL2+14), -(SL3+10)};  Plane Surface(SS+10) = {SS+10};
  Curve Loop(SS+11) = {SL1+15, SL3+15, -(SL2+15), -(SL3+11)};  Plane Surface(SS+11) = {SS+11};
  Curve Loop(SS+12) = {SL1+16, SL3+16, -(SL2+16), -(SL3+12)};  Plane Surface(SS+12) = {SS+12};
EndFor

// Define x surfaces to be transfinite
Printf("==== SETTING X SURFACES AS TRANSFINITE ====");
For i_f In{ 1 : 2 }
  If(i_f == 1)
    SS = SURF_X_LOW;   SP1 = POINT_FLOOR;  SP2 = POINT_CUBE;
  EndIf
  If(i_f == 2)
    SS = SURF_X_HIGH;  SP1 = POINT_CUBE;   SP2 = POINT_TOP;
  EndIf
  For j In{ 1 : 4 }
    For i In{ 1 : 4 }
      FP = (j-1) * 4 + i;
      If(j < 4)
        SS = SS + 1;
        Printf("First point for X surfaces %g is %g", SS, FP);
        Transfinite Surface {SS} = {SP1+FP, SP1+FP+4, SP2+FP+4, SP2+FP};
      EndIf
    EndFor
  EndFor
EndFor

//------------------------
// Recombine all surfaces
//------------------------
Recombine Surface "*";

//==============================================================================
//
// VOLUMES
//
//==============================================================================

SV = 0;

Printf("==== DEFINING VOLUMES ====");
For i_f In{ 1 : 2 }
  If(i_f == 1)
    SSZ1 = SURF_Z_FLOOR;  SSZ2 = SURF_Z_CUBE;  SSY  = SURF_Y_LOW;  SSX  = SURF_X_LOW;
  EndIf
  If(i_f == 2)
    SSZ1 = SURF_Z_CUBE; SSZ2 = SURF_Z_TOP;  SSY  = SURF_Y_HIGH;  SSX  = SURF_X_HIGH;
  EndIf
  For j In{ 1 : 3 }
    For i In{ 1 : 3 }
      SV = SV + 1;
      FF  = (j-1) * 3 + i;    // first face in z direction
      FFX = (j-1) * 4 + i;    // first face in x direction
      Printf("Volume: %g %g %g %g %g %g",  SSZ1+FF, SSZ2+FF, FF+SSY, FF+SSY+3,     FFX+SSX, FFX+SSX+1);
      Surface Loop(SV) = {SSZ1+FF,
                          SSZ2+FF,
                          SSY +FF,
                          SSY +FF  + 3,
                          SSX +FFX,
                          SSX +FFX + 1};
      Volume(SV) = {SV};
    EndFor
  EndFor
EndFor

SV = 0;
For i_f In{ 1 : 2 }
  If(i_f == 1)  SP = POINT_FLOOR;  EndIf
  If(i_f == 2)  SP = POINT_CUBE;   EndIf
  For j In{ 1 : 4 }
    For i In{ 1 : 4 }
      FP = (j-1) * 4 + i;
      If(i < 4) If(j < 4)
        Printf("First point for transfinite volumes %g", FP);
        SV = SV + 1;
//      Transfinite Surface {SS} = {SP+FP, SP+FP+1, SP+FP+5, SP+FP+4};
        Transfinite Volume{SV} = {SP+FP,
                                  SP+FP+1,
                                  SP+FP+5,
                                  SP+FP+4,
                                  POINT_DELTA+SP+FP,
                                  POINT_DELTA+SP+FP+1,
                                  POINT_DELTA+SP+FP+5,
                                  POINT_DELTA+SP+FP+4};
      EndIf EndIf
    EndFor
  EndFor
EndFor

// Make the cube
Recursive Delete {
  Volume{5}; 
}

//==============================================================================
//
// BOUNDARY CONDITONS
//
//==============================================================================

Physical Surface("periodic_x")
  = {Surface In BoundingBox{   -TINY, -HUGE, -HUGE,
                               +TINY, +HUGE, +HUGE},
     Surface In BoundingBox{4*H-TINY, -HUGE, -HUGE,
                            4*H+TINY, +HUGE, +HUGE}};

Physical Surface("periodic_y")
  = {Surface In BoundingBox{-HUGE,    -TINY,  -HUGE,
                            +HUGE,    +TINY,  +HUGE},
     Surface In BoundingBox{-HUGE, 4*H-TINY,  -HUGE,
                            +HUGE, 4*H+TINY,  +HUGE}};

Physical Surface("lower_wall")
  = {Surface In BoundingBox{-HUGE, -HUGE,    -TINY,
                            +HUGE, +HUGE,    +TINY}};
Physical Surface("top_wall")
  = {Surface In BoundingBox{-HUGE, -HUGE,   h-TINY,
                            +HUGE, +HUGE,   h+TINY}};

Physical Surface("cube_walls")
  = {Surface In BoundingBox{2*H-H/2-TINY, 2*H-H/2-TINY,    -TINY,
                            2*H+H/2+TINY, 2*H+H/2+TINY,   H+TINY}};

Physical Volume("FLUID", 613) = {1, 2, 3, 4, 6, 7, 8, 9,
                                 10, 11, 12, 13, 14, 15, 16, 17, 18};
