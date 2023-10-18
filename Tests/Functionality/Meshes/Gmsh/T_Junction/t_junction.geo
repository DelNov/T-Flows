// Gmsh project created on Wed Oct  9 18:07:24 2019

//------------------
// Enable mesh copy
//------------------
Geometry.CopyMeshingMethod = 1;

//--------------------
//
// Set some constants
//
//--------------------

PI   = 3.14159265359;
PI_4 = PI/4.0;

R  =  1;     // outer (true) radius
LJ =  2;     // length of the junction itself
LO = 20;     // lenght of the outlet leg
LI =  6;     // lenght of the inlet legs
r  =  0.7;   // inner radius
NC = 11;     // number of nodes in the core
NB =  7;     // number of nodes in the bundary layer
NJ = 16;     // number of nodes longitudinally
NO = 60;     // number of cell layers towards the outlet
NI = 15;     // number of cell layers towards the outlet

DELTA_MIN =   0.5 * (LJ / (NJ-1) + (LJ-R) / (NJ-1));
MAX_ITERS = 120;

TINY = 1.0e-9;
HUGE = 1.0e+9;

//---------------
//
// Define points
//
//---------------
p = 0;

// Center point at z = 0
Point(p++) = {0, 0, 0};  //  0

// Points defining outer surface at z = 0
Point(p++) = { 0,           -R,           0          };  //  1
Point(p++) = { R*Cos(PI_4), -R*Cos(PI_4), 0          };
Point(p++) = { R,            0,           0          };  //  3
Point(p++) = { R*Cos(PI_4),  R*Cos(PI_4), R*Cos(PI_4)};
Point(p++) = { 0,            R,           R          };  //  5

// Points defining inner surface at z = 0
Point(p++) = { 0,           -r,           0          };  //  6
Point(p++) = { r*Cos(PI_4), -r*Cos(PI_4), 0          };
Point(p++) = { r,            0,           0          };  //  8
Point(p++) = { r*Cos(PI_4),  r*Cos(PI_4), r*Cos(PI_4)};
Point(p++) = { 0,            r,           r          };  // 10

// Center point at z = LJ
Point(p++) = {0, 0, LJ};  // 11

// Points defining outer surface at z = LJ
Point(p++) = { 0,           -R,           LJ          };  // 12
Point(p++) = { R*Cos(PI_4), -R*Cos(PI_4), LJ          };
Point(p++) = { R,            0,           LJ          };  // 14
Point(p++) = { R*Cos(PI_4),  R*Cos(PI_4), LJ          };
Point(p++) = { 0,            R,           LJ          };  // 16

// Points defining inner surface at z = LJ
Point(p++) = { 0,           -r,           LJ          };  // 17
Point(p++) = { r*Cos(PI_4), -r*Cos(PI_4), LJ          };
Point(p++) = { r,            0,           LJ          };  // 19
Point(p++) = { r*Cos(PI_4),  r*Cos(PI_4), LJ          };
Point(p++) = { 0,            r,           LJ          };  // 21

//-------------------------------------
//
// Define connections and their meshes
//
//-------------------------------------
c = 1;

//----------------------
// Connections at z = 0
//----------------------
For n In{ 1 : 4 }
  Ellipse(c) = {n, 0, 11, n+1};
  Transfinite Curve {c++} = NC Using Progression 1;
EndFor

// Lines enclosing the core
For n In{ 6 : 9 }
  Line(c) = {n, n+1};
  Transfinite Curve {c++} = NC Using Progression 1;
EndFor

// Lines from center to the core
For n In{ 6 : 10 : 2 }
  Line(c) = {0, n};
  Transfinite Curve {c++} = NC Using Progression 1;
EndFor

// Lines from core to the outer surface
For n In{ 6 : 10 }
  Line(c) = {n, n-5};
  Transfinite Curve {c++} = NB Using Progression 1;
EndFor

//----------------------
// Connections at z = LJ
//----------------------
For n In{ 12 : 15 }
  Ellipse(c) = {n, 11, 0, n+1};
  Transfinite Curve {c++} = NC Using Progression 1;
EndFor

// Lines enclosing the core
For n In{ 17 : 20 }
  Line(c) = {n, n+1};
  Transfinite Curve {c++} = NC Using Progression 1;
EndFor

// Lines from center to the core
For n In{ 17 : 21 : 2 }
  Line(c) = {11, n};
  Transfinite Curve {c++} = NC Using Progression 1;
EndFor

// Lines from core to the outer surface
For n In{ 17 : 21 }
  Line(c) = {n, n-5};
  Transfinite Curve {c++} = NB Using Progression 1;
EndFor

//-------------------------------------
// Connections between z = 0 and z = LJ
//-------------------------------------
For n In{ 0 : 10 }
  Line(c) = {n, n+11};
  Transfinite Curve {c++} = NJ Using Progression 1;
EndFor

Printf("Defined %g connections", c-1);

//----------------------------------
//
// Create surfaces and their meshes
//
//----------------------------------
s = 1;

//----------
// At z = 0
//----------

// Outer surfaces
Curve Loop(s) = {1, -13, -5, 12};  Plane Surface(s) = {s};  Transfinite Surface {s++};
Curve Loop(s) = {2, -14, -6, 13};  Plane Surface(s) = {s};  Transfinite Surface {s++};
Curve Loop(s) = {3, -15, -7, 14};  Plane Surface(s) = {s};  Transfinite Surface {s++};
Curve Loop(s) = {4, -16, -8, 15};  Plane Surface(s) = {s};  Transfinite Surface {s++};

// In the core
Curve Loop(s) = {5, 6, -10,  9};  Plane Surface(s) = {s};  Transfinite Surface {s++};
Curve Loop(s) = {7, 8, -11, 10};  Plane Surface(s) = {s};  Transfinite Surface {s++};

//----------
// At z = LJ
//----------

// Outer surfaces
Curve Loop(s) = {1+16, -13-16, -5-16, 12+16};  Plane Surface(s) = {s};  Transfinite Surface {s++};
Curve Loop(s) = {2+16, -14-16, -6-16, 13+16};  Plane Surface(s) = {s};  Transfinite Surface {s++};
Curve Loop(s) = {3+16, -15-16, -7-16, 14+16};  Plane Surface(s) = {s};  Transfinite Surface {s++};
Curve Loop(s) = {4+16, -16-16, -8-16, 15+16};  Plane Surface(s) = {s};  Transfinite Surface {s++};

// In the core
Curve Loop(s) = {5+16, 6+16, -10-16,  9+16};  Plane Surface(s) = {s};  Transfinite Surface {s++};
Curve Loop(s) = {7+16, 8+16, -11-16, 10+16};  Plane Surface(s) = {s};  Transfinite Surface {s++};

//----------------------------------
// Surfaces between z = 0 and z = LJ
//----------------------------------

// On the outer radius
Curve Loop(s) = {1, 35, -17, -34};  Surface(s) = {s};  Transfinite Surface {s++};
Curve Loop(s) = {2, 36, -18, -35};  Surface(s) = {s};  Transfinite Surface {s++};
Curve Loop(s) = {3, 37, -19, -36};  Surface(s) = {s};  Transfinite Surface {s++};
Curve Loop(s) = {4, 38, -20, -37};  Surface(s) = {s};  Transfinite Surface {s++};

// Between the core and outer radius
Curve Loop(s) = {12, 34, -28, -39};  Plane Surface(s) = {s};  Transfinite Surface {s++};
Curve Loop(s) = {13, 35, -29, -40};  Plane Surface(s) = {s};  Transfinite Surface {s++};
Curve Loop(s) = {14, 36, -30, -41};  Plane Surface(s) = {s};  Transfinite Surface {s++};
Curve Loop(s) = {15, 37, -31, -42};  Plane Surface(s) = {s};  Transfinite Surface {s++};
Curve Loop(s) = {16, 38, -32, -43};  Plane Surface(s) = {s};  Transfinite Surface {s++};

// Core's outer mantle
Curve Loop(s) = {5, 40, -21, -39};  Plane Surface(s) = {s};  Transfinite Surface {s++};
Curve Loop(s) = {6, 41, -22, -40};  Plane Surface(s) = {s};  Transfinite Surface {s++};
Curve Loop(s) = {7, 42, -23, -41};  Plane Surface(s) = {s};  Transfinite Surface {s++};
Curve Loop(s) = {8, 43, -24, -42};  Plane Surface(s) = {s};  Transfinite Surface {s++};

// Core's inner mantle
Curve Loop(s) = { 9, 39, -25, -33};  Plane Surface(s) = {s};  Transfinite Surface {s++};
Curve Loop(s) = {10, 41, -26, -33};  Plane Surface(s) = {s};  Transfinite Surface {s++};
Curve Loop(s) = {11, 43, -27, -33};  Plane Surface(s) = {s};  Transfinite Surface {s++};

Printf("Defined %g surfaces", s-1);
Recombine Surface { 1 : s-1 };

//----------------
//
// Define volumes
//
//----------------

Surface Loop(1) = {1, 13,  7, 17, 22, 18};  Volume(1) = {1};
Surface Loop(2) = {2, 14,  8, 19, 23, 18};  Volume(2) = {2};
Surface Loop(3) = {5, 26, 11, 27, 22, 23};  Volume(3) = {3};
Surface Loop(4) = {3, 15,  9, 20, 24, 19};  Volume(4) = {4};
Surface Loop(5) = {4, 16, 21, 10, 25, 20};  Volume(5) = {5};
Surface Loop(6) = {6, 28, 12, 25, 27, 24};  Volume(6) = {6};

Transfinite Volume{ 1 : 6 };

//----------------------------
// Copy the basic unit around
//----------------------------

Dilate {{0, 0, 0}, {-1, 1, 1}} {
  Duplicata { Volume{1}; Volume{2}; Volume{3}; }
}

Dilate {{0, 0, 0}, {-1, 1, 1}} {
  Duplicata { Volume{4}; Volume{5}; Volume{6}; }
}

Symmetry {0, -1, 1, 0} {
  Duplicata { Volume{4}; Volume{5}; Volume{6};  
              Volume{123}; Volume{154}; Volume{185}; }
}

Symmetry {0, 0, 1, 0} {
  Duplicata {
              Volume{1};
              Volume{3};
              Volume{2};
              Volume{4};
              Volume{5};
              Volume{6};
              Volume{44};
              Volume{75};
              Volume{106};
              Volume{123};
              Volume{154};
              Volume{185};
              Volume{197};
              Volume{228};
              Volume{259};
              Volume{290};
              Volume{321};
              Volume{352};
            }
}

//-----------------------------------------------
// Progressivelly expand more towards the outlet
//-----------------------------------------------
p = 0.99;  // initial progression
d = 0.05;  // initial increment in progression

L_TARGET = LO-LJ;

// Initial length
l = DELTA_MIN * (1 - p^(NO)) / (1 - p);
Printf("Length with progression %g is %g", p, l);

// Length through iterations
For iter In{ 1 : MAX_ITERS }
  If(l > L_TARGET)  // if current length is larger than the target ...
    p = p - d;      // ... reduce the progression
    l = DELTA_MIN * (1 - p^(NO)) / (1 - p);
    If(l < L_TARGET) d = d * 0.5; EndIf
  EndIf
  If(l < L_TARGET)  // if current length is smaller than the target ...
    p = p + d;      // ... increase the progression
    l = DELTA_MIN * (1 - p^(NO)) / (1 - p);
    If(l > L_TARGET) d = d * 0.5; EndIf
  EndIf
  Printf("Iteration %g; Length with progression %g is %g", iter, p, l);
EndFor

For lay In{ 1 : NO }
  If(lay == 1)
    delta_cur = DELTA_MIN;  // current delta
    z_cur     = LJ;
  Else
    z_cur     += delta_cur;
    delta_cur *= p;
  EndIf
  Extrude {0, 0, delta_cur} {
    Surface{Surface In BoundingBox{-HUGE, -HUGE, z_cur-TINY,
                                    HUGE,  HUGE, z_cur+TINY}};
    Layers {1};
    Recombine;
  }
  Printf("Layer %g; placed at %g", lay, z_cur + delta_cur);
EndFor

//-----------------------------------------------
// Progressivelly expand more towards the inlets
//-----------------------------------------------
p = 0.99;  // initial progression
d = 0.05;  // initial increment in progression

L_TARGET = LI-LJ;

// Initial length
l = DELTA_MIN * (1 - p^(NI)) / (1 - p);
Printf("Length with progression %g is %g", p, l);

// Length through iterations
For iter In{ 1 : MAX_ITERS }
  If(l > L_TARGET)  // if current length is larger than the target ...
    p = p - d;      // ... reduce the progression
    l = DELTA_MIN * (1 - p^(NI)) / (1 - p);
    If(l < L_TARGET) d = d * 0.5; EndIf
  EndIf
  If(l < L_TARGET)  // if current length is smaller than the target ...
    p = p + d;      // ... increase the progression
    l = DELTA_MIN * (1 - p^(NI)) / (1 - p);
    If(l > L_TARGET) d = d * 0.5; EndIf
  EndIf
  Printf("Iteration %g; Length with progression %g is %g", iter, p, l);
EndFor

For lay In{ 1 : NI }
  If(lay == 1)
    delta_cur = DELTA_MIN;  // current delta
    z_cur     = LJ;
  Else
    z_cur     += delta_cur;
    delta_cur *= p;
  EndIf
  Extrude {0, 0, -delta_cur} {
    Surface{Surface In BoundingBox{-HUGE, -HUGE, -z_cur-TINY,
                                   +HUGE, +HUGE, -z_cur+TINY}};
    Layers {1};
    Recombine;
  }
  Printf("Layer %g; placed at %g", lay, - z_cur - delta_cur);
EndFor

For lay In{ 1 : NI }
  If(lay == 1)
    delta_cur = DELTA_MIN;  // current delta
    z_cur     = LJ;
  Else
    z_cur     += delta_cur;
    delta_cur *= p;
  EndIf
  Extrude {0, delta_cur, 0} {
    Surface{Surface In BoundingBox{-HUGE, z_cur-TINY, -HUGE,
                                   +HUGE, z_cur+TINY, +HUGE}};
    Layers {1};
    Recombine;
  }
  Printf("Layer %g; placed at %g", lay, z_cur + delta_cur);
EndFor

//-------------------------
//
// Set boundary conditions
//
//-------------------------

//-------------------
// Walls long z axis
//-------------------
For a In{ 0 : 7 }
  angle_1 =  a    * PI_4;
  angle_2 = (a+1) * PI_4;
  x1  = R * Cos(angle_1);    y1  = R * Sin(angle_1);
  x2  = R * Cos(angle_2);    y2  = R * Sin(angle_2);
  x_1 = Min(x1, x2) - TINY;  y_1 = Min(y1, y2) - TINY;
  x_2 = Max(x1, x2) + TINY;  y_2 = Max(y1, y2) + TINY;
  If(a == 0)
    Physical Surface("outer_wall")
    =  {Surface In BoundingBox{x_1, y_1, -HUGE, x_2, y_2, HUGE}};
  Else
    Physical Surface("outer_wall")
    += {Surface In BoundingBox{x_1, y_1, -HUGE, x_2, y_2, HUGE}};
  EndIf
EndFor

//-------------------
// Walls long y axis
//-------------------
For a In{ 0 : 7 }
  angle_1 =  a    * PI_4;
  angle_2 = (a+1) * PI_4;
  x1  = R * Cos(angle_1);    z1  = R * Sin(angle_1);
  x2  = R * Cos(angle_2);    z2  = R * Sin(angle_2);
  x_1 = Min(x1, x2) - TINY;  z_1 = Min(z1, z2) - TINY;
  x_2 = Max(x1, x2) + TINY;  z_2 = Max(z1, z2) + TINY;
  // Printf("Coords %g %g %g %g", x_1, z_1, x_2, z_2);
  If(a == 0)
    Physical Surface("outer_wall")
    += {Surface In BoundingBox{x_1, -HUGE, z_1, x_2, HUGE, z_2}};
  Else
    Physical Surface("outer_wall")
    += {Surface In BoundingBox{x_1, -HUGE, z_1, x_2, HUGE, z_2}};
  EndIf
EndFor

//--------------------
// Inlets and outlets
//--------------------
Physical Surface("inlet_z") = {Surface In BoundingBox{-HUGE, -HUGE, -LI-TINY,
                                                      +HUGE, +HUGE, -LI+TINY}};
Physical Surface("inlet_y") = {Surface In BoundingBox{-HUGE, LI-TINY, -HUGE,
                                                      +HUGE, LI+TINY, +HUGE}};
Physical Surface("outlet_z") = {Surface In BoundingBox{-HUGE, -HUGE, LO-TINY,
                                                       +HUGE, +HUGE, LO+TINY}};
Physical Volume("fluid") = {Volume In BoundingBox{-HUGE, -HUGE, -HUGE,
                                                  +HUGE, +HUGE, +HUGE}};
