//==============================================================================
//
//
// Definition of constants
//
//
//------------------------------------------------------------------------------
LI =  4.0;
LO = 12.0;
R(1) = 0.6;  // core radius
R(2) = 0.9;  // boundary layer radius
R(3) = 1.0;  // outer radius

// Number of nodes (resolutions for lines)
N_CORE_ARC    = 11;      // number of nodes in the core
N_STREAM_I    = 41;      // number of nodes in the stream-wise direction
N_STREAM_O    = 61;      // number of nodes in the stream-wise direction
N_MIDDLE      =  9;
N_BLAYER      =  7;

// Progression for lines
P_CORE_RAD    = 1.1;
P_STREAM_I    = 1.05;
P_STREAM_O    = 1.05;
P_MIDDLE      = 1.15;
P_BLAYER      = 1.3;

// Connector types (will be used to define resolution and progression)
TYPE_ARC      = 111;  // all arches and ellipses
TYPE_BLAYER   = 222;  // lines in the buffer layer
TYPE_CORE_RAD = 333;  // lines defining the core
TYPE_CORE_TAN = 444;  // lines defining the core
TYPE_MIDDLE   = 555;  // lines connecting core and boundary layer
TYPE_STREAM_I = 666;  // lines in the streamwise directions in inlet
TYPE_STREAM_O = 777;  // lines in the streamwise directions in outlet

// Tolerance for merging nodes
NODE_MERGE_TOL = 1.0e-6;

TINY = 1.0e-9;
HUGE = 1.0e+9;

//==============================================================================
//
//
// Definition of macros
//
//
//------------------------------------------------------------------------------

Macro MergePoints         // stores coppies of the nodes
  cmp = 0;
  For p In{ 1 : pc }
    pnt_copy[p] = p;      // initialize by setting the point as its own copy
    For q In{ 1 : p-1 }
      c[] = Point{p};       // coordinates of point p
      d[] = Point{q};       // coordinates of point q
      If(pnt_copy[p] == p)  // if it hasn't found a copy yet
        If(Abs(c[0]-d[0]) + Abs(c[1]-d[1]) + Abs(c[2]-d[2]) < NODE_MERGE_TOL)
          Printf("Nodes %g and %g are matching", p, q);
          pnt_copy[p] = q;
          cmp++;
        EndIf
      EndIf
      // Printf("Coordinates of point %g are %g %g %g", p, c[0], c[1], c[2]);
    EndFor
  EndFor
  Printf("Number of merged nodes %g!", cmp);
Return

Macro NPointsAroundX  // place eight points in a circle in plane x
  For n In{ nf : nl }
    a = (n-1) * (2.0 * Pi / (nl - nf + 1));
    pc++;  Point(pc) = {x, r * Cos(a), r * Sin(a)};
  EndFor
Return

Macro NPointsMinusX  // place eight points in a circle in plane x
  For n In{ nf : nl }
    a = (n-1) * (2.0 * Pi / (nl - nf + 1));
    pc++;  Point(pc) = {-r * Sin(a), r * Cos(a), r * Sin(a)};
  EndFor
Return

Macro NPointsPlusX  // place eight points in a circle in plane x
  For n In{ nf : nl }
    a = (n-1) * (2.0 * Pi / (nl - nf + 1));
    pc++;  Point(pc) = {+r * Sin(a), r * Cos(a), r * Sin(a)};
  EndFor
Return

Macro NPointsAroundZ  // place eight points in a circle in plane z
  For n In{ nf : nl }
    a = Pi / 2.0 + (n-1) * (2.0 * Pi / (nl - nf + 1));
    pc++;  Point(pc) = {r * Cos(a), r * Sin(a), z};
  EndFor
Return

Macro NLinesAround  // connect eight lines around a circle
  For n In{ nf : nl }
    lc++;
    pnt_1 = pnt_1st + (n - 1);
    pnt_2 = pnt_1 + 1;
    If(n == nl)
      pnt_2 = pnt_1st;
    EndIf
    Printf("Creating line %g around from points: %g, %g", lc, pnt_1, pnt_2);
    Line(lc) = {pnt_copy(pnt_1), pnt_copy(pnt_2)};
    lin_type[lc] = type;
  EndFor
Return

Macro NArcsAround  // connect eight circles around core
  For n In{ nf : nl }
    lc++;
    pnt_1 = pnt_1st + (n - 1);
    pnt_2 = pnt_1 + 1;
    If(n == nl)
      pnt_2 = pnt_1st;
    EndIf
    Printf("Creating arc %g from points: %g, %g", lc, pnt_1, pnt_2);
    Circle(lc) = {pnt_copy(pnt_1),
                  pnt_copy(pnt_cent),
                  pnt_copy(pnt_2)};
    lin_type[lc] = type;
  EndFor
Return

Macro NEllipseArcsAround  // connect eight circles around core
  For n In{ nf : nl }
    lc++;
    pnt_1 = pnt_1st + (n - 1);
    pnt_2 = pnt_1 + 1;
    If(n == nl)
      pnt_2 = pnt_1st;
    EndIf
    Printf("Creating ellips arc %g from points: %g, %g", lc, pnt_1, pnt_2);
    Ellipse(lc) = {pnt_copy(pnt_1),
                   pnt_copy(pnt_cent), pnt_copy(pnt_axis),
                   pnt_copy(pnt_2)};
    lin_type[lc] = type;
  EndFor
Return

Macro NLinesFromCenter  // connect four lines from a central point
  For n In{ nf : nl : 2 }
    lc++;
    pnt_1 = pnt_1st + (n - 1);
    Printf("Creating line %g from from points: %g, %g", lc, pnt_1, pnt_cent);
    Line(lc) = {pnt_copy(pnt_1), pnt_copy(pnt_cent)};
    lin_type[lc] = type;
  EndFor
Return

Macro NLinesBetweenLayers  // connect eight lines from a central point
  For n In{ nf : nl }
    lc++;
    pnt_1 = pnt_1st  + (n - 1);
    pnt_2 = pnt_2nd + (n - 1);
    Printf("Creating line %g between layers from points: %g, %g",
           lc, pnt_1, pnt_2);
    Line(lc) = {pnt_copy(pnt_1), pnt_copy(pnt_2)};
    lin_type[lc] = type;
  EndFor
Return

Macro NLinesBetweenLayersNeg  // connect eight lines from a central point
  For n In{ nf : nl }
    lc++;
    pnt_1 = pnt_1st  - (n - 1);
    pnt_2 = pnt_2nd + (n - 1);
    Printf("Creating line %g between layers from points: %g, %g",
           lc, pnt_1, pnt_2);
    Line(lc) = {pnt_copy(pnt_1), pnt_copy(pnt_2)};
    lin_type[lc] = type;
  EndFor
Return

Macro SurfaceCoreOrthog
  For c In{ cf : cl }
    sc++;                         // increase surface counter
    con_1 = con_first + (c - 1);
    con_2 = con_1 + 1;
    If(c == cl)
      con_2 = con_first;
    EndIf
    con_4 = con_second + (c - 1) * 2;
    con_3 = con_4 + 1;
    Printf("Creating core surface %g from lines: %g, %g %g, %g", 
            sc, con_1, -con_2, -con_3, -con_4);
    Curve Loop(sc) = {con_1, -con_2, -con_3, -con_4};
    Plane Surface(sc) = {sc};
  EndFor
Return

Macro SurfaceOrthog
  For c In{ cf : cl }
    sc++;                         // increase surface counter
    con_1 = con_first + (c - 1);
    con_3 = con_1 + 1;            // I hope this 96 is general enough :-(
    If(c == cl)
      con_3 = con_first;
    EndIf
    con_2 = con_1 - 48;
    con_4 = con_1 + 48;
    Printf("Creating middle surface %g from lines: %g, %g %g, %g", 
            sc, con_1, con_2, -con_3, -con_4);
    Curve Loop(sc) = {con_1, con_2, -con_3, -con_4};
    Plane Surface(sc) = {sc};
  EndFor
Return

//==============================================================================
//
//
// Start of the algorithm
//
//
//------------------------------------------------------------------------------

// Initialize counters
pc = 0;  // point count
lc = 0;  // line count
sc = 0;  // surface count
vc = 0;  // volume count

// Initialize list of line types and point coppies
lin_type[] = {};
pnt_copy[] = {};

//--------
//
// Points
//
//--------

// Create points along the main axes (1-4)
pc++;  Point(pc) = {-LI, 0, 0};
pc++;  Point(pc) = {+LO, 0, 0};
pc++;  Point(pc) = {0,   0, LI};
pc++;  Point(pc) = {0,   0, 0};

// Points in the core
nf = 1;  nl = 8;
x = -LI;  r = R(1);  Call NPointsAroundX; // point  5 on
x = +LO;  r = R(1);  Call NPointsAroundX; // point 13 on
z = +LI;  r = R(1);  Call NPointsAroundZ; // point 21 on
x = 0.0;  r = R(1);  Call NPointsAroundX; // point 29 on
r = R(1);  Call NPointsMinusX;            // point 37 on
r = R(1);  Call NPointsPlusX;             // point 45 on

// Points in the boundary layer
nf = 1;  nl = 8;
x = -LI;  r = R(2);  Call NPointsAroundX;  // point 53 on
x = +LO;  r = R(2);  Call NPointsAroundX;  // point 61 on
z = +LI;  r = R(2);  Call NPointsAroundZ;  // point 69 on
x = 0.0;  r = R(2);  Call NPointsAroundX;  // point 77 on
r = R(2);  Call NPointsMinusX;             // point 85 on
r = R(2);  Call NPointsPlusX;              // point 93 on

// Points at the cylinder wall
nf = 1;  nl = 8;
x = -LI;  r = R(3);  Call NPointsAroundX;  // point 101 on
x = +LO;  r = R(3);  Call NPointsAroundX;  // point 109 on
z = +LI;  r = R(3);  Call NPointsAroundZ;  // point 117 on
x = 0.0;  r = R(3);  Call NPointsAroundX;  // point 125 on
r = R(3);  Call NPointsMinusX;             // point 133 on
r = R(3);  Call NPointsPlusX;              // point 141 on

//-------------------------------------------------------------
// Try merging all the points after defining their coordinates
//-------------------------------------------------------------
Call MergePoints;

//------------
//
// Connectors
//
//------------

// Connectors are created in the following way:
// 1 - center in the core in the radial direction
// 2 - connectors wrapping the core (on the border with middle layer)
// 3 - middle radial connectors
// 4 - connectors wrapping the bundary layers (arches and ellipses)
// 5 - boundary layer radial connectors
// 6 - connectors for the outer walls (arches and ellipses)
//
// Planes are ordered in the following way
// 1 - at x min
// 2 - at x max
// 3 - at z max
// 4 - at x = 0
// 5 - at x = 0 rotated by 45 degrees anticlockwise
// 6 - at x = 0 rotated by 45 degrees clockwise

//--------------------------------------------
// Center in the core in the radial direction
//--------------------------------------------

// Lines in radial direction of the cores
type = TYPE_CORE_RAD;

nf = 1;  nl = 8;
pnt_cent = 1;  pnt_1st =  5;  Call NLinesFromCenter;  // lines  1 -  4
pnt_cent = 2;  pnt_1st = 13;  Call NLinesFromCenter;  // lines  5 -  8
pnt_cent = 3;  pnt_1st = 21;  Call NLinesFromCenter;  // lines  9 - 12
pnt_cent = 4;  pnt_1st = 29;  Call NLinesFromCenter;  // lines 13 - 16
pnt_cent = 4;  pnt_1st = 37;  Call NLinesFromCenter;  // lines 17 - 20
pnt_cent = 4;  pnt_1st = 45;  Call NLinesFromCenter;  // lines 21 - 24

//----------------------------------------------------------------
// Connectors wrapping the core (on the border with middle layer)
//----------------------------------------------------------------
type = TYPE_CORE_TAN;

nf = 1;  nl = 8;
pnt_1st =  5;  Call NLinesAround;  // lines 25 - 32
pnt_1st = 13;  Call NLinesAround;  // lines 33 - 40
pnt_1st = 21;  Call NLinesAround;  // lines 41 - 48
pnt_1st = 29;  Call NLinesAround;  // lines 49 - 56
pnt_1st = 37;  Call NLinesAround;  // lines 57 - 64
pnt_1st = 45;  Call NLinesAround;  // lines 65 - 72

//--------------------------
// Middle radial connectors
//--------------------------
type = TYPE_MIDDLE;

nf = 1; nl = 8;
pnt_1st = 53;  pnt_2nd =  5;  Call NLinesBetweenLayers;
pnt_1st = 61;  pnt_2nd = 13;  Call NLinesBetweenLayers;
pnt_1st = 69;  pnt_2nd = 21;  Call NLinesBetweenLayers;
pnt_1st = 77;  pnt_2nd = 29;  Call NLinesBetweenLayers;
pnt_1st = 85;  pnt_2nd = 37;  Call NLinesBetweenLayers;
pnt_1st = 93;  pnt_2nd = 45;  Call NLinesBetweenLayers;

//--------------------------------------------------------------
// Connectors wrapping the bundary layers (arches and ellipses)
//--------------------------------------------------------------
type = TYPE_ARC;

nf = 1;  nl = 8;
pnt_1st = 53;  pnt_cent = 1;  Call NArcsAround;
pnt_1st = 61;  pnt_cent = 2;  Call NArcsAround;
pnt_1st = 69;  pnt_cent = 3;  Call NArcsAround;
pnt_1st = 77;  pnt_cent = 4;  Call NArcsAround;
pnt_1st = 85;  pnt_cent = 4;  pnt_axis = 87;  Call NEllipseArcsAround;
pnt_1st = 93;  pnt_cent = 4;  pnt_axis = 95;  Call NEllipseArcsAround;

//----------------------------------
// Boundary layer radial connectors
//----------------------------------
type = TYPE_BLAYER;

nf = 1; nl = 8;
pnt_1st = 101;  pnt_2nd = 53;  Call NLinesBetweenLayers;
pnt_1st = 109;  pnt_2nd = 61;  Call NLinesBetweenLayers;
pnt_1st = 117;  pnt_2nd = 69;  Call NLinesBetweenLayers;
pnt_1st = 125;  pnt_2nd = 77;  Call NLinesBetweenLayers;
pnt_1st = 133;  pnt_2nd = 85;  Call NLinesBetweenLayers;
pnt_1st = 141;  pnt_2nd = 93;  Call NLinesBetweenLayers;

//--------------------------------------------------------------
// Connectors wrapping the bundary layers (arches and ellipses)
//--------------------------------------------------------------
type = TYPE_ARC;

nf = 1;  nl = 8;
pnt_1st = 101;  pnt_cent = 1;  Call NArcsAround;
pnt_1st = 109;  pnt_cent = 2;  Call NArcsAround;
pnt_1st = 117;  pnt_cent = 3;  Call NArcsAround;
pnt_1st = 125;  pnt_cent = 4;  Call NArcsAround;
pnt_1st = 133;  pnt_cent = 4;  pnt_axis = 135;  Call NEllipseArcsAround;
pnt_1st = 141;  pnt_cent = 4;  pnt_axis = 143;  Call NEllipseArcsAround;

//------------------------------------
// Lines in the streamwise directions
//------------------------------------

// Towards xmin: from axis, through core, boundary layers and pipe walls
type = TYPE_STREAM_I;
lc++;  Line(lc) = {4, 1};  lin_type(lc) = type;  // axis
nf = 1; nl = 4;  pnt_1st = 37;  pnt_2nd =  5;  Call NLinesBetweenLayers;
nf = 5; nl = 8;  pnt_1st = 29;  pnt_2nd =  5;  Call NLinesBetweenLayers;
nf = 1; nl = 4;  pnt_1st = 85;  pnt_2nd = 53;  Call NLinesBetweenLayers;
nf = 5; nl = 8;  pnt_1st = 77;  pnt_2nd = 53;  Call NLinesBetweenLayers;
nf = 1; nl = 4;  pnt_1st =133;  pnt_2nd =101;  Call NLinesBetweenLayers;
nf = 5; nl = 8;  pnt_1st =125;  pnt_2nd =101;  Call NLinesBetweenLayers;

// Towards xmax: from axis, through core, boundary layers and pipe walls
type = TYPE_STREAM_O;
lc++;  Line(lc) = {4, 2};  lin_type(lc) = type;  // axis
nf = 1; nl = 4;  pnt_1st = 45;  pnt_2nd = 13;  Call NLinesBetweenLayers;
nf = 5; nl = 8;  pnt_1st = 29;  pnt_2nd = 13;  Call NLinesBetweenLayers;
nf = 1; nl = 4;  pnt_1st = 93;  pnt_2nd = 61;  Call NLinesBetweenLayers;
nf = 5; nl = 8;  pnt_1st = 77;  pnt_2nd = 61;  Call NLinesBetweenLayers;
nf = 1; nl = 4;  pnt_1st =141;  pnt_2nd =109;  Call NLinesBetweenLayers;
nf = 5; nl = 8;  pnt_1st =125;  pnt_2nd =109;  Call NLinesBetweenLayers;

// Towards zmax
type = TYPE_STREAM_I;
lc++;  Line(lc) = {4, 3};  lin_type(lc) = type;  // axis
nf = 1; nl = 4;  pnt_1st = 37;  pnt_2nd = 21;  Call NLinesBetweenLayers;
nf = 5; nl = 8;  pnt_1st = 53;  pnt_2nd = 21;  Call NLinesBetweenLayersNeg;
nf = 1; nl = 4;  pnt_1st = 85;  pnt_2nd = 69;  Call NLinesBetweenLayers;
nf = 5; nl = 8;  pnt_1st =101;  pnt_2nd = 69;  Call NLinesBetweenLayersNeg;
nf = 1; nl = 4;  pnt_1st =133;  pnt_2nd =117;  Call NLinesBetweenLayers;
nf = 5; nl = 8;  pnt_1st =149;  pnt_2nd =117;  Call NLinesBetweenLayersNeg;

//-------------------------------------------
// Set all connectors to be transfinite with
//     proper resolution and progression
//-------------------------------------------
For l In{ 1 : lc }
  If(lin_type(l) == TYPE_ARC)
    Transfinite Curve {l} = N_CORE_ARC  Using Progression 1.0;
  EndIf
  If(lin_type(l) == TYPE_BLAYER)
    Transfinite Curve {l} = N_BLAYER    Using Progression P_BLAYER;
  EndIf
  If(lin_type(l) == TYPE_CORE_RAD)
    Transfinite Curve {l} = N_CORE_ARC  Using Progression P_CORE_RAD;
  EndIf
  If(lin_type(l) == TYPE_CORE_TAN)
    Transfinite Curve {l} = N_CORE_ARC  Using Progression 1.0;
  EndIf
  If(lin_type(l) == TYPE_MIDDLE)
    Transfinite Curve {l} = N_MIDDLE    Using Progression P_MIDDLE;
  EndIf
  If(lin_type(l) == TYPE_STREAM_I)
    Transfinite Curve {l} = N_STREAM_I  Using Progression P_STREAM_I;
  EndIf
  If(lin_type(l) == TYPE_STREAM_O)
    Transfinite Curve {l} = N_STREAM_O  Using Progression P_STREAM_O;
  EndIf
EndFor

//----------
//
// Surfaces
//
//----------

//-----------------------------
// Surfaces orthogonal to axes
//-----------------------------

// Surfaces in the core center
cf = 1;  cl = 4;
con_first =  1;  con_second = 25;  Call SurfaceCoreOrthog;  // surfs  1 -  4
con_first =  5;  con_second = 33;  Call SurfaceCoreOrthog;  // surfs  5 -  8
con_first =  9;  con_second = 41;  Call SurfaceCoreOrthog;  // surfs  9 - 12
con_first = 13;  con_second = 49;  Call SurfaceCoreOrthog;  // surfs 13 - 16
con_first = 17;  con_second = 57;  Call SurfaceCoreOrthog;  // surfs 17 - 20
con_first = 21;  con_second = 65;  Call SurfaceCoreOrthog;  // surfs 21 - 24

// Surfaces in the middle
cf = 1;  cl = 8;
con_first =  73;  Call SurfaceOrthog;  // surfs 25 - 32
con_first =  81;  Call SurfaceOrthog;  // surfs 33 - 40
con_first =  89;  Call SurfaceOrthog;  // surfs 41 - 48
con_first =  97;  Call SurfaceOrthog;  // surfs 49 - 56
con_first = 105;  Call SurfaceOrthog;  // surfs 57 - 64
con_first = 113;  Call SurfaceOrthog;  // surfs 65 - 72

// Surfaces in the boundary layers
con_first = 169;  Call SurfaceOrthog;  // surfs  73 -  80
con_first = 177;  Call SurfaceOrthog;  // surfs  81 -  88
con_first = 185;  Call SurfaceOrthog;  // surfs  89 -  96
con_first = 193;  Call SurfaceOrthog;  // surfs  97 - 104
con_first = 201;  Call SurfaceOrthog;  // surfs 105 - 112
con_first = 209;  Call SurfaceOrthog;  // surfs 113 - 120

//---------------------
// Streamwise surfaces
//---------------------

// Surfaces 121 - 124; from axis towards core wrap at xmin
// 265 is axis connecting center to xmin
sc++;  Curve Loop(sc) = {1, -265, -17, 266};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {2, -265, -18, 268};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {3, -265, -15, 270};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {4, -265, -16, 272};  Plane Surface(sc) = {sc};

// Surfaces 125 - 128; from axis towards core wrap at xmax
// 290 is axis connecting center to xmin
sc++;  Curve Loop(sc) = {5, -290, -21, 291};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {6, -290, -22, 293};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {7, -290, -15, 295};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {8, -290, -16, 297};  Plane Surface(sc) = {sc};

// Surfaces 129 - 132; from axis towards core wrap at zmax
// 315 is axis connecting center to xmin
sc++;  Curve Loop(sc) = { 9, -315, -17, 316};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {10, -315, -18, 318};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {11, -315, -23, 320};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {12, -315, -22, 322};  Plane Surface(sc) = {sc};

// Core wrapping surfaces at xmin
sc++;  Curve Loop(sc) = {25, -267, -57, 266};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {26, -268, -58, 267};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {27, -269, -59, 268};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {28, -270, -60, 269};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {29, -271, -53, 270};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {30, -272, -54, 271};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {31, -273, -55, 272};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {32, -266, -56, 273};  Plane Surface(sc) = {sc};

// Core wrapping surfaces at xmax
sc++;  Curve Loop(sc) = {33, -292, -65, 291};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {34, -293, -66, 292};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {35, -294, -67, 293};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {36, -295, -68, 294};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {37, -296, -53, 295};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {38, -297, -54, 296};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {39, -298, -55, 297};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {40, -291, -56, 298};  Plane Surface(sc) = {sc};

// Core wrapping surfaces at zmax
sc++;  Curve Loop(sc) = {41, -317, -57, 316};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {42, -318, -58, 317};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {43, -319, -59, 318};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {44, -320, -60, 319};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {45, -321,  68, 320};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {46, -322,  67, 321};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {47, -323,  66, 322};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {48, -316,  65, 323};  Plane Surface(sc) = {sc};

// Surfaces in the middle at xmin leg
sc++;  Curve Loop(sc) = {73, -266, -105, 274};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {74, -267, -106, 275};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {75, -268, -107, 276};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {76, -269, -108, 277};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {77, -270, -101, 278};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {78, -271, -102, 279};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {79, -272, -103, 280};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {80, -273, -104, 281};  Plane Surface(sc) = {sc};

// Surfaces in the middle at xmin leg
sc++;  Curve Loop(sc) = {81, -291, -113, 299};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {82, -292, -114, 300};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {83, -293, -115, 301};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {84, -294, -116, 302};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {85, -295, -109, 303};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {86, -296, -102, 304};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {87, -297, -103, 305};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {88, -298, -104, 306};  Plane Surface(sc) = {sc};

// Surfaces in the middle at zmax leg
sc++;  Curve Loop(sc) = {89, -316, -105, 324};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {90, -317, -106, 325};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {91, -318, -107, 326};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {92, -319, -108, 327};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {93, -320, -101, 328};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {94, -321, -116, 329};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {95, -322, -115, 330};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {96, -323, -114, 331};  Plane Surface(sc) = {sc};

// Wrapping the boundary layer at xmin
sc++;  Curve Loop(sc) = {121, -275, -153, 274};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {122, -276, -154, 275};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {123, -277, -155, 276};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {124, -278, -156, 277};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {125, -279, -149, 278};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {126, -280, -150, 279};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {127, -281, -151, 280};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {128, -274, -152, 281};  Surface(sc) = {sc};

// Wrapping the boundary layer at xmax
sc++;  Curve Loop(sc) = {129, -300, -161, 299};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {130, -301, -162, 300};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {131, -302, -163, 301};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {132, -303, -164, 302};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {133, -304, -149, 303};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {134, -305, -150, 304};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {135, -306, -151, 305};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {136, -299, -152, 306};  Surface(sc) = {sc};

// Wrapping the boundary layer at zmax
sc++;  Curve Loop(sc) = {137, -325, -153, 324};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {138, -326, -154, 325};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {139, -327, -155, 326};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {140, -328, -156, 327};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {141, -329,  164, 328};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {142, -330,  163, 329};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {143, -331,  162, 330};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {144, -324,  161, 331};  Surface(sc) = {sc};

// In boundary layers at xmin
sc++;  Curve Loop(sc) = {169, -274, -193, 282};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {170, -275, -202, 283};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {171, -276, -203, 284};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {172, -277, -204, 285};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {173, -278, -197, 286};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {174, -279, -198, 287};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {175, -280, -199, 288};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {176, -281, -200, 289};  Plane Surface(sc) = {sc};

// In boundary layers at xmax
sc++;  Curve Loop(sc) = {177, -299, -193, 307};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {178, -300, -210, 308};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {179, -301, -211, 309};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {180, -302, -212, 310};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {181, -303, -197, 311};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {182, -304, -198, 312};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {183, -305, -199, 313};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {184, -306, -200, 314};  Plane Surface(sc) = {sc};

// In boundary layers at zmax
sc++;  Curve Loop(sc) = {185, -324, -193, 332};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {186, -325, -202, 333};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {187, -326, -203, 334};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {188, -327, -204, 335};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {189, -328, -197, 336};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {190, -329, -212, 337};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {191, -330, -211, 338};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {192, -331, -210, 339};  Plane Surface(sc) = {sc};

Printf("Started defining outer walls with surface: %g", sc+1);
wall_first = sc + 1;

// Outer walls at xmin
sc++;  Curve Loop(sc) = {217, -283, -249, 282};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {218, -284, -250, 283};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {219, -285, -251, 284};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {220, -286, -252, 285};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {221, -287, -245, 286};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {222, -288, -246, 287};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {223, -289, -247, 288};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {224, -282, -248, 289};  Surface(sc) = {sc};

// Outer walls at xmax
sc++;  Curve Loop(sc) = {225, -308, -257, 307};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {226, -309, -258, 308};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {227, -310, -259, 309};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {228, -311, -260, 310};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {229, -312, -245, 311};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {230, -313, -246, 312};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {231, -314, -247, 313};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {232, -307, -248, 314};  Surface(sc) = {sc};

// Outer walls at zmax
sc++;  Curve Loop(sc) = {233, -333, -249, 332};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {234, -334, -250, 333};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {235, -335, -251, 334};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {236, -336, -252, 335};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {237, -337,  260, 336};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {238, -338,  259, 337};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {239, -339,  258, 338};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {240, -332,  257, 339};  Surface(sc) = {sc};

Printf("Finished defining outer walls with surface: %g", sc);
wall_last = sc;

//-------------------------------------------------------
// Set all surfaces to be transfinite
// (Recombine them at once to prevent nodes from moving
//-------------------------------------------------------
For s In{ 1 : sc }
  Transfinite Surface {s};
  Recombine Surface{s};
EndFor

//-----------------------------------------------------------
// These surfaces were created automatically
// here I might be brave enough to delete them
//-----------------------------------------------------------
Recursive Delete {
  Surface{13}; Surface{14};  // core surfaces in the central orthogonal set
  Surface{19}; Surface{20};  // core surfaces in the central anticlockw set
  Surface{23}; Surface{24};  // core surfaces in the central clockwise set
  Surface{49}; Surface{50}; Surface{51}; Surface{52};  // middle orthogonal set
  Surface{61}; Surface{62}; Surface{63}; Surface{64};  // middle anticlockwise
  Surface{72}; Surface{71}; Surface{70}; Surface{69};  // middle clockwise
  Surface{97}; Surface{98}; Surface{99}; Surface{100}; // blayer orthogonal
  Surface{109}; Surface{110}; Surface{111}; Surface{112};  // bl anticlock
  Surface{120}; Surface{119}; Surface{118}; Surface{117};  // bl clockwise
}

//---------
//
// Volumes
//
//---------

// Volumes in the core at xmin
vc++;  Surface Loop(vc) = { 1, 17, 121, 133, 134, 122};  Volume(vc) = {vc};  Transfinite Volume{vc} = {1, 4, 29, 5, 7, 39, 38, 6};
vc++;  Surface Loop(vc) = { 2, 18, 122, 135, 136, 123};  Volume(vc) = {vc};  Transfinite Volume{vc} = {9, 33, 4, 1, 8, 40, 39, 7};
vc++;  Surface Loop(vc) = { 3, 15, 123, 137, 138, 124};  Volume(vc) = {vc};  Transfinite Volume{vc} = {10, 34, 35, 11, 9, 33, 4, 1};
vc++;  Surface Loop(vc) = { 4, 16, 124, 139, 140, 121};  Volume(vc) = {vc};  Transfinite Volume{vc} = {11, 35, 36, 12, 1, 4, 29, 5};

// Volumes in the core at xmax
vc++;  Surface Loop(vc) = { 5, 21, 125, 141, 142, 126};  Volume(vc) = {vc};  Transfinite Volume{vc} = {4, 2, 13, 29, 47, 15, 14, 46};
vc++;  Surface Loop(vc) = { 6, 22, 126, 143, 144, 127};  Volume(vc) = {vc};  Transfinite Volume{vc} = {33, 17, 2, 4, 48, 16, 15, 47};
vc++;  Surface Loop(vc) = { 7, 15, 127, 145, 146, 128};  Volume(vc) = {vc};  Transfinite Volume{vc} = {34, 18, 19, 35, 33, 17, 2, 4};
vc++;  Surface Loop(vc) = { 8, 16, 128, 147, 148, 125};  Volume(vc) = {vc};  Transfinite Volume{vc} = {35, 19, 20, 36, 4, 2, 13, 29};

// Volumes in the core at zmax
vc++;  Surface Loop(vc) = { 9, 17, 129, 149, 150, 130};  Volume(vc) = {vc};  Transfinite Volume{vc} = {39, 4, 29, 38, 23, 3, 21, 22};
vc++;  Surface Loop(vc) = {10, 18, 130, 151, 152, 131};  Volume(vc) = {vc};  Transfinite Volume{vc} = {40, 33, 4, 39, 24, 25, 3, 23};
vc++;  Surface Loop(vc) = {11, 22, 131, 153, 154, 132};  Volume(vc) = {vc};  Transfinite Volume{vc} = {33, 48, 47, 4, 25, 26, 27, 3};
vc++;  Surface Loop(vc) = {12, 21, 132, 155, 156, 129};  Volume(vc) = {vc};  Transfinite Volume{vc} = {4, 47, 46, 29, 3, 27, 28, 21};

Printf("Core finished with volume %g", vc);

// Volumes in the middle at xmin
vc++;  Surface Loop(vc) = {25, 57, 133, 157, 181, 158};  Volume(vc) = {vc};  Transfinite Volume{vc} = {5, 29, 77, 53, 6, 38, 86, 54};
vc++;  Surface Loop(vc) = {26, 58, 134, 158, 182, 159};  Volume(vc) = {vc};  Transfinite Volume{vc} = {7, 39, 38, 6, 55, 87, 86, 54};
vc++;  Surface Loop(vc) = {27, 59, 135, 159, 183, 160};  Volume(vc) = {vc};  Transfinite Volume{vc} = {8, 40, 39, 7, 56, 88, 87, 55};
vc++;  Surface Loop(vc) = {28, 60, 136, 160, 184, 161};  Volume(vc) = {vc};  Transfinite Volume{vc} = {57, 81, 33, 9, 56, 88, 40, 8};
vc++;  Surface Loop(vc) = {29, 53, 137, 161, 185, 162};  Volume(vc) = {vc};  Transfinite Volume{vc} = {58, 82, 34, 10, 57, 81, 33, 9};
vc++;  Surface Loop(vc) = {30, 54, 138, 162, 186, 163};  Volume(vc) = {vc};  Transfinite Volume{vc} = {58, 82, 83, 59, 10, 34, 35, 11};
vc++;  Surface Loop(vc) = {31, 55, 139, 163, 187, 164};  Volume(vc) = {vc};  Transfinite Volume{vc} = {59, 83, 84, 60, 11, 35, 36, 12};
vc++;  Surface Loop(vc) = {32, 56, 140, 164, 188, 157};  Volume(vc) = {vc};  Transfinite Volume{vc} = {12, 36, 84, 60, 5, 29, 77, 53};

// Volumes in the middle at xmax
vc++;  Surface Loop(vc) = {33, 65, 165, 189, 166, 141};  Volume(vc) = {vc};  Transfinite Volume{vc} = {29, 13, 61, 77, 46, 14, 62, 94};
vc++;  Surface Loop(vc) = {34, 66, 166, 190, 167, 142};  Volume(vc) = {vc};  Transfinite Volume{vc} = {47, 15, 14, 46, 95, 63, 62, 94};
vc++;  Surface Loop(vc) = {35, 67, 167, 191, 168, 143};  Volume(vc) = {vc};  Transfinite Volume{vc} = {48, 16, 15, 47, 96, 64, 63, 95};
vc++;  Surface Loop(vc) = {36, 68, 168, 192, 169, 144};  Volume(vc) = {vc};  Transfinite Volume{vc} = {81, 65, 17, 33, 96, 64, 16, 48};
vc++;  Surface Loop(vc) = {37, 53, 169, 193, 170, 145};  Volume(vc) = {vc};  Transfinite Volume{vc} = {82, 66, 18, 34, 81, 65, 17, 33};
vc++;  Surface Loop(vc) = {38, 54, 170, 194, 171, 146};  Volume(vc) = {vc};  Transfinite Volume{vc} = {82, 66, 18, 34, 83, 67, 19, 35};
vc++;  Surface Loop(vc) = {39, 55, 171, 195, 172, 147};  Volume(vc) = {vc};  Transfinite Volume{vc} = {83, 67, 68, 84, 35, 19, 20, 36};
vc++;  Surface Loop(vc) = {40, 56, 172, 196, 165, 148};  Volume(vc) = {vc};  Transfinite Volume{vc} = {36, 20, 68, 84, 29, 13, 61, 77};

// Volumes in the middle at zmax
vc++;  Surface Loop(vc) = {41, 57, 149, 173, 197, 174};  Volume(vc) = {vc};  Transfinite Volume{vc} = {38, 29, 77, 86, 22, 21, 69, 70};
vc++;  Surface Loop(vc) = {42, 58, 150, 174, 198, 175};  Volume(vc) = {vc};  Transfinite Volume{vc} = {87, 39, 38, 86, 71, 23, 22, 70};
vc++;  Surface Loop(vc) = {43, 59, 151, 175, 199, 176};  Volume(vc) = {vc};  Transfinite Volume{vc} = {88, 40, 39, 87, 72, 24, 23, 71};
vc++;  Surface Loop(vc) = {44, 60, 152, 176, 200, 177};  Volume(vc) = {vc};  Transfinite Volume{vc} = {72, 73, 25, 24, 88, 81, 33, 40};
vc++;  Surface Loop(vc) = {45, 68, 153, 177, 201, 178};  Volume(vc) = {vc};  Transfinite Volume{vc} = {81, 96, 48, 33, 73, 74, 26, 25};
vc++;  Surface Loop(vc) = {46, 67, 154, 178, 202, 179};  Volume(vc) = {vc};  Transfinite Volume{vc} = {96, 95, 47, 48, 74, 75, 27, 26};
vc++;  Surface Loop(vc) = {47, 66, 155, 179, 203, 180};  Volume(vc) = {vc};  Transfinite Volume{vc} = {47, 95, 94, 46, 27, 75, 76, 28};
vc++;  Surface Loop(vc) = {48, 65, 156, 180, 204, 173};  Volume(vc) = {vc};  Transfinite Volume{vc} = {29, 46, 94, 77, 21, 28, 76, 69};

Printf("Middle finished with volume %g", vc);

// Volumes in the boundary layer at xmin (like middle at xmin +52 except second column in the lower part?)
vc++;  Surface Loop(vc) = {73, 105, 181, 205, 229, 206};  Volume(vc) = {vc};  Transfinite Volume{vc} = {53, 77, 125, 101, 54, 86, 134, 102};
vc++;  Surface Loop(vc) = {74, 106, 182, 206, 230, 207};  Volume(vc) = {vc};  Transfinite Volume{vc} = {55, 87, 86, 54, 103, 135, 134, 102};
vc++;  Surface Loop(vc) = {75, 107, 183, 207, 231, 208};  Volume(vc) = {vc};  Transfinite Volume{vc} = {56, 88, 87, 55, 104, 136, 135, 103};
vc++;  Surface Loop(vc) = {76, 108, 184, 208, 232, 209};  Volume(vc) = {vc};  Transfinite Volume{vc} = {57, 81, 88, 56, 105, 129, 136, 104};
vc++;  Surface Loop(vc) = {77, 101, 185, 209, 233, 210};  Volume(vc) = {vc};  Transfinite Volume{vc} = {106, 130, 82, 58, 105, 129, 81, 57};
vc++;  Surface Loop(vc) = {78, 102, 186, 210, 234, 211};  Volume(vc) = {vc};  Transfinite Volume{vc} = {106, 130, 131, 107, 58, 82, 83, 59};
vc++;  Surface Loop(vc) = {79, 103, 187, 211, 235, 212};  Volume(vc) = {vc};  Transfinite Volume{vc} = {107, 131, 132, 108, 59, 83, 84, 60};
vc++;  Surface Loop(vc) = {80, 104, 188, 212, 236, 205};  Volume(vc) = {vc};  Transfinite Volume{vc} = {60, 84, 132, 108, 53, 77, 125, 101};

// Volumes in the boundary layer at xmax
vc++;  Surface Loop(vc) = {81, 113, 189, 213, 237, 214};  Volume(vc) = {vc};  Transfinite Volume{vc} = {77, 61, 109, 125, 94, 62, 110, 142};
vc++;  Surface Loop(vc) = {82, 114, 190, 214, 238, 215};  Volume(vc) = {vc};  Transfinite Volume{vc} = {95, 63, 62, 94, 143, 111, 110, 142};
vc++;  Surface Loop(vc) = {83, 115, 191, 215, 239, 216};  Volume(vc) = {vc};  Transfinite Volume{vc} = {96, 64, 63, 95, 144, 112, 111, 143};
vc++;  Surface Loop(vc) = {84, 116, 192, 216, 240, 217};  Volume(vc) = {vc};  Transfinite Volume{vc} = {129, 113, 65, 81, 144, 112, 64, 96};
vc++;  Surface Loop(vc) = {85, 101, 193, 217, 241, 218};  Volume(vc) = {vc};  Transfinite Volume{vc} = {130, 114, 66, 82, 129, 113, 65, 81};
vc++;  Surface Loop(vc) = {86, 102, 194, 218, 242, 219};  Volume(vc) = {vc};  Transfinite Volume{vc} = {130, 114, 115, 131, 82, 66, 67, 83};
vc++;  Surface Loop(vc) = {87, 103, 195, 219, 243, 220};  Volume(vc) = {vc};  Transfinite Volume{vc} = {131, 115, 116, 132, 83, 67, 68, 84};
vc++;  Surface Loop(vc) = {88, 104, 196, 220, 244, 213};  Volume(vc) = {vc};  Transfinite Volume{vc} = {132, 116, 109, 125, 84, 68, 61, 77};
//+
// Volumes in the boundary layer at zmax
vc++;  Surface Loop(vc) = {89, 105, 197, 221, 245, 222};  Volume(vc) = {vc};  Transfinite Volume{vc} = {86, 77, 125, 134, 70, 69, 117, 118};
vc++;  Surface Loop(vc) = {90, 106, 198, 222, 246, 223};  Volume(vc) = {vc};  Transfinite Volume{vc} = {135, 87, 86, 134, 119, 71, 70, 118};
vc++;  Surface Loop(vc) = {91, 107, 199, 223, 247, 224};  Volume(vc) = {vc};  Transfinite Volume{vc} = {136, 88, 87, 135, 120, 72, 71, 119};
vc++;  Surface Loop(vc) = {92, 108, 200, 224, 248, 225};  Volume(vc) = {vc};  Transfinite Volume{vc} = {136, 129, 81, 88, 120, 121, 73, 72};
vc++;  Surface Loop(vc) = {93, 116, 201, 225, 249, 226};  Volume(vc) = {vc};  Transfinite Volume{vc} = {129, 144, 96, 81, 121, 122, 74, 73};
vc++;  Surface Loop(vc) = {94, 115, 202, 226, 250, 227};  Volume(vc) = {vc};  Transfinite Volume{vc} = {96, 144, 143, 95, 74, 122, 123, 75};
vc++;  Surface Loop(vc) = {95, 114, 203, 227, 251, 228};  Volume(vc) = {vc};  Transfinite Volume{vc} = {95, 143, 142, 94, 75, 123, 124, 76};
vc++;  Surface Loop(vc) = {96, 113, 204, 228, 252, 221};  Volume(vc) = {vc};  Transfinite Volume{vc} = {77, 94, 142, 125, 69, 76, 124, 117};

Printf("Boundary layer finished with volume %g", vc);

Physical Volume("FLUID") = { 1 : vc };

//-------------------------
//
// Set boundary conditions
//
//-------------------------

Physical Surface("PIPE_WALL") = { wall_first : wall_last };


Physical Surface("X_MIN") = {Surface In BoundingBox{-LI-TINY, -HUGE, -HUGE,
                                                    -LI+TINY, +HUGE, +HUGE}};

Physical Surface("X_MAX") = {Surface In BoundingBox{+LO-TINY, -HUGE, -HUGE,
                                                    +LO+TINY, +HUGE, +HUGE}};

Physical Surface("Z_MAX") = {Surface In BoundingBox{-HUGE, -HUGE, LI-TINY,
                                                    +HUGE, +HUGE, LI+TINY}};

//--------------------------------------------------------------------
//
// Very important: call Cohrence just before the end - for without it
// the parallel version will not work properly because the buffers,
// being based on nodes, will not work if some nodes are duplicated
//
//--------------------------------------------------------------------
Coherence;

