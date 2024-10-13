//==============================================================================
//
//
// Handling command line options
//
//
//------------------------------------------------------------------------------
If(!Exists(FLUID) || !Exists(SOLID))
  Printf("Variables FLUID and/or SOLID are not defined from command line.");
  Printf("The proper invokation of Gmsh for this scrip is, for example:");
  Printf("");
  Printf("gmsh -3 t_junction.geo -setnumber FLUID 1 -setnumber SOLID 0 -o fluid.msh");
  Printf("");
  Printf("to mesh only fluid domain, or:");
  Printf("");
  Printf("gmsh -3 t_junction.geo -setnumber FLUID 0 -setnumber SOLID 1 -o solid.msh");
  Printf("");
  Printf("to mesh only solid domain.");
  Abort;
EndIf

//==============================================================================
//
//
// Definition of constants
//
//
//------------------------------------------------------------------------------
LI =  4.0;
LO = 12.0;
R(1) = 0.6;   // core radius
R(2) = 0.9;  // boundary layer radius
R(3) = 1.0;  // outer radius
R(4) = 1.5;  // outer radius

Printf("FLUID: %g", FLUID);

// Number of nodes (resolutions for lines)
N_CORE_ARC    = 11;      // number of nodes in the core
N_STREAM_I    = 41;      // number of nodes in the stream-wise direction
N_STREAM_O    = 61;      // number of nodes in the stream-wise direction
N_MIDDLE      =  9;
N_BLAYER      =  7;
N_SLAYER      =  5;

// Progression for lines
P_CORE_RAD    = 1.1;
P_STREAM_I    = 1.05;
P_STREAM_O    = 1.05;
P_MIDDLE      = 1.15;
P_BLAYER      = 1.3;
P_SLAYER      = 1.3;

// Connector types (will be used to define resolution and progression)
TYPE_ARC      = 111;  // all arches and ellipses
TYPE_BLAYER   = 222;  // lines in the buffer layer
TYPE_CORE_RAD = 333;  // lines defining the core
TYPE_CORE_TAN = 444;  // lines defining the core
TYPE_MIDDLE   = 555;  // lines connecting core and boundary layer
TYPE_STREAM_I = 666;  // lines in the streamwise directions in inlet
TYPE_STREAM_O = 777;  // lines in the streamwise directions in outlet
TYPE_SLAYER   = 888;  // lines in the solid layer

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
    pnt_1 = pnt_1st + (n - 1);
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

// Points at the wet cylinder wall
nf = 1;  nl = 8;
x = -LI;  r = R(3);  Call NPointsAroundX;  // point 101 on
x = +LO;  r = R(3);  Call NPointsAroundX;  // point 109 on
z = +LI;  r = R(3);  Call NPointsAroundZ;  // point 117 on
x = 0.0;  r = R(3);  Call NPointsAroundX;  // point 125 on
r = R(3);  Call NPointsMinusX;             // point 133 on
r = R(3);  Call NPointsPlusX;              // point 141 on

// Points at the dry cylinder wall
nf = 1;  nl = 8;
x = -LI;  r = R(4);  Call NPointsAroundX;  // point 149 on
x = +LO;  r = R(4);  Call NPointsAroundX;  // point 157 on
z = +LI;  r = R(4);  Call NPointsAroundZ;  // point 165 on
x = 0.0;  r = R(4);  Call NPointsAroundX;  // point 173 on
r = R(4);  Call NPointsMinusX;             // point 181 on
r = R(4);  Call NPointsPlusX;              // point 189 to 196

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
// 7 - boundary layer radial connectors
// 8 - connectors for the outer walls (arches and ellipses)
//
// Planes are ordered in the following way
// 1 - at xmin
// 2 - at xmax
// 3 - at zmax
// 4 - at x = 0
// 5 - at x = 0 rotated by 45 degrees backslash
// 6 - at x = 0 rotated by 45 degrees slash

//------------------------------------------------
// 1 - center in the core in the radial direction
//------------------------------------------------

// Lines in radial direction of the cores
type = TYPE_CORE_RAD;

nf = 1;  nl = 8;
pnt_cent = 1;  pnt_1st =  5;  Call NLinesFromCenter;  // lines  1 -  4
pnt_cent = 2;  pnt_1st = 13;  Call NLinesFromCenter;  // lines  5 -  8
pnt_cent = 3;  pnt_1st = 21;  Call NLinesFromCenter;  // lines  9 - 12
pnt_cent = 4;  pnt_1st = 29;  Call NLinesFromCenter;  // lines 13 - 16
pnt_cent = 4;  pnt_1st = 37;  Call NLinesFromCenter;  // lines 17 - 20
pnt_cent = 4;  pnt_1st = 45;  Call NLinesFromCenter;  // lines 21 - 24

//--------------------------------------------------------------------
// 2 - connectors wrapping the core (on the border with middle layer)
//--------------------------------------------------------------------
type = TYPE_CORE_TAN;

nf = 1;  nl = 8;
pnt_1st =  5;  Call NLinesAround;  // lines 25 - 32
pnt_1st = 13;  Call NLinesAround;  // lines 33 - 40
pnt_1st = 21;  Call NLinesAround;  // lines 41 - 48
pnt_1st = 29;  Call NLinesAround;  // lines 49 - 56
pnt_1st = 37;  Call NLinesAround;  // lines 57 - 64
pnt_1st = 45;  Call NLinesAround;  // lines 65 - 72

//------------------------------
// 3 - middle radial connectors
//------------------------------
type = TYPE_MIDDLE;

nf = 1; nl = 8;
pnt_1st = 53;  pnt_2nd =  5;  Call NLinesBetweenLayers;
pnt_1st = 61;  pnt_2nd = 13;  Call NLinesBetweenLayers;
pnt_1st = 69;  pnt_2nd = 21;  Call NLinesBetweenLayers;
pnt_1st = 77;  pnt_2nd = 29;  Call NLinesBetweenLayers;
pnt_1st = 85;  pnt_2nd = 37;  Call NLinesBetweenLayers;
pnt_1st = 93;  pnt_2nd = 45;  Call NLinesBetweenLayers;

//------------------------------------------------------------------
// 4 - connectors wrapping the bundary layers (arches and ellipses)
//------------------------------------------------------------------
type = TYPE_ARC;

nf = 1;  nl = 8;
pnt_1st = 53;  pnt_cent = 1;  Call NArcsAround;
pnt_1st = 61;  pnt_cent = 2;  Call NArcsAround;
pnt_1st = 69;  pnt_cent = 3;  Call NArcsAround;
pnt_1st = 77;  pnt_cent = 4;  Call NArcsAround;
pnt_1st = 85;  pnt_cent = 4;  pnt_axis = 87;  Call NEllipseArcsAround;
pnt_1st = 93;  pnt_cent = 4;  pnt_axis = 95;  Call NEllipseArcsAround;

//--------------------------------------
// 5 - boundary layer radial connectors
//--------------------------------------
type = TYPE_BLAYER;

nf = 1; nl = 8;
pnt_1st = 101;  pnt_2nd = 53;  Call NLinesBetweenLayers;
pnt_1st = 109;  pnt_2nd = 61;  Call NLinesBetweenLayers;
pnt_1st = 117;  pnt_2nd = 69;  Call NLinesBetweenLayers;
pnt_1st = 125;  pnt_2nd = 77;  Call NLinesBetweenLayers;
pnt_1st = 133;  pnt_2nd = 85;  Call NLinesBetweenLayers;
pnt_1st = 141;  pnt_2nd = 93;  Call NLinesBetweenLayers;

//------------------------------------------------------------------
// 6 - connectors wrapping the bundary layers (arches and ellipses)
//------------------------------------------------------------------
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

//-----------------------------------
// 7 - solid layer radial connectors
//-----------------------------------
type = TYPE_SLAYER;
nf = 1; nl = 8;
pnt_1st = 101;  pnt_2nd = 149;  Call NLinesBetweenLayers;
pnt_1st = 109;  pnt_2nd = 157;  Call NLinesBetweenLayers;
pnt_1st = 117;  pnt_2nd = 165;  Call NLinesBetweenLayers;  // at zmax
pnt_1st = 125;  pnt_2nd = 173;  Call NLinesBetweenLayers;
pnt_1st = 133;  pnt_2nd = 181;  Call NLinesBetweenLayers;
pnt_1st = 141;  pnt_2nd = 189;  Call NLinesBetweenLayers;

//----------------------------------------------------------------
// 8 - connectors wrapping the solid layers (arches and ellipses)
//----------------------------------------------------------------
type = TYPE_ARC;
nf = 1;  nl = 8;
pnt_1st = 149;  pnt_cent = 1;  Call NArcsAround;
pnt_1st = 157;  pnt_cent = 2;  Call NArcsAround;
pnt_1st = 165;  pnt_cent = 3;  Call NArcsAround;
pnt_1st = 173;  pnt_cent = 4;  Call NArcsAround;
pnt_1st = 181;  pnt_cent = 4;  pnt_axis = 135;  Call NEllipseArcsAround;
pnt_1st = 189;  pnt_cent = 4;  pnt_axis = 143;  Call NEllipseArcsAround;

// Towards xmin: from axis, through core, boundary layers and pipe walls
type = TYPE_STREAM_I;
nf = 1; nl = 4;  pnt_1st =181;  pnt_2nd =149;  Call NLinesBetweenLayers;
nf = 5; nl = 8;  pnt_1st =173;  pnt_2nd =149;  Call NLinesBetweenLayers;

// Towards xmax: from axis, through core, boundary layers and pipe walls
type = TYPE_STREAM_O;
nf = 1; nl = 4;  pnt_1st =189;  pnt_2nd =157;  Call NLinesBetweenLayers;
nf = 5; nl = 8;  pnt_1st =173;  pnt_2nd =157;  Call NLinesBetweenLayers;

// Towards zmax
type = TYPE_STREAM_I;
nf = 1; nl = 4;  pnt_1st =181;  pnt_2nd =165;  Call NLinesBetweenLayers;
nf = 5; nl = 8;  pnt_1st =197;  pnt_2nd =165;  Call NLinesBetweenLayersNeg;

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
  If(lin_type(l) == TYPE_SLAYER)
    Transfinite Curve {l} = N_SLAYER    Using Progression P_SLAYER;
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

// Surfaces in the middle at xmin leg ("orthogonal" to walls)
sc++;  Curve Loop(sc) = {73, -266, -105, 274};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {74, -267, -106, 275};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {75, -268, -107, 276};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {76, -269, -108, 277};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {77, -270, -101, 278};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {78, -271, -102, 279};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {79, -272, -103, 280};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {80, -273, -104, 281};  Plane Surface(sc) = {sc};

// Surfaces in the middle at xmin leg ("orthogonal" to walls)
sc++;  Curve Loop(sc) = {81, -291, -113, 299};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {82, -292, -114, 300};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {83, -293, -115, 301};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {84, -294, -116, 302};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {85, -295, -109, 303};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {86, -296, -102, 304};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {87, -297, -103, 305};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {88, -298, -104, 306};  Plane Surface(sc) = {sc};

// Surfaces in the middle at zmax leg ("orthogonal" to walls)
sc++;  Curve Loop(sc) = {89, -316, -105, 324};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {90, -317, -106, 325};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {91, -318, -107, 326};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {92, -319, -108, 327};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {93, -320, -101, 328};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {94, -321, -116, 329};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {95, -322, -115, 330};  Plane Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {96, -323, -114, 331};  Plane Surface(sc) = {sc};

// Wrapping the boundary layer at xmin ("parallel" to walls)
sc++;  Curve Loop(sc) = {121, -275, -153, 274};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {122, -276, -154, 275};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {123, -277, -155, 276};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {124, -278, -156, 277};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {125, -279, -149, 278};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {126, -280, -150, 279};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {127, -281, -151, 280};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {128, -274, -152, 281};  Surface(sc) = {sc};

// Wrapping the boundary layer at xmax ("parallel" to walls)
sc++;  Curve Loop(sc) = {129, -300, -161, 299};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {130, -301, -162, 300};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {131, -302, -163, 301};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {132, -303, -164, 302};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {133, -304, -149, 303};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {134, -305, -150, 304};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {135, -306, -151, 305};  Surface(sc) = {sc};
sc++;  Curve Loop(sc) = {136, -299, -152, 306};  Surface(sc) = {sc};

// Wrapping the boundary layer at zmax ("parallel" to walls)
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

Printf("Started defining wetted walls with surface: %g", sc+1);
wet_wall_first = sc + 1;

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

Printf("Finished defining wetted walls with surface: %g", sc);
wet_wall_last = sc;

//-----------------------------------------------------------------
// Surfaces orthogonal to axes (increment factors are 124 and 172)
//-----------------------------------------------------------------

// at xmin
n = 217;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124-8, -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;

// at xmax
n = 225;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124-8, -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;

// at zmax
n = 233;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124-8, -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;

// at core
n = 241;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124-8, -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
n = 249;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124-8, -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
n = 257;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124,   -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+124-8, -n-172+1, -n-124+1};  Plane Surface(sc) = {sc};  n++;

//---------------------
// Streamwise surfaces
//---------------------

// Orthogonal to walls at xmin
n = 282;
sc++;  Curve Loop(sc) = {n, n+58, -n-154, -n-90+8};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+58, -n-154, -n-90};    Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+58, -n-154, -n-90};    Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+58, -n-154, -n-90};    Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+58, -n-154, -n-90};    Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+58, -n-154, -n-90+8};  Plane Surface(sc) = {sc};  n++;  // 287
sc++;  Curve Loop(sc) = {n, n+58, -n-154, -n-90+8};  Plane Surface(sc) = {sc};  n++;  // 288
sc++;  Curve Loop(sc) = {n, n+58, -n-154, -n-90+8};  Plane Surface(sc) = {sc};  n++;

n = 307;
sc++;  Curve Loop(sc) = {n, n+41, -n-137, -n-57};     Plane Surface(sc) = {sc};  n++;  // 307
sc++;  Curve Loop(sc) = {n, n+41, -n-137, -n-57-16};  Plane Surface(sc) = {sc};  n++;  // 308
sc++;  Curve Loop(sc) = {n, n+41, -n-137, -n-57-16};  Plane Surface(sc) = {sc};  n++;  // 309
sc++;  Curve Loop(sc) = {n, n+41, -n-137, -n-57-16};  Plane Surface(sc) = {sc};  n++;  // 310
sc++;  Curve Loop(sc) = {n, n+41, -n-137, -n-57-8};   Plane Surface(sc) = {sc};  n++;  // 311
sc++;  Curve Loop(sc) = {n, n+41, -n-137, -n-57};     Plane Surface(sc) = {sc};  n++;  // 312
sc++;  Curve Loop(sc) = {n, n+41, -n-137, -n-57};     Plane Surface(sc) = {sc};  n++;  // 313
sc++;  Curve Loop(sc) = {n, n+41, -n-137, -n-57};     Plane Surface(sc) = {sc};  n++;  // 314

// Orthogonal to walls at z min
n = 332;
sc++;  Curve Loop(sc) = {n, n+24, -n-120, -n-40+8};  Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+24, -n-120, -n-40};    Plane Surface(sc) = {sc};  n++;  // 333
sc++;  Curve Loop(sc) = {n, n+24, -n-120, -n-40};    Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+24, -n-120, -n-40};    Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+24, -n-120, -n-40+8};  Plane Surface(sc) = {sc};  n++;  // 336
sc++;  Curve Loop(sc) = {n, n+24, -n-120, -383};     Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+24, -n-120, -382};     Plane Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n+24, -n-120, -381};     Plane Surface(sc) = {sc};  n++;  // 339

Printf("Started defining dry walls with surface: %g", sc+1);
dry_wall_first = sc + 1;

// Outer solid walls wrappings at xmin
n = 436;
sc++;  Curve Loop(sc) = {n, n-48, -n-1,   -n+16};  Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n-48, -n-1,   -n+16};  Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n-48, -n-1,   -n+16};  Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n-48, -n-1,   -n+16};  Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n-48, -n-1,   -n+24};  Surface(sc) = {sc};  n++;  // 440
sc++;  Curve Loop(sc) = {n, n-48, -n-1,   -n+24};  Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n-48, -n-1,   -n+24};  Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n-48, -n-1+8, -n+24};  Surface(sc) = {sc};  n++;

// Outer solid walls wrappings at xmax
n = 444;
sc++;  Curve Loop(sc) = {n, n-48, -n-1,   -n+16};  Surface(sc) = {sc};  n++;  // 444
sc++;  Curve Loop(sc) = {n, n-48, -n-1,   -n+16};  Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n-48, -n-1,   -n+16};  Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n-48, -n-1,   -n+16};  Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n-48, -n-1,   -n+32};  Surface(sc) = {sc};  n++;  // 448
sc++;  Curve Loop(sc) = {n, n-48, -n-1,   -n+32};  Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n-48, -n-1,   -n+32};  Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n-48, -n-1+8, -n+32};  Surface(sc) = {sc};  n++;

// Outer solid walls wrappings at xmax
n = 452;
sc++;  Curve Loop(sc) = {n, n-48, -n-1,   -n+32};  Surface(sc) = {sc};  n++;  // 452
sc++;  Curve Loop(sc) = {n, n-48, -n-1,   -n+32};  Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n-48, -n-1,   -n+32};  Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n-48, -n-1,   -n+32};  Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n-48, -n-1,    431};   Surface(sc) = {sc};  n++;  // 456
sc++;  Curve Loop(sc) = {n, n-48, -n-1,    430};   Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n-48, -n-1,    429};   Surface(sc) = {sc};  n++;
sc++;  Curve Loop(sc) = {n, n-48, -n-1+8,  428};   Surface(sc) = {sc};  n++;

Printf("Finished defining dry walls with surface: %g", sc);
dry_wall_last = sc;

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
  Surface{ 13}; Surface{ 14};  // core surfaces in the central orthogonal set
  Surface{ 19}; Surface{ 20};  // core surfaces in the central anticlockw set
  Surface{ 23}; Surface{ 24};  // core surfaces in the central clockwise set
  Surface{ 49}; Surface{ 50}; Surface{ 51}; Surface{ 52};  // middle orthogonal set
  Surface{ 61}; Surface{ 62}; Surface{ 63}; Surface{ 64};  // middle backslash
  Surface{ 72}; Surface{ 71}; Surface{ 70}; Surface{ 69};  // middle slash
  Surface{ 97}; Surface{ 98}; Surface{ 99}; Surface{100};  // blayer orthogonal
  Surface{109}; Surface{110}; Surface{111}; Surface{112};  // blayer backslash
  Surface{120}; Surface{119}; Surface{118}; Surface{117};  // blayer slash
  Surface{277}; Surface{278}; Surface{279}; Surface{280};  // slayer orthogonal
  Surface{289}; Surface{290}; Surface{291}; Surface{292};  // slayer backslash
  Surface{297}; Surface{298}; Surface{299}; Surface{300};  // slayer slash
}

Coherence;

//---------
//
// Volumes
//
//---------

wet_volume_first = vc + 1;

// Volumes in the core at xmin
vc++;  Surface Loop(vc) = { 1, 17, 121, 133, 134, 122};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = { 2, 18, 122, 135, 136, 123};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = { 3, 15, 123, 137, 138, 124};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = { 4, 16, 124, 139, 140, 121};  Volume(vc) = {vc};
Transfinite Volume{vc};

// Volumes in the core at xmax
vc++;  Surface Loop(vc) = { 5, 21, 125, 141, 142, 126};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = { 6, 22, 126, 143, 144, 127};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = { 7, 15, 127, 145, 146, 128};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = { 8, 16, 128, 147, 148, 125};  Volume(vc) = {vc};
Transfinite Volume{vc};

// Volumes in the core at zmax
vc++;  Surface Loop(vc) = { 9, 17, 129, 149, 150, 130};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {10, 18, 130, 151, 152, 131};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {11, 22, 131, 153, 154, 132};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {12, 21, 132, 155, 156, 129};  Volume(vc) = {vc};
Transfinite Volume{vc};

Printf("Core finished with volume %g", vc);

// Volumes in the middle at xmin
vc++;  Surface Loop(vc) = {25, 57, 133, 157, 181, 158};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {26, 58, 134, 158, 182, 159};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {27, 59, 135, 159, 183, 160};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {28, 60, 136, 160, 184, 161};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {29, 53, 137, 161, 185, 162};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {30, 54, 138, 162, 186, 163};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {31, 55, 139, 163, 187, 164};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {32, 56, 140, 164, 188, 157};  Volume(vc) = {vc};
Transfinite Volume{vc};

// Volumes in the middle at xmax
vc++;  Surface Loop(vc) = {33, 65, 165, 189, 166, 141};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {34, 66, 166, 190, 167, 142};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {35, 67, 167, 191, 168, 143};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {36, 68, 168, 192, 169, 144};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {37, 53, 169, 193, 170, 145};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {38, 54, 170, 194, 171, 146};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {39, 55, 171, 195, 172, 147};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {40, 56, 172, 196, 165, 148};  Volume(vc) = {vc};
Transfinite Volume{vc};

// Volumes in the middle at zmax
vc++;  Surface Loop(vc) = {41, 57, 149, 173, 197, 174};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {42, 58, 150, 174, 198, 175};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {43, 59, 151, 175, 199, 176};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {44, 60, 152, 176, 200, 177};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {45, 68, 153, 177, 201, 178};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {46, 67, 154, 178, 202, 179};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {47, 66, 155, 179, 203, 180};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {48, 65, 156, 180, 204, 173};  Volume(vc) = {vc};
Transfinite Volume{vc};

Printf("Middle finished with volume %g", vc);

// Volumes in the boundary layer at xmin (like middle
// at xmin +52 except second column in the lower part?)
vc++;  Surface Loop(vc) = {73, 105, 181, 205, 229, 206};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {74, 106, 182, 206, 230, 207};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {75, 107, 183, 207, 231, 208};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {76, 108, 184, 208, 232, 209};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {77, 101, 185, 209, 233, 210};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {78, 102, 186, 210, 234, 211};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {79, 103, 187, 211, 235, 212};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {80, 104, 188, 212, 236, 205};  Volume(vc) = {vc};
Transfinite Volume{vc};

// Volumes in the boundary layer at xmax
vc++;  Surface Loop(vc) = {81, 113, 189, 213, 237, 214};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {82, 114, 190, 214, 238, 215};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {83, 115, 191, 215, 239, 216};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {84, 116, 192, 216, 240, 217};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {85, 101, 193, 217, 241, 218};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {86, 102, 194, 218, 242, 219};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {87, 103, 195, 219, 243, 220};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {88, 104, 196, 220, 244, 213};  Volume(vc) = {vc};
Transfinite Volume{vc};

// Volumes in the boundary layer at zmax
vc++;  Surface Loop(vc) = {89, 105, 197, 221, 245, 222};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {90, 106, 198, 222, 246, 223};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {91, 107, 199, 223, 247, 224};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {92, 108, 200, 224, 248, 225};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {93, 116, 201, 225, 249, 226};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {94, 115, 202, 226, 250, 227};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {95, 114, 203, 227, 251, 228};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {96, 113, 204, 228, 252, 221};  Volume(vc) = {vc};
Transfinite Volume{vc};

Printf("Boundary layer finished with volume %g", vc);

wet_volume_last = vc;

//--------------------------------------------------------------------
//
// Very important: call Cohrence just before the end - for without it
// the parallel version will not work properly because the buffers,
// being based on nodes, will not work if some nodes are duplicated
//
//--------------------------------------------------------------------
Coherence;

// Hide "*";
// Show {
//  Point{101}; Point{102}; Point{103}; Point{104}; Point{105}; Point{106}; Point{107}; Point{108}; Point{109}; Point{110}; Point{111}; Point{112}; Point{113}; Point{114}; Point{115}; Point{116}; Point{117}; Point{118}; Point{119}; Point{120}; Point{121}; Point{122}; Point{123}; Point{124}; Point{125}; Point{129}; Point{130}; Point{131}; Point{132}; Point{134}; Point{135}; Point{136}; Point{142}; Point{143}; Point{144}; Point{149}; Point{150}; Point{151}; Point{152}; Point{153}; Point{154}; Point{155}; Point{156}; Point{157}; Point{158}; Point{159}; Point{160}; Point{161}; Point{162}; Point{163}; Point{164}; Point{165}; Point{166}; Point{167}; Point{168}; Point{169}; Point{170}; Point{171}; Point{172}; Point{173}; Point{177}; Point{178}; Point{179}; Point{180}; Point{182}; Point{183}; Point{184}; Point{190}; Point{191}; Point{192}; Curve{217}; Curve{218}; Curve{219}; Curve{220}; Curve{221}; Curve{222}; Curve{223}; Curve{224}; Curve{225}; Curve{226}; Curve{227}; Curve{228}; Curve{229}; Curve{230}; Curve{231}; Curve{232}; Curve{233}; Curve{234}; Curve{235}; Curve{236}; Curve{237}; Curve{238}; Curve{239}; Curve{240}; Curve{245}; Curve{246}; Curve{247}; Curve{248}; Curve{249}; Curve{250}; Curve{251}; Curve{252}; Curve{257}; Curve{258}; Curve{259}; Curve{260}; Curve{282}; Curve{283}; Curve{284}; Curve{285}; Curve{286}; Curve{287}; Curve{288}; Curve{289}; Curve{307}; Curve{308}; Curve{309}; Curve{310}; Curve{311}; Curve{312}; Curve{313}; Curve{314}; Curve{332}; Curve{333}; Curve{334}; Curve{335}; Curve{336}; Curve{337}; Curve{338}; Curve{339}; Curve{340}; Curve{341}; Curve{342}; Curve{343}; Curve{344}; Curve{345}; Curve{346}; Curve{347}; Curve{348}; Curve{349}; Curve{350}; Curve{351}; Curve{352}; Curve{353}; Curve{354}; Curve{355}; Curve{356}; Curve{357}; Curve{358}; Curve{359}; Curve{360}; Curve{361}; Curve{362}; Curve{363}; Curve{364}; Curve{368}; Curve{369}; Curve{370}; Curve{371}; Curve{373}; Curve{374}; Curve{375}; Curve{381}; Curve{382}; Curve{383}; Curve{388}; Curve{389}; Curve{390}; Curve{391}; Curve{392}; Curve{393}; Curve{394}; Curve{395}; Curve{396}; Curve{397}; Curve{398}; Curve{399}; Curve{400}; Curve{401}; Curve{402}; Curve{403}; Curve{404}; Curve{405}; Curve{406}; Curve{407}; Curve{408}; Curve{409}; Curve{410}; Curve{411}; Curve{416}; Curve{417}; Curve{418}; Curve{419}; Curve{420}; Curve{421}; Curve{422}; Curve{423}; Curve{428}; Curve{429}; Curve{430}; Curve{431}; Curve{436}; Curve{437}; Curve{438}; Curve{439}; Curve{440}; Curve{441}; Curve{442}; Curve{443}; Curve{444}; Curve{445}; Curve{446}; Curve{447}; Curve{448}; Curve{449}; Curve{450}; Curve{451}; Curve{452}; Curve{453}; Curve{454}; Curve{455}; Curve{456}; Curve{457}; Curve{458}; Curve{459}; Surface{229}; Surface{230}; Surface{231}; Surface{232}; Surface{233}; Surface{234}; Surface{235}; Surface{236}; Surface{237}; Surface{238}; Surface{239}; Surface{240}; Surface{241}; Surface{242}; Surface{243}; Surface{244}; Surface{245}; Surface{246}; Surface{247}; Surface{248}; Surface{249}; Surface{250}; Surface{251}; Surface{252}; Surface{253}; Surface{254}; Surface{255}; Surface{256}; Surface{257}; Surface{258}; Surface{259}; Surface{260}; Surface{261}; Surface{262}; Surface{263}; Surface{264}; Surface{265}; Surface{266}; Surface{267}; Surface{268}; Surface{269}; Surface{270}; Surface{271}; Surface{272}; Surface{273}; Surface{274}; Surface{275}; Surface{276}; Surface{281}; Surface{282}; Surface{283}; Surface{284}; Surface{285}; Surface{286}; Surface{287}; Surface{288}; Surface{293}; Surface{294}; Surface{295}; Surface{296}; Surface{301}; Surface{302}; Surface{303}; Surface{304}; Surface{305}; Surface{306}; Surface{307}; Surface{308}; Surface{309}; Surface{310}; Surface{311}; Surface{312}; Surface{313}; Surface{314}; Surface{315}; Surface{316}; Surface{317}; Surface{318}; Surface{319}; Surface{320}; Surface{321}; Surface{322}; Surface{323}; Surface{324}; Surface{325}; Surface{326}; Surface{327}; Surface{328}; Surface{329}; Surface{330}; Surface{331}; Surface{332}; Surface{333}; Surface{334}; Surface{335}; Surface{336}; Surface{337}; Surface{338}; Surface{339}; Surface{340}; Surface{341}; Surface{342}; Surface{343}; Surface{344}; Surface{345}; Surface{346}; Surface{347}; Surface{348}; 
// }

dry_volume_first = vc + 1;

// Volumes in the solid at xmin
vc++;  Surface Loop(vc) = {253, 325, 285, 301, 229, 302};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {254, 326, 230, 286, 303, 302};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {255, 327, 231, 287, 304, 303};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {256, 328, 288, 304, 232, 305};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {257, 329, 305, 281, 306, 233};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {258, 330, 282, 307, 234, 306};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {259, 331, 283, 235, 308, 307};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {260, 332, 284, 236, 301, 308};  Volume(vc) = {vc};
Transfinite Volume{vc};

// Volumes in the solid at xmax
vc++;  Surface Loop(vc) = {261, 333, 293, 309, 237, 310};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {262, 334, 294, 238, 310, 311};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {263, 335, 295, 239, 312, 311};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {264, 336, 296, 312, 313, 240};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {265, 337, 281, 241, 313, 314};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {266, 338, 282, 242, 315, 314};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {267, 339, 283, 243, 316, 315};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {268, 340, 284, 316, 309, 244};  Volume(vc) = {vc};
Transfinite Volume{vc};

// Volumes in the solid at zmax
vc++;  Surface Loop(vc) = {341, 269, 318, 317, 285, 245};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {318, 270, 342, 246, 319, 286};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {271, 343, 287, 247, 320, 319};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {272, 344, 288, 321, 248, 320};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {273, 345, 321, 322, 296, 249};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {274, 346, 250, 322, 323, 295};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {275, 347, 324, 323, 251, 294};  Volume(vc) = {vc};
Transfinite Volume{vc};
vc++;  Surface Loop(vc) = {276, 348, 317, 252, 293, 324};  Volume(vc) = {vc};
Transfinite Volume{vc};

dry_volume_last = vc;

Printf("Solid finished with volume %g", vc);

Coherence;

If(FLUID == 1)
  Physical Volume("FLUID") = { wet_volume_first : wet_volume_last };
Else
  For v In { wet_volume_first : wet_volume_last }
    Recursive Delete {
      Volume{v};
    }
  EndFor
EndIf

If(SOLID == 1)
  Physical Volume("SOLID") = { dry_volume_first : dry_volume_last };
Else
  For v In { dry_volume_first : dry_volume_last }
    Recursive Delete {
      Volume{v};
    }
  EndFor
EndIf

//-------------------------
//
// Set boundary conditions
//
//-------------------------

If(FLUID == 1)
  Physical Surface("WET_WALL") = { wet_wall_first : wet_wall_last };
EndIf

If(SOLID == 1)
  Physical Surface("DRY_WALL") = { dry_wall_first : dry_wall_last };
EndIf

Physical Surface("X_MIN") = {Surface In BoundingBox{-LI-TINY, -HUGE, -HUGE,
                                                    -LI+TINY, +HUGE, +HUGE}};

Physical Surface("X_MAX") = {Surface In BoundingBox{+LO-TINY, -HUGE, -HUGE,
                                                    +LO+TINY, +HUGE, +HUGE}};

Physical Surface("Z_MAX") = {Surface In BoundingBox{-HUGE, -HUGE, LI-TINY,
                                                    +HUGE, +HUGE, LI+TINY}};

