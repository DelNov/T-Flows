//==============================================================================
//
// Definition of constants
//
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//  Geometry (m)
//------------------------------------------------------------------------------

Rad = 0.0025; // radius of the pipe
L   = 0.1;   // length of the pipe

//------------------------------------------------------------------------------
//  Meshing 
//------------------------------------------------------------------------------

// layers
R(1) = Rad/2;   // core radius
R(2) = Rad*8/9; // middle layer radius
R(3) = Rad;     // boundary layer radius

//  Options
side           =  8;  // number of sides in the core (must be a multiple of 4 to avoid issues with ellipse generation)
coreRad_skip   =  2;  // alternance in radial connections for the core
middleRad_skip =  1;  // alternance in radial connections for the middle
blayerRad_skip =  1;  // alternance in radial connections for the bounday layer
wallRad_skip   =  1;  // alternance in radial connections for the wall

// Number of nodes (resolutions for lines)
N_CORE                 =  10;      // number of nodes in the radial direction inside the core 
N_MIDDLE               =  12;      // number of nodes in the radial direction inside the middle layer 
N_BLAYER               =   7;      // number of nodes in the radial direction inside the boundary layer 
N_STREAM               =  110;      // number of nodes in the inlet stream-wise direction 
N_ARC                  =  10;      // number of nodes in the arches and ellipses 

// Progression for lines
P_CORE                 = 1.00;
P_MIDDLE               = 1.10;
P_BLAYER               = 1.10;
P_STREAM               = 1.00;        
P_ARC                  = 1.0;

// Connector types (will be used to define resolution and progression)
TYPE_CORE              = 1;  // lines defining the core
TYPE_MIDDLE            = 2;  // lines connecting core and boundary layer
TYPE_BLAYER            = 3;  // lines in the buffer layer
TYPE_STREAM            = 4;  // streamwise lines  
TYPE_ARC               = 5;  // all arches and ellipses

// Tolerance for merging nodes
NODE_MERGE_TOL = 1.0e-6;

//==============================================================================
//
// Definition of macros
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

// place eight points in a circle in plane x
Macro NPointsAroundX  
  For n In{ from : to }
    a = (n-1) * (2.0 * Pi / (to - from + 1));
    pc++;  Point(pc) = {x, r * Cos(a), r * Sin(a)};
  EndFor
Return

//------------------------------------------------------------------------------

// connection from central points to core points
Macro LinesCoreRad  
  For c In { 0 : repetition - 1 : step}
    pnt_cent_c = pnt_cent + c;
    pnt_1st_c = pnt_1st + side * c;
    For n In{ from : to : skip }
      lc++;
      pnt_1 = pnt_1st_c + (n - 1);
      Printf("Creating radial line %g from points: %g, %g", lc, pnt_1, pnt_cent_c);
      Line(lc) = {pnt_copy(pnt_1), pnt_copy(pnt_cent_c)};
      lin_type[lc] = type;
    EndFor
  EndFor
Return

// connection of consecutive points in a circle
Macro LinesCoreTan 
  For c In { 0 : repetition - 1 : step }
    pnt_1st_c = pnt_1st + side * c;
    For n In{ from : to }
      lc++;
      pnt_1 = pnt_1st_c + (n - 1);
      pnt_2 = pnt_1 + 1;
      If(n == to)
        pnt_2 = pnt_1st_c;
      EndIf
      Printf("Creating tangential line %g around from points: %g, %g", lc, pnt_1, pnt_2);
      Line(lc) = {pnt_copy(pnt_1), pnt_copy(pnt_2)};
      lin_type[lc] = type;
    EndFor
  EndFor
Return

// connection between points in 2 concentric circles (radially)
Macro LinesBtwLayers  
  For c In { 0 : repetition - 1 : step }
    pnt_1st_c = pnt_1st + side * c;
    pnt_2nd_c = pnt_2nd + side * c;
    For n In{ from : to }
      lc++;
      pnt_1 = pnt_1st_c  + (n - 1); 
      pnt_2 = pnt_2nd_c + (n - 1);
      Printf("Creating connecting line %g from points: %g, %g",
           lc, pnt_1, pnt_2);
      Line(lc) = {pnt_copy(pnt_1), pnt_copy(pnt_2)};
      lin_type[lc] = type;
    EndFor
  EndFor
Return

// connection of consecutive point around a central one with circular arcs
Macro NArcsAround 
  For c In { 0 : repetition - 1 : step}
    pnt_cent_c = pnt_cent + c;
    pnt_1st_c = pnt_1st + side * c; 
    For n In{ from : to }
      lc++;
      pnt_1 = pnt_1st_c + (n - 1);
      pnt_2 = pnt_1 + 1;
      If(n == to)
        pnt_2 = pnt_1st_c;
      EndIf
      Printf("Creating arc %g from points: %g, %g", lc, pnt_1, pnt_2);
      Circle(lc) = {pnt_copy(pnt_1),
                    pnt_copy(pnt_cent_c),
                    pnt_copy(pnt_2)};
      lin_type[lc] = type;
    EndFor
  EndFor
Return

//------------------------------------------------------------------------------

// orthogonal surfaces inside the core
Macro SurfaceCoreRadial
  For c In { 0 : repetition - 1 : step }
    con_1st_c = con_1st + side/skip * c;
    con_2nd_c = con_2nd + side * c;
    For n In{ from : to }
      sc++;                         
      con_1 = con_1st_c + (n - 1);
      con_2 = con_1 + 1;
      If(n == to)
        con_2 = con_1st_c;
      EndIf
      con_4 = con_2nd_c+ (n - 1) * skip;
      con_3 = con_4 + 1;
      Printf("Creating radial surface %g from lines: %g, %g %g, %g",
              sc, con_1, -con_2, -con_3, -con_4);
      Curve Loop(sc) = {con_1, -con_2, -con_3, -con_4};
      Surface(sc) = {sc};
    EndFor
  EndFor
Return

// orthogonal surfaces in the circual layers
Macro SurfaceRadial
  For c In { 0 : repetition - 1 : step }
    con_1st_c = con_1st + side/skip * c;
    con_2nd_c = con_2nd + side * c;
    con_3rd_c = con_3rd + side * c;
    For n In{ from : to }
      sc++;                         // increase surface counter
      con_1 = con_1st_c + (n - 1);
      con_3 = con_1 + 1;            
      If(n == to)
        con_3 = con_1st_c;
      EndIf
      con_2 = con_2nd_c + (n - 1);
      con_4 = con_3rd_c + (n - 1);
      Printf("Creating radial surface %g from lines: %g, %g %g, %g",
            sc, con_1, con_2, -con_3, -con_4);
      Curve Loop(sc) = {con_1, con_2, -con_3, -con_4};
      Surface(sc) = {sc};
    EndFor
  EndFor
Return

// trasversal surfaces crossing an axis
Macro SurfaceCoreAxial
  For n In{ from : to }
    sc++;
    con_1 = con_1st + (n - 1); 
    con_2 = 1;
    con_3 = con_2nd + (n - 1);
    con_4 = con_3rd + (n - 1) * skip;
    Printf("Creating axial surface %g from lines: %g, %g %g, %g",
            sc, con_1, -con_2, -con_3, con_4);
    Curve Loop(sc) = {con_1, con_2, -con_3, con_4};
    Surface(sc) = {sc};
  EndFor
Return    

// trasversal surfaces 
Macro SurfaceAxial
  For n In{ from : to }
    sc++;
    con_1 = con_1st + (n - 1); 
    con_2 = con_2nd + (n - 1);
    con_3 = con_3rd + (n - 1);
    con_4 = con_4th + (n - 1);   
    Printf("Creating axial surface %g from lines: %g, %g %g, %g",
            sc, con_1, -con_2, -con_3, con_4);
    Curve Loop(sc) = {con_1, -con_2, -con_3, con_4};
    Surface(sc) = {sc};
  EndFor
Return 

// trasversal surfaces around an axis
Macro SurfaceAround
  For n In{ from : to }
    sc++;
    con_1 = con_1st + (n - 1); 
    con_2 = con_2nd + (n - 1);
    con_3 = con_rot + (n - 1);
    con_4 = con_3 + 1;
    If(n == to)
      con_4 = con_rot;
    EndIf
    Printf("Creating tangential surface %g from lines: %g, %g %g, %g",
            sc, con_1, -con_2, con_3, -con_4);
    Curve Loop(sc) = {con_1, -con_2, con_3, -con_4};
    Surface(sc) = {sc};
  EndFor
Return

//------------------------------------------------------------------------------

// volumes inside the core
Macro VolumeCore
  For n In{ from : to }
    vc++;
    sur_1 = sur_rad1 + (n - 1);             // radial start
    sur_2 = sur_rad2 + (n - 1);             // radial end
    sur_3 = sur_tan + (n - 1) * skip;       // tangential start
    sur_4 = sur_3 + 1;                      // tangential end
    sur_5 = sur_axi  + (n - 1);             // axial start
    sur_6 = sur_5 + 1;                      // axial end
    If(n == to)
      sur_6 = sur_axi;
    EndIf
    Printf("Creating core volume %g from surfaces: %g, %g %g, %g, %g, %g",
            vc, sur_1, sur_2, sur_3, sur_4, sur_5, sur_6);
    Surface Loop(vc) = {sur_1, sur_2, sur_3, sur_4, sur_5, sur_6};
    Volume(vc) = {vc};
    Transfinite Volume{vc};
  EndFor
Return 

// volumes outside the core
Macro VolumeOuter
  For n In{ from : to }
    vc++;
    sur_1 = sur_rad1 + (n - 1);             // radial start
    sur_2 = sur_rad2 + (n - 1);             // radial end
    sur_3 = sur_tan1 + (n - 1) * skip;      // tangential start
    sur_4 = sur_tan2 + (n - 1) * skip;      // tangential end
    sur_5 = sur_axi  + (n - 1);             // axial start
    sur_6 = sur_5 + 1;                      // axial end
    If(n == to)
      sur_6 = sur_axi;
    EndIf
    Printf("Creating core volume %g from surfaces: %g, %g %g, %g, %g, %g",
            vc, sur_1, sur_2, sur_3, sur_4, sur_5, sur_6);
    Surface Loop(vc) = {sur_1, sur_2, sur_3, sur_4, sur_5, sur_6};
    Volume(vc) = {vc};
    Transfinite Volume{vc};
  EndFor
Return

//------------------------------------------------------------------------------

// Delete surfaces
Macro SurfaceDelete
  For n In { 0 : repetition - 1 }
    s = start + n;
    Recursive Delete { Surface{s}; }
  EndFor
Return

//==============================================================================
//
// Inizialization
//
//------------------------------------------------------------------------------

// Initialize counters
pc  = 0;  // point count
lc  = 0;  // line count
sc  = 0;  // surface count
vc  = 0;  // volume count

// Initialize list of line types and point coppies
lin_type[] = {};
pnt_copy[] = {};

//==============================================================================
//
// Points
//
//------------------------------------------------------------------------------

// Create points along the main axes
Printf("Creating center points");
center_pnt_start = pc + 1;
pc++;  Point(pc) = {-L/2, 0, 0 };  // inlet center
pc++;  Point(pc) = {+L/2, 0, 0 };  // outlet center

// Points in the core 
core_pnt_start = pc + 1;
Printf("Creating core points from point %g", core_pnt_start);
from = 1;  to = side;
x = -L/2;  r = R(1);  Call NPointsAroundX; 
x = +L/2;  r = R(1);  Call NPointsAroundX; 

// Points in the boundary layer
middle_pnt_start = pc + 1; 
Printf("Creating middlet points from point %g", middle_pnt_start);
from = 1;  to = side;
x = -L/2;  r = R(2);  Call NPointsAroundX;
x = +L/2;  r = R(2);  Call NPointsAroundX;


// Points at the internal wall
blayer_pnt_start = pc + 1;
Printf("Creating blayer points from point %g", blayer_pnt_start);
from = 1;  to = side;
x = -L/2;  r = R(3);  Call NPointsAroundX;
x = +L/2;  r = R(3);  Call NPointsAroundX;

// Merging the points in the same location
Printf("Finishing point creation with point %g", pc);
Call MergePoints;

//==============================================================================
//
// Lines (Order: Radial, Tangential, Axial)
//
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Axis
//------------------------------------------------------------------------------

Printf("Creating axis lines from line %g", lc + 1);
lc++;  Line(lc) = {1, 2};  lin_type(lc) = TYPE_STREAM;

//------------------------------------------------------------------------------
// Core
//------------------------------------------------------------------------------

Printf("Creating core lines from line %g", lc + 1);

// Lines in radial direction 
type = TYPE_CORE; skip = coreRad_skip;
coreRad_lin_start = lc + 1;    
from = 1;  to = side;  repetition = 2;  step = 1;  //function parameters
pnt_cent = center_pnt_start;  pnt_1st = core_pnt_start;   Call LinesCoreRad;  

// Lines in tangential direction 
type = TYPE_ARC; 
coreTan_lin_start = lc + 1;    
from = 1;  to = side;  repetition = 2;  step = 1;  //function parameters
pnt_1st = core_pnt_start;   Call LinesCoreTan;

// Lines in axial direction 
coreAxi_lin_start = lc + 1; 
from = 1;  to = side;  repetition = 1; //function parameters

type = TYPE_STREAM;      
pnt_1st = core_pnt_start + side;  pnt_2nd = core_pnt_start;  Call LinesBtwLayers;

//------------------------------------------------------------------------------
// Middle
//------------------------------------------------------------------------------

Printf("Creating middle lines from line %g", lc + 1);

// Lines in radial direction 
type = TYPE_MIDDLE; skip = middleRad_skip;

middleRad_lin_start = lc + 1;    
from = 1;  to = side;  repetition = 2;  step = 1;  //function parameters
pnt_1st = middle_pnt_start;  pnt_2nd = core_pnt_start;  Call LinesBtwLayers; 

// Lines in tangential direction 
type = TYPE_ARC; 

middleTan_lin_start = lc + 1;    

from = 1;  to = side;  repetition = 2;  step = 1;  //function parameters
pnt_cent = center_pnt_start;  pnt_1st = middle_pnt_start;   Call NArcsAround;  // circle
 
// Lines in axial direction 
middleAxi_lin_start = lc + 1; 
from = 1;  to = side;  repetition = 1; //function parameters

type = TYPE_STREAM;      
pnt_1st = middle_pnt_start + side;  pnt_2nd = middle_pnt_start;  Call LinesBtwLayers;

//------------------------------------------------------------------------------
// Boundary Layer
//------------------------------------------------------------------------------

Printf("Creating boundary layer lines from line %g", lc + 1);

// Lines in radial direction 
type = TYPE_BLAYER; skip = blayerRad_skip;
blayerRad_lin_start = lc + 1;    
from = 1;  to = side;  repetition = 2;  step = 1;  //function parameters
pnt_1st = blayer_pnt_start;  pnt_2nd = middle_pnt_start;  Call LinesBtwLayers;  

// Lines in tangential direction 
type = TYPE_ARC; 
blayerTan_lin_start = lc + 1;    

from = 1;  to = side;  repetition = 2;  step = 1;  //function parameters
pnt_cent = center_pnt_start;  pnt_1st = blayer_pnt_start;  Call NArcsAround;  // circle
 
// Lines in axial direction 
blayerAxi_lin_start = lc + 1; 
from = 1;  to = side;  repetition = 1; //function parameters

type = TYPE_STREAM;      
pnt_1st = blayer_pnt_start + side;  pnt_2nd = blayer_pnt_start;  Call LinesBtwLayers;

//==============================================================================
//
// Set all connectors to be transfinite with
//     proper resolution and progression
//
//------------------------------------------------------------------------------

Printf("Finishing lines creation with line %g", lc);

For l In{ 1 : lc }
  If(lin_type(l) == TYPE_CORE)
    Transfinite Curve {l} = N_CORE      Using Progression P_CORE;
  EndIf
  If(lin_type(l) == TYPE_MIDDLE)
    Transfinite Curve {l} = N_MIDDLE    Using Progression P_MIDDLE;
  EndIf
  If(lin_type(l) == TYPE_BLAYER)
    Transfinite Curve {l} = N_BLAYER    Using Progression P_BLAYER;
  EndIf
  If(lin_type(l) == TYPE_STREAM)
    Transfinite Curve {l} = N_STREAM    Using Progression P_STREAM;
  EndIf
  If(lin_type(l) == TYPE_ARC)
    Transfinite Curve {l} = N_ARC       Using Progression P_ARC;
  EndIf
EndFor

//==============================================================================
//
// Surfaces (Order: Radial, Tangential, Axial)
//
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Core
//------------------------------------------------------------------------------

Printf("Creating core surfaces from surface %g", sc + 1);
skip = coreRad_skip;

// Radial surfaces
coreRad_sur_start = sc + 1; 
from = 1;  to = side/skip;  repetition = 2;  step = 1;  // function parameters
con_1st = coreRad_lin_start;  con_2nd = coreTan_lin_start;  Call SurfaceCoreRadial; 

// Tangential surfaces
coreTan_sur_start = sc + 1; 
from = 1;  to = side;  // function parameters
 
con_1st = coreTan_lin_start;  con_2nd = coreTan_lin_start + side;  
con_rot = coreAxi_lin_start;  Call SurfaceAround; 

// Axial surfaces
coreAxi_sur_start = sc + 1; 
from = 1;  to = side/skip;   

con_1st = coreRad_lin_start;  con_2nd = coreRad_lin_start + side/skip; 
con_3rd = coreAxi_lin_start;  Call SurfaceCoreAxial; 

//------------------------------------------------------------------------------
// Middle
//------------------------------------------------------------------------------

Printf("Creating middle surfaces from surface %g", sc + 1);
skip = middleRad_skip;

// Radial surfaces
middleRad_sur_start = sc + 1;  
from = 1;  to = side/skip;  repetition = 2;  step = 1;  // function parameters
con_1st = middleRad_lin_start;  con_2nd = coreTan_lin_start; 
con_3rd = middleTan_lin_start;  Call SurfaceRadial; 

// Tangential surfaces
middleTan_sur_start = sc + 1; 
from = 1;  to = side;  // function parameters
 
con_1st = middleTan_lin_start;  con_2nd = middleTan_lin_start + side;               
con_rot = middleAxi_lin_start;  Call SurfaceAround; 

// Axial surfaces
middleAxi_sur_start = sc + 1;
from = 1;  to = side/skip;   

con_1st = middleRad_lin_start;  con_2nd = middleRad_lin_start + side/skip;        
con_3rd = coreAxi_lin_start;    con_4th = middleAxi_lin_start ;
Call SurfaceAxial; 

//------------------------------------------------------------------------------
// Boundary Layer
//------------------------------------------------------------------------------

Printf("Creating boundary layer surfaces from surface %g", sc + 1);
skip = blayerRad_skip;

// Radial surfaces
blayerRad_sur_start = sc + 1; 
from = 1;  to = side/skip;  repetition = 2;  step = 1;  // function parameters
con_1st = blayerRad_lin_start;  con_2nd = middleTan_lin_start; 
con_3rd = blayerTan_lin_start;  Call SurfaceRadial;  

// Tangential surfaces
blayerTan_sur_start = sc + 1; 
from = 1;  to = side;  // function parameters
 
con_1st = blayerTan_lin_start;  con_2nd = blayerTan_lin_start + side;               
con_rot = blayerAxi_lin_start;  Call SurfaceAround; 

// Axial surface
blayerAxi_sur_start = sc + 1;
from = 1;  to = side/skip;   

con_1st = blayerRad_lin_start;  con_2nd = blayerRad_lin_start + side/skip;        
con_3rd = middleAxi_lin_start;  con_4th = blayerAxi_lin_start;
Call SurfaceAxial;

//==============================================================================
//
// Set all surfaces to be transfinite
// (Recombine them at once to prevent nodes from moving)
//
//------------------------------------------------------------------------------

Printf("Finishing surfaces creation with surface %g", sc);

For s In{ 1 : sc }
  Transfinite Surface {s};
  Recombine Surface{s};
EndFor

//==============================================================================
//
// Volumes Generation
//
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Core
//------------------------------------------------------------------------------

Printf("Creating core volumes from volume %g", vc + 1);  
skip = coreRad_skip;  core_vol_start = vc + 1;
from = 1;  to = side/skip;  
 
sur_rad1 = coreRad_sur_start;  sur_rad2 = coreRad_sur_start + side/skip;
sur_tan  = coreTan_sur_start;  sur_axi  = coreAxi_sur_start;  
Call VolumeCore; 

//------------------------------------------------------------------------------
// Middle
//------------------------------------------------------------------------------

Printf("Creating middle volumes from volume %g", vc + 1); 
skip = middleRad_skip;  middle_vol_start = vc + 1; 
from = 1;  to = side/skip;  

sur_rad1 = middleRad_sur_start;  sur_rad2 = middleRad_sur_start + side/skip;
sur_tan1 = coreTan_sur_start;    sur_tan2 = middleTan_sur_start;
sur_axi  = middleAxi_sur_start;  
Call VolumeOuter;  

//------------------------------------------------------------------------------
// Boundary Layer
//------------------------------------------------------------------------------

Printf("Creating boundary layer volumes from volume %g", vc + 1); 
skip = blayerRad_skip;  blayer_vol_start = vc + 1;
from = 1;  to = side/skip;  

sur_rad1 = blayerRad_sur_start;  sur_rad2 = blayerRad_sur_start + side/skip;
sur_tan1 = middleTan_sur_start;  sur_tan2 = blayerTan_sur_start;
sur_axi  = blayerAxi_sur_start;  
Call VolumeOuter;  

Coherence;

//==============================================================================
//
// Set physical volumes and surfaces
//
//------------------------------------------------------------------------------

Printf("Solid finished with volume %g", vc);

Physical Volume("FLUID") = {1 : vc};
  
Physical Surface("INTERNAL_WALL") = { blayerTan_sur_start : blayerAxi_sur_start - 1 };
  
Physical Surface("INLET")  = { coreRad_sur_start   : coreRad_sur_start   + side/2 - 1,
                               middleRad_sur_start : middleRad_sur_start + side   - 1,
                               blayerRad_sur_start : blayerRad_sur_start + side   - 1 };
                               
Physical Surface("OUTLET")  = {coreRad_sur_start   + side/2 : coreRad_sur_start   + 2*side/2 - 1,
                               middleRad_sur_start + side   : middleRad_sur_start + 2*side   - 1,
                               blayerRad_sur_start + side   : blayerRad_sur_start + 2*side   - 1 };
                                                


