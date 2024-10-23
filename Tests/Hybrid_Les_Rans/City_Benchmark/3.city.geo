/*******************************************************************************
*                                                                              *
*  Benchmark city mesh                                                         *
*                                                                              *
*  This script reads file "3_building.geo"                                     *
*                                                                              *
*******************************************************************************/

//------------------------------------------------------------------------------
//
// Parameters defining domain extents, buildings and mesh resolution
// (This should be fiddled with, but with care)
//
//------------------------------------------------------------------------------

PERIODIC  =   0;  // or 0
ANGLE_DEG =   0.0;

// Number of layers
N_LAYERS     = 70;
N_SKY_LAYERS = 35;

// Building characteristic dimension and the pitch between them
A = 50.0;
P = 30.0;

// Height of the volume
SKY_HIGH    = 550.0;      // max height of the domain
URBAN_HIGH  = 175.0;      // height 25.0 above the buildings

DELTA_H_MIN = URBAN_HIGH / N_LAYERS;
L_TARGET    = SKY_HIGH - URBAN_HIGH;

// Coordinates of the problem domain (the whole piece of simulated land)
GROUND_X_MIN =  -500.0;
GROUND_X_MAX =  1500.0;
GROUND_Y_MIN =  -500.0;
GROUND_Y_MAX =   500.0;

// Coordinates of the city (where buildings will reside)
CITY_X_MIN = -275.0;
CITY_X_MAX =  725.0;
CITY_Y_MIN = -250.0;
CITY_Y_MAX =  250.0;

// Resolutions in the city (min) and country side (max)
DELTA_MIN      = 10.0;
DELTA_MAX      = 25.0;
DELTA_BUILDING = 10.0;

// Transition between fine and coarse mesh, ...
// ... between city and the rest of the domain
CITY_LIMIT_WIDTH = 40.0;

//------------------------------------------------------------------------------
//
// Parameters for the problem definition algorithms
// (These are better kept untouched)
//
//------------------------------------------------------------------------------
GROUND_LOOP         = 100;                // ground loop number
HOLE_LOOP_START     = GROUND_LOOP + 100;  // hole definitions
BASE_LOOP_START     = GROUND_LOOP + 500;  // buildings' base definitions
GROUND_SURF         = 100;                // ground surface number
MAX_BUILDING_HEIGHT = 300;
MAXN                =   8;                // max nodes per building
TINY                =   1.0e-9;
HUGE                =   1.0e+9;
ANGLE_RAD           = ANGLE_DEG * Pi / 180.0;

//------------------------------------------------------------------------------
//
// Points and lines defining the extents of the ground domain
//
//------------------------------------------------------------------------------

//SetFactory("OpenCASCADE");

Printf("Defining ground");
Point(1) = {GROUND_X_MIN, GROUND_Y_MIN, 0};
Point(2) = {GROUND_X_MAX, GROUND_Y_MIN, 0};
Point(3) = {GROUND_X_MAX, GROUND_Y_MAX, 0};
Point(4) = {GROUND_X_MIN, GROUND_Y_MAX, 0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(GROUND_LOOP) = {1, 2, 3, 4};

//---------------------------------
// Define all buildings' positions
//---------------------------------
Printf("Including file 3_buildings.geo");
Include "3_buildings.geo";

// Rotate them
For b In { 1 : n_buildings }
  For n In{ 1 : node_b(b) }
    i = MAXN*b + n;
    x(i) = xo(i) * Cos(ANGLE_RAD) - yo(i) * Sin(ANGLE_RAD);
    y(i) = xo(i) * Sin(ANGLE_RAD) + yo(i) * Cos(ANGLE_RAD);
  EndFor
EndFor

//-----------------------------------------------------------------------
// Browse through all buildings to define points and lines defining them
//-----------------------------------------------------------------------
Printf("Defining all buildings");
For b In { 1 : n_buildings }

  // Browse through nodes of the building "b"
  For n In{ 1 : node_b(b) }
    delta = DELTA_MIN;
    If( d(MAXN*b + n)/3.0 < delta )
      delta = d(MAXN*b + n)/6.0;
    EndIf
    Point(MAXN*b + n) = {x(MAXN*b + n),
                         y(MAXN*b + n),
                         0,
                         delta};
  EndFor

  // Lines defining the building
  For n In{ 1 : node_b(b)-1 }
    Line(MAXN*b + n) = {MAXN*b + n, MAXN*b + n + 1};
  EndFor
  Line(MAXN*b + node_b(b)) = {MAXN*b + node_b(b), MAXN*b + 1};

  // Anticlockwise curves (holes)
  Curve Loop(HOLE_LOOP_START + b) = {MAXN*b + 1 : MAXN*b + node_b(b)};

  // Clockwise curves (buidlings' bases)
  Curve Loop(BASE_LOOP_START + b) = {MAXN*b + 1 : MAXN*b + node_b(b)};
EndFor


//------------------------------------------------------------------------------
//
// Define surfaces
//
//------------------------------------------------------------------------------


//----------------------------------------
// Define ground surface witout buildings
//----------------------------------------
Printf("Defining ground surface %g", GROUND_SURF);
If(n_buildings > 0)
  Plane Surface(GROUND_SURF) = {GROUND_LOOP,
                                HOLE_LOOP_START + 1:
                                HOLE_LOOP_START + n_buildings};
Else
  Plane Surface(GROUND_SURF) = {GROUND_LOOP};
EndIf

//-----------------------------------
// Add individual buildings surfaces
//-----------------------------------
For b In { 1 : n_buildings }
  Printf("Defining building surface %g", GROUND_SURF + b);
  Plane Surface(GROUND_SURF + b) = {BASE_LOOP_START + b};
EndFor


//------------------------------------------------------------------------------
//
// Define mesh
//
//------------------------------------------------------------------------------


//----------------------------
// First the spacing function
//----------------------------
Printf("Defining spacing function with DELTA_MIN - DELTA_MAX: %g - %g",
       DELTA_MIN, DELTA_MAX);
Field[1] = MathEval;
Field[1].F = Sprintf("  (%5.2g)
                      + (%5.2g) * (1.0-(0.5*(  tanh((x-(%5.2g))/(%5.2g))
                                             - tanh((x-(%5.2g))/(%5.2g))))
                                      *(0.5*(  tanh((y-(%5.2g))/(%5.2g))
                                             - tanh((y-(%5.2g))/(%5.2g)))) )",
  DELTA_MIN,
  DELTA_MAX - DELTA_MIN,
  CITY_X_MIN, CITY_LIMIT_WIDTH,
  CITY_X_MAX, CITY_LIMIT_WIDTH,
  CITY_Y_MIN, CITY_LIMIT_WIDTH,
  CITY_Y_MAX, CITY_LIMIT_WIDTH);
Background Field = 1;


//---------------------------------------------------
// Generate mesh on the ground, and on all the bases
//---------------------------------------------------
Recombine Surface{GROUND_SURF};
For b In { 1 : n_buildings }
  Recombine Surface{GROUND_SURF + b};
EndFor

// This experimental algorithm could give better results
// Mesh.Algorithm = 8;

// Uncomment the following line to try the full-quad algorithm:
// Mesh.RecombinationAlgorithm = 2; // or 3

//---------------
// Create volume
//---------------

Extrude {0,0, URBAN_HIGH} {
  Surface {GROUND_SURF:GROUND_SURF+n_buildings};
  Layers {N_LAYERS};
  Recombine;
}

//--------------------------------------------
// Progressivelly expand more towards the sky
//--------------------------------------------
p           = 0.99;  // initial progression
d           = 0.05;  // initial increment in progression

// Initial length
l = DELTA_H_MIN * (1 - p^(N_SKY_LAYERS)) / (1 - p);
Printf("Length with progression %g is %g", p, l);

// Length through iterations
For iter In{1:40}
  If(l > L_TARGET)  // if current length is larger than the target ...
    p = p - d;      // ... reduce the progression
    l = DELTA_H_MIN * (1 - p^(N_SKY_LAYERS)) / (1 - p);
    If(l < L_TARGET) d = d * 0.5; EndIf
  EndIf
  If(l < L_TARGET)  // if current length is smaller than the target ...
    p = p + d;      // ... increase the progression
    l = DELTA_H_MIN * (1 - p^(N_SKY_LAYERS)) / (1 - p);
    If(l > L_TARGET) d = d * 0.5; EndIf
  EndIf
  Printf("Iteration %g; Length with progression %g is %g", iter, p, l);
EndFor

For lay In{ 1 : N_SKY_LAYERS }
  If(lay == 1)
    delta_cur = DELTA_H_MIN;  // current delta
    z_cur     = URBAN_HIGH;
  Else
    z_cur     += delta_cur;
    delta_cur *= p;
  EndIf
  Printf("Layer %g; height %g", lay, z_cur);
  Extrude {0, 0, delta_cur} {
    Surface{Surface In BoundingBox{-HUGE, -HUGE, z_cur-TINY,
                                   +HUGE, +HUGE, z_cur+TINY}};
    Layers {1};
    Recombine;
  }
EndFor

//------------------------------------------------------------------------------
//
// Define boundary conditions
//
//------------------------------------------------------------------------------

//--------
// Ground
//--------
Physical Surface("ground") = {GROUND_SURF};

//-----------
// Buildings
//-----------

For h In { 0 : MAX_BUILDING_HEIGHT }  // browse through all heights
  exists = 0;

  For b In { 1 : n_buildings }        // browse through all buildings

    // Add the buildings whose height matches "h"
    If(height_b(b) == h)  // found a building with matching height
      If(exists == 0)     // is it the first building at this height
        exists = 1;       // not any more
        Physical Surface(Sprintf("building_%03g_%03g",h,b)) = {GROUND_SURF+b};
      Else                // not the first building at this height, append
        Physical Surface(Sprintf("building_%03g_%03g",h,b)) += {GROUND_SURF+b};
      EndIf
    EndIf

  EndFor
EndFor

//---------------------------------------------
// East, west, south, north and eventually top
//---------------------------------------------

If( PERIODIC == 0 )
  Physical Surface("west")
    = {Surface In BoundingBox{GROUND_X_MIN-TINY, -HUGE, -HUGE,
                              GROUND_X_MIN+TINY, +HUGE, +HUGE}};
  Physical Surface("east")
    = {Surface In BoundingBox{GROUND_X_MAX-TINY, -HUGE, -HUGE,
                              GROUND_X_MAX+TINY, +HUGE, +HUGE}};
  Physical Surface("south")
    = {Surface In BoundingBox{-HUGE, GROUND_Y_MIN-TINY, -HUGE,
                              +HUGE, GROUND_Y_MIN+TINY, +HUGE}};
  Physical Surface("north")
    = {Surface In BoundingBox{-HUGE, GROUND_Y_MAX-TINY, -HUGE,
                              +HUGE, GROUND_Y_MAX+TINY, +HUGE}};
Else
  Physical Surface("east-west")
    = {Surface In BoundingBox{GROUND_X_MIN-TINY, -HUGE, -HUGE,
                              GROUND_X_MIN+TINY, +HUGE, +HUGE}};
    + {Surface In BoundingBox{GROUND_X_MAX-TINY, -HUGE, -HUGE,
                              GROUND_X_MAX+TINY, +HUGE, +HUGE}};
  Physical Surface("north-south")
    = {Surface In BoundingBox{-HUGE, GROUND_Y_MIN-TINY, -HUGE,
                              +HUGE, GROUND_Y_MIN+TINY, +HUGE}};
    + {Surface In BoundingBox{-HUGE, GROUND_Y_MAX-TINY, -HUGE,
                              +HUGE, GROUND_Y_MAX+TINY, +HUGE}};
EndIf
Physical Surface("top")
  = {Surface In BoundingBox{-HUGE, -HUGE, SKY_HIGH-TINY,
                            +HUGE, +HUGE, SKY_HIGH+TINY}};

//------------------------------------------------------------------------------
//
// Create pysical volume
//
//------------------------------------------------------------------------------
Physical Volume("interior") = { 1 : 1 + n_buildings + 10000 };
