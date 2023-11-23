/*******************************************************************************
*                                                                              *
*  Benchmark city mesh                                                         *
*                                                                              *
*  This script reads file "5_obstacles.geo"                                    *
*                                                                              *
*******************************************************************************/

//------------------------------------------------------------------------------
//
// Parameters defining domain extents, obstacles and mesh resolution
// (This should be fiddled with, but with care)
//
//------------------------------------------------------------------------------

PERIODIC  =   0;  // or 0
ANGLE_DEG =   0.0;

// Two different types of obstacles (buildings or porous zones)
SOLID  = 1;
POROUS = 0;

// Number of layers
N_LAYERS     = 70;
N_SKY_LAYERS = 35;

// Building characteristic dimension and the pitch between them
A = 50.0;
P = 30.0;

// Height of the volume
SKY_HIGH    = 500.0;      // max height of the domain
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
DELTA_BUILDING = DELTA_MIN;
DELTA_OBSTACLE = DELTA_MAX;

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
BASE_LOOP_START     = GROUND_LOOP + 500;  // obstacles' base definitions
GROUND_SURF         = 100;                // ground surface number
MAX_BUILDING_HEIGHT = 300;
MAXN                =   8;                // max nodes per obstacle
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
// Define all obstacles' positions
//---------------------------------
Printf("Including file 5_obstacles.geo");
Include "5_obstacles.geo";

// Rotate them
For o In { 1 : n_obstacles }
  For n In{ 1 : node_b(o) }
    i = MAXN*o + n;
    x(i) = xo(i) * Cos(ANGLE_RAD) - yo(i) * Sin(ANGLE_RAD);
    y(i) = xo(i) * Sin(ANGLE_RAD) + yo(i) * Cos(ANGLE_RAD);
  EndFor
EndFor

//-----------------------------------------------------------------------
// Browse through all obstacles to define points and lines defining them
//-----------------------------------------------------------------------
Printf("Defining all obstacles");
For o In { 1 : n_obstacles }

  // Browse through nodes of the obstacle "o"
  For n In{ 1 : node_b(o) }
    delta = 0.0;  // ?DELTA_MAX; I hope it neglects it if it's zero
    If(type_o(o) == SOLID)
      delta = DELTA_MIN;
      If( d(MAXN*o + n)/3.0 < delta )
        delta = d(MAXN*o + n)/6.0;
      EndIf
    EndIf
    Point(MAXN*o + n) = {x(MAXN*o + n),
                         y(MAXN*o + n),
                         0,
                         delta};
  EndFor

  // Lines defining the obstacle
  For n In{ 1 : node_b(o)-1 }
    Line(MAXN*o + n) = {MAXN*o + n, MAXN*o + n + 1};
  EndFor
  Line(MAXN*o + node_b(o)) = {MAXN*o + node_b(o), MAXN*o + 1};

  // Anticlockwise curves (holes)
  Curve Loop(HOLE_LOOP_START + o) = {MAXN*o + 1 : MAXN*o + node_b(o)};

  // Clockwise curves (buidlings' bases)
  Curve Loop(BASE_LOOP_START + o) = {MAXN*o + 1 : MAXN*o + node_b(o)};
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
If(n_obstacles > 0)
  Plane Surface(GROUND_SURF) = {GROUND_LOOP,
                                HOLE_LOOP_START + 1:
                                HOLE_LOOP_START + n_obstacles};
Else
  Plane Surface(GROUND_SURF) = {GROUND_LOOP};
EndIf

//------------------------------------
// Add individual obstacles' surfaces
//------------------------------------
For o In { 1 : n_obstacles }
  Printf("Defining obstacle surface %g", GROUND_SURF + o);
  Plane Surface(GROUND_SURF + o) = {BASE_LOOP_START + o};
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
For o In { 1 : n_obstacles }
  Recombine Surface{GROUND_SURF + o};
EndFor

// This experimental algorithm could give better results
Mesh.Algorithm = 8;

// Uncomment the following line to try the full-quad algorithm:
// Mesh.RecombinationAlgorithm = 2; // or 3

//---------------
// Create volume
//---------------

Extrude {0,0, URBAN_HIGH} {
  Surface {GROUND_SURF:GROUND_SURF+n_obstacles};
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
// Obstacles
//-----------

nb = 0;
np = 0;
For o In { 1 : n_obstacles }            // browse through all buildings
  For h In { 0 : MAX_BUILDING_HEIGHT }  // browse through all heights

    // Add the obstacles whose height matches "h"
    If(height_o(o) == h)  // found an obstacle with matching height
      If(type_o(o) == SOLID)
        nb++;
        Physical Surface(Sprintf("building_%03g_%03g",h,nb)) = {GROUND_SURF+o};
      Else
        np++;
        Physical Surface(Sprintf("porosity_%03g_%03g",h,np)) = {GROUND_SURF+o};
      EndIf
    EndIf

  EndFor    // through all heights
EndFor  // through all obstacles

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
Physical Volume("interior") = { 1 : 1 + n_obstacles + 10000 };
