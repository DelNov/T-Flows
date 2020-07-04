/*******************************************************************************
*                                                                              *
*  Cylinder mesh                                                               *
*                                                                              *
*******************************************************************************/

//------------------------------------------------------------------------------
//
// Parameters defining domain extents, buildings and mesh resolution
// (This should be fiddled with, but with care)
//
//------------------------------------------------------------------------------

n_box_layers =   6;                // number of mesh layers in height
n_bnd_layers =  11;                // number of node layers in boundary layer
depth_z      =   1.0;              // depth of the domain

// Coordinates of the problem domain (the whole piece of simulated land)
box_x_min = 0.0;
box_x_max = 4.0;
box_y_min = 0.0;
box_y_max = 4.0;

cyl_x[1] = (box_x_max - box_x_min) * 0.33333333;
cyl_y[1] = (box_y_max - box_y_min) * 0.33333333;
cyl_x[2] = (box_x_max - box_x_min) * 0.66666667;
cyl_y[2] = (box_y_max - box_y_min) * 0.66666667;
cyl_r = 0.5;
cyl_b = 0.3;

// Transition between fine and coarse mesh, ...
// ... between city and the rest of the domain
delta_min = 0.03;   // resolution in the wake
delta_max = 0.10;   // resolution in the box
wake_limit_width = 0.5;

//------------------------------------------------------------------------------
//
// Parameters for the problem definition algorithms
// (These are better kept untouched)
//
//------------------------------------------------------------------------------
COS_45       = 0.70710678118;
BOX_SURF     = 100;                // ground surface number
LAYER_SURF   = 200;                // ground surface number
TINY         =   1.0e-9;
HUGE         =   1.0e+9;
CYL_POINT    = 100;       // starting point for cylinder definition
LAY_POINT    = 200;       // starting point for boundary layer
CEN_POINT    = 301;       // central point for cylinder
BOX_LINE     = 100;
LAY_LINE     = 200;
CYL_CIRC     = 300;
LAY_CIRC     = 400;
C            =   2;       // number of circles
N            =   3;       // number of circle segments

// Define a closed chain over N segment
For n In { 1 : N }
  NEXT(n) = n + 1;
EndFor
NEXT(N) = 1;

//------------------------------------------------------------
// Points and lines defining the extents of the ground domain
//------------------------------------------------------------
Point(1) = {box_x_min, box_y_min, 0};
Point(2) = {box_x_max, box_y_min, 0};
Point(3) = {box_x_max, box_y_max, 0};
Point(4) = {box_x_min, box_y_max, 0};
Line(BOX_LINE+1) = {1, 2};
Line(BOX_LINE+2) = {2, 3};
Line(BOX_LINE+3) = {3, 4};
Line(BOX_LINE+4) = {4, 1};

// Browse through all circles
For c In{ 1 : C }

  alpha_o =       Pi / N;
  alpha_s = 2.0 * Pi / N;

  // Inner circle points
  For n In{ 1 : N }
    alpha = alpha_o + (n-1) * alpha_s;
    Point(CYL_POINT+(c-1)*N+n) = {cyl_x[c] - cyl_r * Cos(alpha),
                                  cyl_y[c] - cyl_r * Sin(alpha), 0};
  EndFor

  // Boundary layer circle points
  For n In{ 1 : N }
    alpha = alpha_o + (n-1) * alpha_s;
    Point(LAY_POINT+(c-1)*N+n) = {cyl_x[c] - (cyl_r+cyl_b)*Cos(alpha),
                                  cyl_y[c] - (cyl_r+cyl_b)*Sin(alpha), 0};
  EndFor

  // Center of the circle
  Point(CEN_POINT + c) = {cyl_x[c], cyl_y[c], 0};

  // Define cylinder surface
  For n In{ 1 : N }
    Circle(CYL_CIRC+(c-1)*N+n) = {CYL_POINT+(c-1)*N+n,
                                  CEN_POINT+c,
                                  CYL_POINT+(c-1)*N+NEXT(n)};
  EndFor

  // Define outer edge of boundary layer
  For n In{ 1 : N }
    Circle(LAY_CIRC+(c-1)*N+n) = {LAY_POINT+(c-1)*N+n,
                                  CEN_POINT+c,
                                  LAY_POINT+(c-1)*N+NEXT(n)};
  EndFor

  // Lines in the boundary layer
  For n In{ 1 : N }
    Line(LAY_LINE+(c-1)*N+n) = {CYL_POINT+(c-1)*N+n, LAY_POINT+(c-1)*N+n};
  EndFor
EndFor

//------------------------------------------------------------------------------
//
// Define surfaces
//
//------------------------------------------------------------------------------

// Box surface
Curve Loop(1) = {BOX_LINE+1:BOX_LINE+4};
Curve Loop(2) = {LAY_CIRC+1:LAY_CIRC+C*N};
Plane Surface(BOX_SURF) = {1, 2};

// Boundary layer surfaces; Browse through all circles
For c In{ 1 : C }
  For n In{ 1 : N }
    Curve Loop(2+(c-1)*N+n) = {(CYL_CIRC+(c-1)*N+n),
                               (LAY_LINE+(c-1)*N+NEXT(n)),
                              -(LAY_CIRC+(c-1)*N+n),
                              -(LAY_LINE+(c-1)*N+n)};
    Plane Surface(LAYER_SURF+(c-1)*N+n) = {2+(c-1)*N+n};
  EndFor
EndFor

//----------------------------------------
// Define ground surface witout buildings
//----------------------------------------
Printf("Defined ground surface %g", BOX_SURF);

//------------------------------------------------------------------------------
//
// Define mesh
//
//------------------------------------------------------------------------------

//----------------------------
// First the spacing function
//----------------------------

// Option 2: with Cylinder
For c In{ 1 : C }
  Field[c]         = Cylinder;
  Field[c].VIn     = delta_min;
  Field[c].VOut    = delta_max;
  Field[c].Radius  = cyl_r + wake_limit_width;
  Field[c].XCenter = cyl_x[c];
  Field[c].YCenter = cyl_y[c];
EndFor

Background Field = 2;
Printf("Spacing function with delta_min - delta_max: %g - %g defined!",
       delta_min, delta_max);

//---------------------------------------------------
// Generate mesh on the ground, and on all the bases
//---------------------------------------------------
Recombine Surface{BOX_SURF};

// This experimental algorithm could give better results
Mesh.Algorithm = 8;

// Uncomment the following line to try the full-quad algorithm:
// Mesh.RecombinationAlgorithm = 2; // or 3

p = 0.99;  // initial progression
d = 0.05;  // initial increment
l = delta_min * (1 - (1.0/p)^(n_bnd_layers)) / (1 - (1.0/p));
Printf("Length with progression %g is %g", p, l);

For iter In{1:20}
  If(l > cyl_b)
    p = p + d;
    l = delta_min * (1 - (1.0/p)^(n_bnd_layers)) / (1 - (1.0/p));
    Printf("Length with progression %g is %g", p, l);
    If(l < cyl_b) d = d * 0.5; EndIf
  EndIf
  If(l < cyl_b)
    p = p - d;
    l = delta_min * (1 - (1.0/p)^(n_bnd_layers)) / (1 - (1.0/p));
    Printf("Length with progression %g is %g", p, l);
    If(l > cyl_b) d = d * 0.5; EndIf
  EndIf
EndFor

For c In{ 1 : C }
  Transfinite Curve {LAY_LINE+(c-1)*N+1:LAY_LINE+(c-1)*N+N}
   = n_bnd_layers Using Progression p;
EndFor

  s = 2.0 * (cyl_r+cyl_b) * Pi / 4.0;
  Printf("Length of one segment: %g", s);
  n = s / delta_min + 1.0;
  Printf("Number of cells in one segment: %g", n);

For c In{ 1 : C }
  Transfinite Curve {LAY_CIRC+(c-1)*N+1 : LAY_CIRC+(c-1)*N+N}
   = n Using Progression 1;
  Transfinite Curve {CYL_CIRC+(c-1)*N+1 : CYL_CIRC+(c-1)*N+N}
   = n Using Progression 1;
EndFor

For c In{ 1 : C }
  For n In{ 1 : N }
    Transfinite Surface {LAYER_SURF+(c-1)*N+n} = {CYL_POINT+(c-1)*N+n,
                                                  LAY_POINT+(c-1)*N+n,
                                                  LAY_POINT+(c-1)*N+NEXT(n),
                                                  CYL_POINT+(c-1)*N+NEXT(n)};
    Recombine Surface {LAYER_SURF+(c-1)*N+n};
  EndFor
EndFor

//------------------------------------------------------------------------------
//
// Create volume
//
//------------------------------------------------------------------------------
Extrude {0, 0, depth_z} {
   Surface{BOX_SURF, LAYER_SURF+1:LAYER_SURF+C*N};
   Layers{n_box_layers};
   Recombine;
}

//-----------------------------------------------------------------------------t
//
// Define boundary conditions
//
//------------------------------------------------------------------------------

Physical Surface("periodic_x")
  = {Surface In BoundingBox{box_x_min-TINY, -HUGE, -HUGE,
                            box_x_min+TINY, +HUGE, +HUGE},
     Surface In BoundingBox{box_x_max-TINY, -HUGE, -HUGE,
                            box_x_max+TINY, +HUGE, +HUGE}};

Physical Surface("bottom")
  = {Surface In BoundingBox{-HUGE, box_y_min-TINY, -HUGE,
                            +HUGE, box_y_min+TINY, +HUGE}};

Physical Surface("top")
  = {Surface In BoundingBox{-HUGE, box_y_max-TINY, -HUGE,
                            +HUGE, box_y_max+TINY, +HUGE}};

Physical Surface("periodic_z")
  = {Surface In BoundingBox{-HUGE, -HUGE, depth_z-TINY,
                            +HUGE, +HUGE, depth_z+TINY},
     Surface In BoundingBox{-HUGE, -HUGE, 0.0    -TINY,
                            +HUGE, +HUGE, 0.0    +TINY}};

If(N == 3)
  Physical Surface("cylinder") = {489, 467, 511, 533, 555, 577};
EndIf

//------------------------------------------------------------------------------
//
// Create pysical volume
//
//------------------------------------------------------------------------------
Physical Volume("interior") = { 1:1+C*N };
