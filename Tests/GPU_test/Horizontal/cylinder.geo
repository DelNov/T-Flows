//------------------------------------------------------------------------------
//
// Parameters defining domain extents, buildings and mesh resolution
// (This should be fiddled with, but with care)
//
//------------------------------------------------------------------------------

n_box_layers =   4;                // number of mesh layers in depth
n_bnd_layers =  21;                // number of node layers in boundary layer
box_z_max    =   0.4;              // depth of the domain

// Coordinates of the problem domain (the whole piece of simulated land)
box_x_min =  0.0;
box_x_max =  2.2;
box_y_min =  0.0;
box_y_max =  0.41;

cyl_x = 0.2;
cyl_y = cyl_x;  // asymmetric on purpose!  (box_y_max - box_y_min) / 2.0;
cyl_r = 0.05;
cyl_b = 0.02;

// Transition between fine and coarse mesh, ...
// ... between city and the rest of the domain
delta_min = 0.003; // resolution in the wake
delta_max = 0.015; // resolution in the box
wake_limit_width = 0.05;

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
CYL_POINT    =   5;       // starting point for cylinder definition
LAY_POINT    =  50;       // starting point for boundary layer
CEN_POINT    = 101;       // central point for cylinder
BOX_LINE     =   0;
LAY_LINE     =  50;
CYL_CIRC     = 100;
LAY_CIRC     = 150;
NEXT(1) = 2; NEXT(2) = 3; NEXT(3) = 4; NEXT(4) = 1;

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

// Inner circle points
Point(CYL_POINT+1) = {cyl_x - cyl_r*COS_45, cyl_y - cyl_r*COS_45, 0};
Point(CYL_POINT+2) = {cyl_x + cyl_r*COS_45, cyl_y - cyl_r*COS_45, 0};
Point(CYL_POINT+3) = {cyl_x + cyl_r*COS_45, cyl_y + cyl_r*COS_45, 0};
Point(CYL_POINT+4) = {cyl_x - cyl_r*COS_45, cyl_y + cyl_r*COS_45, 0};

// Boundary layer circle points
Point(LAY_POINT+1) = {cyl_x - (cyl_r+cyl_b)*COS_45,
                      cyl_y - (cyl_r+cyl_b)*COS_45, 0};
Point(LAY_POINT+2) = {cyl_x + (cyl_r+cyl_b)*COS_45,
                      cyl_y - (cyl_r+cyl_b)*COS_45, 0};
Point(LAY_POINT+3) = {cyl_x + (cyl_r+cyl_b)*COS_45,
                      cyl_y + (cyl_r+cyl_b)*COS_45, 0};
Point(LAY_POINT+4) = {cyl_x - (cyl_r+cyl_b)*COS_45,
                      cyl_y + (cyl_r+cyl_b)*COS_45, 0};

// Center of the circle
Point(CEN_POINT) = {cyl_x, cyl_y, 0};

// Define cylinder surface
For n In{ 1 : 4 }
  Circle(CYL_CIRC+n) = {CYL_POINT+n, CEN_POINT, CYL_POINT+NEXT(n)};
EndFor

// Define outer edge of boundary layer
For n In{ 1 : 4 }
  Circle(LAY_CIRC+n) = {LAY_POINT+n, CEN_POINT, LAY_POINT+NEXT(n)};
EndFor

// Lines in the boundary layer
For n In{ 1 : 4 }
  Line(LAY_LINE+n) = {CYL_POINT+n, LAY_POINT+n};
EndFor

//------------------------------------------------------------------------------
//
// Define surfaces
//
//------------------------------------------------------------------------------

// Box surface
Curve Loop(1) = {BOX_LINE+1:BOX_LINE+4};
Curve Loop(2) = {LAY_CIRC+1:LAY_CIRC+4};
Plane Surface(BOX_SURF) = {1, 2};

// Boundary layer surfaces
For n In{ 1 : 4 }
  Curve Loop(2+n) = {(CYL_CIRC+n),
                     (LAY_LINE+NEXT(n)),
                    -(LAY_CIRC+n),
                    -(LAY_LINE+n)};
  Plane Surface(LAYER_SURF+n) = {2+n};
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

Field[2]         = Cylinder;
Field[2].VIn     = delta_min;
Field[2].VOut    = delta_max;
Field[2].Radius  = cyl_r + wake_limit_width;
Field[2].XCenter = cyl_x;
Field[2].YCenter = cyl_y;

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

Transfinite Curve {LAY_LINE+1:LAY_LINE+4} = n_bnd_layers Using Progression p;

s = 2.0 * (cyl_r+cyl_b) * Pi / 4.0;
Printf("Length of one segment: %g", s);
n = s / delta_min + 1.0;
Printf("Number of cells in one segment: %g", n);
Transfinite Curve {LAY_CIRC+1:LAY_CIRC+4} = n Using Progression 1;
Transfinite Curve {CYL_CIRC+1:CYL_CIRC+4} = n Using Progression 1;

For n In{ 1 : 4 }
  Transfinite Surface {LAYER_SURF+n} = {CYL_POINT+n,
                                        LAY_POINT+n,
                                        LAY_POINT+NEXT(n),
                                        CYL_POINT+NEXT(n)};

  Recombine Surface {LAYER_SURF+n};

EndFor

//------------------------------------------------------------------------------
//
// Create volume for cylinder
//
//------------------------------------------------------------------------------
Extrude {0, 0, box_z_max} {
  Surface{BOX_SURF, LAYER_SURF+1:LAYER_SURF+4};
  Layers{n_box_layers};
  Recombine;
}

//------------------------------------------------------------------------------
//
// Define boundary conditions
//
//------------------------------------------------------------------------------
Physical Surface("cylinder_wall") = {321, 299, 255, 277};
Physical Surface("in")
  = {Surface In BoundingBox{box_x_min-TINY, -HUGE, -HUGE,
                            box_x_min+TINY, +HUGE, +HUGE}};

Physical Surface("out")
  = {Surface In BoundingBox{box_x_max-TINY, -HUGE, -HUGE,
                            box_x_max+TINY, +HUGE, +HUGE}};

Physical Surface("flat_walls")
  = {Surface In BoundingBox{-HUGE, box_y_min-TINY, -HUGE,
                            +HUGE, box_y_min+TINY, +HUGE},
     Surface In BoundingBox{-HUGE, box_y_max-TINY, -HUGE,
                            +HUGE, box_y_max+TINY, +HUGE}};

Physical Surface("periodic_z")
  = {Surface In BoundingBox{-HUGE, -HUGE, box_z_max-TINY,
                            +HUGE, +HUGE, box_z_max+TINY},
     Surface In BoundingBox{-HUGE, -HUGE, 0.0      -TINY,
                            +HUGE, +HUGE, 0.0      +TINY}};

//------------------------------------------------------------------------------
//
// Create pysical volumes
//
//------------------------------------------------------------------------------
Physical Volume("cylinder") = { 1:5 };

