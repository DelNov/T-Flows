radius = 1.0;
length = 1.8;
width  = 0.3;
inner  = 0.8 * radius;

//------------------------------------------------------------------------------
// Variants of the geometry:
//
// 12 - upper domain (with two legs)
// 21 - lower domain with one outlet leg
// 22 - lower domain with two outlet legs
// 30 - membrane (without legs)
//------------------------------------------------------------------------------
VARIANT = 30;

n_width = 11;    // number of cells per width

// Work out resolutions from n_width
delta    = width / n_width;
angle    = 4.0 * Asin(width / radius);
fraction = angle / (2.0*Pi - angle);
n_circ   = n_width / fraction;

If(VARIANT == 12)
  z_min  =  0.01;
  z_max  = +0.15;
  n_lay  = 11;
ElseIf(VARIANT == 21 || VARIANT == 22)
  z_min  = -0.25;
  z_max  = -0.01;
  n_lay  = 21;
ElseIf(VARIANT == 30)
  z_min  = -0.01;
  z_max  = +0.01;
  n_lay  =  5;
EndIf

PNT_CENT  =  20;
PNT_DOM   =  40; // starting point for domain
PNT_INN   =  50; // starting point for domain
PNT_LEGS  =  60; // starting point for legs
CIRC_OUT  =  80;
CIRC_INN  =  90;
LINE_LEG  = 100;
LINE_INN  = 110;
SURF_INN  = 120;
SURF_LEG  = 140;
TINY      =  1.0e-9;
HUGE      =  1.0e+9;

//----------------------
//
// Definition of points
//
//----------------------

q = (radius^2 + (width/2.0)^2)^0.5;

// Center point
Point(PNT_CENT) = { 0, 0, z_min};

// At maximum x
Point(PNT_DOM+1) = { q, -width/2.0, z_min};
Point(PNT_DOM+2) = { q,  width/2.0, z_min};

// At maximum y
Point(PNT_DOM+3) = { width/2.0, q, z_min};
Point(PNT_DOM+4) = {-width/2.0, q, z_min};

// At minimum x
Point(PNT_DOM+5) = {-q,  width/2.0, z_min};
Point(PNT_DOM+6) = {-q, -width/2.0, z_min};

// At minimum y
Point(PNT_DOM+7) = {-width/2.0, -q, z_min};
Point(PNT_DOM+8) = { width/2.0, -q, z_min};

// Leg at maximum x
Point(PNT_LEGS+1) = { length, -width/2.0, z_min};
Point(PNT_LEGS+2) = { length,  width/2.0, z_min};

// Leg at maximum y
Point(PNT_LEGS+3) = { width/2.0, length, z_min};
Point(PNT_LEGS+4) = {-width/2.0, length, z_min};

// Leg at minimum x
Point(PNT_LEGS+5) = {-length,  width/2.0, z_min};
Point(PNT_LEGS+6) = {-length, -width/2.0, z_min};

// Leg at minimum y
Point(PNT_LEGS+7) = {-width/2.0, -length, z_min};
Point(PNT_LEGS+8) = { width/2.0, -length, z_min};

width  = width  * 0.8;
q = (inner^2 + (width/2.0)^2)^0.5;

// At maximum x
Point(PNT_INN+1) = { q, -width/2.0, z_min};
Point(PNT_INN+2) = { q,  width/2.0, z_min};

// At maximum y
Point(PNT_INN+3) = { width/2.0, q, z_min};
Point(PNT_INN+4) = {-width/2.0, q, z_min};

// At minimum x
Point(PNT_INN+5) = {-q,  width/2.0, z_min};
Point(PNT_INN+6) = {-q, -width/2.0, z_min};

// At minimum y
Point(PNT_INN+7) = {-width/2.0, -q, z_min};
Point(PNT_INN+8) = { width/2.0, -q, z_min};

//------------------------------------------------------------------------------
//
// Lines
//
//------------------------------------------------------------------------------

// For the domains
For n In{ 1 : 7 }
  Circle(CIRC_OUT + n) = {PNT_DOM + n, PNT_CENT, PNT_DOM + n + 1};
EndFor
Circle(CIRC_OUT + 8) = {PNT_DOM + 8, PNT_CENT, PNT_DOM + 1};

For n In{ 1 : 7 }
  Circle(CIRC_INN + n) = {PNT_INN + n, PNT_CENT, PNT_INN + n + 1};
EndFor
Circle(CIRC_INN + 8) = {PNT_INN + 8, PNT_CENT, PNT_INN + 1};

// Legs at maximum and minimum x
If(VARIANT == 12)
  Line(LINE_LEG+ 1) = {PNT_LEGS+1, PNT_LEGS+2};
  Line(LINE_LEG+ 2) = {PNT_DOM +1, PNT_LEGS+1};
  Line(LINE_LEG+ 3) = {PNT_DOM +2, PNT_LEGS+2};
  Line(LINE_LEG+ 4) = {PNT_LEGS+5, PNT_LEGS+6};
  Line(LINE_LEG+ 5) = {PNT_DOM +5, PNT_LEGS+5};
  Line(LINE_LEG+ 6) = {PNT_DOM +6, PNT_LEGS+6};
// Legs at maximum and minimum y
ElseIf(VARIANT == 22)
  Line(LINE_LEG+ 1) = {PNT_LEGS+3, PNT_LEGS+4};
  Line(LINE_LEG+ 2) = {PNT_DOM +3, PNT_LEGS+3};
  Line(LINE_LEG+ 3) = {PNT_DOM +4, PNT_LEGS+4};
  Line(LINE_LEG+ 4) = {PNT_LEGS+7, PNT_LEGS+8};
  Line(LINE_LEG+ 5) = {PNT_DOM +7, PNT_LEGS+7};
  Line(LINE_LEG+ 6) = {PNT_DOM +8, PNT_LEGS+8};
ElseIf(VARIANT == 21)
  Line(LINE_LEG+ 1) = {PNT_LEGS+3, PNT_LEGS+4};
  Line(LINE_LEG+ 2) = {PNT_DOM +3, PNT_LEGS+3};
  Line(LINE_LEG+ 3) = {PNT_DOM +4, PNT_LEGS+4};
EndIf

// Lines connecting inner and outer circle
For n In{ 1 : 8 }
  Line(LINE_INN+n) = {PNT_INN+n, PNT_DOM+n};
EndFor

// Set resolutions on the lines connecting inner and outer circle
Transfinite Curve {LINE_INN+1:LINE_INN+8} = n_width Using Progression 1;
For n In{ 1 : 7 : 2 }
  Transfinite Curve {CIRC_OUT+n} = n_width Using Progression 1;
  Transfinite Curve {CIRC_INN+n} = n_width Using Progression 1;
EndFor
For n In{ 2 : 8 : 2 }
  Transfinite Curve {CIRC_OUT+n} = n_circ Using Progression 1;
  Transfinite Curve {CIRC_INN+n} = n_circ Using Progression 1;
EndFor

//------------------------------------------------------------------------------
//
// Surfaces
//
//------------------------------------------------------------------------------

// Create surfaces on the boundaries of the domain
For n In{ 1 : 8 }
  If(n < 8)
    Curve Loop(n) = {LINE_INN+n, CIRC_OUT+n, -(LINE_INN+n+1), -(CIRC_INN+n)};
    Plane Surface(SURF_INN+n) = {n};
  Else
    Curve Loop(n) = {LINE_INN+n, CIRC_OUT+n, -(LINE_INN+1), -(CIRC_INN+n)};
    Plane Surface(SURF_INN+n) = {n};
  EndIf
EndFor
For n In{ 1 : 8 }
  If(n < 8)
    Transfinite Surface{SURF_INN+n} = {PNT_INN+n,   PNT_DOM+n,
                                       PNT_DOM+n+1, PNT_INN+n+1};
    Recombine Surface{SURF_INN+n};
  Else
    Transfinite Surface{SURF_INN+n} = {PNT_INN+n, PNT_DOM+n,
                                       PNT_DOM+1, PNT_INN+1};
    Recombine Surface{SURF_INN+n};
  EndIf
EndFor

Curve Loop(9) = {CIRC_INN+1:CIRC_INN+8};
Plane Surface(SURF_INN+9) = {9};
Recombine Surface{SURF_INN+9};

// Surface for legs in x VARIANT
If(VARIANT == 12)
  Curve Loop(10) = {CIRC_OUT+1, LINE_LEG+3, -(LINE_LEG+1), -(LINE_LEG+2)};
  Plane Surface(SURF_LEG+1) = {10};
  Recombine Surface{SURF_LEG+1};
  Curve Loop(11) = {CIRC_OUT+5, LINE_LEG+6, -(LINE_LEG+4), -(LINE_LEG+5)};
  Plane Surface(SURF_LEG+2) = {11};
  Recombine Surface{SURF_LEG+2};
// Surface for legs in y VARIANT
ElseIf(VARIANT == 22)
  Curve Loop(10) = {CIRC_OUT+3, LINE_LEG+3, -(LINE_LEG+1), -(LINE_LEG+2)};
  Plane Surface(SURF_LEG+1) = {10};
  Recombine Surface{SURF_LEG+1};
  Curve Loop(11) = {CIRC_OUT+7, LINE_LEG+6, -(LINE_LEG+4), -(LINE_LEG+5)};
  Plane Surface(SURF_LEG+2) = {11};
  Recombine Surface{SURF_LEG+2};
ElseIf(VARIANT == 21)
  Curve Loop(10) = {CIRC_OUT+3, LINE_LEG+3, -(LINE_LEG+1), -(LINE_LEG+2)};
  Plane Surface(SURF_LEG+1) = {10};
  Recombine Surface{SURF_LEG+1};
EndIf

//------------------------------------------------------------------------------
//
// Mesh
//
//------------------------------------------------------------------------------
Field[1] = MathEval;
Field[1].F = Sprintf("(%5.2g)", delta);
Background Field = 1;

// This experimental algorithm could give better results
Mesh.Algorithm = 8;

// Uncomment the following line to try the full-quad algorithm:
Mesh.RecombinationAlgorithm = 2; // or 3

//------------------------------------------------------------------------------
//
// Volume
//
//------------------------------------------------------------------------------
If(VARIANT == 12 || VARIANT == 22)
  Extrude {0, 0, (z_max-z_min)} {
    Surface{SURF_INN+1:SURF_INN+9, SURF_LEG+1:SURF_LEG+2};
    Layers{n_lay};
    Recombine;
  }
ElseIf(VARIANT == 30)
  Extrude {0, 0, (z_max-z_min)} {
    Surface{SURF_INN+1:SURF_INN+9};
    Layers{n_lay};
    Recombine;
  }
ElseIf(VARIANT == 21)
  Extrude {0, 0, (z_max-z_min)} {
    Surface{SURF_INN+1:SURF_INN+9, SURF_LEG+1};
    Layers{n_lay};
    Recombine;
  }
EndIf

//------------------------------------------------------------------------------
//
// Boundary conditions
//
//------------------------------------------------------------------------------
Printf("Radius is %g", radius);
Physical Surface("top")
  = {Surface In BoundingBox{-(radius+0.2), -(radius+0.2), z_max-TINY,
                            +(radius+0.2), +(radius+0.2), z_max+TINY}};
Physical Surface("bottom")
  = {Surface In BoundingBox{-(radius+0.2), -(radius+0.2), z_min-TINY,
                            +(radius+0.2), +(radius+0.2), z_min+TINY}};

If(VARIANT == 12)
  Physical Surface("cylinder") = {177, 199, 309, 221, 265,
                                  287, 373, 381, 395, 403,
                                  404, 142, 382, 141};
  Physical Surface("west") = {399};
  Physical Surface("east") = {377};
ElseIf(VARIANT == 22)
  Physical Surface("cylinder") = {155, 177, 221, 243, 265,
                                  309, 373, 381, 395, 403,
                                  142, 404, 382, 141};
  Physical Surface("south") = {399};
  Physical Surface("north") = {377};
ElseIf(VARIANT == 30)
  Physical Surface("cylinder") = {142, 164, 186, 208,
                                  230, 252, 274, 296};
EndIf

If(VARIANT == 21)
  Physical Surface("cylinder") = {141, 381, 286, 308, 264, 154,
                                  242, 220, 176, 380, 372};
  Physical Surface("north")    = {376};
EndIf

//------------------------------------------------------------------------------
//
// Volume physical entity
//
//------------------------------------------------------------------------------
If(VARIANT == 12 || VARIANT == 22)
  Physical Volume("interior") = {1:11};
ElseIf(VARIANT == 21)
  Physical Volume("interior") = {1:10};
ElseIf(VARIANT == 30)
  Physical Volume("interior") = {1:9};
EndIf
