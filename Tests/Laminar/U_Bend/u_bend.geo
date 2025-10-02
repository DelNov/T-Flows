/*******************************************************************************
*                                                                              *
*  U-bend and precursor domain                                                 *
*                                                                              *
*  This script can create either U-bend or its precursor, or both,             *
*  depending if variables PRECURSOR and U_BEND are set to 1 here.              *
*                                                                              *
*  If PRECURSOR is set to 1 and U_BEND to 0, you should run it with:           *
*  gmsh -3 u_bend.geo -o precursor.msh                                         *
*                                                                              *
*  Otherwise, if U_BEND is set to 1 and PRECURSOR to 0, you run it with:       *
*  gmsh -3 u_bend.geo -o u_bend.msh                                            *
*                                                                              *
*  This version adjusts progression parameters towards the bend.               *
*  It should give higher quality meshes, but it might prove to be unstable.    *                                                                             *
*******************************************************************************/

//-------------------------------
//
// Handling command line options
//
//-------------------------------
If(!Exists(U_BEND) || !Exists(PRECURSOR))
  Printf("Variables U_BEND and/or PRECURSOR are not defined from command line.");
  Printf("The proper invokation of Gmsh for this scrip is, for example:");
  Printf("");
  Printf("gmsh -3 u_bend.geo -setnumber U_BEND 1 -setnumber PRECURSOR 0 -o u_bend.msh");
  Printf("");
  Printf("to mesh only the U-bend domain, or:");
  Printf("");
  Printf("gmsh -3 u_bend.geo -setnumber U_BEND 0 -setnumber PRECURSOR 1 -o precursor.msh");
  Printf("");
  Printf("to mesh only the precursor domain.");
  Abort;
EndIf

// Parameters related to geometrical quantities
L_UPPER     = 0.38;     // length of the upper leg
L_LOWER     = 0.4572;   // length of the lower leg
L_PRECURSOR = 0.2;      // length of the precursor domain
D_PRECURSOR = 0.02;     // distance between precursor and U-bend domain
W_BOTH      = 0.02;     // width of both domains
R_INSIDE    = 0.01905;  // inside radius of the U-bend
R_OUTSIDE   = 0.05715;  // outside radius of the U-bend
R_DELTA     = R_OUTSIDE - R_INSIDE;

// Parameters related to grid resolution
N_UPPER     =  61;  // number of nodes in the upper leg
N_LOWER     = 121;  // number of nodes in the lower leg
N_ARC       = 100;  // number of nodes in the arc (if odd, gmsh complains)
N_ACROSS    =  51;  // number of nodes across the height of both domains
N_PRECURSOR =  13;  // number of nodes in the streamwise direction of precursor
N_WIDTH     =   4;  // number of nodes in the spanwise direction

// Parameter related to grid stretching/clustering
B_ACROSS    =   0.2;   // "bump" across the height of both domains

// Parameters related to grid stretching/clustering
p_upper_leg =   0.99;  // progression for cells in the upper leg
p_lower_leg =   0.99;  // progression for cells in the lower leg

//-----------------------------
//
// Compute progression factors
//
//-----------------------------

//-------------------
// For the upper leg
//-------------------
p = p_upper_leg;  // set p to initial progression
delta_min = (R_INSIDE+R_OUTSIDE)*0.5/(N_ARC-1)*Pi;
Printf("Initial delta_min based on average u-bend radius is %g", delta_min);

d = (R_INSIDE+R_OUTSIDE)*0.5/(N_ARC-1)*Pi;  // initial increment
l = delta_min * (1 - (1.0/p)^(N_UPPER)) / (1 - (1.0/p));
Printf("Length with progression %g is %g", p, l);
For iter In{1:30}
  If(l > L_UPPER)
    p = p + d;
    l = delta_min * (1 - (1.0/p)^(N_UPPER)) / (1 - (1.0/p));
    If(iter % 10 == 0)
      Printf("  Iteration %g; Length with progression %g is %g", iter, p, l);
    EndIf
    If(l < L_UPPER) d = d * 0.5; EndIf
  EndIf
  If(l < L_UPPER)
    p = p - d;
    l = delta_min * (1 - (1.0/p)^(N_UPPER)) / (1 - (1.0/p));
    If(iter % 10 == 0)
      Printf("  Iteration %g; Length with progression %g is %g", iter, p, l);
    EndIf
    If(l > L_UPPER) d = d * 0.5; EndIf
  EndIf
EndFor

delta_min_fin = L_UPPER * (1 - (1.0 / p)) / (1 - (1.0 / p)^(N_UPPER));
Printf("Final delta_min based on final p %g is %g", p, delta_min_fin);
p_upper_leg = p;

//-------------------
// For the lower leg
//-------------------
p = p_lower_leg;  // set p to initial progression
delta_min = (R_INSIDE+R_OUTSIDE)*0.5/(N_ARC-1)*Pi;
Printf("Initial delta_min based on average u-bend radius is %g", delta_min);

d = (R_INSIDE+R_OUTSIDE)*0.5/(N_ARC-1)*Pi;  // initial increment
l = delta_min * (1 - (1.0/p)^(N_LOWER)) / (1 - (1.0/p));
Printf("Length with progression %g is %g", p, l);
For iter In{1:30}
  If(l > L_LOWER)
    p = p + d;
    l = delta_min * (1 - (1.0/p)^(N_LOWER)) / (1 - (1.0/p));
    If(iter % 10 == 0)
      Printf("  Iteration %g; Length with progression %g is %g", iter, p, l);
    EndIf
    If(l < L_LOWER) d = d * 0.5; EndIf
  EndIf
  If(l < L_LOWER)
    p = p - d;
    l = delta_min * (1 - (1.0/p)^(N_LOWER)) / (1 - (1.0/p));
    If(iter % 10 == 0)
      Printf("  Iteration %g; Length with progression %g is %g", iter, p, l);
    EndIf
    If(l > L_LOWER) d = d * 0.5; EndIf
  EndIf
EndFor

delta_min_fin = L_LOWER * (1 - (1.0 / p)) / (1 - (1.0 / p)^(N_LOWER));
Printf("Final delta_min based on final p %g is %g", p, delta_min_fin);
p_lower_leg = p;

//--------
//
// Points
//
//--------

// Central point
Point(1) = {0, 0, -W_BOTH/2};

// Points defining inner wall
Point(2) = {-L_UPPER,  R_INSIDE, -W_BOTH/2};
Point(3) = {0,         R_INSIDE, -W_BOTH/2};
Point(4) = {0,        -R_INSIDE, -W_BOTH/2};
Point(5) = {-L_LOWER, -R_INSIDE, -W_BOTH/2};

// Points defining outer wall
Point(6) = {-L_UPPER,  R_OUTSIDE, -W_BOTH/2};
Point(7) = {0,         R_OUTSIDE, -W_BOTH/2};
Point(8) = {0,        -R_OUTSIDE, -W_BOTH/2};
Point(9) = {-L_LOWER, -R_OUTSIDE, -W_BOTH/2};

// Points for precursor domain
Point(10) = {-L_UPPER -D_PRECURSOR,                R_INSIDE,  -W_BOTH/2};
Point(11) = {-L_UPPER -D_PRECURSOR,                R_OUTSIDE, -W_BOTH/2};
Point(12) = {-L_UPPER -D_PRECURSOR - L_PRECURSOR,  R_INSIDE,  -W_BOTH/2};
Point(13) = {-L_UPPER -D_PRECURSOR - L_PRECURSOR,  R_OUTSIDE, -W_BOTH/2};

//----------------
//
// Lines and arcs
//
//----------------

// Lines on the top
Line(1) = {2, 3};
Line(2) = {6, 7};
Transfinite Curve {1, 2} = N_UPPER Using Progression p_upper_leg;

// Lines on the bottom
Line(3) = {9, 8};
Line(4) = {5, 4};
Transfinite Curve {3, 4} = N_LOWER Using Progression p_lower_leg;

// Arcs
Circle(5) = {8, 1, 7};
Circle(6) = {4, 1, 3};
Transfinite Curve {6, 5} = N_ARC Using Progression 1;

// Across the U
Line(7) = {5, 9};
Line(8) = {4, 8};
Line(9) = {3, 7};
Line(10) = {2, 6};
Transfinite Curve {7, 8, 9, 10} = N_ACROSS Using Bump B_ACROSS;

// Lines along the precursor
Line(11) = {12, 10};
Line(12) = {13, 11};
Transfinite Curve {11, 12} = N_PRECURSOR Using Progression 1;

// Across the pre-cursor
Line(13) = {10, 11};
Line(14) = {12, 13};
Transfinite Curve {14, 13} = N_ACROSS Using Bump B_ACROSS;

//----------
//
// Surfaces
//
//----------

// U-bend
Curve Loop(1) = {3, -8, -4, 7};
Plane Surface(1) = {1};
Curve Loop(2) = {5, -9, -6, 8};
Plane Surface(2) = {2};
Curve Loop(3) = {1, 9, -2, -10};
Plane Surface(3) = {3};
Transfinite Surface {1} = {9, 8, 4, 5};
Transfinite Surface {2} = {8, 7, 3, 4};
Transfinite Surface {3} = {3, 7, 6, 2};
Recombine Surface{1, 2, 3};

// Precursor
Curve Loop(4) = {11, 13, -12, -14};
Plane Surface(4) = {4};
Transfinite Surface {4} = {12, 10, 11, 13};
Recombine Surface{4};

//--------
//
// Volume
//
//--------
Extrude {0, 0, W_BOTH} {
  Surface{4}; Surface{3}; Surface{2}; Surface{1}; Layers {N_WIDTH-1}; Recombine;
}

//------------------------------------
//
// Set boundary and volume conditions
//
//------------------------------------

If(PRECURSOR == 0)
  Recursive Delete {
    Volume{1};
  }
  Physical Surface("PERIODIC_Z", 103) = {58, 3, 80, 2, 102, 1};
  Physical Surface("U_INLET", 104) = {57};
  Physical Surface("U_OUTLET", 105) = {101};
  Physical Surface("WALLS_OUTER", 106) = {53, 67, 89};
  Physical Surface("WALLS_INNER", 107) = {45, 75, 97};
  Physical Volume("FLUID", 108) = {2, 3, 4};
EndIf

If(U_BEND == 0)
  Recursive Delete {
    Volume{4}; Volume{3}; Volume{2};
  }
  Physical Surface("PERIODIC_X", 102) = {27, 35};
  Physical Surface("PERIODIC_Z", 103) = {36, 4};
  Physical Surface("WALLS", 104) = {23, 31};
  Physical Volume("FLUID", 105) = {1};
EndIf

Printf("Cross-sectional area is  %g", R_DELTA * W_BOTH);
Printf("Total s outer coordinate %g", L_UPPER + L_LOWER + R_OUTSIDE * Pi);
Printf("Total s inner coordinate %g", L_UPPER + L_LOWER + R_INSIDE  * Pi);
