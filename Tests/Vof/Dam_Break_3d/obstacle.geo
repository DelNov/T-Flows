// Define domain and step dimensions
L_DOM_X  = 3.22;
L_DOM_Y  = 1.0;
L_BUF_Y  = 0.06;  // length of buffers in y direction of the domain
L_DOM_Z  = 1.5;
L_OBST_X = 0.161;
L_OBST_Y = 0.403;
L_OBST_Z = 0.161;  // height of the obstacle
S_OBST_X = 0.6635;  // shift for the obstacle (distance from the left wall)

// Define domain and step resolution
N_OBST_X = 17;
N_OBST_Y = 41;
N_OBST_Z = 17;

N_DOM_X_LEFT  = 51;  // from the obstacle towards left (towards xmin)
N_DOM_X_RIGHT = 121;  // from the obstacle towards left (towards xmax)
N_DOM_Y_LEFT  = 41;  // spanwise resolution on the left (at xmin)
N_DOM_Y_RIGHT = 21;  // spanwise resolution on the right (at xmax)
N_BUF_Y       =  4;

N_ABOVE_OBST_Z = 134;

TINY = 1.0e-9;
HUGE = 1.0e+9;

//--------
// Points
//--------

// Points defining the domain
Point( 1) = {0,                         L_BUF_Y,         0};
Point( 2) = {S_OBST_X + L_OBST_X / 2.0, L_BUF_Y,         0};
Point( 3) = {L_DOM_X,                   L_BUF_Y,         0};
Point( 4) = {L_DOM_X,                   L_DOM_Y-L_BUF_Y, 0};
Point( 5) = {S_OBST_X + L_OBST_X / 2.0, L_DOM_Y-L_BUF_Y, 0};
Point( 6) = {0,                         L_DOM_Y-L_BUF_Y, 0};

// Points defining the obstacle
Point( 7) = {S_OBST_X,           L_DOM_Y/2 - L_OBST_Y/2.0, 0};
Point( 8) = {S_OBST_X+L_OBST_X,  L_DOM_Y/2 - L_OBST_Y/2.0, 0};
Point( 9) = {S_OBST_X+L_OBST_X,  L_DOM_Y/2 + L_OBST_Y/2.0, 0};
Point(10) = {S_OBST_X,           L_DOM_Y/2 + L_OBST_Y/2.0, 0};

//-------
// Lines
//-------

// Lines defining the domain
Line(1) = {2, 1};
Line(2) = {5, 6};
Line(3) = {2, 3};
Line(4) = {5, 4};
Line(5) = {1, 6};
Line(6) = {3, 4};

// Lines defining the obstacle
Line( 7) = {7, 8};
Line( 8) = {8, 9};
Line( 9) = {9, 10};
Line(10) = {10, 7};

// Along the length (short side) of the obstacle)
Transfinite Curve {7, 9}  = N_OBST_X Using Progression 1;

// Along the width (longer side) of the obstacle)
Transfinite Curve {8, 10} = N_OBST_Y Using Progression 1;

// On the boundary, from the obstacle center towards the beginning
Transfinite Curve {1, 2}  = N_DOM_X_LEFT Using Progression 1.01;

// On the boundary, from the obstacle towards the right
Transfinite Curve {3, 4}  = N_DOM_X_RIGHT Using Progression 1.01;

// Along the width of the domain, left side
Transfinite Curve {5}  = N_DOM_Y_LEFT Using Progression 1;

// Along the width of the domain, right side
Transfinite Curve {6}  = N_DOM_Y_RIGHT Using Progression 1;

//----------
// Surfaces
//----------

Curve Loop(1) = {1, 5, -2, 4, -6, -3};
Curve Loop(2) = {10, 7, 8, 9};
Plane Surface(1) = {1, 2};
Recombine Surface {1};

// Extrude the lower surface
Extrude {0, 0, L_OBST_Z} {
  Surface {1};
  Layers{N_OBST_Z};
  Recombine;
}

// Along the width (longer side) of the obstacle)
Transfinite Curve {8, 10} = N_OBST_Y Using Progression 1;

// Define additional surface on top of the obstacle
Curve Loop(3) = {20, 21, 18, 19};
Plane Surface(63) = {3};
Transfinite Surface(63);
Recombine Surface {63};

// Extrude the surfaces on top of the obstacle
Extrude {0, 0, L_DOM_Z - L_OBST_Z} {
  Surface {62, 63};
  Layers{N_ABOVE_OBST_Z};
  Recombine;
}

// Extrude buffer layers on bottom and top
Extrude {0, -L_BUF_Y, 0} {
  Surface{78}; Surface{98}; Surface{25}; Surface{45}; Layers{N_BUF_Y};
  Recombine;
}
Extrude {0, L_BUF_Y, 0} {
  Surface{90}; Surface{37}; Surface{33}; Surface{86}; Layers{N_BUF_Y};
  Recombine;
}

//-----------------------------------
// Define volumes as physical spaces
//-----------------------------------
Physical Volume("FLUID") = {1:11};

//----------------------------
// Define boundary conditions
//----------------------------
Physical Surface("wall")
   = {Surface In BoundingBox{-TINY,        -HUGE,        -HUGE,
                             +TINY,        +HUGE,        +HUGE},
      Surface In BoundingBox{L_DOM_X-TINY, -HUGE,        -HUGE,
                             L_DOM_X+TINY, +HUGE,        +HUGE},
      Surface In BoundingBox{-HUGE,        -TINY,        -HUGE,
                             +HUGE,        +TINY,        +HUGE},
      Surface In BoundingBox{-HUGE,        L_DOM_Y-TINY, -HUGE,
                             +HUGE,        L_DOM_Y+TINY, +HUGE},
      Surface In BoundingBox{-HUGE,        -HUGE,        -TINY,
                             +HUGE,        +HUGE,        +TINY}};
Physical Surface("top")
   = {Surface In BoundingBox{-HUGE,        -HUGE,        L_DOM_Z-TINY,
                             +HUGE,        +HUGE,        L_DOM_Z+TINY}};
Physical Surface("step") = {49, 57, 61, 53, 63};
