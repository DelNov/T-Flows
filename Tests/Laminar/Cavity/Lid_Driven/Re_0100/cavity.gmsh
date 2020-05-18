//------------------------------------------------------------------------------
// Parameters defining domain extents, buildings and mesh resolution
// (This should be fiddled with, but with care)
//------------------------------------------------------------------------------
lx     =  1.0;
ly     =  0.2;
lz     =  1.0;
n_half = 30;       // number of cells along half of the domain
n_lay  =  5;       // number of cells along half of the domain
prog   =  1.05;    // progression for cell clustering

//------------------------------------------------------------------------------
// Parameters for the problem definition algorithms
// (These are better kept untouched)
//------------------------------------------------------------------------------
TINY = 1.0e-9;
HUGE = 1.0e+9;

//------------------------------------------------------------------------------
//
//     4-----6-----7-----5----3
//     |           |          |
//     |           |          |
//     7          11          4
//     |           |          |
//     |           |          |
//     8----12-----9----10----6
//     |           |          |
//     |           |          |
//     8           9          3
//     |           |          |
//     |           |          |
//     1-----1-----5-----2----2
//
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Points
//------------------------------------------------------------------------------
Point(1) = {0,     0,     0,     1.0};
Point(2) = {lx,    0,     0,     1.0};
Point(3) = {lx,    0,     lz,    1.0};
Point(4) = {0,     0,     lz,    1.0};
Point(5) = {lx/2,  0,     0,     1.0};
Point(6) = {lx,    0,     lz/2,  1.0};
Point(7) = {lx/2,  0,     lz,    1.0};
Point(8) = {0,     0,     lz/2,  1.0};
Point(9) = {lx/2,  0,     lz/2,  1.0};

//------------------------------------------------------------------------------
// Lines
//------------------------------------------------------------------------------
Line (1) = {1, 5};  Line (2) = {2, 5};  Line (3) = {2, 6};  Line (4) = {3, 6};
Line (5) = {3, 7};  Line (6) = {4, 7};  Line (7) = {4, 8};  Line (8) = {1, 8};
Line (9) = {5, 9};  Line(10) = {6, 9};  Line(11) = {7, 9};  Line(12) = {8, 9};

// Define all lines as transfinite
Transfinite Line "*" = n_half+1 Using Progression prog;

//------------------------------------------------------------------------------
// Surfaces
//------------------------------------------------------------------------------
Line Loop(1) = {  1,  9, -12, -8};  Plane Surface(1) = {1};
Line Loop(2) = { -2,  3,  10, -9};  Plane Surface(2) = {2};
Line Loop(3) = {-10, -4,   5, 11};  Plane Surface(3) = {3};
Line Loop(4) = { 12,-11,  -6,  7};  Plane Surface(4) = {4};

// Define all surfaces as transfinite
Recombine Surface    "*";
Transfinite Surface  "*";

//------------------------------------------------------------------------------
// Volume
//------------------------------------------------------------------------------
Extrude {0, ly, 0} {
  Surface {1, 2, 3, 4};
  Layers{n_lay};
  Recombine;
}

//------------------------------------------------------------------------------
// Boundary conditions
//------------------------------------------------------------------------------
Physical Surface("moving_wall")
  = {Surface In BoundingBox{-HUGE, -HUGE, lz-TINY,
                            +HUGE, +HUGE, lz+TINY}};

Physical Surface("static_wall")
  = {Surface In BoundingBox{  -TINY, -HUGE, -HUGE,
                              +TINY, +HUGE, +HUGE},
     Surface In BoundingBox{lx-TINY, -HUGE, -HUGE,
                            lx+TINY, +HUGE, +HUGE},
     Surface In BoundingBox{  -HUGE, -HUGE, -TINY,
                              +HUGE, +HUGE, +TINY}};

Physical Surface("periodic_y")
  = {Surface In BoundingBox{-HUGE, ly-TINY, -HUGE,
                            +HUGE, ly+TINY, +HUGE},
     Surface In BoundingBox{-HUGE,   -TINY, -HUGE,
                            +HUGE,   +TINY, +HUGE}};

//------------------------------------------------------------------------------
// Pysical volume
//------------------------------------------------------------------------------
Physical Volume("interior") = {1:4};
