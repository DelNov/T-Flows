//------------------------------------------------------------------------------
// Parameters defining domain extents, buildings and mesh resolution
// (This should be fiddled with, but with care)
//------------------------------------------------------------------------------
lx     =   1.0;
ly     =   0.2;
lz     =   1.0;
n_half =  80;       // number of nodes along half of the domain
                    // values I would like to use here are 20, 40, 80, 160
n_lay  =   3;       // number of cells along half of the domain
bump   =   1.0;     // progression for cell clustering
angle  =  45.0;     // cavity skewness - it should be the same as the angle
                    // defined in extract_*.py Paraview macros

//------------------------------------------------------------------------------
// Parameters for the problem definition algorithms
// (These are better kept untouched)
//------------------------------------------------------------------------------
TINY = 1.0e-9;
HUGE = 1.0e+9;

//------------------------------------------------------------------------------
// Points
//------------------------------------------------------------------------------
sin = Sin(angle * Pi / 180.0);
cos = Cos(angle * Pi / 180.0);

Point(1) = {0,               0,   0,          1.0};
Point(2) = {lx,              0,   0,          1.0};
Point(3) = {lx + lx * cos,   0,   lz * sin,   1.0};
Point(4) = {     lx * cos,   0,   lz * sin,   1.0};

//------------------------------------------------------------------------------
// Lines
//------------------------------------------------------------------------------
Line (1) = {1, 2};  Line (2) = {2, 3};  Line (3) = {3, 4};  Line (4) = {4, 1};

// Define all lines as transfinite
Transfinite Line "*" = n_half Using Bump bump;

//------------------------------------------------------------------------------
// Surfaces
//------------------------------------------------------------------------------
Line Loop(1) = {  1,  2,   3,  4};  Plane Surface(1) = {1};

// Define all surfaces as transfinite
Recombine Surface    "*";
Transfinite Surface  "*";

//------------------------------------------------------------------------------
// Volume
//------------------------------------------------------------------------------
Extrude {0, ly, 0} {
  Surface {1};
  Layers{n_lay};
  Recombine;
}

//------------------------------------------------------------------------------
// Boundary conditions
//------------------------------------------------------------------------------
Physical Surface("periodic_y", 27) = {1, 26};
Physical Surface("top_wall", 28) = {21};
Physical Surface("bottom_wall", 29) = {13};
Physical Surface("west_wall", 30) = {25};
Physical Surface("east_wall", 31) = {17};

//------------------------------------------------------------------------------
// Pysical volume
//------------------------------------------------------------------------------
Physical Volume("interior") = {1};
