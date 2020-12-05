//------------------------------------------------
// Follow these steps in GMSH to create the mesh:
//
// Mesh -> 2D
// Mesh -> Recombine 2D
// Mesh -> 3D
//------------------------------------------------

// Handy constants
PI = Acos(-1);

//--------------------------------------
// Set parameters defining the geometry
//--------------------------------------
D   =   1.0;     // cylinder diameter
P   =   1.25;    // pitch
H   =   D * PI;  // cylinder height
N_C = 120;       // number of cells in circumferential direction
N_H = 120;       // number of cells in z direction (height)

// Computed parameters
delta = D * PI / N_C;

//--------
// Points
//--------
Point(1) = {D/2,   0,     0, delta};
Point(2) = {P-D/2, 0,     0, delta};
Point(3) = {P,     D/2,   0, delta};
Point(4) = {P,     P-D/2, 0, delta};
Point(5) = {P-D/2, P,     0, delta};
Point(6) = {D/2,   P,     0, delta};
Point(7) = {0,     P-D/2, 0, delta};
Point(8) = {0,     D/2,   0, delta};

Point( 9) = {0,     0,     0, delta};
Point(10) = {P,     0,     0, delta};
Point(11) = {P,     P,     0, delta};
Point(12) = {0,     P,     0, delta};

//------------
// Connectors
//------------
Circle(1) = { 8,  9,  1};
Circle(2) = { 2, 10,  3};
Circle(3) = { 4, 11,  5};
Circle(4) = { 6, 12,  7};
Line  (5) = {1, 2};
Line  (6) = {3, 4};
Line  (7) = {5, 6};
Line  (8) = {7, 8};

//------------
// Surface(s)
//------------
Curve Loop(1) = {1, 5, 2, 6, 3, 7, 4, 8};
Plane Surface(1) = {1};

//-----------
// Volume(s)
//-----------
Extrude {0, 0, H} {
  Surface{1}; Layers{N_H}; Recombine;
}

//-------------------
// Physical entities
//-------------------
Physical Surface("wall") = {29, 33, 37, 25, 21, 49, 41, 45, 50, 1};
Physical Volume("interior") = {1};
