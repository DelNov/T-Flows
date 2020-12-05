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
H   =   D * PI;  // cylinder height
N_C =  80;       // number of cells in circumferential direction
N_H =  80;       // number of cells in z direction (height)

// Computed parameters
delta = D * PI / N_C;

//--------
// Points
//--------
Point(1) = {0,   -D/2, 0, delta};
Point(2) = {D/2,  0,   0, delta};
Point(3) = {0.0,  D/2, 0, delta};
Point(4) = {-D/2, 0.0, 0, delta};
Point(5) = {0.0,  0.0, 0, delta};

//------------
// Connectors
//------------
Circle(1) = {2, 5, 3};
Circle(2) = {3, 5, 4};
Circle(3) = {4, 5, 1};
Circle(4) = {1, 5, 2};

//------------
// Surface(s)
//------------
Curve Loop(1) = {3, 4, 1, 2};
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
Physical Surface("wall") = {26, 1, 17, 13, 21, 25};
Physical Volume("interior") = {1};
