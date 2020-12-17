// Gmsh project created on Wed Oct  9 18:07:24 2019

//------------------
// Enable mesh copy
//------------------
Geometry.CopyMeshingMethod = 1;

//--------------------
//
// Set some constants
//
//--------------------

PI   = 3.14159265359;
PI_4 = PI/4.0;

R  =  1;     // outer (true) radius
L  =  6;     // length of one leg
r  =  0.7;   // inner radius
d  =  0.1;   // mesh resolution
NC = 11;     // number of nodes in the core
NB =  7;     // number of nodes in the bundary layer
NL = 41;     // number of nodes longitudinally

//---------------
//
// Define points
//
//---------------

// Center point at z = 0
Point(1) = {0, 0, 0, d};

// Points defining outer surface at z = 0
Point(2) = {R,           0,           0,             d};
Point(3) = {R*Cos(PI_4), R*Cos(PI_4), R*Cos(PI_4),   d};
Point(4) = {0,           R,           R,             d};

// Points defining inner surface at z = 0
Point(5) = {r,           0,           0,             d};
Point(6) = {r*Cos(PI_4), r*Cos(PI_4), r*Cos(PI_4),   d};
Point(7) = {0,           r,           r,             d};

// Center point at z = L
Point(8) = {0, 0, L, d};

// Points defining outer surface at z = L
Point( 9) = {R,           0,           L,   d};
Point(10) = {R*Cos(PI_4), R*Cos(PI_4), L,   d};
Point(11) = {0,           R,           L,   d};

// Points defining inner surface at z = 0
Point(12) = {r,           0,           L,   d};
Point(13) = {r*Cos(PI_4), r*Cos(PI_4), L,   d};
Point(14) = {0,           r,           L,   d};

//--------------------
//
// Define connections
//
//--------------------

// Connections at z = 0
Ellipse(1) = {2, 1, 4, 3};
Ellipse(2) = {3, 1, 4, 4};

Line(3) = {1, 5};
Line(4) = {5, 2};
Line(5) = {1, 7};
Line(6) = {7, 4};
Line(7) = {5, 6};
Line(8) = {7, 6};
Line(9) = {6, 3};

// Connections at z = L
Ellipse(10) = {9, 8, 11, 10};
Ellipse(11) = {10, 8, 11, 11};

Line(12) = { 8, 12};
Line(13) = {12,  9};
Line(14) = { 8, 14};
Line(15) = {14, 11};
Line(16) = {12, 13};
Line(17) = {14, 13};
Line(18) = {13, 10};

// Connections in between
Line(19) = { 1,  8};
Line(20) = { 2,  9};
Line(21) = { 3, 10};
Line(22) = { 4, 11};
Line(23) = { 5, 12};
Line(24) = { 6, 13};
Line(25) = { 7, 14};

//-----------------
//
// Create surfaces
//
//-----------------

// Resolution in the core
Transfinite Curve {12, 16, 17, 14, 10, 11, 3, 7, 8, 5, 1, 2} = NC Using Progression 1;

// Resolution in the boundary layers
Transfinite Curve {13, 18, 15, 4, 9, 6} = NB Using Progression 1;

// Longitudinal resolution
Transfinite Curve {19, 23, 20, 24, 21, 25, 22} = NL Using Progression 1;

// Create faces at z = 0
Curve Loop(1) = { 3,  7, -8, -5};
Curve Loop(2) = { 4,  1, -9, -7};
Curve Loop(3) = { 8,  9,  2, -6};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};

// Create faces at z = L
Curve Loop(4) = { 12,  16, -17, -14};
Curve Loop(5) = { 16,  18, -10, -13};
Curve Loop(6) = { 17,  18,  11, -15};

Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};

// Create faces in between
Curve Loop( 7) = {12, -23, -3, 19};
Curve Loop( 8) = {16, -24, -7, 23};
Curve Loop( 9) = {17, -24, -8, 25};
Curve Loop(10) = {14, -25, -5, 19};
Curve Loop(11) = {13, -20, -4, 23};
Curve Loop(12) = {18, -21, -9, 24};
Curve Loop(13) = {15, -22, -6, 25};
Curve Loop(14) = {10, -21, -1, 20};
Curve Loop(15) = {11, -22, -2, 21};

Plane Surface( 7) = { 7};
Plane Surface( 8) = { 8};
Plane Surface( 9) = { 9};
Plane Surface(10) = {10};

Plane Surface(11) = {11};
Plane Surface(12) = {12};
Plane Surface(13) = {13};

// These are part of the cylinder
Surface(14) = {14};
Surface(15) = {15};

//----------------
//
// Define volumes
//
//----------------

Surface Loop(1) = {4, 7, 1, 10, 8, 9};
Surface Loop(2) = {5, 14, 2, 11, 12, 8};
Surface Loop(3) = {6, 15, 13, 3, 12, 9};
Volume(1) = {1};
Volume(2) = {2};
Volume(3) = {3};

//-----------------------
//
// Define surface meshes
//
//-----------------------

Transfinite Surface {1, 2, 3};
Recombine Surface   {1, 2, 3};
Transfinite Surface {4, 5, 6};
Recombine Surface   {4, 5, 6};
Transfinite Surface {7, 8, 9, 10};
Recombine Surface   {7, 8, 9, 10};
Transfinite Surface {11, 12, 13};
Recombine Surface   {11, 12, 13};
Transfinite Surface {14, 15};
Recombine Surface   {14, 15};

//----------------------
//
// Define volume meshes
//
//----------------------

Transfinite Volume{1} = {1, 5, 6, 7, 8, 12, 13, 14};
Transfinite Volume{2} = {5, 2, 3, 6, 12, 9, 10, 13};
Transfinite Volume{3} = {7, 6, 3, 4, 14, 13, 10, 11};

//----------------------------
// Copy the basic unit around
//----------------------------

Dilate {{0, 0, 0}, {-1, 1, 1}} {
  Duplicata { Volume{1}; Volume{2}; Volume{3}; }
}

Symmetry {0, -1, 1, 0} {
  Duplicata { Volume  {3};  Volume {88};  Volume {57};
              Volume {26};  Volume  {1};  Volume  {2};}
}

Dilate {{0, 0, 0}, {1, 1, -1}} {
  Duplicata { Volume  {2};  Volume  {1};  Volume {26};
              Volume {57};  Volume {88};  Volume  {3};
              Volume{105};  Volume{136};  Volume{260};
              Volume{229};  Volume{198};  Volume{167}; }
}

Dilate {{0, 0, 0}, {1, -1, 1}} {
  Duplicata { Volume  {2};  Volume  {1};  Volume {26};
              Volume {57};  Volume  {3};  Volume {88};
              Volume{105};  Volume{136};  Volume{260};
              Volume{229};  Volume{198};  Volume{167};
              Volume{525};  Volume{556};  Volume{463};
              Volume{587};  Volume{494};  Volume{618};
              Volume{277};  Volume{308};  Volume{432};
              Volume{339};  Volume{401};  Volume{370}; }
}

Physical Volume("INTERIOR") = {  88,   57,    3,   26,    1,    2,
                                656,  625,  687,  718,  780,  749,
                                401,  432,  277,  308,  339,  370,
                               1338, 1276, 1214, 1183, 1245, 1307,
                                966,  842,  935, 1152, 1090,  904,
                                811,  873, 1121, 1028,  997, 1059,
                                494,  463,  618,  587,  556,  525,
                                167,  198,  229,  260,  136,  105};

//----------------------------
//
// Define boundary conditions
//
//----------------------------
Physical Surface("WALLS") = { 624,  500,  173,  469,  111,  266,  142,  531,
                               15,   94,   63,   14,  631,  755,  786,  724,
                              438,  283, 1189,  407,  376, 1313, 1251, 1344,
                             1065,  848,  879,  817,  972, 1127, 1158, 1003};
Physical Surface("Y_MIN") = { 843,  812,  874,  998, 1060, 1122,
                             1153,  967,  936,  905, 1029, 1091};
Physical Surface("Y_MAX") = { 168,  137,  106,  261,  526,  464,
                              495,  619,  199,  230,  557,  588};
Physical Surface("Z_MIN") = { 433,  402,  371,  278, 1184, 1246,
                             1308, 1339,  309, 1215, 1277,  340};
Physical Surface("Z_MAX") = {  89,    6,   58,    5,  719,  626,
                              750,  781,   27,  688,    4,  657};
