// Specify mesh size
DELTA = 0.06;

//--------
//
// Points
//
//--------

// Points at the bottom
Point(1) = {-0.4,  0,    0,  DELTA};
Point(2) = { 0.4,  0,    0,  DELTA};
Point(3) = { 0.8,  0.6,  0,  DELTA};
Point(4) = {-0.8,  0.6,  0,  DELTA};

// Bottom arc's center point
Point(5)  = { 0.0, -0.1, 0};

// Points at the top
Point(6) = {-0.3,  -0.3,   0.5,  DELTA};
Point(7) = { 0.3,  -0.3,   0.5,  DELTA};
Point(8) = { 0.6,   0.9,   0.5,  DELTA};
Point(9) = {-0.6,   0.9,   0.5,  DELTA};

// Top arc's center point
Point(10) = { 0.0, -0.375, 0.5};

//------------
//
// Connectors
//
//------------

Circle(1) = {1, 5, 2};
Line(2)   = {2, 3};
Line(3)   = {3, 4};
Line(4)   = {4, 1};

Circle(5) = {6, 10, 7};
Line(6)   = {7, 8};
Line(7)   = {8, 9};
Line(8)   = {9, 6};

// Lines connecting bottom and top
Line(9) = {2, 7};
Line(10) = {3, 8};
Line(11) = {4, 9};
Line(12) = {1, 6};

//----------
//
// Surfaces
//
//----------

// Circular dent
Curve Loop(1) = {1, 9, -5, -12};  Surface(1) = {1};

// Right and left
Curve Loop(2) = {2, 10, -6, -9};  Surface(2) = {2};
Curve Loop(3) = {11, 8, -12, -4};  Surface(3) = {3};

// Behind
Curve Loop(4) = {3, 11, -7, -10};  Surface(4) = {4};

// Bottom and top
Curve Loop(5) = {1, 2, 3, 4};  Surface(5) = {5};
Curve Loop(6) = {5, 6, 7, 8};  Surface(6) = {6};

Surface Loop(1) = {5, 1, 2, 4, 3, 6};
Volume(1) = {1};

//-----------------
//
// Physical groups
//
//-----------------
Physical Surface("bottom") = {5};
Physical Surface("top") = {6};
Physical Surface("right") = {2};
Physical Surface("left") = {3};
Physical Surface("front") = {1};
Physical Surface("back") = {4};
Physical Volume("wedge") = {1};
