X_A =  0.0;
X_B =  2.0;
X_C = 16.0;
Y_A =  0.0;
Y_B =  2.0;
Z_A = -1.0;
Z_B =  0.0;
Z_C =  1.0;

//---------------
// Define points
//---------------

Point( 1) = {X_A,   Y_A,   Z_B};
Point( 2) = {X_B,   Y_A,   Z_B};
Point( 3) = {X_A,   Y_B,   Z_B};
Point( 4) = {X_B,   Y_B,   Z_B};
Point( 5) = {X_A,   Y_A,   Z_C};
Point( 6) = {X_B,   Y_A,   Z_C};
Point( 7) = {X_A,   Y_B,   Z_C};
Point( 8) = {X_B,   Y_B,   Z_C};
Point( 9) = {X_B,   Y_A,   Z_A};
Point(10) = {X_C,   Y_A,   Z_A};
Point(11) = {X_B,   Y_B,   Z_A};
Point(12) = {X_C,   Y_B,   Z_A};
Point(13) = {X_C,   Y_A,   Z_B};
Point(14) = {X_C,   Y_B,   Z_B};
Point(15) = {X_C,   Y_A,   Z_C};
Point(16) = {X_C,   Y_B,   Z_C};

//--------------
// Define lines
//--------------

// Lines at the inlet
Line(1) = {1, 2};
Line(2) = {3, 4};
Line(3) = {5, 6};
Line(4) = {7, 8};

// Lines at the outlet
Line(5) = {9, 10};
Line(6) = {2, 13};
Line(7) = {6, 15};
Line(8) = {11, 12};
Line(9) = {4, 14};
Line(10) = {8, 16};

// Lines normal to the walls
Line(11) = {1, 5};
Line(12) = {2, 6};
Line(13) = {3, 7};
Line(14) = {4, 8};
Line(15) = {9, 2};
Line(16) = {10, 13};
Line(17) = {11, 4};
Line(18) = {12, 14};
Line(19) = {13, 15};
Line(20) = {14, 16};

// Lines in spanwise direction
Line(21) = {1, 3};
Line(22) = {2, 4};
Line(23) = {9, 11};
Line(24) = {10, 12};
Line(25) = {13, 14};
Line(26) = {15, 16};
Line(27) = {6, 8};
Line(28) = {5, 7};

//------------------------
// Set lines' resolutions
//------------------------

// At the inlet
Transfinite Curve {1, 2, 3, 4} = 41 Using Progression 1;

// At the outlet
Transfinite Curve {5, 6, 7, 8, 9, 10} = 161 Using Progression 1.005;

// Normal to the walls
Transfinite Curve {11, 12, 13, 14, 15, 16, 17, 18, 19, 20} = 21 Using Progression 1;

// In spanwise direction
Transfinite Curve {21, 22, 23, 24, 25, 26, 27, 28} = 41 Using Progression 1;

//-----------------
// Define surfaces
//-----------------

// In y direction
Curve Loop(1) = {1, 12, -3, -11};   Plane Surface(1) = {1};
Curve Loop(2) = {2, 14, -4, -13};   Plane Surface(2) = {2};
Curve Loop(3) = {5, 16, -6, -15};   Plane Surface(3) = {3};
Curve Loop(4) = {8, 18, -9, -17};   Plane Surface(4) = {4};
Curve Loop(5) = {6, 19, -7, -12};   Plane Surface(5) = {5};
Curve Loop(6) = {9, 20, -10, -14};  Plane Surface(6) = {6};

// In x direction
Curve Loop(7) = {21, 13, -28, -11};   Plane Surface(7) = {7};
Curve Loop(8) = {22, 14, -27, -12};   Plane Surface(8) = {8};
Curve Loop(9) = {23, 17, -22, -15};   Plane Surface(9) = {9};
Curve Loop(10) = {24, 18, -25, -16};  Plane Surface(10) = {10};
Curve Loop(11) = {25, 20, -26, -19};  Plane Surface(11) = {11};

// In z direction
Curve Loop(12) = {1, 22, -2, -21};   Plane Surface(12) = {12};
Curve Loop(13) = {5, 24, -8, -23};   Plane Surface(13) = {13};
Curve Loop(14) = {6, 25, -9, -22};   Plane Surface(14) = {14};
Curve Loop(15) = {3, 27, -4, -28};   Plane Surface(15) = {15};
Curve Loop(16) = {7, 26, -10, -27};  Plane Surface(16) = {16};

Transfinite Surface {1} = {1, 2, 6, 5};
Transfinite Surface {3} = {9, 10, 13, 2};
Transfinite Surface {5} = {2, 13, 15, 6};
Transfinite Surface {2} = {3, 4, 8, 7};
Transfinite Surface {4} = {11, 12, 14, 4};
Transfinite Surface {6} = {4, 14, 16, 8};
Transfinite Surface {7} = {1, 3, 7, 5};
Transfinite Surface {8} = {2, 4, 8, 6};
Transfinite Surface {9} = {9, 11, 4, 2};
Transfinite Surface {11} = {13, 14, 16, 15};
Transfinite Surface {10} = {10, 12, 14, 13};
Transfinite Surface {12} = {1, 2, 4, 3};
Transfinite Surface {13} = {9, 10, 12, 11};
Transfinite Surface {14} = {2, 13, 14, 4};
Transfinite Surface {15} = {5, 6, 8, 7};
Transfinite Surface {16} = {6, 15, 16, 8};

Recombine Surface "*";

//----------------
// Define volumes
//----------------
Surface Loop(1) = {1, 12, 2, 15, 7, 8};   Volume(1) = {1};
Surface Loop(2) = {13, 3, 10, 4, 9, 14};  Volume(2) = {2};
Surface Loop(3) = {16, 5, 11, 6, 14, 8};  Volume(3) = {3};

Transfinite Volume{1} = {1, 2, 4, 3, 5, 6, 8, 7};
Transfinite Volume{2} = {9, 10, 12, 11, 2, 13, 14, 4};
Transfinite Volume{3} = {2, 13, 14, 4, 6, 15, 16, 8};

//--------------------------------
// Boundary and volume conditions
//--------------------------------
Physical Surface("in", 29) = {7};
Physical Surface("out", 30) = {10, 11};
Physical Surface("low_wall", 31) = {12, 9, 13};
Physical Surface("top_wall", 32) = {15, 16};
Physical Surface("periodic", 33) = {1, 2, 5, 6, 3, 4};
Physical Volume("fluid", 34) = {1, 2, 3};
