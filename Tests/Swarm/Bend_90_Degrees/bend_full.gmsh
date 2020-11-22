// Gmsh project created on Wed Aug 12 17:50:13 2020
SetFactory("OpenCASCADE");

//++++++++++++++++++++++++++++++++//
// User's input for customization //
//++++++++++++++++++++++++++++++++//
//                                //
//    Geometry & mesh controls    //
//                                //
//++++++++++++++++++++++++++++++++//

// Diameters parametrizing the geometry
r_c    = 0.04;          // radius of curvature
dia_1  = 0.02;          // pipe outer diameter (diameter 1)
rad_1  = dia_1/2.0;     // outer radius (r1)
rad_2  = rad_1 * 0.9 ;  // inner radius (r2)
rad_3  = rad_1 * 0.6 ;  // quadron radius (r3)
r_a    = Pi/2;          // Angle of rotation (for the bend)

arm_i  = 5 * dia_1  ;   // inner pipe arm length 
arm_o  = 5 * dia_1  ;   // outer pipe arm length 

len_i  = arm_i + r_c;   // Length of inlet pipe (from origin)
len_o  = arm_o + r_c;   // Length of inlet pipe (from origin)

// Number of nodes on each part of the ciscumference (all)
n_c  = 10;      // number of nodes on each arc
n_bl =  6;      // number of nodes from inner to outer diameter (BL)
n_ci = 10;      // number of nodes from octagon core to inner diameter
n_cc = n_c;     // number of nodes from center to octagon core

// Number of layers (for refining after surface extrusion)
la_i = 100;     // number of layers in the inlet pipe
la_o = la_i;    // number of layers in the outlet pipe
la_b =  50;     // number of layers in the bend section

// Progression (mesh clustering)
p_bl = -0.8;    // boundary layer progression
p_ci =  1.0;    // progression core to inner circum
p_cc =  1.0;    // progression center to octagon core


//================//------------------------------------------------------------
// Geo. algorithm //      Please don't change nothing from this point
//================//------------------------------------------------------------

// 5 points defining the central wire
Point(1) = {0.0, r_c, 0.0, 1.0};
Point(2) = {r_c, 0.0, 0.0, 1.0};
Point(3) = {0.0, len_i, 0.0, 1.0};
Point(4) = {len_o, 0.0, 0.0, 1.0};
Point(5) = {r_c, r_c, 0.0, 1.0};

// Center-line wire (needed for extrusion of the 2d plane y=0.14 at top)
Line(1) = {3, 1};
Line(2) = {2, 4};
Circle(3) = {1, 5, 2};

// Top-circle points (to draw the circumference)
Point(6) = {rad_1,  len_i, 0.0, 1.0};
Point(7) = {0.0,   len_i, -rad_1, 1.0};
Point(8) = {-rad_1, len_i, 0.0, 1.0};
Point(9) = {0.0,   len_i, rad_1, 1.0};

//or1= rad_1 * cos(45);

// Points dividing the outer radius quarters 
Point(10) = { 0.0070710676812  , len_i,  -0.0070710676812, 1.0};
Point(11) = {-0.0070710676812  , len_i,  -0.0070710676812, 1.0};
Point(12) = {-0.0070710676812  , len_i,   0.0070710676812, 1.0};
Point(13) = { 0.0070710676812  , len_i,   0.0070710676812, 1.0};

//or2= rad_2 * cos(45);

// Inner radius -> boundary layer refinement
Point(14) = {rad_2, len_i, 0.0, 1.0};
Point(15) = {-rad_2, len_i, 0.0, 1.0};
Point(16) = {0.0, len_i, -rad_2, 1.0};
Point(17) = {0.0, len_i, rad_2, 1.0};
Point(18) = {0.006363961031,  len_i, -0.006363961031, 1.0};
Point(19) = {-0.006363961031, len_i, 0.006363961031, 1.0};
Point(20) = {-0.006363961031, len_i, -0.006363961031, 1.0};
Point(21) = {0.006363961031,  len_i, 0.006363961031, 1.0};

// Inner circle parts
Circle(8) = {14, 3, 18};
Circle(9) = {18, 3, 16};
Circle(10) = {16, 3, 20};
Circle(11) = {20, 3, 15};
Circle(12) = {15, 3, 19};
Circle(13) = {19, 3, 17};
Circle(14) = {17, 3, 21};
Circle(15) = {21, 3, 14};

// Outer circle parts
Circle(16) = {8, 3, 12};
Circle(17) = {12, 3, 9};
Circle(18) = {9, 3, 13};
Circle(19) = {13, 3, 6};
Circle(20) = {6, 3, 10};
Circle(21) = {10, 3, 7};
Circle(22) = {7, 3, 11};
Circle(23) = {11, 3, 8};

//or3= rad_3 * cos(45);

// Quadron points (for core lines construction)
Point(22) = {rad_3,  len_i, 0.0, 1.0};
Point(23) = {0.0,    len_i, -rad_3, 1.0};
Point(24) = {-rad_3, len_i, 0.0, 1.0};
Point(25) = {0.0,    len_i, rad_3, 1.0};
Point(26) = {0.004242640687,  len_i, -0.004242640687, 1.0};
Point(27) = {-0.004242640687, len_i, -0.004242640687, 1.0};
Point(28) = {-0.004242640687, len_i, 0.004242640687, 1.0};
Point(29) = {0.004242640687,  len_i, 0.004242640687, 1.0};

// Quadron lines
Line(24) = {22, 26};
Line(25) = {26, 23};
Line(26) = {23, 27};
Line(27) = {27, 24};
Line(28) = {24, 28};
Line(29) = {28, 25};
Line(30) = {25, 29};
Line(31) = {29, 22};

//=====================//
// Center to boundries //
//=====================//

// Inner to outer lines 
Line(32) = {14, 6};
Line(33) = {18, 10};
Line(34) = {16, 7};
Line(35) = {20, 11};
Line(36) = {15, 8};
Line(37) = {19, 12};
Line(38) = {17, 9};
Line(39) = {21, 13};

// Quadron to inner lines
Line(40) = {14, 22};
Line(41) = {18, 26};
Line(42) = {16, 23};
Line(43) = {20, 27};
Line(44) = {15, 24};
Line(45) = {19, 28};
Line(46) = {17, 25};
Line(47) = {21, 29};

// Center to quadron lines
Line(48) = {14, 3};
Line(49) = {18, 3};
Line(50) = {16, 3};
Line(51) = {20, 3};
Line(52) = {15, 3};
Line(53) = {19, 3};
Line(54) = {17, 3};
Line(55) = {21, 3};

//================//
// Plane surfaces //
//================//
Curve Loop(3) = {8, 41, -24, -40};
Plane Surface(2) = {3};

Curve Loop(4) = {41, 25, -42, -9};
Plane Surface(3) = {4};

Curve Loop(5) = {42, 26, -43, -10};
Plane Surface(4) = {5};

Curve Loop(6) = {43, 27, -44, -11};
Plane Surface(5) = {6};
 
Curve Loop(7) = {12, 45, -28, -44};
Plane Surface(6) = {7};

Curve Loop(8) = {45, 29, -46, -13};
Plane Surface(7) = {8};

Curve Loop(9) = {30, -47, -14, 46};
Plane Surface(8) = {9};

Curve Loop(11) = {14, 39, -18, -38};
Plane Surface(10) = {11};

Curve Loop(12) = {32, -19, -39, 15};
Plane Surface(11) = {12};

Curve Loop(13) = {20, -33, -8, 32};
Plane Surface(12) = {13};

Curve Loop(14) = {21, -34, -9, 33};
Plane Surface(13) = {14};

Curve Loop(15) = {22, -35, -10, 34};
Plane Surface(14) = {15};

Curve Loop(16) = {11, 36, -23, -35};
Plane Surface(15) = {16};

Curve Loop(17) = {12, 37, -16, -36};
Plane Surface(16) = {17};

Curve Loop(18) = {17, -38, -13, 37};
Plane Surface(17) = {18};

Curve Loop(19) = {40, -31, -47, 15};
Plane Surface(18) = {19};

//============//
// Core lines //
//============//
Line(56) = {3, 29};
Line(57) = {3, 22};
Line(58) = {26, 3};
Line(59) = {23, 3};
Line(60) = {27, 3};
Line(61) = {24, 3};
Line(62) = {28, 3};
Line(63) = {25, 3};

// Partitioning into number of nodes
Transfinite Curve{20}  = n_c Using Progression 1;
Transfinite Curve{8}   = n_c Using Progression 1;
Transfinite Line{24}   = n_c Using Progression 1;
Transfinite Line{27}   = n_c Using Progression 1;
Transfinite Curve{11}  = n_c Using Progression 1;
Transfinite Curve{23}  = n_c Using Progression 1;

Transfinite Curve{19}  = n_c Using Progression 1;
Transfinite Curve{15}  = n_c Using Progression 1;
Transfinite Line{31}   = n_c Using Progression 1;
Transfinite Line{28}   = n_c Using Progression 1;
Transfinite Curve{12}  = n_c Using Progression 1;
Transfinite Curve{16}  = n_c Using Progression 1;

Transfinite Curve{21}  = n_c Using Progression 1;
Transfinite Curve{9}   = n_c Using Progression 1;
Transfinite Line{25}   = n_c Using Progression 1;
Transfinite Line{30}   = n_c Using Progression 1;
Transfinite Curve{14}  = n_c Using Progression 1;
Transfinite Curve{18}  = n_c Using Progression 1;

Transfinite Curve{22}  = n_c Using Progression 1;
Transfinite Curve{10}  = n_c Using Progression 1;
Transfinite Line{26}   = n_c Using Progression 1;
Transfinite Line{29}   = n_c Using Progression 1;
Transfinite Curve{13}  = n_c Using Progression 1;
Transfinite Curve{17}  = n_c Using Progression 1;

//=======================//
// Triangles at the core //
//=======================//
// -=-> Line Loop(20) = {56, 57, 31};
// -=-> Plane Surface(19) = {20};
// -=-> 
// -=-> Line Loop(22) = {57, 58, 24};
// -=-> Plane Surface(20) = {22};
// -=-> 
// -=-> Line Loop(24) = {58, 59, 25};
// -=-> Plane Surface(21) = {24};
// -=-> 
// -=-> Line Loop(26) = {59, 60, 26};
// -=-> Plane Surface(22) = {26};
// -=-> 
// -=-> Line Loop(28) = {60, 61, 27};
// -=-> Plane Surface(23) = {28};
// -=-> 
// -=-> Line Loop(30) = {61, 62, 28};
// -=-> Plane Surface(24) = {30};
// -=-> 
// -=-> Line Loop(32) = {62, 63, 29};
// -=-> Plane Surface(25) = {32};
// -=-> 
// -=-> Line Loop(34) = {63, 56, 30};
// -=-> Plane Surface(26) = {34};
//+++++++++++++++++++++++++++++++++++

Line Loop(20) = {57, 24, 25, 59};
Plane Surface(19) = {20};

Line Loop(22) = {59, 26, 27, 61};
Plane Surface(20) = {22};

Line Loop(24) = {61, 28, 29, 63};
Plane Surface(21) = {24};

Line Loop(26) = {63, 30, 31, 57};
Plane Surface(22) = {26};

//=================//
// 2D surface mesh //
//=================//
Transfinite Line{32,33,34,35,36,37,38,39}    = n_bl Using Progression p_bl;
Transfinite Line{40,41,42,43,44,45,46,47}    = n_ci Using Progression p_ci;
Transfinite Line{57,59,61,63}                = n_cc Using Progression p_cc;
Transfinite Surface{2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22};
Recombine Surface{2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22};

// Wire of extrusion 
Wire(29) = {1, 3, 2};

//===================//
// Surface extrusion //
//===================//
Extrude {0,-arm_i,0} {Surface{2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22}; Layers{la_i}; Recombine; }
Extrude{{0,0,1}, {r_c,r_c,0}, r_a} {Surface{81,86,85,83,27,78,51,47,43,39,35,31,62,59,55,76,74,71,68,65}; Layers{la_b}; Recombine; }
Extrude{arm_o,0,0} {Surface{102,91,95,99,121,124,126,106,109,112,115,118,145,148,150,130,133,136,139,142}; Layers{la_o}; Recombine; }
Extrude { Line{56,58,60,62};    } Using Wire {29}

 // Not used extrusions but maybe needed in the future!
 //Extrude { Surface{2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22}; } Using Wire {29}
 //Extrude { Curve{8,9,10,11,12,13,14,15};   } Using Wire {29}
 //Extrude { Curve{16,17,18,19,20,21,22,23}; } Using Wire {29}
 //Extrude { Line{40,41,42,43,44,45,46,47};  } Using Wire {29}
 //Extrude { Line{24,25,26,27,28,29,30,31};  } Using Wire {29}
 //Extrude { Line{56,57,58,59,60,61,62,63};  } Using Wire {29}
 //Extrude { Line{32,33,34,35,36,37,38,39};  } Using Wire {29}

//=====================//
// Boundary conditions //
//=====================//
Physical Surface("bend_inlet") = {21, 20, 5, 4, 14, 13, 12, 2, 19, 22, 18, 11, 10, 8, 7, 17, 16, 6, 15, 3};
Physical Surface("bend_walls") = {63, 75, 53, 57, 60, 66, 70, 73, 127, 149, 146, 143, 140, 137, 134, 131, 195, 198, 201, 204, 207, 210, 213, 192};
Physical Surface("bend_outlet") = {173, 176, 188, 212, 190, 214, 166, 163, 159, 155, 170, 194, 197, 200, 179, 182, 206, 185, 209, 203};
Physical Volume("fluid") = {12, 3, 13, 4, 14, 5, 19, 6, 15, 7, 8, 20, 16, 9, 10, 1, 17, 18, 2, 11, 40, 39, 31, 32, 33, 25, 21, 24, 30, 38, 37, 29, 23, 22, 26, 34, 27, 35, 28, 36, 55, 47, 54, 46, 48, 56, 42, 41, 45, 53, 43, 57, 49, 44, 52, 60, 51, 59, 58, 50};
