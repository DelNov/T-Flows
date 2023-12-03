// Gmsh project created on Thu Dec  5 21:07:29 2423
SetFactory("OpenCASCADE");

Mesh.Algorithm = 8;

//----------------------------------------
//
//   POINTS & LINES
//
//----------------------------------------

//Outer
Point(1) = {-72.5, -75.5, 0, 5};
Point(2) = {127.5, -75.5, 0, 5};
Point(3) = {127.5, 74.5, 0, 5};
Point(4) = {-72.5, 74.5, 0, 5};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};

//Boundary Control
Point(1000) = {-57.5, -60.5, 0, 5};
Point(2000) = {112.5, -60.5, 0, 5};
Point(3000) = {112.5, 59.5, 0, 5};
Point(4000) = {-57.5, 59.5, 0, 5};

Line(1000) = {1000, 2000};
Line(2000) = {2000, 3000};
Line(3000) = {3000, 4000};
Line(4000) = {4000, 1000};

Curve Loop(1000) = {1000, 2000, 3000, 4000};

//Refinement NW
Point(13) = {0, -1.5, 0, 1.0};
Point(14) = {1.5, -1.5, 0, 1.0};
Point(15) = {3, -1.5, 0, 1.0};
Point(16) = {5, -1.5, 0, 1.0};
Point(17) = {5.5, -1, 0, 1.0};
Point(18) = {5.5, 0, 0, 1.0};
Point(19) = {5, 0.5, 0, 1.0};
Point(20) = {3, 0.5, 0, 1.0};
Point(21) = {1.5, 0.5, 0, 1.0};
Point(22) = {0, 0.5, 0, 1.0};
Point(23) = {-0.5, 0, 0, 1.0};
Point(24) = {-0.5, -1, 0, 1.0};
Point(25) = {0, -1, 0, 1.0};
Point(28) = {5, -1, 0, 1.0};
Point(29) = {5, 0, 0, 1.0};
Point(32) = {0, 0, 0, 1.0};

Line(13) = {13, 14};
Line(14) = {14, 15};
Line(15) = {15, 16};
Circle(16) = {16, 28, 17};
Line(17) = {17, 18};
Circle(18) = {18, 29, 19};
Line(19) = {19, 20};
Line(20) = {20, 21};
Line(21) = {21, 22};
Circle(22) = {22, 32, 23};
Line(23) = {23, 24};
Circle(24) = {24, 25, 13};

Curve Loop(3) = {13:24};

//Walls
Point(26) = {1.5, -1, 0, 1.0};
Point(27) = {3, -1, 0, 1.0};
Point(30) = {3, 0, 0, 1.0};
Point(31) = {1.5, 0, 0, 1.0};

Line(25) = {25, 26};
Line(26) = {26, 27};
Line(27) = {27, 28};
Line(28) = {28, 29};
Line(29) = {29, 30};
Line(30) = {30, 31};
Line(31) = {31, 32};
Line(32) = {32, 25};

Line(33) = {25, 13};
Line(34) = {26, 14};
Line(35) = {27, 15};
Line(36) = {28, 16};
Line(37) = {28, 17};
Line(38) = {29, 18};
Line(39) = {29, 19};
Line(40) = {30, 20};
Line(41) = {31, 21};
Line(42) = {32, 22};
Line(43) = {32, 23};
Line(44) = {25, 24};

//Refinement NW2
Point(33) = {0, -2.5, 0, 1.0};
Point(34) = {1.5, -2.5, 0, 1.0};
Point(35) = {3, -2.5, 0, 1.0};
Point(36) = {5, -2.5, 0, 1.0};
Point(37) = {6.5, -1, 0, 1.0};
Point(38) = {6.5, 0, 0, 1.0};
Point(39) = {5, 1.5, 0, 1.0};
Point(40) = {3, 1.5, 0, 1.0};
Point(41) = {1.5, 1.5, 0, 1.0};
Point(42) = {0, 1.5, 0, 1.0};
Point(43) = {-1.5, 0, 0, 1.0};
Point(44) = {-1.5, -1, 0, 1.0};

Line(53) = {33, 34};
Line(54) = {34, 35};
Line(55) = {35, 36};
Circle(56) = {36, 28, 37};
Line(57) = {37, 38};
Circle(58) = {38, 29, 39};
Line(59) = {39, 40};
Line(60) = {40, 41};
Line(61) = {41, 42};
Circle(62) = {42, 32, 43};
Line(63) = {43, 44};
Circle(64) = {44, 25, 33};

Curve Loop(4) = {53:64};

Line(73) = {13, 33};
Line(74) = {14, 34};
Line(75) = {15, 35};
Line(76) = {16, 36};
Line(77) = {17, 37};
Line(78) = {18, 38};
Line(79) = {19, 39};
Line(80) = {20, 40};
Line(81) = {21, 41};
Line(82) = {22, 42};
Line(83) = {23, 43};
Line(84) = {24, 44};

//-------------------------------------------
//
//   CURVE LOOPS & SURFACES
//
//----------------------------------------

//Outer
Plane Surface(1) = {1, 1000};

//Inner
Plane Surface(1000) = {1000, 4};

//Refinement NW
Curve Loop(5) = {13, -34, -25, 33};
Plane Surface(4) = {5};
Curve Loop(6) = {14, -35, -26, 34};
Plane Surface(5) = {6};
Curve Loop(7) = {15, -36, -27, 35};
Plane Surface(6) = {7};
Curve Loop(8) = {36, 16, -37};
Plane Surface(7) = {8};
Curve Loop(9) = {37, 17, -38, -28};
Plane Surface(8) = {9};
Curve Loop(10) = {38, 18, -39};
Plane Surface(9) = {10};
Curve Loop(11) = {-29, 39, 19, -40};
Plane Surface(10) = {11};
Curve Loop(912) = {-30, 40, 20, -41};
Plane Surface(11) = {912};
Curve Loop(13) = {-31, 41, 21, -42};
Plane Surface(12) = {13};
Curve Loop(14) = {-43, 42, 22};
Plane Surface(13) = {14};
Curve Loop(15) = {-44, -32, 43, 23};
Plane Surface(14) = {15};
Curve Loop(16) = {-33, 44, 24};
Plane Surface(15) = {16};

//Refinement NW2
Curve Loop(25) = {64, -73, -24, 84};
Plane Surface(24) = {25};
Curve Loop(26) = {53, -74, -13, 73};
Plane Surface(25) = {26};
Curve Loop(27) = {54, -75, -14, 74};
Plane Surface(26) = {27};
Curve Loop(28) = {55, -76, -15, 75};
Plane Surface(27) = {28};
Curve Loop(29) = {56, -77, -16, 76};
Plane Surface(28) = {29};
Curve Loop(30) = {57, -78, -17, 77};
Plane Surface(29) = {30};
Curve Loop(31) = {58, -79, -18, 78};
Plane Surface(30) = {31};
Curve Loop(32) = {59, -80, -19, 79};
Plane Surface(31) = {32};
Curve Loop(33) = {60, -81, -20, 80};
Plane Surface(32) = {33};
Curve Loop(34) = {61, -82, -21, 81};
Plane Surface(33) = {34};
Curve Loop(35) = {62, -83, -22, 82};
Plane Surface(34) = {35};
Curve Loop(36) = {63, -84, -23, 83};
Plane Surface(35) = {36};

//-------------------------------------------
//
//   TRANSFINITE SURFACES
//
//----------------------------------------
Transfinite Surface {4:6, 8, 10:12, 14, 24:35};

//-------------------------------------------
//
//   TRANSFINITE CURVES
//
//----------------------------------------

//Radial (Wall Normal)
Transfinite Curve {33:44} = 27 Using Progression 1.08; //Inner
Transfinite Curve {73:84} = 14 Using Progression 1.08; //Outer

//Round Corners
Transfinite Curve {16, 18, 22, 24} = 27; //Inner
Transfinite Curve {56, 58, 62, 64} = 27 Using Bump 0.475; //Outer

//Cylinder Orthogonal Walls
Transfinite Curve {63, 23, 32, 28, 17} = 15 Using Bump 0.15;
Transfinite Curve {57} = 15;

	//Leading Edge
		//Upper
		Transfinite Curve {31, 21, 61} = 28 Using Progression 1/1.06;
		//Lower
		Transfinite Curve {25, 13, 53} = 28 Using Progression 1.06;
		
	//Center
		//Upper
		Transfinite Curve {30, 20, 60} = 12 Using Progression 1/1.02;
		//Lower
		Transfinite Curve {26, 14, 54} = 12 Using Progression 1.02;
		
	//Trailing Edge
		//Upper
		Transfinite Curve {29, 19, 59} = 28 Using Progression 1.06;
		//Lower
		Transfinite Curve {27, 15, 55} = 28 Using Progression 1/1.06;
	
	

//-------------------------------------------
//
//   SURFACE RECOMBINATION
//
//----------------------------------------
Recombine Surface {1, 4:15, 24:35, 1000};

//-------------------------------------------
//
//  REFINEMENT BOXES
//
//----------------------------------------
Field[1] = Box;
Field[1].VIn =  0.3;
Field[1].VOut = 5;
Field[1].XMin = -10;
Field[1].XMax = 60;
Field[1].YMin = -8.5;
Field[1].YMax = 7.5;
Field[1].Thickness = 10;

Field[2] = Box;
Field[2].VIn =  0.15;
Field[2].VOut = 5;
Field[2].XMin = -5;
Field[2].XMax = 30;
Field[2].YMin = -3.5;
Field[2].YMax = 2.5;
Field[2].Thickness = 10;

Field[3] = Box;
Field[3].VIn = 0.06;
Field[3].VOut = 5;
Field[3].XMin = 5;
Field[3].XMax = 5.1;
Field[3].YMin = -1.3;
Field[3].YMax = 0.3;
Field[3].Thickness = 0.1;

Field[4] = Box;
Field[4].VIn = 0.028;
Field[4].VOut = 5;
Field[4].XMin = -1;
Field[4].XMax = 6;
Field[4].YMin = -2;
Field[4].YMax = 1;
Field[4].Thickness = 0.1;

Field[9] = Min;
Field[9].FieldsList = {1, 2, 3, 4};
Background Field = 9;

//-------------------------------------------
//
//   SURFACE EXTRUSION
//
//----------------------------------------
Extrude {0, 0, 20} {
  Surface{1, 4:15, 24:35, 1000};
  Layers{148};
  Recombine;
}

//-------------------------------------------
//
//   PHYSICAL BOUNDARIES
//
//----------------------------------------
Physical Surface("IN") = {1004};
Physical Surface("OUT") = {1002};
Physical Surface("WALL") = {1012, 1017, 1021, 1028, 1033, 1037, 1041, 1050};
Physical Surface("TOP") = {1003};
Physical Surface("BOTTOM") = {1001};
Physical Surface("PERIODIC") = {1, 4:15, 24:35, 1000, 1009, 1090, 1057, 1060, 1063, 1066, 1069, 1072, 1075, 1078, 1081, 1084, 1087, 1089, 1014, 1018, 1022, 1025, 1029, 1032, 1036, 1040, 1044, 1047, 1051, 1053};
Physical Volume("FLUID") = {1:26};
