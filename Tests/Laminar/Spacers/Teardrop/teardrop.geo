//------------------------------------------------
//
// Characteristic dimensions and grid resolutions
//
// (See figure teardrop.fig)
//
//------------------------------------------------
L =  5;       // length of the computational domain
H =  1.5;     // height of the computational domain
R1 = 0.5;     // leading edge radius
R2 = R1/4;    // trailing edge radis
P  = R1*1.5;  // pitch between radii

NH_QUARTER = 10;  // number of cells along quarter-height
NS         = 14;  // number of cells covering spacers
NB_HALF    = 40;  // number of cells between spacers

//-------------------------
//
// Basic settings for GMSH
//
//-------------------------
SetFactory("OpenCASCADE");
Geometry.Tolerance = 1e-8; // adjust value here for correct merge result
Coherence Mesh;

//-----------------
//
// Define geometry
//
//-----------------

//----------------------
// Create airfoil shape
//----------------------

// https://en.wikipedia.org/wiki/Tangent_lines_to_circles#Tangent_lines_to_two_circles
N = 8;  // must be even

x1 = 0.0;
y1 = 0.0;

x2 = P;
y2 = 0.0;   // must be zero, algorithm is not general enough

// Calculate angles
gama = -ArcTan((y2-y1)/(x2-x1));
beta = ArcSin((R1-R2) / Sqrt((x2-x1)^2 + (y2-y1)^2));
alfa = gama - beta;

// I find it more intuitive to have alfa as positive
alfa = -alfa;

x3 = x1 + R1*Sin(alfa);
y3 = y1 + R1*Cos(alfa);

x4 = x2 + R2*Sin(alfa);
y4 = y2 + R2*Cos(alfa);

Printf("gama = %f", gama);
Printf("beta = %f", beta * 180/Pi);
Printf("alfa = %f", alfa * 180/Pi);

// Initialize point counter
pc = 0;

// Start from the trailing edge
angle0 = 0;
For n In{ 0 : N }
  d = (Pi/2 - alfa) / N;
  angle = angle0 + d*n;
  Printf("angle = %f", angle * 180/Pi);
  x = x2 + R2*Cos(angle);
  y = y2 + R2*Sin(angle);
  pc++;
  Point(pc) = {x, 0, y};
EndFor

// Straight line
dx = (x4 - x3) / N;
dy = (y4 - y3) / N;
For n In{ 1 : N-1 }
  x = x4 - dx * n;
  y = y4 - dy * n;
  pc++;
  Point(pc) = {x, 0, y};
EndFor

// End with a leading edge
angle0 = (Pi/2 - alfa);
For n In{ 0 : N }
  d = (Pi/2 + alfa) / N;
  angle = angle0 + d*n;
  x = x1 + R1*Cos(angle);
  y = y1 + R1*Sin(angle);
  pc++;
  Point(pc) = {x, 0, y};
EndFor

//---------------------------------
// Box defining the overall domain
//---------------------------------
Box(1) = {0, 0, 0,   L, L, H};

// Cut the domin with rectangle in two parts
Rectangle(7) = {0, 0, H/2, L, L, 0};
BooleanFragments{Volume{1}; Surface{7}; Delete; }{}

//--------------------
// Now define spacers
//--------------------
Spline(1) = {1:N*3+1};
Line(2) = {N*3+1, 1};
Recursive Delete {Point{1:N*3+1};}

Curve Loop(18) = {1, 2}; Plane Surface(18) = {18};
Dilate {{0, 0, 0}, {1, 1, -1}} {Duplicata { Surface{18}; }}
Extrude {0, L+H, 0} {Surface{18}; Surface{19};}

// Put them all in place
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/2} {Duplicata { Volume{4}; Volume{3}; }}
Translate {L/2, -H/2, H/4} {Volume{4}; Volume{3};}
Translate {L+H/2, L/2, 3*H/4} {Volume{5}; Volume{6};}

BooleanUnion{ Volume{4}; Delete; }{ Volume{3}; Delete; }
BooleanUnion{ Volume{5}; Delete; }{ Volume{6}; Delete; }
BooleanUnion{ Volume{7}; Delete; }{ Volume{8}; Delete; }

// Cut the spacers from the domain
BooleanDifference{Volume{2}; Delete;}{ Volume{3};}
BooleanDifference{Volume{1}; Delete;}{ Volume{3};}

// Delete spacers - they are void, not part of a simulation
Recursive Delete {Volume{3};}

Coherence;

//--------------------------------------------------
// Further cuts are needed to make hexahedral grids
//--------------------------------------------------

Extrude {0, 0, -H} {Point{5}; Point{15}; Point{14}; Point{20}; }
Extrude {0, 0, +H} {Point{8}; Point{9}; Point{27}; Point{28}; }

Line(78) = {42, 44};  Line(79) = {43, 45};
Line(80) = {5, 14};   Line(81) = {15, 20};
Line(82) = {46, 48};  Line(83) = {47, 49};
Line(84) = {8, 27};   Line(85) = {9, 28};

Curve Loop(33) = {70, 78, -72, -80};  Plane Surface(33) = {33};
Curve Loop(34) = {71, 79, -73, -81};  Plane Surface(34) = {34};
Curve Loop(35) = {74, 82, -76, -84};  Plane Surface(35) = {35};
Curve Loop(36) = {83, -77, -85, 75};  Plane Surface(36) = {36};

Coherence;

BooleanFragments{ Volume{7}; Delete; }{ Surface{33}; Delete; }
BooleanFragments{ Volume{8}; Delete; }{ Surface{34}; Delete; }
BooleanFragments{ Volume{6}; Delete; }{ Surface{49}; Delete; }
BooleanFragments{ Volume{11}; Delete; }{ Surface{61}; Delete; }
BooleanFragments{ Volume{4}; Delete; }{ Surface{35}; Delete; }
BooleanFragments{ Volume{14}; Delete; }{ Surface{36}; Delete; }
BooleanFragments{ Volume{5}; Delete; }{ Surface{15}; Delete; }
BooleanFragments{ Volume{5}; Delete; }{ Surface{99}; Delete; }
BooleanFragments{ Volume{17}; Delete; }{ Surface{111}; Delete; }

Recursive Delete {Surface{74, 86, 124, 136}; }

Rectangle(201) = {0, 0,   H/4, L, L, 0};
Rectangle(202) = {0, 0, 3*H/4, L, L, 0};

BooleanFragments{ Volume{7}; Delete; }{ Surface{201}; }
BooleanFragments{ Volume{8}; Delete; }{ Surface{201}; }
BooleanFragments{ Volume{9}; Delete; }{ Surface{201}; }
BooleanFragments{ Volume{10}; Delete; }{ Surface{201}; }
BooleanFragments{ Volume{11}; Delete; }{ Surface{201}; }
BooleanFragments{ Volume{12}; Delete; }{ Surface{201}; }

BooleanFragments{ Volume{13}; Delete; }{ Surface{202}; }
BooleanFragments{ Volume{14}; Delete; }{ Surface{202}; }
BooleanFragments{ Volume{15}; Delete; }{ Surface{202}; }
BooleanFragments{ Volume{16}; Delete; }{ Surface{202}; }
BooleanFragments{ Volume{17}; Delete; }{ Surface{202}; }
BooleanFragments{ Volume{18}; Delete; }{ Surface{202}; }

Recursive Delete {
  Surface{201}; Surface{202}; 
  Surface{213}; Surface{281}; 
  Surface{225}; Surface{292}; 
  Surface{236}; Surface{304}; 
  Surface{247}; Surface{315}; 
  Surface{259}; Surface{326}; 
  Surface{270}; Surface{338}; 
}

// Additional intervention is needed to have continuous splines
Recursive Delete {Volume{27}; }
Translate {-L/2, 0, -3*H/4} {Duplicata { Volume{35}; }}
Rotate {{0, 1, 0}, {0, 0, 0}, Pi} {Volume{43}; }
Rotate {{0, 0, 1}, {0, 0, 0}, -Pi/2} {Volume{43}; }
Translate {0, L/2, H/4} {Volume{43};}

Coherence;

//-------------
//
// Define mesh
//
//-------------

//--------
// Curves
//--------

// In z direction; layer by layer; from bottom up
Transfinite Curve {81, 87, 45, 43, 83, 90, 53, 48, 102, 106, 69, 66, 120, 123, 
                   78, 76} = NH_QUARTER Using Progression 1;
Transfinite Curve {93, 96, 56, 55, 95, 100, 58, 57, 111, 118, 64, 60, 109, 113, 
                   73, 71} = NH_QUARTER Using Progression 1;
Transfinite Curve {125, 128, 143, 141, 127, 132, 148, 145, 163, 179, 166, 181, 
                   165, 170, 183, 186} = NH_QUARTER Using Progression 1;
Transfinite Curve {133, 136, 151, 149, 135, 140, 156, 153, 171, 175, 187, 188, 
                   173, 178, 191, 193} = NH_QUARTER Using Progression 1;

// In x direction, without leading and trailing lines
Transfinite Curve {176, 174, 138, 137, 168, 167, 130, 129, 114, 117, 98, 97, 
                   122, 107, 91, 86} = NB_HALF Using Progression 1;
Transfinite Curve {194, 189, 154, 150, 184, 180, 147, 142, 39, 28, 17, 10, 
                   79, 68, 54, 46} = NB_HALF Using Progression 1;

// In x direction, leading edges (these might need compression)
Transfinite Curve {115, 104, 88, 85} = NB_HALF Using Progression 1;

// In x direction, trailing edges (these might need compression)
Transfinite Curve {44, 51, 63, 74} = NB_HALF Using Progression 1;

// In y direction, without leading and trailing lines
Transfinite Curve {134, 139, 155, 152, 94, 99, 21, 20, 82, 89, 50, 47, 84, 92, 
                   52, 49} = NB_HALF Using Progression 1;
Transfinite Curve {172, 177, 190, 192, 112, 119, 42, 41, 110, 116, 75, 72, 121, 
                   124, 80, 77} = NB_HALF Using Progression 1;

// In y direction, leading edges (these might need compression)
Transfinite Curve {126, 146, 131, 144} = NB_HALF Using Progression 1;

// In y direction, trailing edges (these might need compression)
Transfinite Curve {164, 169, 182, 185} = NB_HALF Using Progression 1;

// On top of spacers in all directions
Transfinite Curve {201, 202, 62, 59, 101, 105, 65, 61, 103, 108, 70, 67, 161, 
                   162, 199, 200, 157, 159, 196, 198, 158, 160, 195, 197} = NS Using Progression 1;

//----------
// Surfaces
//----------

// From bottom up; in z direction
Transfinite Surface {49} = {45, 50, 52, 48};
Transfinite Surface {22} = {32, 29, 34, 36};
Transfinite Surface {59} = {48, 52, 60, 58};
Transfinite Surface {35} = {36, 34, 39, 40};
Transfinite Surface {69} = {58, 60, 68, 67};
Transfinite Surface {43} = {40, 39, 43, 44};
Transfinite Surface {47} = {46, 49, 51, 47};
Transfinite Surface {20} = {31, 30, 33, 35};
Transfinite Surface {31} = {35, 33, 37, 38};
Transfinite Surface {63} = {57, 59, 64, 61};
Transfinite Surface {39} = {38, 37, 41, 42};
Transfinite Surface {53} = {53, 55, 56, 54};
Transfinite Surface {93} = {55, 9, 13, 56};
Transfinite Surface {9} = {9, 6, 14, 13};
Transfinite Surface {126} = {54, 56, 66, 63};
Transfinite Surface {29} = {13, 14, 20, 22};
Transfinite Surface {65} = {63, 66, 65, 62};
Transfinite Surface {119} = {66, 22, 28, 65};
// Transfinite Surface {119} = {66, 65, 28, 22};
Transfinite Surface {17} = {22, 20, 26, 28};
Transfinite Surface {57} = {47, 51, 59, 57};
Transfinite Surface {73} = {69, 71, 72, 70};
Transfinite Surface {92} = {71, 78, 80, 72};
Transfinite Surface {83} = {78, 77, 79, 80};
Transfinite Surface {100} = {85, 87, 88, 86};
Transfinite Surface {120} = {87, 93, 95, 88};
Transfinite Surface {110} = {93, 94, 96, 95};
Transfinite Surface {78} = {73, 75, 76, 74};
Transfinite Surface {96} = {75, 82, 84, 76};
Transfinite Surface {88} = {82, 81, 83, 84};
Transfinite Surface {105} = {89, 91, 92, 90};
Transfinite Surface {123} = {91, 97, 99, 92};
Transfinite Surface {116} = {97, 98, 100, 99};

// In y direction
Transfinite Surface {48} = {48, 52, 51, 47};
Transfinite Surface {23} = {36, 34, 33, 35};
Transfinite Surface {54} = {47, 51, 56, 54};
Transfinite Surface {27} = {35, 33, 14, 13};
Transfinite Surface {74} = {54, 56, 72, 70};
Transfinite Surface {85} = {13, 14, 79, 80};
Transfinite Surface {79} = {70, 72, 76, 74};
Transfinite Surface {97} = {72, 80, 84, 76};
Transfinite Surface {94} = {56, 13, 80, 72};
Transfinite Surface {90} = {84, 80, 79, 83};
Transfinite Surface {58} = {58, 60, 59, 57};
Transfinite Surface {34} = {40, 39, 37, 38};
Transfinite Surface {30} = {38, 37, 20, 22};
Transfinite Surface {64} = {57, 59, 66, 63};
Transfinite Surface {99} = {63, 66, 87, 85};
Transfinite Surface {118} = {66, 22, 93, 87};
Transfinite Surface {108} = {22, 20, 94, 93};
Transfinite Surface {104} = {85, 87, 91, 89};
Transfinite Surface {113} = {93, 94, 98, 97};
Transfinite Surface {68} = {67, 68, 64, 61};
Transfinite Surface {42} = {44, 43, 41, 42};
Transfinite Surface {112} = {28, 26, 96, 95};
Transfinite Surface {38} = {42, 41, 26, 28};
Transfinite Surface {62} = {61, 64, 65, 62};
Transfinite Surface {101} = {62, 65, 88, 86};
Transfinite Surface {121} = {65, 28, 95, 88};
Transfinite Surface {106} = {86, 88, 92, 90};
Transfinite Surface {124} = {88, 95, 99, 92};
Transfinite Surface {117} = {95, 96, 100, 99};
Transfinite Surface {46} = {45, 50, 49, 46};
Transfinite Surface {18} = {32, 29, 30, 31};
Transfinite Surface {24} = {31, 30, 6, 9};
Transfinite Surface {52} = {46, 49, 55, 53};
Transfinite Surface {72} = {53, 55, 71, 69};
Transfinite Surface {81} = {9, 6, 77, 78};
Transfinite Surface {91} = {55, 9, 78, 71};
Transfinite Surface {77} = {69, 71, 75, 73};
Transfinite Surface {95} = {71, 78, 82, 75};
Transfinite Surface {86} = {78, 77, 81, 82};
Transfinite Surface {122} = {87, 93, 97, 91};

// In x direction
Transfinite Surface {67} = {67, 58, 57, 61};
Transfinite Surface {56} = {58, 48, 47, 57};
Transfinite Surface {45} = {48, 45, 46, 47};
Transfinite Surface {61} = {61, 57, 63, 62};
Transfinite Surface {125} = {57, 47, 54, 63};
Transfinite Surface {51} = {47, 46, 53, 54};
Transfinite Surface {98} = {62, 63, 85, 86};
Transfinite Surface {71} = {54, 53, 69, 70};
Transfinite Surface {76} = {70, 69, 73, 74};
Transfinite Surface {103} = {86, 85, 89, 90};
Transfinite Surface {50} = {52, 50, 49, 51};
Transfinite Surface {66} = {64, 59, 66, 65};
Transfinite Surface {127} = {59, 51, 56, 66};
Transfinite Surface {55} = {51, 49, 55, 56};
Transfinite Surface {44} = {44, 40, 38, 42};
Transfinite Surface {36} = {40, 36, 35, 38};
Transfinite Surface {21} = {36, 32, 31, 35};
Transfinite Surface {40} = {28, 22, 38, 42};
Transfinite Surface {32} = {38, 35, 13, 22};
Transfinite Surface {26} = {13, 9, 31, 35};
Transfinite Surface {109} = {28, 22, 93, 95};
Transfinite Surface {84} = {13, 9, 78, 80};
Transfinite Surface {114} = {95, 93, 97, 99};
Transfinite Surface {89} = {80, 78, 82, 84};
Transfinite Surface {41} = {43, 39, 37, 41};
Transfinite Surface {33} = {39, 34, 33, 37};
Transfinite Surface {37} = {41, 37, 20, 26};
Transfinite Surface {28} = {37, 33, 14, 20};
Transfinite Surface {25} = {33, 30, 6, 14};
Transfinite Surface {19} = {34, 29, 30, 33};
Transfinite Surface {111} = {26, 20, 94, 96};
Transfinite Surface {115} = {96, 94, 98, 100};
Transfinite Surface {87} = {79, 77, 81, 83};
Transfinite Surface {82} = {14, 6, 77, 79};
Transfinite Surface {60} = {51, 59, 60, 52};
Transfinite Surface {70} = {59, 64, 68, 60};
Transfinite Surface {102} = {65, 66, 87, 88};
Transfinite Surface {107} = {88, 87, 91, 92};
Transfinite Surface {75} = {56, 55, 71, 72};
Transfinite Surface {80} = {72, 71, 75, 76};

Recombine Surface "*";

//---------
// Volumes
//---------

Transfinite Volume{25} = {45, 50, 52, 48, 46, 49, 51, 47};
Transfinite Volume{26} = {46, 49, 51, 47, 53, 55, 56, 54};
Transfinite Volume{31} = {53, 55, 56, 54, 69, 71, 72, 70};
Transfinite Volume{32} = {69, 71, 72, 70, 73, 75, 76, 74};
Transfinite Volume{35} = {55, 9, 13, 56, 71, 78, 80, 72};
Transfinite Volume{36} = {71, 78, 80, 72, 75, 82, 84, 76};
Transfinite Volume{19} = {32, 29, 34, 36, 31, 30, 33, 35};
Transfinite Volume{20} = {31, 30, 33, 35, 9, 6, 14, 13};
Transfinite Volume{33} = {9, 6, 14, 13, 78, 77, 79, 80};
Transfinite Volume{34} = {78, 77, 79, 80, 82, 81, 83, 84};
Transfinite Volume{28} = {48, 52, 60, 58, 47, 51, 59, 57};
Transfinite Volume{43} = {47, 51, 59, 57, 54, 56, 66, 63};
Transfinite Volume{22} = {36, 34, 39, 40, 35, 33, 37, 38};
Transfinite Volume{21} = {35, 33, 37, 38, 13, 14, 20, 22};
Transfinite Volume{24} = {44, 40, 39, 43, 42, 38, 37, 41};
Transfinite Volume{23} = {42, 38, 37, 41, 28, 22, 20, 26};
Transfinite Volume{39} = {28, 22, 20, 26, 95, 93, 94, 96};
Transfinite Volume{40} = {95, 93, 94, 96, 99, 97, 98, 100};
Transfinite Volume{30} = {67, 58, 60, 68, 61, 57, 59, 64};
Transfinite Volume{29} = {61, 57, 59, 64, 62, 63, 66, 65};
Transfinite Volume{37} = {62, 63, 66, 65, 86, 85, 87, 88};
Transfinite Volume{38} = {86, 85, 87, 88, 90, 89, 91, 92};
Transfinite Volume{41} = {66, 22, 28, 65, 87, 93, 95, 88};
Transfinite Volume{42} = {88, 87, 93, 95, 92, 91, 97, 99};

//---------------------------------------
//
// Define boundary and volume conditions
//
//---------------------------------------
TINY = 1.0e-3;
HUGE = 1.0e+3;

Physical Surface("periodic_x")
  = {Surface In BoundingBox{ -TINY, -HUGE, -HUGE,
                             +TINY, +HUGE, +HUGE},
     Surface In BoundingBox{L-TINY, -HUGE, -HUGE,
                            L+TINY, +HUGE, +HUGE}};

Physical Surface("periodic_y")
  = {Surface In BoundingBox{-HUGE,  -TINY, -HUGE,
                            +HUGE,  +TINY, +HUGE},
     Surface In BoundingBox{-HUGE, L-TINY, -HUGE,
                            +HUGE, L+TINY, +HUGE}};

Physical Surface("bottom")
  = {Surface In BoundingBox{-HUGE, -HUGE, -TINY,
                            +HUGE, +HUGE, +TINY}};

Physical Surface("top")
  = {Surface In BoundingBox{-HUGE, -HUGE, H-TINY,
                            +HUGE, +HUGE, H+TINY}};

Physical Surface("spacer_walls") = {119, 93, 29, 126, 21, 36, 44, 
                                    26, 32, 40, 50, 60, 70, 55, 
                                    127, 66, 74, 94, 85, 79, 97, 
                                    90, 99, 118, 108, 104, 122, 113};

Physical Volume("fluid") = {25, 26, 31, 32, 35, 36, 19, 20, 33, 34, 28, 43, 
                            22, 21, 30, 29, 37, 38, 41, 42, 24, 23, 39, 40};

