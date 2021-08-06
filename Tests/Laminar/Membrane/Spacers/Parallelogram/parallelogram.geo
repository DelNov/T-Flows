//------------------------------------------------
//
// Characteristic dimensions and grid resolutions
//
// (See figure parallelogram.fig)
//
//------------------------------------------------
L =  5;    // length of the computational domain
H =  1.5;  // height of the computational domain
S =  0.6;  // spacer width
T =  0.2;  // thinner spacer dimension

NH_HALF = 20;  // number of cells along half-height
NS      = 12;  // number of cells covering spacers
NB_HALF = 40;  // number of cells between spacers

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

//---------------------------------
// Box defining the overall domain
//---------------------------------
Box(1) = {0, 0, 0,   L, L, H};

// Cut the domin with rectangle in two parts
Rectangle(7) = {0, 0, H/2, L, L, 0};
BooleanFragments{Volume{1}; Surface{7}; Delete; }{}

//-----------------------------
// Paralleilograms for spacers
//-----------------------------
X1 = L/2 - S/2 - T/2;
X2 = X1 + S;
X3 = X1 + S + T;
X4 = X1 + T;
Point(21) = {X1,  -T,  0,    1.0};
Point(22) = {X2,  -T,  0,    1.0};
Point(23) = {X3,  -T,  H/2,  1.0};
Point(24) = {X4,  -T,  H/2,  1.0};

Line(33) = {21, 22};  Line(34) = {22, 23};
Line(35) = {23, 24};  Line(36) = {24, 21};
Curve Loop(18) = {34, 35, 36, 33};  Plane Surface(18) = {18};

Point(25) = {-T,  X4,  H/2,  1.0};
Point(26) = {-T,  X3,  H/2,  1.0};
Point(27) = {-T,  X2,  H,    1.0};
Point(28) = {-T,  X1,  H,    1.0};

Line(37) = {25, 26};  Line(38) = {26, 27};
Line(39) = {27, 28};  Line(40) = {28, 25};
Curve Loop(19) = {37, 38, 39, 40};  Plane Surface(19) = {19};

// Create volumes from those parallelograms ...
Extrude {0, L+2*T, 0} {Surface{18}; Layers{5}; Recombine;}
Extrude {L+2*T, 0, 0} {Surface{19}; Layers{5}; Recombine;}

// ... and unite them
BooleanUnion{Volume{3}; Delete;}{Volume{4}; Delete;}

// Cut the spacers from the domain
BooleanDifference{Volume{2}; Delete;}{ Volume{3};}
BooleanDifference{Volume{1}; Delete;}{ Volume{3};}

// Delete spacers - they are void, not part of a simulation
Recursive Delete {Volume{3};}

//---------------------------------------------------
// Cut domains further to be able to make a hex mesh
//---------------------------------------------------

Translate {0, 0,-H/2} {Duplicata{Point{48}; Point{53}; Point{43};Point{50};}}
Translate {0, 0, H/2} {Duplicata{Point{46}; Point{45}; Point{55};Point{56};}}

Line(133) = {48, 73}; Line(134) = {53, 74};
Line(135) = {43, 75}; Line(136) = {50, 76};
Line(137) = {46, 77}; Line(138) = {45, 78};
Line(139) = {55, 79}; Line(140) = {56, 80};

Line(141) = {73, 75}; Line(142) = {74, 76};
Line(143) = {77, 79}; Line(144) = {78, 80};

Curve Loop(64) = {135, -141, -133, -77, -44, -76};  Plane Surface(64) = {64};
Curve Loop(65) = {136, -142, -134, -88, -48, -89};  Plane Surface(65) = {65};
Curve Loop(66) = {137, 143, -139, -98, -35, -81};  Plane Surface(66) = {66};
Curve Loop(67) = {144, -140, -96, -53, -79, 138};  Plane Surface(67) = {67};

BooleanFragments{ Volume{4}; Delete; }{ Surface{67}; Delete; }
BooleanFragments{ Volume{9}; Delete; }{ Surface{66}; Delete; }
BooleanFragments{ Volume{5}; Delete; }{ Surface{79}; Delete; }
BooleanFragments{ Volume{12}; Delete; }{ Surface{90}; Delete; }
BooleanFragments{ Volume{7}; Delete; }{ Surface{64}; Delete; }
BooleanFragments{ Volume{15}; Delete; }{ Surface{65}; Delete; }
BooleanFragments{ Volume{6}; Delete; }{ Surface{127}; Delete; }
BooleanFragments{ Volume{18}; Delete; }{ Surface{138}; Delete; }

//-------------------------------------------------------------------------
//
// Super important - call coherence to avoid duplicate geometrical entries
//
//-------------------------------------------------------------------------
Coherence;

//-------------
//
// Define mesh
//
//-------------

//--------
// Curves
//--------

// X direction
Transfinite Curve {367, 368, 373, 301, 378, 305, 326, 300, 329, 303, 325, 327} = NB_HALF+1 Using Progression 1;
Transfinite Curve {348, 352, 356, 310, 361, 317, 340, 308, 345, 315, 338, 341} = NB_HALF+1 Using Progression 1;
Transfinite Curve {279, 295, 319, 44, 48, 318, 337, 320, 332, 333} = NS+1 Using Progression 1;

// Y direction
Transfinite Curve {350, 351, 369, 313, 365, 316, 306, 311, 299, 314, 304, 297} = NB_HALF+1 Using Progression 1;
Transfinite Curve {359, 362, 379, 344, 376, 335, 331, 342, 324, 334, 328, 322} = NB_HALF+1 Using Progression 1;
Transfinite Curve {355, 357, 374, 353, 372, 35, 53, 370, 215, 200} = NS+1 Using Progression 1;

// Z direction
Transfinite Curve {307, 312, 339, 309, 302, 343, 175, 202, 153, 186, 296, 336, 
                   330, 298, 321, 323} = NH_HALF+1 Using Progression 1;
Transfinite Curve {346, 349, 354, 347, 366, 358, 225, 244, 266, 284, 363, 360, 
                   377, 364, 371, 375} = NH_HALF+1 Using Progression 1;

//----------
// Surfaces
//----------
Recursive Delete {Surface{103}; Surface{114}; Surface{151}; Surface{162}; }

// No particular order here :-(
Transfinite Surface {212} = {161, 163, 127, 162};
Transfinite Surface {198} = {156, 155, 157, 109};
Transfinite Surface {217} = {162, 127, 130, 164};
Transfinite Surface {203} = {109, 157, 158, 117};
Transfinite Surface {221} = {164, 130, 166, 165};
Transfinite Surface {207} = {117, 158, 159, 160};
Transfinite Surface {167} = {131, 136, 30, 134};
Transfinite Surface {177} = {136, 140, 24, 30};
Transfinite Surface {215} = {134, 30, 33, 144};
Transfinite Surface {162} = {130, 117, 25, 33};
Transfinite Surface {201} = {24, 142, 152, 25};
Transfinite Surface {183} = {144, 33, 148, 146};
Transfinite Surface {188} = {148, 33, 25, 150};
Transfinite Surface {193} = {150, 25, 152, 154};
Transfinite Surface {165} = {133, 132, 135, 84};
Transfinite Surface {176} = {84, 135, 139, 92};
Transfinite Surface {171} = {139, 138, 141, 92};
Transfinite Surface {181} = {143, 96, 147, 145};
Transfinite Surface {186} = {147, 96, 99, 149};
Transfinite Surface {191} = {149, 99, 151, 153};
Transfinite Surface {210} = {161, 163, 136, 131};
Transfinite Surface {195} = {156, 155, 137, 140};
Transfinite Surface {211} = {162, 127, 30, 134};
Transfinite Surface {151} = {30, 127, 109, 24};
Transfinite Surface {169} = {140, 137, 138, 139};
Transfinite Surface {164} = {131, 136, 135, 132};
Transfinite Surface {175} = {136, 140, 139, 135};
Transfinite Surface {220} = {165, 146, 148, 166};
Transfinite Surface {206} = {160, 150, 154, 159};
Transfinite Surface {182} = {146, 145, 147, 148};
Transfinite Surface {189} = {147, 149, 150, 148};
Transfinite Surface {194} = {149, 153, 154, 150};
Transfinite Surface {180} = {144, 33, 96, 143};
Transfinite Surface {190} = {25, 152, 151, 99};
Transfinite Surface {185} = {96, 33, 25, 99};
Transfinite Surface {173} = {140, 137, 142, 24};
Transfinite Surface {196} = {155, 157, 142, 137};
Transfinite Surface {170} = {137, 142, 141, 138};
Transfinite Surface {200} = {157, 158, 152, 142};
Transfinite Surface {205} = {158, 159, 154, 152};
Transfinite Surface {192} = {152, 154, 153, 151};
Transfinite Surface {197} = {156, 109, 24, 140};
Transfinite Surface {204} = {109, 117, 25, 24};
Transfinite Surface {208} = {117, 160, 150, 25};
Transfinite Surface {172} = {140, 24, 92, 139};
Transfinite Surface {114} = {24, 25, 99, 92};
Transfinite Surface {187} = {149, 99, 25, 150};
Transfinite Surface {168} = {135, 136, 30, 84};
Transfinite Surface {213} = {136, 163, 127, 30};
Transfinite Surface {219} = {165, 164, 144, 146};
Transfinite Surface {214} = {164, 162, 134, 144};
Transfinite Surface {209} = {162, 161, 131, 134};
Transfinite Surface {179} = {146, 144, 143, 145};
Transfinite Surface {163} = {134, 131, 132, 133};
Transfinite Surface {222} = {148, 33, 130, 166};
Transfinite Surface {218} = {33, 30, 127, 130};
Transfinite Surface {166} = {133, 84, 30, 134};
Transfinite Surface {178} = {84, 92, 24, 30};
Transfinite Surface {174} = {92, 141, 142, 24};
Transfinite Surface {199} = {157, 109, 24, 142};
Transfinite Surface {202} = {158, 117, 25, 152};
Transfinite Surface {184} = {96, 147, 148, 33};
Transfinite Surface {216} = {164, 130, 33, 144};

Recombine Surface "*";

//---------
// Volumes
//---------
Transfinite Volume{17} = {161, 163, 127, 162, 131, 136, 30, 134};
Transfinite Volume{14} = {156, 155, 157, 109, 140, 137, 142, 24};
Transfinite Volume{18} = {162, 127, 130, 164, 134, 30, 33, 144};
Transfinite Volume{15} = {109, 157, 158, 117, 24, 142, 152, 25};
Transfinite Volume{19} = {164, 130, 166, 165, 144, 33, 148, 146};
Transfinite Volume{16} = {117, 158, 159, 160, 25, 152, 154, 150};
Transfinite Volume{8} = {131, 136, 30, 134, 132, 135, 84, 133};
Transfinite Volume{10} = {136, 140, 24, 30, 135, 139, 92, 84};
Transfinite Volume{9} = {140, 137, 142, 24, 139, 138, 141, 92};
Transfinite Volume{11} = {144, 33, 148, 146, 143, 96, 147, 145};
Transfinite Volume{12} = {33, 25, 150, 148, 96, 99, 149, 147};
Transfinite Volume{13} = {25, 152, 154, 150, 99, 151, 153, 149};

//---------------------------------------
//
// Define boundary and volume conditions
//
//---------------------------------------
TINY = 1.0e-6;
HUGE = 1.0e+6;

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


Physical Surface("spacer_walls") = {201, 190, 174, 180, 166, 215, 185, 178, 
                                    197, 204, 208, 222, 218, 213, 177, 188};

Physical Volume("fluid") = {9, 10, 8, 14, 17, 15, 18, 13, 12, 11, 19, 16};
