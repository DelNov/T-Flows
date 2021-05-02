A =  2.0;
N = 20;
gridsize = A / N;

// All numbering counterclockwise from bottom-left corner
Point(1) = {-A/2, -A/2,  0,  gridsize};
Point(2) = { A/2, -A/2,  0,  gridsize};
Point(3) = { A/2,  A/2,  0,  gridsize};
Point(4) = {-A/2,  A/2,  0,  gridsize};

Line(1) = {1, 2};                // bottom line
Line(2) = {2, 3};                // right line
Line(3) = {3, 4};                // top line
Line(4) = {4, 1};                // left line
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

surfaceVector[] = Extrude {0, 0, A} {
 Surface{6};
 Layers{N};

 Recombine;
};
Transfinite Surface {6};
Recombine Surface {6};

Physical Surface("bottom") = {15};
Physical Surface("top") = {23};
Physical Surface("side") = {19, 28, 27, 6};
Physical Volume("volume") = {1};
