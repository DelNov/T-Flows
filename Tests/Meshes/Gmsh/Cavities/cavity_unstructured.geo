 Point(1) = {0.0, 0.0, 0.0}; 
 Point(2) = {1.0, 0.0, 0.0};
 Point(3) = {1.0, 0.0, 1.0};
 Point(4) = {0.0, 0.0, 1.0};
 Line(1) = {1, 2};
 Line(2) = {2, 3};
 Line(3) = {3, 4};
 Line(4) = {4, 1};
 Line Loop(13) = {1, 2, 3, 4};
 Plane Surface(14) = {13};
 Transfinite Line {1, 2, 3, 4} = 21;
 Transfinite Volume "*";
 //+
 Extrude {0, 0.1, 0} {
   Surface{14}; 
   Layers{6};
   Recombine;
 }
 Transfinite Line {16, 17, 18, 19} = 21;
 Transfinite Line {21, 22, 26, 30} = 6;
 Transfinite Surface {35, 31, 27, 23};
//+
Physical Surface("inlet") = {35};
//+
Physical Surface("outlet") = {27};
//+
Physical Surface("wall") = {31, 23};
//+
Physical Surface("periodic") = {14, 36};
//+
Physical Volume("fluid") = {1};
