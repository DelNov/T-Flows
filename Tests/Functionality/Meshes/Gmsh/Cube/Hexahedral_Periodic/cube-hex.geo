Point(1) = {0.0, 0.0, 0.0}; 
Point(2) = {1.0, 0.0, 0.0};
Point(3) = {0.0, 1.0, 0.0};
Point(4) = {1.0, 1.0, 0.0};
Point(5) = {0.0, 0.0, 1.0};
Point(6) = {1.0, 0.0, 1.0};
Point(7) = {0.0, 1.0, 1.0};
Point(8) = {1.0, 1.0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {4, 3};
Line(4) = {3, 1};
Line(5) = {5, 6};
Line(6) = {6, 8};
Line(7) = {8, 7};
Line(8) = {7, 5};
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

Line Loop(13) = { 1,  2,   3,  4};  Plane Surface(14) = {13};
Line Loop(15) = { 5,  6,   7,  8};  Plane Surface(16) = {15};
Line Loop(17) = { 9,  5, -10, -1};  Plane Surface(18) = {17};
Line Loop(19) = {11, -7, -12,  3};  Plane Surface(20) = {19};
Line Loop(21) = { 9, -8, -11,  4};  Plane Surface(22) = {21};
Line Loop(23) = {10,  6, -12, -2};  Plane Surface(24) = {23};

Surface Loop(25) = {14, 22, 20, 18, 16, 24};
Volume(26) = {25};

Physical Surface("z_per") = {14, 16};
Physical Surface("y_per") = {18, 20};
Physical Surface("x_per") = {22, 24};
Physical Volume("cube") = {26};

Transfinite Line {1, 3, 5, 7} = 11;
Transfinite Line {9, 10, 11, 12} = 11;
Transfinite Line {2, 4, 8, 6} = 11;
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";
