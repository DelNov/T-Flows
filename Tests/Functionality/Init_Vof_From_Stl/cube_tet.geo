
SetFactory("OpenCASCADE");
Box(1) = {-0.5, -0.5, -0.5, 1, 1, 1};

MeshSize {2, 4, 6, 8} = 0.05;
MeshSize {1, 3, 5, 7} = 0.1;

Physical Surface("all_walls", 13) = {1:6};
Physical Volume("inside", 14) = {1};
