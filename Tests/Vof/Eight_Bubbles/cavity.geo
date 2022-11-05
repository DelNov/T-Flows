SetFactory("OpenCASCADE");

Box(1) = {-0.5, -0.5, -0.5, 1, 1, 1};

MeshSize {1:8} = 0.02;
//+
Physical Surface("WALL", 13) = {1, 2, 3, 4, 5, 6};
//+
Physical Volume("FLUID", 14) = {1};
