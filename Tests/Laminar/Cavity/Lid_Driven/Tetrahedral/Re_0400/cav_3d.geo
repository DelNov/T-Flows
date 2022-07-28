SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 1, 0.2, 1};

MeshSize {1, 2, 3, 4, 5, 6, 7, 8} = 0.025;

Physical Surface("symmetry", 13) = {3, 4};
Physical Surface("lid", 14) = {6};
Physical Surface("wall", 15) = {1, 5, 2};
Physical Volume("fluid", 16) = {1};
