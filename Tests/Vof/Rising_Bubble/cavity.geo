
SetFactory("OpenCASCADE");
Box(1) = {-0.5, -0.5, -0.5, 1, 1, 2};

MeshSize {1:8} = 0.015625;

Physical Surface("wall", 13) = {1:6};
Physical Volume("fluid", 14) = {1};
