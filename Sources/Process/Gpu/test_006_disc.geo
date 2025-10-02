SetFactory("OpenCASCADE");
Circle(1) = {0, -0, 0, 0.5, 0, 2*Pi};

Curve Loop(1) = {1};

MeshSize {1} = 0.03;

Plane Surface(1) = {1};
Recombine Surface{1};

Extrude {0, 0, 0.1} {
  Surface{1}; Layers {5}; Recombine;
}

Physical Surface("top_and_bottom", 4) = {1, 3};
Physical Surface("perimeter", 5) = {2};
Physical Volume("disc", 6) = {1};
