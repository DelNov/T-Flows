    Depth    = 0.5;           // x dimension
    Width    = 2.0;           // y dimension
    Height   = Depth * 2.0;   // z dimension

    N = 32;

    gridsize = Depth / N;    // Mijail used 64.0

    // All numbering counterclockwise from bottom-left corner
    Point(1) = {0.0,   0.0, 0.0,    gridsize};
    Point(2) = {Depth, 0.0, 0.0,    gridsize};
    Point(3) = {Depth, 0.0, Width, gridsize};
    Point(4) = {0.0,   0.0, Width, gridsize};
    Line(1) = {1, 2};				// bottom line
    Line(2) = {2, 3};				// right line
    Line(3) = {3, 4};				// top line
    Line(4) = {4, 1};				// left line
    Line Loop(5) = {1, 2, 3, 4};
    // the order of lines in Line Loop is used again in surfaceVector[]
    Plane Surface(6) = {5};

    surfaceVector[] = Extrude {0.0, Height, 0.0} {
     Surface{6};
     Layers{N * 2};

     Recombine;
    };
    Transfinite Surface {6};
    Recombine Surface {6};
    /* surfaceVector contains in the following order:
    [0] - front surface (opposed to source surface)
    [1] - extruded volume
    [2] - bottom surface (belonging to 1st line in "Line Loop (6)")
    [3] - right surface (belonging to 2nd line in "Line Loop (6)")
    [4] - top surface (belonging to 3rd line in "Line Loop (6)")
    [5] - left surface (belonging to 4th line in "Line Loop (6)") */
    Physical Volume("internal") = surfaceVector[1];
    Physical Surface("wall") = {surfaceVector[2], surfaceVector[4], 6, surfaceVector[5], surfaceVector[0]};
    Physical Surface("symmetry") = {surfaceVector[3]};
