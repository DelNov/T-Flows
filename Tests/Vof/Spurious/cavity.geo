    Width = 2.0; 
    Height = 2.0;
    R = 0.4;
    R_dx = 12.8;
    gridsize = R / R_dx;
    meshThickness = gridsize * 3.0; 
 
    // All numbering counterclockwise from bottom-left corner
    Point(1) = {0.0, 0.0, 0.0, gridsize};
    Point(2) = {Width, 0.0, 0.0, gridsize};
    Point(3) = {Width, 0.0, Height, gridsize};
    Point(4) = {0.0, 0.0, Height, gridsize};
    Line(1) = {1, 2};				// bottom line
    Line(2) = {2, 3};				// right line
    Line(3) = {3, 4};				// top line
    Line(4) = {4, 1};				// left line
    Line Loop(5) = {1, 2, 3, 4}; 	
    // the order of lines in Line Loop is used again in surfaceVector[]
    Plane Surface(6) = {5};
    
    surfaceVector[] = Extrude {0.0, meshThickness, 0.0} {
     Surface{6};
     Layers{3};

     Recombine;
    };
    Transfinite Surface {6};
    Recombine Surface {6};
    /* surfaceVector contains in the following order:
    [0]	- front surface (opposed to source surface)
    [1] - extruded volume
    [2] - bottom surface (belonging to 1st line in "Line Loop (6)")
    [3] - right surface (belonging to 2nd line in "Line Loop (6)")
    [4] - top surface (belonging to 3rd line in "Line Loop (6)")
    [5] - left surface (belonging to 4th line in "Line Loop (6)") */
    Physical Volume("internal") = surfaceVector[1];
    Physical Surface("wall") = {surfaceVector[4], surfaceVector[3], surfaceVector[5], surfaceVector[2]};
    Physical Surface("periodic") = {surfaceVector[0], 6}; // from Plane Surface (6) 
