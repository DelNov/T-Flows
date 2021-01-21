//----------------------------------
// Parameters defining the geometry
//----------------------------------
WIDTH  = 1.0;
HEIGHT = 2.0;
DELTA  = WIDTH / 80.0; // Controls grid size (WIDTH / 40.0 ,
                       //                     WIDTH / 80.0 ,
                       //                     WIDTH / 160.0,
                       //                     WIDTH / 320.0)
LENGTH = DELTA * 3.0;

//--------
// Points
//--------

// All numbering counterclockwise from bottom-left corner
Point(1) = {0.0,   0.0, 0.0,    DELTA};
Point(2) = {WIDTH, 0.0, 0.0,    DELTA};
Point(3) = {WIDTH, 0.0, HEIGHT, DELTA};
Point(4) = {0.0,   0.0, HEIGHT, DELTA};

//-------
// Lines
//-------

Line(1) = {1, 2};        // bottom line
Line(2) = {2, 3};        // right line
Line(3) = {3, 4};        // top line
Line(4) = {4, 1};        // left line

//---------
// Surface
//---------
Line Loop(5) = {1, 2, 3, 4};
// the order of lines in Line Loop is used again in surfaceVector[]
Plane Surface(6) = {5};

//--------
// Volume
//--------
surfaceVector[] = Extrude {0.0, LENGTH, 0.0} {
 Surface{6};
 Layers{3};

 Recombine;
};
Transfinite Surface {6};
Recombine Surface {6};

//---------------------
// Boundary conditions
//---------------------

/* surfaceVector contains in the following order:
[0] - front surface (opposed to source surface)
[1] - extruded volume
[2] - bottom surface (belonging to 1st line in "Line Loop (6)")
[3] - right surface  (belonging to 2nd line in "Line Loop (6)")
[4] - top surface    (belonging to 3rd line in "Line Loop (6)")
[5] - left surface   (belonging to 4th line in "Line Loop (6)") */

Physical Volume ("internal") = surfaceVector[1];
Physical Surface("wall")     = {surfaceVector[2], surfaceVector[4]};
Physical Surface("symmetry") = {surfaceVector[3], surfaceVector[5]};
Physical Surface("periodic") = {surfaceVector[0], 6}; // from Plane Surface (6) 
