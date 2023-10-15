//--------------------------------
// Definition of buildings begins
//--------------------------------
n_obstacles = 0;    // number of obstacles

// Obstacle with ID: 1
n_obstacles++;
o         =    n_obstacles;
node_b(o) =    4;
xo(MAXN*o+1) = -A/2-A-P;  yo(MAXN*o+1) = -A/2;  d(MAXN*o+1) = DELTA_BUILDING;
xo(MAXN*o+2) = -A/2-P;    yo(MAXN*o+2) = -A/2;  d(MAXN*o+2) = DELTA_BUILDING;
xo(MAXN*o+3) = -A/2-P;    yo(MAXN*o+3) =  A/2;  d(MAXN*o+3) = DELTA_BUILDING;
xo(MAXN*o+4) = -A/2-A-P;  yo(MAXN*o+4) =  A/2;  d(MAXN*o+4) = DELTA_BUILDING;
height_o(o)  =  A;
type_o(o)    = SOLID;

// Obstacle with ID: 2
n_obstacles++;
o         =    n_obstacles;
node_b(o) =    4;
xo(MAXN*o+1) = -A/2;      yo(MAXN*o+1) = -A/2;  d(MAXN*o+1) = DELTA_BUILDING;
xo(MAXN*o+2) =  A/2;      yo(MAXN*o+2) = -A/2;  d(MAXN*o+2) = DELTA_BUILDING;
xo(MAXN*o+3) =  A/2;      yo(MAXN*o+3) =  A/2;  d(MAXN*o+3) = DELTA_BUILDING;
xo(MAXN*o+4) = -A/2;      yo(MAXN*o+4) =  A/2;  d(MAXN*o+4) = DELTA_BUILDING;
height_o(o)  =  A*2;
type_o(o)    = SOLID;

// Obstacle with ID: 3
n_obstacles++;
o         =    n_obstacles;
node_b(o) =    4;
xo(MAXN*o+1) =  A/2+P;    yo(MAXN*o+1) = -A/2;  d(MAXN*o+1) = DELTA_BUILDING;
xo(MAXN*o+2) =  A/2+A+P;  yo(MAXN*o+2) = -A/2;  d(MAXN*o+2) = DELTA_BUILDING;
xo(MAXN*o+3) =  A/2+A+P;  yo(MAXN*o+3) =  A/2;  d(MAXN*o+3) = DELTA_BUILDING;
xo(MAXN*o+4) =  A/2+P;    yo(MAXN*o+4) =  A/2;  d(MAXN*o+4) = DELTA_BUILDING;
height_o(o)  =  A*3;
type_o(o)    = SOLID;

// Obstacle with ID: 4
n_obstacles++;
o         =    n_obstacles;
node_b(o) =    4;
xo(MAXN*o+1) = -A*8;  yo(MAXN*o+1) = -A*4;  d(MAXN*o+1) = DELTA_OBSTACLE;
xo(MAXN*o+2) = -A*4;  yo(MAXN*o+2) = -A*4;  d(MAXN*o+2) = DELTA_OBSTACLE;
xo(MAXN*o+3) = -A*4;  yo(MAXN*o+3) =  A*4;  d(MAXN*o+3) = DELTA_OBSTACLE;
xo(MAXN*o+4) = -A*8;  yo(MAXN*o+4) =  A*4;  d(MAXN*o+4) = DELTA_OBSTACLE;
height_o(o)  =  A/2;
type_o(o)    = POROUS;

// Obstacle with ID: 5
n_obstacles++;
o         =    n_obstacles;
node_b(o) =    8;
xo(MAXN*o+1) =  A* 7;  yo(MAXN*o+1) = -A*3;  d(MAXN*o+1) = DELTA_OBSTACLE;
xo(MAXN*o+2) =  A* 9;  yo(MAXN*o+2) = -A*5;  d(MAXN*o+2) = DELTA_OBSTACLE;
xo(MAXN*o+3) =  A*15;  yo(MAXN*o+3) = -A*5;  d(MAXN*o+3) = DELTA_OBSTACLE;
xo(MAXN*o+4) =  A*17;  yo(MAXN*o+4) = -A*3;  d(MAXN*o+4) = DELTA_OBSTACLE;
xo(MAXN*o+5) =  A*17;  yo(MAXN*o+5) =  A*3;  d(MAXN*o+5) = DELTA_OBSTACLE;
xo(MAXN*o+6) =  A*15;  yo(MAXN*o+6) =  A*5;  d(MAXN*o+6) = DELTA_OBSTACLE;
xo(MAXN*o+7) =  A* 9;  yo(MAXN*o+7) =  A*5;  d(MAXN*o+7) = DELTA_OBSTACLE;
xo(MAXN*o+8) =  A* 7;  yo(MAXN*o+8) =  A*3;  d(MAXN*o+8) = DELTA_OBSTACLE;
height_o(o)  =  A/2;
type_o(o)    = POROUS;

