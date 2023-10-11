//--------------------------------
// Definition of buildings begins
//--------------------------------
n_buildings =  1;  // number of buildings

// Building with ID:      2
b         =    1;
node_b(b) =    4;
xo(MAXN*b+1) = -A/2;      yo(MAXN*b+1) = -A/2;  d(MAXN*b+1) = DELTA_BUILDING;
xo(MAXN*b+2) =  A/2;      yo(MAXN*b+2) = -A/2;  d(MAXN*b+2) = DELTA_BUILDING;
xo(MAXN*b+3) =  A/2;      yo(MAXN*b+3) =  A/2;  d(MAXN*b+3) = DELTA_BUILDING;
xo(MAXN*b+4) = -A/2;      yo(MAXN*b+4) =  A/2;  d(MAXN*b+4) = DELTA_BUILDING;
height_b(b)  =  A*2;
