h = 20.0;

L = 9.6 * h;
W = 4.8 * h;
H = 1.0 * h;

SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, L, W, H};
