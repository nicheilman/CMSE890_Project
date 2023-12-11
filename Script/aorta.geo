SetFactory("OpenCASCADE");

lc = 0.15;
r = 0.8;
L = 10; 

Point(1) = {r, 0, 0, lc};
Point(2) = {r, 0, 2, lc};
Point(3) = {3, 0, 2.5, lc};
Point(4) = {4, 0, 3, lc};
Point(5) = {4, 0, 5, lc};
Point(6) = {3, 0, 5.5, lc};
Point(7) = {r, 0, 6, lc};
Point(8) = {r, 0, L, lc};
Point(9) = {0, 2, 0, lc};
Point(10) = {4, 0, 4, lc};
Point(11) = {0, 0, 0, lc};
Point(12) = {0, 0, L, lc};

Line(101) = {1, 2};
Line(102) = {4, 5};
Line(103) = {7, 8};
BSpline(104) = {2, 3, 4, 10};
BSpline(105) = {10, 5, 6, 7};
Line(106) = {8, 12};
Line(107) = {12, 11}; 
Line(108) = {11, 1};

Curve Loop(110) = {101, 104, 105, 103, 106, 107, 108};
Plane Surface(201) = {110};

//Physical Curve("test") = {101, 102, 103, 104, 105, 106, 107};

lineVect1[] = Extrude{{0, 0, 1}, {0, 0, 0}, Pi/2}{Surface{201}; Recombine;};
lineVect2[] = Extrude{{0, 0, 1}, {0, 0, 0}, Pi/2}{Surface{208}; Recombine;};
lineVect3[] = Extrude{{0, 0, 1}, {0, 0, 0}, Pi/2}{Surface{215}; Recombine;};
lineVect4[] = Extrude{{0, 0, 1}, {0, 0, 0}, Pi/2}{Surface{222}; Recombine;};

Surface Loop(250) = {lineVect1[5], lineVect2[5], lineVect3[5], lineVect4[5]};
Surface Loop(251) = {lineVect1[7], lineVect2[7], lineVect3[7], lineVect4[7]};
Surface Loop(253) = {lineVect1[2], lineVect1[3], lineVect1[4], lineVect1[5], 
                     lineVect2[2], lineVect2[3], lineVect2[4], lineVect2[5], 
                     lineVect3[2], lineVect3[3], lineVect3[4], lineVect3[5], 
                     lineVect4[2], lineVect4[3], lineVect4[4], lineVect4[5]};

Physical Surface("Inlet") = {250};
Physical Surface("Outlet") = {251};
Physical Volume("Vol") = {lineVect1[1], lineVect2[1], lineVect3[1], lineVect4[1]};

Coherence;
