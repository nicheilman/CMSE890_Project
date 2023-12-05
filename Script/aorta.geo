

Point(1) = {2, 0, 0};
Point(2) = {2, 0, 2};
Point(3) = {3, 0, 2.5};
Point(4) = {4, 0, 3};
Point(5) = {4, 0, 5};
Point(6) = {3, 0, 5.5};
Point(7) = {2, 0, 6};
Point(8) = {2, 0, 8};
Point(9) = {0, 2, 0};
Point(10) = {0, 2, 8};
Point(11) = {0, 0, 0};
Point(12) = {0, 0, 8};

Line(101) = {1, 2};
Line(102) = {4, 5};
Line(103) = {7, 8};
BSpline(104) = {2, 3, 4};
BSpline(105) = {5, 6, 7};
Circle(106) = {9, 11, 1};
Circle(107) = {8, 12, 10};

Physical Curve("test") = {101, 102, 103, 104, 105, 106, 107};

lineVect1[] = Extrude{{0, 0, 1}, {0, 0, 0}, Pi/2}{101, 102, 103, 104, 105, 106, 107};


