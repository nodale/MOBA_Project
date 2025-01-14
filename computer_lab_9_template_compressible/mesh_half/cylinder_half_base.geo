// Cylinder center.
Point(1) = {0, 0, 0, 1.0};

r = 0.5;

// middle
h_cylinder = 0.2;
h_wake     = 0.4;
h_inner    = 0.5;
h_outer    = 2.0;

// Cylinder "egdes".
Point(2) = {r*Cos(45*Pi/180),  r*Sin(45*Pi/180),  0, h_cylinder};
Point(3) = {r*Cos(135*Pi/180), r*Sin(135*Pi/180), 0, h_cylinder};
Point(4) = {r*Cos(  0*Pi/180), r*Sin(  0*Pi/180), 0, h_cylinder};
Point(5) = {r*Cos(180*Pi/180), r*Sin(180*Pi/180), 0, h_cylinder};

// Cylinder lines.
Circle(1) = {4, 1, 2};
Circle(2) = {2, 1, 3};
Circle(3) = {3, 1, 5};

// Domain corners.
Point(6) = {-16, 15, 0, h_outer};
Point(7) = {30, 15, 0,  h_outer};

// Mesh structure points square around cylinder.
Point(10) = {-5, 5, 0, h_inner};
Point(13) = {-5, 0, 0, h_inner};

// Wake
Point(11) = {5, 5, 0, h_wake};
Point(12) = {5, 0, 0, h_wake};
Point(18) = {30, 0, 0, h_wake}; // Mesh structure points outlet.
Point(19) = {30, 5, 0, h_wake};

// Mesh structure points inlet.
Point(14) = {-16, 5, 0, h_outer};
Point(15) = {-16, 0, 0, h_outer};

// Mesh structure points top.
Point(16) = {-5, 15, 0, h_outer};
Point(17) = { 5, 15, 0, h_outer};

Line(4) = {15, 14};
Line(5) = {14, 6};
Line(6) = {6, 16};
Line(7) = {16, 17};
Line(8) = {17, 7};
Line(9) = {7, 19};
Line(10) = {19, 18};
Line(11) = {18, 12};
Line(12) = {12, 4};
Line(13) = {5, 13};
Line(14) = {13, 15};
Line(15) = {10, 14};
Line(16) = {10, 13};
Line(17) = {10, 16};
Line(18) = {10, 11};
Line(19) = {11, 12};
Line(20) = {11, 17};
Line(21) = {11, 19};

Curve Loop(1) = {15, 5, 6, -17};
Plane Surface(1) = {1};
Curve Loop(2) = {18, 20, -7, -17};
Plane Surface(2) = {2};
Curve Loop(3) = {21, -9, -8, -20};
Plane Surface(3) = {3};
Curve Loop(4) = {14, 4, -15, 16};
Plane Surface(4) = {4};
Curve Loop(5) = {13, -16, 18, 19, 12, 1, 2, 3};
Plane Surface(5) = {5};
Curve Loop(6) = {11, -19, 21, 10};
Plane Surface(6) = {6};

Physical Curve("inlet")    = {4, 5};
Physical Curve("outlet")   = {10, 9};
Physical Curve("symmtery") = {14, 13, 12, 11, 6, 7, 8};
Physical Curve("cylinder") = {3, 2, 1};
Physical Surface("domain") = {1, 4, 5, 2, 6, 3};
