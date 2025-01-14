// Cylinder center.
Point(1) = {0, 0, 0, 1.0};

r = 0.5;

// coarse
/*h_cylinder = 0.2;
h_wake     = 0.6;
h_inner    = 0.8;
h_outer    = 1.0;*/

// middle
h_cylinder = 0.2;
h_wake     = 0.4;
h_inner    = 0.5;
h_outer    = 2.0;

// fine
/*h_cylinder = 0.1;
h_wake     = 0.15;
h_inner    = 0.2;
h_outer    = 1.5;*/

// Cylinder "egdes".
Point(2) = {r*Cos(45*Pi/180),  r*Sin(45*Pi/180),  0, h_cylinder};
Point(3) = {r*Cos(135*Pi/180), r*Sin(135*Pi/180), 0, h_cylinder};
Point(4) = {r*Cos(225*Pi/180), r*Sin(225*Pi/180), 0, h_cylinder};
Point(5) = {r*Cos(315*Pi/180), r*Sin(315*Pi/180), 0, h_cylinder};

// Cylinder lines.
Circle(1) = {3, 1, 2};
Circle(2) = {2, 1, 5};
Circle(3) = {5, 1, 4};
Circle(4) = {4, 1, 3};

// Domain corners.
Point(6) = {-16, 15, 0, h_outer};
Point(7) = {30, 15, 0,  h_outer};
Point(8) = {-16, -15, 0, h_outer};
Point(9) = {30, -15, 0, h_outer};

// Mesh structure points square around cylinder.
Point(10) = {-5, 5, 0, h_inner};
Point(13) = {-5, -5, 0, h_inner};

// Wake
Point(11) = {5, 5, 0, h_wake};
Point(12) = {5, -5, 0, h_wake};
Point(18) = {30, -5, 0, h_wake}; // Mesh structure points outlet.
Point(19) = {30, 5, 0, h_wake};

// Mesh structure points inlet.
Point(14) = {-16, 5, 0, h_outer};
Point(15) = {-16, -5, 0, h_outer};

// Mesh structure points top.
Point(16) = {-5, 15, 0, h_outer};
Point(17) = { 5, 15, 0, h_outer};

// Mesh structure points bottom.
Point(20) = {-5, -15, 0, h_outer};
Point(21) = {5, -15, 0, h_outer};

// Domain edges.
Line(5) = {8, 15};	// Inlet.
Line(6) = {15, 14};
Line(7) = {14, 6};
Line(8) = {6, 16};	// Top.
Line(9) = {16, 17};
Line(10) = {17, 7};
Line(11) = {7, 19};	// Outlet.
Line(12) = {19, 18};
Line(13) = {18, 9};
Line(14) = {9, 21};	// Bottom.
Line(15) = {21, 20};
Line(16) = {20, 8};

// Inner lines.
Line(17) = {13, 10};
Line(18) = {10, 11};
Line(19) = {11, 12};
Line(20) = {12, 13};
Line(21) = {13, 15};
Line(22) = {10, 14};
Line(23) = {10, 16};
Line(24) = {11, 17};
Line(25) = {11, 19};
Line(26) = {12, 18};
Line(27) = {12, 21};
Line(28) = {13, 20};
Line(29) = {4, 13};
Line(30) = {3, 10};
Line(31) = {2, 11};
Line(32) = {5, 12};

// Definition of surfaces.
Curve Loop(1) = {7, 8, -23, 22};
Plane Surface(1) = {1};
Curve Loop(2) = {23, 9, -24, -18};
Plane Surface(2) = {2};
Curve Loop(3) = {24, 10, 11, -25};
Plane Surface(3) = {3};
Curve Loop(4) = {19, 26, -12, -25};
Plane Surface(4) = {4};
Curve Loop(5) = {27, -14, -13, -26};
Plane Surface(5) = {5};
Curve Loop(6) = {28, -15, -27, 20};
Plane Surface(6) = {6};
Curve Loop(7) = {5, -21, 28, 16};
Plane Surface(7) = {7};
Curve Loop(8) = {6, -22, -17, 21};
Plane Surface(8) = {8};
Curve Loop(9) = {17, -30, -4, 29};
Plane Surface(9) = {9};
Curve Loop(10) = {30, 18, -31, -1};
Plane Surface(10) = {10};
Curve Loop(11) = {2, 32, -19, -31};
Plane Surface(11) = {11};
Curve Loop(12) = {29, -20, -32, 3};
Plane Surface(12) = {12};

// Boundaries
Physical Line("inlet")     = {5, 6, 7};
Physical Line("outlet")    = {11, 12, 13};
Physical Line("cylinder")  = {3, 4, 1, 2};
Physical Surface("domain") = {8, 1, 2, 10, 9, 12, 11, 4, 3, 5, 6, 7};
Physical Curve("symmetry") = {8, 9, 10, 14, 15, 16};
