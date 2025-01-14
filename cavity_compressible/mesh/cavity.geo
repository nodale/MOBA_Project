Li = 3;		// Length of upstream domain.
Lo = 10;	// Length of downstream domain.
H  = 1; 	// Height of channel.

L = 6; 		// Length of cavity.
D = 1; 		// Reynolds number based on cavity depth. Must be 1.

h_cavity = 0.07;
h_in     = 0.1;
h_out    = 0.1;
h_shear  = 0.02;

Point(1)  = {0,     0, 0, h_shear};
Point(2)  = {-Li,   0, 0, h_in};
Point(3)  = {-Li,   H, 0, h_in};
Point(4)  = { L+Lo, H, 0, h_out};
Point(5)  = { L+Lo, 0, 0, h_out};
Point(6)  = { L,    0, 0, h_shear};
Point(7)  = { L,   -D, 0, h_cavity};
Point(8)  = { 0,   -D, 0, h_cavity};
Point(9)  = { 0,    H, 0, h_cavity};
Point(10) = { L,    H, 0, h_cavity};

Line(1) = {8, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 9};
Line(5) = {9, 1};
Line(6) = {1, 6};
Line(7) = {6, 10};
Line(8) = {10, 9};
Line(9) = {8, 7};
Line(10) = {7, 6};
Line(11) = {10, 4};
Line(12) = {4, 5};
Line(13) = {5, 6};

Curve Loop(1) = {2, 3, 4, 5};
Plane Surface(1) = {1};
Curve Loop(2) = {6, 7, 8, 5};
Plane Surface(2) = {2};
Curve Loop(3) = {13, 7, 11, 12};
Plane Surface(3) = {3};
Curve Loop(4) = {9, 10, -6, -1};
Plane Surface(4) = {4};

Physical Curve("inlet")  = {3};
Physical Curve("walls")  = {2, 8, 4, 11, 13, 10, 9, 1};
Physical Curve("outlet") = {12};

Physical Surface("domain") = {1, 2, 3, 4};
