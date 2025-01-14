Li    = 5;
Lo    = 60;
H     = 2;    // Reynolds number is based on half of height, so H/2.
Gamma = 0.5;

h_in   = 0.1;
h_step = 0.05;
h_3    = 0.1;
h_2_3  = 0.2;
h_out  = 0.3;

Point(1)  = { 0,      0,       0, h_step};
Point(2)  = { 0,      H*Gamma, 0, h_step};
Point(3)  = {-Li,     H*Gamma, 0, h_in};
Point(4)  = {-Li,     H,       0, h_in};
Point(5)  = { Lo,     H,       0, h_out};
Point(6)  = { Lo,     0,       0, h_out}; 
Point(7)  = {  0,     H,       0, h_step};
Point(8)  = { Lo/3,   0,       0, h_3};
Point(9)  = { Lo/3,   H,       0, h_3};
Point(10) = { 2*Lo/3, 0,       0, h_2_3};
Point(11) = { 2*Lo/3, H,       0, h_2_3}; 

Line(1)  = {1, 2};
Line(2)  = {2, 3};
Line(3)  = {3, 4};
Line(4)  = {4, 7};
Line(5)  = {7, 9};
Line(6)  = {9, 11};
Line(7)  = {11, 5};
Line(8)  = {5, 6};
Line(9)  = {6, 10};
Line(10) = {10, 8};
Line(11) = {8, 1};
Line(12) = {2, 7};
Line(13) = {9, 8};
Line(14) = {11, 10};

Curve Loop(1) = {2, 3, 4, -12};
Plane Surface(1) = {1};
Curve Loop(2) = {11, 1, 12, 5, 13};
Plane Surface(2) = {2};
Curve Loop(3) = {10, -13, 6, 14};
Plane Surface(3) = {3};
Curve Loop(4) = {9, -14, 7, 8};
Plane Surface(4) = {4};

Physical Curve("inlet")  = {3};
Physical Curve("walls")  = {4, 2, 1, 11, 5, 10, 6, 9, 7};
Physical Curve("outlet") = {8};

Physical Surface("domain") = {1, 2, 3, 4};
