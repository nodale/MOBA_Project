Li    = 6;
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
Point(12) = { Lo/3,   H*Gamma, 0, h_3};
Point(13) = { Lo,     H*Gamma, 0, h_out};
Point(14) = { 2*Lo/3, H*Gamma, 0, h_2_3};

Line(1)  = {1, 2};
Line(2)  = {2, 3};
Line(3)  = {3, 4};
Line(4)  = {4, 7};
Line(5)  = {7, 9};
Line(6)  = {9, 11};
Line(7)  = {11, 5};
Line(9)  = {6, 10};
Line(10) = {10, 8};
Line(11) = {8, 1};
Line(12) = {2, 7};
Line(14) = {2, 12};
Line(15) = {9, 12};
Line(16) = {12, 8};
Line(17) = {12, 14};
Line(18) = {14, 13};
Line(19) = {11, 14};
Line(20) = {14, 10};
Line(21) = {5, 13};
Line(22) = {13, 6};

Curve Loop(1) = {3, 4, -12, 2};
Plane Surface(1) = {1};
Curve Loop(2) = {5, 15, -14, 12};
Plane Surface(2) = {2};
Curve Loop(3) = {1, 14, 16, 11};
Plane Surface(3) = {3};
Curve Loop(4) = {10, -16, 17, 20};
Plane Surface(4) = {4};
Curve Loop(5) = {6, 19, -17, -15};
Plane Surface(5) = {5};
Curve Loop(6) = {7, 21, -18, -19};
Plane Surface(6) = {6};
Curve Loop(7) = {22, 9, -20, 18};
Plane Surface(7) = {7};

Physical Curve("inlet")  = {3};
Physical Curve("walls")  = {4, 2, 1, 11, 5, 10, 6, 9, 7};
Physical Curve("outlet") = {21, 22};

Physical Surface("domain") = {1, 2, 3, 4, 5, 6, 7};

Transfinite Surface {1} = {3, 4, 7, 2};
Transfinite Surface {3} = {1, 2, 12, 8};
Transfinite Surface {2} = {2, 7, 9, 12};
Transfinite Surface {4} = {8, 12, 14, 10};
Transfinite Surface {5} = {12, 9, 11, 14};
Transfinite Surface {6} = {14, 11, 5, 13};
Transfinite Surface {7} = {10, 14, 13, 6};

Transfinite Curve {5, 11, 14} = 100 Using Progression 1;
Transfinite Curve {4, 2} = 80 Using Progression 1;
Transfinite Curve {6, 10} = 80 Using Progression 1;
Transfinite Curve {7, 9, 18} = 60 Using Progression 1;
Transfinite Curve {17} = 80 Using Progression 1;
Transfinite Curve {3, 12, 1, 15, 16, 19, 20, 21, 22} = 30 Using Progression 1;

//+
Transfinite Curve {3, 12, 1, 15, 16, 20, 19, 21, 22} = 10 Using Progression 1;
//+
Transfinite Curve {7, 18, 9} = 25 Using Progression 1;
//+
Transfinite Curve {6, 17, 10} = 40 Using Progression 1;
//+
Transfinite Curve {5, 14, 11} = 60 Using Progression 1;
//+
Transfinite Curve {4, 2} = 20 Using Progression 1;
