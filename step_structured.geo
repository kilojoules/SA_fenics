C = 2;
h1 = 1. / 1. / C;
h2 = 1. / 30. / C;
Point(newp) = {0.0,0.0,0.0,h2};
Point(newp) = {0.0,2.5,0.0,h2};
Point(newp) = {10.0,2.5,0.0,h1};
Point(newp) = {10.0,0.0,0.0,h1};

Point(newp) = {0.0,2.5,0.0,h2};
Point(newp) = {-3.0,2.5,0.0,h1};
Point(newp) = {-3.0,1.0,0.0,h2};
Point(newp) = {0.0,1.0,0.0,h2};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+

//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 5};
//+

Line Loop(6) = {1, 2, 3, 4};
Line Loop(8) = {5, 6, 7, 8};
//+
Plane Surface(7) = {6};
//+
Transfinite Surface {6};
//+
Plane Surface(9) = {8};
//+
Transfinite Surface {8};
//+
////Transfinite Line {1, 2, 3, 4, 5, 6} = 10 Using Progression 1;
//+
Recombine Surface {7, 9};

