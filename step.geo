h = 1. / 30.;
Point(newp) = {0.0,0.0,0.0,h};
Point(newp) = {3.0,0.0,0.0,h};
Point(newp) = {3.0,2.5,0.0,h};
Point(newp) = {-3.0,2.5,0.0,h};
Point(newp) = {-3.0,1.0,0.0,h};
Point(newp) = {0.0,1.0,0.0,h};
//+
Line(1) = {6, 5};
//+
Line(2) = {5, 4};
//+
Line(3) = {4, 3};
//+
Line(4) = {3, 2};
//+
Line(5) = {2, 1};
//+
Line(6) = {1, 6};
//+

Line Loop(1) = {1, 2, 3, 4, 5, 6};
//+
Plane Surface(1) = {1};
//+
///Transfinite Surface {1};
//+
///Transfinite Curve {1, 2, 3, 4} = 10 Using Progression 1;
//+
///Recombine Surface {1};

