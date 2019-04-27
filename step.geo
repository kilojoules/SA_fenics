C = 2;
h1 = 0.01 / C;
h2 = 0.1 / C;
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, 1, 0, 1.0};
//+
Point(3) = {-0.2, 1, 0, 1.0};
//+
Point(4) = {-0.2, 2.5, 0, 1.0};
//+
Point(5) = {10, 2.5, 0, 1.0};
//+
Point(6) = {10, 0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 1};
//+
Line Loop(1) = {4, 5, 6, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Field[1] = Box;
//+
Field[1].VIn = h1;
//+
Field[1].VOut = h2;
//+
Field[1].XMax = 0.2;
//+
Field[1].XMin = -3.2;
//+
Field[1].YMax = 1.2;
//+
Field[1].YMin = -0.2;
//+
Background Field = 1;
