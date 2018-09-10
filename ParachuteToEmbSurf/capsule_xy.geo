// This is the mesh size
cl_capsule = 0.1;

//The top circle of the capsule is centered at (0, y=ya), with radius -xa,
//Point (xa, ya) is on the top circle, modify it to change the size of the capsule. 
//xa = -0.31142370413265841
xa = 0.31142370413265841;
ya = 0.0;


xb = -xa;
x4 = 2.25;
x3 = x4 - 1.012*Tan(36.9*Pi/180);
x2 = x3 - 0.506*Tan(59*Pi/180);
x1 = x2 - 0.5*Tan(33.95*Pi/180);

y1 = 0;
y2 = -0.5;
y3 = -0.5-0.506;
y4 = -0.5-0.506-1.012;
o_x = 0.0;
o_y = y4 - 2.25*Tan(20*Pi/180) + 0.5* 2.25/Cos(20*Pi/180);
nose_x_a = o_x + 1.125*Cos(-70*Pi/180);
nose_y_a = o_y + 1.125*Sin(-70*Pi/180);
nose_x_b = o_x + 1.125*Cos(-110*Pi/180);
nose_y_b = o_y + 1.125*Sin(-110*Pi/180);

Printf("o_y before scale is %g", o_y);

Printf("nose_y_b before scale is %g", nose_y_b);


scale = (xb - xa)/(2*x1);
Printf("scale is %g",scale);
x4 *= scale;
x3 *= scale;
x2 *= scale;
x1 *= scale;
y1 *= scale;
y2 *= scale;
y3 *= scale;
y4 *= scale;
o_x *= scale;
o_y *= scale;

nose_x_a *= scale;
nose_x_b *= scale;
nose_y_a *= scale;
nose_y_b *= scale;



Point(6) = {xa,ya,0,cl_capsule};
Point(7) = {0 + xa + x1,0 +ya - y1,0,cl_capsule};
Point(8) = {x1 + xa + x1,y1 +ya - y1,0,cl_capsule};
Point(9) = {x2 + xa + x1,y2 +ya - y1,0,cl_capsule};
Point(10) = {x3 + xa + x1,y3 +ya - y1,0,cl_capsule};
Point(11) = {x4 + xa + x1,y4 +ya - y1,0,cl_capsule};
Point(12) = {o_x + xa + x1,o_y +ya - y1,0,cl_capsule};
Point(13) = {nose_x_a + xa + x1,nose_y_a +ya - y1,0,cl_capsule};
Point(14) = {nose_x_b + xa + x1,nose_y_b +ya - y1,0,cl_capsule};
Point(15) = {-x4 + xa + x1,y4 +ya - y1,0,cl_capsule};
Point(16) = {-x3 + xa + x1,y3 +ya - y1,0,cl_capsule};
Point(17) = {-x2 + xa + x1,y2 +ya - y1,0,cl_capsule};

Printf("y_min is %g",nose_y_a +ya - y1);

Printf("scale is %g", scale);

Printf("y_a - y_1 after scale is %g", ya - y1);

Printf("o_y after scale is %g", o_y);

Printf("nose_y_b after scale is %g", nose_y_b);

//+
Line(5) = {6, 17};
//+
Line(6) = {17, 16};
//+
Line(7) = {16, 15};
//+
Line(8) = {15, 14};
//+
Circle(9) = {14, 12, 13};

//+
Line(10) = {13, 11};
//+
Line(11) = {11, 10};
//+
Line(12) = {10, 9};
//+
Line(13) = {9, 8};
//+
Line(14) = {8, 7};
//+
Line(15) = {7, 6};


//+
Line Loop(16) = {5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};



Extrude {{0, 1, 0}, {0, 0, 0}, Pi/2} {
  Line{7, 6, 5, 15, 14, 13, 12, 11, 10, 9, 8};
}

Extrude {{0, 1, 0}, {0, 0, 0}, Pi/2} {
  Line{32, 35, 39, 43, 47, 51, 55, 17, 21, 25, 29};
}
//+
Physical Surface("CapusleSurface") = {20, 58, 85, 50, 77, 81, 54, 89, 73, 46, 42, 69, 93, 38, 97, 65, 24, 28, 100, 34, 61, 31};
