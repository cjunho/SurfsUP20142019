//24*(4*8)
//for Mpi parellel n4
// rate number of elemenet between hight and longitudinal: 32:12.  longitudinal should be more than hight
//8 phisical surfaces. success 

ls=.5;

Lcm=12;
Lcp=12;
trho=4;
max=-7.002421210670288;
h1=8;
h2=4;

Point(1) = {-Lcm,max-h1,0.0, ls}; 
Point(2) = {-Lcm,max-h1-5*h2,0.0, ls}; 
Point(3) ={-Lcm,max-h1-6*h2,0.0, ls}; 
Point(8)= {Lcp,max-h1-6*h2,0.0, ls}; 
Point(4)= {Lcp,max-h1-5*h2,0.0, ls}; 
Point(5) = {Lcp,max-h1,0.0, ls}; 
Point(6) = {Lcp,max,0.0, ls}; 
Point(7) = {-Lcm,max,0.0, ls}; 

Point(9) = {-Lcm,max-0.5*h1,0.0, ls}; 

Point(10) = {Lcp,max-0.5*h1,0.0, ls}; 

Point(11) = {-Lcm,max-h1-h2,0.0, ls}; 
Point(12) = {-Lcm,max-h1-2*h2,0.0, ls}; 
Point(13) = {-Lcm,max-h1-3*h2,0.0, ls}; 
Point(14) = {-Lcm,max-h1-4*h2,0.0, ls}; 

Point(15) = {Lcp,max-h1-h2,0.0, ls}; 
Point(16) = {Lcp,max-h1-2*h2,0.0, ls}; 
Point(17) = {Lcp,max-h1-3*h2,0.0, ls}; 
Point(18) = {Lcp,max-h1-4*h2,0.0, ls}; 





//+
Line(1) = {7, 9};
//+
Line(2) = {9, 1};
//+
Line(3) = {1, 11};
//+
Line(4) = {11, 12};
//+
Line(5) = {12, 13};
//+
Line(6) = {13, 14};
//+
Line(7) = {14, 2};
//+
Line(8) = {2, 3};
//+
Line(9) = {3, 8};
//+
Line(10) = {8, 4};
//+
Line(11) = {4, 18};
//+
Line(12) = {18, 17};
//+
Line(13) = {17, 16};
//+
Line(14) = {16, 15};
//+
Line(15) = {15, 5};
//+
Line(16) = {5, 10};
//+
Line(17) = {10, 6};
//+
Line(18) = {6, 7};
//+
Line(19) = {9, 10};
//+
Line(20) = {1, 5};
//+
Line(21) = {11, 15};
//+
Line(22) = {12, 16};
//+
Line(23) = {13, 17};
//+
Line(24) = {14, 18};
//+
Line(25) = {2, 4};
//+
Curve Loop(1) = {1, 19, 17, 18};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {2, 20, 16, -19};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {3, 21, 15, -20};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {4, 22, 14, -21};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {5, 23, 13, -22};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {6, 24, 12, -23};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {7, 25, 11, -24};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {8, 9, 10, -25};
//+
Plane Surface(8) = {8};
//+
Physical Surface(1) = {1};
//+
Physical Surface(2) = {2};
//+
Physical Surface(3) = {3};
//+
Physical Surface(4) = {4};
//+
Physical Surface(5) = {5};
//+
Physical Surface(6) = {6};
//+
Physical Surface(7) = {7};
//+
Physical Surface(8) = {8};

Transfinite Line{18,19,20,21,22,23,24,25,9}=12*(Lcm+Lcp)+1;
Transfinite Line{1,17}=32*h2+1;
Transfinite Line{2,16}=32*h2+1;
Transfinite Line{3,15}=32*h2+1;
Transfinite Line{4,14}=32*h2+1;
Transfinite Line{5,13}=32*h2+1;
Transfinite Line{6,12}=32*h2+1;
Transfinite Line{7,11}=32*h2+1;
Transfinite Line{8,10}=32*h2+1;

Transfinite Surface{1,2,3,4,5,6,7,8};
Recombine Surface{1,2,3,4,5,6,7,8};

