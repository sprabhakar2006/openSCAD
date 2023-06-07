include<dependencies.scad>

// explanation:
// radius of the torus = a*d
// section radius of the torus = d
// p in number of cycles of the wrapping coil over torus
// q in the number of turns of the wrapping coil over torus

a=5;
d=5;
p=3;
q=35;

step=.25;
path=[for(t=[0:step:360-step])[d*cos(p*t)*(a+cos(q*t)),
d*sin(p*t)*(a+cos(q*t)),
-d*sin(q*t)]];

r=1;
sec=circle(r);

sol=path_extrudec(sec,path);

swp_c(flip(sol));


sec1=circle(d-r);
path1=c2t3(circle(a*d));

sol1=path_extrudec(sec1,path1);
%swp_c(flip(sol1));
