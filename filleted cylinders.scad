include<dependencies.scad>

cyl1=trns([-50,0,0],q_rot(["y90"],cylinder(r=10.1,h=100,s=150)));
cyl2=q_rot(["y90","z90"],cylinder(r=10,h=50,s=150));
cyl3=q_rot(["y90","z270"],cylinder(r=10,h=50,s=150));
fillet1=ipf(cyl1,cyl2,4,1);
fillet2=ipf(cyl1,cyl3,4,1);

swp_c(flip(fillet1));
swp_c(flip(fillet2));
%swp(cyl1);
%swp(cyl2);
%swp(cyl3);
