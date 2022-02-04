include<dependencies.scad>

cyl1=trns([-50,0,0],q_rot(["y90"],cyl(r=10.1,h=100)));
cyl2=q_rot(["y90","z90"],cyl(r=10,h=50));
cyl3=q_rot(["y90","z270"],cyl(r=10,h=50));
fillet1=ipf(cyl1,cyl2,3,1);
fillet2=ipf(cyl1,cyl3,3,1);

swp_c(flip(fillet1));
swp_c(flip(fillet2));
%swp(cyl1);
%swp(cyl2);
%swp(cyl3);
