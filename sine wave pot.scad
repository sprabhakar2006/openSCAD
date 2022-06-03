include<dependencies.scad>

//sine wave function
sec=circle(15);
path=[for(i=[0:30])[5*sin(12*i),i]];
prism=prism(sec,path);
prism1=prism(f_offset(sec,-1),path);

//function line x=sin(y)
linear_extrude(.5)p_lineo(path);

swp_prism_h(prism,prism1);


sec1=cr(pts1([[-15,-15,8],[30,0,8],[0,30,8],[-30,0,8]]),10);
prism2=prism(sec1,path);
prism3=prism(f_offset(sec1,-1),path);
translate([50,0,0])
swp_prism_h(prism2,prism3);

translate([0,50,0])
swp(prism);

translate([50,50,0])
swp(prism(sec1,path));
