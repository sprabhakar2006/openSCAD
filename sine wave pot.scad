include<dependencies.scad>

//sine wave function
sec=cir(15);
path=[for(i=[0:30])[5*sin(12*i),i]];
prism=prism(sec,path);
//function line x=sin(y)
linear_extrude(.5)p_lineo(path);

surf_extrudec(sec,path,2,2);

sec1=cr(pts1([[-15,-15,6],[30,0,6],[0,30,6],[-30,0,6]]),10);
translate([50,0,0])
surf_extrudec(sec1,path,1,2);

translate([0,50,0])
swp(prism);

translate([50,50,0])
swp(prism(sec1,path));