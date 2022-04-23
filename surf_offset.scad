include<dependencies.scad>

sec=cir(10);
path=cr(pts1([[2,0],[-2,0,2],[-1,10,2],[-4,0]]),5);
prism=prism(sec,path);
//swp(prism);
prism1= surf_offset(prism,-1);
swp_prism_h(prism,prism1);