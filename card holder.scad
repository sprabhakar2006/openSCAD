include<dependencies.scad>

sec=cr(pts1([[0,0,1],[95,0,1],[0,10,1],[-95,0,1]]),10);
prism=l_extrude(sec,50);
prism1=l_extrude(offset(sec,1),50);


prism2=l_extrude(offset(sec,1),1);


sec1=cr(pts1([[10,0,5],[50,0,5],[10,30],[-70,0]]),10);

prism3=trns([12.5,2,21],q_rot(["x90"],l_extrude(sec1,5)));

difference(){
union(){
swp_prism_h(prism1,prism);
swp(prism2);}
swp(prism3);}