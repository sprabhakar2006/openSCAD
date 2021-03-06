include<dependencies.scad>

sec=circle(40,s=100);
path=cr(pts1([[-5,0],[5,0,5],[3,15,10],[-3,15,5],[-15,0]]),20);
prism=q_rot(["z50"],prism(sec,path));

swp(prism);

sec1=circle(7);
path1=cr(pts1([[0,0],[-2,35,2],[-2,0]]),5);
prism1=trns([34.5,0,20],prism(sec1,path1));

swp(prism1);

prism3=ipf(prism,prism1,3,1);

swp_c(flip(prism3));
