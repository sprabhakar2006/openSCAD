include<dependencies.scad>

sec=circle(10);
path=cr(pts1([[-2,0],[2,0,2],[0,20,2],[-2,0]]),10);
path1=cr(pts1([[-2,0],[-1,0,.1],[-1,1,.1],[0,18,.1],[1,1,.1],[1,0]]),5);
prism=prism(sec,path);
prism1=prism(sec,path1);




swp_prism_h(prism,prism1);
