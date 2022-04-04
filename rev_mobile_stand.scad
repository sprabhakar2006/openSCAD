include<dependencies.scad>

sec=cr(pts1([[0,0,1],[100,0,2],[15*cos(120),15*sin(120),1],[-5,0,1],[10*cos(-60),10*sin(-60),2],[-15,0,5],[50*cos(120),50*sin(120),1],[-5,0,1],[50*cos(-60),50*sin(-60),3],[-69.5,0,1]]),10);

path=cr(pts1([[-.5,0],[.5,0,.5],[0,150,.5],[-.5,0]]),5);

prism=q_rot(["x90"],prism(sec,path));



sketch=cr(pts1([[-50,0],[50,20,100],[50,-10,100],[50,20,100],[50,-30,100],[50,20]]),10);
path1=[[0,0,21],[100,0,21]];

surf=surf_extrude(sketch,path1);

difference(){
swp(prism);
surf_base(surf,70);}

v=q_rot(["x90"],c2t3([[1,0]*rm(210)])).x;

plane1=trns([72+2.5+3,-190,0],plane(v,500));
path2=sort_points(surf[0],ip(plane1,surf));
//%swp(plane1);
//%surf_extrude(sketch,path1);
intersection(){
p_line3d(path2,3,$fn=20);
swp(prism);
}
//points(surf[0],2);