include<dependencies.scad>

sec=cr(pts1([[-15,0,2.5],[0,15,3],[30,0,3],[0,-15,2.5],[5,0,2.5],[0,20,7],[-40,0,7],[0,-20,2.5]]),10);

//points(sec,.2);

path=cr(pts1([[-2,0],[2,0,2],[0,20,2],[-2,0]]),10);

prism=prism(sec,path,1);

//for(p=prism)points(p,.2);

swp(prism);