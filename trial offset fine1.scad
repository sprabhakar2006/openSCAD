include<dependencies.scad>

sec=cr(pts1([[0,0,.2],[8,3,3],[5,7,1],[-8,0,2],[-5,20,1]]),20);
//sec=cr(pts1([[0,0,.5],[7,5,2],[5,7,3],[-5,7,5],[-7,5,5]]),20);

path=cr(pts1([[2,0],[-2,0,2],[0,7,4],[-4,0]]),20);
//path=cr(pts1([[2,0],[-2,0,2],[0,7,5],[-5,0]]),20);

prism=prism1(sec,path);

//for(p=prism)p_line3dc(p,.1);

swp(prism);