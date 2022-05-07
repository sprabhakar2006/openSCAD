include<dependencies.scad>


 sec=cr(pts1([[0,0,1],[8,3,3],[5,7,1],[-8,0,2],[-5,20,1]]),20);
 prism=l_extrude(sec,30);
 plane1=plane([0,0,1],50);
 %swp(prism);
 %swp(plane1);
 prism1=ipe(plane1,prism,r=2,option=1,s=10);

 swp(prism1);

 plane2=trns([0,0,20],plane([0,0,1],50));
 prism2=ipe(plane2,flip(prism),r=1,option=1,s=10);
 //%swp(plane2);
 swp(prism2);