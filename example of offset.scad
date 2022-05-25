include<dependencies.scad>

sec=cr(pts1([[0,0,.2],[8,3,3],[5,7,1],[-8,0,2],[-5,20,1]]),20);
s=[for(i=[-4:.5:20])offset(sec,i)];
for(p=s)p_line(p,.1);

color("green")
p_line(sec,.1);
