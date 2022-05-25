include<dependencies.scad>

sec=cr(pts1([[0,0,.5],[7,5,2],[5,7,3],[-5,7,5],[-7,5,5]]),20);
s=[for(i=[-5:.5:20])offset(sec,i)];
for(p=s)p_line(p,.1);

color("green")
p_line(sec,.1);