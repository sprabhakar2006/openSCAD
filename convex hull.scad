include<dependencies.scad>
a=rands(2,20,100);
b=rands(0,15,100);

pnts=[for(i=[0:len(a)-1])[a[i],b[i]]];

points(pnts,.2);

p_line(c_hull(pnts),.1);
