include<dependencies.scad>

 a=rands(0,10,100);
 b=rands(0,7,100);
 pnts=[for(i=[0:len(a)-1])[a[i],b[i]]];
 points(pnts,.3);
 c_hull=c_hull(pnts);
 color("green")
 p_line(c_hull,.1);
