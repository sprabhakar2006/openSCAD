include<dependencies.scad>


a=rands(0,10,100);
 b=rands(0,7,100);
 pnts=[for(i=[0:len(a)-1])[a[i],b[i]]];
 points(pnts,.3);
 cn_hull=concave_hull(pnts,1);
 color("green")
 p_line(cn_hull,.1);
// color("magenta")
// p_line(c_hull(pnts),.05);
 
echo(s_int2(seg(cn_hull)));
 

 
 
