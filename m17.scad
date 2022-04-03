include<dependencies.scad>



arc2=2cir_tarc(20,10,[0,0],[52.5,0],45);
arc0=2cir_filleto(20,10,[0,0],[52.5,0],18.75)[0];
arc1=
let(
v1=arc0[len(arc0)-1]-[52.5,0],
v2=arc2[0]-[52.5,0],
a1=ang(v1.x,v1.y),a2=ang(v2.x,v2.y)+360,
)arc(10,a1,a2,[52.5,0]);

arc3=let(
v1=arc2[len(arc2)-1]-[0,0],
v2=arc0[0]-[0,0],
a1=ang(v1.x,v1.y),
a2=ang(v2.x,v2.y)
)arc(20,a1,a2,[0,0]);

sec=concat(arc0,arc1,arc2,arc3);

rotate([0,0,30])
linear_extrude(10)
difference(){
polygon(sec);
polygon(cir(10));
polygon(cir(5,[52.5,0]));}

echo(cw(sec));