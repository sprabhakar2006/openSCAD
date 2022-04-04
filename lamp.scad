include<dependencies.scad>


path=[[42,0],[62,50],[25,100],[25,180]];

path1=cytz(bez(path,.01));
    
path2=[for(i=[0:len(path1)-1])
    q([0,0,1],path1[i],i/(len(path1)-1)*(180-2))];
    
prism=q_rot(["z2"],[for(i=[0:len(path2)-1])let(theta=ang(path2[i].x,path2[i].y))trns(path2[i],q_rot([str("z",theta)],f_offset(sqr([6,6],true),i/len(path2)*-1.5)))]);

translate([0,0,2])
for(i=[0:360/15:360])rotate([0,0,i])
{

swp(prism);

mirror([0,1,0])
swp(prism);
}

linear_extrude(2)
difference(){
circle(r=45);
circle(r=30);
}

