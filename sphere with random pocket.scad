include<dependencies.scad>


seca=cir(20,s=5);
secb=c3t2(q_rot(["z36"],cir(5,s=5)));

r=2;
secc=cr([for(i=[0:len(seca)-1])each[[seca[i].x,seca[i].y,r],[secb[i].x,secb[i].y,r]]],5);

prism=l_extrude(m_points(secc,1),50);
difference(){
swp(spr(20));
swp(spr(18));
swp(prism);
}

path1=ip(spr(20),prism);
path2=ip(spr(18),prism);

//points(path1,.5);
//points(path2,.5);

fillet=let(
fillet=[for(i=[0:len(path1)-1])
let(
i_plus=i<len(path1)-1?i+1:0,
p0=path1[i],
p2=path2[i],
cp=p0+(p2-p0)/2,
axis=path1[i_plus]-path1[i],
p1=cp+q(axis,p0-cp,-90),
fillet=3p_3d_arc([p0,p1,p2],s=10)
)fillet]
)[for(i=[0:len(fillet)-1])each
i<len(fillet)-1?[fillet[i]]:[fillet[i],fillet[0]]

];


swp_c(flip(fillet));