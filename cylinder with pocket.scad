include<dependencies.scad>

sec=cir(20);
path=[[0,0,0],[0,50,0]];

d=25;
prism=trns([0,25,0],l_extrude(m_points(cr(pts1([[-d/2,-d/2,4],[d,0,4],[0,d,4],[-d,0,4]]),10),1),50));


path1=ip(q_rot(["x-90"],cyl(r=20,h=50)),prism);
path2=ip(q_rot(["x-90"],cyl(r=18,h=50)),prism);

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

rotate([90,0,0]){
difference(){
swp_h(sec,path,-2);
swp(prism);}
intersection(){
swp_h(sec,path,-2);
swp_c(fillet);
}}
