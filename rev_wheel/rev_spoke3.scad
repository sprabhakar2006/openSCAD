include<dependencies.scad>

function spk3()=let(
sec=cr(pts1([[0,0],[2.5,0,2.5],[0,45,2.5],[-5,0,2.5],[0,-45,2.5]]),5),
prism=trns([30,0,205],q_rot(["x90","z90"],l_extrude(m_points_sc(sec,30,.5),300))),

sec4=cir(50),
path4=cr(pts1([[-5,200],[5,0,5],[0,35,5],[-25,20,5],[-5,0]]),10),
prism4=q_rot(["z90"],prism(sec4,path4)),

p1=ip(prism4,prism),
p2=ip(prism4,rsz3dc(prism,bb(prism)+[4,4,4])),
p3=ip(rsz3dc(prism4,bb(prism4)+[4,4,4]),prism),
prism5=[for(i=[0:len(p1)-1])each i<len(p1)-1?[3p_3d_fillet(p3[i],p1[i],p2[i],3p_3d_r([p3[i],p1[i],p2[i]]))]:[3p_3d_fillet(p3[i],p1[i],p2[i],3p_3d_r([p3[i],p1[i],p2[i]]))
 ,3p_3d_fillet(p3[0],p1[0],p2[0],3p_3d_r([p3[0],p1[0],p2[0]]))
]],
//prism5=ipe(prism4,prism,2,1),
path3=cr(pts1([[0,0],[0,5],[-13,20,5],[0,20,5],[-12,20,5],[0,78,5],[12,20,5],[0,64,5],[13,20,1],[0,5]]),5),
sec3=cir(245,s=72),
prism3=q_rot(["z90"],prism(sec3,path3)),

p4=ip(prism3,prism),
p5=ip(prism3,rsz3dc(prism,bb(prism)+[4,4,4])),
p6=ip(rsz3dc(prism3,bb(prism3)-[4,4,0]),prism),
prism6=[for(i=[0:len(p4)-1])each i<len(p4)-1?[3p_3d_fillet(p5[i],p4[i],p6[i],3p_3d_r([p5[i],p4[i],p6[i]]))]:[3p_3d_fillet(p5[i],p4[i],p6[i],3p_3d_r([p5[i],p4[i],p6[i]])),3p_3d_fillet(p5[0],p4[0],p6[0],3p_3d_r([p5[0],p4[0],p6[0]]))]],
prism7=rsz3dc(prism3,bb(prism3)+[1,1,0])
)
[prism,prism3,prism4,prism5,prism6,prism7];
//  0   1       2       3       4   5
c=spk3();

module spk3(){
render()scale(1.002){
intersct(){
swp(c[5]);
swp(c[0]);
swp_c(c[3]);
swp_c(c[4]);
}}

module intersct(){
for(i=[1:$children-1])
intersection(){
children(i);
children(0);
}}
//swp(c[2]);


}

for(i=[0:360/12:360])rotate([0,0,i])
spk3();
