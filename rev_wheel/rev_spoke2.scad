include<dependencies.scad>


function spk2()=let(
path=cr(pts1([[0,0],[0,5],[-13,20],[0,20],[-12,20],[0,78],[12,20],[0,64],[13,20],[0,5]]),5),
sec=cir(245,s=72),
prism=q_rot(["z90"],prism(sec,path)),
prism1=q_rot(["z90"],prism(f_offset(sec,5),path)),
//swp_prism_h(prism,prism1),

sketch=cr(pts1([[-100,0],[100,30,500],[100,-30]]),20),
path1=trns([0,0,212],q_rot(["y7"],cytz(cr(pts1([[0,0],[130,30,1000],[260,-10,0]]),20)))),
arc1=arc(300,-15,15,s=30),
sec1=f_offset(cr([[0,0,10],c2t3([arc1[0]]).x+[0,0,10],each loop(arc1,2,28),c2t3([arc1[30]]).x+[0,0,10] ],10),-2.5),
sec2=trns([50,0,0],f_offset(cr([[0,0,10],c2t3([arc1[0]]).x+[0,0,10],each loop(arc1,2,28),c2t3([arc1[30]]).x+[0,0,10] ],20),-2.5)),
sec3=cr(pts1([[3,0],[0,3,1.5],[-3,0,1.5],[0,-3]]),5),
prism2=l_extrude(m_points_sc(sec1,10,2),300),
prism3=l_extrude(m_points_sc(sec2,10,2),300),
//prism4=cyl(r=10,h=300,cp=[100,0,0],s=30),
surf=surf_extrude(sketch,path1),
surf1=trns([0,0,-5],surf),
fillet1=ipf(surf,prism2,3),
fillet2=ipf(surf1,flip(prism2),3,1),
path2=ip(surf,prism3),
//path3=ip(surf,prism4),
prism5=p_extrudec(sec3,path2),
//prism6=p_extrudec(sec3,path3),

sec4=cir(50),
path4=cr(pts1([[-5,200],[5,0,5],[0,40,5],[-25,15,5],[-5,0]]),10),
prism7=prism(sec4,path4),
arc2=arc(50,-15,15,s=10),
sec5=cr([[0,0,15],c2t3([arc2[0]]).x+[0,0,1],each loop(arc2,1,9),c2t3([arc2[10]]).x+[0,0,1] ],5),
prism8=l_extrude(sec5,300),
fillet3=ipf(surf,prism8,3,1),
fillet4=ipf(surf1,flip(prism8),3),
prism9=q_rot(["z90"],prism(cir(400,s=72),path)),

//prism10=ipe(surf1,prism4,2,1),
//prism11=[c2t3(c3t2(prism10[0])),each prism10],

prism12=ipe(surf1,prism3,2,1),
prism13=[c2t3(c3t2(prism12[0])),each prism12])
[sketch,path1,prism13,prism2,prism,prism3,fillet1,prism5,fillet2,fillet3,fillet4];
// 0    1       2       3       4       5   6       7       8       9       10

b=spk2();

module spk2(){
render()scale(1.002){
difference(){
intersection(){
intersection(){
difference(){
surf_extrude(b[0],b[1],-5);
//swp(prism11);
swp(b[2]);
}
//swp(b[3]);
swp(b[4]);
}
swp(b[3]);}
swp(b[5]);
//swp(prism4);
}

intersct()
{
swp(b[4]);
swp_c(flip(b[6]));
swp_c(flip(b[7]));
//swp_c(prism6);
swp_c(b[8]);
}

swp_c(b[9]);
swp_c(flip(b[10]));}
//swp(prism7);
module intersct(){
for(i=[1:$children-1])
intersection(){
children(i);
children(0);
}}}

for(i=[0:360/6:300])
rotate([0,0,i])
rotate([0,0,-15])
//render(){
spk2();
//}
