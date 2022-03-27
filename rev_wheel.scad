include<dependencies.scad>

//spoke 1

function spk1()=let(
path=cr(pts1([[0,0],[0,5],[-13,20],[0,20],[-12,20],[0,78],[12,20],[0,64],[13,20],[0,5]]),5),
sec=cir(245,s=72),
prism=q_rot(["z90"],prism(sec,path)),
prism1=q_rot(["z90"],prism(f_offset(sec,5),path)),
//swp_prism_h(prism,prism1),

sketch=cr(pts1([[-100,0],[100,30,500],[100,-30]]),20),
path1=trns([0,0,202],cytz(cr(pts1([[0,0],[130,30,1000],[260,-10,0]]),20))),
arc1=arc(300,-15,15,s=30),
sec1=f_offset(cr([[0,0,10],c2t3([arc1[0]]).x+[0,0,10],each loop(arc1,2,28),c2t3([arc1[30]]).x+[0,0,10] ],10),-2.5),
sec2=trns([50,0,0],f_offset(cr([[0,0,30],c2t3([arc1[0]]).x+[0,0,10],each loop(arc1,2,28),c2t3([arc1[30]]).x+[0,0,10] ],20),-2.5)),
sec3=cr(pts1([[3,0],[0,3,1.5],[-3,0,1.5],[0,-3]]),5),
prism2=l_extrude(m_points_sc(sec1,10,2),300),
prism3=l_extrude(m_points_sc(sec2,10,2),300),
prism4=cyl(r=10,h=300,cp=[100,0,0],s=30),
surf=surf_extrude(sketch,path1),
surf1=trns([0,0,-5],surf),
fillet1=ipf(surf,prism2,3),
fillet2=ipf(surf1,flip(prism2),3,1),
path2=ip(surf,prism3),
path3=ip(surf,prism4),
prism5=p_extrudec(sec3,path2),
prism6=p_extrudec(sec3,path3),

sec4=cir(50),
path4=cr(pts1([[-5,200],[5,0,5],[0,40,5],[-25,15,5],[-5,0]]),10),
prism7=prism(sec4,path4),
arc2=arc(50,-15,15,s=10),
sec5=cr([[0,0,10],c2t3([arc2[0]]).x+[0,0,1],each loop(arc2,1,9),c2t3([arc2[10]]).x+[0,0,1] ],5),
prism8=l_extrude(sec5,300),
fillet3=ipf(surf,prism8,3,1),
fillet4=ipf(surf1,flip(prism8),3),
prism9=q_rot(["z90"],prism(cir(400,s=72),path)),

prism10=ipe(surf1,prism4,2,1),
prism11=[c2t3(c3t2(prism10[0])),each prism10],

prism12=ipe(surf1,prism3,2,1),
prism13=[c2t3(c3t2(prism12[0])),each prism12])
[sketch,path1,prism11,prism13,prism2,prism,prism3,prism4,fillet1,prism5,prism6,
// 0    1       2       3       4       5   6       7       8       9       10
fillet2,fillet3,fillet4];
// 11   12      13

a=spk1();

module spk1(){
render(){
difference(){
intersection(){
difference(){
surf_extrude(a[0],a[1],-5);
swp(a[2]);
swp(a[3]);
}
swp(a[4]);
swp(a[5]);}
swp(a[6]);
swp(a[7]);
}

intersct(a[5])
{
swp_c(a[8]);
swp_c(flip(a[9]));
swp_c(a[10]);
swp_c(a[11]);
}


swp_c(a[12]);
swp_c(a[13]);}

module intersct(prism){
for(i=[0:$children-1])
intersection(){
children(i);
swp(prism);
}}}

// spoke2
function spk2()=let(
path=cr(pts1([[0,0],[0,5],[-13,20],[0,20],[-12,20],[0,78],[12,20],[0,64],[13,20],[0,5]]),5),
sec=cir(245,s=72),
prism=q_rot(["z90"],prism(sec,path)),
prism1=q_rot(["z90"],prism(f_offset(sec,5),path)),
//swp_prism_h(prism,prism1),

sketch=cr(pts1([[-100,0],[100,30,500],[100,-30]]),20),
path1=trns([0,0,207],q_rot(["y7"],cytz(cr(pts1([[0,0],[130,30,1000],[260,-10,0]]),20)))),
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
sec5=cr([[0,0,10],c2t3([arc2[0]]).x+[0,0,1],each loop(arc2,1,9),c2t3([arc2[10]]).x+[0,0,1] ],5),
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

difference(){
intersection(){
difference(){
surf_extrude(b[0],b[1],-5);
//swp(prism11);
swp(b[2]);
}
swp(b[3]);
swp(b[4]);}
swp(b[5]);
//swp(prism4);
}

intersct(b[4])
{
swp_c(b[6]);
swp_c(b[7]);
//swp_c(prism6);
swp_c(flip(b[8]));
}


swp_c(b[9]);
swp_c(b[10]);
//swp(prism7);
module intersct(prism){
for(i=[0:$children-1])
intersection(){
children(i);
swp(prism);
}}}

// spoke3

function spk3()=let(
sec=cr(pts1([[0,0],[2.5,0,2.5],[0,45,2.5],[-5,0,2.5],[0,-45,2.5]]),5),
prism=trns([5,0,205],q_rot(["x90","z90"],l_extrude(m_points_sc(sec,30,.5),300))),

sec4=cir(50),
path4=cr(pts1([[-5,200],[5,0,5],[0,35,5],[-25,20,5],[-5,0]]),10),
prism4=q_rot(["z90"],prism(sec4,path4)),

p1=ip(prism4,prism),
p2=ip(prism4,rsz3dc(prism,bb(prism)+[5,5,5])),
p3=ip(rsz3dc(prism4,bb(prism4)+[5,5,5]),prism),
prism5=[for(i=[0:len(p1)-1])each i<len(p1)-1?[3p_3d_fillet(p3[i],p1[i],p2[i],3p_3d_r([p3[i],p1[i],p2[i]]))]:[3p_3d_fillet(p3[i],p1[i],p2[i],3p_3d_r([p3[i],p1[i],p2[i]]))
 ,3p_3d_fillet(p3[0],p1[0],p2[0],3p_3d_r([p3[0],p1[0],p2[0]]))
]],

path3=cr(pts1([[0,0],[0,5],[-13,20],[0,20],[-12,20],[0,78],[12,20],[0,64],[13,20],[0,5]]),5),
sec3=cir(245,s=72),
prism3=q_rot(["z90"],prism(sec3,path3)),

p4=ip(prism3,prism),
p5=ip(prism3,rsz3dc(prism,bb(prism)+[4,4,4])),
p6=ip(rsz3dc(prism3,bb(prism3)-[4,4,0]),prism),
prism6=[for(i=[0:len(p4)-1])each i<len(p4)-1?[3p_3d_fillet(p5[i],p4[i],p6[i],3p_3d_r([p5[i],p4[i],p6[i]]))]:[3p_3d_fillet(p5[i],p4[i],p6[i],3p_3d_r([p5[i],p4[i],p6[i]])),3p_3d_fillet(p5[0],p4[0],p6[0],3p_3d_r([p5[0],p4[0],p6[0]]))]])
[prism,prism3,prism4,prism5,prism6];
//  0   1       2       3       4
c=spk3();

module spk3(){
difference(){
intersection(){
swp(c[0]);
swp(c[1]);
}
swp(c[2]);
}
swp_c(c[3]);
swp_c(c[4]);
}

// Hub and Rim

function hub_rim()=let(
path=cr(pts1([[0,0],[0,5],[-13,20,5],[0,20,5],[-12,20,5],[0,78,5],[12,20,5],[0,64,1],[13,20,1],[0,5]]),5),
sec=cir(245,s=72),
prism=q_rot(["z90"],prism(sec,path)),
prism1=q_rot(["z90"],prism(f_offset(sec,5),path)),
sec4=cir(50),
path4=cr(pts1([[-5,200],[5,0,5],[0,35,5],[-25,20,5],[-5,0]]),10),
prism7=prism(sec4,path4),
prism8=cyl(r=20,h=300)

)[prism,prism1,prism7,prism8];

d=hub_rim();
module hub_rim(){
swp_prism_h(d[0],d[1]);
difference(){
swp(d[2]);
swp(d[3]);}}

// rendering the wheel

for(i=[0:360/6:360])
rotate([0,0,i])
rotate([0,0,15])
render(){
spk1();}

for(i=[0:360/6:360])
rotate([0,0,i])
rotate([0,0,-15])
render(){
spk2();}

for(i=[0:360/12:360])rotate([0,0,i])
spk3();

hub_rim();
