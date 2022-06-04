include<dependencies.scad>
r=3;

sec1=circle(245,s=72);
path1=cr(pts1([[0,0],[0,5,5],[-13,20,5],[0,20,5],[-12,20,5],[0,78,5],[12,20,5],[0,64,5],[13,20,5],[0,5]]),5);
prism1=q_rot(["z50"],prism(sec1,path1));
prism7=q_rot(["z50"],prism(offset(sec1,-r),path1));

sec2=circle(45);
path2=cr(pts1([[0,200],[5,0,5],[0,40,5],[-25,15,5],[-24.5,0]]),5);
prism2=q_rot(["z30"],prism(sec2,path2));

sec3=cr(pts1([[0,0],[2.5,0,2.5],[0,45,2.5],[-5,0,2.5],[0,-45,2.5]]),5);
prism3=trns([20,0,205],q_rot(["x90","z90"],l_extrude(remove_extra_points(m_points(sec3,1)),300)));
prism4=trns([20,0,205],q_rot(["x90","z90"],l_extrude(remove_extra_points(m_points(offset(sec3,r),1)),300)));

swp(prism2);


p1=ip(prism2,prism4);
p2=ip(prism2,prism3);
p3=ip(surf_offset(prism2,r),prism3);


prism5=[for(p=cpo([p3,p2,p1]))3p_3d_fillet(p.x,p.y,p.z,3p_3d_r(p)*1.5)];
prism6=[each prism5,prism5[0]];


for(i=[0:360/12:360])rotate([0,0,i]){
swp_c(prism6);

intersection(){
swp(prism3);
swp(prism1);}
}







