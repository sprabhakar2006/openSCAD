include<dependencies.scad>

sec=cir(60,s=72);
path=cr(pts1([[-5,0],[5,0,5],[3,20,10],[-3,20,5],[-15,0]]),5);
prism=q_rot(["z50"],prism(sec,path));

fillet_radius=2;
sec1=m_points_sc(cr(pts1([[-15,0,2.5],[0,15,3],[30,0,3],[0,-15,2.5],[5,0,2.5],[0,20,7],[-40,0,7],[0,-20,2.5]]),10),10);

prism1=trns([0,40,10],q_rot(["x90","z230"],[c2t3(sec1),trns([0,0,130],scl2d_c(sec1,.6))]));

//sec2=m_points(f_offset(cr(pts1([[-15,0,2.5],[0,15,3],[30,0,3],[0,-15,2.5],[5,0,2.5],[0,20,7],[-40,0,7],[0,-20,2.5]]),10),fillet_radius),1);
//
//prism2=trns([0,40,10],q_rot(["x90","z230"],[c2t3(sec2),trns([0,0,130],scl2d_c(sec2,.6))]));
//
//p1=ip(prism,prism1);
//p2=sort_points(p1,ip(q_rot(["z50"],prism(sec,f_offset(path,fillet_radius))),prism1));
//p3=sort_points(p1,ip(prism,prism2));
//prism3=fillet(p1,p2,p3,fillet_radius);

prism3=ipf(prism,prism1,fillet_radius,1);
swp(prism);
swp(prism1);
swp_c(prism3);
