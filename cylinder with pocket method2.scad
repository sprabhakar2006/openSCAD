include<dependencies.scad>

prism=l_extrude(circle(20,s=100),50);
prism1=l_extrude(circle(18,s=100),50);



sec=m_points_sc(cr(pts1([[-15,-15,5],[30,0,5],[0,30,5],[-30,0,5]]),10),10);
prism2=trns([0,0,25],q_rot(["y90"],l_extrude(sec,50)));

//%swp(prism2);

prism3=[trns([10,0,0],ipe(prism,flip(prism2),1)[0]),each ipe(prism,flip(prism2),1)];
prism4=[trns([-10,0,0],ipe(prism1,prism2,1,1)[0]),each ipe(prism1,prism2,1,1)];
prism5=[each prism3,each flip(prism4)];

//swp(prism3);

difference(){
swp_prism_h(prism,prism1);
swp(prism5);
}
