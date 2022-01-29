include<dependencies.scad>

//spoke3

module spoke3a(){
sec7=cir(245,s=72);

path7=cr(pts1([[0,0],[0,5],[-13,20],[0,20],[-12,20],[0,78],[12,20],[0,64],[13,20],[0,5]]),5);
sec8=cir(45);
path8=cr(pts1([[0,200],[5,0,5],[0,40,5],[-25,15,5],[-24.5,0]]),5);
prism8=q_rot(["z30"],prism(sec8,path8));
sec15=cr(pts1([[0,0],[2.5,0,2.5],[0,45,2.5],[-5,0,2.5],[0,-45,2.5]]),5);
path15=[[0,0,205],[300,0,205]];

prism9=q_rot(["z50"],prism(sec7,path7));

intersection(){
p_extrude(sec15,path15);
swp(prism9);
}
prism1=p_extrude(m_points(sec15,1),path15);
prism2=ipf(prism8,prism1,1,1);
swp_c(flip(prism2));


prism3=ipf(prism9,flip(prism1),1);
swp_c(prism3);
}