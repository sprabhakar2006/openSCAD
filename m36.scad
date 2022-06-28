include<dependencies.scad>


cir=circle(20,s=150);
arc1=2p_arc([0,0],[4,28],58,-1,30);
arc2=2p_arc([-2.5,2.5],[4-2.5,28],58-2.5,-1,30);
path=cr([[-5,0],[0,0,5],each loop(arc1,3,27),[4,28,1],[3,28]],20);
path1=cr([[-5,2.5],[-2.5,2.5,2.5],each loop(arc2,3,27),[4-2.5,28,1],[4-1.5,28]],20);
prism=q_rot(["z30"],prism1(cir,path));
prism1=q_rot(["z30"],prism1(cir,path1));
swp_prism_h(prism,prism1);

swp(l_extrude(circle(15),2.5));
point=[20,4];
cp=[24+6,21-3];
arc3=arc(7.75,-60,155,cp,s=70);
path2=c2t3([point,each arc3]);
path3=c2t3([point,arc3[0]]);
sec2=cr(pts1([[-1.25,-2,1.25],[2.5,0,1.25],[0,4,1.25],[-2.5,0,1.25]]),20);
rotate([90,0,0])
p_extrude(sec2,path2);

prism2=q_rot(["y46.4"],l_extrude(sec2,50));
prism3=q_rot(["x90"],p_extrude(sec2,path3));

fillet1=ipf(prism,prism2,3,1,s=10);
fillet2=ipf(prism,prism3,3,s=10);
swp_c(fillet1);
swp_c(fillet2);
