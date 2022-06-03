include<dependencies.scad>

module bend(){

sketch=m_points_so([[-20,0],[20,0]],10);
path=cytz(m_points_so(cr(pts1([[-10,20],[55,-.1,1],[60*cos(45),60*sin(45)]]),5),10,m=.5));
surf=surf_extrude(sketch,path);

point1=[45,-15];
point2=[45,15];
cir=circle(r=7.5,[45+22.5,0]);
p1=p_cir_t(point1,cir);
p2=cir_p_t(cir,point2);
a1=ang_v(p1-[45+22.5,0]);
a2=ang_v(p2-[45+22.5,0]);

sec1=trns([0,0,-5],l_extrude([point1,each arc(7.5,a1,a2+360,[45+22.5,0]),point2],5));
//p_line(cir,.2);
sec2=trns([45,0,0],q_rot(["y-45"],trns([-45,0,0],sec1)));

//swp(sec2);

arc1=3d_3p_fillet([0,-15,0],[45,-15,0],sec2[1][1],1);
arc2=3d_3p_fillet(sec2[1][len(sec2[1])-2],[45,15,0],[0,15,0],1);

prism1=[[0,-15,0],each arc1, each loop(sec2[1],1,len(sec2[1])-2),each arc2,[0,15,0]];

arc3=3d_3p_fillet([0,-15,-5],[45,-15,-5],sec2[0][0],6);
arc4=3d_3p_fillet(sec2[0][len(sec2[1])-1],[45,15,-5],[0,15,-5],6);

prism2=[[0,-15,-5],each arc3,each loop(sec2[0],1,len(sec2[0])-2),each arc4,[0,15,-5]];

prism3=trns([30,0,14],[prism2,prism1]);
difference(){
swp(prism3);
swp(trns([75,0,14],q_rot(["y-45"],trns([22.5,0,-11],cylinder(r=3.75,h=12)))));}}

sec=cr(pts1([[0,-15],[45,0,15],[0,30,15],[-45,0]]),10);

prism=l_extrude(sec,22.5);

sec1=sqr([30,9]);
prism1=trns([-.1,-4.5,-0.1],l_extrude(sec1,22.7));

sec4=circle(6);
path4=pts([[0,-31/2],[0,5],[-3,0],[0,21],[3,0],[0,5]]);
prism7=trns([12,0,11.25],q_rot(["x90"],prism(sec4,path4)));

sec5=cr(pts1([[0,0],[9,0],[0,20,4.5],[-9,0,4.5]]),10);
prism8=trns([55.5,-15.5,0],l_extrude(sec5,20));

difference(){
union(){
swp(prism);
bend();}

swp(prism1);
swp(trns([30,0,-.1],cylinder(r=7.5,h=22.7)));
swp(prism7);
swp(prism8);
}

