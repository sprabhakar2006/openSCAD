include<dependencies.scad>

c1=cir(10,[-35,0]);
c2=cir(10,[35,0]);
c3=cir(10,[0,44.06]);

c4=cir(7.5,[-35,0]);
c5=cir(7.5,[35,0]);
c6=cir(7.5,[0,44.06]);
//p_line(c1,.2);
//p_line(c2,.2);
//p_line(c3,.2);

p0=[-45,0];
p1=cir_p_t(c1,[0,0]);
arc1=2p_arc(p0,p1,10,-1,10);
p2=[0,0,27.5];
p3=p_cir_t([0,0],c2);
p4=cir_p_t(c2,[45,35]);
arc2=2p_arc(p3,p4,10,-1,10);
p5=[45,35,10];
p6=p_cir_t([45,35],c3);
p7=cir_p_t(c3,[-45,35]);
arc3=2p_arc(p6,p7,10,-1,10);
p8=[-45,35,10];

sec=cr([each arc1,p2,each arc2,p5,each arc3,p8],10);

p10=[-45+2.5,0];
p11=cir_p_t(c4,[0,2.5]);
arc11=2p_long_arc(p10,p11,7.5,1,30);
p12=[0,2.5,32.5];
p13=p_cir_t([0,2.5],c5);
p14=cir_p_t(c5,[45-2.5,25]);
arc12=2p_long_arc(p13,p14,7.5,1,30);
p15=[45-2.5,28.5];
p16=2cir_tangent(7.5,7.5,[35,28.5],[0,44.06]);
arc13=2p_arc(p15,p16[2],7.5,-1,10);
p17=2cir_tangent(7.5,7.5,[-35,28.5],[0,44.06]);
arc14=2p_long_arc(p16[3],p17[0],7.5,1,30);
p18=[-45+2.5,28.5];
arc15=2p_arc(p17[1],p18,7.5,-1);

sec1=cr([c2t3([arc11[0]]).x+[0,0,2],each loop(arc11,5,25),c2t3([arc11[29]]).x+[0,0,2],p12,c2t3([arc12[0]]).x+[0,0,2],each loop(arc12,5,25),c2t3([arc12[29]]).x+[0,0,2], each arc13,c2t3([arc14[0]]).x+[0,0,2], each loop(arc14,4,26),c2t3([arc14[30]]).x+[0,0,2], each arc15],10);
//p_line(sec,.2);
//p_line(sec1,.2);

path1=cr(pts1([[-1,0],[1,0,1],[0,10,1],[-1,0]]),5);
path2=cr(pts1([[-1,3],[1,0,1],[0,7,1],[1,0]]),5);
prism1=prism(sec,path1);
prism2=prism(sec1,path2);

c7=cir(5,[-35,0]);
c8=cir(5,[35,0]);
c9=cir(5,[0,44.06]);
path3=cr(pts1([[.5,0],[-.5,.5],[0,9],[.5,.5]]),5);
prism3=prism(c7,path3);
prism4=prism(c8,path3);
prism5=prism(c9,path3);
render(){
difference(){
swp(prism1);
swp(prism2);
swp(prism3);
swp(prism4);
swp(prism5);
}}

