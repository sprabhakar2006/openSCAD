include<dependencies.scad>

sec1=2cir_tangent(17.5,7.5,[0,0],[14.5,0]);
sec2=2cir_tangent(17.5,7.5,[0,0],[-14.5,0]);

prism1=l_extrude(sec1,5.5);
prism2=l_extrude(sec2,5.5);
prism3=cyl(r=17.5,h=5.5);
prism4=cyl(r=7.5,h=5.5,cp=[14.5,0]);
prism5=cyl(r=7.5,h=5.5,cp=[-14.5,0]);

sec=cir(10);
path=c2t3(cr(pts1([[0,0],[0,29,17.5],[-44,0]]),10));

prism6=trns([-40,0,29],q_rot(["y-90"],cyl(r=39/2,h=4)));

rotate([90,0,0])
swp_h(sec,path,-3.5);

prism7=trns([0.001,0,.001],q_rot(["x90"],p_extrude(sec,path)));
fillet1=ipf(prism3,prism7,1.25);
fillet2=ipf(prism6,flip(prism7),1.25,1);


difference(){
union(){
swp(prism1);
swp(prism4);
swp(prism2);
swp(prism5);
swp(prism3);
}
swp(cyl(d=13,h=5.5));
swp(trns([0,0,-.5],cyl(d=6,h=7,cp=[14.5,0])));
swp(trns([0,0,-.5],cyl(d=6,h=7,cp=[-14.5,0])));
}


difference(){
swp(prism6);
translate([-40,0,29])
rotate([0,-90,0]){
for(i=[0:90:360])rotate([0,0,i])
swp(cyl(d=5,h=4,cp=[15,0]));
swp(cyl(d=13,h=4));}
}
swp_c(fillet1);
swp_c(fillet2);