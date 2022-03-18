include<dependencies.scad>

sec=cr(pts1([[0,-15],[45,0,15],[0,30,15],[-45,0]]),10);

prism=l_extrude(sec,22.5);

sec1=sqr([30,9]);
prism1=trns([-.1,-4.5,-0.1],l_extrude(sec1,22.7));

l1=[[0,0],[75,0]];l2=pts([[75,0],[30*cos(45),30*sin(45)]]);
l3=offst_l(l1,5);l4=offst_l(l2,5);
v=q_rot(["x90"],l2);
ip=i_p2d(l3,l4);

sec2=[each l1,l2[1],l4[1],ip,l3[0]];
prism4=trns([0,15,0],q_rot(["x90"],l_extrude(sec2,30)));

sec3=[[0,-15],p_cir_t([0,-15],cir(7.5,[22.5,0])),cir_p_t(cir(7.5,[22.5,0]),[0,15]),[0,15]];
prism2=trns([0,0,-10.5],l_extrude(sec3,11));
prism3=trns([0,0,-10.5],l_extrude(cir(7.5,[22.5,0]),11));
prism6=trns([0,0,-10.5],l_extrude(cir(3.75,[22.5,0]),11));

prism5=trns([0,-35/2,-10],cub([35,35,10]));

sec4=cir(6);
path4=pts([[0,-31/2],[0,5],[-3,0],[0,21],[3,0],[0,5]]);
prism7=trns([12,0,11.25],q_rot(["x90"],prism(sec4,path4)));

sec5=cr(pts1([[0,0],[9,0],[0,20,4.5],[-9,0,4.5]]),10);
prism8=trns([55.5,-15.5,0],l_extrude(sec5,20));

difference(){
union(){
swp(prism);
translate([0,0,9+5])
difference(){
swp(prism4);
translate([75,0,0])
rotate([0,-45,0]){
difference(){
swp(prism5);
difference(){
swp(prism2);
swp(trns([22.5,0,-11],cyl(r=3.75,h=12)));
}
swp_prism_h(prism3,prism6);}

}}}

swp(prism1);
swp(trns([30,0,-.1],cyl(r=7.5,h=22.7)));
swp(prism7);
swp(prism8);
}



