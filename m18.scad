include<dependencies.scad>

sec=cr(pts1([[-55/2,0],[55/2,5,200],[55/2,-5],[-1.5,3,0],[-52/2,5,200],[-52/2,-5]]),10);
sec1=f_offset(sec,1);
sec2=f_offset(sec1,-.5);
sec3=f_offset(sec,-.5);
sec4=f_offset(sec1,-1.1);
sec5=cr(pts1([[-5,0],[8,0],[0,3],[-3,0,1],[0,3,.5],[-2,0,.5],[-3,-3,.5]]),5);
prism=q_rot(["x90"],l_extrude(sec1,10));
prism1=q_rot(["x90"],l_extrude(sec,45));
prism2=trns([0,-1,-.6],q_rot(["x90"],l_extrude(sec2,8)));
prism3=trns([0,-1,-.6],q_rot(["x90"],l_extrude(sec3,43)));
prism4=trns([0,-8,-1.6],q_rot(["x90"],l_extrude(sec4,5)));
rotate([-90,0,0])
translate([-35,0,0])
union(){
difference(){
    union(){
swp(prism);
swp(prism1);}
//swp(prism2);
//swp(prism3);
//swp(prism4);

//swp(trns([-9+55/2-10,-44,0],l_extrude(sqr([9,10]),15)));
//swp(trns([-5+55/2-10-2,-35,0],l_extrude(sqr([5,4]),15)));

//swp(trns([-9+55/2-10-17-9,-44,0],l_extrude(sqr([9,10]),15)));
//swp(trns([-5+55/2-10-2-21-5,-35,0],l_extrude(sqr([5,4]),15)));

swp(trns([-7/2,-44,0],l_extrude(sqr([7,20]),15)));
}
translate([-7/2,-44,2.5]){
translate([-.25,0,0])
linear_extrude(3)p_line(sqr([7.5,20.5]),.5);
swp(trns([1,25.5,0],q_rot(["x90"],cub([5,3,20]))));
swp(trns([1+5,5+4,3],q_rot(["x90","z90","y180"],l_extrude(sec5,5))));
    }

sec6=trns([-3.5,0],sqr([7,2]));
points(sec6,.2);

path6=cytz(cr(pts1([[0,0],[.1,2,2],[7,3,10],[7,0.25]]),5));

swp(trns([13,-34,4.4],q_rot(["x180","z-90"],p_extrude(sec6,path6))));

swp(trns([-13,-34,4.4],q_rot(["x180","z-90"],p_extrude(sec6,path6))));}

//sec7=[[0,0],[3,0],[0,-3,]];
//translate([-2.5,5.5,18.5])
//rotate([0,90,0])
//linear_extrude(5)polygon(sec7);

rotate([-90,0,0])
translate([35,0,0])
union(){
difference(){
    union(){
swp(prism);
swp(prism1);}
//swp(prism2);
//swp(prism3);
//swp(prism4);

//swp(trns([-9+55/2-10,-44,0],l_extrude(sqr([9,10]),15)));
//swp(trns([-5+55/2-10-2,-35,0],l_extrude(sqr([5,4]),15)));

//swp(trns([-9+55/2-10-17-9,-44,0],l_extrude(sqr([9,10]),15)));
//swp(trns([-5+55/2-10-2-21-5,-35,0],l_extrude(sqr([5,4]),15)));

swp(trns([-7/2,-44,0],l_extrude(sqr([7,20]),15)));}
translate([-7/2,-44,2.5]){
translate([-.25,0,0])
linear_extrude(3)p_line(sqr([7.5,20.5]),.5);
swp(trns([1,25.5,0],q_rot(["x90"],cub([5,3,20]))));
swp(trns([1+5,5+4,3],q_rot(["x90","z90","y180"],l_extrude(sec5,5))));
    }

sec6=trns([-3.5,0],sqr([7,2]));
points(sec6,.2);

path6=cytz(cr(pts1([[0,0],[.1,2,2],[7,3,10],[7,.25]]),5));

swp(trns([13,-34,4.4],q_rot(["x180","z-90"],p_extrude(sec6,path6))));

swp(trns([-13,-34,4.4],q_rot(["x180","z-90"],p_extrude(sec6,path6))));}
