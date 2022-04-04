include<dependencies.scad>

function mobile_stand()=let(
sec=cr(pts1([[-75,0,2],[150,0,2],[0,10,2],[-75,7,500],[-75,-7,2]]),10),
path=cr(pts1([[-1,0],[1,0,1],[0,5,1],[-1,0]]),5),
prism=trns([0,0,3],q_rot(["x120"],prism(sec,path))),
sec1=cr(pts1([[-75,0,2],[150,0,2],[0,35,2],[-75,20,500],[-75,-20,2]]),10),
prism1=trns([0,-19,3],q_rot(["x120"],prism(sec1,path))),
//swp(trns([0,-6,0],q_rot(["x"],[120],prism(sec1,path))));
sec2=cr(pts1([[-75,0,2],[150,0,2],[0,85,2],[-150,0,2]]),5),
prism2=trns([0,-85,0],prism(sec2,path)),
sec3=cr(pts1([[0,7],[7*cos(240),7*sin(240),4],[19,0,10],[20*cos(60),20*sin(60)]]),5),
sec4=cr(pts1([[0,8],[8*cos(240),8*sin(240),4],[14,0]]),5),
sec01=f_offset(sec1,-10),
prism01=trns([0,-18.5,3],q_rot(["x120"],l_extrude(sec01,7))),
sec5=cr(pts1([[-70,-30,2],[137,0,2],[0,-47,2],[-105,0,2]]),5))
[prism,prism1,prism2,sec3,sec4,prism01,sec5,];
//0     1       2    3      4   5       6

module mobile_stand(a){
//difference(){
union(){
swp(a[0]);
swp(a[1]);
swp(a[2]);
translate([74,-6.5,1.5])
rotate([90,0,0])
rotate([0,-90,0])
linear_extrude(148)
p_lineo(a[3],4);
translate([74,-6.5-19,1.5])
rotate([90,0,0])
rotate([0,-90,0])
linear_extrude(148)
p_lineo(a[4],4);

/*translate([0,-78,4.5])
linear_extrude(1)
rotate([0,0,180])
text("SDMC  Satbari",5,font="tahoma: style=Bold",);*/

}

//swp(a[5]);
//translate([0,0,-1])
//linear_extrude(10)polygon(a[6]);}
}

data=mobile_stand();




mobile_stand(data);
//translate([0,-90,0])
//mobile_stand(data);
