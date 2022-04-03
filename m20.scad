include<dependencies.scad>

//echo(175*sin(45));

sec=cr(pts1([[0,0],[124,0,3],[27*cos(45),27*sin(45),1],[15*cos(135),15*sin(135),1],[5*cos(180+45),5*sin(180+45),1],[10*cos(-45),10*sin(-45),1],[17*cos(180+45),17*sin(180+45),1],[10*cos(135),10*sin(135),1],[8*cos(225),8*sin(225),1],[20*cos(135),20*sin(135),1],[8*cos(45),8*sin(45),1],[105*cos(135),105*sin(135),1],[8*cos(225),8*sin(225),1],[20*cos(135),20*sin(135),1],[8*cos(45),8*sin(45),1],[10*cos(135),10*sin(135),1],[17*cos(45),17*sin(45),1],[10*cos(-45),10*sin(-45),1],[5*cos(45),5*sin(45),1],[15*cos(135),15*sin(135),1],[17*cos(225),17*sin(225),1],[0,30,2],[-7,0,1]]),5);

sec1=cr(pts1([[0,0,2],[17,0,2],[0,85,2],[-17,17,2]]),5);
path1=[[0,0],[0,6]];
prism=trns([123.5,7.7,30],q_rot(["x90","z45"],prism(sec1,path1)));

prism1=trns([5,129,130],q_rot(["x90","z45"],prism(sec1,path1)));

sec2=cr(pts1([[0,0,4],[10,0,4],[0,40,4],[-10,0,4]]),5);
path2=[[0,0],[0,8]];
prism2=trns([-0.25,140,40],q_rot(["x90","z90"],prism(sec2,path2)));
prism3=trns([-0.25,140,170],q_rot(["x90","z90"],prism(sec2,path2)));

sec4=cr(pts1([[0,0,0],[27,0,5],[0,90,5],[-27,0,0]]),5);
path4=[[0,0],[0,3]];
prism4=trns([90,25,0],q_rot(["z45"],prism(sec4,path4)));

sec5=cr(pts1([[0,0],[120,0,0],[32*cos(135),32*sin(135),1],[8*cos(45),8*sin(45),1],[105*cos(135),105*sin(135),1],[8*cos(225),8*sin(225),1],[33*cos(135),33*sin(135),0]]),5);

sec6=cr(pts1([[0,0,5],[90,0,5],[0,110,5],[-90,90,5]]),5);
path6=[[0,0],[0,7]];
prism6=trns([15,6,25],q_rot(["x90"],prism(sec6,path6)));
prism7=trns([-1,15,25],q_rot(["x90","z90"],prism(sec6,path6)));
prism8=trns([95,30,25],q_rot(["x90","z135"],prism(sec6,path6)));

rotate([0,-90,0])rotate([0,0,-90])
difference(){
union(){
difference(){
linear_extrude(250)polygon(sec);
swp(prism);
swp(prism1);
swp(prism2);
swp(prism3);
}
swp(prism4);}
translate([0,0,-.5])linear_extrude(251)offset(r=-4)polygon(sec5);
swp(prism6);
swp(prism7);
swp(prism8);
}
