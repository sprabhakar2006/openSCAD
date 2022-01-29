include<dependencies.scad>


module hub_rim(){
sec7=cir(245,s=144);

path7=cr(pts1([[0,0],[0,5],[-13,20],[0,20],[-12,20],[0,78],[12,20],[0,64],[13,20],[0,5]]),5);
sec8=cir(45);
path8=cr(pts1([[0,200],[5,0,5],[0,40,5],[-25,15,5],[-24.5,0]]),5);
prism8=q_rot(["z30"],prism(sec8,path8));

sec9=[each path7,each flip([for(p=path7)p+[5,0]])];
prism9=q_rot(["z30"],prism(sec7,path7));

swp(prism8);

rotate_extrude($fn=72) 
translate([245,0,0])polygon(sec9);}