include<dependencies.scad>

module anchor(){
sec=cr(pts1([[0,0,1],[10,0,1],[0,5,1],[-10,0,1]]),5);
sec1=cr(pts1([[0,0],[10,0],[0,9,1],[-1,1,1],[-8,0,1],[-1,-1,1]]),5);

prism=l_extrude(sec,10);
prism1=trns([0,7,0],q_rot(["x90"],l_extrude(sec1,10)));

prism2=trns([5,-2,7],q_rot(["x-90"],cylinder(d=5,h=10)));

sec2=cr(pts1([[0,0],[2,0],[0,10],[-2,0]]),5);
prism3=trns([4,-2,7],l_extrude(sec2,5));
prism4=ipe(plane([0,0,1],100),prism,1,1);
swp(prism4);
difference(){
intersection(){
swp(prism);
swp(prism1);}
swp(prism3);
swp(prism2);}}

for(i=[0,12,24,36,48])translate([i+5,-2.5,1.9])anchor();

sec3=cr(pts1([[0,-5,2],[65,0,2],[0,10,2],[-65,0,2]]),5);
prism4=l_extrude(sec3,2);

swp(prism4);