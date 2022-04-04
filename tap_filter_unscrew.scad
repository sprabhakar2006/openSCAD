include<dependencies.scad>

sec=cir(26/2);


line=trns([0,12,0],[[-10,0],[10,0]]);

linear_extrude(8){
p_line(sec,2);
intersection(){
p_lineo(line,1);
circle(r=13);}
mirror([0,1,0])
intersection(){
p_lineo(line,1);
circle(r=13);}

line1=[[0,-25],[0,25]];
difference(){
p_lineo(line1,3);
circle(r=13);}}
