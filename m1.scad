include<dependencies.scad>

sec1=cir(10);
sec2=c3t2(q_rot(["z30"],trns([52.5,0],cir(5))));

sec3=2cir_tarc(20,10,[0,0],[52.5,0]*rm(30),45);
sec4=2cir_filleto(20,10,[0,0],[52.5,0]*rm(30),18.75)[0];
sec5=2p_arc(sec3[len(sec3)-1],sec4[0],r=20,cw=-1);
sec6=2p_arc(sec4[len(sec4)-1],sec3[0],10,-1);


sec7=concat(sec5,sec4,sec6,sec3);

linear_extrude(10)
difference(){
polygon(sec7);
polygon(sec1);
polygon(sec2);}

