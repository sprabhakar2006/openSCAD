include<dependencies.scad>
module spoke1a(){
sec=cr(pts1([[-100,0],[100,30,500],[100,-30],[0,.01],[-100,30,500],[-100,-30]]),20);
sec3=cr(pts1([[0,0],[0,3,1.5],[-3,0,1.5],[0,-3]]),5);

sec5=cr(pts1([[-100,-5],[100,30,500],[100,-30],[0,5],[-100,30,500],[-100,-30]]),20);

path7=cr(pts1([[0,0],[0,5],[-13,20],[0,20],[-12,20],[0,78],[12,20],[0,64],[13,20],[0,5]]),5);
sec7=cir(245,s=72);

path10=trns([0,0,185],q_rot(["y-5","z-15"],cytz(cr(pts1([[0,0],[130,30,1000],[260,-10,0]]),20))));
a=arc(300,-30,0,s=20);
sec11=cr([[0,0,5],[a[0].x,a[0].y,5],[a[len(a)-1].x,a[len(a)-1].y,5]],10);
sec12=trns([50*cos(-15),50*sin(-15),0],cr([[0,0,25],[a[0].x,a[0].y,5],[a[len(a)-1].x,a[len(a)-1].y,5]],30));

path11=ip(p_extrude(sec,path10),l_extrude(m_points(sec12,10),300));
path12=[for(i=[0:len(path11)-2])if(norm(path11[i]-path11[i+1])>.1)path11[i]];
 
path13=ip(p_extrude(sec,path10),cyl(r=10,h=300,cp=[85*cos(-15),85*sin(-15)]));
path14=[for(i=[0:len(path13)-2])if(norm(path13[i]-path13[i+1])>.1)path13[i]];
 
sec15=cr([[0,0,5],each arc(50,-30,0,s=10)]);

render(){
intersection(){
difference(){
intersection(){
p_extrude(sec5,path10);
swp(l_extrude(sec11,300));
}
swp(trns([85*cos(-15),85*sin(-15)],cyl(r=10,h=300)));
swp(l_extrude(sec12,300));}

swp(prism(sec7,path7));}}


render(){
intersection(){
union(){
p_extrudec(sec3,path12);
swp_c(ipf(p_extrude(sec,path10),l_extrude(m_points(f_offset(sec11,-2.5),10),300),3));}
swp(prism(sec7,path7));}}
p_extrudec(sec3,path14);
swp_c(ipf(p_extrude(sec,path10),l_extrude(sec15,300),3,1));
}