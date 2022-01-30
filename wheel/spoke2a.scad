include<dependencies.scad>

//spoke 2

function spoke2a()=let(
sec=cr(pts1([[-100,0],[100,30,500],[100,-30],[0,.01],[-100,30,500],[-100,-30]]),20),
path=trns([0,0,220],q_rot(["y10","z15"],cytz(cr(pts1([[0,0],[130,30,1000],[260,-10,0]]),20)))),

a=arc(300,0,30,s=20),
sec1=cr([[0,0,5],[a[0].x,a[0].y,10],[a[len(a)-1].x,a[len(a)-1].y,10]],10),

sec2=trns([50*cos(15),50*sin(15),0],cr([[0,0,10],[a[0].x,a[0].y,5],[a[len(a)-1].x,a[len(a)-1].y,5]],20)),

path3=ip(p_extrude(sec,path),l_extrude(m_points(sec2,10),300)),
sec3=cr(pts1([[0,0],[0,3,1.5],[-3,0,1.5],[0,-3]]),5),

path4=[for(i=[0:len(path3)-2])if(norm(path3[i]-path3[i+1])>.1)path3[i]],

sec5=cr(pts1([[-100,-5,1],[100,30,500],[100,-30,1],[0,5,1],[-100,30,500],[-100,-30,1]]),20),

sec6=cr([[0,0],each arc(50,0,30,s=10)],5),
 
path7=cr(pts1([[0,0],[0,5],[-13,20],[0,20],[-12,20],[0,78],[12,20],[0,64],[13,20],[0,5]]),5),
sec7=cir(245,s=72),
b=arc(50,0,30,s=10),
sec10=cr([[0,0,5],[b[0].x,b[0].y,2],each loop(b,1,len(b)-2),[b[len(b)-1].x,b[len(b)-1].y,5]],5))
[sec,sec1,sec2,sec3,sec5,sec6,sec7,path,path3,path4,path7,sec10];
//0    1     2   3   4    5      6   7   8       9   10   11

module spoke2a(b){
difference(){
intersection(){
p_extrude(b[4],b[7]);
swp(l_extrude(b[1],300));
swp(prism(b[6],b[10]));
}
swp(l_extrude(b[2],300));}

intersection(){
p_extrudec(b[3],b[9]);
swp(prism(b[6],b[10]));}

intersection(){
swp_c(ipf(p_extrude(b[0],b[7]),l_extrude(m_points(f_offset(b[1],-2.5),10),300),3));
swp(prism(b[6],b[10]));}




swp_c(ipf(p_extrude(b[0],b[7]),l_extrude(b[11],300),3,1));
}

b=spoke2a();
for(i=[0:60:300])rotate([0,0,i])
spoke2a(b);


