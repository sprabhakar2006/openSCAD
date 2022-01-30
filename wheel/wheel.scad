include<dependencies.scad>

function spoke1a()=let(

sec=cr(pts1([[-100,0],[100,30,500],[100,-30],[0,.01],[-100,30,500],[-100,-30]]),20),
sec3=cr(pts1([[0,0],[0,3,1.5],[-3,0,1.5],[0,-3]]),5),

sec5=cr(pts1([[-100,-5,1],[100,30,500],[100,-30,1],[0,5,1],[-100,30,500],[-100,-30,1]]),20),

path7=cr(pts1([[0,0],[0,5],[-13,20],[0,20],[-12,20],[0,78],[12,20],[0,64],[13,20],[0,5]]),5),
sec7=cir(245,s=72),

path10=trns([0,0,185],q_rot(["y-5","z-15"],cytz(cr(pts1([[0,0],[130,30,1000],[260,-10,0]]),20)))),
a=arc(300,-30,0,s=20),
sec11=cr([[0,0,5],[a[0].x,a[0].y,5],[a[len(a)-1].x,a[len(a)-1].y,5]],10),
sec12=trns([50*cos(-15),50*sin(-15),0],cr([[0,0,25],[a[0].x,a[0].y,5],[a[len(a)-1].x,a[len(a)-1].y,5]],30)),

path11=ip(p_extrude(sec,path10),l_extrude(m_points(sec12,10),300)),
path12=[for(i=[0:len(path11)-2])if(norm(path11[i]-path11[i+1])>.1)path11[i]],
 
path13=ip(p_extrude(sec,path10),cyl(r=10,h=300,cp=[85*cos(-15),85*sin(-15)])),
path14=[for(i=[0:len(path13)-2])if(norm(path13[i]-path13[i+1])>.1)path13[i]],
b=arc(50,-30,0,s=10),
sec15=cr([[0,0,5],[b[0].x,b[0].y,2],each loop(b,1,len(b)-2),[b[len(b)-1].x,b[len(b)-1].y,2]]),
sec16=cr(pts1([[-100,0],[100,30,500],[100,-30]]),20))
[sec,sec3,sec5,sec7,sec11,sec12,sec15,path7,path10,path12,path14,sec16];
// 0   1    2     3   4       5   6       7   8       9       10    11
a=spoke1a();

module spoke1a(a){
intersection(){
difference(){
intersection(){
p_extrude1(a[2],a[8]);
swp(l_extrude(a[4],300));
}
swp(trns([85*cos(-15),85*sin(-15)],cyl(r=10,h=300)));
swp(l_extrude(a[5],300));}

swp(prism(a[3],a[7]));}

intersection(){
union(){
p_extrudec(a[1],a[9]);
swp_c(ipf(p_extrude(a[0],a[8]),l_extrude(m_points(f_offset(a[4],-2.5),10),300),3));}
swp(prism(a[3],a[7]));}
p_extrudec(a[1],a[10]);
swp_c(flip(ipf(p_extrude(a[0],a[8]),l_extrude(a[6],300),3,1)));
}

for(i=[0:60:300])rotate([0,0,i])
spoke1a(a);

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

//spoke3

function spoke3a()=let(
sec7=cir(245,s=72),
path7=cr(pts1([[0,0],[0,5],[-13,20],[0,20],[-12,20],[0,78],[12,20],[0,64],[13,20],[0,5]]),5),
sec8=cir(45),
path8=cr(pts1([[0,200],[5,0,5],[0,40,5],[-25,15,5],[-24.5,0]]),5),
prism8=q_rot(["z30"],prism(sec8,path8)),
sec15=cr(pts1([[0,0],[2.5,0,2.5],[0,45,2.5],[-5,0,2.5],[0,-45,2.5]]),5),
path15=[[0,0,205],[300,0,205]],
prism9=q_rot(["z50"],prism(sec7,path7)),
prism1=p_extrude(m_points(sec15,1),path15),
prism2=ipf(prism8,prism1,1,1),
prism3=ipf(prism9,flip(prism1),1),
)[sec15,path15,prism9,prism2,prism3];

module spoke3a(a){
intersection(){
p_extrude(a[0],a[1]);
swp(a[2]);
}
swp_c(flip(a[3]));
swp_c(a[4]);
}

c=spoke3a();
for(i=[0:30:330])rotate([0,0,i])
spoke3a(c);

//hub and rim
function hub_rim()=let(
sec7=cir(245,s=144),

path7=cr(pts1([[0,0],[0,5],[-13,20],[0,20],[-12,20],[0,78],[12,20],[0,64],[13,20],[0,5]]),5),
sec8=cir(45),
path8=cr(pts1([[0,200],[5,0,5],[0,40,5],[-25,15,5],[-24.5,0]]),5),
prism8=q_rot(["z30"],prism(sec8,path8)),

sec9=[each path7,each flip([for(p=path7)p+[5,0]])],
prism9=q_rot(["z30"],prism(sec7,path7)))
[prism8,sec9];

module hub_rim(d){
swp(d[0]);
rotate_extrude($fn=72) 
translate([245,0,0])polygon(d[1]);}

d=hub_rim();

hub_rim(d);
