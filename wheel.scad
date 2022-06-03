include<dependencies.scad>

function spoke1a()=let(

sec=cr(pts1([[-100,0],[100,30,500],[100,-30]]),20),
sec3=cr(pts1([[3,0],[0,3,1.5],[-3,0,1.5],[0,-3]]),5),

path7=cr(pts1([[0,0],[0,5],[-13,20],[0,20],[-12,20],[0,78],[12,20],[0,64],[13,20],[0,5]]),5),
sec7=circle(245,s=72),

path10=trns([0,0,185],q_rot(["y-5","z-15"],cytz(cr(pts1([[0,0],[130,30,1000],[260,-10,0]]),20)))),
a=arc(300,-30,0,s=20),
sec11=cr([[0,0,10],[a[0].x,a[0].y,10],[a[len(a)-1].x,a[len(a)-1].y,10]],10),
sec12=trns([50*cos(-15),50*sin(-15),0],cr([[0,0,25],[a[0].x,a[0].y,5],[a[len(a)-1].x,a[len(a)-1].y,5]],30)),

path11=ip(surf_extrude(sec,path10),l_extrude(m_points(sec12,10),300)),

path13=ip(surf_extrude(sec,path10),cylinder(r=10,h=300,cp=[85*cos(-15),85*sin(-15)])),
prism=ipe(trns([0,0,-5],surf_extrude(sec,path10)),l_extrude(circle(10,[85*cos(-15),85*sin(-15)]),300),2,1),
prism1=[c2t3(c3t2(prism[0])),each prism],
prism2=ipe(trns([0,0,-5],surf_extrude(sec,path10)),l_extrude(m_points(sec12,10),300),2,1),
prism3=[c2t3(c3t2(prism2[0])),each prism2],
b=arc(50,-30,0,s=10),
sec15=cr([[0,0,5],[b[0].x,b[0].y,2],each loop(b,1,len(b)-2),[b[len(b)-1].x,b[len(b)-1].y,2]]),
sec16=cr(pts1([[-100,0],[100,30,500],[100,-30]]),20))
[sec,sec3,sec7,sec11,sec12,sec15,path7,path10,path11,path13,sec16,prism1, prism3];
// 0   1    2     3   4       5   6       7   8       9       10    11  12    
a=spoke1a();

module spoke1a(a){
intersection(){
difference(){
intersection(){
surf_extrude(a[0],a[7],-5);
swp(l_extrude(a[3],300));
}
swp(trns([85*cos(-15),85*sin(-15)],cylinder(r=10,h=300)));
swp(l_extrude(a[4],300));
swp(a[11]);
swp(a[12]);
}

swp(prism(a[2],a[6]));}

intersection(){
union(){
p_extrudec(a[1],a[8]);
swp_c(ipf(surf_extrude(a[0],a[7]),l_extrude(m_points(f_offset(a[3],-2.5),10),300),3));
swp_c(ipf(trns([0,0,-5],surf_extrude(a[0],a[7])),flip(l_extrude(m_points(f_offset(a[3],-2.5),10),300)),3,1));
}
swp(prism(a[2],a[6]));}
p_extrudec(a[1],a[9]);
swp_c(ipf(surf_extrude(a[0],a[7]),l_extrude(a[5],300),3,1));
swp_c(ipf(trns([0,0,-5],surf_extrude(a[0],a[7])),flip(l_extrude(a[5],300)),3));
}


//spoke 2

function spoke2a()=let(
sec=cr(pts1([[-100,0],[100,30,500],[100,-30]]),20),
path=trns([0,0,217],q_rot(["y10","z15"],cytz(cr(pts1([[0,0],[130,30,1000],[260,-10,0]]),20)))),

a=arc(300,0,30,s=20),
sec1=cr([[0,0,10],[a[0].x,a[0].y,10],[a[len(a)-1].x,a[len(a)-1].y,10]],10),

sec2=trns([50*cos(15),50*sin(15),0],cr([[0,0,10],[a[0].x,a[0].y,5],[a[len(a)-1].x,a[len(a)-1].y,5]],20)),

path3=ip(p_extrude(sec,path),l_extrude(m_points(sec2,10),300)),
prism=ipe(trns([0,0,-5],surf_extrude(sec,path)),l_extrude(m_points(sec2,10),300),2,1),
prism1=[c2t3(c3t2(prism[0])),each prism],
sec3=cr(pts1([[3,0],[0,3,1.5],[-3,0,1.5],[0,-3]]),5),

path4=[for(i=[0:len(path3)-2])if(norm(path3[i]-path3[i+1])>.1)path3[i]],

sec5=cr(pts1([[-100,-5,1],[100,30,500],[100,-30,1],[0,5,1],[-100,30,500],[-100,-30,1]]),20),

sec6=cr([[0,0],each arc(50,0,30,s=10)],5),
 
path7=cr(pts1([[0,0],[0,5],[-13,20],[0,20],[-12,20],[0,78],[12,20],[0,64],[13,20],[0,5]]),5),
sec7=circle(245,s=72),
b=arc(50,0,30,s=10),
sec10=cr([[0,0,5],[b[0].x,b[0].y,2],each loop(b,1,len(b)-2),[b[len(b)-1].x,b[len(b)-1].y,5]],5))
[sec,sec1,sec2,sec3,sec5,sec6,sec7,path,path3,path4,path7,sec10, prism1];
//0    1     2   3   4    5      6   7   8       9   10   11        12

b=spoke2a();

module spoke2a(b){
difference(){
intersection(){
surf_extrude(b[0],b[7],-5);
swp(l_extrude(b[1],300));
swp(prism(b[6],b[10]));
}
swp(l_extrude(b[2],300));
swp(b[12]);
}

intersection(){
p_extrudec(b[3],b[9]);
swp(prism(b[6],b[10]));}

intersection(){
swp_c(ipf(surf_extrude(b[0],b[7]),l_extrude(m_points(f_offset(b[1],-2.5),10),300),3));
swp(prism(b[6],b[10]));}

intersection(){
swp_c(ipf(trns([0,0,-5],surf_extrude(b[0],b[7])),flip(l_extrude(m_points(f_offset(b[1],-2.5),10),300)),3,1));
swp(prism(b[6],b[10]));}

swp_c(ipf(surf_extrude(b[0],b[7]),l_extrude(b[11],300),3,1));
swp_c(ipf(trns([0,0,-5],surf_extrude(b[0],b[7])),flip(l_extrude(b[11],300)),3));
}

//spoke3

function s3()=let(
sec1=circle(245,s=72),
path1=cr(pts1([[0,0],[0,5,5],[-13,20,5],[0,20,5],[-12,20,5],[0,78,5],[12,20,5],[0,64,5],[13,20,5],[0,5]]),5),
prism1=q_rot(["z50"],prism(sec1,path1)),

sec2=circle(45),
path2=cr(pts1([[0,200],[5,0,5],[0,40,5],[-25,15,5],[-24.5,0]]),5),
prism2=q_rot(["z30"],prism(sec2,path2)),

sec3=cr(pts1([[0,0],[2.5,0,2.5],[0,45,2.5],[-5,0,2.5],[0,-45,2.5]]),5),
path3=[[0,0],[0,300]],
prism3=trns([20,0,205],q_rot(["x90","z90"],prism(sec3,path3,20))),

prism4=ipf(prism1,flip(prism3),2),
p5=ip(prism2,prism3),
p6=sort_points(p5,ip(q_rot(["z30"],prism(sec2,path_offset(path2,2))),prism3)),
p7=sort_points(p5,ip(prism2,trns([20,0,205],q_rot(["x90","z90"],prism(f_offset(sec3,2),path3,20))))),
p8=[for(i=[0:len(p5)-1])each i<len(p5)-1?[3p_3d_fillet(p6[i],p5[i],p7[i],2)]:[3p_3d_fillet(p6[i],p5[i],p7[i],2),3p_3d_fillet(p6[0],p5[0],p7[0],2)]]
)[prism4,p8,prism3,prism1];

c=s3();

module s3(data){
swp_c(c[1]);
swp_c(c[0]);
intersection(){
swp(c[2]);
swp(c[3]);
}}



//hub and rim

let(

sec=cr(pts1([[0,0],[0,5],[-13,20,5],[0,20,5],[-12,20,5],[0,78,5],[12,20,5],[0,64,5],[13,20,3],[0,5]]),5),
sec1=concat(sec,c3t2(trns([5,0,0],flip(sec)))),

sec2=cr(pts1([[0,197],[50,0,5],[0,43,5],[-25,15,5],[-25,0]]),5),

){
difference(){
rotate_extrude($fn=72) polygon(sec2);
cylinder(r=20,h=300);}
rotate_extrude($fn=144) translate([245,0,0])polygon(sec1);
}
render(){
difference(){
for(i=[0:60:300])rotate([0,0,i]){
//scale(1.002)
spoke1a(a);
spoke2a(b);}
cylinder(r=20,h=300);}

for(i=[0:30:330])rotate([0,0,i])
s3(c);}
