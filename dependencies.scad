

function prism(sec,path,m_points=100)=[for(p=path)[for(p1=sort_points(m_points(sec,m_points),m_points(f_offset(sec,round(p.x*100)/100),m_points)))[p1.x,p1.y,p.y]]];
    
function surf(sec,path)=[for(p=path)[for(p1=sec)[p.x,p1.y,p.y]]];
             
function ang(x,y)= x>=0&&y>=0?atan(y/x):x<0&&y>=0?180-abs(atan(y/x)):x<0&&y<0?180+abs(atan(y/x)):360-abs(atan(y/x));
            
function ang3d(v1,v2)=let(  
    u1=v1/norm(v1),
    u2=v2/norm(v2))asin(norm(u2-u1)/2);

            
function q(vector=[1,0,0],point=[0,5,0],theta=0)=

let(t=theta,
v=vector/norm(vector),
p=[cos(t/2),v*sin(t/2)],
p1=[p.x,-p.y],
q=[0,len(point)==2?[point.x,point.y,0]:point],
pq=[p.x*q.x-p.y*q.y,p.x*q.y+p.y*q.x+cross(p.y,q.y)],
pqp1=[pq.x*p1.x-pq.y*p1.y,pq.x*p1.y+pq.y*p1.x+cross(pq.y,p1.y)],
transformation=pqp1.y
)
//assert(!is_undef(transformation),str(v,theta,p.y,q.y))
transformation
;

function qmr1(s,r,pl,n=0)= n==len(s)?pl:
qmr1(s,r,
    let(
    v1=s[n]=="x"?[1,0,0]:s[n]=="y"?[0,1,0]:[0,0,1],
    r1=r[n]==undef?0:r[n])
    [for(p=pl)q(v1,p,r1)],n+1);

function qmr2(s,r,pl,n=0)= n==len(s)?pl:qmr2(s,r,let(
    v1=s[n]=="x"?[1,0,0]:s[n]=="y"?[0,1,0]:[0,0,1],
    r1=r[n]==undef?0:r[n])
[for(i=[0:len(pl)-1])[for(p=pl[i])q(v1,p,r1)]],n+1);
    
function q_rot(s,pl)= is_num(pl[0][0])?qmr1([for(p=s)cvar(p)[0]],[for(p=s)cvar(p)[1]],pl):qmr2([for(p=s)cvar(p)[0]],[for(p=s)cvar(p)[1]],pl);

function sort(list,n=0)=
let(
list1=[for(i=[0:len(list)-1])[list[i]+i*.0000000001,i]],
a=lookup(min(list1*[1,0]),list1),
list2=[for(i=[0:len(list1)-1])if (lookup(list1[i].x,list1)!=a)list1[i]]
)n==0?min(list1*[1,0]):sort(list2*[1,0],n-1);

function sort_list(list)=[for(i=[0:len(list)-1])sort(list,i)];
    
function sortv(vector_list,m=.2)=
let(
list=[for(i=[0:len(vector_list)-2])[norm(vector_list[i+1]-vector_list[i]),i]],
list1=[for(i=[0:len(list)-1])if(list[i].x>=m)list[i].y]    

)[for(i=list1)vector_list[i]];
    

module sec_hull(sec,t=.01,options=0) {
for(i=[0:len(sec)-2])
    for(j=[0:len(sec[i])-2])
        let(
    i_plus=i+1,
    j_plus=j+1,
    p=options==1?[0,0,t/2*1.1]:options==2?[0,0,-t/2*1.1]:options==3?[0,t/2*1.1,0]:options==4?[0,-t/2*1.1,0]:options==5?[t/2*1.1,0,0]:options==6?[-t/2*1.1,0,0]:[0,0,0]){
       
    hull(){
    translate(sec[i][j])translate(p)cube(t,true);
    translate(sec[i][j_plus])translate(p)cube(t,true);
    translate(sec[i_plus][j])translate(p)cube(t,true);}
    
    hull(){
    translate(sec[i_plus][j_plus])translate(p)cube(t,true);
    translate(sec[i][j_plus])translate(p)cube(t,true);
    translate(sec[i_plus][j])translate(p)cube(t,true);}}
    
    }

module sec_hull_c(sec,t=.01) {
for(i=[0:len(sec)-1])
    for(j=[0:len(sec[i])-1])
        let(
    i_plus=i<len(sec)-1?i+1:0,
    j_plus=j<len(sec[i])-1?j+1:0,
    p=options==1?[0,0,t/2]:options==2?[0,0,-t/2]:options==3?[0,t/2,0]:options==4?[0,-t/2,0]:options==5?[t/2,0,0]:options==6?[-t/2,0,0]:[0,0,0]){
       
    hull(){
    translate(sec[i][j])translate(p)cube(t,true);
    translate(sec[i][j_plus])translate(p)cube(t,true);
    translate(sec[i_plus][j])translate(p)cube(t,true);}
    
    hull(){
    translate(sec[i_plus][j_plus])translate(p)cube(t,true);
    translate(sec[i][j_plus])translate(p)cube(t,true);
    translate(sec[i_plus][j])translate(p)cube(t,true);}}
    
    }
    
function surf_extrude(sec,path)=[

    for(i=[0:len(path)-2])
    let(
    p0=path[i],
    p1=path[i+1],
    v=p1-p0,
    a1=ang(v.x,v.y),
    a2=ang(sqrt(v.x^2+v.y^2),v.z),
    sec1=trns(p0,q_rot(["x90","z-90",str("y",-a2),str("z",a1)],sec)),
    sec2=trns(p1,q_rot(["x90","z-90",str("y",-a2),str("z",a1)],sec))

    )i<len(path)-2?sec1:sec2];
    
 module surf_extrude(sec,path,t=.01,o=1){
     if(o==1){
     surf=surf_extrude(sec,path);
     surf1=trns([0,0,t],surf);
     for(i=[0:len(surf)-2])
         for(j=[0:len(surf[i])-2])
           if(t>0)
             swp([[surf1[i][j],surf1[i+1][j],surf1[i+1][j+1],surf1[i][j+1]],[surf[i][j],surf[i+1][j],surf[i+1][j+1],surf[i][j+1]]]);
         else
             swp([[surf[i][j],surf[i+1][j],surf[i+1][j+1],surf[i][j+1]],[surf1[i][j],surf1[i+1][j],surf1[i+1][j+1],surf1[i][j+1]]]);}
             else if(o==2){ 
                 s1=prism(sec,path);
                 s2=prism(f_offset(sec,t),path);
                 for(i=[0:len(s1)-2])
                     for(j=[0:len(s2[i])-2])
                 if(t>0)       
                 swp(flip([[s1[i][j],s2[i][j],s2[i][j+1],s1[i][j+1]],[s1[i+1][j],s2[i+1][j],s2[i+1][j+1],s1[i+1][j+1]]]));
                  else
                 swp(flip([[s1[i+1][j],s2[i+1][j],s2[i+1][j+1],s1[i+1][j+1]],[s1[i][j],s2[i][j],s2[i][j+1],s1[i][j+1]]]));
                 
                 }
     }
    
 module surf_extrudec(sec,path,t=.01,o=1){
     if(o==1){
     surf=surf_extrude(sec,path);
     surf1=trns([0,0,t],surf);
     for(i=[0:len(surf)-2])
         for(j=[0:len(surf[i])-2])
           if(t>0)
             swp([[surf1[i][j],surf1[i+1][j],surf1[i+1][j+1],surf1[i][j+1]],[surf[i][j],surf[i+1][j],surf[i+1][j+1],surf[i][j+1]]]);
         else
             swp([[surf[i][j],surf[i+1][j],surf[i+1][j+1],surf[i][j+1]],[surf1[i][j],surf1[i+1][j],surf1[i+1][j+1],surf1[i][j+1]]]);}
             else if(o==2){ 
                 s1=prism(sec,path);
                 s2=prism(f_offset(sec,t),path);
                 for(i=[0:len(s1)-2])
                     for(j=[0:len(s1[i])-1])let(j_plus=j<len(s1[i])-1?j+1:0)
                 if(t>0)       
                 swp([[s1[i][j],s2[i][j],s2[i][j_plus],s1[i][j_plus]],[s1[i+1][j],s2[i+1][j],s2[i+1][j_plus],s1[i+1][j_plus]]]);
                  else
                 swp([[s1[i+1][j],s2[i+1][j],s2[i+1][j_plus],s1[i+1][j_plus]],[s1[i][j],s2[i][j],s2[i][j_plus],s1[i][j_plus]]]);
                 
                 }
     }
     
//function avg_v(vector)=
//let( 
// x=len(vector[0])==3?vector*[1,0,0]:vector*[1,0],
// mx=[for(i=[0:len(x)-1])1],
//avg_x=(x*mx)/len(x),
// 
// y=len(vector[0])==3?vector*[0,1,0]:vector*[0,1],
// my=[for(i=[0:len(y)-1])1],
//avg_y=(y*my)/len(y),
// 
// z=len(vector[0])==3?vector*[0,0,1]:vector*[0,0],
// mz=[for(i=[0:len(z)-1])1],
//avg_z=(z*mz)/len(z)
// )len(vector[0])==3?[avg_x,avg_y,avg_z]:[avg_x,avg_y];
 
// function avg_v3d(vector)=let(
//v=[for(i=[0:len(vector)-1])avg_v(vector[i])],
//max_z=max(v*[0,0,1]),
//min_z=min(v*[0,0,1]))[avg_v(v).x,avg_v(v).y,(max_z-min_z)/2];


function near(sec,p)=let(
    a=[for(i=[0:len(sec)-1])[norm(p-[sec[i].x,sec[i].y]),i]],
    b=min(a*[1,0]),
    c=lookup(b,a))sec[c];

    
function cytz(path)=[for(p=path)[p.x,0,p.y]];
    
function sort_p(sec,path)=[for(p=sec)near(path,p)];
    
function shift_s(sec,dist)=[for(i=[0:len(sec)-1])
    [for(p=sec[i])p+dist]];
        
function cir(r,cp=[0,0],s=50)=[for(i=[0:360/s:360-360/s])[cp.x+r*cos(i),cp.y+r*sin(i)]];

module p_line(path,size=.5){
    for(i=[0:len(path)-1])
        let(p0=path[i],p1=i<len(path)-1?path[i+1]:path[0])
    
    hull(){
    translate(p0)circle(size/2,$fn=20);
    translate(p1)circle(size/2,$fn=20);}}
    
module p_lineo(path,size=.5){
    for(i=[0:len(path)-2])
        let(p0=path[i],p1=path[i+1])
    
    hull(){
    translate(p0)circle(size/2,$fn=20);
    translate(p1)circle(size/2,$fn=20);}}
    
module rd_line(path,size=.5){
    for(i=[0:len(path)-1])
        let(p0=path[i],p1=i<len(path)-1?path[i+1]:path[0])
    
    hull(){
    translate(p0)sphere(size,true,$fn=30);
    translate(p1)sphere(size,true,$fn=30);}}
    
function list_ang(sec)=[for(i=[0:len(sec)-1])
    let(
p0=sec[i],p1=i<len(sec)-1?sec[i+1]:sec[0],p2=i<len(sec)-2?sec[i+2]:i<len(sec)-1?sec[0]:sec[1],
v1=p1-p0, v2=p2-p1,
angle1=ang(v1.x,v1.y),angle2=ang(v2.x,v2.y),
angle=angle2-angle1

)if(is_num(angle))angle];

function cw(sec)=let(p=mode_sign(list_ang(sec)))
p[0]>p[1]?-1:1;

function i_p2d(l1,l2)=let(
p0=l1[0],p1=l1[1],
p2=l2[0],p3=l2[1],
v1=p1-p0,
v2=p3-p2,
//p0+v1*t1=p2+v2*t2
//v1*t1-v2*t2=p2-p0
//[[v1.x,-v2.x],[v1.y,-v2.y]]*[t1,t2]=[p2.x-p0.x,p2.y-p0.y]
t1=(i_m2d([[v1.x,-v2.x],[v1.y,-v2.y]])*[p2.x-p0.x,p2.y-p0.y])[0],
pi=p0+v1*t1

)pi;

function i_p3d(l1,l2)=
let(

v1=l1[1]-l1[0],u1=v1/norm(v1),
v2=l2[1]-l2[0],u2=v2/norm(v2),
v3=l2[0]-l1[0],
//l1[0]+v1*t1=l2[0]+v2*t2
//v1*t1-v2*t2=v3 where v3=l2[0]-l1[0]
t1=(i_m3d([[v1.x,-v2.x,1],
           [v1.y,-v2.y,1],
           [v1.z,-v2.z,1]])*
      [v3.x,v3.y,v3.z])[0],
ip=l1[0]+v1*t1
)ip;

function mode_sign(p,plus=0,minus=0,n=0)=n==len(p)?[plus,minus]:mode_sign(p,p[n]>0?plus+1:plus+0,p[n]<0?minus+1:minus+0,n+1);

module poly_h(prism)
let(
sec1=[for(i=[0:len(prism)-1])each prism[i]],
points_list=[for(i=[0:len(prism)-1])
if(i==0)[for(j=[len(prism[0])-1:-1:0])j]
    else if(i==len(prism)-1)[for(j=[0:len(prism[0])-1])j+i*len(prism[0])]    

],

intermediate_points_list=[for(i=[0:len(prism)-2])
    for(j=[0:len(prism[0])-1])
    j<len(prism[0])-1?
    [
    j+1+i*len(prism[0])
    ,j+1+(i+1)*len(prism[0])
    ,j+(i+1)*len(prism[0])
    ,j+i*len(prism[0])  
    ]:
    [
    j-len(prism[0])+1+(i)*len(prism[0])
    ,j-len(prism[0])+1+(i+1)*len(prism[0])
    ,j+(i+1)*len(prism[0])
    ,j+i*len(prism[0])
    ]
    ])
    polyhedron(sec1,[each points_list,each intermediate_points_list],convexity=10);
    
 module poly_h1(prism)
let(
sec1=[for(i=[0:len(prism)-1])each prism[i]],
points_list=[for(i=[0:len(prism)-1])
if(i==0)[for(j=[0:len(prism[0])-1])j]
    else if(i==len(prism)-1)[for(j=[len(prism[0])-1:-1:0])j+i*len(prism[0])]    

],

intermediate_points_list=[for(i=[0:len(prism)-2])
    for(j=[0:len(prism[0])-1])
    j<len(prism[0])-1?
    [
    j+i*len(prism[0]) 
    ,j+(i+1)*len(prism[0])
    ,j+1+(i+1)*len(prism[0])
    ,j+1+i*len(prism[0])
    
    
     
    ]:
    [
    j+i*len(prism[0])
    ,j+(i+1)*len(prism[0])
    ,j-len(prism[0])+1+(i+1)*len(prism[0])
    ,j-len(prism[0])+1+(i)*len(prism[0])
    
    
    
    ]
    ])
    polyhedron(sec1,[each points_list,each intermediate_points_list],convexity=10);
    
function arc(radius,ang1=0,ang2=355,cp=[0,0],s=20)=[for(i=[ang1:(ang2-ang1)/s:ang2])cp+[radius*cos(i),radius*sin(i)]];
    

    
function 2cir_fillet(r1=10,r2=10,c1=[0,0],c2=[20,0],r=10)=
let(
l1=norm(c2-c1),l2=r1+r,l3=r2+r,
t=(l1^2+l2^2-l3^2)/(2*l1),
h=sqrt(l2^2-t^2),
v=c2-c1,u=v/norm(v),
p1=c1+u*t+u*[[0,1],[-1,0]]*h,
a1=ang((c1-p1).x,(c1-p1).y),
a2=ang((c2-p1).x,(c2-p1).y),
p2=c1+u*t+u*[[0,-1],[1,0]]*h,
a3=ang((c2-p2).x,(c2-p2).y),
a4=ang((c1-p2).x,(c1-p2).y),
a5=ang((p1-c1).x,(p1-c1).y),
a6=ang((p2-c1).x,(p2-c1).y),
a7=ang((p1-c2).x,(p1-c2).y),
a8=ang((p2-c2).x,(p2-c2).y),


arc1=arc(r,a2<a1?360+a2:a2,a1,p1),
arc2=arc(r,a4<a3?360+a4:a4,a3,p2),
arc3=arc(r2,a7<a8?a7+360:a7,a8,c2),
arc4=arc(r1,a5,a6<a5?a6+360:a6,c1)
)
concat(arc2,arc1);

function 2cir_filleto(r1=10,r2=10,c1=[0,0],c2=[20,0],r=10)=
let(
l1=norm(c2-c1),l2=r1+r,l3=r2+r,
t=(l1^2+l2^2-l3^2)/(2*l1),
h=sqrt(l2^2-t^2),
v=c2-c1,u=v/norm(v),
p1=c1+u*t+u*[[0,1],[-1,0]]*h,
a1=ang((c1-p1).x,(c1-p1).y),
a2=ang((c2-p1).x,(c2-p1).y),
p2=c1+u*t+u*[[0,-1],[1,0]]*h,
a3=ang((c2-p2).x,(c2-p2).y),
a4=ang((c1-p2).x,(c1-p2).y),
a5=ang((p1-c1).x,(p1-c1).y),
a6=ang((p2-c1).x,(p2-c1).y),
a7=ang((p1-c2).x,(p1-c2).y),
a8=ang((p2-c2).x,(p2-c2).y),


arc1=arc(r,a2<a1?360+a2:a2,a1,p1),
arc2=arc(r,a4<a3?360+a4:a4,a3,p2),
arc3=arc(r2,a7<a8?a7+360:a7,a8,c2),
arc4=arc(r1,a5,a6<a5?a6+360:a6,c1)
)
[arc2,arc1];

function rm(theta)=[[cos(theta),sin(theta)],[-sin(theta),cos(theta)]];

function 2df(p1,p2,p3,r0,r1,r2,theta0,theta1,theta2,u2,u3,s)=let(
l1=norm(p1-p2),
l2=r0*tan(theta0)+r1*tan(theta1),
l3=norm(p3-p2),
l4=r1*tan(theta1)+r2*tan(theta2),
rf1=l1>l2?r1:l1/l2*r1,
rf2=l3>l4?r1:l3/l4*r1,
rf=min(rf1,rf2),

p=p2+u2*rf*tan(theta1),
cp=cw([p1,p2,p3])==-1?p-u2*rm(90)*rf:p-u2*rm(-90)*rf,
a1=ang((p-cp).x,(p-cp).y),
a2=cw([p1,p2,p3])==-1?a1+2*theta1:a1-2*theta1,
arc=arc(rf,a1,a2,cp,s)

)
r1==0 || r1==undef||norm(u2-u3)<.2?[p2]:arc;
    


function 2dfillet(pl,rl,s)=[for(i=[0:len(pl)-1])let(ep=[.0001,.0001],
p0=i==0?pl[len(pl)-2]:i==1?pl[len(pl)-1]:pl[i-2],
p1=i==0?pl[len(pl)-1]:pl[i-1],
p2=pl[i],
p3=i<len(pl)-1?pl[i+1]:pl[0],
p4=i<len(pl)-2?pl[i+2]:i<len(pl)-1?pl[0]:pl[1],
r0=i==0?rl[len(rl)-1]:rl[i-1],
r1=rl[i],
r2=i<len(rl)-1?rl[i+1]:rl[0],
v0=p0-p1,u0=v0/norm(v0),
v1=p2-p1,u1=v1/norm(v1),
v2=p1-p2,u2=v2/norm(v2),
v3=p3-p2,u3=v3/norm(v3),
v4=p2-p3,u4=v4/norm(v4),
v5=p4-p3,u5=v5/norm(v5),
ang0=ang(u0.x,u0.y),
ang1=ang(u1.x,u1.y),
ang2=ang(u2.x,u2.y),
ang3=ang(u3.x,u3.y),
ang4=ang(u4.x,u4.y),
ang5=ang(u5.x,u5.y),


theta0=abs(180-((ang0<ang1?ang0+360:ang0)-ang1))/2,
theta1=abs(180-((ang2<ang3?ang2+360:ang2)-ang3))/2,
theta2=abs(180-((ang4<ang5?ang4+360:ang4)-ang5))/2        


)each 
2df(p1,p2,p3,r0,r1,r2,theta0,theta1,theta2,u2,u3,s)];


function cr1(pl,s=20)=let(
pl1=[for(i=[0:len(pl)-1])[pl[i].x,pl[i].y]],
rl=[for(i=[0:len(pl)-1])pl[i].z==undef?0:pl[i].z]

)2dfillet(pl1,rl,s);

function cr(pl,s=20)=let(
sec=cr1(pl,s),
sec1=[for(i=[0:len(sec)-1])if(norm(sec[i<len(sec)-1?i+1:0]-sec[i])>.01)sec[i]]
//r=min_r(sec1),
//sec01=[for(i=[0:len(sec1)-1])
//let(
//p0=sec1[i],
//p1=i<len(sec1)-1?sec1[i+1]:sec1[0],
//list=norm(p1-p0)>abs(r)?l([p0,p1],round(norm(p1-p0)/r)):[p0],
//)each list    
//],
//sec0=[for(i=[0:len(sec01)-1])if(i<len(sec01)-1&&norm(sec01[i]-sec01[i+1])>.1)sec01[i] else sec01[i]],

)sec1;

function tr1(tm,sec)=[for(p=sec)(len(tm)==2?[tm.x,tm.y,0]:tm) + (len(p)==2?[p.x,p.y,0]:p)];
    
function tr2(tm,sec)=[for(i=[0:len(sec)-1])[for(p=sec[i])(len(tm)==2?[tm.x,tm.y,0]:tm) + (len(p)==2?[p.x,p.y,0]:p)]];
    
function trns(m,sec)=is_num(sec[0][0])?tr1(m,sec):tr2(m,sec);

function scl2d(sec,sl)=let(
sec1=len(sec[0])==2?[for(p=sec)[p.x,p.y,0]]:sec,
xmin=min(sec1*[1,0,0]),
xmax=max(sec1*[1,0,0]),
ymin=min(sec1*[0,1,0]),
ymax=max(sec1*[0,1,0]),

cp=[xmin,ymin]+[(xmax-xmin)/2,(ymax-ymin)/2,0],
x_rev=[for(p=sec1)let(
p0=[cp.x,0,0],p1=[p.x,0,0],
p2=p0+(p1-p0)*sl
)p2],
y_rev=[for(p=sec1)let(
p0=[0,cp.y,0],p1=[0,p.y,0],
p2=p0+(p1-p0)*sl
)p2],
rev_sec=[for(i=[0:len(sec1)-1])x_rev[i]+y_rev[i]],
len1=norm([0,cp.y,0]-[0,ymin,0]),
len2=len1*sl,
len3=len2-len1,
rev_pl=[for(p=rev_sec)p+[0,len3,0]]

)rev_pl;

function scl3d(sec,sl)=[for(i=[0:len(sec)-1])let(
sec0=len(sec[0][0])==2?[for(p=sec[0])[p.x,p.y,0]]:sec[0],
x0min=min(sec0*[1,0,0]),
x0max=max(sec0*[1,0,0]),
y0min=min(sec0*[0,1,0]),
y0max=max(sec0*[0,1,0]),
z0min=min(sec0*[0,0,1]),
z0max=max(sec0*[0,0,1]),   

cp0=[x0min,y0min,z0min]+[(x0max-x0min)/2,(y0max-y0min)/2,(z0max-z0min)/2],

sec1=len(sec[i][0])==2?[for(p=sec[i])[p.x,p.y,0]]:sec[i],
   
xmin=min(sec1*[1,0,0]),
xmax=max(sec1*[1,0,0]),
ymin=min(sec1*[0,1,0]),
ymax=max(sec1*[0,1,0]),
zmin=min(sec1*[0,0,1]),
zmax=max(sec1*[0,0,1]),

cp=[xmin,ymin,zmin]+[(xmax-xmin)/2,(ymax-ymin)/2,(zmax-zmin)/2],

x_rev=[for(p=sec1)let(
p0=[cp.x,0,0],p1=[p.x,0,0],
p2=p0+(p1-p0)*sl
)p2],

y_rev=[for(p=sec1)let(
p0=[0,cp.y,0],p1=[0,p.y,0],
p2=p0+(p1-p0)*sl
)p2],

z_rev=[for(p=sec1)let(
p0=[0,0,cp.z],p1=[0,0,p.z],
p2=p0+(p1-p0)*sl
)p2],

rev_sec=[for(i=[0:len(sec1)-1])x_rev[i]+y_rev[i]+z_rev[i]],
len1=norm([0,0,cp.z]-[0,0,zmin]),
len2=len1*sl,
len3=len2-len1,
rev_pl=[for(p=rev_sec)p+[0,0,len3]],
vector=cp0-cp,
vector1=-vector*sl
)trns(vector1,trns(vector,rev_pl))];

function ipa(prism,prism1)=
[for(i=[0:len(prism1)-2])
    for(j=[0:len(prism1[i])-1])
        for(k=[0:len(prism)-2])
            for(l=[0:len(prism[k])-2])
                let(
            k_plus=k+1,l_plus=l<len(prism[k])-1?l+1:0,
            pa=prism[k][l],pb=prism[k][l_plus],pc=prism[k_plus][l],pd=prism[k_plus][l_plus],
            p0=prism1[i][j],p1=prism1[i+1][j],
            v1=p1-p0,v2=pb-pa,v3=pc-pa,
            v4=p1-p0,v5=pb-pd,v6=pc-pd,
//            p0+v1*t1=pa+v2*t2+v3*t3
//            p0-pa=-v1*t1+v2*t2+v3*t3           
            t1=(p0-pa)*cross(v2,v3)/(-v1*cross(v2,v3)),
            t2=(p0-pa)*cross(v3,-v1)/(-v1*cross(v2,v3)),
            t3=(p0-pa)*cross(-v1,v2)/(-v1*cross(v2,v3)),
            t4=(p0-pd)*cross(v5,v6)/(-v4*cross(v5,v6)),
            t5=(p0-pd)*cross(v6,-v4)/(-v4*cross(v5,v6)),
            t6=(p0-pd)*cross(-v4,v5)/(-v4*cross(v5,v6))
            
            )if(lim(t1,0,1)&&lim(t2,0,1)&&lim(t3,0,1)&&lim(t2+t3,0,1))p0+v1*t1
            else if(lim(t4,0,1)&&lim(t5,0,1)&&lim(t6,0,1)&&lim(t5+t6,0,1))p0+v4*t4];
            
function ip(prism,prism1)=let(sec=ipa(prism,prism1))[for(i=[0:len(sec)-1])let(i_plus=i<len(sec)-1?i+1:0)if(norm(sec[i]-sec[i_plus])>.1)sec[i]];
            
function ip1(prism,prism1)=
[for(i=[0:len(prism1)-2])
    
        for(k=[0:len(prism)-2])
            for(l=[0:len(prism[k])-1])
                let(
            k_plus=k+1,l_plus=l<len(prism[k])-1?l+1:0,
            pa=prism[k][l],pb=prism[k][l_plus],pc=prism[k_plus][l],pd=prism[k_plus][l_plus],
            p0=prism1[i],p1=prism1[i+1],
            
            v1=p1-p0,v2=pb-pa,v3=pc-pa,
            v4=p1-p0,v5=pb-pd,v6=pc-pd,
//            p0+v1*t1=pa+v2*t2+v3*t3
//            p0-pa=-v1*t1+v2*t2+v3*t3           
            t1=(p0-pa)*cross(v2,v3)/(-v1*cross(v2,v3)),
            t2=(p0-pa)*cross(v3,-v1)/(-v1*cross(v2,v3)),
            t3=(p0-pa)*cross(-v1,v2)/(-v1*cross(v2,v3)),
            t4=(p0-pd)*cross(v5,v6)/(-v4*cross(v5,v6)),
            t5=(p0-pd)*cross(v6,-v4)/(-v4*cross(v5,v6)),
            t6=(p0-pd)*cross(-v4,v5)/(-v4*cross(v5,v6))

            )if(lim(t1,0,1)&&lim(t2,0,1)&&lim(t3,0,1)&&lim(t2+t3,0,1))p0+v1*t1
            else if(lim(t4,0,1)&&lim(t5,0,1)&&lim(t6,0,1)&&lim(t5+t6,0,1))p0+v4*t4];
            
function ip2(prism,prism1)=
[for(i=[0:len(prism1)-2])
    
        
                let(
            pa=prism[0],pb=prism[1],pc=prism[2],
            p0=prism1[i],p1=prism1[i+1],
            
            v1=p1-p0,v2=pb-pa,v3=pc-pa,
            //p0+v1*t1=pa+v2*t2+v3*t3
            //p0-pa=-v1*t1+v2*t2+v3*t3
            t1= cross(v2,v3)*(p0-pa)/(-v1*cross(v2,v3)),
            t2=cross(v3,-v1)*(p0-pa)/(-v1*cross(v2,v3)),
            t3=cross(-v1,v2)*(p0-pa)/(-v1*cross(v2,v3))
            
            )if(lim(t1,0,1))p0+v1*t1
];
            
function nv3d(v)=[
v.x==0&&v.y==0&&v.z>0?[-1,0,0]:
v.x==0&&v.y==0&&v.z<0?[1,0,0]:
v.x==0&&v.z==0&&v.y>0?[-1,0,0]:
v.x==0&&v.z==0&&v.y<0?[1,0,0]:
v.y==0&&v.z==0&&v.x>0?[0,1,0]:
v.y==0&&v.z==0&&v.x<0?[0,-1,0]:
let(v1=[v.x,v.y]*[[0,1],[-1,0]])
[v1.x,v1.y,0]

].x;
            
function c3t2(sec)=[for(p=sec)[p.x,p.y]];
function sum_v2d(v,n=0,s=[0,0])=n==len(v)-1?s:sum_v2d(v,n+1,s+v[n]);
//function avg_v2d(sec)=sum_v2d(sec)/len(sec);

function 2ctp(r1,r2,cp1,cp2)=
let(
v1=cp2-cp1,
u1=v1/norm(v1),
ang1=asin((r2-r1)/norm(cp2-cp1)),

t1=cp1+u1*r1*rm(90+ang1),
t2=cp2+u1*r2*rm(90+ang1),

t3=cp1+u1*r1*rm(-90-ang1),
t4=cp2+u1*r2*rm(-90-ang1))[t1,t2];

function 2ctpf(r1,r2,cp1,cp2)=
let(
v1=cp2-cp1,
u1=v1/norm(v1),
ang1=asin((r2-r1)/norm(cp2-cp1)),

t1=cp1+u1*r1*rm(90+ang1),
t2=cp2+u1*r2*rm(90+ang1),

t3=cp1+u1*r1*rm(-90-ang1),
t4=cp2+u1*r2*rm(-90-ang1))[t1,t2,t4,t3];
    
module p_line3d(path,r,rec=0){
    for(i=[0:len(path)-2])
        
    hull(){
    translate(path[i])if(rec==0)sphere(r); else cube(r*2,true);
    translate(path[i+1])if(rec==0)sphere(r);else cube(r*2,true);
    }}
module p_line3dc(path,r,rec=0){
    for(i=[0:len(path)-1])
        let(
    i_plus=i<len(path)-1?i+1:0
    )
    hull(){
    translate(path[i])if(rec==0)sphere(r); else cube(r*2,true);
    translate(path[i_plus])if(rec==0)sphere(r);else cube(r*2,true);
    }}
    
function ipw(prism,prism1,r)=
[for(i=[0:len(prism1)-2])
    for(j=[0:len(prism1[i])-1])
        for(k=[0:len(prism)-2])
            for(l=[0:len(prism[k])-2])
                let(ep=[.0001,.0001,.0001],
            k_plus=k+1,l_plus=l<len(prism[k])-1?l+1:0,
            pa=prism[k][l],pb=prism[k][l_plus],pc=prism[k_plus][l],pd=prism[k_plus][l_plus],

            p0=prism1[i][j],p1=prism1[i+1][j],
            p2=prism1[i][j+1],p3=prism1[i+1][j+1],
            v1=p1-p0,v2=pb-pa,v3=pc-pa,
            v4=p1-p0,v5=pb-pd,v6=pc-pd,
//            p0+v1*t1=pa+v2*t2+v3*t3
//            p0-pa=-v1*t1+v2*t2+v3*t3           
            t1=(p0-pa)*cross(v2,v3)/(-v1*cross(v2,v3)),
            t2=(p0-pa)*cross(v3,-v1)/(-v1*cross(v2,v3)),
            t3=(p0-pa)*cross(-v1,v2)/(-v1*cross(v2,v3)),
            t4=(p0-pd)*cross(v5,v6)/(-v4*cross(v5,v6)),
            t5=(p0-pd)*cross(v6,-v4)/(-v4*cross(v5,v6)),
            t6=(p0-pd)*cross(-v4,v5)/(-v4*cross(v5,v6))
            )
            if(lim(t1,0,1)&&lim(t2,0,1)&&lim(t3,0,1)&&lim(t2+t3,0,1)) [p0+v1*t1,p0+v1*t1+(p1-p0)/norm(p1-p0)*r,pa,pb,pc]
            else if(lim(t4,0,1)&&lim(t5,0,1)&&lim(t6,0,1)&&lim(t5+t6,0,1))[p0+v4*t4,p0+v4*t4+(p1-p0)/norm(p1-p0)*r,pd,pb,pc]
            ];
            
           
 function ipr(prism,prism1,r,option=0,s=5)=let(list=ipw(prism,prism1,r),
            p1=[for(i=[0:len(list)-1])list[i][0]],
            p2=[for(i=[0:len(list)-1])list[i][1]],
            p3=[for(i=[0:len(list)-1])list[i][2]],
            p4=[for(i=[0:len(list)-1])list[i][3]],
            p5=[for(i=[0:len(list)-1])list[i][4]]
            //p6=[for(i=[0:len(p1)-1])i<len(p1)-1?p1[i+1]:p1[0]]
            
            
            
            )[for(i=[0:len(p1)-1])
            let(i_plus=i<len(p1)-1?i+1:0,
            v1=p1[i_plus]-p1[i],
            //v2=p4[i]-p3[i],u2=v2/norm(v2),
            //v3=p5[i]-p3[i],u3=v3/norm(v3),
            cir=option==0?[for(j=[0:-20:-180])if(norm(v1)>.01) p1[i]+q(v1,p2[i]-p1[i],j)]:[for(j=[0:20:180])if(norm(v1)>.01)p1[i]+q(v1,p2[i]-p1[i],j)],
            p7=norm(v1)>.01?ip2([p3[i],p4[i],p5[i]],cir):[]
            
            
            ) //[v1,p2[i]-p1[i],5]
            if(! is_undef(p7[0]))3p_3d_fillet(p2[i],p1[i],p7[0],r,s)
            
            ];
            
 function ipr1(prism,prism1,r,option=0,s=5)=let(list=ipw(prism,prism1,r),
            p1=[for(i=[0:len(list)-1])list[i][0]],
            p2=[for(i=[0:len(list)-1])list[i][1]],
            p3=[for(i=[0:len(list)-1])list[i][2]],
            p4=[for(i=[0:len(list)-1])list[i][3]],
            p5=[for(i=[0:len(list)-1])list[i][4]]
            //p6=[for(i=[0:len(p1)-1])i<len(p1)-1?p1[i+1]:p1[0]]
            
            
            
            )[for(i=[0:len(p1)-1])
            let(i_plus=i<len(p1)-1?i+1:0,
            v1=p1[i_plus]-p1[i],
            //v2=p4[i]-p3[i],u2=v2/norm(v2),
            //v3=p5[i]-p3[i],u3=v3/norm(v3),
            cir=option==0?[for(j=[0:-20:-180])if(norm(v1)>.01) p1[i]+q(v1,p2[i]-p1[i],j)]:[for(j=[0:20:180])if(norm(v1)>.01)p1[i]+q(v1,p2[i]-p1[i],j)],
            p7=norm(v1)>.01?ip2([p3[i],p4[i],p5[i]],cir):[]
            
            
            ) //[v1,p2[i]-p1[i],5]
            if(! is_undef(p7[0]))3p_3d_fillet_wo_pivot(p2[i],p1[i],p7[0],r,s)
            
            ];
 
 
  function ipf(prism,prism1,r,option=0,s=5)=let(sec=ipr(prism,prism1,r,option,s=s))
            [for(i=[0:len(sec)-1])each i<len(sec)-1?[sec[i]]:[sec[i],sec[0]]];
            
 function ipe(prism,prism1,r,option=0,s=5)=let(sec=ipr1(prism,prism1,r,option,s=s))
            [for(i=[0:len(sec)-1])each i<len(sec)-1?[sec[i]]:[sec[i],sec[0]]];
                
function cyl(r1=1,r2=1,h=1,cp=[0,0],s=50,r,d,d1,d2,center=false)=let(
     ra=is_num(r)?r:is_num(d)?d/2:is_num(d1)?d1/2:r1,
     rb=is_num(r)?r:is_num(d)?d/2:is_num(d2)?d2/2:r2,
     sec=cir(ra,cp,s),
     
path=pts([[-ra+.1,0],[ra-.1,0],[rb-ra,h],[-rb+.1,0]]),
    prism=center==true?trns([0,0,-h/2],prism(sec,path)):prism(sec,path))
    prism;
            
function flip(sec)=[for(i=[len(sec)-1:-1:0])sec[i]];
     
function l_extrude(sec,h=1,a=0,steps=1)=[for(i=[0:a==0?1:(a-0)/steps:a==0?1:a])
    trns([0,0,a==0?h*i:h/a*i],q_rot([str("z",a==0?0:i)],sec))];
 
function sqr(s,center=false)=
let(
m=is_num(s)?s:s.x,
n=is_num(s)?s:s.y,
sec=[[0,0],[m,0],[m,n],[0,n]],
sec1=center==true?[for(p=sec)p-[m/2,n/2]]:sec)
    sec1;

function cub(p,center=false)=
let(
m=is_num(p)?p:p.x,
n=is_num(p)?p:p.y,
o=is_num(p)?p:p.z,

path=pts([[-m/2,0],[m/2,0],[0,o],[-m/2,0]]),
prism=center==true?trns([-m/2,-n/2,-o/2],rsz3d(prism(sqr(m),path),[m,n,o])):rsz3d(prism(sqr(m),path),[m,n,o])
)
prism;

function spr(r,cp=[0,0,0],s=50)=let(
path=arc(r,-90,90,s=s),
prism=[for(p=path)trns([0,0,p.y]+cp,cir(p.x,s=s))])
    prism;

function add_p(p,p1=[0,0],n,i=0)= n==0?p1:add_p(p,[p[i].x+p1.x,p[i].y+p1.y],n-1,i+1);
function pts(p)=[for(n=[1:len(p)])add_p(p=p,p1=[0,0],n=n,i=0)];
    
function add_p1(p,p1=[0,0,0],n,i=0)= n==0?p1:add_p1(p,[p[i].x+p1.x,p[i].y+p1.y,p[i].z],n-1,i+1);
function pts1(p)=[for(n=[1:len(p)])add_p1(p=p,p1=[0,0,0],n=n,i=0)];
    
function add_p2(p,p1=[0,0,0,0],n,i=0)= n==0?p1:add_p2(p,[p[i][0]+p1[0],p[i][1]+p1[1],p[i][2]+p1[2],p[i][3]],n-1,i+1);
function pts2(p)=[for(n=[1:len(p)])add_p2(p=p,p1=[0,0,0,0],n=n,i=0)];

module points(p,d=.5){
    for(i=p)translate(i)cube(size=d,center=true);
    
    }
 
 function offst(sec,r)=let(
//    rev_r=r<0?(min_r(sec)>abs(r)?r:-(min_r(sec)-.1)):r,
rev_r=r,
sec1=[for(i=[0:len(sec)-1])
    let(
p0=i==0?sec[len(sec)-1]:sec[i-1],
p1=sec[i],
p2=i<len(sec)-1?sec[i+1]:sec[0],
v1=p0-p1,u1=v1/norm(v1),
v2=p2-p1,u2=v2/norm(v2),
a1=ang(u1.x,u1.y),
a2=ang(u2.x,u2.y),
a=cw(sec)==-1?(a1>a2?a1:a1+360)-a2:(a2>a1?a2:a2+360)-a1,
p3=p1+u2*rev_r/cos((180-a)/2)*cw(sec)*rm(a/2)
)p3]
    

)sec1;   

function sec_r(sec)=let(

sec1=[for(i=[0:len(sec)-1])
let(
p0=i==0?sec[len(sec)-1]:sec[i-1],
p1=sec[i],
p2=i<len(sec)-1?sec[i+1]: sec[0],
v1=p0-p1,u1=v1/norm(v1),
v2=p2-p1,u2=v2/norm(v2),
p3=p0-u1*rm(90),
p4=p2-u2*rm(-90),
p5=i_p2d([p0,p3],[p2,p4]),
length=norm(p5-p0)
)length
],
minimum_r=min(sec1)


)minimum_r;

function det3d(m)=let(
m11=m[0][0],m12=m[0][1],m13=m[0][2],
m21=m[1][0],m22=m[1][1],m23=m[1][2],
m31=m[2][0],m32=m[2][1],m33=m[2][2],

s11=m22*m33-m32*m23,s12=-(m21*m33-m31*m23),s13=m21*m32-m31*m22,
s21=-(m12*m33-m32*m13),s22=m11*m33-m31*m13,s23=-(m11*m32-m31*m12),
s31=m12*m23-m22*m13,s32=-(m11*m23-m21*m13),s33=m11*m22-m21*m12,

d=m11*s11+m12*s12+m13*s13
)d;



function det2d(m)=let(
m11=m[0][0],m12=m[0][1],
m21=m[1][0],m22=m[1][1],

s11=m22,s12=-m21,
s21=-m12,s22=m11,

d=m11*m22-m21*m12
)d;



function i_m3d(m)=let(
m11=m[0][0],m12=m[0][1],m13=m[0][2],
m21=m[1][0],m22=m[1][1],m23=m[1][2],
m31=m[2][0],m32=m[2][1],m33=m[2][2],

s11=m22*m33-m32*m23,s12=-(m21*m33-m31*m23),s13=m21*m32-m31*m22,
s21=-(m12*m33-m32*m13),s22=m11*m33-m31*m13,s23=-(m11*m32-m31*m12),
s31=m12*m23-m22*m13,s32=-(m11*m23-m21*m13),s33=m11*m22-m21*m12,

d=m11*s11+m12*s12+m13*s13
) 1/d*[[s11,s21,s31],[s12,s22,s32],[s13,s23,s33]];


function i_m2d(m)=let(
m11=m[0][0],m12=m[0][1],
m21=m[1][0],m22=m[1][1],

s11=m22,s12=-m21,
s21=-m12,s22=m11,

d=m11*m22-m21*m12
)1/d*[[s11,s21],[s12,s22]];

function add_v(v,s=[0,0],n=0)=n==len(v)?s:add_v(v,s+v[n],n+1);
function fact(n,m=1)=n==0?m:fact(n-1,m*n);
function comb(n,i)=fact(n)/(fact(i)*fact(n-i));
function bez(p,s=.1)=[for(t=[0:s:1])
    let(n=len(p)-1)add_v([for(i=[0:n])comb(n,i)*(1-t)^(n-i)*t^i*p[i]])];

function 2cir_tarc(r1,r2,cp1,cp2,r)=
let(
l1=norm(cp2-cp1),
l2=r-r1,
l3=r-r2,
//l2^2-x^2=l3^2-(l1-x)^2
//l2^2-x^2=l3^2-l1^2+2*l1*x-x^2
x=(l2^2-l3^2+l1^2)/(2*l1),
h=sqrt(l2^2-x^2),
v1=cp2-cp1,u1=v1/norm(v1),
p0=cp1+u1*x,
cp3=p0-u1*h*rm(90),
v2=cp2-cp3,u2=v2/norm(v2),
v3=cp1-cp3,u3=v3/norm(v3),
ang1=ang(u2.x,u2.y),
ang2=ang(u3.x,u3.y)

)arc(r,ang1,ang2,cp3);
    
function 2p_arc(l,r,cw=1,s=20)=let(
p1=l[0],p2=l[1],
p3=p1+(p2-p1)/2,
d=norm(p3-p1),
l=sqrt(r^2-d^2),
v=p1-p3,u=v/norm(v),
cp=p3+u*l*rm(cw==-1?-90:90),
v1=p1-cp,v2=p2-cp,
a1=ang(v1.x,v1.y),a2=ang(v2.x,v2.y),
a3=cw==-1?(a2<a1?a2+360:a2):(a2<a1?a2:a2-360)

)arc(r,a1,a3,cp,s);

function 2p_arc_cp(l,r,cw=1,s=20)=let(
p1=l[0],p2=l[1],
p3=p1+(p2-p1)/2,
d=norm(p3-p1),
l=sqrt(r^2-d^2),
v=p1-p3,u=v/norm(v),
cp=p3+u*l*rm(cw==-1?-90:90),
v1=p1-cp,v2=p2-cp,
a1=ang(v1.x,v1.y),a2=ang(v2.x,v2.y),
a3=cw==-1?(a2<a1?a2+360:a2):(a2<a1?a2:a2-360)

)cp;


function 2r(l,r,cw=1,s=20)=let(
p1=l[0],p2=l[1],
p3=p1+(p2-p1)/2,
d=norm(p3-p1),
l=sqrt(r^2-d^2),
v=p1-p3,u=v/norm(v),
cp=p3+u*l*rm(cw==-1?90:-90),
v1=p1-cp,v2=p2-cp,
a1=ang(v1.x,v1.y),a2=ang(v2.x,v2.y),
a3=cw==-1?(a2<a1?a2+360:a2):(a2<a1?a2:a2-360)

)arc(r,a1,a3,cp,s);

function 3p_arc(l,s=30)=

let(
p1=l[0],p2=l[1],p3=l[2],
p4=p1+(p2-p1)/2,
p5=p2+(p3-p2)/2,
v1=p2-p4,u1=v1/norm(v1),
v2=p3-p5,u2=v2/norm(v2),
p6=p4+u1*rm(90),
p7=p5+u2*rm(90),
cp=i_p2d([p4,p6],[p5,p7]),
r=norm(p1-cp),
v3=p1-cp,v4=p2-cp,v5=p3-cp,
a1=ang(v3.x,v3.y),
a2=ang(v4.x,v4.y),
a3=ang(v5.x,v5.y),
a4=cw([p1,p2,p3])==-1?(a3<a1?a3+360:a3):(a3<a1?a3:a3-360)

)arc(r,a1,a4,cp,s);

function ellipse(r1,r2,cp,s=30)=
let(
sec=[for(i=[0:360/s:360-360/s])cp+[r1*cos(i),r2*sin(i)]]
)sec;

function l_cir_fillet(line,r1,r2,cp)=let(
p1=line[0],p2=line[1],
v1=p2-p1,u1=v1/norm(v1),
p12=p1+u1*cp.y,
cp1=p12+u1*cp.x*rm(-90),
theta=acos((norm(cp1-p12)-r2)/(r1+r2)),
p121=p12-u1*(r1+r2)*sin(theta),
cp2=p121+u1*r2*rm(-90),
p122=p12+u1*(r1+r2)*sin(theta),
cp3=p122+u1*r2*rm(-90),
v2=p121-cp2,v3=cp1-cp2,v4=cp2-cp1,v5=cp3-cp1,v6=cp1-cp3,v7=p122-cp3,
a1=ang(v2.x,v2.y),
a2=ang(v3.x,v3.y)>a1?ang(v3.x,v3.y)-360:ang(v3.x,v3.y),
a3=ang(v4.x,v4.y),
a4=ang(v5.x,v5.y)<a3?ang(v5.x,v5.y)+360:ang(v5.x,v5.y),
a5=ang(v6.x,v6.y),
a6=ang(v7.x,v7.y)>a5?ang(v7.x,v7.y)-360:ang(v7.x,v7.y)
)
[p1,each arc(r2,a1,a2,cp2),each arc(r1,a3,a4,cp1),each arc(r2,a5,a6,cp3),p2];

function avg_v(vector)=let(
x=len(vector[0])==3?sum(vector*[1,0,0])/len(vector):sum(vector*[1,0])/len(vector),
y=len(vector[0])==3?sum(vector*[0,1,0])/len(vector):sum(vector*[0,1])/len(vector),
z=len(vector[0])==3?sum(vector*[0,0,1])/len(vector):[],
)len(vector[0])==3?[x,y,z]:[x,y];


function rsz3d(prism,rsz=[1,1,1])=
let(
rev_p_list=[for(p=prism) each[for(p1=p)p1]],
max_x=max(rev_p_list*[1,0,0]),
max_y=max(rev_p_list*[0,1,0]),
max_z=max(rev_p_list*[0,0,1]),
min_x=min(rev_p_list*[1,0,0]),
min_y=min(rev_p_list*[0,1,0]),
min_z=min(rev_p_list*[0,0,1]),
avg=avg_v(rev_p_list),

r_x=rsz.x/(max_x-min_x),
r_y=rsz.y/(max_y-min_y),
r_z=rsz.z/(max_z-min_z)
)[for(i=[0:len(prism)-1])
    [for(p=prism[i])avg+[r_x*(p.x-avg.x),r_y*(p.y-avg.y),r_z*(p.z-avg.z)+(avg.z-min_z)*r_z-(avg.z-min_z)]]];
        
 function rsz(sec,rsz=[1,1,1])=
let(
avg=avg_v(sec),
max_x=len(rsz)==2?max(sec*[1,0]):max(sec*[1,0,0]),
min_x=len(rsz)==2?min(sec*[1,0]):min(sec*[1,0,0]),
max_y=len(rsz)==2?max(sec*[0,1]):max(sec*[0,1,0]),
min_y=len(rsz)==2?min(sec*[0,1]):min(sec*[0,1,0]),

r_x=rsz.x/(max_x-min_x),
r_y=rsz.y/(max_y-min_y)


)[for(i=[0:len(sec)-1])let(
p=avg+[r_x*(sec[i].x-avg.x),r_y*(sec[i].y-avg.y)-((min_y-avg.y)*r_y-(min_y-avg.y))]

)p];

function l(l,s=20)=
let(
p0=l[0],p1=l[1],
v=p1-p0,u=v/norm(v),
length=norm(v)
)[for(i=[0:length/s:length])p0+u*i];
    
function l1(l,s=20)=
let(
p0=l[0],p1=l[1],
v=p1-p0,u=v/norm(v),
length=norm(v)
)[for(i=[0:length/s:length])p0+u*i];
    
function 3p_r(p1,p2,p3)=

let(
p4=p1+(p2-p1)/2,
p5=p2+(p3-p2)/2,
v1=p2-p4,u1=v1/norm(v1),
v2=p3-p5,u2=v2/norm(v2),
p6=p4+u1*rm(90),
p7=p5+u2*rm(90),
cp=i_p2d([p4,p6],[p5,p7]),
r=norm(p1-cp)

)r;

function min_r(sec)=
min([for(i=[0:len(sec)-1])3p_r(sec[i==0?len(sec)-1:i-1],sec[i],sec[i<len(sec)-1?i+1:0])]);
    
function m_points(sec,sl=20)=
[for(i=[0:len(sec)-1])let(
p0=sec[i],
p1=sec[i<len(sec)-1?i+1:0],
lnth=norm(p1-p0),
sec1=lnth>sl?l1([p0,p1],lnth/sl):[p0],
sec2=[for(i=[0:len(sec1)-1])if(sec1[i]!=sec1[i<len(sec1)?i+1:0])sec1[i]])
each sec2];

function cum_sum(list,list1,n,s=1)=n==0?list1:cum_sum(list,[for(i=[0:s])list[i]]*[for(i=[0:s])1],n-1,s+1);
function add_paths(path1,path2)=
let(
sec1=[for(i=[0:len(path2)-2])norm(path2[i+1]-path2[i])],
sec2=[0,each [for(i=[0:len(sec1)-1])cum_sum(sec1,sec1[0],n=i)]],
sec3=[for(i=[0:len(sec2)-1])[sec2[i],0]],
length=sec2[len(sec2)-1],
height=max(path1*[0,1])-min(path1*[0,1]),
path3=rsz(path1,[length,height]),
d=path3[0]-path1[0],
path4=[for(p=path3)p-d],

path5=[for(i=[0:len(sec3)-1])each [for(j=[0:len(path4)-2])let(
l1=[sec3[i],sec3[i]+[0,1]],
l2=[path4[j],path4[j+1]],
ip=i_p2d(l1,l2),
v1=l2[1]-l2[0],u1=v1/norm(v1),
v2=ip-l2[0],u2=v2/norm(v2),
lnth1=norm(v1),lnth2=norm(v2)
)if(i==0&&j==0)l2[0] else if(l2[0].x<=ip.x &&l2[1].x>=ip.x)ip]],

path6=[for(i=[0:len(path2)-1])[path2[i].x,path2[i].y,path5[i].y]]    

)path6;


//module swp1(prism) let(
//n=len(prism[0]),
//points=[for(p=prism)each [for(p1=p)p1]],
//faces1=[for(j=[0:n:len(points)-2*n])[for(i=[j:j+n-1])i]],
//faces2=[for(j=[len(points)-n])[for(i=[j+n-1:-1:j])i]],
//faces3=[for(j=[0:n:len(points)-2*n])[for(i=[j:j+n-2])[i,i+n,i+n+1,i+1]]],
//faces4=[for(i=[0:n:len(points)-n-1])[i,i+n-1,i+2*n-1,i+n]]
//)polyhedron(points,[each faces1,each faces2,each each faces3,each faces4],convexity=10);

module swp(surf1)

let(l=len(surf1[0]),
p0=[for(j=[0:len(surf1)-1])each surf1[j]],
p1=[each [for(j=[0:len(surf1)-1])if(j==0)[for(i=[0:l-1])i+j*l]],
each [for(j=[0:len(surf1)-2])each [for(i=[0:l-1])let(i_plus=i<l-1?i+1:0)[i+l*j,i+l+l*j,i_plus+l+l*j,i_plus+l*j]]],
each [for(j=[0:len(surf1)-1])if(j==len(surf1)-1)[for(i=[l-1:-1:0])i+l*j]]
    ]
)
polyhedron(p0,p1,convexity=10);

module swp_c(surf1)

let(l=len(surf1[0]),
p0=[for(j=[0:len(surf1)-1])each surf1[j]],
p1=[//each [for(j=[0:len(surf1)-1])if(j==0)[for(i=[0:l-1])i+j*l]],
each [for(j=[0:len(surf1)-2])each [for(i=[0:l-1])let(i_plus=i<l-1?i+1:0)[i+l*j,i+l+l*j,i_plus+l+l*j,i_plus+l*j]]]//,
//each [for(j=[0:len(surf1)-1])if(j==len(surf1)-1)[for(i=[l-1:-1:0])i+l*j]]
    ]
)
polyhedron(p0,p1,convexity=10);

function 2cyl_fillet(r1,r2,cp1,cp2,r,path)=[for(p=path)trns([0,0,p.y],2cir_fillet(r1+p.x,r2+p.x,cp1,cp2,r))];
    
function 2cir_fillet1(r1,r2,c1,c2,r)=
let(
l1=norm(c2-c1),l2=r1+r,l3=r2+r,
t=(l1^2+l2^2-l3^2)/(2*l1),
h=sqrt(l2^2-t^2),
v=c2-c1,u=v/norm(v),
p1=c1+u*t+u*[[0,1],[-1,0]]*h,
a1=ang((c1-p1).x,(c1-p1).y),
a2=ang((c2-p1).x,(c2-p1).y),
p2=c1+u*t+u*[[0,-1],[1,0]]*h,
a3=ang((c2-p2).x,(c2-p2).y),
a4=ang((c1-p2).x,(c1-p2).y),
a5=ang((p1-c1).x,(p1-c1).y),
a6=ang((p2-c1).x,(p2-c1).y),
a7=ang((p1-c2).x,(p1-c2).y),
a8=ang((p2-c2).x,(p2-c2).y),
a9=((a4<a3?360+a4:a4)+a3)/2,
p10=p2+[r/cos(a9-a3)*cos(a9),r/cos(a9-a3)*sin(a9)],

arc1=arc(r,a2<a1?360+a2:a2,a1,p1),
arc2=arc(r,a4<a3?360+a4:a4,a3,p2),
arc3=arc(r2,a7<a8?a7+360:a7,a8,c2),
arc4=arc(r1,a5,a6<a5?a6+360:a6,c1)
)
[p10,each arc2];

function 2spr_fillet(r1,r2,cp1,cp2,r)=
let(
v=cp2-cp1,u=v/norm(v),
l=norm(cp2-cp1),
c1=[0,0],
c2=[l,0],
a1=u==[0,0,1]?90:u==[0,0,-1]?-90:ang(sqrt(v.x^2+v.y^2),v.z),
a2=u==[0,0,1]||u==[0,0,-1]?0:ang(v.x,v.y),


sec=2cir_fillet1(r1,r2,c1,c2,r),
prism=trns(cp1,q_rot([str("y",-a1),str("z",a2)],[for(i=[0:5:360])[for(p=sec)q([1,0,0],[p.x,p.y,0],i)]]))
    
)flip(prism);

function l_cir_ip(cir,line)=
let(
p0=line[0],p1=line[1],
v=p1-p0,u=v/norm(v),
ip=[for(i=[0:len(cir)-1])let(
p2=cir[i],p3=i<len(cir)-1?cir[i+1]:cir[0],
int_p=i_p2d([p2,p3],[p0,p1]),
v1=p3-p2,u1=v1/norm(v1),
v2=int_p-p2,u2=v2/norm(v2),
l1=norm(p3-p2),l2=norm(int_p-p2)

)//if(u1==u2&&l1>=l2)int_p
   if(l1>=l2&&norm(u1-u2)<.01)int_p ]

)ip;
   
function offst_l(l,d)=
let(
v=l[1]-l[0],u=v/norm(v),
p0=l[0]+u*d*rm(-90),
p1=l[1]+u*d*rm(-90)
)[p0,p1];
   
function perp(line,point)=
let(
v1=line[1]-line[0],
slope1=v1.y/v1.x,
slope2=-1/slope1,
line1=[point,point+[1,slope2]],
ip=i_p2d(line,line1)

)[point,ip];
   
function 2cir_tangent(r1,r2,cp1,cp2)=
let(
v=cp2-cp1,u=v/norm(v),
theta=ang(v.x,v.y),
theta1=atan((r1-r2)/norm(v)),
p1=cp1+u*r1*rm(90-theta1),
p0=cp2+u*r2*rm(90-theta1),
p2=cp1+u*r1*rm(-(90-theta1)),
p3=cp2+u*r2*rm(-(90-theta1))

)[p0,p1,p2,p3];

function cvar(a)=
let(
text=a[0],
b=search("-",a)!=[]&&search(".",a)!=[]?[for(i=[2:len(a)-1])if(a[i]!=".")a[i]]:
search("-",a)!=[]&&search(".",a)==[]?[for(i=[2:len(a)-1])a[i]]:
search("-",a)==[]&&search(".",a)!=[]?[for(i=[1:len(a)-1])if(a[i]!=".")a[i]]:[for(i=[1:len(a)-1])a[i]],
n=[["0"],["1"],["2"],["3"],["4"],["5"],["6"],["7"],["8"],["9"]],
l=len(b),
c=[for(i=[0:l-1])each search(b[i],n)],
d=[for(i=[0:len(c)-1])10^(l-i-1)],
e=c*d,
f=search(".",a)!=[]&&search("-",a)!=[]?e*10^-(l-(search(".",a)[0]-2)):search(".",a)!=[]&&search("-",a)==[]?e*10^-(l-(search(".",a)[0]-1)):e,
g=search("-",a)!=[]?f*-1:f
)[text,g];

function o_set1(sec,d=1)=let(

sec1=[for(i=[0:len(sec)-1])let(
p0=i==0?sec[len(sec)-1]:sec[i-1],
p1=sec[i],
p2=(i<len(sec)-1?sec[i+1]:sec[0]),
v1=p0-p1,u1=v1/norm(v1),
v2=p2-p1,u2=v2/norm(v2),
theta=acos(u1*u2),theta1=(180-theta)/2,
start_angle=ang(u2.x,u2.y)
) d<0?arc(abs(d)/cos(theta1),start_angle,start_angle+theta,sec[i],2)[1]:arc(abs(d)/cos(theta1),start_angle-360,start_angle+theta,sec[i],2)[1]],
sec2=[for(i=[0:len(sec1)-1])let(
p0=sec1[i],p1=i<len(sec1)-1?sec1[i+1]:sec1[0],
v1=p1-p0,u1=v1/norm(v1),
ip=[for(j=[0:len(sec1)-1])let(p2=sec1[j],p3=j<len(sec1)-1?sec1[j+1]:sec1[0])if(j!=i)i_p2d([p0,p1],[p2,p3])],
l1=norm(p1-p0),
ipf=[for(p=ip)let(u2=(p-p0)/norm(p-p0))if(norm(p-p0)<l1 && sign(u1.x)==sign(u2.x) && sign(u1.y)==sign(u2.y))p]

)if (len(ipf)>0)each ipf else sec1[i]],
 sec3=[for(p=sec2)if(min([for(p1=m_points(sec,1))norm(p-p1)])>abs(d))p]

)sort_points(sec,remove_extra_points(sec3));

 
function o_set(sec,d=1)=abs(d)<=sec_r(sec)?offst(sec,d):o_set1(sec,d);

//function offst(sec,d)=d<=0?(abs(d)>sec_r(sec)?o_set(sec,d):offst1(sec,d)):offst1(sec,d);

function rot(axis,ang,sec)=[for(p=sec)let(point=len(p)==2?[p.x,p.y,0]:p)q(axis,point,ang)];
 
function s_pnt(sec)=let(
y_min=min(sec*[0,1]),
loc=search(y_min,sec,0,1)[0]
)sec[loc];

function reduced_list(sec,list)=[for(p=sec)each [for(p1=list)if(norm(p-p1)>.01|| p1==[])p]];
 
function list_of_points_to_omit(sec,point)=let(
list=[for(i=[0:len(sec)-1])if(norm(sec[i]-point)<.001)i],
list1=len(list)>1?[for(i=[1:len(list)-1])list[i]]:[]

)list1;

function revised_list(sec,index_list)=let(
a=[for(i=[0:len(sec)-1]) if(search(0,[for(j=index_list)i-j],0)==[])i],
sec1=[for(i=a)sec[i]]
)sec1;

function remove_extra_points(sec,n=0)=
n==len(sec)?sec:remove_extra_points(
let(
a=list_of_points_to_omit(sec,sec[n]),
b=revised_list(sec,a)
)b,n+1
);

function sort_points(sec,list)=[if(list!=[])let(
a=[for(p=sec)min([for(i=[0:len(list)-1])norm(list[i]-p)])],
b=[for(p=sec)[for(i=[0:len(list)-1])norm(list[i]-p)]],
c=[for(i=[0:len(sec)-1])each search(a[i],b[i])],
d=[for(i=c)list[i]]
)d][0];


module swp_h(sec,path,t=-.5){
prism1=p_extrude(sec,path);
prism2=p_extrude(f_offset(sec,t),path);
prism3=[each each prism1,each each prism2];
let(
n=len(sec),
p=len(path),
faces1=[for(i=[0:n-1])i<n-1?[i,i+1,i+1+n*p,i+n*p]:[i,i+1-n,n*p,i+n*p]],
faces2=[for(i=[0:n*p-n-1])(i+1)%n==0?[i,i+n,i+1,i+1-n]:[i,i+n,i+1+n,i+1]],
faces3=[for(i=[n*p:2*n*p-1-n])(i+1)%n==0?[i,i+1-n,i+1,i+n]:[i,i+1,i+1+n,i+n]],
faces4=[for(i=[n*p-n:n*p-1])i<n*p-1?[i,i+n*p,i+1+n*p,i+1]:[i,i+n*p,i+1-n+n*p,i+1-n]]
)polyhedron(prism3,[each faces1,each faces2, each faces3, each faces4]);
}

module swp_prism_h(prism,prism1){

prism2=[each each prism,each each prism1];
let(
n=len(prism[0]),
p=len(prism),
faces1=[for(i=[0:n-1])i<n-1?[i,i+1,i+1+n*p,i+n*p]:[i,i+1-n,n*p,i+n*p]],
faces2=[for(i=[0:n*p-n-1])(i+1)%n==0?[i,i+n,i+1,i+1-n]:[i,i+n,i+1+n,i+1]],
faces3=[for(i=[n*p:2*n*p-1-n])(i+1)%n==0?[i,i+1-n,i+1,i+n]:[i,i+1,i+1+n,i+n]],
faces4=[for(i=[n*p-n:n*p-1])i<n*p-1?[i,i+n*p,i+1+n*p,i+1]:[i,i+n*p,i+1-n+n*p,i+1-n]]
)polyhedron(prism2,[each faces1,each faces2, each faces3, each faces4]);
}

function outer_offset(sec1,d)=d==0?(cw(sec)==1?flip(sec1):sec1):
let(
sec=cw(sec1)==1?flip(sec1):sec1,
r=abs(d),
op=[for(i=[0:len(sec)-1])let(
p0=i==0?sec[len(sec)-1]:sec[i-1],
p1=sec[i],
p2=i<len(sec)-1?sec[i+1]:sec[0],
v1=p0-p1, u1=v1/norm(v1),
v2=p2-p1, u2=v2/norm(v2),
theta=acos (u1*u2),
alpha=180-theta,
pa=p1+u1*r*tan(alpha/2),
pb=p1+u2*r*tan(alpha/2),
cp=2p_arc_cp([pa,pb],r,1),
pc=p1+u1*r*rm(90),
pd=p1+u2*r*rm(-90)
) cw([p0,p1,p2])==-1?2p_arc([pc, pd],r,-1,s=norm(pc-pd)<1?0:5):[cp]],

op01=[for(i=[0:len(sec)-1])let(
p0=i==0?sec[len(sec)-1]:sec[i-1],
p1=sec[i],
p2=i<len(sec)-1?sec[i+1]:sec[0],
radius=3p_r(p0,p1,p2)
) if((radius>=r)||(cw( [p0,p1,p2])==-1))each op[i]],
op02=[for(i=[0:len(op01)-1])let(
p0=op01[i],p1=i<len (op01)-1?op01[i+1]:op01[0],
v1=p1-p0, u1=v1/norm(v1),
ip=[for(j=i==0?[len(op01)-2,i+2]:i==1?[len(op01)-1,i+2]:i==len(op01)-1?[i-2,1]:i==len(op01)-2?[i-2,0]:[i-2,i+2])let(p2=op01[j],p3=j==len(op01)-1?op01[0]:op01[j+1])i_p2d([p0,p1],[p2,p3])],
l1=norm(p1-p0),
ipf=[for(p=ip)let(u2=(p-p0)/norm(p-p0)) if(norm(p-p0)<l1 && sign(u1.x)==sign(u2.x) && sign(u1.y)==sign(u2.y))p]
)if (len(ipf)>0) each ipf else op01[i]],
op03=[for(p=op02) if(min([for(p1=m_points(sec,r))norm(p-p1)])>=abs(d)-.1)p]
)sort_points(sec, remove_extra_points(op03));


function inner_offset(sec1,d)=d==0?(cw(sec1)==1?flip(sec1):sec1):
let(
sec=cw(sec1)==1?flip(sec1):sec1,
r=abs(d),
op=[for(i=[0:len(sec)-1])let(
p0=i==0?sec[len(sec)-1]:sec[i-1],
p1=sec[i],
p2=i<len(sec)-1?sec[i+1]:sec[0],
v1=p0-p1, u1=v1/norm(v1),
v2=p2-p1, u2=v2/norm(v2),
theta=acos (u1*u2),
alpha=180-theta,
pa=p1+u1*r*tan(alpha/2),
pb=p1+u2*r*tan(alpha/2),
cp=2p_arc_cp([pa,pb],r,-1),
pc=p1+u1*r*rm(-90),
pd=p1+u2*r*rm(90)
) cw( [p0, p1, p2])==-1?[cp]:2p_arc([pc, pd],r,1,s=norm(pc-pd)<1?0:5)],
op01=[for(i=[0:len(sec)-1])let(
p0=i==0?sec[len(sec)-1]:sec[i-1],
p1=sec[i],
p2=i<len(sec)-1?sec[i+1]:sec[0],
radius=3p_r(p0, p1, p2)
) if((radius>=r)||(cw( [p0, p1, p2])==1)) each op[i]],
op02=[for(i=[0:len(op01)-1])let(
p0=op01[i],p1=i<len(op01) -1?op01[i+1]: op01[0],
v1=p1-p0, u1=v1/norm(v1),
ip=[for(j=i==0?[len(op01)-2,i+2]:i==1?[len(op01)-1,i+2]:i==len(op01)-1?[i-2,1]:i==len(op01)-2?[i-2,0]:[i-2,i+2])let(p2=op01[j],p3=j==len(op01)-1?op01[0]:op01[j+1])i_p2d([p0,p1],[p2,p3])],
l1=norm(p1-p0),
ipf=[for(p=ip) let(u2=(p-p0)/norm(p-p0))if(norm(p-p0)<l1 && sign(u1.x) ==sign(u2.x) && sign(u1.y)==sign(u2.y))p])
if (len(ipf)>0)each ipf else op01[i]],
op03=[for(p=op02) if(min([for(p1=m_points (sec,r))norm(p-p1)])>=abs(d)-.001)p]
) sort_points (sec, remove_extra_points (op03));

function f_offset(sec,d)=d<=0?inner_offset(sec,d):outer_offset(sec,d);


function v_sec_extrude(sec,path,o)=[for(i=[0:len(path)-2])let(
off=o/(len(path)-1),
sec=f_offset(sec,off*i),
p0=path[i],
p1=path[i+1],
v=p1-p0,
v1=[v.x,v.y,0],
u=[v.x,v.y]/norm([v.x,v.y]),
u1=v/norm(v),
u2=v1/norm(v1),
theta=!is_num(u.x)?0:(u.y<0?360-acos([1,0]*u):acos([1,0]*u)),
a=u1==u2?0:u1.z<0?360-acos(u1*u2):acos(u1*u2),
alpha=a-90,
rev_sec=q_rot(["x90","z90",str("y",-a),str("z",theta)],sec)
)each i<len(path)-2?[trns(p0,rev_sec)]:[trns(p0,rev_sec),trns(p1,rev_sec)]];

function p_extrudec1(sec,path)=[for(i=[0:len(path)-1])let(
p0=path[i],
p1=i<len(path)-1?path[i+1]:path[0],
v=p1-p0,
u=[v.x,v.y,0]/norm([v.x,v.y,0]),
u1=v/norm(v),
theta=!is_num(u)?0:(u.y<0?360-acos([1,0,0]*u):acos([1,0,0]*u)),
alpha=u1.z<0?360-acos([1,0,0]*u1):acos([1,0,0]*u1)

)trns(p0,q_rot(["x90","z90",str("y",-alpha),str("z",theta)],sec))];

module p_extrude1(sec,path) swp([for(i=[0:len(path)-2])let(
p0=path[i],
p1=path[i+1],
v=p1-p0,
v1=[v.x,v.y,0],
u=[v.x,v.y]/norm([v.x,v.y]),
u1=v/norm(v),
u2=v1/norm(v1),
theta=!is_num(u.x)?0:(u.y<0?360-acos([1,0]*u):acos([1,0]*u)),
a=u1.z<0?360-acos(u1*u2):acos(u1*u2),
alpha=a-90,
rev_sec=q_rot(["x90","z90",str("y",-a),str("z",theta)],sec)
)each i<len(path)-2?[trns(p0,rev_sec)]:[trns(p0,rev_sec),trns(p1,rev_sec)]]);

module p_extrudec1(sec,path) swp([for(i=[0:len(path)-1])let(
p0=path[i],
p1=i<len(path)-1?path[i+1]:path[0],
v=p1-p0,
u=[v.x,v.y,0]/norm([v.x,v.y,0]),
u1=v/norm(v),
theta=!is_num(u)?0:(u.y<0?360-acos([1,0,0]*u):acos([1,0,0]*u)),
alpha=u1.z<0?360-acos([1,0,0]*u1):acos([1,0,0]*u1)

)trns(p0,q_rot(["x90","z90",str("y",-alpha),str("z",theta)],sec))]);

module p_extrudec(sec,path) swp_c([for(i=[0:len(path)-1])let(
p0=path[i],
p1=i<len(path)-1?path[i+1]:path[0]-(path[1]-path[0])*.01,
v=p1-p0,
u=v/norm(v),
theta=u.y<0?360-acos([1,0,0]*u):acos([1,0,0]*u),
prism=q_rot(["x90","z90"],sec),

p2=path[0],
p3=path[1],
v1=p3-p2,
u1=v1/norm(v1),
theta1=u1.y<0?360-acos([1,0,0]*u1):acos([1,0,0]*u1),
prism1=q_rot(["x90","z90"],sec)

)each i<len(path)-1?[trns(p0,q_rot([str("z",theta)],prism))]:[trns(p0,q_rot([str("z",theta)],prism)),trns(p1,q_rot([str("z",theta)],prism)),trns(p2,q_rot([str("z",theta1)],prism1))]]);

module p_extrude(sec,path) swp([for(i=[0:len(path)-2])let(
p0=path[i],
p1=path[i+1],
v=p1-p0,
v1=[v.x,v.y,0],
u=[v.x,v.y]/norm([v.x,v.y]),
u1=v/norm(v),
u2=v1/norm(v1),
theta=!is_num(u.x)?0:(u.y<0?360-acos([1,0]*u):acos([1,0]*u)),
a=u1==u2?0:u1.z<0?360-acos(u1*u2):acos(u1*u2),
alpha=a-90,
rev_sec=q_rot(["x90","z90",str("y",-a),str("z",theta)],sec)
)each i<len(path)-2?[trns(p0,rev_sec)]:[trns(p0,rev_sec),trns(p1,rev_sec)]]);

function p_extrudec(sec,path)= [for(i=[0:len(path)-1])let(
p0=path[i],
p1=i<len(path)-1?path[i+1]:path[0]-(path[1]-path[0])*.01,
v=p1-p0,
u=v/norm(v),
theta=u.y<0?360-acos([1,0,0]*u):acos([1,0,0]*u),
prism=q_rot(["x90","z90"],sec),

p2=path[0],
p3=path[1],
v1=p3-p2,
u1=v1/norm(v1),
theta1=u1.y<0?360-acos([1,0,0]*u1):acos([1,0,0]*u1),
prism1=q_rot(["x90","z90"],sec)

)each i<len(path)-1?[trns(p0,q_rot([str("z",theta)],prism))]:[trns(p0,q_rot([str("z",theta)],prism)),trns(p1,q_rot([str("z",theta)],prism)),trns(p2,q_rot([str("z",theta1)],prism1))]];

function p_extrude(sec,path)= [for(i=[0:len(path)-2])let(
p0=path[i],
p1=path[i+1],
v=p1-p0,
v1=[v.x,v.y,0],
u=[v.x,v.y]/norm([v.x,v.y]),
u1=v/norm(v),
u2=v1/norm(v1),
theta=!is_num(u.x)?0:(u.y<0?360-acos([1,0]*u):acos([1,0]*u)),
a=u1==u2?0:u1.z<0?360-acos(u1*u2):acos(u1*u2),
alpha=a-90,
rev_sec=q_rot(["x90","z90",str("y",-a),str("z",theta)],sec)
)each i<len(path)-2?[trns(p0,rev_sec)]:[trns(p0,rev_sec),trns(p1,rev_sec)]];

function p_ex(sec,path)= [for(i=[0:len(path)-2])let(
p0=path[i],
p1=path[i+1],
v=p1-p0,
v1=[v.x,v.y,0],
u=[v.x,v.y]/norm([v.x,v.y]),
u1=v/norm(v),
u2=v1/norm(v1),
theta=!is_num(u.x)?0:(u.y<0?360-acos([1,0]*u):acos([1,0]*u)),
a=u1==u2?0:u1.z<0?360-acos(u1*u2):acos(u1*u2),
alpha=a-90,
rev_sec=q_rot(["x90","z90",str("y",-a),str("z",theta)],sec)
) trns(p0,rev_sec)];

function 3p_3d_fillet(p0,p1,p2,r=1, s=5)=
let(
v1=p0-p1, u1=v1/norm(v1),
v2=p2-p1, u2=v2/norm(v2),
n=cross (u1, u2),
theta=acos (u1*u2),
alpha= (180-theta)/2,
pa=r*tan (alpha) *u1,
pb=r*tan (alpha) *u2,
pap=pa+q(n, u1,90),
pbp=pb+q(n, u2,-90),
l1=[pa, pap],
l2=[pb, pbp],
cp=i_p3d (l1,l2),
arc=trns(p1+cp,[for(i=[0:alpha*2/s:alpha*2])q(n,pb-cp,i)])

) [p1,each arc];

function 3p_3d_fillet_wo_pivot(p0,p1,p2,r=1, s=5)=[
let(
v1=p0-p1, u1=v1/norm(v1),
v2=p2-p1, u2=v2/norm(v2),
n=cross (u1, u2)==[0,0,0]?nv3d(u2):cross (u1, u2),
theta=acos (u1*u2),
alpha= (180-theta)/2,
pa=r*tan (alpha) *u1,
pb=r*tan (alpha) *u2,
pap=pa+q(n, u1,90),
pbp=pb+q(n, u2,-90),
l1=[pa, pap],
l2=[pb, pbp],
cp=i_p3d (l1,l2),
arc=trns(p1+cp,[for(i=[0:alpha*2/s:alpha*2])q(n,pb-cp,i)])
)each r==0?[p1]:arc];

function 3p_3d_arc(points=[p0, p1,p2], s=5)=
let(
v1=points[0]-points[1], u1=v1/norm(v1),
v2=points[2]-points[1], u2=v2/norm(v2),
n=cross(u1, u2),
alpha=acos(u1*u2),
pa=v1/2,
pb=v2/2,
pap=pa+q(n,u1,90),
pbp=pb+q(n,u2,-90),
l1=[pa, pap],
l2=[pb, pbp],
cp=i_p3d (l1,l2),
v3=points[0]-(points[1]+cp),u3=v3/norm(v3),
v4=points[2]-(points[1]+cp),u4=v4/norm(v4),
theta=alpha<90?360-acos(u3*u4):acos(u3*u4),
radius=norm(pa-cp),
arc=trns(points[1]+cp,[for(i=[0:theta/s:theta])q(n,points[0]-(points[1]+cp),-i)])
)arc;

function 3p_3d_r(points=[p0, p1,p2])=
let(
v1=points[0]-points[1], u1=v1/norm(v1),
v2=points[2]-points[1], u2=v2/norm(v2),
n=cross(u1, u2),
alpha=acos(u1*u2),
pa=v1/2,
pb=v2/2,
pap=pa+q(n,u1,90),
pbp=pb+q(n,u2,-90),
l1=[pa, pap],
l2=[pb, pbp],
cp=i_p3d (l1,l2),
v3=points[0]-(points[1]+cp),u3=v3/norm(v3),
v4=points[2]-(points[1]+cp),u4=v4/norm(v4),
theta=alpha<90?360-acos(u3*u4):acos(u3*u4),
radius=norm(pa-cp))
radius;

function 3d_arc(v, r, theta1=0, theta2=180, cw=-1,s=50)=
let(
v=v+[0,0,.0001],
u=v/norm (v),
v1=[r,0,0], u1=v1/norm (v1),
n=cross (v, v1),
theta=90-acos (u*u1),
alpha=u.y<0?360-acos ([1,0] * [u.x, u.y]): acos ([1,0]*[u
.x,u.y]),
v2=q(v,q(n, v1, theta),alpha),
arc=[for (i=[theta1: (theta2-theta1)/s: theta2]) q(v,
v2,-i*cw)])arc;

function c2t3(sec)=trns([0,0,0],sec);

 function lim(t,s=0,e=1)=t>=s&&t<=e;
 
function t(m)=[[m.x.x,m.y.x,m.z.x],[m.x.y,m.y.y,m.z.y],[m.x.z,m.y.z,m.z.z]];

function loop(sec,a,b)=[for(i=[a:b])sec[i]];

function cr3d(l,s=5)=let(
p=[for(i=[0:len(l)-1])[l[i].x,l[i].y,l[i].z]],

r=[for(i=[0:len(l)-1])is_undef(l[i][3])?0:l[i][3]],

theta=[for(i=[0:len(p)-1])let(
i_minus=i==0?len(p)-1:i-1,i_plus=i<len(p)-1?i+1:0,
p0=p[i_minus],p1=p[i],p2=p[i_plus],
v1=p0-p1,v2=p2-p1,
alpha=acos(uv(v1)*uv(v2))
)(180-alpha)/2],

l1=[for(i=[0:len(p)-1])let(
i_minus=i==0?len(p)-1:i-1,
p0=p[i_minus],p1=p[i]
)norm(p0-p1)],

l2=[for(i=[0:len(p)-1])let(
i_minus=i==0?len(p)-1:i-1)r[i_minus]*tan(theta[i_minus])+r[i]*tan(theta[i])],
compare=[for(i=[0:len(p)-1])l1[i]>=l2[i]],

arcs=[for(i=[0:len(p)-1])each assert(compare[i],"radius too big")let(
i_minus=i==0?len(p)-1:i-1,i_plus=i<len(p)-1?i+1:0)3p_3d_fillet_wo_pivot(p[i_plus],p[i],p[i_minus],r[i],s=s)]
)arcs;


function uv(v)=v/norm(v);

function sum(list)=[for(i=[0:len(list)-1])1]*list;

function cumsum(list)=[for(i=[0:len(list)-1])sum([for(j=[0:i])list[j]])];

function add_p3(p,p1=[0,0,0,0,list],n,i=0)= n==0?p1:add_p3(p,[p[i][0]+p1[0],p[i][1]+p1[1],p[i][2]+p1[2],p[i][3],p[i][4]],n-1,i+1);
function pts3(p)=[for(n=[1:len(p)])add_p3(p=p,p1=[0,0,0,0],n=n,i=0)];

function rsz3d_enc(prism,rsz=[1,1,1])=
let(
rev_points_list=[for(p=prism)each[for(p1=p)p1]],
x=max(rev_points_list*[1,0,0])-min(rev_points_list*[1,0,0])+rsz.x,
y=max(rev_points_list*[0,1,0])-min(rev_points_list*[0,1,0])+rsz.y,
z=max(rev_points_list*[0,0,1])-min(rev_points_list*[0,0,1])+rsz.z,
rev_prism=rsz3d(prism,[x,y,z]),
//avg of prism
x1=(max(rev_points_list*[1,0,0])+min(rev_points_list*[1,0,0]))/2,
y1=(max(rev_points_list*[0,1,0])+min(rev_points_list*[0,1,0]))/2,
z1=(max(rev_points_list*[0,0,1])+min(rev_points_list*[0,0,1]))/2,
avg=[x1,y1,z1],
//avg of rev_prism
rev_p_list=[for(p=rev_prism) each[for(p1=p)p1]],
x2=(max(rev_p_list*[1,0,0])+min(rev_p_list*[1,0,0]))/2,
y2=(max(rev_p_list*[0,1,0])+min(rev_p_list*[0,1,0]))/2,
z2=(max(rev_p_list*[0,0,1])+min(rev_p_list*[0,0,1]))/2,
avg_rev=[x2,y2,z2],
)trns(avg-avg_rev,rev_prism)
    ;
    
 function rsz3d_offset(prism,d=1)=
let(
rev_prism=[for(i=[0:len(prism)-2])[for(j=[0:len(prism[i])-1])
let(
j_plus=j<len(prism[i])-1?j+1:0,
p0=prism[i][j],p1=prism[i][j_plus],p2=prism[i+1][j],
v1=p1-p0,v2=p2-p0,
v3=cross(uv(v1),uv(v2))*d,
rev_p=p0+v3
)rev_p
]]
)rev_prism;

function bb(prism)=
let(
p=[for(p=prism)each[for(p1=p) p1]],
bb=[max(p*[1,0,0])-min(p*[1,0,0]),max(p*[0,1,0])-min(p*[0,1,0]),max(p*[0,0,1])-min(p*[0,0,1])])bb;

function avg_prism(prism)=
let(
rev_points_list=[for(p=prism)each each[for(p1=p)p1]],
avg_x=sum(rev_points_list*[1,0,0])/len(rev_points_list),
avg_y=sum(rev_points_list*[0,1,0])/len(rev_points_list),
avg_z=sum(rev_points_list*[0,0,1])/len(rev_points_list)

)[avg_x,avg_y,avg_z];
