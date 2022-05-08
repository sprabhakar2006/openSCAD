
//function to make a prism with combination of 2d section and 2d path
// Example:
// sec=cir(10);
// path=cr(pts1([[2,0],[-2,0,2],[0,10,3],[-3,0]]),5);
// prism=prism(sec,path);
// swp(prism);

function prism(sec,path,m_points=0)=[for(p=path)[for(p1=sort_points(m_points_sc(sec,m_points),m_points_sc(f_offset(sec,round(p.x*100)/100),m_points)))[p1.x,p1.y,p.y]]];

// high quality prism, takes slightly longer

function prism1(sec,path)=[for(p=path)trns([0,0,p.y],offset(sec,rnd(p.x,3)))];

//function to calculate angle of a 2d vector starting from origin and end point with x and y co-ordinates
// example:
// p1=[3,4];p2=[-3,2];
// v=p2-p1;
// p_lineo([p1,p2],.2);
// ang= ang(v.x,v.y);
// echo(ang);

function ang(x,y)= x>=0&&y>=0?atan(y/x):x<0&&y>=0?180-abs(atan(y/x)):x<0&&y<0?180+abs(atan(y/x)):360-abs(atan(y/x));

//function to rotate a point around a vector(axis) with angle theta           
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
// function is input to another function q_rot
function qmr1(s,r,pl,n=0)= n==len(s)?pl:
qmr1(s,r,
    let(
    v1=s[n]=="x"?[1,0,0]:s[n]=="y"?[0,1,0]:[0,0,1],
    r1=r[n]==undef?0:r[n])
    [for(p=pl)q(v1,p,r1)],n+1);
//function is input to another function q_rot
function qmr2(s,r,pl,n=0)= n==len(s)?pl:qmr2(s,r,let(
    v1=s[n]=="x"?[1,0,0]:s[n]=="y"?[0,1,0]:[0,0,1],
    r1=r[n]==undef?0:r[n])
[for(i=[0:len(pl)-1])[for(p=pl[i])q(v1,p,r1)]],n+1);
  
//function to rotate a group of points "pl" around a series of axis with defined angles e.g q_rot(s=["z20","x40","y80"],pl=[[2,0],[10,2]])=> will rotate the line first around z axis by 20 deg then around x axis by 40 degrees and then around y axis by 80 degrees.  
function q_rot(s,pl)= is_num(pl[0][0])?qmr1([for(p=s)cvar(p)[0]],[for(p=s)cvar(p)[1]],pl):qmr2([for(p=s)cvar(p)[0]],[for(p=s)cvar(p)[1]],pl);

//function used as input to sort
function sort1(list,n=0)=
let(
list1=[for(i=[0:len(list)-1])[list[i]+i*.0000000001,i]],
a=lookup(min(list1*[1,0]),list1),
list2=[for(i=[0:len(list1)-1])if (lookup(list1[i].x,list1)!=a)list1[i]]
)n==0?min(list1*[1,0]):sort1(list2*[1,0],n-1);

// function to sort a list of real numbers in ascending order
function sort(list)=[for(i=[0:len(list)-1])sort1(list,i)];
       
//function to make surface with a polyline 2d sketch and a 3d path(there is no render here but points can be visualised with following command for(p=surf_extrude(sec,path))points(p,.2);) 
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

    )each i<len(path)-2?[sec1]:[sec1,sec2]];
    
//module to render surface with a polyline 2d sketch and a 3d path. thickness of the surface can be set with parameter "t". positive and negative value creates thickness towards +z and -z directions respectively 
    
 module surf_extrude(sec,path,t=.01){
     
     surf=surf_extrude(sec,path);
     surf1=trns([0,0,t],surf);
     for(i=[0:len(surf)-2])
         for(j=[0:len(surf[i])-2])
           if(t>0)
             swp([[surf1[i][j],surf1[i+1][j],surf1[i+1][j+1],surf1[i][j+1]],[surf[i][j],surf[i+1][j],surf[i+1][j+1],surf[i][j+1]]]);
         else
             swp([[surf[i][j],surf[i+1][j],surf[i+1][j+1],surf[i][j+1]],[surf1[i][j],surf1[i+1][j],surf1[i+1][j+1],surf1[i][j+1]]]);}
              
    
//function to convert the y co-ordinates to z co-ordinates e.g.[x,y]=>[x,0,y]. 2d to 3d coordinate system
    
function cytz(path)=[for(p=path)[p.x,0,p.y]];

//function for creating points in circle with radius "r", center point "cp" and number of segments "s"           
function cir(r,cp=[0,0],s=50)=[for(i=[0:360/s:360-360/s])[cp.x+r*cos(i),cp.y+r*sin(i)]];

//module for drawing a closed 2d polyline from a group of points "path" and width of the polyline is defined by parameter "size".
module p_line(path,size=.5){
    for(i=[0:len(path)-1])
        let(p0=path[i],p1=i<len(path)-1?path[i+1]:path[0])
    
    hull(){
    translate(p0)circle(size/2,$fn=20);
    translate(p1)circle(size/2,$fn=20);}}
    
//module for drawing an open 2d polyline from a group of points "path" and width of the polyline is defined by parameter "size".
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
//function to identify whether the section is clockwise or counter clockwise. cw(sec)==1 means clockwise and -1 means counterclockwise. e.g. echo(cw([[0,0],[4,0],[0,4],[-4,0]]));// -1

function cw(sec)=let(p=mode_sign(list_ang(sec)))
p[0]>p[1]?-1:1;

//function to calculate the intersection point between 2 lines e.g. echo(i_p2d(l1=[[0,0],[1,4]],l2=[[10,0],[7,2]])); => //ECHO: [1.42857, 5.71429]
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

//function to calculate intersection point between 2 lines in 3d space (mostly if these lines lie on the same plane)
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

//function to draw points in circular arc with radius, start angle "ang1" , end angle "ang2", center point of the arc "cp" and number of segments required in the arc "s". e.g. following code will draw an arc of radius 5 from 0 to 90 degrees centered at [0,0] with 20 segments in the arc: p_lineo(arc(radius=5,ang1=0,ang2=90,cp=[0,0],s=20),.1);
    
function arc(radius,ang1=0,ang2=355,cp=[0,0],s=20)=[for(i=[ang1:(ang2-ang1)/s:ang2])cp+[radius*cos(i),radius*sin(i)]];
    
//function to create 2d fillet between 2 circles, where r1,r2 and c1,c2 are radiuses and enter points of the 2 circles respectively. r-> fillet radius
////example:
//%p_line(cir(5),.2);
//%p_line(cir(3,[7,0]),.2);
//fillet=2cir_fillet(r1=5,r2=3,c1=[0,0],c2=[7,0],r=1);
//p_line(fillet,.2);
    
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

//function to draw the fillet radius "r" between the 2 circle with radiuses "r1" and "r2" centered at "c1" and "c2" respectively.This function gives an additional flexibility for drawing fillet only one side. e.g try following example
//fillet=2cir_filleto(r1=10,r2=10,c1=[0,0],c2=[20,0],r=10);
//p_lineo(fillet[0],.1);

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

//function to rotate a vector by "theta" degrees e.g. try following code:
//line=[[0,0],[5,3]];
//line1=line*rm(30);
//
//p_lineo(line,.1);
//p_lineo(line1,.1);

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

//function to create section with corner radiuses. e.g. following code has 3 points at [0,0],[10,0] and [7,15] and radiuses of 0.5,2 and 1 respectively,s=5 represent the number of segments at each corner radius.
//sec=cr(pl=[[0,0,.5],[10,0,2],[7,15,1]],s=5);
//p_line(sec,.1);

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
    
//function to translate a group of points "sl" by "m" distance defined in [x,y,z].e.g. try following code:
//sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5);
//p_line3dc(trns([2,5,10],sec),.1);
 
function trns(m,sec)=is_num(sec[0][0])?tr1(m,sec):tr2(m,sec);

//function to scale a 2d section by an amount "sl" which has to be >0 (keeps the y-coordinates same). e.g.following code scales the section by 0.7 (70% of the original shape)
//sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5);
//p_line(sec,.1);
//p_line(scl2d(sec,.7),.1);

function scl2d(sec,sl)=let(
 cp=avg_v(sec),
 rev=[for(p=sec)cp+(p-cp)*sl],
 y1=cp-[0,each min(sec*[0,1])],
 y2=cp-[0,each min(rev*[0,1])],
 d=y2-y1
 )c3t2(trns(d,rev));

// //function to scale a 2d section by an amount "sl" which has to be >0 (keeps the revised section in center). e.g.following code scales the section by 0.7 (70% of the original shape)
//sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5);
//p_line(sec,.1);
//p_line(scl2d_c(sec,.7),.1);
 
function scl2d_c(sec,sl)=let(
 cp=avg_v(sec),
 rev=[for(p=sec)cp+(p-cp)*sl]
 )rev;

// function to scale a 3d prism keeping the base z-coordinate same. takes 2 arguments "prism" to scale and the scaling factor "s". scale factor can take any real number negative values will scale the prism and turn the prism upside down.
// try the following code to understand better:
// sec=cir(10);
// path=cr(pts1([[2,0],[-2,0,2],[0,10,3],[-3,0]]),5);
// prism=prism(sec,path);
// %swp(prism);
// swp(scl3d(prism,.7));

function scl3d(prism,s=1)=let(
cp=avg_v(prism),
rev=[for(p=prism)[for(p1=p)cp+(p1-cp)*s]],
flat_p1=[for(p=prism)each p],
flat_p2=[for(p1=rev)each p1],
z1=min(flat_p1*[0,0,1]),
z2=min(flat_p2*[0,0,1]),
d=z1-z2
)trns([0,0,d],rev);

// function to scale a 3d prism keeping the prism centered. takes 2 arguments "prism" to scale and the scaling factor "s". scale factor can take any real number negative values will scale the prism and turn the prism upside down.
// try the following code to understand better:
// sec=cir(10);
// path=cr(pts1([[2,0],[-2,0,2],[0,10,3],[-3,0]]),5);
// prism=prism(sec,path);
// %swp(prism);
// swp(scl3d_c(prism,.7));

function scl3d_c(prism,s=1)=let(
cp=avg_v(prism),
rev=[for(p=prism)[for(p1=p)cp+(p1-cp)*s]]
)rev;

// used as input for function ip()
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

// function to calculate intersection point between two 3d prisms. "prism" is the 3d object which is intersected with "prism1".
// try below code for better understanding:
// sec=cir(10);
// path=cr(pts1([[2,0],[-2,0,2],[0,10,3],[-9.9,0]]),5);
// prism=prism(sec,path);
// prism1=q_rot(["y40"],cyl(r=3,h=15,s=30));
//
// %swp(prism);
// %swp(prism1);
// ip=ip(prism,prism1);
// points(ip,.2);
          
function ip(prism,prism1)=let(sec=ipa(prism,prism1))[for(i=[0:len(sec)-1])let(i_plus=i<len(sec)-1?i+1:0)if(norm(sec[i]-sec[i_plus])>.1)sec[i]];

//This function is not very often used and may get removed
          
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
            
// used as input to another function

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
    
// function to draw normal vector to a given vector "v". Not used very often and may be removed    
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
 
//function to convert 3d to 2d, it just removes the z-coordinate from the points list 
function c3t2(sec)=is_undef(sec.x.x)?[sec.x,sec.y]:[for(p=sec)[p.x,p.y]];

// function to calculate the cumulative sum of all the points of a 2d or 3d points list.
// e.g.
// sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5);
// echo(sum_v(sec)); //[95.9558, 82.3332]

function sum_v(prism)=let(
decision=is_num(prism.x.x)?0:1,
sum=decision==0?
[sum(sec*[1,0]),sum(sec*[0,1])]:
let(cg=[for(p=prism)[sum(p*[1,0,0]),sum(p*[0,1,0]),sum(p*[0,0,1])]])[sum(cg*[1,0,0]),sum(cg*[0,1,0]),sum(cg*[0,0,1])]
)sum;

//function to draw tangent line joining 2 circles with radiuses "r1" and "r2" with center points "cp1" and "cp2" respectively. This function draws tangent line only one side
// e.g. try this code below:
// sec=2ctp(r1=10,r2=5,cp1=[0,0],cp2=[15,6]);
// p_line(cir(10),.1);
// p_line(cir(5,[15,6]),.1);
// p_line(sec,.1);

function 2ctp(r1,r2,cp1,cp2)=
let(
v1=cp2-cp1,
u1=v1/norm(v1),
ang1=asin((r2-r1)/norm(cp2-cp1)),

t1=cp1+u1*r1*rm(90+ang1),
t2=cp2+u1*r2*rm(90+ang1),

t3=cp1+u1*r1*rm(-90-ang1),
t4=cp2+u1*r2*rm(-90-ang1))[t1,t2];

//function to draw tangent line joining 2 circles with radiuses "r1" and "r2" with center points "cp1" and "cp2" respectively. This function draws tangent line on both the sides
// e.g. try this code below:
// sec=2ctpf(r1=10,r2=5,cp1=[0,0],cp2=[15,6]);
// p_line(cir(10),.1);
// p_line(cir(5,[15,6]),.1);
// p_line(sec,.1);

function 2ctpf(r1,r2,cp1,cp2)=
let(
v1=cp2-cp1,
u1=v1/norm(v1),
ang1=asin((r2-r1)/norm(cp2-cp1)),

t1=cp1+u1*r1*rm(90+ang1),
t2=cp2+u1*r2*rm(90+ang1),

t3=cp1+u1*r1*rm(-90-ang1),
t4=cp2+u1*r2*rm(-90-ang1))[t1,t2,t4,t3];

//module to draw a polyline in 3d space (loop not closed)
// e.g. try following code:
// sec=trns([5,10,6],q_rot(["x45"],cir(10)));
// p_line3d(sec,.2);
    
module p_line3d(path,r,rec=0){
    for(i=[0:len(path)-2])
        
    hull(){
    translate(path[i])if(rec==0)sphere(r); else cube(r*2,true);
    translate(path[i+1])if(rec==0)sphere(r);else cube(r*2,true);
    }}

//module to draw a polyline in 3d space (loop closed)
// e.g. try following code:
// sec=trns([5,10,6],q_rot(["x45"],cir(10)));
// p_line3dc(sec,.2);    

module p_line3dc(path,r,rec=0){
    for(i=[0:len(path)-1])
        let(
    i_plus=i<len(path)-1?i+1:0
    )
    hull(){
    translate(path[i])if(rec==0)sphere(r); else cube(r*2,true);
    translate(path[i_plus])if(rec==0)sphere(r);else cube(r*2,true);
    }}
  
//function used as input to another function
  
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
            
//function used as input to another function 
 
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
            
// function used as input to another function
            
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
            if(! is_undef(p7[0]))3p_3d_fillet_wo_pivot(p7[0],p1[i],p2[i],r,s)
            
            ];
 
// function for creating fillet: this function first finds the intersection point between prism and prism1 and then calculates the fillet with radius "r". option "0" and "1" creates fillet either outside or inside.parameter "s" is for number of segments in the fillet
// an example below will be more clear (try changing option from 1 =>0 or flip the direction of prism1 by flip(prism1))
// try below code for better understanding:
// sec=cir(10);
// path=cr(pts1([[2,0],[-2,0,2],[0,10,3],[-9.9,0]]),5);
// prism=prism(sec,path);
// prism1=q_rot(["y40"],cyl(r=3,h=15,s=30));
//
// %swp(prism);
// %swp(prism1);
// fillet=ipf(prism,prism1,r=1,option=1,s=5);
// swp_c(fillet);

  function ipf(prism,prism1,r,option=0,s=5)=let(sec=ipr(prism,prism1,r,option,s=s))
            [for(i=[0:len(sec)-1])each i<len(sec)-1?[sec[i]]:[sec[i],sec[0]]];
            
 //function for creating a fillet by intersection between a plane and a prism
//example: 
// sec=cr(pts1([[0,0,1],[8,3,3],[5,7,1],[-8,0,2],[-5,20,1]]),20);
// prism=l_extrude(sec,30);
// plane1=plane([0,0,1],60);
// %swp(prism);
// %swp(plane1);
// prism1=ipe(plane1,prism,r=2,option=1,s=10);
//
// swp(prism1);
//
// plane2=trns([0,0,20],plane([0,0,1],60));
// prism2=ipe(plane2,flip(prism),r=1,option=1,s=10);
// %swp(plane2);
// swp(prism2);
 
 function ipe(prism,prism1,r,option=0,s=5)=
 let(
 sec=ipr1(prism,prism1,r,option,s=s),
 sec1=[for(i=[0:len(sec[0])-1])[for(p=sec)p[i]]]
 )sec1;

// draws a cylinder try swp(cyl(r=5,h=15)); 

function cyl(r1=1,r2=1,h=1,cp=[0,0],s=50,r,d,d1,d2,center=false)=let(
     ra=is_num(r)?r:is_num(d)?d/2:is_num(d1)?d1/2:r1,
     rb=is_num(r)?r:is_num(d)?d/2:is_num(d2)?d2/2:r2,
     sec=cir(ra,cp,s),
     
path=pts([[-ra+.1,0],[ra-.1,0],[rb-ra,h],[-rb+.1,0]]),
    prism=center==true?trns([0,0,-h/2],prism(sec,path)):prism(sec,path))
    prism;
    
// function flips the direction of points of 2d section or 3d prism
            
function flip(sec)=[for(i=[len(sec)-1:-1:0])sec[i]];
     
// function for linear extrude a section by height "h", also the section can be rotated by an angle "a" in number of steps "steps"
// try following code for better understanding (also try changing "a" and "steps"):
// sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5);
// prism=l_extrude(sec,h=15,a=0,steps=1);
// swp(prism);

function l_extrude(sec,h=1,a=0,steps=1)=[for(i=[0:a==0?1:(a-0)/steps:a==0?1:a])
    trns([0,0,a==0?h*i:h/a*i],q_rot([str("z",a==0?0:i)],sec))];
 
// function to draw a rectangle
// e.g. p_line(sqr([10,5]),.1); or polygon(sqr([10,5]));
 
function sqr(s,center=false)=
let(
m=is_num(s)?s:s.x,
n=is_num(s)?s:s.y,
sec=[[0,0],[m,0],[m,n],[0,n]],
sec1=center==true?[for(p=sec)p-[m/2,n/2]]:sec)
    sec1;

// function to draw cube
// swp(cub(p=[10,5,4]));
    
function cub(p,center=false)=
let(
m=is_num(p)?p:p.x,
n=is_num(p)?p:p.y,
o=is_num(p)?p:p.z,

path=pts([[-m/2,0],[m/2,0],[0,o],[-m/2,0]]),
prism=center==true?trns([-m/2,-n/2,-o/2],rsz3d(prism(sqr(m),path),[m,n,o])):rsz3d(prism(sqr(m),path),[m,n,o])
)
prism;

// function for creating sphere with radius "r", center point "cp" and number of segments "s".
// try following code:
// swp(spr(r=3,cp=[4,5,6],s=30));

function spr(r,cp=[0,0,0],s=50)=let(
path=arc(r,-90,90,s=s),
prism=[for(p=path)trns([0,0,p.y]+cp,cir(p.x,s=s))])
    prism;

//function is used as input to another function    
    
function add_p(p,p1=[0,0],n,i=0)= n==0?p1:add_p(p,[p[i].x+p1.x,p[i].y+p1.y],n-1,i+1);

// function is used like a turtle move to create 2d shapes.
// following example will create a rectangle with sides 10 x 5:
// sec=pts([[0,0],[10,0],[0,5],[-10,0]]); // starts at [0,0] then moves 10 units to +x direction then moves 5 units towards +y direction and then moves 10 units to -x direction
// p_line(sec,.1);

function pts(p)=[for(n=[1:len(p)])add_p(p=p,p1=[0,0],n=n,i=0)];

// function is used as input to another function
    
function add_p1(p,p1=[0,0,0],n,i=0)= n==0?p1:add_p1(p,[p[i].x+p1.x,p[i].y+p1.y,p[i].z],n-1,i+1);

//same as pts(p) with only difference that it keeps the z value unchanged
// for example:
// sec=pts1([[0,0,1],[10,0,1],[0,5,1],[-10,0,1]]); // starts at [0,0] then moves 10 units to +x direction then moves 5 units towards +y direction and then moves 10 units to -x direction
//  echo(sec); // ECHO: [[0, 0, 1], [10, 0, 1], [10, 5, 1], [0, 5, 1]]
//  this function is mainly used with function cr(pl,s) (please see the example of function cr(pl,s))

function pts1(p)=[for(n=[1:len(p)])add_p1(p=p,p1=[0,0,0],n=n,i=0)];
 
// function is used as input to another function
 
function add_p2(p,p1=[0,0,0,0],n,i=0)= n==0?p1:add_p2(p,[p[i][0]+p1[0],p[i][1]+p1[1],p[i][2]+p1[2],p[i][3]],n-1,i+1);

// same as pts and pts1 and is used as a turtle movement in 3d space and keeps the 4th point same to be used as radius for rounding in function cr3d() (check example fo function cr3d()

function pts2(p)=[for(n=[1:len(p)])add_p2(p=p,p1=[0,0,0,0],n=n,i=0)];

// module for rendering points along the various shapes 2d or 3d. parameter "d" is the size of cube which is used as point. a list has to be provided for parameter "p"
// try following code:
// sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5);
// prism=l_extrude(sec,h=15,a=90,steps=20);
// %swp(prism);
// for(p=prism) points(p,.2);
 
module points(p,d=.5){
    for(i=p)translate(i)cube(size=d,center=true);
    
    }
 
// function to get the minimum radius for a defined section
// example:
// sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5);
// echo(sec_r(sec)); //=>ECHO: 0.5


function sec_r(sec)=let(

list_of_radius=[for(i=[0:len(sec)-1])
  let(
  i_minus=i==0?len(sec)-1:i-1,
  i_plus=i<len(sec)-1?i+1:0,
  p0=sec[i_minus],p1=sec[i],p2=sec[i_plus],
  r=3p_r(p0,p1,p2))r]
)
min(list_of_radius);

// math function to calculate the determinant of a 3 x 3 matrix

function det3d(m)=let(
m11=m[0][0],m12=m[0][1],m13=m[0][2],
m21=m[1][0],m22=m[1][1],m23=m[1][2],
m31=m[2][0],m32=m[2][1],m33=m[2][2],

s11=m22*m33-m32*m23,s12=-(m21*m33-m31*m23),s13=m21*m32-m31*m22,
s21=-(m12*m33-m32*m13),s22=m11*m33-m31*m13,s23=-(m11*m32-m31*m12),
s31=m12*m23-m22*m13,s32=-(m11*m23-m21*m13),s33=m11*m22-m21*m12,

d=m11*s11+m12*s12+m13*s13
)d;

// math function to calculate the determinant of a 2 x 2 matrix

function det2d(m)=let(
m11=m[0][0],m12=m[0][1],
m21=m[1][0],m22=m[1][1],

s11=m22,s12=-m21,
s21=-m12,s22=m11,

d=m11*m22-m21*m12
)d;

// math function to calculate the inverse of a 3 x 3 matrix
// example:
// v1=[2,3,4];
// v2=[3,4,1];
// v3=[4,5,6];
// echo(i_m3d(t([v1,v2,v3])));// =>ECHO: [[-2.375, 1.75, 0.125], [-0.25, 0.5, -0.25], [1.625, -1.25, 0.125]]

function i_m3d(m)=let(
m11=m[0][0],m12=m[0][1],m13=m[0][2],
m21=m[1][0],m22=m[1][1],m23=m[1][2],
m31=m[2][0],m32=m[2][1],m33=m[2][2],

s11=m22*m33-m32*m23,s12=-(m21*m33-m31*m23),s13=m21*m32-m31*m22,
s21=-(m12*m33-m32*m13),s22=m11*m33-m31*m13,s23=-(m11*m32-m31*m12),
s31=m12*m23-m22*m13,s32=-(m11*m23-m21*m13),s33=m11*m22-m21*m12,

d=m11*s11+m12*s12+m13*s13
) 1/d*[[s11,s21,s31],[s12,s22,s32],[s13,s23,s33]];

// math function to calculate the inverse of a 2 x 2 matrix
// example:
// v1=[2,3];
// v2=[3,4];
// echo(i_m2d(t([v1,v2]))); //=> ECHO: [[-4, 3], [3, -2]]

function i_m2d(m)=let(
m11=m[0][0],m12=m[0][1],
m21=m[1][0],m22=m[1][1],

s11=m22,s12=-m21,
s21=-m12,s22=m11,

d=m11*m22-m21*m12
)1/d*[[s11,s21],[s12,s22]];

// function is used as input to bezier curve function

function add_v(v,s=[0,0],n=0)=n==len(v)?s:add_v(v,s+v[n],n+1);

// math function to calculate factorial of a number

function fact(n,m=1)=n==0?m:fact(n-1,m*n);

// math function to calculate number of possible combinations for "n" items with "i" selected items

function comb(n,i)=fact(n)/(fact(i)*fact(n-i));

//function for calculating bezier curve with control points "p" and with number of segments 1/s
// example:
// p=[[0,0],[10,5],[0,15],[12,20]]; 
// b=bez(p,.1); 
// points(b,.5);
// //control points
// color("green")
// points(p,.5);

function bez(p,s=10)=[for(t=[0:1/s:1])
    let(n=len(p)-1)add_v([for(i=[0:n])comb(n,i)*(1-t)^(n-i)*t^i*p[i]])];

// function for creating arc which is tangent to 2 circles
// try this code as an example:
// sec=2cir_tarc(10,5,[0,0],[20,5],20);
// p_lineo(sec,.2);
// p_line(cir(10),.2);
// p_line(cir(5,[20,5]),.2);
    
function 2cir_tarc(r1,r2,cp1,cp2,r)=
assert(r>=(r1+r2+norm(cp2-cp1))/2,str("arc radius : ",r," is smaller than the minimum required radius of ",(r1+r2+norm(cp2-cp1))/2))
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
    
// function creates a shortest 2d arc with 2 points with a radius "r" and number of segments "s". parameter cw(clockwise=1 and counter clockwise=-1) defines the order of arc
//try this example for better understanding:
// sec=2p_arc(p1=[2,3],p2=[6,5],r=2.25,cw=-1,s=20);
// p_lineo(sec,.2);
    
function 2p_arc(p1,p2,r,cw=1,s=20)=
assert(r>=norm(p2-p1)/2,str("radius : ",r," is smaller than ",norm(p2-p1)/2))
let(
p3=p1+(p2-p1)/2,
d=norm(p3-p1),
l=sqrt(r^2-d^2),
v=p1-p3,u=v/norm(v),
cp=p3+u*l*rm(cw==-1?-90:90),
v1=p1-cp,v2=p2-cp,
a1=ang(v1.x,v1.y),a2=ang(v2.x,v2.y),
a3=cw==-1?(a2<a1?a2+360:a2):(a2<a1?a2:a2-360)

)arc(r,a1,a3,cp,s);

// function to calculate the center point for arc where 2 points "p1" and "p2" and radius "r" are known (clockwise and counter clockwise will have different center points
// example:
// pnt=2p_arc_cp(p1=[2,3],p2=[6,5],r=5,cw=-1);
// points([pnt],.5);

function 2p_arc_cp(p1,p2,r,cw=1)=
assert(r>=norm(p2-p1)/2,str("radius : ",r," is smaller than ",norm(p2-p1)/2))
let(
p3=p1+(p2-p1)/2,
d=norm(p3-p1),
l=sqrt(r^2-d^2),
v=p1-p3,u=v/norm(v),
cp=p3+u*l*rm(cw==-1?-90:90),
v1=p1-cp,v2=p2-cp,
a1=ang(v1.x,v1.y),a2=ang(v2.x,v2.y),
a3=cw==-1?(a2<a1?a2+360:a2):(a2<a1?a2:a2-360)
)cp;

// function creates a longest 2d arc with 2 points with a radius "r" and number of segments "s". parameter cw(clockwise=1 and counter clockwise=-1) defines the order of arc
//try this example for better understanding:
// sec=2r(p1=[2,3],p2=[6,5],r=3,cw=-1,s=20);
// p_lineo(sec,.2);

function 2r(p1,p2,r,cw=1,s=20)=let(
p3=p1+(p2-p1)/2,
d=norm(p3-p1),
l=sqrt(r^2-d^2),
v=p1-p3,u=v/norm(v),
cp=p3+u*l*rm(cw==-1?90:-90),
v1=p1-cp,v2=p2-cp,
a1=ang(v1.x,v1.y),a2=ang(v2.x,v2.y),
a3=cw==-1?(a2<a1?a2+360:a2):(a2<a1?a2:a2-360)

)arc(r,a1,a3,cp,s);

// function to create arc with 3 points in 2d
// example:
// sec=3p_arc([1,2],[3,7],[7,3]);
// p_lineo(sec,.2);
// points([[1,2],[3,7],[7,3]],.5);

function 3p_arc(p1,p2,p3,s=30)=

let(
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


function 3p_cir(p1,p2,p3,s=30)=

let(
p4=p1+(p2-p1)/2,
p5=p2+(p3-p2)/2,
v1=p2-p4,u1=v1/norm(v1),
v2=p3-p5,u2=v2/norm(v2),
p6=p4+u1*rm(90),
p7=p5+u2*rm(90),
cp=i_p2d([p4,p6],[p5,p7]),
r=norm(p1-cp),
//v3=p1-cp,v4=p2-cp,v5=p3-cp,
//a1=ang(v3.x,v3.y),
//a2=ang(v4.x,v4.y),
//a3=ang(v5.x,v5.y),
//a4=cw([p1,p2,p3])==-1?(a3<a1?a3+360:a3):(a3<a1?a3:a3-360)

)arc(r,0,360,cp,s);
// function to draw an ellipse with semi-major and semi-minor axis "r1" and "r2" respectively and with center "cp" and number of segment "s"
// example:
// sec=ellipse(r1=5,r2=3,cp=[2,3],s=30);
// p_line(sec,.2);

function ellipse(r1,r2,cp,s=30)=
let(
sec=[for(i=[0:360/s:360-360/s])cp+[r1*cos(i),r2*sin(i)]]
)sec;

// experimental function
// example:
// sec=l_cir_fillet(line=[[0,0],[0,20]],r1=5,r2=1,cp=[5,10]);
// p_lineo(sec,.2);

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

// function to calculate average of a group of points either 2d or 3d
// example:
// sec=cir(10);
// path=cr(pts1([[2,0],[-2,0,2],[0,10,3],[-3,0]]),5);
// prism=prism(sec,path);
// %swp(prism);
// avg=avg_v(prism);
// echo(avg);
// points([avg],.5);
 
function avg_v(prism)=let(
decision=is_num(prism.x.x)?0:1,
cp=decision==0?
[sum(prism*[1,0])/len(prism),sum(prism*[0,1])/len(prism)]:
let(cg=[for(p=prism)[sum(p*[1,0,0])/len(p),sum(p*[0,1,0])/len(p),sum(p*[0,0,1])/len(p)]])[sum(cg*[1,0,0])/len(cg),sum(cg*[0,1,0])/len(cg),sum(cg*[0,0,1])/len(cg)]
)cp;

// function to calculate the resized prism
//example:
// sec=cir(10);
// path=cr(pts1([[2,0],[-2,0,2],[0,10,3],[-3,0]]),5);
// prism=prism(sec,path);
// %swp(prism);
// resized_prism=rsz3d(prism,[5,5,7]);
// swp(resized_prism);

function rsz3d(prism,rsz=[1,1,1])=
let(
rev_p_list=[for(p=prism) each[for(p1=p)p1]],
max_x=max(rev_p_list*[1,0,0]),
max_y=max(rev_p_list*[0,1,0]),
max_z=max(rev_p_list*[0,0,1]),
min_x=min(rev_p_list*[1,0,0]),
min_y=min(rev_p_list*[0,1,0]),
min_z=min(rev_p_list*[0,0,1]),
avg=avg_v(prism),

r_x=rsz.x/(max_x-min_x),
r_y=rsz.y/(max_y-min_y),
r_z=rsz.z/(max_z-min_z),
rev_prism=[for(i=[0:len(prism)-1])
    [for(p=prism[i])avg+[r_x*(p.x-avg.x),r_y*(p.y-avg.y),r_z*(p.z-avg.z)]]],
t=(bb(rev_prism)-bb(prism))/2
    )trns(t,rev_prism);
 
// function to calculate the resized prism- centered
//example:
// sec=cir(10);
// path=cr(pts1([[2,0],[-2,0,2],[0,10,3],[-3,0]]),5);
// prism=prism(sec,path);
// %swp(prism);
// resized_prism=rsz3dc(prism,[5,5,7]);
// swp(resized_prism);
 
function rsz3dc(prism,rsz=[1,1,1])=
let(
rev_p_list=[for(p=prism) each[for(p1=p)p1]],
max_x=max(rev_p_list*[1,0,0]),
max_y=max(rev_p_list*[0,1,0]),
max_z=max(rev_p_list*[0,0,1]),
min_x=min(rev_p_list*[1,0,0]),
min_y=min(rev_p_list*[0,1,0]),
min_z=min(rev_p_list*[0,0,1]),
avg=avg_v(prism),

r_x=rsz.x/(max_x-min_x),
r_y=rsz.y/(max_y-min_y),
r_z=rsz.z/(max_z-min_z),
rev_prism=[for(i=[0:len(prism)-1])
    [for(p=prism[i])avg+[r_x*(p.x-avg.x),r_y*(p.y-avg.y),r_z*(p.z-avg.z)]]],
t=(bb(rev_prism)-bb(prism))/2
    )rev_prism;

// function to calculate 2d resized section - placed at minimum y value
// example:
// sec=cir(10);
// rsz_sec=rsz(sec,[5,3]);
// %p_line(sec,.2);
// p_line(rsz_sec,.2);
    
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

// function to calculate 2d resized section - centered
// example:
// sec=cir(10);
// rsz_sec=rsz_c(sec,[5,3]);
// %p_line(sec,.2);
// p_line(rsz_sec,.2);
  
function rsz_c(sec,rsz=[1,1,1])=
  let(
  avg=avg_v(sec),
  max_x=len(rsz)==2?max(sec*[1,0]):max(sec*[1,0,0]),
  min_x=len(rsz)==2?min(sec*[1,0]):min(sec*[1,0,0]),
  max_y=len(rsz)==2?max(sec*[0,1]):max(sec*[0,1,0]),
  min_y=len(rsz)==2?min(sec*[0,1]):min(sec*[0,1,0]),

  r_x=rsz.x/(max_x-min_x),
  r_y=rsz.y/(max_y-min_y)

  )[for(i=[0:len(sec)-1])let(
  p=avg+[r_x*(sec[i].x-avg.x),r_y*(sec[i].y-avg.y)]
  )p];

// function to create a line with number of segments "s"
// example:
// line=l([[0,0],[4,3]],10);
// points(line,.2);  
  
function l(l,s=20)=
let(
p0=l[0],p1=l[1],
v=p1-p0,u=v/norm(v),
length=norm(v)
)[for(i=[0:length/s:length])p0+u*i];
    
// function to find radius of arc with 3 known points in 2d
// example:
// radius=3p_r([1,2],[3,7],[7,3]);
// echo(radius); //=> ECHO: 3.30892
    
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

// function to get the minimum radius for a defined section
// example:
// sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5);
// echo(min_r(sec)); //=>ECHO: 0.5
  
  
function min_r(sec)=
min([for(i=[0:len(sec)-1])3p_r(sec[i==0?len(sec)-1:i-1],sec[i],sec[i<len(sec)-1?i+1:0])]);
    
// function for placing multiple points on the straight line segments of a closed loop section. parameter "sl" is for placing points with pitch distance defined by "sl"
// example:
// sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5);
// points(sec,.2);
// 
// translate([15,0])
// points(m_points(sec,2),.2);// segment length=> 2 units  
    
function m_points(sec,sl=20)=
[for(i=[0:len(sec)-1])let(
p0=sec[i],
p1=sec[i<len(sec)-1?i+1:0],
lnth=norm(p1-p0),
sec1=lnth>sl?l([p0,p1],lnth/sl):[p0],
sec2=[for(i=[0:len(sec1)-1])if(sec1[i]!=sec1[i<len(sec1)?i+1:0])sec1[i]])
each sec2];

// function for placing multiple points on the straight line segments of an open section. parameter "sl" is for placing points with pitch distance defined by "sl"
// example:
// sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5);
// points(sec,.2);
// 
// translate([15,0])
// points(m_points_o(sec,2),.2);// segment length=> 2 units

function m_points_o(sec,sl=20)=
[for(i=[0:len(sec)-2])let(
p0=sec[i],
p1=sec[i+1],
lnth=norm(p1-p0),
sec1=lnth>sl?l([p0,p1],lnth/sl):[p0],
sec2=[for(i=[0:len(sec1)-1])if(sec1[i]!=sec1[i<len(sec1)?i+1:0])sec1[i]])
each sec2];

// function for calculating multiple points on the straight line segments of a closed section. sec-> closed section; s -> number of segments for each straight line segment of closed section; m-> minimum segment length, if the derived segment length < m, then it is omitted
// sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5);
// points(sec,.2);
// 
// translate([15,0])
// points(m_points_sc(sec,s=5,m=.5),.2);// number of segments=> 5

function m_points_sc(sec1,s,m=.5)=
let(
l=[for(i=[0:len(sec1)-1])
let(
i_plus=i<len(sec1)-1?i+1:0,
l=norm(sec1[i_plus]-sec1[i]),
u=uv(sec1[i_plus]-sec1[i])
)each l/s>=m?[for(j=[0:l/s:l])sec1[i]+j*u]:[sec1[i]]]
)l;

// function for calculating multiple points on the straight line segments of an open section. sec-> closed section; s -> number of segments for each straight line segment of closed section; m-> minimum segment length, if the derived segment length < m, then it is omitted
// sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5);
// points(sec,.2);
// 
// translate([15,0])
// points(m_points_so(sec,s=5,m=.5),.2);// number of segments=> 5

function m_points_so(sec1,s,m=.5)=
let(
l=[for(i=[0:len(sec1)-2])
let(
i_plus=i+1,
l=norm(sec1[i_plus]-sec1[i]),
u=uv(sec1[i_plus]-sec1[i])
)each l/s>=m?[for(j=[0:l/s:l])sec1[i]+j*u]:[sec1[i]]]
)l;

// used as input to another function

function cum_sum(list,list1,n,s=1)=n==0?list1:cum_sum(list,[for(i=[0:s])list[i]]*[for(i=[0:s])1],n-1,s+1);

// experimental and need more work

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

// module for rendering various 3d prism 
// //example1:
// sec=cir(10);
// path=cr(pts1([[2,0],[-2,0,2],[0,10,3],[-3,0]]),5);
// prism=prism(sec,path);
// swp(prism);
// //example2:
// prism1=l_extrude(sqr([10,6]),15);
// translate([13,0,0])
// swp(prism1);
// //example3:
// sec2=cr(pts1([[0,0,1],[5,0,1],[-2.5,4,1]]),5);
// path2=[for(i=[0:5:360*5])[10*cos(i),10*sin(i),i/360*5]]; 
// prism2=p_extrude(sec2,path2);
// translate([35,0,0])
// swp(prism2);

module swp(surf1)

let(l=len(surf1[0]),
p0=[for(j=[0:len(surf1)-1])each surf1[j]],
p1=[each [for(j=[0:len(surf1)-1])if(j==0)[for(i=[0:l-1])i+j*l]],
each [for(j=[0:len(surf1)-2])each [for(i=[0:l-1])let(i_plus=i<l-1?i+1:0)[i+l*j,i+l+l*j,i_plus+l+l*j,i_plus+l*j]]],
each [for(j=[0:len(surf1)-1])if(j==len(surf1)-1)[for(i=[l-1:-1:0])i+l*j]]
    ]
)
polyhedron(p0,p1,convexity=10);

// module for rendering 3d prisms with closed section
// example:
// sec=cir(10);
// path=cr(pts1([[2,0],[-2,0,2],[0,10,3],[-9.9,0]]),5);
// prism=prism(sec,path);
// prism1=q_rot(["y40"],cyl(r=3,h=15,s=30));
//
// %swp(prism);
// %swp(prism1);
// fillet=ipf(prism,prism1,r=1,option=1,s=5);
// swp_c(fillet);

module swp_c(surf1)

let(l=len(surf1[0]),
p0=[for(j=[0:len(surf1)-1])each surf1[j]],
p1=[
each [for(j=[0:len(surf1)-2])each [for(i=[0:l-1])let(i_plus=i<l-1?i+1:0)[i+l*j,i+l+l*j,i_plus+l+l*j,i_plus+l*j]]]
    ]
)
polyhedron(p0,p1,convexity=10);

//function for creating fillet between 2 cylinders. r1, r2 and cp1,cp2 are the radiuses and center points of 2 cylinders respectively. r -> is the fillet radius. path -> is given for rounding the cylinder edges
//// example
// path=[[0,0],[0,10]];
// %swp(cyl(r=5,h=15));
// %swp(cyl(r=3,h=10,cp=[7,0]));
// swp(2cyl_fillet(5,3,[0,0],[7,0],1,path));

function 2cyl_fillet(r1,r2,cp1,cp2,r,path)=[for(p=path)trns([0,0,p.y],2cir_fillet(r1+p.x,r2+p.x,cp1,cp2,r))];

//function to create 2d fillet between 2 circles (creates fillet only one side), where r1,r2 and c1,c2 are radiuses and enter points of the 2 circles respectively. r-> fillet radius
////example:
//%p_line(cir(5),.2);
//%p_line(cir(3,[7,0]),.2);
//fillet=2cir_fillet1(r1=5,r2=3,c1=[0,0],c2=[7,0],r=1);
//p_line(fillet,.2);
    
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

// function for creating fillet between 2 spheres. r1, r2 and cp1,cp2 are the radiuses and center points of 2 spheres. r-> fillet radius.
//// example:
// swp(spr(r=5,cp=[0,0,0]));
// swp(spr(r=3,cp=[7,0,0]));
//
// fillet=2spr_fillet(r1=5,r2=3,cp1=[0,0,0],cp2=[7,0,0],1);
// swp(fillet);

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

// // function to get intersection point between a line and circle
// // example
//  line=[[0,0],[3,5]];
//  cir=cir(5);
//  %p_line(line,.2);
//  %p_line(cir,.2);
//
//  pnt=l_cir_ip(line,cir);
//  color("green")
//  points(pnt,.5);
//  echo(pnt);

function l_cir_ip(line,cir)=
let(
ip=[for(i=[0:len(cir)-1])let(
i_plus=i<len(cir)-1?i+1:0,
p0=cir[i],p1=cir[i_plus],
pa=line[0],pb=line[1],
v1=p1-p0,v2=pb-pa,

//p0+v1*t1=pa+v2*t2
//v1*t1-v2*t2=pa-p0
t1=(i_m2d(t([v1,-v2]))*(pa-p0))[0],
ip=p0+v1*t1,
u1=uv(v2),u2=uv(ip-pa)
)if(lim(t1,0,1)&&u1==u2)ip]
)ip;
 
// function to offset a line "l" by distance "d" 
//  example
//  line=[[0,0],[3,5]];
//  %p_line(line,.2);
//  p_line(offst_l(line,2),.2);

function offst_l(l,d)=
let(
v=l[1]-l[0],u=v/norm(v),
p0=l[0]+u*d*rm(-90),
p1=l[1]+u*d*rm(-90)
)[p0,p1];
   
// function to find intersection point at a shortest distance between a point and a line
// example
// line=[[0,0],[3,5]];
// point=[-3,5];
//
// %p_line(line,.2);
// %points([point],.3);
//
// p=perp(line,point);
// points([p],.5);
// echo(p);
   
function perp(line,point)=
let(
v1=line[1]-line[0],
u=uv(v1),
u1=u*rm(90),

line1=[point,point+u1],
ip=i_p2d(line,line1)

)ip;
  
// function to create tangent between 2 circle where r1, r2 and cp1, cp2 are the radiuses and center points of the 2 circles respectively. 
// example:
//  tangent=2cir_tangent(5,3,[0,0],[7,0]);
//  %p_line(cir(5),.2);
//  %p_line(cir(3,[7,0]),.2);
//  p_line(tangent,.2);
 
function 2cir_tangent(r1,r2,cp1,cp2)=
let(
v=cp2-cp1,u=v/norm(v),
theta=ang(v.x,v.y),
theta1=asin((r1-r2)/norm(v)),
p1=cp1+u*r1*rm(90-theta1),
p0=cp2+u*r2*rm(90-theta1),
p2=cp1+u*r1*rm(-(90-theta1)),
p3=cp2+u*r2*rm(-(90-theta1))

)[p0,p1,p2,p3];

// function used as input for another function

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

// function to rotate an object around any arbitrary axis
// example
//  sec=cir(10);
//  path=cr(pts1([[2,0],[-2,0,2],[-1,5,3],[-4,0]]),5);
//  prism=trns([15,0],prism(sec,path));
//  prism1=rot([3,4,7],prism,180);
//  swp(prism);
//  swp(prism1);
//  p_line([[0,0,0],[3,4,7]*10],.2);

function rot(axis,prism,ang)=let(
decision=is_num(prism.x.x)?0:1,
rot=decision==0?
[for(p=prism)q(axis,p,ang)]:
[for(p=prism)[for(p1=p)let(pc=len(p1)==2?[p1.x,p1.y,0]:p1)q(axis,pc,ang)]]
)rot;
 
// function used as input to another function c_hull
 
function s_pnt(sec)=let(
y=sec*[0,1],
loc=[for(i=[0:len(y)-1])if(abs(min(y)-y[i])<.001)i],
x=[for(i=loc)sec[i]],
x_min=min(x*[1,0]),
i=search(x_min,x,0,0)[0],
pnt=x[i]

)pnt;//sec[loc];

// function to subtract points from a list of points
// example:
// list=[[1,2,2],[3,4,5],[10,2,9],[11,1,9]];
// list1=[[1,2,2],[10,2,9]];
// revised_list=reduced_list(list,list1);
// echo(revised_list); //=> ECHO: [[3, 4, 5], [11, 1, 9]]

function reduced_list(list,list1)=revised_list(list,[for(p=list1)each each search([p],list,0)]);
 
 // function used as input to another function
 
function list_of_points_to_omit(sec,point)=let(
list=[for(i=[0:len(sec)-1])if(norm(sec[i]-point)<.001)i],
list1=len(list)>1?[for(i=[1:len(list)-1])list[i]]:[]

)list1;

// function used as input to another function

function revised_list(sec,index_list)=let(
a=[for(i=[0:len(sec)-1]) if(search(0,[for(j=index_list)i-j],0)==[])i],
sec1=[for(i=a)sec[i]]
)sec1;

// function used as input to another function

function remove_extra_points(sec,n=0)=
n==len(sec)?sec:remove_extra_points(
let(
a=list_of_points_to_omit(sec,sec[n]),
b=revised_list(sec,a)
)b,n+1
);

// function to match the number of points between 2 sections
// example:
// sec=cr(pts1([[-2.5,-2.5,1],[5,0,1],[0,5,1],[-5,0,1]]),5);
// cir=cir(5);
// echo(len(sec), len(cir));
//
// sec1=sort_points(cir,sec);
// echo(len(sec1),len(cir));
//
// points(sec1,.2);
// points(cir,.2);

function sort_points(sec,list)=[if(list!=[])let(
a=[for(p=sec)min([for(i=[0:len(list)-1])norm(list[i]-p)])],
b=[for(p=sec)[for(i=[0:len(list)-1])norm(list[i]-p)]],
c=[for(i=[0:len(sec)-1])each search(a[i],b[i])],
d=[for(i=c)list[i]]
)d][0];

// experimental and not yet mature

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

// module to create a hollow prism with defined 2 prisms. first prism is the outer one and the second one is the inner one. number of points of the outer and inner prism should be exactly same.
// example:
// sec=cir(10);
// prism=l_extrude(sec,15);
// prism1=l_extrude(f_offset(sec,-1),15);
//
// swp_prism_h(prism,prism1);

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

// function is used as input to another function

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
cp=2p_arc_cp(pa,pb,r,1),
pc=p1+u1*r*rm(90),
pd=p1+u2*r*rm(-90)
) cw([p0,p1,p2])==-1?2p_arc(pc, pd,r,-1,s=norm(pc-pd)<1?0:5):[cp]],

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
op03=[for(p=op02) if(min([for(p1=m_points_sc(sec,10,.2))norm(p-p1)])>=abs(d)-.1)p]
)sort_points(sec, remove_extra_points(op03));

// function is used as input for another function


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
cp=2p_arc_cp(pa,pb,r,-1),
pc=p1+u1*r*rm(-90),
pd=p1+u2*r*rm(90)
) cw( [p0, p1, p2])==-1?[cp]:2p_arc(pc, pd,r,1,s=norm(pc-pd)<1?0:5)],

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

op03=[for(p=op02) if(min([for(p1=m_points_sc(sec,10,.2))norm(p-p1)])>=abs(d)-.1)p]

) sort_points (sec, remove_extra_points (op03));

// function for creating offset of a defined section
// example:
// sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5);
// p_line(sec,.2);
//
// sec1=f_offset(sec,-1);
// p_line(sec1,.2);

function f_offset(sec,d)=d<=0?inner_offset(sec,d):outer_offset(sec,d);

// function for offset for very simple use cases where there are no corner radiuses
// example:
//sec=pts([[-15,0],[0,15],[30,0],[0,-15],[5,0],[0,20],[-40,0],[0,-20]]);
//p_line(sec,.2);
//p_line(offst(sec,2),.2);

function offst(sec,r)=let(
    rev_r=r<0?(min_r(sec)>abs(r)?r:-(min_r(sec)-.1)):r,

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

// function to extrude a section along a path by varying section defined by offset "o"
// example:
// sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5);
//
// path=c2t3(arc(20,0,120,s=20));
//
// p_line3d(path,.2);
//
// prism=v_sec_extrude(sec,path,-2);
//
// swp(prism);


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
rev_sec=q_rot(["x90","z-90",str("y",-a),str("z",theta)],sec)
)each i<len(path)-2?[trns(p0,rev_sec)]:[trns(p0,rev_sec),trns(p1,rev_sec)]];

// module to extrude a section along a closed loop path. 2d section "sec" and a 3d path "path" are the 2 arguments to be filled.
// example
// sec=cr([[0,0,.5],[10,0,2],[7,15,1]],10);
//
// path=c2t3(arc(20,0,355,s=72));
//
// p_line3d(path,.2);
//
// p_extrudec(sec,path);

module p_extrudec(sec,path) 
swp_c(p_extrudec(sec,path));

// module to extrude a section along a open loop path. 2d section "sec" and a 3d path "path" are the 2 arguments to be filled.
// example
// sec=cr([[0,0,.5],[10,0,2],[7,15,1]],10);
//
// path=c2t3(arc(20,0,355,s=72));
//
// p_line3d(path,.2);
//
// p_extrude(sec,path);


module p_extrude(sec,path) swp(p_extrude(sec,path));

// function to extrude a section along a closed loop path. 2d section "sec" and a 3d path "path" are the 2 arguments to be filled.
// example
// sec=cr([[0,0,.5],[10,0,2],[7,15,1]],10);
//
// path=c2t3(arc(20,0,355,s=72));
//
// p_line3d(path,.2);
// 
// prism=p_extrudec(sec,path);
//
// swp_c(prism);

function p_extrudec(sec,path)= let(
prism=[for(i=[0:len(path)-1])let(
p0=path[i],
p1=i<len(path)-1?path[i+1]:path[0],
v=p1-p0,
u=uv([v.x,v.y]),
theta=u.y<0?360-acos([1,0]*u):acos([1,0]*u),
prism=trns(p0,q_rot(["x90","z-90",str("z",theta)],sec))
)prism],
prism1=[for(i=[0:len(path)-1])each i<len(path)-1?[prism[i]]:[prism[i],prism[0]]]
)prism1;

// function to extrude a section along a open loop path. 2d section "sec" and a 3d path "path" are the 2 arguments to be filled.
// example
// sec=cr([[0,0,.5],[10,0,2],[7,15,1]],10);
//
// path=c2t3(arc(20,0,355,s=72));
//
// p_line3d(path,.2);
//
// prism=p_extrude(sec,path);
//
// swp(prism);

function p_extrude(sec,path)= [for(i=[0:len(path)-2])let(
p0=path[i],
p1=path[i+1],
v=p1-p0,
u=uv([v.x,v.y]),
theta=u.y<0?360-acos([1,0]*u):acos([1,0]*u),
prism=i<len(path)-2?[trns(p0,q_rot(["x90","z-90",str("z",theta)],sec))]:[trns(p0,q_rot(["x90","z-90",str("z",theta)],sec)),trns(p1,q_rot(["x90","z-90",str("z",theta)],sec))]

)each prism];

// experimental

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
rev_sec=q_rot(["x90","z-90",str("y",-a),str("z",theta)],sec)
) trns(p0,rev_sec)];

// function to create a fillet with 3 known points with radius "r" and number of segments "s"
// example
// p0=[2,3,5];
// p1=[3,7,2];
// p2=[5,8,3];
// 
// r=2;
// s=10;
// 
// fillet=3p_3d_fillet(p0,p1,p2,r,s);
// $fn=20;
// p_line3dc(fillet,.1);

function 3p_3d_fillet(p0,p1,p2,r=1, s=5)=
let(
n=nv([p0,p1,p2]),
theta=(180-acos(uv(p0-p1)*uv(p2-p1)))/2,
alpha=acos(uv(p0-p1)*uv(p2-p1)),
l=r*tan(theta),
cp=p1+q(n,uv(p0-p1)*r/cos(theta),alpha/2),
pa=p1+uv(p0-p1)*l,
arc=[for(i=[0:theta*2/s:theta*2])cp+q(n,pa-cp,-i)],
a=arc[0],b=loop(arc,1,s-1),c=arc[s]
)[p1,each arc];//cr3d([p1,[a.x,a.y,a.z,.01],each b,[c.x,c.y,c.z,.01]],5);

// function to calculate center point for a fillet

function 3p_3d_fillet_cp(p0,p1,p2,r=1, s=5)=
let(
n=nv([p0,p1,p2]),
theta=(180-acos(uv(p0-p1)*uv(p2-p1)))/2,
alpha=acos(uv(p0-p1)*uv(p2-p1)),
l=r*tan(theta),
cp=p1+q(n,uv(p0-p1)*r/cos(theta),alpha/2),
//pa=p1+uv(p0-p1)*l,
//arc=[for(i=[0:theta*2/s:theta*2])cp+q(n,pa-cp,-i)],
//a=arc[0],b=loop(arc,1,s-1),c=arc[s]
)cp;
// function to create a fillet with 3 known points with radius "r" and number of segments "s". point p1 is omitted while drawing the arc
// example
// p0=[2,3,5];
// p1=[3,7,2];
// p2=[5,8,3];
// 
// r=2;
// s=10;
// 
// fillet=3p_3d_fillet_wo_pivot(p0,p1,p2,r,s);
// $fn=20;
// p_line3d(fillet,.1);

function 3p_3d_fillet_wo_pivot(p0,p1,p2,r=1, s=5)=
let(
n=nv([p0,p1,p2]),
theta=(180-acos(uv(p0-p1)*uv(p2-p1)))/2,
alpha=acos(uv(p0-p1)*uv(p2-p1)),
l=r*tan(theta),
cp=p1+q(n,uv(p0-p1)*r/cos(theta),alpha/2),
pa=p1+uv(p0-p1)*l,
arc=[for(i=[0:theta*2/s:theta*2])cp+q(n,pa-cp,-i)]
)arc;

// function for creating 3d arc with 3 known points.
// example:
// p0=[2,3,5];
// p1=[3,7,2];
// p2=[5,8,3];
// points([p0,p1,p2],.3);
// arc=3p_3d_arc([p0,p1,p2],s=20);
// $fn=20;
// p_line3d(arc,.1);

function 3p_3d_arc(points, s=5)=
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

// function to find the radius with 3 known points in 3d space.
// example:
// p0=[2,3,5];
// p1=[3,7,2];
// p2=[5,8,3];
// echo(3p_3d_r([p0,p1,p2])); //=> ECHO: 1.89252

function 3p_3d_r(points)=
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

// function to draw a 3d arc on a plane defined by a normal vector "n" with radius "r" from angle "theta1" to "theta2". Rotation of the arc can be defined as clockwise (cw=1) or counter clockwise (cw=-1). Number of segments of the arc can be defined with "s".
// Example:
// nv=[3,7,5];
// arc=3d_arc(v=nv,r=10,theta1=0,theta2=180,cw=-1,s=50);
// p_line3d(arc,.2);
// p_line3d([o(),nv],.2);

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

// function to convert a 2d section to 3d
// example:
// sec=cr([[0,0,.5],[10,0,2],[7,15,1]],10);
// path=c2t3(arc(20,0,355,s=72));
//
// p_line3d(path,.2);
//
// prism=p_extrude(sec,path);
//
// swp(prism);

function c2t3(sec)=is_undef(sec.x.x)?[sec.x,sec.y,0]:trns([0,0,0],sec);

// Boolean function which returns "true" ot "false" if the value of a variable "t" is between "s" and "e".
// example:
// t=.5;
// echo(lim(t,0,1)); // => true
// echo(lim(t,10,20)); // => false

function lim(t,s=0,e=1)=t>=s&&t<=e;
 
// function to transpose a 3 x 3 matrix
// example:
// v1=[2,3,5];
// v2=[7,8,9];
// v3=[10,11,12];
// echo(t([v1,v2,v3])); // => ECHO: [[2, 7, 10], [3, 8, 11], [5, 9, 12]]
 
function t(m)=[[m.x.x,m.y.x,m.z.x],[m.x.y,m.y.y,m.z.y],[m.x.z,m.y.z,m.z.z]];

// function to select in between points of a section
// example:
// sec=arc(10,0,70,s=10);
// %points(sec,.5);
// points(loop(sec,1,9),.3);

function loop(sec,a,b)=[for(i=[a:b])sec[i]];

// function to create 3d path
// example:
// path=cr3d(pts2([[0,0,0],[5,3,2,1],[3,3,8,2],[-7,4,1]]),10);
// p_line3d(path,.2);

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
i_minus=i==0?len(p)-1:i-1,i_plus=i<len(p)-1?i+1:0)3d_3p_fillet(p[i_minus],p[i],p[i_plus],r[i],s=s)]
)arcs;

// function to calculate a unit vector for a given vector
// example: 
// echo(uv([2,3,4])); // => ECHO: [0.371391, 0.557086, 0.742781]

function uv(v)=v/norm(v);

// function to find sum of a list of numbers
// example:
// echo(sum([1,3,2,5,7])); //=> echo: 18


function sum(list)=[for(i=[0:len(list)-1])1]*list;

// function to find cumsum of a list of numbers
// example:
// echo(cumsum([1,3,2,5,7])); //=> echo: [1, 4, 6, 11, 18]

function cumsum(list)=[for(i=[0:len(list)-1])sum([for(j=[0:i])list[j]])];

// experimental

function add_p3(p,p1=[0,0,0,0,list],n,i=0)= n==0?p1:add_p3(p,[p[i][0]+p1[0],p[i][1]+p1[1],p[i][2]+p1[2],p[i][3],p[i][4]],n-1,i+1);

// experimental

function pts3(p)=[for(n=[1:len(p)])add_p3(p=p,p1=[0,0,0,0],n=n,i=0)];

    
// function to calculate the bounding box dimensions of a prism
// example:
// echo(bb(rsz3d(spr(4),[5,5,8]))); // => ECHO: [5, 5, 8]

function bb(prism)=
let(
p=[for(p=prism)each[for(p1=p) p1]],
bb=[max(p*[1,0,0])-min(p*[1,0,0]),max(p*[0,1,0])-min(p*[0,1,0]),max(p*[0,0,1])-min(p*[0,0,1])])bb;

// function is used as input to another function

function align_xy(sec,nv)=
let(
v1=[nv.x,nv.y],
u1=uv(v1),
theta=u1.y<0?360-acos([1,0]*u1):acos([1,0]*u1),
v2=q([0,0,1],nv,-theta),u2=uv(v2),
theta1=u2.x>0?360-acos([0,0,1]*u2):acos([0,0,1]*u2),
v3=q([0,1,0],v2,theta1),
aligned_sec=c3t2(q_rot([str("z",-theta),str("y",theta1)],sec))
)aligned_sec;

// function is used as input to another function

function 3d_offset_input(sec,nv,o)=
let(

v1=[nv.x,nv.y],
u1=uv(v1),
theta=u1.y<0?360-acos([1,0]*u1):acos([1,0]*u1),
//nv1=v1==[0,0]?nv+[.001,.001,0]:nv,
v2=q([0,0,1],nv,-theta),u2=uv(v2),
theta1=u2.x>0?360-acos([0,0,1]*u2):acos([0,0,1]*u2),
v3=q([0,1,0],v2,theta1),
z_value=q_rot([str("z",-theta),str("y",theta1)],sec)[0].z,
a_sec=f_offset(c3t2(q_rot([str("z",-theta),str("y",theta1)],sec)),o),
rev_align=trns([0,0,z_value],a_sec),
back=q_rot([str("y",-theta1),str("z",theta)],rev_align)

)back;

// function to calculate an offset of sectio  in 3d space
// example:
// sec=trns([7,8,20],align_v([2,3,5],cir(5)));
// p_line3dc(sec,.2);
// p_line3dc(3d_offset(sec,nv(sec),1),.2);

function 3d_offset(sec,nv,o=1)=3d_offset_input(sec,nv+[0.00001,.00001,0],o);

// function to calculate the normal vector of a known section
// example:
// sec=trns([7,8,20],align_v([2,3,5],cir(5)));
// echo(nv(sec)); // =>ECHO: [-0.0160329, -0.0240482, -0.0400802]

function nv(sec)= cross(sec[0]-sec[1],sec[2]-sec[1]);

/*function rnd(sec1)=cr([for(i=[0:len(sec1)-1])
let(
i_minus=i==0?len(sec1)-1:i-1,i_plus=i<len(sec1)-1?i+1:0,
p0=sec1[i_minus],p1=sec1[i],p2=sec1[i_plus],
v1=p0-p1,v2=p2-p1,
u1=uv(v1),u2=uv(v2),
theta=acos(u1*u2)
)[p1.x,p1.y,.1]],5);*/

// function to offset a given 2d path
// example:
// path=cr(pts1([[2,0],[-2,0,2],[-1,10,2],[-2,0]]),5);
// p_lineo(path,.2);
// p_lineo(path_offset(path,1),.2);


function path_offset(path,d)=[for(i=[0:len(path)-2])let(p0=path[i],p1=path[i+1],line=[p0,p1],rev_point=offst_l(line,d))each i<len(path)-2?[rev_point[0]]:rev_point];

function fillet(p1,p2,p3,r)=[for(i=[0:len(p1)-1])each i<len(p1)-1?[3p_3d_fillet(p3[i],p1[i],p2[i],r)]:[3p_3d_fillet(p3[i],p1[i],p2[i],r),3p_3d_fillet(p3[0],p1[0],p2[0],r)]];

// function to create a fillet with 3 known points
// example:
// p0=[2,3,5];
// p1=[3,7,2];
// p2=[5,8,3];
//
// arc=3d_3p_fillet(p0,p1,p2,r=2,s=10);
// p_line3d(arc,.1);
// points([p0,p1,p2],.2);

function 3d_3p_fillet(p0,p1,p2,r,s=5)=
let(
n=nv([p0,p1,p2]),
theta=(180-acos(uv(p0-p1)*uv(p2-p1)))/2,
alpha=acos(uv(p0-p1)*uv(p2-p1)),
l=r*tan(theta),
cp=assert(l<=norm(p0-p1)&&l<=norm(p2-p1),str("radius :",r," is too big"))p1+q(n,uv(p0-p1)*r/cos(theta),alpha/2),
pa=p1+uv(p0-p1)*l,
arc=[for(i=[0:theta*2/s:theta*2])cp+q(n,pa-cp,-i)]
)arc;

function cir_v(r,cp,v)=
let(r=3,
u=uv(v),
p0=[0,0],
p1=p0+u*100000,
p=perp([p0,p1],cp),
a=norm(p-cp)<.1?(u.y<0?360-acos([1,0]*u)+90:acos([1,0]*u)+90):(uv(p-cp).y<0?360-acos([1,0]*uv(p-cp)):acos([1,0]*uv(p-cp))),

tp=cp+[r,0]*rm(a),
p2=tp+v,
tp1=cp+[r,0]*rm(a+180),
p3=tp1+v

)[[tp,p2],[tp1,p3]];

// function to create arc with 2 points and center. parameter "s" is to define number of segments in the arc
// example
// p0=[2,3,5];
// p1=[7,8,9];
// cp=(p0+p1)/2+[0.001,0,0];
//
// points([p0,p1,cp],.3);
//
// arc=2pnc_arc(p0,p1,cp,20);
// p_line3d(arc,.2);

function 2pnc_arc(p0,p1,cp,cw=-1,s=20)=let(
n=uv(nv(len(p0)==2?c2t3([p0,cp,p1]):[p0,cp,p1])),
theta=acos(uv(p0-cp)*uv(p1-cp)),
r1=norm(p0-cp),r2=norm(p1-cp),
arc=assert(abs(norm((p0-cp))-norm((p1-cp)))<.1,str("radiuses ",r1," and ",r2," are unequal"))cw==1?[for(i=[0:theta/s:theta])cp+q(n,p0-cp,i)]:[for(i=[0:(360-theta)/s:(360-theta)])cp+q(n,p0-cp,-i)]
)arc;

// function used as input to function c_hull

function n_pnt(list,s_pnt,a=0)=let(
a1=a==0||a==360?0:a,
r_list=reduced_list(list,[s_pnt]),
n_pnt=[for(p=r_list)let(v=uv(p-s_pnt),ang=v.y<0?360-acos(v*[1,0]):acos(v*[1,0]),ang1=ang==360||ang<.001?0:ang)round(ang1*1000)/1000],
n1=[for(i=[0:len(n_pnt)-1])if(n_pnt[i]>=a1)round(n_pnt[i]*1000)/1000],
n2=search(min(n1),n_pnt,0),
n3=[for(i=n2)norm(s_pnt-r_list[i])],
n4=search(min(n3),n3,0)[0],
point=r_list[n2[n4]],
//n2=search(n11[0],n_pnt,0)[0],
//point=r_list[n2],
v=point-s_pnt,
ang=round(ang(v.x,v.y)*1000)/1000,
ang1=ang==360||ang<.001?0:ang

)[point,ang1];

// function used as input to function c_hull

function c_hull1(list,s_pnt,n_pnt,revised_list)= 
n_pnt.x==s_pnt?revised_list:c_hull1(
list,s_pnt,
n_pnt=n_pnt(list,n_pnt.x,n_pnt.y-.1),
revised_list=concat(revised_list,[n_pnt.x])
);

// function to create a convex hull of a group of points
// example:
// a=rands(0,10,30);
// b=rands(0,7,30);
// pnts=[for(i=[0:len(a)-1])[a[i],b[i]]];
// points(pnts,.3);
// c_hull=c_hull(pnts);
// color("green")
// p_line(c_hull,.2);

function c_hull(list)=
c_hull1(list=list,
s_pnt=s_pnt(list),
n_pnt=n_pnt(list,s_pnt(list)),
revised_list=[s_pnt(list),n_pnt(list,s_pnt(list)).x]);

function f_surf(list,list1)=let(
index=[for(p=list1)each each search([p],c3t2(list),0)]
)[for(i=index)list[i]];

module partial_surf(surf,t){
     
     surf1=trns([0,0,t],surf);
     for(i=[0:len(surf)-2])
         for(j=[0:len(surf[i])-1])
         let(j_plus=j<len(surf[i])-1?j+1:0){
           if(t>0)
             swp([[surf1[i][j],surf1[i+1][j],surf1[i+1][j_plus],surf1[i][j_plus]],[surf[i][j],surf[i+1][j],surf[i+1][j_plus],surf[i][j_plus]]]);
         else
             swp([[surf[i][j],surf[i+1][j],surf[i+1][j_plus],surf[i][j_plus]],[surf1[i][j],surf1[i+1][j],surf1[i+1][j_plus],surf1[i][j_plus]]]);}}
             
function resurf1(list,c_hull,revised_list)=
    len(list)<=2?revised_list:
    resurf1(list=reduced_list(list,c_hull),
    c_hull=c_hull(list),
    revised_list=concat(revised_list,[c_hull]));
    

// function to reorganise a set of random points
// example:
// sketch=cr(pts1([[-25,0],[25,20,100],[25,-20]]),20);
// path=cytz(cr(pts1([[0,-5],[50,30,50],[20,-25]]),20));
// surf=surf_extrude(sketch,path);
//
// sec=cr(pts1([[10,-20,20],[60,0,20],[0,40,20],[-60,0,20]]),30);
//
// p_surf=[for(p=surf)each [for(p1=p)[p1.x,p1.y]]];
// p_pnts=pies(p_surf,sec);
//
// //points(p_surf,.3);
//
// //%p_line(sec,.2);
// color("green")
// points(p_pnts,.5);
//
// resurf=resurf(p_pnts);
// for(p=resurf)p_line(p,.2);

function resurf(list)=resurf1(list=list,c_hull=c_hull(list),revised_list=[c_hull(list)]);

//function: intersection between section and point
function ibsap(sec,pnt)=let( 
ip=[for(i=[0:len(sec)-1])let(ep=[0,.00001],
i_plus=i<len(sec)-1?i+1:0,
p0=sec[i],p1=sec[i_plus],
p2=pnt,p3=p2+[1,0],
v1=p1-p0+ep,v2=p3-p2-ep,u1=uv(v1),u2=uv(v2),
//p0+v1*t1=p2+v2*t2
//v1*t1-v2*t2=p2-p0
t1=(i_m2d(t([v1,-v2]))*(p2-p0))[0],
ip=p0+v1*t1,
v3=ip-p2,u3=uv(v3)
)if((lim(t1,0,1)&&sign(u2.x)==sign(u3.x)))ip]
)ip;

//function: points inside enclosed section
// example:
// sketch=cr(pts1([[-25,0],[25,20,100],[25,-20]]),20);
// path=cytz(cr(pts1([[0,-5],[50,30,50],[20,-25]]),20));
// surf=surf_extrude(sketch,path);
//
// sec=cr(pts1([[10,-20,20],[60,0,20],[0,40,20],[-60,0,20]]),30);
//
// p_surf=[for(p=surf)each [for(p1=p)[p1.x,p1.y]]];
// p_pnts=pies(p_surf,sec);
//
// points(p_surf,.3);
//
// p_line(sec,.2);
// color("green")
// points(p_pnts,.5);

function pies(pnts,sec)=let(
pwir=[for(p=pnts)let(
ip=ibsap(sec,p)

)if(ip!=[]&&len(ip)%2==1)p]

)pwir;

function flat(dia=10,cp=[0,0,0])=trns(cp,[cir(.001),cir(dia/2)]);

// function to draw a helix with diameter "dia", pitch "pitch" and number of turns "turns"
// example:
// helix=helix(dia=20,pitch=5,turns=7);
// p_line3d(helix,.2);

function helix(dia=10,pitch=3,turns=5)=[for(i=[0:5:360*turns])[dia/2*cos(i),dia/2*sin(i),i/360*pitch]];

// function to define a plane with normal vector "nv" and diameter of the surface "dia"
// example:
// plane= plane(nv=[2,3,5],dia=20);
// swp(plane);
//
// example 2:
// prism=l_extrude(cir(5,s=50),50);
// p1=ipe(trns([0,0,0],plane([0,0,1],50)),prism,1);
// p2=ipe(trns([0,0,50],plane([0,0,1],50)),flip(prism),1,1);
// swp([each p1,each flip(p2)]);


function plane(nv, dia)=let(
sec1=3d_arc(nv,.01,0,360,-1),
sec2=3d_arc(nv,dia/2,0,360,-1),
plane=[sec1,sec2]
)plane;

// function to define origin
// example:
// v=[2,3,5];
// p_line3d([o(),v],.2,$fn=20);

function o()=[0,0,0];

// function to align any shape "prism" with a vector "v"
// example:
// v=[20,30,50];
// prism=l_extrude(cir(1),50);
// aligned_prism=align_v(v,prism);
// %swp(aligned_prism);
// p_line3d([o(),v],.2);

function align_v(v,prism)=let(
v=v+[.0001,0,0],
theta1=ang(v.x,v.y),
theta2=ang(norm([v.x,v.y]),v.z),
avg=avg_v(prism),
t=[avg.x,avg.y],
rev_prism=trns(t-o(),q_rot([str("y",90-theta2),str("z",theta1)],trns(o()-t,prism)))
)rev_prism;

// module to create a solid with base on x-y plane for a surface, produced with function surf_extrude(). Parameter "h" gives distance of base from x-y plane e.g. -ve value of "h" meansthe base is below the x-y plane and +ve value means it is above the x-y plane.
// example:
// sketch=cr(pts1([[-25,0],[25,20,100],[25,-20]]),20);
// path=cytz(cr(pts1([[0,-5],[50,30,50],[20,-25]]),20));
// surf=surf_extrude(sketch,path);
// surf_base(surf,h=-10);

 module surf_base(surf,h=0){
     surf1=trns([0,0,h],[for(p=surf)[for(p1=p)[p1.x,p1.y]]]);

         for(i=[0:len(surf)-2])
         for(j=[0:len(surf[i])-2])
           if(surf1.x.x.z>surf.x.x.z)
             swp([[surf1[i][j],surf1[i+1][j],surf1[i+1][j+1],surf1[i][j+1]],[surf[i][j],surf[i+1][j],surf[i+1][j+1],surf[i][j+1]]]);
         else
             swp([[surf[i][j],surf[i+1][j],surf[i+1][j+1],surf[i][j+1]],[surf1[i][j],surf1[i+1][j],surf1[i+1][j+1],surf1[i][j+1]]]);}
   
function remove_duplicate(path)=[for(i=[0:len(path)-1])let(
p0=path[i],p1=i<len(path)-1?path[i+1]:path[i]+path[i]*100
)if(norm(p1-p0)>.01)p0];

// function to find the tanget to a circle from a point outside the circle
// example:
// point=[10,0];
// cir=cir(r=7.5,cp=[22.5,15]);
// p_line([point,p_cir_t(point,cir)],.2);
// p_line(cir,.2);

function p_cir_t(pnt,cir)=let(
ang=[for(i=[0:len(cir)-1])let(
i_plus=i<len(cir)-1?i+1:0,
p0=cir[i],p1=cir[i_plus],
v=p1-p0,
a1=ang_v(v),
v1=p0-pnt,
a2=ang_v(v1),
a=360/len(cir)/2)
abs(a1-a2)],

ang1=[for(i=[0:len(cir)-1])let(
i_plus=i<len(cir)-1?i+1:0,
p0=cir[i],p1=cir[i_plus],
v=p1-p0,
a1=ang_v(v),
v1=p0-pnt,
a2=ang_v(v1),
a=360/len(cir)/2)
if(abs(a1-a2)<a)abs(a1-a2)][0],

i=search(ang1,ang,0)[0],

point=cir[i+1]
)point;
// function to find the tanget from a circle to a point outside the circle
// example:
// point=[10,0];
// cir=cir(r=7.5,cp=[22.5,15]);
// p_line([cir_p_t(cir,point),point],.2);
// p_line(cir,.2);

function cir_p_t(cir,pnt)=let(
ang=[for(i=[0:len(cir)-1])let(
i_plus=i<len(cir)-1?i+1:0,
p0=cir[i],p1=cir[i_plus],
v=p1-p0,
a1=ang_v(v),
v1=pnt-p0,
a2=ang_v(v1),
a=360/len(cir)/2)
abs(a1-a2)],

ang1=[for(i=[0:len(cir)-1])let(
i_plus=i<len(cir)-1?i+1:0,
p0=cir[i],p1=cir[i_plus],
v=p1-p0,
a1=ang_v(v),
v1=pnt-p0,
a2=ang_v(v1),
a=360/len(cir)/2)
if(abs(a1-a2)<a)abs(a1-a2)][0],

i=search(ang1,ang,0)[0],

point=cir[i]
)point;

// function to find the angle of a 2d vector with [1,0]
// example
//  point=[10,0];
//  cir=cir(r=7.5,cp=[22.5,15]);
//  tangent_point=p_cir_t(point,cir);
//  v=tangent_point-point;
//  ang=ang_v(v);
//  echo(ang); // ECHO: 27.6865

function ang_v(v)=ang(v.x,v.y);

//function to change the orientation of points of a prism. for example check prism and prism1
//sec=cr(pts1([[5,20],[20,10,50],[20,-7]]),20);
//sec1=trns([0,-.05,0],sec);
//prism=[for(i=[0:5:355])rot([1,0,0],sec,i)];
//prism1=cpo(prism);
//
//translate([-60,0,0])
//for(p=prism)points(p,.5);
//for(p=prism)p_line3d(p,.2);
//translate([60,0,0])
//for(p=prism1)p_line3d(p,.2);
//
//translate([-60,0,30])rotate([90,0,0])text("points");
//translate([0,0,30])rotate([90,0,0])text("prism");
//translate([60,0,30])rotate([90,0,0])text("prism1");

function cpo(prism)=[for(i=[0:len(prism[0])-1])[for(p=prism)p[i]]];

// function to draw an offset for a 3d prism
// example:
// sec=cir(10);
// path=cr(pts1([[2,0],[-2,0,2],[-1,10,2],[-4,0]]),5);
// prism=prism(sec,path);
// //swp(prism);
// prism1= surf_offset(prism,-1);
// swp_prism_h(prism,prism1);

function surf_offset(prism,d)=
[for(i=[0:len(prism)-1])[for(j=[0:len(prism[0])-1])
let(
j_plus=j<len(prism[0])-1?j+1:0,
p0=prism[i][j],
p1=i<len(prism)-1?prism[i][j_plus]:prism[i-1][j],
p2=i<len(prism)-1?prism[i+1][j]:prism[i][j_plus],
v1=p1-p0,v2=p2-p0,
u1=uv(v1),u2=uv(v2),

p3=cross(u1,u2)*d

)p0+p3
]];

// function to round a vector, "v" to "n" decimal points
// example:
// echo(rnd_v(v=[2.3456,3.27598,5.876921],n=3)); // ECHO: [2.346, 3.276, 5.877]

function rnd_v(v,n)=[for(p=v)round(p*10^n)/10^n];

// function to round a value "v" to "n" decimal points
// example:
// echo(rnd(v=7.9816523,n=2)); // ECHO: 7.98

function rnd(v,n)=round(v*10^n)/10^n;

function udef(s)=[if(s!=undef)s else []].x;

// function to round a list of points to "n" decimal points
// example:
// x=rands(0,10,50);
// y=rands(3,10,50);
// p=[for(i=[0:len(x)-1])[x[i],y[i]]];
// echo(rnd_list(p,3));

function rnd_list(list,n)=[for(p=list)rnd_v(p,n)];

// input to offset function

function io(s1,r)=let(
s=cw(s1)==1?flip(s1):s1,
sec=convert_sec(s,abs(r)),

sec1=[for(i=[0:len(sec)-1])let(
i_plus=i<len(sec)-1?i+1:0,
p0=sec[i],
p1=sec[i_plus],
l=[p0,p1]
)offst_l(l,r)],

sec2=[for(i=[0:len(sec1)-1])let(
i_minus=i==0?len(sec)-1:i-1,
i_plus=i<len(sec)-1?i+1:0,
p0=sec[i_minus],
p1=sec[i],
p2=sec[i_plus],
cw=cw([p0,p1,p2]),
p3=i_p2d(sec1[i_minus],sec1[i])
)each if(cw==1)sec1[i]],

//sec3=[for(p=sec1)each p],

//sec4=[for(i=[0:len(sec3)-1])let(
//i_plus=i<len(sec3)-1?i+1:0,
//p0=sec3[i],
//p1=sec3[i_plus],
//l=[p0,p1]
//)l],

sec5=[for(p=sec1)each remove_extra_points([for(p1=sec1)let(
l1=p,l2=p1,
v1=l1.y-l1.x,v2=l2.y-l2.x,
im=i_m2d(t([v1,-v2])),
t=(im*(l2.x-l1.x)).x,
u=(im*(l2.x-l1.x)).y,
)if(p!=p1&&lim(t,0,1)&&lim(u,0,1))l1.x+t*v1])],

sec6=[for(i=[0:len(sec)-1])let(
i_plus=i<len(sec)-1?i+1:0,
p0=sec[i],p1=sec[i_plus]
)[p0,p1]],

sec7=remove_extra_points([for(p=[each sec5, each sec2])if(min([for(l=sec6)

let(
v1=l.y-l.x,
v2=p-l.x,
u1=uv(v1),
d=v1*v2/norm(v1),
t=rnd(d/norm(v1),3),
p1=l.x+u1*d
)lim(t,0,1)?rnd(norm(p-p1),3):10^5])==abs(r))p])
)sort_points(s,sec7);

function oo(s,r)=let(
sec=[for(i=[0:len(s)-1])let(
i_minus=i==0?len(s)-1:i-1,
i_plus=i<len(s)-1?i+1:0,
p0=s[i_minus],p1=s[i],p2=s[i_plus],
cw=cw([p0,p1,p2])
)each if(cw==-1)offst_l([p1,p2],r)],

sec1=[for(i=[0:len(s)-1])let(
i_minus=i==0?len(s)-1:i-1,
i_plus=i<len(s)-1?i+1:0,
p0=s[i_minus],p1=s[i],p2=s[i_plus],
cw=cw([p0,p1,p2])
)each offst_l([p1,p2],r)],

sec2=[for(i=[0:len(sec1)-1])let(
i_plus=i<len(sec1)-1?i+1:0,
p0=sec1[i],p1=sec1[i_plus]
)[p0,p1]],

sec3=[for(p=sec2)each [for(p1=sec2)let(
l1=p,l2=p1,
v1=l1.y-l1.x,v2=l2.y-l2.x,
im=i_m2d(t([v1,-v2])),
t=(im*(l2.x-l1.x)).x,
u=(im*(l2.x-l1.x)).y,
)if(p!=p1&&lim(t,0,1)&&lim(u,0,1))l1.x+t*v1]],

sec4=reduced_list(sec3,sec),

sec5=[each sec4,each sec],

sec6=[for(i=[0:len(s)-1])let(
i_plus=i<len(s)-1?i+1:0,
p0=s[i],p1=s[i_plus]
)[p0,p1]],

sec7=remove_extra_points([for(p=sec5)if(min([for(l=sec6)let(
v1=l.y-l.x,
v2=p-l.x,
u1=uv(v1),
d=v1*v2/norm(v1),
t=rnd(d/norm(v1),3),
p1=l.x+u1*d
)lim(t,0,1)?rnd(norm(p-p1),3):10^5])==abs(r))p])
)sort_points(s,sec7);


// function for drawing a offset to a section. This is a finer quality and takes longer than the f_offset function
// example
// sec=cr(pts1([[0,0,.5],[7,5,2],[5,7,3],[-5,7,5],[-7,5,5]]),10);
//
// path=cr(pts1([[2,0],[-2,0,2],[0,7,5],[-5,0]]),20);
// prism=prism1(sec,path);
//  
// swp(prism);


function offset(s,r)=r==0?s:r<0?(sec_r(s)>=abs(r)?f_offset(s,r):io(s,r)):outer_offset(s,r);

function sec_d(s,p,r)=[if(min([for(i=[0:len(s)-1])let(
i_minus=i==0?len(s)-1:i-1,
i_plus=i<len(s)-1?i+1:0,
p0=s[i_minus],
p1=s[i],
p2=s[i_plus],
l1=[p1,p2],
p3=perp(l1,p),
d=rnd(norm(p3-p),3),
d1=rnd(norm(p3-p1),3),d2=rnd(norm(p2-p1),3),
u1=rnd_v(uv(p3-p1),3),u2=rnd_v(uv(p2-p1),3),
)if(d1<=d2&&u1==u2)d else 10^5])==abs(r))p];

function sec_d_min(s,p,r)=[min([for(i=[0:len(s)-1])let(
i_minus=i==0?len(s)-1:i-1,
i_plus=i<len(s)-1?i+1:0,
p0=s[i_minus],
p1=s[i],
p2=s[i_plus],
l1=[p1,p2],
p3=perp(l1,p),
d=rnd(norm(p3-p),3),
d1=rnd(norm(p3-p1),3),d2=rnd(norm(p2-p1),3),
u1=rnd_v(uv(p3-p1),3),u2=rnd_v(uv(p2-p1),3),
)if(d1<=d2&&u1==u2)d else 10^5])];


// function to calculate center point of a triangle "t".
// example:
// t=[[0,0],[10,0],[5,10]];
// p_line(t,.1);
// cp=triangle_cp(t);
// points([cp],0.3);

function triangle_cp(t)=let(
p=cw(t)==1?flip(t):t,
c1=(t[0]+t[1])/2,c2=(t[1]+t[2])/2,c3=(t[2]+t[0])/2,
u1=uv(t[1]-t[0]),u2=uv(t[2]-t[1]),u3=uv(t[0]-t[2]),
l1=[c1,c1+u1*rm(90)], l2=[c2,c2+u2*rm(90)], l3=[c3,c3+u3*rm(90)],
ip=i_p2d(l1,l2)
)ip;

// function to find bottom left point from a group of points
function bl_pnt(sec)=let(
y=sec*[0,1],
loc=[for(i=[0:len(y)-1])if(abs(min(y)-y[i])<.001)i],
x=[for(i=loc)sec[i]],
x_min=min(x*[1,0]),
i=search(x_min,x,0,0)[0],
pnt=x[i]
)pnt;

// function to find bottom right point from a group of points
function br_pnt(sec)=let(
y=sec*[0,1],
loc=[for(i=[0:len(y)-1])if(abs(min(y)-y[i])<.001)i],
x=[for(i=loc)sec[i]],
x_max=max(x*[1,0]),
i=search(x_max,x,0,0)[0],
pnt=x[i]
)pnt;

// function to find top left point from a group of points

function tl_pnt(sec)=let(
y=sec*[0,1],
loc=[for(i=[0:len(y)-1])if(abs(max(y)-y[i])<.001)i],
x=[for(i=loc)sec[i]],
x_min=min(x*[1,0]),
i=search(x_min,x,0,0)[0],
pnt=x[i]
)pnt;

// function to find top right point from a group of points

function tr_pnt(sec)=let(
y=sec*[0,1],
loc=[for(i=[0:len(y)-1])if(abs(max(y)-y[i])<.001)i],
x=[for(i=loc)sec[i]],
x_max=max(x*[1,0]),
i=search(x_max,x,0,0)[0],
pnt=x[i]
)pnt;

// function to find the bounding box dimensions for group of 2d points
function bb2d(p)=[max(p*[1,0])-min(p*[1,0]),max(p*[0,1])-min(p*[0,1])];

// function to find the left most point from a list of points
function l_m(p)=p[search(min(p*[1,0]),p,0,0).x];

// function to find the right most point from a list of points
function r_m(p)=p[search(max(p*[1,0]),p,0,0).x];

// function to find the top most point from a list of points
function t_m(p)=p[search(max(p*[0,1]),p,0,1).x];

// function to find the bottom most point from a list of points
function b_m(p)=p[search(min(p*[0,1]),p,0,1).x];

// function to sort the list and remove equal numbers.
function unique_sort(list)=let(
a=sort(list),
b=[for(i=[0:len(a)-1])let(
i_plus=i<len(a)-1?i+1:0,
)if(rnd(a[i],3)!=rnd(a[i_plus],3))a[i]]
)b;

// function to sort the points in lexicographic order
// example:
// sec=m_points_so([[-10,0],[10,0]],10);
// path=m_points_so([[0.001,0,0.001],[20,0,0.001]],10);
// p=surf_extrude(sec,path);
// p1=c3t2([for(n=p)each n]);
// sec1=ellipse(10,7,[10,0],s=100);
// sec3=rnd_list(pies(p1,sec1),3);
// p2=lexicographic_sort(sec3);
// for(i=[0:len(p2)-1])translate(p2[i])text(str(i),.5);

function lexicographic_sort(list)=[let(
a=list*[1,0],
b=unique_sort(a),
c=[for(i=[0:len(b)-1])search(rnd(b[i],3),rnd_list(list,3),0,0)],
d=[for(i=c)let(
e=[for(j=i)list[j]]
)e],
f=[for(p=d)let(
g=sort(p*[0,1]),
h=[for(n=g)p[search(rnd(n,3),rnd_list(p,3),0,1).x]]
)each h]

)f ].x;

// function to convert a section with corner radius to without radius for a given radius. for eaxmple 
// sec=cr(pts1([[0,0,.2],[8,3,3],[5,7,1],[-8,0,2],[-5,20,1]]),20);
// //sec=cr(pts1([[0,0,.5],[7,5,2],[5,7,3],[-5,7,5],[-7,5,5]]),20);
// %p_line(sec,.1);
// sec1=convert_sec(sec,4); // in case the corner radius of a section is < 4, the corner radius reduced to 0.
// p_line(sec1,.1);

function convert_sec(sec,d)=let(
r=[for(i=[0:len(sec)-1])let(
i_2minus=i==0?len(sec)-2:i==1?len(sec)-1:i-2,
i_minus=i==0?len(sec)-1:i-1,
i_plus=i<len(sec)-1?i+1:0,
i_2plus=i<len(sec)-2?i+2:i==len(sec)-2?0:1,
pi_2minus=sec[i_2minus],
pi_minus=sec[i_minus],
pi=sec[i],
pi_plus=sec[i_plus],
pi_2plus=sec[i_2plus],
v1=pi_minus-pi_2minus,
v2=pi-pi_minus,
v3=pi_plus-pi,
v4=pi_2plus-pi_plus,
l1=rnd(norm(v1),3),
l2=rnd(norm(v2),3),
l3=rnd(norm(v3),3),
l4=rnd(norm(v4),3),
r1=rnd(3p_r(pi_2minus,pi_minus,pi),3),
r2=rnd(3p_r(pi_minus,pi,pi_plus),3),
r3=rnd(3p_r(pi,pi_plus,pi_2plus),3)
)if(l2!=l3&&(r1!=r2 || r2!=r3))0 else r2],

sec1=[for(i=[0:len(r)-1])let(
i_minus=i==0?len(sec)-1:i-1,
i_plus=i<len(sec)-1?i+1:0,

p0=sec[i_minus],
p1=sec[i],
p2=sec[i_plus],
cw=cw([p0,p1,p2])
)if((r[i]==0||r[i]>=d)||cw==1)sec[i]],

sec2=[for(p=sec1)search([p],sec,0).x.x],

sec3=[for(i=[0:len(sec2)-1])let(
i_minus=i==0?len(sec2)-1:i-1,
i_plus=i<len(sec2)-1?i+1:0,
i_2plus=i<len(sec2)-2?i+2:i==len(sec2)-2?0:1,
)sec2[i_plus]-sec2[i]>1?i_p2d([sec[sec2[i_minus]],sec[sec2[i]]],[sec[sec2[i_plus]],sec[sec2[i_2plus]]]):sec[sec2[i]]],

sec4=[for(i=[0:len(sec3)-1])let(
i_minus=i==0?len(sec3)-1:i-1,
i_plus=i<len(sec3)-1?i+1:0,
v1=sec3[i]-sec3[i_minus],
v2=sec3[i_plus]-sec3[i_minus],
u1=rnd_v(uv(v1),3),
u2=rnd_v(uv(v2),3)
)if(u1!=u2)sec3[i]]
)
sec4;

function top_bottom_sort(list)=[let(
a=list*[0,1],
b=flip(unique_sort(a)),
c=[for(i=[0:len(b)-1])search(rnd(b[i],3),rnd_list(list,3),0,1)],
d=[for(i=c)let(
e=[for(j=i)list[j]]
)e],
f=[for(p=d)let(
g=sort(p*[1,0]),
h=[for(n=g)p[search(rnd(n,3),rnd_list(p,3),0,0).x]]
)each h]

)f ].x;
