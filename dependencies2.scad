 

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
    
    
// module for rendering points along the various shapes 2d or 3d. parameter "d" is the size of cube which is used as point. a list has to be provided for parameter "p"
// try following code:
// sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5);
// prism=l_extrude(sec,h=15,a=90,steps=20);
// %swp(prism);
// for(p=prism) points(p,.2);
 
module points(p,d=.5){
    for(i=p)translate(i)cube(size=d,center=true);
    
    }
    
//module to draw a polyline in 3d space (loop not closed)
// e.g. try following code:
// sec=trns([5,10,6],q_rot(["x45"],circle(10)));
// p_line3d(sec,.2);
    
module p_line3d(path,d,rec=0,$fn=20){
    for(i=[0:len(path)-2])
        
    hull(){
    translate(path[i])if(rec==0)sphere(d/2); else cube(d,true);
    translate(path[i+1])if(rec==0)sphere(d/2);else cube(d,true);
    }}

//module to draw a polyline in 3d space (loop closed)
// e.g. try following code:
// sec=trns([5,10,6],q_rot(["x45"],circle(10)));
// p_line3dc(sec,.2);    

module p_line3dc(path,d,rec=0,$fn=20){
    for(i=[0:len(path)-1])
        let(
    i_plus=i<len(path)-1?i+1:0
    )
    hull(){
    translate(path[i])if(rec==0)sphere(d/2); else cube(d,true);
    translate(path[i_plus])if(rec==0)sphere(d/2);else cube(d,true);
    }}
    
function faces(sol)=

//    calculate the faces for the vertices with shape l x m with first and the last end closed
    let(
    l=len(sol),
    m=len(sol[0]),
    n1=[for(i=[0:m-1])i],
    n2=[for(i=[0:l-2]) each ([ for(j=[0:m-1])
    each
    j<m-1?[[(j+1)+i*m,j+i*m,j+(i+1)*m],[(j+1)+i*m,j+(i+1)*m,(j+1)+(i+1)*m]]:
    [[0+i*m,j+i*m,j+(i+1)*m],[0+i*m,j+(i+1)*m,0+(i+1)*m]]
    ])],
    n3=[for(i=[0:m-1])i+(l-1)*m],
    n4=[for(i=[len(n3)-1:-1:0])n3[i]],
    n=[n1,each (n2),n4]
    )n;
    

function faces_1(sol)=

//    calculate the faces for the vertices with shape l x m with first and the last end open
    let(
    l=len(sol),
    m=len(sol[0]),
    n2=[for(i=[0:l-2])each([ for(j=[0:m-1])
    each
    j<m-1?[[(j+1)+i*m,j+i*m,j+(i+1)*m],[(j+1)+i*m,j+(i+1)*m,(j+1)+(i+1)*m]]:
    [[0+i*m,j+i*m,j+(i+1)*m],[0+i*m,j+(i+1)*m,0+(i+1)*m]]
    ])]
    
    )n2; 
   
 function vertices(sol)=
[each for (p=sol)p];

// module for rendering the polyhedron with ends closed
module swp(sol){
let(
v1=vertices(sol),
f1=faces(sol)
)
polyhedron(v1,f1,convexity=10);

}

// module for rendering polyhedron with ends open (mainly for closed polyhedron)
module swp_c(sol){
let(
v1=vertices(sol),
f1=faces_1(sol)
)
polyhedron(v1,f1,convexity=10);

}
// function to select specific section of points from a list of points
function loop ( list , a , b ) = [ for (i = [ a : b ] ) list[i] ];

function faces_surf(sol)=

//    calculate the faces for the vertices with shape l x m with first and the last end open and considering the ends are not closed. like 2 straight lines surface
    let(
    l=len(sol),
    m=len(sol[0]),
    n2=[for(i=[0:l-2])each [for(j=[0:m-2]) each
[[m*i+j,m*(i+1)+j,m*i+(j+1)],
 [m*i+(j+1),m*(i+1)+j,m*(i+1)+(j+1)]] ]]
    
    )n2;
    
    // module for rendering polyhedron with ends open (mainly for open surfaces like 2 straight lines)
module swp_surf(sol){
let(
v1=vertices(sol),
f1=faces_surf(sol)
)
polyhedron(v1,f1,convexity=10);

}
// function to create slices in solid

function slice_sol(sol_1,n=10)=
cpo([for(p=cpo(sol_1)) m_points_so(p,n)]);


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

// function to calculate a unit vector for a given vector
// example: 
// echo(uv([2,3,4])); // => ECHO: [0.371391, 0.557086, 0.742781]

function uv(v)=v/norm(v);

    // module for rendering polyhedron with a closed polygon
module swp_sec(sec){
let(
v1=sec,
f1=[for(i=[0:len(sec)-1]) i]
)
polyhedron(v1,[f1],convexity=10);

}

// function to convert a list of points in a section to list of line segments

function seg(sec)=[for(i=[0:len(sec)-1])let(
 i_plus=i<len(sec)-1?i+1:0,
 p0=sec[i],
 p1=sec[i_plus],
 l=[p0,p1]
 )l];

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

function tr1(tm,sec)=[for(p=sec)(len(tm)==2?[tm.x,tm.y,0]:tm) + (len(p)==2?[p.x,p.y,0]:p)];
    
function tr2(tm,sec)=[for(i=[0:len(sec)-1])[for(p=sec[i])(len(tm)==2?[tm.x,tm.y,0]:tm) + (len(p)==2?[p.x,p.y,0]:p)]];
    
//function to translate a group of points "sl" by "m" distance defined in [x,y,z].e.g. try following code:
//sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5);
//p_line3dc(trns([2,5,10],sec),.1);
 
function trns(m,sec)=is_num(sec[0][0])?tr1(m,sec):tr2(m,sec);

//function to calculate angle of a 2d vector starting from origin and end point with x and y co-ordinates
// example:
// p1=[3,4];p2=[-3,2];
// v=p2-p1;
// p_lineo([p1,p2],.2);
// ang= ang(v.x,v.y);
// echo(ang);

function ang(x,y)= x>=0&&y>=0?atan(y/x):x<0&&y>=0?180-abs(atan(y/x)):x<0&&y<0?180+abs(atan(y/x)):360-abs(atan(y/x));

function sec2vector(v1,sec)=
let(
    theta_y=v1[2]==0?0:ang((v1[0]^2+v1[1]^2)^.5,v1[2]),
    theta_z=ang(v1[0],v1[1])
    )
    q_rot(["x90","z-90",str("y",-theta_y),str("z",theta_z)],sec);
    
    
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

//module to draw a polyline in 3d space (loop not closed) with an arrow to show the direction of polyline
// e.g. try following code:
// sec=circle(20);
// p_line3da(sec,.1);

module p_line3da(l1,d=0.1,rec=0,$fn=20){
p_line3d(l1,d,rec);
l1=c2t3(l1);
l2=seg(l1);
for(n=[0:len(l1)-2])
if (l2[n][0]!=l2[n][1])
let(
p0=l2[n][0],
p1=l2[n][1],
l=norm(p1-p0),
a=l*.05,
p2=[-a/2,0,0],
p3=[a/2,0,0],
p4=[0,sqrt(3)/2*a,0],
cp1=(p2+p3+p4)/3,
px=[for (i=[p2,p3,p4]) i-cp1],
py=sec2vector(p1-p0,px)
){

hull(){
for (i=py)translate(p0+i+(p1-p0)*.9)if(rec==0)sphere(d/2);else cube(d,center=true);
translate(p1) if(rec==0)sphere(d/2,$fn=20); else cube(d,center=true);
}
}
}

//module to draw a polyline in 3d space (loop closed) with an arrow to show the direction of polyline
// e.g. try following code:
// sec=circle(20);
// p_line3dca(sec,.1);

module p_line3dca(l1,d=0.1,rec=0,$fn=20){
p_line3dc(l1,d,rec);
l1=c2t3(l1);
l2=seg(l1);
for(n=[0:len(l1)-1])
if (l2[n][0]!=l2[n][1])
let(
p0=l2[n][0],
p1=l2[n][1],
l=norm(p1-p0),
a=l*.05,
p2=[-a/2,0,0],
p3=[a/2,0,0],
p4=[0,sqrt(3)/2*a,0],
cp1=(p2+p3+p4)/3,
px=[for (i=[p2,p3,p4]) i-cp1],
py=sec2vector(p1-p0,px)
){

hull(){
for (i=py)translate(p0+i+(p1-p0)*.9)if(rec==0)sphere(d/2);else cube(d,center=true);
translate(p1) if(rec==0)sphere(d/2,$fn=20); else cube(d,center=true);
}
}
}

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