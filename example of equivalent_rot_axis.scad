

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


function equivalent_rot_axis(r1=[])=
let(
vz=[0,0,1],
v1=q_rot(r1,[vz])[0],
v2=cross(vz,v1),
theta=acos(vz*v1/norm(v1))
)
[theta,v2];



//original openscad method of rotating the objects

rotate([0,45,0])
rotate([70,0,0])
rotate([0,10,0])
rotate([30,40,100])
cylinder(h=50,$fn=30);


//original openscad method of reverse rotating the objects
rotate([-30,0,0])
rotate([0,-40,0])
rotate([0,0,-100])
rotate([0,-10,0])
rotate([-70,0,0])
rotate([0,-45,0])


rotate([0,45,0])
rotate([70,0,0])
rotate([0,10,0])
rotate([30,40,100])

cylinder(h=50,$fn=30);

// equivalent rotation axis method of rotation
r1=["x30","y40","z100","y10","x70","y45"];
a=equivalent_rot_axis(r1);
color("magenta")
rotate(a.x,a.y)
cylinder(h=50,$fn=30);


// equivalent rotation axis method of reverse rotation

color("blue")
rotate(-a.x,a.y)
rotate(a.x,a.y)
cylinder(h=50,$fn=30);



