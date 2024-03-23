 

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
