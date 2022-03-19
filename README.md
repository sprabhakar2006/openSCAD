# openSCAD
dependencies.scad has functions which help in making complex 3d objects.
examples:

<img width="200" alt="Screenshot 2022-02-03 at 10 23 30 PM" src="https://user-images.githubusercontent.com/55306937/152389738-6b938653-9dd1-442e-995a-18f48796e8bc.png" title="Sine wave pots"> <img width="200" alt="Screenshot 2022-02-04 at 8 40 00 PM" src="https://user-images.githubusercontent.com/55306937/152553171-694931f9-aeca-4330-bc90-bde2db84b6ee.png" title="filleted cylinders"> <img width="200" alt="Screenshot 2022-02-05 at 12 23 04 PM" src="https://user-images.githubusercontent.com/55306937/152632165-6f31fa7e-d76c-444e-a963-a1524acdf38a.png" title="cylinder with pocket"> <img width="200" alt="Screenshot 2022-02-05 at 1 04 33 PM" src="https://user-images.githubusercontent.com/55306937/152633155-0f888318-3081-42bf-9599-96c444ed89e7.png" title="cylinder with random pocket"> <img width="200" alt="Screenshot 2022-02-05 at 3 33 37 PM" src="https://user-images.githubusercontent.com/55306937/152637301-7f30f85c-3c95-4cba-914a-816577e34135.png" title="sphere with random pocket"> <img width="200" alt="Screenshot 2022-02-08 at 8 56 34 PM" src="https://user-images.githubusercontent.com/55306937/153019249-7dc470d5-2cec-4c73-a010-e7103244dc76.png" title="sphere plane fillet"> <img width="200" alt="Screenshot 2022-02-05 at 6 48 42 PM" src="https://user-images.githubusercontent.com/55306937/152643819-28a4f867-7ba5-482a-aba6-eb56c19f6818.png" title="variable rounding"> <img width="200" alt="Screenshot 2022-02-05 at 8 41 41 PM" src="https://user-images.githubusercontent.com/55306937/152647519-1a3cd910-8bf9-4825-8c14-22090f773745.png" title="hub with post"> <img width="200" alt="Screenshot 2022-02-07 at 8 59 35 PM" src="https://user-images.githubusercontent.com/55306937/152818899-26bcf636-667d-496f-b803-b52a2bd73858.png" title="surface cylinder fillet"> <img width="200" alt="Screenshot 2022-02-12 at 8 49 57 PM" src="https://user-images.githubusercontent.com/55306937/153717282-156a1b18-b681-4a76-aff8-381b974685c5.png" title="rounded c section"> <img width="200" alt="Screenshot 2022-02-12 at 9 24 28 PM" src="https://user-images.githubusercontent.com/55306937/153718365-656046d9-5a63-4d48-9519-82ddba260b90.png" title="swp prism hollow example">
<img width="200" alt="Screenshot 2022-02-17 at 9 38 42 PM" src="https://user-images.githubusercontent.com/55306937/154522816-e4e4e5ba-7ec4-484d-9460-2c1fb1be4a83.png" title="wheel"> <img width="200" alt="Screenshot 2022-02-18 at 7 03 48 AM" src="https://user-images.githubusercontent.com/55306937/154600866-2b630b2b-fe7b-4f56-a4fc-0f0ba1440606.png" title="hub section fillet"> 
<img width="200" alt="Screenshot 2022-02-24 at 11 09 21 PM" src="https://user-images.githubusercontent.com/55306937/155578133-df3f7889-0be7-4b7a-a171-fd73e2199fd3.png" title="m10"> 
<img width="200" alt="Screenshot 2022-03-12 at 7 50 44 PM" src="https://user-images.githubusercontent.com/55306937/158021710-be8622a2-7a38-4c99-a650-dfd8eb9a6340.png" title="surface"> 
<img width="200" alt="Screenshot 2022-03-12 at 8 07 58 PM" src="https://user-images.githubusercontent.com/55306937/158022323-d251ab86-49a2-4f12-9b3f-fe41573e97ac.png" title="surf_base example"> 
<img width="200" alt="Screenshot 2022-03-17 at 8 40 46 PM" src="https://user-images.githubusercontent.com/55306937/158833600-38717028-7c62-46fd-89e3-0297e70da4b1.png" title="convex hull"> 
<img width="200" alt="Screenshot 2022-03-17 at 8 58 39 PM" src="https://user-images.githubusercontent.com/55306937/158834677-47baa02d-1773-4eab-b318-07e77e92def6.png" title="surface method3"> 
<img width="200" alt="Screenshot 2022-03-19 at 11 16 35 AM" src="https://user-images.githubusercontent.com/55306937/159109020-2a161626-95fe-470b-93b7-e44c0ca8959e.png" title="m35">

Explanation of functions in file dependencies.scad:

functions	Brief explanation
2cir_fillet(r1=10,r2=10,c1=[0,0],c2=[20,0],r=10)	"function to create 2d fillet between 2 circles, where r1,r2 and c1,c2 are radiuses and enter points of the 2 circles respectively. r-> fillet radius
//example:
%p_line(cir(5),.2);
%p_line(cir(3,[7,0]),.2);
fillet=2cir_fillet(r1=5,r2=3,c1=[0,0],c2=[7,0],r=1);
p_line(fillet,.2);"
2cir_fillet1(r1,r2,c1,c2,r)	"function to create 2d fillet between 2 circles (creates fillet only one side), where r1,r2 and c1,c2 are radiuses and enter points of the 2 circles respectively. r-> fillet radius
//example:
%p_line(cir(5),.2);
%p_line(cir(3,[7,0]),.2);
fillet=2cir_fillet1(r1=5,r2=3,c1=[0,0],c2=[7,0],r=1);
p_line(fillet,.2);"
2cir_filleto(r1=10,r2=10,c1=[0,0],c2=[20,0],r=10)	"function to draw the fillet radius ""r"" between the 2 circle with radiuses ""r1"" and ""r2"" centered at ""c1"" and ""c2"" respectively.This function gives an additional flexibility for drawing fillet only one side. e.g try following example
fillet=2cir_filleto(r1=10,r2=10,c1=[0,0],c2=[20,0],r=10);
p_lineo(fillet[0],.1);"
2cir_tangent(r1,r2,cp1,cp2)	" function to create tangent between 2 circle where r1, r2 and cp1, cp2 are the radiuses and center points of the 2 circles respectively. 
 example:
  tangent=2cir_tangent(5,3,[0,0],[7,0]);
  %p_line(cir(5),.2);
  %p_line(cir(3,[7,0]),.2);
  p_line(tangent,.2);"
2cir_tarc(r1,r2,cp1,cp2,r)	" function for creating arc which is tangent to 2 circles
 try this code as an example:
 sec=2cir_tarc(10,5,[0,0],[20,5],20);
 p_lineo(sec,.2);
 p_line(cir(10),.2);
 p_line(cir(5,[20,5]),.2);"
2ctp(r1,r2,cp1,cp2)	"function to draw tangent line joining 2 circles with radiuses ""r1"" and ""r2"" with center points ""cp1"" and ""cp2"" respectively. This function draws tangent line only one side
 e.g. try this code below:
 sec=2ctp(r1=10,r2=5,cp1=[0,0],cp2=[15,6]);
 p_line(cir(10),.1);
 p_line(cir(5,[15,6]),.1);
 p_line(sec,.1);"
2ctpf(r1,r2,cp1,cp2)	"function to draw tangent line joining 2 circles with radiuses ""r1"" and ""r2"" with center points ""cp1"" and ""cp2"" respectively. This function draws tangent line on both the sides
 e.g. try this code below:
 sec=2ctpf(r1=10,r2=5,cp1=[0,0],cp2=[15,6]);
 p_line(cir(10),.1);
 p_line(cir(5,[15,6]),.1);
 p_line(sec,.1);"
2cyl_fillet(r1,r2,cp1,cp2,r,path)	"function for creating fillet between 2 cylinders. r1, r2 and cp1,cp2 are the radiuses and center points of 2 cylinders respectively. r -> is the fillet radius. path -> is given for rounding the cylinder edges
// example
 path=[[0,0],[0,10]];
 %swp(cyl(r=5,h=15));
 %swp(cyl(r=3,h=10,cp=[7,0]));
 swp(2cyl_fillet(5,3,[0,0],[7,0],1,path));"
2p_arc_cp(p1,p2,r,cw=1)	" function to calculate the center point for arc where 2 points ""p1"" and ""p2"" and radius ""r"" are known (clockwise and counter clockwise will have different center points
 example:
 pnt=2p_arc_cp(p1=[2,3],p2=[6,5],r=5,cw=-1);
 points([pnt],.5);"
2p_arc(p1,p2,r,cw=1,s=20)	" function creates a shortest 2d arc with 2 points with a radius ""r"" and number of segments ""s"". parameter cw(clockwise=1 and counter clockwise=-1) defines the order of arc
try this example for better understanding:
 sec=2p_arc(p1=[2,3],p2=[6,5],r=2.25,cw=-1,s=20);
 p_lineo(sec,.2);"
2pnc_arc(p0,p1,cp,s)	" function to create arc with 2 points and center. parameter ""s"" is to define number of segments in the arc
 example
 p0=[2,3,5];
 p1=[7,8,9];
 cp=(p0+p1)/2+[0.001,0,0];

 points([p0,p1,cp],.3);

 arc=2pnc_arc(p0,p1,cp,20);
 p_line3d(arc,.2);"
2r(p1,p2,r,cw=1,s=20)	" function creates a longest 2d arc with 2 points with a radius ""r"" and number of segments ""s"". parameter cw(clockwise=1 and counter clockwise=-1) defines the order of arc
try this example for better understanding:
 sec=2r(p1=[2,3],p2=[6,5],r=3,cw=-1,s=20);
 p_lineo(sec,.2);"
2spr_fillet(r1,r2,cp1,cp2,r)	" function for creating fillet between 2 spheres. r1, r2 and cp1,cp2 are the radiuses and center points of 2 spheres. r-> fillet radius.
// example:
 swp(spr(r=5,cp=[0,0,0]));
 swp(spr(r=3,cp=[7,0,0]));

 fillet=2spr_fillet(r1=5,r2=3,cp1=[0,0,0],cp2=[7,0,0],1);
 swp(fillet);"
3d_3p_fillet(p0,p1,p2,r,s=5)	" function to create a fillet with 3 known points
 example:
 p0=[2,3,5];
 p1=[3,7,2];
 p2=[5,8,3];

 arc=3d_3p_fillet(p0,p1,p2,r=2,s=10);
 p_line3d(arc,.1);
 points([p0,p1,p2],.2);"
3d_arc(v, r, theta1=0, theta2=180, cw=-1,s=50)	" function to draw a 3d arc on a plane defined by a normal vector ""n"" with radius ""r"" from angle ""theta1"" to ""theta2"". Rotation of the arc can be defined as clockwise (cw=1) or counter clockwise (cw=-1). Number of segments of the arc can be defined with ""s"".
 Example:
 nv=[3,7,5];
 arc=3d_arc(v=nv,r=10,theta1=0,theta2=180,cw=-1,s=50);
 p_line3d(arc,.2);
 p_line3d([o(),nv],.2);"
3d_offset(sec,nv,o=1)	" function to calculate an offset of sectio  in 3d space
 example:
 sec=trns([7,8,20],align_v([2,3,5],cir(5)));
 p_line3dc(sec,.2);
 p_line3dc(3d_offset(sec,nv(sec),1),.2);"
3p_3d_arc(points, s=5)	" function for creating 3d arc with 3 known points.
 example:
 p0=[2,3,5];
 p1=[3,7,2];
 p2=[5,8,3];
 points([p0,p1,p2],.3);
 arc=3p_3d_arc([p0,p1,p2],s=20);
 $fn=20;
 p_line3d(arc,.1);"
3p_3d_fillet_wo_pivot(p0,p1,p2,r=1, s=5)	" function to create a fillet with 3 known points with radius ""r"" and number of segments ""s"". point p1 is omitted while drawing the arc
 example
 p0=[2,3,5];
 p1=[3,7,2];
 p2=[5,8,3];
 
 r=2;
 s=10;
 
 fillet=3p_3d_fillet_wo_pivot(p0,p1,p2,r,s);
 $fn=20;
 p_line3d(fillet,.1);"
3p_3d_fillet(p0,p1,p2,r=1, s=5)	" function to create a fillet with 3 known points with radius ""r"" and number of segments ""s""
 example
 p0=[2,3,5];
 p1=[3,7,2];
 p2=[5,8,3];
 
 r=2;
 s=10;
 
 fillet=3p_3d_fillet(p0,p1,p2,r,s);
 $fn=20;
 p_line3dc(fillet,.1);"
3p_3d_r(points)	" function to find the radius with 3 known points in 3d space.
 example:
 p0=[2,3,5];
 p1=[3,7,2];
 p2=[5,8,3];
 echo(3p_3d_r([p0,p1,p2])); //=> ECHO: 1.89252"
3p_arc(p1,p2,p3,s=30)	" function to create arc with 3 points in 2d
 example:
 sec=3p_arc([1,2],[3,7],[7,3]);
 p_lineo(sec,.2);
 points([[1,2],[3,7],[7,3]],.5);"
3p_r(p1,p2,p3	" function to find radius of arc with 3 known points in 2d
 example:
 radius=3p_r([1,2],[3,7],[7,3]);
 echo(radius); //=> ECHO: 3.30892"
align_v(v,prism)	" function to align any shape ""prism"" with a vector ""v""
 example:
 v=[20,30,50];
 prism=l_extrude(cir(1),50);
 aligned_prism=align_v(v,prism);
 %swp(aligned_prism);
 p_line3d([o(),v],.2);"
ang_v(v)	" function to find the angle of a 2d vector with [1,0]
 example
  point=[10,0];
  cir=cir(r=7.5,cp=[22.5,15]);
  tangent_point=p_cir_t(point,cir);
  v=tangent_point-point;
  ang=ang_v(v);
  echo(ang); // ECHO: 27.6865"
ang(x,y)	"function to calculate angle of a 2d vector starting from origin and end point with x and y co-ordinates
 example:
 p1=[3,4];p2=[-3,2];
 v=p2-p1;
 p_lineo([p1,p2],.2);
 ang= ang(v.x,v.y);
 echo(ang);"
arc(radius,ang1=0,ang2=355,cp=[0,0],s=20)	function to draw points in circular arc with radius, start angle "ang1" , end angle "ang2", center point of the arc "cp" and number of segments required in the arc "s". e.g. following code will draw an arc of radius 5 from 0 to 90 degrees centered at [0,0] with 20 segments in the arc: p_lineo(arc(radius=5,ang1=0,ang2=90,cp=[0,0],s=20),.1);
avg_v(prism	" function to calculate average of a group of points either 2d or 3d
 example:
 sec=cir(10);
 path=cr(pts1([[2,0],[-2,0,2],[0,10,3],[-3,0]]),5);
 prism=prism(sec,path);
 %swp(prism);
 avg=avg_v(prism);
 echo(avg);
 points([avg],.5);"
bb(prism)	" function to calculate the bounding box dimensions of a prism
 example:
 echo(bb(rsz3d(spr(4),[5,5,8]))); // => ECHO: [5, 5, 8]"
bez(p,s=.1)	"function for calculating bezier curve with control points ""p"" and with number of segments 1/s
 example:
 p=[[0,0],[10,5],[0,15],[12,20]]; 
 b=bez(p,.1); 
 points(b,.5);
 //control points
 color(""green"")
 points(p,.5);"
c_hull(list)	" function to create a convex hull of a group of points
 example:
 a=rands(0,10,30);
 b=rands(0,7,30);
 pnts=[for(i=[0:len(a)-1])[a[i],b[i]]];
 points(pnts,.3);
 c_hull=c_hull(pnts);
 color(""green"")
 p_line(c_hull,.2);"
c2t3(sec)	" function to convert a 2d section to 3d
 example:
 sec=cr([[0,0,.5],[10,0,2],[7,15,1]],10);
 path=c2t3(arc(20,0,355,s=72));

 p_line3d(path,.2);

 prism=p_extrude(sec,path);

 swp(prism);"
c3t2(sec)	function to convert 3d to 2d, it just removes the z-coordinate from the points list
cir_p_t(cir,pnt)	" function to find the tanget from a circle to a point outside the circle
 example:
 point=[10,0];
 cir=cir(r=7.5,cp=[22.5,15]);
 p_line([cir_p_t(cir,point),point],.2);
 p_line(cir,.2);"
cir(r,cp=[0,0],s=50)	function for creating points in circle with radius "r", center point "cp" and number of segments "s"
comb(n,i)	math function to calculate number of possible combinations for "n" items with "i" selected items
cr(pl,s=20)	"function to create section with corner radiuses. e.g. following code has 3 points at [0,0],[10,0] and [7,15] and radiuses of 0.5,2 and 1 respectively,s=5 represent the number of segments at each corner radius.
sec=cr(pl=[[0,0,.5],[10,0,2],[7,15,1]],s=5);
p_line(sec,.1);"
cr3d(l,s=5)	" function to create 3d path
 example:
 path=cr3d(pts2([[0,0,0],[5,3,2,1],[3,3,8,2],[-7,4,1]]),10);
 p_line3d(path,.2);"
cub(p,center=false)	" function to draw cube
 swp(cub(p=[10,5,4]));"
cumsum(list)	" function to find cumsum of a list of numbers
 example:
 echo(cumsum([1,3,2,5,7])); //=> echo: [1, 4, 6, 11, 18]"
cw(sec)	function to identify whether the section is clockwise or counter clockwise. cw(sec)==1 means clockwise and -1 means counterclockwise. e.g. echo(cw([[0,0],[4,0],[0,4],[-4,0]]));// -1
cyl(r1=1,r2=1,h=1,cp=[0,0],s=50,r,d,d1,d2,center=false)	draws a cylinder try swp(cyl(r=5,h=15)); 
cytz(path)	function to convert the y co-ordinates to z co-ordinates e.g.[x,y]=>[x,0,y]. 2d to 3d coordinate system
det2d(m)	math function to calculate the determinant of a 2 x 2 matrix
det3d(m)	math function to calculate the determinant of a 3 x 3 matrix
ellipse(r1,r2,cp,s=30)	" function to draw an ellipse with semi-major and semi-minor axis ""r1"" and ""r2"" respectively and with center ""cp"" and number of segment ""s""
 example:
 sec=ellipse(r1=5,r2=3,cp=[2,3],s=30);
 p_line(sec,.2);"
f_offset(sec,d)	" function for creating offset of a defined section
 example:
 sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5);
 p_line(sec,.2);

 sec1=f_offset(sec,-1);
 p_line(sec1,.2);"
fact(n,m=1)	math function to calculate factorial of a number
flip(sec)	function flips the direction of points of 2d section or 3d prism
helix(dia=10,pitch=3,turns=5)	" function to draw a helix with diameter ""dia"", pitch ""pitch"" and number of turns ""turns""
 example:
 helix=helix(dia=20,pitch=5,turns=7);
 p_line3d(helix,.2);"
i_m2d(m)	" math function to calculate the inverse of a 2 x 2 matrix
 example:
 v1=[2,3];
 v2=[3,4];
 echo(i_m2d(t([v1,v2]))); //=> ECHO: [[-4, 3], [3, -2]]"
i_m3d(m)	" math function to calculate the inverse of a 3 x 3 matrix
 example:
 v1=[2,3,4];
 v2=[3,4,1];
 v3=[4,5,6];
 echo(i_m3d(t([v1,v2,v3])));// =>ECHO: [[-2.375, 1.75, 0.125], [-0.25, 0.5, -0.25], [1.625, -1.25, 0.125]]"
i_p2d(l1,l2)	function to calculate the intersection point between 2 lines e.g. echo(i_p2d(l1=[[0,0],[1,4]],l2=[[10,0],[7,2]])); => //ECHO: [1.42857, 5.71429]
i_p3d(l1,l2)	function to calculate intersection point between 2 lines in 3d space (mostly if these lines lie on the same plane)
ip(prism,prism1)	" function to calculate intersection point between two 3d prisms. ""prism"" is the 3d object which is intersected with ""prism1"".
 try below code for better understanding:
 sec=cir(10);
 path=cr(pts1([[2,0],[-2,0,2],[0,10,3],[-9.9,0]]),5);
 prism=prism(sec,path);
 prism1=q_rot([""y40""],cyl(r=3,h=15,s=30));

 %swp(prism);
 %swp(prism1);
 ip=ip(prism,prism1);
 points(ip,.2);"
ipe(prism,prism1,r,option=0,s=5)	
ipf(prism,prism1,r,option=0,s=5)	" function for creating fillet: this function first finds the intersection point between prism and prism1 and then calculates the fillet with radius ""r"". option ""0"" and ""1"" creates fillet either outside or inside.parameter ""s"" is for number of segments in the fillet
 an example below will be more clear (try changing option from 1 =>0 or flip the direction of prism1 by flip(prism1))
 try below code for better understanding:
 sec=cir(10);
 path=cr(pts1([[2,0],[-2,0,2],[0,10,3],[-9.9,0]]),5);
 prism=prism(sec,path);
 prism1=q_rot([""y40""],cyl(r=3,h=15,s=30));

 %swp(prism);
 %swp(prism1);
 fillet=ipf(prism,prism1,r=1,option=1,s=5);
 swp_c(fillet);"
l_cir_fillet(line,r1,r2,cp)	" experimental function
 example:
 sec=l_cir_fillet(line=[[0,0],[0,20]],r1=5,r2=1,cp=[5,10]);
 p_lineo(sec,.2);"
l_cir_ip(line,cir)	" // function to get intersection point between a line and circle
 // example
  line=[[0,0],[3,5]];
  cir=cir(5);
  %p_line(line,.2);
  %p_line(cir,.2);

  pnt=l_cir_ip(line,cir);
  color(""green"")
  points(pnt,.5);
  echo(pnt);"
l_extrude(sec,h=1,a=0,steps=1)	" function for linear extrude a section by height ""h"", also the section can be rotated by an angle ""a"" in number of steps ""steps""
 try following code for better understanding (also try changing ""a"" and ""steps""):
 sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5);
 prism=l_extrude(sec,h=15,a=0,steps=1);
 swp(prism);"
l(l,s=20)	" function to create a line with number of segments ""s""
 example:
 line=l([[0,0],[4,3]],10);
 points(line,.2);  "
lim(t,s=0,e=1)	" Boolean function which returns ""true"" ot ""false"" if the value of a variable ""t"" is between ""s"" and ""e"".
 example:
 t=.5;
 echo(lim(t,0,1)); // => true
 echo(lim(t,10,20)); // => false"
list_ang(sec)	
loop(sec,a,b)	" function to select in between points of a section
 example:
 sec=arc(10,0,70,s=10);
 %points(sec,.5);
 points(loop(sec,1,9),.3);"
m_points_o(sec,sl=20)	" function for placing multiple points on the straight line segments of an open section. parameter ""sl"" is for placing points with pitch distance defined by ""sl""
 example:
 sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5);
 points(sec,.2);
 
 translate([15,0])
 points(m_points_o(sec,2),.2);// segment length=> 2 units"
m_points_sc(sec1,s,m=.5)	" function for calculating multiple points on the straight line segments of a closed section. sec-> closed section; s -> number of segments for each straight line segment of closed section; m-> minimum segment length, if the derived segment length < m, then it is omitted
 sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5);
 points(sec,.2);
 
 translate([15,0])
 points(m_points_sc(sec,s=5,m=.5),.2);// number of segments=> 5"
m_points_so(sec1,s,m=.5)	" function for calculating multiple points on the straight line segments of an open section. sec-> closed section; s -> number of segments for each straight line segment of closed section; m-> minimum segment length, if the derived segment length < m, then it is omitted
 sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5);
 points(sec,.2);
 
 translate([15,0])
 points(m_points_so(sec,s=5,m=.5),.2);// number of segments=> 5"
m_points(sec,sl=20)	" function for placing multiple points on the straight line segments of a closed loop section. parameter ""sl"" is for placing points with pitch distance defined by ""sl""
 example:
 sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5);
 points(sec,.2);
 
 translate([15,0])
 points(m_points(sec,2),.2);// segment length=> 2 units  "
min_r(sec)	" function to get the minimum radius for a defined section
 example:
 sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5);
 echo(min_r(sec)); //=>ECHO: 0.5"
mode_sign(p,plus=0,minus=0,n=0)	
nv(sec)	" function to calculate the normal vector of a known section
 example:
 sec=trns([7,8,20],align_v([2,3,5],cir(5)));
 echo(nv(sec)); // =>ECHO: [-0.0160329, -0.0240482, -0.0400802]"
o()	" function to define origin
 example:
 v=[2,3,5];
 p_line3d([o(),v],.2,$fn=20);"
offst_l(l,d)	" function to offset a line ""l"" by distance ""d"" 
  example
  line=[[0,0],[3,5]];
  %p_line(line,.2);
  p_line(offst_l(line,2),.2);"
p_cir_t(pnt,cir)	" function to find the tanget to a circle from a point outside the circle
 example:
 point=[10,0];
 cir=cir(r=7.5,cp=[22.5,15]);
 p_line([point,p_cir_t(point,cir)],.2);
 p_line(cir,.2);"
p_extrude(sec,path)	" function to extrude a section along a open loop path. 2d section ""sec"" and a 3d path ""path"" are the 2 arguments to be filled.
 example
 sec=cr([[0,0,.5],[10,0,2],[7,15,1]],10);

 path=c2t3(arc(20,0,355,s=72));

 p_line3d(path,.2);

 prism=p_extrude(sec,path);

 swp(prism);"
p_extrudec(sec,path)	" function to extrude a section along a closed loop path. 2d section ""sec"" and a 3d path ""path"" are the 2 arguments to be filled.
 example
 sec=cr([[0,0,.5],[10,0,2],[7,15,1]],10);

 path=c2t3(arc(20,0,355,s=72));

 p_line3d(path,.2);
 
 prism=p_extrudec(sec,path);

 swp_c(prism);"
path_offset(path,d)	" function to offset a given 2d path
 example:
 path=cr(pts1([[2,0],[-2,0,2],[-1,10,2],[-2,0]]),5);
 p_lineo(path,.2);
 p_lineo(path_offset(path,1),.2);"
perp(line,point)	" function to find intersection point at a shortest distance between a point and a line
 example
 line=[[0,0],[3,5]];
 point=[-3,5];

 %p_line(line,.2);
 %points([point],.3);

 p=perp(line,point);
 points([p],.5);
 echo(p);"
pies(pnts,sec)	"function: points inside enclosed section
 example:
 sketch=cr(pts1([[-25,0],[25,20,100],[25,-20]]),20);
 path=cytz(cr(pts1([[0,-5],[50,30,50],[20,-25]]),20));
 surf=surf_extrude(sketch,path);

 sec=cr(pts1([[10,-20,20],[60,0,20],[0,40,20],[-60,0,20]]),30);

 p_surf=[for(p=surf)each [for(p1=p)[p1.x,p1.y]]];
 p_pnts=pies(p_surf,sec);

 points(p_surf,.3);

 p_line(sec,.2);
 color(""green"")
 points(p_pnts,.5);"
plane(nv, dia)	" function to define a plane with normal vector ""nv"" and diameter of the surface ""dia""
 example:
 plane= plane(nv=[2,3,5],dia=20);
 swp(plane);

 example 2:
 prism=l_extrude(cir(5,s=50),50);
 p1=ipe(trns([0,0,0],plane([0,0,1],50)),prism,1);
 p2=ipe(trns([0,0,50],plane([0,0,1],50)),flip(prism),1,1);
 swp([each p1,each flip(p2)]);"
prism(sec,path,m_points=0)	"function to make a prism with combination of 2d section and 2d path
 Example:
 sec=cir(10);
 path=cr(pts1([[2,0],[-2,0,2],[0,10,3],[-3,0]]),5);
 prism=prism(sec,path);
 swp(prism);"
pts(p)	" function is used like a turtle move to create 2d shapes.
 following example will create a rectangle with sides 10 x 5:
 sec=pts([[0,0],[10,0],[0,5],[-10,0]]); // starts at [0,0] then moves 10 units to +x direction then moves 5 units towards +y direction and then moves 10 units to -x direction
 p_line(sec,.1);"
pts1(p)	"same as pts(p) with only difference that it keeps the z value unchanged
 for example:
 sec=pts1([[0,0,1],[10,0,1],[0,5,1],[-10,0,1]]); // starts at [0,0] then moves 10 units to +x direction then moves 5 units towards +y direction and then moves 10 units to -x direction
  echo(sec); // ECHO: [[0, 0, 1], [10, 0, 1], [10, 5, 1], [0, 5, 1]]
  this function is mainly used with function cr(pl,s) (please see the example of function cr(pl,s))"
pts2(p)	same as pts and pts1 and is used as a turtle movement in 3d space and keeps the 4th point same to be used as radius for rounding in function cr3d() (check example fo function cr3d()
q_rot(s,pl)	function to rotate a group of points "pl" around a series of axis with defined angles e.g q_rot(s=["z20","x40","y80"],pl=[[2,0],[10,2]])=> will rotate the line first around z axis by 20 deg then around x axis by 40 degrees and then around y axis by 80 degrees.
q(vector=[1,0,0],point=[0,5,0],theta=0)	function to rotate a point around a vector(axis) with angle theta
reduced_list(list,list1)	" function to subtract points from a list of points
 example:
 list=[[1,2,2],[3,4,5],[10,2,9],[11,1,9]];
 list1=[[1,2,2],[10,2,9]];
 revised_list=reduced_list(list,list1);
 echo(revised_list); //=> ECHO: [[3, 4, 5], [11, 1, 9]]"
remove_duplicate(path)	
resurf(list)	" function to reorganise a set of random points
 example:
 sketch=cr(pts1([[-25,0],[25,20,100],[25,-20]]),20);
 path=cytz(cr(pts1([[0,-5],[50,30,50],[20,-25]]),20));
 surf=surf_extrude(sketch,path);

 sec=cr(pts1([[10,-20,20],[60,0,20],[0,40,20],[-60,0,20]]),30);

 p_surf=[for(p=surf)each [for(p1=p)[p1.x,p1.y]]];
 p_pnts=pies(p_surf,sec);

 //points(p_surf,.3);

 //%p_line(sec,.2);
 color(""green"")
 points(p_pnts,.5);

 resurf=resurf(p_pnts);
 for(p=resurf)p_line(p,.2);"
rm(theta)	"function to rotate a vector by ""theta"" degrees e.g. try following code:
line=[[0,0],[5,3]];
line1=line*rm(30);

p_lineo(line,.1);
p_lineo(line1,.1);"
rot(axis,prism,ang)	" function to rotate an object around any arbitrary axis
 example
  sec=cir(10);
  path=cr(pts1([[2,0],[-2,0,2],[-1,5,3],[-4,0]]),5);
  prism=trns([15,0],prism(sec,path));
  prism1=rot([3,4,7],prism,180);
  swp(prism);
  swp(prism1);
  p_line([[0,0,0],[3,4,7]*10],.2);"
rsz_c(sec,rsz=[1,1,1])	" function to calculate 2d resized section - centered
 example:
 sec=cir(10);
 rsz_sec=rsz_c(sec,[5,3]);
 %p_line(sec,.2);
 p_line(rsz_sec,.2);"
rsz(sec,rsz=[1,1,1])	" function to calculate 2d resized section - placed at minimum y value
 example:
 sec=cir(10);
 rsz_sec=rsz(sec,[5,3]);
 %p_line(sec,.2);
 p_line(rsz_sec,.2);"
rsz3d(prism,rsz=[1,1,1])	" function to calculate the resized prism
example:
 sec=cir(10);
 path=cr(pts1([[2,0],[-2,0,2],[0,10,3],[-3,0]]),5);
 prism=prism(sec,path);
 %swp(prism);
 resized_prism=rsz3d(prism,[5,5,7]);
 swp(resized_prism);"
rsz3dc(prism,rsz=[1,1,1])	" function to calculate the resized prism- centered
example:
 sec=cir(10);
 path=cr(pts1([[2,0],[-2,0,2],[0,10,3],[-3,0]]),5);
 prism=prism(sec,path);
 %swp(prism);
 resized_prism=rsz3dc(prism,[5,5,7]);
 swp(resized_prism);"
s_pnt(sec)	
scl2d_c(sec,sl)	"function to scale a 2d section by an amount ""sl"" which has to be >0 (keeps the revised section in center). e.g.following code scales the section by 0.7 (70% of the original shape)
sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5);
p_line(sec,.1);
p_line(scl2d_c(sec,.7),.1);"
scl2d(sec,sl)	"function to scale a 2d section by an amount ""sl"" which has to be >0 (keeps the y-coordinates same). e.g.following code scales the section by 0.7 (70% of the original shape)
sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5);
p_line(sec,.1);
p_line(scl2d(sec,.7),.1);"
scl3d_c(prism,s=1)	" function to scale a 3d prism keeping the prism centered. takes 2 arguments ""prism"" to scale and the scaling factor ""s"". scale factor can take any real number negative values will scale the prism and turn the prism upside down.
 try the following code to understand better:
 sec=cir(10);
 path=cr(pts1([[2,0],[-2,0,2],[0,10,3],[-3,0]]),5);
 prism=prism(sec,path);
 %swp(prism);
 swp(scl3d_c(prism,.7));"
scl3d(prism,s=1)	" function to scale a 3d prism keeping the base z-coordinate same. takes 2 arguments ""prism"" to scale and the scaling factor ""s"". scale factor can take any real number negative values will scale the prism and turn the prism upside down.
 try the following code to understand better:
 sec=cir(10);
 path=cr(pts1([[2,0],[-2,0,2],[0,10,3],[-3,0]]),5);
 prism=prism(sec,path);
 %swp(prism);
 swp(scl3d(prism,.7));"
sec_r(sec)	" function to get the minimum radius for a defined section
 example:
 sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5);
 echo(sec_r(sec)); //=>ECHO: 0.5"
sort_points(sec,list)	" function to match the number of points between 2 sections
 example:
 sec=cr(pts1([[-2.5,-2.5,1],[5,0,1],[0,5,1],[-5,0,1]]),5);
 cir=cir(5);
 echo(len(sec), len(cir));

 sec1=sort_points(cir,sec);
 echo(len(sec1),len(cir));

 points(sec1,.2);
 points(cir,.2);"
sort(list)	function to sort a list of real numbers in ascending order
spr(r,cp=[0,0,0],s=50)	" function for creating sphere with radius ""r"", center point ""cp"" and number of segments ""s"".
 try following code:
 swp(spr(r=3,cp=[4,5,6],s=30));"
sqr(s,center=false)	" function to draw a rectangle
 e.g. p_line(sqr([10,5]),.1); or polygon(sqr([10,5]));"
sum_v(prism)	" function to calculate the cumulative sum of all the points of a 2d or 3d points list.
 e.g.
 sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5);
 echo(sum_v(sec)); //[95.9558, 82.3332]"
sum(list)	" function to find sum of a list of numbers
 example:
 echo(sum([1,3,2,5,7])); //=> echo: 18"
surf_extrude(sec,path)	function to make surface with a polyline 2d sketch and a 3d path(there is no render here but points can be visualised with following command for(p=surf_extrude(sec,path))points(p,.2); 
swp_c(surf1)	" module for rendering 3d prisms with closed section
 example:
 sec=cir(10);
 path=cr(pts1([[2,0],[-2,0,2],[0,10,3],[-9.9,0]]),5);
 prism=prism(sec,path);
 prism1=q_rot([""y40""],cyl(r=3,h=15,s=30));

 %swp(prism);
 %swp(prism1);
 fillet=ipf(prism,prism1,r=1,option=1,s=5);
 swp_c(fillet);"
swp_prism_h(prism,prism1)	" module to create a hollow prism with defined 2 prisms. first prism is the outer one and the second one is the inner one. number of points of the outer and inner prism should be exactly same.
 example:
 sec=cir(10);
 prism=l_extrude(sec,15);
 prism1=l_extrude(f_offset(sec,-1),15);

 swp_prism_h(prism,prism1);"
swp(surf1)	" module for rendering various 3d prism 
 //example1:
 sec=cir(10);
 path=cr(pts1([[2,0],[-2,0,2],[0,10,3],[-3,0]]),5);
 prism=prism(sec,path);
 swp(prism);
 //example2:
 prism1=l_extrude(sqr([10,6]),15);
 translate([13,0,0])
 swp(prism1);
 //example3:
 sec2=cr(pts1([[0,0,1],[5,0,1],[-2.5,4,1]]),5);
 path2=[for(i=[0:5:360*5])[10*cos(i),10*sin(i),i/360*5]]; 
 prism2=p_extrude(sec2,path2);
 translate([35,0,0])
 swp(prism2);"
t(m)	" function to transpose a 3 x 3 matrix
 example:
 v1=[2,3,5];
 v2=[7,8,9];
 v3=[10,11,12];
 echo(t([v1,v2,v3])); // => ECHO: [[2, 7, 10], [3, 8, 11], [5, 9, 12]]"
trns(m,sec)	"function to translate a group of points ""sl"" by ""m"" distance defined in [x,y,z].e.g. try following code:
sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5);
p_line3dc(trns([2,5,10],sec),.1);"
uv(v)	" function to calculate a unit vector for a given vector
 example: 
 echo(uv([2,3,4])); // => ECHO: [0.371391, 0.557086, 0.742781]"
v_sec_extrude(sec,path,o)	" function to extrude a section along a path by varying section defined by offset ""o""
 example:
 sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5);

 path=c2t3(arc(20,0,120,s=20));

 p_line3d(path,.2);

 prism=v_sec_extrude(sec,path,-2);

 swp(prism);"
![image](https://user-images.githubusercontent.com/55306937/159120168-64551408-30cc-4bf5-a860-ae7c9b1c3152.png)
