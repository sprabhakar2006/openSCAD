 

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
    
module p_line3d(path,r,rec=0,$fn=20){
    for(i=[0:len(path)-2])
        
    hull(){
    translate(path[i])if(rec==0)sphere(r); else cube(r*2,true);
    translate(path[i+1])if(rec==0)sphere(r);else cube(r*2,true);
    }}

//module to draw a polyline in 3d space (loop closed)
// e.g. try following code:
// sec=trns([5,10,6],q_rot(["x45"],circle(10)));
// p_line3dc(sec,.2);    

module p_line3dc(path,r,rec=0,$fn=20){
    for(i=[0:len(path)-1])
        let(
    i_plus=i<len(path)-1?i+1:0
    )
    hull(){
    translate(path[i])if(rec==0)sphere(r); else cube(r*2,true);
    translate(path[i_plus])if(rec==0)sphere(r);else cube(r*2,true);
    }}