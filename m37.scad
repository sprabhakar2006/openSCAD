include<dependencies.scad>

sec=[for(i=[-20:20])[i,3*sin(i*18)]];

//p_lineo(sec,.1);

path=[for(i=[-30:30])[i,0,20+3*sin(i*12)]];

//p_line3d(path,.1);
surf=surf_extrude(sec,path);
//%surf_extrude(sec,path);

sec1=cr(pts1([[0,0,1.9],[20,0,1.9],[0,20,1.9]]));

intersection(){
surf_base(surf);
swp(l_extrude(sec1,100));

}
