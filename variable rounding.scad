include<dependencies.scad>

sec=[[-20,0],[20,0]];
path=cytz(m_points([[-40,-5],[40,30]],2));
//p_line3d(sec,.2);
prism=surf_extrude(sec,path);
//surf_extrude(sec,path);
sec1=cr(pts1([[-5,-7.5,1],[10,0,4],[0,15,1],[-10,0,4]]),20);
path1=cr(pts1([[-2,0],[2,0,1],[0,2]]),15);
prism1=l_extrude(m_points(sec1,1),100);
//
fillet1=ipe(prism,flip(prism1),1,1,s=15);

fillet2=ipe(plane([0,0,1],100),prism1,1,0,s=15);

prism2=[each fillet2, each flip(fillet1)];

swp(prism2);
