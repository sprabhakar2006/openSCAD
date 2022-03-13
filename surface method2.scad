include<dependencies.scad>
sketch=cr(pts1([[-25,0],[25,20,100],[25,-20]]),20);
path=cytz(cr(pts1([[0,-5],[50,30,50],[20,-25]]),20));
surf=surf_extrude(sketch,path);

sec=cr(pts1([[10,-20,20],[60,0,20],[0,40,20],[-60,0,20]]),30);
prism=[c2t3(m_points_sc(sec,10,.5)),trns([0,0,25],m_points_sc(f_offset(sec,5),10,.5))];
prism1=ipe(plane([0,0,1],200),prism,2);
prism2=ipe(surf,flip(prism),2,1,s=10);
x=57;
prism3=[for(x=[10:5:55])ip(surf,rsz3dc(prism,bb(prism)-[x,x/50*40,0]))];
prism4=[each prism1, each flip(prism2),each prism3];
swp(prism4);
