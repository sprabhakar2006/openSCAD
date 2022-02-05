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
//
//for(p=fillet1)points(p,.2);
//swp(fillet1);

prism2=[for(i=[len(fillet1[0])-1:-1:0])sort_points(m_points(sec1,1),[for(p=fillet1)p[i]])];

prism3=prism(sec1,path1);
prism5=[for(p=prism3)sort_points(m_points(sec1,1),m_points(p,1))];
prism6=[each prism5,each prism2];

//for(p=prism5)points(p,.2);
//for(p=prism2)points(p,.2);

swp(prism6);