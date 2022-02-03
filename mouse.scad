include<dependencies.scad>

sec=2p_arc([[-45,0],[45,0]],60,s=30);
path=cytz(cr(pts1([[0,0],[70,30,60],[50,0,40],[10,-30]]),20));


prism=surf_extrude(sec,path);


sec1=cr(pts1([[0,-40,40],[147,0,40],[0,80,40],[-147,0,40]]),30);
prism1=l_extrude(m_points(sec1,5),100);


path1=ip(prism,flip(prism1));


prism2=ipe(prism,flip(prism1),10,1);
prism3=[for(i=[len(prism2[0])-1:-1:0])[for(p=prism2)p[i]]];

prism4=[for(p=prism3)sort_points(m_points(sec1,5),p)];

prism5=[for(i=[8:2:40])m_points(f_offset(sec1,-i),5)];
prism6=[for(p=prism5)sort_points(m_points(sec1,5),ip(prism,l_extrude(p,100)))];

path2=cr(pts1([[-3,0],[3,0,3],[3,3]]),5);
prism7=prism(f_offset(sec1,-10),path2);
prism71=[for(p=prism7)m_points(p,5)];

prism8=[each prism71, each prism4, each prism6];

swp(prism8);