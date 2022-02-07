include <dependencies.scad>

prism=cyl(r=20,h=50,s=72);
%swp(cyl(r=20,h=50,s=72));

sec1=[[-30,0],[30,0]];
    
path1=m_points([[-30,0,15],[30,0,30]],10);

%surf_extrude(sec1,path1);
prism1=surf_extrude(sec1,path1);

prism2=ipf(prism1,prism,2,1);
prism3=ipf(prism1,flip(prism),2);

swp_c(prism2);
swp_c(prism3);

