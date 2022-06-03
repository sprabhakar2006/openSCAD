include<dependencies.scad>

sec=m_points(square(60,true),5);
sec1=circle(.1,s=72);

sec2=sort_points(sec,sec1);

prism=q_rot(["y20"],[sec2,sec]);

prism1=sphere(20);
prism2=q_rot(["y20"],l_extrude(sec,.01));
prism3=ipf(prism1,prism,2,s=10);

%swp(prism1);
%swp(prism2);
swp_c(prism3);

