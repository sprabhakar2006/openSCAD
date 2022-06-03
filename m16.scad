include<dependencies.scad>
   
sec=cr([[20,0,3],[40,0,3],each arc(40,5,35,s=10),[40*cos(40),40*sin(40),3],[20*cos(40),20*sin(40),3],each arc(20,35,5,s=10)],5);
        

prism=trns([0,0,.001],l_extrude(m_points_sc(sec,5,.5),100));


prism1=cylinder(r1=41.5,r2=18,h=10);  

plane=trns([0,0,20],plane([0,0,1],200));

//for(p=prism)points(p,.2);

ip1=ipe(prism1,prism,1,1);
sec1=c3t2(ip1[len(ip1)-1]);
prism2=l_extrude(sec1,100);
ip2=flip(ipe(plane,flip(prism2),1,1));

prism3=[each ip1,each ip2];

//swp(ip1);
//swp(ip2);
swp(prism3);
