include<dependencies.scad>

    
sec1=arc(20,0,359,s=72);
//sec1=cir(20,s=72);    
path1=[[0,0],[-5,25]];

surf1=prism(sec1,path1);

swp(surf1);

sec2=cr(pts1([[-25,0],[10,5,5],[10,-3,10],[10,5,5],[10,-8,7],[10,1]]),10);
    
path2=cytz(cr(pts1([[-55,5,0],[10,8,20],[20,-5,10],[20,8,20],[10,-9,20],[10,1,0]]),10));

surf2=surf_extrude(sec2,path2);

surf_extrude(sec2,path2);

prism=ipf(surf2,surf1,2,1);

swp_c(prism);
