include<dependencies.scad>
sketch=[for(i=[-25:25])[i,2*sin(360/25*i)]];
path=cytz([for(i=[-35:35])[i,10+4*sin(360/35*i)]]);
surf_base(surf_extrude(sketch,path));