include<dependencies.scad>

sec=cr(pts1([[-10,0],[0,2,1],[20,0,1],[0,-2]]),10);
sec1=cr(pts1([[-10,2],[0,-2,1],[20,0,1],[0,2]]),10);
path=cytz(cr(pts1([[20,10],[15,0.001,4],[0.001,20,.1],[30,0.001,.1],[0.001,-5]]),10));
path2=cytz(cr(pts1([[20,0],[47,.01,2],[0.01,7]]),10));
sec2=cr(pts1([[0,-10,10],[30,0],[0,20],[-30,0,10]]),30);
path1=cr(pts1([[-2,0],[2,0,1],[0,12,1],[-2,0]]),10);
prism1=prism(sec2,path1);
surf=surf_extrude(sec,path);
surf1=surf_extrude(sec1,path2);

render(){
intersection(){
surf_base(surf);
surf_base(surf1,50);
}
swp(prism1);}