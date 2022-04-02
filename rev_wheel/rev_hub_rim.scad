include<dependencies.scad>


path=cr(pts1([[245,0],[0,5],[-13,20,5],[0,20,5],[-12,20,5],[0,78,5],[12,20,5],[0,64,5],[13,20,2],[0,5]]),5);
path1=flip([for(p=path) [5,0]+p]);
path2=[each path,each path1];
rotate_extrude($fn=100)
polygon(path2);
sec4=cir(50);
path4=cr(pts1([[-5,200],[5,0,5],[0,35,5],[-25,20,5],[-5,0]]),10);
prism7=prism(sec4,path4);
prism8=cyl(r=18,h=300);

difference(){
swp(prism7);
swp(prism8);}