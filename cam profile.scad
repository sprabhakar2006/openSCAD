include<dependencies.scad>

sec=remove_extra_points(m_points_so(cr(pts1([[0,10],[0,5],[12,0,5],[10,-5],[23,0]]),5),40,.1));
max_x=max(sec*[1,0]);
//p_line(sec,.1);
echo(max_x);

sec1=[for(i=[0:len(sec)-1])[13*cos(sec[i].x/max_x*360),13*sin(sec[i].x/max_x*360),sec[i].y]];
sec2=[for(i=[0:len(sec)-1])[10*cos(sec[i].x/max_x*360),10*sin(sec[i].x/max_x*360),sec[i].y]];

sec3=cpo([sec1,sec2]);

sec4=loop(sec,44,127);

//p_lineo(sec4);

path=[for(i=[0:len(sec4)-1])[11.5*cos(sec4[i].x/max_x*360),11.5*sin(sec4[i].x/max_x*360),sec4[i].y]];
sec5=circle(.8,s=50);

prism=p_extrude(sec5,path);

render()
difference(){
surf_base(flip(sec3));
swp(flip(prism));}



