include<dependencies.scad>

sec=m_points_so(cr(pts1([[0,10],[0,5],[12,0,5],[10,-5],[23,0]]),5),40,.1);
max_x=max(sec*[1,0]);
//p_line(sec,.1);
echo(max_x);

sec1=[for(i=[0:len(sec)-1])[13*cos(sec[i].x/max_x*360),13*sin(sec[i].x/max_x*360),sec[i].y]];
sec2=[for(i=[0:len(sec)-1])[10*cos(sec[i].x/max_x*360),10*sin(sec[i].x/max_x*360),sec[i].y]];

sec3=cpo([sec1,sec2]);

sec4=loop(sec,44,131);

//p_lineo(sec4);

path=[for(i=[0:len(sec4)-1])[11.5*cos(sec4[i].x/max_x*360),11.5*sin(sec4[i].x/max_x*360),sec4[i].y]];
sec5=cir(.8,s=30);

prism=[for(i=[0:len(path)-2])let(
i_plus=i+1,
p0=path[i],
p1=path[i_plus],
v1=p1-p0,
v2=[p0.x,p0.y],
a1=ang(norm([v1.x,v1.y]),v1.z),
a2=ang(v2.x,v2.y)+90,
ang1=is_num(a1)?a1:0,
ang2=is_num(a2)?a2:0,
sec=q_rot(["x90","z-90"],sec5),
a_sec=trns(p0,q_rot([str("y",ang1),str("z",ang2)],sec))
)a_sec];


difference(){
surf_base(flip(sec3));
swp(flip(prism));}


