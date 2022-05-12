include <dependencies.scad>

stages =200;

stage_height = 1.25;

rad = 50;

f1 = 25;

f2 = 25;

phase1 = 0;

phase2 = 180;

height_depth=5;

depth1 = 20;

depth2 = 20;

thickness = 2;

bottom_thickness = 3;

myslices = 5;

angle_step=0.1;

// generate outer points

points_base = [for (i = [0:1:stages]) let(f = pow(sin(i/stages * 120),2) * 7 + 1, var = 1 , a=((sin((i/stages*360*f)%360) * 0.5 + 0.5) * (var * height_depth))) [for(j = [0:angle_step:360-angle_step]) [sin(j) * (rad+a+(pauw(sin(j *f1+phase1),0.5)*0.5+0.5)*depth1*i/stages+(pauw(sin(j *f2+phase2),0.5)*0.5+0.5)*depth2*(1-i/stages)), cos(j) *(rad+a+(pauw(sin(j *f1+phase1),0.5)*0.5+0.5)*depth1*i/stages+(pauw(sin(j *f2+phase2),0.5)*0.5+0.5)*depth2*(1-i/stages)),i*stage_height]]];

points_base_e=surf_offset(points_base,2);

swp_prism_h(points_base,points_base_e);

swp([points_base[0],points_base[1],points_base[2],points_base[3],points_base[4]]);

function pauw(x,p)=sign(x)*abs(x)^p;
