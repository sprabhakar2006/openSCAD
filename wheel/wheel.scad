include<dependencies.scad>
include<spoke1a.scad>
include<spoke2a.scad>
include<spoke3a.scad>
include<hub and rim.scad>


//spoke1a();
//spoke2a();
//spoke3a();
//hub_rim();

difference(){
union(){
import("hub and rim.stl");
for(i=[0:30:330])rotate([0,0,i])
import("spoke3.stl");
for(i=[0:60:300])rotate([0,0,i]){
import("spoke2.stl");
import("spoke1.stl");}}

swp(cyl(r=20,h=300));}