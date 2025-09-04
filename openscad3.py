from numpy import *
from numpy.linalg import *
import time
from scipy.spatial import cKDTree, Delaunay, ConvexHull
import sympy as sp
import math
# from numba import njit, jit
# from skimage import measure
import open3d as o3d
# from stl import mesh

def arc(radius=0,start_angle=0,end_angle=0,cp=[0,0],s=20):
    """
function for calculating 2d arc
'cp': center point of the arc
's': number of segments in the arc

example:
a1=arc(radius=10,start_angle=180,end_angle=270,cp=[10,10],s=20)
fo(f'''
color("blue") p_line3d({a1},.3);
''')
    """
    cp=array(cp)
    r=linspace(start_angle,end_angle,s+1)
    x=radius*cos(pi/180*r)
    y=radius*sin(pi/180*r)
    c=(cp+array([x,y]).swapaxes(0,1))
    return c.tolist()        

def pts(p):
    """
calculates the cumulative sum of 2d list of points 'p'
example:
pts([[0,0],[4,0],[2,3],[5,-8]]) will produce following output
[[0, 0], [4, 0], [6, 3], [11, -5]]
    """
    return array(p)[:,0:2].cumsum(axis=0).tolist()

turtle2d=pts


def pts1(p):
    """
'p' is a list of points
function calculates the cumulative sum of x,y values in the list while z value remains the same.
this is mainly used in function corner_radius(pl,s).
example:
a=pts1([[0,0,1],[10,0,1],[0,5,1],[-10,0,1]])
# if used with function corner_radius(pl,s)
b=corner_radius(a,5) # This is a rounded rectangle of dim 10 x 5 with corner radius of 1 at each corner
fo(f'''
color("blue") p_line3dc({b},.2);
''')
    """
    
    p=[[a[0],a[1],0] if len(a)==2 else a for a in p]
    b=array(p)[:,0:2].cumsum(axis=0)
    c=array([array(p)[:,2].tolist()])
    return concatenate((b,c.T),1).tolist()


    
def cw(sec):
    """
function to identify if an enclosed section is clockwise(cw) or counterclockwise(ccw)
this returns 1 if section is clockwise and -1 if it is counterclockwise
example:
cw([[0,0],[10,0],[5,5]]) => -1
cw([[0,0],[10,0],[5,-5]]) => 1
    """
    sec=remove_extra_points(sec)
    sec=rot2d(0.001,sec)
    if len(sec)==3:
        p=array(sec)
        return -1 if cross(c23(p[1]-p[0]),c23(p[2]-p[0]))[-1]>0 else 1
    else:
        cp1=array(sec).mean(0)
        p0=[array(sec)[:,0].min()-1,cp1[1]]
        v1=array([[1,0]]*len(sec))

        p1=array(sec)
        p2=array(sec[1:]+[sec[0]])
        v2=p2-p1
        iim=array([v1,-v2]).transpose(1,0,2).transpose(0,2,1)
        im=inv(iim)
        p=p1-p0
        t=einsum('ijk,ik->ij',im,p)

        t1=t[:,0][(t[:,0]>0)&(t[:,1]>=0)&(t[:,1]<=1)]
        ip1=(array(p0)+array([[1,0]])*t1[:,None])
        a=ip1[:,0]
        b=sort(ip1[:,0])
        ip2=ip1[a==b[0]].tolist()+ip1[a==b[1]].tolist()
        v3=ip2[1]-p1
        u3=(v3/norm(v3,axis=1).reshape(-1,1)).round(3)
        u2=(v2/norm(v2,axis=1).reshape(-1,1)).round(3)
        p3=array(p1)[(norm(v3,axis=1)<=norm(v2,axis=1))&(u2==u3).all(1)].tolist()[0]
        n=arange(len(sec))[(p1==p3).all(1)][0]
        p4=sec[n+1 if n<len(sec)-1 else 0]

        cw=-1 if cross(c23(array(p4)-array(p3)),c23(array(ip2[0])-array(p3)))[-1]>0 else 1
        return cw

check_orientation=cw

def cwv(sec):
    """
function to identify whether each point in a section is clockwise or counter clockwise. 
cw(sec)==1 means clockwise and -1 means counterclockwise. 
example:
cwv(pts([[0,0],[4,0],[0,4],[2,0],[0,2],[-6,0]])) => [-1, -1, 1, -1, -1, -1]
    """
    p=sec
    p0=[p[-1]]+p[:-1]
    p1=p
    p2=p[1:]+[p[0]]
    p0,p1,p2=array([p0,p1,p2])
    p=array([p0,p1,p2]).transpose(1,0,2)
    return [ -1 if cross(c23(p1[1]-p1[0]),c23(p1[2]-p1[0]))[-1]>=0 else 1 for p1 in p]

check_each_point_orientation=cwv

def ang(x,y):
    """
function to calculate angle of a 2d vector starting from origin and end point with x and y co-ordinates
 example:
 p1,p2=array([[3,4],[-3,2]])
 v=p2-p1
 ang= ang(v[0],v[1])
 
    """
    if x>=0 and y>=0:
        return l_(arctan(y/(0.000001 if x==0 else x))*180/pi)
    elif x<0 and y>=0:
        return l_(180-abs(arctan(y/x))*180/pi)
    elif  x<0 and y<0:
        return l_(180+abs(arctan(y/x))*180/pi)
    else:
        return l_(360-abs(arctan(y/(0.000001 if x==0 else x)))*180/pi)


def flip(sec):
    """
function to flip the sequence of a list or a list of points
example:
list=[1,2,3,4,5]
flipped_list=flip(list) => [5, 4, 3, 2, 1]

list=[[1,2,3],[4,5,6],[7,8,9]]
flipped_list=flip(list) => [[7, 8, 9], [4, 5, 6], [1, 2, 3]]
    """
    return sec[::-1]
    

def max_r(sec):
    """
function calculates the maximum radius in a given closed section
example:
sec=cr_c(pts1([[0,0,.2],[8,3,3],[5,7,1],[-8,0,2],[-5,20,1]]),20)
max_r(sec) => 3.0
    
    """
    c=[]
    for i in range(len(sec)):
        i_2minus=len(sec)-2 if i==0 else len(sec)-1 if i==1 else i-2
        i_minus=len(sec)-1 if i==0 else i-1
        i_plus=i+1 if i<len(sec)-1 else 0
        i_2plus=i+2 if i<len(sec)-2 else 0 if i<len(sec)-1 else 1
        pi_2minus=sec[i_2minus]
        pi_minus=sec[i_minus]
        pi=sec[i]
        pi_plus=sec[i_plus]
        pi_2plus=sec[i_2plus]
        v1=subtract(pi_minus,pi_2minus)
        v2=subtract(pi,pi_minus)
        v3=subtract(pi_plus,pi)
        v4=subtract(pi_2plus,pi_plus)
        l1=norm(v1).round(3)
        l2=norm(v2).round(3)
        l3=norm(v3).round(3)
        l4=norm(v4).round(3)
        r1=r_3p([pi_2minus,pi_minus,pi]).round(3)
        r2=r_3p([pi_minus,pi,pi_plus]).round(3)
        r3=r_3p([pi,pi_plus,pi_2plus]).round(3)
        c.append(0 if l2!=l3 and (r1!=r2 or r2!=r3) else r2)
    return l_(max(c))
        
max_radius=max_r

def offset_l(l,d):
    u=uv(subtract(l[1],l[0]))
    p0=add(l[0],dot(u,multiply(d,rm(-90)))).tolist()
    p1=add(l[1],dot(u,multiply(d,rm(-90)))).tolist()
    return [p0,p1]

def seg(sec):
    """
function to create a segment from a list of points or a list
example:
list=[1,2,3,4,5,6]
seg(list)=> [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 1]]

list=[[1,2,3],[4,5,6],[7,8,9]]
seg(list) => [[[1, 2, 3], [4, 5, 6]], [[4, 5, 6], [7, 8, 9]], [[7, 8, 9], [1, 2, 3]]]
    """
    c=[]
    for i in range(len(sec)):
        i_plus=i+1 if i<len(sec)-1 else 0
        p0=sec[i]
        p1=sec[i_plus]
        l=[p0,p1]
        c.append(l)
    return c


def offset_segv(sec,d):
    """
function makes the segments of the original section and offset each segment by a distance 'd'
    
    """
    s=sec
    s1=s[1:]+[s[0]]
    x=(array(s1)-array(s))
    y=norm(x,axis=1)
    u=x/y.reshape(-1,1)
    p0=array(s)+u@array(rm(-90))*d
    p1=array(s1)+u@array(rm(-90))*d
    return swapaxes([p0,p1],0,1).tolist()

def offset_points(sec,r):
    """
function to calculate offset of a list of 2d points
in defining sections, providing corner radius is a must
Not an important function to study
    """
    return array(offset_segv(sec,r))[:,0].tolist()


def offset_seg_cw(sec,r):
    """
function offsets the segment only when the point is clockwise
    
    """
    c=[]
    for i in range(len(sec)):
        i_minus=len(sec)-1 if i==0 else i-1
        i_plus=i+1 if i<len(sec)-1 else 0
        p0=sec[i_minus]
        p1=sec[i]
        p2=sec[i_plus]
        clock=cw([p0,p1,p2])
        if clock==1:
            c.append(offset_l([p1,p2],r))
    d=[]
    for a in c:
        for b in a:
            d.append(b)
    return d


def remove_extra_points(points_list):
    """
function removes all the duplicates from a 
example:
list=[9,5,1,2,3,4,5,2,4,6,9]
remove_extra_points(list) => [9, 5, 1, 2, 3, 4, 6]

list=[[7,8,9],[1,2,3],[10,11,12],[4,5,6],[7,8,9],[10,11,12]]
remove_extra_points(list) => [[7, 8, 9], [1, 2, 3], [10, 11, 12], [4, 5, 6]]
    """
    return array(points_list)[sort(unique(points_list,axis=0,return_index=True)[1])].tolist()

remove_duplicates = remove_extra_points
    


def list_r(sec):
    """
function list the corner radiuses of a given section (only where the radius is specified)
example:
sec=corner_radius(pts1([[0,0,.1],[7,5,2],[5,7,3],[-5,7,5],[-7,5,5]]),5)
list_r(sec) => 
array([0.   , 0.1  , 0.1  , 0.1  , 0.1  , 0.   , 0.   , 2.   , 2.   ,
   2.   , 2.   , 0.   , 0.   , 3.   , 3.   , 3.   , 3.   , 0.   ,
   0.   , 4.077, 4.077, 4.077, 4.077, 4.077, 4.077, 4.077, 4.077,
   4.077, 0.   , 0.   ])
    """
    pi_2minus=sec[-2:]+sec[:-2]
    pi_minus=[sec[-1]]+sec[:-1]
    p_i=sec
    pi_plus=sec[1:]+[sec[0]]
    pi_2plus=sec[2:]+sec[:2]

    v1=array(pi_minus)-array(pi_2minus)
    v2=array(p_i)-array(pi_minus)
    v3=array(pi_plus)-array(p_i)
    v4=array(pi_2plus)-array(pi_plus)

    l1=norm(v1,axis=1).round(3)
    l2=norm(v2,axis=1).round(3)
    l3=norm(v3,axis=1).round(3)
    l4=norm(v4,axis=1).round(3)

    p4=array(pi_2minus)+(array(pi_minus)-array(pi_2minus))/2
    p5=array(pi_minus)+(array(p_i)-array(pi_minus))/2

    u1=(array(pi_minus)-p4)/norm(array(pi_minus)-p4,axis=1).reshape(-1,1)
    u2=(array(p_i)-p5)/norm(array(p_i)-p5).reshape(-1,1)

    v5=array(pi_minus)-p4
    v6=(v5/norm(v5,axis=1).reshape(-1,1))
    r1=r_3pv(array(pi_2minus),array(pi_minus),array(p_i)).round(3)
    r2=r_3pv(array(pi_minus),array(p_i),array(pi_plus)).round(3)
    r3=r_3pv(array(p_i),array(pi_plus),array(pi_2plus)).round(3)
    r=where((l2!=l3) & ((r1!=r2) | (r2!=r3)),0,r2)
    return r

def list_ra(sec):
    """
calculates list of radiuses for all the points of a given section
example:
sec=corner_radius(pts1([[0,0,.1],[7,5,2],[5,7,3],[-5,7,5],[-7,5,5]]),5)
list_ra(sec) =>
array([  0.,   0.,   0.,   0.,   0.,  19., 124.,   2.,   2.,   2.,   2.,
    95.,  28.,   3.,   3.,   3.,   3.,  26.,  92.,   4.,   4.,   4.,
     4.,   4.,   4.,   4.,   4.,   4.,  40.,   8.])

    """
    pi_2minus=sec[-2:]+sec[:-2]
    pi_minus=[sec[-1]]+sec[:-1]
    p_i=sec
    pi_plus=sec[1:]+[sec[0]]
    pi_2plus=sec[2:]+sec[:2]

    v1=array(pi_minus)-array(pi_2minus)
    v2=array(p_i)-array(pi_minus)
    v3=array(pi_plus)-array(p_i)
    v4=array(pi_2plus)-array(pi_plus)

    l1=norm(v1,axis=1).round(3)
    l2=norm(v2,axis=1).round(3)
    l3=norm(v3,axis=1).round(3)
    l4=norm(v4,axis=1).round(3)

    p4=array(pi_2minus)+(array(pi_minus)-array(pi_2minus))/2
    p5=array(pi_minus)+(array(p_i)-array(pi_minus))/2

    u1=(array(pi_minus)-p4)/norm(array(pi_minus)-p4,axis=1).reshape(-1,1)
    u2=(array(p_i)-p5)/norm(array(p_i)-p5).reshape(-1,1)

    v5=array(pi_minus)-p4
    v6=(v5/norm(v5,axis=1).reshape(-1,1))
    r1=r_3pv(array(pi_2minus),array(pi_minus),array(p_i)).round(3)
    r2=r_3pv(array(pi_minus),array(p_i),array(pi_plus)).round(3)
    r3=r_3pv(array(p_i),array(pi_plus),array(pi_2plus)).round(3)
    r=where((l2!=l3) & ((r1!=r2) | (r2!=r3)),0,r2)
    return r2

list_radiuses=list_ra
    
def rm(theta):
    """
function to rotate a vector by "theta" degrees e.g. try following code:
example:
line=[[0,0],[5,0]]
line1=array(line)@rm(30)
line1=line1.tolist()
fo(f'''
color("blue") p_line3d({line},.2);
color("magenta") p_line3d({line1},.2);
''')
    """
    
    pi=3.141592653589793
    return [[cos(theta * pi/180),sin(theta * pi/180)],[-sin(theta * pi/180),cos(theta * pi/180)]]

def max_rv(sec):
    """
function calculates the maximum radius in a given closed section
example:
sec=corner_radius(pts1([[0,0,.2],[8,3,3],[5,7,1],[-8,0,2],[-5,20,1]]),20)
max_rv(sec) => 3.0
    
    """
    pi_2minus=sec[-2:]+sec[:-2]
    pi_minus=[sec[-1]]+sec[:-1]
    p_i=sec
    pi_plus=sec[1:]+[sec[0]]
    pi_2plus=sec[2:]+sec[:2]

    v1=array(pi_minus)-array(pi_2minus)
    v2=array(p_i)-array(pi_minus)
    v3=array(pi_plus)-array(p_i)
    v4=array(pi_2plus)-array(pi_plus)

    l1=norm(v1,axis=1).round(3)
    l2=norm(v2,axis=1).round(3)
    l3=norm(v3,axis=1).round(3)
    l4=norm(v4,axis=1).round(3)

    p4=array(pi_2minus)+(array(pi_minus)-array(pi_2minus))/2
    p5=array(pi_minus)+(array(p_i)-array(pi_minus))/2

    u1=(array(pi_minus)-p4)/norm(array(pi_minus)-p4,axis=1).reshape(-1,1)
    u2=(array(p_i)-p5)/norm(array(p_i)-p5).reshape(-1,1)

    v5=array(pi_minus)-p4
    v6=(v5/norm(v5,axis=1).reshape(-1,1))
    r1=r_3pv(array(pi_2minus),array(pi_minus),array(p_i)).round(3)
    r2=r_3pv(array(pi_minus),array(p_i),array(pi_plus)).round(3)
    r3=r_3pv(array(p_i),array(pi_plus),array(pi_2plus)).round(3)
    return max(where((l2!=l3) & ((r1!=r2) | (r2!=r3)),0,r2))

def r_3p(p):
    """
function calculates radius of the circle drawn with 3 points 'p1','p2','p3'
example:
p1,p2,p3=[3,0],[0,0],[0,3]
radius=r_3p([p1,p2,p3]) => 2.1213203435596424
    """
    p4=add(p[0],divide(subtract(p[1],p[0]),2)).tolist()
    p5=add(p[1],divide(subtract(p[2],p[1]),2)).tolist()
    u1=uv(subtract(p[1],p4))
    u2=uv(subtract(p[2],p5))
    p6=add(p4,dot(u1,rm(90))).tolist()
    p7=add(p5,dot(u2,rm(90))).tolist()
    cp=i_p2d([p4,p6],[p5,p7])
    r=norm(subtract(p[0],cp))
    return r
    
radius_3points=r_3p




def s_int(sec1):
    """
calulates the self intersection points of a list of line segments 's'
it also picks the points in case the 2 lines are just connected at 1 point and are not crossing
 function s_int1 is applicable in most of the calculations
    """
    n=len(sec1)
    a=array(sec1)[comb_list(n)]
    p0=a[:,0][:,0]
    p1=a[:,0][:,1]
    p2=a[:,1][:,0]
    p3=a[:,1][:,1]
    v1=p1-p0
    v2=p3-p2
    iim=array([v1,-v2+.00001]).transpose(1,0,2).transpose(0,2,1)
    im=inv(iim)
    p=p2-p0

    t=einsum('ijk,ik->ij',im,p)
    dcn=(t[:,0].round(4)>=0)&(t[:,0].round(4)<=1)&(t[:,1].round(4)>=0)&(t[:,1].round(4)<=1)
    i_p1=p0+einsum('ij,i->ij',v1,t[:,0])
    i_p1=i_p1[dcn].tolist()
    return i_p1


def r_3pv(p1,p2,p3):
    p4=p1+(p2-p1)/2
    p5=p2+(p3-p2)/2
    u1=(p2-p4)/norm(p2-p4,axis=1).reshape(-1,1)
    u2=(p3-p5)/norm(p3-p5,axis=1).reshape(-1,1)
    p6=p4+u1@array([[0,1],[-1,0]])
    p7=p5+u2@array([[0,1],[-1,0]])
    cp=i_p2dv(p4,p6,p5,p7)
    r=norm(p1-cp,axis=1)
    return r


def sort_points(sec,list1):
    """
function picks the nearest point of a section from a reference section and matches the length of points for the 2 compared sections
example:
sec=circle(10)
sec1=corner_radius(pts1([[-2.5,-2.5,1],[5,0,1],[0,5,1],[-5,0,1]]))
sec2=sort_points(sec,sec1)
sec3=cpo([sec2,sec])
fo(f'''
p_line({sec},.2);
p_line({sec2},.2);
color("blue") for(p={sec3})p_line(p,.1);
    ''')
    """
    return array(list1)[cKDTree(list1).query(sec)[1]].tolist()
            

def m_points(sec,pitch=1):
    """
multiple points within straight lines of a closed section 'sec' with equal segment length 'pitch' in the straight line segments
example:
l1=[[0,0],[10,0],[10,10],[0,10]]
l2=translate_2d([20,0],m_points(l1,pitch=1))
fo(f'''
color("blue") points({l1},.2);
color("magenta") points({l2},.2);
''')
    """
    sec1=[]
    for p in seg(sec):
        n=1 if int(round(l_len(p)/pitch,0))==0 else int(round(l_len(p)/pitch,0))
        sec1.append(ls(p,n))
    return concatenate(sec1).tolist()



def m_points_o(sec,pitch=1):
    """
multiple points within straight lines of a open section 'sec' with equal segment length 'pitch' in the straight line segments
l1=[[0,0],[10,0],[10,10],[0,10]]
l2=translate_2d([20,0],m_points_o(l1,pitch=1))
fo(f'''
color("blue") points({l1},.2);
color("magenta") points({l2},.2);
''')
    """
    sec1=[]
    for p in seg(sec)[:-1]:
        n=1 if int(round(l_len(p)/pitch,0))==0 else int(round(l_len(p)/pitch,0))
        sec1.append(ls(p,n))
    return concatenate(sec1).tolist()+[sec[len(sec)-1]]


def ls(line,n):
    """
function to draw number of points 'n' in a line 'line'
example:
line=[[0,0],[10,0]]
line1=ls(line,5) => [[0.0, 0.0], [2.0, 0.0], [4.0, 0.0], [6.0, 0.0], [8.0, 0.0], [10.0, 0.0]]
    """
    p0,p1=array(line)
    v1=p1-p0
    return array([p0+v1/n*i for i in range(n)]).tolist()


def l_len(l):
    """
calculates length of a line 'l'
example:
line=[[0,0],[10,0]]
l_len(line) =>10
    """
    p0,p1=array(l[0]),array(l[1])
    v=p1-p0
    u=[v/(norm(v)+.00001)]
    length=norm(v)
    return length.tolist()

line_length=l_len
    
def arc_2p(p1,p2,r,cw=1,s=20):
    """
arc with 2 points 'p1,p2' with radius 'r' and with orientation clockwise (1) or counterclock wise(-1)
example:
p1=[4,0]
p2=[0,4]
arc1=arc_2p(p1,p2,3,1,20)
arc2=arc_2p(p1,p2,3,-1,20)
fo(f'''
// arc clockwise
color("blue")
p_line3d({arc1},.1);

// arc counter clockwise
color("magenta")
p_line3d({arc2},.1);

color("cyan")points({[p1,p2]},.3);
    
''')
    """
    p1,p2=array([p1,p2])
    p3=p1+(p2-p1)/2
    d=norm(p3-p1)
    l=sqrt(abs(r**2-d**2))
    v=p1-p3
    u=v/norm(v)
    cp=p3+(u*l)@rm(-90 if cw==-1 else 90)
    v1,v2=p1-cp,p2-cp
    a1,a2=ang(v1[0],v1[1]),ang(v2[0],v2[1])
    a3= (a2+360 if a2<a1 else a2) if cw==-1 else (a2 if a2<a1 else a2-360)
    return arc(r,a1,a3,cp,s)

arc_2points_short=arc_2p
    
def arc_long_2p(p1,p2,r,cw=1,s=20):
    """
arc long with 2 points 'p1,p2' with radius 'r' and with orientation clockwise (1) or counterclock wise(-1)
example:
p1=[4,0]
p2=[0,4]
arc1=arc_long_2p(p1,p2,3,1,20)
arc2=arc_long_2p(p1,p2,3,-1,20)
fo(f'''
// arc clockwise
color("blue")
p_line3d({arc1},.1);

// arc counter clockwise
color("magenta")
p_line3d({arc2},.1);

color("cyan")points({[p1,p2]},.3);
''')
    """
    p1,p2=array([p1,p2])
    p3=p1+(p2-p1)/2
    d=norm(p3-p1)
    l=sqrt(abs(r**2-d**2))
    v=p1-p3
    u=v/norm(v)
    cp=p3+(u*l)@rm(90 if cw==-1 else -90)
    v1,v2=p1-cp,p2-cp
    a1,a2=ang(v1[0],v1[1]),ang(v2[0],v2[1])
    a3=(a2+360 if a2<a1 else a2) if cw==-1 else (a2 if a2<a1 else a2-360)
    return arc(r,a1,a3,cp,s)

arc_2points_long=arc_long_2p
    
def cir_2p(p1,p2,r,cw=1,s=20):
    """
circle with 2 points 'p1,p2' with radius 'r' and with orientation clockwise (1) or counterclock wise(-1)
example:
p1=[5,2]
p2=[7,8]
r=10
c1=cir_2p(p1,p2,r,cw=-1,s=50)
cp1=cp_arc(c1)
fo(f'''
color("cyan")points({[p1,p2]},.5);
color("blue")p_line3dc({c1},.1);
color("magenta")points({[cp1]},.5);
''')
    """
    p1,p2=array([p1,p2])
    p3=p1+(p2-p1)/2
    d=norm(p3-p1)
    l=sqrt(abs(r**2-d**2))
    v=p1-p3
    u=v/norm(v)
    cp=p3+(u*l)@rm(-90 if cw==-1 else 90)
    v1,v2=p1-cp,p2-cp
    a1,a2=ang(v1[0],v1[1]),ang(v2[0],v2[1])
    a3= (a2+360 if a2<a1 else a2) if cw==-1 else (a2 if a2<a1 else a2-360)
    return arc(r,a1,a1+360,cp,s)

circle_2points=cir_2p


def arc_2p_cp(p1,p2,r,cw=-1):
    """
center point of an arc with 2 points 'p1,p2' with radius 'r' and with orientation clockwise (1) or counterclock wise(-1).
example:
p1=[2,0]
p2=[0,2]
arc1=arc_2p(p1,p2,r=2,cw=1,s=20)
cp1=arc_2p_cp(p1,p2,r=2,cw=1)
 
fo(f'''
// arc clockwise with center point of the arc
color("blue")
{{p_line3d({arc1},.05);
points({[cp1]},.2);}}
''') 
    """
    p1,p2=array([p1,p2])
    p3=p1+(p2-p1)/2
    d=norm(p3-p1)
    l=sqrt(abs(r**2-d**2))
    v=p1-p3
    u=v/norm(v)
    cp=p3+(u*l)@rm(-90 if cw==-1 else 90)
    return cp.tolist()

center_arc_2points=arc_2p_cp
    
def offset_2(sec,r):
    """
calculates offset for a section 'sec' by amount 'r'
    """

    if convex(sec):
        if r <0:
            return inner_concave_offset(sec,r)
        elif r >0:
            return outer_convex_offset(sec,r)
        elif r==0:
            return sec
    else:
        if r<0:
            return inner_concave_offset(sec,r)
        elif r>0:
            return outer_concave_offset(sec,r)
        elif r==0:
            return sec


def prism(sec,path,type=1):
    """
function to make a prism with combination of 2d section and 2d path
type: can be set to either "1" or "2" which means prism calculation is done by offset method "1" or offset method "2", by defualt it is "1"
example:
sec=corner_radius(pts1([[-15,0,2.4],[0,15,3],[30,0,3],[0,-15,2.4],[5,0,2.4],[0,20,7],[-40,0,7],[0,-20,2.4]]),10)
path=corner_radius(pts1([[2,0],[-2,0,2],[0,7,2],[-2,0]]),10)
sol=prism(sec,path)
fo(f'''
{swp(sol)}
''') 
    """
    s1=flip(sec) if cw(sec)==1 else sec
    return [array(translate([0,0,y],offset(s1,x,type))).tolist() for (x,y) in path]



def f_prism(sec,path):
    """
function to make a prism with combination of 2d section and 2d path.
use prism function instead of this
    """
    s1=flip(sec) if cw(sec)==1 else sec
    return [translate([0,0,y],oset(s1,x)) for (x,y) in path]




def translate(p,sec):#translates a prism or section by [x,y,z] distance
    """
function to translate a group of points "sec" by "p" distance defined in [x,y,z].
example:
sec=corner_radius([[0,0,.5],[10,0,2],[7,15,1]],5)
sec1=translate(p=[10,5,3],sec=sec)
fo(f'''
color("blue")p_line3dc({sec},.1); 
color("magenta")p_line3dc({sec1},.1);
''') 
    """
    return (array(c2t3(sec))+c2t3(p)).tolist()

def translate_2d(p,sec):#translates a 2d section by [x,y] distance
    """
function to translate a group of points "sec" by "p" distance defined in [x,y].
example:
sec=corner_radius([[0,0,.5],[10,0,2],[7,15,1]],5)
sec1=translate_2d(p=[10,5],sec=sec)
fo(f'''
color("blue")p_line3dc({sec},.1); 
color("magenta")p_line3dc({sec1},.1);
''') 
    """
    return c3t2((array(c2t3(sec))+c2t3(p)))

    
def offset_points_cw(sec,r):
    """
function to offset only those points which are clockwise 
    """
    s=seg(sec)
    c=[]
    for i in range(len(sec)):
        i_minus=len(sec)-1 if i==0 else i-1
        i_plus=i+1 if i<len(sec)-1 else 0
        p0=sec[i_minus]
        p1=sec[i]
        p2=sec[i_plus]
        if cw([p0,p1,p2])==1:
            c.append(offset_l([p1,p2],r)[0])
    return c

def offset_points_ccw(sec,r):
    """
function to offset only those points which are counter clockwise
    """
    s=seg(sec)
    c=[]
    for i in range(len(sec)):
        i_minus=len(sec)-1 if i==0 else i-1
        i_plus=i+1 if i<len(sec)-1 else 0
        p0=sec[i_minus]
        p1=sec[i]
        p2=sec[i_plus]
        if cw([p0,p1,p2])==-1:
            c.append(offset_l([p1,p2],r)[0])
    return c

def cytz(path):# converts 'y' points to 'z' points in a 2d list of points
    """
function to convert the y co-ordinates to z co-ordinates e.g.[x,y]=>[x,0,y]. 2d to 3d coordinate system
example:
a=circle(20)
b=cytz(a)
fo(f'''
color("blue") p_line3dc({a},.3);
color("magenta") p_line3dc({b},.3);
''')
    """
    return [[p[0],0,p[1]] for p in path]

convert_y_to_z_coordinate=cytz

def surf_extrude(sec,path):# extrudes an open section 'sec' to a 'path' to create surface
    """
function to make surface with a polyline 2d sketch and a 3d path
use path_extrude_open instead of this
    """
    p0=path
    p1=p0[1:]+[p0[0]]
    p0,p1=array(p0),array(p1)
    v=p1-p0
    a1=vectorize(ang)(v[:,0],v[:,1])
    b=sqrt(v[:,0]**2+v[:,1]**2)
    a2=vectorize(ang)(b,v[:,2])
    c=[]
    for i in range(len(path)-1):
        sec1=translate(p0[i],rot(f'x90z-90y{-a2[i]}z{a1[i]}',sec))
        sec2=translate(p1[i],rot(f'x90z-90y{-a2[i]}z{a1[i]}',sec))
        if i<len(path)-2:
            c.append([sec1])
        else:
            c.append([sec1,sec2])
    return concatenate(c).tolist()

def cpo(prism): # changes the orientation of points of a prism
    """
function to change the orientation of the points of the prism
a=sphere(20,s=20)
b=cpo(sphere(20,[45,0,0],s=20))
fo(f'''
color("blue") for(p={a}) p_line3dc(p,.3);
color("magenta") for(p={b}) p_line3d(p,.3);
''')
    """
    return swapaxes(array(prism),0,1).tolist()

change_prism_points_orientation=cpo

def c2t3(p):# converts 2d list to 3d
    """
function to convert 2d to 3d, it just adds the z-coordinate to the points list 
example:
list=c2t3([[1,2],[3,4],[6,7]])
output=> [[1, 2, 0], [3, 4, 0], [6, 7, 0]]
    """
    return (array(p)@[[1,0,0],[0,1,0]]).tolist() if array(p).shape[-1]==2 else array(p).tolist()

convert2d_to_3d=c2t3
    
def c3t2(a): # converts 3d list to 2d list 
    """
function to convert 3d to 2d, it just removes the z-coordinate from the points list 
example:
list=c3t2([[1,2,3],[3,4,5],[6,7,8]])
output=> [[1, 2], [3, 4], [6, 7]]
    """
    return (array(a)@[[1,0],[0,1],[0,0]]).tolist() if array(a).shape[-1]==3 else array(a).tolist()

convert3d_to_2d=c3t2

def nv(p):# normal vector to the plane 'p' with atleast 3 known points
    """
given 3 points ['p1','p2',p3] function calculates unit normal vector
example:
p1,p2,p3=[1,0,0],[0,10,0],[-5,0,0]
nv([p1,p2,p3]) => [0.0, 0.0, -1.0]
    """
    l1=len(p)
    p0,p1,p2=array(translate([0,0,0],[p[0],p[int(l1/3)],p[int(l1*2/3)]]))
    nv=cross(p0-p1,p2-p1)
    m=1/norm(nv) if norm(nv)>0 else 1e5
    return (nv*m).tolist()
    
def nv1(p):# normal vector to the plane 'p' with atleast 3 known points
    """
given 3 points ['p1','p2',p3] function calculates normal vector
example:
p1,p2,p3=[1,0,0],[0,10,0],[-5,0,0]
nv1([p1,p2,p3]) => [0.0, 0.0, -60.0]
    """
    l1=len(p)
    p0,p1,p2=array(translate([0,0,0],[p[0],p[int(l1/3)],p[int(l1*2/3)]]))
    nv=cross(p0-p1,p2-p1)
    return nv.tolist()

def fillet_3p_3d(p0,p1,p2,r,s):# fillet with 3 known points 'p0,p1,p2' in 3d space. 'r' is the radius of fillet and 's' is the number of segments in the fillet
    """
function to create fillet given 3 points 'p1','p2','p3' 
r: radius of the fillet
s: number of segments in the fillet
example:
p1,p2,p3=[[3,0,0],[0,0,0],[0,3,3]]
fillet=fillet_3p_3d(p1,p2,p3,3,10)
fo(f'''
color("blue")p_line3dc({fillet},.05);
color("cyan")points({[p1,p2,p3]},.2);
    ''')
    """
    p0,p1,p2=array(translate([0,0,0],[p0,p1,p2]))
    n=array(nv([p0,p1,p2]))
    u1=(p0-p1)/(norm(p0-p1)+.00001)
    u2=(p2-p1)/(norm(p2-p1)+.00001)
    theta=(180-arccos(u1@u2)*180/pi)/2
    alpha=arccos(u1@u2)*180/pi
    l=r*tan(theta*pi/180)
    cp=p1+axis_rot(n,u1*r/cos(theta*pi/180),alpha/2)
    pa=p1+u1*l
    arc=[ cp+axis_rot(n,pa-cp,-i) for i in linspace(0,theta*2,s)]
    a,b,c=arc[0],arc[1:s-1],arc[s-1]
    return concatenate([[p1],arc]).tolist()

def fillet_3p_3d_cp(p0,p1,p2,r):# center point 'cp' of the fillet with 3 known points 'p0,p1,p2' in 3d space. 'r' is the radius of fillet
    """
function to find the center point of the fillet created by given 3 points 'p1','p2','p3' 
r: radius of the fillet
use function cp_cir_3d instead of this
    """
    p0,p1,p2=array(translate([0,0,0],[p0,p1,p2]))
    n=array(nv([p0,p1,p2]))
    u1=(p0-p1)/(norm(p0-p1)+.00001)
    u2=(p2-p1)/(norm(p2-p1)+.00001)
    theta=(180-arccos(u1@u2)*180/pi)/2
    alpha=arccos(u1@u2)*180/pi
    l=r*tan(theta*pi/180)
    cp=p1+axis_rot(n,u1*r/cos(theta*pi/180),alpha/2)
    return cp.tolist()
    

def arc_3p_3d(points,s=20):
    """
draws an arc through the 3 points list
's' is the number of segments of the circle
example:
p1,p2,p3=[[3,0,0],[0,1,0],[5,3,2]]
arc1=arc_3p_3d([p1,p2,p3],50)
fo(f'''
color("magenta")p_line3d({arc1},.05);
color("blue")points({[p1,p2,p3]},.2);
    ''')
    """
    n1=array(nv(points))+[.000001,.000001,0]
    a1=cross(n1,[0,0,-1])
    t1=r2d(arccos(n1@[0,0,-1]))
    sec1=translate(-array(points).mean(0),points)
    sec2=c3t2(axis_rot(a1,sec1,t1))
    l1=len(sec2)
    p0,p1,p2=[sec2[0],sec2[int(l1/3)],sec2[int(l1*2/3)]]
    arc1=arc_3p(p0,p1,p2,s=s)
    arc1=translate(array(points).mean(0),axis_rot(a1,arc1,-t1))
    return arc1

arc_3points_3d=arc_3p_3d

def r_3p_3d(points):
    """
calculates the radius of circle made by 3 points in 3d space
example:
p1,p2,p3=[[3,0,0],[0,1,0],[5,3,2]]
arc1=arc_3p_3d([p1,p2,p3],50)
r_3p_3d(arc1) => 2.95147084082169
    """
    n1=array(nv(points))+[.000001,.000001,0]
    a1=cross(n1,[0,0,-1])
    t1=r2d(arccos(n1@[0,0,-1]))
    sec1=translate(-array(points).mean(0),points)
    sec2=c3t2(axis_rot(a1,sec1,t1))
    l1=len(sec2)
    p0,p1,p2=[sec2[0],sec2[int(l1/3)],sec2[int(l1*2/3)]]
    
    return l_(r_3p([p0,p1,p2]))

radius_arc_3p_3d=r_3p_3d

def scl2d(sec,sl):# scale the 2d section 'sec' by a scaling factor 'sl'. this places the scaled section in the bottom center of the original section
    """
function to scale a 2d section by an amount "sl" which has to be >0 (keeps the y-coordinates same). 
following code scales the section by 0.7 (70% of the original shape)
example:
sec=corner_radius([[0,0,.5],[10,0,2],[7,15,1]],5)
sec1=scl2d(sec,.7)
fo(f'''
color("blue") for(p={[sec]}) p_line3dc(p,.3);
color("magenta") for(p={[sec1]}) p_line3dc(p,.3);
''')
    """
    s1=array(translate([0,0,0],sec))
    cp=array(s1).mean(axis=0)
    rev=array(s1).mean(axis=0)+(array(s1)-array(s1).mean(axis=0))*sl
    y1=cp-array([0,array(s1)[:,1].min(),0])
    y2=cp-array([0,rev[:,1].min(),0])
    d=y2-y1
    return c3t2(translate(d,rev))

scale2d=scl2d
    
def scl2d_c(sec,sl):# scale the 2d section 'sec' with scaling factor 'sl'. this places the scaled section in the center of original section or the center of both original and scaled section remains the same.
    """
function to scale a 2d section by an amount "sl" which has to be >0 (keeps the revised section in center). 
following code scales the section by 0.7 (70% of the original shape)
example:
sec=corner_radius([[0,0,.5],[10,0,2],[7,15,1]],5)
sec1=scl2d_c(sec,.7)
fo(f'''
color("blue") for(p={[sec]}) p_line3dc(p,.3);
color("magenta") for(p={[sec1]}) p_line3dc(p,.3);
''')
    """
    s1=array(translate([0,0,0],sec))
    cp=array(s1).mean(axis=0)
    rev=array(s1).mean(axis=0)+(array(s1)-array(s1).mean(axis=0))*sl
    return c3t2(rev)

scale2d_centered=scl2d_c

def scl3d(p,s):# scale 3d prism 'p' with scaling factor 's'. This places the scaled prism at the same bottom of the original prism
    """
function to scale a 3d prism keeping the base z-coordinate same. 
takes 2 arguments "p" to scale and the scaling factor "s". 
scale factor can take any real number negative values will scale the prism and turn the prism upside down.
try the following code to understand better:
example:
sec=circle(10);
path=corner_radius(pts1([[2,0],[-2,0,2],[0,10,3],[-3,0]]),5)
sol=prism(sec,path)
sol1=scl3d(sol,.7)
fo(f'''
%{swp(sol)}
{swp(sol1)}
''')
    """
    p=array(p)
    cp=p.reshape(-1,3).mean(axis=0)
    rev=cp+(p-cp)*s
    z1=p.reshape(-1,3)[:,2].min()
    z2=rev.reshape(-1,3)[:,2].min()
    d=z1-z2
    return translate([0,0,d],rev)
    
scale3d=scl3d

def scl3dc(p,s):# scale a 3d prism 'p' with scaling factor 's'. This places the scaled prism in the center of the original prism or the center of both the prism is same
    """
function to scale a 3d prism keeping the prism centered. takes 2 arguments "p" to scale and 
 the scaling factor "s". 
scale factor can take any real number negative values will scale the prism and turn the prism upside down.
try the following code to understand better:
example:
sec=circle(10);
path=corner_radius(pts1([[2,0],[-2,0,2],[0,10,3],[-3,0]]),5)
sol=prism(sec,path)
sol1=scl3dc(sol,.7)
fo(f'''
%{swp(sol)}
{swp(sol1)}
''')
    """
    p=array(p)
    cp=p.reshape(-1,3).mean(axis=0)
    rev=cp+(p-cp)*s
    return rev.tolist()

scale3d_centered=scl3dc
    
def m_points1(sec,s,d=.25):# multiple points with in the straight lines in the closed section 'sec'. 's' is the number of segments between each straight line
    """
adds 's' number of points in each straight line segment of a section 'sec'
'd' is the minimum segment length where multipe points to be added
example:
s1=square(10);
s2=translate_2d([12,0,0],m_points1(s1,10))
fo(f'''
color("blue",.2) p_line3dc({s1},.3);
color("magenta") points({s1},.3);

color("blue",.2) p_line3dc({s2},.3);
color("magenta") points({s2},.3);
''')
    """
    c=[]
    for i in range(len(sec)):
        i_plus=i+1 if i<len(sec)-1 else 0
        if l_len([sec[i],sec[i_plus]])/s>=d:
            c.append(ls([sec[i],sec[i_plus]],s))
        else:
            c.append([sec[i],sec[i_plus]])
    return remove_extra_points(concatenate(c))

divide_each_segment_of_closed_loop=m_points1
    

def m_points1_o(sec,s,d=.25):# multiple points with in the straight lines in the open section 'sec'. 's' is the number of segments between each straight line
    """
adds 's' number of points in each straight line segment of an open loop section 'sec'
'd' is the minimum segment length where multipe points to be added
example:
s1=square(10);
s2=translate_2d([12,0,0],m_points1_o(s1,10))
fo(f'''
color("blue",.2) p_line3d({s1},.3);
color("magenta") points({s1},.3);

color("blue",.2) p_line3d({s2},.3);
color("magenta") points({s2},.3);
''')
    """
    c=[]
    for i in range(len(sec)-1):
        i_plus=i+1 if i<len(sec)-1 else 0
        if l_len([sec[i],sec[i_plus]])/s>=d:
            c.append(ls([sec[i],sec[i_plus]],s))
        else:
            c.append([sec[i],sec[i_plus]])
    return remove_extra_points(concatenate(c))+[sec[-1]]

divide_each_segment_of_open_loop=m_points1_o


def ibsap(sec,pnt):# intersection between section and a point. used to find whether the poin is inside the section or outside the section
    p0=array(pnt)
    p2=sec
    p3=sec[1:]+[sec[0]]
    p2,p3=array([p2,p3])
    v1=[1,0]
    v2=(p3-p2)+[0,.00001]
    im=pinv(array([[v1]*len(v2),-v2]).transpose(1,0,2).transpose(0,2,1))
    p=p2-p0
    t1=einsum('ijk,ik->ij',im,p)[:,0]
    t2=einsum('ijk,ik->ij',im,p)[:,1]
    c1=(t2>=0)&(t2<=1)
    c2=t1>=0
    t=t1[c1&c2]
    p4=p0[None,:]+array(v1)*t.reshape(-1,1)
    return p4.tolist()



def circle(r,cp=[0,0],s=50): # circle with radius r and center point cp, s is the number of segments in the circle
    """
function for creating points in circle with radius "r", center point "cp" and number of segments "s" 
    """
    return array([ [cp[0]+r*cos(i*pi/180),cp[1]+r*sin(i*pi/180)] for i in linspace(0,360,s)][0:-1]).tolist()

    
def linear_extrude(sec,h=1,a=0,steps=1):
    """
function to linear extrude a section where
sec: section to extrude
h: height of the extrusion
a: angle of twist while extruding
steps: number of steps in each angular extrusion
example:
sec=corner_radius(pts1([[-5,-2.5,1],[10,0,1],[0,5,1],[-10,0,1]]),20)
sol=linear_extrude(sec,h=50,a=360,steps=200)
fo(f'''
color("magenta") for(p={sol})p_line3dc(p,.05,1);
{swp(sol)} 
''')
    """
    s=2 if a==0 else steps
    return [translate([0,0,h*i if a==0 else h/a*i],rot(f"z{0 if a==0 else i}",sec)) for i in linspace(0,1 if a==0 else a,s)]

def cylinder(r1=1,r2=1,h=1,s=50,r=0,d=0,d1=0,d2=0,center=False):
    """
function for making a cylinder
r1 or r: radius of circle at the bottom
r2 or r: radius of circle at the top
d1 or d: diameter of circle at the bottom
d2 or d: diameter of circle at the top
h: height of the cylinder
    """
    ra=r if r>0 else d/2 if d>0 else d1/2 if d1>0 else r1
    rb=r if r>0 else d/2 if d>0 else d2/2 if d2>0 else r2
    sec=c2t3(circle(ra,s=s))
    sec1=translate([0,0,h],circle(rb,s=s))
    sol=[sec,sec1]
    if center==True:
        return translate([0,0,-h/2],sol)
    else:
        return sol

def square(s=0,r=0,seg=10,center=False):
    """
function to draw a square
example:
s1=square(10);
s2=translate_2d([12,0],square([10,5]))
fo(f'''
color("blue",.2) p_line3dc({s1},.3);
color("blue",.2) p_line3dc({s2},.3);
''')
    """
    m= s if type(s)==int or type(s)==float else s[0]
    n= s if type(s)==int or type(s)==float else s[1]
    if r==0:
        sec=turtle2d([[0,0],[m,0],[0,n],[-m,0]])
    else:
        sec=cr2dt([[0,0,.01],[m,0,.01],[0,n,.01],[-m,0,.01]],seg)
    sec1= [[p[0]-m/2,p[1]-n/2] for p in sec] if center==True else sec
    return sec1

def rsz3d(prism,rsz):
    """
function to resize a 'prism' to dimensions 'rsz'
bottom left corner of both the prisms would be same
example:
sec=corner_radius(pts1([[0,0,1],[10,0,1],[0,5,1],[-10,0,1]]),5)
sol=translate([-5,0,0],linear_extrude(sec,10))
sol1=rsz3d(sol,[5,6,7])
fo(f'''
%{swp(sol)}
{swp(sol1)} 
''') 
    """
    prism1=array(prism).reshape(-1,3)
    max_x=prism1[:,0].max()
    max_y=prism1[:,1].max()
    max_z=prism1[:,2].max()
    min_x=prism1[:,0].min()
    min_y=prism1[:,1].min()
    min_z=prism1[:,2].min()
    avg=prism1.mean(axis=0)
    
    r_x=rsz[0]/(max_x-min_x)
    r_y=rsz[1]/(max_y-min_y)
    r_z=rsz[2]/(max_z-min_z)
    
    rev_prism=[[[avg[0]+r_x*(p[0]-avg[0]),avg[1]+r_y*(p[1]-avg[1]),avg[2]+r_z*(p[2]-avg[2])] for p in prism[i]] 
               for i in range(len(prism))]
    t=((array(bb(rev_prism))-array(bb(prism)))/2).tolist()
    return l_(translate(t,rev_prism))

resize_solid=rsz3d

def rsz3dc(prism,rsz):
    """
function to resize a 'prism' to dimensions 'rsz'
resized prism will be placed in the center of the original prism or center point of both the prisms will be same
example:
sec=corner_radius(pts1([[0,0,1],[10,0,1],[0,5,1],[-10,0,1]]),5)
sol=translate([-5,0,0],linear_extrude(sec,10))
sol1=rsz3dc(sol,[5,6,7])
fo(f'''
%{swp(sol)}
{swp(sol1)} 
''') 
    """
    prism1=array(prism).reshape(-1,3)
    max_x=prism1[:,0].max()
    max_y=prism1[:,1].max()
    max_z=prism1[:,2].max()
    min_x=prism1[:,0].min()
    min_y=prism1[:,1].min()
    min_z=prism1[:,2].min()
    avg=prism1.mean(axis=0)
    
    r_x=rsz[0]/(max_x-min_x)
    r_y=rsz[1]/(max_y-min_y)
    r_z=rsz[2]/(max_z-min_z)
    
    rev_prism=[[[avg[0]+r_x*(p[0]-avg[0]),avg[1]+r_y*(p[1]-avg[1]),avg[2]+r_z*(p[2]-avg[2])] for p in prism[i]] 
               for i in range(len(prism))]
    return l_(rev_prism)

resize_solid_centered=rsz3dc

def bb(prism):
    """
function to find the bounding box dimensions of a prism
example:
s1=square(10);
sol=linear_extrude(s1,10)
bb(sol) => [10.0, 10.0, 10.0]
    """
    prism1=array(prism).reshape(-1,3)
    max_x=prism1[:,0].max()
    max_y=prism1[:,1].max()
    max_z=prism1[:,2].max()
    min_x=prism1[:,0].min()
    min_y=prism1[:,1].min()
    min_z=prism1[:,2].min()
    return l_([max_x-min_x,max_y-min_y,max_z-min_z])

bounding_box_3d=bb

def cube(size=1,center=False):
    """
function to draw cube with size 'size'
example:
sol1=cube(10)
sol2=translate([12,0,0],cube([10,5,3]))
fo(f'''
{swp(sol1)}
{swp(sol2)}
''')
    """
    if type(size)==list:
        i,j,k=size[0],size[1],size[2]
        sec=pts([[0,0],[i,0],[0,j],[-i,0]])
        sol=linear_extrude(sec,k)
        if center==1:
            return translate([-i/2,-j/2,-k/2],sol)
        else:
            return sol
    elif type(size)==int:
        i,j,k=size,size,size
        sec=pts([[0,0],[i,0],[0,j],[-i,0]])
        sol=linear_extrude(sec,k)
        if center==1:
            return translate([-i/2,-j/2,-k/2],sol)
        else:
            return sol


def sphere(r=0,cp=[0,0,0],s=50):
   """
function to draw sphere with radius 'r' , center point 'cp' and number of segments 's'
example:
sol1=sphere(10)
fo(f'''
{swp(sol1)}
''')
   """
   path=arc(r,-90,90,s=int(s/2))
   p=[ translate([cp[0],cp[1],p[1]+cp[2]],circle(p[0],s=s)) for p in path]
   return array(p).tolist()


def rsz2d(sec,rsz):
    """
function to resize a 2d section to dimensions 'rsz'
resized section will be placed on bottom center of the original section
example:
s1=square(10)
s2=rsz2d(s1,[5,5])
fo(f'''
color("blue",.2) p_line3dc({s1},.3);
color("magenta",.2) p_line3dc({s2},.3);
''')
    """
    avg=array(sec).mean(axis=0)
    max_x=array(sec)[:,0].max()
    min_x=array(sec)[:,0].min()
    max_y=array(sec)[:,1].max()
    min_y=array(sec)[:,1].min()
    r_x=rsz[0]/(max_x-min_x)
    r_y=rsz[1]/(max_y-min_y)
    s=array([ avg+array([r_x*(sec[i][0]-avg[0]),r_y*(sec[i][1]-avg[1])-((min_y-avg[1])*r_y-(min_y-avg[1]))]) for i in range(len(sec))]).round(4)
    return s[sort(unique(s,axis=0,return_index=True)[1])].tolist()

resize2d=rsz2d

def rsz2dc(sec,rsz):
    """
function to resize a 2d section to dimensions 'rsz'
resized section will be placed in center of the original section
s1=square(10)
s2=rsz2dc(s1,[5,5])
fo(f'''
color("blue",.2) p_line3dc({s1},.3);
color("magenta",.2) p_line3dc({s2},.3);
''')
    """
    avg=array(sec).mean(axis=0)
    max_x=array(sec)[:,0].max()
    min_x=array(sec)[:,0].min()
    max_y=array(sec)[:,1].max()
    min_y=array(sec)[:,1].min()
    r_x=rsz[0]/(max_x-min_x)
    r_y=rsz[1]/(max_y-min_y)
    s=array([ avg+array([r_x*(sec[i][0]-avg[0]),r_y*(sec[i][1]-avg[1])]) for i in range(len(sec))]).round(4)
    return s[sort(unique(s,axis=0,return_index=True)[1])].tolist()

resize2d_centered=rsz2dc


def ip(sol1,sol2):
    """
function to calculate intersection point between two 3d prisms. 
"sol1" is the 3d object which is intersected with "sol2".
use function ip_sol2sol instead of this
    """
    line=array([ seg(p)[:-1] for p in cpo(sol2)])
    v,f1=vnf2(sol1)
    tri=array(v)[array(f1)]
    line=array([ seg(p)[:-1] for p in cpo(sol2)])
    tri.shape,line.shape
    la,lb=line[:,:,0],line[:,:,1]
    p0,p1,p2=tri[:,0],tri[:,1],tri[:,2]
    lab=lb-la
    p01,p02=p1-p0,p2-p0
    t=einsum('kl,ijkl->ijk',cross(p01,p02),la[:,:,None]-p0)/(einsum('ijl,kl->ijk',(-lab),cross(p01,p02))+.00001)
    u=einsum('ijkl,ijkl->ijk',cross(p02[None,None,:,:],(-lab)[:,:,None,:]),(la[:,:,None,:]-p0[None,None,:,:]))/(einsum('ijl,kl->ijk',(-lab),cross(p01,p02))+.00001)
    v=einsum('ijkl,ijkl->ijk',cross((-lab)[:,:,None,:],p01[None,None,:,:]),(la[:,:,None,:]-p0[None,None,:,:]))/(einsum('ijl,kl->ijk',(-lab),cross(p01,p02))+.00001)
    condition=(t>=0)&(t<=1)&(u>=-0.0001)&(u<=1)&(v>=-.0001)&(v<=1)&(u+v<1)
    i_p=(array([la]*len(p0)).transpose(1,2,0,3)+einsum('ijl,ijk->ijkl',lab,t))[condition].tolist()
    return i_p



def s_int1(sec1):
    """
calulates the self intersection points of a list of line segments 's'
it picks the intersection points only if the 2 lines are crossing each other
example:
sec=[[0,0],[10,0],[10,15]]
sec1=offset_segv(sec,-2)
p0=s_int1(sec1)
fo(f'''
color("blue",.2) for(p={[sec]}) p_line3dc(p,.3);
color("cyan",.2) for(p={sec1}) p_line3dc(p,.3);
color("magenta") points({p0},.5);
''')
    """
    n=len(sec1)
    a=array(sec1)[comb_list(n)]
    p0=a[:,0][:,0]
    p1=a[:,0][:,1]
    p2=a[:,1][:,0]
    p3=a[:,1][:,1]
    v1=a_(rot2d(0.00001,p1-p0))
    v2=p3-p2
    iim=array([v1,-v2]).transpose(1,0,2).transpose(0,2,1)
    im=inv(iim)
    p=p2-p0

    t=einsum('ijk,ik->ij',im,p)
    dcn=(t[:,0].round(4)>0)&(t[:,0].round(4)<1)&(t[:,1].round(4)>0)&(t[:,1].round(4)<1)
    i_p1=p0+einsum('ij,i->ij',v1,t[:,0])
    i_p1=i_p1[dcn].tolist()
    return i_p1

    
def self_intersections(sec1):
    """
self intersections based on Bentley-Ottmann line sweep
avoid using this, use s_int1 instead
    """
    sec2=lexicographic_seg_sort_xy(sec1)
    s,s1,b=[],[],[]
    for i in range(len(sec2)):
        s=s+[sec2[i]]
        s1.append(s_int1(s))
        if i>1:
            for j in range(len(s)):
                if s[-1][0][0]>=s[j][1][0]:
                    b=exclude_seg(s,[s[j]])
            s=b if b!=[] else s
            b=[]

    s1=[p for p in s1 if p!=[]]
    s1=concatenate(s1).tolist()
    return s1

def comb(n,i): 
    """
calculates number of possible combinations for "n" items with "i" selected items
example:
comb(8,2) => 28
    """
    return int(math.factorial(n)/(math.factorial(i)*math.factorial(n-i)))

combinations=comb

def bezier(pl,s=20):
    """
bezier curve defined by points 'pl' and number of segments 's'
example:
sec=[[0,0],[10,0],[10,15],[20,10]]
sec1=bezier(sec,20)
fo(f'''
color("blue",.2) for(p={[sec]}) p_line3d(p,.3);
color("cyan",.2) for(p={[sec1]}) p_line3d(p,.3);
''')
    """
    def B(n,i,t):
        """
        Bernstein function
        """
        if i>n or i<0:
            return 0
        else:
            return comb(n,i)*t**i*(1-t)**(n-i)
            
    b=[l_(a_([a_(pl[i])*B(len(pl)-1,i,t) for i in range(len(pl))]).sum(0)) 
       for t in linspace(0,1,s)]
    return b

def arc_3d(v=[0,0,1],r=1,theta1=0,theta2=360,cw=-1,s=50):
    """
3d arc defined by normal vector 'v', radius 'r1', start angle 'theta1', 
end angle 'theta2' , clockwise(1) or counter clockwise(-1) and number of segments 's'
example:
v1=[1,2,2]
a1=arc_3d(v=v1,r=5,theta1=0,theta2=270,cw=-1,s=50)
fo(f'''
color("blue",.2) for(p={[a1]}) p_line3d(p,.3);
color("magenta",.2) p_line3d({point_vector([0,0,0],v1)},.3);
''')
    """

    if uv(v)==[0,0,1]:
        arc1=arc(r,theta1,theta2,[0,0],s) if cw==-1 else flip(arc(r,theta1,theta2,[0,0],s))
        return c2t3(arc1)
    elif uv(v)==[0,0,-1]:
        arc1=arc(r,theta1,theta2,[0,0],s) if cw==-1 else flip(arc(r,theta1,theta2,[0,0],s))
        arc1=rot('y180',arc1)
        return arc1
    else:
        sec=arc(r,theta1,theta2,[0,0],s) if cw==-1 else flip(arc(r,theta1,theta2,[0,0],s))
        s=rot('x90z-90',sec)
        v1=array(v)+array([0,0,0.00001])
        va=[v1[0],v1[1],0]
        u1=array(uv(v1))
        ua=array(uv(va))
        v2=cross(va,v1)
        a1=arccos(u1@ua)*180/pi
        a2=ang(v1[0],v1[1])
        s1=rot(f'z{a2}',s)
        sec1=[axis_rot(v2,p,a1) for p in s1]
        return sec1



def l_cir_ip(line,cir):
    """
line circle intersection point
example:
c1=circle(10)
l1=point_vector([-5,-5],[3,2])
p0=l_cir_ip(l1,c1)
fo(f'''
color("blue",.2) for(p={[c1,l1]}) p_line3d(p,.3);
color("magenta") points({p0},.5);
''')
    """
    
    p0=line[0]
    p1=line[1]
    p0,p1=array(p0),array(p1)
    v1=p1-p0
    cp=cp_3p(cir[0],cir[int(len(cir)/2)],cir[int(len(cir)*2/3)])
    r=r_3p([cir[0],cir[int(len(cir)/2)],cir[int(len(cir)*2/3)]])

    a=v1[0]**2+v1[1]**2
    b=2*p0[0]*v1[0]-2*v1[0]*cp[0]+2*p0[1]*v1[1]-2*v1[1]*cp[1]
    c=p0[0]**2+p0[1]**2+cp[0]**2+cp[1]**2-2*p0[0]*cp[0]-2*p0[1]*cp[1] -r**2
    t1=(-b-sqrt(b**2-4*a*c))/(2*a)
    t2=(-b+sqrt(b**2-4*a*c))/(2*a)
    p2=p0+v1*t1
    p3=p0+v1*t2
    return [p2.tolist(),p3.tolist()]

line_circle_ip=l_cir_ip
    
def s_pnt(pnt): # starting point for calculating convex hull (bottom left point)
    """
starting point for calculating convex hull (bottom left point)
    """
    pnt=array(pnt)
    c1=pnt[:,1]==pnt[:,1].min()
    s1=pnt[c1]
    c2=s1[:,0]==s1[:,0].min()
    return s1[c2][0].tolist()

def n_pnt(pnt,sp,an):
    pnt,sp=array(pnt),array(sp)
    pnt=pnt[(pnt!=sp).all(1)]
    a=pnt-sp
    a1=vectorize(ang)(a[:,0],a[:,1])
    n_pnt=pnt[a1==a1[a1>=an].min()][0].tolist()
    return [n_pnt,a1[a1>=an].min().tolist()]

def c_hull(pnt): # convex hull for an array of points
    """
function to calculate convex hull for a list of points 'pnt'
    """
    c=[]
    np=n_pnt(pnt,s_pnt(pnt),0)
    for i in range(len(pnt)):
        c.append(np[0])
        np=n_pnt(pnt,np[0],np[1])
        if np[0]==s_pnt(pnt):
            break
    return [s_pnt(pnt)]+c

def convex(sec):
    """
function to check whether a section is convex or not
example:
sec1=corner_radius(pts1([[0,0,.2],[8,3,3],[5,7,1],[-8,0,2],[-5,20,1]]),20)
sec2=corner_radius(pts1([[0,0,.1],[7,5,2],[5,7,3],[-5,7,5],[-7,5,5]]),20)
convex(sec1),convex(sec2) => (False, True)

refer file "example of various functions" for application example
    """
    sec= c3t2(sec) if array(sec).shape[-1]==3 else sec
    return (array(cwv(sec))==-1).all()|(array(cwv(sec))==1).all()

def oo_convex(sec,r): #outer offset of a convex section
    s=flip(sec) if cw(sec)==1 else sec
    return offset_points(sec,r)

def cir_p_t(cir,p):
    """
    circle to point tangent line (point should be outside the circle)
example:
c1=circle(10)
p0=[15,15]
p1=cir_p_t(c1,p0)
fo(f'''
color("blue",.2) for(p={[c1]}) p_line3d(p,.3);
color("magenta") points({[p0,p1]},.5);
color("cyan",.2) p_line3d({[p1,p0]},.3);
''')
    """
    cp=cp_3p(cir[0],cir[int(len(cir)/3)],cir[int(len(cir)*2/3)])
    r=r_3p([cir[0],cir[int(len(cir)/3)],cir[int(len(cir)*2/3)]])
    if array(r).round(4)==norm(array(p)-array(cp)).round(4):
        tp=array(p).tolist()
    else:
        l1=l_len([p,cp])
        v1=array(p)-array(cp)
        theta1=ang(v1[0],v1[1])
        theta2=arccos(r/l1)*180/pi
        theta3=(theta1-theta2)*pi/180
        tp=array([r*cos(theta3),r*sin(theta3)])+array(cp)
        tp=tp.tolist()
    return tp

circle_to_point_tangent=cir_p_t

def p_cir_t(p,cir): # point to circle tangent line (point should be outside the circle)
    """
point to circle tangent line (point should be outside the circle)
example:
c1=circle(10)
p0=[15,15]
p1=p_cir_t(p0,c1)
fo(f'''
color("blue",.2) for(p={[c1]}) p_line3d(p,.3);
color("magenta") points({[p0,p1]},.5);
color("cyan",.2) p_line3d({[p0,p1]},.3);
''')
    """
    cp=cp_3p(cir[0],cir[int(len(cir)/3)],cir[int(len(cir)*2/3)])
    r=r_3p([cir[0],cir[int(len(cir)/3)],cir[int(len(cir)*2/3)]])
    if array(r).round(4)==norm(array(p)-array(cp)).round(4):
        tp=array(p).tolist()
    else:
        l1=l_len([p,cp])
        v1=array(p)-array(cp)
        theta1=ang(v1[0],v1[1])
        theta2=arccos(r/l1)*180/pi
        theta3=(theta1+theta2)*pi/180
        tp=array([r*cos(theta3),r*sin(theta3)])+array(cp)
        tp=tp.tolist()
    return tp

point_to_circle_tangent=p_cir_t

def v_sec_extrude(sec,path,o):
    """
extrude a section 'sec' through a path 'path' 
section will vary from start to end such that at the end the section will be offset by 'o' distance
example:
sec=corner_radius(pts1([[0,0,1],[10,0,1],[0,5,1],[-10,0,1]]),10)
path=helix(20,15,1,5)
sol=v_sec_extrude(sec,path,-2)
fo(f'''
{swp(sol)}
    ''')
    """
    sec0=[offset(sec,i) for i in linspace(0,o,len(path))]
    p1=path[:-1]
    p2=path[1:]
    p1,p2=array([p1,p2])
    v1=p2-p1
    u1=v1/norm(v1,axis=1).reshape(-1,1)
    v2=concatenate([[u1[0]],(u1[1:]+u1[:-1])/2,[u1[-1]]])
    sec2=[]
    for i in range(len(path)):
        sec1=translate(path[i],sec2vector(v2[i],sec0[i]))
        sec2.append(sec1)
    return sec2


def two_cir_tarc(c1,c2,r,side=0,s=50): #two circle tangent arc
    """
function draws a arc which is tangent to 2 circles 'c1' and 'c2'    's' is the number of segments of the tangent arc
'r' is the radius of the tangent arc 
'side' there are 2 sides of the circles where the arc could be created defined by '0' and '1'
example:
c1=circle(10,[20,10])
c2=circle(10)
a1=two_cir_tarc(c1,c2,5)
a2=two_cir_tarc(c2,c1,5)
fo(f'''
color("blue",.2) for(p={[c1,c2]}) p_line3d(p,.3);
color("magenta") for(p={[a1]}) p_line3d(p,.3);
color("cyan") for(p={[a2]}) p_line3d(p,.3);
''')
    """
    def filleto_2cir(r1,r2,c1,c2,r,s=50): # fillet between 2 circles with radius 'r1' and 'r2' and center points 'c1' and 'c2' and 'r' is the radius of the fillet. This is an open fillet where first or the second fillet can be called based on requirement
        """
function to draw the fillet radius "r" between the 2 circle with radiuses "r1" and "r2" centered at "c1" and "c2" respectively.
This function gives an additional flexibility for drawing fillet only one side. e.g 
fillet=filleto_2cir(r1=10,r2=10,c1=[0,0],c2=[20,0],r=10)
fillet[0] will calculate fillet on one side
refer to the file "example of various functions" to see the application
        """
        
        c1,c2=array([c1,c2])
        l1=norm(c2-c1)
        l2=r1+r
        l3=r2+r
        t=(l1**2+l2**2-l3**2)/(2*l1)
        h=sqrt(l2**2-t**2)
        v=c2-c1
        u=v/norm(v)
        p1=c1+u*t+(u@rm(90))*h
        a1=ang((c1-p1)[0],(c1-p1)[1])
        a2=ang((c2-p1)[0],(c2-p1)[1])
        p2=c1+u*t+u@rm(-90)*h
        a3=ang((c2-p2)[0],(c2-p2)[1])
        a4=ang((c1-p2)[0],(c1-p2)[1])
        a5=ang((p1-c1)[0],(p1-c1)[1])
        a6=ang((p2-c1)[0] ,(p2-c1)[1])
        a7=ang((p1-c2)[0] ,(p1-c2)[1])
        a8=ang((p2-c2)[0] ,(p2-c2)[1])
    
        arc1=arc(r,360+a2 if a2<a1 else a2,a1,p1,s=s)
        arc2=arc(r,360+a4 if a4<a3 else a4,a3,p2,s=s)
        arc3=arc(r2,360+a7 if a7<a8 else a7,a8,c2,s=s)
        arc4=arc(r1,a5,360+a6 if a6<a5 else a6,c1,s=s)
    
        return [arc2,arc1]


    
    r1,r2,cp1,cp2=r_arc(c1),r_arc(c2),cp_arc(c1),cp_arc(c2)
    cp1,cp2=array([cp1,cp2])
    l1=norm(cp2-cp1)
    if r>(r1+r2+l1)/2:
        l2=r-r1
        l3=r-r2
        x=(l2**2-l3**2+l1**2)/(2*l1)
        h=sqrt(l2**2-x**2)
        v1=cp2-cp1
        u1=v1/norm(v1)
        p0=cp1+u1*x
        if side==1:
            cp3=p0-(u1@rm(90))*h
        elif side==0:
            cp3=p0+(u1@rm(90))*h

        v2=cp2-cp3
        u2=v2/norm(v2)
        v3=cp1-cp3
        u3=v3/norm(v3)
        p1=cp2+u2*r2
        p2=cp1+u3*r1

        if side==1:
            arc1=arc_2p(p1,p2,r,-1,s=s)
        elif side==0:
            arc1=arc_2p(p2,p1,r,-1,s=s)



    else:
        if side==1:
            arc1=filleto_2cir(r1,r2,cp1,cp2,r,s=s)[1]
        elif side==0:
            arc1=filleto_2cir(r1,r2,cp1,cp2,r,s=s)[0]
            
    return arc1

two_circles_tangent_arc=two_cir_tarc




def arc_3p(p1,p2,p3,s=30):
    """ 
function to draw arc with 3 known points 'p1','p2','p3' 
's' is the number of segments of the arc
example:
p0,p1,p2=[[0,0],[10,0],[3,7]]
a1=arc_3p(p0,p1,p2)
fo(f'''
color("blue",.2) for(p={[a1]}) p_line3d(p,.3);
color("magenta") points({[p0,p1,p2]},.5);
''')
    """
    p1,p2,p3=array([p1,p2,p3])
    p4=p1+(p2-p1)/2
    p5=p2+(p3-p2)/2
    v1=p2-p4
    u1=v1/norm(v1)
    v2=p3-p5
    u2=v2/norm(v2)
    p6=p4+u1@rm(90)
    p7=p5+u2@rm(90)
    cp=i_p2d([p4,p6],[p5,p7])
    r=norm(p1-cp)
    v3=p1-cp
    v4=p2-cp
    v5=p3-cp
    a1=ang(v3[0],v3[1])
    a2=ang(v4[0],v4[1])
    a3=ang(v5[0],v5[1])
    a4=(a3+360 if a3<a1 else a3) if cw([p1,p2,p3])==-1 else (a3 if a3<a1 else a3-360)
    return arc(r,a1,a4,cp,s)

arc_3points=arc_3p

def cir_3p(p1,p2,p3,s=30):
    """ 
function to draw circle with 3 known points 'p1','p2','p3' 
's' is the number of segments of the circle
example:
p0,p1,p2=[[0,0],[10,0],[3,7]]
a1=cir_3p(p0,p1,p2)
fo(f'''
color("blue",.2) for(p={[a1]}) p_line3d(p,.3);
color("magenta") points({[p0,p1,p2]},.5);
''')
    """
    p1,p2,p3=array([p1,p2,p3])
    p4=p1+(p2-p1)/2
    p5=p2+(p3-p2)/2
    v1=p2-p4
    u1=v1/norm(v1)
    v2=p3-p5
    u2=v2/norm(v2)
    p6=p4+u1@rm(90)
    p7=p5+u2@rm(90)
    cp=i_p2d([p4,p6],[p5,p7])
    r=norm(p1-cp)

    return circle(r,cp,s)

circle_with_3points=cir_3p

def cp_3p(p1,p2,p3):
    """
function to calculate center point of a circle created from 3 known points 'p1','p2','p3'
example:
p0,p1,p2=[[0,0],[10,0],[3,7]]
a1=cir_3p(p0,p1,p2)
c1=cp_3p(p0,p1,p2)
fo(f'''
color("blue",.2) for(p={[a1]}) p_line3d(p,.3);
color("magenta") points({[p0,p1,p2,c1]},.5);
''')
    """
    p1,p2,p3=array([p1,p2,p3])
    p4=p1+(p2-p1)/2
    p5=p2+(p3-p2)/2
    v1=p2-p4
    u1=v1/norm(v1)
    v2=p3-p5
    u2=v2/norm(v2)
    p6=p4+u1@rm(90)
    p7=p5+u2@rm(90)
    cp=i_p2d([p4,p6],[p5,p7])
    return array(cp).tolist()

center_circle_3points=cp_3p

def ip_surf(surf2,surf1):
    """
function to calculate intersection point between two 3d prisms or between surface and solid. 
 "surf2" is the 3d object which is intersected with "surf1".
try below code for better understanding:
example:
sec=circle(10)
path=corner_radius(pts1([[2,0],[-2,0,2],[0,10,3],[-9.9,0]]),5)
p1=prism(sec,path)
p2=rot("y40",cylinder(r=3,h=15,s=30))
fo(f'''
%swp({p1});
%swp({p2});
ip={ip_surf(c_(p1),p2)};
color("magenta")points(ip,.2);
''')
    """
    i,j,_=array(surf2).shape
    a=surf2
    b=surf1
    p1=array([[[[a[i][j],a[i+1][j],a[i][j+1]],[a[i+1][j+1],a[i][j+1],a[i+1][j]]] 
            for j in range(j-1)]  for i in range(i-1)]).reshape(-1,3,3)
    p2=array([[[b[i][j],b[i+1][j]] for j in range(len(b[i]))] for i in range(len(b)-1)]).reshape(-1,2,3)
    pm=p1[:,0]
    pn=p1[:,1]
    po=p1[:,2]
    px=p2[:,0]
    py=p2[:,1]
    v1,v2,v3=py-px,pn-pm,po-pm
    t1=einsum('ijk,jk->ij',px[:,None]-pm,cross(v2,v3))/(einsum('ik,jk->ij',-v1,cross(v2,v3))+.00001)
    t2=einsum('ijk,ijk->ij',px[:,None]-pm,cross(v3,-v1[:,None]))/(einsum('ik,jk->ij',-v1,cross(v2,v3))+.00001)
    t3=einsum('ijk,ijk->ij',px[:,None]-pm,cross(-v1[:,None],v2))/(einsum('ik,jk->ij',-v1,cross(v2,v3))+.00001)
    p=px[:,None]+einsum('ik,ij->ijk',v1,t1)
    condition=(t1>=0)&(t1<=1)&(t2>=0)&(t2<=1)&(t3>=0)&(t3<=1)&((t2+t3)>=0)&((t2+t3)<=1)
    return p[condition].tolist()

def perp(sec,point,radius):
    sec=array(seg(sec))
    p0=sec[:,0]
    p1=sec[:,1]
    v1=p1-p0
    u1=v1/(norm(v1,axis=1).reshape(-1,1)+.00001)
    v2=array(point)-p0
    v1norm=norm(v1,axis=1)
    v2norm=norm(v2,axis=1)
    v2cost=einsum('ij,ij->i',u1,v2)
    cond1=v2cost>=0
    cond2=v2cost<=v1norm
    d=sqrt(v2norm**2-v2cost**2)
    d=min(d[(cond1)&(cond2)]).round(4)
    cond3=d==round(abs(radius),3)
    return point if cond3 else []


def near_points(points,s_p,n):
    l=array([ norm(array(p)-array(s_p)) for p in points])
    l1=sort(l)[0:n+1]
    index=array([[i for i in range(len(l)) if p==l[i]]for p in l1]).reshape(-1)
    p1=array(points)[index].tolist()
    return p1[1:]

def next_point(points,s_p):
    a1=[270+(360-ang((array(p)-array(s_p))[0],(array(p)-array(s_p))[1]))
        if ang((array(p)-array(s_p))[0],(array(p)-array(s_p))[1])>270 else
        270-ang((array(p)-array(s_p))[0],(array(p)-array(s_p))[1])
        for p in points]
    n_p=array(points)[a1==max(a1)][0].tolist()
    return n_p


def exclude_points(list1,list_to_exclude):
    """
exclude a list of points from a defined list
example:
sec=circle(10)
path=corner_radius(pts1([[2,0],[-2,0,2],[0,10,3],[-9.9,0]]),5)
p1=prism(sec,path)
p2=rot("y40",cylinder(r=3,h=15,s=30))
ip1=ip_surf(c_(p1),p2)
ip2=ip1[:4] # points to exclude
ip3=exclude_points(ip1,ip2)
fo(f'''
//original points list 
color("magenta") points({ip1},.2);
//excluded points list
%color("blue",.2) points({ip3},.3);
''')
    """
    la,lb=array(list1).round(5),array(list_to_exclude).round(5)
    return a_(list1)[~(lb==la[:,None]).all(2).any(1)].tolist()
    
def exclude_seg(list,list_to_exclude):
    """
similar like exclude_points. This function excludes segments from a segment list
example:
a=cr2dt([[0,0],[10,0],[0,10],[-5,5],[-5,-5]])
b=round(seg(a),3)
c=round(seg(a)[2:4],3)
exclude_seg(b,c) => 
[[[0.0, 0.0], [10.0, 0.0]],
[[10.0, 0.0], [10.0, 10.0]],
[[-0.0, 10.0], [0.0, 0.0]]]
    """
    return array(list)[~ (array(list)==array(list_to_exclude)[:,None]).all(2).all(2).transpose(1,0).any(1)].tolist()

        
def rev_pnts(sec,pnts):
    s8,s4=[sec,pnts]
    p0=array(s4)
    p2=s8
    p3=s8[1:]+[s8[0]]
    p2,p3=array([p2,p3])
    v2=(p3-p2)
    con1=v2.round(3)[:,1]==0
    y_list=p2[con1][:,1].round(3)
    p0=p0[(p0[:,1][:,None]!=y_list).all(1)]
    return p0.tolist()
    
def pies1(section,pnts):
    """
function to find points which are inside the enclosed section
example:
a=random.random(10000)*(20-(-5))+(-5)
b=random.random(10000)*(25-(-2))+(-2)
points=array([a,b]).transpose(1,0).tolist()
sec=corner_radius(pts1([[2,1,.1],[7,5,2],[5,7,3],[-5,7,2],[-7,5,3]]),20)
# sec=corner_radius(pts1([[8,0],[11,10,10],[0,10,5],[-10,0,5],[-1,-6,0.3],[-1,6,5],[-10,0,5],[0,-10,10]]))
fo(f'''
color("blue") p_line({sec},.05);
points({points},.2);
color("magenta")points({pies1(sec,points)},.2);
    ''')
    """
    s1=section
    v1=[1,0.00001]
    v2=a_([line_as_vector(p) for p in seg(s1)])+[0,.000001]
    p0=a_(pnts)
    p1=a_(s1)
    # p0+v1*t1=p1+v2*t2
    # v1*t1-v2*t2=p1-p0
    b=[]
    for i in range(len(p0)):
        v3=a_([v1]*len(s1))
        iim=a_([v3,-v2]).transpose(1,0,2).transpose(0,2,1)
        im=inv(iim)
        p=(p1-p0[i][None,:])
        p.shape,im.shape
        t1,t2=einsum('ijk,ik->ij',im,p).transpose(1,0)
        dec=(t1>0)&(t2>0)&(t2<1)
        if l_(dec.sum()%2!=0):
            b.append(pnts[i])
    return b

def r_sec(r1,r2,cp1,cp2,s=20):
    """
creates a rounded section around a line defined by points 'cp1' and 'cp2'
radius around 'cp1' is 'r1' and radius around 'cp2' is 'r2'
example:
line=[[1,2],[5,-6]]
r1,r2=5,1
sec1=r_sec(r1,r2,line[0],line[1])
fo(f'''
color("blue")p_line3dc({sec1},.1);
color("magenta")p_line3d({line},.1);
color("cyan")points({line},.2);
''')
    """
    sec=tctpf(r2,r1,cp2,cp1)
    a1=arc_long_2p(sec[1],sec[2],r1,-1,s=s) if r1>r2 else arc_2p(sec[1],sec[2],r1,-1,s=s)
    a2=arc_long_2p(sec[3],sec[0],r2,-1,s=s) if r2>r1 else arc_2p(sec[3],sec[0],r2,-1,s=s)
    sec1=a1+a2
    return sec1


def cs1(sec,d):
    """
creates a cleaning section for removing excess points for offseting a section 'sec' with offset distance 'd'
example:
sec=corner_radius(pts1([[0,0,1],[8,3,3],[5,7,1],[-8,0,2],[-5,20,1]]),10)
fo(f'''
color("blue")p_line({sec},.1); 
color("magenta")for(p={cs1(sec,-2)})p_line(p,.08); 
''')
    """
    r=abs(d)
    a=seg(sec)
    cs=[r_sec(r-r/1000,r-r/1000,p2[0],p2[1]) for p2 in a ]
    return cs
    
def cs2(sec,d):
    r=abs(d)
    a=seg(sec)
    cs=[r_sec(r-r/1000,r-r/1000,p2[0],p2[1]) for p2 in a if l_len(p2)>.5]
    return cs


def swp(bead2):
    """
function to render various 3d shapes
example:
swp(cylinder(d=10,h=20)) will render a cylinder with dia 10 and height 20
    """
    n1=arange(len(bead2[0])).tolist()
    n2=array([[[[(j+1)+i*len(bead2[0]),j+i*len(bead2[0]),j+(i+1)*len(bead2[0])],[(j+1)+i*len(bead2[0]),j+(i+1)*len(bead2[0]),(j+1)+(i+1)*len(bead2[0])]] \
             if j<len(bead2[0])-1 else \
             [[0+i*len(bead2[0]),j+i*len(bead2[0]),j+(i+1)*len(bead2[0])],[0+i*len(bead2[0]),j+(i+1)*len(bead2[0]),0+(i+1)*len(bead2[0])]] \
                 for j in range(len(bead2[0]))] for i in range(len(bead2)-1)]).reshape(-1,3).tolist()
    n3=(array(flip(arange(len(bead2[0]))))+(len(bead2)-1)*len(bead2[0])).tolist()
    n=[n1]+n2+[n3]
    pnt=array(bead2).reshape(-1,3).round(4).tolist()
    return f'polyhedron({pnt},{n},convexity=10);'

poly_h=swp
    

def swp_c(bead2):
    """
function to render various polyhedron with closed loop shapes e.g. fillets
example:
sec=circle(5)
path=c23(circle(20))
sol=path_extrude_closed(sec,path)
fo(f'''
{swp_c(sol)}
''')
    """
    n1=arange(len(bead2[0])).tolist()
    n2=array([[[[(j+1)+i*len(bead2[0]),j+i*len(bead2[0]),j+(i+1)*len(bead2[0])],[(j+1)+i*len(bead2[0]),j+(i+1)*len(bead2[0]),(j+1)+(i+1)*len(bead2[0])]] \
             if j<len(bead2[0])-1 else \
             [[0+i*len(bead2[0]),j+i*len(bead2[0]),j+(i+1)*len(bead2[0])],[0+i*len(bead2[0]),j+(i+1)*len(bead2[0]),0+(i+1)*len(bead2[0])]] \
                 for j in range(len(bead2[0]))] for i in range(len(bead2)-1)]).reshape(-1,3).tolist()
    n3=(array(flip(arange(len(bead2[0]))))+(len(bead2)-1)*len(bead2[0])).tolist()
    n=[n1]+n2+[n3]
    pnt=array(bead2).reshape(-1,3).round(4).tolist()
    return f'polyhedron({pnt},{n2},convexity=10);'

poly_h_closed_loop=swp_c
    
def mesh_vf(bead2):
    """
finds the vertices and faces of a solid
example:
sec=square(10)
sol=linear_extrude(sec,10)
v,f=mesh_vf(sol) # 'v' are vertices and 'f' are faces
v,f =>
(array([[ 0.,  0.,  0.],
        [10.,  0.,  0.],
        [10., 10.,  0.],
        [ 0., 10.,  0.],
        [ 0.,  0., 10.],
        [10.,  0., 10.],
        [10., 10., 10.],
        [ 0., 10., 10.]]),
 array([[1, 0, 4],
        [1, 4, 5],
        [2, 1, 5],
        [2, 5, 6],
        [3, 2, 6],
        [3, 6, 7],
        [0, 3, 7],
        [0, 7, 4]]))

    """
    n1=arange(len(bead2[0])).tolist()
    n2=array([[[[(j+1)+i*len(bead2[0]),j+i*len(bead2[0]),j+(i+1)*len(bead2[0])],[(j+1)+i*len(bead2[0]),j+(i+1)*len(bead2[0]),(j+1)+(i+1)*len(bead2[0])]] \
             if j<len(bead2[0])-1 else \
             [[0+i*len(bead2[0]),j+i*len(bead2[0]),j+(i+1)*len(bead2[0])],[0+i*len(bead2[0]),j+(i+1)*len(bead2[0]),0+(i+1)*len(bead2[0])]] \
                 for j in range(len(bead2[0]))] for i in range(len(bead2)-1)]).reshape(-1,3)
    pnt=array(bead2).reshape(-1,3).round(4)
    return [pnt,n2]


def swp_prism_h(prism_big,prism_small):
    """
creats a hollow prism with 2 similar prisms (1 big and 1 smaller)
example:
sec=square(10)
sol=linear_extrude(sec,10)
sol1=offset_sol(sol,-1)
sol2=swp_prism_h(sol,sol1)
fo(f'''
{swp_c(sol2)}
''')
    """
    
    p1=prism_big
    p2=flip(prism_small)
    p3=p1+p2+[p1[0]]
    return p3
    
create_tubular_solids=swp_prism_h
    

def surf_base(surf,h=0):
    """
creates a solid from any surface, 'h' is the height of the base of the surface
example:
sec2=corner_radius(pts1([[-25,0],[10,5,5],[10,-3,10],[10,5,5],[10,-8,7],[10,1]]),10)  
path2=cytz(corner_radius(pts1([[-35,5,0],[10,8,20],[20,-5,10],[20,8,20],[10,-9,20],[10,1,0]]),10))
surf2=path_extrude_open(sec2,path2)
sol=surf_base(surf2,0)
fo(f'''
{swp(sol)}
''') 
    """
    s=cpo(surf)
    s1=translate([0,0,h],c2t3(c3t2([flip(p) for p in s])))
    s2=array([s,s1]).transpose(1,0,2,3)
    
    i,j,k,l=s2.shape
    s2=s2.reshape(i,j*k,l).tolist()
    t=array(surf).reshape(-1,3).mean(0)[2]
    return s2 if h>t else flip(s2)


def helix(radius=10,pitch=10, number_of_coils=1, step_angle=1):
    """
creates a helix with radius, pitch and number of coils parameters
example:
a=helix(radius=10,pitch=10,number_of_coils=5,step_angle=5)
fo(f'''
color("blue") p_line3d({a},.5);
''')
    """
    return l_(a_([[radius*cos(d2r(i)),radius*sin(d2r(i)),i/360*pitch] for i in arange(0,360*number_of_coils,step_angle)]))

def multiple_sec_extrude(path_points=[],radiuses_list=[],sections_list=[],option=0,s=10):
    """
explanation of the function 'multiple_sec_extrude'
path_points: are the points at which sections needs to be placed,
radiuses: radius required at each path_point. this can be '0' in case no radius required in the path
sections_list= list of sections required at each path_points. same section can be provided for various path_points as well
option: can be '0' in case the number of points in each section do not match or '1' in case number of points for each section are same
s: in case value of radiuses is provided 's' is the number of segments in that path curve
example:
sec1=circle(10)
sec2=square(5,True)
sec3=corner_radius(pts1([[-3.5,-2.5,2.49],[7,0,2.49],[0,5,2.49],[-7,0,2.49]]))
sections=[sec1,sec2,sec3,sec1,sec3]
path=array([[0,0,0],[20,2,5],[-7,25,3],[5,10,30],[-30,15,3]]).cumsum(0)
r=[0,5,7,12,0]
sol=multiple_sec_extrude(path,r,sections,0,20)
sol=cpo([ bezier(p,200)  for p in cpo(sol)])
fo(f'''
//color("blue")for(p={multiple_sec_extrude(path,r,sections,0,20)})p_line3dc(p,.1);
{swp(sol)}
    ''')
    """
    p=array(path_points)
    r=radiuses_list
    if option==0:
        sections=[sections_list[0]]+[sort_points(sections_list[0],p) for p in sections_list[1:]]
    else:
        sections=sections_list
        
    s1=[]
    for i in range(len(p)):
        if r[i]==0 and i<len(p)-1:
            p0=p[i].tolist()
            p1=(p0+(p[i+1]-p[i])*.01).tolist()
            s1.append([p0,p1])
        elif r[i]==0 and i==len(p)-1:
            p0=p[i].tolist()
            p1=(p0+(p[i]-p[i-1])*.01).tolist()
            s1.append([p0,p1])
        else:
            s1.append(fillet_3p_3d(p[i-1],p[i],p[i+1],r[i],s)[1:])
    
    s1=[remove_extra_points(p) for p in s1]

    s4=[]
    for i in range(len(s1)):
        for j in range(len(s1[i])-1):
            p0,p1=array([s1[i][j],s1[i][j+1]])
            v1=p1-p0
            va=[v1[0],v1[1]+.00001,0]
            u1=array(uv(v1))
            ua=array(uv(va))
            v2=cross(v1,va)
            a1=arccos(u1@ua)*180/pi
            a2=ang(v1[0],v1[1])
            s2=rot(f'x90z-90z{a2}',sections[i])
            s3=translate(p0,flip([axis_rot(v2,p,-a1) for p in s2]))
            s4.append(s3)
    return s4

def pntsnfaces(bead2):
    """
    function returns points and faces of a prism
    """
    n1=arange(len(bead2[0])).tolist()
    n2=array([[[[(j+1)+i*len(bead2[0]),j+i*len(bead2[0]),j+(i+1)*len(bead2[0])],[(j+1)+i*len(bead2[0]),j+(i+1)*len(bead2[0]),(j+1)+(i+1)*len(bead2[0])]] \
             if j<len(bead2[0])-1 else \
             [[0+i*len(bead2[0]),j+i*len(bead2[0]),j+(i+1)*len(bead2[0])],[0+i*len(bead2[0]),j+(i+1)*len(bead2[0]),0+(i+1)*len(bead2[0])]] \
                 for j in range(len(bead2[0]))] for i in range(len(bead2)-1)]).reshape(-1,3).tolist()
    n3=(array(flip(arange(len(bead2[0]))))+(len(bead2)-1)*len(bead2[0])).tolist()
    n=[n1]+n2+[n3]
    pnt=array(bead2).reshape(-1,3).round(4).tolist()
    return [pnt,n]



def path_offset(path,d):
    """
function to offset a 2d path by distance 'd'
example:
r=-1
path=corner_radius(pts1([[-5,0],[5,0,2],[-1,10,2],[5,0]]),20)
# path=circle(10)
# path=corner_radius(pts1([[0,0],[10,0,.1],[-10,5,.1],[-10,-5]]),20)
path1=path_offset(path,r)
fo(f'''
// orginal path
color("blue")p_lineo({path},.1);
//offset path
color("magenta")p_lineo({path1},.1);
    ''')
    """
    if d==0:
        return path
    elif len(rationalise_path(path))==2:
        return offset_segv(path,d)[:-1][0]
    else:
        a=offset_segv(path,d)[:-1]
        b=[a[0][0]]+intersections(a)[1:]+[a[-1][-1]]
        c=s_int1(seg(b))
        c=b+c if c!=[] else b
        d=cs1(path,abs(d))[:-1]
        e=[ pies1(p,c) for p in d]
        e=[p for p in e if p!=[]]
        f=remove_extra_points(concatenate(e)) if e!=[] else []
        g=exclude_points(c,f) if f!=[] else c
        g=sort_points(path,g)
        return g


def bb2d(sec):
    return [array(sec)[:,0].max()-array(sec)[:,0].min(),array(sec)[:,1].max()-array(sec)[:,1].min()]

bounding_dim2d=bb2d


    
def inner_concave_offset(sec,r):
    sec=flip(sec) if cw(sec)==1 else sec
    r=round(r,3)
    sec1=offset_segv(sec,r)
    s=intersections(sec1)
    a=s_int1(seg(s))
    if a!=[]:
        sec2=a+s
        sec2=pies1(sec,sec2)
        sec2=array(sec2)
        clean=cs1(sec,abs(r)-.01)
        clean1=[p[1:]+[p[0]] for p in clean]
        m,n,_=array(clean).shape
        o,_=sec2.shape
        v1=array([[[1,0]]*n]*m)
        v2=array(clean1)-array(clean)
        iim=array([v1,-v2]).transpose(1,2,0,3).transpose(0,1,3,2)+[0,.00001]
        im=array([pinv(iim)]*o)
        p=(array(clean)[:,:,None]-sec2).transpose(2,0,1,3)
        t=einsum('ijklm,ijkm->ijkl',im,p)
        decision1=((t[:,:,:,0]>=0)&(t[:,:,:,1]>=0)&(t[:,:,:,1]<=1))
        sec3=sec2[(decision1.sum(2)==1).any(1)]
        sec4=sort_points(sec,exclude_points(sec2,sec3))
    else:
        sec4=s
    return sec4

    
def outer_concave_offset(sec,r):
    sec=flip(sec) if cw(sec)==1 else sec
    r=round(r,3)
    sec1=offset_segv(sec,r)
    s=intersections(sec1)
    a=s_int1(seg(s))
    if a!=[]:
        sec2=a+s
        #sec2=pies1(sec,sec2)
        sec2=array(sec2)
        clean=cs1(sec,abs(r)-.01)
        clean1=[p[1:]+[p[0]] for p in clean]
        m,n,_=array(clean).shape
        o,_=sec2.shape
        v1=array([[[1,0]]*n]*m)
        v2=array(clean1)-array(clean)
        iim=array([v1,-v2]).transpose(1,2,0,3).transpose(0,1,3,2)+[0,.00001]
        im=array([pinv(iim)]*o)
        p=(array(clean)[:,:,None]-sec2).transpose(2,0,1,3)
        t=einsum('ijklm,ijkm->ijkl',im,p)
        decision1=((t[:,:,:,0]>=0)&(t[:,:,:,1]>=0)&(t[:,:,:,1]<=1))
        sec3=sec2[(decision1.sum(2)==1).any(1)]
        sec4=sort_points(sec,exclude_points(sec2,sec3))
    else:
        sec4=s
    return sec4
    
    
def outer_convex_offset(sec,d):
    segments=offset_segv(sec,d)
    return intersections(segments)


def intersections(segments):
    """
    calculates the intersections of adjacent line segments only
    """
    a=[segments[len(segments)-1]]+segments[:-1]
    b=segments
    a,b=array([a,b])
    p0,p1,p2,p3=a[:,0],a[:,1],b[:,0],b[:,1]
    v1,v2=p1-p0,p3-p2+.00001
    #     v1t1-v2t2=p2-p0
    im=inv(array([v1,-v2]).transpose(1,0,2).transpose(0,2,1))
    p=p2-p0  
    t=einsum('ijk,ik->ij',im,p)[:,0]
    p0.shape,v1.shape,t.shape
    points=(p0+einsum('ij,i->ij',v1,t)).tolist()
    return points




def c2ro(sol,s=1):#circular to rectangulat orientation
    """
change the orientation of points of a cylinder from circular to rectangular orientation
'sol': is a cylindrical type 3d shape
's': number of segments required between each straight line segments
example:
cyl=linear_extrude(circle(5,s=31),10)
cyl1=translate([15,0,0],c2ro(cyl,1))
cyl2=translate([30,0,0],cpo(c2ro(cyl,5)))
fo(f'''
color("blue")for(p={cyl})p_line3dc(p,.05);
color("blue")for(p={cyl1})p_line3dc(p,.05);
color("blue")for(p={cyl2})p_line3d(p,.05);
swp({cyl});
for(p={cyl1})swp_sec(p);
    ''')
    """
    # angle=360/len(sol[0])/2
    sol=cpo(sol)
    return [m_points1(sol[i]+flip(sol[len(sol)-1-i]),s,.001) for i in range(int(len(sol)/2))]

change_orientation_type2=c2ro



def pts2(path):
    """
returns the cumulative sum of points
path=[[0,0,1],[0,5,10],[10,3,20]]
pts2(path)=> [[0, 0, 1], [0, 5, 11], [10, 8, 31]]
example:
a=turtle3d([[0,0,0],[10,0,0],[0,0,10],[0,10,0],[-10,0,0],[0,-10,0]])
fo(f'''
color("blue") p_line3d({a},.5);
''')
    """
    return array(path).cumsum(0).tolist()

turtle3d=pts2

def axis_rot(axis,solid,angle):
    """
rotate a solid around an axis
example:
c1=cylinder(r=5,h=20)
c2=axis_rot([1,1,0],c1,90)
l1=vector2line([1,1,0],10)
fo(f'''
//axis to rotate
color("magenta") p_line3d({l1},.3);
//original cylinder
%{swp(c1)}
//rotated cylinder
color("blue",.5){swp(c2)}
''')
    """

    return (c2t3(solid)@arot(axis,angle)).tolist()

rotate_around_axis=axis_rot



def end_cap(sol,r,s=20,t=1):
    """
create a rounded edge instead of sharp edge for a solid created with linear_extrude 
or path_extrude_open function
"t" is the type of offset to be used, in most of the cases it will be default 1.
in case 1 does not work use 2 instead
example:
c1=cylinder(r=5,h=20)
e1=end_cap(c1,1)
fo(f'''
difference(){{
{swp(c1)}
for(p={e1}) swp_c(p);
}}
''')
    """
    l1=sol[0]
    l2=offset_3d(l1,-r,t)
    l3=i_p_p(sol,l1,r)
    f1=cpo(convert_3lines2fillet(l2,l3,l1,s))[:-1]
    f11=surface_offset(f1,-r)
    s1=f11+flip(f1)+[f11[0]]

    l1=sol[-1]
    l2=offset_3d(l1,-r,t)
    l3=i_p_p(sol,l1,-r)
    f1=cpo(convert_3lines2fillet(l3,l2,l1,s))[:-1]
    f11=surface_offset(f1,-r)
    s2=f11+flip(f1)+[f11[0]]

    return [s1,s2]

def end_cap_1(sol,r,s=20,t=1):
    """
create a rounded edge instead of sharp edge for a hole created with linear_extrude 
or path_extrude_open function
 "t" is the type of offset to be used, in most of the cases it will be default 1.
in case 1 does not work use 2 instead
example:
s1=linear_extrude(square(20,center=True),20)
c1=cylinder(r=5,h=20)
e1=end_cap_1(c1,2)
fo(f'''
difference(){{
{swp(s1)}
{swp(c1)}
for(p={e1}) swp(p);
}}
''')
    """
    l1=sol[0]
    l2=offset_3d(l1,r,t)
    l3=i_p_p(sol,l1,r)
    f1=cpo(convert_3lines2fillet(l2,l3,l1,s))[:-1]
    f11=surface_offset(f1,r)[0]
    s1=[f11]+f1

    l1=sol[-1]
    l2=offset_3d(l1,r,t)
    l3=i_p_p(sol,l1,-r)
    f1=cpo(convert_3lines2fillet(l3,l2,l1,s))[:-1]
    f11=surface_offset(f1,r)[-1]
    s2=f1+[f11]

    return [s1,s2]
        
def d2r(d):
    """
    converts degrees to radians
    
    """
    return l_(radians(d))
    
def r2d(r):
    """
    converts radians to degrees
    """
    return l_(rad2deg(r))
    


def convert_3lines2fillet(pnt3,pnt2,pnt1,s=10,orientation=0,r=10,style=2):
    """
Develops a fillet with 3 list of points in 3d space
s: number of segments in the fillet, increase the segments in case finer model is required
f: higher number of factor 'f' reduces the concavity, very high number like >10 will be like chamfer
example:
l1=sinewave(50,3,2,100)
l2=path_offset(l1,-5)
l3=translate([0,0,5],l1)
l1,l2=c23([l1,l2])
f1=convert_3lines2fillet(l2,l3,l1)
fo(f'''
color("blue") for(p={[l1,l2,l3]}) p_line3d(p,.3);
{swp_c(f1)}
''') 
    """
    l4=i_p_t_o(pnt1)
    r=1e5 if r==0 else r
    sol=l_(array([pnt3,pnt1,pnt2]).transpose(1,0,2))
    sol1=[]
    for i in range(len(sol)):
        p0,p1,p2=sol[i]
        la,lb,lc=a_([p0,p1,p2])
        lba=lb-la
        lbc=lb-lc
        v1=cross(lba,lbc)
        v1=v1/norm(v1)
        clock=1 if norm(v1-l4[i])>1 else -1
        if style==1:
            sol1.append(arc_2p_3d(l4[i],p0,p2,r,clock,s)+[p1])
        elif style==2:
            sol1.append(bezier([p0,p1,p2],s)+[p1])

    return sol1 if orientation==0 else cpo(sol1)[:-1]
    
def min_d_points(sec,min_d=.1):
    """ 
rationalises the number points in a section based on the minimum distance between 2 points
i.e. all the points which are less than the defined minimum distance "min_d" will be omitted from the section "sec" 
example:
l1=circle(10,s=100)
l2=translate_2d([25,0],min_d_points(l1,min_d=5))
l3=translate_2d([50,0],homogenise(l1,pitch=5,closed_loop=1))
fo(f'''
color("blue") for(p={[l1,l2,l3]}) points(p,.3);
''')
    """
    b=sec[0]
    c=[sec[0]]
    for i in range(1,len(sec)):
        if l_len([b,sec[i]])>=min_d:
            c.append(sec[i])
            b=sec[i]
            
    return c
    
def wrap_around(sec,path):
    """
function to wrap a line or section around a 3d path. This only works with path in y-direction and section/ polyline strictly in +y quadrants
example:
c1=translate([0,20.1,0],circle(20))
path=rot('y90',circle(40.2/(2*pi)+.2))
c2=wrap_around(c1,path)
fo(f'''
color("blue") p_line3d({c1},.2);
color("cyan") p_line3d({path},.2);
color("magenta") p_line3d({c2},.2);
''')    
    """
    line1=array(sec) if array(sec).shape[-1]==3 else array(c2t3(sec))
    c_dist_path=array([0]+[l_len(p) for p in seg(path)[:-1]]).cumsum()
    v1_path=array([array(p[0])-array(p[1]) for p in seg(path)[:-1]])
    v1_path=v1_path/norm(v1_path,axis=1).reshape(-1,1)
    c_dist_line=line1[:,1]
    sec1=[]
    for i in range(len(line1)):
        if c_dist_line[i]==0:
            a=c_dist_path[1]-c_dist_line[i]
            b=array(path[1])+v1_path[0]*a
            c=cross(v1_path[0],[1,0,0])
            c=c/norm(c)
            c=b+c*line1[i][2]+array([line1[i][0],0,0])
        else:
            j=arange(len(path))[c_dist_line[i]>c_dist_path][-1]+1
            a=c_dist_path[j]-c_dist_line[i]
            b=array(path[j])+v1_path[j-1]*a
            c=cross(v1_path[j-1],[1,0,0])
            c=c/norm(c)
            c=b+c*line1[i][2]+array([line1[i][0],0,0])
        sec1.append(c.tolist())  
    return sec1

wrap_section_around_path=wrap_around
    
def align_sec(sec1,sec2,ang=10):
    """
function to align 2 3d sections to obtain the non twisted optimised solid
ang: is the resolution for the angle of rotation, 1 degree will have much higher resolution and hence will take longer to compute
example:
l1=c23(square(4))
c1=translate([2,2,5],circle(5,s=5))
c2=align_sec(l1,c1,ang=1)[1]
fo(f'''
//original lines
color("blue") for(p={[l1,c1]}) p_line3d(p,.3);
// aligned section
color("magenta") p_line3d({c2},.3);
''')
    """
    nv1=nv(sec2)
    cp1=array(sec2).mean(0)
    sec2=translate(-cp1,sec2)
    i=arange(0,360,ang)
    area1=[norm((array(axis_rot(nv1,sec2,j))+cp1)-array(sec1),axis=1).sum() for j in i]
    sec3=(array(axis_rot(nv1,sec2,array(area1).argmin()*ang))+cp1).tolist()
    sol2=[sec1]+[sec3]
    return sol2
    

def s2v(v1,sec):
    sec=rot('x90z-90',sec)
    vx=[1,0,0]
    v2=uv(c23(v1[:2]))
    theta1=ang_v(v1[:2])
    theta2=r2d(arccos(a_(v2)@a_(v1)))
    ax1=cross(v2,v1)
    b=axis_rot(ax1,rot(f'z{theta1}',sec),theta2)
    return b
    
def sec2vector(v1=[1,0,0],sec=[]):
    """
function to align a section 'sec' with a vector 'v1'
example:
l1=c23(square(4))
v1=[1,-1,0.5]
l2=sec2vector(v1,l1)
l3=vector2line(v1,5)
fo(f'''
//original section
color("blue") for(p={[l1]}) p_line3dc(p,.3);
// vector to align section
color("magenta") p_line3d({l3},.3);
//aligned sec
color("cyan") p_line3dc({l2},.3);
''')
    """
    vz=[0,0,-1]
    vz,v1=array([vz,v1])
    # v1=v1+[.0001,.0001,0]
    nvzv1=cross(vz,v1)
    u1=v1/norm(v1)
    theta=r2d(arccos(u1@vz))
    sec=c2t3(flip(sec))

    sec1=sec@arot([1,0,0],theta)

    theta1=ang(v1[0],v1[1])
    sec1=sec1@zrot(-90)@zrot(theta1)
    return sec1.tolist()

align_a_section_to_vector=sec2vector
    
def sec2vector1(v1,sec):
    """
same as sec2vector but simpler method
example:
l1=c23(square(4))
v1=[1,-1,0.5]
l2=sec2vector1(v1,l1)
l3=vector2line(v1,5)
fo(f'''
//original section
color("blue") for(p={[l1]}) p_line3dc(p,.3);
// vector to align section
color("magenta") p_line3d({l3},.3);
//aligned sec
color("cyan") p_line3dc({l2},.3);
''')
    """
    theta_y=ang((v1[0]**2+v1[1]**2)**.5,v1[2])
    theta_z=ang(v1[0],v1[1])
    return rot(f'x90z-90y{-theta_y}z{theta_z}',sec)
    return (sec@xrot(90)@zrot(-90)@yrot(-theta_y)@zrot(theta_z)).tolist()


def cut_plane(nv=[0,0,1],size=[5,5],thickness=10,trns1=0,trns2=0,trns3=0,theta=[0,0,0]): #oriented solid
    """
function for defining a solid (cutting plane) oriented as per the defined normal vector
nv: normal vector for defining plane orientation of the section
thickness: thickness or height of the cutting plane
trns1: translate the solid in the direction of normal vector 'nv'
trns2: translate the solid in the direction 'right' to the normal vector 'nv'
trns3: translate the solid in the direction 'up' to the normal vector 'nv'
'-ve' values given to the trns1,trns2,trns3 will translate the solid in the reverse direction
    """
    sec=square(size,center=True)
    plane1=sec2vector(nv,sec)
    v1=array(nv)
    u1=v1/norm(v1)
    ua=array([0,0,1]) if u1[2]==0 else array([0,-1,0]) if (u1==[0,0,1]).all() else array([0,1,0]) if (u1==[0,0,-1]).all() else array([0,0,1])
    v2=cross(u1,ua) if u1[2]>0 else cross(u1,ua)
    u2=v2/norm(v2)
    v3=cross(u2,u1) if u1[2]==0 else cross(u2,u1)
    u3=v3/norm(v3)
    plane2=translate(u1*thickness,plane1)
    sol=[plane1]+[plane2]
    sol=axis_rot(u3,axis_rot(-u2,axis_rot(nv,sol,theta[0]),theta[1]),theta[2])
    sol=translate(u1*trns1,sol)
    sol=translate(u2*trns2,sol)
    sol=translate(u3*trns3,sol)
    return sol
    
    
def slice_sol(sol,n=10):
    """
function to slice a solid with 'n' intermediate steps.
this creats n steps for each turn in sol
example:
l1=c23(square(4))
l2=translate([0,0,5],circle(5,[2,2],s=5))
l3=translate([0,0,10],l1)
sol=[l1,l2,l3]
sol1=translate([10,0,0],sol)
sol1=slice_sol(sol1,10)
fo(f'''
// original solid
color("blue") for(p={sol}) p_line3dc(p,.1);
{swp(sol)}
// sliced solid
color("magenta") for(p={sol1}) p_line3dc(p,.1);
{swp(sol1)}
''')
    """
    a=cpo(sol)
    sol1=[[ls(p,n)+[p[1]] for p in seg(a[i])[:-1]] for i in range(len(a))]
    sol2=array(sol1).transpose(1,2,0,3)
    b,c,d,e=sol2.shape
    sol2=sol2.reshape(b*c,d,e).tolist()
    return sol2

def slice_sol_1(sol_1,n=10):
    """
function to slice a solid with 'n' intermediate steps
example:
l1=c23(square(4))
l2=translate([0,0,5],circle(5,[2,2],s=5))
l3=translate([0,0,10],l1)
sol=[l1,l2,l3]
sol1=translate([10,0,0],sol)
sol1=slice_sol_1(sol1,10)
fo(f'''
// original solid
color("blue") for(p={sol}) p_line3dc(p,.1);
{swp(sol)}
// sliced solid
color("magenta") for(p={sol1}) p_line3dc(p,.1);
{swp(sol1)}
''')
    """
    return cpo([equidistant_path(p,n) for p in cpo(sol_1)])

    
def cp_arc(arc1):
    """
function returns the center point of a given circle or arc
example:
a1=circle(10,[5,6])
cp1=cp_arc(a1)
fo(f'''
// arc
color("blue") for(p={[a1]}) p_line3dc(p,.1);

// center point calculated
color("magenta") points({[cp1]},.3);
''')
cp1 # calculated center point
    """
    n=int(len(arc1)/360*120)
    p0=arc1[0]
    p1=arc1[n]
    p2=arc1[n*2]
    return cp_3p(p0,p1,p2)

center_arc=cp_arc
    
def r_arc(arc1):
    """
function returns the radius of a given circle or arc
example:
a1=circle(10,[5,6])
r1=r_arc(a1)
txt1=dim_radial(a1[:5])
fo(f'''
// arc
color("blue") for(p={[a1]}) p_line3dc(p,.1);
// calculated radius
{txt1}
''')
r1 # calculated radius
    """
    n=int(len(arc1)/360*120)
    p0=arc1[0]
    p1=arc1[n]
    p2=arc1[n*2]
    return l_(r_3p([p0,p1,p2]))

radius_arc2d=r_arc
    

def o_solid(nv=[0,0,1],sec=[],thickness=10,trns1=0,trns2=0,trns3=0, theta=[0,0,0]): #oriented solid
    """
function for defining a solid with any defined section. solid gets oriented as per the defined normal vector
nv: normal vector for defining plane orientation of the section
sec: cross section of the solid
thickness: thickness or height of the solid
trns1: translate the solid in the direction of normal vector 'nv'
trns2: translate the solid in the direction 'right' to the normal vector 'nv'
trns3: translate the solid in the direction 'up' to the normal vector 'nv'
'-ve' values given to the trns1,trns2,trns3 will translate the solid in the reverse direction 
theta: rotate the section around axis  fox example if nv is [1,0,0] or x-axis, the sequence of rotation will be x, y ,z axis
example:
sol=o_solid(nv=[1,0,0],sec=circle(2.5),thickness=20,trns1=-5,trns2=10,trns3=20)
l1=turtle3d([[0,0,20],[-5,0,0],[0,-10,0],[0.001,0,-20]])
txt1=dim_linear(l1[:2])
txt2=dim_linear(l1[1:3])
txt3=dim_linear(l1[2:4])
fo(f'''
// oriented_solid
{swp(sol)}
{txt1}{txt2}{txt3}
''')
    """
    plane1=sec2vector(nv,sec)
    v1=array(nv)
    u1=v1/norm(v1)
    ua=array([0,0,1]) if u1[2]==0 else array([0,-1,0]) if (u1==[0,0,1]).all() else array([0,1,0]) if (u1==[0,0,-1]).all() else array([0,0,1])
    v2=cross(u1,ua) if u1[2]>0 else cross(u1,ua)
    u2=v2/norm(v2)
    v3=cross(u2,u1) if u1[2]==0 else cross(u2,u1)
    u3=v3/norm(v3)

    plane2=translate(u1*thickness,plane1)
    sol=[plane1]+[plane2]
    sol=axis_rot(u3,axis_rot(-u2,axis_rot(nv,sol,theta[0]),theta[1]),theta[2])
    sol=translate(u1*trns1,sol)
    sol=translate(u2*trns2,sol)
    sol=translate(u3*trns3,sol)
    return sol

oriented_solid=o_solid

def ppesec(p0,sec): #point's projection on an enclosed 3d section
    """
function to find projected points of a given point list 'p0' on a 3d sec which is on 1 plane
"""
    
    v1=array(nv(sec))
    u1=v1/norm(v1)
    sec1=array(sec)
    v2=seg(array([sec1[0]-p for p in sec1[1:]]).tolist())[:-1]
    iim=array([[u1.tolist()]+p for p in v2]).transpose(0,2,1)
    im=array([inv(iim)]*len(p0))
    p=array([sec1[0]-array(p0)]*len(v2)).transpose(1,0,2)
    t=einsum('ijkl,ijl->ijk',im,p)
    decision=(t[:,:,1]>=0)&(t[:,:,1]<=1)&(t[:,:,2]>=0)&(t[:,:,2]<=1)&((t[:,:,1]+t[:,:,2])<=1)
    t1=t[decision][:,0]
    p1=array([p0]*len(v2)).transpose(1,0,2)
    p1=p1[decision]
    ip1=(p1+u1*t1[:,None]).tolist()
    ip2=p1.tolist()
    return [ip1,ip2]

def ppplane(p0,v1,loc):#point's projection on a plane
    """
function to find projected points of a given list of points 'p0' on a plane defined by normal'v1' and location 'loc'
example:
p0=[20,0,0]
v1=[2,3,4]
loc=[0,10,0]
ppplane([p0],v1,loc) => [19.310359216374945, -1.034461175437585, -1.3792815672501133]
    
    """
    p0=array(p0)
    sec=pts([[-50/2,-50/2],[50,0],[0,50],[-50,0]])
    plane1=translate(loc,o_solid(v1,sec,.001))
    v1=array(v1)
    u1=v1/norm(v1)
    sec1=array(plane1[0])
    v2=seg(array([sec1[0]-p for p in sec1[1:]]).tolist())[:-1]
    iim=array([[u1.tolist()]+p for p in v2]).transpose(0,2,1)
    im=array([inv(iim)]*len(p0))
    p=sec1[0]-array(p0)
    t=einsum('ijkl,il->ijk',im,p)[:,:,0][:,0]
    ip1=(p0+u1*t[:,None]).tolist()
    return ip1

def honeycomb(r,n1,n2):
    """
function to draw a honeycomb structure with radius 'r' 
n1: number of hexagons in 1 layer
n2: number of layers
example:
l1=honeycomb(5,3,2)
fo(f'''
color("blue") for(p={l1}) p_line3dc(p,.2);
''')
    """
    cir1=circle(r,s=7)
    cir2=c3t2(rot('z30',cir1))
    sec=[translate([i,0,0],cir1) for i in arange(0,3*n1*r,3*r)]
    sec1=[translate([i,r*sin(d2r(60)),0],cir1) for i in arange(r*1.5,3*n1*r,3*r)]
    sec2=array([sec,sec1]).transpose(1,0,2,3)
    a,b,c,d=sec2.shape
    sec3=sec2.reshape(a*b,c,d)
    sec3=array([translate([0,i,0],sec3) for i in arange(0,r*sin(d2r(60))*(n2*2),2*r*sin(d2r(60)))])
    a,b,c,d=sec3.shape
    sec4=c3t2(sec3.reshape(a*b,c,d))
    return sec4

    
    
def path_extrude2msec(sec_list,path):
    """
function to extrude multiple sections 'sec_list' along an open path 'path'
number of sections in the 'sec_list' >= len(path)
example:
sec=square(4)
p1=cr2dt([[0,0],[2,0],[3,0],[-6,0]])
sec_list=c32(prism(sec,p1))
path=cr3dt([[0,0,0],[20,0,0],[0,0,20],[30,0,0]])
sol=path_extrude2msec(sec_list,path)
fo(f'''
{swp(sol)}
''')
    """
    p1=path[:-1]
    p2=path[1:]
    p1,p2=array([p1,p2])
    v1=p2-p1
    u1=v1/norm(v1,axis=1).reshape(-1,1)
    v2=concatenate([[u1[0]],(u1[1:]+u1[:-1])/2,[u1[-1]]])
    sec2=[]
    for i in range(len(path)):
        sec=sec_list[i]
        sec1=translate(path[i],sec2vector(v2[i],sec))
        sec2.append(sec1)
    return sec2
    
def sol2vector(v1=[],sol=[],loc=[0,0,0]):
    sol1=translate(loc,[sec2vector(-array(v1),p) for p in sol])
    return sol1
    

    
def ip_sol2line(sol,line):# when line has more than 2 points
    """
function to calculate intersection point between a 3d solid and a line. 
 "sol" is the 3d object which is intersected with a "line".
 try below code for better understanding:
example:
sec=circle(10)
path=cr2d(pts1([[-10+.1,0],[12,0],[-2,0,2],[0,10,3],[-10,0]]),5)
sol=prism(sec,path)
line=[[0,0,-1],[20,20,10]]
ip1=ip_sol2line(sol,line)
fo(f'''
color("blue") p_line3d({line},.2);
color("magenta") points({ip1},.5);
%{swp(sol)}
''')
    """


    pa=sol
    p1=array([[ [[pa[i][j],pa[i][j+1],pa[i+1][j]],[pa[i+1][j+1],pa[i+1][j],pa[i][j+1]]] if j<len(pa[i])-1 
     else [[pa[i][j],pa[i][0],pa[i+1][j]],[pa[i+1][0],pa[i+1][j],pa[i][0]]] 
     for j in range(len(pa[i]))] 
              for i in range(len(pa)-1)]).reshape(-1,3,3)
    pm=p1[:,0]
    pn=p1[:,1]
    po=p1[:,2]
    px=array(line[:-1])
    py=array(line[1:])
    v1,v2,v3=py-px,pn-pm,po-pm
    a,_=v1.shape
    b,_=v2.shape
    v1=array([v1]*b)
    v2=-array([v2]*a).transpose(1,0,2)
    v3=-array([v3]*a).transpose(1,0,2)
    iim=array([v1,v2,v3]).transpose(1,2,0,3).transpose(0,1,3,2)+.00001
    im=inv(iim)
    p=array([pm]*a).transpose(1,0,2)-array([px]*b)
    t=einsum('ijkl,ijl->ijk',im,p)
    condition=(t[:,:,0]>=0)&(t[:,:,0]<=1)&(t[:,:,1]>=0)&(t[:,:,1]<=1)&(t[:,:,2]>=0)&(t[:,:,2]<=1)&((t[:,:,1]+t[:,:,2])<=1)
    t1=t[:,:,0][condition]
    i_p1=array([px]*b)[condition]+einsum('ij,i->ij',v1[condition],t1)
    i_p2=i_p1[argsort([norm(p-px[0]) for p in i_p1])]

    return i_p2.tolist()
    


def align_sol(sol,ang=10):
    """
function to straighten the twists in the path_extruded sections for better alignments
example:
l1=c23(square(4))
c1=translate([2,2,5],circle(5,s=5))
sol=[l1,c1]
sol1=align_sol(sol,ang=1)
fo(f'''
//original solid
{swp(sol)}
// aligned solid
translate([10,0,0]) {swp(sol1)}
''')
    """
    sol1=[sol[0]]
    for i in range(1,len(sol)):
        a=align_sec(sol1[i-1],sol[i],ang=ang)
        sol1.append(a[1])
    return sol1

def extrude_sol2path(sec,path1,path2):
    """
extrude a solid to a different path
"sec" and "path1" defines the original solid
"path2" defines the path where the shape of the original solid to be extruded
use function sol2path instead
    """
    min_l=array([l_len(p) for p in seg(path1)[:-1]]).min()
    path1=m_points_o(path1,min_l)
    l1=array([l_len(p) for p in seg(path2)[:-1]]).sum()
    l2=l1/len(path1)
    path3=m_points_o(path2,l2)[:len(path1)]

    sec_list=[offset(sec,x) for (x,y) in path1]
    sol=path_extrude2msec(sec_list,path3)
    return sol
    


def ip_normal_sol2line(sol,line):
    """
    function to find the normal from intersection points between a 3d solid and a line. 
     "sol" is the 3d object which is intersected with a "line".
    """


    pa=sol
    p1=array([[ [[pa[i][j],pa[i][j+1],pa[i+1][j]],[pa[i+1][j+1],pa[i+1][j],pa[i][j+1]]] if j<len(pa[i])-1 
     else [[pa[i][j],pa[i][0],pa[i+1][j]],[pa[i+1][0],pa[i+1][j],pa[i][0]]] 
     for j in range(len(pa[i]))] 
              for i in range(len(pa)-1)]).reshape(-1,3,3)
    pm=p1[:,0]
    pn=p1[:,1]
    po=p1[:,2]
    px=array(line[:-1])
    py=array(line[1:])
    v1,v2,v3=py-px,pn-pm,po-pm
    a,_=v1.shape
    b,_=v2.shape
    v1=array([v1]*b)
    v2=-array([v2]*a).transpose(1,0,2)
    v3=-array([v3]*a).transpose(1,0,2)
    iim=array([v1,v2,v3]).transpose(1,2,0,3).transpose(0,1,3,2)+.00001
    im=inv(iim)
    p=array([pm]*a).transpose(1,0,2)-array([px]*b)
    t=einsum('ijkl,ijl->ijk',im,p)
    condition=(t[:,:,0]>=0)&(t[:,:,0]<=1)&(t[:,:,1]>=0)&(t[:,:,1]<=1)&(t[:,:,2]>=0)&(t[:,:,2]<=1)&((t[:,:,1]+t[:,:,2])<=1)
    t1=t[:,:,0][condition]
    i_p1=array([px]*b)[condition]+einsum('ij,i->ij',v1[condition],t1)
    i_p2=i_p1[argsort([norm(p-px[0]) for p in i_p1])]
    s_planes=array([p1]*a).transpose(1,0,2,3)[condition][argsort([norm(p-px[0]) for p in i_p1])]
    nv1=[nv(p) for p in s_planes]

    i_p2,s_planes,nv1=i_p2.tolist(),s_planes.tolist(),array(nv1).tolist()
    un1=array(nv1)/norm(array(nv1),axis=1).reshape(-1,1)
    return un1.tolist()
    
def pntsnfaces_c(bead2):
    """
    returns vertices and faces for a closed loop solid
    
    """
    n1=arange(len(bead2[0])).tolist()
    n2=array([[[[(j+1)+i*len(bead2[0]),j+i*len(bead2[0]),j+(i+1)*len(bead2[0])],[(j+1)+i*len(bead2[0]),j+(i+1)*len(bead2[0]),(j+1)+(i+1)*len(bead2[0])]] \
             if j<len(bead2[0])-1 else \
             [[0+i*len(bead2[0]),j+i*len(bead2[0]),j+(i+1)*len(bead2[0])],[0+i*len(bead2[0]),j+(i+1)*len(bead2[0]),0+(i+1)*len(bead2[0])]] \
                 for j in range(len(bead2[0]))] for i in range(len(bead2)-1)]).reshape(-1,3).tolist()
    n3=(array(flip(arange(len(bead2[0]))))+(len(bead2)-1)*len(bead2[0])).tolist()
    n=[n1]+n2+[n3]
    pnt=array(bead2).reshape(-1,3).round(4).tolist()
    return [pnt,n2]

def vnf1(surf):
    """
    function calculates vertices and faces for a given surface
    refer file 'example of various functions.ipynb' for application examples
    """
    n1,n2,_=array(surf).shape
    v=array(surf).reshape(-1,3)    
    f1=array([[[[i*n2+j,i*n2+j+1,(i+1)*n2+j],[(i+1)*n2+j,i*n2+j+1,(i+1)*n2+j+1]] for j in range(n2-1)] for i in range(n1-1)]).reshape(-1,3)
    return [v.tolist(),f1.tolist()]
    

    
def equidistant_path(path,s=10,pitch=[]):
    """
divides a path in to equally spaced points
example:
l1=c23(square(20))
l2=equidistant_path(l1,20)
fo(f'''
//original line
color("blue",.2) p_line3d({l1},.3);
color("cyan") points({l1},1);
// equidistant_path line

translate([30,0,0]) {{
color("magenta",.2) p_line3d({l2},.3);
color("cyan") points({l2},1);
}}
''')
    """
    s= l_lenv_o(path)/pitch if a_(pitch).size>0 else s
    v=[p[1]-p[0] for p in array(seg(path)[:-1])]
    l=[l_len(p) for p in seg(path)[:-1]]
    c=array(l).cumsum().tolist()
    l1=c[-1]/s
    d=[l1*i for i in arange(1,s+1)]
    p_rev=[]
    for i in range(len(c)):
        for j in range(len(d)):
            if c[i]>d[j]:
                t=d[j]/l[i] if i==0 else (d[j]-c[i-1])/l[i]
                px=array(path[i])+array(v[i])*t
                p_rev.append(px.tolist())
                d[j]=c[-1]+1
    p_rev=[path[0]]+p_rev+[path[-1]]
    return p_rev[:int(s)+1] if s%1==0 else p_rev[:int(s)+1]+[path[-1]]

def equidistant_pathc(path,s=10,pitch=[]):
    """
divides a closed path in to equally spaced points
example:
l1=c23(square(20))
l2=equidistant_pathc(l1,20)
fo(f'''
//original line
color("blue",.2) p_line3dc({l1},.3);
color("cyan") points({l1},1);
// equidistant_path line

translate([30,0,0]) {{
color("magenta",.2) p_line3dc({l2},.3);
color("cyan") points({l2},1);
}}
''')
    """
    s= l_lenv(path)/pitch if a_(pitch).size>0 else s
    v=[p[1]-p[0] for p in array(seg(path))]
    l=[l_len(p) for p in seg(path)]
    c=array(l).cumsum().tolist()
    l1=c[-1]/s
    d=[l1*i for i in arange(1,s+1)]
    p_rev=[]
    for i in range(len(c)):
        for j in range(len(d)):
            if c[i]>d[j]:
                t=d[j]/l[i] if i==0 else (d[j]-c[i-1])/l[i]
                px=array(path[i])+array(v[i])*t
                p_rev.append(px.tolist())
                d[j]=c[-1]+1
    p_rev=[path[0]]+p_rev
    return p_rev[:int(s)+1] if a_(pitch).size>0 else p_rev[:s]


    
def ang_2lineccw(p0,p1,p2):
    """
ccw angle of the line p0p2 from base line p0p1
example:
p0=[0,0]
p1=[5,2]
p2=[7,8]
txt1=dim_angular([p0,p1],[p1,p2],text_size=.5)
txt2=dim_angular([p1,p0],[p0,p2],text_size=.5)
fo(f'''

// case 1: counter-clockwise angle between line p1p0 and p1p2 is 229.76 degrees
color("blue")points({[p0,p1,p2]},.2);
color("cyan")p_line3d({[p0,p1,p2]},.1);
color("magenta")translate([0,-.5]){{
translate({p0})text("p0",.5);
translate({p1})text("p1",.5);
translate({p2})text("p2",.5);
}}
{txt1}
// case 2: counter-clockwise angle between line p0p1 and p0p2 is 27.01 degrees

translate([10,0,0]){{
color("blue")points({[p0,p1,p2]},.2);
color("cyan")p_line3d({[p1,p0,p2]},.1);

color("magenta")translate([0,-.5]){{
translate({p0})text("p0",.5);
translate({p1})text("p1",.5);
translate({p2})text("p2",.5);
}}
{txt2}
}}

''')
ang_2lineccw(p1,p0,p2),ang_2lineccw(p0,p1,p2)
    """
    p0,p1,p2=array([p0,p1,p2])
    v1,v2=p1-p0,p2-p0
    a1=ang(v1[0],v1[1])
    a2=ang(v2[0],v2[1])
    return 360 if a1-a2==0 else 360-(a1-a2) if a2<a1 else a2-a1 


def ang_2linecw(p0,p1,p2):
    """
cw angle of the line p0p2 from the base line p0p1
example:
p0=[0,0]
p1=[5,2]
p2=[7,8]
fo(f'''
// case 1: clockwise angle between line p1p0 and p1p2 is 130.23 degrees
color("blue")points({[p0,p1,p2]},.2);
color("cyan")p_line3d({[p0,p1,p2]},.1);
color("magenta")translate([0,-.5]){{
translate({p0})text("p0",.5);
translate({p1})text("p1",.5);
translate({p2})text("p2",.5);
}}

// case 2: clockwise angle between line p0p1 and p0p2 is 332.98 degrees
translate([10,0,0]){{
color("blue")points({[p0,p1,p2]},.2);
color("cyan")p_line3d({[p1,p0,p2]},.1);

color("magenta")translate([0,-.5]){{
translate({p0})text("p0",.5);
translate({p1})text("p1",.5);
translate({p2})text("p2",.5);
}}
}}
''')
l_(ang_2linecw(p1,p0,p2)),l_(ang_2linecw(p0,p1,p2))
    """
    p0,p1,p2=array([p0,p1,p2])
    v1,v2=p1-p0,p2-p0
    a=ang(v1[0],v1[1])
    b=ang(v2[0],v2[1])
    return l_(abs(a-b) if b<=a else (360-b)+a)

def l_lenv(l):
    """
calculates sum of lengths of all the segments in a line 'l' considering the section is closed
example:
l_lenv(square(5)) => 20.0
    """
    return l_(array([l_len(p) for p in seg(l)]).sum())

length_closed_section=l_lenv

def l_lenv_o(l):
    """
    calculates sum of lengths of all the segments in a line 'l' considering the section is open
example:
l_lenv_o(square(5)) => 15.0
    """
    return l_(array([l_len(p) for p in seg(l)[:-1]]).sum())

length_open_section=l_lenv_o

def a_3seg(s):
    """
area of the triangle enclosed with in 3 vertices 's'
example:
a_3seg([[0,0,0],[10,0,0],[0,0,10]]) => 50.0
    """
    return l_(norm(cross(array(s[1])-array(s[0]),array(s[2])-array(s[0])))/2)



def offset_sol(sol,d,type=1):
    """
function to calculate offset of a 3d object by distance 'd'
example:
c1=cylinder(r=10,h=40)
c2=offset_sol(c1,-3)
c3=offset_sol(c1,3)
fo(f'''
//original cylinder
color("blue",.2) swp_c({c1});
// offset inwards by 3mm
color("magenta",.2) swp_c({c2});
// offset outwards 3mm
color("cyan",.2) swp_c({c3});
''')
    """
    
    sol=[offset_3d(p,d,type) for p in sol]
    return sol
    


def ip_sol2sol(sol1,sol2,n=0):
    """
function finds the intersection point between 2 solids
example:
s1=sphere(20)
c1=cylinder(r=5,h=50)
l1=ip_sol2sol(s1,c1)
fo(f'''
%{swp(s1)}
%{swp(c1)}
color("blue") points({l1},.3);
''')
    
    """
    line=array([ seg(p)[:-1] for p in cpo(sol2)])
    v,f1=vnf2(sol1)
    tri=array(v)[array(f1)]
    line=array([ seg(p)[:-1] for p in cpo(sol2)])
    tri.shape,line.shape
    la,lb=line[:,:,0],line[:,:,1]
    p0,p1,p2=tri[:,0],tri[:,1],tri[:,2]
    lab=lb-la
    p01,p02=p1-p0,p2-p0
    t=einsum('kl,ijkl->ijk',cross(p01,p02),la[:,:,None]-p0)/(einsum('ijl,kl->ijk',(-lab),cross(p01,p02))+.00000)
    u=einsum('ijkl,ijkl->ijk',cross(p02[None,None,:,:],(-lab)[:,:,None,:]),(la[:,:,None,:]-p0[None,None,:,:]))/(einsum('ijl,kl->ijk',(-lab),cross(p01,p02))+.00000)
    v=einsum('ijkl,ijkl->ijk',cross((-lab)[:,:,None,:],p01[None,None,:,:]),(la[:,:,None,:]-p0[None,None,:,:]))/(einsum('ijl,kl->ijk',(-lab),cross(p01,p02))+.00000)
    condition=(t>=0)&(t<=1)&(u>=0)&(u<=1)&(v>=0)&(v<=1)&(u+v<1)

    a=(la[:,None,:,None,:]+lab[:,None,:,None,:]*t[:,None,:,:,None])
    b=condition[:,None,:,:]
    c=[]
    for i in range(len(a)):
        c.append(a[i][b[i]].tolist())

    return [p[n] for p in c if p!=[]]

def ip_surf2sol(sol1,sol2,n=0):
    """
intersection points between surface to surface or surface to solid
example:
l0=translate([-25,0,0],rot('x90',cosinewave(50,3,2,100)))
s1=surface_line_vector(l0,[0,50,0],1)
l1=translate([-10,0,-10],sinewave(20,2,2,50))
s2=surface_line_vector(l1,[5,5,50])
l1=ip_surf2sol(s1,s2)
fo(f'''
%{swp_surf(s1)}
%{swp_surf(s2)}
color("blue") points({l1},.3);
''')
    """
    line=array([ seg(p)[:-1] for p in cpo(sol2)])
    v,f1=vnf1(sol1)
    tri=array(v)[array(f1)]
    line=array([ seg(p)[:-1] for p in cpo(sol2)])
    tri.shape,line.shape
    la,lb=line[:,:,0],line[:,:,1]
    p0,p1,p2=tri[:,0],tri[:,1],tri[:,2]
    lab=lb-la
    p01,p02=p1-p0,p2-p0
    t=einsum('kl,ijkl->ijk',cross(p01,p02),la[:,:,None]-p0)/(einsum('ijl,kl->ijk',(-lab),cross(p01,p02))+.00000)
    u=einsum('ijkl,ijkl->ijk',cross(p02[None,None,:,:],(-lab)[:,:,None,:]),(la[:,:,None,:]-p0[None,None,:,:]))/(einsum('ijl,kl->ijk',(-lab),cross(p01,p02))+.00000)
    v=einsum('ijkl,ijkl->ijk',cross((-lab)[:,:,None,:],p01[None,None,:,:]),(la[:,:,None,:]-p0[None,None,:,:]))/(einsum('ijl,kl->ijk',(-lab),cross(p01,p02))+.00000)
    condition=(t>=0)&(t<=1)&(u>=0)&(u<=1)&(v>=0)&(v<=1)&(u+v<1)

    a=(la[:,None,:,None,:]+lab[:,None,:,None,:]*t[:,None,:,:,None])
    b=condition[:,None,:,:]
    c=[]
    for i in range(len(a)):
        c.append(a[i][b[i]].tolist())

    return [p[n] for p in c if p!=[]]


    
def vnf2(bead2):
    """
    function returns vertices and faces of 3d shapes with first and the last section triangulated. only works with convex sections
    
    """
    n1=arange(len(bead2[0])).tolist()
    cp1=array(bead2[0]).mean(0).tolist()
    n1=[[0,p[0]+1,p[1]+1] for p in seg(n1)]
    n2=array([[[[(j+1)+i*len(bead2[0]),j+i*len(bead2[0]),j+(i+1)*len(bead2[0])],[(j+1)+i*len(bead2[0]),j+(i+1)*len(bead2[0]),(j+1)+(i+1)*len(bead2[0])]] \
             if j<len(bead2[0])-1 else \
             [[0+i*len(bead2[0]),j+i*len(bead2[0]),j+(i+1)*len(bead2[0])],[0+i*len(bead2[0]),j+(i+1)*len(bead2[0]),0+(i+1)*len(bead2[0])]] \
                 for j in range(len(bead2[0]))] for i in range(len(bead2)-1)]).reshape(-1,3)
    n2=(n2+1).tolist()
    n3=(array(flip(arange(len(bead2[0]))))+(len(bead2)-1)*len(bead2[0]))
    n3=n3+1
    cp2=array(bead2[-1]).mean(0).tolist()
    n3=[[len(bead2)*len(bead2[0])+1,p[0],p[1]] for p in seg(n3)]

    n=n2
    pnt=[cp1]+array(bead2).reshape(-1,3).round(4).tolist()+[cp2]
    return [pnt,n]
    
    
def convex_hull(pnts):
    """
calculates convex hull for a list of points 'sec'
example:
a=translate_2d([10,-5],square(10))
b=circle(10)
c=convex_hull(a+b)
fo(f'''
color("blue") for(p={[a,b]}) p_line3dc(p,.2);
color("magenta") p_line3dc({c},.21);
''')
    """
    pnts=remove_extra_points(pnts)
    p_x=[]
    a=array(pnts)
    b=a[a[:,1]==a[:,1].min()]
    if len(b)>1:
        c=b[b[:,0].argmin()]

    s1=c if len(b)>1 else b[0]
    p_x.append(s1.tolist())
    a=array(exclude_points(a,s1))

    b1=a[array([l_len([s1,p]) for p in a]).argsort()]
    c=b1[array([ang_2linecw(s1,[s1[0],-1e5],p) for p in b1]).argmax()]
    p_x.append(c.tolist())
    a=array(exclude_points(a,c)+[s1])

    while (c.tolist()!=s1.tolist() ):
        b=a[array([l_len([p_x[-1],p]) for p in a ]).argsort()]
        c=b[array([ang_2linecw(p_x[-1],p_x[-2],p) for p in b]).argmax()]
        p_x.append(c.tolist())
        a=array(exclude_points(a,c))
    return p_x[:-1]

def lexicographic_sort_xy(p):
    """
function sorts the points list 'p' first with x and then with y smallest to largest 
example:
a=random.random(50)*(20-5)+5
b=random.random(50)*(30-10)+10
p=array([a.round(1),b.round(1)]).transpose(1,0)
fo(f'''
p2={lexicographic_sort_xy(p)};
color("blue")for(i=[0:len(p2)-1])translate(p2[i])text(str(i),.3);
p_lineo(p2,.05);
color("magenta")points(p2,.2);
    ''')
    """
    p=array(p)
    p1=p[p[:,0].argsort()]
    pux=unique(p1[:,0])
    p2=concatenate([p1[p1[:,0]==p][p1[p1[:,0]==p][:,1].argsort()] for p in pux])
    p2=p2.tolist()
    return p2

def lexicographic_sort_yx(p):
    """
function sorts the points list 'p' first with y and then with x smallest to largest 
example:
a=random.random(50)*(20-5)+5
b=random.random(50)*(30-10)+10
p=array([a.round(1),b.round(1)]).transpose(1,0)
fo(f'''
p2={lexicographic_sort_yx(p)};
color("blue")for(i=[0:len(p2)-1])translate(p2[i])text(str(i),.3);
p_lineo(p2,.05);
color("magenta")points(p2,.2);
    ''')
    """
    p=array(p)
    p1=p[p[:,1].argsort()]
    puy=unique(p1[:,1])
    p2=concatenate([p1[p1[:,1]==p][p1[p1[:,1]==p][:,0].argsort()] for p in puy])
    p2=p2.tolist()
    return p2
    
def lexicographic_seg_sort_xy(sec1):
    """
sorts segments in lexicographic order
see the animation below
example:
a=offset_segv(circle(10,s=5),2)
b=lexicographic_seg_sort_xy(a)
for i in range(12):
    c=b[i%4]
    fo(f'''
    color("blue") for(p={a}) p_line3dc(p,.2);
    color("magenta") p_line3d({c},.5);
    ''')
    time.sleep(1)
    
    """

    a=lexicographic_sort_xy(array(sec1)[:,0])
    sec2=[[[p1,p2] if p1[0]<p2[0] else [p2,p1] for (p1,p2) in sec1 if (p2==p) or (p1==p)] for p in a]
    sec2=remove_extra_points(concatenate(sec2))
    return sec2



def rationalise_path(path):
    """
removes all the points which are in straight line
example:
a=m_points1_o(square(10),10)
b=translate([15,0], rationalise_path(a))
fo(f'''
color("blue") for(p={[a]}) points(p,.5);
color("magenta") points({b},.5);
''')
    """
    a=seg(path)[:-1]
    b=[line_as_unit_vector(p) for p in a]
    c=[path[0]]
    for i in range(1,len(a)):
        if norm(a_(b[i])-a_(b[i-1]))>.001:
            c.append(path[i])
    return c+[path[-1]]

def cir_3p_3d(points,s=20):
    """
draws a circle through the 3 points list
's' is the number of segments of the circle
example:
a=[[0,0,0],[10,0,0],[2,5,4]]
b=cir_3p_3d(a,s=50)
fo(f'''
color("blue") for(p={[a]}) points(p,.5);
color("magenta") p_line3dc({b},.2);
''')
    """
    points=rot('x.0001',points)
    n1=array(nv(points))
    a1=cross(n1,[0,0,-1])
    t1=r2d(arccos(n1@[0,0,-1]))
    sec1=translate(-array(points).mean(0),points)
    sec2=c3t2(axis_rot(a1,sec1,t1))
    l1=len(sec2)
    p0,p1,p2=[sec2[0],sec2[int(l1/3)],sec2[int(l1*2/3)]]
    cir1=cir_3p(p0,p1,p2,s=s)
    cir1=translate(array(points).mean(0),axis_rot(a1,cir1,-t1))
    return cir1

circle_3points_3d=cir_3p_3d

def cp_cir_3d(cir):
    """
center point of circle with atleast 3 known list of 'points' in 3d space
a=[[0,0,0],[10,0,0],[2,5,4]]
b=cir_3p_3d(a,s=50)
c=cp_cir_3d(b)
fo(f'''
color("blue") for(p={[a]}) points(p,.5);
color("magenta") p_line3dc({b},.2);
color("cyan") points({[c]},.5);
''')
    """
    n1=array(nv(cir))
    a1=cross(n1,[0.0000001,0.0000001,-1])
    t1=r2d(arccos(n1@[0.0000001,0.0000001,-1]))
    sec1=translate(-array(cir).mean(0),cir)
    sec2=c3t2(axis_rot(a1,sec1,t1))
    l1=len(sec2)
    p0,p1,p2=[sec2[0],sec2[int(l1/3)],sec2[int(l1*2/3)]]
    cp=cp_3p(p0,p1,p2)
    cp=translate(array(cir).mean(0),axis_rot(a1,[cp],-t1))[0]
    return cp

center_circle3d=cp_cir_3d
center_arc3d=cp_cir_3d


def tangents_along_path(path,scale=1):
    p1=array(seg(path))
    p2=array(path)
    v1=array([(p[1]-p[0])/norm(p[1]-p[0]) for p in p1])
    t_v=array([ (v1[-1]+v1[i])/2 if i==0 else
         (v1[i-1]+v1[i])/2
        for i in range(len(p1))])

    n_v=array([ cross(p2[i]-p2[-1],p2[i+1]-p2[i]) if i==0 else
         cross(p2[i]-p2[i-1],p2[i+1]-p2[i]) if i<len(p2)-1 else
         cross(p2[i]-p2[i-1],p2[0]-p2[i])
        for i in range(len(p2))])
    o_v=array([cross(n_v[i],t_v[i]) for i in range(len(t_v))])

    t_v=t_v/norm(t_v,axis=1).reshape(-1,1)
    n_v=n_v/norm(n_v,axis=1).reshape(-1,1)
    o_v=o_v/norm(o_v,axis=1).reshape(-1,1)


    t_v1=array([array([p2[i],p2[i]+t_v[i]*scale]).tolist() for i in range(len(p2))])
    n_v1=array([array([p2[i],p2[i]+n_v[i]*scale]).tolist() for i in range(len(p2))])
    o_v1=array([array([p2[i],p2[i]+o_v[i]*scale]).tolist() for i in range(len(p2))])
    
    return t_v1.tolist()

def normals_along_path(path,scale=1):
    p1=array(seg(path))
    p2=array(path)
    v1=array([(p[1]-p[0])/norm(p[1]-p[0]) for p in p1])
    t_v=array([ (v1[-1]+v1[i])/2 if i==0 else
         (v1[i-1]+v1[i])/2
        for i in range(len(p1))])

    n_v=array([ cross(p2[i]-p2[-1],p2[i+1]-p2[i]) if i==0 else
         cross(p2[i]-p2[i-1],p2[i+1]-p2[i]) if i<len(p2)-1 else
         cross(p2[i]-p2[i-1],p2[0]-p2[i])
        for i in range(len(p2))])
    o_v=array([cross(n_v[i],t_v[i]) for i in range(len(t_v))])

    t_v=t_v/norm(t_v,axis=1).reshape(-1,1)
    n_v=n_v/norm(n_v,axis=1).reshape(-1,1)
    o_v=o_v/norm(o_v,axis=1).reshape(-1,1)


    t_v1=array([array([p2[i],p2[i]+t_v[i]*scale]).tolist() for i in range(len(p2))])
    n_v1=array([array([p2[i],p2[i]+n_v[i]*scale]).tolist() for i in range(len(p2))])
    o_v1=array([array([p2[i],p2[i]+o_v[i]*scale]).tolist() for i in range(len(p2))])
    
    return n_v1.tolist()

def orthos_along_path(path,scale=1):
    p1=array(seg(path))
    p2=array(path)
    v1=array([(p[1]-p[0])/norm(p[1]-p[0]) for p in p1])
    t_v=array([ (v1[-1]+v1[i])/2 if i==0 else
         (v1[i-1]+v1[i])/2
        for i in range(len(p1))])

    n_v=array([ cross(p2[i]-p2[-1],p2[i+1]-p2[i]) if i==0 else
         cross(p2[i]-p2[i-1],p2[i+1]-p2[i]) if i<len(p2)-1 else
         cross(p2[i]-p2[i-1],p2[0]-p2[i])
        for i in range(len(p2))])
    o_v=array([cross(n_v[i],t_v[i]) for i in range(len(t_v))])

    t_v=t_v/norm(t_v,axis=1).reshape(-1,1)
    n_v=n_v/norm(n_v,axis=1).reshape(-1,1)
    o_v=o_v/norm(o_v,axis=1).reshape(-1,1)


    t_v1=array([array([p2[i],p2[i]+t_v[i]*scale]).tolist() for i in range(len(p2))])
    n_v1=array([array([p2[i],p2[i]+n_v[i]*scale]).tolist() for i in range(len(p2))])
    o_v1=array([array([p2[i],p2[i]+o_v[i]*scale]).tolist() for i in range(len(p2))])
    
    return o_v1.tolist()


def l_sec_ip(line,sec):
    """
calculates the intersection point between a line and a closed loop section
example:
a=square(10)
b=point_vector([-5,-5],[4,8])
c=l_sec_ip(b,a)
fo(f'''
color("blue") p_line3dc({a},.2);
color("magenta") p_line3dc({b},.2);
color("cyan") points({c},.5);
''')
    """
    l1=array(line)
    s1=array(seg(sec))
    v1=l1[1]-l1[0]
    p_l=[]
    for p in s1:
        v2=p[1]-p[0]

        iim=array([v1,-v2]).transpose(1,0)+.00001
        im=inv(iim)
        px=p[0]-l1[0]
        t2=(im@px)[1]
        pnts=p[0]+v2*t2
        if 0<=t2<=1:
            p_l.append(pnts.tolist())
    try:
        p_l=l_(a_(p_l)[cKDTree([line[0]]).query(p_l)[0].argsort()])
    except:
        p_l=[]
    return p_l

line_section_ip=l_sec_ip


def l_sec_ip_3d(sec,line):
    """
    same as function l_sec_ip but for 3d coordinates.
    line ans section should be in the same plane
    """
    n1=array(nv(sec))
    a1=cross(n1,[0,0,-1])
    t1=r2d(arccos(n1@[0,0,-1]))
    sec1=translate(-array(sec).mean(0),sec)
    sec2=c3t2(axis_rot(a1,sec1,t1))
    line1=translate(-array(sec).mean(0),line)
    line2=c3t2(axis_rot(a1,line1,t1))
    l1=len(sec2)
    p0,p1,p2=[sec2[0],sec2[int(l1/3)],sec2[int(l1*2/3)]]
    pnts=line_section_ip(line2,sec2)
    pnts=translate(array(sec).mean(0),axis_rot(a1,pnts,-t1)) if pnts!=[] else []
    return pnts

line_section_ip_3d=l_sec_ip_3d

def path_offset_n(sec,r):
    sec=flip(sec) if cw(sec)==1 else sec
    r=round(r,3)
    sec1=offset_segv(sec,r)[:-1]
    s=offset_points(sec,r)[:-1]+[sec1[-1][1]]
    a=s_int1(sec1)
    if a!=[]:
        sec2=a+s
        sec2=array(sec2)
        clean=cs1(sec,abs(r)-.01)[:-1]
        clean1=[p[1:]+[p[0]] for p in clean]
        m,n,_=array(clean).shape
        o,_=sec2.shape
        v1=array([[[1,0]]*n]*m)
        v2=array(clean1)-array(clean)
        iim=array([v1,-v2]).transpose(1,2,0,3).transpose(0,1,3,2)+[0,.00001]
        im=array([pinv(iim)]*o)
        p=(array(clean)[:,:,None]-sec2).transpose(2,0,1,3)
        t=einsum('ijklm,ijkm->ijkl',im,p)
        decision1=((t[:,:,:,0]>=0)&(t[:,:,:,1]>=0)&(t[:,:,:,1]<=1))
        sec3=sec2[(decision1.sum(2)==1).any(1)]
        sec4=sort_points(sec,exclude_points(sec2,sec3))
    else:
        sec4=s
    return sec4

path_offset_new=path_offset_n

def faces(l:int,m:int):
    """
    calculate the faces for the vertices with shape l x m with first and the last end closed
    """
    n1=arange(m,dtype=int)
    n2=array([[[[(j+1)+i*m,j+i*m,j+(i+1)*m],[(j+1)+i*m,j+(i+1)*m,(j+1)+(i+1)*m]] \
             if j<m-1 else \
             [[0+i*m,j+i*m,j+(i+1)*m],[0+i*m,j+(i+1)*m,0+(i+1)*m]] \
                 for j in range(m)] for i in range(l-1)],dtype=int).reshape(-1,3)
    n3=(array(flip(arange(m)),dtype=int)+(l-1)*m)
    n=[n1.tolist()]+n2.tolist()+[n3.tolist()]
    return n

def faces_1(l,m):
    """
    calculate the faces for the vertices with shape l x m 
    """
    n=array([[[[(j+1)+i*m,j+i*m,j+(i+1)*m],[(j+1)+i*m,j+(i+1)*m,(j+1)+(i+1)*m]] \
             if j<m-1 else \
             [[0+i*m,j+i*m,j+(i+1)*m],[0+i*m,j+(i+1)*m,0+(i+1)*m]] \
                 for j in range(m)] for i in range(l-1)]).reshape(-1,3).tolist()

    return n



def faces_2(l,m):
    """
    returns the faces for the vertices with shape l x m with first and the last end open
    """
    return concatenate(faces(l,m)[1:-1]).tolist()

def faces_3(l,m):
    """
    calculate the faces for the vertices with shape l x m with first and the last end open and faces flipped
    """
    n=array([[[[(j+1)+i*m,j+i*m,j+(i+1)*m],[(j+1)+i*m,j+(i+1)*m,(j+1)+(i+1)*m]] \
             if j<m-1 else \
             [[0+i*m,j+i*m,j+(i+1)*m],[0+i*m,j+(i+1)*m,0+(i+1)*m]] \
                 for j in range(m-1)] for i in range(l-1)]).reshape(-1,3).tolist()

    return [flip(p) for p in n]



def prism_center(sol):
    """
calculates the center of the prism or solid object, may not be the mean.
This calculates the center of the bounding box for a solid.
    
    """
    x_max=array(sol)[:,:,0].max()
    x_min=array(sol)[:,:,0].min()

    y_max=array(sol)[:,:,1].max()
    y_min=array(sol)[:,:,1].min()

    z_max=array(sol)[:,:,2].max()
    z_min=array(sol)[:,:,2].min()

    return l_([array([x_max,x_min]).mean(),array([y_max,y_min]).mean(),array([z_max,z_min]).mean()])
    

    
def align_sol_1(sol):
    """
function to straighten the twists in the path_extruded sections for better alignments
example:
example:
a=circle(5)
b=rot('x90',circle(20))
c=path_extrude_closed(a,b)
d=translate([70,0,0],align_sol_1(c))
fo(f'''
// misaligned solid
{swp_c(c)}
// aligned solid
{swp_c(d)}
''')
    """
    sol1=[sol[0]]
    for i in range(1,len(sol)):
        a=align_sec_1(sol1[i-1],sol[i])
        sol1.append(a[1])
    return sol1
    
def align_sec_1(sec1,sec2):
    """
    function to align 2 3d sections to obtain the non twisted optimised solid
    """
    
    area1=[ norm(array(sec2[i:]+sec2[:i])-array(sec1),axis=1).sum() for i in range(len(sec2)) ]
    i=array(area1).argmin()
    sol2=[sec1]+[sec2[i:]+sec2[:i]]
    return sol2


def convert_3lines2fillet_closed(pnt3,pnt2,pnt1,s=10, orientation=0,r=10,style=2):
    """
Develops a fillet with 3 list of points in 3d space
s: number of segments in the fillet, increase the segments in case finer model is required
example:
l1=c23(circle(10))
l2=offset_3d(l1,3)
l3=translate([0,0,3],l1)
f1=convert_3lines2fillet_closed(l2,l3,l1)
fo(f'''
color("blue") for(p={[l1,l2,l3]}) p_line3d(p,.3);
{swp_c(f1)}
''')
    """
    sol1=convert_3lines2fillet(pnt3,pnt2,pnt1,s=s, orientation=orientation,r=r,style=style)
    sol1=sol1+[sol1[0]]
    return sol1 if orientation==0 else cpo(sol1)[:-1]


    
def gcd(a,b):
    """
    calculates the greatest common divisor of 2 numbers 'a','b'
    """
    for _ in range(max([a,b])):
        if a>b:
            a=a-b
        elif b>a:
            b=b-a
        elif a==b:
            break
    return a

def lcm(a,b):
    """
    calculates the least common multiple of 2 numbers 'a','b'
    """
    return a*b/gcd(a,b)


def axis_rot_o(axis,solid,angle):
    """
    rotate a solid around an axis considering the solid is centered at origin
    """
    s_1=len(array(solid).shape)
    cp1=array(solid).mean(0) if s_1==2 else array(prism_center(solid))
    solid=translate(-cp1,solid)
    return translate(cp1,solid@arot(axis,angle))


    
def edges(l,m):
    return array([[[i+j*m,i+(j+1)*m] for j in range(l-1)] for i in range(m)]).reshape(-1,2)




def i_p_n(px,sol1):
    tri=array(ip_triangle(px,sol1))
    p0,p1,p2=tri[:,0],tri[:,1],tri[:,2]
    p01,p02=p1-p0,p2-p0
    v3=cross(p01,p02)
    v3=v3/norm(v3,axis=1).reshape(-1,1)
    return v3


def arc_2p_3d(n1,p0,p1,r,cw=1,s=20):
    """
draws an arc through 2 points 
n1: normal vector to define plane on which the arc will be drawn
r: radius of the arc
cw: '1' stands for clockwise and '-1'stands for counter-clockwise
's' is the number of segments of the circle
example:
n1=[0,1,0]
p0,p1=[[0,0,0],[10,0,0]]
a=arc_2p_3d(n1,p0,p1,r=8,cw=1,s=50)
fo(f'''
color("blue") p_line3d({a},.2);
color("magenta") p_line3d({vector2line(n1,5)},.2);
%{swp_surf(plane(n1,[25,25]))}
''')
    """
    n1=array(n1)
    a1=cross(n1,[0,0,-1])+[.000001,.0000001,0]
    t1=r2d(arccos(n1@[0,0,-1]))
    sec1=translate(-array([p0,p1]).mean(0),[p0,p1])
    sec2=c3t2(axis_rot(a1,sec1,t1))
    pa,pb=sec2
    arc1=arc_2p(pa,pb,r,cw,s=s)
    arc1=translate(array([p0,p1]).mean(0),axis_rot(a1,arc1,-t1))
    return arc1

def arc_long_2p_3d(n1,p0,p1,r,cw=1,s=20):
    """
draws a long arc through 2 points 
n1: normal vector to define plane on which the arc will be drawn
r: radius of the arc
cw: '1' stands for clockwise and '-1'stands for counter-clockwise
's' is the number of segments of the circle
example:
n1=[0,1,0]
p0,p1=[[0,0,0],[10,0,0]]
a=arc_long_2p_3d(n1,p1,p0,r=8,cw=1,s=50)
fo(f'''
color("blue") p_line3d({a},.2);
color("magenta") p_line3d({vector2line(n1,5)},.2);
%{swp_surf(plane(n1,[25,25]))}
''')
    """
    n1=array(n1)
    a1=cross(n1,[0,0,-1])
    t1=r2d(arccos(n1@[0,0,-1]))
    sec1=translate(-array([p0,p1]).mean(0),[p0,p1])
    sec2=c3t2(axis_rot(a1,sec1,t1))
    pa,pb=sec2
    arc1=arc_long_2p(pa,pb,r,cw,s=s)
    arc1=translate(array([p0,p1]).mean(0),axis_rot(a1,arc1,-t1))
    return arc1
    
def arc_2p_3d_cp(n1,p0,p1,r,cw=1):
    """
calculates the center point of the circle drawn through 2 points 
n1: normal vector to define plane on which the arc/ circle drawn
r: radius of the arc/ circle
cw: '1' stands for clockwise and '-1'stands for counter-clockwise
use function cp_cir_3d instead
    """
    n1=array(n1)
    a1=cross(n1,[0,0,-1])
    t1=r2d(arccos(n1@[0,0,-1]))
    sec1=translate(-array([p0,p1]).mean(0),[p0,p1])
    sec2=c3t2(axis_rot(a1,sec1,t1))
    pa,pb=sec2
    cp=arc_2p_cp(pa,pb,r,cw)
    cp=translate(array([p0,p1]).mean(0),axis_rot(a1,[cp],-t1))[0]
    return cp

    
def axis_rot_1(sol,ax1,loc1,theta):
    """
rotate a solid on any pivot point 'loc1' with axis of rotation 'ax1' by an angle 'theta'
example:
c1=cylinder(r=5,h=20)
ax1=[1,1,0]
loc=[0,0,20]
l1=line2length(point_vector(loc,ax1),10)
c2=axis_rot_1(c1,ax1,loc,-90)
fo(f'''
//axis to rotate
color("magenta") p_line3d({l1},.3);
//original cylinder
%{swp(c1)}
//rotated cylinder
color("blue",.5){swp(c2)}
''')
    """
    s_1=len(array(sol).shape)
    c1=array(sol).mean(0) if s_1==2 else array(prism_center(sol))
    loc1=array(loc1)
    c2=c1-loc1
    s1=translate(-c1,sol)
    s1=translate(c2,s1)
    s2=s1@arot(ax1,theta)
    s2=translate(-c2,s2)
    s2=translate(c1,s2)
    return s2


    
def path2path1(path1,path):
    """
function to match the points of path1 with path
i.e. path1 is independent variable and path is dependent variable
example:
p1=c23(cr2dt([[-5,-5,2],[10,0,2],[0,10,2],[-10,0,2]],10))
p2=translate([0,0,10],cr2dt([[-5,-5],[10,0],[0,10]]))
p2=path2path1(p1,p2)
fo(f'''
color("blue",.2) p_line3d({p1},.2);
color("magenta",.2) p_line3d({p2},.2);
color("cyan") for(p={[p1,p2]}) points(p,.4);
''')
    """
    v=[p[1]-p[0] for p in array(seg(path)[:-1])]
    l=[l_len(p) for p in seg(path)[:-1]]
    c=array(l).cumsum().tolist()
    e=[l_len(p) for p in seg(path1)[:-1]]
    f=array(e).cumsum()
    d=[l_lenv_o(path)/f[-1]*p for p in f[:-1]]
    p_rev=[]
    for i in range(len(c)):
        for j in range(len(d)):
            if c[i]>d[j]:
                t=d[j]/l[i] if i==0 else (d[j]-c[i-1])/l[i]
                px=array(path[i])+array(v[i])*t
                p_rev.append(px.tolist())
                d[j]=c[-1]+1
    p_rev=[path[0]]+p_rev+[path[-1]]
    return p_rev

match_points_of_a_path_to_reference_path_open=path2path1
    
def path2path1_closed(path1,path):
    """
function to match the points of path1 with path. Both the paths are closed loop
i.e. path1 is independent variable and path is dependent variable
example:
p1=c23(cr2dt([[-5,-5,2],[10,0,2],[0,10,2],[-10,0,2]],10))
p2=translate([0,0,10],circle(5))
p2=path2path1_closed(p1,p2)
fo(f'''
color("blue",.2) p_line3dc({p1},.2);
color("magenta",.2) p_line3dc({p2},.2);
color("cyan") for(p={[p1,p2]}) points(p,.4);
''')
    """
    v=[p[1]-p[0] for p in array(seg(path))]
    l=[l_len(p) for p in seg(path)]
    c=array(l).cumsum().tolist()
    e=[l_len(p) for p in seg(path1)]
    f=array(e).cumsum()
    d=[l_lenv(path)/f[-1]*p for p in f]
    p_rev=[]
    for i in range(len(c)):
        for j in range(len(d)):
            if c[i]>d[j]:
                t=d[j]/l[i] if i==0 else (d[j]-c[i-1])/l[i]
                px=array(path[i])+array(v[i])*t
                p_rev.append(px.tolist())
                d[j]=c[-1]+1
    p_rev=[path[0]]+p_rev
    return p_rev[:len(path1)]

match_points_of_a_path_to_reference_path_closed=path2path1_closed

def ip_triangle(ip,sol1):
    """
function to find the triangles on the solid 'sol1' where the intersection points list 'ip' lies
example:
sec=pts([[0,0],[5,0],[0,5],[-5,0]])
sec=m_points1(sec,7)
path=corner_radius(pts1([[-2.49,0],[2.49,0,0.25],[0,5,0.25],[-2.49,0]]),10)
sol1=prism(sec,path)
v1=[2,0,8]
sec1=axis_rot([0,0,1],circle(1.5,s=100),90)
sol2=o_solid(v1,sec1,10,-1.5,-2.5,0)
p1=ip_sol2sol(sol1,sol2)
p2=ip_sol2sol(sol1,sol2,-1)
p3=flip(p1)+p2
tri_1=ip_triangle(p3,sol1)
fo(f'''
%{swp(sol1)}
%{swp(sol2)}
// intersection points
color("magenta")points({p3},.05);
// triangles where the intersection points fall
color("blue")for(p={tri_1})p_line3dc(p,.04,rec=1);
''')
    """
    v,f1=vnf2(sol1)
    tri=array(v)[array(f1)]
    p0,p1,p2=tri[:,0],tri[:,1],tri[:,2]
    p01,p02=p1-p0,p2-p0
    n1=cross(p01,p02)
    n1=n1/(norm(n1,axis=1).reshape(-1,1)+.00001)
    tri=tri[~((n1==[0,0,0]).all(1))]
    n1=n1[~((n1==[0,0,0]).all(1))]
    p0,p1,p2=tri[:,0],tri[:,1],tri[:,2]
    p01,p02=p1-p0,p2-p0
    la=array(ip)
    lab=n1

    iim=array([lab,-p01,-p02]).transpose(1,0,2).transpose(0,2,1)
    im=inv(iim)

    x=einsum('jkl,ijl->ijk',im,(p0-la[:,None]))
    t=x[:,:,0]
    u=x[:,:,1]
    v=x[:,:,2]
    decision=(t>=-0.01)&(t<=1.01)&(u>=-0.01)&(u<=1.01)&(v>=-0.01)&(v<=1.01)&((u+v)<=1.01)
    tri_1=array([tri[decision[i]][0] for i in range(len(ip))]).tolist()

    return tri_1



def o_3d(i_p,sol,r,o=0,f=1,closed=0):
    """
function to offset the intersection points 'i_p' on a solid 'sol' by distance 'r'. option 'o' can have values '0' or '1' and changes the direction of offset.
for closed loop path set closed=1
example:
sec=pts([[0,0],[5,0],[0,5],[-5,0]])
sec=m_points1(sec,7)
path=corner_radius(pts1([[-2.49,0],[2.49,0,0.25],[0,5,0.25],[-2.49,0]]),10)
sol1=prism(sec,path)
v1=[2,0,8]
sec1=axis_rot([0,0,1],circle(1.5,s=100),90)
sol2=o_solid(v1,sec1,10,-1.5,-2.5,0)
p1=ip_sol2sol(sol1,sol2)
p2=ip_sol2sol(sol1,sol2,-1)
p3=flip(p1)+p2
tri_1=ip_triangle(p3,sol1)
p4=o_3d(p3,sol1,-.5)
p5=o_3d(p3,sol2,.5)

fo(f'''
%{swp(sol1)}
%{swp(sol2)}
// intersection points
color("magenta")points({p3},.05);
color("magenta") p_line3dc({p3},.03);
// offset 3d of intersection points on sol1
color("cyan") points({p4},.05);
color("cyan") p_line3dc({p4},.03);
// offset 3d of intersection points on sol2
color("green") points({p5},.05);
color("green") p_line3dc({p5},.03);
''')
    """
    a=i_p_n(i_p,sol)
    if closed==0:
        b=i_p_t_o(i_p)
    elif closed==1:
        b=i_p_t(ip)
    if o==0:
        c=array(i_p)+cross(b,a)*r
    elif o==1:
        c=array(i_p)+cross(a,b)*r
    s=array([c+a*r*f,c-a*r*f])
    i_p1=ip_sol2sol(sol,s)
    # i_p1=[p[0] for p in i_p1]
    return i_p1

offset_intersection_line_on_solid=o_3d

def o_3d_rev(i_p,sol,r,o=0,closed=0,type=0,dist=0,vx=[],edges_closed=1,cg=0):
    """
function to offset the intersection points 'i_p' on a solid 'sol' by distance 'r'. option 'o' can have values '0' or '1' and changes the direction of offset.
for closed loop path set closed=1
o: set '1' to shift the line on other side
type: set '1' if prism_center lies outside the solid. type can also be set to '2'.
for using type=2, refer function definition psos_v_1 and also checkout the video explanation
prism_center is the center of the bounding box of the solid
example:
sec=pts([[0,0],[5,0],[0,5],[-5,0]])
sec=m_points1(sec,7)
path=corner_radius(pts1([[-2.49,0],[2.49,0,0.25],[0,5,0.25],[-2.49,0]]),10)
sol1=prism(sec,path)
v1=[2,0,8]
sec1=axis_rot([0,0,1],circle(1.5,s=100),90)
sol2=o_solid(v1,sec1,10,-1.5,-2.5,0)
p1=ip_sol2sol(sol1,sol2)
p2=ip_sol2sol(sol1,sol2,-1)
p3=flip(p1)+p2
tri_1=ip_triangle(p3,sol1)
p4=o_3d_rev(p3,sol1,-.5)
p5=o_3d_rev(p3,sol2,.5)

fo(f'''
%{swp(sol1)}
%{swp(sol2)}
// intersection points
color("magenta")points({p3},.05);
color("magenta") p_line3dc({p3},.03);
// offset 3d of intersection points on sol1
color("cyan") points({p4},.05);
color("cyan") p_line3dc({p4},.03);
// offset 3d of intersection points on sol2
color("green") points({p5},.05);
color("green") p_line3dc({p5},.03);
''')
    """
    if edges_closed==1:
        sol=cpo(cpo(sol)+[cpo(sol)[0]])
    elif edges_closed==0:
        sol=sol
    a=i_p_n(i_p,sol)
    if closed==0:
        b=i_p_t_o(i_p)
    elif closed==1:
        b=i_p_t(ip)
    if o==0:
        c=array(i_p)+cross(b,a)*r
    elif o==1:
        c=array(i_p)+cross(a,b)*r

    if type==0:
        i_p1=psos_v(sol,[l_(c)],prism_center(sol),dist=(abs(r) if dist==0 else dist))[0]
    elif type==1:
        s=[i_p,l_(c)]
        i_p1=psos_n_b(sol,s,dist=(abs(r) if dist==0 else dist))[-1]
    elif type==2:
        if cg==0:
            i_p1=psos_v_1(sol,[l_(c)],prism_center(sol),vx,dist=(abs(r) if dist==0 else dist))[0]
        elif cg==1:
            i_p1=psos_v_1(sol,[l_(c)],cog(sol),vx,dist=(abs(r) if dist==0 else dist))[0]
            
    return i_p1

offset_intersection_line_on_solid_rev=o_3d_rev

def o_3d_surf(i_p,sol,r,o=0,f=1,closed=0):
    """
function to offset the intersection points 'i_p' on a solid 'sol' by distance 'r'. option 'o' can have values '0' or '1' and changes the direction of offset
for closed loop path, set closed=1
example:
l1=rot('x90',sinewave(50,3,2,50))
l2=rot('x90z90',cosinewave(50,3,2,50))
s1=surface_from_2_waves(l1,l2,2)
c1=translate([0,0,-10],linear_extrude(circle(5,[25,25]),30))
l3=ip_surf2sol(s1,c1)
l4=o_3d_surf(l3,s1,-3)
l5=o_3d(l3,c_(c1),3)
fo(f'''
%{swp_surf(s1)}
%{swp(c1)}
color("blue") p_line3dc({l3},.3);
color("magenta") p_line3dc({l4},.3);
color("cyan") p_line3dc({l5},.3);
''')
    """
    a=i_p_n_surf(i_p,sol)
    if closed==0:
        b=i_p_t_o(i_p)
    elif closed==1:
        b=i_p_t(ip)
    if o==0:
        c=array(i_p)+cross(b,a)*r
    elif o==1:
        c=array(i_p)+cross(a,b)*r
    i_p1=psos_v_2(sol,[c],a)[0]

    return i_p1

offset_intersection_line_on_surface=o_3d_surf
    



def ip_fillet(sol1,sol2,r1,r2,s=20,o=0,type=1,dist=0,vx=[],style=2,edges_closed=1,c=0):
    """
calculates a fillet at the intersection of 2 solids.
r1 and r2 would be same in most of the cases, but the signs can be different depending on which side the fillet is required
r1 is the distance by which intersection line offsets on sol1 and similarly r2 is on sol2 
for type, dist and vx parameters refer function o_3d_rev
example:
s1=sphere(20)
l1=translate([-10,0,0],sinewave(20,2,2,50))
s2=surface_line_vector(l1,[5,5,50])
f1=ip_fillet(s1,s2,2,2)
fo(f'''
%{swp(s1)}
%{swp_surf(s2)}
{swp(f1)}
''')
    """
    sol1=cpo(cpo(sol1)+[cpo(sol1)[0]])
    p1=ip_sol2sol(sol1,sol2,o)
    if style==1:
        p2=i_p_p(sol2,p1,r2/2)
        p3= o_3d_rev(p1,sol1,r1/2,type=type, dist=dist,vx=vx, edges_closed=edges_closed, cg=c)
        fillet1=convert_3lines2fillet(p3,p2,p1,s=s,r=max([r1,r2]),style=1)
    elif style==2:
        p2=i_p_p(sol2,p1,r2)
        p3= o_3d_rev(p1,sol1,r1,type=type, dist=dist,vx=vx, edges_closed=edges_closed, cg=c)
        fillet1=convert_3lines2fillet(p3,p2,p1,s=s,style=2)
    
    return fillet1

def ip_fillet_closed(sol1,sol2,r1,r2,s=20,o=0,type=1,dist=0,vx=[],style=2,edges_closed=1,c=0):
    """
calculates a fillet at the intersection of 2 solids.
r1 and r2 would be same in most of the cases, but the signs can be different depending on which side the fillet is required
r1 is the distance by which intersection line offsets on sol1 and similarly r2 is on sol2 
for type, dist and vx parameters refer function o_3d_rev
example:
s1=sphere(20)
c1=rot('y45',cylinder(r=5,h=50,s=50))
f1=ip_fillet_closed(s1,c1,-2,2)
fo(f'''
{swp(s1)}
{swp(c1)}
{swp_c(f1)}
''')
    """
    sol1=cpo(cpo(sol1)+[cpo(sol1)[0]])
    p1=ip_sol2sol(sol1,sol2,o)
    p2=i_p_p(sol2,p1,r2)
    p3=o_3d_rev(p1,sol1,r1,type=type,dist=dist,vx=vx,edges_closed=edges_closed,cg=c)
    fillet1=convert_3lines2fillet_closed(p3,p2,p1,s=s,style=style)
    
    return fillet1


def ip_fillet_surf(surf,sol,r1,r2,s=20,type=1,dist=0,vx=[],style=2,edges_closed=0,c=0):
    """
calculates a fillet at the intersection of surface with solid.
r1 and r2 would be same in most of the cases, but the signs can be different depending on which side the fillet is required
r1 is the distance by which intersection line offsets on surf and similarly r2 is on sol 
example:
l0=translate([-25,0,0],rot('x90',cosinewave(50,3,2,100)))
s1=surface_line_vector(l0,[0,50,0],1)
l1=translate([-10,0,-10],sinewave(20,2,2,50))
s2=surface_line_vector(l1,[5,5,50])
f1=ip_fillet(s1,s2,2,2)
fo(f'''
%{swp_surf(s1)}
%{swp_surf(s2)}
{swp(f1)}
''')
    """
        
    p1=ip_surf2sol(surf,sol)
    p2=i_p_p(sol,p1,r2)
    p3=o_3d_rev(p1,surf,r1,type=type,dist=dist,vx=vx,edges_closed=edges_closed,cg=c)
    fillet1=convert_3lines2fillet(p3,p2,p1,s=s,style=style)

    return fillet1

def ip_fillet_surf_closed(surf,sol,r1,r2,s=20,type=1,dist=0,vx=[],style=2,f=1,edges_closed=0,c=0):
    """
calculates a fillet at the intersection of surface with solid.
r1 and r2 would be same in most of the cases, but the signs can be different depending on which side the fillet is required
r1 is the distance by which intersection line offsets on surf and similarly r2 is on sol 
example:
l0=translate([-25,0,0],rot('x90',cosinewave(50,3,2,100)))
s1=surface_line_vector(l0,[0,50,0],1)
s2=translate([0,0,-10],cylinder(r=10,h=50))
f1=ip_fillet_closed(s1,s2,-2,2)
fo(f'''
%{swp_surf(s1)}
%{swp_surf(s2)}
{swp(f1)}
''')
    """
        
    p1=ip_surf2sol(surf,sol)
    p2=i_p_p(sol,p1,r2)
    p3=o_3d_rev(p1,surf,r1,type=type,dist=dist,vx=vx,edges_closed=edges_closed,cg=c)
    fillet1=convert_3lines2fillet_closed(p3,p2,p1,s=s,style=style)

    return fillet1


def ip_random(sol1,sol2):
    v,f1=vnf2(sol1)
    a=array(v)[array(f1)]
    p0,p1,p2=a[:,0],a[:,1],a[:,2]
    p01,p02=p1-p0,p2-p0
    v,f1=vnf2(sol2)
    b=array(v)[array(f1)]
    c=b[:,0]
    d=b[:,1]
    e=b[:,2]
    lcd,lde,lec=d-c,e-d,c-e
    x=zeros(len(p01)*3).reshape(len(p01),3)
    y=zeros(len(lcd)*3).reshape(len(lcd),3)

    t=einsum('jk,ijk->ij',cross(p01,p02),c[:,None]-p0[None,:])/(einsum('ik,jk->ij',-lcd,cross(p01,p02))+.00001)
    u=einsum('ijk,ijk->ij',cross(p02[None,:],-lcd[:,None]),c[:,None]-p0[None,:])/(einsum('ik,jk->ij',-lcd,cross(p01,p02))+.00001)
    v=einsum('ijk,ijk->ij',cross(-lcd[:,None],p01[None,:]),c[:,None]-p0[None,:])/(einsum('ik,jk->ij',-lcd,cross(p01,p02))+.00001)
    dcn=(t>=0)&(t<=1)&(u>=0)&(u<=1)&(v>=0)&(v<=1)&(u+v<=1)
    i_p=((c[:,None]+x[None,:])+einsum('ijk,ij->ijk',(lcd[:,None]+x[None,:]),t))[dcn]
    
    return i_p.tolist()
    
def ellipse(a,b,s=50):
    """
function to form ellipse
example:
l1=ellipse(a=20,b=10,s=100)
fo(f'''
color("blue") p_line3dc({l1},.3);
''')
    """
    return l_([[a*cos(d2r(i)),b*sin(d2r(i))]  for i in linspace(0,360,s)[:-1]])
    

def sol2path(sol,path):
    """
function to extrude a solid to path
example:
sec=circle(5)
path1=corner_radius(pts1([[-3,0],[3,0,3],[0,15,2],[-2,5,3],[2,5,2],[0,20,3],[-3,0]]),10)
path1=equidistant_path(path1,200)
sol=prism(sec,path1)
path=bspline_open(cr3dt([[10,0,0,0],[10,15,10,8],[0,0,15,6],[20,0,0,0]],20),2,50)
sol2=sol2path(sol,path)
fo(f'''
%{swp(sol)}
%{swp(sol2)}
color("cyan")p_line3d({path},.5);
''')
    """
    sol1=c3t2(sol)
    zpath=[[0,0,p[0][2]] for p in sol]
    path2=path2path1(zpath,path)
    sol2=align_sol_1(path_extrude2msec(sol1,path2))
    return sol2

extrude_a_solid_to_defined_path=sol2path

def comb_list(n):
    """
finds the combination in pairs e.g. there are 4 teams, how many matches should be played so that everyone plays with everyone else.
example:
comb_list(4) => 
array([[0, 1],
       [0, 2],
       [0, 3],
       [1, 2],
       [1, 3],
       [2, 3]])
    """
    n=arange(n)
    a=array([n]*len(n)).transpose(1,0)
    b=array([n]*len(n))
    c=array([a,b]).transpose(1,0,2).transpose(0,2,1)
    d=concatenate([c[i][i+1:] for i in range(len(c))])
    return d

def ip_nv_sol2sol(sol1,sol2):
    """
function finds the intersection point and the normal vector to that intersection point between 2 solids
sol1: solid on who's surface intersection points needs to be found
sol2: solid which intersects sol1
    
    """
    line=array([ seg(p)[:-1] for p in cpo(sol2)])
    v,f1=vnf2(sol1)
    tri=array(v)[array(f1)]
    line=array([ seg(p)[:-1] for p in cpo(sol2)])
    tri.shape,line.shape
    la,lb=line[:,:,0],line[:,:,1]
    p0,p1,p2=tri[:,0],tri[:,1],tri[:,2]
    lab=lb-la
    p01,p02=p1-p0,p2-p0
    t=einsum('kl,ijkl->ijk',cross(p01,p02),la[:,:,None]-p0)/(einsum('ijl,kl->ijk',(-lab),cross(p01,p02))+.00000)
    u=einsum('ijkl,ijkl->ijk',cross(p02[None,None,:,:],(-lab)[:,:,None,:]),(la[:,:,None,:]-p0[None,None,:,:]))/(einsum('ijl,kl->ijk',(-lab),cross(p01,p02))+.00000)
    v=einsum('ijkl,ijkl->ijk',cross((-lab)[:,:,None,:],p01[None,None,:,:]),(la[:,:,None,:]-p0[None,None,:,:]))/(einsum('ijl,kl->ijk',(-lab),cross(p01,p02))+.00000)
    condition=(t>=0)&(t<=1)&(u>=0)&(u<=1)&(v>=0)&(v<=1)&(u+v<1)

    a=(la[:,None,:,None,:]+lab[:,None,:,None,:]*t[:,None,:,:,None])
    n1=-array([cross(p01,p02)/norm(cross(p01,p02),axis=1).reshape(-1,1)]*len(line))[:,None,None,:]
    b=condition[:,None,:,:]
    c,d=[],[]
    for i in range(len(a)):
        c.append(a[i][b[i]].tolist())
        d.append(n1[i][b[i]].tolist())

    p=[p for p in c if p!=[]]
    n1=[d[i] for i in range(len(c)) if c[i]!=[]]
    return [p,n1]


def corner_radius(sec,s=20):
    """
function to create section with corner radiuses. e.g. 
following code has 3 points at [0,0],[10,0] and [7,15] and radiuses of 0.5,2 and 1 respectively,
s=5 represent the number of segments at each corner radius.
example:
sec1=corner_radius(sec=[[0,0,.5],[10,0,2],[7,15,1]],s=5)
fo(f'''
color("blue") p_line3dc({sec1},.3);
''')
    
    """
    r_l=[0 if len(p)==2 else p[2] for p in sec]
    sec=[ [p[0],p[1]] for p in sec]

    p0=[sec[-1]]+sec[:-1]
    p1=sec
    p2=sec[1:]+[sec[0]]

    p0,p1,p2=array([p0,p1,p2])
    v1,v2=p0-p1,p2-p1
    u1,u2=v1/norm(v1,axis=1).reshape(-1,1),v2/norm(v2,axis=1).reshape(-1,1)
    p_o=cwv(p1.tolist())

    theta=[(180-ang_2lineccw(p1[i],p0[i],p2[i]))/2  if p_o[i]==1 else (180-ang_2linecw(p1[i],p0[i],p2[i]))/2 for i in range(len(p0))]
    th0=[theta[-1]]+theta[:-1]
    th1=theta
    c_p1=array([p1[i] if r_l[i]==0 else  
                array(p1[i])+(u1[i]*r_l[i]/cos(d2r(theta[i])))@[[cos(d2r((180-2*theta[i])/2)),sin(d2r((180-2*theta[i])/2))],[-sin(d2r((180-2*theta[i])/2)),cos(d2r((180-2*theta[i])/2))]]
                if p_o[i]==1 else 
                array(p1[i])+(u1[i]*r_l[i]/cos(d2r(theta[i])))@[[cos(d2r((180-2*theta[i])/2)),-sin(d2r((180-2*theta[i])/2))],[sin(d2r((180-2*theta[i])/2)),cos(d2r((180-2*theta[i])/2))]]
                for i in range(len(p0))]).tolist()

    cir_1=[ 
        [p1[i].tolist()]
        if r_l[i]==0 else
        circle(r_l[i],cp=c_p1[i])
        for i in range(len(p0))]

    # radiuses
    r0=[r_l[-1]]+r_l[:-1]
    r1=r_l
    # circles
    c0=[cir_1[-1]]+cir_1[:-1]
    c1=cir_1
    # center points
    cp0=[c_p1[-1]]+c_p1[:-1]
    cp1=c_p1
    # orientations
    p_o_0=[p_o[-1]]+p_o[:-1]
    p_o_1=p_o
    a=[]
    for i in range(len(p0)):
        if r0[i]==0 and r1[i]==0:
            a.append([p1[i].tolist()])
        elif p_o_0[i]==-1 and r0[i]>0 and r1[i]==0:
            a.append([cir_p_t(c0[i],p1[i]),p1[i].tolist()])
        elif p_o_0[i]==1 and r0[i]>0 and r1[i]==0:
            a.append([p_cir_t(p1[i],c0[i]),p1[i].tolist()])
        elif p_o_1[i]==-1 and r0[i]==0 and r1[i]>0:
            a.append([p_cir_t(p0[i],c1[i])])
        elif p_o_1[i]==1 and r0[i]==0 and r1[i]>0:
            a.append([cir_p_t(c1[i],p0[i])])
        elif p_o_0[i]==1 and p_o_1[i]==1 and r0[i]>0 and r1[i]>0:
            if (r0[i]*tan(d2r(th0[i]))+r1[i]*tan(d2r(th1[i])))>norm(p1[i]-p0[i]):
                raise ValueError('radiuses more than acceptable limit')
            else:
                a.append(flip(tctpf(r0[i],r1[i],cp0[i],cp1[i])[2:]))
        elif p_o_0[i]==1 and p_o_1[i]==-1 and r0[i]>0 and r1[i]>0:
            a.append(tcct(r0[i],r1[i],cp0[i],cp1[i],1))
        elif p_o_0[i]==-1 and p_o_1[i]==1 and r0[i]>0 and r1[i]>0:
            a.append(tcct(r0[i],r1[i],cp0[i],cp1[i],-1))
        elif p_o_0[i]==-1 and p_o_1[i]==-1 and r0[i]>0 and r1[i]>0:
            if (r0[i]*tan(d2r(th0[i]))+r1[i]*tan(d2r(th1[i])))>norm(p1[i]-p0[i]):
                raise ValueError('radiuses more than acceptable limit')
            else:
                a.append(tctpf(r0[i],r1[i],cp0[i],cp1[i])[:2])

    b=[0]
    for i in range(len(p0)):
        if r_l[i]>0:
            b.append(b[-1]+2)
        else:
            b.append(b[-1]+1)

    b=array(b[1:])-1


    c=concatenate(a).tolist()
    c=c if r_l[-1]==0 else c[1:]+[c[0]]
    d=[]
    for i in range(len(p0)):
        if r_l[i]==0:
            d.append([c[b[i]]])
        else:
            d.append( arc_2p(c[b[i]-1],c[b[i]],r_l[i],p_o[i],s))
            

    d=(concatenate(d).round(8)).tolist()
    d=min_d_points(d,.0001)
    d=c3t2(rot('z.0001',d))
    return d

cr2d=corner_radius

def surround(path,r,s=20):
    """
function to surround a path to create a rounded section
example:
sec1=corner_radius(sec=[[0,0,.5],[10,0,2],[7,15,1]],s=5)
l1=surround(homogenise(sec1,.25),.5)
fo(f'''
color("blue") p_line3d({sec1},.3);
color("magenta") p_line3dc({l1},.3);
''')
    
    """
    a=path_offset_n(path,r)
    b=path_offset_n(path,-r)
    b=flip(b)
    arc1=arc_2p(a[-1],b[0],r,-1,s)
    arc2=arc_2p(b[-1],a[0],r,-1,s)
    sec=a[1:-1]+arc1+b[1:-1]+arc2
    return sec


def oset(sec,r):
    sec0=intersections(offset_segv(sec,r))
    i_p1=s_int1(seg(sec0))
    d=int_seg_list(seg(sec0))
    sec2=array(sec0)
    for i in range(len(d)):
        a=arange(d[i][0]+1)
        b=arange(d[i][0],d[i][1]+1)
        c=arange(d[i][1],len(sec2))
        if (len(a)+len(c))<len(b):
            sec2[a]=array([i_p1[i]]*len(a))
        else:
            sec2[b]=array([i_p1[i]]*len(b))
    sec2=array(sec2).tolist()
            
    return sec2

def arot(v,theta):
    """
    rotation matrix for rotating objects around any arbitrary axis defined by vector 'v'
    follows right hand thumb rule for rotation
    """
    u=v/norm(v)
    return array([
    [cos(d2r(theta))+u[0]**2*(1-cos(d2r(theta))),u[1]*u[0]*(1-cos(d2r(theta)))+u[2]*sin(d2r(theta)),u[2]*u[0]*(1-cos(d2r(theta)))-u[1]*sin(d2r(theta))],
    [u[0]*u[1]*(1-cos(d2r(theta)))-u[2]*sin(d2r(theta)),cos(d2r(theta))+u[1]**2*(1-cos(d2r(theta))),u[2]*u[1]*(1-cos(d2r(theta)))+u[0]*sin(d2r(theta))],
    [u[0]*u[2]*(1-cos(d2r(theta)))+u[1]*sin(d2r(theta)),u[1]*u[2]*(1-cos(d2r(theta)))-u[0]*sin(d2r(theta)),cos(d2r(theta))+u[2]**2*(1-cos(d2r(theta)))]
    ])

def xrot(theta):
    """
    rotation matrix to rotate objects around x-axis
    follows right hand thumb rule for rotation
    """
    return array([
        [1,0,0],
        [0,cos(d2r(theta)),sin(d2r(theta))],
        [0,-sin(d2r(theta)),cos(d2r(theta))]
    ])

def yrot(theta):
    """
    rotation matrix to rotate objects around y-axis
    follows right hand thumb rule for rotation
    """
    return array([
        [cos(d2r(theta)),0,-sin(d2r(theta))],
        [0,1,0],
        [sin(d2r(theta)),0,cos(d2r(theta))]
    ])

def zrot(theta):
    """
    rotation matrix to rotate objects around z-axis
    follows right hand thumb rule for rotation
    """
    return array([[cos(d2r(theta)),sin(d2r(theta)),0],
                 [-sin(d2r(theta)),cos(d2r(theta)),0],
                  [0,0,1]
                 ])

def surface_from_2_waves(p0,p1,amplitude=1):
    """
function to draw surface based on 2 waves perpendicular to each other.
waves are multiplied
example:
l1=rot('x90',sinewave(50,3,2,50))
l2=rot('x90z90',cosinewave(50,3,2,50))
s1=surface_from_2_waves(l1,l2,2)
fo(f'''
{swp_surf(s1)}
''')
    """
    p2=array([[[i[0],j[1],(i@j)]  for j in array(p1)]  for i in array(p0)])
    a=p2[:,:,2].max()
    p2=array([[[j[0],j[1],j[2]/a*amplitude]  for j in i] for i in p2]).tolist()
    return p2

def surface_from_2_waves_add(p0,p1,amplitude=1):
    """
function to draw surface based on 2 waves perpendicular to each other.
waves are added

example:
l1=rot('x90',sinewave(50,3,2,50))
l2=rot('x90z90',cosinewave(50,3,2,50))
s1=surface_from_2_waves_add(l1,l2,2)
fo(f'''
{swp_surf(s1)}
''')
    """
    p2=array([[[i[0],j[1],i[2]+j[2]]  for j in array(p1)]  for i in array(p0)])
    a=p2[:,:,2].max()
    p2=array([[[j[0],j[1],j[2]/a*amplitude]  for j in i] for i in p2]).tolist()
    return p2

def surface_from_2_waves_min(p0,p1,amplitude=1):
    """
function to draw surface based on 2 waves perpendicular to each other.
maximum point in the 2 waves will be considered
example:
l1=rot('x90',sinewave(50,3,2,50))
l2=rot('x90z90',cosinewave(50,3,2,50))
s1=surface_from_2_waves_min(l1,l2,2)
fo(f'''
{swp_surf(s1)}
''')
    """
    p2=array([[[i[0],j[1],a_([i[2],j[2]]).min()]  for j in array(p1)]  for i in array(p0)])
    a=p2[:,:,2].max()
    p2=array([[[j[0],j[1],j[2]/a*amplitude]  for j in i] for i in p2]).tolist()
    return p2

def surface_from_2_waves_max(p0,p1,amplitude=1):
    """
function to draw surface based on 2 waves perpendicular to each other.
maximum point in the 2 waves will be considered
example:
l1=rot('x90',sinewave(50,3,2,50))
l2=rot('x90z90',cosinewave(50,3,2,50))
s1=surface_from_2_waves_max(l1,l2,2)
fo(f'''
{swp_surf(s1)}
''')
    """
    p2=array([[[i[0],j[1],max([i[2],j[2]])]  for j in array(p1)]  for i in array(p0)])
    a=p2[:,:,2].max()
    p2=array([[[j[0],j[1],j[2]/a*amplitude]  for j in i] for i in p2]).tolist()
    return p2

def surface_from_2_waves_norm(p0,p1,amplitude=1):
    """
function to draw surface based on 2 waves perpendicular to each other.
norm of the 2 waves will be considered
example:
l1=rot('x90',sinewave(50,3,2,50))
l2=rot('x90z90',cosinewave(50,3,2,50))
s1=surface_from_2_waves_norm(l1,l2,2)
fo(f'''
{swp_surf(s1)}
''')
    """
    p2=array([[[i[0],j[1],norm([i[2],j[2]])]  for j in array(p1)]  for i in array(p0)])
    a=p2[:,:,2].max()
    p2=array([[[j[0],j[1],j[2]/a*amplitude]  for j in i] for i in p2]).tolist()
    return p2


def fillet_intersection_lines(l1,l2,r,s=10,side=0):
    """
function calculates the fillet at intersection between 2 lines
'l1' and 'l2'
r: radius of fillet
s: segments of fillet
example:
l1=point_vector([-5,-5],[10,10])
l2=point_vector([-5,20],[10,-10])
l3=fillet_intersection_lines(l1,l2,r=3)
fo(f'''
color("blue") p_line3d({l1},.3);
color("cyan") p_line3d({l2},.3);
color("magenta") p_line3d({l3},.3);
''')
    """
    p0=i_p2d(l1,l2)
    # p1=l1[a_([l_len([p0,p]) for p in l1]).argsort()[-1]]
    # p2=l2[a_([l_len([p0,p]) for p in l2]).argsort()[-1]]
    p1=l1[0]
    p2=l2[-1]
    clock=cw([p1,p0,p2])
    if side==0:
        p3=path_offset([p1,p0,p2],-r)[1] if clock==-1 else \
        path_offset([p1,p0,p2],r)[1]
        p4=vcost1(l1,p3)
        p5=vcost1(l2,p3)
        arc1=arc_2p(p4,p5,r,clock,s)
    if side==1:
        p3=path_offset([p1,p0,p2],r)[1] if clock==-1 else \
        path_offset([p1,p0,p2],-r)[1]
        p4=vcost1(l1,p3)
        p5=vcost1(l2,p3)
        arc1=arc_2p(p4,p5,r,clock,s)
    return arc1
    
def cir_line_tangent(c1,l1,side=0):
    """
function to draw a line tangent to a circle
c1 is circle
l1 is line
side=0, draws tangent at one side and 1 draws tangent on the other side
example:
c1=circle(10,[-30,0])
l1=point_vector([0,0],[10,3])
l2=cir_line_tangent(c1,l1,side=0)
l3=cir_line_tangent(c1,l1,side=1)
fo(f'''
color("blue") p_line3d({c1},.3);
color("cyan") p_line3d({l1},.3);
color("magenta") p_line3d({l2},.3);
color("grey") p_line3d({l3},.3);
''')
    """
    r=r_arc(c1)
    cp=cp_arc(c1)
    v1=array(l1[1])-array(l1[0])
    theta=ang(v1[0],v1[1]) if side==0 else ang(v1[0],v1[1])+180
    l=l_len(l1)
    p0=[r*cos(d2r(270+theta)),r*sin(d2r(270+theta))]
    p1=(array(p0)+rot2d((theta if side==0 else theta-180),[l,0])).tolist()
    p0,p1= (array([p0,p1])+cp).tolist()
    return [p0,p1]

def spiral_poly(r=1,d=.3,n=4,t=100):
    """
create a spiral polygon
r: initial length of the line
d: increment length every iteration
n: number of sides of the polygon
t: number of turns
example:
l1=spiral_poly(r=1,d=.2,n=5,t=100)
fo(f'''
color("blue") p_line3d({l1},.3);
''')
    """
    theta=360/n
    sec=[[r,0]]
    for i in range(1,t):
        r=r+d
        a=array([r,0])
        sec.append(array(sec[-1])+rot2d(i*theta,a))

    sec=array(sec).tolist()
    return sec

def equate_points(sec,sec1):
    """
function to make the points in 2 sections equal without changing the location of points
    """
    c=array(sec1).shape[-1]
    b=array(sec).shape[-1]
    a=lcm(len(sec),len(sec1))
    sec=array([sec]*int(a/len(sec))).transpose(1,0,2).reshape(-1,b).tolist()
    sec1=array([sec1]*int(a/len(sec1))).transpose(1,0,2).reshape(-1,c).tolist()
    return [sec,sec1]

def pies2(sec,pnts):
    """
    function to find 3d points 'pnts' which are inside an enclosed 2d section 'sec'
    """
    pnts=rev_pnts(sec,pnts)
    if pnts!=[]:
        s8,s4=[sec,pnts]
        p0=array(s4)
        p2=s8
        p3=s8[1:]+[s8[0]]
        p2,p3=array([p2,p3])
        # v1=array([[[1,0]]*len(p2)]*len(p0))
        v1=array([ones(len(p2)),zeros(len(p2))]).transpose(1,0)
        v2=(p3-p2)+.000001
        p=p2-p0[:,:2][:,None]
        im=pinv(array([v1,-v2]).transpose(1,0,2).transpose(0,2,1))
        im=array([im]*len(p0))
        t=einsum('ijkl,ijl->ijk',im,p)

        s10=[p0[i] for i in range(len(p0)) if \
                t[i][(t[i][:,0]>=0)&(t[i][:,1]>=0)&(t[i][:,1]<=1)].shape[0]%2 \
             ==1]
        return array(s10).tolist()

def sinewave(l,n,a,p):
    """
creates a sinewave with length 'l', number of cycles 'n'
amplitude 'a' and number of points 'p'
example:
l1=sinewave(l=50,n=3,a=2,p=50)
fo(f'''
color("blue") p_line3d({l1},.3);
''')
    """

    w1=[[i,a*sin(d2r(n*i*360/l))]  for i in linspace(0,l,p)]
    return l_(w1)

def cosinewave(l,n,a,p):
    """
creates a cosinewave with length 'l', number of cycles 'n'
amplitude 'a' and number of points 'p'
l1=cosinewave(l=50,n=3,a=2,p=50)
fo(f'''
color("blue") p_line3d({l1},.3);
''')
    """

    w1=[[i,a*cos(d2r(n*i*360/l))]  for i in linspace(0,l,p)]
    return l_(w1)
    
def mod(a,b):
    """
    function to calculate remainder of a divison of numbers
    example:
    mod(6,2) => 0
    """
    return round(a-sign(a)/sign(b)*b*floor(round(abs(a/b),10)),10);

def e_wave(l=50,a=1,w=0.1,t=100):
    """
create a graph of exponential function a*e^-(wt) where
w: omega
t: time steps
a: amplitude
l: length of time
example:
l1=e_wave(l=20,a=10,w=.2,t=100)
fo(f'''
color("blue") p_line3d({l1},.3);
''')
    """
    return l_([[i,a*exp(-i*w)]  for i in linspace(0,l,t)])

def waves_2d_multiply(w1,w2,a=1):
    """
multiply 2 waves 'w1' and 'w2' and amplitude of the resultant wave is set to 'a'
example:
l1=e_wave(l=20,a=10,w=.2,t=100)
l2=sinewave(20,3,2,100)
l3=waves_2d_multiply(l1,l2,10)
fo(f'''
color("blue") p_line3d({l3},.3);
''')
    """
    w3=[ [w1[i][0],(w1[i][1]*w2[i][1]) ] for i in range(len(w1))]
    a_max=array(w3)[:,1].max()
    w3=[ [p[0],p[1]/a_max*a ] for p in w3]
    return l_(a_(w3))

def waves_2d_add(w1,w2,a=1):
    """
add 2 waves 'w1' and 'w2' and amplitude of the resultant wave is set to 'a'
example:
l1=e_wave(l=20,a=10,w=.2,t=100)
l2=sinewave(20,3,2,100)
l3=waves_2d_add(l1,l2,10)
fo(f'''
color("blue") p_line3d({l3},.3);
''')
    """
    w3=[ [w1[i][0],(w1[i][1]+w2[i][1]) ] for i in range(len(w1))]
    a_max=array(w3)[:,1].max()
    w3=[ [p[0],p[1]/a_max*a ] for p in w3]
    return l_(a_(w3))

def waves_2d_max(w1,w2,a=1):
    """
max of 2 waves 'w1' and 'w2' and amplitude of the resultant wave is set to 'a'
example:
l1=e_wave(l=20,a=10,w=.2,t=100)
l2=sinewave(20,3,2,100)
l3=waves_2d_max(l1,l2,10)
fo(f'''
color("blue") p_line3d({l3},.3);
''')
    """
    w3=[ [w1[i][0],a_([w1[i][1],w2[i][1]]).max() ] for i in range(len(w1))]
    a_max=array(w3)[:,1].max()
    w3=[ [p[0],p[1]/a_max*a ] for p in w3]
    return l_(a_(w3))

def waves_2d_min(w1,w2,a=1):
    """
min of 2 waves 'w1' and 'w2' and amplitude of the resultant wave is set to 'a'
example:
l1=e_wave(l=20,a=10,w=.2,t=100)
l2=sinewave(20,3,2,100)
l3=waves_2d_min(l1,l2,10)
fo(f'''
color("blue") p_line3d({l3},.3);
''')
    """
    w3=[ [w1[i][0],a_([w1[i][1],w2[i][1]]).min() ] for i in range(len(w1))]
    a_max=array(w3)[:,1].max()
    w3=[ [p[0],p[1]/a_max*a ] for p in w3]
    return l_(a_(w3))

def waves_2d_norm(w1,w2,a=1):
    """
norm 2 waves 'w1' and 'w2' and amplitude of the resultant wave is set to 'a'
example:
l1=e_wave(l=20,a=10,w=.2,t=100)
l2=sinewave(20,3,2,100)
l3=waves_2d_norm(l1,l2,10)
fo(f'''
color("blue") p_line3d({l3},.3);
''')
    """
    w3=[ [w1[i][0],norm([w1[i][1],w2[i][1]]) ] for i in range(len(w1))]
    a_max=array(w3)[:,1].max()
    w3=[ [p[0],p[1]/a_max*a ] for p in w3]
    return l_(a_(w3))

def fit_curve2d(curve,pnt):
    """
    function to transform a 2d shape to a curve
    an example can make this clear
    """
    x_1,y_1=array(curve)[:,0],array(curve)[:,1]
    b=array([l_lenv_o(curve[:i])  for i in range(1,len(curve)+1)])
    x_2=norm([pnt[0],pnt[1]])
    n2=arange(len(b))[b>=x_2][0]
    n1=n2-1
    d_y=b[n2]-b[n1]
    d1=x_2-b[n1]
    d_x=x_1[n2]-x_1[n1]
    x_2=x_1[n1]+d_x/d_y*d1
    d_y=y_1[n2]-y_1[n1]
    d_x=x_1[n2]-x_1[n1]
    d_x1=x_2-x_1[n1]
    e_p=[x_2/norm([pnt[0],pnt[1]])*pnt[0],x_2/norm([pnt[0],pnt[1]])*pnt[1],y_1[n1]+d_y/d_x*d_x1]
    return e_p

def extrude_wave2path(w1,c1):
    """
function to extrude a wave 'w1' to any defined path 'c1'
example:
l1=c23(sinewave(50,5,4,100))
p1=c23(arc(5,0,360,s=99))
l2=extrude_wave2path(l1,p1)
fo(f'''
// wave
color("magenta",.2) p_line3d({l1},.2);
//path
color("cyan",.2) p_line3d({p1},.2);
// final extruded sketch
color("blue") p_line3dc({l2},.3);
''')
    """
    w2=[[0,p[1],p[2]] for p in w1]
    c2=array(c1)
    c3=array(c1[1:]+c1[:1])
    v1=c3-c2
    w3=[]
    for i in range(len(w2)):
        a1=ang(v1[i][0],v1[i][1])
        a2=ang(norm(v1[i][:2]),v1[i][2])
        p3=rot(f'y{-a2}z{a1}',w2[i])
        p3=c2[i]+p3
        w3.append(p3.tolist())
    w3=w3[:-1]
    return w3

def x_fit(curve,pnt):
    """
    fit a point's 'z' coordinate as per a curve in the x-z plane
    """
    a=array(curve)[:,0]
    b=array(curve)[:,2]
    c=pnt[0]
    
    n2=arange(len(a))[a>c][0]
    n1=n2-1
    dy=b[n2]-b[n1]
    dx=a[n2]-a[n1]
    dybdx=dy/dx
    dx1=c-a[n1]
    d=b[n1]+dybdx*dx1
    pnt1=[c,pnt[1],d]
    return pnt1

def convert_3lines2surface(l1,l2,l3,s=50):
    """
convert 3 lines to surface
's' is the number of segments on each surface line
example:
l1=translate([10,0,0],rot('x90y-90',sinewave(100,3,2,100)))
l2,l3=[rot(f'z{i}',point_vector([10,0,0],[0,0,100])) for i in [-90,90]]
s1=convert_3lines2surface(l2,l1,l3)
fo(f'''
color("blue") for(p={[l1,l2,l3]}) p_line3d(p,.3);
color("magenta",.2) for(p={s1}) p_line3d(p,.3);
{swp_surf(s1)}
''')
    """
    l1=path2path1(l2,l1)
    l3=path2path1(l2,l3)
    l1,l2,l3=array([l1,l2,l3])
    surf_1=[arc_3p_3d([l1[i]+[0,0,.00001],l2[i],l3[i]],s)  for i in range(len(l1))]
    return surf_1

def convert_lines2surface_spline(lines,s=50):
    """
    create surface with lines, method used is bspline
    """
    lines=[lines[0]]+[path2path1(lines[0],lines[i])  for i in range(1,len(lines))]
    surf_1=[bspline(p,2,s) for p in cpo(lines)]
    return surf_1



def SurfaceFrom3LinesInDifferentPlanes(l1,l2,l3,s1=20,s2=20):
    """
create surface with 3 lines in different plane.
example:
l1=rot('x90',arc_2p([0,0],[50,0],40))
l2=translate([0,30,0],rot('z20',l1))
l3=rot('x45',sinewave(50,3,2,50))
l3=fit_pline2line(l3,[l1[0],l2[0]])
s1=SurfaceFrom3LinesInDifferentPlanes(l1,l2,l3,50,50)
fo(f'''
color("blue") for(p={[l1,l2,l3]}) p_line3d(p,.3);
color("magenta",.2) for(p={s1}) p_line3d(p,.3);
{swp_surf(s1)}
''')
    
    """
    l1=equidistant_path(l1,s1)
    l2=equidistant_path(l2,s1)
    l3=equidistant_path(l3,s2)
    
    
    l4=fit_pline2line(l3,[l1[-1],l2[-1]])
    surf_1=surface_from_4_lines(l1,l2,l3,l4)
    return surf_1

def mid_point(w1):
    """
finds the mid point of a line or polyline
example:
l1=[[0,0],[10,0]]
p0=mid_point(l1)
l2=translate_2d([0,5],square(10))
p1=mid_point(l2)
fo(f'''
color("blue") p_line3d({l1},.3);
color("magenta") points({[p0]},.5);
color("cyan") p_line3d({l2},.3);
color("brown") points({[p1]},.5);
''')
    """
    return equidistant_path(w1,2)[1]



def cw_3p_3d(points):
    """
    finds the orientation of 3 points in 3d.
    1 means clockwise
    -1 means counter clockwise
    """
    n1=array(nv(points))+[.000001,.000001,0]
    a1=cross(n1,[0,0,-1])
    t1=r2d(arccos(n1@[0,0,-1]))
    sec1=translate(-array(points).mean(0),points)
    sec2=c3t2(axis_rot(a1,sec1,t1))
    l1=len(sec2)
    p0,p1,p2=[sec2[0],sec2[int(l1/3)],sec2[int(l1*2/3)]]
    
    return cw([p0,p1,p2])

def smoothen_3d(p0,s=50):
    """
    draw smooth curves with random points 'p0'
    """
    r0=[r_3p_3d([p0[i],p0[i+1],p0[i+2]])  for i in range(len(p0)-2)]
    n0=[nv([p0[i],p0[i+1],p0[i+2]])  for i in range(len(p0)-2)]
    c0=[cw_3p_3d([p0[i],p0[i+1],p0[i+2]])  for i in range(len(p0)-2)]
    
    arc_1=arc_2p_3d(n0[0],p0[0],p0[1],r0[0],c0[0],50)
    arc_2=arc_2p_3d(n0[-1],p0[-2],p0[-1],r0[-1],c0[-1],50)
    arc_m=[
        [arc_2p_3d(n0[i],p0[i+1],p0[i+2],r0[i],c0[i],50),
        arc_2p_3d(n0[i+1],p0[i+1],p0[i+2],r0[i+1],c0[i+1],50)]
        for i in range(len(p0)-3)]
    pnts=[equidistant_path(p[0],10)[:4]+equidistant_path(p[1],10)[6:]  for p in arc_m]
    arc_n=concatenate([bezier(p,50)  for p in pnts]).tolist()
    return equidistant_path(arc_1+arc_n+arc_2,s)

def smoothen_2d(p0,s=50):
    """
    draw smooth curves with random points 'p0'
    """
    r0=[r_3p([p0[i],p0[i+1],p0[i+2]])  for i in range(len(p0)-2)]
    c0=[cw([p0[i],p0[i+1],p0[i+2]])  for i in range(len(p0)-2)]
    
    arc_1=arc_2p(p0[0],p0[1],r0[0],c0[0],50)
    arc_2=arc_2p(p0[-2],p0[-1],r0[-1],c0[-1],50)
    arc_m=[
        [arc_2p(p0[i+1],p0[i+2],r0[i],c0[i],50),
        arc_2p(p0[i+1],p0[i+2],r0[i+1],c0[i+1],50)]
        for i in range(len(p0)-3)]
    pnts=[equidistant_path(p[0],10)[:4]+equidistant_path(p[1],10)[6:]  for p in arc_m]
    arc_n=concatenate([bezier(p,50)  for p in pnts]).tolist()
    return equidistant_path(arc_1+arc_n+arc_2,s)

def curve_2d(p0,s=50):
    """
    draw smooth curves with random points 'p0'
    """
    r0=[r_3p([p0[i],p0[i+1],p0[i+2]])  for i in range(len(p0)-2)]
    c0=[cw([p0[i],p0[i+1],p0[i+2]])  for i in range(len(p0)-2)]
    
    arc_m=[
        [arc_2p(p0[i+1],p0[i+2],r0[i],c0[i],50),
        arc_2p(p0[i+1],p0[i+2],r0[i+1],c0[i+1],50)]
        for i in range(len(p0)-3)]
    pnts=[[mid_point([mid_point(p[0]),mid_point(p[1])]),
           mid_point([mid_point(p[0]),mid_point(p[1])]),
           p[0][-1]]  for p in arc_m]
    pnts=[p0[0],p0[1]]+concatenate(pnts).tolist()+[p0[-1]]
    pnts=array(pnts).reshape(-1,3,2).tolist()
    arc_n=[ arc_3p(p[0],p[1],p[2]) for p in pnts]
    return concatenate(arc_n).tolist()

def faces_surface(l,m):
    """
    calculate the faces for an open surface
    """
    return array([[
            [
            [m*i+j,m*(i+1)+j,m*i+(j+1)],
            [m*i+(j+1),m*(i+1)+j,m*(i+1)+(j+1)]
            ] 
        for j in range(m-1)
        ] 
           for i in range(l-1)
      ]).reshape(-1,3).tolist()



def surface_offset(surf,d=1):
    """
calculates an offset surface by distance 'd'
example:
l1=rot('x90',sinewave(20,1,2,40))
l2=rot('x90z90',sinewave(20,2,2,40))
s1=surface_from_2_waves(l1,l2,2)
s2=surface_offset(s1,1)
fo(f'''
{swp_surf(s1)}
color("cyan"){swp_surf(s2)}
''')
    """
    f=faces_surface(len(surf),len(surf[0]))
    v=a_(surf).reshape(-1,3)
    tri=v[f]
    p0,p1,p2=tri[:,0],tri[:,1],tri[:,2]
    v1,v2=p1-p0,p2-p0
    n1=cross(v1,v2)
    n1=a_([n1,n1,n1]).transpose(1,0,2).reshape(-1,3)
    g=concatenate(f)
    n1.shape,p0.shape,g.shape
    n2=a_([n1[arange(len(g))[(g==i)]].mean(0) for i in range(len(v))])+.00001
    n2=n2/norm(n2,axis=1).reshape(-1,1)
    v3=v+n2*d
    v3=l_(v3.reshape(len(surf),len(surf[0]),3))
    return v3

def swp_surf(surf_1):
    l,m=len(surf_1),len(surf_1[0])
    f_1=faces_surface(l,m)
    v_1=array(surf_1).reshape(-1,3).tolist()
    return f'polyhedron({v_1},{f_1},convexity=10);'

def surface_thicken(surf_1,d=1):
    """
thicken the surface by amount 'd'
example:
sec2=corner_radius(pts1([[-25,0],[10,5,5],[10,-3,10],[10,5,5],[10,-8,7],[10,1]]),10)  
path2=cytz(corner_radius(pts1([[-35,5,0],[10,8,20],[20,-5,10],[20,8,20],[10,-9,20],[10,1,0]]),10))
surf2=path_extrude_open(sec2,path2)
surf3=surface_thicken(surf2,-1)
fo(f'''
{swp(surf3)}
//color("blue")for(p={surf3})p_line3dc(p,.1,rec=1);
''') 
    """
    surf_2=surface_offset(surf_1,d)
    sol=[surf_1[i]+flip(surf_2[i])  for i in range(len(surf_1))]
    return sol

def surface_thicken_1(surf,d=1):
    """
function same as surface_thicken but with different orientation
example:
sec1=circle(10)
path=corner_radius(pts1([[2,0],[-2,0,2],[0,10,3],[-3,0]]),5)
sol=prism(sec1,path)
sol2=surface_thicken_1(sol,1)
fo(f'''
{swp_c(sol2)} 
''') 
    """
    a=surface_offset(surf,d)
    b=swp_prism_h(surf,a) if d<0 else swp_prism_h(a,surf)
    return b

def boundary_edges_sol(sol):
    """
    finds the boundary edges of the solid
    """
    l,m=len(sol),len(sol[0])
    f_1=faces(l,m)[1:-1]
    f_2=array([ seg(p)  for p in f_1]).reshape(-1,2)
    f_3=[]
    for p in f_2:
     if (f_2==flip(p)).all(1).any()==0:
         f_3.append(p)
    
    f_3=array(f_3)
    f_4=[f_3[0]]
    f_3=array(exclude_points(f_3,f_4))
    while ((f_3[:,0]==f_4[-1][1]).any()):
        f_4.append(f_3[f_3[:,0]==f_4[-1][1]][0])
        f_3=array(exclude_points(f_3,[f_4[-1]]))
    
    if f_3.tolist()==[]:
        a=f_4.tolist()
        return flip(array(a)[0][:,1]).tolist()
    else:
        f_5=[f_3[0]]
        f_3=array(exclude_points(f_3,f_5))
        while (f_3.tolist()!=[]):
            f_5.append(f_3[f_3[:,0]==f_5[-1][1]][0])
            f_3=array(exclude_points(f_3,[f_5[-1]]))
        a=[f_4,f_5]
        return array([flip(array(a)[0][:,1]),flip(array(a)[1][:,0])]).tolist()

def boundary_edges_surf(surf):
    """
    finds the boundary edges of a surface
    """
    l,m=len(surf),len(surf[0])
    f_1=faces_surface(l,m)
    f_2=array([ seg(p)  for p in f_1]).reshape(-1,2)
    f_3=[]
    for p in f_2:
     if (f_2==flip(p)).all(1).any()==0:
         f_3.append(p)
    
    f_3=array(f_3)
    f_4=[f_3[0]]
    f_3=array(exclude_points(f_3,f_4))
    while (f_3.tolist()!=[]):
        f_4.append(f_3[f_3[:,0]==f_4[-1][1]][0])
        f_3=array(exclude_points(f_3,[f_4[-1]]))
    
    a=array(f_4).tolist()
    return flip(array(a)[:,1]).tolist()

def extend_arc2d(a,theta=0,s=20,both=0):
    """
extend a 2d arc by theta degrees
example:
l1=arc_2p([0,0],[10,0],5)
l2=extend_arc2d(l1,90,s=50)
fo(f'''
color("blue",.2) p_line3d({l1},.3);
color("magenta") p_line3d({l2},.28);
''')
    """
    p1,p2=array([a[0],a[-1]])
    r=r_arc(a)
    cp=cp_arc(a)
    v1,v2=p1-cp,p2-cp
    a1,a2=ang(v1[0],v1[1]),ang(v2[0],v2[1])
    a3= (a2+360 if a2<a1 else a2) if cw(a)==-1 else (a2 if a2<a1 else a2-360)
    if both==0:
        return arc(r,a1,a3+theta,cp,s) if cw(a)==-1 else arc(r,a1,a3-theta,cp,s)
    elif both==1:
        return arc(r,a1-theta,a3+theta,cp,s) if cw(a)==-1 else arc(r,a1+theta,a3-theta,cp,s)

def extend_arc3d(a,theta=0,s=20,both=0):
    """
extend a 3d arc by theta degrees
example:
l1=arc_2p_3d([0,1,-.3],[0,0,0],[10,0,0],5)
l2=extend_arc3d(l1,90,s=50)
fo(f'''
color("blue",.2) p_line3d({l1},.3);
color("magenta") p_line3d({l2},.28);
''')
    """
    n1=array(nv(a))
    p0,p1=array([a[0],a[-1]])
    r=r_3p_3d(a)
    a1=cross(n1,[0,0,-1])+[.000001,.0000001,0]
    t1=r2d(arccos(n1@[0,0,-1]))
    sec1=translate(-array([p0,p1]).mean(0),[p0,p1])
    sec2=c3t2(axis_rot(a1,sec1,t1))
    sec3=c3t2(axis_rot(a1,a,t1))
    
    
    pa,pb=sec2
    arc1=arc_2p(pa,pb,r,cw(sec3),s=s)
    arc1=extend_arc2d(arc1,theta,s,both)
    arc1=translate(array([p0,p1]).mean(0),axis_rot(a1,arc1,-t1))
    return arc1

def vnf_surf(surf_1):
    l,m=len(surf_1),len(surf_1[0])
    f_1=faces_surface(l,m)
    v_1=array(surf_1).reshape(-1,3).tolist()
    return [v_1,f_1]

def boundary_edges(faces_list):
    """
    finds the boundary edges of a surface
    """
    f_1=faces_list
    f_2=array([ seg(p)  for p in f_1]).reshape(-1,2)
    f_3=[]
    for p in f_2:
     if (f_2==flip(p)).all(1).any()==0:
         f_3.append(p)
    
    f_3=array(f_3)
    f_4=[f_3[0]]
    f_3=array(exclude_points(f_3,f_4))
    while (f_3.tolist()!=[]):
        f_4.append(f_3[f_3[:,0]==f_4[-1][1]][0])
        f_3=array(exclude_points(f_3,[f_4[-1]]))
    
    a=array(f_4).tolist()
    return flip(array(a)[:,1]).tolist()

def create_mesh(v,f):
    """ function to create stl_mesh from vertices and faces"""
    v_1,f_1=array(v),array(f)
    sol = mesh.Mesh(zeros(f_1.shape[0], dtype=mesh.Mesh.dtype))
    sol.vectors=v_1[f_1]
    return sol

def arc_with_start_pt_and_cp(start_point=[],center_point=[],theta=90,segments=30):
    """
    function to draw an arc with known center_point and start_point
    """
    center_point,start_point=array([center_point,start_point])
    v1=start_point-center_point
    arc_1=center_point+[rot2d(i,v1) for i in linspace(0,theta,segments)]
    arc_1=arc_1.tolist()
    return arc_1

def fillet_line_circle(l1,c1,r2,cw=-1,option=0,s=50):
    """
function to draw a fillet between a line and a circle
option can be '0' or '1' to flip the fillet from one side to another
's' is the number of segments in the arc
example:
h=15
line=[[-10,h],[30,h]]
cir1=circle(10,[10,10])
r2=5
s=20
fillet1=fillet_line_circle(line,cir1,r2,1)
fillet2=fillet_line_circle(line,cir1,r2,2)
fillet3=fillet_line_circle(line,cir1,r2,3)
fillet4=fillet_line_circle(line,cir1,r2,4)
fillet5=fillet_line_circle_internal(line,cir1,1,1)
fillet6=fillet_line_circle_internal(line,cir1,r2,2)
fillet7=fillet_line_circle_internal(line,cir1,1,3)
fillet8=fillet_line_circle_internal(line,cir1,r2,4)

fo(f'''
 
color("blue",.1)p_line({line},.3);
color("violet",.2)p_line({cir1},.3);
color("cyan")p_lineo({fillet1},.3);
color("blue")p_lineo({fillet2},.3);
color("magenta")p_lineo({fillet3},.3);
color("green")p_lineo({fillet4},.3);
color("blue")p_lineo({fillet5},.3);
color("magenta")p_lineo({fillet6},.3);
color("cyan")p_lineo({fillet7},.3);
color("green")p_lineo({fillet8},.3);

''')
    """
    v1=array(l1[1])-array(l1[0])
    u1=v1/norm(v1)
    cp1=cp_arc(c1)
    v2=array(cp1)-array(l1[0])
    u2=v2/norm(v2)
    l_1=norm(cross(c23(v1),c23(v2)))/norm(v1)
    p3=array(l1[0])+u1*(u1@v2)
    r1=r_arc(c1)
    d=sqrt((r1+r2)**2-(l_1+r2)**2) if option==0 else sqrt((r1+r2)**2-(l_1-r2)**2)
    p2=p3-u1*d
    v3=array(cp1)-p3
    u3=v3/norm(v3)
    cp2=p2-u3*r2 if option==0 else p2+u3*r2
    v4=array(cp1)-cp2
    u4=v4/norm(v4)
    p4=cp2+u4*r2
    return arc_2p(p2,p4,r2,cw=cw,s=s)

def fillet_intersection_lines_3d(l1,l2,r,s=10,side=0):
    """
function calculates the fillet at intersection between 2 3d lines in 1 plane
'l1' and 'l2'
r: radius of fillet
s: segments of fillet
example:
l1=point_vector([-12,0,0],[5,0,5])
l2=flip(point_vector([0,-12,0],[0,5,5]))
f1=fillet_intersection_lines_3d(l1,l2,5,20)
fo(f'''
color("blue") p_line3d({l1},.3);
color("cyan") p_line3d({l2},.3);
color("magenta") p_line3d({f1},.2);
''')
    """
    p0=i_p3d(l1,l2)
    # p1=l1[a_([norm([p0,p]) for p in l1]).argsort()[-1]]
    # p2=l2[a_([norm([p0,p]) for p in l2]).argsort()[-1]]
    p1=l1[0]
    p2=l2[-1]
    sec=[p1,p0,p2]
    sec1=rot_sec2xy_plane(sec)
    sec2=c32(sec1)
    l1=[sec2[0],sec2[1]]
    l2=[sec2[1],sec2[2]]
    sec3=fillet_intersection_lines(l1,l2,r,s,side)
    a=a_(sec[0])
    b=a_(c23(c32(sec)[0]))
    v1=uv(best_fit_plane(sec)[0])
    if (l_(a_(v1).round(4))==[0,0,-1]) | (l_(a_(v1).round(4))==[0,0,1]):
        sec5=translate(a-b,sec3)
    else:
        v2=[0,0,-1]
        theta=l_(r2d(arccos(a_(v1)@v2)))
        ax_1=cross(v1,v2)
        sec4=axis_rot_1(c23(sec3),ax_1,c23(sec2)[0],-theta)
        sec5=translate(a-b,sec4)
    return sec5

def fillet_line_circle_internal_3d(line,cir,r=1,o=1,s=20):
    """
creates fillet between a line and an arc in 3d space.
fillet is created internal
The line and arc should lie in same plane
'r' is the radius of the fillet
'o' is option of fillet, means which side of the line and arc it should be created.
the value of 'o' can be set from 1 to 4, as there could be only 4 sides a fillet can be created
's' is the number of segments in the fillet
example:
h=15
line=[[-10,h],[30,h]]
cir1=circle(10,[10,10])
r2=5
s=20
fillet1=fillet_line_circle(line,cir1,r2,1)
fillet2=fillet_line_circle(line,cir1,r2,2)
fillet3=fillet_line_circle(line,cir1,r2,3)
fillet4=fillet_line_circle(line,cir1,r2,4)
fillet5=fillet_line_circle_internal(line,cir1,1,1)
fillet6=fillet_line_circle_internal(line,cir1,r2,2)
fillet7=fillet_line_circle_internal(line,cir1,1,3)
fillet8=fillet_line_circle_internal(line,cir1,r2,4)

fo(f'''
 
color("blue",.1)p_line({line},.3);
color("violet",.2)p_line({cir1},.3);
color("cyan")p_lineo({fillet1},.3);
color("blue")p_lineo({fillet2},.3);
color("magenta")p_lineo({fillet3},.3);
color("green")p_lineo({fillet4},.3);
color("blue")p_lineo({fillet5},.3);
color("magenta")p_lineo({fillet6},.3);
color("cyan")p_lineo({fillet7},.3);
color("green")p_lineo({fillet8},.3);

''')
    """
    a=line
    b=cir
    n1,incpt=best_fit_plane(b)
    n1=uv(n1)
    c=rot_sec2xy_plane(a+b)
    d=c32(c)
    e=d[:2]
    f=d[2:]
    a1=fillet_line_circle_internal(e,f,r,o,s)
    n2=uv(cross(n1,[0,0,-1]))
    theta=r2d(arccos(a_(n1)@a_([0,0,-1])))
    g,h,a2=[axis_rot_1(p,n2,f[0],-theta) for p in [e,f,a1]]
    d1=a_(b[0])-a_(h[0])
    i,j,a3=[translate(d1,p) for p in [g,h,a2]]
    return a3



def fillet_line_circle_3d(line,cir,r=1,o=1,s=20):
    """
creates fillet between a line and an arc in 3d space.
fillet is created external
The line and arc should lie in same plane
'r' is the radius of the fillet
'o' is option of fillet, means which side of the line and arc it should be created.
the value of 'o' can be set from 1 to 4, as there could be only 4 sides a fillet can be created
's' is the number of segments in the fillet
    """
    a=line
    b=cir
    n1,incpt=best_fit_plane(b)
    n1=uv(n1)
    c=rot_sec2xy_plane(a+b)
    d=c32(c)
    e=d[:2]
    f=d[2:]
    a1=fillet_line_circle(e,f,r,o,s)
    n2=uv(cross(n1,[0,0,-1]))
    theta=r2d(arccos(a_(n1)@a_([0,0,-1])))
    g,h,a2=[axis_rot_1(p,n2,f[0],-theta) for p in [e,f,a1]]
    d1=a_(b[0])-a_(h[0])
    i,j,a3=[translate(d1,p) for p in [g,h,a2]]
    return a3

def mirror_line(p1,n1,loc):
    """
function to mirror the points list 'p0' defined by mirroring plane 'n1' with location 'loc'
example:
l1=circle(5,[10,10])
mirroring_plane=[1,1,0]
passing_through_point=[0,0,0]
pl1=plane(mirroring_plane,size=[50,50],intercept=passing_through_point)
l2=mirror_line(l1,mirroring_plane,passing_through_point)
fo(f'''
//original line
color("blue") p_line3d({l1},.3);
// mirrored along defined plane
color("magenta") p_line3d({l2},.28);
%{swp_surf(pl1)}
''')
    """
    p0=c23(p1)
    a=ppplane(p0,n1,loc)
    b=[]
    for i in range(len(a)):
        v1=array(a[i])-p0[i]
        b.append((array(a[i])+v1).tolist())
    return c32(b) if a_(p1).shape[-1]==2 else b


def solid_from_fillet(fillet_1,d):
    """
    creates a cutting edge from fillet
    """
    fillet_1=cpo(fillet_1)[:-1]
    fillet_2=surface_offset(fillet_1,d)
    sol=fillet_1+flip(fillet_2)
    return cpo(sol)

def solid_from_fillet_closed(fillet_1,d):
    """
    creates a cutting edge from fillet with closed section lines
    """
    fillet_1=cpo(fillet_1)[:-1]
    fillet_2=surface_offset(fillet_1,d)
    sol=fillet_1+flip(fillet_2)
    sol=cpo(sol)
    sol=sol+[sol[0]]
    return sol


def surface_normal(s1,length=1):
    """
    calculates the normal of a surface (average)
    """
    a=faces_1(len(s1),len(s1[0]))
    b=array(s1).reshape(-1,3)
    c=b[a]
    pa,pb,pc=c[:,0],c[:,1],c[:,2]
    v1,v2=pb-pa,pc-pa
    v3=cross(v1,v2)
    u3=v3/norm(v3,axis=1).reshape(-1,1)
    u3=u3.mean(0)
    u3=u3/norm(u3)*length
    return u3.tolist()
  


# def path_offset_3d(sec,d):
#     """
#     path_offsets an enclosed section in 3d space, in case the section is in 1 plane
#     sec: section in 3d space
#     d: offset distance -ve sign means inner offset and +ve sign is outer offset
#     refer to the file"example of various functions" for application examples
    
#     """
#     sec0=remove_extra_points(sec)
#     sec0=rot('z.00001',sec0)
#     avg1=array(sec0).mean(0)
#     sec1=translate(-avg1,sec0)

#     # nv1=-array(nv(sec1))
#     nv1=-a_(best_fit_plane(sec1)[0])
#     nz=[0,0,1]
#     nr=cross(nv1,nz) if abs(nv1).tolist()!=[0,0,1] else nv1
#     theta=r2d(arccos(nv1@array(nz)))
#     sec1=axis_rot(nr,sec1,theta)
#     z_values=array(sec1)[:,2]-avg1[2]
#     sec1=ppplane(sec1,[0,0,1],[0,0,0])
#     sec1=c3t2(sec1)
#     x_values=array([l_len([[0,0],p])  for p in sec1])
#     sec2=path_offset_n(sec1,d)
#     x1_values=array([l_len([[0,0],p])  for p in sec2])
#     z1_values=z_values/x_values*x1_values
#     z1_values=array([[0,0,p] for p in z1_values])
#     sec2=array(c2t3(sec2))
#     sec2=axis_rot(nr,sec2,-theta)
#     sec2=translate(array(sec).mean(0),sec2)
#     return sort_points(sec,sec2)

def surface_line_vector(line=[[0,0,0],[10,0,0]],vector=[0,0,1],both_sides=0):
    """
draw a surface base on a line and a vector.
if surface is required both sides of the line, both_sides option should be marked as '1' else default for only one side is '0'
example:
l1=rot('x90',sinewave(50,3,2,50))
s1=surface_line_vector(line=l1,vector=[0,20,0],both_sides=1)
fo(f'''
color("blue") p_line3d({l1},.3);
%{swp_surf(s1)}
''')
    """
    l_1=translate(array(vector),line)
    l_2=translate(-array(vector),line) if both_sides==1 else line
    return [l_2,l_1]


def mirror_surface(surf_1,n1,loc=[0,0,0]):
    """
function to mirror a solid or surface base on a mirroring plane given by vector 'n1'
passing through a point 'loc'
example:
l1=rot('x90',sinewave(50,3,2,50))
s1=surface_line_vector(line=l1,vector=[0,20,0],both_sides=1)
s2=mirror_surface(s1,[1,1,0],[0,-30,0])
fo(f'''
color("blue") p_line3d({l1},.3);
// original surface
color("blue",.3){swp_surf(s1)}
//mirrored surface
color("cyan",.3){swp_surf(s2)}
// mirroring plane
%{swp_surf(plane([1,1,0],[100,100],[0,-30,0]))}
''')
    """
    surf_1_1=[mirror_line(surf_1[i],n1,loc) for i in range(len(surf_1))]
    return surf_1_1

def solid_from_2surfaces(surf_1,surf_2):
    """
function to make a solid from 2 surfaces 
example:
l1=rot('x90',sinewave(20,1,2,40))
l2=rot('x90z90',sinewave(20,2,2,40))
s1=surface_from_2_waves(l1,l2,1)
s2=surface_offset(s1,1)
sol=solid_from_2surfaces(s1,s2)
fo(f'''
{swp(sol)}
''')
    """
    sol=[surf_1[i]+flip(surf_2[i]) for i in range(len(surf_1))]
    # sol=cpo(sol)
    return sol

def sec2surface(surf_1,s=1):
    """
function to convert a closed polygon to lines
example:
a=c2t3(circle(10))
b=sec2surface(a)
fo(f'''
color("blue")p_line3d({a},.2,1);
color("magenta")for(p={b})p_line3d(p,.2,1);
{swp_surf(b)}''')
    """
    return [equidistant_path([surf_1[i],surf_1[-i-1]] ,s)
            for i in range(int(len(surf_1)/2))]



def rot_sec2xy_plane(sec):
    """
function to rotate any section open or closed parallel to x-y plane
example:
a=cylinder(r=10,h=50,s=100)
s1=surface_line_vector([[-20,-20,5],[-20,20,5]],[40,0,30])
l1=ip_sol2sol(s1,a)
l2=rot_sec2xy_plane(l1)
fo(f'''
%{swp(a)}
%{swp_surf(s1)}
//original line
color("blue") p_line3dc({l1},.3);
// line rotated parallel to x-y plane
color("magenta") p_line3dc({l2},.3);
''')
    """
    # n1=nv(sec)
    n1=uv(best_fit_plane(sec)[0])
    if (array(n1).round(5).tolist()==[0,0,1]) | (array(n1).round(5).tolist()==[0,0,-1]) :
        return sec
    else:
        v1=cross(n1,[0,0,-1])
        t1=r2d(arccos(array(n1)@[0,0,-1]))

        l_3=axis_rot_1([sec],v1,sec[0],t1)[0]

        return l_3
    
def surround_3d(path,r,s=20):
    """
    function to surround a path to create a rounded section on a 3d path
    
    """
    n1=nv(path)
    if (array(n1).round(5).tolist()==[0,0,1]) | (array(n1).round(5).tolist()==[0,0,-1]) :
        path1= c3t2(path)
    else:
        v1=cross(n1,[0,0,-1])
        t1=r2d(arccos(array(n1)@[0,0,-1]))
        l_1=[[0,0,0],(array(n1)*5).tolist()]
        l_2=[[0,0,0],(array(v1)*5).tolist()]
        l_3=axis_rot_1([path],v1,path[0],t1)[0]
        path1=c3t2(l_3)
    path1=remove_extra_points(array(path1).round(5))
    a=path_offset_n(path1,r)
    b=path_offset_n(path1,-r)
    b=flip(b)
    arc1=arc_2p(a[-1],b[0],r,-1,s)
    arc2=arc_2p(b[-1],a[0],r,-1,s)
    sec=a[1:-1]+arc1+b[1:-1]+arc2
    sec=translate([0,0,path[0][2]],sec)
    if (array(n1).round(5).tolist()==[0,0,1]) | (array(n1).round(5).tolist()==[0,0,-1]) :
        sec=sec
    else:
        sec=axis_rot_1([sec],v1,path[0],-t1)[0]
    return equidistant_pathc(sec,s)

def plane(nv,size=[100,100],intercept=[0,0,0]):
    """
plane defined by normal vector 'nv' with size as defined
example:
n1=[1,1,0] # normal vector to define plane
pl1=plane(nv=n1,size=[100,100],intercept=[0,0,0])
l1=vector2line(n1,30)
fo(f'''
%{swp_surf(pl1)}
color("blue") p_line3d({l1},.5);
''')
    """
    d1,d2=(size/2,size/2) if array(size).shape==() else (size[0]/2,size[1]/2)
    n1=array(nv)/norm(nv)
    if array(n1).round(5).tolist()==[0,0,1]:
        v1=[1,0,0]
    else:
        v1=c2t3(rot2d(90,c3t2(n1)))
    v1=array(v1)/norm(v1)
    v2=axis_rot(v1,[n1],90)[0]
    v1=array([array(v1)*-d1,array(v1)*d1])
    v2=array(v2)*d2
    s1=surface_line_vector(v1,v2,1)
    s1=translate(intercept,s1)
    return s1

def align_sec_2(sec1):
    """
    function to align points of a section to obtain the non twisted optimised surface after using
    function sec2surface(sec1)
    """
    
    area1=[array([l_lenv_o(p) for p in sec2surface(sec1[i:]+sec1[:i])]).sum() for i in range(len(sec1)) ]
    i=array(area1).argmin()
    sol2=sec1[i:]+sec1[:i]
    return sol2

def s_int1_list(sec1):
    """
calulates the self intersection list numbers of segment lists 'sec1'
it picks the intersection points only if the 2 lines are crossing each other
example:
a=circle(10,s=7)
b=point_vector([-10,-10],[20,5])
segments=seg(b)[:-1]+seg(a)
n=s_int1_list(segments)
p0=s_int1(segments)
tx=[label_linear(segments[i],f's-{i}') for i in range(len(segments))]
fo(f'''
color("blue") for(p={[a,b]}) p_line3dc(p,.3);
color("magenta") points({p0},1);
{tx[0]}{tx[1]}{tx[2]}{tx[3]}{tx[4]}{tx[5]}{tx[6]}
''')
n # list of paired segment numbers where intersection is observed
    """
    n=len(sec1)
    a=array(sec1)[comb_list(n)]
    p0=a[:,0][:,0]
    p1=a[:,0][:,1]
    p2=a[:,1][:,0]
    p3=a[:,1][:,1]
    v1=a_(rot2d(0.00001,p1-p0))
    v2=p3-p2
    iim=array([v1,-v2]).transpose(1,0,2).transpose(0,2,1)
    im=inv(iim)
    p=p2-p0

    t=einsum('ijk,ik->ij',im,p)
    dcn=(t[:,0].round(4)>0)&(t[:,0].round(4)<1)&(t[:,1].round(4)>0)&(t[:,1].round(4)<1)

    return comb_list(n)[dcn]



def points_inside_offset_surround(sec,sec2,r):
    """
    finds all the points, in a list of points 'sec2' which are inside the offset surround of an enclosed section 'sec'
    """
    sec2=array(sec2)
    clean=cs1(sec,abs(r))
    clean1=[p[1:]+[p[0]] for p in clean]
    m,n,_=array(clean).shape
    o,_=sec2.shape
    v1=array([[[1,0]]*n]*m)
    v2=array(clean1)-array(clean)
    iim=array([v1,-v2]).transpose(1,2,0,3).transpose(0,1,3,2)+[0,.00001]
    im=array([pinv(iim)]*o)
    p=(array(clean)[:,:,None]-sec2).transpose(2,0,1,3)
    t=einsum('ijklm,ijkm->ijkl',im,p)
    decision1=((t[:,:,:,0]>=0)&(t[:,:,:,1]>=0)&(t[:,:,:,1]<=1))
    sec3=sec2[(decision1.sum(2)==1).any(1)].tolist()
    return sec3



def list_r1(sec):
    """
    function to list the radius at each point of a section
    """
    x_6=list_r(sec)
    x_6[0]=x_6[1] if x_6[1]>0 else x_6[0]
    x_6[-1]=x_6[-2] if x_6[-2]>0 else x_6[-1]
    for i in range(1,len(x_6)-1):
        x_6[i]= x_6[i+1] if x_6[i+1]>0 else x_6[i-1] if x_6[i-1]>0 else x_6[i]
    return x_6

def exclude_numbers(a,b):
    """
    exclude list of numbers 'b' from 'a'
    example:
    a=[1,2,3,4,5]
    b=[4,5,6,7]
    exclude_numbers(a,b) => array([1, 2, 3])
    """
    a,b=array(a),array(b)
    return a[~(a[:,None]==b[None,:]).any(1)]

def subset(b,a):
    """
returns True if 'b' is subset of 'a'
example:
a=arange(10)
b=arange(10)[2:-2]
subset(b,a)
    """
    return l_((array(b)[:,None]==array(a)[None,:]).any(1).all())

def points_inside_offset_surround_list(sec,sec2,r):
    """
finds all the points list, in a list of points 'sec2' which are inside the offset surround of an enclosed section 'sec'
example:
sec=corner_radius(pts1([[0,0,1],[7,5,2],[5,7,3],[-5,7,2],[-7,5,3]]),10)
a=(random.random(5000)*(30-0)+(-10))
b=(random.random(5000)*(30-0)+(-5))
c=array([a,b]).transpose(1,0).tolist()
r=3
sec2=cs1(sec,r-.01)
p_0=points_inside_offset_surround(sec,c,r-.01)
fo(f'''
    color("blue")p_line3dc({sec},.1,1);
    color("cyan")points({c},.2);
    color("green",.1)for(p={sec2})p_line3dc(p,.1,1);
    color("magenta")points({p_0},.2);
    
    ''')
    """
    sec2=array(sec2)
    clean=cs1(sec,abs(r))
    clean1=[p[1:]+[p[0]] for p in clean]
    m,n,_=array(clean).shape
    o,_=sec2.shape
    v1=array([[[1,0]]*n]*m)
    v2=array(clean1)-array(clean)
    iim=array([v1,-v2]).transpose(1,2,0,3).transpose(0,1,3,2)+[0,.00001]
    im=array([pinv(iim)]*o)
    p=(array(clean)[:,:,None]-sec2).transpose(2,0,1,3)
    t=einsum('ijklm,ijkm->ijkl',im,p)
    decision1=((t[:,:,:,0]>=0)&(t[:,:,:,1]>=0)&(t[:,:,:,1]<=1))
    return arange(len(sec2))[(decision1.sum(2)==1).any(1)]

def offset_1(sec,r):
    """
    function returns offset of a enclosed section 'sec' by a distance 'r'
    """
    if r==0:
        return sec
    else:
        sec1=sec
        sec=sec1 if cw(sec1)==-1 else flip(sec1)
        s1=offset_segv(sec,r)
        s2=intersections(s1)
        s3=s_int1(seg(s2))
        if s3!=[]:
            n1=s_int1_list(seg(s2))
            s4=seg(s2)
            for i in range(len(n1)):
                s4[n1[i][0]]=[s4[n1[i][0]][0],s3[i] ,s4[n1[i][0]][1]]
                s4[n1[i][1]]=[s4[n1[i][1]][0],s3[i] ,s4[n1[i][1]][1]]

            s5=remove_extra_points(concatenate(s4))
            s6=points_inside_offset_surround(sec,s5,abs(r)-.01)
            s7=points_inside_offset_surround_list(sec,s5,abs(r)-.01)
            s8=arange(len(s5))[(array(s5)[:,None]==array(s3)[None,:]).any(1).all(1)]
            a=~(arange(len(s5))[:,None]==s7[None,:]).any(1)
            s9=[array(s5)[a][0]] if a[0]==0 else [array(s5)[0]]
            for i in range(1,len(a)):
                if a[i]==0:
                    s9.append(s9[-1])
                else:
                    s9.append(array(s5)[i])
            s10=array(s9)[~(arange(len(s9))[:,None]==s8[None,:]).any(1)].tolist()   
        else:
            s10=s2
        return s10 if cw(sec1)==-1 else flip(s10)
    
def vertices(sol_1):
    """
    returns the vertices of a solid 'sol_1'
    """
    return array(sol_1).reshape(-1,3)

def coil(r1,r2,n1=1):
    """
function to draw a coil with initial radius 'r1', final_radius 'r2' and numbe of coils 'n1'
example:
coil=array([i/360*array([cos(d2r(i)),sin(d2r(i)),-1]) for i in linspace(0,3600,720)]).tolist()
cyl1=translate([0,0,-10],cylinder(r1=10,r2=0.1,h=10))
fo(f'''
translate([0,0,10]){{  
color("blue")p_line3d({coil},.05);
%{swp(cyl1)}
}}
    ''')
    """
    r_l=[i for i in linspace(r1,r2,360*n1)]
    theta_1=[i for i in linspace(0,360*n1,360*n1)]
    c1=[[r_l[i]*cos(d2r(theta_1[i])),r_l[i]*sin(d2r(theta_1[i]))] for i in range(len(r_l))]
    return c1

def corner_n_radius_list(p0,r_l,n=10):
    """
corner list 'p0' and radius list 'r_l' will create a smothened section
'n' is the number of segments in each filleted corner
example:
sec=[[0,0],[10,0],[10,10],[0,10]]
r_l=[1,2,3,4]
sec=corner_n_radius_list(sec,r_l,10)
txt1=dim_radial(sec[2:8],.2,'magenta',.5)
txt2=dim_radial(sec[15:20],.2,'magenta',.5)
txt3=dim_radial(sec[25:30],.2,'magenta',.5)
txt4=dim_radial(sec[35:40],.2,'magenta',.5)

fo(f'''
color("blue") p_line3dc({sec},.3);
{txt1}{txt2}{txt3}{txt4}
''')
    """
    # r_l=[.01 if i==0 else i for i in r_l]
    p1=seg(p0)
    p2=[p1[-1]]+p1[:-1]
    p3=p1
    s1=[fillet_intersection_lines(p2[i],p3[i],r_l[i],s=n) if r_l[i]!=0 else [p0[i]] for i in range(len(p1))]
    s2=concatenate(s1).tolist()
    return remove_extra_points(array(s2).round(5))

def corner_n_radius_list_3d(p0,r_l,n=10):
    """
corner list 'p0' and radius list 'r_l' will create a smothened 3d path or closed section
'n' is the number of segments in each filleted corner
example:
sec=[[0,0,0],[10,0,0],[10,10,3],[0,10,0]]
r_l=[1,2,3,4]
sec=corner_n_radius_list_3d(sec,r_l,10)
txt1=dim_radial(sec[2:8],.2,'magenta',.5)
txt2=dim_radial(sec[12:20],.2,'magenta',.5)
txt3=dim_radial(sec[22:30],.2,'magenta',.5)
txt4=dim_radial(sec[32:40],.2,'magenta',.5)

fo(f'''
color("blue") p_line3dc({sec},.3);
{txt1}{txt2}{txt3}{txt4}
''')
    """

    p1=[p0[-1]]+p0[:-1]
    s1=[]
    for i in range(len(p0)):
        if i<len(p0)-1 and r_l[i]>0:
            a=fillet_3p_3d(p1[i],p0[i],p0[i+1],r_l[i],n)[1:]
        elif i<len(p0)-1 and r_l[i]==0:
            a=[p0[i]]
        elif i==len(p0)-1 and r_l[i]==0:
            a=[p0[i]]
        elif i==len(p0)-1 and r_l[i]>0:
            a=fillet_3p_3d(p1[i],p0[i],p0[0],r_l[i],n)[1:]
        s1.append(a)
    
    s1=remove_extra_points(concatenate(s1).round(5))
    return s1

def path_extrude_closed(sec_1,path,twist=0):
    """
function to extrude a closed section to a closed path
twist=0 for simple path extrudes
set twist=1 for complex path extrudes
example:
sec=corner_radius(pts1([[0,0,.2],[3,0,.2],[0,2,1.49],[-3,0,1.49]]),10)
path=c2t3(circle(10))
sol=path_extrude_closed(sec,path)
fo(f'''
{swp_c(sol)}
color("blue")p_line3d({path},.1);
    ''')
    """
    if twist==0:
        p1=path
        p2=path[1:]+[path[0]]
        p1,p2=array([p1,p2])
        v1=p2-p1
        u1=v1/norm(v1,axis=1).reshape(-1,1)
        v2=concatenate([[(u1[-1]+u1[0])/2], (u1[1:]+u1[:-1])/2])
        sec2=[]
        for i in range(len(path)):
            sec1=translate(path[i],sec2vector(v2[i],sec_1))
            sec2.append(sec1)
        sec2=sec2+[sec2[0]]
        return sec2
    elif twist==1:
        sec=[[1,0],[0,1],[-1,0],[0,-1]]
        t_l=tangents_along_path(path)
        v1=array([array(p[1])-array(p[0]) for p in t_l]).tolist()
        sol_1=[translate(path[i],sec2vector(v1[i],sec))  for i in range(len(path))]
        sol_2=align_sol(sol_1,360/len(t_l)/2)
        cp_1=array(sec_1).mean(0)
        sec_2=translate_2d(-cp_1,sec_1)
        sol_3=[translate(path[i],sec2vector(v1[i],sec_2))  for i in range(len(path))]
        sol_4=align_sol(sol_3,360/len(sol_3)/2)
        sol_5=[]
        for i in range(len(path)):
            a=mid_point([sol_1[i][0],sol_1[i][2]])
            v2=array(sol_2[i][0])-array(a)
            v3=array(sol_2[i][1])-array(a)
            sol_5.append(translate(-v3*cp_1[0],translate(-v2*cp_1[1],sol_4[i])))
        return sol_5+[sol_5[0]]

def path_extrude_open(sec_1,path,twist=0):
    """
function to extrude a closed section to an open path
example:
sec=corner_radius(pts1([[0,0,.2],[3,0,.2],[0,2,1],[-3,0,1]]),10)
path=helix(10,2.5,4,5.01)
sol=path_extrude_open(sec,path)
fo(f'''
{swp(sol)}
color("blue")p_line3d({path},.1);
    ''')
    """
    if twist==0:
        p1=path[:-1]
        p2=path[1:]
        p1,p2=array([p1,p2])
        v1=p2-p1
        u1=v1/norm(v1,axis=1).reshape(-1,1)
        v2=concatenate([[u1[0]],(u1[1:]+u1[:-1])/2,[u1[-1]]])
        sec2=[]
        for i in range(len(path)):
            sec1=translate(path[i],sec2vector(v2[i],sec_1))
            sec2.append(sec1)
        return sec2
    elif twist==1:
        sec=[[1,0],[0,1],[-1,0],[0,-1]]
        t_l=tangents_along_path(path)[1:-1]
        v1=array([array(path)[1]-array(path)[0]]).tolist()+ \
        array([array(p[1])-array(p[0]) for p in t_l]).tolist()+ \
        array([array(path)[-1]-array(path)[-2]]).tolist()
        sol_1=[translate(path[i],sec2vector(v1[i],sec))  for i in range(len(path))]
        sol_2=align_sol(sol_1,360/len(t_l)/2)
        cp_1=array(sec_1).mean(0)
        sec_2=translate_2d(-cp_1,sec_1)
        sol_3=[translate(path[i],sec2vector(v1[i],sec_2))  for i in range(len(path))]
        sol_4=align_sol(sol_3,360/len(sol_3)/2)
        sol_5=[]
        for i in range(len(path)):
            a=mid_point([sol_1[i][0],sol_1[i][2]])
            v2=array(sol_2[i][0])-array(a)
            v3=array(sol_2[i][1])-array(a)
            sol_5.append(translate(-v3*cp_1[0],translate(-v2*cp_1[1],sol_4[i])))
        return sol_5

def offset(sec,r,type=1):
    """
function to offset a 2d closed loop section by distance 'r'
type 1 nad type 2 are 2 different types of offsets available and can be tried optionally if one type fails to perform. In most of the cases default type 1 should work
example:
sec=corner_radius(pts1([[-15,0,2.4],[0,15,3],[30,0,3],[0,-15,2.4],[5,0,2.4],[0,20,7],[-40,0,7],[0,-20,2.4]]),15)
os=linspace(-2.5,50,50)
sec1=[offset(sec,i) for i in os] #
fo(f'''
color("magenta")for(p={sec1})p_line3dc(p,.1,1);
color("blue")p_line3dc({sec},.1,1);
''') 
    """
    if type==1:
        return offset_1(sec,r)
    elif type==2:
        return offset_2(sec,r)

def sort_random_points(l_1,n_1,k=3):
    """
function to arrange random points in order
l_1: list of random points in space
n_1: is the normal vector to a plane from which all the points can be distinctly seen
k: is a factor which can have values >2 , default is 2, if the result is not satisfactory the values can be changed to see if the result is better. It has to be an integer
    
    """
    avg_1=array(l_1).mean(0).tolist()
    l_3=rot_sec2xy_plane(ppplane(l_1,n_1,avg_1))
    l_4=c3t2(l_3)
    l_5=concave_hull(l_4,k)
    l_6=array(l_1)[cKDTree(l_4).query(l_5)[1]].tolist()
    return l_6

def t_vec(path):
    """
    find the array of tangent vectors to a given path
    """
    p1=array(seg(path))
    p2=array(path)
    v1=array([(p[1]-p[0])/norm(p[1]-p[0]) for p in p1])
    t_v=array([ (v1[-1]+v1[i])/2 if i==0 else
         (v1[i-1]+v1[i])/2
        for i in range(len(p1))])
    t_v=t_v/norm(t_v,axis=1).reshape(-1,1)
    return t_v

def o_vec(path,n_v):
    """
    finds the array of orthogal vectors to a given path.
    normal vector at each point needs to be defined for this calculation
    """
    t_v=t_vec(path)
    o_v=array([cross(n_v[i],t_v[i]) for i in range(len(t_v))])
    o_v=o_v/norm(o_v,axis=1).reshape(-1,1)
    return o_v

def swp_sec(sec):
    """
    function to render a sec
    """
    n1=arange(len(sec)).tolist()
    
    return f'polyhedron({sec},{[n1]},convexity=10);'


def arc_with_start_pt_and_cp_3d(n1,start_point=[],center_point=[],theta=90,segments=30):
    """
    function to draw an arc with known center_point and start_point in 3d space
    n1 is the normal which defines the plane
    """
    center_point,start_point=array([center_point,start_point])
    v1=start_point-center_point
    arc_1=center_point+[axis_rot(n1,v1,i) for i in linspace(0,theta,segments)]
    arc_1=arc_1.tolist()
    return arc_1

def offset_3d(sec,d,type_of_offset=1):
    """
offsets an enclosed section in 3d space, in case the section is in 1 plane
sec: section in 3d space
d: offset distance -ve sign means inner offset and +ve sign is outer offset
refer to the file"example of various functions" for application examples
type_of_offset: offset type default is '1' in case of any issue in offset, try with '2'
example:
sec=corner_radius(pts1([[3,2,1],[8,3,3],[5,7,1],[-8,0,2],[-5,20,1]]),30)
sec1=o_solid([2,3,5],sec,1,10)[0]
sec2=offset_3d(sec1,-2)
fo(f'''
//original section
color("blue")p_line3dc({sec1},.1);
//offset section
color("magenta")p_line3dc({sec2},.1);
''')
    """
    sec1=rot_sec2xy_plane(sec)
    sec2=c32(sec1)
    sec3=offset(sec2,d,type_of_offset)
    a=a_(sec[0])
    b=a_(c23(c32(sec)[0]))
    v1=uv(best_fit_plane(sec)[0])
    if (l_(a_(v1).round(4))==[0,0,-1]) | (l_(a_(v1).round(4))==[0,0,1]):
        sec5=translate(a-b,sec3)
    else:
        v2=[0,0,-1]
        theta=l_(r2d(arccos(a_(v1)@v2)))
        ax_1=cross(v1,v2)
        sec4=axis_rot_1(c23(sec3),ax_1,c23(sec2)[0],-theta)
        sec5=translate(a-b,sec4)
    return sec5

def intersection_between_2_sketches(s1,s2):
    """
    finds the intersection between 2 2d sketches
    """
    a=seg(s1)
    b=seg(s2)
    c=a+b
    ip_1=s_int1(c)
    return ip_1

def pol(p1,l1):
    """
    point on line
    finds whether a point is on a line or not
    returns True or False
    """
    v1=l1[:-1]
    v2=l1[1:]
    v3=array(v2)-array(v1)
    v4=array(p1)-array(v1)
    u3=(v3/norm(v3,axis=1).reshape(-1,1)).round(4)
    u4=(v4/norm(v4,axis=1).reshape(-1,1)).round(4)
    dec=(array(p1).round(4)==array(l1[0]).round(4)).all()
    return True if dec==True else (u3[:,None,:]==u4[None,:]).all(2).any()



def i_p_p(surf_1,i_p_l,d=1.):
    """
function to project the intersection point on the cutting lines based on the distance 'r'
example:
a=cylinder(r=10,h=50,s=100)
b=plane([-1,0,1],[100,100],[0,0,15])
l1=ip_sol2sol(b,a)
l2=i_p_p(a,l1,5)
fo(f'''
%{swp_surf(b)}
%{swp(a)}
color("blue") p_line3dc({l1},.3);
color("magenta") p_line3dc({l2},.3);
''')
    """
    surf_1=cpo(surf_1) if d>0. else(cpo(flip(surf_1)))
    r_ipl=[]
    for j in range(len(surf_1)):
        for i in range(len(i_p_l)):
            if pol(i_p_l[i],surf_1[j]):
                a=cKDTree(surf_1[j]).query(i_p_l[i],2)[1]
                if pol(i_p_l[i],a_(surf_1[j])[a]) and a[0]<a[1]:
                    p1=equidistant_path([i_p_l[i]]+surf_1[j][a[1]:],pitch=abs(d))[1]
                elif pol(i_p_l[i],a_(surf_1[j])[a]) and a[1]<a[0]:
                    p1=equidistant_path([i_p_l[i]]+surf_1[j][a[0]:],pitch=abs(d))[1]
                elif pol(i_p_l[i],a_(surf_1[j])[[a[0],a[0]+1]]):
                    p1=equidistant_path([i_p_l[i]]+surf_1[j][a[0]+1:],pitch=abs(d))[1]
                elif pol(i_p_l[i],a_(surf_1[j])[[a[0],a[0]-1]]):
                    p1=equidistant_path([i_p_l[i]]+surf_1[j][a[0]:],pitch=abs(d))[1]
        
                r_ipl.append(p1)
    return path2path1(i_p_l,r_ipl) if len(r_ipl)!=len(i_p_l) else r_ipl


def best_fit_plane(pnts):
    """
returns a best fit plane approximation and intercept
example:
a=cylinder(r=10,h=50,s=100)
b=plane([-1,0,1],[100,100],[0,0,15])
l1=ip_sol2sol(b,a)
l2=i_p_p(a,l1,5)
n1,intcpt=best_fit_plane(l2)
pl1=plane(n1,[100,100],intcpt)
fo(f'''
//best fit plane
%{swp_surf(pl1)}
color("magenta") p_line3dc({l2},.3);
color("blue") p_line3d({vector2line(n1,10,intcpt)},.5);
''')
    """
    pnts=axis_rot_o([1,1,0],pnts,.00001)
    x2=(a_(pnts)[:,0]**2).sum()
    xy=(a_(pnts)[:,0]*a_(pnts)[:,1]).sum()
    x=(a_(pnts)[:,0]).sum()
    y2=(a_(pnts)[:,1]**2).sum()
    y=(a_(pnts)[:,1]).sum()
    hx=(a_(pnts)[:,2]*a_(pnts)[:,0]).sum()
    hy=(a_(pnts)[:,2]*a_(pnts)[:,1]).sum()
    h=(a_(pnts)[:,2]).sum()
    m1=l_([[x2,xy,x],[xy,y2,y],[x,y,len(pnts)]])
    m11=l_(a_(m1).transpose(1,0))
    m2=l_([hx,hy,h])
    a0,a1,a2=l_(inv(m1)@m2)
    return [[a0,a1,-1],l_(a_(pnts).mean(0))]

def path_extrude_over_multiple_sec_open(sec_1,path,twist=0):
    """
function to extrude multiple closed sections to an open path
    """
    if twist==0:
        p1=path[:-1]
        p2=path[1:]
        p1,p2=array([p1,p2])
        v1=p2-p1
        u1=v1/norm(v1,axis=1).reshape(-1,1)
        v2=concatenate([[u1[0]],(u1[1:]+u1[:-1])/2,[u1[-1]]])
        sec2=[]
        for i in range(len(path)):
            sec1=translate(path[i],sec2vector(v2[i],sec_1[i]))
            sec2.append(sec1)
        return sec2
    elif twist==1:
        sec=[[1,0],[0,1],[-1,0],[0,-1]]
        t_l=tangents_along_path(path)[1:-1]
        v1=array([array(path)[1]-array(path)[0]]).tolist()+ \
        array([array(p[1])-array(p[0]) for p in t_l]).tolist()+ \
        array([array(path)[-1]-array(path)[-2]]).tolist()
        sol_1=[translate(path[i],sec2vector(v1[i],sec))  for i in range(len(path))]
        sol_2=align_sol(sol_1,360/len(t_l)/2)
        sol_3=[]
        for i in range(len(path)):
            cp_1=array(sec_1[i]).mean(0)
            sec_2=translate_2d(-cp_1,sec_1[i])
            sol_3.append(translate(path[i],sec2vector(v1[i],sec_2)))
            
        sol_4=align_sol(sol_3,360/len(sol_3)/2)
        sol_5=[]
        for i in range(len(path)):
            a=mid_point([sol_1[i][0],sol_1[i][2]])
            v2=array(sol_2[i][0])-array(a)
            v3=array(sol_2[i][1])-array(a)
            sol_5.append(translate(-v3*cp_1[0],translate(-v2*cp_1[1],sol_4[i])))
        return sol_5

def nv2v_xy(v1):
    """
    returns a normal vector to any vector in x-y plane
    """
    u1=v1/norm(v1)
    ua=array([0,0,-1]) if u1[2]==0 else array([0,-1,0]) if (u1==[0,0,1]).all() else array([-1,0,0]) if (u1==[0,0,-1]).all() else array([u1[0],u1[1],0])
    v2=cross(u1,ua)
    u2=v2/norm(v2)
    u3=axis_rot(u2,u1,-90)
    u1,u2,u3=array([u1,u2,u3]).tolist()
    return u2

def nv2v_z(v1):
    """
    returns a normal vector to any vector in +/- z direction
    """
    u1=v1/norm(v1)
    ua=array([0,0,-1]) if u1[2]==0 else array([0,-1,0]) if (u1==[0,0,1]).all() else array([-1,0,0]) if (u1==[0,0,-1]).all() else array([u1[0],u1[1],0])
    v2=cross(u1,ua)
    u2=v2/norm(v2)
    u3=cross(u2,u1)
    u1,u2,u3=array([u1,u2,u3]).tolist()
    return u3

def s_p(p_l): # starting point
        """
        find the starting point for a convex hull
        bottom left point
        """
        l_1=array(p_l).round(5)
        a=l_1[l_1[:,1].argsort()]
        if len(a[a[:,1]==a[:,1].min()])>1:
            b=a[a[:,1]==a[:,1].min()]
            s_pnt=b[b[:,0].argsort()][0]
        else:
            s_pnt=a[0]
        return s_pnt.tolist()


def sec2surface_1(sec1,s=1):
    """
create an aligned surface from a section 'sec1'
considers the min sum method
example:
a=cr2dt([[0,0,2],[10,0,2],[0,10,2],[-10,0,2]],11)
s1=sec2surface_1(a)
fo(f'''
color("magenta") p_line3dc({a},.1);
color("blue") for(p={s1}) p_line3d(p,.05,1);
{swp_surf(c23(s1))}
''')
    """
    sec2=[sec2surface(sec1[i:]+sec1[:i],s) for i in range(len(sec1))]
    sec2=[sum([l_len(p1) for p1 in p]) for p in sec2]
    n2=array(sec2).argmin()
    sec3=sec2surface(sec1[n2:]+sec1[:n2],s)
    return sec3

def sec2surface_2(sec1,s=1):
    """
create an aligned surface from a section 'sec1'
considers the max sum method
example:
a=cr2dt([[0,0,2],[10,0,2],[0,10,2],[-10,0,2]],11)
s1=sec2surface_2(a)
fo(f'''
color("magenta") p_line3dc({a},.1);
color("blue") for(p={s1}) p_line3d(p,.05,1);
{swp_surf(c23(s1))}
''')
    """
    sec2=[sec2surface(sec1[i:]+sec1[:i],s) for i in range(len(sec1))]
    sec2=[sum([l_len(p1) for p1 in p]) for p in sec2]
    n2=array(sec2).argmax()
    sec3=sec2surface(sec1[n2:]+sec1[:n2],s)
    return sec3
    
def vector2length(v,l=10):
    """
    draw a defined vector to length 'l'
    example:
    vector2length([1,0,0],10) => [10.0,0.0,0.0]
    """
    u=array(v)/norm(v)
    v1=(u*l).tolist()
    return v1

def tangent_on_cir_from_pnt(c,p,l=1,side=0):
    """
function to draw a tangent on a circle 'c' from any given
point 'p' on the circle. length of the tangent 'l'
side can be set to either '0' or '1' to draw tangents in
2 different directions
example:
a=circle(20)
b=tangent_on_cir_from_pnt(a,a[10],15)
c=tangent_on_cir_from_pnt(a,a[10],15,side=1)
fo(f'''
color("magenta") p_line3dc({a},.3);
color("blue") p_line3d({b},.3);
color("orange") p_line3d({c},.3);
color("cyan") points({[a[10]]},1);
''')
    """
    v1=array(p)-array(cp_arc(c))
    theta1=r2d(arctan(l/norm(v1)))
    theta2=ang(v1[0],v1[1])
    l2=norm(v1)/cos(d2r(theta1))
    pa=translate_2d(cp_arc(c),rot2d(theta1+theta2,[l2,0])) if side==0 \
    else translate_2d(cp_arc(c),rot2d(-theta1+theta2,[l2,0]))
    return [p,pa]

def l_(a):
    """
    convert an array to list
    """
    return array(a).tolist()

def a_(l):
    """
    convert a list to array
    """
    return array(l)

def convert_to_triangles(sol):
    """
    convert a solid to a triangular mesh
    """
    f1=faces_1(len(sol),len(sol[0]))
    v1=vertices(sol)
    return l_(v1[f1])

def match_2_points_list(s0,s1):
    """
    match 2 sets of list of points or sections in space without loosing any point from both the lists
    
    """
    i=cKDTree(s1).query(s0)[1]
    j=cKDTree(s0).query(s1)[1]
    l_i=a_([arange(len(i)),i]).transpose(1,0)
    l_j=a_([j,arange(len(j))]).transpose(1,0)
    l_j,l_i
    s2=a_(exclude_points(l_j,l_i)+l_(l_i))
    
    s3=a_(sorted(s2,key=lambda x:(x[0],x[1]))).transpose(1,0)
    s4=l_(a_([a_(s0)[s3[0]],a_(s1)[s3[1]]]))
    return s4

def convert_to_triangles_surface(surf):
    """
    convert a surface to a triangular mesh
    """
    f1=faces_surface(len(surf),len(surf[0]))
    v1=vertices(surf)
    return l_(v1[f1])


def lexico(pnts=[],seq=[0,1,2],ord=[1,1,1]):
    """
lexicographic ordering of a points list
seq: defines the seduence in which the points needs to be ordered
e.g. [0,1,2] means first on x then on y and lastly on z coordinates
order: means asceding or descending order '1' means ascending and '-1' means descending order
e.g. [1,1,1] means all the coordinates should be in ascending order
    """
    if len(pnts[0])==2:        
        return sorted(pnts,key=lambda x:(ord[0]*x[seq[0]],ord[1]*x[seq[1]]))
    elif len(pnts[0])==3:
        return sorted(pnts,key=lambda x:(ord[0]*x[seq[0]],ord[1]*x[seq[1]],ord[2]*x[seq[2]]))

def fillet3points(sec,r=1,s=20):
    '''
    fillet from a list of 3 points
    sec: list of 3 points
    r: radius
    s: number of segments in the fillet
    '''
    p0=seg(sec)
    l1=p0[0]
    l2=p0[1]
    l3=path_offset(l1,abs(r)) if cw(sec)==1 else path_offset(l1,-abs(r))
    l4=path_offset(l2,abs(r)) if cw(sec)==1 else path_offset(l2,-abs(r))
    px=i_p2d(l3,l4)
    px1=vcost1(l1,px)
    px2=vcost1(l2,px)
    a1=arc_2p(px1,px2,r,cw(sec),s)
    return a1

def surface2sec(surf):
    """
reverse of sec2surface_1 function
    """
    a,b,c=array(surf).transpose(1,0,2).shape
    d=l_(array(surf).transpose(1,0,2))
    e=[d[i]+flip(d[-i-1]) for i in range(int(a/2))]
    return e


def bspline_open(pl,deg=3,s=100):
    """
draws bspline curve for control points list 'pl'
degree of curve 'deg' which should lie 0<=deg<len(p0)
s: number of points in the resultant curve
example:
p0=pts2([[0,0,0],[10,-2,2],[0,10,3],[5,-3,-6],[-5,5,1],[0,7,-2],[0,0,5]])
c1=bezier(p0,100)
s1=bspline_open(p0,2,100)
fo(f'''
color("blue")points({p0},.5);
color("cyan")p_line3d({c1},.05,rec=1);
color("green")p_line3d({s1},.05,rec=1);
''')
    """
    def t(i,k,n):
        return 0 if i<=k else i-k if (i>k and i<=n) else n-k+1

    def N(i,k,u,n,ak):

        if k>0:
            a1=(u-t(i,ak,n))
            b1=(t(i+k,ak,n)-t(i,ak,n))
            a2=(t(i+k+1,ak,n)-u)
            b2=(t(i+k+1,ak,n)-t(i+1,ak,n))
            a=(a1/b1 if b1!=0 else 0)*N(i,k-1,u,n,ak)+  \
               (a2/b2 if b2!=0 else 0)*N(i+1,k-1,u,n,ak)     
            return a
        elif k==0 and u>t(i,ak,n) and u<=t(i+1,ak,n):
            return 1
        else:
            return 0
    # tr=a_(pl[0])
    # p0=a_(translate(-tr,pl))
    p0=a_(pl)
    n=len(p0)-1
    k=deg
    p1=l_(a_([a_([p0[i]*N(i,k,u,n,ak=deg) for i in range(len(p0))]).sum(0) for u in linspace(0,n-k+1,s)]))
    return [pl[0]]+p1[1:]

def bspline_closed(pl,deg=3,s=100):
    """
draws bspline closed loop curve for control points list 'pl'
degree of curve 'deg' which should lie 0<=deg<len(p0)
s: number of points in the resultant curve
example:
px=pts2([[0,0,0],[5,2,3],[5,8,5],[-15,-3,1],[-3,5,-5],[3,5,6],[-13,-10,10]])
# b=bspline_open(px,2,100)
c=bspline_closed(px,2,100)
fo(f'''
color("blue")points({px},.5);
color("cyan")p_line3d({px},.1);
color("blue")p_line3d({c},.1);
''')
    """
    def t(i,k,n):
        return i-k if i<=k else i-k if (i>k and i<=n) else i-k

    def N(i,k,u,n,ak):

        if k>0:
            a1=(u-t(i,ak,n))
            b1=(t(i+k,ak,n)-t(i,ak,n))
            a2=(t(i+k+1,ak,n)-u)
            b2=(t(i+k+1,ak,n)-t(i+1,ak,n))
            a=(a1/b1 if b1!=0 else 0)*N(i,k-1,u,n,ak)+  \
               (a2/b2 if b2!=0 else 0)*N(i+1,k-1,u,n,ak)     
            return a
        elif k==0 and u>t(i,ak,n) and u<=t(i+1,ak,n):
            return 1
        else:
            return 0
    pl=pl+pl[:deg]
    # tr=a_(pl[0])
    # p0=a_(translate(-tr,pl))
    p0=a_(pl)
    n=len(p0)-1
    k=deg
    p1=l_(a_([a_([p0[i]*N(i,k,u,n,ak=deg) for i in range(len(p0))]).sum(0) for u in linspace(0,n-k+1,s)]))
    return p1[:-1]



def bspline_surface(pl,deg1=3,deg2=3,s1=100,s2=100,a=[1,1]):
    """
draws bspline surface from 2 control points list 'pl1' and 'pl2'
degree of curves are 'deg1' and 'deg2' in 2 directions 
which should lie 0<=deg1<len(pl[0]) for 'deg1' and
0<=deg2<=len(pl) for 'deg2'
s: number of points in the resultant curve
a: if a[0]==1 means set 1st shape of the surface to closed else to open and similarly for a[1]. open here means function bspline_open and closed means bspline_closed
example:
l1=cr2dt([[0,0],[10,0],[0,10],[10,0],[0,-10],[10,0],[0,10],[10,0],[0,-10],[10,0]])
l1=rot('x90',l1)
l2=rot('z90',l1)
s1=surface_from_2_waves(l1,l2,10)
s2=bspline_surface(s1,deg1=3,deg2=3,s1=50,s2=50,a=[0,0])
fo(f'''
color("magenta") p_line3d({l1},.3);
color("blue") p_line3d({l2},.3);
//color("orange") for(p={s1}) p_line3d(p,.1);
//{swp_surf(s1)}
{swp_surf(s2)}
''')
    """
    p1=[bspline_closed(p,deg1,s1) if a[0]==1 else bspline_open(p,deg1,s1) for p in pl]
    p2=cpo([bspline_closed(p,deg2,s2) if a[1]==1 else bspline_open(p,deg2,s2) for p in cpo(p1)])
    return p2

def bspline_sol(pl,deg1=3,deg2=3,s1=100,s2=100,a=[1,1]):
    """
    draws bspline surface from 2 control points list 'pl1' and 'pl2'
    degree of curves are 'deg1' and 'deg2' in 2 directions 
    which should lie 0<=deg1<len(pl[0]) for 'deg1' and
    0<=deg2<=len(pl) for 'deg2'
    s: number of points in the resultant curve
    a: if a[0]==1 means set 1st shape of the surface to closed else to open and similarly for a[1]. open here means function bspline_open and closed means bspline_closed
    """
    p1=[bspline_closed(p,deg1,s1) if a[0]==1 else bspline_open(p,deg1,s1) for p in pl]
    p2=cso([bspline_closed(p,deg2,s2) if a[1]==1 else bspline_open(p,deg2,s2) for p in cso(p1)])
    return p2

def polp(l,t):# point on line parameteric
    """
    find a point on line 'l' at parameter 't' where 0<=t<=1
    """
    def pol1(l,t):
        return l_(a_(l[0])*(1-t)+a_(l[1])*t)
    a=l_lenv_o(l)*t
    b=a_([0]+[l_len(p) for p in seg(l)[:-1]]).cumsum()
    if a<b[-1]:
        n=arange(len(b))[b<=a][-1]
        d1=a-b[n]
        d2=l_len(b[n:n+2])
        t1=d1/d2
        return pol1(l[n:n+2],t1)
    elif a==l_lenv_o(l):
        return l[-1]
    elif len(l)==2:
        return pol1(l,t)
    else:
        raise ValueError('value of t should be between 0 - 1')

def interpolate(p0,s=100):
    """
    create a smooth line joining all the points in the list p0
    """
    def l_basis(i,n,u):
        t=[j/(n-1) for j in range(n)]
        a=[]
        for j in range(n):
            a.append((u-t[j])/(t[i]-t[j]) if (t[i]-t[j])!=0 else 1)
        return prod(a).round(5)
    n=len(p0)
    p1=l_(a_([a_([a_(p0[i])*l_basis(i,n,u)  for i in range(n)]).sum(0) 
           for u in linspace(0,1,s)]))
    return p1


def lineFromPointTillEnd(l1,pnt,dist=.01):
    """
    draws line from point 'pnt' till end of line 'l1'
    the point with in distance 'dist' will be picked up
    deletes the rest of the line
    """
    def lfpte(l1,t):
        """
        removes the line till parameter t where 0<=t<=1
        keeps the rest of the line
        """
        p0=polp(l1,t)
        a=l_lenv_o(l1)*t
        b=a_([0]+[l_len(p) for p in seg(l1)[:-1]]).cumsum()
        if a<b[-1]:
            n=arange(len(b))[b<=a][-1]
        if n==len(l1)-2:
            l2=[p0,l1[n+1]]
        else:
            l2=[p0]+l1[n+1:]
        return l2
    return lfpte(l1,timeToReachPoint(pnt,l1,dist=dist))

def lineFromStartTillPoint(l1,pnt,dist=.01):
    """
    draws line from start till point 'pnt'.
    the point with in distance 'dist' will be picked up
    deletes the rest of the line
    """
    t=timeToReachPoint(pnt,l1,dist)    
    p0=polp(l1,t)
    a=l_lenv_o(l1)*t
    b=a_([0]+[l_len(p) for p in seg(l1)[:-1]]).cumsum()
    if a<b[-1]:
        n=arange(len(b))[b<=a][-1]
    if n==0:
        l2=[l1[0],p0]
    else:
        l2=l1[:n+1]+[p0]
    return l2

def vcost(l1,p0,dist=.2):
    """
    finds the projection of the point 'p0' on line 'l1' which is within distance 'dist' from the line
    """
    v1=a_(l1[1])-a_(l1[0])
    u1=v1/norm(v1)
    d=norm(v1)
    v2=a_(p0)-a_(l1[0])
    d1,d2=v1@v2/norm(v1),norm(cross(c23(v1),c23(v2)))/norm(v1)
    if d1>=0 and d1<=d and abs(d2)<=dist:
        p1=l_(a_(l1[0])+u1*d1)
        t1=d1/d
    else:
        p1=[]
    return p1


def timeToReachPoint(p0,l1,dist=.1):
    """
    p0: point
    l1: line
    dist: point within distance 'dist' of line will be picked
    it return time where 0<=time<=1
    """
    sec1=seg(l1)[:-1]
    n=0
    for i in range(len(sec1)):
        if vcost(sec1[i],p0,dist):
            n=i
    a=l_lenv_o(l1[:n+1])
    b=l_len([l1[n],p0])
    c=l_lenv_o(l1)
    return (a+b)/c


def movePointOnLine(l1,p0,d):
    """
    move a point 'p0' on line 'l1' by a distance 'd'
    """
    if p0!=l1[0]:
        t1=d/l_lenv_o(l1)
        t2=timeToReachPoint(p0,l1)
        t3=t1+t2
        p1=polp(l1,t3)
        return p1
    elif p0==l1[0]:
        t1=d/l_lenv_o(l1)
        return polp(l1,t1)

def ip_triangle_surf(ip,sol1):
    """
    function to find the triangles on the solid 'sol1' where the intersection points list 'ip' lies
    """
    v,f1=vnf1(sol1)
    tri=array(v)[array(f1)]
    p0,p1,p2=tri[:,0],tri[:,1],tri[:,2]
    p01,p02=p1-p0,p2-p0
    n1=cross(p01,p02)
    n1=n1/(norm(n1,axis=1).reshape(-1,1)+0)
    tri=tri[~((n1==[0,0,0]).all(1))]
    n1=n1[~((n1==[0,0,0]).all(1))]
    p0,p1,p2=tri[:,0],tri[:,1],tri[:,2]
    p01,p02=p1-p0,p2-p0
    la=array(ip)
    lab=n1

    iim=array([lab,-p01,-p02]).transpose(1,0,2).transpose(0,2,1)
    im=inv(iim)

    x=einsum('jkl,ijl->ijk',im,(p0-la[:,None]))
    t=x[:,:,0].round(3)
    u=x[:,:,1].round(3)
    v=x[:,:,2].round(3)
    decision=(t>=-0.01)&(t<=1.01)&(u>=-0.01)&(u<=1.01)&(v>=-0.01)&(v<=1.01)&((u+v)<=1.01)
    tri_1=array([tri[decision[i]][0] for i in range(len(ip))]).tolist()

    return tri_1

def i_p_n_surf(px,sol1):
    tri=array(ip_triangle_surf(px,sol1))
    p0,p1,p2=tri[:,0],tri[:,1],tri[:,2]
    p01,p02=p1-p0,p2-p0
    v3=cross(p01,p02)
    v3=v3/norm(v3,axis=1).reshape(-1,1)
    return v3

def c23(pl):
    """
convert a points list from 2d to 3d coordinate system
example:
c23([[0,0],[10,0]]) => [[0, 0, 0], [10, 0, 0]]
    """
    return c2t3(pl)

def c32(pl):
    """
convert a points list from 3d to 2d coordinate system
example:
c32([[0,0,0],[10,0,0]]) => [[0, 0], [10, 0]]
    """
    return c3t2(pl)

def prism2cpo(s1):
    """
change the orientation of points on the prism in such a way that all the open holes if any are plugged
example:
c1=circle(10)
p1=cr2dt([[3,0],[-3,0,3],[0,6,2],[-2,0]])
sol=prism(c1,p1)
sol1=translate([40,0,0],prism2cpo(sol))
fo(f'''
{swp_c(sol)}
{swp_c(sol1)}
''')
    """
    s2=[cpo(s1)[i]+flip(cpo(s1)[-(i+1)]) for i in range(int(len(s1[0])/2))]
    a=[mid_point([s2[0][i],s2[0][-(i+1)]]) for i in range(int(len(s2[0])/2))]
    b=[mid_point([s2[-1][i],s2[-1][-(i+1)]]) for i in range(int(len(s2[-1])/2))]
    
    a=a+flip(a)
    b=b+flip(b)
    s2=[a]+s2+[b]
    return s2


def psos(s2,s3,v1,dist=100000,unidirection=1):
    """
project a surface on to another without loosing the original points
surface 's3' will be projected on surface 's2'
'v1' is vector for projection
'dist' is the maximum distance through which projection can happen
unidirection: if the projection is to be in both direction set parameter
unidirection to '1' else to '0'
example:
s1=sphere(20)
s2=c23(a_(grid2d(circle(10),2)).reshape(11,11,2))
v1=[1,0,1]
s3=psos(c_(s1),s2,v1)

fo(f'''
//surface to project upon
%{swp(s1)}
// surface to project
{swp_surf(s2)}
//projected surface
color("cyan"){swp_surf(s3)}
// vector to project the surface
color("blue") p_line3d({vector2line(v1,10)},.3);
''')
    """
    p0=a_(s3).reshape(-1,3)
    f=faces_surface(len(s2),len(s2[0]))
    v=a_(s2).reshape(-1,3)
    tri=v[f]
    p2,p3,p4=tri[:,0],tri[:,1],tri[:,2]
    n1=a_([v1]*len(p2))
    v2,v3=p3-p2,p4-p2
    
    iim=a_([n1,-v2,-v3+.0000001]).transpose(1,0,2).transpose(0,2,1)+.000001
    im=inv(iim)
    px=[]
    for i in range(len(p0)):
        t=(im@(p2-p0[i][None,:])[:,:,None]).reshape(-1,3)
        t1,t2,t3=t[:,0],t[:,1],t[:,2]
        if unidirection==0:
            dec=(t2>=0)&(t2<=1)&(t3>=0)&(t3<=1)&((t2+t3)<=1)
        elif unidirection==1:
            dec=(t1>=0)&(t2>=0)&(t2<=1)&(t3>=0)&(t3<=1)&((t2+t3)<=1)
        if dec.any()==1 and norm(a_(v1)*sorted(t1[arange(len(p2))[dec]],key=abs)[0])<=dist:
            px.append(p0[i]+a_(v1)*sorted(t1[arange(len(p2))[dec]],key=abs)[0])
        elif dec.any()==0 or norm(a_(v1)*sorted(t1[arange(len(p2))[dec]],key=abs)[0])>dist:
            px.append(p0[i])

    px=l_(a_(px).reshape(len(s3),len(s3[0]),3))
    return px


def psos_v(s2,s3,v1,dist=100000,unidirection=0):
    """
project a surface on to another without loosing the original points
surface 's3' will be projected on surface 's2'
'v1' is vector for projection. this is a focal vector 
from where the rays are emitted for projection
example:
s1=sphere(20)
s2=c23(a_(grid2d(circle(10),2)).reshape(11,11,2))
v1=[1,0,1]
s3=psos(c_(s1),s2,v1)

fo(f'''
//surface to project upon
%{swp(s1)}
// surface to project
{swp_surf(s2)}
//projected surface
color("cyan"){swp_surf(s3)}
// vector to project the surface
color("blue") p_line3d({vector2line(v1,10)},.3);
''')
    """
    p0=a_(s3).reshape(-1,3)
    f=faces_surface(len(s2),len(s2[0]))
    v=a_(s2).reshape(-1,3)
    tri=v[f]
    p2,p3,p4=tri[:,0],tri[:,1],tri[:,2]
    
    px=[]
    for i in range(len(p0)):
        n1=a_([uv(p0[i]-v1)]*len(p2))
        v2,v3=p3-p2,p4-p2
        iim=a_([n1,-v2,-v3+.0000001]).transpose(1,0,2).transpose(0,2,1)+.000001
        im=inv(iim)
        # im.shape,p0[198].shape
        t=(im@(p2-a_(p0[i])[None,:])[:,:,None]).reshape(-1,3)
        t1,t2,t3=t[:,0],t[:,1],t[:,2]
        if unidirection==0:
            dec=(t2>=0)&(t2<=1)&(t3>=0)&(t3<=1)&((t2+t3)<=1)
        elif unidirection==1:
            dec=(t1>=0)&(t2>=0)&(t2<=1)&(t3>=0)&(t3<=1)&((t2+t3)<=1)

        if dec.any()==1 and norm(a_(n1[0])*sorted(t1[arange(len(p2))[dec]],key=abs)[0])<=dist:
            px.append(a_(p0[i])+a_(n1[0])*sorted(t1[arange(len(p2))[dec]],key=abs)[0])
        elif dec.any()==0 or norm(a_(n1[0])*sorted(t1[arange(len(p2))[dec]],key=abs)[0])>dist:
            px.append(p0[i])
    
    px=l_(a_(px).reshape(len(s3),len(s3[0]),3))
    return px

# def psos_v_1(s2,s3,v1,vx,dist=100000,unidirection=0):
#     """
#     project a surface on to another without loosing the original points
#     surface 's3' will be projected on surface 's2'
#     'v1' is vector for projection. this is a focal vector
#     'vx' is a vector to define plane through which each vector works for projecting
#     from where the rays are emitted for projection
#     """
#     p0=a_(s3).reshape(-1,3)
#     f=faces_surface(len(s2),len(s2[0]))
#     v=a_(s2).reshape(-1,3)
#     tri=v[f]
#     p2,p3,p4=tri[:,0],tri[:,1],tri[:,2]
    
#     px=[]
#     for i in range(len(p0)):
#         v4=ppplane([v1],vx,p0[i])[0]
#         n1=a_([uv(v4-p0[i])]*len(p2))
#         v2,v3=p3-p2,p4-p2
#         iim=a_([n1,-v2,-v3+.0000001]).transpose(1,0,2).transpose(0,2,1)+.000001
#         im=inv(iim)

#         t=(im@(p2-a_(p0[i])[None,:])[:,:,None]).reshape(-1,3)
#         t1,t2,t3=t[:,0],t[:,1],t[:,2]
#         if unidirection==0:
#             dec=(t2>=0)&(t2<=1)&(t3>=0)&(t3<=1)&((t2+t3)<=1)
#         elif unidirection==1:
#             dec=(t1>=0)&(t2>=0)&(t2<=1)&(t3>=0)&(t3<=1)&((t2+t3)<=1)

#         if dec.any()==1 and norm(a_(n1[0])*sorted(t1[arange(len(p2))[dec]],key=abs)[0])<=dist:
#             px.append(a_(p0[i])+a_(n1[0])*sorted(t1[arange(len(p2))[dec]],key=abs)[0])
#         elif dec.any()==0 or norm(a_(n1[0])*sorted(t1[arange(len(p2))[dec]],key=abs)[0])>dist:
#             px.append(p0[i])
    
#     px=l_(a_(px).reshape(len(s3),len(s3[0]),3))
#     return px

def psos_v_1(s2,s3,l1,dist=100000,unidirection=0):
    """
project a surface on to another without loosing the original points
surface 's3' will be projected on surface 's2'
'l1' is a line for projection
example:
s1=sphere(20)
s2=c23(a_(grid2d(circle(10),2)).reshape(11,11,2))
l1=point_vector([-10,0,-10],[20,0,0]) # point from where the rays are emitting for projection
s3=psos_v_1(c_(s1),s2,l1,unidirection=1)
gd=a_(s3).reshape(-1,3)
vx=[vcost1(l1,p) for p in gd]
vx=cpo([vx,gd])
fo(f'''
//surface to project upon
%{swp(s1)}
// surface to project
{swp_surf(s2)}
//projected surface
color("cyan"){swp_surf(s3)}
// rays projecting out of source
color("blue") for(p={vx}) p_line3d(p,.03);
//light source
color("magenta") p_line3d({l1},.5);
''')
    """
    p0=a_(s3).reshape(-1,3)
    f=faces_surface(len(s2),len(s2[0]))
    v=a_(s2).reshape(-1,3)
    tri=v[f]
    p2,p3,p4=tri[:,0],tri[:,1],tri[:,2]
    p1=a_([vcost1(l1,p) for p in p0])
    v1=p0-p1
    px=[]
    for i in range(len(p0)):
        n1=a_([uv(v1[i])]*len(p2))
        v2,v3=p3-p2,p4-p2
        iim=a_([n1,-v2,-v3+.0000001]).transpose(1,0,2).transpose(0,2,1)+.000001
        im=inv(iim)
        # im.shape,p0[198].shape
        t=(im@(p2-a_(p0[i])[None,:])[:,:,None]).reshape(-1,3)
        t1,t2,t3=t[:,0],t[:,1],t[:,2]
        if unidirection==0:
            dec=(t2>=0)&(t2<=1)&(t3>=0)&(t3<=1)&((t2+t3)<=1)
        elif unidirection==1:
            dec=(t1>=0)&(t2>=0)&(t2<=1)&(t3>=0)&(t3<=1)&((t2+t3)<=1)

        if dec.any()==1 and norm(a_(n1[0])*sorted(t1[arange(len(p2))[dec]],key=abs)[0])<=dist:
            px.append(a_(p0[i])+a_(n1[0])*sorted(t1[arange(len(p2))[dec]],key=abs)[0])
        elif dec.any()==0 or norm(a_(n1[0])*sorted(t1[arange(len(p2))[dec]],key=abs)[0])>dist:
            px.append(p0[i])
    
    px=l_(a_(px).reshape(len(s3),len(s3[0]),3))
    return px


def normal_vector(sec):
    """
    calculates the normal vector of a given sec
    """
    sec=c23(sec)
    if len(sec)<3:
        raise ValueError("To calculate Normal number of points should be atleast 3 ")
    elif len(sec)==3:
        d=a_(sec)
        n3=l_(cross(d[1]-d[0],d[2]-d[0]))
    else:
        e=sec2surface_1(sec,1)
        v=a_(e).reshape(-1,3)
        f=faces_surface(len(e),len(e[0]))
        tri=v[f]
        
        p0,p1,p2=tri[:,0],tri[:,1],tri[:,2]
        v1,v2=p1-p0,p2-p0
        n1=cross(v1,v2)
        n1=a_([n1,n1,n1]).transpose(1,0,2).reshape(-1,3)
        g=concatenate(f)
        n1.shape,p0.shape,g.shape
        n2=a_([n1[arange(len(g))[(g==i)]].mean(0) for i in range(len(v))])
        n2=n2/norm(n2,axis=1).reshape(-1,1)
        n3=n2.mean(0)
        n3=l_(n3/norm(n3))
    return n3


def surface_normals(surf,direction=1):
    """
    calculates normals from each point on the surface
    normals may need to be reshaped to the shape of the 'surf'
    """
    f=faces_surface(len(surf),len(surf[0]))
    v=a_(surf).reshape(-1,3)
    tri=v[f]
    p0,p1,p2=tri[:,0],tri[:,1],tri[:,2]
    v1,v2=p1-p0,p2-p0
    n1=cross(v1,v2)
    n1=a_([n1,n1,n1]).transpose(1,0,2).reshape(-1,3)
    g=concatenate(f)
    n1.shape,p0.shape,g.shape
    n2=a_([n1[arange(len(g))[(g==i)]].mean(0) for i in range(len(v))])+.00001
    n2=n2/norm(n2,axis=1).reshape(-1,1)*direction

    return n2

def psos_n(s2,s3,direction=1,dist=100000):
    """
project a surface on to another without loosing the original points
surface 's3' will be projected on surface 's2'
the projection is based on the normal to the surface s3 and the direction of normals can be changed from '1' to '-1'
    """
    p0=a_(s3).reshape(-1,3)
    f=faces_surface(len(s2),len(s2[0]))
    v=a_(s2).reshape(-1,3)
    tri=v[f]
    p2,p3,p4=tri[:,0],tri[:,1],tri[:,2]
    v1=l_(surface_normals(s3,direction))
    px=[]
    for i in range(len(p0)):
        n1=a_([v1[i]]*len(p2))
        v2,v3=p3-p2,p4-p2
        iim=a_([n1,-v2,-v3+.0000001]).transpose(1,0,2).transpose(0,2,1)+.000001
        im=inv(iim)
        # im.shape,p0[198].shape
        t=(im@(p2-a_(p0[i])[None,:])[:,:,None]).reshape(-1,3)
        t1,t2,t3=t[:,0],t[:,1],t[:,2]
        dec=(t1>=0)&(t2>=0)&(t2<=1)&(t3>=0)&(t3<=1)&((t2+t3)<=1)

        if dec.any()==1 and norm(a_(n1[0])*sorted(t1[arange(len(p2))[dec]],key=abs)[0])<=dist:
            px.append(a_(p0[i])+a_(n1[0])*sorted(t1[arange(len(p2))[dec]],key=abs)[0])
        elif dec.any()==0 or norm(a_(n1[0])*sorted(t1[arange(len(p2))[dec]],key=abs)[0])>dist:
            px.append(p0[i])
    
    px=l_(a_(px).reshape(len(s3),len(s3[0]),3))
    return px

def psos_n_b(s2,s3,dist=100000):
    """
project a surface on to another without loosing the original points
surface 's3' will be projected on surface 's2'
the projection is based on the normal to the surface s3 and the direction of normals can be changed from '1' to '-1'
    """
    p0=a_(s3).reshape(-1,3)
    f=faces_surface(len(s2),len(s2[0]))
    v=a_(s2).reshape(-1,3)
    tri=v[f]
    p2,p3,p4=tri[:,0],tri[:,1],tri[:,2]
    v1=l_(surface_normals(s3,1))
    px=[]
    for i in range(len(p0)):
        n1=a_([v1[i]]*len(p2))
        v2,v3=p3-p2,p4-p2
        iim=a_([n1,-v2,-v3+.0000001]).transpose(1,0,2).transpose(0,2,1)+.000001
        im=inv(iim)

        t=(im@(p2-a_(p0[i])[None,:])[:,:,None]).reshape(-1,3)
        t1,t2,t3=t[:,0],t[:,1],t[:,2]
        dec=(t2>=0)&(t2<=1)&(t3>=0)&(t3<=1)&((t2+t3)<=1)

        if dec.any()==1 and norm(a_(n1[0])*sorted(t1[arange(len(p2))[dec]],key=abs)[0])<=dist:
            px.append(a_(p0[i])+a_(n1[0])*sorted(t1[arange(len(p2))[dec]],key=abs)[0])
        elif dec.any()==0 or norm(a_(n1[0])*sorted(t1[arange(len(p2))[dec]],key=abs)[0])>dist:
            px.append(p0[i])
    
    px=l_(a_(px).reshape(len(s3),len(s3[0]),3))
    return px

def cog(sol):
    """
    calculate the center of gravity
    """
    return l_(a_(sol).reshape(-1,3).mean(0))



def interpolation_bspline_closed(p0,s=50,f=3.425):
    """
Interpolates through points
p0: original points
s: number of points in the resultant curve
f: factor of smoothness, lower values may result in too wavyness and larger values
are as good as connecting points in straight line
example:
l1=cr2dt([[0,0],[10,0],[-3,7],[15,-5],[0,10],[-10,10],[-15,-5]])
# l2=interpolation_bspline_open(l1,50)
l3=interpolation_bspline_closed(l1,50)
fo(f'''
color("blue") points({l1},.5);
//color("magenta") p_line3d({l2},.2);
color("cyan") p_line3d({l3},.2);
''')
    """
    p2=[ mid_point(p) for p in seg(p0)]
    l1=[ l_len(p)/f for p in seg(p0)]
    l2=[a_(p).mean() for p in seg(l1)]
    p3=p0+[p0[0]]
    pa=a_(p3[:-1])
    pb=a_(p3[1:])
    v1=pb-pa
    v2=[ a_(uv(a_(p).mean(0))) for p in seg(v1)]
    p3=a_([ [a_(p3[i+1])-v2[i]*l2[i],a_(p3[i+1])+v2[i]*l2[i]] for i in range(len(v2))])
    p4=l_(concatenate(p3))
    p5=bspline_closed(p4,2,s)
    return p5[:-1]

def interpolation_bspline_open(p0,s=50,f=3.425):
    """
Interpolates through points
p0: original points
s: number of points in the resultant curve
f: factor of smoothness, lower values may result in too wavyness and larger values
are as good as connecting points in straight line
example:
l1=cr2dt([[0,0],[10,0],[-3,7],[15,-5],[0,10],[-10,10],[-15,-5]])
l2=interpolation_bspline_open(l1,50)
#l3=interpolation_bspline_closed(l1,50)
fo(f'''
color("blue") points({l1},.5);
color("magenta") p_line3d({l2},.2);
//color("cyan") p_line3d({l3},.2);
''')
    
    """
    p2=[ mid_point(p) for p in seg(p0)[:-1]]
    l1=[ l_len(p)/f for p in seg(p0)[:-1]]
    l2=[a_(p).mean() for p in seg(l1)[:-1]]
    pa=a_(p0[:-1])
    pb=a_(p0[1:])
    v1=pb-pa
    v2=[ a_(uv(a_(p).mean(0))) for p in seg(v1)[:-1]]
    p3=a_([ [a_(p0[i+1])-v2[i]*l2[i],a_(p0[i+1])+v2[i]*l2[i]] for i in range(len(v2))])
    p4=[p0[0]]+l_(concatenate(p3))+[p0[-1]]
    p5=bspline_open(p4,2,s)
    return p5

def interpolation_surface_open(pl,f1=3.425,f2=3.425,s1=100,s2=100):
    """
draws bspline interpolation surface from 2 control points list 'pl1' and 'pl2'
'f1' and 'f2' are factors of smoothness in 2 directions 
s1 and s2: number of points in 2 direction in the resultant curve
example:
l1=cr2dt([[0,0],[10,0],[3,10],[10,0],[3,-10],[10,0],[3,10],[10,0],[3,-10],[10,0]])
l1=rot('x90',l1)
l2=rot('z90',l1)
s1=surface_from_2_waves(l1,l2,10)
s2=interpolation_surface_open(s1,f1=10,f2=10,s1=50,s2=50)
fo(f'''
color("magenta") p_line3d({l1},.3);
color("blue") p_line3d({l2},.3);
//color("orange") for(p={s1}) p_line3d(p,.1);
//{swp_surf(s1)}
{swp_surf(s2)}
''')
    """
    p1=[interpolation_bspline_open(p,s1,f1) for p in pl]
    p2=cpo([interpolation_bspline_open(p,s2,f2) for p in cpo(p1)])
    return p2

def bezier_surface(pl,s1=100,s2=100):
    """
draws bezier surface from control points list 'pl'
s1 and s2: number of points in the resultant curve in 2 direction
example:
l1=cr2dt([[0,0],[10,0],[3,10],[10,0],[3,-10],[10,0],[3,10],[10,0],[3,-10],[10,0]])
l1=m_points1_o(rot('x90',l1),10)
l2=m_points1_o(rot('z90',l1),10)
s1=surface_from_2_waves(l1,l2,10)
s2=bezier_surface(s1,s1=50,s2=50)
fo(f'''
color("magenta") p_line3d({l1},.3);
color("blue") p_line3d({l2},.3);
//color("orange") for(p={s1}) p_line3d(p,.1);
//{swp_surf(s1)}
{swp_surf(s2)}
''')
    """
    p1=[bezier(p,s1) for p in pl]
    p2=cpo([bezier(p,s2) for p in cpo(p1)])
    return p2



def il_fillet(il,sol1,sol2,r1,r2,s=20,o=0,type=1,dist=0,vx=[],style=2,f=1):
    """
    calculates a fillet at the intersection of 2 solids.
    r1 and r2 would be same in most of the cases, but the signs can be different depending on which side the fillet is required
    r1 is the distance by which intersection line offsets on sol1 and similarly r2 is on sol2 
    for type, dist and vx parameters refer function o_3d_rev
    """
        
    p1=il
    p2=o_3d_rev(p1,sol2,r2,type=type,dist=dist,vx=vx,edges_closed=1)
    p3=o_3d_rev(p1,sol1,r1,type=type,dist=dist,vx=vx,edges_closed=1)
    fillet1=convert_3lines2fillet(p3,p2,p1,s=s,style=style,f=f)
    
    return fillet1

def il_fillet_surf(il,surf1,sol2,r1,r2,s=20,o=0,type=1,dist=0,vx=[],style=2,f=1):
    """
    calculates a fillet at the intersection of surface with solid.
    r1 and r2 would be same in most of the cases, but the signs can be different depending on which side the fillet is required
    r1 is the distance by which intersection line offsets on sol1 and similarly r2 is on sol2 
    for type, dist and vx parameters refer function o_3d_rev
    """
        
    p1=il
    p2=o_3d_rev(p1,sol2,r2,type=type,dist=dist,vx=vx,edges_closed=1)
    p3=o_3d_rev(p1,surf1,r1,type=type,dist=dist,vx=vx,edges_closed=0)
    fillet1=convert_3lines2fillet(p3,p2,p1,s=s,style=style,f=f)
    
    return fillet1

def line2length(l1,length=1):
    """
change the length of a line to the parameter distance 'length' keeping the initial point same
    """
    v1=a_(l1[1])-a_(l1[0])
    u1=v1/norm(v1)
    return [l1[0],l_(l1[0]+u1*length)]

def ip_surf2line(surf,line):# when line has more than 2 points
    """
    function to calculate intersection point between a 3d surface and a line. 
     "surf" is the 3d object which is intersected with a "line".
    
    """


    pa=surf
    p1=array([[ [[pa[i][j],pa[i][j+1],pa[i+1][j]],[pa[i+1][j+1],pa[i+1][j],pa[i][j+1]]]
     for j in range(len(pa[i])-1)] 
              for i in range(len(pa)-1)]).reshape(-1,3,3)
    pm=p1[:,0]
    pn=p1[:,1]
    po=p1[:,2]
    px=array(line[:-1])
    py=array(line[1:])
    v1,v2,v3=py-px,pn-pm,po-pm
    a,_=v1.shape
    b,_=v2.shape
    v1=array([v1]*b)
    v2=-array([v2]*a).transpose(1,0,2)
    v3=-array([v3]*a).transpose(1,0,2)
    iim=array([v1,v2,v3]).transpose(1,2,0,3).transpose(0,1,3,2)+.00001
    im=inv(iim)
    p=array([pm]*a).transpose(1,0,2)-array([px]*b)
    t=einsum('ijkl,ijl->ijk',im,p)
    condition=(t[:,:,0]>=0)&(t[:,:,0]<=1)&(t[:,:,1]>=0)&(t[:,:,1]<=1)&(t[:,:,2]>=0)&(t[:,:,2]<=1)&((t[:,:,1]+t[:,:,2])<=1)
    t1=t[:,:,0][condition]
    i_p1=array([px]*b)[condition]+einsum('ij,i->ij',v1[condition],t1)
    t1=a_([timeToReachPoint(p,line) for p in i_p1]).argsort()
    i_p1=l_(a_(i_p1)[t1])
    return i_p1

def line_as_axis(l1):
    """
    convert line to an axis vector
    """
    v1=l_(a_(l1[1])-a_(l1[0]))
    return v1


def psos_v_2(s2,s3,v1,dist=100000,unidirection=0):
    """
    project a surface on to another without loosing the original points
    surface 's3' will be projected on surface 's2'
    'v1' are the vectors for projection. v1 are the list of vectors
    """
    p0=a_(s3).reshape(-1,3)
    f=faces_surface(len(s2),len(s2[0]))
    v=a_(s2).reshape(-1,3)
    tri=v[f]
    p2,p3,p4=tri[:,0],tri[:,1],tri[:,2]
    
    px=[]
    for i in range(len(p0)):
        n1=a_([uv(v1[i])]*len(p2))
        v2,v3=p3-p2,p4-p2
        iim=a_([n1,-v2,-v3+.0000001]).transpose(1,0,2).transpose(0,2,1)+.000001
        im=inv(iim)
        # im.shape,p0[198].shape
        t=(im@(p2-a_(p0[i])[None,:])[:,:,None]).reshape(-1,3)
        t1,t2,t3=t[:,0],t[:,1],t[:,2]
        if unidirection==0:
            dec=(t2>=0)&(t2<=1)&(t3>=0)&(t3<=1)&((t2+t3)<=1)
        elif unidirection==1:
            dec=(t1>=0)&(t2>=0)&(t2<=1)&(t3>=0)&(t3<=1)&((t2+t3)<=1)

        if dec.any()==1 and norm(a_(n1[0])*sorted(t1[arange(len(p2))[dec]],key=abs)[0])<=dist:
            px.append(a_(p0[i])+a_(n1[0])*sorted(t1[arange(len(p2))[dec]],key=abs)[0])
        elif dec.any()==0 or norm(a_(n1[0])*sorted(t1[arange(len(p2))[dec]],key=abs)[0])>dist:
            px.append(p0[i])
    
    px=l_(a_(px).reshape(len(s3),len(s3[0]),3))
    return px

def rot(a,sol):
    """
    function to rotate an object by axis and angles defined by a string 'a'
    example:
    rot(a='z30x70y-90',sol=[[0,0],[10,0]])
    will rotate line [[0,0],[10,0]] first through z-axis by 30 deg 
    then through x-axis by 70 deg and then through y-axis by -90 deg 
    """
    b=[i for i,p in enumerate(a) if p=='x']
    c=[i for i,p in enumerate(a) if p=='y']
    d=[i for i,p in enumerate(a) if p=='z']
    
    e=l_(sort(b+c+d))+[len(a)]
    th=[l_(a_(a[p[0]+1:p[1]]).astype(float)) for p in seg(e)[:-1]]
    ax=[[1,0,0] if a[i]=='x' else [0,1,0] if a[i]=='y' else [0,0,1] for i in e[:-1]]
    for i in range(len(ax)):
        sol=axis_rot(ax[i],sol,th[i])
    return sol

def rot2d(a,sec):
    """
    function to rotate a 2d section by a defined angle 'a'
    """
    return c32(axis_rot([0,0,1],sec,a))


def sub_d_2(b,a=[1/4,3/4]):
    """
    input to function fillet_by_subdivison
    """
    c=[[[polp(p1,a[0]),polp(p1,a[1])] for p1 in seg(p)[:-1]] for p in cpo(b)]
    n1,n2,n3,n4=a_(c).shape
    c=[b[0]]+cpo(a_(c).reshape(n1,n2*n3,n4))+[b[-1]]
    
    return c

def fillet_by_subdivison(sol,a=[1/4,3/4]):
    """
    produce fillet by subdividing an object
    """
    for _ in range(4):
        sol=sub_d_2(sol,a)
    return sol

def cso(sol):
    """
    change the orientation of points in solid
    """
    a=cpo(sol)
    b=[ a[i]+flip(a[-(i+1)]) for i in range(int(len(a)/2))]
    return b


def offset_3d_sec(line1,d=1):
    l1=line1
    s1=sec2surface(l1)
    s2=surface_offset(s1,1)
    l2=surface2sec(s2)[0]
    s3=cpo([l1,l2])
    s4=surface_offset(s3,d)
    s5=cpo(s4)[0]
    return s5


def remove_points_within_dist_closed(path,line,dist=1):
    """
    function to remove points which are with in defined distance of a given path
    path : is the original path
    line: is the line from where the points needs to be removed
    dist: is the distance with in which all the points to be excluded
    """

    l2=line
    p0=a_(path)
    p1=a_(path[1:]+path[:1])
    v1=p1-p0
    l3=[]
    for i in range(len(l2)):
        p2=a_(l2[i])
        v2=p2-p0
        a,b=einsum('ij,ij->i',v1,v2),norm(v1,axis=1)
        a=a/b
        c=norm(cross(v1,v2),axis=1)/b
        # d=norm(v2,axis=1)
        if ~(c[arange(len(path))[((a>=0) & (a<=b))]]<dist).any():
            l3.append(l2[i])
    return l3

def list_remove_points_within_dist_closed(path,line,dist=1):
    """
    function to find list of points to remove which are with in defined distance of a given path
    path : is the original path
    line: is the line from where the points needs to be removed
    dist: is the distance with in which all the points to be excluded
    """

    l2=line
    p0=a_(path)
    p1=a_(path[1:]+path[:1])
    v1=p1-p0
    l3=[]
    for i in range(len(l2)):
        p2=a_(l2[i])
        v2=p2-p0
        a,b=einsum('ij,ij->i',v1,v2),norm(v1,axis=1)
        a=a/b
        c=norm(cross(v1,v2),axis=1)/b
        # d=norm(v2,axis=1)
        if ~(c[arange(len(path))[((a>=0) & (a<=b))]]<dist).any():
            l3.append(i)
    return l3

def remove_points_within_dist_open(path,line,dist=1):
    """
    function to remove points which are with in defined distance of a given path
    path : is the original path
    line: is the line from where the points needs to be removed
    dist: is the distance with in which all the points to be excluded
    """

    l2=line
    p0=a_(path)[:-1]
    p1=a_(path)[1:]
    v1=p1-p0
    l3=[]
    for i in range(len(l2)):
        p2=a_(l2[i])
        v2=p2-p0
        a,b=einsum('ij,ij->i',v1,v2),norm(v1,axis=1)
        a=a/b
        c=norm(cross(v1,v2),axis=1)/b
        # d=norm(v2,axis=1)
        if ~(c[arange(len(path)-1)[((a>=0) & (a<=b))]]<dist).any():
            l3.append(l2[i])
    return l3

def list_remove_points_within_dist_open(path,line,dist=1):
    """
    function to find list of points to remove which are with in defined distance of a given path
    path : is the original path
    line: is the line from where the points needs to be removed
    dist: is the distance with in which all the points to be excluded
    """

    l2=line
    p0=a_(path)[:-1]
    p1=a_(path)[1:]
    v1=p1-p0
    l3=[]
    for i in range(len(l2)):
        p2=a_(l2[i])
        v2=p2-p0
        a,b=einsum('ij,ij->i',v1,v2),norm(v1,axis=1)
        a=a/b
        c=norm(cross(v1,v2),axis=1)/b
        # d=norm(v2,axis=1)
        if ~(c[arange(len(path)-1)[((a>=0) & (a<=b))]]<dist).any():
            l3.append(i)
    return l3

def ppos(surf,pnt,vect,unidirection=1):#project point on surface
    """
    function to project a point on a surface
    'vect' is the vector in which direction the line would be projected.
    'unidirection' to be set to '1' if the projection is required in the 
    direction of the vector and set to '0' in case it is required in either of 
    the direction. 
    Projections happens at the nearest location on the surface from a line
    """
    return psos(surf,[[pnt]],vect,unidirection=unidirection)[0][0]

def plos(surf,line,vect,unidirection=1):#project line on surface
    """
    function to project a line on a surface
    'vect' is the vector in which direction the line would be projected.
    'unidirection' to be set to '1' if the projection is required in the 
    direction of the vector and set to '0' in case it is required in either of 
    the direction. 
    Projections happens at the nearest location on the surface from a line
    """
    return psos(surf,[line],vect,unidirection=unidirection)[0]

def barycentric_normals(surf,edges_closed=1):
    """
    finds normals for each triangle mesh of a surface
    if it is a closed loop surface like cylinder, 
    parameter edges_closed to be set to '1' else set to '0'
    """
    if edges_closed==1:
        f=faces_1(len(surf),len(surf[0]))
    elif edges_closed==0:
        f=faces_surface(len(surf),len(surf[0]))
    v=a_(surf).reshape(-1,3)
    triangles=v[f]
    p0=triangles[:,0]
    p1=triangles[:,1]
    p2=triangles[:,2]
    v1=p1-p0
    v2=p2-p0
    n1=cross(v1,v2)
    n1=n1/norm(n1,axis=1).reshape(-1,1)
    return n1

def barycenter(surf,edges_closed=1):
    """
    finds center for each triangle mesh of a surface
    if it is a closed loop surface like cylinder, 
    parameter edges_closed to be set to '1' else set to '0'
    """
    if edges_closed==1:
        f=faces_1(len(surf),len(surf[0]))
    elif edges_closed==0:
        f=faces_surface(len(surf),len(surf[0]))
    v=a_(surf).reshape(-1,3)
    triangles=v[f]
    cp=triangles.mean(1)
    return cp



def smoothening_by_subdivison(sec,iterations=4,closed=0):
    """
    smoothen a polyline by subdivision method.
    if the polyline is closed loop, parameter 'closed' should be set at '1'
    else if the polyline is not a closed loop parameter 'closed' to be set to '0'
    """
    for i in range(iterations):
        if closed==1:
            sec=l_(concatenate([[polp(p1,1/4),polp(p1,3/4)] for p1 in seg(sec)]))
        elif closed==0:
            sec=[sec[0]]+l_(concatenate([[polp(p1,1/4),polp(p1,3/4)] for p1 in seg(sec)[:-1]]))+[sec[-1]]
    return sec



def smoothening_by_subdivison_surf(sol,iterations=4,o=[0,0]):
    """
    smoothen a solid with sub-divison method
    iterations: number of iterations
    o: if o[0]==0 (means loop remains open for the 1st orientation, and '1' means the loop remains closed for the 1st orientation of the shape) 
    similarly o[1]==0 (means loop remains open for the 2nd orientation, and '1' means the loop remains closed for the 2n orientation of the shape)
    an example can make this more clear
    """
    sol=[smoothening_by_subdivison(p,iterations,o[0]) for p in sol]
    sol=cpo([smoothening_by_subdivison(p,iterations,o[1]) for p in cpo(sol)])
    return sol

def wrap_x(l1,path):
    """
same as function wrap_around but in x-direction
example:
    c1=circle(10,[10.1,0],s=100)
p1=rot('x90',sinewave(25,3,1,100))
c2=wrap_x(c1,p1)
fo(f'''
// original section
color("magenta",.2) p_line3d({c1},.2);
//path
color("cyan",.2) p_line3d({p1},.2);
// section wrapped over path
color("blue") p_line3dc({c2},.3);
''')
    
    """
    
    l2=[l_len(p) for p in seg(path)[:-1]]
    l2=a_([0]+l_(a_(l2).cumsum()))
    l3=a_(l1)[:,0]
    l4=[[0,p[1],p[2]] for p in c23(l1)]
    l5=[[0,p[1],0] for p in c23(l1)]
    t2=[]
    for i in range(len(l3)):
        n=arange(len(l2))[l2<l3[i]][-1]
        t1=(l3[i]-l2[n])/l_len(seg(path)[n])
        p0=seg(path)[n][0]
        p1=seg(path)[n][1]
        p2=a_(p0)*(1-t1)+a_(p1)*t1
        v1=line_as_axis([p0,p1])
        u1=v1/norm(v1)
        u2=[0,sign(l4[i][2])*1,0]
        d1=cross(u1,u2)*abs(l4[i][2])
        t2.append(translate(p2+d1,l5[i]))

    return t2

def wrap_y(l1,path):
    
    l2=[l_len(p) for p in seg(path)[:-1]]
    l2=a_([0]+l_(a_(l2).cumsum()))
    l3=a_(l1)[:,1]
    l4=[[p[0],0,p[2]] for p in c23(l1)]
    l5=[[p[0],0,0] for p in c23(l1)]
    t2=[]
    for i in range(len(l3)):
        n=arange(len(l2))[l2<l3[i]][-1]
        t1=(l3[i]-l2[n])/l_len(seg(path)[n])
        p0=seg(path)[n][0]
        p1=seg(path)[n][1]
        p2=a_(p0)*(1-t1)+a_(p1)*t1
        v1=line_as_axis([p0,p1])
        u1=v1/norm(v1)
        u2=[sign(l4[i][2])*-1,0,0]
        d1=cross(u1,u2)*abs(l4[i][2])
        t2.append(translate(p2+d1,l5[i]))

    return t2

def thicken_surface(surf,t=1):
    """
    gives thickness to a surface by amount 't'
    """
    a=surf
    b=surface_offset(surf,t)
    sol=[a[i]+flip(b[i]) for i in range(len(a))]
    return sol

def find_points_beyond_distance_of_surface(surf,pnts,d=1,edges_closed=1):
    c=surf
    if edges_closed==1:
        f1=faces_1(len(c),len(c[0]))
    elif edges_closed==0:
        f1=faces_surface(len(c),len(c[0]))
    vert1=a_(c).reshape(-1,3)
    trngl=vert1[f1]
    p0,p1,p2=trngl[:,0],trngl[:,1],trngl[:,2]
    v1,v2=p1-p0,p2-p0
    v0=barycentric_normals(c)
    l0=a_(pnts)
    # l0+v0*t0=p0+v1*t1+v2*t2
    iim=a_([v0,-v1,-v2]).transpose(1,0,2).transpose(0,2,1)
    im=inv(iim)
    pa=[]
    # pb=[]
    for i in range(len(l0)):
        p=p0-l0[i]
        t0,t1,t2=einsum('ijk,ik->ij',im,p).transpose(1,0)
        dec=(t1>=0)&(t1<=1)&(t2>=0)&(t2<=1)&((t1+t2)<=1)
        px=norm(einsum('ij,i->ij',v0,t0),axis=1)
        if (px[dec]>d).all() and l_(px[dec])!=[]:
            pa.append(pnts[i])
            # pb.append(px[dec].min())
    return pa

def points_inside_solid(pnts,surf,edges_closed=1):
    """
    finds all the points which are inside the solid
    """
    c=surf
    if edges_closed==1:
        f1=faces_1(len(c),len(c[0]))
    elif edges_closed==0:
        f1=faces_surface(len(c),len(c[0]))
    a1,a2,a3=a_([[1,0,0],[0,1,0],[0,0,1]])
    l0=a_(pnts)
    # f1=faces_1(len(c),len(c[0]))
    vert1=a_(c).reshape(-1,3)
    trngl=vert1[f1]
    p0,p1,p2=trngl[:,0],trngl[:,1],trngl[:,2]
    v1,v2=p1-p0,p2-p0
    # l0+a1*t0=p0+v1*t1+v2*t2
    # l0+a2*t0=p0+v1*t1+v2*t2
    # l0+a3*t0=p0+v1*t1+v2*t2
    iim1=a_([a_([a1]*len(v1)),-v1,-v2+.000001]).transpose(1,0,2).transpose(0,2,1)+.000001
    im1=inv(iim1)
    iim2=a_([a_([a2]*len(v1))+.000001,-v1,-v2]).transpose(1,0,2).transpose(0,2,1)+.000001
    im2=inv(iim2)
    iim3=a_([a_([a3]*len(v1))+.000001,-v1,-v2]).transpose(1,0,2).transpose(0,2,1)+.000001
    im3=inv(iim3)
    pb=[]
    for i in range(len(l0)):
        
        p=p0-l0[i]
        t0,t1,t2=einsum('ijk,ik->ij',im1,p).transpose(1,0)
        dec=(t0>=0)&(t1>=0)&(t1<=1)&(t2>=0)&(t2<=1)&((t1+t2)<=1)
        x1=ones(len(v1))[dec]
        
        t0,t1,t2=einsum('ijk,ik->ij',im2,p).transpose(1,0)
        dec=(t0>=0)&(t1>=0)&(t1<=1)&(t2>=0)&(t2<=1)&((t1+t2)<=1)
        x2=ones(len(v1))[dec]
        
        t0,t1,t2=einsum('ijk,ik->ij',im3,p).transpose(1,0)
        dec=(t0>=0)&(t1>=0)&(t1<=1)&(t2>=0)&(t2<=1)&((t1+t2)<=1)
        x3=ones(len(v1))[dec]
        if (a_([x1.sum()%2,x2.sum()%2,x3.sum()%2])==1).sum()>=2:
            pb.append(pnts[i])
    return pb   


def swp_triangles(sol):
    """
    function to render triangular mesh
    """
    vert1=l_(a_(sol).reshape(-1,3))
    f1=l_(arange(len(sol)*len(sol[0])).reshape(-1,3))
    return f'polyhedron({vert1},{f1},convexity=10);'

def surface_split(sol,v1):
    """
    function to split a 3d mesh surface, based on vector 
    example, if you are looking at the object from the top use vector [0,0,1]
    if from right [1,0,0] from left [-1,0,0], from front [0,1,0] , from back [0,-1,0] and so on...
    """
    p0,p1,p2=a_(sol)[:,0],a_(sol)[:,1],a_(sol)[:,2]
    p01=p1-p0
    p02=p2-p0
    a=sol.mean(1)
    b=[]
    for i in range(len(sol)):
        la=a[i]+a_(v1)*.00001
        lb=a[i]+a_(v1)
        lab=lb-la
        x1=cross(p01,p02)
        x2=cross(p02,-lab)
        x3=cross(-lab,p01)
        t=einsum('ij,ij->i',x1,la-p0)/einsum('j,ij->i',-lab,x1)
        u=einsum('ij,ij->i',x2,la-p0)/einsum('j,ij->i',-lab,x1)
        v=einsum('ij,ij->i',x3,la-p0)/einsum('j,ij->i',-lab,x1)
        dec=(t>=0)&(u>=0)&(u<=1)&(v>=0)&(v<=1)&((u+v)<=1)
        if ~dec.any():
            b.append(i)
    
    sol1=sol[a_(b)]
    return sol1

def surface_reshape_with_2lines_and_pnt(l1,l2,pnt,s=20):
    """
reshape a surface with a point 'pnt'
example:
l1=rot('x90',arc_2p([0,0],[50,0],50,s=50))
l2=translate([0,50,0],rot('z30',l1))
p0=[25,25,-5]
surf1=surface_reshape_with_2lines_and_pnt(l1,l2,p0,50)
fo(f'''
color("blue") for(p={[l1,l2]}) p_line3d(p,.3);
{swp_surf(surf1)}
color("magenta") points({[p0]},1);
''')
    """
    p0,p1=l1[cKDTree(l1).query(pnt)[1]],l2[cKDTree(l1).query(pnt)[1]]
    l3=[p0,pnt,p1]
    d1=len(l1)
    d2=cKDTree(l1).query(pnt)[1]
    a1=linspace(0,d2,d2)*90/d2
    a2=linspace(0,d1-d2,(d1-d2))*90/(d1-d2)
    a=a_(l_(sin(d2r(a1)))+l_(cos(d2r(a2))))
    # a=linspace(0,len(l1),len(l1))*180/len(l1)
    pnt1=pnt-a_(l1[cKDTree(l1).query(pnt)[1]])
    b=[]
    for i in range(len(l1)):
        b.append([l1[i],l_(a_(l1[i])+[a[i]*pnt1[0],a[i]*pnt1[1],a[i]*pnt1[2]]),l2[i]])
    b=[bspline_open(p,2,s) for p in b]
    return b

def lineFromPointToPointOnLine(l1,p0,p1,dist=.01):
    """
    Draw a line from a defined point to another point on a line.
    points with in distance 'dist' will be picked up
    Rest of the line will be deleted
    
    """
    if p0==l1[0]:
        l2=[]
    else:
        l2=lineFromStartTillPoint(l1,p0,dist)
    if p1==l1[-1]:
        l3=l1
    else:
        l3=lineFromStartTillPoint(l1,p1,dist)
    l4=[l2[-1]]+(exclude_points(l3,l2) if p0!=l1[0] else l3)
    return l4

def distanceOfPointFromLine(p0,l1):
    """
    calculates the perpendicular distance of a point from a polyline
    """
    px=[vcost(p,p0,1e5) for p in seg(l1)[:-1] if vcost(p,p0,1e5)!=[]]
    return l_len([px[0],p0])  if px!=[] else []

def perpendicularProjectionOfPointOnLine(p0,l1):
    """
    finds the perpendicular projection of a point on a polyline, if one exists
    """
    px=[vcost(p,p0,1e5) for p in seg(l1)[:-1] if vcost(p,p0,1e5)!=[]]
    return px[0]  if px!=[] else []

def switch_orientation(sol1):
    """
switch the orientation of the solid
example:
sol=cylinder(r=10,h=20)
sol1=switch_orientation(sol)
fo(f'''
{swp_surf(sol1)}
''')
    """
    if len(sol1[0])%2==1:
        sol1=[equidistant_path(p,len(sol1[0])) for p in sol1]
    sol2=a_(sol1).transpose(1,0,2)
    x,y,z=a_(sol2).shape
    c=a_(sol2).reshape(2,int(x/2),y,3)[0]
    d=a_(sol2).reshape(2,int(x/2),y,3)[1][::-1]
    d=d[::-1].transpose(1,0,2)[::-1].transpose(1,0,2)[::-1]
    return l_(a_([c,d]).transpose(1,0,2,3).reshape(-1,y*2,z))

def i_line_fillet(intersection_line,solid_1,solid_2,distance_1=1,distance_2=1,segments=20):
    """
    if the intersection line is defined between 2 solids. fillet can be drawn with the
    information
    """
    l1,sol1,sol2,r1,r2=intersection_line,solid_1,solid_2,distance_1,distance_2
    l2=o_3d_rev(l1,sol1,r1)
    l3=o_3d_rev(l1,sol2,r2)
    f1=convert_3lines2fillet(l2,l3,l1,s=segments)
    return f1

def i_line_fillet_closed(intersection_line,solid_1,solid_2,distance_1=1,distance_2=1,segments=20):
    """
    if the intersection line is defined between 2 solids. fillet (Closed loop) can be drawn with the
    information
    """
    l1,sol1,sol2,r1,r2=intersection_line,solid_1,solid_2,distance_1,distance_2
    l2=o_3d_rev(l1,sol1,r1)
    l3=o_3d_rev(l1,sol2,r2)
    f1=convert_3lines2fillet_closed(l2,l3,l1,s=segments)
    return f1

def mid_line(l1,l2):
    """
    draw a line in center of 2 lines 'l1' and 'l2'
example:
# complex case study
l1=cr3dt([[-10,-5,0,5],[0,10,0,5],[20,0,0,5],[0,-10,0,5],[0,0,10,5],[0,10,0,5],
         [-20,0,0,5],[0,-10,0,5]],20)
l1=homogenise(l1,1,1)
l2=point_vector([0,10,5],[1,0,0])
l2=[vcost1(l2,p) for p in l1]
l3=mid_line(l1,l2)
fo(f'''
//original line
color("blue") p_line3d({l1},.3);
// line 2
color("magenta") p_line3d({l2},.3);
//mid line
color("cyan") p_line3d({l3},.2);
''')
    """
    return [mid_point(p) for p in cpo([l1,l2])]

def list_ang(l1):
    """
    list of angles for each point for a closed loop section
    """
    a=[]
    for i in range(len(l1)):
        if i==0:
            p1=l1[i]
            p0=l1[-1]
            p2=l1[i+1]
        elif i==len(l1)-1:
            p1=l1[i]
            p0=l1[i-1]
            p2=l1[0]
        else:
            p1=l1[i]
            p0=l1[i-1]
            p2=l1[i+1]
        p0,p1,p2=a_([p0,p1,p2])
        p1p2=p2-p1
        p0p1=p1-p0
        u1=p1p2/norm(p1p2)
        u2=p0p1/norm(p0p1)
        a.append(l_(r2d(arccos(u1@u2))))
    return a

def surface_correction_after_offset_closed(surf_original,surf_off,dist=1):
    """
    remove all the points from the offset surface which are less than the offset distance 'dist'
    """
    sol2=[]
    for i in range(len(surf_original)):
        d=remove_points_within_dist_closed(surf_original[i],surf_off[i],dist-.01)
        if len(d)==len(surf_original[i]):
            sol2.append(d)
        else:
            d=path2path1(surf_original[i],d)
            sol2.append(d)
    return sol2

def surface_correction_after_offset_open(surf_original,surf_off,dist=1):
    """
    remove all the points from the offset surface which are less than the offset distance 'dist'
    """
    sol2=[]
    for i in range(len(surf_original)):
        d=remove_points_within_dist_open(surf_original[i],surf_off[i],dist-.01)
        if len(d)==len(surf_original[i]):
            sol2.append(d)
        else:
            d=path2path1(surf_original[i],d)
            sol2.append(d)
    return sol2

def plos_v_2(s2,l1,v1,dist=100000,unidirection=0):
    """
project a line on to a surface without loosing the original points
line 'l1' will be projected on surface 's2'
'v1' are the vectors for projection. v1 are the list of vectors
example:
l1=cr3dt([[-10,-5,0,5],[0,10,0,5],[20,0,0,5],[0,-10,0,5],[0,0,10,5],[0,10,0,5],
         [-20,0,0,5],[0,-10,0,5]],20)
l1=homogenise(l1,1,1)
l2=scl3dc(l1,.7)
v1=lines2vectors(cpo([l2,l1]))
s1=sphere(15)
l3=plos_v_2(c_(s1),l1,v1,unidirection=1)
l4=cpo([l2,l3])
fo(f'''
%{swp(s1)}
//original line
color("blue") p_line3d({l1},.3);
// projected line
color("magenta") p_line3d({l3},.3);
//source of light
color("cyan") p_line3d({l2},.2);
// rays for projection
color("grey") for(p={l4}) p_line3d(p,.05);
''')
    """
    return psos_v_2(s2,[l1],v1,dist=dist,unidirection=unidirection)[0]

def plos_v_1(s2,l1,l2,dist=100000,unidirection=0):
    """
project a line on to a surface without loosing the original points
line 'l1' will be projected on surface 's2'
'l2' is the line from where the rays are emitting to project  
example:
s1=sphere(10)
l1=m_points1_o(point_vector([0,1,0],[5,5,0]),10)
l2=point_vector([0,-5,-2],[5,0,0])
l3=plos_v_1(s1,l1,l2)
x1=[ vcost1(l2,p) for p in l1]
fo(f'''
%{swp(c_(s1))}
// line to project
color("blue") p_line3d({l1},.2);
// light source
color("cyan") p_line3d({l2},.2);
// projected line
color("magenta") p_line3d({l3},.2);
//rays of projection
color("grey") for(p={cpo([x1,l3])}) p_line3d(p,.03);
''')
    """
    return psos_v_1(s2,[l1],l2,dist=dist,unidirection=unidirection)[0]


def plos_v(s2,l1,v1,dist=100000,unidirection=0):
    """
project a line on to a surface without loosing the original points
line 'l1' will be projected on surface 's2'
'v1' is vector for projection. this is a focal vector 
from where the rays are emitted for projection
example:
s1=sphere(10)
l1=m_points1_o(point_vector([0,1,0],[5,5,0]),10)
v1=[0,0,-5]
l2=plos_v(s1,l1,v1)
fo(f'''
%{swp(c_(s1))}
color("blue") p_line3d({l1},.2);
color("cyan") p_line3d({l2},.2);
color("magenta") points({[v1]},.5);
//rays of projection
color("grey") for(p={l2}) p_line3d([{v1},p],.03);
''')
    """
    return psos_v(s2,[l1],v1,dist=dist,unidirection=unidirection)[0]

def lines2vectors(lines):
    """
convert lines to vectors
example:
c1=circle(10)
l1=lines2vectors(seg(c1))
fo(f'''
color("blue") p_line3d({c1},.2);
color("cyan") for(p={l1})p_line3d([[0,0],p],.02);
''')
    """
    return [line_as_axis(p) for p in lines]

def mirror_point(pnt,n1,loc):
    """
function to mirror a point 'pnt' defined by mirroring plane 'n1' passing through intercept 'loc'
example:
p0=[10,5,5]
p1=mirror_point(p0,n1=[0,1,0],loc=[0,0,0])
p2=mirror_point(p0,n1=[1,0,0],loc=[0,0,0])
p3=mirror_point(p0,n1=[0,0,1],loc=[0,0,0])
fo(f'''
color("blue") points({[p0]},.5);
color("cyan") points({[p1]},.5);
color("magenta") points({[p2]},.5);
color("grey") points({[p3]},.5);
''')
    """
    return mirror_line([pnt],n1,loc)[0]

def corner_radius3d(pnts,s=5): # Corner radius 3d where 'pnts' are the list of points with 4th coordinate in each point is radius 'rds' and 's' is number of segments for each arc
    """
Corner radius 3d where 'pnts' are the list of points with 4th coordinate in each point is radius and 's' is number of segments for each radial arc
example:
# cr3d is the short name of this function
l1=cr3d([[0,0,0,5],[0,10,0,5],[20,10,0,5],[20,0,0,5],[20,0,10,5],[20,10,10,5],
         [0,10,10,5],[0,0,10,5]],20)
fo(f'''
color("blue") p_line3d({l1},.3);
''')
    """
    rds=[pnts[i][3] if len(pnts[i])==4 else 0 for i in range(len(pnts))]
    pnts=[pnts[i][:3] for i in range(len(pnts))]
    c=[]
    for i in range(len(pnts)):
        if i==0:
            p0=pnts[len(pnts)-1]
            p1=pnts[i]
            p2=pnts[i+1]
        elif i<len(pnts)-1:
            p0=pnts[i-1]
            p1=pnts[i]
            p2=pnts[i+1]
        else:
            p0=pnts[i-1]
            p1=pnts[i]
            p2=pnts[0]
        c.append(fillet_3p_3d(p0,p1,p2,rds[i],s)[1:])
    c=array(c).reshape(-1,3).tolist()
    return remove_extra_points(array(c).round(5))

cr3d=corner_radius3d

def twoCircleCrossTangent(c1,c2,cw=-1): # two circle cross tangent
    """
function to draw cross tangent between 2 circles
example:
cir1=circle(10)
cir2=circle(5,[15,6])
p0=twoCircleCrossTangent(cir1,cir2)
fo(f'''
color("blue") for(p={[cir1,cir2]}) p_line3d(p,.3);
color("magenta") points({p0},.5);
%color("cyan",.2) p_line3d({p0},.3);
''')
    """
    r1=r_arc(c1)
    r2=r_arc(c2)
    cp1=cp_arc(c1)
    cp2=cp_arc(c2)
    v1=[1,1]
    v2=[-r2,r1]
    cp1,cp2=array([cp1,cp2])
    d=norm(cp2-cp1)
    d1=(inv(array([v1,v2]).T)@array([d,0]))[0]
    d2=(inv(array([v1,v2]).T)@array([d,0]))[1]
    a=arcsin(r1/d1)*180/pi
    v3=cp2-cp1
    u3=v3/norm(v3)
    b=arccos(u3@array([1,0]))*180/pi
    if cw==-1:
        if v3[0]>0 and v3[1]<=0:
            theta1=270+a-b
            theta2=90+a-b
        elif v3[0]>=0 and v3[1]>0:
            theta1=270+a+b
            theta2=90+a+b
        elif v3[0]<0 and v3[1]>=0:
            theta1=270+a+b
            theta2=90+a+b
        else:
            theta1=270+a-b
            theta2=90+a-b
    else:
        if v3[0]>0 and v3[1]<=0:
            theta2=270-a-b
            theta1=90-a-b
        elif v3[0]>=0 and v3[1]>0:
            theta2=270-a+b
            theta1=90-a+b
        elif v3[0]<0 and v3[1]>=0:
            theta2=270-a+b
            theta1=90-a+b
        else:
            theta2=270-a-b
            theta1=90-a-b
        
    p0=(cp1+array([r1*cos(theta1*pi/180),r1*sin(theta1*pi/180)])).tolist()
    p1=(cp2+array([r2*cos(theta2*pi/180),r2*sin(theta2*pi/180)])).tolist()
    return [p0,p1]

two_circle_cross_tangent=twoCircleCrossTangent

def twoCircleTangentPoints(c1,c2,side=0): #2 circle tangent point full (both the sides)
    """
function to draw tangent line joining 2 circles 'c1' and 'c2'.
this works counter-clockwise
This function draws tangent line on both the sides
example:
cir1=circle(10)
cir2=circle(5,[15,6])
p0=twoCircleTangentPoints(cir1,cir2)
fo(f'''
color("blue") for(p={[cir1,cir2]}) p_line3d(p,.3);
color("magenta") points({p0},.5);
%color("cyan",.2) p_line3d({p0},.3);
''')
    """
    r1=r_arc(c1)
    r2=r_arc(c2)
    cp1=cp_arc(c1)
    cp2=cp_arc(c2)
    cp1,cp2=array([cp1,cp2])
    v1=cp2-cp1,
    u1=v1/norm(v1)
    ang1=arcsin((r2-r1)/norm(cp2-cp1))*180/pi

    t1=cp1+u1@rm(90+ang1)*r1
    t2=cp2+u1@rm(90+ang1)*r2

    t3=cp1+u1@rm(-90-ang1)*r1
    t4=cp2+u1@rm(-90-ang1)*r2
    t5=[t3[0].tolist(),t4[0].tolist(),t2[0].tolist(),t1[0].tolist()]
    return t5[:2] if side==0 else t5[2:]

two_circle_tangent=twoCircleTangentPoints

def pts3(pl):
    """
this functionis required for 3d turtle movement and works with function 'corner_radius3d'
example:
a=pts3([[0,0,0],[10,7,5,5],[0,0,20,5],[0,10,0]]) => 
[[0, 0, 0, 0], [10, 7, 5, 5], [10, 7, 25, 5], [10, 17, 25, 0]]
    """
    a=pl
    b=[p if len(p)==4 else p+[0] for p in a]
    c=a_(b)[:,:3].cumsum(0)
    d=a_([a_(b)[:,3]])
    return l_(concatenate([c,d.T],axis=1))

def line_as_vector(line):
    """
convert a line to vector
example:
line_as_vector([[0,0],[10,0]]) => [10,0]
    """
    return line_as_axis(line)

def line_as_unit_vector(line):
    """
convert a line to unit vector
example:
line_as_unit_vector([[0,0],[10,0]]) => [1.0,0.0]
    """
    v1=line_as_vector(line)

    return l_(v1/norm(v1))


def derivative2d(line):
    b=seg(line)[:-1]
    c=[line_as_unit_vector(p) for p in b]
    d=[p[1]/p[0] for p in c]
    x=l_(a_(line)[:,0])[:-1]
    e=cpo([x,d])
    return e

def derivative3d(line):
    b=seg(line)[:-1]
    c=[line_as_unit_vector(p) for p in b]
    d=l_([p[2]/norm(p[:2]) for p in c])
    x=l_(a_(line)[:,:2][:-1])
    e=[x[i]+[d[i]] for i in range(len(x))]
    return e

def points_to_meshes(pnts,voxel_size=.1,iso_level_percentile=2,flip=0):
    """
create triangle meshes from a list of points using marching cube method
example:
# for this to work library "skimage" needs to be imported
from skimage import measure
l1=c23(circle(5))
fo(f'''
color("blue") p_line3d({l1},.3);
{points_to_meshes(l1,.5,10)}
''')
    """
    pnts=a_(pnts)
    mins=min(pnts,axis=0)-3*voxel_size
    maxs=max(pnts,axis=0)+3*voxel_size
    x = arange(mins[0], maxs[0], voxel_size)
    y = arange(mins[1], maxs[1], voxel_size)
    z = arange(mins[2], maxs[2], voxel_size)
    x, y, z = meshgrid(x, y, z, indexing='ij')
    grid_points = vstack([x.ravel(), y.ravel(), z.ravel()]).T
    distances=cKDTree(pnts).query(grid_points)[0]
    scalar_field=distances.reshape(x.shape)
    iso_level = percentile(distances, iso_level_percentile)
    verts, faces, _, _ = measure.marching_cubes(scalar_field, level=iso_level, allow_degenerate=False)
    if flip==0:
        faces=faces
    elif flip==1:
        faces=a_([p[::-1] for p in faces])
    verts = verts * voxel_size + mins
    return f'polyhedron({l_(verts)},{l_(faces)},convexity=10);'

def iso_surfaces(pnts,level_size=1):
    """
    create various iso levels to divide the points
    """
    ss=level_size
    m1=a_(pnts).min(axis=0)-1
    m2=a_(pnts).max(axis=0)+1
    ns=int(round((m2[2]-m1[2])/ss,0))
    a=[[m1,[m1[0],m2[1],m1[2]]],[[m1[0],m1[1],m2[2]],[m1[0],m2[1],m2[2]]]]
    a=m_points1_o(a,ns)
    b=translate([m2[0]-m1[0],0,0],a)
    s1=cpo([a,b])
    return s1

def iso_surfaces3d(pnts,level_size=1):
    """
    create various iso levels to divide the points
    """
    ss=level_size
    m1=a_(pnts).min(axis=0)
    m2=a_(pnts).max(axis=0)
    # ns=int(round((m2[2]-m1[2])/ss,0))
    x0,y0,z0=l_(m1)
    x1,y1,z1=l_(m2)
    ns1=int(round((x1-x0)/ss,0))
    ns2=int(round((y1-y0)/ss,0))
    ns3=int(round((z1-z0)/ss,0))
    
    l1=m_points1_o([[x0,y0,z0],[x1,y0,z0]],ns1,.00001)
    l2=m_points1_o([[x0,y1,z0],[x1,y1,z0]],ns1,.00001)
    s1=m_points1_o([l1,l2],ns2,.00001)
    l1=m_points1_o([[x0,y0,z1],[x1,y0,z1]],ns1,.00001)
    l2=m_points1_o([[x0,y1,z1],[x1,y1,z1]],ns1,.00001)
    s2=m_points1_o([l1,l2],ns2,.00001)
    surface=m_points1_o([s1,s2],ns3,.00001)
    return surface

def ip_sol2sol_each_line(sol1,sol2,n=0):
    line=array([ seg(p)[:-1] for p in cpo(sol2)])
    v,f1=vnf2(sol1)
    tri=array(v)[array(f1)]
    line=array([ seg(p)[:-1] for p in cpo(sol2)])
    tri.shape,line.shape
    la,lb=line[:,:,0],line[:,:,1]
    p0,p1,p2=tri[:,0],tri[:,1],tri[:,2]
    lab=lb-la
    p01,p02=p1-p0,p2-p0
    t=einsum('kl,ijkl->ijk',cross(p01,p02),la[:,:,None]-p0)/(einsum('ijl,kl->ijk',(-lab),cross(p01,p02))+.00000)
    u=einsum('ijkl,ijkl->ijk',cross(p02[None,None,:,:],(-lab)[:,:,None,:]),(la[:,:,None,:]-p0[None,None,:,:]))/(einsum('ijl,kl->ijk',(-lab),cross(p01,p02))+.00000)
    v=einsum('ijkl,ijkl->ijk',cross((-lab)[:,:,None,:],p01[None,None,:,:]),(la[:,:,None,:]-p0[None,None,:,:]))/(einsum('ijl,kl->ijk',(-lab),cross(p01,p02))+.00000)
    condition=(t>=0)&(t<=1)&(u>=0)&(u<=1)&(v>=0)&(v<=1)&(u+v<1)

    a=(la[:,None,:,None,:]+lab[:,None,:,None,:]*t[:,None,:,:,None])
    b=condition[:,None,:,:]
    c=[]
    for i in range(len(a)):
        c.append(a[i][b[i]].tolist())

    return c

def grid2d(pnts,resolution=1,offset=0):
    """
creates a grid of points to cover a 2d shape
example:
l1=circle(5)
l2=grid2d(l1,1,offset=2)
fo(f'''
color("blue") p_line3d({l1},.3);
color("magenta") points({l2},.1);
//%color("cyan",.2) p_line3d({l2},.05);
''')
    """
    ss=resolution
    m1=a_(pnts).min(axis=0)-offset
    m2=a_(pnts).max(axis=0)+offset
    
    x0,y0=l_(m1)
    x1,y1=l_(m2)

    ns1=int(round((x1-x0)/ss,0))
    ns2=int(round((y1-y0)/ss,0))
    l1=m_points1_o([[x0,y0],[x1,y0]],ns1,.00001)
    l2=m_points1_o([[x0,y1],[x1,y1]],ns1,.00001)
    s1=m_points1_o([l1,l2],ns2,.00001)
    
    return l_(a_(s1).reshape(-1,2))

def ang_v(v):
    """
finds angle of a vector
example:
ang_v([4,4]) => 45.0
    """
    return ang(v[0],v[1])


def sterguss(n):
    return 3 if n<=200 else round(1+log2(n))

def h_lines(sec,n=10,o=.1):
    """
horizontal lines are drawn covering the bounding box of a sketch
example:
l1=circle(5)
l2=h_lines(l1,n=10)
fo(f'''
color("blue") p_line3d({l1},.3);
color("magenta") for(p={l2}) p_line3d(p,.1);
''')
    """
    m1=a_(sec).min(axis=0)+[-o,o]
    m2=a_(sec).max(axis=0)+[o,-o]
    x0,y0=l_(m1)
    x1,y1=l_(m2)
    b=equidistant_path([[[x0,y0],[x1,y0]],[[x0,y1],[x1,y1]]],n-1)
    return b

def h_lines_sec(sec,n=10,o=.1):
    """
horizontal lines are drawn covering the closed loop secton
example:
l1=circle(5)
l2=h_lines_sec(l1,n=10)
fo(f'''
color("blue") p_line3d({l1},.3);
color("magenta") for(p={l2}) p_line3d(p,.1);
''')
    """
    d=rot2d(.001,h_lines(sec,n,o))
    e=l_(a_(s_int1(d+seg(sec))).round(4))
    g=seg(lexicographic_sort_yx(e))
    g=l_(a_(g)[(a_([line_as_unit_vector(p) for p in g]).round(4)[:,0]>0)])
    m_p=[mid_point(p) for p in g]
    m_p1=pies1(sec,m_p)
    la,lb=array(m_p).round(5),array(m_p1).round(5)
    h=a_(g)[(lb==la[:,None]).all(2).any(1)].tolist()
    return h

def v_lines(sec,n=10,o=.1):
    """
verticle lines are drawn covering the bounding box of a sketch
example:
l1=circle(5)
l2=v_lines(l1,n=10)
fo(f'''
color("blue") p_line3d({l1},.3);
color("magenta") for(p={l2}) p_line3d(p,.1);
''')
    """
    m1=a_(sec).min(axis=0)+[o,-o]
    m2=a_(sec).max(axis=0)+[-o,o]
    x0,y0=l_(m1)
    x1,y1=l_(m2)
    b=equidistant_path([[[x0,y0],[x0,y1]],[[x1,y0],[x1,y1]]],n-1)
    return b

def v_lines_sec(sec,n=10,o=.1):
    """
verticle lines are drawn to fill the closed section
example:
l1=cr2dt([[0,0,2],[10,0,2],[0,10,2],[-10,0,2]],10)
l2=v_lines_sec(l1,n=10)
fo(f'''
color("blue") p_line3d({l1},.3);
color("magenta") for(p={l2}) p_line3d(p,.1);
''')
    """
    d=rot2d(-.001,v_lines(sec,n,o))
    e=l_(a_(s_int1(d+seg(sec))).round(4))
    g=seg(lexicographic_sort_xy(e))
    g=l_(a_(g)[(a_([line_as_unit_vector(p) for p in g]).round(4)[:,1]>0)])
    m_p=[mid_point(p) for p in g]
    m_p1=pies1(sec,m_p)
    la,lb=array(m_p).round(5),array(m_p1).round(5)
    h=a_(g)[(lb==la[:,None]).all(2).any(1)].tolist()
    return h

def corner_radius_with_turtle(sec,s=20):
    """
function to create section with corner radiuses.
s represent the number of segments at each corner radius.
control points can be given as turtle movement e.g. [[0,0],[10,0],[0,10]] means the polyline starts at [0,0] and then moves 10 units in x-direction and then moves 10 units in y-direction
example:
# cr2dt is the short name of this function. below is square with rounded corners with radius 2 mm at each corner
l1=cr2dt([[0,0,2],[10,0,2],[0,10,2],[-10,0,2]],10)
fo(f'''
color("blue") p_line3d({l1},.3);
''')
    """
    sec=pts1(sec)
    r_l=array(sec)[:,2].tolist()
    sec=array(sec)[:,:2].tolist()

    p0=[sec[-1]]+sec[:-1]
    p1=sec
    p2=sec[1:]+[sec[0]]

    p0,p1,p2=array([p0,p1,p2])
    v1,v2=p0-p1,p2-p1
    u1,u2=v1/norm(v1,axis=1).reshape(-1,1),v2/norm(v2,axis=1).reshape(-1,1)
    p_o=cwv(p1.tolist())

    theta=[(180-ang_2lineccw(p1[i],p0[i],p2[i]))/2  if p_o[i]==1 else (180-ang_2linecw(p1[i],p0[i],p2[i]))/2 for i in range(len(p0))]
    th0=[theta[-1]]+theta[:-1]
    th1=theta
    c_p1=array([p1[i] if r_l[i]==0 else  
                array(p1[i])+(u1[i]*r_l[i]/cos(d2r(theta[i])))@[[cos(d2r((180-2*theta[i])/2)),sin(d2r((180-2*theta[i])/2))],[-sin(d2r((180-2*theta[i])/2)),cos(d2r((180-2*theta[i])/2))]]
                if p_o[i]==1 else 
                array(p1[i])+(u1[i]*r_l[i]/cos(d2r(theta[i])))@[[cos(d2r((180-2*theta[i])/2)),-sin(d2r((180-2*theta[i])/2))],[sin(d2r((180-2*theta[i])/2)),cos(d2r((180-2*theta[i])/2))]]
                for i in range(len(p0))]).tolist()

    cir_1=[ 
        [p1[i].tolist()]
        if r_l[i]==0 else
        circle(r_l[i],cp=c_p1[i])
        for i in range(len(p0))]

    # radiuses
    r0=[r_l[-1]]+r_l[:-1]
    r1=r_l
    # circles
    c0=[cir_1[-1]]+cir_1[:-1]
    c1=cir_1
    # center points
    cp0=[c_p1[-1]]+c_p1[:-1]
    cp1=c_p1
    # orientations
    p_o_0=[p_o[-1]]+p_o[:-1]
    p_o_1=p_o
    a=[]
    for i in range(len(p0)):
        if r0[i]==0 and r1[i]==0:
            a.append([p1[i].tolist()])
        elif p_o_0[i]==-1 and r0[i]>0 and r1[i]==0:
            a.append([cir_p_t(c0[i],p1[i]),p1[i].tolist()])
        elif p_o_0[i]==1 and r0[i]>0 and r1[i]==0:
            a.append([p_cir_t(p1[i],c0[i]),p1[i].tolist()])
        elif p_o_1[i]==-1 and r0[i]==0 and r1[i]>0:
            a.append([p_cir_t(p0[i],c1[i])])
        elif p_o_1[i]==1 and r0[i]==0 and r1[i]>0:
            a.append([cir_p_t(c1[i],p0[i])])
        elif p_o_0[i]==1 and p_o_1[i]==1 and r0[i]>0 and r1[i]>0:
            if (r0[i]*tan(d2r(th0[i]))+r1[i]*tan(d2r(th1[i])))>norm(p1[i]-p0[i]):
                raise ValueError('radiuses more than acceptable limit')
            else:
                a.append(flip(tctpf(r0[i],r1[i],cp0[i],cp1[i])[2:]))
        elif p_o_0[i]==1 and p_o_1[i]==-1 and r0[i]>0 and r1[i]>0:
            a.append(tcct(r0[i],r1[i],cp0[i],cp1[i],1))
        elif p_o_0[i]==-1 and p_o_1[i]==1 and r0[i]>0 and r1[i]>0:
            a.append(tcct(r0[i],r1[i],cp0[i],cp1[i],-1))
        elif p_o_0[i]==-1 and p_o_1[i]==-1 and r0[i]>0 and r1[i]>0:
            if (r0[i]*tan(d2r(th0[i]))+r1[i]*tan(d2r(th1[i])))>norm(p1[i]-p0[i]):
                raise ValueError('radiuses more than acceptable limit')
            else:
                a.append(tctpf(r0[i],r1[i],cp0[i],cp1[i])[:2])

    b=[0]
    for i in range(len(p0)):
        if r_l[i]>0:
            b.append(b[-1]+2)
        else:
            b.append(b[-1]+1)

    b=array(b[1:])-1


    c=concatenate(a).tolist()
    c=c if r_l[-1]==0 else c[1:]+[c[0]]
    d=[]
    for i in range(len(p0)):
        if r_l[i]==0:
            d.append([c[b[i]]])
        else:
            d.append( arc_2p(c[b[i]-1],c[b[i]],r_l[i],p_o[i],s))
            

    d=(concatenate(d).round(8)).tolist()
    d=min_d_points(d,.0001)
    d=c3t2(rot('z.0001',d))
    return d

cr2dt=corner_radius_with_turtle

def corner_radius3d_with_turtle(pnts,s=5): # Corner radius 3d where 'pnts' are the list of points with 4th coordinate in each point is radius 'rds' and 's' is number of segments for each arc
    """
    Corner radius 3d where 'pnts' are the list of points with 4th coordinate in each point is radius and 's' is number of segments for each radial arc
example:
# cr3dt is the short name of this function
l1=cr3dt([[0,0,0,5],[0,10,0,5],[20,0,0,5],[0,-10,0,5],[0,0,10,5],[0,10,0,5],
         [-20,0,0,5],[0,-10,0,5]],20)
fo(f'''
color("blue") p_line3d({l1},.3);
''')
    """
    pnts=pts3(pnts)
    rds=[pnts[i][3] if len(pnts[i])==4 else 0 for i in range(len(pnts))]
    pnts=[pnts[i][:3] for i in range(len(pnts))]
    c=[]
    for i in range(len(pnts)):
        if i==0:
            p0=pnts[len(pnts)-1]
            p1=pnts[i]
            p2=pnts[i+1]
        elif i<len(pnts)-1:
            p0=pnts[i-1]
            p1=pnts[i]
            p2=pnts[i+1]
        else:
            p0=pnts[i-1]
            p1=pnts[i]
            p2=pnts[0]
        c.append(fillet_3p_3d(p0,p1,p2,rds[i],s)[1:])
    c=array(c).reshape(-1,3).tolist()
    return remove_extra_points(array(c).round(5))

cr3dt=corner_radius3d_with_turtle

def trim_sec_ip(sec,p0,p1,side=0,dist=.1):
    """
trim any closed loop section given 2 points (p0,p1) on the sec
dist: points can be approximately with in a defined distance "dist" from the section 
example:
c1=circle(10)
l1=point_vector([0,0],[15,10])
l2=rot2d(100,l1)
p0,p1=s_int1([l1,l2]+seg(c1))
l3=trim_sec_ip(c1,p0,p1)
fo(f'''
color("blue") p_line3d({c1},.2);
color("cyan") p_line3d({l1},.2);
color("grey") p_line3d({l2},.2);
color("green") points({[p0,p1]},.5);
color("magenta") p_line3d({l3},.3);
''')
    """
    sec=sec+[sec[0]]
    c2=lineFromPointToPointOnLine(sec,p0,p1,dist) if side==0 else lineFromPointToPointOnLine(sec,p1,p0,dist)
    c2=lineFromPointTillEnd(sec,c2[0],dist)+lineFromStartTillPoint(sec,c2[1],dist) if len(c2)==2 else c2
    return remove_duplicates(c2) if side==0 else flip(remove_duplicates(c2))

def homogenise_points(a=[],pitch=1,closed_loop=0):
    """
    function to homogenise the points of a open or closed section, 
    meaning the pitch between 2 consecutive points of the section remains same,
    it can be multi-dimensional list of array.
    here:
    "a" is the list which needs to be homogenise
    "pitch" is the pitch between 2 consecutive points
    "closed_loop" to be set to "1" if the section is closed loop or else set it to "0"
    """
    try:
        b=a_(a)
    except:
        b=[ equidistant_path(p,10) for p in a]
    if len(a_(b).shape)>2:
        if closed_loop==1:
            return [equidistant_pathc(p,pitch=pitch) for p in a]
        elif closed_loop==0:
            return [ equidistant_path(p,pitch=pitch) for p in a]
    else:
        if closed_loop==1:
            return equidistant_pathc(a,pitch=pitch)
        elif closed_loop==0:
            return equidistant_path(a,pitch=pitch)

def homogenise(a=[],pitch=1,closed_loop=0):
    """
refer function homogenise_points
In addition to the function homogenise_points, it mixes all the points together in one level.
example:
c1=circle(10,[5,0])
c2=[rot2d(i,c1) for i in linspace(0,360,5)[:-1]]
c3=homogenise(c2,pitch=1,closed_loop=1)
fo(f'''
color("blue") points({c3},.2);
''')
    """
    try:
        if len(a_(a).shape)==2:
            return homogenise_points(a,pitch,closed_loop)
            
        elif len(a)>1 and len(a_(a[0]).shape)==2:
            return l_(concatenate(homogenise_points(a,pitch,closed_loop)))
    except:
        if len(a)>1 and len(a_(a[0]).shape)==2:
            return l_(concatenate(homogenise_points(a,pitch,closed_loop)))

def concave_hull(points,n=3,engaging_angle=270):
    
    def k_nearest_points(points_list,point,n):
        """
        find n nearest points from a points_list near to point.
        points will be sorted w.r.t. distance near to far
        """
        return l_(a_(points_list)[cKDTree(points_list).query(point,n)[1]])

    def spoa(points_list,p0,p1):
        """
        sort the points_list w.r.t. angle each point makes with line p0p1.
        where p0 is the center for rotation. e.g. it is the ange between line p0p1 and p0p2
        where consider p2 as one of the point in 
        """
        return l_(a_(points_list)[a_([ang_2linecw(p0,p1,p)  for p in points_list]).argsort()[::-1]])

    def spaci(points_list,line,engaging_angle):
        """
        sort points on angle and check self intersection
        """
        x1=a_([ang_2linecw(line[-1],line[-2],p) for p in points_list]).round(2)
        x1[x1>engaging_angle]=0
        a=l_(a_(points_list)[x1.argsort()[::-1]])
        for i in range(len(a)):
            if s_int1(seg(line+[a[i]])[:-1])==[]:
                line=line+[a[i]]
                break
        return line

    def ffp(points,n):
        """
        first four points calculation
        """
        p0=s_pnt(points)
        pnts=exclude_points(points,p0)
        l1=k_nearest_points(pnts,p0,n)
        p1=spoa(l1,p0,translate_2d([0,-1],p0))[0]
        pnts=exclude_points(pnts,p1)
        l1=k_nearest_points(pnts,p1,n)
        p2=spoa(l1,p1,p0)[0]
        pnts=exclude_points(pnts,p2)
        l1=k_nearest_points(pnts,p2,n)
        p3=spoa(l1,p2,p1)[0]
        pnts=exclude_points(pnts,p3)+[p0]
        l2=[p0,p1,p2,p3]
        return [l2,pnts]

    
    l2,pnts=ffp(points,n)
    while(len(pnts)>2):
        if l2[-1]==l2[0]:
            break
        l1=k_nearest_points(pnts,l2[-1],n)
        l3=spaci(l1,l2,engaging_angle)
        if l2==l3:
            n=n+1
            l2,pnts=ffp(points,n)
        else:
            l2=l3
        pnts=exclude_points(pnts,l2[-1])
    return l2[:-1]

def fillet_line_circle(l1,c1,r=1,o=1,s=10):
    """
fillet between a line and a circle
'o' is option which can be set from 1 to 4 for fillets in 4 different direction.
's' is the number of segments of the fillet
example:
h=15
line=[[-10,h],[30,h]]
cir1=circle(10,[10,10])
r2=5
s=20
fillet1=fillet_line_circle(line,cir1,r2,1)
fillet2=fillet_line_circle(line,cir1,r2,2)
fillet3=fillet_line_circle(line,cir1,r2,3)
fillet4=fillet_line_circle(line,cir1,r2,4)
fo(f'''
color("blue",.1)p_line({line},.3);
color("violet",.2)p_line({cir1},.3);
color("cyan")p_lineo({fillet1},.3);
color("blue")p_lineo({fillet2},.3);
color("magenta")p_lineo({fillet3},.3);
color("green")p_lineo({fillet4},.3);
''')
    """
    def vcost1(l1,p0):
        """
        finds the projection of the point 'p0' on line 'l1'
        """
        v1=a_(l1[1])-a_(l1[0])
        u1=v1/norm(v1)
        v2=a_(p0)-a_(l1[0])
        d1,d2=v1@v2/norm(v1),norm(cross(c23(v1),c23(v2)))/norm(v1)
        p1=l_(a_(l1[0])+u1*d1)
        return p1
    l1=rot2d(.0001,l1)
    r1=r_arc(c1)
    r2=r
    if o==1 or o==3:
        l2=path_offset(l1,-r2)
    elif o==2 or o==4:
        l2=path_offset(l1,r2)
    cp=cp_arc(c1)
    p0=vcost1(l2,cp)
    u1=a_(line_as_unit_vector(l2))
    if o==1 or o==2:
        p1=l_(a_(p0)+u1*sqrt((r1+r2)**2-l_len([cp,p0])**2))
    elif o==3 or o==4:
        p1=l_(a_(p0)-u1*sqrt((r1+r2)**2-l_len([cp,p0])**2))
        
    u2=a_(line_as_unit_vector([cp,p1]))
    p2=l_(a_(cp)+u2*r1)
    p3=vcost1(l1,p1)
    if o==1 or o==4:
        a1=arc_2p(p2,p3,r2,-1,s=s)
    elif o==2 or o==3:
        a1=arc_2p(p2,p3,r2,1,s=s)
    return a1

def fillet_line_circle_internal(l1,c1,r=1,o=1,s=10):
    """
fillet between a line and a circle. This is inside the circle fillet
'o' is option which can be set from 1 to 4 for fillets in 4 different direction.
's' is the number of segments of the fillet
example:
h=15
line=[[-10,h],[30,h]]
cir1=circle(10,[10,10])
r2=5
s=20
fillet5=fillet_line_circle_internal(line,cir1,1,1)
fillet6=fillet_line_circle_internal(line,cir1,r2,2)
fillet7=fillet_line_circle_internal(line,cir1,1,3)
fillet8=fillet_line_circle_internal(line,cir1,r2,4)
fo(f'''
color("blue",.1)p_line({line},.3);
color("violet",.2)p_line({cir1},.3);
color("blue")p_lineo({fillet5},.3);
color("magenta")p_lineo({fillet6},.3);
color("cyan")p_lineo({fillet7},.3);
color("green")p_lineo({fillet8},.3);
''')
    """
    def vcost1(l1,p0):
            """
            finds the projection of the point 'p0' on line 'l1'
            """
            v1=a_(l1[1])-a_(l1[0])
            u1=v1/norm(v1)
            v2=a_(p0)-a_(l1[0])
            d1,d2=v1@v2/norm(v1),norm(cross(c23(v1),c23(v2)))/norm(v1)
            p1=l_(a_(l1[0])+u1*d1)
            return p1
    l1=rot2d(.0001,l1)
    r1=r_arc(c1)
    r2=r
    if o==1 or o==3:
        l2=path_offset(l1,-r2)
    elif o==2 or o==4:
        l2=path_offset(l1,r2)
    cp=cp_arc(c1)
    p0=vcost1(l2,cp)
    u1=a_(line_as_unit_vector(l2))
    if o==1 or o==2:
        p1=l_(a_(p0)+u1*sqrt((r1-r2)**2-l_len([cp,p0])**2))
    elif o==3 or o==4:
        p1=l_(a_(p0)-u1*sqrt((r1-r2)**2-l_len([cp,p0])**2))
        
    u2=a_(line_as_unit_vector([cp,p1]))
    p2=l_(a_(cp)+u2*r1)
    p3=vcost1(l1,p1)
    
    if o==1 or o==4:
        a1=arc_2p(p2,p3,r2,1,s=s)
    elif o==2 or o==3:
        a1=arc_2p(p2,p3,r2,-1,s=s)
    return a1

def extend_line(line,sec):
    """
extend a line to an intersecting section
example:
c1=circle(10)
l1=point_vector([-10,-5],[10,0])
l2=extend_line(l1,c1)
txt1=label_linear(l1,"original line",3,text_size=.5)
txt2=label_linear(l2,"extended line",-3, text_size=.5)
fo(f'''
color("blue") for(p={[c1]}) p_line3d(p,.3);
%color("cyan",.2) p_line3d({l1},.3);
color("magenta") p_line3d({l2},.1);
{txt1}{txt2}
''')
    """
    l1=line[-2:]
    v1=line_as_unit_vector(l1)
    v2=lines2vectors(seg(sec))
    v1=[v1]*len(v2)
    iim=a_([a_(v1),-a_(v2)]).transpose(1,0,2).transpose(0,2,1)
    im=inv(iim)
    p=a_(sec)-a_(l1[1])
    tx=einsum('ijk,ik->ij',im,p)
    ta,tb=tx.transpose(1,0)
    d1=(ta>0.001)&(tb>=0)&(tb<1)
    p0=l_(a_(l1[1])+a_(v1[0])[None,:]*ta[d1][:,None])
    p0=l_(a_(p0)[(cKDTree([l1[1]]).query(p0))[0].argsort()][0])
    return line[:-1]+[p0]



def line_multi_sections_ip(line,sections=[]):
    """
intersection points between line and multiple sections.
intersection points sorted w.r.t. distance from the line's starting point
example:
c1=circle(10)
c2=[translate_2d([i,0,0],c1) for i in linspace(0,60,5)]
l1=point_vector([-10,-5],[80,0])
p0=line_multi_sections_ip(l1,c2)
fo(f'''
color("blue") for(p={c2}) p_line3d(p,.3);
color("cyan") p_line3d({l1},.3);
color("magenta") points({p0},.5);
''')
    """
    l1=line[-2:]
    s1=l_(concatenate([seg(p) for p in sections]))
    v1=line_as_unit_vector(l1)
    v2=lines2vectors(s1)
    v1=[v1]*len(v2)
    iim=a_([a_(v1),-a_(v2)]).transpose(1,0,2).transpose(0,2,1)
    im=inv(iim)
    p=a_(s1)[:,0]-a_(l1[0])
    tx=einsum('ijk,ik->ij',im,p)
    ta,tb=tx.transpose(1,0)
    d1=(tb>=0)&(tb<1)
    p0=l_(a_(l1[0])+a_(v1[0])[None,:]*ta[d1][:,None])
    try:
        p0=l_(a_(p0)[(cKDTree([l1[0]]).query(p0))[0].argsort()])
    except:
        p0=[]

    return p0


def dim_radial(a1,cross_hair_size=2,text_color="blue",text_size=1,line_color="blue",arc_color="magenta",outside=0):
    """
radial dimensions with defined arc or circle 'a1'
example:
c1=circle(10)
txt1=dim_radial(c1)
fo(f'''
color("blue") for(p={[c1]}) p_line3d(p,.1);
{txt1}
''')
    """
    def point_vector(point,vector):
        """
        draw a line by defining a point and vector
        """
        p0=a_(point)
        v0=a_(vector)
        p1=p0+v0
        return l_(a_([p0,p1]))
    a1=c23(a1)
    c1=center_circle3d(a1)
    p0=mid_point(a1)
    vector=line_as_vector([c1,p0])
    l1=point_vector(c1,vector)
    chs=cross_hair_size
    tc=text_color
    ts=text_size
    lc=line_color
    ac=arc_color
    pl1=plane(line_as_axis(l1),[chs,chs],l1[0])
    l2=[mid_point(pl1[0]),mid_point(pl1[1])]
    l3=mid_line(pl1[0],pl1[1])
    l4,l5=translate(vector,[l2,l3])
    if outside==1:
        l1=mirror_line(l1,vector,[l4[0]])
        l2,l3=mirror_line(l2,vector,l4[0]),mirror_line(l3,vector,l4[0])
        
    p0=mid_point(l1)
    
    txt=f"""
    color("{ac}") p_line3d({a1},.1);
    color("{lc}")for(p={[l1,l2,l3,l4,l5]})p_line3d(p,.1);
    color("{tc}")translate({p0})linear_extrude(.2)text(str("R",{round(l_len([c1,a1[0]]),2)}),{ts});"""
    return txt

def dim_angular(l1,l2,text_color="blue",text_size=1,line_color="blue",arc_color="magenta"):
    """
angular dimension between 2 lines 'l1' and 'l2'
example:
p0,p1,p2=[[10,0],[0,0],rot2d(120,[10,0])]
txt1=dim_angular([p0,p1],[p1,p2])
fo(f'''
color("blue") for(p={[[p0,p1,p2]]}) p_line3d(p,.1);
color("magenta") for(p={[[p0,p1,p2]]}) points(p,.3);
{txt1}
''')
    """
    l1,l2=c23([l1,l2])
    tc=text_color
    ts=text_size
    lc=line_color
    ac=arc_color
    v1,v2=[line_as_unit_vector(p) for p in [l1,l2]]
    v1=l_(a_(v1)*-1)
    a1=l_(r2d(arcsin(round(norm(cross(v1,v2))/(norm(v1)*norm(v2)),4))))
    p2=l_(v1)
    p3=vcost1([p2,[0,0,0]],l_(v2))
    d1=l_len([p2,p3])
    a1=180-a1 if d1>1 else a1
    p4=movePointOnLine(flip(l1),l1[-1],(l_len(l1)+l_len(l2))/2/4)
    p5=movePointOnLine(l2,l2[0],(l_len(l1)+l_len(l2))/2/4)
    d2=l_len([p4,p5])
    a2=arc_2p_3d(nv([v1,[0,0,0],v2]),p4,p5,d2)
    p6=mid_point(a2)
    txt=f"""
    color("{lc}")for(p={[l1,l2]})p_line3d(p,.1);
    color("{ac}")p_line3d({a2},.1);
    color("{tc}")translate({p6})linear_extrude(.2)text(str({round(a1,2)},"deg"),{ts});
    """
    return txt

def vcost1(l1,p0):
    """
finds the projection of the point 'p0' on line 'l1'
example:
l1=helix(10,10,5,20)
l2=point_vector([0,0,0],[0,0,1])
p0=[ vcost1(l2,p) for p in l1]
l3=cpo([l1,p0])
fo(f'''
color("blue") points({l1},.3);
color("magenta") points({p0},.3);
color("cyan") p_line3d({l2},.3);
color("grey") for(p={l3}) p_line3d(p,.1);
''')
    """
    v1=a_(l1[1])-a_(l1[0])
    u1=v1/norm(v1)
    v2=a_(p0)-a_(l1[0])
    d1,d2=v1@v2/norm(v1),norm(cross(c23(v1),c23(v2)))/norm(v1)
    p1=l_(a_(l1[0])+u1*d1)
    return p1

def dim_linear(l1,gap=2,cross_hair_size=2,text_color="blue",text_size=1,line_color="blue"):
    """
linear dimensions with defined line l1
example:
l1=point_vector([-5,-5],[10,4])
txt1=dim_linear(l1)
fo(f'''
color("blue") for(p={[l1]}) p_line3d(p,.3);
{txt1}
''')
    """
    def point_vector(point,vector):
        """
        draw a line by defining a point and vector
        """
        p0=a_(point)
        v0=a_(vector)
        p1=p0+v0
        return l_(a_([p0,p1]))
    l1=c23(l1)
    p0=l1[0]
    vector=line_as_axis(l1)
    l1=point_vector(p0,vector)
    chs=cross_hair_size
    tc=text_color
    ts=text_size
    lc=line_color
    pl1=plane(vector,[gap,gap],l1[0])
    l2=[mid_point(pl1[0]),mid_point(pl1[1])]
    l3=mid_line(pl1[0],pl1[1])
    p0=l3[0]
    p1=translate(vector,p0)
    pl1=plane(vector,[chs,chs],p0)
    l2=[mid_point(pl1[0]),mid_point(pl1[1])]
    l3=mid_line(pl1[0],pl1[1])
    l4,l5=translate(vector,[l2,l3])
    p2=mid_point([p0,p1])
    txt=f"""
    color("{lc}")for(p={[[p0,p1],l2,l3,l4,l5]})p_line3d(p,.1);
    color("{tc}")translate({p2})linear_extrude(.2)text(str({l_len(l1)}),{ts});"""
    return txt

def point_vector(point,vector):
    """
draw a line by defining a point and vector
example:
a line [[10,0],[10,10]] can be defined as point_vector(point=[10,0],vector=[0,10])
    """
    p0=a_(point)
    v0=a_(vector)
    p1=p0+v0
    return l_(a_([p0,p1]))

def wrap_surface_around_path(surf,path):
    """
wrap a surface around a specified 3d path
example:
c1=translate_2d([0,20.1],circle(20))
s1=h_lines_sec(c1,100)
path=rot('y90',circle(40.2/(2*pi)+.2))
c2=wrap_around(c1,path)
s2=[wrap_around(p,path) for p in s1]
fo(f'''
color("blue") p_line3d({c1},.2);
color("cyan") p_line3d({path},.2);
color("magenta") p_line3d({c2},.2);

color("blue") for(p={s1}) p_line3d(p,.1,1);
color("magenta") for(p={s2}) p_line3d(p,.1,1);
''')
    """
    surf=translate([0,.001,0],surf)
    return [ wrap_around(p,path) for p in surf]

def s_int1_first(sec1):
    """
similar to 's_int1' function but only calculates intersection with the first segment with all the other segment
example:
l1=point_vector([-20,-20],[50,40])
c1=circle(20)
p0=s_int1_first([l1]+seg(c1))
fo(f'''
color("blue") for(p={[l1]}) p_line3d(p,.3);
color("cyan") p_line3d({c1},.3);
color("magenta") points({p0},.5);
''')
    """
    n=len(sec1)
    a=array(sec1)[comb_list(n)[comb_list(n)[:,0]==0]]
    p0=a[:,0][:,0]
    p1=a[:,0][:,1]
    p2=a[:,1][:,0]
    p3=a[:,1][:,1]
    v1=a_(rot2d(0.00001,p1-p0))
    v2=p3-p2
    iim=array([v1,-v2]).transpose(1,0,2).transpose(0,2,1)
    im=inv(iim)
    p=p2-p0

    t=einsum('ijk,ik->ij',im,p)
    dcn=(t[:,0].round(4)>0)&(t[:,0].round(4)<1)&(t[:,1].round(4)>0)&(t[:,1].round(4)<1)
    i_p1=p0+einsum('ij,i->ij',v1,t[:,0])
    i_p1=i_p1[dcn].tolist()
    return i_p1

def points_projection_on_line_within(line,points):
    """
finds the perpendicular projection of points on a line.
only those points will be projected which are within the bounds of the line.
example:
a=helix(10,10,5,20)
l1=point_vector([0,0,10],[0,0,30])
b=points_projection_on_line_within(l1,a)
l2=[[p,vcost_within(l1,p)] for p in a if vcost_within(l1,p)!=[]]
fo(f'''
color("blue") for(p={[a]}) points(p,.3);
color("cyan") p_line3d({l1},.1);
color("magenta") for(p={[b]}) points(p,.2);
color("grey") for(p={l2}) p_line3d(p,.1);
''')
    """
    def vcost_within(l1,p0):
        """
        finds the projection of the point 'p0' on line 'l1'.
        projects only those points which are strictly projected within the line
        """
        v1=a_(l1[1])-a_(l1[0])
        u1=v1/norm(v1)
        d=norm(v1)
        v2=a_(p0)-a_(l1[0])
        d1,d2=v1@v2/norm(v1),norm(cross(c23(v1),c23(v2)))/norm(v1)
        if d1>=0 and d1<=d :
            p1=l_(a_(l1[0])+u1*d1)
            t1=d1/d
        else:
            p1=[]
        return p1
    a=[ vcost_within(line,p) for p in points]
    b=[p for p in a if p!=[]]
    return b

def vcost_within(l1,p0):
    """
    finds the projection of the point 'p0' on line 'l1'.
    projects only those points which are strictly projected within the line
    """
    v1=a_(l1[1])-a_(l1[0])
    u1=v1/norm(v1)
    d=norm(v1)
    v2=a_(p0)-a_(l1[0])
    d1,d2=v1@v2/norm(v1),norm(cross(c23(v1),c23(v2)))/norm(v1)
    if d1>=0 and d1<=d :
        p1=l_(a_(l1[0])+u1*d1)
        t1=d1/d
    else:
        p1=[]
    return p1

def projected_points_on_line_within(line,points):
    def vcost_within_point(l1,p0):
        """
        finds the projected point, refer function vcost_within
        """
        v1=a_(l1[1])-a_(l1[0])
        u1=v1/norm(v1)
        d=norm(v1)
        v2=a_(p0)-a_(l1[0])
        d1,d2=v1@v2/norm(v1),norm(cross(c23(v1),c23(v2)))/norm(v1)
        if d1>=0 and d1<=d :
            p1=p0
            
        else:
            p1=[]
        return p1
    a=[ vcost_within_point(line,p) for p in points]
    b=[ p for p in a if p!=[]]
    return b

def line_projection_sorting(line,points):
    """
    sort the projected points in the order in which they are projected 
    from line's first point to the last point
    """
    p0=projected_points_on_line_within(line,points)
    p1=points_projection_on_line_within(line,points)
    p2=l_(a_(p0)[cKDTree(p1).query(line[0],len(p1))[1]])
    return p2

def sort_seg_on_line_length(segments):
    l1=segments
    l2=a_([l_len(p) for p in l1]).argsort()[::-1]
    l1=l_(a_(l1)[l2])
    return l1

def ang3points(p0,p1,p2):
    """
angle of 3 consecutive points p0,p1,p2
function shows the angle p0p1p2 only in 2d coordinate system
example:
p0,p1,p2=[[10,0],[0,0],rot2d(120,[10,0])]
txt1=dim_angular([p0,p1],[p1,p2])
fo(f'''
color("blue") for(p={[[p0,p1,p2]]}) p_line3d(p,.1);
color("magenta") for(p={[[p0,p1,p2]]}) points(p,.3);
{txt1}
''')
ang3points(p0,p1,p2)
    """
    if cw([p0,p1,p2])==1:
       a1= ang_2lineccw(p1,p0,p2)
    elif cw([p0,p1,p2])==-1:
        a1=ang_2linecw(p1,p0,p2)
    return a1

def fit_pline2line(polyline,line):
    """
fit a polyline between 2 defined points or a line
example:
l1=turtle3d([[0,10,0],[0,0,20]])
l2=sinewave(50,3,2,50)
l3=fit_pline2line(l2,l1)
fo(f'''
color("blue") for(p={[l1]}) p_line3d(p,.1);
color("magenta") for(p={[l2]}) p_line3d(p,.3);
color("cyan") for(p={[l3]}) p_line3d(p,.3);
''')
    """
    def pline2length(line,length):
        a=line
        d=length
        l1=l_len([a[0],a[-1]])
        rt=d/l1
        b=[[0,0,0]]+turtle3d([line_as_vector( c23(line2length(p,l_len(p)*rt))) for p in seg(a)[:-1]])
        return c32(b) if a_(a).shape[1]==2 else b
    a=c23(polyline)
    l1=line
    d=l_len(l1)
    b=pline2length(a,d)
    t1=a_(l1[0])-a_(b[0])
    b=translate(t1,b)
    u1=line_as_unit_vector(l1)
    u2=line_as_unit_vector([b[0],b[-1]])
    a1=cross(u1,u2)
    theta1=-l_(r2d(arccos(a_(u1)@a_(u2))))
    b=axis_rot_1(b,a1,b[0],theta1)
    return b

def surface_from_4_lines(l1,l2,l3,l4):
    """
create a surface by 4 lines enclosure
l1 and l2 are opposite side
l3 and l4 are opposite side
l1 and l3 starts from the same point
l2 and l4 ends at the same point
example:
l1=rot('x90',arc_2p([0,0],[100,0],100,s=100))
l2=translate([0,100,0],rot('z30',l1))
l3=rot('x90',sinewave(100,4,2,100))
l3=fit_pline2line(l3,[l1[0],l2[0]])
l4=rot('x90',cosinewave(100,4,2,100))
l4=fit_pline2line(l4,[l1[-1],l2[-1]])
surf1=surface_from_4_lines(l1,l2,l3,l4)
txt1=label_linear([l1[0],l1[-1]],"L1",10,text_size=3)
txt2=label_linear([l2[0],l2[-1]],"L2",-10,text_size=3)
txt3=label_linear([l3[0],l3[-1]],"L3",-10,text_size=3)
txt4=label_linear([l4[0],l4[-1]],"L4",10,text_size=3)
fo(f'''
color("blue") for(p={[l1,l2,l3,l4]}) p_line3d(p,.3);
{txt1}{txt2}{txt3}{txt4}
%{swp_surf(surf1)}
''')
    """
    l2=path2path1(l1,l2)
    l4=path2path1(l3,l4)
    s1=slice_sol([l1,l2],len(l3)-1)
    s2=cpo([l3,l4])
    s3=[ fit_pline2line(s1[i],s2[i]) for i in range(len(s1))]
    return s3

def s_int1_3d(sec1):
    """
    same as function s_int1, applicable to 3d coordinates
    """
    n=len(sec1)
    a=array(sec1)[comb_list(n)]
    p0=a[:,0][:,0]
    p1=a[:,0][:,1]
    p2=a[:,1][:,0]
    p3=a[:,1][:,1]
    v1=a_(rot('z.00001',p1-p0))
    v2=a_(rot('x.00001',p3-p2))
    iim=a_([v1[:,:2],-v2[:,:2]]).transpose(1,0,2).transpose(0,2,1)
    im=inv(iim)
    p=(p2-p0)[:,:2]
    t=einsum('ijk,ik->ij',im,p)
    
    iim1=a_([v1[:,1:],-v2[:,1:]]).transpose(1,0,2).transpose(0,2,1)
    im1=inv(iim1)
    px=(p2-p0)[:,1:]
    t1=einsum('ijk,ik->ij',im1,px)
    
    dcn=(t[:,0].round(4)>0)&(t[:,0].round(4)<1)&(t[:,1].round(4)>0)&(t[:,1].round(4)<1) 
    dcn1=(t1[:,0].round(4)>0)&(t1[:,0].round(4)<1)&(t1[:,1].round(4)>0)&(t1[:,1].round(4)<1)
    dcn2=dcn&dcn1
    i_p1=p0+einsum('ij,i->ij',v1,t[:,0])
    i_p1=i_p1[dcn2].tolist()
    return i_p1

def sec_start_pos(sec,n=0):
    """
    function to change the starting position/ point of a closed loop section.
    'n' can be any number between 0 to length of sec
    """
    return sec[n:]+sec[:n]

def lines_fillets_solid(l1,l2,l3,l4,l5,s=20,o=1):
    """
convert 5 lines to solid for removing sharp edges from a shape
example:
sol=linear_extrude(circle(20),10)
l1=offset_3d(sol[0],-5)
l2=sol[0]
l3=mid_line(sol[0],sol[-1])
l4=sol[-1]
l5=offset_3d(sol[-1],-5)
f1=lines_fillets_solid(l1,l2,l3,l4,l5,s=20,o=5)
f1=f1+[f1[0]]
txt1=label_radial(l1,"L1",text_size=2)
txt2=label_radial(l2[10:],"L2",text_size=2)
txt3=label_radial(l3[20:],"L3",text_size=2)
txt4=label_radial(l4[30:],"L4",text_size=2)
txt5=label_radial(l5[40:],"L5",text_size=2)
fo(f'''
color("blue") for(p={[l1,l2,l3,l4,l5]}) p_line3dc(p,.3);
{txt1}{txt2}{txt3}{txt4}{txt5}
difference(){{
//{swp(sol)}
#{swp_c(f1)}
}}
''')
    """
    f1=cpo(cpo(convert_3lines2fillet(l1,l3,l2,s))[:-1])
    f2=cpo(cpo(convert_3lines2fillet(l3,l5,l4,s))[:-1])
    f3=[ f1[i]+f2[i] for i in range(len(f1))]
    f4=surface_offset(f3,o)
    sol1=solid_from_2surfaces(f3,f4)
    return sol1

def fo(txt='',fname='trial.scad'):
   with open(fname,'w+') as f:
       f.write(txt)
       f.write(f'''
       //module for drawing a closed 2d polyline from a group of points "path" and width of the polyline is defined by parameter "size".
module p_line(path,size=.5){{
    for(i=[0:len(path)-1])
        let(p0=path[i],p1=i<len(path)-1?path[i+1]:path[0])
    
    hull(){{
    translate(p0)circle(size/2,$fn=20);
    translate(p1)circle(size/2,$fn=20);}}
}}
    
//module for drawing an open 2d polyline from a group of points "path" and width of the polyline is defined by parameter "size".
module p_lineo(path,size=.5){{
    for(i=[0:len(path)-2])
        let(p0=path[i],p1=path[i+1])
    
    hull(){{
    translate(p0)circle(size/2,$fn=20);
    translate(p1)circle(size/2,$fn=20);}}
          }}
    
    
// module for rendering points along the various shapes 2d or 3d. parameter "d" is the size of cube which is used as point. a list has to be provided for parameter "p"
// try following code:
// sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5);
// prism=l_extrude(sec,h=15,a=90,steps=20);
// %swp(prism);
// for(p=prism) points(p,.2);
 
module points(p,d=.5){{
    for(i=p)translate(i)cube(size=d,center=true);
    
    }}
    
//module to draw a polyline in 3d space (loop not closed)

    
module p_line3d(path,d,rec=0,$fn=20){{
    for(i=[0:len(path)-2])
        
    hull(){{
    translate(path[i])if(rec==0)sphere(d/2); else cube(d,true);
    translate(path[i+1])if(rec==0)sphere(d/2);else cube(d,true);
    }}
}}

//module to draw a polyline in 3d space (loop closed)
  

module p_line3dc(path,d,rec=0,$fn=20){{
    for(i=[0:len(path)-1])
        let(
    i_plus=i<len(path)-1?i+1:0
    )
    hull(){{
    translate(path[i])if(rec==0)sphere(d/2); else cube(d,true);
    translate(path[i_plus])if(rec==0)sphere(d/2);else cube(d,true);
    }}
}}
    
function faces(sol)=

//    calculate the faces for the vertices with shape l x m with first and the last end closed
    let(
    l=len(sol),
    m=len(sol[0]),
    n1=[for(i=[0:m-1])i],
    n2=[for(i=[0:l-2]) each ([ for(j=[0:m-1])
    each
    j<m-1?[[(j+1)+i*m,j+i*m,j+(i+1)*m],[(j+1)+i*m,j+(i+1)*m,(j+1)+(i+1)*m]]:
    [[0+i*m,j+i*m,j+(i+1)*m],[0+i*m,j+(i+1)*m,0+(i+1)*m]]
    ])],
    n3=[for(i=[0:m-1])i+(l-1)*m],
    n4=[for(i=[len(n3)-1:-1:0])n3[i]],
    n=[n1,each (n2),n4]
    )n;
    

function faces_1(sol)=

//    calculate the faces for the vertices with shape l x m with first and the last end open
    let(
    l=len(sol),
    m=len(sol[0]),
    n2=[for(i=[0:l-2])each([ for(j=[0:m-1])
    each
    j<m-1?[[(j+1)+i*m,j+i*m,j+(i+1)*m],[(j+1)+i*m,j+(i+1)*m,(j+1)+(i+1)*m]]:
    [[0+i*m,j+i*m,j+(i+1)*m],[0+i*m,j+(i+1)*m,0+(i+1)*m]]
    ])]
    
    )n2; 
   
 function vertices(sol)=
[each for (p=sol)p];

// module for rendering the polyhedron with ends closed
module swp(sol){{
let(
v1=vertices(sol),
f1=faces(sol)
)
polyhedron(v1,f1,convexity=10);

}}

// module for rendering polyhedron with ends open (mainly for closed polyhedron)
module swp_c(sol){{
let(
v1=vertices(sol),
f1=faces_1(sol)
)
polyhedron(v1,f1,convexity=10);

}}
// function to select specific section of points from a list of points
function loop ( list , a , b ) = [ for (i = [ a : b ] ) list[i] ];

function faces_surf(sol)=

//    calculate the faces for the vertices with shape l x m with first and the last end open and considering the ends are not closed. like 2 straight lines surface
    let(
    l=len(sol),
    m=len(sol[0]),
    n2=[for(i=[0:l-2])each [for(j=[0:m-2]) each
[[m*i+j,m*(i+1)+j,m*i+(j+1)],
 [m*i+(j+1),m*(i+1)+j,m*(i+1)+(j+1)]] ]]
    
    )n2;
    
    // module for rendering polyhedron with ends open (mainly for open surfaces like 2 straight lines)
module swp_surf(sol){{
let(
v1=vertices(sol),
f1=faces_surf(sol)
)
polyhedron(v1,f1,convexity=10);

}}
// function to create slices in solid

function slice_sol(sol_1,n=10)=
cpo([for(p=cpo(sol_1)) m_points_so(p,n)]);


//function to change the orientation of points of a prism. for example check prism and prism1
//sec=cr(pts1([[5,20],[20,10,50],[20,-7]]),20);
//sec1=trns([0,-.05,0],sec);
//prism=[for(i=[0:5:355])rot([1,0,0],sec,i)];
//prism1=cpo(prism);
//
//translate([-60,0,0])
//for(p=prism)points(p,.5);
//for(p=prism)p_line3d(p,.2);
//translate([60,0,0])
//for(p=prism1)p_line3d(p,.2);
//
//translate([-60,0,30])rotate([90,0,0])text("points");
//translate([0,0,30])rotate([90,0,0])text("prism");
//translate([60,0,30])rotate([90,0,0])text("prism1");

function cpo(prism)=[for(i=[0:len(prism[0])-1])[for(p=prism)p[i]]];

// function for calculating multiple points on the straight line segments of an open section. sec-> closed section; s -> number of segments for each straight line segment of closed section; m-> minimum segment length, if the derived segment length < m, then it is omitted
// sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5);
// points(sec,.2);
// 
// translate([15,0])
// points(m_points_so(sec,s=5,m=.5),.2);// number of segments=> 5

function m_points_so(sec1,s,m=.5)=
let(
l=[for(i=[0:len(sec1)-2])
let(
i_plus=i+1,
l=norm(sec1[i_plus]-sec1[i]),
u=uv(sec1[i_plus]-sec1[i])
)each l/s>=m?[for(j=[0:l/s:l])sec1[i]+j*u]:[sec1[i]]]
)l;

// function to calculate a unit vector for a given vector
// example: 
// echo(uv([2,3,4])); // => ECHO: [0.371391, 0.557086, 0.742781]

function uv(v)=v/norm(v);

    // module for rendering polyhedron with a closed polygon
module swp_sec(sec){{
let(
v1=sec,
f1=[for(i=[0:len(sec)-1]) i]
)
polyhedron(v1,[f1],convexity=10);

}}

       ''')

def i_p2d(l1,l2):
    """
    function to calculate the intersection point between 2 lines in 2d space
    example:
    i_p2d(l1=[[0,0],[1,4]],l2=[[10,0],[7,2]]) =>  [1.42857, 5.71429]
    """
    l1,l2=array([l1,l2])
    v1=l1[1]-l1[0]
    v2=l2[1]-l2[0]
    t1,t2=inv(array([v1,-v2]).transpose(1,0))@(l2[0]-l1[0])
    return (l1[0]+v1*t1).tolist()

def uv(v):
    """
    function to calculate unit vector of a given vector
    example:
    vector=[2,3,5]
    unit_vector=uv(vector) => [0.3244428422615251, 0.48666426339228763, 0.8111071056538127]
    """
    v=array(v)
    return (v/norm(v)).tolist()

unit_vector=uv

def tcct(r1,r2,cp1,cp2,cw=-1): # two circle cross tangent
    """
    function to draw cross tangent between 2 circles
    refer to the file "example of various functions " for application examples
    """
    v1=[1,1]
    v2=[-r2,r1]
    cp1,cp2=array([cp1,cp2])
    d=norm(cp2-cp1)
    d1=(inv(array([v1,v2]).T)@array([d,0]))[0]
    d2=(inv(array([v1,v2]).T)@array([d,0]))[1]
    a=arcsin(r1/d1)*180/pi
    v3=cp2-cp1
    u3=v3/norm(v3)
    b=arccos(u3@array([1,0]))*180/pi
    if cw==-1:
        if v3[0]>0 and v3[1]<=0:
            theta1=270+a-b
            theta2=90+a-b
        elif v3[0]>=0 and v3[1]>0:
            theta1=270+a+b
            theta2=90+a+b
        elif v3[0]<0 and v3[1]>=0:
            theta1=270+a+b
            theta2=90+a+b
        else:
            theta1=270+a-b
            theta2=90+a-b
    else:
        if v3[0]>0 and v3[1]<=0:
            theta2=270-a-b
            theta1=90-a-b
        elif v3[0]>=0 and v3[1]>0:
            theta2=270-a+b
            theta1=90-a+b
        elif v3[0]<0 and v3[1]>=0:
            theta2=270-a+b
            theta1=90-a+b
        else:
            theta2=270-a-b
            theta1=90-a-b
        
    p0=(cp1+array([r1*cos(theta1*pi/180),r1*sin(theta1*pi/180)])).tolist()
    p1=(cp2+array([r2*cos(theta2*pi/180),r2*sin(theta2*pi/180)])).tolist()
    return [p0,p1]

def tctpf(r1,r2,cp1,cp2): #2 circle tangent point full (both the sides)
    """
    function to draw tangent line joining 2 circles with radiuses "r1" and "r2" with center points "cp1" and "cp2" respectively. 
    This function draws tangent line on both the sides
    example:
    cir1=circle(10)
    cir2=circle(5,[15,6])
    sec=tctpf(r1=10,r2=5,cp1=[0,0],cp2=[15,6])
    
    
    """
    cp1,cp2=array([cp1,cp2])
    v1=cp2-cp1,
    u1=v1/norm(v1)
    ang1=arcsin((r2-r1)/norm(cp2-cp1))*180/pi

    t1=cp1+u1@rm(90+ang1)*r1
    t2=cp2+u1@rm(90+ang1)*r2

    t3=cp1+u1@rm(-90-ang1)*r1
    t4=cp2+u1@rm(-90-ang1)*r2

    return [t3[0].tolist(),t4[0].tolist(),t2[0].tolist(),t1[0].tolist()]

def reorient_sec(sec):
    """
    re-orient section to create a better surface through 'prism' function
    """
    sec1=sec2surface_1(sec)
    sec2=[p[0] for p in sec1]+flip([p[-1] for p in sec1])
    return sec2

def i_p_t(path):
    """
    function to calculate tangent vectors to a given path
    """
    p1=array(seg(path))
    p2=array(path)
    v1=array([(p[1]-p[0])/norm(p[1]-p[0]) for p in p1])
    t_v=array([ (v1[-1]+v1[i])/2 if i==0 else
         (v1[i-1]+v1[i])/2
        for i in range(len(p1))])

    t_v=t_v/norm(t_v,axis=1).reshape(-1,1)

    return t_v

def i_p_t_o(path):
    """
    function to calculate tangent vectors to a given open path
    """
    p1=array(seg(path)[:-1])
    p2=array(path)
    v1=array([(p[1]-p[0])/norm(p[1]-p[0]) for p in p1])
    t_v=[(v1[i-1]+v1[i])/2
        for i in range(1,len(p1))]
    t_v=a_([t_v[0]]+t_v+[t_v[-1]])

    t_v=t_v/norm(t_v,axis=1).reshape(-1,1)

    return t_v

def o_p_p(sol,i_p,d):
    """
    calculates projected points on the surface of a solid 
    sol: solid on which the points to be projected
    i_p: list of points in 3d space near the solid
    d: approximate distance of the points from the surface, specifying too big distance 
    may create multiple projection of the same point on the solid
    """
    l,m,_=array(sol).shape
    f1=faces(l,m)[1:-1]
    v=array(sol).reshape(-1,3)
    tri=v[f1]
    v2,v3=tri[:,1]-tri[:,0], tri[:,2]-tri[:,0]
    v1=cross(v2,v3)
    v1=v1/norm(v1,axis=1).reshape(-1,1)
    p0=array(i_p)
    p2=tri[:,0]
    n1,n2=len(p0),len(p2)
    v1=array([v1]*n1)
    v2=array([v2]*n1)
    v3=array([v3]*n1)
    p0=array([p0]*n2).transpose(1,0,2)
    p2=array([p2]*n1)
    iim=array([v1,-v2,-v3]).transpose(1,2,0,3).transpose(0,1,3,2)
    im=inv(iim)
    t=einsum('ijkl,ijl->ijk',im,p2-p0)
    d1=(t[:,:,0]>-d)&(t[:,:,0]<d)&(t[:,:,1]>=0)&(t[:,:,1]<=1) \
    &(t[:,:,2]>=0)&(t[:,:,2]<=1)&(t[:,:,1]+t[:,:,2]<1)

    p0.shape,v1.shape,t[:,:,0].shape
    nx=p0+einsum('ijk,ij->ijk',v1,t[:,:,0])
    nx=nx[d1].tolist()
    return nx

def int_seg_list(sec1):
    n=len(sec1)
    a=array(sec1)[comb_list(n)]
    p0=a[:,0][:,0]
    p1=a[:,0][:,1]
    p2=a[:,1][:,0]
    p3=a[:,1][:,1]
    v1=p1-p0
    v2=p3-p2
    iim=array([v1,-v2+.00001]).transpose(1,0,2).transpose(0,2,1)
    im=inv(iim)
    p=p2-p0

    t=einsum('ijk,ik->ij',im,p)
    dcn=(t[:,0].round(4)>0)&(t[:,0].round(4)<1)&(t[:,1].round(4)>0)&(t[:,1].round(4)<1)

    d=comb_list(n)[dcn]
    return d

def surface_for_fillet(sol1=[],sol2=[],factor1=50,factor2=20,factor3=4,factor4=25,dia=8):
    """
    sol1: Solid on which the surface needs to be created
    sol2: Intersecting solid
    factor1: number of segments in the circle
    factor2: number of layers or slices of surface
    factor3: decides the size of the surface lower value means bigger size. value can be set between 1 to any number
    factor4: any high number should be ok like maybe 100 or greater, basically greater than the bounding box dimension of the "sol1"
    dia: diameter around the solid 2 where surfavce needs to be created
    """
    v,f1=partial_surface(sol1,prism_center(sol2),dia)
    tri=array(v)[array(f1)]
    s1=shield(sol1,sol2,factor1,factor2,factor3,factor4)
    p0,p1,p2=tri[:,0],tri[:,1],tri[:,2]
    s2=array(s1).transpose(0,2,1,3)
    la,lb=s2[:,:,0],s2[:,:,1]
    p01,p02,lab=p1-p0,p2-p0,lb-la

    p0.shape,array(s1).shape,s2.shape,la.shape,p01.shape,lab.shape
    t=einsum('ijklm,ijklm->ijkl',cross(p01,p02)[None,None,None,:,:],(la[:,:,None,None,:]-p0[None,None,None,:,:]))/(einsum('ijklm,ijklm->ijkl',(-lab)[:,:,None,None,:],cross(p01,p02)[None,None,None,:,:])+.00001)
    u=einsum('ijklm,ijklm->ijkl',cross(p02[None,None,None,:,:],(-lab)[:,:,None,None,:]),(la[:,:,None,None,:]-p0[None,None,None,:,:]))/(einsum('ijklm,ijklm->ijkl',(-lab)[:,:,None,None,:],cross(p01,p02)[None,None,None,:,:])+.00001)
    v=einsum('ijklm,ijklm->ijkl',cross((-lab)[:,:,None,None,:],p01[None,None,None,:,:]),(la[:,:,None,None,:]-p0[None,None,None,:,:]))/(einsum('ijklm,ijklm->ijkl',(-lab)[:,:,None,None,:],cross(p01,p02)[None,None,None,:,:])+.00001)
    condition=(t>=0)&(t<=1)&(u>=0)&(u<=1)&(v>=0)&(v<=1)&(u+v<=1)
    i_p=la[:,:,None,None,:]+einsum('ijklm,ijkl->ijklm',lab[:,:,None,None,:],t)
    i_p1=[]
    for i in range(1,len(i_p)):
        i_p1.append(i_p[i][condition[i]].tolist())
    return align_sol_1([equidistant_pathc(p,factor1) for p in i_p1])

def shield(sol1=[],sol2=[],factor1=50,factor2=10,factor3=1,factor4=100):
    """
    function to check function surface_for_fillet.
    sol1: Solid on which the surface needs to be created
    sol2: Intersecting solid
    factor1: number of segments in the circle
    factor2: number of layers or slices of surface
    factor3: decides the size of the surface lower value means bigger size. value can be set between 1 to any number
    factor4: any high number should be ok like maybe 100 or greater, basically greater than the bounding box dimension of the "sol1"
    """
    p0= array(prism_center(sol1))
    p1= array(prism_center(sol2))

    v1=p1-p0
    u1=v1/norm(v1)
    p3=p0-u1*factor4
    p3=ip_sol2line(sol1,[p0,p3])[0]
    cir1=circle(1,s=factor1)
    sur1=sol2vector(v1,cpo([ls([[0,0],p],factor2) for p in cir1]),u1*factor3)
    lines1=[[[[0,0,0],(array(p1)/norm(p1)*factor4).tolist()] for p1 in p] for p in sur1]
    sol3=[translate(p3,cpo(p)) for p in lines1]
    return sol3

def partial_surface(sol,cp1,dia):
    """
    return the part of the surface of a given solid 'sol' in terms of vertices and faces
    cp1: center point around which the surface is needed
    dia: dia of the sphere inside which the surface needs tobe extracted
    """
    v,f1=vnf2(sol)
    v,f1=array(v),array(f1)
    bc1=v[f1].mean(1)
    f2=f1[(sqrt((bc1[:,0]-cp1[0])**2+(bc1[:,1]-cp1[1])**2+(bc1[:,2]-cp1[2])**2)<=dia)]
    return [v.tolist(),f2.tolist()]

def i_p2dv(p0,p1,p2,p3):
    v1=p1-p0
    v2=p3-p2
    a=pinv(swapaxes(transpose(array([v1,-v2])),0,1))
    b=p2-p0
    t=einsum('ijk,ik->ij',a,b)[:,0]
    return p0+einsum('ij,i->ij',v1,t)

def extend_line3d(line,sec):
    """
    same as function extend_line, but for 3d
    """
    l1=sec
    l2=line
    l3=l1+l2
    n1,incpt=best_fit_plane(l3)
    n1=uv(n1)
    n2=[0,0,-1]
    n3=cross(n1,n2)
    theta=r2d(arccos(a_(n1)@a_(n2)))
    l4=c32(rot_sec2xy_plane(l3))
    l5=l4[:len(l1)]
    l6=l4[len(l1):]
    l7=extend_line(l6,l5)
    l8=c23(l5+l7)
    d1=a_(l1[0])-a_(l8[0])
    l8=translate(d1,axis_rot_1(l8,n3,l8[0],-theta))
    l9=l8[:len(l1)]
    l10=l8[len(l1):]
    return l10

def vcost2(line,point):
    """
    perpendicular projection of a point on a polyline
    projection will be returned if it lies within the line
    """
    for p in seg(line)[:-1]:
        px=vcost1(p,point)
        d1=l_len([p[0],px])
        d2=l_len(p)
        v1=line_as_unit_vector(p)
        v2=line_as_unit_vector([p[0],px])
        if (d1<=d2) & (l_((a_(v1).round(4)==a_(v2).round(4)).all())):
            p5=px
    return p5

# def path_offset3d(line,d=1,nv=[]):
#     """
#     same as function path_offset, but for 3d
#     in case normal vector is known use parameter 'nv' to define the same
#     """
#     l1=line
    
#     n1,incpt=best_fit_plane(l1)
#     n1=uv(n1) if nv==[] else uv(nv)
#     n2=[0,0,-1]
#     n3=cross(n1,n2)
#     theta=r2d(arccos(a_(n1)@a_(n2)))
#     l2=c32(rot_sec2xy_plane(l1))
#     l3=path_offset(l2,d)
#     l4=c23(l3)
#     d1=a_(l1[0])-a_(l4[0])
#     l5=translate(d1,axis_rot_1(l4,n3,l4[0],-theta))
#     return l5

def path_offset3d(path,d=1,nv=[]):
    """
function to offset a 3d path. path has to be in the same plane
Incase the normal vector 'nv' is known can be specified
example:
a=cr3dt([[0,0,0],[10,0,0],[0,10,3],[-5,7,5],[-15,-1,-10]])
pl1=plane([1,1,20],[100,100],[0,0,10])
b=plos(pl1,a,[1,1,20],unidirection=0)
c=path_offset3d(b,1)
fo(f'''
%{swp_surf(pl1)}
//color("blue") p_line3d({a},.3);
color("magenta") p_line3d({b},.3);
color("cyan") p_line3d({c},.3);
''')    
    """
    b=path
    n1=uv(best_fit_plane(b)[0]) if nv==[] else nv
    n2=[0,0,-1]
    n3=uv(cross(n1,n2))
    theta=r2d(arccos(a_(n1)@a_(n2)))
    b1=rot_sec2xy_plane(b)
    b2=c32(b1)
    tr1=a_(b)[0]-a_(c23(b2))[0]
    b3=path_offset(b2,d)
    b4=axis_rot_1(c23(b3),n3,c23(b2)[0],-theta)
    b5=translate(tr1,b4)
    return b5


def ilo(il,s1,s2,r,t=1):
    """
intersection line offset
type "t" can be set to 1,2 or 3
il: intersection line between 2 surfaces s1 and s2
s1: surface which is intersected
s2: surface which intersects the surface s1
r: distance to offset the intersecting line
example:
s1=sphere(20)
c1=rot('y45',cylinder(r=5,h=50))
l1=ip_sol2sol(s1,c1)
l2=ilo(l1,s1,c1,-3,1)
fo(f'''
%{swp(s1)}
%{swp(c1)}
color("blue") p_line3dc({l1},.3);
color("magenta") p_line3dc({l2},.3);
''')
    """
    a=i_p_p(s2,il,.5)
    v1=i_p_n(il,s1)
    v2=i_p_t(il)
    v3=il+cross(v2,v1)*r
    b=surface_offset([il,a],r)[0]
    v4=line_as_unit_vector([il[0],a[0]])
    v5=concatenate(s2).mean(0)
    if t==1:
        c=plos(c_(s1),b,v4,0)
    elif t==2:
        c=plos_v_2(c_(s1),v3,v1,unidirection=0)
    elif t==3:
        c=plos_v(c_(s1),v3,v5,unidirection=0)
    return c

def c_(sol):
    """
close loop a solid. In the example below the ends of the surface are closed by applying the function
example:
c1=cylinder(r=10,h=50,s=10)
c2=translate([40,0,0],c1)
c2=c_(c2)
fo(f'''
{swp_surf(c1)}
{swp_surf(c2)}
''')
    """
    return cpo(cpo(sol)+[cpo(sol)[0]])

def union(a=[],n=3,pitch=.5,closed_loop=1):
    """
creating a union of various individual elements
example:
c1=circle(20,[15,0])
c1,c2,c3,c4=[ rot2d(i,c1) for i in linspace(0,360,5)[:-1]]
c4=union([c1,c2,c3,c4],5,1)
fo(f'''
color("blue",.2) for(p={[c1,c2,c3,c4]}) p_line3dc(p,.3);
color("magenta") p_line3dc({c4},.32);
''')
    """
    sx=homogenise(a,pitch=pitch,closed_loop=closed_loop)
    return concave_hull(sx,n)

def fillet_2spheres(sp1,sp2,r,s1=20,s2=50):
    """
create a fillet between 2 spheres
r is the radius of the fillet
s1 and s2 are the number of sections in the fillet 
and number of points in each section
example:
s1=sphere(10)
s2=sphere(7,[15,15,0])
f1=fillet_2spheres(s1,s2,7,s1=10,s2=40)
fo(f'''
{swp(s1)}
{swp(s2)}
{swp(f1)}
''')
    """
    cp1=cp_cir_3d( cpo(sp1)[0])
    cp2=cp_cir_3d( cpo(sp2)[0])
    l1=[cp1,cp2]
    u1=line_as_unit_vector(l1)
    u2=[0,0,-1]
    u3=cross(u1,u2)
    u4=cross(u1,u3)
    c1=translate(cp1,sec2vector(u4,circle(.5)))
    c2=translate(cp2,sec2vector(u4,circle(.5)))
    c1=plos_v_1(sp1,c1,point_vector(cp1,u4))
    c2=plos_v_1(sp2,c2,point_vector(cp2,u4))
    
    a1=two_cir_tarc3d(c1,c2,r,s=s1)
    f1=cpo([ axis_rot_1(a1,u1,cp1,i)  for i in linspace(0,360,s2+1)[:-1]])
    return f1

def two_cir_tarc3d(c1,c2,r,side=0,s=50):
    """
    it is similar to function two_cir_tarc but in 3d coordinates.
    both the circles should lie in the same plane
    """
    n1,intcpt=best_fit_plane(c1)
    n1=uv(n1)
    u1=[0,0,-1]
    u2=cross(n1,u1)
    theta=r2d(arccos(a_(n1)@a_(u1)))
    c3=c1+c2
    c4=c32(rot_sec2xy_plane(c3))
    c5,c6=c4[:len(c1)],c4[len(c1):]
    a1=two_cir_tarc(c5,c6,r,side,s)
    tr1=a_(c3[0])-a_(c23(c4)[0])
    a2=translate(tr1,a1)
    a3=axis_rot_1(a2,u2,c3[0],-theta)
    return a3

def bezier_c(sec,s=30):
    """
    approximates a bezier closed loop
example:
a=cr3dt([[0,0,0],[10,0,0],[0,10,3],[-5,7,5],[-15,-1,-10]])
a=m_points1(a,3)
b=bezier_c(a,100)
fo(f'''
// control points
color("blue") for(p={[a]}) points(p,.3);
color("blue",.2) p_line3dc({a},.2);
// bezier closed
color("magenta") for(p={[b]}) p_line3d(p,.2);
''')
    """
    a=sec
    b=[mid_point(bezier(sec_start_pos(a,i),s)) for i in range(len(a))]
    c=equidistant_pathc(bspline_closed(b,3,s),s)
    return c

def rationalise_sec(sec):
    """
    removes all the points which are in straight line in a closed loop section or polyline
    """
    a=rationalise_path(sec)
    if pol(a[-1],[a[-2],a[0]]):
        return a[:-1]
    else:
        return a

def label_linear(l1,s1,gap=2,text_color="blue",text_size=1,line_color="blue"):
    """
    label linear elements with defined line l1
    """
    def point_vector(point,vector):
        """
        draw a line by defining a point and vector
        """
        p0=a_(point)
        v0=a_(vector)
        p1=p0+v0
        return l_(a_([p0,p1]))
    l1=c23(l1)
    p0=l1[0]
    vector=line_as_axis(l1)
    l1=point_vector(p0,vector)
    # chs=cross_hair_size
    tc=text_color
    ts=text_size
    lc=line_color
    pl1=plane(vector,[gap,gap],l1[0])
    l2=[mid_point(pl1[0]),mid_point(pl1[1])]
    l3=mid_line(pl1[0],pl1[1])
    p0=l3[0]
    p1=translate(vector,p0)
    # pl1=plane(vector,[chs,chs],p0)
    # l2=[mid_point(pl1[0]),mid_point(pl1[1])]
    # l3=mid_line(pl1[0],pl1[1])
    # l4,l5=translate(vector,[l2,l3])
    p2=mid_point([p0,p1])
    txt=f"""
    color("{lc}")for(p={[[p0,p1]]})p_line3d(p,.1);
    color("{tc}")translate({p2})linear_extrude(.2)text("{s1}",{ts});"""
    return txt


def label_radial(a1,s1,text_color="blue",text_size=1,cross_hair_size=1,line_color="blue",arc_color="magenta",outside=0):
    """
    label radial element with defined arc or circle 'a1'
    """
    def point_vector(point,vector):
        """
        draw a line by defining a point and vector
        """
        p0=a_(point)
        v0=a_(vector)
        p1=p0+v0
        return l_(a_([p0,p1]))
    a1=c23(a1)
    c1=center_circle3d(a1)
    p0=mid_point(a1)
    vector=line_as_vector([c1,p0])
    l1=point_vector(c1,vector)
    chs=cross_hair_size
    tc=text_color
    ts=text_size
    lc=line_color
    ac=arc_color
    pl1=plane(line_as_axis(l1),[chs,chs],l1[0])
    l2=[mid_point(pl1[0]),mid_point(pl1[1])]
    l3=mid_line(pl1[0],pl1[1])
    l4,l5=translate(vector,[l2,l3])
    if outside==1:
        l1=mirror_line(l1,vector,[l4[0]])
        l2,l3=mirror_line(l2,vector,l4[0]),mirror_line(l3,vector,l4[0])
        
    p0=mid_point(l1)
    
    txt=f"""
    color("{ac}") p_line3d({a1},.1);
    color("{lc}")for(p={[l1]})p_line3d(p,.1);
    color("{tc}")translate({p0})linear_extrude(.2)text("{s1}",{ts});"""
    return txt

def vector2line(v1,length=1,intercept=[0,0,0]):
    return translate(intercept,point_vector([0,0,0],a_(c23(v1))*length))

def two_cir_tarc_internal(c1,c2,r,s=20):
    """
function to draw an internal arc between 2 intersecting circles or arcs
example:
c1=circle(10)
c2=circle(8,[13,-4])
a1=two_cir_tarc_internal(c1,c2,r=2,s=20)
a2=two_cir_tarc_internal(c2,c1,r=1,s=20)
fo(f'''
color("blue") for(p={[c1,c2]}) p_line3dc(p,.1);
color("magenta") p_line3d({a1},.1);
color("cyan") p_line3d({a2},.1);
''')
    """
    # a,b,c,d=sp.symbols('a b c d') # where a=r1,b=r2,c=l
    # e1=a**2-d**2-b**2+(c-d)**2
    # sp.solve(e1,d) # to find equation for d
    r1,r2,p0,p1=r_arc(c1),r_arc(c2),cp_arc(c1),cp_arc(c2)
    l1=[p0,p1]
    l=l_len(l1)
    d=((r1-r)**2-(r2-r)**2+l**2)/(2*l)
    h=sqrt((r1-r)**2-d**2)
    theta1=r2d(arctan(h/d))
    theta2=r2d(arctan(h/(l-d)))
    u1=a_(line_as_unit_vector(l1))
    p2=c32(axis_rot_1([p0+sign(theta1)*u1*r1],[0,0,1],p0,theta1)[0])
    p3=c32(axis_rot_1([p1-sign(theta2)*u1*r2],[0,0,1],p1,-theta2)[0])
    return arc_2p(p2,p3,r,-1,s)

def i_p3d(l1,l2):
    """
    find intersection point (if exists) between 2 lines in 3d 
    """
    sec1=[l1,l2]
    n=len(sec1)
    a=array(sec1)[comb_list(n)]
    p0=a[:,0][:,0]
    p1=a[:,0][:,1]
    p2=a[:,1][:,0]
    p3=a[:,1][:,1]
    v1=a_(rot('z.00001',p1-p0))
    v2=a_(rot('x.00001',p3-p2))
    iim=a_([v1[:,:2],-v2[:,:2]]).transpose(1,0,2).transpose(0,2,1)
    im=inv(iim)
    p=(p2-p0)[:,:2]
    t=einsum('ijk,ik->ij',im,p)
    
    iim1=a_([v1[:,1:],-v2[:,1:]]).transpose(1,0,2).transpose(0,2,1)
    im1=inv(iim1)
    px=(p2-p0)[:,1:]
    t1=einsum('ijk,ik->ij',im1,px)
    
    i_p1=p0+einsum('ij,i->ij',v1,t[:,0])
    i_p2=p0+einsum('ij,i->ij',v1,t1[:,0])
    i_p1=l_(i_p1[0]) if l_(i_p1.round(4))==l_(i_p2.round(4)) else []
    return i_p1

def concave_hull_3d(pnts,n=3,engaging_angle=270):
    # pnts=rot('y.001',pnts)
    sec1=rot_sec2xy_plane(pnts)
    sec2=c32(sec1)
    sec3=concave_hull(sec2,n=n,engaging_angle=engaging_angle)
    a=a_(pnts[0])
    b=a_(c23(c32(pnts)[0]))
    v1=uv(best_fit_plane(pnts)[0])
    if (l_(a_(v1).round(4))==[0,0,-1]) | (l_(a_(v1).round(4))==[0,0,1]):
        sec5=translate(a-b,sec3)
    else:
        v2=[0,0,-1]
        theta=l_(r2d(arccos(a_(v1)@v2)))
        ax_1=cross(v1,v2)
        sec4=axis_rot_1(c23(sec3),ax_1,c23(sec2)[0],-theta)
        sec5=translate(a-b,sec4)
    return sec5

def convex_hull_3d(pnts):
    # pnts=rot('y.001',pnts)
    sec1=rot_sec2xy_plane(pnts)
    sec2=c32(sec1)
    sec3=convex_hull(sec2)
    a=a_(pnts[0])
    b=a_(c23(c32(pnts)[0]))
    v1=uv(best_fit_plane(pnts)[0])
    if (l_(a_(v1).round(4))==[0,0,-1]) | (l_(a_(v1).round(4))==[0,0,1]):
        sec5=translate(a-b,sec3)
    else:
        v2=[0,0,-1]
        theta=l_(r2d(arccos(a_(v1)@v2)))
        ax_1=cross(v1,v2)
        sec4=axis_rot_1(c23(sec3),ax_1,c23(sec2)[0],-theta)
        sec5=translate(a-b,sec4)
    return sec5

def ip_unordered(sol1,sol2):
    '''
    calculates the maximum number of intersection points between 2 solids in random order
    '''
    def unordered_ip(sol1,sol2):
        a=sol1
        b=sol2
        v,f1=vnf1(b)
        x=concatenate([seg(p) for p in f1]).min(1)
        y=concatenate([seg(p) for p in f1]).max(1)
        c=remove_duplicates(a_([x,y]).transpose(1,0))
        l1=a_(v)[c]
        la=l1[:,0]
        lb=l1[:,1]
        lab=lb-la
        v,f1=vnf1(a)
        tr=a_(v)[f1]
        p0,p1,p2=tr[:,0],tr[:,1],tr[:,2]
        p01=p1-p0
        p02=p2-p0
        # la+lab*t1=p0+p0p1*t2+p0p2*t3
        # lab*t1-p0p1*t2-p0p2*t3=p0-la
        t=einsum('jk,ijk->ij',cross(p01,p02),(la[:,None]-p0[None,:]))/ \
        einsum('ik,jk->ij',(-lab),cross(p01,p02))
        u=einsum('ijk,ijk->ij',cross(p02[None,:],-lab[:,None]),(la[:,None]-p0[None,:]))/ \
        einsum('ik,jk->ij',(-lab),cross(p01,p02))
        v=einsum('ijk,ijk->ij',cross(-lab[:,None],p01[None,:]),(la[:,None]-p0[None,:]))/ \
        einsum('ik,jk->ij',(-lab),cross(p01,p02))
        condition=(t>=0)&(t<=1)&(u>=0)&(u<=1)&(v>=0)&(v<=1)&(u+v<1)
        l2=l_((la[:,None]+lab[:,None]*t[:,:,None])[condition])
        return l2
    
    a=unordered_ip(sol1,sol2)
    b=unordered_ip(sol2,sol1)
    return a+b

def ip_unsorted(sol1,sol2):
    a=ip(sol1,sol2)
    b=ip(sol2,sol1)
    return a+b

def xyc(pnts):
    '''
    makes z-coordinate 0 of the 3d points list
    '''
    s1=a_(pnts).shape
    s2=a_(pnts).reshape(-1,s1[-1])
    a=s2[:,0]
    b=s2[:,1]
    c=zeros(len(s2))
    return l_(a_([a,b,c]).transpose(1,0).reshape(s1))

def xzc(pnts):
    '''
    makes y-coordinate 0 of the 3d points list
    '''
    s1=a_(pnts).shape
    s2=a_(pnts).reshape(-1,s1[-1])
    a=s2[:,0]
    b=zeros(len(s2))
    c=s2[:,2]
    return l_(a_([a,b,c]).transpose(1,0).reshape(s1))

def yzc(pnts):
    '''
    makes x-coordinate 0 of the 3d points list
    '''
    s1=a_(pnts).shape
    s2=a_(pnts).reshape(-1,s1[-1])
    a=zeros(len(s2))
    b=s2[:,1]
    c=s2[:,2]
    return l_(a_([a,b,c]).transpose(1,0).reshape(s1))

def bp_mesh_from_points(pnts,normals,n=10000,f=3):
    '''
    creates a ball pivoting mesh from a list of points and normals
    'n' is the max number of triangles required in the resultant mesh
    'f' is a factor which can be increased to increase the ball size defining the mesh
    '''
    pcd=o3d.geometry.PointCloud()
    pcd.points=o3d.utility.Vector3dVector(pnts)
    pcd.normals=o3d.utility.Vector3dVector(normals)
    distances=pcd.compute_nearest_neighbor_distance()
    avg_dist=mean(distances)
    radius=f*avg_dist
    bp_mesh=o3d.geometry.TriangleMesh.create_from_point_cloud_ball_pivoting(pcd,
                o3d.utility.DoubleVector([radius,radius*2]))
    dc_mesh=bp_mesh.simplify_quadric_decimation(n)
    tri_m=l_(a_(dc_mesh.vertices)[a_(dc_mesh.triangles)])
    return tri_m


def poisson_mesh_from_points(pnts,normals,n=10000):
    '''
    creates poisson's mesh from points and normals
    'n' is the max number of triangles required in the resultant mesh
    '''
    pcd=o3d.geometry.PointCloud()
    pcd.points=o3d.utility.Vector3dVector(pnts)
    pcd.normals=o3d.utility.Vector3dVector(normals)
    bp_mesh=o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(pcd,depth=8,width=0,scale=1.1,linear_fit=False)
    dc_mesh=bp_mesh[0].simplify_quadric_decimation(n)
    tri_m=l_(a_(dc_mesh.vertices)[a_(dc_mesh.triangles)])
    return tri_m

def alpha_mesh_from_points(pnts,normals,alpha=.5,n=10000):
    '''
    creates alpha shape mesh from pointa and normals.
    'alpha' is a factor need to be defined based on the mesh quality, 
    it can be any number upto around 20
    'n' is the max number of triangles required in the resultant mesh
    '''
    pcd=o3d.geometry.PointCloud()
    pcd.points=o3d.utility.Vector3dVector(pnts)
    pcd.normals=o3d.utility.Vector3dVector(normals)
    bp_mesh=o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(pcd,alpha)
    dc_mesh=bp_mesh.simplify_quadric_decimation(n)
    tri_m=l_(a_(dc_mesh.vertices)[a_(dc_mesh.triangles)])
    return tri_m

def points_near_line(list_of_points, line, dist=.2):
    '''
    find all the points which are at < distance 'dist' from the line from the list_of_points
    '''
    a=a_(list_of_points)
    b=a_(seg(line)[:-1])
    # a=concatenate([l1,l2,l3])
    # b=a_(seg(cpo(sol2)[0])[:-1])
    v1=a_(lines2vectors(b))
    d1=norm(v1,axis=1)
    v2=a[None,:]-b[:,0][:,None]
    x1=einsum('ik,ijk->ij',v1,v2)
    x2=einsum('ij,i->ij',x1,1/d1)
    x3=norm(cross(v1[:,None],v2),axis=2)/d1[:,None]
    
    pnts=l_(a_([a]*len(b))[(x2<=d1[:,None])&(x2>=0)&(x3<=dist)])
    pnts=l_(a_(pnts)[a_([timeToReachPoint(p,line,dist) for p in pnts]).argsort()])
    return pnts

def sol2path_closed(sol,path):
    """
function to extrude a solid to a closed loop path
example:
c1=circle(5)
p1=rot2d(90,sinewave(100,5,1,50))
s3=prism(c1,p1)
p2=c23(arc(20,0,360,[50,50],51))
s3=sol2path_closed(s3,p2)
fo(f'''
{swp(s3)}
''')
    """
    sol1=c3t2(sol)
    zpath=[[0,0,p[0][2]] for p in sol]
    path2=path2path1(zpath,path)
    sol2=align_sol_1(path_extrude2msec_closed(sol1,path2))
    return sol2

def path_extrude2msec_closed(sec_list,path):
    p1=path
    p2=path[1:]+[path[0]]
    p1,p2=array([p1,p2])
    v1=p2-p1
    u1=v1/norm(v1,axis=1).reshape(-1,1)
    v2=concatenate([[(u1[-1]+u1[0])/2], (u1[1:]+u1[:-1])/2])
    sec2=[]
    for i in range(len(path)):
        sec1=translate(path[i],sec2vector(v2[i],sec_list[i]))
        sec2.append(sec1)
    sec2=sec2+[sec2[0]]
    return sec2