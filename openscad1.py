from numpy import *
from numpy.linalg import *
import matplotlib.pyplot as plt
import time
from scipy.spatial import cKDTree,Delaunay, Voronoi
import pandas as pd
# import igl



def arc(radius=0,start_angle=0,end_angle=0,cp=[0,0],s=20):
    '''
    function for calculating 2d arc
    'cp': center point of the arc
    's': number of segments in the arc
    refer file "example of various functions" for application example
    '''
    cp=array(cp)
    r=linspace(start_angle,end_angle,s+1)
    x=radius*cos(pi/180*r)
    y=radius*sin(pi/180*r)
    c=(cp+array([x,y]).swapaxes(0,1))
    return c.tolist()        

def pts(p):
    '''
    calculates the cumulative sum of 2d list of points 'p'
    e.g.
    pts([[0,0],[4,0],[2,3],[5,-8]]) will produce following output
    [[0, 0], [4, 0], [6, 3], [11, -5]]
    '''
    return array(p)[:,0:2].cumsum(axis=0).tolist()




def pts1(p):
    '''
    'p' is a list of points
    function calculates the cumulative sum of x,y values in the list while z value remains the same.
    this is mainly used in function cr(pl,s).
    example:
    pts1([[0,0,1],[10,0,1],[0,5,1],[-10,0,1]]) => [[0, 0, 1], [10, 0, 1], [10, 5, 1], [0, 5, 1]]
    
    if used with function cr(pl,s)
    cr(pts1([[0,0,1],[10,0,1],[0,5,1],[-10,0,1]]),5) => This is a rounded rectangle of dim 10 x 5 with corner radius of 1 at each corner
    refer to file 'example of various functions' for application example
    '''
    
    p=[[a[0],a[1],0] if len(a)==2 else a for a in p]
    b=array(p)[:,0:2].cumsum(axis=0)
    c=array([array(p)[:,2].tolist()])
    return concatenate((b,c.T),1).tolist()


    
def cw(sec):
    '''
    function to identify if an enclosed section is clockwise(cw) or counterclockwise(ccw)
    this returns 1 if section is clockwise and -1 if it is counterclockwise
    '''
    if len(sec)==3:
        p=array(sec)
        return -1 if cross(p[1]-p[0],p[2]-p[0])>0 else 1
    else:
        cp1=array(sec).mean(0)
        p0=[array(sec)[:,0].min()-1,cp1[1]]
        v1=array([[1,0]]*len(sec))

        p1=array(sec)
        p2=array(sec[1:]+[sec[0]])
        v2=p2-p1+.0000001
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

        cw=-1 if cross(array(p4)-array(p3),array(ip2[0])-array(p3))>0 else 1
        return cw


def cwv(sec):
    '''
    function to identify whether each point in a section is clockwise or counter clockwise. 
    cw(sec)==1 means clockwise and -1 means counterclockwise. 
    e.g.
    cwv(pts([[0,0],[4,0],[0,4],[2,0],[0,2],[-6,0]])) => [-1, -1, 1, -1, -1, -1]
    '''
    p=sec
    p0=[p[-1]]+p[:-1]
    p1=p
    p2=p[1:]+[p[0]]
    p0,p1,p2=array([p0,p1,p2])
    p=array([p0,p1,p2]).transpose(1,0,2)
    return [ -1 if cross(p1[1]-p1[0],p1[2]-p1[0])>=0 else 1 for p1 in p]


def ang(x,y):
    '''
function to calculate angle of a 2d vector starting from origin and end point with x and y co-ordinates
 example:
 p1,p2=array([[3,4],[-3,2]])
 v=p2-p1
 ang= ang(v[0],v[1])
 
    '''
    if x>=0 and y>=0:
        return arctan(y/(0.000001 if x==0 else x))*180/pi
    elif x<0 and y>=0:
        return 180-abs(arctan(y/x))*180/pi
    elif  x<0 and y<0:
        return 180+abs(arctan(y/x))*180/pi
    else:
        return 360-abs(arctan(y/(0.000001 if x==0 else x)))*180/pi



def q(vector=[1,0,0],point=[0,5,0],theta=0):
    '''
    function to rotate a point around a vector(axis) with angle theta
    example:
    q(vector=[1,0,0],point=[0,5,0],theta=90)
    output: [0,0,5]
    '''

    t=theta
    v=vector/(norm(vector))
    a=t/2*pi/180
    p=[cos(a),multiply(v,sin(a))]
    p1=[p[0],-p[1]]
    q=[0,[point[0],point[1],0] if len(point)==2 else point]
    pq=[p[0]*q[0]-p[1]@q[1],multiply(p[0],q[1])+p[1]*q[0]+cross(p[1],q[1])]
    pqp1=[pq[0]*p1[0]-pq[1]@p1[1],pq[0]*p1[1]+pq[1]*p1[0]+cross(pq[1],p1[1])]
    transformation=pqp1[1].tolist()
    return transformation

def uv(v):
    '''
    function to calculate unit vector of a given vector
    example:
    vector=[2,3,5]
    unit_vector=uv(vector) => [0.3244428422615251, 0.48666426339228763, 0.8111071056538127]
    '''
    v=array(v)
    return (v/norm(v)).tolist()

# def norm(v):
#     return norm(v)

def fillet2d(pl,rl,s):
    p0=array(array(pl)[len(pl)-2:len(pl)].tolist()+array(pl)[0:len(pl)-2].tolist())
    p1=array([array(pl)[len(pl)-1].tolist()]+array(pl)[0:len(pl)-1].tolist())
    p2=array(pl)
    p3=array(array(pl)[1:len(pl)].tolist()+[array(pl)[0].tolist()])
    p4=array(array(pl)[2:len(pl)].tolist()+array(pl)[0:2].tolist())
    r0=array([array(rl)[len(rl)-1].tolist()]+array(rl)[0:len(rl)-1].tolist())
    r1=array(rl)
    r2=array(array(rl)[1:len(rl)].tolist()+[array(rl)[0].tolist()])
    u0=(p0-p1)/(norm(p0-p1,axis=1)).reshape(-1,1)
    u1=(p2-p1)/(norm(p2-p1,axis=1)).reshape(-1,1)
    u2=(p1-p2)/(norm(p1-p2,axis=1)).reshape(-1,1)
    u3=(p3-p2)/(norm(p3-p2,axis=1)).reshape(-1,1)
    u4=(p2-p3)/(norm(p2-p3,axis=1)).reshape(-1,1)
    u5=(p4-p3)/(norm(p4-p3,axis=1)).reshape(-1,1)
    theta0= (180-arccos(einsum('ij,ij->i',u0,u1))*180/pi)/2
    theta1= (180-arccos(einsum('ij,ij->i',u2,u3))*180/pi)/2
    theta2= (180-arccos(einsum('ij,ij->i',u4,u5))*180/pi)/2
    return f2d(p1,p2,p3,r0,r1,r2,theta0,theta1,theta2,u2,u3,s)
    


def each(a):
    c=[]
    for p in a:
        for p1 in p:
            c.append(p1)
    return c

def cr1(pl,s=20):
    pl1=array(pl)[:,0:2].tolist()
    rl=[0 if len(p)==2 else p[2] for p in pl]
    return fillet2d(pl1,rl,s)

def cr(pl,s=20):
    '''
    function to create section with corner radiuses. e.g. 
    following code has 3 points at [0,0],[10,0] and [7,15] and radiuses of 0.5,2 and 1 respectively,
    s=5 represent the number of segments at each corner radius.
    sec=cr(pl=[[0,0,.5],[10,0,2],[7,15,1]],s=5)
    
    refer file "example of various functions" for application
    '''
    sec=array(cr1(pl,s)).round(8)
    s1=sec[sort(unique(sec,axis=0,return_index=True)[1])].tolist()
    return s1

def cr_c(pl,s=20):
    sec=array(cr1(pl,s)).round(8)
    s1=sec[sort(unique(sec,axis=0,return_index=True)[1])].tolist()
    p0,p1=array([s1[len(s1)-1],s1[0]])
    v=p1-p0
    p=(p0+v*.999).tolist()
    
    return s1+[p]

def f2d(p1,p2,p3,r0,r1,r2,theta0,theta1,theta2,u2,u3,s):
    l1=norm(p1-p2,axis=1)
    l2=r0*tan(theta0*pi/180)+r1*tan(theta1*pi/180)
    l3=norm(p3-p2,axis=1)
    l4=r1*tan(theta1*pi/180)+r2*tan(theta2*pi/180)
    rf1=[r1[i] if l1[i]>l2[i] else 0 if l2[i]==0 else l1[i]/l2[i]*r1[i] for i in range(len(l1))]
    rf2=[r1[i] if l3[i]>l4[i] else 0 if l4[i]==0 else l3[i]/l4[i]*r1[i] for i in range(len(l3))]
    rf=swapaxes([rf1,rf2],0,1).min(axis=1)
    p=p2+u2*(rf*tan(theta1*pi/180)).reshape(-1,1)
    q=array([p1,p2,p3]).transpose(1,0,2)

    r=array([-1 if cross(p[1]-p[0],p[2]-p[0])>0 else 1 for p in q])
    n=r==-1
    n1=p-u2@array(rm(90))*rf.reshape(-1,1)
    n2=p-u2@array(rm(-90))*rf.reshape(-1,1)
    cp=[]
    for i in range(len(n)):
        if n[i]==True:
            cp.append(n1[i])
        else:
            cp.append(n2[i])

    cp=array(cp)
    a1=[]
#     alpha=(p-cp)/norm(p-cp,axis=1).reshape(-1,1)
    alpha=[ [0,0] if norm(p[i]-cp[i])==0 else (p[i]-cp[i])/norm(p[i]-cp[i]) for i in range(len(p))]
    for i in range(len(alpha)):
        a1.append(ang(alpha[i][0],alpha[i][1]))
    a1=array(a1)
    boo=[]
    for i in range(len(p1)):
        boo.append(cw([p1[i],p2[i],p3[i]]))
    boo=array(boo)   
    a2=where(boo==-1,a1+2*theta1,a1-2*theta1)
    ar=[]
    for i in range(len(rf)):
        ar.append(arc(rf[i],a1[i],a2[i],cp[i],s))
    ar=array(ar)
    c1=r1==0
    c2=norm(u2-u3,axis=1)<.2
    d=[]
    for i in range(len(c1)):
        if c1[i] or c2[i]:
            d.append([p2[i].tolist()])
        else:
            d.append(ar[i].tolist())
    return concatenate(d).tolist()


def flip(sec):
    '''
    function to flip the sequence of a list or a list of points
    example:
    list=[1,2,3,4,5]
    flipped_list=flip(list) => [5, 4, 3, 2, 1]
    
    list=[[1,2,3],[4,5,6],[7,8,9]]
    flipped_list=flip(list) => [[7, 8, 9], [4, 5, 6], [1, 2, 3]]
    '''
    return sec[::-1]
    

def max_r(sec):
    '''
    function calculates the maximum radius in a given closed section
    example:
    sec=cr_c(pts1([[0,0,.2],[8,3,3],[5,7,1],[-8,0,2],[-5,20,1]]),20)
    max_r(sec) => 3.0
    
    '''
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
    return max(c)
        

def offset_l(l,d):
    u=uv(subtract(l[1],l[0]))
    p0=add(l[0],dot(u,multiply(d,rm(-90)))).tolist()
    p1=add(l[1],dot(u,multiply(d,rm(-90)))).tolist()
    return [p0,p1]

def seg(sec):
    '''
    function to create a segment from a list of points or a list
    example:
    list=[1,2,3,4,5,6]
    seg(list)=> [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 1]]
    
    list=[[1,2,3],[4,5,6],[7,8,9]]
    seg(list) => [[[1, 2, 3], [4, 5, 6]], [[4, 5, 6], [7, 8, 9]], [[7, 8, 9], [1, 2, 3]]]
    '''
    c=[]
    for i in range(len(sec)):
        i_plus=i+1 if i<len(sec)-1 else 0
        p0=sec[i]
        p1=sec[i_plus]
        l=[p0,p1]
        c.append(l)
    return c


def offset_segv(sec,d):
    '''
    function makes the segments of the original section and offset each segment by a distance 'd'
    refer the file "example of various functions" for application examples
    
    '''
    s=sec
    s1=s[1:]+[s[0]]
    x=(array(s1)-array(s))
    y=norm(x,axis=1)
    u=x/y.reshape(-1,1)
    p0=array(s)+u@array(rm(-90))*d
    p1=array(s1)+u@array(rm(-90))*d
    return swapaxes([p0,p1],0,1).tolist()

def offset_points(sec,r):
    '''
    function to calculate offset of a list of 2d points
    in defining sections, providing corner radius is a must
    e.g. pts([[0,0],[10,0],[0,5],[-10,0]]) will fail
    while cr(pts1([[0,0,.1],[10,0,.1],[0,5,.1],[-10,0,.1]])) will work perfectly
    refer the file "example of various functions" for application examples
    '''
    return array(offset_segv(sec,r))[:,0].tolist()


def offset_pointsv(sec,r):
    '''
    function to calculate offset of a list of 2d points
    in defining sections, providing corner radius is a must
    e.g. pts([[0,0],[10,0],[0,5],[-10,0]]) will fail
    while cr(pts1([[0,0,.1],[10,0,.1],[0,5,.1],[-10,0,.1]])) will work perfectly
    refer the file "example of various functions" for application examples
    '''
    return array(offset_segv(sec,r))[:,0].tolist()

def offset_seg_cw(sec,r):
    '''
    function offsets the segment only when the point is clockwise
    refer to file 'example of various functions' for application example
    '''
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
    '''
    function removes all the duplicates from a 
    example:
    list=[9,5,1,2,3,4,5,2,4,6,9]
    remove_extra_points(list) => [9, 5, 1, 2, 3, 4, 6]
    
    list=[[7,8,9],[1,2,3],[10,11,12],[4,5,6],[7,8,9],[10,11,12]]
    remove_extra_points(list) => [[7, 8, 9], [1, 2, 3], [10, 11, 12], [4, 5, 6]]
    '''
    return array(points_list)[sort(unique(points_list,axis=0,return_index=True)[1])].tolist()

def convert_secv(sec):
    '''
    function removes all the radiuses from the section 'sec' where points are ccw
    example:
    sec=cr_c(pts1([[0,0,.1],[7,5,2],[5,7,3],[-5,7,5],[-7,5,5]]),20)
    sec1=convert_secv(sec)
    sec1 will remove all the radius in 'sec' where points are ccw
    refer to file "examples of various functions" for application example
    '''
    a=list_r(sec)
    if (a==a[0]).all() or list_r(sec).max()>array(bb2d(sec)).max():
        return sec
    else:
        d=max_rv(sec)+1
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
        arr=swapaxes([pi_minus,p_i,pi_plus],0,1)
        clock=array([-1 if cross(p[1]-p[0],p[2]-p[0])>0 else 1 for p in arr])
        c1=where(r==0,True,False)
        c2=where(r>=d,True,False)
        c3=where(clock==1,True,False)
        p=array(sec)[c1 | c2 | c3 ]
        p1=cKDTree(array(sec)).query(p)[1].tolist()
        p2=[p1[len(p1)-1]]+p1[0:len(p1)-1]
        p3=p1[1:len(p1)]+[p1[0]]
        p4=p1[2:len(p1)]+p1[0:2]
        a=i_p2dv(array(sec)[p2],array(sec)[p1],array(sec)[p3],array(sec)[p4])
        b=array(sec)[p1]
        c=array(p3)-array(p1)>1
        d=[]
        for i in range(len(c)):
            if c[i]==True:
                d.append(a[i].tolist())
            else:
                d.append(b[i].tolist())
        d_minus=[d[len(d)-1]]+d[0:len(d)-1]
        d_plus=d[1:len(d)]+[d[0]]
        va=array(d)-array(d_minus)
        vb=array(d_plus)-array(d_minus)
        normva=1/norm(va,axis=1)
        normvb=1/norm(vb,axis=1)
        ua=einsum('ij,i->ij',va,normva)
        ub=einsum('ij,i->ij',vb,normvb)
        sec1=array(d)[(ua!=ub).all(axis=1)].tolist()
        a=[sec1[len(sec1)-1]]+sec1[:-1]
        b=sec1
        c=sec1[1:]+[sec1[0]]
        a,b,c=array([a,b,c])
        v1,v2=b-a,c-a
        n1=1/einsum('ij,ij->i',v1,v1)**.5
        n2=1/einsum('ij,ij->i',v2,v2)**.5
        u1=einsum('ij,i->ij',v1,n1).round(3)
        u2=einsum('ij,i->ij',v2,n2).round(3)
        decision=~(u1==u2).all(1)
        return b[decision].tolist()


def convert_secv1(sec):
    '''
    function removes all the radiuses from the section 'sec'
    example:
    sec=cr_c(pts1([[0,0,.1],[7,5,2],[5,7,3],[-5,7,5],[-7,5,5]]),20)
    sec1=convert_secv1(sec)
    sec1 will remove all the radius in 'sec' 
    refer to file "examples of various functions" for application example
    '''
    a=list_r(sec)
    if (a==a[0]).all() or list_r(sec).max()>array(bb2d(sec)).max():
        return sec
    else:
        d=max_rv(sec)+1
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
        arr=swapaxes([pi_minus,p_i,pi_plus],0,1)
        clock=array([-1 if cross(p[1]-p[0],p[2]-p[0])>0 else 1 for p in arr])
        c1=where(r==0,True,False)
        c2=where(r>=d,True,False)
        c3=where(clock==-1,True,False)
        p=array(sec)[c1 | c2 ]
        p1=cKDTree(array(sec)).query(p)[1].tolist()
        p2=[p1[len(p1)-1]]+p1[0:len(p1)-1]
        p3=p1[1:len(p1)]+[p1[0]]
        p4=p1[2:len(p1)]+p1[0:2]
        a=i_p2dv(array(sec)[p2],array(sec)[p1],array(sec)[p3],array(sec)[p4])
        b=array(sec)[p1]
        c=array(p3)-array(p1)>1
        d=[]
        for i in range(len(c)):
            if c[i]==True:
                d.append(a[i].tolist())
            else:
                d.append(b[i].tolist())
        d_minus=[d[len(d)-1]]+d[0:len(d)-1]
        d_plus=d[1:len(d)]+[d[0]]
        va=array(d)-array(d_minus)
        vb=array(d_plus)-array(d_minus)
        normva=1/norm(va,axis=1)
        normvb=1/norm(vb,axis=1)
        ua=einsum('ij,i->ij',va,normva)
        ub=einsum('ij,i->ij',vb,normvb)
        sec1=array(d)[(ua!=ub).all(axis=1)].tolist()
        a=[sec1[len(sec1)-1]]+sec1[:-1]
        b=sec1
        c=sec1[1:]+[sec1[0]]
        a,b,c=array([a,b,c])
        v1,v2=b-a,c-a
        n1=1/einsum('ij,ij->i',v1,v1)**.5
        n2=1/einsum('ij,ij->i',v2,v2)**.5
        u1=einsum('ij,i->ij',v1,n1).round(3)
        u2=einsum('ij,i->ij',v2,n2).round(3)
        decision=~(u1==u2).all(1)
        return b[decision].tolist()


def convert_secv2(sec,d):
    '''
    function removes all the radiuses from the section 'sec' where points are cw
    example:
    sec=cr_c(pts1([[0,0,.1],[7,5,2],[5,7,3],[-5,7,5],[-7,5,5]]),20)
    sec1=convert_secv2(sec,5)
    sec1 will remove all the radius in 'sec' where points are ccw
    refer to file "examples of various functions" for application example
    '''
    a=list_r(sec)
    if (a==a[0]).all() or list_r(sec).max()>array(bb2d(sec)).max():
        return sec
    else:
        
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
        arr=swapaxes([pi_minus,p_i,pi_plus],0,1)
        clock=array([-1 if cross(p[1]-p[0],p[2]-p[0])>0 else 1 for p in arr])
        c1=where(r==0,True,False)
        c2=where(r>=d,True,False)
        c3=where(clock==-1,True,False)
        p=array(sec)[c1 | c2 | c3 ]
        p1=cKDTree(array(sec)).query(p)[1].tolist()
        p2=[p1[len(p1)-1]]+p1[0:len(p1)-1]
        p3=p1[1:len(p1)]+[p1[0]]
        p4=p1[2:len(p1)]+p1[0:2]
        a=i_p2dv(array(sec)[p2],array(sec)[p1],array(sec)[p3],array(sec)[p4])
        b=array(sec)[p1]
        c=array(p3)-array(p1)>1
        d=[]
        for i in range(len(c)):
            if c[i]==True:
                d.append(a[i].tolist())
            else:
                d.append(b[i].tolist())
        d_minus=[d[len(d)-1]]+d[0:len(d)-1]
        d_plus=d[1:len(d)]+[d[0]]
        va=array(d)-array(d_minus)
        vb=array(d_plus)-array(d_minus)
        normva=1/norm(va,axis=1)
        normvb=1/norm(vb,axis=1)
        ua=einsum('ij,i->ij',va,normva)
        ub=einsum('ij,i->ij',vb,normvb)
        sec1=array(d)[(ua!=ub).all(axis=1)].tolist()
        a=[sec1[len(sec1)-1]]+sec1[:-1]
        b=sec1
        c=sec1[1:]+[sec1[0]]
        a,b,c=array([a,b,c])
        v1,v2=b-a,c-a
        n1=1/einsum('ij,ij->i',v1,v1)**.5
        n2=1/einsum('ij,ij->i',v2,v2)**.5
        u1=einsum('ij,i->ij',v1,n1).round(3)
        u2=einsum('ij,i->ij',v2,n2).round(3)
        decision=~(u1==u2).all(1)
        return b[decision].tolist()



def list_r(sec):
    '''
    function list the corner radiuses of a given section (only where the radius is specified)
    example:
    sec=cr_c(pts1([[0,0,.1],[7,5,2],[5,7,3],[-5,7,5],[-7,5,5]]),5)
    list_r(sec) => 
    array([0.   , 0.1  , 0.1  , 0.1  , 0.1  , 0.   , 0.   , 2.   , 2.   ,
       2.   , 2.   , 0.   , 0.   , 3.   , 3.   , 3.   , 3.   , 0.   ,
       0.   , 4.077, 4.077, 4.077, 4.077, 4.077, 4.077, 4.077, 4.077,
       4.077, 0.   , 0.   ])
    '''
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
    '''
    calculates list of radiuses for all the points of a given section
    sec=cr_c(pts1([[0,0,.1],[7,5,2],[5,7,3],[-5,7,5],[-7,5,5]]),5)
    list_ra(sec) =>
    array([  0.,   0.,   0.,   0.,   0.,  19., 124.,   2.,   2.,   2.,   2.,
        95.,  28.,   3.,   3.,   3.,   3.,  26.,  92.,   4.,   4.,   4.,
         4.,   4.,   4.,   4.,   4.,   4.,  40.,   8.])
    
    '''
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

def rm(theta):
    '''
    function to rotate a vector by "theta" degrees e.g. try following code:
    line=[[0,0],[5,3]]
    line1=array(line)@rm(30)
    line1=line1.tolist()
    refer file "examples of various functions" for application
    '''
    pi=3.141592653589793
    return [[cos(theta * pi/180),sin(theta * pi/180)],[-sin(theta * pi/180),cos(theta * pi/180)]]

def max_rv(sec):
    '''
    function calculates the maximum radius in a given closed section
    example:
    sec=cr_c(pts1([[0,0,.2],[8,3,3],[5,7,1],[-8,0,2],[-5,20,1]]),20)
    max_rv(sec) => 3.0
    
    '''
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
    '''
    function calculates radius of the circle drawn with 3 points 'p1','p2','p3'
    example
    p1,p2,p3=[3,0],[0,0],[0,3]
    radius=r_3p([p1,p2,p3]) => 2.1213203435596424
    '''
    p4=add(p[0],divide(subtract(p[1],p[0]),2)).tolist()
    p5=add(p[1],divide(subtract(p[2],p[1]),2)).tolist()
    u1=uv(subtract(p[1],p4))
    u2=uv(subtract(p[2],p5))
    p6=add(p4,dot(u1,rm(90))).tolist()
    p7=add(p5,dot(u2,rm(90))).tolist()
    cp=i_p2d([p4,p6],[p5,p7])
    r=norm(subtract(p[0],cp))
    return r


def i_p2d(l1,l2):
    '''
    function to calculate the intersection point between 2 lines in 2d space
    e.g. i_p2d(l1=[[0,0],[1,4]],l2=[[10,0],[7,2]]) =>  [1.42857, 5.71429]
    '''
    p0,p1,p2,p3=l1[0],l1[1],l2[0],l2[1]
    p0,p1,p2,p3=array([p0,p1,p2,p3])
    v1=p1-p0
    v2=p3-p2
    im=pinv(array([v1,-v2]).T)
    t1=(im@(p2-p0))[0]
    ip=p0+v1*t1
    
    return ip.tolist()

def s_int(s):
    '''
    calulates the self intersection points of a list of line segments 's'
    it also picks the points in case the 2 lines are just connected at 1 point and are not crossing
    refer to file 'example of various functions' for application example
    '''
    c=[]
    for i in range(len(s)):
        p0=array([s[i]]*len(s))[:,0]
        p1=array([s[i]]*len(s))[:,1]
        v1=p1-p0
        p2=array(s)[:,0]
        p3=array(s)[:,1]
        v2=p3-p2
        m=swapaxes([swapaxes([v1.T[0],-v2.T[0]],0,1),swapaxes([v1.T[1],-v2.T[1]],0,1)],0,1)
        n=m[where(det(m)!=0)]
        pa=p0[where(det(m)!=0)]
        pb=p2[where(det(m)!=0)]
        v=v1[where(det(m)!=0)]
        A=inv(n)
        B=pb-pa
        def mul(a,b):
            return a@b
        t=einsum('ijk,ik->ij',A,B)[:,0].round(4)
        u=einsum('ijk,ik->ij',A,B)[:,1].round(4)
        t1=where(t>=0,where(t<=1,True,False),False)
        u1=where(u>=0,where(u<=1,True,False),False)
        d=(pa+v*t.reshape(-1,1))[where(t1&u1==True)].tolist()
        if d!=[]:
            c=c+d
    return c


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

def i_p2dv(p0,p1,p2,p3):
    v1=p1-p0
    v2=p3-p2
    a=pinv(swapaxes(transpose(array([v1,-v2])),0,1))
    b=p2-p0
    t=einsum('ijk,ik->ij',a,b)[:,0]
    return p0+einsum('ij,i->ij',v1,t)

def sort_points(sec,list1):
    '''
    function picks the nearest point of a section from a reference section and matches the length of points for the 2 compared sections
    refer file "example of various functions" for application example
    
    '''
    return array(list1)[cKDTree(list1).query(sec)[1]].tolist()
            
def sort_pointsv(sec,list1):
    '''
    function picks the nearest point of a section from a reference section and matches the length of points for the 2 compared sections
    refer file "example of various functions" for application example
    '''
    return array(list1)[cKDTree(list1).query(sec)[1]].tolist()



def m_points(sec,sl=20):
    '''
    multiple points within straight lines of a closed section 'sec' with equal segment length 'sl' in the straight line segments
    refer file "example of various functions" for application example
    '''
    sec1=[]
    for p in seg(sec):
        n=1 if int(round(l_len(p)/sl,0))==0 else int(round(l_len(p)/sl,0))
        sec1.append(ls(p,n))
    return concatenate(sec1).tolist()
    

def m_points_o(sec,sl=20):
    '''
    multiple points within straight lines of a open section 'sec' with equal segment length 'sl' in the straight line segments
    refer file "example of various functions" for application example
    '''
    sec1=[]
    for p in seg(sec)[:-1]:
        n=1 if int(round(l_len(p)/sl,0))==0 else int(round(l_len(p)/sl,0))
        sec1.append(ls(p,n))
    return concatenate(sec1).tolist()+[sec[len(sec)-1]]


def ls(line,n):
    '''
    function to draw number of points 'n' in a line 'line'
    example:
    line=[[0,0],[10,0]]
    line1=ls(line,5) => [[0.0, 0.0], [2.0, 0.0], [4.0, 0.0], [6.0, 0.0], [8.0, 0.0], [10.0, 0.0]]
    '''
    p0,p1=array(line)
    v1=p1-p0
    return array([p0+v1/n*i for i in range(n)]).tolist()


def l_len(l):
    '''
    calculates length of a line 'l'
    example:
    line=[[0,0],[10,0]]
    l_len(line) =>10
    '''
    p0,p1=array(l[0]),array(l[1])
    v=p1-p0
    u=[v/(norm(v)+.00001)]
    length=norm(v)
    return length.tolist()

def arc_2p(p1,p2,r,cw=1,s=20):
    '''
    arc with 2 points 'p1,p2' with radius 'r' and with orientation clockwise (1) or counterclock wise(-1)
    refer file "example of various functions" for application example
    '''
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

def arc_long_2p(p1,p2,r,cw=1,s=20):
    '''
    long arc with 2 points 'p1,p2' with radius 'r' and with orientation clockwise (1) or counterclock wise(-1)
    
    refer file "example of various functions" for application example
    '''
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
    
    
def cir_2p(p1,p2,r,cw=1,s=20):
    '''
    circle with 2 points 'p1,p2' with radius 'r' and with orientation clockwise (1) or counterclock wise(-1)
    refer file "example of various functions" for application example
    '''
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




def arc_2p_cp(p1,p2,r,cw=-1):
    '''
    center point of an arc with 2 points 'p1,p2' with radius 'r' and with orientation clockwise (1) or counterclock wise(-1)
    
    refer file "example of various functions" for application example
    
    '''
    p1,p2=array([p1,p2])
    p3=p1+(p2-p1)/2
    d=norm(p3-p1)
    l=sqrt(abs(r**2-d**2))
    v=p1-p3
    u=v/norm(v)
    cp=p3+(u*l)@rm(-90 if cw==-1 else 90)
    return cp.tolist()

def offset(sec,r):
    '''
    calculates offset for a section 'sec' by amount 'r'
    refer file "example of various functions" for application example
    '''
#     return io(sec,r) if r<0 else sec if r==0 else oo_convex(sec,r) if convex(sec)==True else outer_offset(sec,r)
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


def prism(sec,path):
    '''
function to make a prism with combination of 2d section and 2d path
Example:
sec=circle(10)
path=cr(pts1([[2,0],[-2,0,2],[0,10,3],[-3,0]]),5)
prism=prism(sec,path)
    '''
    s1=flip(sec) if cw(sec)==1 else sec
    return [array(translate([0,0,y],offset(s1,round(x,3)))).tolist() for (x,y) in path]



def translate(p,sec):#translates a prism or section by [x,y,z] distance
    '''
    function to translate a group of points "sec" by "p" distance defined in [x,y,z].e.g. try following code:
    sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5)
    sec1=translate(p=[2,5,3],sec=sec)
    
    refer to file "example of various functions " for application
    '''
    p=[p[0],p[1],0] if len(p)==2 else p
    if len(array(sec).shape)==3:
        l,m,n=array(sec).shape
        a=array(sec).reshape(-1,n)
        if n==3:
            return array([a[:,0]+p[0],a[:,1]+p[1],a[:,2]+p[2]]).transpose(1,0).reshape(l,m,n).tolist()
        elif n==2:
            return array([a[:,0]+p[0],a[:,1]+p[1],p[2]]).transpose(1,0).reshape(l,m,n+1).tolist()
            
        
    elif len(array(sec).shape)==2:
        m,n=array(sec).shape
        a=array(sec).reshape(-1,n)
        if n==3:
            return array([a[:,0]+p[0],a[:,1]+p[1],a[:,2]+p[2]]).transpose(1,0).reshape(m,n).tolist()
        elif n==2:
            return array([a[:,0]+p[0],a[:,1]+p[1],array([p[2]]*m)]).transpose(1,0).reshape(m,n+1).tolist()



def prism1(sec,path,s=100,sp=0):
    '''
    produces a prism with each section divided in 's' number of points
    sp: this gives a flexibility of defining the starting point of the section
    '''
    sec0=equidistant_pathc(sec,s)
    sol=prism(sec0,path)
    sol=[p[sp:]+p[:sp] for p in sol]
    
    return sol

def prism2(sec,path,s=100,sp=0):
    '''
    similar to prism1 only order of offset is changed
    '''
    sol=prism(sec,path)
    sol=[equidistant_pathc(p,s) for p in sol]
    sol=[p[sp:]+p[:sp] for p in sol]
    
    return sol

    
def offset_points_cw(sec,r):
    '''
    function to offset only those points which are clockwise 
    refer file "example of various functions" for application example
    
    '''
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
    '''
    function to offset only those points which are counter clockwise 
    refer file "example of various functions" for application example
    
    '''
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
    '''
    function to convert the y co-ordinates to z co-ordinates e.g.[x,y]=>[x,0,y]. 2d to 3d coordinate system
    '''
    return [[p[0],0,p[1]] for p in path]

def surf_extrude(sec,path):# extrudes an open section 'sec' to a 'path' to create surface
    '''
    function to make surface with a polyline 2d sketch and a 3d path
    (there is no render here but points can be visualised with following command:
    for(p=surf_extrude(sec,path))points(p,.2);)
    example:
    sec2=cr(pts1([[-25,0],[10,5,5],[10,-3,10],[10,5,5],[10,-8,7],[10,1]]),10)  
    path2=cytz(cr(pts1([[-35,5,0],[10,8,20],[20,-5,10],[20,8,20],[10,-9,20],[10,1,0]]),10))
    surf2=surf_extrude(sec2,path2)
    
    refer file "example of various functions"
    '''
    p0=path
    p1=p0[1:]+[p0[0]]
    p0,p1=array(p0),array(p1)
    v=p1-p0
    a1=vectorize(ang)(v[:,0],v[:,1])
    b=sqrt(v[:,0]**2+v[:,1]**2)
    a2=vectorize(ang)(b,v[:,2])
    c=[]
    for i in range(len(path)-1):
        sec1=translate(p0[i],q_rot(['x90','z-90',f'y{-a2[i]}',f'z{a1[i]}'],sec))
        sec2=translate(p1[i],q_rot(['x90','z-90',f'y{-a2[i]}',f'z{a1[i]}'],sec))
        if i<len(path)-2:
            c.append([sec1])
        else:
            c.append([sec1,sec2])
    return concatenate(c).tolist()

def cpo(prism): # changes the orientation of points of a prism
    '''
    function to change the orientation of the points of the prism
    refer to the file "example of various functions" for application example
    
    
    '''
    return swapaxes(array(prism),0,1).tolist()

def c2t3(p):# converts 2d list to 3d
    '''
    function to convert 2d to 3d, it just adds the z-coordinate to the points list 
    example:
    list=c2t3([[1,2],[3,4],[6,7]])
    output=> [[1, 2, 0], [3, 4, 0], [6, 7, 0]]
    '''
    if len(array(p).shape)>2:
        return [translate([0,0,0],p) for p in p]
    else:
        return translate([0,0,0],p)

def c3t2(a): # converts 3d list to 2d list 
    '''
    function to convert 3d to 2d, it just removes the z-coordinate from the points list 
    example:
    list=c3t2([[1,2,3],[3,4,5],[6,7,8]])
    output=> [[1, 2], [3, 4], [6, 7]]
    '''
    if len(array(a).shape)==3:
        return array([ swapaxes([p[:,0],p[:,1]],0,1) for p in array(a)]).tolist()
    else:
        p=array(a)
        return swapaxes([p[:,0],p[:,1]],0,1).tolist()

def nv(p):# normal vector to the plane 'p' with atleast 3 known points
    '''
    given 3 points ['p1','p2',p3] function calculates unit normal vector
    example:
    p1,p2,p3=[1,0,0],[0,10,0],[-5,0,0]
    nv([p1,p2,p3]) => [0.0, 0.0, -1.0]
    '''
    l1=len(p)
    p0,p1,p2=array(translate([0,0,0],[p[0],p[int(l1/3)],p[int(l1*2/3)]]))
    nv=cross(p0-p1,p2-p1)
    m=1/norm(nv) if norm(nv)>0 else 1e5
    return (nv*m).tolist()

def fillet_3p_3d(p0,p1,p2,r,s):# fillet with 3 known points 'p0,p1,p2' in 3d space. 'r' is the radius of fillet and 's' is the number of segments in the fillet
    '''
    function to create fillet given 3 points 'p1','p2','p3' 
    r: radius of the fillet
    s: number of segments in the fillet
    refer file "example of various functions" for application example
    '''
    p0,p1,p2=array(translate([0,0,0],[p0,p1,p2]))
    n=array(nv([p0,p1,p2]))
    u1=(p0-p1)/(norm(p0-p1)+.00001)
    u2=(p2-p1)/(norm(p2-p1)+.00001)
    theta=(180-arccos(u1@u2)*180/pi)/2
    alpha=arccos(u1@u2)*180/pi
    l=r*tan(theta*pi/180)
    cp=p1+q(n,u1*r/cos(theta*pi/180),alpha/2)
    pa=p1+u1*l
    arc=[ cp+q(n,pa-cp,-i) for i in linspace(0,theta*2,s)]
    a,b,c=arc[0],arc[1:s-1],arc[s-1]
    return concatenate([[p1],arc]).tolist()

def fillet_3p_3d_cp(p0,p1,p2,r):# center point 'cp' of the fillet with 3 known points 'p0,p1,p2' in 3d space. 'r' is the radius of fillet
    '''
    function to find the center point of the fillet created by given 3 points 'p1','p2','p3' 
    r: radius of the fillet
    
    refer file "example of various functions" for application example
    '''
    p0,p1,p2=array(translate([0,0,0],[p0,p1,p2]))
    n=array(nv([p0,p1,p2]))
    u1=(p0-p1)/(norm(p0-p1)+.00001)
    u2=(p2-p1)/(norm(p2-p1)+.00001)
    theta=(180-arccos(u1@u2)*180/pi)/2
    alpha=arccos(u1@u2)*180/pi
    l=r*tan(theta*pi/180)
    cp=p1+q(n,u1*r/cos(theta*pi/180),alpha/2)
    return cp.tolist()

def i_p3d(l1,l2): # intersection point between 2 lines 'l1' and 'l2' in 3d space where both the lines are in the same plane
    '''
    function to calculate intersection point between 2 lines in 3d space 
    (only if these lines lie on the same plane)
    function is similar to i_p2d
    '''
    l1,l2=array(l1),array(l2)
    v1=l1[1]-l1[0]
    v2=l2[1]-l2[0]
    u1=v1/(norm(v1)+.00001)
    u2=v2/(norm(v2)+.00001)
    v3=l2[0]-l1[0]
    t1= (pinv(array([v1,-v2,[1,1,1]]).T)@array(v3))[0]
    ip=l1[0]+v1*t1
    return ip.tolist()

def arc_3p_3d(points,s=20):
    '''
    draws an arc through the 3 points list
    's' is the number of segments of the circle
    '''
    n1=array(nv(points))
    a1=cross(n1,[0,0,-1])
    t1=r2d(arccos(n1@[0,0,-1]))
    sec1=translate(-array(points).mean(0),points)
    sec2=c3t2(axis_rot(a1,sec1,t1))
    l1=len(sec2)
    p0,p1,p2=[sec2[0],sec2[int(l1/3)],sec2[int(l1*2/3)]]
    arc1=arc_3p(p0,p1,p2,s=s)
    arc1=translate(array(points).mean(0),axis_rot(a1,arc1,-t1))
    return arc1

def r_3p_3d(points):# radius of the circle with 3 known list of 'points' in 3d space
    '''
    function to find the radius of a circle created by 3 given points 'p1','p2','p3' in 3d space
    example:
    p1,p2,p3=[[3,0,0],[0,0,0],[0,3,2]]
    r_3p_3d([p1,p2,p3])=>1.8027906380190175
    '''
    points=array(points)
    v1=points[0]-points[1]
    v2=points[2]-points[1]
    u1=v1/(norm(v1)+.00001)
    u2=v2/(norm(v2)+.00001)
    n=cross(u1,u2)
    alpha=arccos(u1@u2)*180/pi
    pa=v1/2
    pb=v2/2
    pap=pa+q(n,u1,90)
    pbp=pb+q(n,u2,-90)
    l1=[pa,pap]
    l2=[pb,pbp]
    cp=i_p3d(l1,l2)
    v3=points[0]-(points[1]+cp)
    u3=v3/(norm(v3)+.00001)
    v4=points[2]-(points[1]+cp)
    u4=v4/(norm(v4)+.00001)
    theta= 360-arccos(u3@u4)*180/pi if alpha<90 else arccos(u3@u4)*180/pi
    radius=norm(pa-cp)
    return radius



def scl2d(sec,sl):# scale the 2d section 'sec' by a scaling factor 'sl'. this places the scaled section in the bottom center of the original section
    '''
    function to scale a 2d section by an amount "sl" which has to be >0 (keeps the y-coordinates same). 
    e.g.following code scales the section by 0.7 (70% of the original shape)
    sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5)
    sec1=scl2d(sec,.7)
    
    refer file "example of various functions" for application
    '''
    s1=array(translate([0,0,0],sec))
    cp=array(s1).mean(axis=0)
    rev=array(s1).mean(axis=0)+(array(s1)-array(s1).mean(axis=0))*sl
    y1=cp-array([0,array(s1)[:,1].min(),0])
    y2=cp-array([0,rev[:,1].min(),0])
    d=y2-y1
    return c3t2(translate(d,rev))

def scl2d_c(sec,sl):# scale the 2d section 'sec' with scaling factor 'sl'. this places the scaled section in the center of original section or the center of both original and scaled section remains the same.
    '''
    function to scale a 2d section by an amount "sl" which has to be >0 (keeps the revised section in center). 
    e.g.following code scales the section by 0.7 (70% of the original shape)
    sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5)
    sec1=scl2d_c(sec,.7)
    
    refer file "example of various functions" for application
    '''
    s1=array(translate([0,0,0],sec))
    cp=array(s1).mean(axis=0)
    rev=array(s1).mean(axis=0)+(array(s1)-array(s1).mean(axis=0))*sl
    return c3t2(rev)

def scl3d(p,s):# scale 3d prism 'p' with scaling factor 's'. This places the scaled prism at the same bottom of the original prism
    '''
    function to scale a 3d prism keeping the base z-coordinate same. 
    takes 2 arguments "p" to scale and the scaling factor "s". 
    scale factor can take any real number negative values will scale the prism and turn the prism upside down.
    try the following code to understand better:
    sec=circle(10);
    path=cr(pts1([[2,0],[-2,0,2],[0,10,3],[-3,0]]),5)
    sol=prism(sec,path)
    sol1=scl3d(sol,.7)
    refer file "example of various functions" for application
    '''
    p=array(p)
    cp=p.reshape(-1,3).mean(axis=0)
    rev=cp+(p-cp)*s
    z1=p.reshape(-1,3)[:,2].min()
    z2=rev.reshape(-1,3)[:,2].min()
    d=z1-z2
    return translate([0,0,d],rev)

def scl3dc(p,s):# scale a 3d prism 'p' with scaling factor 's'. This places the scaled prism in the center of the original prism or the center of both the prism is same
    '''
     function to scale a 3d prism keeping the prism centered. takes 2 arguments "p" to scale and 
     the scaling factor "s". 
     scale factor can take any real number negative values will scale the prism and turn the prism upside down.
     try the following code to understand better:
     sec=circle(10)
     path=cr(pts1([[2,0],[-2,0,2],[0,10,3],[-3,0]]),5)
     sol=prism(sec,path)
     sol1=scl3dc(p,.7)
     
    refer file "example of various functions" for application
    '''
    p=array(p)
    cp=p.reshape(-1,3).mean(axis=0)
    rev=cp+(p-cp)*s
    return rev.tolist()


    
def m_points1(sec,s,d=.25):# multiple points with in the straight lines in the closed section 'sec'. 's' is the number of segments between each straight line
    '''
    adds 's' number of points in each straight line segment of a section 'sec'
    'd' is the minimum segment length where multipe points to be added
    refer to the file "example of various functions" for application example
    '''
    c=[]
    for i in range(len(sec)):
        i_plus=i+1 if i<len(sec)-1 else 0
        if l_len([sec[i],sec[i_plus]])>=d:
            c.append(ls([sec[i],sec[i_plus]],s))
        else:
            c.append([sec[i],sec[i_plus]])
    return remove_extra_points(concatenate(c))


def m_points1_o(sec,s,d=.25):# multiple points with in the straight lines in the open section 'sec'. 's' is the number of segments between each straight line
    '''
    adds 's' number of points in each straight line segment of an open section 'sec'
    'd' is the minimum segment length where multipe points to be added
    refer to the file "example of various functions" for application example
    '''
    c=[]
    for i in range(len(sec)-1):
        i_plus=i+1 if i<len(sec)-1 else 0
        if l_len([sec[i],sec[i_plus]])>=d:
            c.append(ls([sec[i],sec[i_plus]],s))
        else:
            c.append([sec[i],sec[i_plus]])
    return remove_extra_points(concatenate(c))


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

def sec_clean(sec,sec1,r):
    sec1=array([p for p in sec1 if len(ibsap(sec,p))%2==1])
    p0=sec
    p1=sec[1:]+[sec[0]]
    sec6=swapaxes(array([p0,p1]),0,1)
    p0,p1=array([p0,p1])
    v1=p1-p0
    v2=sec1[:,None]-p0
    v3=sec1[:,None]-p1
    u1=v1/norm(v1,axis=1).reshape(-1,1)
    n=1/norm(v1,axis=1)
    u1.shape,v2.shape
    d=einsum('jk,ijk->ij',u1,v2)
    t=einsum('ij,j->ij',d,n).round(3)
    u1.shape,d.shape
    n1=einsum('jk,ij->ijk',u1,d)
    p1=p0+n1
    sec1.shape,p1.shape
    n2=sec1[:,None]-p1
    n3=sqrt(einsum('ijk,ijk->ij',n2,n2)).round(3)
    n4=where((t>=0)&(t<=1),n3,1e5).min(axis=1)
    m=sec1[(n4>=abs(r)-.02)&(n4<=abs(r)+.02)].tolist()
    return array(m)[cKDTree(m).query(sec)[1]].tolist()


def sec_clean1(sec,sec1,r):
#     sec1=array([p for p in sec1 if len(ibsap(sec,p))%2==1])
    p0=sec
    p1=sec[1:]+[sec[0]]
    sec6=swapaxes(array([p0,p1]),0,1)
    p0,p1=array([p0,p1])
    v1=p1-p0
    v2=sec1[:,None]-p0
    v3=sec1[:,None]-p1
    u1=v1/norm(v1,axis=1).reshape(-1,1)
    n=1/norm(v1,axis=1)
    u1.shape,v2.shape
    d=einsum('jk,ijk->ij',u1,v2)
    t=einsum('ij,j->ij',d,n).round(3)
    u1.shape,d.shape
    n1=einsum('jk,ij->ijk',u1,d)
    p1=p0+n1
    sec1.shape,p1.shape
    n2=sec1[:,None]-p1
    n3=sqrt(einsum('ijk,ijk->ij',n2,n2)).round(3)
    n4=where((t>=0)&(t<=1),n3,1e5).min(axis=1)
    m=sec1[(n4>=abs(r)-.02)&(n4<=abs(r)+.02)].tolist()
    return array(m)[cKDTree(m).query(sec)[1]].tolist()



def fillet_2cir(r1,r2,c1,c2,r,s=50): # fillet between 2 circles with radius 'r1' and 'r2' and center points 'c1' and 'c2' and 'r' is the radius of the fillet
    '''
    function to create 2d fillet between 2 circles, where r1,r2 and c1,c2 are radiuses and enter points of the 2 circles respectively. r-> fillet radius
    example:
    fillet=fillet_2cir(r1=5,r2=3,c1=[0,0],c2=[7,0],r=1)
    
    refer to file "examples of various functions"
   
    '''
    
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

    return arc2+arc1

def filleto_2cir(r1,r2,c1,c2,r,s=50): # fillet between 2 circles with radius 'r1' and 'r2' and center points 'c1' and 'c2' and 'r' is the radius of the fillet. This is an open fillet where first or the second fillet can be called based on requirement
    '''
    function to draw the fillet radius "r" between the 2 circle with radiuses "r1" and "r2" centered at "c1" and "c2" respectively.
    This function gives an additional flexibility for drawing fillet only one side. e.g 
    fillet=filleto_2cir(r1=10,r2=10,c1=[0,0],c2=[20,0],r=10)
    fillet[0] will calculate fillet on one side
    refer to the file "example of various functions" to see the application
    '''
    
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

def tctp(r1,r2,cp1,cp2): # 2 circle tangent points (one side) r1 and r2 are the radius of 2 circles and cp1 and cp2 are the center points
    '''
    function to draw tangent line joining 2 circles with radiuses "r1" and "r2" with center points "cp1" and "cp2" respectively. 
    This function draws tangent line on only one side
     e.g. try this code below:
     sec=tctp(r1=10,r2=5,cp1=[0,0],cp2=[15,6]);
     
     refer to file "example of various functions" for application
 
    '''
    cp1,cp2=array([cp1,cp2])
    v1=cp2-cp1,
    u1=v1/norm(v1)
    ang1=arcsin((r2-r1)/norm(cp2-cp1))*180/pi

    t1=cp1+u1@rm(90+ang1)*r1
    t2=cp2+u1@rm(90+ang1)*r2

    t3=cp1+u1@rm(-90-ang1)*r1
    t4=cp2+u1@rm(-90-ang1)*r2
    return [t1[0].tolist(),t2[0].tolist()]

def tctpf(r1,r2,cp1,cp2): #2 circle tangent point full (both the sides)
    '''
    function to draw tangent line joining 2 circles with radiuses "r1" and "r2" with center points "cp1" and "cp2" respectively. 
    This function draws tangent line on both the sides
    example:
    cir1=circle(10)
    cir2=circle(5,[15,6])
    sec=tctpf(r1=10,r2=5,cp1=[0,0],cp2=[15,6])
    
    refer file "example of various functions" for application
    '''
    cp1,cp2=array([cp1,cp2])
    v1=cp2-cp1,
    u1=v1/norm(v1)
    ang1=arcsin((r2-r1)/norm(cp2-cp1))*180/pi

    t1=cp1+u1@rm(90+ang1)*r1
    t2=cp2+u1@rm(90+ang1)*r2

    t3=cp1+u1@rm(-90-ang1)*r1
    t4=cp2+u1@rm(-90-ang1)*r2
    return [t1[0].tolist(),t2[0].tolist(),t4[0].tolist(),t3[0].tolist()]

def circle(r,cp=[0,0],s=50): # circle with radius r and center point cp, s is the number of segments in the circle
    '''
    function for creating points in circle with radius "r", center point "cp" and number of segments "s" 
    '''
    return array([ [cp[0]+r*cos(i*pi/180),cp[1]+r*sin(i*pi/180)] for i in linspace(0,360,s)][0:-1]).tolist()

def circle_c(r,cp=[0,0],s=50):
    c=array([ [cp[0]+r*cos(i*pi/180),cp[1]+r*sin(i*pi/180)] for i in linspace(0,360,s)][0:-1]).tolist()
    p0,p1=array([c[len(c)-1],c[0]])
    v=p1-p0
    p=(p0+v*.999).tolist()
    return c+[p]

def qmr1(s,r,pl):
    for i in range(len(s)):
        a=[1,0,0] if s[i]=='x' else [0,1,0] if s[i]=='y' else [0,0,1]
        b=r[i]
        pl=[q(a,p,b) for p in pl]
    return pl

def qmr2(s,r,pl):
    for i in range(len(s)):
        a=[1,0,0] if s[i]=='x' else [0,1,0] if s[i]=='y' else [0,0,1]
        b=r[i]
        pl=[[q(a,p1,b) for p1 in p]for p in pl]
    return pl

def q_rot(s,pl):
    '''
    function to rotate a group of points "pl" around a series of axis with defined angles 
    example:
    q_rot(s=["z20","x40","y80"],pl=[[2,0],[10,2]])
    => 
    will rotate the line first around z axis by 20 deg then around x axis by 40 degrees and then around y axis by 80 degrees.
    '''
    if len(array(pl).shape)==2:
        return qmr1([p[0] for p in s],[0 if len(p)==1 else float(p[1:]) for p in s],pl)
    else:
        return qmr2([p[0] for p in s],[0 if len(p)==1 else float(p[1:]) for p in s],pl)
    
def linear_extrude(sec,h=1,a=0,steps=1):
    '''
    function to linear extrude a section where
    sec: section to extrude
    h: height of the extrusion
    a: angle of twist while extruding
    steps: number of steps in each angular extrusion
    refer to the file ' example of various functions' for application example
    '''
    s=2 if a==0 else steps
    return [translate([0,0,h*i if a==0 else h/a*i],q_rot([f"z{0 if a==0 else i}"],sec)) for i in linspace(0,1 if a==0 else a,s)]

def cylinder(r1=1,r2=1,h=1,cp=[0,0],s=50,r=0,d=0,d1=0,d2=0,center=False):
    '''
    function for making a cylinder
    r1 or r: radius of circle at the bottom
    r2 or r: radius of circle at the top
    d1 or d: diameter of circle at the bottom
    d2 or d: diameter of circle at the top
    h: height of the cylinder
    '''
    ra=r if r>0 else d/2 if d>0 else d1/2 if d1>0 else r1
    rb=r if r>0 else d/2 if d>0 else d2/2 if d2>0 else r2
    sec=circle(ra,cp,s)
    
    path=pts([[-ra+.1,0],[ra-.1,0],[rb-ra,h],[-rb+.1,0]])
    p= translate([0,0,-h/2],prism(sec,path)) if center==True else prism(sec,path)
    return p

def square(s=0,center=False):
    m= s if type(s)==int or type(s)==float else s[0]
    n= s if type(s)==int or type(s)==float else s[1]
    sec=cr(pts1([[0,0,.01],[m,0,.01],[0,n,.01],[-m,0,.01]]),10)
    sec1= [[p[0]-m/2,p[1]-n/2] for p in sec] if center==True else sec
    return sec1

def rsz3d(prism,rsz):
    '''
    function to resize a 'prism' to dimensions 'rsz'
    bottom left corner of both the prisms would be same
    refer to file 'example of various functions' for application example
    '''
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
    return translate(t,rev_prism)

def rsz3dc(prism,rsz):
    '''
    function to resize a 'prism' to dimensions 'rsz'
    resized prism will be placed in the center of the original prism or center point of both the prisms will be same
    refer to file 'example of various functions' for application example
    '''
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
    return rev_prism


def bb(prism):
    '''
    function to find the bounding box dimensions of a prism
    refer to the file "example of various functions " for application example
    '''
    prism1=array(prism).reshape(-1,3)
    max_x=prism1[:,0].max()
    max_y=prism1[:,1].max()
    max_z=prism1[:,2].max()
    min_x=prism1[:,0].min()
    min_y=prism1[:,1].min()
    min_z=prism1[:,2].min()
    return [max_x-min_x,max_y-min_y,max_z-min_z]



def cube(s,center=False):
    '''
    function to draw cube with size 's'
    refer to the file "example of various functions " for application example
    
    '''
    if center==False:
        return linear_extrude(square([s[0],s[1]]),s[2])
    elif center==True:
        return translate([0,0,-s[2]/2],linear_extrude(square([s[0],s[1]],True),s[2]))


def sphere(r=0,cp=[0,0,0],s=50):
   '''
   function to draw sphere with radius 'r' , center point 'cp' and number of segments 's'
   refer to the file "example of various functions " for application example
   
   '''
   path=arc(r,-90,90,s=int(s/2))
   p=[ translate([cp[0],cp[1],p[1]+cp[2]],circle(p[0],s=s)) for p in path]
   return array(p).tolist()

# def sphere(r=1,c=[0,0,0],s=20):
#     '''
#     function to draw sphere with radius 'r' , center point 'c' and number of segments 's'
#     refer to the file "example of various functions " for application example
    
#     '''
#     p_l=[]
#     for i in linspace(0,r,int(s/2)+1):
#         i=r*cos(d2r(180/(r+.0001)*(i+.0001)))
#         a=sqrt(r**2-i**2)
#         for j in linspace(a,-a,s)[:-1]:
#             j1=a*cos(d2r(180/(a+.0001)*(j+.0001)))
#             k=sqrt((r**2-i**2-j1**2).round(5))*sign(j)
#             p_l.append([j1+c[0],k+c[1],i+c[2]])
#     return array(p_l).reshape(-1,s-1,3).tolist()

def rsz2d(sec,rsz):
    '''
    function to resize a 2d section to dimensions 'rsz'
    resized section will be placed on bottom center of the original section
    refer the file "example of various functions" for application example
    '''
    avg=array(sec).mean(axis=0)
    max_x=array(sec)[:,0].max()
    min_x=array(sec)[:,0].min()
    max_y=array(sec)[:,1].max()
    min_y=array(sec)[:,1].min()
    r_x=rsz[0]/(max_x-min_x)
    r_y=rsz[1]/(max_y-min_y)
    s=array([ avg+array([r_x*(sec[i][0]-avg[0]),r_y*(sec[i][1]-avg[1])-((min_y-avg[1])*r_y-(min_y-avg[1]))]) for i in range(len(sec))]).round(4)
    return s[sort(unique(s,axis=0,return_index=True)[1])].tolist()
    
def rsz2dc(sec,rsz):
    '''
    function to resize a 2d section to dimensions 'rsz'
    resized section will be placed in center of the original section
    refer the file "example of various functions" for application example
    '''
    avg=array(sec).mean(axis=0)
    max_x=array(sec)[:,0].max()
    min_x=array(sec)[:,0].min()
    max_y=array(sec)[:,1].max()
    min_y=array(sec)[:,1].min()
    r_x=rsz[0]/(max_x-min_x)
    r_y=rsz[1]/(max_y-min_y)
    s=array([ avg+array([r_x*(sec[i][0]-avg[0]),r_y*(sec[i][1]-avg[1])]) for i in range(len(sec))]).round(4)
    return s[sort(unique(s,axis=0,return_index=True)[1])].tolist()

# def ip(prism,prism1):
#     '''
#     function to calculate intersection point between two 3d prisms. 
#      "prism" is the 3d object which is intersected with "prism1".
#      try below code for better understanding:
#     sec=circle(10)
#     path=cr(pts1([[2,0],[-2,0,2],[0,10,3],[-9.9,0]]),5)
#     p=prism(sec,path)
#     p1=cylinder(r=3,h=15,s=30)
#     ip1=ip(p,p1)
    
#     refer to file "example of various functions" for application
#     '''
#     pa=prism
#     pb=prism1
#     p1=array([[ [[pa[i][j],pa[i][j+1],pa[i+1][j]],[pa[i+1][j+1],pa[i+1][j],pa[i][j+1]]] if j<len(pa[i])-1 
#      else [[pa[i][j],pa[i][0],pa[i+1][j]],[pa[i+1][0],pa[i+1][j],pa[i][0]]] 
#      for j in range(len(pa[i]))] 
#               for i in range(len(pa)-1)]).reshape(-1,3,3)
#     p2=array([[[pb[i][j],pb[i+1][j]] for j in range(len(pb[i]))] for i in range(len(pb)-1)]).reshape(-1,2,3)
#     pm=p1[:,0]
#     pn=p1[:,1]
#     po=p1[:,2]
#     px=p2[:,0]
#     py=p2[:,1]
#     v1,v2,v3=py-px,pn-pm,po-pm
#     t1=einsum('ijk,jk->ij',px[:,None]-pm,cross(v2,v3))/einsum('ik,jk->ij',-v1,cross(v2,v3)+[.00001,.00001,.00001])
#     t2=einsum('ijk,ijk->ij',px[:,None]-pm,cross(v3,-v1[:,None]))/einsum('ik,jk->ij',-v1,cross(v2,v3)+[.00001,.00001,.00001])
#     t3=einsum('ijk,ijk->ij',px[:,None]-pm,cross(-v1[:,None],v2))/einsum('ik,jk->ij',-v1,cross(v2,v3)+[.00001,.00001,.00001])
#     p=px[:,None]+einsum('ik,ij->ijk',v1,t1)
#     condition=(t1>=0)&(t1<=1)&(t2>=0)&(t2<=1)&(t3>=0)&(t3<=1)&((t2+t3)>=0)&((t2+t3)<=1)
#     p=p[condition]
# #     p=p[unique(p,return_index=True)[1]]
#     return p.tolist()
    
def ip(prism,prism1,side=-1):
    '''
    function to calculate intersection point between two 3d prisms. 
     "prism" is the 3d object which is intersected with "prism1".
     side: when a ray intersects a solid it can intersect at 2 locations, if the ray is travelling from outside, in that case if '0' is given meaning only the first intersection point is considered, and in case '-1' is given meaning the last intersection point will be considered.
     try below code for better understanding:
    sec=circle(10)
    path=cr(pts1([[2,0],[-2,0,2],[0,10,3],[-9.9,0]]),5)
    p=prism(sec,path)
    p1=cylinder(r=3,h=15,s=30)
    ip1=ip(p,p1)
    
    refer to file "example of various functions" for application
    '''
    pa=prism
    pb=prism1
    p1=array([[ [[pa[i][j],pa[i][j+1],pa[i+1][j]],[pa[i+1][j+1],pa[i+1][j],pa[i][j+1]]] if j<len(pa[i])-1 
     else [[pa[i][j],pa[i][0],pa[i+1][j]],[pa[i+1][0],pa[i+1][j],pa[i][0]]] 
     for j in range(len(pa[i]))] 
              for i in range(len(pa)-1)]).reshape(-1,3,3)
    p2=array([[[pb[i][j],pb[i+1][j]] for j in range(len(pb[i]))] for i in range(len(pb)-1)]).reshape(-1,2,3)
    pm=p1[:,0]
    pn=p1[:,1]
    po=p1[:,2]
    px=p2[:,0]
    py=p2[:,1]
    v1,v2,v3=py-px,pn-pm,po-pm
    t1=einsum('ijk,jk->ij',px[:,None]-pm,cross(v2,v3))/einsum('ik,jk->ij',-v1,cross(v2,v3)+[.00001,.00001,.00001])
    t2=einsum('ijk,ijk->ij',px[:,None]-pm,cross(v3,-v1[:,None]))/einsum('ik,jk->ij',-v1,cross(v2,v3)+[.00001,.00001,.00001])
    t3=einsum('ijk,ijk->ij',px[:,None]-pm,cross(-v1[:,None],v2))/einsum('ik,jk->ij',-v1,cross(v2,v3)+[.00001,.00001,.00001])
    p=px[:,None]+einsum('ik,ij->ijk',v1,t1)
    condition=(t1>=0)&(t1<=1)&(t2>=0)&(t2<=1)&(t3>=0)&(t3<=1)&((t2+t3)>=0)&((t2+t3)<=1)
#     p=p[condition]
#     p=p[unique(p,return_index=True)[1]]
    p=array([[p[i][condition[i]],i%len(pb[0])] for i in range(len(p))],dtype=object)
    n=array([p1[1] for p1 in p if p1[0].tolist()!=[]])
    p=concatenate([concatenate([[[p2,i] for p2 in p1[0]] for p1 in p if (p1[0].tolist()!=[])&(p1[1]==i)],dtype=object) for i in n],dtype=object)
    p=array([p]*len(n),dtype=object)[p[:,1]==unique(p[:,1])[:,None]]
    p=array([[array([p1[0] for p1 in p if p1[1]==i],dtype=object),i] for i in n],dtype=object)
    if side=='all':
        p=concatenate([a[array([l_len([p2[:,0][b],p1]) for p1 in a],dtype=object).argsort()] for (a,b) in p ],dtype=object)
    else:
        p=array([a[array([l_len([p2[:,0][b],p1]) for p1 in a],dtype=object).argsort()[side]] for (a,b) in p if a.tolist()!=[]],dtype=object)
        
    return p.tolist()



def ipf(prism,prism1,r,s,o=0):
    '''
    function to calculate fillet at the intersection point of 2 solids
    'prism': solid 1 or surface 1
    'prism1': solid 2
    'r': radius of the fillet
    's': number of segments in the fillet, more number of segments will give finer finish
    'o': option '0' produces fillet in outer side of the intersection and '1' in the inner side of the intersections
    refer to the file "example of various functions" for application
    '''
    pa=prism
    pb=prism1
    p1=array([[ [[pa[i][j],pa[i][j+1],pa[i+1][j]],[pa[i+1][j+1],pa[i+1][j],pa[i][j+1]]] if j<len(pa[i])-1 
     else [[pa[i][j],pa[i][0],pa[i+1][j]],[pa[i+1][0],pa[i+1][j],pa[i][0]]] 
     for j in range(len(pa[i])-1)] 
              for i in range(len(pa)-1)]).reshape(-1,3,3)
#    segment=seg(pa[0][1:])[:-1]
#    seg1=[[pa[0][0],p1[0],p1[1]] for p1 in segment]
    
#    segment=seg(pa[len(pa)-1][1:])[:-1]
#    seg2=[[pa[len(pa)-1][0],p1[0],p1[1]] for p1 in segment]
#    p1=array(seg1+p1.tolist()+seg2)



    p2=array([[[pb[i][j],pb[i+1][j]] for j in range(len(pb[i]))] for i in range(len(pb)-1)]).reshape(-1,2,3)
    pm=p1[:,0]
    pn=p1[:,1]
    po=p1[:,2]
    px=p2[:,0]
    py=p2[:,1]
    v1,v2,v3=py-px,pn-pm,po-pm
#     px+v1*t1=pm+v2*t2+v3*t3
#     v1*t1-v2*t2-v3*t3=pm-px
    u1=v1/(norm(v1,axis=1).reshape(-1,1)+.0001)
    t1=einsum('ijk,jk->ij',px[:,None]-pm,cross(v2,v3))/(einsum('ik,jk->ij',-v1,cross(v2,v3))+.0001)
    t2=einsum('ijk,ijk->ij',px[:,None]-pm,cross(v3,-v1[:,None]))/(einsum('ik,jk->ij',-v1,cross(v2,v3))+.0001)
    t3=einsum('ijk,ijk->ij',px[:,None]-pm,cross(-v1[:,None],v2))/(einsum('ik,jk->ij',-v1,cross(v2,v3))+.0001)
    p=px[:,None]+einsum('ik,ij->ijk',v1,t1)
    pq=p+(u1*r)[:,None]
    p.shape,pq.shape
    condition=(t1>=0)&(t1<=1)&(t2>=0)&(t2<=1)&(t3>=0)&(t3<=1)&((t2+t3)>=0)&((t2+t3)<=1)
    p=p[condition].tolist()
    pp=p[1:]+[p[0]]
    pq=pq[condition].tolist()
    v4=array(pp)-array(p)
    pnt=array(pq)-array(p)
    n=cross(v4,pnt)
    n=n/(norm(n,axis=1).reshape(-1,1)+.0001)*r
    pnt=n
    cir=[[(p[i]+array(q(v4[i],pnt[i],t))).tolist() for t in linspace(-90,90,10)]for i in range(len(v4))] if o==0 else \
    [[(p[i]+array(q(v4[i],pnt[i],-t))).tolist() for t in linspace(90,270,10)]for i in range(len(v4))] 
    p2=[[ [cir[i][j],cir[i][j+1]] for j in range(len(cir[i])-1)] for i in range(len(cir))]
    p2=array(p2).reshape(-1,2,3)
    px=p2[:,0]
    py=p2[:,1]
    v1,v2,v3=py-px,pn-pm,po-pm
    u1=v1/(norm(v1,axis=1).reshape(-1,1)+.0001)
    t1=einsum('ijk,jk->ij',px[:,None]-pm,cross(v2,v3))/(einsum('ik,jk->ij',-v1,cross(v2,v3))+.0001)
    t2=einsum('ijk,ijk->ij',px[:,None]-pm,cross(v3,-v1[:,None]))/(einsum('ik,jk->ij',-v1,cross(v2,v3))+.0001)
    t3=einsum('ijk,ijk->ij',px[:,None]-pm,cross(-v1[:,None],v2))/(einsum('ik,jk->ij',-v1,cross(v2,v3))+.0001)
    m=px[:,None]+einsum('ik,ij->ijk',v1,t1)
    condition=(t1>=0)&(t1<=1)&(t2>=0)&(t2<=1)&(t3>=0)&(t3<=1)&((t2+t3)>=0)&((t2+t3)<=1)
    m=m[condition]
    m=unique(m,axis=0)[:-1]
    m=m[cKDTree(m).query(p)[1]].tolist()
    p=swapaxes(array([m,p,pq]),0,1)
    p=[[fillet_3p_3d(p2,p1,p0,r_3p_3d([p0,p1,p2])*1.9,s)]for (p0,p1,p2) in p]
    p=array(p).reshape(-1,s+1,3).tolist()
    return p+[p[0]]

def ipf1(p,p1,r,s,o=0):
    pa=[[[[p[i][j],p[i][j+1],p[i+1][j]],[p[i+1][j+1],p[i+1][j],p[i][j+1]]] if j<len(p[0])-1 else \
         [[p[i][j],p[i][0],p[i+1][j]],[p[i+1][0],p[i+1][j],p[i][0]]] \
         for j in range(len(p[0]))] for i in range(len(p)-1)]
    pa=array(pa).reshape(-1,3,3)

    pb=[[[p1[i][j],p1[i+1][j]] for j in range(len(p1[0]))] for i in range(len(p1)-1)]
    pb=array(pb).reshape(-1,2,3)
    
    pm=pa[:,0]
    pn=pa[:,1]
    po=pa[:,2]
    px=pb[:,0]
    py=pb[:,1]
    v1,v2,v3=py-px,pn-pm,po-pm
    v5=array([cross(v2,v3).tolist()]*len(v1))
#     px+v1*t1=pm+v2*t2+v3*t3
#     v1*t1-v2*t2-v3*t3=pm-px
    u1=v1/(norm(v1,axis=1).reshape(-1,1)+.0001)
    t1=einsum('ijk,jk->ij',px[:,None]-pm,cross(v2,v3))/(einsum('ik,jk->ij',-v1,cross(v2,v3))+.0001)
    t2=einsum('ijk,ijk->ij',px[:,None]-pm,cross(v3,-v1[:,None]))/(einsum('ik,jk->ij',-v1,cross(v2,v3))+.0001)
    t3=einsum('ijk,ijk->ij',px[:,None]-pm,cross(-v1[:,None],v2))/(einsum('ik,jk->ij',-v1,cross(v2,v3))+.0001)
    p=px[:,None]+einsum('ik,ij->ijk',v1,t1)
    pq=p+(u1*r)[:,None]
    p.shape,pq.shape
    condition=(t1>=0)&(t1<=1)&(t2>=0)&(t2<=1)&(t3>=0)&(t3<=1)&((t2+t3)>=0)&((t2+t3)<=1)
    p=p[condition].tolist()
    pp=p[1:]+[p[0]]
    pq=pq[condition].tolist()
    v4=array(pp)-array(p)
#     pnt=array(pq)-array(p)
#     n=cross(v4,pnt)
#     n=n/(norm(n,axis=1).reshape(-1,1)+.0001)*r
#     pnt=n
    pnt=v5[condition]
    
    cir=[[(p[i]+array(q(v4[i],pnt[i],t))).tolist() for t in linspace(-90,90,10)]for i in range(len(v4))] if o==0 else \
    [[(p[i]+array(q(v4[i],pnt[i],-t))).tolist() for t in linspace(90,270,10)]for i in range(len(v4))] 
    p2=[[ [cir[i][j],cir[i][j+1]] for j in range(len(cir[i])-1)] for i in range(len(cir))]
    p2=array(p2).reshape(-1,2,3)
    px=p2[:,0]
    py=p2[:,1]
    v1,v2,v3=py-px,pn-pm,po-pm
    u1=v1/(norm(v1,axis=1).reshape(-1,1)+.0001)
    t1=einsum('ijk,jk->ij',px[:,None]-pm,cross(v2,v3))/(einsum('ik,jk->ij',-v1,cross(v2,v3))+.0001)
    t2=einsum('ijk,ijk->ij',px[:,None]-pm,cross(v3,-v1[:,None]))/(einsum('ik,jk->ij',-v1,cross(v2,v3))+.0001)
    t3=einsum('ijk,ijk->ij',px[:,None]-pm,cross(-v1[:,None],v2))/(einsum('ik,jk->ij',-v1,cross(v2,v3))+.0001)
    m=px[:,None]+einsum('ik,ij->ijk',v1,t1)
    condition=(t1>=0)&(t1<=1)&(t2>=0)&(t2<=1)&(t3>=0)&(t3<=1)&((t2+t3)>=0)&((t2+t3)<=1)
    m=m[condition]
    m=unique(m,axis=0)[:-1]
    m=m[cKDTree(m).query(p)[1]].tolist()
    p=swapaxes(array([m,p,pq]),0,1)
    p=[[fillet_3p_3d(p2,p1,p0,r_3p_3d([p0,p1,p2]),s)]for (p0,p1,p2) in p]
    p=array(p).reshape(-1,s+1,3).tolist()
    return p+[p[0]]
    
    
def ipf_co(prism,prism1,r,s,o=0):
    '''
    function to change the orientation of a fillet to create a solid
    prism: solid
    prism1: is another 3d solid which intersects 'prism' to create fillet
    r: radius of the fillet
    s: number of segments in the fillet
    o: options '0' and '1' (refer to the explanation of options in function fillet_sol2sol())
    
    refer to the file "explanation of various functions" for application example
    '''
    a=cpo(ipf(prism,prism1,2,10,0))[1:]
    return a


def ipe(prism,prism1,r,s,o):
    '''
    function to change the orientation of a fillet to create a solid
    prism: solid
    prism1: is another 3d solid which intersects 'prism' to create fillet
    r: radius of the fillet
    s: number of segments in the fillet
    o: options '0' and '1' (refer to the explanation of options in function fillet_sol2sol())
    
    refer to the file "explanation of various functions" for application example
    '''
    a=cpo(ipf(prism,prism1,2,10,0))[1:]
    return a


def s_int1(s): #creates intersection between all the segments of a section which are crossing
    '''
    calulates the self intersection points of a list of line segments 's'
    it picks the intersection points only if the 2 lines are crossing each other
    e.g.
    sec=seg([[0,0],[10,0],[15,7]])
    s_int1(sec) => []
    
    sec=offset_segv([[0,0],[10,0],[15,7]],-1)
    s_int1(sec) => 
    [[9.485, 1.0],
     [4.508, 1.0],
     [9.485266528793264, 0.9998381937190964],
     [11.974266528793265, 4.484438193719097],
     [4.507385465331124, 0.9999168600047348],
     [11.974385465331125, 4.484516860004734]]
     
    refer to file 'example of various functions' for application example
    '''
    p0=array([array(s)[:,0]]*len(s)).transpose(1,0,2)
    p1=array([array(s)[:,1]]*len(s)).transpose(1,0,2)
    v1=p1-p0
    p2=array([array(s)[:,0]]*len(s))
    p3=array([array(s)[:,1]]*len(s))
    v2=p3-p2
    v1.shape,v2.shape
    A=inv(array([v1+[.00001,0],-v2+[.00001,.00001]]).transpose(1,0,2,3).transpose(0,2,1,3).transpose(0,1,3,2))
    B=p2-p0
    t=einsum('ijkl,ijl->ijk',A,B)[:,:,0].round(4)
    u=einsum('ijkl,ijl->ijk',A,B)[:,:,1].round(4)
    condition=(t>0)&(t<1)&(u>0)&(u<1)
    d=(p0+einsum('ijk,ij->ijk',v1,t))[condition].tolist()
    return remove_extra_points(d)

def comb(n,i): 
    '''
    calculates number of possible combinations for "n" items with "i" selected items
    comb(8,2) => 28
    '''
    return int(math.factorial(n)/(math.factorial(i)*math.factorial(n-i)))

def bezier(p,s=10):
    '''
    bezier curve defined by points 'p' and number of segments 's'
    refer file "example of various functions" for application
    '''
    return array([array([ comb((len(p)-1),i)*(1-t)**((len(p)-1)-i)*t**i*array(p[i])  for i in range(len(p))]).sum(0) for t in linspace(0,1,s)]).tolist()

def arc_3d(v=[0,0,1],r=1,theta1=0,theta2=360,cw=-1,s=50):
    '''
    3d arc defined by normal vector 'v', radius 'r1', start angle 'theta1', 
    end angle 'theta2' , clockwise(1) or counter clockwise(-1) and number of segments 's'
    
    refer file "example of various functions" for application example
    '''

    if uv(v)==[0,0,1]:
        arc1=arc(r,theta1,theta2,[0,0],s) if cw==-1 else flip(arc(r,theta1,theta2,[0,0],s))
        return c2t3(arc1)
    elif uv(v)==[0,0,-1]:
        arc1=arc(r,theta1,theta2,[0,0],s) if cw==-1 else flip(arc(r,theta1,theta2,[0,0],s))
        arc1=q_rot(['y180'],arc1)
        return arc1
    else:
        sec=arc(r,theta1,theta2,[0,0],s) if cw==-1 else flip(arc(r,theta1,theta2,[0,0],s))
        s=q_rot(['x90','z-90'],sec)
        v1=array(v)+array([0,0,0.00001])
        va=[v1[0],v1[1],0]
        u1=array(uv(v1))
        ua=array(uv(va))
        v2=cross(va,v1)
        a1=arccos(u1@ua)*180/pi
        a2=ang(v1[0],v1[1])
        s1=q_rot([f'z{a2}'],s)
        sec1=[q(v2,p,a1) for p in s1]
        return sec1

def plane(nv,radius):
    '''
    plane defined by normal 'nv' and 'radius'
    
    refer file "example of various functions" for application example
    '''
    sec1=arc_3d(nv,.0001,0,360,-1)
    sec2=arc_3d(nv,radius,0,360,-1)
    plane=[sec1,sec2]
    return plane

def l_cir_ip(line,cir):
    '''
    line circle intersection point
    '''
    
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


def s_pnt(pnt): # starting point for calculating convex hull (bottom left point)
    '''
    starting point for calculating convex hull (bottom left point)
    '''
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
    '''
    function to calculate convex hull for a list of points 'pnt'
    
    refer file "example of various functions" for application example
    '''
    c=[]
    np=n_pnt(pnt,s_pnt(pnt),0)
    for i in range(len(pnt)):
        c.append(np[0])
        np=n_pnt(pnt,np[0],np[1])
        if np[0]==s_pnt(pnt):
            break
    return [s_pnt(pnt)]+c

def convex(sec):
    '''
    function to check whether a section is convex or not
    example:
    sec1=cr_c(pts1([[0,0,.2],[8,3,3],[5,7,1],[-8,0,2],[-5,20,1]]),20)
    sec2=cr_c(pts1([[0,0,.1],[7,5,2],[5,7,3],[-5,7,5],[-7,5,5]]),20)
    convex(sec1),convex(sec2) => (False, True)
    
    refer file "example of various functions" for application example
    '''

    return (array(cwv(c3t2(sec)))==-1).all()|(array(cwv(c3t2(sec)))==1).all()

def oo_convex(sec,r): #outer offset of a convex section
    s=flip(sec) if cw(sec)==1 else sec
    return offset_points(sec,r)

def cir_p_t(cir,p):
    '''
    circle to point tangent line (point should be outside the circle)
    refer file "example of various functions" for application example
    '''
    cp=cp_3p(cir[0],cir[int(len(cir)/3)],cir[int(len(cir)*2/3)])
    r=r_3p([cir[0],cir[int(len(cir)/3)],cir[int(len(cir)*2/3)]])
    l1=l_len([p,cp])
    v1=array(p)-array(cp)
    theta1=ang(v1[0],v1[1])
    theta2=arccos(r/l1)*180/pi
    theta3=(theta1-theta2)*pi/180
    tp=array([r*cos(theta3),r*sin(theta3)])+array(cp)
    tp=tp.tolist()
    return tp

def p_cir_t(p,cir): # point to circle tangent line (point should be outside the circle)
    '''
    point to circle tangent line (point should be outside the circle)
    refer file "example of various functions" for application example
    '''
    cp=cp_3p(cir[0],cir[int(len(cir)/3)],cir[int(len(cir)*2/3)])
    r=r_3p([cir[0],cir[int(len(cir)/3)],cir[int(len(cir)*2/3)]])
    l1=l_len([p,cp])
    v1=array(p)-array(cp)
    theta1=ang(v1[0],v1[1])
    theta2=arccos(r/l1)*180/pi
    theta3=(theta1+theta2)*pi/180
    tp=array([r*cos(theta3),r*sin(theta3)])+array(cp)
    tp=tp.tolist()
    return tp


def v_sec_extrude(sec,path,o):
    '''
    extrude a section 'sec' through a path 'path' 
    section will vary from start to end such that at the end the section will be offset by 'o' distance
    refer to the file "example of various functions" for application example
    '''
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

def t_cir_tarc(r1,r2,cp1,cp2,r,side=0,s=50): #two circle tangent arc
    '''
    function draws a arc which is tangent to 2 circles defined by radiuses 'r1' and 'r2' and center points 'cp1' and 'cp2'
    's' is the number of segments of the tangent arc
    'r' is the radius of the tangent arc (it should be >= (r1+r2+center distance of 2 circles)/2)
    'side' there are 2 sides of the circles where the arc could be created defined by '0' and '1'
    refer the file "example of various functions " for application examples
    '''
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
        if side==0:
            cp3=p0-(u1@rm(90))*h
        elif side==1:
            cp3=p0+(u1@rm(90))*h

        v2=cp2-cp3
        u2=v2/norm(v2)
        v3=cp1-cp3
        u3=v3/norm(v3)
        p1=cp2+u2*r2
        p2=cp1+u3*r1

        if side==0:
            arc1=arc_2p(p1,p2,r,-1,s=s)
        elif side==1:
            arc1=arc_2p(p2,p1,r,-1,s=s)



    else:
        if side==0:
            arc1=filleto_2cir(r1,r2,cp1,cp2,r,s=s)[1]
        elif side==1:
            arc1=filleto_2cir(r1,r2,cp1,cp2,r,s=s)[0]
            
    return arc1
    
    
def tcct(r1,r2,cp1,cp2,cw=-1): # two circle cross tangent
    '''
    function to draw cross tangent between 2 circles
    refer to the file "example of various functions " for application examples
    '''
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

def arc_3p(p1,p2,p3,s=30):
    ''' 
    function to draw arc with 3 known points 'p1','p2','p3' 
    's' is the number of segments of the arc
    refer to the file "example of various functions " for application examples
    
    
    '''
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

def cir_3p(p1,p2,p3,s=30):
    ''' 
    function to draw circle with 3 known points 'p1','p2','p3' 
    's' is the number of segments of the circle
    refer to the file "example of various functions " for application examples
    
    
    '''
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
#     v3=p1-cp
#     v4=p2-cp
#     v5=p3-cp
#     a1=ang(v3[0],v3[1])
#     a2=ang(v4[0],v4[1])
#     a3=ang(v5[0],v5[1])
#     a4=(a3+360 if a3<a1 else a3) if cw([p1,p2,p3])==-1 else (a3 if a3<a1 else a3-360)
    return circle(r,cp,s)

def cp_3p(p1,p2,p3):
    '''
    function to calculate center point of a circle created from 3 known points 'p1','p2','p3'
    refer to the file "example of various functions " for application examples
    
    
    '''
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



def ip_surf(surf2,surf1):
    '''
     function to calculate intersection point between two 3d prisms or between surface and solid. 
     "surf2" is the 3d object which is intersected with "surf1".
 try below code for better understanding:
 sec=circle(10);
 path=cr(pts1([[2,0],[-2,0,2],[0,10,3],[-9.9,0]]),5);
 prism=prism(sec,path);
 prism1=q_rot(["y40"],cylinder(r=3,h=15,s=30));
 %swp(prism);
 %swp(prism1);
 ip=ip_surf(prism,prism1);
 points(ip,.2);
    '''
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





def perp_points(line,points):
    p0=line[0]
    p1=line[1]
    p0,p1=array([p0,p1])
    v1=p1-p0
    u1=v1/(norm(v1)+.00001)
    v2=array(points)-p0
    v2cost=einsum('j,ij->i',u1,v2)
    cond=(v2cost>=0)&(v2cost<=l_len(line))
    pnts=array(points)[cond]
    v2=pnts-p0
    
    return pnts.tolist()

def perp_points_with_dist(line,points,f=1):
    '''
    function returns all the points which are projected on the line and has a perpendicular distance of line_length/f.
    
    
    '''

    p0=line[0]
    p1=line[1]
    p0,p1=array([p0,p1])
    v1=p1-p0
    u1=v1/(norm(v1)+.00001)
    v2=array(points)-p0
    v2cost=einsum('j,ij->i',u1,v2)
    cond=(v2cost>=0)&(v2cost<=l_len(line))
    pnts=array(points)[cond]
    v2=pnts-p0
    d=cross(v1,v2)/norm(v1)
    cond=d<=l_len(line)/f
    pnts=pnts[cond].tolist()
    return pnts

def perp_points_within_d(line,points,l=1):
    '''
    function returns all the points which are projected on the line and has a perpendicular distance of < l.
    '''

    p0=line[0]
    p1=line[1]
    p0,p1=array([p0,p1])
    v1=p1-p0
    u1=v1/(norm(v1)+.00001)
    v2=array(points)-p0
    v2cost=einsum('j,ij->i',u1,v2)
    cond=(v2cost>=0)&(v2cost<=l_len(line))
    pnts=array(points)[cond]
    v2=pnts-p0
    d=cross(v1,v2)/norm(v1)
    cond=d<l
    pnts=pnts[cond].tolist()
    return pnts


    
def perp_dist(line,point):
    p0=line[0]
    p1=line[1]
    p0,p1=array([p0,p1])
    v1=p1-p0
    u1=v1/(norm(v1)+.00001)
    v2=array(point)-p0
    n1=cross(v1,v2)
    d=n1/norm(v1)
    return d




def sq(d,cp=[0,0]):
    cp=array(cp)-d/2
    cp=[cp[0],cp[1],0]
    return c3t2(translate(cp,[[0,0],[d,0],[d,d],[0,d]]))

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

# def exclude_points(points,pnts):
#     return [p for p in points if p not in pnts]

def exclude_points(list1,list_to_exclude):
    list1,list_to_exclude=array(list1),array(list_to_exclude)
    return list1[~(list_to_exclude==list1[:,None]).all(2).any(1)].tolist()

def i_p2dw(l1,l):
    p0,p1=array(l1)
    p2,p3=array(l)
    v1=p1-p0
    v2=p3-p2
#                     p0+v1*t1=p2+v2*t2
#                     v1*t1-v2*t2=p2-p0
    im=inv(array([v1,-v2]).transpose(1,0)
                  +array([[.000001,.000002],[.000002,.000003]]))
    t=(im@(p2-p0))[0]
    u=(im@(p2-p0))[1]
    return  (p0+v1*t).tolist() if (0<t<1)& (0<u<1) else []



        
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
    

def pies1(sec,pnts):
    '''
    function to find points 'pnts' which are inside an enclosed section 'sec'
    refer to the file "example of various functions " for application examples
    
    
    '''
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
        p=p2-p0[:,None]
        im=pinv(array([v1,-v2]).transpose(1,0,2).transpose(0,2,1))
        im=array([im]*len(p0))
        t=einsum('ijkl,ijl->ijk',im,p)

        s10=[p0[i] for i in range(len(p0)) if \
                t[i][(t[i][:,0]>=0)&(t[i][:,1]>=0)&(t[i][:,1]<=1)].shape[0]%2 \
             ==1]
        return array(s10).tolist()


def rsec(line,radius):
    p0=line[0]
    p1=line[1]
    p0,p1=array([p0,p1])
    v=p1-p0
    a=ang(v[0],v[1])
    return arc(radius,a+90,a+270,p0,int(round(10+log10(radius+1)**6,0)))+arc(radius,a-90,a+90,p1,int(round(10+log10(radius+1)**6,0)))



def cleaning_seg(sec):
    r=-max_r(sec)-1
    s=seg(sec)
    s1=offset_points(sec,r)
    s2=seg(s1)
    u=array([(array(p[1])-array(p[0]))/norm(array(p[1])-array(p[0])) for p in s])
    u1=array([(array(p[1])-array(p[0]))/norm(array(p[1])-array(p[0])) for p in s2])
    s3=array(s)[norm(u-u1,axis=1)<1].tolist()
    return s3

def cleaning_sec_inner(sec,r):
    s=cleaning_seg(sec)
    s1=[rsec(p,abs(r)-.01) for p in s]
    return s1

def cleaning_sec_outer(sec,r):
    s=cleaning_seg(sec)
    s1=[rsec(p,abs(r)-.1) for p in s]
    return s1



#def r_sec(r1,r2,cp1,cp2):
#    '''
#    rounded section around a line with points 'cp1' and  'cp2' with radiuses 'r1' and 'r2' respectively
#    '''
#    l=tctpf(r1,r2,cp1,cp2)
#    l=l[:2]+arc_2p(l[1],l[2],r1)+l[2:]+arc_2p(l[3],l[0],r2)
#    return l

def r_sec(r1,r2,cp1,cp2):
    '''
    creates a rounded section around a line defined by points 'cp1' and 'cp2'
    radius around 'cp1' is 'r1' and radius around 'cp2' is 'r2'
    
    '''
    sec=tctpf(r1,r2,cp1,cp2)
    a1=arc_2p(sec[3],sec[0],r1,1) if r1<r2 else arc_long_2p(sec[3],sec[0],r1,1)
    a2=arc_long_2p(sec[1],sec[2],r2,1) if r1<r2 else arc_2p(sec[1],sec[2],r2,1)
    sec1=a2+a1
    return sec1



def cs(sec,d):
    '''
    please use function cs1 instead of this
    cleaning section for removing excess points from offset
    refer to the file "example of various functions " for application examples
    
    '''
    sec=c3t2(q_rot(['z.0001'],sec))
    r=abs(d)
    a=array(sec)[array(list_r(sec))==0]
    a=seg(a)
    p1=array([a[i] for i in range(len(a)) if i%2!=0]).tolist()
    cs=[r_sec(r,r,p2[0],p2[1]) for p2 in p1]
    return cs

def cs1(sec,d):
    '''
    creates a cleaning section for removing excess points for offseting a section 'sec' with offset distance 'd'
    
    '''
    r=abs(d)
    a=seg(sec)
    cs=[r_sec(r-r/1000,r-r/1000,p2[0],p2[1]) for p2 in a]
    return cs
    
def cs2(sec,d):
    r=abs(d)
    a=seg(sec)
    cs=[r_sec(r-r/1000,r-r/1000,p2[0],p2[1]) for p2 in a if l_len(p2)>.5]
    return cs




def swp(bead2):
    '''
    function to render various 3d shapes
    example:
    swp(cylinde(d=10,h=20)) will render a cylinder with dia 10 and height 20
    refer to the file "example of various functions " for application examples
    
    '''
    n1=arange(len(bead2[0])).tolist()
    n2=array([[[[(j+1)+i*len(bead2[0]),j+i*len(bead2[0]),j+(i+1)*len(bead2[0])],[(j+1)+i*len(bead2[0]),j+(i+1)*len(bead2[0]),(j+1)+(i+1)*len(bead2[0])]] \
             if j<len(bead2[0])-1 else \
             [[0+i*len(bead2[0]),j+i*len(bead2[0]),j+(i+1)*len(bead2[0])],[0+i*len(bead2[0]),j+(i+1)*len(bead2[0]),0+(i+1)*len(bead2[0])]] \
                 for j in range(len(bead2[0]))] for i in range(len(bead2)-1)]).reshape(-1,3).tolist()
    n3=(array(flip(arange(len(bead2[0]))))+(len(bead2)-1)*len(bead2[0])).tolist()
    n=[n1]+n2+[n3]
    pnt=array(bead2).reshape(-1,3).round(4).tolist()
    return f'polyhedron({pnt},{n},convexity=10);'
    
def swp1(bead2):
    '''
    function to render various 3d shapes with first and the last section triangulated. only works with convex sections
    example:
    swp(cylinde(d=10,h=20)) will render a cylinder with dia 10 and height 20
    refer to the file "example of various functions " for application examples
    
    '''
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
    n=n1+n2+n3
    pnt=[cp1]+array(bead2).reshape(-1,3).round(4).tolist()+[cp2]
    return f'polyhedron({pnt},{n},convexity=10);'
    

    

def swp_c(bead2):
    '''
    function to render various polyhedron with closed loop shapes e.g. fillets
    refer to the file "example of various functions " for application examples
    
    '''
    n1=arange(len(bead2[0])).tolist()
    n2=array([[[[(j+1)+i*len(bead2[0]),j+i*len(bead2[0]),j+(i+1)*len(bead2[0])],[(j+1)+i*len(bead2[0]),j+(i+1)*len(bead2[0]),(j+1)+(i+1)*len(bead2[0])]] \
             if j<len(bead2[0])-1 else \
             [[0+i*len(bead2[0]),j+i*len(bead2[0]),j+(i+1)*len(bead2[0])],[0+i*len(bead2[0]),j+(i+1)*len(bead2[0]),0+(i+1)*len(bead2[0])]] \
                 for j in range(len(bead2[0]))] for i in range(len(bead2)-1)]).reshape(-1,3).tolist()
    n3=(array(flip(arange(len(bead2[0]))))+(len(bead2)-1)*len(bead2[0])).tolist()
    n=[n1]+n2+[n3]
    pnt=array(bead2).reshape(-1,3).round(4).tolist()
    return f'polyhedron({pnt},{n2},convexity=10);'
    
def mesh_vf(bead2):
    '''
    function to render various polyhedron with closed loop shapes e.g. fillets
    refer to the file "example of various functions " for application examples
    
    '''
    n1=arange(len(bead2[0])).tolist()
    n2=array([[[[(j+1)+i*len(bead2[0]),j+i*len(bead2[0]),j+(i+1)*len(bead2[0])],[(j+1)+i*len(bead2[0]),j+(i+1)*len(bead2[0]),(j+1)+(i+1)*len(bead2[0])]] \
             if j<len(bead2[0])-1 else \
             [[0+i*len(bead2[0]),j+i*len(bead2[0]),j+(i+1)*len(bead2[0])],[0+i*len(bead2[0]),j+(i+1)*len(bead2[0]),0+(i+1)*len(bead2[0])]] \
                 for j in range(len(bead2[0]))] for i in range(len(bead2)-1)]).reshape(-1,3)
    pnt=array(bead2).reshape(-1,3).round(4)
    return [pnt,n2]

def resurf(surf,f):
    base=array(c3t2(surf)).reshape(-1,2).tolist()
    c=[]
    for i in range(len(surf)):
        if len(base)<=2:
            break
        else:
            base1=concave_hull(base,f)
            base=exclude_points(base,base1)
            c.append(base1)
    base=concave_hull(array(c3t2(surf)).reshape(-1,2).tolist(),f)
    c=[array(p)[cKDTree(p).query(base)[1]].tolist() for p in c]
    base=array(c3t2(surf)).reshape(-1,2).tolist()
    surf=array(surf).reshape(-1,3)
    c= [surf[cKDTree(base).query(p)[1]].tolist() for p in c]
    return c

def surf_extrudef(surf,t=-.05):
    '''
    surface with a polyline 2d sketch and a 3d path. thickness of the surface can be set with parameter "t". 
    positive and negative value creates thickness towards +z and -z directions respectively
    refer file "example of various functions"
    '''
    s=cpo(surf)
    s1=translate([0,0,t],[flip(p) for p in s])
    s2=array([s,s1]).transpose(1,0,2,3)
    
    i,j,k,l=s2.shape
    s2=s2.reshape(i,j*k,l).tolist()
    return s2 if t>0 else flip(s2)



def swp_prism_h(prism_big,prism_small):
    '''
    
    creats a hollow prism with 2 similar prisms (1 big and 1 smaller)
    
    refer the file "example of various functions" for application example
    '''
    
    p1=prism_big
    p2=flip(prism_small)
    p3=p1+p2+[p1[0]]
    return p3
    

    
def pmdp(line,pnts,f=1): #perpendicular minimum distance point
    '''
    function returns a points which projects on the line and has minimum perpendicular distance from a list of points with distances <= line_length/f

    '''
    pnts=remove_extra_points(pnts)
    pnts=perp_points_with_dist(line,pnts,f=f)
    if pnts==[]:
        return line
    else:
        a=array([perp_dist(line,p) for p in pnts])
        if a.tolist()!=[] and (a>0).any():
            b=array(pnts)[a>0][a[a>0].argmin()].tolist()
            return [line[0],b,line[1]]
        else:
            return line


def surf_base(surf,h=0):
    '''
    creates a solid from any surface, 'h' is the height of the base of the surface
    refer the file "example of various functions" for application example
    
    '''
    s=cpo(surf)
    s1=translate([0,0,h],c2t3(c3t2([flip(p) for p in s])))
    s2=array([s,s1]).transpose(1,0,2,3)
    
    i,j,k,l=s2.shape
    s2=s2.reshape(i,j*k,l).tolist()
    t=array(surf).reshape(-1,3).mean(0)[2]
    return s2 if h>t else flip(s2)

def cr_3d(p,s=5): # Corner radius 3d where 'p' are the list of points (turtle movement) and 's' is number of segments for each arc
    pnts=array(p)[:,0:3]
    pnts=pnts.cumsum(0)

    rds=array(p)[:,3]
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
    return remove_extra_points(c) 



def helix(radius=10,pitch=10, number_of_coils=1, step_angle=1):
    '''
    creates a helix with radius, pitch and number of coils parameters
    
    refer to file "example of various functions" for application example
    
    '''
    return [[radius*cos(i*pi/180),radius*sin(i*pi/180),i/360*pitch] for i in arange(0,360*number_of_coils,step_angle)]


def surf_offset(sec,d):
    '''
    function to offset the surface 'sec' by a distance 'd'
    refer to file "example of various functions" for application example
    '''
    sec1=[]
    for i in range(len(sec)):
        for j in range(len(sec[0])):
            if j<len(sec[0])-1 and i<len(sec)-1:
                p0=sec[i][j]
                p1=sec[i][j+1]
                p2=sec[i+1][j]
            elif j==len(sec[0])-1 and i<len(sec)-1:
                p0=sec[i][j]
                p1=sec[i][0]
                p2=sec[i+1][0]
            elif j<len(sec[0])-1 and i==len(sec)-1:
                p0=sec[i][j]
                p1=sec[i-1][j]
                p2=sec[i-1][j+1]
            elif j==len(sec[0])-1 and i==len(sec)-1:
                p0=sec[i][j]
                p1=sec[i-1][j]
                p2=sec[i-1][0]

            p0,p1,p2=array([p0,p1,p2])
            v1,v2=p1-p0,p2-p0
            v3=cross(v1,v2)
            nv3=v3/norm(v3)
            p3=(p0+nv3*d).tolist()
            sec1.append(p3)
    a,b,c=array(sec).shape
    sec1=array(sec1).reshape(a,b,c).tolist()
    return sec1
    
    
def path_to_vectors(path):
    c=[]
    for i in range(len(path)):
        i_plus=i+1 if i<len(path)-1 else 0
        p0=path[i]
        p1=path[i_plus]
        p0,p1=array([p0,p1])
        v1=p1-p0
        c.append(v1.round(4).tolist())
    return vector_correct(c)

def vector_correct(c):
    for i in range(len(c)):
        if i>0 and i<len(c)-1:
            if c[i][0]==0 and abs(c[i-1][0])>0 and abs(c[i+1][0])>0:
                c[i][0]=.001
            if c[i][1]==0 and abs(c[i-1][1])>0 and abs(c[i+1][1])>0:
                c[i][1]=.001
            if c[i][2]==0 and abs(c[i-1][2])>0 and abs(c[i+1][2])>0:
                c[i][2]=.001
    path=c
    return path




def path_extrude(sec,path):
    '''
    function to extrude a section 'sec' along a open path 'path'
    refer to file "example of various functions" for application example
    '''
    p1=path[:-1]
    p2=path[1:]
    p1,p2=array([p1,p2])
    v1=p2-p1
    u1=v1/norm(v1,axis=1).reshape(-1,1)
    v2=concatenate([[u1[0]],(u1[1:]+u1[:-1])/2,[u1[-1]]])
    sec2=[]
    for i in range(len(path)):
        sec1=translate(path[i],sec2vector(v2[i],sec))
        sec2.append(sec1)
    return sec2


def path_extrudec(sec,path):
    '''
    function to extrude a section 'sec' along a closed loop path 'path'
    refer to file "example of various functions" for application example
    '''
    p1=path
    p2=path[1:]+[path[0]]
    p1,p2=array([p1,p2])
    v1=p2-p1
    u1=v1/norm(v1,axis=1).reshape(-1,1)
    v2=concatenate([[(u1[-1]+u1[0])/2], (u1[1:]+u1[:-1])/2])
    sec2=[]
    for i in range(len(path)):
        sec1=translate(path[i],sec2vector(v2[i],sec))
        sec2.append(sec1)
    sec2=sec2+[sec2[0]]
    # sec3=concatenate([align_sec(sec2[i-1],sec2[i]) for i in range(1,len(sec2))]).tolist()
    return sec2


def multiple_sec_extrude(path_points=[],radiuses_list=[],sections_list=[],option=0,s=10):
    '''
    explanation of the function 'multiple_sec_extrude'
    path_points: are the points at which sections needs to be placed,
    radiuses: radius required at each path_point. this can be '0' in case no radius required in the path
    sections_list= list of sections required at each path_points. same section can be provided for various path_points as well
    option: can be '0' in case the number of points in each section do not match or '1' in case number of points for each section are same
    s: in case value of radiuses is provided 's' is the number of segments in that path curve
    
    refer to file "example of various functions" for application example
    '''
    p=array(path_points)
    r=radiuses_list
    if option==0:
        sections=[sections_list[0]]+[sort_pointsv(sections_list[0],p) for p in sections_list[1:]]
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
            s2=q_rot(['x90','z-90',f'z{a2}'],sections[i])
            s3=translate(p0,flip([q(v2,p,-a1) for p in s2]))
            s4.append(s3)
    return s4

def pntsnfaces(bead2):
    '''
    function returns points and faces of a prism
    refer file "example of various functions" for application example
    '''
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
    '''
    function to offset a 'path' by 'd' distance
    example:
    line=[[0,0],[10,0]]
    path_offset(line,-3) => [[0,3],[10,3]]
    
    refer file "example of various functions" for application example
    '''
    p=array([offset_l(p,d) for p in seg(path)[:-1]])
    return p[:,0].tolist()+[p[len(p)-1][1].tolist()]


def fillet_sol2sol(p=[],p1=[],r=1,s=10,o=0,f=1.8):
    ''' 
    function to calculate fillet at the intersection point of 2 solids
    'p': solid 1
    'p1': solid 2
    'r': radius of the fillet
    's': number of segments in the fillet, more number of segments will give finer finish
    'o': option '0' produces fillet in outer side of the intersection and '1' in the inner side of the intersections
    refer file "example of various functions" for application example
    '''
    pa=[[[[p[i][j],p[i][j+1],p[i+1][j]],[p[i+1][j+1],p[i+1][j],p[i][j+1]]] if j<len(p[0])-1 else \
         [[p[i][j],p[i][0],p[i+1][j]],[p[i+1][0],p[i+1][j],p[i][0]]] \
         for j in range(len(p[0]))] for i in range(len(p)-1)]
    pa=array(pa).reshape(-1,3,3)
#    segment=seg(p[0][1:])[:-1]
#    seg1=[[p[0][0],p1[0],p1[1]] for p1 in segment]
    
#    segment=seg(p[len(p)-1][1:])[:-1]
#    seg2=[[p[len(p)-1][0],p1[0],p1[1]] for p1 in segment]
#    pa=array(seg1+pa.tolist()+seg2)

#     pb=[[[p1[i][j],p1[i+1][j]] for j in range(len(p1[0]))] for i in range(len(p1)-1)]
#     pb=array(pb).reshape(-1,2,3)
    p2=cpo(p1)
    pb=[[[p2[i][j],p2[i][j+1]] for j in range(len(p2[0])-1)] for i in range(len(p2))]
    pb=array(pb).reshape(-1,2,3)
    
    p01,p02,p03,p04,p05=pa[:,0],pa[:,1],pa[:,2],pb[:,0],pb[:,1]

    v1,v2,v3=p05-p04,p02-p01,p03-p01
    i,j=len(v1),len(v2)
#     array([(-v1).tolist()]*j).transpose(1,0,2).shape,array([v2]*i).shape,array([cross(v2,v3)]*i).shape,(p04[:,None]-p01).shape
#     cross(v3,-v1[:,None]).shape,cross(-v1[:,None],v2).shape
    a=einsum('ijk,ijk->ij',array([cross(v2,v3)]*i),p04[:,None]-p01)
    b=(1/einsum('ijk,ijk->ij',array([-v1]*j).transpose(1,0,2),array([cross(v2,v3)+.00001]*i)))
    t1=einsum('ij,ij->ij',a,b)

    a=einsum('ijk,ijk->ij',cross(v3,-v1[:,None]),p04[:,None]-p01)
    b=(1/einsum('ijk,ijk->ij',array([-v1]*j).transpose(1,0,2),array([cross(v2,v3)+.00001]*i)))
    t2=einsum('ij,ij->ij',a,b)

    a=einsum('ijk,ijk->ij',cross(-v1[:,None],v2),p04[:,None]-p01)
    b=(1/einsum('ijk,ijk->ij',array([-v1]*j).transpose(1,0,2),array([cross(v2,v3)+.00001]*i)))
    t3=einsum('ij,ij->ij',a,b)

    condition=(t1>=0) & (t1<=1) & (t2>=0) & (t2<=1) & (t3>=0) & (t3<=1) & (t2+t3>=0) & (t2+t3<=1)
#     condition.shape

    pnt1=(p04[:,None]+einsum('ijk,ij->ijk',array([v1]*j).transpose(1,0,2),t1))[condition]
#     pnt1.shape

    uv1=v1/norm(v1,axis=1).reshape(-1,1)
    uv1=array([uv1]*j).transpose(1,0,2)[condition]
#     uv1.shape
#     pnt2=pnt1+uv1*r

    a=cross(v2,v3)
    b=a/(norm(a,axis=1).reshape(-1,1)+.00001)
    b=array([b]*i)[condition]
#     b.shape

    nxt_pnt=array(pnt1[1:].tolist()+[pnt1[0]])
    v_rot=nxt_pnt-pnt1

    if o==0:
        cir=array([[pnt1[i]+array(q(v_rot[i],b[i]*r,t)) for t in linspace(0,180,5)] for i in arange(len(pnt1))]).tolist()
    else:
        cir=array([[pnt1[i]+array(q(v_rot[i],b[i]*r,-t)) for t in linspace(0,180,5)] for i in arange(len(pnt1))]).tolist()

    pc=array([[[cir[i][j],cir[i][j+1]]  for j in arange(len(cir[0])-1)] for i in arange(len(cir))]).reshape(-1,2,3)
#     pc.shape

    p01,p02,p03,p04,p05=pa[:,0],pa[:,1],pa[:,2],pc[:,0],pc[:,1]

    v1,v2,v3=p05-p04,p02-p01,p03-p01
    i,j=len(v1),len(v2)

    a1=einsum('ijk,ijk->ij',array([cross(v2,v3)]*i),p04[:,None]-p01)
    b1=(1/einsum('ijk,ijk->ij',array([-v1]*j).transpose(1,0,2),array([cross(v2,v3)+.00001]*i)))
    t1=einsum('ij,ij->ij',a1,b1)

    a1=einsum('ijk,ijk->ij',cross(v3,-v1[:,None]),p04[:,None]-p01)
    b1=(1/einsum('ijk,ijk->ij',array([-v1]*j).transpose(1,0,2),array([cross(v2,v3)+.00001]*i)))
    t2=einsum('ij,ij->ij',a1,b1)

    a1=einsum('ijk,ijk->ij',cross(-v1[:,None],v2),p04[:,None]-p01)
    b1=(1/einsum('ijk,ijk->ij',array([-v1]*j).transpose(1,0,2),array([cross(v2,v3)+.00001]*i)))
    t3=einsum('ij,ij->ij',a1,b1)

    condition=(t1>=0) & (t1<=1) & (t2>=0) & (t2<=1) & (t3>=0) & (t3<=1) & (t2+t3>=0) & (t2+t3<=1)
#     condition.shape

    pnt3=(p04[:,None]+einsum('ijk,ij->ijk',array([v1]*j).transpose(1,0,2),t1))[condition]
    pnt3=sort_pointsv(pnt1,pnt3) if len(pnt1)!=len(pnt3) else pnt3.tolist()

#     if o==0:
    cir=array([[pnt1[i]+array(q(v_rot[i],b[i]*r,-t)) for t in linspace(-90,90,5)] for i in arange(len(pnt1))]).tolist()
#     else:
#         cir=array([[pnt1[i]+array(q(v_rot[i],b[i]*r,t)) for t in linspace(-90,90,5)] for i in arange(len(pnt1))]).tolist()
    
    pa=[[[[p1[i][j],p1[i][j+1],p1[i+1][j]],[p1[i+1][j+1],p1[i+1][j],p1[i][j+1]]] if j<len(p1[0])-1 else \
         [[p1[i][j],p1[i][0],p1[i+1][j]],[p1[i+1][0],p1[i+1][j],p1[i][0]]] \
         for j in range(len(p1[0]))] for i in range(len(p1)-1)]
    pa=array(pa).reshape(-1,3,3)
    
    pc=array([[[cir[i][j],cir[i][j+1]]  for j in arange(len(cir[0])-1)] for i in arange(len(cir))]).reshape(-1,2,3)

    p01,p02,p03,p04,p05=pa[:,0],pa[:,1],pa[:,2],pc[:,0],pc[:,1]

    v1,v2,v3=p05-p04,p02-p01,p03-p01
    i,j=len(v1),len(v2)

    a1=einsum('ijk,ijk->ij',array([cross(v2,v3)]*i),p04[:,None]-p01)
    b1=(1/einsum('ijk,ijk->ij',array([-v1]*j).transpose(1,0,2),array([cross(v2,v3)+.00001]*i)))
    t1=einsum('ij,ij->ij',a1,b1)

    a1=einsum('ijk,ijk->ij',cross(v3,-v1[:,None]),p04[:,None]-p01)
    b1=(1/einsum('ijk,ijk->ij',array([-v1]*j).transpose(1,0,2),array([cross(v2,v3)+.00001]*i)))
    t2=einsum('ij,ij->ij',a1,b1)

    a1=einsum('ijk,ijk->ij',cross(-v1[:,None],v2),p04[:,None]-p01)
    b1=(1/einsum('ijk,ijk->ij',array([-v1]*j).transpose(1,0,2),array([cross(v2,v3)+.00001]*i)))
    t3=einsum('ij,ij->ij',a1,b1)

    condition=(t1>=0) & (t1<=1) & (t2>=0) & (t2<=1) & (t3>=0) & (t3<=1) & (t2+t3>=0) & (t2+t3<=1)
#     condition.shape

    pnt2=(p04[:,None]+einsum('ijk,ij->ijk',array([v1]*j).transpose(1,0,2),t1))[condition]
    pnt2=sort_pointsv(pnt1,pnt2) if len(pnt2)!=len(pnt1) else pnt2.tolist()

    
    sol=array([pnt3,pnt1,pnt2]).transpose(1,0,2)
    sol=[fillet_3p_3d(p3,p2,p1,r_3p_3d([p1,p2,p3])*f,s) for (p1,p2,p3) in sol]
    sol=sol+[sol[0]]
    return sol
    
def fillet_sol2sol_co(p=[],p1=[],r=1,s=10,o=0,f=1.8):
    '''
    fillet with changed orientation
    many times it is helpful
    see example in file 'examples of various functions'
    
    '''
    sol=fillet_sol2sol(p,p1,r,s,o,f)
    return cpo(sol)[1:]


def fillet_surf2sol(p=[],p1=[],r=1,s=10,o=0,f=1.8):
    '''
    function to calculate fillet at the intersection point of 2 solids
    'p': solid 1
    'p1': solid 2
    'r': radius of the fillet
    's': number of segments in the fillet, more number of segments will give finer finish
    'o': option '0' produces fillet in outer side of the intersection and '1' in the inner side of the intersections
    refer file "example of various functions" for application
    '''
    pa=[[[[p[i][j],p[i][j+1],p[i+1][j]],[p[i+1][j+1],p[i+1][j],p[i][j+1]]] if j<len(p[0])-1 else \
         [[p[i][j],p[i][0],p[i+1][j]],[p[i+1][0],p[i+1][j],p[i][0]]] \
         for j in range(len(p[0])-1)] for i in range(len(p)-1)]
    pa=array(pa).reshape(-1,3,3)

#     pb=[[[p1[i][j],p1[i+1][j]] for j in range(len(p1[0]))] for i in range(len(p1)-1)]
#     pb=array(pb).reshape(-1,2,3)
    p2=cpo(p1)
    pb=[[[p2[i][j],p2[i][j+1]] for j in range(len(p2[0])-1)] for i in range(len(p2))]
    pb=array(pb).reshape(-1,2,3)
    
    p01,p02,p03,p04,p05=pa[:,0],pa[:,1],pa[:,2],pb[:,0],pb[:,1]

    v1,v2,v3=p05-p04,p02-p01,p03-p01
    i,j=len(v1),len(v2)
#     array([(-v1).tolist()]*j).transpose(1,0,2).shape,array([v2]*i).shape,array([cross(v2,v3)]*i).shape,(p04[:,None]-p01).shape
#     cross(v3,-v1[:,None]).shape,cross(-v1[:,None],v2).shape
    a=einsum('ijk,ijk->ij',array([cross(v2,v3)]*i),p04[:,None]-p01)
    b=(1/einsum('ijk,ijk->ij',array([-v1]*j).transpose(1,0,2),array([cross(v2,v3)+.00001]*i)))
    t1=einsum('ij,ij->ij',a,b)

    a=einsum('ijk,ijk->ij',cross(v3,-v1[:,None]),p04[:,None]-p01)
    b=(1/einsum('ijk,ijk->ij',array([-v1]*j).transpose(1,0,2),array([cross(v2,v3)+.00001]*i)))
    t2=einsum('ij,ij->ij',a,b)

    a=einsum('ijk,ijk->ij',cross(-v1[:,None],v2),p04[:,None]-p01)
    b=(1/einsum('ijk,ijk->ij',array([-v1]*j).transpose(1,0,2),array([cross(v2,v3)+.00001]*i)))
    t3=einsum('ij,ij->ij',a,b)

    condition=(t1>=0) & (t1<=1) & (t2>=0) & (t2<=1) & (t3>=0) & (t3<=1) & (t2+t3>=0) & (t2+t3<=1)
#     condition.shape

    pnt1=(p04[:,None]+einsum('ijk,ij->ijk',array([v1]*j).transpose(1,0,2),t1))[condition]
#     pnt1.shape

    uv1=v1/norm(v1,axis=1).reshape(-1,1)
    uv1=array([uv1]*j).transpose(1,0,2)[condition]
#     uv1.shape
#     pnt2=pnt1+uv1*r

    a=cross(v2,v3)
    b=a/(norm(a,axis=1).reshape(-1,1)+.00001)
    b=array([b]*i)[condition]
#     b.shape

    nxt_pnt=array(pnt1[1:].tolist()+[pnt1[0]])
    v_rot=nxt_pnt-pnt1

    if o==0:
        cir=array([[pnt1[i]+array(q(v_rot[i],b[i]*r,t)) for t in linspace(0,180,5)] for i in arange(len(pnt1))]).tolist()
    else:
        cir=array([[pnt1[i]+array(q(v_rot[i],b[i]*r,-t)) for t in linspace(0,180,5)] for i in arange(len(pnt1))]).tolist()

    pc=array([[[cir[i][j],cir[i][j+1]]  for j in arange(len(cir[0])-1)] for i in arange(len(cir))]).reshape(-1,2,3)
#     pc.shape

    p01,p02,p03,p04,p05=pa[:,0],pa[:,1],pa[:,2],pc[:,0],pc[:,1]

    v1,v2,v3=p05-p04,p02-p01,p03-p01
    i,j=len(v1),len(v2)

    a1=einsum('ijk,ijk->ij',array([cross(v2,v3)]*i),p04[:,None]-p01)
    b1=(1/einsum('ijk,ijk->ij',array([-v1]*j).transpose(1,0,2),array([cross(v2,v3)+.00001]*i)))
    t1=einsum('ij,ij->ij',a1,b1)

    a1=einsum('ijk,ijk->ij',cross(v3,-v1[:,None]),p04[:,None]-p01)
    b1=(1/einsum('ijk,ijk->ij',array([-v1]*j).transpose(1,0,2),array([cross(v2,v3)+.00001]*i)))
    t2=einsum('ij,ij->ij',a1,b1)

    a1=einsum('ijk,ijk->ij',cross(-v1[:,None],v2),p04[:,None]-p01)
    b1=(1/einsum('ijk,ijk->ij',array([-v1]*j).transpose(1,0,2),array([cross(v2,v3)+.00001]*i)))
    t3=einsum('ij,ij->ij',a1,b1)

    condition=(t1>=0) & (t1<=1) & (t2>=0) & (t2<=1) & (t3>=0) & (t3<=1) & (t2+t3>=0) & (t2+t3<=1)
#     condition.shape

    pnt3=(p04[:,None]+einsum('ijk,ij->ijk',array([v1]*j).transpose(1,0,2),t1))[condition]
    pnt3=sort_pointsv(pnt1,pnt3) if len(pnt3)!= len(pnt1) else pnt3.tolist()

#     if o==0:
    cir=array([[pnt1[i]+array(q(v_rot[i],b[i]*r,-t)) for t in linspace(-90,90,5)] for i in arange(len(pnt1))]).tolist()
#     else:
#         cir=array([[pnt1[i]+array(q(v_rot[i],b[i]*r,t)) for t in linspace(-90,90,5)] for i in arange(len(pnt1))]).tolist()
    
    pa=[[[[p1[i][j],p1[i][j+1],p1[i+1][j]],[p1[i+1][j+1],p1[i+1][j],p1[i][j+1]]] if j<len(p1[0])-1 else \
         [[p1[i][j],p1[i][0],p1[i+1][j]],[p1[i+1][0],p1[i+1][j],p1[i][0]]] \
         for j in range(len(p1[0]))] for i in range(len(p1)-1)]
    pa=array(pa).reshape(-1,3,3)
    
    pc=array([[[cir[i][j],cir[i][j+1]]  for j in arange(len(cir[0])-1)] for i in arange(len(cir))]).reshape(-1,2,3)

    p01,p02,p03,p04,p05=pa[:,0],pa[:,1],pa[:,2],pc[:,0],pc[:,1]

    v1,v2,v3=p05-p04,p02-p01,p03-p01
    i,j=len(v1),len(v2)

    a1=einsum('ijk,ijk->ij',array([cross(v2,v3)]*i),p04[:,None]-p01)
    b1=(1/einsum('ijk,ijk->ij',array([-v1]*j).transpose(1,0,2),array([cross(v2,v3)+.00001]*i)))
    t1=einsum('ij,ij->ij',a1,b1)

    a1=einsum('ijk,ijk->ij',cross(v3,-v1[:,None]),p04[:,None]-p01)
    b1=(1/einsum('ijk,ijk->ij',array([-v1]*j).transpose(1,0,2),array([cross(v2,v3)+.00001]*i)))
    t2=einsum('ij,ij->ij',a1,b1)

    a1=einsum('ijk,ijk->ij',cross(-v1[:,None],v2),p04[:,None]-p01)
    b1=(1/einsum('ijk,ijk->ij',array([-v1]*j).transpose(1,0,2),array([cross(v2,v3)+.00001]*i)))
    t3=einsum('ij,ij->ij',a1,b1)

    condition=(t1>=0) & (t1<=1) & (t2>=0) & (t2<=1) & (t3>=0) & (t3<=1) & (t2+t3>=0) & (t2+t3<=1)
#     condition.shape

    pnt2=(p04[:,None]+einsum('ijk,ij->ijk',array([v1]*j).transpose(1,0,2),t1))[condition]
    pnt2=sort_pointsv(pnt1,pnt2) if len(pnt2)!= len(pnt1) else pnt2.tolist()

    
    sol=array([pnt3,pnt1,pnt2]).transpose(1,0,2)
    sol=[fillet_3p_3d(p3,p2,p1,r_3p_3d([p1,p2,p3])*f,s) for (p1,p2,p3) in sol]
    sol=sol+[sol[0]]
    return sol
    
def fillet_surf2sol_co(p=[],p1=[],r=1,s=10,o=0,f=1.8):
    '''
    fillet with changed orientation
    many times it is helpful
    see example in file 'examples of various functions'
    
    '''
    sol=fillet_surf2sol(p,p1,r,s,o,f)
    return cpo(sol)[1:]


def sec_radiuses(sec):
    a=list_r(sec)
    if (a==a[0]).all() or list_r(sec).max()>array(bb2d(sec)).max():
        return zeros(len(sec)).tolist()
    else:
        a=list_r(sec)
        b=[a[i] for i in range(len(a)-1) if a[i+1]==0]
        c=[b[i] for i in range(len(b)-1) if b[i+1]==0]+[b[-1]]
        d= [0.0]+c if len(c)!=len(convert_secv1(sec)) else c
        return d

def bb2d(sec):
    return [array(sec)[:,0].max()-array(sec)[:,0].min(),array(sec)[:,1].max()-array(sec)[:,1].min()]




    
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
    '''
    calculates the intersections of adjacent line segments only
    '''
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



def c2ro(sol,s):#circular to rectangulat orientation
    '''
    change the orientation of points of a cylinder from circular to rectangular orientation
    'sol': is a cylindrical type 3d shape
    's': number of segments required between each straight line segments
    refer to the file 'example of various functions' for application examples 
    '''
    angle=360/len(sol[0])/2
    sol=cpo(sol)
    return q_rot([f'z{angle}'],[m_points1(sol[i]+flip(sol[len(sol)-1-i]),s) for i in range(int(len(sol)/2))])


def c2ro_(sol,s):#circular to rectangulat orientation
    '''
    change the orientation of points of a cylinder from circular to rectangular orientation
    'sol': is a cylindrical type 3d shape
    's': number of segments required between each straight line segments
    refer to the file 'example of various functions' for application examples 
    '''
    angle=360/len(sol[0])/2
    sol=cpo(sol)
    return q_rot([f'z{0}'],[m_points1(sol[i]+flip(sol[len(sol)-1-i]),s) for i in range(int(len(sol)/2))])
    
    
def vsp_extrude(sec,extrude_path, shape_path):
    '''
    function variable section and path extrude
    sec: section to extrude
    extrude_path: is the path on which the section needs to be extruded
    shape_path: sculpting path
    
    extrude path should always be a little longer than the sculpting shape
    an example will make this more clear
    refer to the file "example of various functions" for the same
    '''
    path=shape_path
    path1=extrude_path
    path2=path1[:-1]
    path3=path1[1:]
    path,path2,path3=array(path),array(path2),array(path3)
    v1=array([path3-path2]*len(path)).transpose(1,0,2)
    v2=path-path2[:,None]
    v1.shape
    v1norm=sqrt(einsum('ijk,ijk->ij',v1,v1))
    inv_v1norm=1/v1norm
    v1.shape,v1norm.shape
    u1=einsum('ijk,ij->ijk',v1,inv_v1norm)
    u1.shape,v2.shape
    p1=einsum('ijk,ijk->ij',u1,v2)
    a,b=p1.shape
    decision=(zeros(a*b).reshape(a,b)<=p1)&(p1<=v1norm)
    p0=array([path2]*len(path)).transpose(1,0,2)
    p0.shape,u1.shape,p1.shape
    points=p0+einsum('ijk,ij->ijk',u1,p1)
    points=array(sort_points(path,points[decision]))
    os=[norm(path[i]-points[i])-norm(path[0]-points[0]) for i in range(len(path))]

    sections=[offset(sec,p) for p in os]

    p=array(cytz(points))
    s2=[]
    for i in range(len(p)-1):
        s=q_rot(['x90','z-90'],sections[i])
        v1=p[i+1]-p[i]+array([0,0,0.00001])
        va=[v1[0],v1[1],0]
        u1=array(uv(v1))
        ua=array(uv(va))
        v2=cross(va,v1)
        a1=arccos(u1@ua)*180/pi
        a2=ang(v1[0],v1[1])
        s1=q_rot([f'z{a2}'],s)
        if i<len(p)-1:
            s2.append(translate(p[i],[q(v2,p,a1) for p in s1]))
        else:
            s2.append(translate(p[i],[q(v2,p,a1) for p in s1]))
            s2.append(translate(p[i+1],[q(v2,p,a1) for p in s1]))

    s3=flip([[p for p in p1 if ~isnan(p[0])] for p1 in s2])
    s3=[p for p in s3 if p!=[]]
    return s3


def sl_int(sec,line):
    '''
    section-line intersection:
    function to find intersection between an enclosed section in 3d space and a line
    refer to the file "example of various functions" for application example
    '''
    vectors=array([p-array(sec[0]) for p in array(sec[1:])])
    segments=array(seg(vectors))[:-1]
    v1=array([array(line[1])-array(line[0])])
    im=inv(array([concatenate([v1,-p]) for p in segments]))
    point=array([array(sec[0])-array(line[0])]*len(im))
    im.shape,point.shape
    t=einsum('ijk,ij->ik',im,point)
    p0=array(line[0])

    decision=(t[:,0]>=0)&(t[:,0]<=1)&(t[:,1]>=0)&(t[:,1]<=1)&(t[:,2]>=0)&(t[:,2]<=1)&((t[:,1]+t[:,2])<=1)
    p0.shape,v1.shape,t.shape
    v2=array([v1]*len(t)).reshape(-1,3)
    intersection_point=(p0+einsum('ij,i->ij',v2,t[:,0])[decision]).tolist()

    return intersection_point
    
def pts2(path):
    '''
    returns the cumulative sum of points
    example:
    path=[[0,0,1],[0,5,10],[10,3,20]]
    pts2(path)=> [[0, 0, 1], [0, 5, 11], [10, 8, 31]]
    
    '''
    return array(path).cumsum(0).tolist()
    
def axis_rot(axis,solid,angle):
    '''
    rotate a solid around an axis
    '''
    if len(array(solid).shape)==3:
        return[[q(axis,p1,angle) for p1 in p] for p in solid]
    else:
        return [q(axis,p,angle) for p in solid]
    
def end_cap(fillet,f=-1):
    '''
    for giving fillet at the end of a path_extruded solid
    factor 'f' can be -1 or 1 depending on which side fillet needs to be created
    
    '''
    fillet1=cpo(fillet)[1:]
    v1=nv(fillet1[0])
    fillet01=translate(f*array(v1),scl3dc(fillet1,1.5))
    fillet1=swp_prism_h(fillet01,fillet1)
    return fillet1

        
def d2r(d):
    '''
    converts degrees to radians
    
    '''
    return radians(d)
def r2d(r):
    '''
    converts radians to degrees
    '''
    return rad2deg(r)
    
def convert_3lines2fillet(pnt1,pnt2,pnt3,f=1.9,s=10):
    '''
    Develops a fillet with 3 list of points (pnt1,pnt2,pnt3) in 3d space
    f: is a factor which can be reduced to 1.5 in case of self intersection observed
    s: number of segments in the fillet, increase the segments in case finer model is required
    refer to the file "example of various functions" for application examples
    
    '''
    
    
    sol=array([pnt3,pnt1,pnt2]).transpose(1,0,2)
    sol=[fillet_3p_3d(p3,p2,p1,r_3p_3d([p1,p2,p3])*f,s) for (p1,p2,p3) in sol]
    sol=sol
    return sol
    
def min_d_points(sec,min_d=.1):
    ''' 
    rationalises the number points in a section based on the minimum distance between 2 points
    i.e. all the points which are less than the defined minimum distance "min_d" will be omitted from the section "sec" 
    
    '''
    b=sec[0]
    c=[sec[0]]
    for i in range(1,len(sec)):
        if l_len([b,sec[i]])>=min_d:
            c.append(sec[i])
            b=sec[i]
            
    return c
    
def wrap_around(sec,path):
    sec=sec if array(sec).shape[-1]==3 else c2t3(sec)
    path=array(path)
    dy=[0]+[l_len([path[i],path[i+1]]) for i in range(len(path)-1)]
    v1=[path[i+1]-path[i] for i in range(len(path)-1)]
    dy,v1=array(dy),array(v1)
    y=dy.cumsum()
    l1=array([p for p in array(sec)[:,1]])

    sec1=[]
    for n in range(len(sec)):
        m=arange(len(y))[y>l1[n]][0]
        m=0 if m-1<0 else m-1
        l2=l1[n]-y[m]
        l2=l1[n] if l2<0 else l2
        p2=path[m]+v1[m]*l2/dy[m+1]+array([sec[n][0],0,sec[n][2]])
        sec1.append(p2.tolist())
    return sec1
    
def wrap_around1(sec,path):
    sec=sec if array(sec).shape[-1]==3 else c2t3(sec)
    path=array(path)
    dy=[0]+[l_len([path[i],path[i+1]]) for i in range(len(path)-1)]
    v1=[path[i+1]-path[i] for i in range(len(path)-1)]
    dy,v1=array(dy),array(v1)
    y=dy.cumsum()
    l1=array([p for p in array(sec)[:,1]])

    sec1=[]
    for n in range(len(sec)):
        m=arange(len(y))[y>l1[n]][0]
        m=0 if m-1<0 else m-1
        l2=l1[n]-y[m]
        l2=l1[n] if l2<0 else l2
        p2=path[m]+v1[m]*l2/dy[m+1]+array([sec[n][0],0,sec[n][2]])
        sec1.append(p2.tolist())
    return sec1


#def align_sec(sec1,sec2,ang=10):
#    '''
#    function to align 2 3d sections to obtain the non twisted optimised solid
#    ang: is the resolution for the angle of rotation, 1 degree will have much higher resolution and hence will take longer to compute
#    refer file "examples of various functions" for application examples
#    '''
#    nv1=nv(sec2)
#    cp1=array(sec2).mean(0)
#    sec2=translate(-cp1,sec2)
#    i=arange(0,360,ang)
#    area1=[]
#    sol1=[]
#    for j in i:
#        sec3=(array(axis_rot(nv1,sec2,j))+cp1).tolist()
#        sol=[sec1]+[sec3]
#        a1=array([l_len(p) for p in cpo(sol)]).sum()
#        area1.append(a1)
#        sol1.append(sol)
#    sol2=sol1[array(area1).argmin()]
#    return sol2
    
def align_sec(sec1,sec2,ang=10):
    '''
    function to align 2 3d sections to obtain the non twisted optimised solid
    ang: is the resolution for the angle of rotation, 1 degree will have much higher resolution and hence will take longer to compute
    refer file "examples of various functions" for application examples
    '''
    nv1=nv(sec2)
    cp1=array(sec2).mean(0)
    sec2=translate(-cp1,sec2)
    i=arange(0,360,ang)
    area1=[norm((array(axis_rot(nv1,sec2,j))+cp1)-array(sec1),axis=1).sum() for j in i]
    sec3=(array(axis_rot(nv1,sec2,array(area1).argmin()*ang))+cp1).tolist()
    sol2=[sec1]+[sec3]
    return sol2
    
def sec2vector(v1=[1,0,0],sec=[]):
    '''
    function to align a section 'sec' with a vector 'v1'
    refer file "example of various function" for application examples
    '''
    vz=[0,0,-1]
    vz,v1=array([vz,v1])

    nvzv1=cross(vz,v1)
    u1=v1/norm(v1)
    theta=r2d(arccos(u1@vz))
    sec=flip(sec)

#    sec=flip(q_rot(['x180'],sec))
    sec1=axis_rot([1,0,0],sec,theta)

    theta1=ang(v1[0],v1[1])
    sec1=q_rot(['z-90',f'z{theta1}'],sec1)
    return sec1
    
def sec2vector1(v1,sec):
    '''
    same as sec2vector but simpler method
    '''
    theta_y=ang((v1[0]**2+v1[1]**2)**.5,v1[2])
    theta_z=ang(v1[0],v1[1])
    return q_rot(['x90','z-90',f'y{-theta_y}',f'z{theta_z}'],sec)

    
#def sec2vector(v1=[1,0,0],sec=[]):
#    '''
#    function to align a section 'sec' with a vector 'v1'
#    refer file "example of various function" for application examples
#    '''
#    sec=flip(sec)
#    th1=ang(v1[0],v1[1])
#    s1=q_rot(['x90','z-90',f'z{th1}'],sec)
#    a=array([v1[0],v1[1],0])
#    b=array(v1)
#    th2=r2d(arccos((a/norm(a))@(b/norm(b)))) if v1[2]!=0 else 0
#    c=cross(a,b)
#    s1=axis_rot(c,s1,th2)
#    return s1



def cut_plane(nv=[0,0,1],size=[5,5],thickness=10,trns1=0,trns2=0,trns3=0): #oriented solid
    '''
    function for defining a solid (cutting plane) oriented as per the defined normal vector
    nv: normal vector for defining plane orientation of the section
    thickness: thickness or height of the cutting plane
    trns1: translate the solid in the direction of normal vector 'nv'
    trns2: translate the solid in the direction 'right' to the normal vector 'nv'
    trns3: translate the solid in the direction 'up' to the normal vector 'nv'
    '-ve' values given to the trns1,trns2,trns3 will translate the solid in the reverse direction 
    '''
    sec=square(size,center=True)
    plane1=sec2vector(nv,sec)
    v1=array(nv)
    u1=v1/norm(v1)
    ua=array([0,0,1]) if u1[2]==0 else array([0,-1,0]) if (u1==[0,0,1]).all() else array([0,1,0]) if (u1==[0,0,-1]).all() else array([0,0,1])
    v2=cross(u1,ua) if u1[2]>0 else cross(u1,ua)
    u2=v2/norm(v2)
    v3=cross(u2,u1) if u1[2]==0 else cross(u2,u1)
    u3=v3/norm(v3)
#     u1,u2,u3=array([u1,u2,u3]).tolist()
    plane2=translate(u1*thickness,plane1)
    sol=[plane1]+[plane2]
    sol=translate(u1*trns1,sol)
    sol=translate(u2*trns2,sol)
    sol=translate(u3*trns3,sol)
    return sol
    
    
def slice_sol(sol,n=10):
    '''
    function to slice a solid with 'n' intermediate steps
    '''
    a=cpo(sol)
    sol1=[[ls(p,n)+[p[1]] for p in seg(a[i])[:-1]] for i in range(len(a))]
    sol2=array(sol1).transpose(1,2,0,3)
    b,c,d,e=sol2.shape
    sol2=sol2.reshape(b*c,d,e).tolist()
    return sol2
    
def cp_arc(arc1):
    '''
    function returns the center point of a given circle or arc
    
    '''
    n=int(len(arc1)/360*120)
    p0=arc1[0]
    p1=arc1[n]
    p2=arc1[n*2]
    return cp_3p(p0,p1,p2)
    
def r_arc(arc1):
    '''
    function returns the radius of a given circle or arc
    
    '''
    n=int(len(arc1)/360*120)
    p0=arc1[0]
    p1=arc1[n]
    p2=arc1[n*2]
    return r_3p([p0,p1,p2])
    
def fillet_l_cir(line=[],cir1=[],fillet_radius=1,s=20):
    '''
    function to draw fillet between a line and a circle
    '''
    p0,p1=array(line)
    cp=array(cp_arc(cir1))
    r1=r_arc(cir1)
    r2=fillet_radius
    v1=p1-p0
    v2=cp-p0
    u1=v1/norm(v1)
    u2=v2/norm(v2)
    d1=u1@v2
    p2=p0+u1*d1
    v3=p2-cp
    u3=v3/norm(v3)
    h=norm(p2-cp)-r2
    r=r1+r2
    d=sqrt(r**2-h**2)
    cp1=cp+u3*h-u1*d
    cp2=cp+u3*h+u1*d
    p3=p0+u1*(d1-d)
    p4=l_cir_ip([cp1,cp],cir1)[0]
    p5=p0+u1*(d1+d)
    p6=l_cir_ip([cp2,cp],cir1)[0]
    fillet1=arc_2p(p3,p4,r2,cw([p0,p2,cp]),s=s)
    fillet2=arc_2p(p5,p6,r2,cw([p1,p2,cp]),s=s)
    return [fillet1, fillet2]




def o_solid(nv=[0,0,1],sec=[],thickness=10,trns1=0,trns2=0,trns3=0, theta=[0,0,0]): #oriented solid
    '''
    function for defining a solid with any defined section. solid gets oriented as per the defined normal vector
    nv: normal vector for defining plane orientation of the section
    sec: cross section of the solid
    thickness: thickness or height of the solid
    trns1: translate the solid in the direction of normal vector 'nv'
    trns2: translate the solid in the direction 'right' to the normal vector 'nv'
    trns3: translate the solid in the direction 'up' to the normal vector 'nv'
    '-ve' values given to the trns1,trns2,trns3 will translate the solid in the reverse direction 
    theta: rotate the section around axis  fox example if nv is [1,0,0] or x-axis, the sequence of rotation will be x, y ,z axis
    '''
    plane1=sec2vector(nv,sec)
    v1=array(nv)
    u1=v1/norm(v1)
    ua=array([0,0,1]) if u1[2]==0 else array([0,-1,0]) if (u1==[0,0,1]).all() else array([0,1,0]) if (u1==[0,0,-1]).all() else array([0,0,1])
    v2=cross(u1,ua) if u1[2]>0 else cross(u1,ua)
    u2=v2/norm(v2)
    v3=cross(u2,u1) if u1[2]==0 else cross(u2,u1)
    u3=v3/norm(v3)
#     u1,u2,u3=array([u1,u2,u3]).tolist()
    plane2=translate(u1*thickness,plane1)
    sol=[plane1]+[plane2]
    sol=axis_rot(u3,axis_rot(-u2,axis_rot(nv,sol,theta[0]),theta[1]),theta[2])
    sol=translate(u1*trns1,sol)
    sol=translate(u2*trns2,sol)
    sol=translate(u3*trns3,sol)
    return sol
    
def ppesec(p0,sec): #point's projection on an enclosed 3d section
    '''
    function to find projected points of a given point list 'p0' on a 3d sec which is on 1 plane
'''
    
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
    '''
    function to find projected points of a given list of points 'p0' on a plane defined by normal'v1' and location 'loc'
    example:
    p0=[20,0,0]
    v1=[2,3,4]
    loc=[0,10,0]
    ppplane([p0],v1,loc) => [19.310359216374945, -1.034461175437585, -1.3792815672501133]
    
    '''
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
    '''
    function to draw a honeycomb structure with radius 'r' 
    n1: number of hexagons in 1 layer
    n2: number of layers
    '''
    cir1=circle(r,s=7)
    cir2=c3t2(q_rot(['z30'],cir1))
    sec=[translate([i,0,0],cir1) for i in arange(0,3*n1*r,3*r)]
    sec1=[translate([i,r*sin(d2r(60)),0],cir1) for i in arange(r*1.5,3*n1*r,3*r)]
    sec2=array([sec,sec1]).transpose(1,0,2,3)
    a,b,c,d=sec2.shape
    sec3=sec2.reshape(a*b,c,d)
    sec3=array([translate([0,i,0],sec3) for i in arange(0,r*sin(d2r(60))*(n2*2),2*r*sin(d2r(60)))])
    a,b,c,d=sec3.shape
    sec4=c3t2(sec3.reshape(a*b,c,d))
    return sec4

def offset_3d(sec,d):
    '''
    offsets an enclosed section in 3d space, in case the section is in 1 plane
    sec: section in 3d space
    d: offset distance -ve sign means inner offset and +ve sign is outer offset
    refer to the file"example of various functions" for application examples
    
    '''
    sec=remove_extra_points(sec)
    avg1=array(sec).mean(0)
    sec1=translate(-avg1,sec)
#     v1=array([array(p)-avg1 for p in sec]).tolist()
#     v2=v1[1:]+[v1[0]]
#     v1,v2=array([v1,v2])
#     n1=cross(v1,v2)
#     nv1=n1.mean(0)
    nv1=-array(nv(sec1))
    nz=[0,0,1]
    nr=cross(nv1,nz) if abs(nv1).tolist()!=[0,0,1] else nv1
    theta=r2d(arccos(nv1@array(nz)))
    sec1=axis_rot(nr,sec1,theta)
    z_values=array(sec1)[:,2]-avg1[2]
    sec1=ppplane(sec1,[0,0,1],[0,0,0])
    sec1=c3t2(sec1)
    x_values=array([l_len([[0,0],p])  for p in sec1])
    sec2=offset(c3t2(sec1),d)
    x1_values=array([l_len([[0,0],p])  for p in sec2])
    z1_values=z_values/x_values*x1_values
    z1_values=array([[0,0,p] for p in z1_values])
    sec2=array(c2t3(sec2))
    sec2=axis_rot(nr,sec2,-theta)
    sec2=translate(array(sec).mean(0),sec2)
    return sec2
    
    
    
def path_extrude2msec(sec_list,path):
    '''
    function to extrude multiple sections 'sec_list' along an open path 'path'
    number of sections in the 'sec_list' >= len(path)
    refer to file "example of various functions" for application example
    '''
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
    '''
    function to calculate intersection point between a 3d solid and a line. 
     "sol" is the 3d object which is intersected with a "line".
     try below code for better understanding:
    sec=circle(10)
    path=cr(pts1([[-10+.1,0],[12,0],[-2,0,2],[0,10,3],[-10,0]]),5)
    sol=prism(sec,path)

    line=[[0,0,-1],[20,20,10]]

    ip1=ip_sol2line(sol,line)
    
    refer to file "example of various functions" for application
    '''


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
#    s_planes=array([p1]*a).transpose(1,0,2,3)[condition][argsort([norm(p-px[0]) for p in i_p1])]
#    nv1=[nv(p) for p in s_planes]

#    i_p2,s_planes,nv1=i_p2.tolist(),s_planes.tolist(),array(nv1).tolist()
    return i_p2.tolist()
    


def align_sol(sol,ang=10):
    '''
    function to straighten the twists in the path_extruded sections for better alignments
    refer to the file "example of various functions.ipynb" for application examples
    '''
    sol1=[sol[0]]
    for i in range(1,len(sol)):
        a=align_sec(sol1[i-1],sol[i],ang=ang)
        sol1.append(a[1])
    return sol1

def extrude_sol2path(sec,path1,path2):
    '''
    extrude a solid to a different path
    "sec" and "path1" defines the original solid
    "path2" defines the path where the shape of the original solid to be extruded
    refer file "example of various functions.ipynb" for application example
    '''
    min_l=array([l_len(p) for p in seg(path1)[:-1]]).min()
    path1=m_points_o(path1,min_l)
    l1=array([l_len(p) for p in seg(path2)[:-1]]).sum()
    l2=l1/len(path1)
    path3=m_points_o(path2,l2)[:len(path1)]
    #    path3=bezier(path2,len(path1))
    sec_list=[offset(sec,x) for (x,y) in path1]
    sol=path_extrude2msec(sec_list,path3)
    return sol
    


def ip_normal_sol2line(sol,line):
    '''
    function to find the normal from intersection points between a 3d solid and a line. 
     "sol" is the 3d object which is intersected with a "line".
     try below code for better understanding:
    sec=circle(10)
    path=cr(pts1([[-10+.1,0],[12,0],[-2,0,2],[0,10,3],[-10,0]]),5)
    sol=prism(sec,path)

    line=[[0,0,-1],[20,20,10]]

    ip1=ip_normal_sol2line(sol,line)
    
    refer to file "example of various functions" for application
    '''


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
    '''
    returns vertices and faces for a closed loop solid
    
    '''
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
    '''
    function calculates vertices and faces for a given surface
    refer file 'example of various functions.ipynb' for application examples
    '''
    n1,n2,_=array(surf).shape
    v=array(surf).reshape(-1,3)    
    f1=array([[[[i*n2+j,i*n2+j+1,(i+1)*n2+j],[(i+1)*n2+j,i*n2+j+1,(i+1)*n2+j+1]] for j in range(n2-1)] for i in range(n1-1)]).reshape(-1,3)
    return [v.tolist(),f1.tolist()]
    
    
def resurface(v,f1):
    '''
    arranges the vertices in a circular fashion for a surface defined by v- vertices and f1- triangulated faces
    refer file 'example of various functions.ipynb' for application examples
    '''
    v,f1=array(v),array(f1)
    f2=f1
    v3=[]
    while(len(f2)>=1):
        a=igl.boundary_loop(f2)
        v2=v[a]
        if len(v3)>0:
            v3.append(sort_points(v3[-1],v2.tolist()))
        else:
            v3.append(v2.tolist())
        f2=f2[~isin(f2,a).any(1)]
    return v3
    
def equidistant_path(path,s=10):
    '''
    divides a path in to equally spaced points
    refer file 'example of various functions.ipynb' for application examples
    '''
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
    return p_rev

def equidistant_pathc(path,s=10):
    '''
    divides a closed path in to equally spaced points
    refer file 'example of various functions.ipynb' for application examples
    '''
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
    return p_rev[:s]





    
def surface_base(v,f2,h,up=0):
    '''
    makes a solid with flat base keeping the top surface as the defined surface
    v : vertices
    f2: triangulated faces for defining the surface for the purpose of drawing polyhedron
    h: 'z' height of the base e.g. 0 means the base will be drawn at z=0 ans so on.
    up: if the base is below the surface 'up' should be 0 and if the base is above the surface 'up' should be 1 
    
    
    '''
    if up==0:
        f2=[flip(p) for p in f2]
        v1=array(c2t3(c3t2(v)))+[0,0,h]
        v1=v1.tolist()
        v2=v+v1
        f3=array(f2)+array(f2).max()+1
        f3=f3.tolist()
        f3=[flip(p) for p in f3]
        f3p=igl.boundary_loop(array(f3)).tolist()
        f2p=flip(igl.boundary_loop(array(f2))).tolist()
        f4p=f3p+f2p
        n1=len(f3p)
        f5=[ [[f4p[i],f4p[i+n1],f4p[i+1]],[f4p[i+1],f4p[i+n1],f4p[i+1+n1]]] 
            if i < n1-1 else  
            [[f4p[i],f4p[i+n1],f4p[0]],[f4p[0],f4p[i+n1],f4p[n1]]]
            for i in range(n1)]
        f5=array(f5).reshape(-1,3).tolist()
        f4=f2+f3+f5
    elif up==1:
        v1=array(c2t3(c3t2(v)))+[0,0,h]
        v1=v1.tolist()
        v2=v+v1
        f3=array(f2)+array(f2).max()+1
        f3=f3.tolist()
        f3=[flip(p) for p in f3]
        f3p=igl.boundary_loop(array(f3)).tolist()
        f2p=flip(igl.boundary_loop(array(f2))).tolist()
        f4p=f3p+f2p
        n1=len(f3p)
        f5=[ [[f4p[i],f4p[i+n1],f4p[i+1]],[f4p[i+1],f4p[i+n1],f4p[i+1+n1]]] 
            if i < n1-1 else  
            [[f4p[i],f4p[i+n1],f4p[0]],[f4p[0],f4p[i+n1],f4p[n1]]]
            for i in range(n1)]
        f5=array(f5).reshape(-1,3).tolist()
        f4=f2+f3+f5
    return [v2,f4]
    
def surface_offset_sol(v,f2,t=-1):
    '''
    function to draw a solid by offsetting the surface
    v and f2 are the vertices and faces of the surface and t is the thickness of the surface
    't' can be given any number - or +
    refer file 'example of various functions.ipynb' for application examples
    
    '''
    v1=(array(v)+igl.per_vertex_normals(array(v),array(f2))*t).tolist()
    if t<0:
        f2=[flip(p) for p in f2]
        v2=v+v1
        f3=array(f2)+array(f2).max()+1
        f3=f3.tolist()
        f3=[flip(p) for p in f3]
        f3p=igl.boundary_loop(array(f3)).tolist()
        f2p=flip(igl.boundary_loop(array(f2))).tolist()
        f4p=f3p+f2p
        n1=len(f3p)
        f5=[ [[f4p[i],f4p[i+n1],f4p[i+1]],[f4p[i+1],f4p[i+n1],f4p[i+1+n1]]] 
            if i < n1-1 else  
            [[f4p[i],f4p[i+n1],f4p[0]],[f4p[0],f4p[i+n1],f4p[n1]]]
            for i in range(n1)]
        f5=array(f5).reshape(-1,3).tolist()
        f4=f2+f3+f5
    elif t>0:
        v2=v+v1
        f3=array(f2)+array(f2).max()+1
        f3=f3.tolist()
        f3=[flip(p) for p in f3]
        f3p=igl.boundary_loop(array(f3)).tolist()
        f2p=flip(igl.boundary_loop(array(f2))).tolist()
        f4p=f3p+f2p
        n1=len(f3p)
        f5=[ [[f4p[i],f4p[i+n1],f4p[i+1]],[f4p[i+1],f4p[i+n1],f4p[i+1+n1]]] 
            if i < n1-1 else  
            [[f4p[i],f4p[i+n1],f4p[0]],[f4p[0],f4p[i+n1],f4p[n1]]]
            for i in range(n1)]
        f5=array(f5).reshape(-1,3).tolist()
        f4=f2+f3+f5
    return [v2,f4]
    
def ang_2lineccw(p0,p1,p2):
    '''
    ccw angle of the line p0p2 from base line p0p1
    '''
    p0,p1,p2=array([p0,p1,p2])
    v1,v2=p1-p0,p2-p0
    a1=ang(v1[0],v1[1])
    a2=ang(v2[0],v2[1])
    return 360 if a1-a2==0 else 360-(a1-a2) if a2<a1 else a2-a1 

def ang_2linecw(p0,p1,p2):
    '''
    cw angle of the line p0p2 from the base line p0p1
    '''
    p0,p1,p2=array([p0,p1,p2])
    v1,v2=p1-p0,p2-p0
    a1=ang(v1[0],v1[1])
    a2=ang(v2[0],v2[1])
    return 0 if a1-a2==0 else a1-a2 if a2<a1 else 360+(a1-a2) 

def l_lenv(l):
    '''
    calculates sum of lengths of all the segments in a line 'l' considering the section is closed
    '''
    return array([l_len(p) for p in seg(l)]).sum()

def l_lenv_o(l):
    '''
    calculates sum of lengths of all the segments in a line 'l' considering the section is open
    '''
    return array([l_len(p) for p in seg(l)[:-1]]).sum()


def a_3seg(s):
    '''
    area of the triangle enclosed with in 3 vertices 's'
    '''
    return norm(cross(array(s[1])-array(s[0]),array(s[2])-array(s[0])))/2




def points2line_min_d_point(line,points,f=1):
    if len(points)>0:
        line=array(line)
        pnts=array(points)
        v1=line[1]-line[0]
        u1=v1/norm(v1)
        v2=pnts-line[0]
        u2=v2/norm(v2,axis=1).reshape(-1,1)
        u1.shape,v2.shape
        v2cost=einsum('j,ij->i',u1,v2)
        v2sint=cross(v1,v2)/norm(v1)
        decision=(v2cost<=l_len(line))&(v2cost>=0)&(v2sint<l_len(line)/f)
        pnts1=pnts[decision]
        v2=pnts1-line[0]
        u2=v2/norm(v2,axis=1).reshape(-1,1)
        u1.shape,v2.shape
        v2sint=cross(v1,v2)/norm(v1)
        if v2sint.tolist()!=[]:
            d=v2sint.argmin()
            return pnts1[d].tolist()
        else:
            return []
    else:
        return []


def offset_sol(sol,d,o=0,s=100):
    '''
    function to calculate offset of a 3d object by distance 'd'
    option 'o' can be set to '0' or '1' depending on the shape of the object.
    in case the shape of the 3d object is twisted,option should be set to '1'
    in few cases this function may not work 
    
    '''
    sol=q_rot(['x.001'],sol)
    n=array([len(remove_extra_points(p)) for p in sol]).argmax()
    if o==0:
        sol=[sort_points(sol[n],offset_3d(p,d)) for p in sol]
    else:
        sol=[equidistant_pathc(offset_3d(p,d),s) for p in sol]
    return sol
    
def ip_sol2sol(sol,sol1,i=0):
    '''
    function to find the intersection point between 2 solids
    this function is to be used where the cutting lines of sol1 are intersecting sol at more than 1 times.
    sol: solid which is intersected
    sol1: this intersects the solid 'sol'
    i: if the first intersection points of all the cutting lines are to be considered, value of i should be '0'.
    if the last intersection points of all the cutting lines are to be considered, value of 'i' should be set to '-1'
    if all the intersection points are required, value of 'i' should be set to 'all'
    '''
    if i=='all':
        a=[ip_sol2line(sol,p) for p in cpo(sol1)]
        a=array([p for p in a if p!=[]]).reshape(-1,3).tolist()
    else:
        a=[ip_sol2line(sol,p) for p in cpo(sol1)]
        a=[p[i] for p in a if p!=[]]
    
    return a
    
    
def vnf2(bead2):
    '''
    function returns vertices and faces of 3d shapes with first and the last section triangulated. only works with convex sections
    
    '''
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
    n=n1+n2+n3
    pnt=[cp1]+array(bead2).reshape(-1,3).round(4).tolist()+[cp2]
    return [pnt,n]


        
def p_outside(triangle,p):
    p0,p1,p2,p3=triangle+[p]
    t1=[]
    if pies1(cir_3p(p0,p1,p2,100),[p3])==[]:
        t1.append([p0,p1,p2])
    if pies1(cir_3p(p1,p2,p3,100),[p0])==[]:
        t1.append([p1,p2,p3])
    if pies1(cir_3p(p2,p3,p0,100),[p1])==[]:
        t1.append([p2,p3,p0])
    if pies1(cir_3p(p3,p0,p1,100),[p2])==[]:
        t1.append([p3,p0,p1])
    t1=[flip(p) if cw(p)==1 else p for p in t1]
    return t1[:-1] if len(t1)>3 else t1
    

def p_inside(triangle,p):
    t1=[[p,p1[0],p1[1]] for p1 in seg(triangle) ]
    t1=[flip(p) if cw(p)==1 else p for p in t1]
    return t1

def p_online(triangle,p):
    for p1 in seg(triangle):
        v1=array(p1[1])-array(p1[0])
        u1=(v1/norm(v1)).round(3)
        v2=array(p)-array(p1[0])
        u2=(v2/norm(v2)).round(3)
        if u1.tolist()==u2.tolist() and l_len([p1[0],p1[1]])>l_len([p1[0],p]):
            px=exclude_points(triangle,p1)[0]
            t1=[[p1[0],p,px],[p1[1],p,px]]
        if u1.tolist()==u2.tolist() and l_len([p1[0],p1[1]])<l_len([p1[0],p]):
            px=exclude_points(triangle,p1)[0]
            t1=[[p1[0],p1[1],px],[p1[1],p,px]]
        elif u1.tolist()==(-u2).tolist():
            px=exclude_points(triangle,p1)[0]
            t1=[[p1[0],p1[1],px],[p1[0],px,p]]
    t1=[flip(p) if cw(p)==1 else p for p in t1]
    return t1

def triangulate_4p(triangle,p):
    c1=[]
    for p1 in seg(triangle):
        v1=array(p1[1])-array(p1[0])
        u1=(v1/norm(v1)).round(3)
        v2=array(p)-array(p1[0])
        u2=(v2/norm(v2)).round(3)
        if u1.tolist()==u2.tolist() or u1.tolist()==(-u2).tolist():
            c1.append(1)
    if c1!=[]:
        return p_online(triangle,p)
    elif pies1(triangle,[p])==[]:
        return p_outside(triangle,p)
    else:
        return p_inside(triangle,p)
    
    
def convex_hull(pnts):
    '''
    calculates convex hull for a list of points 'sec'
    
    '''
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
    '''
    function sorts the points list 'p' first with x and then with y smallest to largest 
    '''
    p=array(p)
    p1=p[p[:,0].argsort()]
    pux=unique(p1[:,0])
    p2=concatenate([p1[p1[:,0]==p][p1[p1[:,0]==p][:,1].argsort()] for p in pux])
    p2=p2.tolist()
    return p2

def lexicographic_sort_yx(p):
    '''
    function sorts the points list 'p' first with y and then with x smallest to largest 
    '''
    p=array(p)
    p1=p[p[:,1].argsort()]
    puy=unique(p1[:,1])
    p2=concatenate([p1[p1[:,1]==p][p1[p1[:,1]==p][:,0].argsort()] for p in puy])
    p2=p2.tolist()
    return p2
    
def equivalent_rot_axis(r1=[]):
    '''
    function returns an equivalent axis for rotation and angle of rotation for a sequence of rotations given by list 'r1'
    example:
    r1=['x30','y40','z100','y10','x70','y45']
    rotation_axis,theta=equivalent_rot_axis(r1)
    following is returned by the function:
    rotation_axis=> [0.33215577139188024, 0.9203587619795575, -0.0]
    theta=> 78.09373872704292
    
    '''
    vz=array([0,0,1])
    v1=array(q_rot(r1,[vz])[0])
    v2=cross(vz,v1)
    theta=r2d(arccos(vz@v1/norm(v1)))
    return [v2.tolist(),theta]
    
def path_extrude_open(sec,path,twist=0):
    '''
    function to extrude a closed section to an open path
    twist can be set either to '0' or '1' depending on the shape produced
    '''
    if twist==0:
        p1=path[:-1]
        p2=path[1:]
        p1,p2=array([p1,p2])
        v1=p2-p1
        u1=v1/norm(v1,axis=1).reshape(-1,1)
        v2=concatenate([[u1[0]],(u1[1:]+u1[:-1])/2,[u1[-1]]])
        sec2=[]
        for i in range(len(path)):
            sec1=translate(path[i],sec2vector(v2[i],sec))
            sec2.append(sec1)
        return sec2

    if twist==1:
        sec=flip(sec) if cw(sec)==-1 else sec
        p1=array(seg(path))[:-1]
        p2=array(path)
        v1=array([(p[1]-p[0])/norm(p[1]-p[0]) for p in p1])
        t_v=array([ v1[i] if i==0 else
                   (v1[i-1]+v1[i])/2 if i<len(p2)-1 else
                   v1[-1]
            for i in range(len(p2))])

        n_v=array([ cross(p2[i+1]-p2[i],p2[i+2]-p2[i+1]) if i==0 else
             cross(p2[i]-p2[i-1],p2[i+1]-p2[i]) if i<len(p2)-1 else
             cross(p2[i-1]-p2[i-2],p2[i]-p2[i-1])
            for i in range(len(p2))])
        o_v=array([cross(n_v[i],t_v[i]) for i in range(len(t_v))])

        t_v=t_v/norm(t_v,axis=1).reshape(-1,1)
        n_v=n_v/norm(n_v,axis=1).reshape(-1,1)
        o_v=o_v/norm(o_v,axis=1).reshape(-1,1)

        map_v=array([t_v,n_v,o_v]).transpose(1,0,2)
        sec2=[]
        for p in map_v:
            v2=[[0,0,0],[0,0,-1],[0,1,0]]
            a1=cross(v2[1],p[0])
            t1=r2d(arccos(array(v2[1])@p[0]))
            sec1=axis_rot(a1,sec,t1)
            v3=axis_rot(a1,v2,t1)
            a2=cross(v3[2],p[1])
            t2=r2d(arccos(array(v3[2])@p[1]))
            sec1=axis_rot(a2,sec1,t2)
            sec2.append(sec1)
        sol=[ translate(path[i],sec2[i]) for i in range(len(path))]

        return sol
    
def path_extrude_closed(sec,path,twist=0):
    '''
    function to extrude a closed section to a closed path
    closed path means the path provided has it's first and the last point same example a circle
    '''
    if twist==0:
        p1=path
        p2=path[1:]+[path[0]]
        p1,p2=array([p1,p2])
        v1=p2-p1
        u1=v1/norm(v1,axis=1).reshape(-1,1)
        v2=concatenate([[(u1[-1]+u1[0])/2], (u1[1:]+u1[:-1])/2])
        sec2=[]
        for i in range(len(path)):
            sec1=translate(path[i],sec2vector(v2[i],sec))
            sec2.append(sec1)
        sec2=sec2+[sec2[0]]
        # sec3=concatenate([align_sec(sec2[i-1],sec2[i]) for i in range(1,len(sec2))]).tolist()
        return sec2
        
    if twist==1:
        sec=flip(sec) if cw(sec)==-1 else sec
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

        map_v=array([t_v,n_v,o_v]).transpose(1,0,2)
        sec2=[]
        for p in map_v:
            v2=[[0,0,0],[0,0,-1],[0,1,0]]
            a1=cross(v2[1],p[0])
            t1=r2d(arccos(array(v2[1])@p[0]))
            sec1=axis_rot(a1,sec,t1)
            v3=axis_rot(a1,v2,t1)
            a2=cross(v3[2],p[1])
            t2=r2d(arccos(array(v3[2])@p[1]))
            sec1=axis_rot(a2,sec1,t2)
            sec2.append(sec1)
        sol=[ translate(path[i],sec2[i]) for i in range(len(path))]
        sol=sol+[sol[0]]
        return sol
        
def rationalise_path(path,eps=.01):
    p2=array(path)
    p_v=array([ p2[i+1]-p2[i] if i<len(p2)-1 else
               p2[i]-p2[i-1]
        for i in range(len(p2))])
    p_v=p_v/norm(p_v,axis=1).reshape(-1,1)
    p3=p2[1:][(abs(p_v[1:]-p_v[:-1])>eps).any(1)].tolist()
    p3=[p2[0].tolist()]+p3
    return p3

def cir_3p_3d(points,s=20):
    '''
    draws a circle through the 3 points list
    's' is the number of segments of the circle
    '''
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

def cp_cir_3d(cir):
    '''
    center point of circle with atleast 3 known list of 'points' in 3d space
    '''
    n1=array(nv(cir))
    a1=cross(n1,[0,0,-1])
    t1=r2d(arccos(n1@[0,0,-1]))
    sec1=translate(-array(cir).mean(0),cir)
    sec2=c3t2(axis_rot(a1,sec1,t1))
    l1=len(sec2)
    p0,p1,p2=[sec2[0],sec2[int(l1/3)],sec2[int(l1*2/3)]]
    cp=cp_3p(p0,p1,p2)
    cp=translate(array(cir).mean(0),axis_rot(a1,[cp],-t1))[0]
    return cp

def centroid_3p_3d(points):
    '''
    calculates the centroid of a triangle in 3d space
    '''
    n1=array(nv(points))
    a1=cross(n1,[0,0,-1])
    t1=r2d(arccos(n1@[0,0,-1]))
    sec1=translate(-array(points).mean(0),points)
    sec2=c3t2(axis_rot(a1,sec1,t1))
    l1=len(sec2)
    p0,p1,p2=[sec2[0],sec2[int(l1/3)],sec2[int(l1*2/3)]]
    cp=centroid_3p([p0,p1,p2])
    cp=translate(array(points).mean(0),axis_rot(a1,[cp],-t1))[0]
    return cp

def centroid_3p(points):
    '''
    calculates the centroid of a triangle in 2d
    '''
    sec1=seg(points)
    mid=[array(p).mean(0) for p in sec1]
    centroid=i_p2d([points[0],mid[1]],[points[1],mid[2]])
    return centroid

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


def i_line_planes(p1,p2):
    n1=nv(p1)
    n2=nv(p2)

    l1=cross(n1,n2)
    t1=axis_rot(l1,[n1],90)[0]
    t2=axis_rot(l1,[n2],90)[0]

    m1=array(p1).mean(0).tolist()
    m2=array(p2).mean(0).tolist()

    i_line1=((cross(n1,n2)*10+array(p1).mean(0))).tolist()
    i_line2=((cross(n1,n2)*10+array(p2).mean(0))).tolist()


    n_line1=((array(n1)*10+array(p1).mean(0))).tolist()
    n_line2=((array(n2)*10+array(p2).mean(0))).tolist()

    t_line1=((array(t1)*10+array(p1).mean(0))).tolist()
    t_line2=((array(t2)*10+array(p2).mean(0))).tolist()

    px=ppplane([m1,t_line1,m2,t_line2],l1,m1)

    i_p1=i_p3d([px[0],px[1]],[px[2],px[3]])

    l_i=(cross(n1,n2)*10+array(i_p1)).tolist()

    line1=[i_p1,l_i]

    pnts=l_sec_ip_3d(p2,line1)
    pnts1=l_sec_ip_3d(p1,line1)

    if pnts!=[] and pnts1!=[]:
        pa,pb=array(pnts),array(pnts1)
        la,lb=pb[1]-pa[0],pa[1]-pa[0]
        la=abs(la/norm(la)).round(4)
        lb=abs(lb/norm(lb)).round(4)
        d1=la.tolist()==lb.tolist()
        v1=pa[1]-pa[0]
        u1=v1/norm(v1)
        v2=pb[0]-pa[0]
        v3=pb[1]-pa[0]
        t1=u1@v2
        t2=u1@v3
        d2=(t1>=0)&(t1<=norm(v1))
        d3=(t2>=0)&(t2<=norm(v1))
        
        v1=pb[1]-pb[0]
        u1=v1/norm(v1)
        v2=pa[0]-pb[0]
        v3=pa[1]-pb[0]
        t3=u1@v2
        t4=u1@v3
        d4=(t3>=0)&(t3<=norm(v1))
        d5=(t4>=0)&(t4<=norm(v1))
        
        
        decision=d1 & ((d2 | d3) | (d4 |d5))
    else:
        decision=0

#     if decision==1:
    return pnts if decision==1 else []


def l_sec_ip(line,sec):
    l1=array(line)
    s1=array(seg(sec))
    v1=l1[1]-l1[0]
    p_l=[]
    for p in s1:
        v2=p[1]-p[0]
    #     l[0]+v1*t1=p[0]+v2*t2
    #     v1*t1-v2*t2=p[0]-l[0]
        iim=array([v1,-v2]).transpose(1,0)+.00001
        im=inv(iim)
        px=p[0]-l1[0]
        t2=(im@px)[1]
        pnts=p[0]+v2*t2
        if 0<=t2<=1:
            p_l.append(pnts.tolist())
    return p_l

def l_sec_ip_3d(sec,line):
    n1=array(nv(sec))
    a1=cross(n1,[0,0,-1])
    t1=r2d(arccos(n1@[0,0,-1]))
    sec1=translate(-array(sec).mean(0),sec)
    sec2=c3t2(axis_rot(a1,sec1,t1))
    line1=translate(-array(sec).mean(0),line)
    line2=c3t2(axis_rot(a1,line1,t1))
    l1=len(sec2)
    p0,p1,p2=[sec2[0],sec2[int(l1/3)],sec2[int(l1*2/3)]]
    pnts=l_sec_ip(line2,sec2)
    pnts=translate(array(sec).mean(0),axis_rot(a1,pnts,-t1)) if pnts!=[] else []
    return pnts

def path_offset_n(sec,r):
    sec=flip(sec) if cw(sec)==1 else sec
    r=round(r,3)
    sec1=offset_segv(sec,r)[:-1]
    s=offset_points(sec,r)[:-1]+[sec1[-1][1]]
    a=s_int1(sec1)
    if a!=[]:
        sec2=a+s
    #         sec2=pies1(sec,sec2)
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
    
def faces(l,m):
    '''
    calculate the faces for the vertices with shape l x m with first and the last end closed
    '''
    n1=arange(m,dtype=int)
    n2=array([[[[(j+1)+i*m,j+i*m,j+(i+1)*m],[(j+1)+i*m,j+(i+1)*m,(j+1)+(i+1)*m]] \
             if j<m-1 else \
             [[0+i*m,j+i*m,j+(i+1)*m],[0+i*m,j+(i+1)*m,0+(i+1)*m]] \
                 for j in range(m)] for i in range(l-1)],dtype=int).reshape(-1,3)
    n3=(array(flip(arange(m)),dtype=int)+(l-1)*m)
    n=[n1.tolist()]+n2.tolist()+[n3.tolist()]
    return n

def faces_1(l,m):
    '''
    calculate the faces for the vertices with shape l x m 
    '''
    n=array([[[[(j+1)+i*m,j+i*m,j+(i+1)*m],[(j+1)+i*m,j+(i+1)*m,(j+1)+(i+1)*m]] \
             if j<m-1 else \
             [[0+i*m,j+i*m,j+(i+1)*m],[0+i*m,j+(i+1)*m,0+(i+1)*m]] \
                 for j in range(m)] for i in range(l-1)]).reshape(-1,3).tolist()

    return n

# def faces_2(l,m):
#     '''
#     calculate the faces for the vertices with shape l x m with first and the last end open
#     '''
#     n=array([[[[(j+1)+i*m,j+i*m,j+(i+1)*m],[(j+1)+i*m,j+(i+1)*m,(j+1)+(i+1)*m]] \
#              if j<m-1 else \
#              [[0+i*m,j+i*m,j+(i+1)*m],[0+i*m,j+(i+1)*m,0+(i+1)*m]] \
#                  for j in range(m-1)] for i in range(l-1)]).reshape(-1,3).tolist()

#     return n

def faces_2(l,m):
    '''
    returns the faces for the vertices with shape l x m with first and the last end open
    '''
    return concatenate(faces(l,m)[1:-1]).tolist()

def faces_3(l,m):
    '''
    calculate the faces for the vertices with shape l x m with first and the last end open and faces flipped
    '''
    n=array([[[[(j+1)+i*m,j+i*m,j+(i+1)*m],[(j+1)+i*m,j+(i+1)*m,(j+1)+(i+1)*m]] \
             if j<m-1 else \
             [[0+i*m,j+i*m,j+(i+1)*m],[0+i*m,j+(i+1)*m,0+(i+1)*m]] \
                 for j in range(m-1)] for i in range(l-1)]).reshape(-1,3).tolist()

    return [flip(p) for p in n]

# def surface_for_fillet(sol,vector,radius,thickness,forward,right,top,num,s=50):
#     '''
#     sol: the object where the surface needs to be created
#     vector,radius,thickness,forward,right and top are the parameters to define the o_solid which covers the surface
#     num: number of slices in a surface
#     s: number of segments in the circle
#     '''
#     sol3=o_solid(vector,circle(radius,s=s),thickness,forward,right,top)
#     a=cpo([ls([array(sol3[0]).mean(0).tolist(),p],num) for p in sol3[0]])
#     b=cpo([ls([array(sol3[1]).mean(0).tolist(),p],num) for p in sol3[1]])
#     c=[[a[i]]+[b[i]] for i in range(len(a))]
#     sol4=[ip_sol2sol(sol1,p,0) for p in c if ip_sol2sol(sol1,p,0)!=[]]
#     return sol4
    
def surface_for_fillet(sol1=[],sol2=[],factor1=50,factor2=10,factor3=1,factor4=100,dia=40):
    '''
    sol1: Solid on which the surface needs to be created
    sol2: Intersecting solid
    factor1: number of segments in the circle
    factor2: number of layers or slices of surface
    factor3: decides the size of the surface lower value means bigger size. value can be set between 1 to any number
    factor4: any high number should be ok like maybe 100 or greater, basically greater than the bounding box dimension of the "sol1"
    dia: diameter around the solid 2 where surfavce needs to be created
    '''
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
    v,f1=vnf2(sol1)
    v,f1=array(v),array(f1)
    bc1=v[f1].mean(1)
    f2=f1[(sqrt((bc1[:,0]-p1[0])**2+(bc1[:,1]-p1[1])**2+(bc1[:,2]-p1[2])**2)<=dia)]
    solx=v[f2].tolist()
    sur2=[ipx(solx,p) for p in sol3][1:]
    
    return sur2

    
def prism_center(sol):
    '''
    calculates the center of the prism or solid object, may not be the mean.
    This calculates the center of the bounding box for a solid.
    
    '''
    x_max=array(sol)[:,:,0].max()
    x_min=array(sol)[:,:,0].min()

    y_max=array(sol)[:,:,1].max()
    y_min=array(sol)[:,:,1].min()

    z_max=array(sol)[:,:,2].max()
    z_min=array(sol)[:,:,2].min()

    return [array([x_max,x_min]).mean(),array([y_max,y_min]).mean(),array([z_max,z_min]).mean()]
    
    
def ipx(prism,prism1,side=-1):
    '''
    function to calculate intersection point between two 3d prisms. 
     "prism" is the 3d object which is intersected with "prism1".
     side: when a ray intersects a solid it can intersect at 2 locations, if the ray is travelling from outside, in that case if '0' is given meaning only the first intersection point is considered, and in case '-1' is given meaning the last intersection point will be considered.
     try below code for better understanding:
    sec=circle(10)
    path=cr(pts1([[2,0],[-2,0,2],[0,10,3],[-9.9,0]]),5)
    p=prism(sec,path)
    p1=cylinder(r=3,h=15,s=30)
    ip1=ip(p,p1)
    
    refer to file "example of various functions" for application
    '''
    pb=prism1
    p1=array(prism)
    p2=array([[[pb[i][j],pb[i+1][j]] for j in range(len(pb[i]))] for i in range(len(pb)-1)]).reshape(-1,2,3)
    pm=p1[:,0]
    pn=p1[:,1]
    po=p1[:,2]
    px=p2[:,0]
    py=p2[:,1]
    v1,v2,v3=py-px,pn-pm,po-pm
    t1=einsum('ijk,jk->ij',px[:,None]-pm,cross(v2,v3))/einsum('ik,jk->ij',-v1,cross(v2,v3)+[.00001,.00001,.00001])
    t2=einsum('ijk,ijk->ij',px[:,None]-pm,cross(v3,-v1[:,None]))/einsum('ik,jk->ij',-v1,cross(v2,v3)+[.00001,.00001,.00001])
    t3=einsum('ijk,ijk->ij',px[:,None]-pm,cross(-v1[:,None],v2))/einsum('ik,jk->ij',-v1,cross(v2,v3)+[.00001,.00001,.00001])
    p=px[:,None]+einsum('ik,ij->ijk',v1,t1)
    condition=(t1>=0)&(t1<=1)&(t2>=0)&(t2<=1)&(t3>=0)&(t3<=1)&((t2+t3)>=0)&((t2+t3)<=1)
#     p=p[condition]
#     p=p[unique(p,return_index=True)[1]]
    p=array([[p[i][condition[i]],i%len(pb[0])] for i in range(len(p))],dtype=object)
    n=array([p1[1] for p1 in p if p1[0].tolist()!=[]])
    p=concatenate([concatenate([[[p2,i] for p2 in p1[0]] for p1 in p if (p1[0].tolist()!=[])&(p1[1]==i)],dtype=object) for i in n],dtype=object)
    p=array([p]*len(n),dtype=object)[p[:,1]==unique(p[:,1])[:,None]]
    p=array([[array([p1[0] for p1 in p if p1[1]==i],dtype=object),i] for i in n],dtype=object)
    if side=='all':
        p=concatenate([a[array([l_len([p2[:,0][b],p1]) for p1 in a],dtype=object).argsort()] for (a,b) in p ],dtype=object)
    else:
        p=array([a[array([l_len([p2[:,0][b],p1]) for p1 in a],dtype=object).argsort()[side]] for (a,b) in p if a.tolist()!=[]],dtype=object)
        
    return p.tolist()
    
def ipx_sol2sol(sol,sol1,i=0):
    '''
    function to find the intersection point between 2 solids
    this function is to be used where the cutting lines of sol1 are intersecting sol at more than 1 times.
    sol: solid which is intersected
    sol1: this intersects the solid 'sol'
    i: if the first intersection points of all the cutting lines are to be considered, value of i should be '0'.
    if the last intersection points of all the cutting lines are to be considered, value of 'i' should be set to '-1'
    if all the intersection points are required, value of 'i' should be set to 'all'
    '''
    if i=='all':
        a=[ipx_sol2line(sol,p) for p in cpo(sol1) if ipx_sol2line(sol,p)!=[]]
    else:
        a=[ipx_sol2line(sol,p)[i] for p in cpo(sol1) if ipx_sol2line(sol,p)!=[]]
    
    return a
    
def ipx_sol2line(sol,line):# when line has more than 2 points
    '''
    function to calculate intersection point between a 3d solid and a line. 
     "sol" is the 3d tri-mesh which is intersected with a "line".
     try below code for better understanding:
    sec=circle(10)
    path=cr(pts1([[-10+.1,0],[12,0],[-2,0,2],[0,10,3],[-10,0]]),5)
    sol=prism(sec,path)

    line=[[0,0,-1],[20,20,10]]

    ip1=ip_sol2line(sol,line)
    
    refer to file "example of various functions" for application
    '''


    pa=sol
    p1=array(sol)
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
    return i_p2
    
def partial_surface(sol,cp1,dia):
    '''
    return the part of the surface of a given solid 'sol' in terms of vertices and faces
    cp1: center point around which the surface is needed
    dia: dia of the sphere inside which the surface needs tobe extracted
    '''
    v,f1=vnf2(sol)
    v,f1=array(v),array(f1)
    bc1=v[f1].mean(1)
    f2=f1[(sqrt((bc1[:,0]-cp1[0])**2+(bc1[:,1]-cp1[1])**2+(bc1[:,2]-cp1[2])**2)<=dia)]
    return [v.tolist(),f2.tolist()]
    
def align_sol_1(sol):
    '''
    function to straighten the twists in the path_extruded sections for better alignments
    refer to the file "example of various functions.ipynb" for application examples
    '''
    sol1=[sol[0]]
    for i in range(1,len(sol)):
        a=align_sec_1(sol1[i-1],sol[i])
        sol1.append(a[1])
    return sol1
    
def align_sec_1(sec1,sec2):
    '''
    function to align 2 3d sections to obtain the non twisted optimised solid
    refer file "examples of various functions" for application examples
    '''
    
    area1=[ norm(array(sec2[i:]+sec2[:i])-array(sec1),axis=1).sum() for i in range(len(sec2)) ]
    i=array(area1).argmin()
    sol2=[sec1]+[sec2[i:]+sec2[:i]]
    return sol2

def convert_3lines2fillet_closed(pnt1,pnt2,pnt3,f=1.9,s=10):
    '''
    Develops a fillet with 3 list of points (pnt1,pnt2,pnt3) in 3d space
    f: is a factor which can be reduced to 1.5 in case of self intersection observed
    s: number of segments in the fillet, increase the segments in case finer model is required
    refer to the file "example of various functions" for application examples
    
    '''
    m1=min(len(pnt1),len(pnt2),len(pnt3))
    sol=array([pnt3[:m1],pnt1[:m1],pnt2[:m1]]).transpose(1,0,2)
    sol=[fillet_3p_3d(p3,p2,p1,r_3p_3d([p1,p2,p3])*f,s) for (p1,p2,p3) in sol]
    sol=sol+[sol[0]]
    return sol

def ip_s2l(sol,line,side=-1):
    '''
    function to calculate intersection point between solid to line. 
     "sol" is the 3d object which is intersected with "line".
     side: when a ray intersects a solid it can intersect at 2 locations, if the ray is travelling from outside, in that case if '0' is given meaning only the first intersection point is considered, and in case '-1' is given meaning the last intersection point will be considered.
     try below code for better understanding:
    refer to file "example of various functions" for application
    '''
    pa=sol
    p1=array([[ [[pa[i][j],pa[i][j+1],pa[i+1][j]],[pa[i+1][j+1],pa[i+1][j],pa[i][j+1]]] if j<len(pa[i])-1 
     else [[pa[i][j],pa[i][0],pa[i+1][j]],[pa[i+1][0],pa[i+1][j],pa[i][0]]] 
     for j in range(len(pa[i]))] 
              for i in range(len(pa)-1)]).reshape(-1,3,3)
    p2=array(seg(line)[:-1])
    pm=p1[:,0]
    pn=p1[:,1]
    po=p1[:,2]
    px=p2[:,0]
    py=p2[:,1]
    # px+v1*t1=pm+v2*t2+v3*t3
    v1,v2,v3=py-px,pn-pm,po-pm

    iim=array([array([v1]*len(v2)).transpose(1,0,2),array([-v2]*len(v1)),array([-v3]*len(v1))]).transpose(1,2,0,3)
    im=inv(iim+.0001)
    p=array([pm]*len(v1))-array([px]*len(v2)).transpose(1,0,2)
    t=einsum('ijk,ijkl->ijl',p,im)
    condition=(t[:,:,0]>=0)&(t[:,:,0]<=1)&(t[:,:,1]>=0)&(t[:,:,1]<=1)&(t[:,:,2]>=0)&(t[:,:,2]<=1)&(t[:,:,1]+t[:,:,2]<=1)

    p=array([px]*len(v2)).transpose(1,0,2)[condition]+array([v1]*len(v2)).transpose(1,0,2)[condition]*t[condition][:,0][:,None]
    p=p[norm(p-px[0],axis=1).argsort()]
    return p.tolist()
    
def gcd(a,b):
    '''
    calculates the greatest common divisor of 2 numbers 'a','b'
    '''
    for _ in range(max(a,b)):
        if a>b:
            a=a-b
        elif b>a:
            b=b-a
        elif a==b:
            break
    return a

def lcm(a,b):
    '''
    calculates the least common multiple of 2 numbers 'a','b'
    '''
    return a*b/gcd(a,b)

def perp_min_dist_point(line,points):
    '''
    out of all the "points" in the list, first this function selects points which have a perpendicular projection on the line.
    subsequently select the point which is shortest distance from the line.
    In case no point is projected on the line, function returns an empty list.
    '''
    p0=line[0]
    p1=line[1]
    p0,p1=array([p0,p1])
    v1=p1-p0
    u1=v1/(norm(v1)+.00001)
    v2=array(points)-p0
    v2cost=einsum('j,ij->i',u1,v2)
    cond=(v2cost>=0)&(v2cost<=l_len(line))
    pnts=array(points)[cond]
    if pnts.tolist()!=[]:
        dist=array(perp_dist(line,pnts)).argmin()
        pnts=pnts[dist]
        return pnts.tolist()
    else:
        return []

def axis_rot_o(axis,solid,angle):
    '''
    rotate a solid around an axis considering the solid is centered at origin
    '''
    cp1=array(prism_center(solid))
    solid=translate(-cp1,solid)
    return translate(cp1,[[q(axis,p1,angle) for p1 in p] for p in solid])

def edges(l,m):
    return array([[[i+j*m,i+(j+1)*m] for j in range(l-1)] for i in range(m)]).reshape(-1,2)
    
# def i_p_n(sol,i_p):
#     '''
#     calculates normal at the intersection points
#     sol: solid on which the normal is required
#     i_p: list of intersection points between 2 solids
#     '''
#     l,m,_=array(sol).shape
#     f1=faces(l,m)
#     v=array(sol).reshape(-1,3)
#     tri=v[f1[1]]
#     v2,v3=tri[:,1]-tri[:,0], tri[:,2]-tri[:,0]
#     v1=cross(v2,v3)
#     v1=v1/norm(v1,axis=1).reshape(-1,1)
#     p0=array(i_p)
#     p2=tri[:,0]
#     n1,n2=len(p0),len(p2)
#     v1=array([v1]*n1)
#     v2=array([v2]*n1)
#     v3=array([v3]*n1)
#     p0=array([p0]*n2).transpose(1,0,2)
#     p2=array([p2]*n1)
#     iim=array([v1,-v2,-v3]).transpose(1,2,0,3).transpose(0,1,3,2)
#     im=inv(iim)
#     t=einsum('ijkl,ijl->ijk',im,p2-p0)
#     d1=(t[:,:,0]>-.1)&(t[:,:,0]<.1)&(t[:,:,1]>=0)&(t[:,:,1]<=1) \
#     &(t[:,:,2]>=0)&(t[:,:,2]<=1)&(t[:,:,1]+t[:,:,2]<1)
#     nx=v1[d1]
#     return nx

def i_p_n(sol,ip):
    '''
    calculates the unit normal vectors for the intersection points "ip" on a solid "sol"
    '''
    sec2=[ip_triangle(sol,p) for p in ip]
    pa,pb,pc=array(sec2)[:,0],array(sec2)[:,1],array(sec2)[:,2]
    n1=cross(pc-pa,pb-pa)
    n1=n1/norm(n1,axis=1).reshape(-1,1)
    return n1




def i_p_p(sol,i_p,r):
    '''
    function to project the intersection point on the cutting lines based on the distance 'r'
    '''
    sol=array(sol)
    i_p=array(i_p)
    s1=[]
    for j in range(len(sol[0])):
        for i in range(len(sol)-1):
            v1=sol[i+1][j]-sol[i][j]
            u1=v1/norm(v1)
            l1=norm(v1)
            v2=i_p[j]-sol[i][j]
            v2cost=u1@v2
            if v2cost>=0 and v2cost<=l1:
                l2=[i_p[j].tolist()]+array(cpo(sol)[j])[arange(i+1,len(sol))].tolist()
                len_l2=path_length(l2)
                s=len_l2/r
                p1=equidistant_path(l2,s)[1]
                s1.append(p1)
    return s1

def path_length(path):
    '''
    calculates the length of the path
    '''
    v=[p[1]-p[0] for p in array(seg(path)[:-1])]
    l=[l_len(p) for p in seg(path)[:-1]]
    c=array(l).cumsum()
    
    return c[-1]

def i_p_t(path):
    '''
    function to calculate tangent vectors to a given path
    '''
    p1=array(seg(path))
    p2=array(path)
    v1=array([(p[1]-p[0])/norm(p[1]-p[0]) for p in p1])
    t_v=array([ (v1[-1]+v1[i])/2 if i==0 else
         (v1[i-1]+v1[i])/2
        for i in range(len(p1))])

    t_v=t_v/norm(t_v,axis=1).reshape(-1,1)

    return t_v

def o_p_p(sol,i_p,d):
    '''
    calculates projected points on the surface of a solid 
    sol: solid on which the points to be projected
    i_p: list of points in 3d space near the solid
    d: approximate distance of the points from the surface, specifying too big distance 
    may create multiple projection of the same point on the solid
    '''
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
    
def trim_points(o,pl,r):
    px=array(seg(pl))
    o=array(o)
    v1=px[:,1]-px[:,0]
    v1=array([v1]*len(o))
    v1norm=einsum('ijk,ijk->ij',v1,v1)**0.5
    u1=einsum('ijk,ij->ijk',v1,1/v1norm)
    v2=o[:,None]-px[:,0]
    v2cost=einsum('ijk,ijk->ij',u1,v2)
    v2sint=norm(cross(v1,v2),axis=2)/norm(v1,axis=2)
    d1=(v2cost>=0)&(v2cost<=norm(v1,axis=2))&(v2sint<r)
    pnts1=array([o]*len(px)).transpose(1,0,2)[d1]
    pnts1=exclude_points(o,pnts1)
    return pnts1

def arc_2p_3d(n1,p0,p1,r,cw=1,s=20):
    '''
    draws an arc through 2 points 
    n1: normal vector to define plane on which the arc will be drawn
    r: radius of the arc
    cw: '1' stands for clockwise and '-1'stands for counter-clockwise
    's' is the number of segments of the circle
    '''
    n1=array(n1)
    a1=cross(n1,[0,0,-1])
    t1=r2d(arccos(n1@[0,0,-1]))
    sec1=translate(-array([p0,p1]).mean(0),[p0,p1])
    sec2=c3t2(axis_rot(a1,sec1,t1))
    pa,pb=sec2
    arc1=arc_2p(pa,pb,r,cw,s=s)
    arc1=translate(array([p0,p1]).mean(0),axis_rot(a1,arc1,-t1))
    return arc1

def arc_long_2p_3d(n1,p0,p1,r,cw=1,s=20):
    '''
    draws a long arc through 2 points 
    n1: normal vector to define plane on which the arc will be drawn
    r: radius of the arc
    cw: '1' stands for clockwise and '-1'stands for counter-clockwise
    's' is the number of segments of the circle
    '''
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
    '''
    calculates the center point of the circle drawn through 2 points 
    n1: normal vector to define plane on which the arc/ circle drawn
    r: radius of the arc/ circle
    cw: '1' stands for clockwise and '-1'stands for counter-clockwise
    '''
    n1=array(n1)
    a1=cross(n1,[0,0,-1])
    t1=r2d(arccos(n1@[0,0,-1]))
    sec1=translate(-array([p0,p1]).mean(0),[p0,p1])
    sec2=c3t2(axis_rot(a1,sec1,t1))
    pa,pb=sec2
    cp=arc_2p_cp(pa,pb,r,cw)
    cp=translate(array([p0,p1]).mean(0),axis_rot(a1,[cp],-t1))[0]
    return cp

def sec_lines_ip2d(sec,lines):
    # p0+v1*t1=p2+v2*t2
    px=array(lines)
    py=array(seg(sec))
    p0=px[:,0]
    v1=px[:,1]-px[:,0]
    p2=py[:,0]
    v2=py[:,1]-py[:,0]
    n1,n2=len(p0),len(p2)
    p0=array([p0]*n2).transpose(1,0,2)
    p2=array([p2]*n1)
    v1=array([v1]*n2).transpose(1,0,2)
    v2=array([v2]*n1)

    iim=array([v1,-v2]).transpose(1,2,0,3).transpose(0,1,3,2)+.0001
    im=inv(iim)
    t=einsum('ijkl,ijl->ijk',im,p2-p0)
    d=(t[:,:,0]>=0)&(t[:,:,0]<=1)&(t[:,:,1]>=0)&(t[:,:,1]<=1)
    i_p1=p0[d]+einsum('ij,i->ij',v1[d],t[:,:,0][d])
    return i_p1.tolist()
    
def aligned_cut_lines_prism(sec,path,s=100):
    sec_list=[offset(sec,x) for (x,y) in path]
    i=array([path_length(p) for p in sec_list]).argmin()
    j=array([path_length(p) for p in sec_list]).argmax()

    cp=array(sec_list[i]).mean(0)
    r=max([l_len([cp,p]) for p in sec_list[j]])+10
    c1=circle(.1,cp,s)
    c2=circle(r,cp,s)
    p1=cpo([c1,c2])
    sec_list_r=[sec_lines_ip2d(p,p1) for p in sec_list]
    sol3=[translate([0,0,path[i][1]],sec_list_r[i]) for i in range(len(path))]
    return sol3

def ch1(p_l,k=3):
    p_l=array(p_l)
    s_p=p_l[p_l[:,1].argmin()]
    p_l=array(exclude_points(p_l,[s_p]))
    p_l1=p_l[norm(p_l-s_p,axis=1).argsort()[:k]]

    n_p=p_l1[array([ang_2linecw(s_p,s_p-[0,1],p) for p in p_l1]).argmax()]

    p_l=array(exclude_points(p_l,[n_p])+[s_p])
    s1=[s_p,n_p]

    while((array(s1[-1]).tolist()!=s_p.tolist())&(len(p_l)>3)):
        p_l1=p_l[norm(p_l-s1[-1],axis=1).argsort()[:k]]
        n_p=p_l1[array([-ang_2linecw(s1[-1],s1[-2],p) for p in p_l1]).argsort()]
        for p in n_p:
            if s_int1(seg(array(s1).tolist()+[p.tolist()])[:-1])==[]:
                s1.append(p)
                p_l=array(exclude_points(p_l,[p]))
                break
            else:
                continue
            break
        else:
            return array(s1).tolist()
        
        if len(p_l)==3:
            n_p=p_l[array([ang_2linecw(s1[-1],s1[-2],p) for p in p_l]).argsort()]
            s1=array(s1).tolist()+array(n_p).tolist()
    return array(s1).tolist()


def concave_hull(p_l,k=3):
    '''
    calculate the concave hull of a random list of points
    larger number for 'k' will give smoother shape
    '''
    s1=ch1(p_l,k)
    py=exclude_points(p_l,s1)
    x=pies1(s1,py)
    pz=exclude_points(py,x) if x!=[] else py

    while (pz!=[]):
        k=k+1
        s1=ch1(p_l,k)
        py=exclude_points(p_l,s1)
        x=pies1(s1,py)
        pz=exclude_points(py,x) if x!=[] else py

    return s1
    
def axis_rot_1(sol,ax1,loc1,theta):
    '''
    rotate a solid on any pivot point 'loc1' with axis of rotation 'ax1' by an angle 'theta'
    
    '''
    c1=array(prism_center(sol))
    loc1=array(loc1)
    c2=c1-loc1
    s1=translate(-c1,sol)
    s1=translate(c2,s1)
    s2=axis_rot(ax1,s1,theta)
    s2=translate(-c2,s2)
    s2=translate(c1,s2)
    return s2

def oset(sec,r):
    '''
    Simpler algorithm for offset
    Takes a little longer to compute
    '''
    s1=offset_segv(sec,r)
    s2=min_d_points(intersections(s1),.01)
    s3=s_int1(seg(s2))
    if s3==[]:
        s4=s2
    else:

        s4=s2+s3
        c1=cs1(sec,abs(r))
        for p in c1:
            p1=pies1(p,s4)
            if p1!=[]:
                s4=exclude_points(s4,p1)
    s4=sort_points(sec,s4)
    return s4

def pa2pb(path,zval):
    '''
    function to match the path to extrude points to the z-values of a solid to extrude along the path
    '''
    v=[p[1]-p[0] for p in array(seg(path)[:-1])]
    l=[l_len(p) for p in seg(path)[:-1]]
    c=array(l).cumsum().tolist()
    d=[l_lenv_o(path)/zval[-1]*p for p in zval[1:-1]]
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
    
def path2path1(path1,path):
    '''
    function to match the points of path1 with path
    i.e. path1 is independent variable and path is dependent variable
    '''
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

def spline_curve(p0,n=10,f=1):
    '''
    cubic spline curve based on the bezier considering 3 points 
    p0: control points
    n: number of points between each 3 points
    f: is a factor which an integer, lower values will make the spline sharp and higher values will be more curvy
    example:
p0 = pts2([[0,0,0],[10,0,2],[0,10,3],[5,-3,-6],[-5,5,1],[0,7,-2],[0,0,5]])
    spline_curve(p0) => will give the spline curve
    
    '''
    p0=[p0[0]]*3+p0

    for i in range(1,len(p0)):
        if i==1:
            
            c=bezier(p0[:3],n)[:-f]
            s1=c
        elif i==len(p0)-2:
            p1=[s1[-1]]+p0[i:len(p0)]
            c=bezier(p1,n)
            s1=s1+c
        else:
            p1=[s1[-1]]+p0[i:i+2]
            c=bezier(p1,n)[:-f]
            s1=s1+c
    return [p0[0]]+exclude_points(s1,p0[0])
    
def ip_triangle(sol1,p0):
    '''
    finds the triangle where the intersection point lies in a solid
    '''
    l,m,_=array(sol1).shape
    f1=faces(l,m)[1:-1]
    v=array(sol1).reshape(-1,3)
    tri=v[f1]
    pa,pb,pc=tri[:,0],tri[:,1],tri[:,2]
    v1,v2=pb-pa,pc-pa
    v0=cross(v1,v2)
#     v0=array([nv(p) for p in tri])
    tri=tri[~(v0==[0,0,0]).all(1)]
#     v0=array([nv(p) for p in tri])
#     p0=array(p3[89])
    p0=array(p0)
    pa,pb,pc=tri[:,0],tri[:,1],tri[:,2]
    v1,v2=pb-pa,pc-pa
    v0=cross(v1,v2)
    iim=array([v0,-v1,-v2]).transpose(1,0,2)
    im=inv(iim)
    p=pa-p0
    t=einsum('ijk,ij->ik',im,p)
    d=(t[:,0]>=-0.01)&(t[:,0]<=1)&(t[:,1]>=0)&(t[:,1]<=1)&(t[:,2]>=0)&(t[:,2]<=1)&(t[:,1]+t[:,2]<=1)
    sec2=tri[d].tolist()
    return sec2[0] if sec2!=[] else []

