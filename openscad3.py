from openscad2 import *

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

def turtle2D(p):
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


    
def check_closed_section_orientation(sec):
    '''
    function to identify if an enclosed section is clockwise(cw) or counterclockwise(ccw)
    this returns 1 if section is clockwise and -1 if it is counterclockwise
    '''
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


def check_closed_section_orientation_of_each_point(sec):
    '''
    function to identify whether each point in a section is clockwise or counter clockwise. 
    checkClosedSectionOrientation(sec)==1 means clockwise and -1 means counterclockwise. 
    e.g.
    checkClosedSectionOrientationAtEachPoint(pts([[0,0],[4,0],[0,4],[2,0],[0,2],[-6,0]])) => [-1, -1, 1, -1, -1, -1]
    '''
    p=sec
    p0=[p[-1]]+p[:-1]
    p1=p
    p2=p[1:]+[p[0]]
    p0,p1,p2=array([p0,p1,p2])
    p=array([p0,p1,p2]).transpose(1,0,2)
    return [ -1 if cross(c23(p1[1]-p1[0]),c23(p1[2]-p1[0]))[-1]>=0 else 1 for p1 in p]


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


def unit_vector(v):
    '''
    function to calculate unit vector of a given vector
    example:
    vector=[2,3,5]
    unit_vector=unitVector(vector) => [0.3244428422615251, 0.48666426339228763, 0.8111071056538127]
    '''
    v=array(v)
    return (v/norm(v)).tolist()


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
    

def max_radius(sec):
    '''
    function calculates the maximum radius in a given closed section
    example:
    sec=cr_c(pts1([[0,0,.2],[8,3,3],[5,7,1],[-8,0,2],[-5,20,1]]),20)
    maxRadius(sec) => 3.0
    
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
        

def offset_line(l,d):
    u=uv(subtract(l[1],l[0]))
    p0=add(l[0],dot(u,multiply(d,rm(-90)))).tolist()
    p1=add(l[1],dot(u,multiply(d,rm(-90)))).tolist()
    return [p0,p1]

def convert_to_segments(sec):
    '''
    function to create a segment from a list of points or a list
    example:
    list=[1,2,3,4,5,6]
    seg(list)=> [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 1]]
    
    list=[[1,2,3],[4,5,6],[7,8,9]]
    convertToSegments(list) => [[[1, 2, 3], [4, 5, 6]], [[4, 5, 6], [7, 8, 9]], [[7, 8, 9], [1, 2, 3]]]
    '''
    c=[]
    for i in range(len(sec)):
        i_plus=i+1 if i<len(sec)-1 else 0
        p0=sec[i]
        p1=sec[i_plus]
        l=[p0,p1]
        c.append(l)
    return c


def offset_segments(sec,d):
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
    return array(offset_segments(sec,r))[:,0].tolist()


def offset_segments_with_CW_orientation(sec,r):
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


def remove_duplicates(points_list):
    '''
    function removes all the duplicates from a 
    example:
    list=[9,5,1,2,3,4,5,2,4,6,9]
    remove_extra_points(list) => [9, 5, 1, 2, 3, 4, 6]
    
    list=[[7,8,9],[1,2,3],[10,11,12],[4,5,6],[7,8,9],[10,11,12]]
    remove_extra_points(list) => [[7, 8, 9], [1, 2, 3], [10, 11, 12], [4, 5, 6]]
    '''
    return array(points_list)[sort(unique(points_list,axis=0,return_index=True)[1])].tolist()

def remove_radiuses_from_sec_with_CCW_orientation_points(sec):
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
        clock=array([-1 if cross(c23(p[1]-p[0]),c23(p[2]-p[0]))[-1]>0 else 1 for p in arr])
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


def remove_all_radiuses_from_closed_section(sec):
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
        clock=array([-1 if cross(c23(p[1]-p[0]),c23(p[2]-p[0]))[-1]>0 else 1 for p in arr])
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


def remove_radiuses_from_section_with_CW_orientation_points(sec,d):
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
        clock=array([-1 if cross(c23(p[1]-p[0]),c23(p[2]-p[0]))[-1]>0 else 1 for p in arr])
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



def list_radiuses_of_closed_section(sec):
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

def max_radius(sec):
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
    l1,l2=array([l1,l2])
    v1=l1[1]-l1[0]
    v2=l2[1]-l2[0]
    t1,t2=inv(array([v1,-v2]).transpose(1,0))@(l2[0]-l1[0])
    return (l1[0]+v1*t1).tolist()


def self_intersections(sec1):
    '''
    calulates the self intersection points of a list of line segments 's'
    it also picks the points in case the 2 lines are just connected at 1 point and are not crossing
     
    refer to file 'example of various functions' for application example
    '''
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



def multiple_points_at_defined_pitch_closed_section(sec,pitch=2):
    '''
    multiple points within straight lines of a closed section 'sec' with equal segment length 'sl' in the straight line segments
    refer file "example of various functions" for application example
    '''
    sec1=[]
    for p in segments(sec):
        n=1 if int(round(l_len(p)/pitch,0))==0 else int(round(l_len(p)/pitch,0))
        sec1.append(ls(p,n))
    return concatenate(sec1).tolist()
    

def multiple_points_at_defined_pitch_open_section(sec,pitch=2):
    '''
    multiple points within straight lines of a open section 'sec' with equal segment length 'sl' in the straight line segments
    refer file "example of various functions" for application example
    '''
    sec1=[]
    for p in seg(sec)[:-1]:
        n=1 if int(round(l_len(p)/pitch,0))==0 else int(round(l_len(p)/pitch,0))
        sec1.append(ls(p,n))
    return concatenate(sec1).tolist()+[sec[len(sec)-1]]


def line_slices(line,n):
    '''
    function to draw number of points 'n' in a line 'line'
    example:
    line=[[0,0],[10,0]]
    line1=ls(line,5) => [[0.0, 0.0], [2.0, 0.0], [4.0, 0.0], [6.0, 0.0], [8.0, 0.0], [10.0, 0.0]]
    '''
    p0,p1=array(line)
    v1=p1-p0
    return array([p0+v1/n*i for i in range(n)]).tolist()


def line_length(l):
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




def center_arc_2p(p1,p2,r,cw=-1):
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

def offset_2(sec,r):
    '''
    calculates offset for a section 'sec' by amount 'r'
    refer file "example of various functions" for application example
    '''

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
    '''
function to make a prism with combination of 2d section and 2d path
type: can be set to either "1" or "2" which means prism calculation is done by offset method "1" or offset method "2", by defualt it is "1"
Example:
sec=circle(10)
path=cr(pts1([[2,0],[-2,0,2],[0,10,3],[-3,0]]),5)
sol=prism(sec,path)
    '''
    s1=flip(sec) if check_closed_section_orientation(sec)==1 else sec
    return [array(translate([0,0,y],offset(s1,x,type))).tolist() for (x,y) in path]


def translate(p,sec):#translates a prism or section by [x,y,z] distance
    '''
    function to translate a group of points "sec" by "p" distance defined in [x,y,z].e.g. try following code:
    sec=cr([[0,0,.5],[10,0,2],[7,15,1]],5)
    sec1=translate(p=[2,5,3],sec=sec)
    
    refer to file "example of various functions " for application
    '''
    return (array(c23(sec))+c23(p)).tolist()

def translate_2d(p,sec):#translates a 2d section by [x,y] distance
    '''
    function to translate a group of points "sec" by "p" distance defined in [x,y].e.g. try following code:
    sec=corner_radius([[0,0,.5],[10,0,2],[7,15,1]],5)
    sec1=translate_2d(p=[2,5],sec=sec)
    
    refer to file "example of various functions " for application
    '''
    return c32((array(c23(sec))+c23(p)))

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


def change_points_orientation(prism): # changes the orientation of points of a prism
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
    return (array(p)@[[1,0,0],[0,1,0]]).tolist() if array(p).shape[-1]==2 else array(p).tolist()

def c3t2(a): # converts 3d list to 2d list 
    '''
    function to convert 3d to 2d, it just removes the z-coordinate from the points list 
    example:
    list=c3t2([[1,2,3],[3,4,5],[6,7,8]])
    output=> [[1, 2], [3, 4], [6, 7]]
    '''
    return (array(a)@[[1,0],[0,1],[0,0]]).tolist() if array(a).shape[-1]==3 else array(a).tolist()

def unit_normal_vector(p):# normal vector to the plane 'p' with atleast 3 known points
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
    
def normal_vector(p):# normal vector to the plane 'p' with atleast 3 known points
    '''
    given 3 points ['p1','p2',p3] function calculates normal vector
    example:
    p1,p2,p3=[1,0,0],[0,10,0],[-5,0,0]
    nv1([p1,p2,p3]) => [0.0, 0.0, -60.0]
    '''
    l1=len(p)
    p0,p1,p2=array(translate([0,0,0],[p[0],p[int(l1/3)],p[int(l1*2/3)]]))
    nv=cross(p0-p1,p2-p1)
    return nv.tolist()

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
    cp=p1+axis_rot(n,u1*r/cos(theta*pi/180),alpha/2)
    pa=p1+u1*l
    arc=[ cp+axis_rot(n,pa-cp,-i) for i in linspace(0,theta*2,s)]
    a,b,c=arc[0],arc[1:s-1],arc[s-1]
    return concatenate([[p1],arc]).tolist()

def cp_fillet_3p_3d(p0,p1,p2,r):# center point 'cp' of the fillet with 3 known points 'p0,p1,p2' in 3d space. 'r' is the radius of fillet
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
    cp=p1+axis_rot(n,u1*r/cos(theta*pi/180),alpha/2)
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



def radius_3p_3d(points):
    '''
    calculates the radius of circle made by 3 points in 3d space
    '''
    n1=array(nv(points))+[.000001,.000001,0]
    a1=cross(n1,[0,0,-1])
    t1=r2d(arccos(n1@[0,0,-1]))
    sec1=translate(-array(points).mean(0),points)
    sec2=c3t2(axis_rot(a1,sec1,t1))
    l1=len(sec2)
    p0,p1,p2=[sec2[0],sec2[int(l1/3)],sec2[int(l1*2/3)]]
    
    return radius_3p([p0,p1,p2])

def scale2d(sec,sl):# scale the 2d section 'sec' by a scaling factor 'sl'. this places the scaled section in the bottom center of the original section
    '''
    function to scale a 2d section by an amount "sl" which has to be >0 (keeps the y-coordinates same). 
    
    refer file "example of various functions" for application
    '''
    s1=array(translate([0,0,0],sec))
    cp=array(s1).mean(axis=0)
    rev=array(s1).mean(axis=0)+(array(s1)-array(s1).mean(axis=0))*sl
    y1=cp-array([0,array(s1)[:,1].min(),0])
    y2=cp-array([0,rev[:,1].min(),0])
    d=y2-y1
    return c32(translate(d,rev))

def scale2d_centered(sec,sl):# scale the 2d section 'sec' with scaling factor 'sl'. this places the scaled section in the center of original section or the center of both original and scaled section remains the same.
    '''
    function to scale a 2d section by an amount "sl" which has to be >0 (keeps the revised section in center). 
    
    refer file "example of various functions" for application
    '''
    s1=array(translate([0,0,0],sec))
    cp=array(s1).mean(axis=0)
    rev=array(s1).mean(axis=0)+(array(s1)-array(s1).mean(axis=0))*sl
    return c32(rev)

def scale3d(p,s):# scale 3d prism 'p' with scaling factor 's'. This places the scaled prism at the same bottom of the original prism
    '''
    function to scale a 3d prism keeping the base z-coordinate same. 
    takes 2 arguments "p" to scale and the scaling factor "s". 
    scale factor can take any real number negative values will scale the prism and turn the prism upside down.
    refer file "example of various functions" for application
    '''
    p=array(p)
    cp=p.reshape(-1,3).mean(axis=0)
    rev=cp+(p-cp)*s
    z1=p.reshape(-1,3)[:,2].min()
    z2=rev.reshape(-1,3)[:,2].min()
    d=z1-z2
    return translate([0,0,d],rev)

def scale3d_centered(p,s):# scale a 3d prism 'p' with scaling factor 's'. This places the scaled prism in the center of the original prism or the center of both the prism is same
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


    
def multiple_points_in_closed_sec(sec,s,d=.25):# multiple points with in the straight lines in the closed section 'sec'. 's' is the number of segments between each straight line
    '''
    adds 's' number of points in each straight line segment of a section 'sec'
    points will be added only if the segmentlength > 'd'
    refer to the file "example of various functions" for application example
    '''
    c=[]
    for i in range(len(sec)):
        i_plus=i+1 if i<len(sec)-1 else 0
        if l_len([sec[i],sec[i_plus]])>=d:
            c.append(ls([sec[i],sec[i_plus]],s))
        else:
            c.append([sec[i],sec[i_plus]])
    return remove_duplicates(concatenate(c))


def multiple_points_in_open_sec(sec,s,d=.25):# multiple points with in the straight lines in the open section 'sec'. 's' is the number of segments between each straight line
    '''
    adds 's' number of points in each straight line segment of an open section 'sec'
    points will be added only if the segmentlength > 'd'
    refer to the file "example of various functions" for application example
    '''
    c=[]
    for i in range(len(sec)-1):
        i_plus=i+1 if i<len(sec)-1 else 0
        if l_len([sec[i],sec[i_plus]])>=d:
            c.append(ls([sec[i],sec[i_plus]],s))
        else:
            c.append([sec[i],sec[i_plus]])
    return remove_extra_points(concatenate(c))+[sec[-1]]


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



def fillet_between_2circles(circle1,circle2,r,s=50): # fillet between 2 circles and 'r' is the radius of the fillet
    '''
    fillet between 2 circles and 'r' is the radius of the fillet
    refer to file "examples of various functions"
   
    '''
    r1=radius_arc(circle1)
    r2=radius_arc(circle2)
    c1=cp_arc(circle1)
    c2=cp_arc(circle2)
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




def circle(r,cp=[0,0],s=50): # circle with radius r and center point cp, s is the number of segments in the circle
    '''
    function for creating points in circle with radius "r", center point "cp" and number of segments "s" 
    '''
    return array([ [cp[0]+r*cos(i*pi/180),cp[1]+r*sin(i*pi/180)] for i in linspace(0,360,s)][0:-1]).tolist()

    
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

def cylinder(r1=1,r2=1,h=1,s=50,r=0,d=0,d1=0,d2=0,center=False):
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
    sec=c2t3(circle(ra,s=s))
    sec1=translate([0,0,h],circle(rb,s=s))
    sol=[sec,sec1]
    if center==True:
        return translate([0,0,-h/2],sol)
    else:
        return sol

def square(s=0,center=False):
    m= s if type(s)==int or type(s)==float else s[0]
    n= s if type(s)==int or type(s)==float else s[1]
    sec=cr(pts1([[0,0,.01],[m,0,.01],[0,n,.01],[-m,0,.01]]),10)
    sec1= [[p[0]-m/2,p[1]-n/2] for p in sec] if center==True else sec
    return sec1

def resize3d(prism,rsz):
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
    return l_(translate(t,rev_prism))

def resize3d_centered(prism,rsz):
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
    return l_(rev_prism)


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



def cube(size=1,center=False):
    '''
    function to draw cube with size 'size'
    refer to the file "example of various functions " for application example
    
    '''
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
   '''
   function to draw sphere with radius 'r' , center point 'cp' and number of segments 's'
   refer to the file "example of various functions " for application example
   
   '''
   path=arc(r,-90,90,s=int(s/2))
   p=[ translate([cp[0],cp[1],p[1]+cp[2]],circle(p[0],s=s)) for p in path]
   return array(p).tolist()


def resize2d(sec,rsz):
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
    
def resize2d_centered(sec,rsz):
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




def intersection_points(sol1,sol2):
    '''
    function to calculate intersection point between two 3d prisms. 
    "sol1" is the 3d object which is intersected with "sol2".
    '''
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


def self_intersections_of_crossing_segments(sec1):
    '''
    calulates the self intersection points of a list of line segments 's'
    it picks the intersection points only if the 2 lines are crossing each other
    e.g.
    sec=segments([[0,0],[10,0],[15,7]])
    self_intersections_crossed_segments(sec) => []
    
    sec=offset_segments([[0,0],[10,0],[15,7]],-1)
    self_intersections_crossed_segments(sec) => 
    [[9.485381933784767, 1.0],
     [4.5075719388876445, 1.0],
     [11.97429003208061, 4.4844710983213805]]
     
    refer to file 'example of various functions' for application example
    '''
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

    
# def self_intersections(sec1):
#     '''self intersections based on Bentley-Ottmann line sweep '''
#     sec2=lexicographic_seg_sort_xy(sec1)
#     s,s1,b=[],[],[]
#     for i in range(len(sec2)):
#         s=s+[sec2[i]]
#         s1.append(s_int1(s))
#         if i>1:
#             for j in range(len(s)):
#                 if s[-1][0][0]>=s[j][1][0]:
#                     b=exclude_seg(s,[s[j]])
#             s=b if b!=[] else s
#             b=[]

#     s1=[p for p in s1 if p!=[]]
#     s1=concatenate(s1).tolist()
#     return s1

def comb(n,i): 
    '''
    calculates number of possible combinations for "n" items with "i" selected items
    comb(8,2) => 28
    '''
    return int(math.factorial(n)/(math.factorial(i)*math.factorial(n-i)))



def bezier(pl,s=20):
    '''
    bezier curve defined by points 'pl' and number of segments 's'
    refer file "example of various functions" for application
    '''
    def B(n,i,t):
        '''
        Bernstein function
        '''
        if i>n or i<0:
            return 0
        else:
            return comb(n,i)*t**i*(1-t)**(n-i)
            
    b=[l_(a_([a_(pl[i])*B(len(pl)-1,i,t) for i in range(len(pl))]).sum(0)) 
       for t in linspace(0,1,s)]
    return b

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
        sec1=[axis_rot(v2,p,a1) for p in s1]
        return sec1



def line_circle_intersection_points(line,cir):
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
    sec= c3t2(sec) if array(sec).shape[-1]==3 else sec
    return (array(cwv(sec))==-1).all()|(array(cwv(sec))==1).all()

def oo_convex(sec,r): #outer offset of a convex section
    s=flip(sec) if cw(sec)==1 else sec
    return offset_points(sec,r)

def circle_point_tangent(cir,p):
    '''
    circle to point tangent line (point should be outside the circle)
    refer file "example of various functions" for application example
    '''
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


def point_circle_tangent(p,cir): # point to circle tangent line (point should be outside the circle)
    '''
    point to circle tangent line (point should be outside the circle)
    refer file "example of various functions" for application example
    '''
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


def variable_sections_extrude(sec,path,o):
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


def two_circle_tangent_arc(c1,c2,r,side=0,s=50): #two circle tangent arc
    '''
    function draws a arc which is tangent to 2 circles 'c1' and 'c2'    's' is the number of segments of the tangent arc
    'r' is the radius of the tangent arc 
    'side' there are 2 sides of the circles where the arc could be created defined by '0' and '1'
    
    '''
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


def arc_with_3_defined_points(p1,p2,p3,s=30):
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

def circle_with_3_defined_points(p1,p2,p3,s=30):
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

    return circle(r,cp,s)

def center_of_circle_with_3_defined_points(p1,p2,p3):
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



def intersection_points_between_2surfaces(surf2,surf1):
    '''
     function to calculate intersection point between two 3d prisms or between surface and solid. 
     "surf2" is the 3d object which is intersected with "surf1".
 try below code for better understanding:
 sec=circle(10);
 path=corner_radius(pts1([[2,0],[-2,0,2],[0,10,3],[-9.9,0]]),5);
 prism=prism(sec,path);
 prism1=q_rot(["y40"],cylinder(r=3,h=15,s=30));
 %swp(prism);
 %swp(prism1);
 ip=intersection_points_between_2surfaces(prism,prism1);
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


def exclude_points_from_list(list1,list_to_exclude):
    '''
    function removes the points from a list of points based on defined list to exclude
    '''
    list1,list_to_exclude=array(list1).round(5),array(list_to_exclude).round(5)
    return list1[~(list_to_exclude==list1[:,None]).all(2).any(1)].tolist()
    
def exclude_segments(list,list_to_exclude):
    '''
    function removes the list of line segments from a segment list
    '''
    return array(list)[~ (array(list)==array(list_to_exclude)[:,None]).all(2).all(2).transpose(1,0).any(1)].tolist()



def points_inside_enclosed_section(section,pnts):
    '''
    function to find points 'pnts' which are inside an enclosed section 'sec'
    refer to the file "example of various functions " for application examples
    '''
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


def rounded_section(line,r1,r2,s=20):
    '''
    creates a rounded section around a line
    radius around first point of line is 'r1' and radius around 2nd point is 'r2'
    
    '''
    cp1=line[0]
    cp2=line[1]
    sec=tctpf(r2,r1,cp2,cp1)
    a1=arc_long_2p(sec[1],sec[2],r1,-1,s=s) if r1>r2 else arc_2p(sec[1],sec[2],r1,-1,s=s)
    a2=arc_long_2p(sec[3],sec[0],r2,-1,s=s) if r2>r1 else arc_2p(sec[3],sec[0],r2,-1,s=s)
    sec1=a1+a2
    return sec1



def cs1(sec,d):
    '''
    creates a cleaning section for removing excess points for offseting a section 'sec' with offset distance 'd'
    
    '''
    r=abs(d)
    a=seg(sec)
    cs=[r_sec(r-r/1000,r-r/1000,p2[0],p2[1]) for p2 in a ]
    return cs


def swp(bead2):
    '''
    function to render various 3d shapes
    example:
    swp(cylinder(d=10,h=20)) will render a cylinder with dia 10 and height 20
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
    

def swp_prism_hollow(prism_big,prism_small):
    '''
    
    creats a hollow prism with 2 similar prisms (1 big and 1 smaller)
    
    refer the file "example of various functions" for application example
    '''
    
    p1=prism_big
    p2=flip(prism_small)
    p3=p1+p2+[p1[0]]
    return p3



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


def helix(radius=10,pitch=10, number_of_coils=1, step_angle=1):
    '''
    creates a helix with radius, pitch and number of coils parameters
    
    refer to file "example of various functions" for application example
    
    '''
    return l_(a_([[radius*cos(d2r(i)),radius*sin(d2r(i)),i/360*pitch] for i in arange(0,360*number_of_coils,step_angle)]))


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
            s3=translate(p0,flip([axis_rot(v2,p,-a1) for p in s2]))
            s4.append(s3)
    return s4

def points_and_faces(sol1):
    '''
    function returns points and faces of a prism
    refer file "example of various functions" for application example
    '''
    bead2=sol1
    n1=arange(len(bead2[0])).tolist()
    n2=array([[[[(j+1)+i*len(bead2[0]),j+i*len(bead2[0]),j+(i+1)*len(bead2[0])],[(j+1)+i*len(bead2[0]),j+(i+1)*len(bead2[0]),(j+1)+(i+1)*len(bead2[0])]] \
             if j<len(bead2[0])-1 else \
             [[0+i*len(bead2[0]),j+i*len(bead2[0]),j+(i+1)*len(bead2[0])],[0+i*len(bead2[0]),j+(i+1)*len(bead2[0]),0+(i+1)*len(bead2[0])]] \
                 for j in range(len(bead2[0]))] for i in range(len(bead2)-1)]).reshape(-1,3).tolist()
    n3=(array(flip(arange(len(bead2[0]))))+(len(bead2)-1)*len(bead2[0])).tolist()
    n=[n1]+n2+[n3]
    pnt=array(bead2).reshape(-1,3).round(4).tolist()
    return [pnt,n]


def fillet_at_instersection_of_2solids(p=[],p1=[],r=1,s=10,o=0,f=1.8):
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

    p2=cpo(p1)
    pb=[[[p2[i][j],p2[i][j+1]] for j in range(len(p2[0])-1)] for i in range(len(p2))]
    pb=array(pb).reshape(-1,2,3)
    
    p01,p02,p03,p04,p05=pa[:,0],pa[:,1],pa[:,2],pb[:,0],pb[:,1]

    v1,v2,v3=p05-p04,p02-p01,p03-p01
    i,j=len(v1),len(v2)

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


    pnt1=(p04[:,None]+einsum('ijk,ij->ijk',array([v1]*j).transpose(1,0,2),t1))[condition]


    uv1=v1/norm(v1,axis=1).reshape(-1,1)
    uv1=array([uv1]*j).transpose(1,0,2)[condition]


    a=cross(v2,v3)
    b=a/(norm(a,axis=1).reshape(-1,1)+.00001)
    b=array([b]*i)[condition]


    nxt_pnt=array(pnt1[1:].tolist()+[pnt1[0]])
    v_rot=nxt_pnt-pnt1

    if o==0:
        cir=array([[pnt1[i]+array(axis_rot(v_rot[i],b[i]*r,t)) for t in linspace(0,180,5)] for i in arange(len(pnt1))]).tolist()
    else:
        cir=array([[pnt1[i]+array(axis_rot(v_rot[i],b[i]*r,-t)) for t in linspace(0,180,5)] for i in arange(len(pnt1))]).tolist()

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


    pnt3=(p04[:,None]+einsum('ijk,ij->ijk',array([v1]*j).transpose(1,0,2),t1))[condition]
    pnt3=sort_pointsv(pnt1,pnt3) if len(pnt1)!=len(pnt3) else pnt3.tolist()


    cir=array([[pnt1[i]+array(axis_rot(v_rot[i],b[i]*r,-t)) for t in linspace(-90,90,5)] for i in arange(len(pnt1))]).tolist()

    
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

    pnt2=(p04[:,None]+einsum('ijk,ij->ijk',array([v1]*j).transpose(1,0,2),t1))[condition]
    pnt2=sort_pointsv(pnt1,pnt2) if len(pnt2)!=len(pnt1) else pnt2.tolist()

    
    sol=array([pnt3,pnt1,pnt2]).transpose(1,0,2)

    sol=[array(bezier([p1,(p1+p2)/2,((p1+p2)/2+(p2+p3)/2)/2,(p2+p3)/2,p3],s)).tolist()[:s+1]+[p2.tolist()] for (p1,p2,p3) in array(sol)]
    sol=sol+[sol[0]]
    return sol



def fillet_at_intersection_of_surface_and_solid(p=[],p1=[],r=1,s=10,o=0,f=1.8):
    '''
    function to calculate fillet at the intersection point of 2 solids
    'p': surface
    'p1': solid
    'r': radius of the fillet
    's': number of segments in the fillet, more number of segments will give finer finish
    'o': option '0' produces fillet in outer side of the intersection and '1' in the inner side of the intersections
    refer file "example of various functions" for application
    '''
    pa=[[[[p[i][j],p[i][j+1],p[i+1][j]],[p[i+1][j+1],p[i+1][j],p[i][j+1]]] if j<len(p[0])-1 else \
         [[p[i][j],p[i][0],p[i+1][j]],[p[i+1][0],p[i+1][j],p[i][0]]] \
         for j in range(len(p[0])-1)] for i in range(len(p)-1)]
    pa=array(pa).reshape(-1,3,3)

    p2=cpo(p1)
    pb=[[[p2[i][j],p2[i][j+1]] for j in range(len(p2[0])-1)] for i in range(len(p2))]
    pb=array(pb).reshape(-1,2,3)
    
    p01,p02,p03,p04,p05=pa[:,0],pa[:,1],pa[:,2],pb[:,0],pb[:,1]

    v1,v2,v3=p05-p04,p02-p01,p03-p01
    i,j=len(v1),len(v2)

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


    pnt1=(p04[:,None]+einsum('ijk,ij->ijk',array([v1]*j).transpose(1,0,2),t1))[condition]

    uv1=v1/norm(v1,axis=1).reshape(-1,1)
    uv1=array([uv1]*j).transpose(1,0,2)[condition]


    a=cross(v2,v3)
    b=a/(norm(a,axis=1).reshape(-1,1)+.00001)
    b=array([b]*i)[condition]

    nxt_pnt=array(pnt1[1:].tolist()+[pnt1[0]])
    v_rot=nxt_pnt-pnt1

    if o==0:
        cir=array([[pnt1[i]+array(axis_rot(v_rot[i],b[i]*r,t)) for t in linspace(0,180,5)] for i in arange(len(pnt1))]).tolist()
    else:
        cir=array([[pnt1[i]+array(axis_rot(v_rot[i],b[i]*r,-t)) for t in linspace(0,180,5)] for i in arange(len(pnt1))]).tolist()

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


    pnt3=(p04[:,None]+einsum('ijk,ij->ijk',array([v1]*j).transpose(1,0,2),t1))[condition]
    pnt3=sort_pointsv(pnt1,pnt3) if len(pnt3)!= len(pnt1) else pnt3.tolist()


    cir=array([[pnt1[i]+array(axis_rot(v_rot[i],b[i]*r,-t)) for t in linspace(-90,90,5)] for i in arange(len(pnt1))]).tolist()

    
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

    pnt2=(p04[:,None]+einsum('ijk,ij->ijk',array([v1]*j).transpose(1,0,2),t1))[condition]
    pnt2=sort_pointsv(pnt1,pnt2) if len(pnt2)!= len(pnt1) else pnt2.tolist()

    
    sol=array([pnt3,pnt1,pnt2]).transpose(1,0,2)

    sol=[array(bezier([p1,(p1+p2)/2,((p1+p2)/2+(p2+p3)/2)/2,(p2+p3)/2,p3],s)).tolist()[:s+1]+[p2.tolist()] for (p1,p2,p3) in array(sol)]
    sol=sol+[sol[0]]
    return sol



def sec_radiuses(sec):
    '''
    function lists the radiuses of a rounded section
    '''
    a=list_r(sec)
    if (a==a[0]).all() or list_r(sec).max()>array(bb2d(sec)).max():
        return zeros(len(sec)).tolist()
    else:
        a=list_r(sec)
        b=[a[i] for i in range(len(a)-1) if a[i+1]==0]
        c=[b[i] for i in range(len(b)-1) if b[i+1]==0]+[b[-1]]
        d= [0.0]+c if len(c)!=len(convert_secv1(sec)) else c
        return l_(d)

def bounding_box2d(sec):
    '''
    bounding boundary for a 2d 
    '''
    return l_([array(sec)[:,0].max()-array(sec)[:,0].min(),array(sec)[:,1].max()-array(sec)[:,1].min()])




    
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


def intersections_of_adjacent_line_segments(segments):
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




def change_orientation_of_solid_from_circular_to_rectangular(sol,s):#circular to rectangulat orientation
    '''
    change the orientation of points of a cylinder from circular to rectangular orientation
    'sol': is a cylindrical type 3d shape
    's': number of segments required between each straight line segments
    refer to the file 'example of various functions' for application examples 
    '''
    # angle=360/len(sol[0])/2
    sol=cpo(sol)
    return [m_points1(sol[i]+flip(sol[len(sol)-1-i]),s) for i in range(int(len(sol)/2))]
    
    

def section_line_intersection_points(sec,line):
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
    
def section_line_intersection_when_line_outside_of_section(sec,line):
    '''
    section-line intersection:
    function to find intersection between an plane in 3d space and a line
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

    p0.shape,v1.shape,t.shape
    v2=array([v1]*len(t)).reshape(-1,3)
    intersection_point=(p0+einsum('ij,i->ij',v2,t[:,0])).tolist()

    return intersection_point


def turtle3D(path):
    '''
    returns the cumulative sum of points
    example:
    path=[[0,0,1],[0,5,10],[10,3,20]]
    turtle3D(path)=> [[0, 0, 1], [0, 5, 11], [10, 8, 31]]
    
    '''
    return array(path).cumsum(0).tolist()
    
def axis_rot_centered_at_origin(axis,solid,angle):
    '''
    rotate a solid around an axis
    '''

    return (c2t3(solid)@arot(axis,angle)).tolist()

        
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
    


def convert_3lines2fillet(pnt3,pnt2,pnt1,s=10,f=1,orientation=0,style=2):
    '''
    Develops a fillet with 3 list of points in 3d space
    s: number of segments in the fillet, increase the segments in case finer model is required
    f: higher number of factor 'f' reduces the concavity, very high number like >10 will be like chamfer
    refer to the file "example of various functions" for application examples
    
    '''
    sol=l_(array([pnt3,pnt1,pnt2]).transpose(1,0,2))
    sol1=[]
    for i in range(len(sol)):
        p0,p1,p2=sol[i]
        p3=mid_point([p0,p2])
        d=l_len([p0,p3])
        d1=l_len([p3,p1])
        p4=movePointOnLine([p3,p1],p3,d/f) if (d1>d or d1>.5) else p1
        if style==0:
            sol1.append(bspline_open([p0,mid_point([p0,p4]),mid_point([p4,p2]),p2],3,s)+[p1])
        elif style==1:
            sol1.append(bezier([p0,mid_point([p0,p4]),mid_point([p4,p2]),p2],s)+[p1])

        elif style==2:
            sol1.append(bezier([p0,p1,p2],s)+[p1])

    return sol1 if orientation==0 else cpo(sol1)[:-1]
    
def minimum_distance_points(sec,min_d=.1):
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
    
def wrap_section_around_defined_path(sec,path):
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


    
def align_points_of_2_sections(sec1,sec2,ang=10):
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
    
    
    
def align_2D_section_to_defined_vector(v1=[1,0,0],sec=[]):
    '''
    function to align a section 'sec' with a vector 'v1'
    refer file "example of various function" for application examples
    '''
    vz=[0,0,-1]
    vz,v1=array([vz,v1])

    nvzv1=cross(vz,v1)
    u1=v1/norm(v1)
    theta=r2d(arccos(u1@vz))
    sec=c2t3(flip(sec))

    sec1=sec@arot([1,0,0],theta)

    theta1=ang(v1[0],v1[1])
    sec1=sec1@zrot(-90)@zrot(theta1)
    return sec1.tolist()



def cut_plane(nv=[0,0,1],size=[5,5],thickness=10,trns1=0,trns2=0,trns3=0,theta=[0,0,0]): #oriented solid
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
    plane2=translate(u1*thickness,plane1)
    sol=[plane1]+[plane2]
    sol=axis_rot(u3,axis_rot(-u2,axis_rot(nv,sol,theta[0]),theta[1]),theta[2])
    sol=translate(u1*trns1,sol)
    sol=translate(u2*trns2,sol)
    sol=translate(u3*trns3,sol)
    return sol
    
    
def slice_sol(sol,n=10):
    '''
    function to slice a solid with 'n' intermediate steps.
    this creats n steps for each turn in sol
    '''
    a=cpo(sol)
    sol1=[[ls(p,n)+[p[1]] for p in seg(a[i])[:-1]] for i in range(len(a))]
    sol2=array(sol1).transpose(1,2,0,3)
    b,c,d,e=sol2.shape
    sol2=sol2.reshape(b*c,d,e).tolist()
    return sol2

def slice_sol_1(sol_1,n=10):
    '''
    function to slice a solid with 'n' intermediate steps
    '''
    return cpo([equidistant_path(p,n) for p in cpo(sol_1)])

    
def center_arc(arc1):
    '''
    function returns the center point of a given circle or arc
    
    '''
    n=int(len(arc1)/360*120)
    p0=arc1[0]
    p1=arc1[n]
    p2=arc1[n*2]
    return cp_3p(p0,p1,p2)
    
def radius_arc(arc1):
    '''
    function returns the radius of a given circle or arc
    
    '''
    n=int(len(arc1)/360*120)
    p0=arc1[0]
    p1=arc1[n]
    p2=arc1[n*2]
    return r_3p([p0,p1,p2])
    
def fillet_between_line_and_circle(line=[],cir1=[],fillet_radius=1,s=20):
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




def oriented_solid(nv=[0,0,1],sec=[],thickness=10,trns1=0,trns2=0,trns3=0, theta=[0,0,0]): #oriented solid
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

    plane2=translate(u1*thickness,plane1)
    sol=[plane1]+[plane2]
    sol=axis_rot(u3,axis_rot(-u2,axis_rot(nv,sol,theta[0]),theta[1]),theta[2])
    sol=translate(u1*trns1,sol)
    sol=translate(u2*trns2,sol)
    sol=translate(u3*trns3,sol)
    return sol
    
def points_list_projection_on_enclosed_3D_section(p0,sec): #point's projection on an enclosed 3d section
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



def points_list_projection_on_defined_plane(p0,v1,loc):#point's projection on a plane
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

    
    
def path_extrude_to_multiple_sections(sec_list,path):
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
    

    
def intersection_points_between_solid_and_line(sol,line):# when line has more than 2 points
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

    return i_p2.tolist()
    


def align_various_sections_of_a_solid(sol,ang=10):
    '''
    function to straighten the twists in the path_extruded sections for better alignments
    refer to the file "example of various functions.ipynb" for application examples
    '''
    sol1=[sol[0]]
    for i in range(1,len(sol)):
        a=align_sec(sol1[i-1],sol[i],ang=ang)
        sol1.append(a[1])
    return sol1

def extrude_solid_to_defined_path(sec,path1,path2):
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

    sec_list=[offset(sec,x) for (x,y) in path1]
    sol=path_extrude2msec(sec_list,path3)
    return sol

    
def points_and_faces_of_closed_ends_solid(bead2):
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

    
def equidistant_path_open(path,s=10,pitch=[]):
    '''
    divides a path in to equally spaced points
    refer file 'example of various functions.ipynb' for application examples
    '''
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

def equidistant_path_closed(path,s=10,pitch=[]):
    '''
    divides a closed path in to equally spaced points
    refer file 'example of various functions.ipynb' for application examples
    '''
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
    return p_rev[:int(s)]

    
def ang_between_2lines_ccw(p0,p1,p2):
    '''
    ccw angle of the line p0p2 from base line p0p1
    '''
    p0,p1,p2=array([p0,p1,p2])
    v1,v2=p1-p0,p2-p0
    a1=ang(v1[0],v1[1])
    a2=ang(v2[0],v2[1])
    return 360 if a1-a2==0 else 360-(a1-a2) if a2<a1 else a2-a1 

def ang_between_2lines_cw(p0,p1,p2):
    '''
    cw angle of the line p0p2 from the base line p0p1
    '''
    p0,p1,p2=array([p0,p1,p2])
    v1,v2=p1-p0,p2-p0
    a1=ang(v1[0],v1[1])
    a2=ang(v2[0],v2[1])
    return 0 if a1-a2==0 else a1-a2 if a2<a1 else 360+(a1-a2) 

def length_of_closed_loop_line_with_multiple_segments(l):
    '''
    calculates sum of lengths of all the segments in a line 'l' considering the section is closed
    '''
    return array([l_len(p) for p in seg(l)]).sum()

def length_of_open_loop_line_with_multiple_segments(l):
    '''
    calculates sum of lengths of all the segments in a line 'l' considering the section is open
    '''
    return array([l_len(p) for p in seg(l)[:-1]]).sum()


def area_of_triangular_section(s):
    '''
    area of the triangle enclosed with in 3 vertices 's'
    '''
    return norm(cross(array(s[1])-array(s[0]),array(s[2])-array(s[0])))/2



def offset_solid(sol,d,o=0):
    '''
    function to calculate offset of a 3d object by distance 'd'
    option 'o' can be set to '0' or '1' depending on the shape of the object.
    in case the shape of the 3d object is twisted,option should be set to '1'
    in few cases this function may not work 
    
    '''
    
    n=array([len(remove_extra_points(p)) for p in sol]).argmax()
    if o==0:
        sol=[sort_points(sol[n],offset_3d(p,d)) for p in sol]
    else:
        sol=[offset_3d(p,d) for p in sol]
    return sol
    


def intersection_points_between_2solids(sol1,sol2,n=0):
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

def intersection_points_between_surface_and_solid(surface,solid,n=0):
    sol1=surface
    sol2=solid
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
    
def lexicographic_segments_sort_xy(sec1):

    a=lexicographic_sort_xy(array(sec1)[:,0])
    sec2=[[[p1,p2] if p1[0]<p2[0] else [p2,p1] for (p1,p2) in sec1 if (p2==p) or (p1==p)] for p in a]
    sec2=remove_extra_points(concatenate(sec2))
    return sec2
    

def circle_from_3points_3d(points,s=20):
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

def center_of_circle_3d(cir):
    '''
    center point of circle with atleast 3 known list of 'points' in 3d space
    '''
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


def tangents_along_a_path(path,scale=1):
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

def normals_along_a_path(path,scale=1):
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

def orthos_along_a_path(path,scale=1):
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


def line_section_intersection_points(line,sec):
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
    return p_l

def line_section_intersection_points_3d(sec,line):
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

def path_offset_new(sec,r):
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
    
def faces_with_first_and_last_end_closed(sol):
    '''
    calculate the faces for the vertices with shape l x m with first and the last end closed
    '''
    l,m=len(sol),len(sol[0])
    n1=arange(m,dtype=int)
    n2=array([[[[(j+1)+i*m,j+i*m,j+(i+1)*m],[(j+1)+i*m,j+(i+1)*m,(j+1)+(i+1)*m]] \
             if j<m-1 else \
             [[0+i*m,j+i*m,j+(i+1)*m],[0+i*m,j+(i+1)*m,0+(i+1)*m]] \
                 for j in range(m)] for i in range(l-1)],dtype=int).reshape(-1,3)
    n3=(array(flip(arange(m)),dtype=int)+(l-1)*m)
    n=[n1.tolist()]+n2.tolist()+[n3.tolist()]
    return n


def faces_with_first_and_the_last_ends_open(sol):
    '''
    returns the faces for the vertices with shape l x m with first and the last end open
    '''
    l,m=len(sol),len(sol[0])
    return faces(l,m)[1:-1]

    

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

    
def align_various_section_points_of_solid(sol):
    '''
    function to straighten the twists in the path_extruded sections for better alignments
    refer to the file "example of various functions.ipynb" for application examples
    '''
    sol1=[sol[0]]
    for i in range(1,len(sol)):
        a=align_sec_1(sol1[i-1],sol[i])
        sol1.append(a[1])
    return sol1
    
def align_points_of_2sections(sec1,sec2):
    '''
    function to align 2 3d sections to obtain the non twisted optimised solid
    refer file "examples of various functions" for application examples
    '''
    
    area1=[ norm(array(sec2[i:]+sec2[:i])-array(sec1),axis=1).sum() for i in range(len(sec2)) ]
    i=array(area1).argmin()
    sol2=[sec1]+[sec2[i:]+sec2[:i]]
    return sol2


def convert_3lines2fillet_closed(pnt3,pnt2,pnt1,s=10,f=1, orientation=0,style=2):
    '''
    Develops a fillet with 3 list of points in 3d space
    s: number of segments in the fillet, increase the segments in case finer model is required
    refer to the file "example of various functions" for application examples
    
    '''
    sol=l_(array([pnt3,pnt1,pnt2]).transpose(1,0,2))
    sol1=[]
    for i in range(len(sol)):
        p0,p1,p2=sol[i]
        p3=mid_point([p0,p2])
        d=l_len([p0,p3])
        d1=l_len([p3,p1])
        p4=movePointOnLine([p3,p1],p3,d/f) if (d1>d or d1>.5) else p1
        if style==0:
            sol1.append(bspline_open([p0,mid_point([p0,p4]),mid_point([p4,p2]),p2],3,s)+[p1])
        elif style==1:
            sol1.append(bezier([p0,mid_point([p0,p4]),mid_point([p4,p2]),p2],s)+[p1])
        elif style==2:
            sol1.append(bezier([p0,p1,p2],s)+[p1])
            
    sol1=sol1+[sol1[0]]
    return sol1 if orientation==0 else cpo(sol1)[:-1]
    
def gcd(a,b):
    '''
    calculates the greatest common divisor of 2 numbers 'a','b'
    '''
    for _ in range(max([a,b])):
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

def axis_rot_centered_at_objects_center(axis,solid,angle):
    '''
    rotate a solid around an axis considering the solid is centered at origin
    '''
    s_1=len(array(solid).shape)
    cp1=array(solid).mean(0) if s_1==2 else array(prism_center(solid))
    solid=translate(-cp1,solid)
    return translate(cp1,solid@arot(axis,angle))

def edges(l,m):
    return array([[[i+j*m,i+(j+1)*m] for j in range(l-1)] for i in range(m)]).reshape(-1,2)


def path_length(path):
    '''
    calculates the length of the path
    '''
    v=[p[1]-p[0] for p in array(seg(path)[:-1])]
    l=[l_len(p) for p in seg(path)[:-1]]
    c=array(l).cumsum()
    
    return c[-1]


def arc_2p_3d(n1,p0,p1,r,cw=1,s=20):
    '''
    draws an arc through 2 points 
    n1: normal vector to define plane on which the arc will be drawn
    r: radius of the arc
    cw: '1' stands for clockwise and '-1'stands for counter-clockwise
    's' is the number of segments of the circle
    '''
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
    
def center_arc_2p_3d(n1,p0,p1,r,cw=1):
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

    
def axis_rot_centered_at_any_defined_location(sol,ax1,loc1,theta):
    '''
    rotate a solid on any pivot point 'loc1' with axis of rotation 'ax1' by an angle 'theta'
    
    '''
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


    
def match_points_of_a_path_with_original_path_open(original_path,path_to_match):
    '''
    function to match the points of original_path with path
    i.e. original_path is independent variable and path_to_match is dependent variable
    '''
    path1=original_path
    path=path_to_match
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

def match_points_of_a_path_with_original_path_closed(original_path,path_to_match):
    '''
    function to match the points of original_path with path.Both the paths are closed loop
    i.e. original_path is independent variable and path_to_match is dependent variable
    '''
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



def find_triangular_meshes_of_intersecting_points(ip,sol1):
    '''
    function to find the triangles on the solid 'sol1' where the intersection points list 'ip' lies
    '''
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
    decision=(t>=-0.01)&(t<=1)&(u>=-0.01)&(u<=1)&(v>=-0.01)&(v<=1)&((u+v)<=1)
    tri_1=array([tri[decision[i]][0] for i in range(len(ip))]).tolist()

    return tri_1


def find_points_with_defined_distance_from_line(line,pnts,d):
    '''
    finds the points in list 'pnts' which are less than distance 'd' from the 'line'
    '''
    if array(line).shape[-1]==3:
        l1=array([line[0],line[1]])
        v1=l1[1]-l1[0]
        v2=array(pnts)-l1[0]
        v2sint=norm(cross(v1,v2),axis=1)/norm(v1)
        v2cost=einsum('j,ij->i',v1,v2)/norm(v1)
        tx=v2cost/norm(v1)
        d1=(tx>=0) & (tx<=1) & (v2sint<d)
        p7=array(pnts)[d1].tolist()
    elif array(line).shape[-1]==2:
        l1=array([line[0],line[1]])
        v1=l1[1]-l1[0]
        v2=array(pnts)-l1[0]
        v2sint=cross(c23(v1),c23(v2))[:,-1]/norm(v1)
        v2cost=einsum('j,ij->i',v1,v2)/norm(v1)
        tx=v2cost/norm(v1)
        d1=(tx>=0) & (tx<=1) & (v2sint<d)
        # p7=array(pnts)[d1].tolist()
        p7=l_(a_(pnts)[d1])
    return p7
    
def find_distance_of_points_from_line(line,pnts):
    '''
    finds the perpendicular distance of the points in list 'pnts' from the 'line'
    only calculates if the projecton of the point lies within the line
    '''
    if array(line).shape[-1]==3:
        l1=array([line[0],line[1]])
        v1=l1[1]-l1[0]
        v2=array(pnts)-l1[0]
        v2sint=norm(cross(v1,v2),axis=1)/norm(v1)
        v2cost=einsum('j,ij->i',v1,v2)/norm(v1)
        tx=v2cost/norm(v1)
        d1=(tx>=0) & (tx<=1)
        v2sint=array(v2sint)[d1].tolist()
    elif array(line).shape[-1]==2:
        l1=array([line[0],line[1]])
        v1=l1[1]-l1[0]
        v2=array(pnts)-l1[0]
        v2sint=cross(c23(v1),c23(v2))[:,-1]/norm(v1)
        v2cost=einsum('j,ij->i',v1,v2)/norm(v1)
        tx=v2cost/norm(v1)
        d1=(tx>=0) & (tx<=1)
        v2sint=array(v2sint)[d1].tolist()
    return v2sint
    

    
def plane_to_plane_intersection_line(pa,pb):#plane to plane intersection line
    '''
    function to calculate intersection line between 2 planes
    '''
    x,y,z=sympy.symbols('x y z')
    p0,p1,p2=array(pa)
    v1,v2=p1-p0,p2-p0
    n1=cross(v1,v2)

    p3,p4,p5=array(pb)
    v1,v2=p4-p3,p5-p3
    n2=cross(v1,v2)

    eq1=n1[0]*x+n1[1]*y+n1[2]*z-p0@n1
    eq2=n2[0]*x+n2[1]*y+n2[2]*z-p3@n2

    f=sympy.solve([eq1,eq2],x,y)

    v4=cross(n1,n2)
    p6=array([f[x].subs(z,0),f[y].subs(z,0),0])
    p7=array(p6)+v4*10
    p8=array(p6)-v4*10
    line=array([p7,p8]).tolist()
    return array(line).astype(float).tolist()

def offset_an_intersection_line_on_solid( intersection_line,solid, offset_distance, options=0, f=1, closed=0):
    '''
    function to offset the intersection points 'i_p' on a solid 'sol' by distance 'r'. option 'o' can have values '0' or '1' and changes the direction of offset.
    for closed loop path set closed=1
    '''
    i_p=intersection_line
    sol=solid
    r=offset_distance
    o=options
    
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

def  offset_an_intersection_line_on_solid_revised(intersection_line,solid, offset_distance, options=0,closed=0,type=0,dist=0,vx=[],edges_closed=1,cg=0):
    i_p=intersection_line
    sol=solid
    r=offset_distance
    o=options
    '''
    function to offset the intersection points 'i_p' on a solid 'sol' by distance 'r'. option 'o' can have values '0' or '1' and changes the direction of offset.
    for closed loop path set closed=1
    o: set '1' to shift the line on other side
    type: set '1' if prism_center lies outside the solid. type can also be set to '2'.
    for using type=2, refer function definition psos_v_1 and also checkout the video explanation
    prism_center is the center of the bounding box of the solid
    '''
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

def offset_an_intersection_line_on_surface(intersection_line, surface, offset_distance, options=0, f=1, closed=0):
    '''
    function to offset the intersection line on a surface by a defined distance . options can have values '0' or '1' and changes the direction of offset
    for closed loop path, set closed=1
    '''
    i_p=intersection_line
    sol=surface
    r=offset_distance
    o=options
    
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


def fillet_at_intersection_line_between_two_solids (sol1, sol2, offset_distance_on_sol1, offset_distance_on_sol2, number_of_segments_in_fillet=20, o=0,type=1,dist=0,vx=[],style=2,f=1,edges_closed=1,c=0):
    '''
    calculates a fillet at the intersection of 2 solids.
    r1 and r2 would be same in most of the cases, but the signs can be different depending on which side the fillet is required
    r1 is the distance by which intersection line offsets on sol1 and similarly r2 is on sol2 
    for type, dist and vx parameters refer function o_3d_rev
    '''
    r1=offset_distance_on_sol1
    r2=offset_distance_on_sol2
    s=number_of_segments_in_fillet
    
    sol1=cpo(cpo(sol1)+[cpo(sol1)[0]])
    p1=ip_sol2sol(sol1,sol2,o)
    p2=i_p_p(sol2,p1,r2)
    p3=o_3d_rev(p1,sol1,r1,type=type,dist=dist,vx=vx,edges_closed=edges_closed,cg=c)
    fillet1=convert_3lines2fillet(p3,p2,p1,s=s,style=style,f=f)
    
    return fillet1

def fillet_at_intersection_line_between_two_solids_closed_loop( sol1, sol2, offset_distance_on_sol1, offset_distance_on_sol2, number_of_segments_in_fillet=20, o=0,type=1,dist=0,vx=[],style=2,f=1,edges_closed=1,c=0):
    '''
    calculates a fillet at the intersection of 2 solids.
    r1 and r2 would be same in most of the cases, but the signs can be different depending on which side the fillet is required
    r1 is the distance by which intersection line offsets on sol1 and similarly r2 is on sol2 
    for type, dist and vx parameters refer function o_3d_rev
    '''
    r1=offset_distance_on_sol1
    r2=offset_distance_on_sol2
    s=number_of_segments_in_fillet
    
    sol1=cpo(cpo(sol1)+[cpo(sol1)[0]])
    p1=ip_sol2sol(sol1,sol2,o)
    p2=i_p_p(sol2,p1,r2)
    p3=o_3d_rev(p1,sol1,r1,type=type,dist=dist,vx=vx,edges_closed=edges_closed,cg=c)
    fillet1=convert_3lines2fillet_closed(p3,p2,p1,s=s,style=style,f=f)
    
    return fillet1


def fillet_at_intersection_line_between_surface_and_solid( surf, sol, offset_distance_on_surf, offset_distance_on_sol, number_of_segments_in_fillet=20 , type=1,dist=0,vx=[],style=2,f=1,edges_closed=0,c=0):
    '''
    calculates a fillet at the intersection of surface with solid.
    r1 and r2 would be same in most of the cases, but the signs can be different depending on which side the fillet is required
    r1 is the distance by which intersection line offsets on surf and similarly r2 is on sol 
    '''
    r1=offset_distance_on_surf
    r2=offset_distance_on_sol
    s=number_of_segments_in_fillet
    
    p1=ip_surf2sol(surf,sol)
    p2=i_p_p(sol,p1,r2)
    p3=o_3d_rev(p1,surf,r1,type=type,dist=dist,vx=vx,edges_closed=edges_closed,cg=c)
    fillet1=convert_3lines2fillet(p3,p2,p1,s=s,style=style,f=f)

    return fillet1

def fillet_at_intersection_line_between_surface_and_solid_closed_loop ( surf, sol, offset_distance_on_surf, offset_distance_on_sol, number_of_segments_in_fillet=20, type=1, dist=0, vx=[], style=2, f=1, edges_closed=0, c=0):
    '''
    calculates a fillet at the intersection of surface with solid.
    r1 and r2 would be same in most of the cases, but the signs can be different depending on which side the fillet is required
    r1 is the distance by which intersection line offsets on surf and similarly r2 is on sol 
    '''
        
    p1=ip_surf2sol(surf,sol)
    p2=i_p_p(sol,p1,r2)
    p3=o_3d_rev(p1,surf,r1,type=type,dist=dist,vx=vx,edges_closed=edges_closed,cg=c)
    fillet1=convert_3lines2fillet_closed(p3,p2,p1,s=s,style=style,f=f)

    return fillet1

    
def ellipse(semi_major,semi_minor,number_of_segments=50):
    '''
    returns the ellipse with defined semi-major and semi-minor axis
    '''
    a,b,s=semi_major,semi_minor,number_of_segments
    return l_([[a*cos(d2r(i)),b*sin(d2r(i))]  for i in linspace(0,360,s)[:-1]])
    
def arc_tangent_to_3_points_triangle(p0,p1,p2,r,s=20):
    '''
    draws an arc with 3 points e.g. p0,p1,p2
    where p0p1 is line1 and p1p2 is line2
    so it draws an arc which is tangent to p0p1 and p1p2 and not inclusive of 3 points
    
    '''
    if r>0:
        p0,p1,p2=array([p0,p1,p2])
        if cw([p0,p1,p2])==1:
            a1=180-ang_2lineccw(p1,p0,p2)
        else:
            a1=180-ang_2linecw(p1,p0,p2)

        d=r*tan(d2r(a1/2))
        u10=(p0-p1)/norm(p0-p1)
        u12=(p2-p1)/norm(p2-p1)
        pa=p1+u10*d
        pb=p1+u12*d

        arc_1=arc_2p(pa,pb,r,cw([p0,p1,p2]),s)
        return arc_1
    else:
        return [p1]


def extrude_solid_to_path(sol,path):
    '''
    function, extrudes solid with slices to any arbitrary path
    '''
    sol1=c3t2(sol)
    zpath=[[0,0,p[0][2]] for p in sol]
    path2=path2path1(zpath,path)
    sol2=align_sol_1(path_extrude2msec(sol1,path2))
    return sol2

def comb_list(n):
    '''
    returns the list of combination for number of integers or characters where order doesn't matter.
    example:
    comb_list(3) =>
    array([[0, 1],
           [0, 2],
           [1, 2]])
    '''
    n=arange(n)
    a=array([n]*len(n)).transpose(1,0)
    b=array([n]*len(n))
    c=array([a,b]).transpose(1,0,2).transpose(0,2,1)
    d=concatenate([c[i][i+1:] for i in range(len(c))])
    return d


def corner_radius(sec,s=20):
    '''
    function to create section with corner radiuses. e.g. 
    following code has 3 points at [0,0],[10,0] and [7,15] and radiuses of 0.5,2 and 1 respectively,
    s=5 represent the number of segments at each corner radius.
    sec=corner_radius(pl=[[0,0,.5],[10,0,2],[7,15,1]],s=5)
    
    refer file "example of various functions" for application
    '''
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
    d=c3t2(q_rot(['z.0001'],d))
    return d

def corner_radius_with_turtle(sec,s=20):
    '''
    function to create section with corner radiuses. e.g. 
    following code has 3 points at [0,0],[10,0] and [7,15] and radiuses of 0.5,2 and 1 respectively,
    s=5 represent the number of segments at each corner radius.
    sec=corner_radius(pl=[[0,0,.5],[10,0,2],[7,15,1]],s=5)
    
    refer file "example of various functions" for application
    '''
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
    d=c3t2(q_rot(['z.0001'],d))
    return d




def surround(path,r,s=20):
    '''
    function to surround a path to create a rounded section
    
    '''
    a=path_offset_n(path,r)
    b=path_offset_n(path,-r)
    b=flip(b)
    arc1=arc_2p(a[-1],b[0],r,-1,s)
    arc2=arc_2p(b[-1],a[0],r,-1,s)
    sec=a[1:-1]+arc1+b[1:-1]+arc2
    return sec



def tangent_line_to_a_circle(circle,theta=90,length=1):
    '''
    function to draw a line tangent to a circle defined by radius 'r'
    center point 'cp' and line is defined by angle 'theta' and length 'l'
    angle theta is from x-axis
    length can be positive or negative
    In case negative the line is drawn on the opposite side
    
    '''
    r=r_arc(circle)
    cp=cp_arc(circle)
    l=length
    p0=[r*cos(d2r(270+theta)),r*sin(d2r(270+theta))]
    p1=(array(p0)+q_rot2d(theta,[l,0])).tolist()
    p0,p1= (array([p0,p1])+cp).tolist()
    return [p0,p1]

def fillet_between_intersection_of_two_lines(l1,l2,r,s=10):
    '''
    function calculates the fillet at intersection between 2 lines
    'l1' and 'l2'
    r: radius of fillet
    s: segments of fillet
    '''
    p0=i_p2d(l1,l2)
    l2=l2 if l_(a_(p0).round(4))!=l_(a_(l2[0]).round(4)) else flip(l2)

    
    clock=cw([l1[0],p0,l2[0]])
    a1=ang_2lineccw(p0,l1[0],l2[0]) if clock==1 else \
    ang_2linecw(p0,l1[0],l2[0])
    a2=180-a1
    l_1=r*tan(d2r(a2/2))
    v1=array(l1[0])-array(p0)
    v2=array(l2[0])-array(p0)
    u1,u2=v1/norm(v1),v2/norm(v2)
    p1=array(p0)+u1*l_1
    p2=array(p0)+u2*l_1
    arc_1=arc_2p(p1,p2,r,clock,s)
    return arc_1
    
def cir_line_tangent(c1,l1,side=0):
    '''
    function to draw a line tangent to a circle
    c1 is circle
    l1 is line
    side=0, draws tangent at one side and 1 draws tangent on the other side
    
    '''
    r=r_arc(c1)
    cp=cp_arc(c1)
    v1=array(l1[1])-array(l1[0])
    theta=ang(v1[0],v1[1]) if side==0 else ang(v1[0],v1[1])+180
    l=l_len(l1)
    p0=[r*cos(d2r(270+theta)),r*sin(d2r(270+theta))]
    p1=(array(p0)+q_rot2d((theta if side==0 else theta-180),[l,0])).tolist()
    p0,p1= (array([p0,p1])+cp).tolist()
    return [p0,p1]

def spiral_poly(r=1,d=.3,n=4,t=100):
    '''
    create a spiral polygon
    r: initial length of the line
    d: increment length every iteration
    n: number of sides of the polygon
    t: number of turns
    '''
    theta=360/n
    sec=[[r,0]]
    for i in range(1,t):
        r=r+d
        a=array([r,0])
        sec.append(array(sec[-1])+q_rot2d(i*theta,a))

    sec=array(sec).tolist()
    return sec

def equate_points(sec,sec1):
    '''
    function to make the points in 2 sections equal without changing the location of points
    '''
    c=array(sec1).shape[-1]
    b=array(sec).shape[-1]
    a=lcm(len(sec),len(sec1))
    sec=array([sec]*int(a/len(sec))).transpose(1,0,2).reshape(-1,b).tolist()
    sec1=array([sec1]*int(a/len(sec1))).transpose(1,0,2).reshape(-1,c).tolist()
    return [sec,sec1]


def sinewave(l,n,a,p):
    '''
    creates a sinewave with length 'l', number of cycles 'n'
    amplitude 'a' and number of points 'p'
    '''

    w1=[[i,a*sin(d2r(n*i*360/l))]  for i in linspace(0,l,p)]
    return l_(w1)

def cosinewave(l,n,a,p):
    '''
    creates a cosinewave with length 'l', number of cycles 'n'
    amplitude 'a' and number of points 'p'
    '''

    w1=[[i,a*cos(d2r(n*i*360/l))]  for i in linspace(0,l,p)]
    return l_(w1)
    
def mod(a,b):
    '''
    function to calculate remainder of a divison of numbers
    example:
    mod(6,2) => 0
    '''
    return round(a-sign(a)/sign(b)*b*floor(round(abs(a/b),10)),10);

def e_wave(l=50,a=1,w=0.1,t=100):
    '''
    create a graph of exponential function a*e^-(wt) where
    w: omega
    t: time steps
    a: amplitude
    l: length of time
    
    '''
    return l_([[i,a*exp(-i*w)]  for i in linspace(0,l,t)])

def waves_2d_multiply(w1,w2,a=1):
    '''
    multiply 2 waves 'w1' and 'w2' and amplitude of the resultant wave is set to 'a'
    '''
    w3=[ [w1[i][0],(w1[i][1]*w2[i][1]) ] for i in range(len(w1))]
    a_max=array(w3)[:,1].max()
    w3=[ [p[0],p[1]/a_max*a ] for p in w3]
    return w3

def waves_2d_add(w1,w2,a=1):
    '''
    add 2 waves 'w1' and 'w2' and amplitude of the resultant wave is set to 'a'
    '''
    w3=[ [w1[i][0],(w1[i][1]+w2[i][1]) ] for i in range(len(w1))]
    a_max=array(w3)[:,1].max()
    w3=[ [p[0],p[1]/a_max*a ] for p in w3]
    return w3

def waves_2d_max(w1,w2,a=1):
    '''
    max of 2 waves 'w1' and 'w2' and amplitude of the resultant wave is set to 'a'
    '''
    w3=[ [w1[i][0],max(w1[i][1],w2[i][1]) ] for i in range(len(w1))]
    a_max=array(w3)[:,1].max()
    w3=[ [p[0],p[1]/a_max*a ] for p in w3]
    return w3

def waves_2d_min(w1,w2,a=1):
    '''
    min of 2 waves 'w1' and 'w2' and amplitude of the resultant wave is set to 'a'
    '''
    w3=[ [w1[i][0],min(w1[i][1],w2[i][1]) ] for i in range(len(w1))]
    a_max=array(w3)[:,1].max()
    w3=[ [p[0],p[1]/a_max*a ] for p in w3]
    return w3

def waves_2d_norm(w1,w2,a=1):
    '''
    norm 2 waves 'w1' and 'w2' and amplitude of the resultant wave is set to 'a'
    '''
    w3=[ [w1[i][0],norm([w1[i][1],w2[i][1]]) ] for i in range(len(w1))]
    a_max=array(w3)[:,1].max()
    w3=[ [p[0],p[1]/a_max*a ] for p in w3]
    return w3

def fit_curve2d(curve,pnt):
    '''
    function to transform a 2d shape to a curve
    an example can make this clear
    '''
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
    '''
    function to extrude a wave 'w1' to any defined path 'c1'
    '''
    w2=[[0,p[1],p[2]] for p in w1]
    c2=array(c1)
    c3=array(c1[1:]+c1[:1])
    v1=c3-c2
    w3=[]
    for i in range(len(w2)):
        a1=ang(v1[i][0],v1[i][1])
        a2=ang(norm(v1[i][:2]),v1[i][2])
        p3=q_rot([f'y{-a2}' ,f'z{a1}'],w2[i])
        p3=c2[i]+p3
        w3.append(p3.tolist())
    w3=w3[:-1]
    return w3

def convert_3lines2surface(l1,l2,l3,s=50):
    '''
    convert 3 lines to surface
    's' is the number of segments on each surface line
    '''
    l1=path2path1(l2,l1)
    l3=path2path1(l2,l3)
    l1,l2,l3=array([l1,l2,l3])
    surf_1=[arc_3p_3d([l1[i]+[0,0,.00001],l2[i],l3[i]],s)  for i in range(len(l1))]
    return surf_1

def convert_lines2surface_spline(lines,s=50):
    '''
    create surface with lines, method used is bspline
    '''
    lines=[lines[0]]+[path2path1(lines[0],lines[i])  for i in range(1,len(lines))]
    surf_1=[bspline(p,2,s) for p in cpo(lines)]
    return surf_1



def SurfaceFrom3LinesInDifferentPlanes(a3,a4,a5,s=50):
    '''
    create surface with 3 lines in different plane.
    option 'o' needs to be adjusted between '0' and '1' based on the result correctness
    example:
    w1=arc_3p_3d([[0,0,0],[20,0,5],[40,0,0]])
    w2=arc_3p_3d([[0,0,0],[-1,15,3],[-2,30,0]])
    w3=arc_3p_3d([[-2,30,0],[15,35,4],[30,40,0]])
    surf_2=SurfaceFrom3LinesInDifferentPlanes(w1,w2,w3,s=30)
    '''
    a3=equidistant_path(a3,s)
    a4=equidistant_path(a4,s)
    a5=equidistant_path(a5,s)
    
    surf_1=slice_sol([a3,a5],s)
    surf_1=[translate(array(a4[i]),translate(-array(surf_1[i][0]),surf_1[i])) for i in range(len(surf_1))]
    return surf_1
    

def mid_point(w1):
    return equidistant_path(w1,2)[1]

def cw_3p_3d(points):
    '''
    finds the orientation of 3 points in 3d.
    1 means clockwise
    -1 means counter clockwise
    '''
    n1=array(nv(points))+[.000001,.000001,0]
    a1=cross(n1,[0,0,-1])
    t1=r2d(arccos(n1@[0,0,-1]))
    sec1=translate(-array(points).mean(0),points)
    sec2=c3t2(axis_rot(a1,sec1,t1))
    l1=len(sec2)
    p0,p1,p2=[sec2[0],sec2[int(l1/3)],sec2[int(l1*2/3)]]
    
    return cw([p0,p1,p2])

def smoothen_3d(p0,s=50):
    '''
    draw smooth curves with random points 'p0'
    '''
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
    '''
    draw smooth curves with random points 'p0'
    '''
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


def faces_surface(l,m):
    '''
    calculate the faces for an open surface
    '''
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
    '''
    thicken the surface by amount 'd'
    '''
    surf_2=surface_offset(surf_1,d)
    sol=[surf_1[i]+flip(surf_2[i])  for i in range(len(surf_1))]
    return sol

def boundary_edges_sol(sol):
    '''
    finds the boundary edges of the solid
    '''
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
    '''
    finds the boundary edges of a surface
    '''
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
    '''
    extend a 2d arc by theta degrees
    '''
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
    '''
    extend a 3d arc by theta degrees
    '''
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
    '''
    finds the boundary edges of a surface
    '''
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
    ''' function to create stl_mesh from vertices and faces'''
    v_1,f_1=array(v),array(f)
    sol = mesh.Mesh(zeros(f_1.shape[0], dtype=mesh.Mesh.dtype))
    sol.vectors=v_1[f_1]
    return sol

def arc_with_start_point_and_center(start_point=[],center_point=[],theta=90,segments=30):
    '''
    function to draw an arc with known center_point and start_point
    '''
    center_point,start_point=array([center_point,start_point])
    v1=start_point-center_point
    arc_1=center_point+[q_rot2d(i,v1) for i in linspace(0,theta,segments)]
    arc_1=arc_1.tolist()
    return arc_1

def fillet_between_line_and_circle(l1,c1,r2,cw=-1,option=0,s=50):
    '''
    function to draw a fillet between a line and a circle
    option can be '0' or '1' to flip the fillet from one side to another
    's' is the number of segments in the arc
    '''
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

def fillet_at_intersection_of_two_lines_3d(l1,l2,r,s=10):
    '''
    function calculates the fillet at intersection between 2 3d lines in 1 plane
    'l1' and 'l2'
    r: radius of fillet
    s: segments of fillet
    '''
    l1=array(l1)
    l2=array(l2)
    
    n1=nv(remove_extra_points(array([l1[0],l1[1],l2[0],l2[1]]).round(5))[:3])
    n1=n1/norm(n1)
    n2=cross(n1,[0,0,-1])
    theta=r2d(arccos(n1@[0,0,-1]))
    l1,l2=c3t2(axis_rot(n2,[l1,l2],theta))
    p0=i_p2d(l1,l2)
    l2=l2 if l_(a_(p0).round(4))!=l_(a_(l2[0]).round(4)) else flip(l2)

    clock=cw([l1[0],p0,l2[0]])
    a1=ang_2lineccw(p0,l1[0],l2[0]) if clock==1 else \
    ang_2linecw(p0,l1[0],l2[0])
    a2=180-a1
    l_1=r*tan(d2r(a2/2))
    v1=array(l1[0])-array(p0)
    v2=array(l2[0])-array(p0)
    u1,u2=v1/norm(v1),v2/norm(v2)
    p1=array(p0)+u1*l_1
    p2=array(p0)+u2*l_1
    arc_1=arc_2p(p1,p2,r,clock,s)
    return axis_rot(n2,arc_1,-theta)

def fillet_inside_between_line_and_circle_3d(l1,c1,r2,cw=-1,option=0,s=50):
    '''
    function to draw a fillet between a line and a circle in 3d space, where the fillet is drawn inside of the circle
    circle and line should be in same plane
    option can be '0' or '1' to flip the fillet from one side to another
    's' is the number of segments in the arc
    '''
    tr=array(l1+c1).mean(0)
    l1=translate(-tr,l1)
    c1=translate(-tr,c1)
    
    n1=nv(c1)
    n1=n1/norm(n1)
    n2=cross(n1,[0,0,-1])
    theta=r2d(arccos(n1@[0,0,-1]))
    l1=c3t2(axis_rot(n2,l1,theta))
    c1=c3t2(axis_rot(n2,c1,theta))
    
    v1=array(l1[1])-array(l1[0])
    u1=v1/norm(v1)
    cp1=cp_arc(c1)
    v2=array(cp1)-array(l1[0])
    u2=v2/norm(v2)
    l_1=norm(cross(c23(v1),c23(v2)))/norm(v1)
    p3=array(l1[0])+u1*(u1@v2)
    r1=r_arc(c1)
    d=sqrt((r1-r2)**2-(l_1+r2)**2) if option==0 else sqrt((r1-r2)**2-(l_1-r2)**2)
    p2=p3-u1*d
    v3=array(cp1)-p3
    u3=v3/norm(v3)
    cp2=p2-u3*r2 if option==0 else p2+u3*r2
    v4=array(cp1)-cp2
    u4=v4/norm(v4)
    p4=cp2-u4*r2
    p2,p4=translate(tr,axis_rot(n2,[p2,p4],-theta))
    return arc_2p_3d(n1,p2,p4,r2,cw=cw,s=s)


def fillet_inside_between_line_and_circle_2d(l1,c1,r2,cw=-1,option=0,s=50):
    '''
    function to draw a fillet between a line and a circle, where the fillet is drawn inside of the circle
    option can be '0' or '1' to flip the fillet from one side to another
    's' is the number of segments in the arc
    '''
    v1=array(l1[1])-array(l1[0])
    u1=v1/norm(v1)
    cp1=cp_arc(c1)
    v2=array(cp1)-array(l1[0])
    u2=v2/norm(v2)
    l_1=norm(cross(c23(v1),c23(v2)))/norm(v1)
    p3=array(l1[0])+u1*(u1@v2)
    r1=r_arc(c1)
    d=sqrt((r1-r2)**2-(l_1+r2)**2) if option==0 else sqrt((r1-r2)**2-(l_1-r2)**2)
    p2=p3-u1*d
    v3=array(cp1)-p3
    u3=v3/norm(v3)
    cp2=p2-u3*r2 if option==0 else p2+u3*r2
    v4=array(cp1)-cp2
    u4=v4/norm(v4)
    p4=cp2-u4*r2
    return arc_2p(p2,p4,r2,cw=cw,s=s)

def fillet_outside_between_line_and_circle_3d(l1,c1,r2,cw=-1,option=0,s=50):
    '''
    function to draw a fillet between a line and a circle in 3d space.
    line and circle should be in the same plane
    option can be '0' or '1' to flip the fillet from one side to another
    's' is the number of segments in the arc
    '''
    tr=array(l1+c1).mean(0)
    l1=translate(-tr,l1)
    c1=translate(-tr,c1)
    n1=nv(c1)
    n1=n1/norm(n1)
    n2=cross(n1,[0,0,-1])
    theta=r2d(arccos(n1@[0,0,-1]))
    l1=c3t2(axis_rot(n2,l1,theta))
    c1=c3t2(axis_rot(n2,c1,theta))
    
    v1=array(l1[1])-array(l1[0])
    u1=v1/norm(v1)
    cp1=cp_arc(c1)
    v2=array(cp1)-array(l1[0])
    u2=v2/norm(v2)
    l_1=norm(cross(v1,v2))/norm(v1)
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
    return translate(tr,axis_rot(n2,arc_2p(p2,p4,r2,cw=cw,s=s),-theta))

def mirror_line(p0,n1,loc):
    '''
    function to mirror the points list 'p0' defined by mirroring plane 'n1' with location 'loc'
    '''
    a=ppplane(p0,n1,loc)
    b=[]
    for i in range(len(a)):
        v1=array(a[i])-p0[i]
        b.append((array(a[i])+v1).tolist())
    return b

def change_orientation(surf):
    '''
    change orientation of a surface to make it suitable for creating solid
    '''
    b=[surf[0]+[surf[j][-1] for j in range(1,len(surf)-1)]+flip(surf[-1])+ \
                 flip([surf[j][0] for j in range(1,len(surf)-1)])]
    for i in range(1,int(len(surf)/2)):
        
        b.append(surf[i][i:-i]+[surf[j][-(1+i)] for j in range(i+1,len(surf)-(i+1))]+flip(surf[-(i+1)][i:-i])+ \
                 flip([surf[j][i] for j in range(i+1,len(surf)-(i+1))]))
    if mod(len(surf),2)==0:
        return sort_surface(b)
    else:
        a=int(len(surf)/2)
        return sort_surface(b+[surf[a][a:-a]])
    
def sort_surface(l3):
    '''
    sort the points of a surface to make them equal
    '''
    b=[l3[0]]
    for i in range(1,len(l3)):
        b.append(sort_points(b[-1],l3[i]))
    return b


def solid_from_fillet(fillet_1,d):
    '''
    creates a cutting edge from fillet
    '''
    fillet_1=cpo(fillet_1)[:-1]
    fillet_2=surface_offset(fillet_1,d)
    sol=fillet_1+flip(fillet_2)
    return cpo(sol)

def solid_from_fillet_closed(fillet_1,d):
    '''
    creates a cutting edge from fillet with closed section lines
    '''
    fillet_1=cpo(fillet_1)[:-1]
    fillet_2=surface_offset(fillet_1,d)
    sol=fillet_1+flip(fillet_2)
    sol=cpo(sol)
    sol=sol+[sol[0]]
    return sol

def points_projection_on_surface(p_0,surf):
    '''
    project a point on to a surface
    '''
    a=array(faces_surface(len(surf),len(surf[0])))
    b=array(surf).reshape(-1,3)
    c=b[a]

    p0,p1,p2=c[:,0],c[:,1],c[:,2]
    v1,v2=p1-p0,p2-p0
    u1,u2=v1/norm(v1,axis=1).reshape(-1,1),v2/norm(v2,axis=1).reshape(-1,1)
    n1=cross(v1,v2)
    un1=n1/norm(n1,axis=1).reshape(-1,1)
    p_0=array(p_0)
    # p_0+un1*t1=p0+v1*t2+v2*t3
    iim=array([un1,-v1,-v2]).transpose(1,0,2).transpose(0,2,1)
    im=inv(iim)
    p=p0-p_0
    p.shape,im.shape
    t=einsum('ijk,ik->ij',im,p)
    t_1,t_2,t_3=t[:,0],t[:,1],t[:,2]
    dec=(t_2>=0) & (t_2<=1) & (t_3>=0) & (t_3<=1) & ((t_2+t_3)<=1)
    px=(p_0+einsum('ij,i->ij',un1,t_1))[dec].tolist()
    return px

def surface_normal(s1,length=1):
    '''
    calculates the normal of a surface (average)
    '''
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

def projecting_a_surface_on_to_another(s_1,s_2,n_1=[]):
    '''
    projecting surface s_2 on surface s_1
    surfaces should ideally cross each other completely
    '''
    n_1=surface_normal(s_2,1) if n_1==[] else n_1
    s_3=[project_line_on_surface(p,s_1,n_1) for p in s_2]
    return s_3
  
    
def project_line_on_surface(l_2,surf_1,n_1=[]):
    '''
    function for projecting a line on to a surface.
    n_1 is a direction vector for projecting the line
    an example video can be refered for clarity
    '''
    n_1=surface_normal(surf_1,1) if n_1==[] else n_1
    # p1+v1*t1=p2+v2*t2+v3*t3
    f_1=faces_surface(len(surf_1),len(surf_1[0]))
    v_1=array(surf_1).reshape(-1,3)
    tri=v_1[f_1]
    p2,p3,p4=tri[:,0],tri[:,1],tri[:,2]
    v2,v3=p3-p2,p4-p2
    v1=array(n_1)
    p1=array(l_2)
    v1=array([[v1]*len(p2)]*len(p1))
    v2=array([v2]*len(p1))
    v3=array([v3]*len(p1))
    iim=array([v1,-v2,-v3+.000001]).transpose(1,2,0,3).transpose(0,1,3,2)
    im=inv(iim)
    p=p2[None,:,:]-p1[:,None,:]
    t=einsum('ijkl,ijl->ijk',im,p)
    t2,t3=t[:,:,1],t[:,:,2]
    dec=(t2>=0)&(t2<=1)&(t3>=0)&(t3<=1)&((t2+t3)<=1)
    ip_1=p1[:,None,:]+einsum('ijk,ij->ijk',v1,t[:,:,0])
    ip_1=ip_1[dec].tolist()
    return ip_1

def path_offset_3d(sec,d):
    '''
    path_offsets an enclosed section in 3d space, in case the section is in 1 plane
    sec: section in 3d space
    d: offset distance -ve sign means inner offset and +ve sign is outer offset
    refer to the file"example of various functions" for application examples
    
    '''
    sec0=remove_extra_points(sec)
    sec0=q_rot(['z.00001'],sec0)
    avg1=array(sec0).mean(0)
    sec1=translate(-avg1,sec0)

    nv1=-array(nv(sec1))
    nz=[0,0,1]
    nr=cross(nv1,nz) if abs(nv1).tolist()!=[0,0,1] else nv1
    theta=r2d(arccos(nv1@array(nz)))
    sec1=axis_rot(nr,sec1,theta)
    z_values=array(sec1)[:,2]-avg1[2]
    sec1=ppplane(sec1,[0,0,1],[0,0,0])
    sec1=c3t2(sec1)
    x_values=array([l_len([[0,0],p])  for p in sec1])
    sec2=path_offset_n(sec1,d)
    x1_values=array([l_len([[0,0],p])  for p in sec2])
    z1_values=z_values/x_values*x1_values
    z1_values=array([[0,0,p] for p in z1_values])
    sec2=array(c2t3(sec2))
    sec2=axis_rot(nr,sec2,-theta)
    sec2=translate(array(sec).mean(0),sec2)
    return sort_points(sec,sec2)

def surface_from_line_and_vector(line=[[0,0,0],[10,0,0]],vector=[0,0,1],both_sides=0):
    '''
    draw a surface base on a line and a vector.
    if surface is required both sides of the line, both_sides option should be marked as '1' else default for only one side is '0'
    '''
    l_1=translate(array(vector),line)
    l_2=translate(-array(vector),line) if both_sides==1 else line
    return [l_1,l_2]

def slice_surfaces(surf_1,surf_1_1,n=10):
    '''
    function creates multiple surfaces between 2 surfaces
    '''
    surf_x=cpo([slice_sol([surf_1[i],surf_1_1[i]],n)  for i in range(len(surf_1))])
    return surf_x

def mirror_surface(surf_1,n1,loc=[0,0,0]):
    '''
    function to mirror a solid or surface base on a mirroring plane given by vector 'n1'
    passing through a point 'loc'
    '''
    surf_1_1=[mirror_line(surf_1[i],n1,loc) for i in range(len(surf_1))]
    return surf_1_1

def solid_from_2surfaces(surf_1,surf_2):
    '''
    function to make a solid from 2 surfaces 
    '''
    sol=[surf_1[i]+flip(surf_2[i]) for i in range(len(surf_1))]
    # sol=cpo(sol)
    return sol

def convert_closed_section_to_surface(surf_1,s=1):
    '''
    function to convert a closed polygon to lines
    e.g.
    a=c2t3(circle(10))
    b=sec2surface(a)
    
    a and b can be visualised by following commands
    color("blue")p_line3d({a},.2,1);
    color("magenta")for(p={b})p_line3d(p,.2,1);
    {swp_surf(b)}
    '''
    return [equidistant_path([surf_1[i],surf_1[-i-1]] ,s)
            for i in range(int(len(surf_1)/2))]

def sec2lines(sec,n=20,s=10):
    '''
    function to convert a polygon to lines (horizontal lines)
    '''
    a=array(sec)[:,1].min()
    b=array(sec)[:,1].max()
    delta=(b-a)/1000
    p0=array([[0,i] for i in linspace(a+delta,b-delta,n)])
    # p0+v1*t1=p1+v2*t2
    pa=array(sec)
    pb=array(sec[1:]+[sec[0]])
    v1=array([[[1,0]]*len(pa)]*len(p0))
    v2=array([pb-pa]*len(p0))
    iim=array([v1,-v2+.000001]).transpose(1,2,0,3).transpose(0,1,3,2)
    im=inv(iim)
    p=pa[None,:,:]-p0[:,None,:]
    t=einsum('ijkl,ijl->ijk',im,p)
    dec=(t[:,:,1]>=0) & (t[:,:,1]<=1)

    pa.shape,v2.shape,t[:,:,1].shape
    c=array(lexicographic_sort_yx((pa[None,:,:]+einsum('ijk,ij->ijk',v2,t[:,:,1]))[dec])).reshape(-1,2,2).tolist()
    c=[lexicographic_sort_xy(equidistant_path(p,s))  for p in c]
    
    return c

def rotate_3d_sketch_parallel_to_xy_plane(sec):
    '''
    function to rotate any section open or closed parallel to x-y plane
    
    '''
    n1=nv(sec)
    if (array(n1).round(5).tolist()==[0,0,1]) | (array(n1).round(5).tolist()==[0,0,-1]) :
        return sec
    else:
        v1=cross(n1,[0,0,-1])
        t1=r2d(arccos(array(n1)@[0,0,-1]))

        l_3=axis_rot_1([sec],v1,sec[0],t1)[0]

        return l_3
    
def surround_3d(path,r,s=20):
    '''
    function to surround a path to create a rounded section on a 3d path
    
    '''
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
    '''
    plane defined by normal vector 'nv' with size as defined
    '''
    d1,d2=(size/2,size/2) if array(size).shape==() else (size[0]/2,size[1]/2)
    n1=array(nv)/norm(nv)
    if array(n1).round(5).tolist()==[0,0,1]:
        v1=[1,0,0]
    else:
        v1=c2t3(q_rot2d(90,c3t2(n1)))
    v1=array(v1)/norm(v1)
    v2=axis_rot(v1,[n1],90)[0]
    v1=array([array(v1)*-d1,array(v1)*d1])
    v2=array(v2)*d2
    s1=surface_line_vector(v1,v2,1)
    s1=translate(intercept,s1)
    return s1

def align_sec_2(sec1):
    '''
    function to align points of a section to obtain the non twisted optimised surface after using
    function sec2surface(sec1)
    '''
    
    area1=[array([l_lenv_o(p) for p in sec2surface(sec1[i:]+sec1[:i])]).sum() for i in range(len(sec1)) ]
    i=array(area1).argmin()
    sol2=sec1[i:]+sec1[:i]
    return sol2

def list_the_index_of_self_intersecting_segments(segments_list):
    '''
    calulates the self intersection list numbers of segment lists 'sec1'
    it picks the intersection points only if the 2 lines are crossing each other

    '''
    sec1=segments_list
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
    '''
    finds all the points, in a list of points 'sec2' which are inside the offset surround of an enclosed section 'sec'
    '''
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



def list_of_radius_at_each_point_of_a_closed_section(sec):
    '''
    function to list the radius at each point of a section
    '''
    x_6=list_r(sec)
    x_6[0]=x_6[1] if x_6[1]>0 else x_6[0]
    x_6[-1]=x_6[-2] if x_6[-2]>0 else x_6[-1]
    for i in range(1,len(x_6)-1):
        x_6[i]= x_6[i+1] if x_6[i+1]>0 else x_6[i-1] if x_6[i-1]>0 else x_6[i]
    return x_6

def exclude_numbers(a,b):
    '''
    exclude list of numbers 'b' from 'a'
    example:
    a=[1,2,3,4,5]
    b=[4,5,6,7]
    exclude_numbers(a,b) => array([1, 2, 3])
    '''
    a,b=array(a),array(b)
    return a[~(a[:,None]==b[None,:]).any(1)]

def subset(b,a):
    return (array(b)[:,None]==array(a)[None,:]).any(1).all()

def index_of_points_inside_offset_surround(sec,sec2,r):
    '''
    finds all the points list, in a list of points 'sec2' which are inside the offset surround of an enclosed section 'sec'
    '''
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
    '''
    function returns offset of a enclosed section 'sec' by a distance 'r'
    '''
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
    '''
    returns the vertices of a solid 'sol_1'
    '''
    return array(sol_1).reshape(-1,3)

def coil(r1,r2,n1=1):
    '''
    function to draw a coil with initial radius 'r1', final_radius 'r2' and numbe of coils 'n1'
    '''
    r_l=[i for i in linspace(r1,r2,360*n1)]
    theta_1=[i for i in linspace(0,360*n1,360*n1)]
    c1=[[r_l[i]*cos(d2r(theta_1[i])),r_l[i]*sin(d2r(theta_1[i]))] for i in range(len(r_l))]
    return c1

# def corner_n_radius_list(p0,r_l,n=10):
#     '''
#     corner list 'p0' and radius list 'r_l' will create a smothened section
#     'n' is the number of segments in each filleted corner
    
#     '''
#     # r_l=[.01 if i==0 else i for i in r_l]
#     p1=seg(p0)
#     p2=[p1[-1]]+p1[:-1]
#     p3=p1
#     s1=[fillet_intersection_lines(p2[i],p3[i],r_l[i],s=n) if r_l[i]!=0 else [p0[i]] for i in range(len(p1))]
#     s2=concatenate(s1).tolist()
#     return remove_extra_points(array(s2).round(5))

# def corner_n_radius_list_3d(p0,r_l,n=10):
#     '''
#     corner list 'p0' and radius list 'r_l' will create a smothened 3d path or closed section
#     'n' is the number of segments in each filleted corner
    
#     '''

#     p1=[p0[-1]]+p0[:-1]
#     s1=[]
#     for i in range(len(p0)):
#         if i<len(p0)-1 and r_l[i]>0:
#             a=fillet_3p_3d(p1[i],p0[i],p0[i+1],r_l[i],n)[1:]
#         elif i<len(p0)-1 and r_l[i]==0:
#             a=[p0[i]]
#         elif i==len(p0)-1 and r_l[i]==0:
#             a=[p0[i]]
#         elif i==len(p0)-1 and r_l[i]>0:
#             a=fillet_3p_3d(p1[i],p0[i],p0[0],r_l[i],n)[1:]
#         s1.append(a)
    
#     s1=remove_extra_points(concatenate(s1).round(5))
#     return s1

def path_extrude_closed(sec_1,path,twist=0):
    '''
    function to extrude a closed section to a closed path
    twist=0 for simple path extrudes
    set twist=1 for complex path extrudes
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
    '''
    function to extrude a closed section to an open path
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
    if type==1:
        return offset_1(sec,r)
    elif type==2:
        return offset_2(sec,r)

def sort_random_points(l_1,n_1,k=3):
    '''
    function to arrange random points in order
    l_1: list of random points in space
    n_1: is the normal vector to a plane from which all the points can be distinctly seen
    k: is a factor which can have values >2 , default is 2, if the result is not satisfactory the values can be changed to see if the result is better. It has to be an integer
    
    '''
    avg_1=array(l_1).mean(0).tolist()
    l_3=rot_sec2xy_plane(ppplane(l_1,n_1,avg_1))
    l_4=c3t2(l_3)
    l_5=concave_hull(l_4,k)
    l_6=array(l_1)[cKDTree(l_4).query(l_5)[1]].tolist()
    return l_6


def t_vec(path):
    '''
    find the array of tangent vectors to a given path
    '''
    p1=array(seg(path))
    p2=array(path)
    v1=array([(p[1]-p[0])/norm(p[1]-p[0]) for p in p1])
    t_v=array([ (v1[-1]+v1[i])/2 if i==0 else
         (v1[i-1]+v1[i])/2
        for i in range(len(p1))])
    t_v=t_v/norm(t_v,axis=1).reshape(-1,1)
    return t_v

def o_vec(path,n_v):
    '''
    finds the array of orthogal vectors to a given path.
    normal vector at each point needs to be defined for this calculation
    '''
    t_v=t_vec(path)
    o_v=array([cross(n_v[i],t_v[i]) for i in range(len(t_v))])
    o_v=o_v/norm(o_v,axis=1).reshape(-1,1)
    return o_v

def swp_sec(sec):
    '''
    function to render a sec
    '''
    n1=arange(len(sec)).tolist()
    
    return f'polyhedron({sec},{[n1]},convexity=10);'



def arc_with_start_point_and_center_3d(n1,start_point=[],center_point=[],theta=90,segments=30):
    '''
    function to draw an arc with known center_point and start_point in 3d space
    n1 is the normal which defines the plane
    '''
    center_point,start_point=array([center_point,start_point])
    v1=start_point-center_point
    arc_1=center_point+[axis_rot(n1,v1,i) for i in linspace(0,theta,segments)]
    arc_1=arc_1.tolist()
    return arc_1

def offset_3d(sec,d,type=1):
    '''
    offsets an enclosed section in 3d space, in case the section is in 1 plane
    sec: section in 3d space
    d: offset distance -ve sign means inner offset and +ve sign is outer offset
    refer to the file"example of various functions" for application examples
    type: offset type default is '1' in case of any issue in offset, try with '2'
    '''
    # sec=axis_rot_o([1,0,0],sec,.001)
    l_2=rot_sec2xy_plane(sec)
    l_3=c3t2(l_2)
    l_4=offset(l_3,d,type)
    avg_1=array(l_2).mean(0)
    avg_2=array(c2t3(l_3)).mean(0)
    l_5=translate(avg_1-avg_2,l_4)
    n_1=array(nv(sec)).round(5)
    n_2=array([0,0,-1])
    if (l_(n_1)==[0,0,1]) | (l_(n_1)==[0,0,-1]):
        return l_5
    else:
        ax_1=cross(n_1,n_2)
        theta=r2d(arccos(n_1@n_2))
        l_6=axis_rot_o(ax_1,[l_5],-theta)[0]
        l_6_1=axis_rot_o(ax_1,[l_2],-theta)[0]
        avg_1=array(sec).mean(0)
        avg_2=array(l_6_1).mean(0)
        l_7=translate(avg_1-avg_2,l_6)
        return l_7

def intersection_between_2_sketches(s1,s2):
    '''
    finds the intersection between 2 2d sketches
    '''
    a=seg(s1)
    b=seg(s2)
    c=a+b
    ip_1=s_int1(c)
    return ip_1

def pol(p1,l1):
    '''
    point on line
    finds whether a point is on a line or not
    returns True or False
    '''
    v1=l1[:-1]
    v2=l1[1:]
    v3=array(v2)-array(v1)
    v4=array(p1)-array(v1)
    u3=(v3/norm(v3,axis=1).reshape(-1,1)).round(4)
    u4=(v4/norm(v4,axis=1).reshape(-1,1)).round(4)
    dec=(array(p1).round(4)==array(l1[0]).round(4)).all()
    return True if dec==True else (u3[:,None,:]==u4[None,:]).all(2).any()



def intersection_points_projection_on_intersecting_solid(surf_1,i_p_l,d=1.):
    '''
    function to project the intersection point on the cutting lines based on the distance 'r'
    '''
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



def path_extrude_over_multiple_sec_open(sec_1,path,twist=0):
    '''
    function to extrude multiple closed sections to an open path
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

def normal_vector_to_any_vector_in_xy_plane(v1):
    '''
    returns a normal vector to any vector in x-y plane
    '''
    u1=v1/norm(v1)
    ua=array([0,0,-1]) if u1[2]==0 else array([0,-1,0]) if (u1==[0,0,1]).all() else array([-1,0,0]) if (u1==[0,0,-1]).all() else array([u1[0],u1[1],0])
    v2=cross(u1,ua)
    u2=v2/norm(v2)
    u3=axis_rot(u2,u1,-90)
    u1,u2,u3=array([u1,u2,u3]).tolist()
    return u2

def normal_vector_to_any_vector(v1):
    '''
    returns a normal vector to any vector in +/- z direction
    '''
    u1=v1/norm(v1)
    ua=array([0,0,-1]) if u1[2]==0 else array([0,-1,0]) if (u1==[0,0,1]).all() else array([-1,0,0]) if (u1==[0,0,-1]).all() else array([u1[0],u1[1],0])
    v2=cross(u1,ua)
    u2=v2/norm(v2)
    u3=cross(u2,u1)
    u1,u2,u3=array([u1,u2,u3]).tolist()
    return u3

def s_p(p_l): # starting point
        '''
        find the starting point for a convex hull
        bottom left point
        '''
        l_1=array(p_l).round(5)
        a=l_1[l_1[:,1].argsort()]
        if len(a[a[:,1]==a[:,1].min()])>1:
            b=a[a[:,1]==a[:,1].min()]
            s_pnt=b[b[:,0].argsort()][0]
        else:
            s_pnt=a[0]
        return s_pnt.tolist()


def convert_closed_section_to_surface_minimum_sum_method(sec1,s=1):
    '''
    create an aligned surface from a section 'sec1'
    considers the min sum method
    '''
    sec2=[sec2surface(sec1[i:]+sec1[:i],s) for i in range(len(sec1))]
    sec2=[sum([l_len(p1) for p1 in p]) for p in sec2]
    n2=array(sec2).argmin()
    sec3=sec2surface(sec1[n2:]+sec1[:n2],s)
    return sec3

def convert_closed_section_to_surface_maximum_sum_method(sec1,s=1):
    '''
    create an aligned surface from a section 'sec1'
    considers the max sum method
    '''
    sec2=[sec2surface(sec1[i:]+sec1[:i],s) for i in range(len(sec1))]
    sec2=[sum([l_len(p1) for p1 in p]) for p in sec2]
    n2=array(sec2).argmax()
    sec3=sec2surface(sec1[n2:]+sec1[:n2],s)
    return sec3
    
def vector2length(v,l=10):
    '''
    draw a defined vector to length 'l'
    '''
    u=array(v)/norm(v)
    v1=(u*l).tolist()
    return v1

def tangent_on_circle_from_any_point_on_circle(c,p,l=1,side=0):
    '''
    function to draw a tangent on a circle 'c' from any given
    point 'p' on the circle. length of the tangent 'l'
    side can be set to either '0' or '1' to draw tangents in
    2 different directions
    '''
    v1=array(p)-array(cp_arc(c))
    theta1=r2d(arctan(l/norm(v1)))
    theta2=ang(v1[0],v1[1])
    l2=norm(v1)/cos(d2r(theta1))
    pa=translate(cp_arc(c),q_rot2d(theta1+theta2,[l2,0])) if side==0 \
    else translate(cp_arc(c),q_rot2d(-theta1+theta2,[l2,0]))
    return [p,pa]

def l_(a):
    '''
    convert an array to list
    '''
    return array(a).tolist()

def a_(l):
    '''
    convert a list to array
    '''
    return array(l)

def convert_to_triangles(sol):
    '''
    convert a solid to a triangular mesh
    '''
    f1=faces_1(len(sol),len(sol[0]))
    v1=vertices(sol)
    return l_(v1[f1])

def match_2_points_list(s0,s1):
    '''
    match 2 sets of list of points or sections in space without loosing any point from both the lists
    
    '''
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
    '''
    convert a surface to a triangular mesh
    '''
    f1=faces_surface(len(surf),len(surf[0]))
    v1=vertices(surf)
    return l_(v1[f1])

def lb2p(p0,p1):
    '''
    line between 2 points,
    it is a perpendicular bisector of the line between 2 points
    '''
    l1=l_len([p0,p1])
    c1=circle(l1,p0)
    c2=circle(l1,p1)
    l2=s_int1(seg(c1)+seg(c2))
    return l2

def lexicographic_sorting_of_points(pnts=[],seq=[0,1,2],ord=[1,1,1]):
    '''
    lexicographic ordering of a points list
    seq: defines the seduence in which the points needs to be ordered
    e.g. [0,1,2] means first on x then on y and lastly on z coordinates
    order: means asceding or descending order '1' means ascending and '-1' means descending order
    e.g. [1,1,1] means all the coordinates should be in ascending order
    '''
    if len(pnts[0])==2:        
        return sorted(pnts,key=lambda x:(ord[0]*x[seq[0]],ord[1]*x[seq[1]]))
    elif len(pnts[0])==3:
        return sorted(pnts,key=lambda x:(ord[0]*x[seq[0]],ord[1]*x[seq[1]],ord[2]*x[seq[2]]))

def fillet_from_3_points(p1,p2,p3,s=10):
    '''
    creates a fillet with 3 defined points
    '''
    l2=[p1,(p1+p2)/2,((p1+p2)/2+(p2+p3)/2)/2,(p2+p3)/2,p3]
    arc1=bezier(l2,s)
    return arc1

# def convert_surface_to_section(surf):
#     '''
#     reverse of sec2surface_1 function
#     '''
#     a,b,c=array(surf).transpose(1,0,2).shape
#     d=l_(array(surf).transpose(1,0,2))
#     e=[d[i]+flip(d[-i-1]) for i in range(int(a/2))]
#     return e

def convert_surface_to_section(surf):
    return surf[0]+flip(surf[-1])


def join_arcs_closed(pl,rl,ol,ls,s):
    '''
    join various arcs
    pl: points list for the arcs
    rl: radius list for the arcs
    ol: orientation list i.e. clockwise (1) or counter clockwise(-1)
    ls: long(1) or small(0) arcs list
    s: number of segments of the list
    '''
    a=seg(pl)
    a_l=[]
    for i in range(len(a)):
        a_l.append(arc_2p(a[i][0],a[i][1],rl[i],ol[i],s) if ls[i]==0 else
                  arc_long_2p(a[i][0],a[i][1],rl[i],ol[i],s))
    a_l=equidistant_path(l_(concatenate(a_l)),s)
    return a_l

def join_arcs_open(pl,rl,ol,ls,s):
    '''
    join various arcs
    pl: points list for the arcs
    rl: radius list for the arcs
    ol: orientation list i.e. clockwise (1) or counter clockwise(-1)
    ls: long(1) or small(0) arcs list
    s: number of segments of the list
    '''
    a=seg(pl)[:-1]
    a_l=[]
    for i in range(len(a)):
        a_l.append(arc_2p(a[i][0],a[i][1],rl[i],ol[i],s) if ls[i]==0 else
                  arc_long_2p(a[i][0],a[i][1],rl[i],ol[i],s))
    a_l=equidistant_path(l_(concatenate(a_l)),s)
    return a_l

def bspline_open(pl,deg=3,s=100):
    '''
    draws bspline curve for control points list 'pl'
    degree of curve 'deg' which should lie 0<=deg<len(p0)
    s: number of points in the resultant curve
    '''
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
    '''
    draws bspline closed loop curve for control points list 'pl'
    degree of curve 'deg' which should lie 0<=deg<len(p0)
    s: number of points in the resultant curve
    '''
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
    '''
    draws bspline surface from 2 control points list 'pl1' and 'pl2'
    degree of curves are 'deg1' and 'deg2' in 2 directions 
    which should lie 0<=deg1<len(pl[0]) for 'deg1' and
    0<=deg2<=len(pl) for 'deg2'
    s: number of points in the resultant curve
    a: if a[0]==1 means set 1st shape of the surface to closed else to open and similarly for a[1]. open here means function bspline_open and closed means bspline_closed
    '''
    p1=[bspline_closed(p,deg1,s1) if a[0]==1 else bspline_open(p,deg1,s1) for p in pl]
    p2=cpo([bspline_closed(p,deg2,s2) if a[1]==1 else bspline_open(p,deg2,s2) for p in cpo(p1)])
    return p2

def bspline_sol(pl,deg1=3,deg2=3,s1=100,s2=100,a=[1,1]):
    '''
    draws bspline surface from 2 control points list 'pl1' and 'pl2'
    degree of curves are 'deg1' and 'deg2' in 2 directions 
    which should lie 0<=deg1<len(pl[0]) for 'deg1' and
    0<=deg2<=len(pl) for 'deg2'
    s: number of points in the resultant curve
    a: if a[0]==1 means set 1st shape of the surface to closed else to open and similarly for a[1]. open here means function bspline_open and closed means bspline_closed
    '''
    p1=[bspline_closed(p,deg1,s1) if a[0]==1 else bspline_open(p,deg1,s1) for p in pl]
    p2=cso([bspline_closed(p,deg2,s2) if a[1]==1 else bspline_open(p,deg2,s2) for p in cso(p1)])
    return p2

def polp(l,t):# point on line parameteric
    '''
    find a point on line 'l' at parameter 't' where 0<=t<=1
    '''
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


def lineFromPointTillEnd(l1,pnt,dist=.01):
    '''
    draws line from point 'pnt' till end of line 'l1'
    the point with in distance 'dist' will be picked up
    deletes the rest of the line
    '''
    def lfpte(l1,t):
        '''
        removes the line till parameter t where 0<=t<=1
        keeps the rest of the line
        '''
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
    '''
    draws line from start till point 'pnt'.
    the point with in distance 'dist' will be picked up
    deletes the rest of the line
    '''
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
    '''
    finds the projection of the point 'p0' on line 'l1' which is within distance 'dist' from the line
    '''
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

def vcost_without_bounds(l1,p0,dist=.2):
    '''
    finds the projection of the point 'p0' on line 'l1' which is within distance 'dist' from the line
    '''
    v1=a_(l1[1])-a_(l1[0])
    u1=v1/norm(v1)
    d=norm(v1)
    v2=a_(p0)-a_(l1[0])
    d1,d2=v1@v2/norm(v1),norm(cross(c23(v1),c23(v2)))/norm(v1)
    if abs(d2)<=dist:
        p1=l_(a_(l1[0])+u1*d1)
        t1=d1/d
    else:
        p1=[]
    return p1

def timeToReachPoint(p0,l1,dist=.01):
    '''
    p0: point
    l1: line
    dist: poin with in distance 'dist' will be picked up
    if the point is on the line, it gives back time
    where 0<=time<=1
    '''
    l2=seg(l1)[:-1]
    n1=[ i for i in range(len(l2)) if vcost(l2[i],p0,dist)!=[]]
    p0=[] if n1==[] else vcost(l2[n1[0]],p0,dist)
    a=a_(seg(l1)[:-1])
    b,c=a[:,0],a[:,1]
    d=p0-b
    d=d/norm(d,axis=1).reshape(-1,1)
    e=c-b
    e=e/norm(e,axis=1).reshape(-1,1)
    n=arange(len(a))[(d.round(2)==e.round(2)).all(1)][0]
    d1=l_len([p0,l1[n]])+l_lenv_o(l1[:n+1])
    t1=d1/l_lenv_o(l1)
    return t1

def movePointOnLine(l1,p0,d):
    '''
    move a point 'p0' on line 'l1' by a distance 'd'
    '''
    if p0!=l1[0]:
        t1=d/l_lenv_o(l1)
        t2=timeToReachPoint(p0,l1)
        t3=t1+t2
        p1=polp(l1,t3)
        return p1
    elif p0==l1[0]:
        t1=d/l_lenv_o(l1)
        return polp(l1,t1)

def find_triangular_meshes_of_intersecting_points_on_surface(ip,surface):
    '''
    function to find the triangles on the solid 'sol1' where the intersection points list 'ip' lies
    '''
    sol1=surface
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
    decision=(t>=0)&(t<=1)&(u>=0)&(u<=1)&(v>=0)&(v<=1)&((u+v)<=1)
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
    '''
    convert a points list from 2d to 3d coordinate system
    '''
    return c2t3(pl)

def c32(pl):
    '''
    convert a points list from 3d to 2d coordinate system
    '''
    return c3t2(pl)

def prism2cpo(s1):
    '''
    change the orientation of the surface to rectangular
    '''
    s2=[cpo(s1)[i]+flip(cpo(s1)[-(i+1)]) for i in range(int(len(s1[0])/2))]
    a=[mid_point([s2[0][i],s2[0][-(i+1)]]) for i in range(int(len(s2[0])/2))]
    b=[mid_point([s2[-1][i],s2[-1][-(i+1)]]) for i in range(int(len(s2[-1])/2))]
    
    a=a+flip(a)
    b=b+flip(b)
    s2=[a]+s2+[b]
    return s2


def project_surface_on_to_another_surface(surface_on_which_to_project, surface_to_project, vector_for_projection, distance_limit_for_projection=100000, unidirection=1):
    '''
    project a surface on to another without loosing the original points
    surface 's3' will be projected on surface 's2'
    'v1' is vector for projection
    'dist' is the maximum distance through which projection can happen
    unidirection: if the projection is to be in both direction set parameter
    unidirection to '1' else to '0'
    '''
    s2,s3,v1=surface_on_which_to_project, surface_to_project, vector_for_projection
    dist=distance_limit_for_projection
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

def reorient_sec(sec):
    '''
    re-orient section to create a better surface through 'prism' function
    '''
    sec1=sec2surface_1(sec)
    sec2=[p[0] for p in sec1]+flip([p[-1] for p in sec1])
    return sec2

def project_surface_on_to_another_surface_with_focal_point (surface_on_which_to_project, surface_to_project, focal_point, distance_limit_for_projection=100000, unidirection=0):
    '''
    project a surface on to another without loosing the original points
    surface 's3' will be projected on surface 's2'
    'v1' is vector for projection. this is a focal vector 
    from where the rays are emitted for projection
    '''
    s2,s3,v1,dist=surface_on_which_to_project, surface_to_project, focal_point, distance_limit_for_projection
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


def hyd_cyl(bore_dia,stroke_length,piston_thk,plunger_length,plunger_dia,top_bottom_thk,side_thk,int_position):
    a,b,c,d,e,f,g,h=bore_dia,stroke_length,piston_thk,plunger_length,plunger_dia,int_position,top_bottom_thk,side_thk
    
    sec2=corner_radius(pts1([[-(a+2*h)/2,-(a+2*h)/2,2*h],[a+2*h,0,2*h],[0,a+2*h,2*h],[-(a+2*h),0,2*h]]),9)
    sec1=circle(a/2,s=len(sec2)+1)
    path1=pts([[0,g],[0,b+c],[-a/2+e/2,0],[0,g]])
    sol1=prism(sec1,path1)
    sol2=linear_extrude(sec2,c+b+2*g)
    cyl=sol2+flip(sol1)
    cyl=align_sol_1(cyl)
    path2=pts([[0,g],[0,c],[(-a+e)/2,0],[0,d]])
    piston=prism(sec1,path2)
    l1=polp([[0,0,0],[0,0,b]],f)
    piston=translate(l1,piston)
    return [cyl,piston]
    
def sleeve(sleeve_dia,flange_dia,sleeve_length,thickness,flange_thickness):
    a,b,c,d,e=sleeve_dia,flange_dia,sleeve_length,thickness,flange_thickness
    sec=circle(a/2)
    path=corner_radius(pts1([[0,0,e/2],[(b-a)/2,0],[0,e],[-(b-a)/2+d,0,d/2],[0,c-e],[-d,0,d/2]]),10)
    sleeve=prism(sec,path)
    return sleeve+[sleeve[0]]

def normal_vector(sec):
    '''
    calculates the normal vector of a given sec
    '''
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

def extrude_sec2path_3d(sec,path):
    '''
    extrude a 3d section to a 3d path
    '''
    s1=sec
    p2=path
    n1=a_(normal_vector(s1))
    p3=tangents_along_path(p2)[1:-1]
    pa,pb=[p2[0],l_(a_(p2[0])+(a_(p2[1])-a_(p2[0])))], \
    [p2[-1],l_(a_(p2[-1])+(a_(p2[-1])-a_(p2[-2])))]
    p3=a_([pa]+p3+[pb])
    s2=[]
    for i in range(len(p3)):
        n2=p3[i][1]-p3[i][0]
        n2=n2/norm(n2)
        if (n1.round(4)==n2.round(4)).all():
            s2.append(translate(a_(p2[i])-a_(p2[0]),s1))
        else:
            a1=cross(n2,n1)
            theta=180-r2d(arccos(n1@n2))
            s2.append(translate(a_(p2[i])-a_(p2[0]),axis_rot_1(s1,a1,p2[0],theta)))
    

    return s2

def surface_normals(surf,direction=1):
    '''
    calculates normals from each point on the surface
    normals may need to be reshaped to the shape of the 'surf'
    '''
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

def project_surface_on_to_another_with_surface_normals( surface_on_which_to_project, surface_to_project, direction=1, dist=100000):
    '''
    project a surface on to another without loosing the original points
    surface 's3' will be projected on surface 's2'
    the projection is based on the normal to the surface s3 and the direction of normals can be changed from '1' to '-1'
    '''
    s2=surface_on_which_to_project
    s3=surface_to_project
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


def cog(sol):
    '''
    calculate the center of gravity
    '''
    return l_(a_(sol).reshape(-1,3).mean(0))


def interpolation_bspline_closed(p0,s=50,f=3.425):
    '''
    Interpolates through points
    p0: original points
    s: number of points in the resultant curve
    f: factor of smoothness, lower values may result in too wavyness and larger values
    are as good as connecting points in straight line
    
    '''
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
    '''
    Interpolates through points
    p0: original points
    s: number of points in the resultant curve
    f: factor of smoothness, lower values may result in too wavyness and larger values
    are as good as connecting points in straight line
    
    '''
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
    '''
    draws bspline interpolation surface from 2 control points list 'pl1' and 'pl2'
    'f1' and 'f2' are factors of smoothness in 2 directions 
    s1 and s2: number of points in 2 direction in the resultant curve
    '''
    p1=[interpolation_bspline_open(p,s1,f1) for p in pl]
    p2=cpo([interpolation_bspline_open(p,s2,f2) for p in cpo(p1)])
    return p2

def bezier_surface(pl,s1=100,s2=100):
    '''
    draws bezier surface from control points list 'pl'
    s1 and s2: number of points in the resultant curve in 2 direction
    '''
    p1=[bezier(p,s1) for p in pl]
    p2=cpo([bezier(p,s2) for p in cpo(p1)])
    return p2

    
def intersection_line_fillet(il,sol1,sol2,r1,r2,s=20,o=0,type=1,dist=0,vx=[],style=2,f=1):
    '''
    calculates a fillet at the intersection of 2 solids.
    r1 and r2 would be same in most of the cases, but the signs can be different depending on which side the fillet is required
    r1 is the distance by which intersection line offsets on sol1 and similarly r2 is on sol2 
    for type, dist and vx parameters refer function o_3d_rev
    '''
        
    p1=il
    p2=o_3d_rev(p1,sol2,r2,type=type,dist=dist,vx=vx,edges_closed=1)
    p3=o_3d_rev(p1,sol1,r1,type=type,dist=dist,vx=vx,edges_closed=1)
    fillet1=convert_3lines2fillet(p3,p2,p1,s=s,style=style,f=f)
    
    return fillet1

def intersection_line_fillet_surface(il,surf1,sol2,r1,r2,s=20,o=0,type=1,dist=0,vx=[],style=2,f=1):
    '''
    calculates a fillet at the intersection of surface with solid.
    r1 and r2 would be same in most of the cases, but the signs can be different depending on which side the fillet is required
    r1 is the distance by which intersection line offsets on sol1 and similarly r2 is on sol2 
    for type, dist and vx parameters refer function o_3d_rev
    '''
        
    p1=il
    p2=o_3d_rev(p1,sol2,r2,type=type,dist=dist,vx=vx,edges_closed=1)
    p3=o_3d_rev(p1,surf1,r1,type=type,dist=dist,vx=vx,edges_closed=0)
    fillet1=convert_3lines2fillet(p3,p2,p1,s=s,style=style,f=f)
    
    return fillet1

def line2length(l1,length=1):
    '''
    change the length of a line to the parameter distance 'length' keeping the initial point same
    '''
    v1=a_(l1[1])-a_(l1[0])
    u1=v1/norm(v1)
    return [l1[0],l_(l1[0]+u1*length)]

def ip_surf2line(surf,line):# when line has more than 2 points
    '''
    function to calculate intersection point between a 3d surface and a line. 
     "surf" is the 3d object which is intersected with a "line".
    
    '''


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
    i_p2=i_p1[argsort([norm(p-px[0]) for p in i_p1])]

    return i_p2.tolist()

def line_as_axis(l1):
    '''
    convert line to an axis vector
    '''
    v1=l_(a_(l1[1])-a_(l1[0]))
    return v1

def project_surface_on_to_another_with_list_of_vectors ( surface_on_which_to_project ,surface_to_project ,list_of_vectors_for_projection_direction ,dist=100000, unidirection=0):
    '''
    project a surface on to another without loosing the original points
    surface 's3' will be projected on surface 's2'
    'v1' are the vectors for projection. v1 are the list of vectors
    '''
    s2=surface_on_which_to_project
    s3=surface_to_project
    v1=list_of_vectors_for_projection_direction
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
    '''
    function to rotate an object by axis and angles defined by a string 'a'
    example:
    rot(a='z30x70y-90',sol=[[0,0],[10,0]])
    will rotate line [[0,0],[10,0]] first through z-axis by 30 deg 
    then through x-axis by 70 deg and then through y-axis by -90 deg 
    '''
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
    '''
    function to rotate a 2d section by a defined angle 'a'
    '''
    return c32(axis_rot([0,0,1],sec,a))


def sub_d_2(b,a=[1/4,3/4]):
    '''
    input to function fillet_by_subdivison
    '''
    c=[[[polp(p1,a[0]),polp(p1,a[1])] for p1 in seg(p)[:-1]] for p in cpo(b)]
    n1,n2,n3,n4=a_(c).shape
    c=[b[0]]+cpo(a_(c).reshape(n1,n2*n3,n4))+[b[-1]]
    
    return c

def fillet_by_subdivison(sol,a=[1/4,3/4]):
    '''
    produce fillet by subdividing an object
    '''
    for _ in range(4):
        sol=sub_d_2(sol,a)
    return sol

def cso(sol):
    '''
    change the orientation of points in solid
    '''
    a=cpo(sol)
    b=[ a[i]+flip(a[-(i+1)]) for i in range(int(len(a)/2))]
    return b

def offset_3d_sec(path,dist=1,o=0):
    '''
    offsets any random closed 3d path by defined distance
    '''
    path1=equidistant_pathc(path,pitch=.1)
    path2=rationalise_path(path1)
    if o==0:
        s1=sec2surface_1(path2,1)
    elif o==1:
        s1=sec2surface(path2,1)
    l2=o_3d_surf(path1,s1,dist)
    l1=list_remove_points_within_dist_closed(path,l2,dist-.02)
    l4=[l1[cKDTree(a_([l1,zeros(len(l1))]).transpose(1,0)).query([i,0])[1]] for i in range(len(l2))]
    l5=[cKDTree(path1).query(path[i])[1] for i in range(len(path))]
    l6=a_(l4)[l5]
    path1=l_(a_(l2)[l6])
    return path1


def remove_points_within_dist_closed(path,line,dist=1):
    '''
    function to remove points which are with in defined distance of a given path
    path : is the original path
    line: is the line from where the points needs to be removed
    dist: is the distance with in which all the points to be excluded
    '''

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
    '''
    function to find list of points to remove which are with in defined distance of a given path
    path : is the original path
    line: is the line from where the points needs to be removed
    dist: is the distance with in which all the points to be excluded
    '''

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
    '''
    function to remove points which are with in defined distance of a given path
    path : is the original path
    line: is the line from where the points needs to be removed
    dist: is the distance with in which all the points to be excluded
    '''

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
    '''
    function to find list of points to remove which are with in defined distance of a given path
    path : is the original path
    line: is the line from where the points needs to be removed
    dist: is the distance with in which all the points to be excluded
    '''

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

def project_a_point_on_surface(surf,pnt,vect,unidirection=1):#project point on surface
    '''
    function to project a point on a surface
    'vect' is the vector in which direction the line would be projected.
    'unidirection' to be set to '1' if the projection is required in the 
    direction of the vector and set to '0' in case it is required in either of 
    the direction. 
    Projections happens at the nearest location on the surface from a line
    '''
    return psos(surf,[[pnt]],vect,unidirection=unidirection)[0][0]

def project_a_line_on_surface(surf,line,vect,unidirection=1):#project line on surface
    '''
    function to project a line on a surface
    'vect' is the vector in which direction the line would be projected.
    'unidirection' to be set to '1' if the projection is required in the 
    direction of the vector and set to '0' in case it is required in either of 
    the direction. 
    Projections happens at the nearest location on the surface from a line
    '''
    return psos(surf,[line],vect,unidirection=unidirection)[0]

def barycentric_normals(surf,edges_closed=1):
    '''
    finds normals for each triangle mesh of a surface
    if it is a closed loop surface like cylinder, 
    parameter edges_closed to be set to '1' else set to '0'
    '''
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
    '''
    finds center for each triangle mesh of a surface
    if it is a closed loop surface like cylinder, 
    parameter edges_closed to be set to '1' else set to '0'
    '''
    if edges_closed==1:
        f=faces_1(len(surf),len(surf[0]))
    elif edges_closed==0:
        f=faces_surface(len(surf),len(surf[0]))
    v=a_(surf).reshape(-1,3)
    triangles=v[f]
    cp=triangles.mean(1)
    return cp

def horizontal_lines_on_closed_section(sec,n=10):
    '''
    divide a 2d section 'sec' in to horizontal lines
    'n' is the number of lines in the section
    '''
    py=m_points1_o([a_(sec)[:,1].min()+.001,a_(sec)[:,1].max()-.001],n)
    px=zeros(len(py))
    p2=a_([px,py]).transpose(1,0)
    a=seg(sec)
    p0=a_(a)[:,0]
    p1=a_(a)[:,1]
    v1=p1-p0
    v2=a_([[1,0]]*len(v1))
    # p0+v1*t1=p2+v2*t2
    # v1*t1-v2*t2=p2-p0
    iim=a_([v1,-v2+.000001]).transpose(1,0,2).transpose(0,2,1)
    im=inv(iim)
    t=einsum('ijk,hik->hij',im,p2[:,None,:]-p0[None,:])
    p=p0+einsum('ij,hi->hij',v1,t[:,:,0])
    dec=(t[:,:,0]>=0) & (t[:,:,0]<=1)
    p1=l_(a_(seg(lexico(p[dec],[1,0],[1,1])))[::2])
    return p1

def vertical_lines_on_closed_section(sec,n=10):
    '''
    divide a 2d section 'sec' in to verticle lines
    'n' is the number of lines in the section
    '''
    px=m_points1_o([a_(sec)[:,0].min()+.001,a_(sec)[:,0].max()-.001],n)
    py=zeros(len(px))
    p2=a_([px,py]).transpose(1,0)
    a=seg(sec)
    p0=a_(a)[:,0]
    p1=a_(a)[:,1]
    v1=p1-p0
    v2=a_([[0,1]]*len(v1))
    # p0+v1*t1=p2+v2*t2
    # v1*t1-v2*t2=p2-p0
    iim=a_([v1,-v2+.000001]).transpose(1,0,2).transpose(0,2,1)
    im=inv(iim)
    t=einsum('ijk,hik->hij',im,p2[:,None,:]-p0[None,:])
    p=p0+einsum('ij,hi->hij',v1,t[:,:,0])
    dec=(t[:,:,0]>=0) & (t[:,:,0]<=1)
    p1=l_(a_(seg(lexico(p[dec],[0,1],[1,1])))[::2])
    return p1

def smoothening_by_subdivison(sec,iterations=4,closed=0):
    '''
    smoothen a polyline by subdivision method.
    if the polyline is closed loop, parameter 'closed' should be set at '1'
    else if the polyline is not a closed loop parameter 'closed' to be set to '0'
    '''
    for i in range(iterations):
        if closed==1:
            sec=l_(concatenate([[polp(p1,1/4),polp(p1,3/4)] for p1 in seg(sec)]))
        elif closed==0:
            sec=[sec[0]]+l_(concatenate([[polp(p1,1/4),polp(p1,3/4)] for p1 in seg(sec)[:-1]]))+[sec[-1]]
    return sec


def smoothening_by_subdivison_surf(sol,iterations=4,o=[0,0]):
    '''
    smoothen a solid with sub-divison method
    iterations: number of iterations
    o: if o[0]==0 (means loop remains open for the 1st orientation, and '1' means the loop remains closed for the 1st orientation of the shape) 
    similarly o[1]==0 (means loop remains open for the 2nd orientation, and '1' means the loop remains closed for the 2n orientation of the shape)
    an example can make this more clear
    '''
    sol=[smoothening_by_subdivison(p,iterations,o[0]) for p in sol]
    sol=cpo([smoothening_by_subdivison(p,iterations,o[1]) for p in cpo(sol)])
    return sol

def wrap_x(l1,path):
    
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
    '''
    gives thickness to a surface by amount 't'
    '''
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


def swp_triangles(sol):
    '''
    function to render triangular mesh
    '''
    vert1=l_(a_(sol).reshape(-1,3))
    f1=l_(arange(len(sol)*len(sol[0])).reshape(-1,3))
    return f'polyhedron({vert1},{f1},convexity=10);'

def surface_split(sol,v1):
    '''
    function to split a 3d mesh surface, based on vector 
    example, if you are looking at the object from the top use vector [0,0,1]
    if from right [1,0,0] from left [-1,0,0], from front [0,1,0] , from back [0,-1,0] and so on...
    '''
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
    '''
    reshape a surface with a point 'pnt'
    '''
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
    '''
    Draw a line from a defined point to another point on a line.
    points with in distance 'dist' will be picked up
    Rest of the line will be deleted
    
    '''
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
    '''
    calculates the perpendicular distance of a point from a polyline
    '''
    px=[vcost(p,p0,1e5) for p in seg(l1)[:-1] if vcost(p,p0,1e5)!=[]]
    return l_len([px[0],p0])  if px!=[] else []

def perpendicularProjectionOfPointOnLine(p0,l1):
    '''
    finds the perpendicular projection of a point on a polyline, if one exists
    '''
    px=[vcost(p,p0,1e5) for p in seg(l1)[:-1] if vcost(p,p0,1e5)!=[]]
    return px[0]  if px!=[] else []

def switch_orientation(sol1):
    '''
    switch the orientation of the solid
    '''
    if len(sol1[0])%2==1:
        sol1=[equidistant_path(p,len(sol1[0])) for p in sol1]
    sol2=a_(sol1).transpose(1,0,2)
    x,y,z=a_(sol2).shape
    c=a_(sol2).reshape(2,int(x/2),y,3)[0]
    d=a_(sol2).reshape(2,int(x/2),y,3)[1][::-1]
    d=d[::-1].transpose(1,0,2)[::-1].transpose(1,0,2)[::-1]
    return l_(a_([c,d]).transpose(1,0,2,3).reshape(-1,y*2,z))

def intersection_line_fillet (intersection_line,solid_1, solid_2, distance_1=1, distance_2=1, segments=20):
    '''
    if the intersection line is defined between 2 solids. fillet can be drawn with the
    information
    '''
    l1,sol1,sol2,r1,r2=intersection_line,solid_1,solid_2,distance_1,distance_2
    l2=o_3d_rev(l1,sol1,r1)
    l3=o_3d_rev(l1,sol2,r2)
    f1=convert_3lines2fillet(l2,l3,l1,s=segments)
    return f1

def intersection_line_fillet_closed (intersection_line,solid_1, solid_2, distance_1=1, distance_2=1, segments=20):
    '''
    if the intersection line is defined between 2 solids. fillet (Closed loop) can be drawn with the
    information
    '''
    l1,sol1,sol2,r1,r2=intersection_line,solid_1,solid_2,distance_1,distance_2
    l2=o_3d_rev(l1,sol1,r1)
    l3=o_3d_rev(l1,sol2,r2)
    f1=convert_3lines2fillet_closed(l2,l3,l1,s=segments)
    return f1

def mid_line(l1,l2):
    '''
    draw a line in center of 2 lines 'l1' and 'l2'
    '''
    return [mid_point(p) for p in cpo([l1,l2])]

def list_ang(l1):
    '''
    list of angles for each point for a closed loop section
    '''
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
    '''
    remove all the points from the offset surface which are less than the offset distance 'dist'
    '''
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
    '''
    remove all the points from the offset surface which are less than the offset distance 'dist'
    '''
    sol2=[]
    for i in range(len(surf_original)):
        d=remove_points_within_dist_open(surf_original[i],surf_off[i],dist-.01)
        if len(d)==len(surf_original[i]):
            sol2.append(d)
        else:
            d=path2path1(surf_original[i],d)
            sol2.append(d)
    return sol2

def project_a_line_on_to_a_surface_with_list_of_vectors (surface,line, list_of_vectors_for_projection, dist=100000, unidirection=0):
    '''
    project a line on to a surface without loosing the original points
    line 'l1' will be projected on surface 's2'
    'v1' are the vectors for projection. v1 are the list of vectors
    '''
    s2=surface
    l1=line
    v1=list_of_vectors_for_projection
    return psos_v_2(s2,[l1],v1,dist=dist,unidirection=unidirection)[0]

def project_a_line_on_to_a_surface_with_focal_line (surface ,line ,focal_line , dist=100000, unidirection=0):
    '''
    project a line on to a surface without loosing the original points
    'line' will be projected on 'surface'
    focal_line is the line from where the rays are emitted for projection
    '''
    s2,l1,l2=surface ,line ,focal_line
    return project_surface_on_to_another_surface_with_focal_line(s2,[l1],l2,dist=dist,unidirection=unidirection)[0]


def project_a_line_on_to_a_surface_with_focal_point(surface ,line ,focal_point ,dist=100000,unidirection=0):
    '''
    project a line on to a surface without loosing the original points
    line 'l1' will be projected on surface 's2'
    'v1' is vector for projection. this is a focal vector 
    from where the rays are emitted for projection
    '''
    s2,l1,v1=surface, line, focal_point
    return psos_v(s2,[l1],v1,dist=dist,unidirection=unidirection)[0]

def lines2vectors(lines):
    '''
    convert lines to vectors
    '''
    return [line_as_axis(p) for p in lines]

def mirror_point(pnt,n1,loc):
    '''
    function to mirror a point 'pnt' defined by mirroring plane 'n1' passing through intercept 'loc'
    '''
    return mirror_line([pnt],n1,loc)[0]

def corner_radius3d(pnts,s=5): # Corner radius 3d where 'pnts' are the list of points with 4th coordinate in each point is radius 'rds' and 's' is number of segments for each arc
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

def corner_radius3d_with_turtle(pnts,s=5): # Corner radius 3d where 'pnts' are the list of points with 4th coordinate in each point is radius 'rds' and 's' is number of segments for each arc
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



# def corner_and_radius3d(pnts,rds,s=5): # Corner radius 3d where 'pnts' are the list of points, 'rds' are the list of radiuses and 's' is number of segments for each arc

#     c=[]
#     for i in range(len(pnts)):
#         if i==0:
#             p0=pnts[len(pnts)-1]
#             p1=pnts[i]
#             p2=pnts[i+1]
#         elif i<len(pnts)-1:
#             p0=pnts[i-1]
#             p1=pnts[i]
#             p2=pnts[i+1]
#         else:
#             p0=pnts[i-1]
#             p1=pnts[i]
#             p2=pnts[0]
#         c.append(fillet_3p_3d(p0,p1,p2,rds[i],s)[1:])
#     c=array(c).reshape(-1,3).tolist()
#     return remove_extra_points(array(c).round(5))

def twoCircleCrossTangent(c1,c2,cw=-1): # two circle cross tangent
    '''
    function to draw cross tangent between 2 circles
    refer to the file "example of various functions " for application examples
    '''
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

def twoCircleTangentPoints(c1,c2): #2 circle tangent point full (both the sides)
    '''
    function to draw tangent line joining 2 circles 'c1' and 'c2'.
    this works counter-clockwise
    This function draws tangent line on both the sides
    example:
    cir1=circle(10)
    cir2=circle(5,[15,6])
    sec=twoCircleTangentPoints(cir1,cir2)
    
    refer file "example of various functions" for application
    '''
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

    return [t3[0].tolist(),t4[0].tolist(),t2[0].tolist(),t1[0].tolist()]

def pts3(pl):
    '''
    this functionis required for 3d turtle movement and works with function 'corner_radius3d'
    example:
    a=pts3([[0,0,0],[10,7,5,5],[0,0,20,5],[0,10,0]]) => 
    [[0, 0, 0, 0], [10, 7, 5, 5], [10, 7, 25, 5], [10, 17, 25, 0]]
    '''
    a=pl
    b=[p if len(p)==4 else p+[0] for p in a]
    c=a_(b)[:,:3].cumsum(0)
    d=a_([a_(b)[:,3]])
    return l_(concatenate([c,d.T],axis=1))

def project_surface_on_to_another_surface_with_focal_line(surface_on_which_to_project, surface_to_project, focal_line, dist=100000, unidirection=0):
    '''
    project a surface on to another without loosing the original points
    focal_line is the line from which the rays are emitting and creates vectors for projection on to the surface
    '''
    s2,s3,v1=surface_on_which_to_project, surface_to_project, focal_line
    p0=a_(s3).reshape(-1,3)
    f=faces_surface(len(s2),len(s2[0]))
    v=a_(s2).reshape(-1,3)
    tri=v[f]
    p2,p3,p4=tri[:,0],tri[:,1],tri[:,2]
    a=[vcost_without_bounds(l1,p,1e5) for p in p0]
    b=cpo([a,p0])
    v1=[line_as_axis(p) for p in b]
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

def center_arc_2p(p1,p2,r,cw=-1):
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

def convert_4lines_enclosure2surface(l1,l2,l3,l4,n=5):
    '''
    create surface with 4 lines enclosed area
    l1 and l2 are opposite lines and l3 and l4 are opposite lines
    starting point of l1 and l3 should match
    end point of l3 and starting point of l2 should match
    end point of l1 and starting point of l4 should match
    end point of l4 and end point of l2 should match
    counterclockwise points should match as following
    l1[0]l3[0]->l1[-1]l4[0]->l4[-1]l2[-1]->l2[0]l3[-1]
    (l1[-1] means end point of l1)
    'n' is the number of segment lines between l1 and l2
    '''
    l5=slice_sol([l1,l2],n)
    l3=path2path1(cpo(l5)[0],l3)
    l4=path2path1(cpo(l5)[-1],l4)
    l5=cpo([l3]+cpo(l5)+[l4])
    s1=cpo(slice_sol([l3,l4],len(l5[0])))
    s1=[path2path1(s1[i],l5[i]) for i in range(len(l5))]
    return [equidistant_path(p,len(l1)-1) for p in s1]

