"direct translation to x=-1"
"domain : dom9_h4.geo; numbers of elements: (32*32)*(12*(12+12))"
"So For each of 4 process workers, (32*32)*(12*(12+12))/4"



from firedrake import *
import numpy as np
def trans(coords,phi,phi1,ddx,trans_gap):
    xi=np.where(abs(coords.dat.data[:,1]-coords.dat.data[0,1])<10**-2)[0]
    nn0=int(trans_gap/ddx)
    "numbers of points on longditudinal :(12+12)*12*4=288"
    sign1=int((coords.dat.data[xi[280],0]-coords.dat.data[xi[0],0])\
              /abs(coords.dat.data[xi[280],0]-coords.dat.data[xi[0],0]))
    nn1=len(xi)
    "for boundary points [0~2n0-1,... ,n1-2n0~n1-1]"    
    a1=int((nn1-1)*.5+sign1*.5*(nn1-2*nn0+1))
    # print(.5*((sign1+1)*(nn1-1)),(sign1+1),(nn1-1))
    y2=np.where(abs(coords.dat.data[:,0]-coords.dat.data[xi[int(.5*(sign1+1)*(nn1-1))],0])<=10**-2)
    for k in range(0,nn0):        
        y1=np.where(abs(coords.dat.data[:,0]-coords.dat.data[xi[a1+k*sign1],0])<=10**-2)
        phi1.dat.data[y1]=phi.dat.data[y2]
    
    for k in range(0,nn1-nn0):
        y2=np.where(abs(coords.dat.data[:,0]-coords.dat.data[xi[k+int(.5*nn0*(sign1+1))],0])<=10**-2)
        y1=np.where(abs(coords.dat.data[:,0]-coords.dat.data[xi[k-int(.5*nn0*(sign1-1))],0])<=10**-2)
        phi1.dat.data[y1]=phi.dat.data[y2]
    return(phi1)
