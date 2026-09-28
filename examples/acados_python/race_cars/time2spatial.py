#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

# author: Daniel Kloeser

import numpy as np
from tracks.readDataFcn import getTrack



def transformProj2Orig(si,ni,alpha,v,filename='LMS_Track.txt'):
    [sref,xref,yref,psiref,_]=getTrack(filename=filename)
    tracklength=sref[-1]
    si=si%tracklength
    idxmindist=findClosestS(si,sref)
    idxmindist2=findSecondClosestS(si,sref,idxmindist)
    t=(si-sref[idxmindist])/(sref[idxmindist2]-sref[idxmindist])
    x0=(1-t)*xref[idxmindist]+t*xref[idxmindist2]
    y0=(1-t)*yref[idxmindist]+t*yref[idxmindist2]
    psi0=(1-t)*psiref[idxmindist]+t*psiref[idxmindist2]

    x=x0-ni*np.sin(psi0)
    y=y0+ni*np.cos(psi0)
    psi=psi0+alpha
    v=v
    return x,y,psi,v


def findClosestS(si,sref):
    # Get number of elements
    if(np.isscalar(si)):
        N=1
    else:
        N=np.array(si).shape[0]
    mindist=100000*np.ones(N)
    idxmindist=np.zeros(N)
    for i in range(sref.size):
        di=abs(si-sref[i])
        idxmindist = np.where(di < mindist,i, idxmindist)
        mindist = np.where(di < mindist, di, mindist)
    idxmindist = np.where(idxmindist==sref.size,1,idxmindist)
    idxmindist = np.where(idxmindist<1,sref.size-1,idxmindist)
    return idxmindist.astype(int)


def findSecondClosestS(si,sref,idxmindist):
    d1=abs(si-sref[idxmindist-1])               # distance to node before
    d2=abs(si-sref[(idxmindist+1)%sref.size])   # distance to node after
    idxmindist2 = np.where(d1>d2,idxmindist+1,idxmindist-1) # decide which node is closer
    idxmindist2 = np.where(idxmindist2==sref.size,0,idxmindist2)    # if chosen node is too large
    idxmindist2 = np.where(idxmindist2<0,sref.size-1,idxmindist2)   # if chosen node is too small

    return idxmindist2

def transformOrig2Proj(x,y,psi,v,filename='LMS_Track.txt'):
    [sref,xref,yref,psiref,_]=getTrack(filename=filename)
    idxmindist=findClosestPoint(x,y,xref,yref)
    idxmindist2=findClosestNeighbour(x,y,xref,yref,idxmindist)
    t=findProjection(x,y,xref,yref,sref,idxmindist,idxmindist2)
    s0=(1-t)*sref[idxmindist]+t*sref[idxmindist2]
    x0=(1-t)*xref[idxmindist]+t*xref[idxmindist2]
    y0=(1-t)*yref[idxmindist]+t*yref[idxmindist2]
    psi0=(1-t)*psiref[idxmindist]+t*psiref[idxmindist2]

    s=s0
    n=np.cos(psi0)*(y-y0)-np.sin(psi0)*(x-x0)
    alpha=psi-psi0
    v=v
    return s,n,alpha,v

def findProjection(x,y,xref,yref,sref,idxmindist,idxmindist2):
    vabs=abs(sref[idxmindist]-sref[idxmindist2])
    vl=np.empty(2)
    u=np.empty(2)
    vl[0]=xref[idxmindist2]-xref[idxmindist]
    vl[1]=yref[idxmindist2]-yref[idxmindist]
    u[0]=x-xref[idxmindist]
    u[1]=y-yref[idxmindist]
    t=(vl[0]*u[0]+vl[1]*u[1])/vabs/vabs
    return t

def findClosestPoint(x,y,xref,yref):
    mindist=1
    idxmindist=0
    for i in range(xref.size):
        dist=dist2D(x,xref[i],y,yref[i])
        if dist<mindist:
            mindist=dist
            idxmindist=i
    return idxmindist

def findClosestNeighbour(x,y,xref,yref,idxmindist):
    distBefore=dist2D(x,xref[idxmindist-1],y,yref[idxmindist-1])
    distAfter=dist2D(x,xref[idxmindist+1],y,yref[idxmindist+1])
    if(distBefore<distAfter):
        idxmindist2=idxmindist-1
    else:
        idxmindist2=idxmindist+1
    if(idxmindist2<0):
        idxmindist2=xref.size-1
    elif(idxmindist==xref.size):
        idxmindist2=0
    return idxmindist2


def dist2D(x1,x2,y1,y2):
    return np.sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))


