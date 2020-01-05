#!/usr/bin/python

import sys, time
import pickle
import numpy as np
from numpy.linalg import norm

def getCMD(fp):
    while(True):
        s=fp.readline()
        if len(s)==0:
            break
        s=s.strip()
        sl=s.split('#')
        if len(sl) > 0:
            s=sl[0].strip()
        if len(s)==0:
            continue
        sl=s.split()
        yield sl
    return

def parseCMD(cps):
    tmp=list(cps)
    cmd=tmp.pop(0)
    params=list(map(float,tmp))
    return cmd, params

gp=dict()
ap=dict()
def getParams():
    for s in getCMD(sys.stdin):
        param,params=parseCMD(s)
        if not param in gp:
            gp[param]=list()
        gp[param].append(params)
    xmi,xma,ymi,yma,zmi,zma=gp['field'][0]
    step,=gp['step'][0]
    nx=int((xma-xmi)/step)+1
    ny=int((yma-ymi)/step)+1
    nz=int((zma-zmi)/step)+1
    ap["mi"]=np.array([xmi,ymi,zmi])
    ap["ma"]=np.array([xma,yma,zma])
    ap["shape"]=np.array([nx,ny,nz], dtype=int)
    ap["kip"]=(ap["ma"]-ap["mi"])/(ap["shape"]-1)
#-------------------------------------------------
    plast=list()
    for pl in gp["rect"]:
        q,xmi,xma,ymi,yma,zmi,zma = pl
        nx=int((xma-xmi)/step)+1
        ny=int((yma-ymi)/step)+1
        nz=int((zma-zmi)/step)+1
        plast.append(list())
        plast[-1].append(q/(nx*ny*nz))
        plast[-1].append(np.zeros((0,3)))
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    a1=plast[-1][-1]
                    a2=np.array([xmi,ymi,zmi])+np.array([ix,iy,iz])*step
                    print(a1.shape,a2.shape)
                    plast[-1][-1]=np.vstack([a1,a2])
            print(f"getParams: ix={ix} .. {nx}")
    ap['plast']=plast
    return

def emfr1(vec, p):	#array[n:3], array[1:3]
    r=vec-p
#    epl=r/np.linalg.norm(r)**3
    epl=r/norm(r)**3
    exyz=np.sum(epl,axis=0)
    return exyz

def main():
    print('-'*40)
    getParams()
    for name,val in gp.items():
        print(name,':')
        print(val)

    print('-'*40)

    nx,ny,nz=map(int,ap['shape'])
    field=np.zeros(ap['shape'])
    print(f"main: field.shape={field.shape}")
#    apkip=ap['kip']
#    apmi=ap['mi']
    for ix in range(nx):
        print(f"{time.ctime(time.time())} main: ix={ix} .. {nx}")
        for iy in range(ny):
            for iz in range(nz):
                ixyz=np.array([ix,iy,iz])
                xyz=ap['kip']*ixyz+ap['mi']
#                xyz=apkip*ixyz+apmi
                e=np.zeros((1,3))
                for plast in ap['plast']:
                    q,pla=plast
                    e+=emfr1(pla,xyz)*q
#                field[ix,iy,iz]=np.linalg.norm(e)
                field[ix,iy,iz]=norm(e)
    print('-'*40)

    with open("gp.dump", 'wb') as fp:
        pickle.dump(gp, fp, pickle.HIGHEST_PROTOCOL)
    with open("ap.dump", 'wb') as fp:
        pickle.dump(ap, fp, pickle.HIGHEST_PROTOCOL)
    with open("field.dump", 'wb') as fp:
        pickle.dump(field, fp, pickle.HIGHEST_PROTOCOL)

    print('-'*40)

if __name__ == "__main__":
    main()
