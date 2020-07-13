#!/usr/bin/python3.7 -u

#import cProfile
#import emfr
 
import os, sys, time, builtins
import pickle
import numpy as np
from numpy.linalg import norm
import numexpr as ne
import multiprocessing as mp
import pdb
#pdb.set_trace()

def dprint(*args, **kwargs):
    if False:
        print(*args, **kwargs)

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
def getParams(fp):
    for s in getCMD(fp):
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
        dprint(f"pl={pl}")
        q,xmi,xma,ymi,yma,zmi,zma = pl
        nx=int((xma-xmi)/step)+1
        ny=int((yma-ymi)/step)+1
        nz=int((zma-zmi)/step)+1
        dprint(f"getParams(plast): nx={nx} ny={ny} nz={nz}")
        plast.append(list())
        plast[-1].append(q/(nx*ny*nz))
#        plast[-1].append(None)	#np.zeros((0,3)))
        xyzmi=np.array([xmi,ymi,zmi])
        if True:
            lst=[[ix,iy,iz] for ix in range(nx) for iy in range(ny) for iz in range(nz)]
            plast[-1].append(np.array(lst)*step+xyzmi)
#            print(f"plast={plast}")
        else:
            for ix in range(nx):
                for iy in range(ny):
                    for iz in range(nz):
                        a1=np.zeros((0,3))
                        a2=np.array([xmi,ymi,zmi])+np.array([ix,iy,iz])*step
                        print(a1.shape,a2.shape)
                        plast[-1].append(np.vstack([a1,a2]))
                dprint(f"getParams: ix={ix} .. {nx}")
    dprint(f"plast={plast}")
    ap['plast']=plast
    return

def emfr1(vec, p):	#array[n:3], array[1:3]
    dprint(f"vec={vec}")
    dprint(f"p={p}")
    r=vec-p
    dprint(f"r.shape={r.shape}")
    dprint(f"r={r}")
#	hypot(x, y) Return the Euclidean norm, sqrt(x*x + y*y)
    rnorm=norm(r,axis=1).reshape((r.shape[0],1))
    dprint(f"rnorm.shape={rnorm.shape}")
    dprint(f"rnorm={rnorm}")
    e=r/(rnorm**3)
    dprint(f"e.shape={e.shape}")
    dprint(f"e={e}")
    e0=np.sum(e,axis=0)
    dprint(f"e0={e0}")
#    sys.exit()
    return e0
#    tmp=r/norm(r)**3
#    rn3=norm(r)**3
#    tmp=ne.evaluate('r/rn3')
#    tmp=r/rn3
#    return np.sum(tmp,axis=0)

kk=0
def emfr2(p):	#array[1:3]
    global kk
    e=np.zeros((1,3))
    for q,pla in ap['plast']:
        e+=emfr1(pla,p)*q
    kk+=1
    if kk%1000 == 0:
        print(f"{time.ctime(time.time())} {kk}")
    return norm(e)

kk=0
def sun(p):	#array[1:3]
    global kk
    p0=ap['kip']*ap['shape']*0.5+ap['mi']
    dprint("p0=",p0)
    r=p-p0
    r0=norm(r)
    dprint(r0)
    dr=norm(ap['kip'])
    if r0<dr:
        r0=dr
        r=np.array([dr]*3)
    e=r/(r0**3)
    if kk%1000 == 0:
        dprint(f"{time.ctime(time.time())} {kk}")
    return norm(e)

if len(sys.argv)!=3:
    print(f"{sys.argv[0]} data/emfr.dat 5")
    sys.exit()
else:
    datfile=sys.argv[1]
    prefix, fn = os.path.split(datfile)
    with open(datfile) as fp:
        getParams(fp)
    nmp=int(sys.argv[2])

def main_():
    print('-'*40)
    print(f"Treads={nmp}")
    print(f"prefix={prefix}")
    for name,val in gp.items():
        print(f"{name}: {val}")
    print('-'*40)

    nx,ny,nz=map(int,ap['shape'])
    field=np.zeros(ap['shape'])
    print(f"main: field.shape={field.shape}")
    if True:
        print(f"ap['kip']={ap['kip']} ap['mi']={ap['mi']}")
        ixyz=[[ix,iy,iz] for ix in range(nx) for iy in range(ny) for iz in range(nz)]
        dprint(f"ixyz={ixyz}")
        xyz=ap['kip']*ixyz+ap['mi']
        print(f"xyz.shape={xyz.shape}")
        dprint(f"xyz={xyz}")
#        exyz=np.array(list(map(sun,xyz)))
        exyz=np.array(list(map(emfr2,xyz)))
        print(f"exyz.shape={exyz.shape}")
        dprint(f"exyz={exyz}")
        for k in range(len(exyz)):
            ix,iy,iz=ixyz[k]
            field[ix,iy,iz]=exyz[k]
#            field[ixyz[k]]=exyz[k]
        print(f"xyz.shape={xyz.shape}")
        print(f"field.shape=",field.shape)
    else:
        for ix in range(nx):
            print(f" main: ix={ix} .. {nx}")
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

    with open(prefix+"/gp.dump", 'wb') as fp:
        pickle.dump(gp, fp, pickle.HIGHEST_PROTOCOL)
    with open(prefix+"/ap.dump", 'wb') as fp:
        pickle.dump(ap, fp, pickle.HIGHEST_PROTOCOL)
    with open(prefix+"/field.dump", 'wb') as fp:
        pickle.dump(field, fp, pickle.HIGHEST_PROTOCOL)

    print('-'*40)

def inbox(p):
    x,y,z = p
    for pl in gp["rect"]:
        q,xmi,xma,ymi,yma,zmi,zma = pl
        if((x>=xmi and x<=xma) and
           (y>=ymi and y<=yma) and 
           (z>=zmi and z<=zma)):
            return True           
    return False

def emfr3(p):	#array[1:3]
    e=np.zeros((1,3))
    for q,pla in ap['plast']:
        e+=emfr1(pla,p)*q
    return norm(e)

def f(num,lin,lout,qin,qout):
    while not qin.empty():
        lin.acquire()
        try:
            xyz_,ixyz_=qin.get()
        finally:
            lin.release()
        e=emfr3(xyz_)
        lout.acquire()
        try:
            qout.put((ixyz_,e))
        finally:
            lout.release()
    return

if __name__ == "__main__":
#    cProfile.run('emfr.main()')
    print('-'*40)
    print(f"Treads={nmp}")
    print(f"prefix={prefix}")
    for name,val in gp.items():
        print(f"{name}: {val}")
    print('-'*40)
#    sys.exit()
    dumpfromfile=True
    try:
        datfile=sys.argv[1]
        prefix, datfn = os.path.split(datfile)
        with open(prefix+"/field.dump", 'rb') as fp:
            field=pickle.load(fp)
    except:
        dumpfromfile=False

    if not dumpfromfile or (ap['shape']!=field.shape).any():
        field=np.zeros(ap['shape'])
    nx,ny,nz=map(int,ap['shape'])
    ixyz=[[ix,iy,iz] for ix in range(nx) for iy in range(ny) for iz in range(nz)]
    xyz=ap['kip']*ixyz+ap['mi']
        
    for k in range(len(ixyz)):
        ix, iy, iz = ixyz[k]
        if field[ix, iy, iz]==0.0 and inbox(xyz[k]):
            field[ix, iy, iz]=np.nan

    print(f"xyz.shape={xyz.shape}")
    print(f"field.shape=",field.shape)

    qin = mp.Queue()
    qout = mp.Queue()
    lin = mp.Lock()
    lout = mp.Lock()
    procs=list()
    for k in range(nmp):
        proc=mp.Process(target=f, args=(k,lin,lout,qin,qout))
        proc.daemon=True
        procs.append(proc)
    npo=0
    for k in range(len(ixyz)):
        ix, iy, iz = ixyz[k]
        if field[ix, iy, iz]!=0.0:
            continue
        qin.put((xyz[k],ixyz[k]))
        npo+=1
    print(f"Start ({len(procs)}: {npo})")
    for proc in procs:
        proc.start()
    t0=time.time()
    work=True
    try:
        while work:
            if qin.qsize()<2*nmp:
                time.sleep(1)
            pflag=any(proc.is_alive() for proc in procs)
            qflag=qout.empty()
            if not pflag and qflag:
                break
            (ix,iy,iz),e=qout.get()
            field[ix,iy,iz]=e
            kk+=1
            if kk%1000 == 0:
                t=time.time()
                tx=(t-t0)/kk*(npo-kk)+t
                print(f"{time.ctime(tx)} {kk}")
    except KeyboardInterrupt:
        work=False
    time.sleep(1)
    print("qq")

    for proc in procs:
        proc.terminate()

    for proc in procs:
        proc.join()

    with open(prefix+"/gp.dump", 'wb') as fp:
        pickle.dump(gp, fp, pickle.HIGHEST_PROTOCOL)
    with open(prefix+"/ap.dump", 'wb') as fp:
        pickle.dump(ap, fp, pickle.HIGHEST_PROTOCOL)
    with open(prefix+"/field.dump", 'wb') as fp:
        pickle.dump(field, fp, pickle.HIGHEST_PROTOCOL)

    print("qqq")
