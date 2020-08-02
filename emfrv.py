#!/usr/bin/python3.7

#Plotting a 3d surface from a list of tuples in matplotlib
#https://stackoverflow.com/questions/21161884/plotting-a-3d-surface-from-a-list-of-tuples-in-matplotlib

#import tkinter
import pdb
#pdb.set_trace()
import os, sys
import pickle
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt

def dprint(*args, **kwargs):
    if False:
        print(*args, **kwargs)

#aa=np.array(list(range(60))).reshape((4,5,3))
#print(aa.tolist())

def getRectEdges(rect):
    xmi,xma,ymi,yma,zmi,zma=rect
    eds={0,1,2}
    for ed in eds:
        oed=list(eds-{ed})
        res=[None, None, None]
        res[ed]=[rect[2*ed],rect[2*ed+1]]
        for k1 in range(2):
            for k2 in range(2):
                res[oed[0]]=[rect[oed[0]*2+k1],rect[oed[0]*2+k1]]
                res[oed[1]]=[rect[oed[1]*2+k2],rect[oed[1]*2+k2]]
                yield res	#[[xmi,xma],[ymi,ymi],[zmi,zmi]]
    return

if len(sys.argv) not in [3,4]:
    print(f"{sys.argv[0]} data/emfr.dat 0.002 [0.8]")
    sys.exit()
datfile=sys.argv[1]
prefix, datfn = os.path.split(datfile)
with open(prefix+"/gp.dump", 'rb') as fp:
    gp=pickle.load(fp)
with open(prefix+"/ap.dump", 'rb') as fp:
    ap=pickle.load(fp)
with open(prefix+"/field.dump", 'rb') as fp:
    field=pickle.load(fp)
delta=float(sys.argv[2])
try:
    porog=float(sys.argv[3])
except IndexError:
    porog=None
except ValueError:
    porog=None

"""
----------------------------------------
Treads=5
prefix=data/120x120x120
step: [[0.02]]
field: [[-1.2, 1.2, -1.2, 1.2, -1.2, 1.2]]
rect: [[1.0, -1.0, 1.0, -1.0, 1.0, 1.0, 1.1], [-1.0, -1.0, 1.0, -1.0, 1.0, -1.1, -1.0]]
----------------------------------------
"""
def main():
    global porog
    print('-'*40)
    dprint("field=",field);

    print(f"NaNs: {len(np.isnan(field))}")
    q=abs(gp["rect"][1][0]-gp["rect"][0][0])
    field[np.isnan(field)]=np.finfo(float).max
    if not porog:
        nh=20
        fmax=5
        bins=np.append(np.linspace(0,fmax,nh), np.finfo(float).max)
        dprint(bins)
        counts, edges = np.histogram(field, bins)
        dprint(counts)
        dprint(edges)
        maxind=np.where(counts == np.amax(counts))[0][0]
        dprint(maxind)
        porog=(edges[maxind]+edges[maxind+1])/2
        dprint(porog)

    ixs=np.argwhere(np.abs(field-porog)<delta)
    xs=ap['kip']*ixs+ap['mi']
    print("field:", field.shape)
    print("xs:", xs.shape)
    """
    ng=11
    fldl=[[list() for i in range(ng)] for j in range(ng)]
    print("fldl:", len(fldl), len(fldl[0]))
    zmi=np.min(xs[:,2])
    zma=np.max(xs[:,2])
    for k in range(xs.shape[0]):
        ix,iy,iz=ixs[k]
        x,y,z=xs[k]
        fi=np.arctan2(y,x)
        ifi=int((fi+np.pi)/(2*np.pi)*(ng-1))
        iz=int((z-zmi)/(zma-zmi)*(ng-1))
        fldl[iz][ifi].append(xs[k])
    flda=[[np.array([np.nan, np.nan, np.nan]) for i in range(ng)] for j in range(ng)]
    for iz in range(ng):
        for ifi in range(ng):
            if len(fldl[iz][ifi]) > 0:
                flda[iz][ifi]=np.mean(np.array(fldl[iz][ifi]), axis=0).tolist()
    fld=np.array(flda)
    print(f"fld.shape={fld.shape}")
    """

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(20,35)	#degree
    ax.set_xlim3d(gp['field'][0][0:2])
    ax.set_ylim3d(gp['field'][0][2:4])
    ax.set_zlim3d(gp['field'][0][4:6])

    for rect in gp['rect']:
        rt=list(rect)
        q=rt.pop(0)
        color='r' if q>0.0 else 'b'
        for edge in getRectEdges(rt):
            ax.plot(*edge,color)
    if True:
        X, Y, Z = xs[:,0], xs[:,1], xs[:,2]
        ax.scatter(X,Y,Z,marker=".")
    else:
        color='g'
        for row in range(ng):
            X, Y, Z = fld[:,row,0], fld[:,row,1], fld[:,row,2]
            ax.plot(X,Y,Z,color)
        for col in range(ng):
            X, Y, Z = fld[col,:,0], fld[col,:,1], fld[col,:,2]
            ax.plot(X,Y,Z,color)

    plt.title(f"step={gp['step'][0][0]}, delta={delta}, porog={round(porog,3)}")
    dfname,dfext=datfn.split('.')
    plt.savefig(os.path.join(prefix,dfname+'.png'), dpi=300)
    plt.show()
    plt.close()
    print('-'*40)

if __name__ == "__main__":
    main()
