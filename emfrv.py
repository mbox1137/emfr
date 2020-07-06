#!/usr/bin/python3.7

#Plotting a 3d surface from a list of tuples in matplotlib
#https://stackoverflow.com/questions/21161884/plotting-a-3d-surface-from-a-list-of-tuples-in-matplotlib

import pdb
#pdb.set_trace()
import os, sys
import pickle
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt


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

if len(sys.argv)!=4:
    print(f"{sys.argv[0]} data/emfr.dat 0.45 0.02")
    sys.exit()
datfile=sys.argv[1]
prefix, fn = os.path.split(datfile)
with open(prefix+"/gp.dump", 'rb') as fp:
    gp=pickle.load(fp)
with open(prefix+"/ap.dump", 'rb') as fp:
    ap=pickle.load(fp)
with open(prefix+"/field.dump", 'rb') as fp:
    field=pickle.load(fp)
porog=float(sys.argv[2])
delta=float(sys.argv[3])

def main():
    print('-'*40)
    print("field=",field);

    fmi=np.min(field)
    fma=np.max(field)
    print(f"fmi={fmi} fma={fma}");
#    sys.exit()

    df=fma-fmi
#    ixs=np.argwhere(np.abs(field-fma*porog)<0.1*fma)
    ixs=np.argwhere(np.abs((field-fmi)/df-porog)<delta)
    xs=ap['kip']*ixs+ap['mi']
    print("field:", field.size)
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

    plt.title(f"step={gp['step'][0][0]}, porog={porog}, delta={delta}")
    plt.savefig('emfrv.png', dpi=300)
    plt.show()
    plt.close()

    print('-'*40)

if __name__ == "__main__":
    main()
