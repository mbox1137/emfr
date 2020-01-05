#!/usr/bin/python

#Plotting a 3d surface from a list of tuples in matplotlib
#https://stackoverflow.com/questions/21161884/plotting-a-3d-surface-from-a-list-of-tuples-in-matplotlib

import sys
import pickle
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt

try:
    porog=float(sys.argv[1])
except IndexError:
    porog=0.6
except ValueError:
    porog=0.6

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

def main():
    print('-'*40)
    with open("gp.dump", 'rb') as fp:
        gp=pickle.load(fp)
    with open("ap.dump", 'rb') as fp:
        ap=pickle.load(fp)
    with open("field.dump", 'rb') as fp:
        field=pickle.load(fp)

    fma=np.max(field)
    ixs=np.argwhere(np.abs(field-fma*porog)<0.01*fma)
    xs=ap['kip']*ixs+ap['mi']
    print("field:", field.size)
    print("xs:", xs.shape)

    fi=np.arctan2(xs[:,1],xs[:,0])
    teta=np.arctan2(xs[:,2],np.linalg.norm(xs[:,0:2],axis=1))

    ng=20
    fld=[[list() for j in range(ng)] for i in range(ng)]

    for k in range(ixs.shape[0]):
        ifi=((fi[k]+np.pi)/(2*np.pi)*(ng-1)).astype(int)
        iteta=((teta[k]+np.pi/2)/np.pi*(ng-1)).astype(int)
        fld[iteta][ifi].append(xs[k,:])
    for iteta in range(ng):
        for ifi in range(ng):
            if len(fld[iteta][ifi]) == 0:
                fld[iteta][ifi]=None
                continue
            fld[iteta][ifi]=np.average(fld[iteta][ifi], axis=0)
    segments=list()
    for iteta in range(ng):	#range(5,ng)
#        print(fld[iteta])
        ztf=[type(fld[iteta][i]) == type(None) for i in range(len(fld[iteta]))]
#        print(ztf)
        z=list(np.where(ztf)[0])
#        print(f"z={z}")
        z1=[0]+z
#        print(f"z1={z1}")
        z2=z+[len(ztf)]
#        print(f"z2={z2}")
#        print(list(zip(z1, z2)))
        le=[(i+1,j) for i, j in zip(z1, z2) if i+1<j]	#Паузы
        print(f"le={le}")
        for k1,k2 in le:
            if k2-k1>ng//2:
                segments.append(fld[iteta][k1:k2])
    print("segments: ", len(segments))
#fld[iteta][ifi]=np.average(fld[iteta][ifi], axis=0)

    for ifi in range(ng):
        pass

#    ifi=((fi+np.pi)/(2*np.pi)*(ng-1)).astype(int)
#    iteta=((teta+np.pi/2)/np.pi*(ng-1)).astype(int)
#    nx,ny,nz=map(int,ap['shape'])

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(20,35)	#degree

    for rect in gp['rect']:
        rt=list(rect)
        q=rt.pop(0)
        color='r' if q>0.0 else 'b'
        for edge in getRectEdges(rt):
            ax.plot(*edge,color)

    if False:
        for segment in segments:
            tmp=np.array(segment)
            ax.plot(tmp[:,0],tmp[:,1],tmp[:,2],'k')

#    ax.plot_wireframe(ixs[:,0], ixs[:,1], ixs[:,2])
#    surf = ax.plot_trisurf(X, Y, Z, linewidth=0, antialiased=False)
    X, Y, Z = xs[:,0], xs[:,1], xs[:,2]
    ax.scatter(X,Y,Z,'.')
#    ax.grid(True)

    plt.title(f"step={gp['step'][0][0]} porog={porog}*max")
    plt.savefig('emfrv.png', dpi=300)
    plt.show()
    plt.close()

    print('-'*40)

if __name__ == "__main__":
    main()
