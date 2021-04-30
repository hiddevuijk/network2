import numpy as np
import matplotlib.pyplot as plt
from sys import exit


topology = open('topology.txt')
lines = topology.readlines()
topology.close()

Nv = int( lines[0].strip() )
Ne = int( lines[1].strip() )
Nb = int( lines[2].strip() )

xList = np.zeros( Nv )
yList = np.zeros( Nv )
edges = np.asarray( [ [-1,-1] ]*Ne )
bends = np.asarray( [ [-1,-1,-1] ]*Nb )

for vi in range(Nv):
    xList[vi] = lines[3+vi].strip().split('\t')[1]
    yList[vi] = lines[3+vi].strip().split('\t')[2]


for ei in range(Ne):
    edge = lines[3+Nv+ei].strip().split('\t')
    edges[ei][0] = edge[0];
    edges[ei][1] = edge[1];

for bi in range(Nb):
    bend = lines[3+Nv+Ne].strip().split('\t')
    bends[bi][0] = bend[0]
    bends[bi][1] = bend[1]
    bends[bi][2] = bend[2]

plt.scatter(xList, yList,color='black', s=100) 

for ei in range(Ne):
    xfrom = xList[ edges[ei][0] ]
    yfrom = yList[ edges[ei][0] ]
    xto = xList[ edges[ei][1] ]
    yto = yList[ edges[ei][1] ]

    plt.plot( [ xfrom, xto], [yfrom, yto], linewidth=3, color='blue' )


plt.show()
