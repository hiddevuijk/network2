import numpy as np
import matplotlib.pyplot as plt
from sys import exit

topology = open('topology.dat')
lines = topology.readlines()
topology.close()

Lx = 12
Ly = Lx*np.sqrt(3/4.)
gamma = 0
Nv = int( lines[0].strip() )
Ne = int( lines[1].strip() )
Nb = int( lines[2].strip() )

# read positions from r0.dat
rList = open('r0.dat').readlines()
xList = np.zeros(Nv)
yList = np.zeros(Nv)

for vi in range(Nv):
	xy = rList[vi].strip().split('\t')
	xList[vi] = xy[0]
	yList[vi] = xy[1]

x0List = np.zeros( Nv )
y0List = np.zeros( Nv )
edges = np.asarray( [ [-1]*4 ]*Ne )
bends = np.asarray( [ [-1]*7 ]*Nb )

# read positions from topology.dat
for vi in range(Nv):
	x0List[vi] = lines[3+vi].strip().split('\t')[1]
	y0List[vi] = lines[3+vi].strip().split('\t')[2]

for ei in range(Ne):
	edges[ei] = lines[3+Nv+ei].strip().split('\t')
	

for bi in range(Nb):
	bends[bi] = lines[3+Nv+Ne+bi].strip().split('\t');

for ei in range(Ne):
	try:	
		xfrom = xList[ edges[ei][0] ]
		yfrom = yList[ edges[ei][0] ]
		xto = xList[ edges[ei][1] ]
		yto = yList[ edges[ei][1] ]

		xb = edges[ei][2];
		yb = edges[ei][3];

		if( xb == 0 and yb == 0): 
			plt.plot( [ xfrom, xto], [yfrom, yto], linewidth=3, color='grey' )
		else: 
			#continue
			plt.plot( [ xfrom, xto + xb*Lx], [yfrom, yto + yb*Ly], linewidth=3, color='lightgrey' )
			plt.plot( [ xto, xfrom-xb*Lx ], [yto , yfrom - yb*Ly], linewidth=3, color='lightgrey' )
	except:
		print(edges[ei][0], edges[ei][1] )		

#plt.scatter(xList, yList,color='black', s=30, zorder=10)

d = 0.0075
#for i in range(Nv):
#    plt.annotate( i, ( xList[i]-d, yList[i]-d ), color='red', zorder = 11 )

for i in range(Nv):
    plt.scatter(xList[i], yList[i], color='red', zorder=11, s=5)
plt.axis("off")
plt.gca().set_aspect('equal')
plt.tight_layout()

plt.show()
