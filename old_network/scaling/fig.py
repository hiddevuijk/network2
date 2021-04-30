import numpy as np
import matplotlib.pyplot as plt
from sys import exit

n = ""
epsilon = 1e-5
ni = 3
nf = 3


f = 0.75
phi = 2.26
gc = 0.268798
gc2 = 0.12784
gc3 = 0.258774
gc = gc

k = [0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
H = []
K = []

c = ["black",  "red", "green", "blue", "orange", "magenta", "brown"]
data = np.loadtxt("strain0"+n+".dat")
g = data[:,0]
dg = g[1] - g[0]
Hs0 = data[:,1]
Hb0 = data[:,2]
H.append(Hs0+Hb0)


data = np.loadtxt("strain1"+n+".dat")
Hs1 = data[:,1]
Hb1 = data[:,2]
H.append(Hs1+Hb1)


data = np.loadtxt("strain2"+n+".dat")
Hs2 = data[:,1]
Hb2 = data[:,2]
H.append(Hs2+Hb2)

data = np.loadtxt("strain3"+n+".dat")
Hs3 = data[:,1]
Hb3 = data[:,2]
H.append(Hs3+Hb3)

data = np.loadtxt("strain4"+n+".dat")
Hs4 = data[:,1]
Hb4 = data[:,2]
H.append(Hs4+Hb4)

data = np.loadtxt("strain5"+n+".dat")
Hs5 = data[:,1]
Hb5 = data[:,2]
H.append(Hs5+Hb5)

data = np.loadtxt("strain6"+n+".dat")
Hs6 = data[:,1]
Hb6 = data[:,2]
H.append(Hs6+Hb6)




for i in range(7):
	sigma = np.gradient(H[i],g)
	K.append( np.gradient(sigma,g) )


fig, axes = plt.subplots( nrows=1, ncols=3, figsize=(12,4) )

ax = axes[0]
for i in range(1,7):
	ax.scatter(g,K[i],color=c[i],label= k[i])

ax.scatter(g,K[0],label= k[0], color="black")
ax.axvline(gc)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim([ 0.5*dg, max(g)] )
ax.set_ylim([ 1e-7, 1] )
ax.legend();
ax.set_ylabel(r"$ K $", fontsize=10)
ax.set_xlabel(r"$\gamma$", fontsize=10)


ax = axes[1]

for i in range(1,7):
	for gi in range(ni, g.shape[0] - nf):
		delta_g = abs( g[gi] - gc )
		if delta_g < epsilon : continue
		ax.scatter( k[i] / ( delta_g**phi) , K[i][gi]/( delta_g**f) ,color=c[i],label= k[i])
ax.set_xlabel(r"$ \kappa | \Delta \gamma|^{-\phi} $", fontsize=10)
ax.set_ylabel(r"$ K | \Delta \gamma|^{-f}$", fontsize=10)
ax.set_yscale('log')
ax.set_xscale('log')

ax.set_xlim([ 1e-9, 1e4] )
ax.set_ylim([ 1e-6, 10] )

ax.set_title("f = "+str(f) + "   phi = "+str(phi) )

ax = axes[2]
epsilon=0.
for i in range(1,7):
	for gi in range(ni, g.shape[0] ) :
		if g[gi] > gc+epsilon: 
			ax.scatter(g[gi]-gc,K[i][gi],color=c[i])

i = 0
for gi in range(ni, g.shape[0] ) :
	if g[gi] > gc+epsilon: 
		ax.scatter(g[gi]-gc,K[i][gi],color=c[i])

x = np.logspace(0.001, 1,100) - 1
N = 0.15
ax.plot( x, N*(x**f) )

ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim([ 0.005, max(g)] )
ax.set_ylim([ 0.01, 1] )

ax.set_ylabel(r"$ K $", fontsize=10)
ax.set_xlabel(r"$\Delta \gamma$", fontsize=10)
ax.set_title("f="+str(f))


fig.tight_layout()

plt.show()
