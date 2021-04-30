import numpy as np
import matplotlib.pyplot as plt
from sys import exit

nu = 1.38
f = 0.75

gc = 0.28085
ni = 2
nf = 2
w = [ 50, 75, 100, 150, 200]
c = [ 'red' , 'blue', 'orange', 'purple', 'black']
Nw = len(w)
H = []
sigma = []
K = []


data = np.loadtxt("strainW50.dat")
g = data[:,0]
dg = g[1] - g[0]

for wi in range(Nw):
	data = np.loadtxt("strainW"+str(w[wi])+".dat")
	H.append(data[:,1] + data[:,2] )
	sigma.append(  np.gradient(H[wi],g) )
	K.append(  np.gradient(sigma[wi],g) )


fig, axes = plt.subplots( nrows=1, ncols=2, figsize=(10,5) )

ax = axes[0]
for wi in range(Nw):
	ax.scatter(g[ni:-nf],K[wi][ni:-nf], label=w[wi], color=c[wi])

ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim([ 1.e-3, max(g)] )
ax.set_ylim([ 1e-8, 10] )
ax.legend();

ax.set_ylabel(r"$K$", fontsize=20)
ax.set_xlabel(r"$\gamma$", fontsize=20)

ax = axes[1]
for wi in range(Nw):
	ax.scatter(abs(g-gc)*(w[wi]**(1/nu) ),K[wi]*(w[wi]**(f/nu)), color= c[wi])

ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim([1e-2, 1e3] )
ax.set_ylim([ 1e-6, 1e2] )

ax.set_ylabel(r"$K W^{f/\nu}$", fontsize=20)
ax.set_xlabel(r"$|\Delta \gamma| W^{1/\nu}$", fontsize=20)

plt.tight_layout()
plt.show()
