import numpy as np
import matplotlib.pyplot as plt
from sys import exit

data = np.loadtxt("strain6_2.dat")
g = data[:,0]
dg = g[1] - g[0]
Hs = data[:,1]
Hb = data[:,2]



sigma = np.gradient(Hs+Hb,g)
K = np.gradient(sigma,g) 


plt.scatter(g,K)

plt.yscale('log')
plt.xscale('log')
plt.xlim([ 0.5*dg, max(g)] )
plt.ylim([ 1e-6, 10] )
plt.legend();
plt.xlabel(r"$\gamma$", fontsize=20)


plt.tight_layout()
plt.show()
