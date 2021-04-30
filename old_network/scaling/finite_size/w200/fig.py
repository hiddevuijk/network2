import numpy as np
import matplotlib.pyplot as plt
from sys import exit



data = np.loadtxt("strainW200.dat")
N = int( data.shape[0])
g = data[:N,0]
Hs = data[:N,1]
Hb = data[:N,2]
sxx = data[:N,3]
sxy = data[:N,4]
syx = data[:N,5]
syy = data[:N,6]

H = Hs + Hb

dg = g[1] - g[0]

sigma = np.gradient(H,g)
K = np.gradient(sigma,g)


plt.scatter(g,H, color='red', label=r"$H$")
plt.scatter(g,sigma, color='magenta', label=r"$\sigma= \frac{d H}{d \gamma}$")
plt.scatter(g,K, color='blue', label=r"$K = \frac{ d^2 H}{d\gamma^2}$")

plt.yscale('log')
plt.xscale('log')
plt.xlim([ 0.5*dg, max(g)] )
plt.ylim([ 1e-17, 10] )
plt.legend();

plt.xlabel(r"$\gamma$", fontsize=20)

plt.show()
