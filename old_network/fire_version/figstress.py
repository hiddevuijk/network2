import numpy as np
import matplotlib.pyplot as plt
from sys import exit



data = np.loadtxt("strain.dat")
N = int( data.shape[0])
#print(N)
g = data[:,0]
for i in range(N): g[i] = abs(g[i])
Hs = data[:,1]
Hb = data[:,2]
sxx = data[:,3]
sxy = data[:,4]
syx = data[:,5]
syy = data[:,6]

for i in range(N): sxx[i] = abs(sxx[i])
for i in range(N): sxy[i] = abs(sxy[i])
for i in range(N): syx[i] = abs(syx[i])
for i in range(N): syy[i] = abs(syy[i])

H = Hs + Hb

dg = g[1] - g[0]
#print(dg)

sigma = np.gradient(H,g)
K = np.gradient(sigma,g)


plt.scatter(g,H, color='red', label=r"$H$")
plt.scatter(g,sigma, color='magenta', label=r"$\sigma= \frac{d H}{d \gamma}$")
plt.scatter(g,K, color='blue', label=r"$K = \frac{ d^2 H}{d\gamma^2}$")
plt.scatter(g,sxx, label= "sxx")
plt.scatter(g,sxy, label= "sxy")
plt.scatter(g,syx, label= "syx")
plt.scatter(g,syy, label= "syy")

plt.yscale('log')
plt.xscale('log')
plt.xlim([ 0.5*dg, max(g)] )
plt.ylim([ 1e-17, 10] )
plt.legend();

plt.xlabel(r"$\epsilon$", fontsize=20)


plt.tight_layout()
plt.show()
