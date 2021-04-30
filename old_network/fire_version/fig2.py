import numpy as np
import matplotlib.pyplot as plt
from sys import exit



data = np.loadtxt("strain1.dat")
N = int( data.shape[0])
g = data[:N,0]
for i in range(g.shape[0]): g[i] = abs(g[i])
Hs = data[:N,1]
Hb = data[:N,2]

H = Hs + Hb

dg = g[1] - g[0]

sigma = np.gradient(H,g)
K = np.gradient(sigma,g)
dK = np.gradient(K,g)
dK2 = np.gradient(np.log(K),np.log(g))


#plt.scatter(g,H, color='red', label=r"$H$")
#plt.scatter(g,sigma, color='magenta', label=r"$\sigma= \frac{d H}{d \gamma}$")
plt.scatter(g,K, color='blue', label=r"$K = \frac{ d^2 H}{d\gamma^2}$")
#plt.scatter(g,dK, color='brown', label=r"$dK $")
#plt.scatter(g,dK2, color='orange', label=r"$\frac{ d log(K)}{d log(\gamma)} $")

data = np.loadtxt("strain2.dat")
N = int( data.shape[0])
g = data[:N,0]
for i in range(g.shape[0]): g[i] = abs(g[i])
Hs = data[:N,1]
Hb = data[:N,2]

H = Hs + Hb

dg = g[1] - g[0]

sigma = np.gradient(H,g)
K = np.gradient(sigma,g)
dK = np.gradient(K,g)
dK2 = np.gradient(np.log(K),np.log(g))


#plt.scatter(g,H, color='blue',  label=r"$H$")
#plt.scatter(g,sigma, color='magenta', facecolor="none", label=r"$\sigma= \frac{d H}{d \gamma}$")
plt.scatter(g,K, color='red',  label=r"$K = \frac{ d^2 H}{d\gamma^2}$")
#plt.scatter(g,dK, color='brown', label=r"$dK $")
#plt.scatter(g,dK2, color='orange', facecolor="none", label=r"$\frac{ d log(K)}{d log(\gamma)} $")


data = np.loadtxt("strain3.dat")
N = int( data.shape[0])
g = data[:N,0]
for i in range(g.shape[0]): g[i] = abs(g[i])
Hs = data[:N,1]
Hb = data[:N,2]

H = Hs + Hb

dg = g[1] - g[0]

sigma = np.gradient(H,g)
K = np.gradient(sigma,g)
dK = np.gradient(K,g)
dK2 = np.gradient(np.log(K),np.log(g))


#plt.scatter(g,H, color='green',  label=r"$H$")
#plt.scatter(g,sigma, color='magenta', facecolor="none", label=r"$\sigma= \frac{d H}{d \gamma}$")
plt.scatter(g,K, color='orange',  facecolor="none", label=r"$K = \frac{ d^2 H}{d\gamma^2}$")
#plt.scatter(g,dK, color='brown', label=r"$dK $")
#plt.scatter(g,dK2, color='orange', facecolor="none", label=r"$\frac{ d log(K)}{d log(\gamma)} $")



plt.yscale('log')
plt.xscale('log')
plt.xlim([ 0.5*dg, max(g)] )
plt.ylim([ 1e-17, 10] )
plt.legend();

plt.xlabel(r"$\gamma$", fontsize=20)

plt.show()
