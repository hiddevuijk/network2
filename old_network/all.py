import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
from sys import exit

cmap = matplotlib.cm.get_cmap("viridis")

kappas = [2,3,4,5,6,0]
for ki, k in enumerate(kappas):
	
	data = np.loadtxt("strain_m{}.dat".format(k) )
	N = int( data.shape[0])
	g = data[:N,0]
	e = data[:N,1]

	dedg = np.gradient(e,g)
	ddedg = np.gradient(dedg,g)

	if k == 0:
		c = (k - kappas[0])/(kappas[-1]-kappas[0])
		plt.scatter(g,ddedg, label=r"$\kappa = 0 $", color="black", facecolor='none')
	else:
		c = (k - kappas[0])/(kappas[-2]-kappas[0])
		plt.scatter(g,ddedg, label=r"$\kappa = 1e-{} $".format(k), color=cmap(c))

	#emax = max( emax, max(ddedg) )

plt.ylim([ 1e-7, 1])
plt.yscale('log')
plt.xlim([ 1e-2, 10] )
plt.xscale('log')
plt.legend();

plt.xlabel(r"$\gamma$", fontsize=20)

plt.show()
