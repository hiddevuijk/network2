import numpy as np
import matplotlib.pyplot as plt
from sys import exit



it_dt = np.loadtxt("dt.dat")
it = it_dt[:,0]
dt = it_dt[:,1]

plt.plot(it,dt)
plt.show()
