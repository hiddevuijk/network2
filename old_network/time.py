import numpy as np
import matplotlib.pyplot as plt
from sys import exit

T70 = 40*60+48
T60 = 25*60+8
T50 = (20*60 + 55 + 13*60+13)/2.
T40 = 9*60+29
T30 = 3*60+56
T20 = 49
T15 = 16.7
T10 = 3.7

T = [T10,T15, T20, T30, T40,T50,T60, T70]
N = [10,15, 20, 30,40,50,60,70]

deg = 2

fit = np.polyfit(N,T,deg)
print(fit)
def f2(x):
	return fit[0]*(x**deg)

def f(x):
	temp = 0
	for ai,a in enumerate(fit):
		n = deg-ai
		print(n)
		temp += a*(x**n)
	return temp
print(f(150)/(60*60))
x = np.linspace(0,100)
y = f(x)
y2 = f2(x)
plt.plot(x,y)
#plt.plot(x,y2, linestyle=":")

plt.scatter(N,T)
plt.show()

