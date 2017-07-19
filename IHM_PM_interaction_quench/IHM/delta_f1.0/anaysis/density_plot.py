import matplotlib.pyplot as plt
from scipy import *
import numpy as np

U_i = 0.0
U_f = 1.0

color_list = ['r', 'b', 'g','m']



U_list  = [0.5,      1.0,   2.0,   3.000]
nA_list = [0.146,  0.177,  0.252,  0.333]
nB_list = [0.854,  0.825,  0.748,  0.667]
ax = plt.axes()
for i,U in enumerate(U_list):
	dA = loadtxt("densityA_U%s"%U)
	dB = loadtxt("densityB_U%s"%U)
	plt.plot(dA[:,0],dA[:,1],color=color_list[i],label=r'$U=%s$'%U)
	plt.plot(dB[:,0],dB[:,1],color=color_list[i])
	ax.arrow(4.5, nA_list[i], 0.4 , 0.0, head_width=0.03, head_length=0.1, fc=color_list[i], ec=color_list[i])
	ax.arrow(4.5, nB_list[i], 0.4 , 0.0, head_width=0.03, head_length=0.1, fc=color_list[i], ec=color_list[i])


plt.text(1.0,0.5,r'$\Delta_i = 0.0t \rightarrow \Delta_f = 1.0t,\beta=16.0$',size=20)
#plt.text(2.0,0.6,r'$U_i = 0.0t \rightarrow U_f = 2.0t$',size = 20)
plt.legend(loc = 5)
#plt.ylim(0.4,0.5)
plt.xlabel("time")
plt.ylabel("density")
plt.savefig("density.eps",format="eps")
plt.show()

ax = plt.axes()
U_list  = [0.5,      1.0,   2.0,   3.000]
dA_list = [0.018,   0.0222, 0.0314,0.0448]
dB_list = [0.73,    0.672,  0.5285,0.379]
print 
for i,U in enumerate(U_list):
	dA = loadtxt("double-occupancyA_U%s"%U)
	dB = loadtxt("double-occupancyB_U%s"%U)
	plt.plot(dA[:,0],dA[:,1],color=color_list[i],label=r'$U=%s$'%U)
	plt.plot(dB[:,0],dB[:,1],color=color_list[i])
	ax.arrow(4.5, dA_list[i], 0.4 , 0.0, head_width=0.03, head_length=0.1, fc=color_list[i], ec=color_list[i])
	ax.arrow(4.5, dB_list[i], 0.4 , 0.0, head_width=0.03, head_length=0.1, fc=color_list[i], ec=color_list[i])


plt.text(1.0,0.5,r'$\Delta_i = 0.0t \rightarrow \Delta_f = 1.0t,\beta=16.0$',size=20)
plt.legend()
plt.xlabel("time")
plt.ylabel("double-occupancy")
plt.savefig("double-occupancy.eps",format="eps")
plt.show()



