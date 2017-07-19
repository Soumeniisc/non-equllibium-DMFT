import matplotlib.pyplot as plt
from scipy import *
import numpy as np

U_i = 0.0
U_f = 1.0

string = "B_A_init_A_B_solver"
string = ""
d = loadtxt("density")
dA = loadtxt("densityA")
dB = loadtxt("densityB")
#print d[:,1]

ax = plt.axes()
#plt.plot(d[:,0],d[:,1],'-o',label="HM_original")
plt.plot(dA[:,0],dA[:,1],label="n_A")
plt.plot(dB[:,0],dB[:,1],label="n_B")
plt.plot(dB[:,0],(np.array(dB[:,1])+np.array(dA[:,1]))/2.0,label="ntotal")
ax.arrow(4.5, 0.252, 0.4 , 0.0, head_width=0.03, head_length=0.1, fc='r', ec='r')
ax.arrow(4.5, 0.7471, 0.4 , 0.0, head_width=0.03, head_length=0.1, fc='r', ec='r')
plt.text(2.0,0.4,r'$\Delta = 1.0t$',size=20)
plt.text(2.0,0.6,r'$U_i = 0.0t \rightarrow U_f = 2.0t$',size = 20)
plt.legend(loc = 5)
#plt.ylim(0.4,0.5)
plt.xlabel("time")
plt.ylabel("density")
plt.savefig("density%s.eps"%string,format="eps")
plt.show()



d = loadtxt("double-occupancy")
dA = loadtxt("double-occupancyA")
dB = loadtxt("double-occupancyB")
#print d[:,1]
#plt.plot(d[:,0],d[:,1],'-o',label="HM_original")
plt.plot(dA[:,0],dA[:,1],'-o',label="HM_A")
plt.plot(dB[:,0],dB[:,1],'-o',label="HM_B")
plt.legend()
plt.xlabel("time")
plt.ylabel("double-occupancy")
plt.savefig("double-occupancy%s.eps"%string,format="eps")
plt.show()


#d = loadtxt("kinetic-energy")
dA = loadtxt("kinetic-energyA")
dB = loadtxt("kinetic-energyB")
#dB1 = loadtxt("kinetic-energyB1")
#print d[:,1]
#plt.plot(d[:,0],d[:,1],'-o',label="HM_original")
plt.plot(dA[:,0],dA[:,1],'-o',label="HM_A")
plt.plot(dB[:,0],dB[:,1],'-o',label="HM_B")
#plt.plot(dB1[:,0],dB1[:,1],'-o',label="HM_B1")
plt.legend()
plt.xlabel("time")
plt.ylabel("kinetic-energy")
plt.savefig("kinetic-energy%s.eps"%string,format="eps")
plt.show()

#d = loadtxt("interaction-energy")
dA = loadtxt("interaction-energyA")
dB = loadtxt("interaction-energyB")
#print d[:,1]
#plt.plot(d[:,0],d[:,1],'-o',label="HM_original")
plt.plot(dA[:,0],dA[:,1],'-o',label="HM_A")
plt.plot(dB[:,0],dB[:,1],'-o',label="HM_B")
plt.legend()
plt.xlabel("time")
plt.ylabel("interaction-energy")
plt.savefig("interaction-energy%s.eps"%string,format="eps")
plt.show()

#d = loadtxt("total-energy")
dA = loadtxt("total-energyA")
dB = loadtxt("total-energyB")
#print d[:,1]
#plt.plot(d[:,0],d[:,1],'-o',label="HM_original")
plt.plot(dA[:,0],dA[:,1],'-o',label="HM_A")
plt.plot(dB[:,0],dB[:,1],'-o',label="HM_B")
plt.legend()
plt.xlabel("time")
plt.ylabel("total-energy")
plt.savefig("total-energy%s.eps"%string,format="eps")
plt.show()

