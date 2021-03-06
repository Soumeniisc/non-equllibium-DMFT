import matplotlib.pyplot as plt
from scipy import *
import numpy as np

U_i = 0.0
U_f = 1.0
delta = 1.0
string = "B_A_init_A_B_solver"
string = ""
dA = loadtxt("densityA")
dB = loadtxt("densityB")
U  = loadtxt("interaction") 

ax = plt.axes()
#plt.plot(d[:,0],d[:,1],'-o',label="HM_original")
plt.plot(dA[:,0],dA[:,1],label="n_A")
plt.plot(dB[:,0],dB[:,1],label="n_B")
plt.plot(dB[:,0],(np.array(dB[:,1])+np.array(dA[:,1]))/2.0,label="ntotal")
ax.arrow(4.5, 0.177, 0.4 , 0.0, head_width=0.03, head_length=0.1, fc='r', ec='r')
ax.arrow(4.5, 0.823, 0.4 , 0.0, head_width=0.03, head_length=0.1, fc='r', ec='r')
plt.text(2.0,0.3,r'$\Delta = 1.0t$',size=20)
plt.text(2.0,0.7,r'$U_i = 0.0t \rightarrow U_f = 1.0t$',size = 20)
plt.legend(loc = 5)
plt.xlabel("time")
plt.ylabel("density")
plt.savefig("density%s.eps"%string,format="eps")
plt.show()



dAI = loadtxt("double-occupancyA")
dBI = loadtxt("double-occupancyB")
#print d[:,1]
#plt.plot(d[:,0],d[:,1],'-o',label="HM_original")
plt.plot(dAI[:,0],dAI[:,1],'-o',label="HM_A")
plt.plot(dBI[:,0],dBI[:,1],'-o',label="HM_B")
plt.legend()
plt.xlabel("time")
plt.ylabel("double-occupancy")
plt.savefig("double-occupancy%s.eps"%string,format="eps")
plt.show()


#d = loadtxt("kinetic-energy")
dA = loadtxt("kinetic-energyA_")
dB = loadtxt("kinetic-energyB_")
#dB1 = loadtxt("kinetic-energyB1")
#print d[:,1]
#plt.plot(d[:,0],d[:,1],'-o',label="HM_original")
plt.plot(dA[:,0],dA[:,1],'-o',label="HM_A")
plt.plot(dB[:,0],dB[:,1],'-o',label="HM_B")
#plt.plot(dB1[:,0],dB1[:,1],'-o',label="HM_B1")
plt.legend()
plt.xlabel("time")
plt.ylabel("kinetic-energy")
plt.savefig("kinetic-energy%s_new.eps"%string,format="eps")
plt.show()

#d = loadtxt("interaction-energy")
plt.plot(dAI[:,0],np.array(U[:,1])*np.array(dAI[:,1]),'-o',label="HM_A")
plt.plot(dBI[:,0],np.array(U[:,1])*np.array(dBI[:,1]),'-o',label="HM_B")
plt.legend()
plt.xlabel("time")
plt.ylabel("interaction-energy")
plt.savefig("interaction-energy%s.eps"%string,format="eps")
plt.show()

dAd = loadtxt("densityA")
dBd = loadtxt("densityB")
plt.plot(dAI[:,0],delta*np.array(dAd[:,1]) + np.array(dA[:,1])+np.array(U[:,1])*np.array(dAI[:,1]),'-o',label="HM_A")
plt.plot(dBI[:,0],-delta*np.array(dBd[:,1]) + np.array(dB[:,1])+np.array(U[:,1])*np.array(dBI[:,1]),'-o',label="HM_B")
plt.legend()
plt.xlabel("time")
plt.ylabel("total-energy")
plt.savefig("total-energy%s_new.eps"%string,format="eps")
plt.show()

plt.plot(dAI[:,0],-delta*np.array(dBd[:,1]),'-o',label="site_energy")
plt.plot(dBI[:,0], np.array(dB[:,1]),'-o',label="kinetic")
plt.plot(dBI[:,0], np.array(U[:,1])*np.array(dBI[:,1]),'-o',label="potential")
plt.plot(dBI[:,0], -delta*np.array(dBd[:,1]) + np.array(dB[:,1]),'-o',label="site+kinetic")
plt.plot(dBI[:,0],-delta*np.array(dBd[:,1]) + np.array(dB[:,1])+np.array(U[:,1])*np.array(dBI[:,1]),'-o',label="total")
plt.legend()
plt.xlabel("time")
plt.ylabel("total-energy")
plt.savefig("total-energyB%s_new.eps"%string,format="eps")
plt.show()


