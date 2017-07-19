import matplotlib.pyplot as plt
from scipy import *
import numpy as np	


beta_list = [1.0,1.11,1.25,1.43,1.67,2.0, 2.5,3.0, 3.5, 4.0, 5.0]
U=2.0
potential = []
kinetic   = []
for beta in beta_list:
        d  = loadtxt("double-occupancy_U%s_beta%s"%(U,beta))
        dk = loadtxt("kinetic-energy_U%s_beta%s"%(U,beta))
	potential.append(U*d[1])
	kinetic.append(dk[1])

plt.plot( 1.0/np.array(beta_list), np.array(potential)+np.array(kinetic),'-o',label="total" )
plt.plot( [1./beta_list[0],1./beta_list[-1]], [-0.35,-0.35],'-o',label="U=2.0" )
#plt.plot( [1./beta_list[0],1./beta_list[-1]], [-0.59,-0.59],'-o',label="U=1.0" )
#plt.plot( [1./beta_list[0],1./beta_list[-1]], [0.05,0.05],'-o',label="U=3.0" )
plt.legend()
plt.xlabel("T")
plt.ylabel("Energy")
plt.savefig("total_energy_vs_beta_U%s.eps"%U,format="eps")
plt.show()

beta_list = [1.0,1.11,1.25,1.43,1.67,2.0, 2.5,3.0, 3.5, 4.0, 5.0,10.0]
U=1.0
potential = []
kinetic   = []
for beta in beta_list:
        d  = loadtxt("double-occupancy_U%s_beta%s"%(U,beta))
        dk = loadtxt("kinetic-energy_U%s_beta%s"%(U,beta))
	potential.append(U*d[1])
	kinetic.append(dk[1])

plt.plot( 1.0/np.array(beta_list), np.array(potential)+np.array(kinetic),'-o',label="total" )
plt.plot( [1./beta_list[0],1./beta_list[-1]], [-0.59,-0.59],'-o',label="U=1.0" )
plt.legend()
plt.xlabel("T")
plt.ylabel("Energy")
plt.savefig("total_energy_vs_beta_U%s.eps"%U,format="eps")
plt.show()

beta_list = [1.0,1.11,1.25,1.43,1.67,2.0, 2.5,3.0, 3.5, 4.0, 5.0]
U=3.3
potential = []
kinetic   = []
for beta in beta_list:
        d  = loadtxt("double-occupancy_U%s_beta%s"%(U,beta))
        dk = loadtxt("kinetic-energy_U%s_beta%s"%(U,beta))
	potential.append(U*d[1])
	kinetic.append(dk[1])

plt.plot( 1.0/np.array(beta_list), np.array(potential)+np.array(kinetic),'-o',label="total" )
plt.plot( [1./beta_list[0],1./beta_list[-1]], [0.0,0.0],'-o',label="U=%s"%U )
plt.legend()
plt.xlabel("T")
plt.ylabel("Energy")
plt.savefig("total_energy_vs_beta_U%s.eps"%U,format="eps")
plt.show()

beta_list = [1.0,1.11,1.25,1.43,1.67]
U=3.0
potential = []
kinetic   = []
for beta in beta_list:
        d  = loadtxt("double-occupancy_U%s_beta%s"%(U,beta))
        dk = loadtxt("kinetic-energy_U%s_beta%s"%(U,beta))
	potential.append(U*d[1])
	kinetic.append(dk[1])

plt.plot( 1.0/np.array(beta_list), np.array(potential)+np.array(kinetic),'-o',label="total" )
plt.plot( [1./beta_list[0],1./beta_list[-1]], [-0.085,-0.085],'-o',label="U=%s"%U )
plt.legend()
plt.xlabel("T")
plt.ylabel("Energy")
plt.savefig("total_energy_vs_beta_U%s.eps"%U,format="eps")
plt.show()


