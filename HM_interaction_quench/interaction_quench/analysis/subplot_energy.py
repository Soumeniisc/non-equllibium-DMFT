import matplotlib.pyplot as plt
from scipy import *
import numpy as np


plt.figure(1,figsize=(8,12))
ax = [plt.subplot(2,1,i+1) for i in range(2)]
color_list = ["b","r","g","m","y","Purple"]
U_list = [0.5,1.0,1.5,2.0,2.5,3.0]
U_list = [2.0]
plt.subplot(211)
for U in U_list:
        d = loadtxt("double-occupancy_U%s"%U)
	dI = loadtxt("interaction-energy_U%s"%U)
        dk = loadtxt("kinetic-energy_U%s"%U)
        dt = loadtxt("total-energy_U%s"%U)
	ax[0].plot(d[:,0],U*np.array(d[:,1]),'-^',label = r'$E_{pot}$')
        ax[0].plot(dk[:,0],np.array(dk[:,1]),'-H',label = r'$E_{kin}$')
	ax[0].plot(dt[:,0],np.array(dk[:,1])+U*np.array(d[:,1]),'-8',label = r'$E_{tot}$')
        #plt.plot(dt[:,0],np.array(dt[:,1]),'-o',label = r'$total$')
	plt.text(1.0,0.0,r'$U_i = 0.0 \rightarrow U_f=2.0$',size=20)
	plt.legend()
	#plt.ylabel("double_occu(t)")
	plt.xlabel("t")
	plt.xlim(0,3)

plt.ylabel("Energy",size=20)
plt.xlabel("t",size=20)
#plt.savefig("energy_U%s.eps"%U,format="eps")
#plt.show()
plt.subplot(212)
U_list = [5.0]
for U in U_list:
	d = loadtxt("double-occupancy_U%s"%U)
	dI = loadtxt("interaction-energy_U%s"%U)
        dk = loadtxt("kinetic-energy_U%s"%U)
        dt = loadtxt("total-energy_U%s"%U)
	plt.plot(d[:,0],U*np.array(d[:,1]),'-^',label = r'$E_{pot}$')
        plt.plot(dk[:,0],np.array(dk[:,1]),'-8',label = r'$E_{kin}$')
	plt.plot(dt[:,0],np.array(dk[:,1])+U*np.array(d[:,1]),'-h',label = r'$E_{tot}$')
        #plt.plot(dt[:,0],np.array(dt[:,1]),'-o',label = r'$total$')
	plt.text(1.0,0.5,r'$U_i = 0.0 \rightarrow U_f = 5.0$',size=20)
	plt.legend(loc = 2)
	#plt.ylabel("double_occu(t)")
	plt.xlabel("t",size=20)
	plt.xlim(0,3)
#plt.savefig("energy_U%s.eps"%U,format="eps")
#plt.show()


#ax[0].text(0.0,1.0,r'$Ui=0.0$',size=25)	
ax[0].legend(loc=8)
plt.subplot(211)
plt.ylabel(r"$Energy$",size=20)
#plt.ylim(0.09, 0.25)
#ax[0].xlabel("t")

ax[0].tick_params(labelsize=20)



plt.subplot(212)
plt.ylabel(r'$Energy$',size=20)
plt.xlabel("t",size=20)
plt.tick_params(labelsize=20)
plt.savefig("energy.eps", format="eps")
plt.show()
