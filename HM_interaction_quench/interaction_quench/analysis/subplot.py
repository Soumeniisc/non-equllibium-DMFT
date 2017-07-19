import matplotlib.pyplot as plt
from scipy import *
import numpy as np


plt.figure(1,figsize=(8,9))
ax = [plt.subplot(2,1,i+1) for i in range(2)]
color_list = ["b","r","g","m","y","Purple"]
U_list = [0.5,1.0,1.5,2.0,2.5,3.0]
#U_list = [4.0,5.0,6.0,7.0]
for i,U in enumerate(U_list):
	d = loadtxt("double-occupancy_U%s"%U)
	ax[0].plot(d[:,0],np.array(d[:,1]),'--',color =color_list[i], label = r'$Uf=%s$'%U, lw=2)

ax[0].text(0.0,1.0,r'$Ui=0.0$',size=25)	
ax[0].legend(loc=8)
plt.subplot(211)
plt.ylabel("d",size=20)
plt.ylim(0.09, 0.25)
#ax[0].xlabel("t")

ax[0].set_xticklabels([])


for i,U in enumerate(U_list):
	d = loadtxt("dn_U%s"%U)
	ax[1].plot(d[:,0],abs(np.array(d[:,1])+1j*np.array(d[:,2]))**2,'--',color =color_list[i],label="U=%s"%U, lw=2)

plt.subplot(212)
plt.ylabel(r'$\Delta n$',size=20)
plt.xlabel("t",size=20)
plt.savefig("dn_d.eps", format="eps")
plt.show()
