import matplotlib.pyplot as plt
from scipy import *
import numpy as np	

U_list = [0.5,1.0,1.5,2.0,2.5,3.0]

for U in U_list:
	d = loadtxt("dn_U%s"%U)
	plt.plot(d[:,0],abs(np.array(d[:,1])+1j*np.array(d[:,2]))**2,'--',label="U=%s"%U)

plt.ylim(0,1.0)
plt.legend(loc=8)
plt.xlabel("t")
plt.ylabel(r'$\Delta n$')
plt.savefig("nt_Uf.eps",format="eps")
plt.show()
