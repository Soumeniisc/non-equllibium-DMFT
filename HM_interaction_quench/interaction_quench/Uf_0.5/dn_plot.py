import matplotlib.pyplot as plt
from scipy import *
import numpy as np	

d = loadtxt("dn")
d2 = loadtxt("dn_U0.5")
U = 0.5
plt.plot(d[:,0],abs(np.array(d[:,1])+1j*np.array(d[:,2])),'-o')
plt.plot(d2[:,0],abs(np.array(d2[:,1])+1j*np.array(d2[:,2])),'-o')
plt.ylim(0,1.2)
plt.savefig("nt_Uf%s.eps"%U,format="eps")
plt.show()
