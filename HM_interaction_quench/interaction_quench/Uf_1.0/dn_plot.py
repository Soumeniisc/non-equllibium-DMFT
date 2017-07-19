import matplotlib.pyplot as plt
from scipy import *
import numpy as np	

d = loadtxt("dn_vol")
d2 = loadtxt("dn_")
U = 0.5
plt.plot(d[:,0],abs(np.array(d[:,1])+1j*np.array(d[:,2]))**2,'-o')
plt.plot(d2[:,0],abs(np.array(d2[:,1])+1j*np.array(d2[:,2]))**2,'-*')
plt.ylim(0,1.0)
plt.savefig("nt_Uf%s.eps"%U,format="eps")
plt.show()
