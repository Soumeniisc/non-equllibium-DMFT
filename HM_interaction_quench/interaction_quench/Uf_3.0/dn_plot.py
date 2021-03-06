import matplotlib.pyplot as plt
from scipy import *
import numpy as np	

d = loadtxt("dn")
#dv = loadtxt("dn_vol")
U = 0.5
plt.plot(d[:,0],abs(np.array(d[:,1])+1j*np.array(d[:,2]))**2,'-o')
#plt.plot(dv[:,0],abs(np.array(dv[:,1])+1j*np.array(dv[:,2]))**2,'-*')
plt.ylim(0,1.0)
plt.savefig("nt_Uf%s.eps"%U,format="eps")
plt.show()
