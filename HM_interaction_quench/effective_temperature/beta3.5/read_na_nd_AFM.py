import matplotlib.pyplot as plt
from scipy import *
import numpy as np	


U_list = [0.0,0.2,0.4,0.6,0.8,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.1]
nA = []
nB = []
nAd = []
nBd = []
f = open("na_nb.dat",'w')
print >> f, "#U, na,nad,nb,nbd, 0.5*(na-nad-nb+ndd),0.5*(-na-nad+nb+ndd),(na+nad+nb+ndd)"
for i,U in enumerate(U_list):
	with open("density_%s"%U, 'r') as verFile:
    		for line in verFile:
			
		 	text = line.split()
			na = float(text[1])
			nb = float(text[2])

	with open("densityd_%s"%U, 'r') as verFile:
    		for line in verFile:
			
		 	text = line.split()
			nad = float(text[1])
			nbd = float(text[2])
	print na,nb
	print >> f,U, na,nad,nb,nbd, 0.5*(na-nad-nb+nbd),0.5*(-na-nad+nb+nbd),(na+nad+nb+nbd)
	nA.append(na)
	nB.append(nb)
	nAd.append(nad)
	nBd.append(nbd)
f.close()

data = loadtxt("na_nb.dat")
plt.plot(data[:,0],data[:,6],'-o')
plt.ylabel("dn")
plt.xlabel("U")
plt.savefig("dnvsU.eps",format="eps")
plt.show()

plt.plot(data[:,0],data[:,5],'-o')
plt.ylabel("ms")
plt.xlabel("U")
plt.savefig("msvsU.eps",format="eps")
plt.show()




