import matplotlib.pyplot as plt
from scipy import *
import numpy as np	

t2=0.0
t=0.5
U_list = [4.0,5.0,6.0,7.0]
U_list = [0.0]#,3.0,3.5,4.0,4.5,5.0,5.5]#,2.0,2.5,3.0,3.5]
symbol=['-o','-^','-8','-h','-s','-','--','.','-*','-*','-*','-*',]
for i,U in enumerate(U_list):
	dA = loadtxt("GfA_%s"%U)
	dB = loadtxt("GfB_%s"%U)
	plt.plot(dA[:,0],np.array(dA[:,2]),symbol[i],label = r'$UA=%s$'%U)
	plt.plot(dB[:,0],np.array(dB[:,2]),symbol[i],label = r'$UB=%s$'%U)
        plt.plot(dA[:,0],np.array(dA[:,1]),symbol[i],label = r'$UA=%s$'%U)
	plt.plot(dB[:,0],np.array(dB[:,1]),symbol[i],label = r'$UB=%s$'%U)
	#plt.text(0.0,1.0,r'$Ui=0.0$',size=25)
	plt.legend()
	plt.ylabel("ImG")
	plt.xlabel("wn")
	plt.xlim(0,10)
	#plt.ylim(-3 ,0)
	plt.savefig("IMGfA_%s.eps"%U, format = "eps")
	plt.show()

