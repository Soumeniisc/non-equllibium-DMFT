import matplotlib.pyplot as plt
from scipy import *
import numpy as np	

t2=0.0
t=0.5
U_list = [0.5,1.0,1.5,2.0,2.5,3.0]
#U_list = [4.0,5.0,6.0,7.0]
for U in U_list:
	d = loadtxt("double-occupancy_U%s"%U)
	plt.plot(d[:,0],np.array(d[:,1]),'--',label = r'$Uf=%s$'%U)
	plt.text(0.0,1.0,r'$Ui=0.0$',size=25)
	plt.legend()
	plt.ylabel("double_occu(t)")
	plt.xlabel("time_t")
	#plt.xlim(-3,3)
plt.legend(loc=8)
plt.xlabel("t")
plt.ylabel(r'$d$')
#plt.savefig("double_occu_large_U.pdf")
plt.savefig("double_occu_small_U.pdf")
plt.show()

U_list = [2.0]
for U in U_list:
        d = loadtxt("double-occupancy_U%s"%U)
	dI = loadtxt("interaction-energy_U%s"%U)
        dk = loadtxt("kinetic-energy_U%s"%U)
        dt = loadtxt("total-energy_U%s"%U)
	plt.plot(d[:,0],U*np.array(d[:,1]),'-o',label = r'$interaction$')
        plt.plot(dk[:,0],np.array(dk[:,1]),'-o',label = r'$kinetic$')
	plt.plot(dt[:,0],np.array(dk[:,1])+U*np.array(d[:,1]),'-o',label = r'$total$')
        #plt.plot(dt[:,0],np.array(dt[:,1]),'-o',label = r'$total$')
	plt.text(3.0,-0.5,r'$U_i = 0.0 \rightarrow U_f=2.0$',size=20)
	plt.legend()
	#plt.ylabel("double_occu(t)")
	plt.xlabel("time_t")
	#plt.xlim(-3,3)
plt.savefig("energy_U%s.eps"%U,format="eps")
plt.show()

U_list = [5.0]
for U in U_list:
	d = loadtxt("double-occupancy_U%s"%U)
	dI = loadtxt("interaction-energy_U%s"%U)
        dk = loadtxt("kinetic-energy_U%s"%U)
        dt = loadtxt("total-energy_U%s"%U)
	plt.plot(d[:,0],U*np.array(d[:,1]),'-o',label = r'$interaction$')
        plt.plot(dk[:,0],np.array(dk[:,1]),'-o',label = r'$kinetic$')
	plt.plot(dt[:,0],np.array(dk[:,1])+U*np.array(d[:,1]),'-o',label = r'$total$')
        #plt.plot(dt[:,0],np.array(dt[:,1]),'-o',label = r'$total$')
	plt.text(3.0,0.5,r'$U_i = 0.0 \rightarrow U_f = 5.0$',size=20)
	plt.legend(loc = 2)
	#plt.ylabel("double_occu(t)")
	plt.xlabel("time_t")
	#plt.xlim(-3,3)
plt.savefig("energy_U%s.eps"%U,format="eps")
plt.show()

U_list = [1.0]
for U in U_list:
	d = loadtxt("double-occupancy_U%s"%U)
	dI = loadtxt("interaction-energy_U%s"%U)
        dk = loadtxt("kinetic-energy_U%s"%U)
        dt = loadtxt("total-energy_U%s"%U)
	plt.plot(d[:,0],U*np.array(d[:,1]),'-o',label = r'$interaction$')
        plt.plot(dk[:,0],np.array(dk[:,1]),'-o',label = r'$kinetic$')
	plt.plot(dt[:,0],np.array(dk[:,1])+U*np.array(d[:,1]),'-o',label = r'$total$')
        #plt.plot(dt[:,0],np.array(dt[:,1]),'-o',label = r'$total$')
	plt.text(3.0,0.5,r'$U_i = 0.0 \rightarrow U_f = %s$'%U,size=20)
	plt.legend(loc = 2)
	#plt.ylabel("double_occu(t)")
	plt.xlabel("time_t")
	#plt.xlim(-3,3)
plt.savefig("energy_U%s.eps"%U,format="eps")
plt.show()

U_list = [3.3]
for U in U_list:
	d = loadtxt("double-occupancy_U%s"%U)
	dI = loadtxt("interaction-energy_U%s"%U)
        dk = loadtxt("kinetic-energy_U%s"%U)
        dt = loadtxt("total-energy_U%s"%U)
	plt.plot(d[:,0],U*np.array(d[:,1]),'-o',label = r'$interaction$')
        plt.plot(dk[:,0],np.array(dk[:,1]),'-o',label = r'$kinetic$')
	plt.plot(dt[:,0],np.array(dk[:,1])+U*np.array(d[:,1]),'-o',label = r'$total$')
        #plt.plot(dt[:,0],np.array(dt[:,1]),'-o',label = r'$total$')
	plt.text(3.0,0.5,r'$U_i = 0.0 \rightarrow U_f = %s$'%U,size=20)
	plt.legend(loc = 2)
	#plt.ylabel("double_occu(t)")
	plt.xlabel("time_t")
	#plt.xlim(-3,3)
plt.savefig("energy_U%s.eps"%U,format="eps")
plt.show()

U_list = [3.0]
for U in U_list:
	d = loadtxt("double-occupancy_U%s"%U)
        dk = loadtxt("kinetic-energy_U%s"%U)
	plt.plot(d[:,0],U*np.array(d[:,1]),'-o',label = r'$interaction$')
        plt.plot(dk[:,0],np.array(dk[:,1]),'-o',label = r'$kinetic$')
	plt.plot(dt[:,0],np.array(dk[:,1])+U*np.array(d[:,1]),'-o',label = r'$total$')
        #plt.plot(dt[:,0],np.array(dt[:,1]),'-o',label = r'$total$')
	plt.text(3.0,0.5,r'$U_i = 0.0 \rightarrow U_f = %s$'%U,size=20)
	plt.legend(loc = 2)
	#plt.ylabel("double_occu(t)")
	plt.xlabel("time_t")
	#plt.xlim(-3,3)
plt.savefig("energy_U%s.eps"%U,format="eps")
plt.show()
