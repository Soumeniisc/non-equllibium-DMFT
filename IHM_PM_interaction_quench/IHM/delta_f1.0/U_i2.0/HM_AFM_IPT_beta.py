'''
to run this code u need following thing to make remeber
1> check the init para. if it is 1 then excutable read GB_old and GBd_old and density_old and densityd_old of previous U. other wise its read density_old and densityd_old  and  create GB and GBd of non interacting Green's function in HM on bethe lattice
2> if init=1 then turn on intialise(U,staggere). there U is the previous U. that GB_%U, GBd_%U, density_%U, densityd_%U would be copied to GB_old, GBd_old, density_old, densityd_old
3> chose U_i, for which U_i you want have calculation.
'''

from scipy import * 
import os,sys,subprocess
import numpy as np




def replace(old_file_,new_file_,stagger): #its put intial small staggerd magnetization for given U for A sublattice
	with open(new_file_, 'w') as new_file:
		with open(old_file_, 'r') as old_file:
			for i, line in enumerate(old_file):
				text = line.split()
				print text
				if i == 0: new_file.write("%d	%f	%f"%(int(text[0]),float(text[1])+stagger,float(text[2])-stagger))
				
def replaced(old_file_,new_file_,stagger): #its put intial small staggerd magnetization for given U for B sublattice
	with open(new_file_, 'w') as new_file:
		with open(old_file_, 'r') as old_file:
			for i, line in enumerate(old_file):
				text = line.split()
				if i == 0: new_file.write("%d	%f	%f"%(int(text[0]),float(text[1])-stagger,float(text[2])+stagger))
		#old_file.close()



def copy(input_file, copy_to):
	cmd = "cp %s  %s"%(input_file, copy_to)
	print os.popen(cmd).read()

def intialise(U,stagger):
	#copy("density_%s"%U, "density_old")
	#copy("densityd_%s"%U, "densityd_old")
	replace("density_%s"%U, "density_old",stagger)
	replaced("densityd_%s"%U, "densityd_old",stagger)
	copy("GfB_%s"%U, "GfB_old")
	copy("GfBd_%s"%U, "GfBd_old")

# Parameters of the nonequilibrium DMFT simulation for the single-band Hubbard model
dos="semicircular"  # density of states
beta=16.0         # inverse temperature
t_max=5.0        # duration of time evolution
U_i=2.0        # initial value of the interaction
U_f=2.0           # final value of the interaction after quench
N_tau=1024        # number of imaginary-time steps
N_t=200           # number of real-time steps
N_e=1000          # number of energy steps
N_iter=4000       # maximum iteration of DMFT self-consistency
tolerance=0.0001  # tolerance for DMFT convergence
solver="IPT"        # impurity solver
init=0            #0 when G0 is non-interacting GF.1 when G0 is interacting GF of previous U value
mix=0.5		 # G = (1-mix)*G_new + G except first iteration
delta_i=0.0	# orbital potential
delta_f=1.0	# orbital potential


#stagger = 0.001
#intialise(1.7,stagger)

#/eq-dmft_IHM_AFM $dos $beta $t_max $U_i $U_f $N_tau $N_t $N_e $N_iter $tolerance $solver $init $mix
#./eq-dmft_IHM_AFM semicircular 100 1 3 2 2048 40 1000 400 0 IPT
cmd_ex = "./noneq-dmft_IHM2 %s %f %f %f %f %i %i %i %i %f %s %i %f %f %f"%(dos, beta, t_max, U_i, U_f, N_tau, N_t, N_e, N_iter, tolerance, solver,init, mix, delta_i,delta_f)
print cmd_ex
IPTA_info = open( "IPT_U%s_mix%s.dat"%(U_i,mix),'w')		
subprocess.call(cmd_ex,shell=True,stdout=IPTA_info,stderr=IPTA_info)
IPTA_info.flush()


copy("GfA", "GfA_%s"%U_i)
copy("GfB", "GfB_%s"%U_i)
copy("WeissA", "WeissA_%s"%U_i)
copy("WeissB", "WeissB_%s"%U_i)
copy("SigA", "SigA_%s"%U_i)
copy("SigB", "SigB_%s"%U_i)
copy("densityA", "densityA_U%s"%U_i)
copy("densityB", "densityB_U%s"%U_i)
copy("double-occupancyA", "double-occupancyA_U%s"%U_i)
copy("double-occupancyB", "double-occupancyB_U%s"%U_i)
copy("kinetic-energyA", "kinetic-energyA_U%s"%U_i)
copy("kinetic-energyB", "kinetic-energyB_U%s"%U_i)
copy("interaction-energyA", "interaction-energyA_U%s"%U_i)
copy("interaction-energyB", "interaction-energyB_U%s"%U_f)
copy("total-energyA", "total-energyA_U%s"%U_i)
copy("total-energyB", "total-energyB_U%s"%U_i)




