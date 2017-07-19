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
	copy("density_U%s_beta%s"%(U,beta), "density_old")
	#copy("densityd_%s"%U, "densityd_old")
	#replace("density_%s"%U, "density_old",stagger)
	#replaced("densityd_%s"%U, "densityd_old",stagger)
	copy("Gf_U%s_beta%s"%(U,beta), "Gf_old")
	#copy("GfBd_%s"%U, "GfBd_old")

U_list = [0.0,0.4,0.7,1.0,1.4,1.8,2.0,2.2,2.6,3.0,3.3]
#beta_list = [3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.5,4.6,4.7,4.8,5.0]
for i, U in enumerate(U_list):
	# Parameters of the nonequilibrium DMFT simulation for the single-band Hubbard model
	dos="semicircular"  # density of states
	beta=1.0       # inverse temperature
	t_max=20        # duration of time evolution
	U_i=U        # initial value of the interaction
	U_f=2.0           # final value of the interaction after quench
	N_tau=2048        # number of imaginary-time steps
	N_t=2           # number of real-time steps
	N_e=1000          # number of energy steps
	N_iter=4000       # maximum iteration of DMFT self-consistency
	tolerance=0.0001  # tolerance for DMFT convergence
	solver="IPT"        # impurity solver
	mix=0.5		 # G = (1-mix)*G_new + G except first iteration
	delta=0.0	# orbital potential

	if i==0: 
		init=   0         #0 when G0 is non-interacting GF.1 when G0 is interacting GF of previous U value
		stagger = 0.00
		#intialise(3.0,stagger)
	if i!=0: 
		init=1            #0 when G0 is non-interacting GF.1 when G0 is interacting GF of previous U value
		stagger = 0.00
		intialise(U_list[i-1],stagger)

	cmd_ex = "./eq-dmft_HM %s %f %f %f %f %i %i %i %i %f %s %i %f %f"%(dos, beta, t_max, U_i, U_f, N_tau, N_t, N_e, N_iter, tolerance, solver,init, mix, delta)
	print cmd_ex

	IPTA_info = open( "IPT_U%s_delta%s.dat"%(U_i,delta),'w')		
	subprocess.call(cmd_ex,shell=True,stdout=IPTA_info,stderr=IPTA_info)
	IPTA_info.flush()

	# coping for up spin
	copy("Gf", "Gf_U%s_beta%s"%(U_i,beta))
	#copy("Sig", "Sig_%s"%U_i)
	copy("density", "density_U%s_beta%s"%(U_i,beta))
	copy("double-occupancy", "double-occupancy_U%s_beta%s"%(U_i,beta))
	copy("kinetic-energy", "kinetic-energy_U%s_beta%s"%(U_i,beta))
