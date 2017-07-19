#!/bin/bash

# Parameters of the nonequilibrium DMFT simulation for the single-band Hubbard model
dos=semicircular  # density of states
beta=16.0         # inverse temperature
t_max=5.0         # duration of time evolution
U_i=0.0           # initial value of the interaction
U_f=5.0           # final value of the interaction after quench
N_tau=1024        # number of imaginary-time steps
N_t=200           # number of real-time steps
N_e=1000          # number of energy steps
N_iter=100        # maximum iteration of DMFT self-consistency
tolerance=0.0001  # tolerance for DMFT convergence
solver=IPT        # impurity solver

./a.out $dos $beta $t_max $U_i $U_f $N_tau $N_t $N_e $N_iter $tolerance $solver
