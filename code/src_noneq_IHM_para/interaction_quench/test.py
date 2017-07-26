# this two term are the asymptotic behaviour of host freen's function.

import numpy as np
beta = 50

sum_1 = 0.0
sum_2 = 0.0
tau = 14
for i in range(2024):

	w= (2*i + 1)*np.pi/beta
	
	sum_1 = sum_1 + (1/w)*2*np.sin(w*tau)
	sum_2 = sum_2 + (1/w**3)*2*np.sin(w*tau)

print sum_1/beta    # its is 0.5
print (sum_2/beta),0.25*tau*(beta-tau)  #0.25*tau*(parm_%beta-tau)

