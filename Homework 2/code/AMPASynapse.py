import numpy as np
from AMPASynapse_step import AMPASynapse_step

def AMPASynapse(t,NT,params):
	'''
	AMPA Synapse taking Vpre as input
	'''
	gmax = params.gmax

	dt = t[1]-t[0]
	ddt = dt * 1e3 # change to millis

	s = np.zeros_like(t)	# state
	cond = np.zeros_like(t) # conductance

	for tt in range(len(t)):
	    s[tt], cond[tt] = AMPASynapse_step(	ddt, 
	    									NT[tt-1],
									    	gmax,
									    	s[tt-1])

	return s, cond

class params:
	gmax = []

if __name__ == '__main__': # run self-test
	t = np.arange(0, 0.2, 1e-5)
	NT = np.ones_like(t)
	params = params()
	params.gmax = 1
	s, cond = AMPASynapse(t,NT,params)
	print(s)
	print(cond)