import numpy as np
from NMDASynapse_step import NMDASynapse_step

def NMDASynapse(t,NT,Vpost,Mg,params):
	'''
	NMDA Synapse is affected by Mg2+ block which is dependent on voltage of
	postsynaptic neuron
	'''
	gmax = params.gmax

	dt = t[1]-t[0]
	ddt = dt * 1e3 # change to millis

	s = np.zeros_like(t)   	# state
	cond = np.zeros_like(t) # conductance

	for tt in range(1, len(t)):
	    s[tt], cond[tt] = NMDASynapse_step(
	        ddt,
	        NT[tt-1],
	        Vpost[tt-1],
	        Mg[tt-1],
	        gmax,
	        s[tt-1])

	return s, cond

class params:
	gmax = []

if __name__ == '__main__': # run self-test
	t = np.arange(0, 0.2, 1e-5)
	NT = np.ones_like(t)
	Vpost = np.zeros_like(t)
	Mg = np.ones_like(t)
	params = params()
	params.gmax = 1
	s, cond = NMDASynapse(t,NT,Vpost,Mg,params)
	print(s)
	print(cond)